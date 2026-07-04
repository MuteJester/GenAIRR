"""P-nucleotide length estimator sub-mixin (behavior-preserving). The
estimator method is moved verbatim from the original single-file module."""
from __future__ import annotations

from typing import Any, Dict, List

from ._common import _load_rearrangements
from ..reference_models import (
    EmpiricalDistributionSpec,
    P_NUCLEOTIDE_END_KEYS,
    P_NUCLEOTIDE_END_KEYS_VJ,
    ReferenceEmpiricalModels,
)


class _PNucleotideLengthEstimatorMixin:
    __slots__ = ()

    def estimate_p_nucleotide_lengths(
        self,
        rearrangements: Any,
        *,
        min_count: int = 1,
        pseudocount: float = 0.0,
        replace: bool = True,
    ) -> "ReferenceCartridgeBuilder":
        """Estimate per-end P-nucleotide length distributions
        from observed AIRR rearrangement records.

        Writes ``EmpiricalDistributionSpec`` instances into
        ``self._reference_models.p_nucleotide_lengths``
        keyed by ``"V_3"`` and ``"J_5"`` (always) plus
        ``"D_5"`` and ``"D_3"`` (VDJ only). VJ cartridges
        produce only the V_3 and J_5 keys; the typed-plane
        validator chain-type-rejects D-end keys on VJ at
        attach time.

        **Provenance warning.** This estimator requires
        AIRR-like records that ALREADY carry GenAIRR's
        P-length fields (``p_v_3_length`` /
        ``p_d_5_length`` / ``p_d_3_length`` /
        ``p_j_5_length``). It does NOT infer P-lengths
        from generic AIRR junction sequences, NP strings,
        or trim arithmetic — that inference problem is
        out of scope for v1. External AIRR tools (IgBLAST,
        MiXCR, …) do not model P-nucleotide additions, so
        their output will either omit these columns
        entirely (rejection storm) or populate them as
        zero (degenerate ``[(0, 1.0)]`` distribution).
        The estimator emits a stage-level warning per key
        when ≥ 95% of contributing rows reported zero —
        diagnostic of P-naïve input. See
        ``docs/p_nucleotide_length_estimation_design.md``
        §5 for the realistic-input enumeration.

        AIRR columns consumed:

        - ``p_v_3_length`` (always).
        - ``p_j_5_length`` (always).
        - ``p_d_5_length`` (VDJ only — VJ rows with a
          non-zero value raise a one-time warning per
          column and the contribution is skipped).
        - ``p_d_3_length`` (same).

        Columns deliberately ignored: every other AIRR
        column. v1 does NOT derive P-lengths from
        junction-arithmetic, NP strings, trim fields, or
        end-loss fields.

        Per-row validation is **field-local**: a malformed
        or missing value in one P-length column drops only
        that key's contribution. Per-row structured
        entries land in ``report.rejected`` with reasons
        ``"missing_required_column"`` /
        ``"malformed_length_value"`` /
        ``"negative_length_value"``.

        ``min_count`` (int, default 1) drops length values
        whose observed count is strictly below the
        threshold before normalisation.

        ``pseudocount`` (float, default 0.0) adds a uniform
        prior to every **observed** length value before
        normalisation (no support expansion). Applied
        AFTER the ``min_count`` filter.

        ``replace`` (default ``True``) controls idempotency.
        When ``False`` AND a prior typed-plane
        ``p_nucleotide_lengths`` is already attached,
        the call raises :class:`ValueError` before
        consuming any records.
        """
        if isinstance(min_count, bool) or not isinstance(min_count, int):
            raise ValueError(
                f"min_count must be an int, got {min_count!r}"
            )
        if min_count < 0:
            raise ValueError(
                f"min_count must be non-negative, got {min_count!r}"
            )
        if isinstance(pseudocount, bool) or not isinstance(
            pseudocount, (int, float)
        ):
            raise ValueError(
                f"pseudocount must be a non-negative number, got {pseudocount!r}"
            )
        if pseudocount < 0.0:
            raise ValueError(
                f"pseudocount must be non-negative, got {pseudocount!r}"
            )

        if not replace:
            existing = (
                self._reference_models.p_nucleotide_lengths
                if self._reference_models is not None
                else None
            )
            if existing:
                raise ValueError(
                    "p_nucleotide_lengths model already attached to "
                    "this cartridge; pass replace=True to overwrite"
                )
        previously_estimated = any(
            entry.get("stage") == "estimate_p_nucleotide_lengths"
            for entry in self._report.stages
        )

        records, source_label = _load_rearrangements(rearrangements)
        chain_has_d = self._chain_type.has_d
        active_keys = (
            P_NUCLEOTIDE_END_KEYS if chain_has_d
            else P_NUCLEOTIDE_END_KEYS_VJ
        )

        key_to_column = {
            "V_3": "p_v_3_length",
            "D_5": "p_d_5_length",
            "D_3": "p_d_3_length",
            "J_5": "p_j_5_length",
        }
        vj_ignored_columns = ("p_d_5_length", "p_d_3_length")

        counters: Dict[str, Dict[int, int]] = {k: {} for k in active_keys}
        contributing_counts: Dict[str, int] = {k: 0 for k in active_keys}
        zero_counts: Dict[str, int] = {k: 0 for k in active_keys}
        skipped = {
            "missing_required_column": {
                key_to_column[k]: 0 for k in active_keys
            },
            "malformed_length_value": {
                key_to_column[k]: 0 for k in active_keys
            },
            "negative_length_value": {
                key_to_column[k]: 0 for k in active_keys
            },
        }
        # VJ chains track nonzero `p_d_*_length` columns as
        # dropped columns with one warning per column.
        dropped_columns: Dict[str, int] = {}
        if not chain_has_d:
            for col in vj_ignored_columns:
                dropped_columns[col] = 0
        rejected_entries: List[Dict[str, Any]] = []
        warnings: List[str] = []
        vj_ignored_warned = {col: False for col in vj_ignored_columns}

        for row_idx, row in enumerate(records):
            # VJ: surface nonzero D-end columns as dropped + warn once per column.
            if not chain_has_d:
                for col in vj_ignored_columns:
                    raw = row.get(col)
                    if raw in (None, "", "0"):
                        continue
                    try:
                        value = int(str(raw).strip())
                    except (TypeError, ValueError):
                        continue
                    if value != 0:
                        dropped_columns[col] += 1
                        if not vj_ignored_warned[col]:
                            warnings.append(
                                f"VJ cartridge: {col} non-zero in "
                                f"input — contribution ignored (no "
                                f"D segment on a VJ chain)"
                            )
                            vj_ignored_warned[col] = True

            # Field-local validation per active key.
            for key in active_keys:
                col = key_to_column[key]
                raw = row.get(col)
                if raw is None or str(raw).strip() == "":
                    skipped["missing_required_column"][col] += 1
                    rejected_entries.append(
                        {
                            "stage": "estimate_p_nucleotide_lengths",
                            "row_index": row_idx,
                            "column": col,
                            "key": key,
                            "reason": "missing_required_column",
                        }
                    )
                    continue
                try:
                    value = int(str(raw).strip())
                except (TypeError, ValueError):
                    skipped["malformed_length_value"][col] += 1
                    rejected_entries.append(
                        {
                            "stage": "estimate_p_nucleotide_lengths",
                            "row_index": row_idx,
                            "column": col,
                            "key": key,
                            "value": raw,
                            "reason": "malformed_length_value",
                        }
                    )
                    continue
                if value < 0:
                    skipped["negative_length_value"][col] += 1
                    rejected_entries.append(
                        {
                            "stage": "estimate_p_nucleotide_lengths",
                            "row_index": row_idx,
                            "column": col,
                            "key": key,
                            "value": value,
                            "reason": "negative_length_value",
                        }
                    )
                    continue
                counters[key][value] = counters[key].get(value, 0) + 1
                contributing_counts[key] += 1
                if value == 0:
                    zero_counts[key] += 1

        # min_count + pseudocount + per-key normalisation.
        below_min: Dict[str, int] = {k: 0 for k in active_keys}
        inferred_pairs: Dict[str, List[tuple]] = {k: [] for k in active_keys}
        zero_fraction: Dict[str, float] = {}

        for key in active_keys:
            raw_counts = counters[key]
            kept = {v: c for v, c in raw_counts.items() if c >= min_count}
            below_min[key] = len(raw_counts) - len(kept)
            if pseudocount > 0.0:
                kept = {v: (c + pseudocount) for v, c in kept.items()}
            total = float(sum(kept.values()))
            if total <= 0.0 or not kept:
                inferred_pairs[key] = []
            else:
                inferred_pairs[key] = [
                    (v, kept[v] / total) for v in sorted(kept.keys())
                ]
            # zero_fraction tracks observed-row provenance before
            # normalisation — useful diagnostic for P-naïve input.
            n_contrib = contributing_counts[key]
            if n_contrib > 0:
                zero_fraction[key] = zero_counts[key] / n_contrib
            else:
                zero_fraction[key] = 0.0

        if any(below_min[k] for k in active_keys):
            warnings.append(
                f"P-nucleotide length values below min_count={min_count} "
                f"dropped — "
                + ", ".join(f"{k}={below_min[k]}" for k in active_keys)
            )

        # Per-key provenance auto-warning (audit §5.2).
        for key in active_keys:
            if (
                contributing_counts[key] > 0
                and zero_fraction[key] >= 0.95
            ):
                warnings.append(
                    f"p_nucleotide_lengths[{key}] is >=95% zero; "
                    f"input may be P-naive or lack P annotations"
                )

        # Build per-key EmpiricalDistributionSpec instances.
        new_p_lengths: Dict[str, EmpiricalDistributionSpec] = {}
        for key, pairs in inferred_pairs.items():
            if not pairs:
                continue
            spec = EmpiricalDistributionSpec(pairs)
            spec.validate(name=f"p_nucleotide_lengths[{key}]")
            new_p_lengths[key] = spec

        # Attach to existing reference_models (creating if absent);
        # other typed planes preserved.
        if self._reference_models is None:
            self._reference_models = ReferenceEmpiricalModels()
        self._reference_models = ReferenceEmpiricalModels(
            np_lengths=self._reference_models.np_lengths,
            trims=self._reference_models.trims,
            np_bases=self._reference_models.np_bases,
            p_nucleotide_lengths=new_p_lengths,
            allele_usage=self._reference_models.allele_usage,
        )
        chain_label = "vdj" if chain_has_d else "vj"
        self._reference_models.validate(chain_type=chain_label)

        self._report.stages.append(
            {
                "stage": "estimate_p_nucleotide_lengths",
                "inputs": {
                    "record_count": len(records),
                    "min_count": int(min_count),
                    "pseudocount": float(pseudocount),
                    "source": source_label,
                    "replaced": previously_estimated,
                },
                "inferred": {
                    "V_3": inferred_pairs.get("V_3", []),
                    "D_5": inferred_pairs.get("D_5", []),
                    "D_3": inferred_pairs.get("D_3", []),
                    "J_5": inferred_pairs.get("J_5", []),
                    "skipped": skipped,
                    "below_min_count": {
                        k: below_min.get(k, 0) for k in P_NUCLEOTIDE_END_KEYS
                    },
                    "zero_fraction": {
                        k: zero_fraction.get(k, 0.0)
                        for k in P_NUCLEOTIDE_END_KEYS
                    },
                    "dropped_columns": dropped_columns,
                },
                "warnings": warnings,
            }
        )
        self._report.rejected.extend(rejected_entries)
        return self
