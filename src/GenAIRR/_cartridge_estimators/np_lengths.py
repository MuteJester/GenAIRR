"""NP-length distribution estimator sub-mixin (behavior-preserving). The
estimator method is moved verbatim from the original single-file module."""
from __future__ import annotations

from typing import Any, Dict, List

from ._common import _load_rearrangements
from ..reference_models import (
    EmpiricalDistributionSpec,
    NP_KEYS,
    ReferenceEmpiricalModels,
)


class _NpLengthEstimatorMixin:
    __slots__ = ()

    def estimate_np_length_distributions(
        self,
        rearrangements: Any,
        *,
        min_count: int = 1,
        pseudocount: float = 0.0,
        replace: bool = True,
    ) -> "ReferenceCartridgeBuilder":
        """Estimate per-key NP length distributions from
        observed AIRR rearrangement records.

        Writes ``EmpiricalDistributionSpec`` instances into
        ``self._reference_models.np_lengths`` keyed by
        ``"NP1"`` and (on VDJ cartridges) ``"NP2"`` (see
        :data:`GenAIRR.reference_models.NP_KEYS`). VJ
        cartridges produce only the ``NP1`` key; the
        cartridge has no NP2 region. The boundary is
        enforced at estimation time because the typed-plane
        validator does NOT chain-type-reject ``NP2`` on
        VJ at attach time (see audit §2.2 of
        ``docs/np_length_estimation_design.md``).

        AIRR columns consumed (per audit §1.2):

        - ``np1_length`` (always).
        - ``np2_length`` (VDJ only — VJ rows with a
          non-zero ``np2_length`` raise a one-time warning
          and the contribution is skipped).

        Columns deliberately ignored: ``np1`` / ``np2``
        (sequence-derived length is sensitive to
        post-claim reabsorption — see audit §5.3),
        ``p_v_3_length`` / ``p_d_5_length`` /
        ``p_d_3_length`` / ``p_j_5_length`` (separate
        P-nucleotide biology — see
        ``docs/p_nucleotide_design.md``),
        and ``junction_length`` (aggregate arithmetic too
        fragile across simulators — audit §5.2).

        Per-row validation is **field-local**: a malformed
        or missing value in one NP column drops only that
        key's contribution; the row's other column still
        feeds. Per-row structured entries land in
        ``report.rejected`` with reasons
        ``"missing_required_column"`` /
        ``"malformed_length_value"`` /
        ``"negative_length_value"``, each carrying the
        AIRR column name.

        ``min_count`` (int, default 1) drops length values
        whose observed integer count is strictly below
        the threshold before normalisation.

        ``pseudocount`` (float, default 0.0) adds a uniform
        prior to every **observed** length value before
        normalisation (no support expansion). Applied
        AFTER the ``min_count`` filter.

        ``replace`` (default ``True``) controls idempotency.
        When ``False`` AND a prior typed-plane
        ``np_lengths`` is attached to
        ``self._reference_models``, the call raises
        :class:`ValueError` before consuming any records.
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
                self._reference_models.np_lengths
                if self._reference_models is not None
                else None
            )
            if existing:
                raise ValueError(
                    "np_length distributions already attached to this "
                    "cartridge; pass replace=True to overwrite"
                )
        previously_estimated = any(
            entry.get("stage") == "estimate_np_length_distributions"
            for entry in self._report.stages
        )

        records, source_label = _load_rearrangements(rearrangements)
        chain_has_d = self._chain_type.has_d
        active_keys = NP_KEYS if chain_has_d else ("NP1",)

        key_to_column = {"NP1": "np1_length", "NP2": "np2_length"}

        counters: Dict[str, Dict[int, int]] = {k: {} for k in active_keys}
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
        # VJ chains track `np2_length` contributions as a dropped
        # column with a single one-time warning.
        dropped_columns: Dict[str, int] = {}
        if not chain_has_d:
            dropped_columns["np2_length"] = 0
        rejected_entries: List[Dict[str, Any]] = []
        warnings: List[str] = []
        vj_np2_warned = False

        for row_idx, row in enumerate(records):
            # On VJ, surface non-zero `np2_length` as a dropped
            # column with one warning across the dataset.
            if not chain_has_d:
                raw_np2 = row.get("np2_length")
                if raw_np2 not in (None, "", "0"):
                    try:
                        if int(str(raw_np2).strip()) != 0:
                            dropped_columns["np2_length"] += 1
                            if not vj_np2_warned:
                                warnings.append(
                                    "VJ cartridge: np2_length column "
                                    "present in records; contribution "
                                    "ignored (no NP2 region on a VJ "
                                    "chain)"
                                )
                                vj_np2_warned = True
                    except (TypeError, ValueError):
                        pass

            # Field-local validation per active key.
            for key in active_keys:
                col = key_to_column[key]
                raw = row.get(col)
                if raw is None or str(raw).strip() == "":
                    skipped["missing_required_column"][col] += 1
                    rejected_entries.append(
                        {
                            "stage": "estimate_np_length_distributions",
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
                            "stage": "estimate_np_length_distributions",
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
                            "stage": "estimate_np_length_distributions",
                            "row_index": row_idx,
                            "column": col,
                            "key": key,
                            "value": value,
                            "reason": "negative_length_value",
                        }
                    )
                    continue
                counters[key][value] = counters[key].get(value, 0) + 1

        # min_count + pseudocount + per-key normalisation.
        below_min: Dict[str, int] = {k: 0 for k in active_keys}
        inferred_pairs: Dict[str, List[tuple]] = {k: [] for k in active_keys}

        for key in active_keys:
            raw_counts = counters[key]
            kept = {v: c for v, c in raw_counts.items() if c >= min_count}
            below_min[key] = len(raw_counts) - len(kept)
            if pseudocount > 0.0:
                kept = {v: (c + pseudocount) for v, c in kept.items()}
            total = float(sum(kept.values()))
            if total <= 0.0 or not kept:
                inferred_pairs[key] = []
                continue
            normalised = [(v, kept[v] / total) for v in sorted(kept.keys())]
            inferred_pairs[key] = normalised

        if any(below_min[k] for k in active_keys):
            warnings.append(
                f"NP-length values below min_count={min_count} dropped — "
                + ", ".join(f"{k}={below_min[k]}" for k in active_keys)
            )

        # Build per-key EmpiricalDistributionSpec instances.
        new_np_lengths: Dict[str, EmpiricalDistributionSpec] = {}
        for key, pairs in inferred_pairs.items():
            if not pairs:
                continue
            spec = EmpiricalDistributionSpec(pairs)
            spec.validate(name=f"np_lengths[{key}]")
            new_np_lengths[key] = spec

        # Attach to existing reference_models (creating if absent).
        # Other typed planes are preserved.
        if self._reference_models is None:
            self._reference_models = ReferenceEmpiricalModels()
        self._reference_models = ReferenceEmpiricalModels(
            np_lengths=new_np_lengths,
            trims=self._reference_models.trims,
            np_bases=self._reference_models.np_bases,
            p_nucleotide_lengths=self._reference_models.p_nucleotide_lengths,
            allele_usage=self._reference_models.allele_usage,
        )
        chain_label = "vdj" if chain_has_d else "vj"
        self._reference_models.validate(chain_type=chain_label)

        self._report.stages.append(
            {
                "stage": "estimate_np_length_distributions",
                "inputs": {
                    "record_count": len(records),
                    "min_count": int(min_count),
                    "pseudocount": float(pseudocount),
                    "source": source_label,
                    "replaced": previously_estimated,
                },
                "inferred": {
                    "NP1": inferred_pairs.get("NP1", []),
                    "NP2": inferred_pairs.get("NP2", []),
                    "skipped": skipped,
                    "below_min_count": {k: below_min.get(k, 0) for k in NP_KEYS},
                    "dropped_columns": dropped_columns,
                },
                "warnings": warnings,
            }
        )
        self._report.rejected.extend(rejected_entries)
        return self
