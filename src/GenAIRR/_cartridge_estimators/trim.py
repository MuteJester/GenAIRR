"""Trim-distribution estimator sub-mixin (behavior-preserving). The
estimator method is moved verbatim from the original single-file module."""
from __future__ import annotations

from typing import Any, Dict, List

from ._common import _load_rearrangements
from ..reference_models import (
    EmpiricalDistributionSpec,
    ReferenceEmpiricalModels,
    TRIM_KEYS,
    TRIM_KEYS_VJ,
)


class _TrimEstimatorMixin:
    __slots__ = ()

    def estimate_trim_distributions(
        self,
        rearrangements: Any,
        *,
        min_count: int = 1,
        pseudocount: float = 0.0,
        replace: bool = True,
    ) -> "ReferenceCartridgeBuilder":
        """Estimate per-key trim distributions from observed
        AIRR rearrangement records.

        Writes ``EmpiricalDistributionSpec`` instances into
        ``self._reference_models.trims`` keyed by ``"V_3"`` /
        ``"D_5"`` / ``"D_3"`` / ``"J_5"`` (see
        :data:`GenAIRR.reference_models.TRIM_KEYS`). VJ cartridges
        only produce the ``V_3`` and ``J_5`` keys (see
        :data:`GenAIRR.reference_models.TRIM_KEYS_VJ`); the D-end
        keys are skipped silently because the cartridge has no D
        segment to trim.

        ``rearrangements`` accepts the same shapes as
        :meth:`estimate_allele_usage`:

        - a list of dicts (each row is one AIRR record),
        - a path-like to an AIRR TSV (parsed via ``csv.DictReader``
          with tab delimiter), or
        - an open text file handle pointing at AIRR TSV.

        AIRR columns consumed (per audit §1.3):

        - **VJ:** ``v_trim_3``, ``j_trim_5``.
        - **VDJ:** ``v_trim_3``, ``d_trim_5``, ``d_trim_3``, ``j_trim_5``.

        Columns deliberately ignored: ``v_trim_5`` / ``j_trim_3``
        (no engine pass — hard-zero in projection), and the
        observation-stage ``end_loss_5_length`` / ``end_loss_3_length``
        (separate corruption surface — see
        ``docs/primer_trim_end_loss_audit.md``).

        Per-row validation is **field-local**: a malformed or
        missing value in one trim column drops that field's
        contribution only — the row still feeds its other
        well-formed fields. Per-row structured entries land in
        ``report.rejected`` with reasons
        ``"missing_required_column"`` / ``"malformed_trim_value"``
        / ``"negative_trim_value"``, each carrying the AIRR
        column name.

        ``min_count`` (int, default 1) drops trim values whose
        observed integer count is strictly below the threshold
        before normalisation.

        ``pseudocount`` (float, default 0.0) adds a uniform prior
        to every **observed** trim value before normalisation
        (no support expansion; values never observed stay
        unobserved). Applied AFTER the ``min_count`` filter so
        the filter looks at raw observations.

        ``replace`` (default ``True``) controls idempotency. When
        ``True``, calling this method twice overwrites the
        previous specs and writes ``replaced=True`` on the new
        stage entry. When ``False`` AND a prior typed-plane
        ``trims`` is already attached to ``self._reference_models``,
        the call raises :class:`ValueError` before consuming any
        records.
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
            existing_trims = (
                self._reference_models.trims
                if self._reference_models is not None
                else None
            )
            if existing_trims:
                raise ValueError(
                    "trim distributions already attached to this cartridge; "
                    "pass replace=True to overwrite"
                )
        previously_estimated = any(
            entry.get("stage") == "estimate_trim_distributions"
            for entry in self._report.stages
        )

        records, source_label = _load_rearrangements(rearrangements)
        chain_has_d = self._chain_type.has_d
        active_keys = TRIM_KEYS if chain_has_d else TRIM_KEYS_VJ

        # Map plane key → AIRR column.
        key_to_column = {
            "V_3": "v_trim_3",
            "D_5": "d_trim_5",
            "D_3": "d_trim_3",
            "J_5": "j_trim_5",
        }
        # AIRR columns the estimator deliberately does NOT consume.
        ignored_columns = ("v_trim_5", "j_trim_3")

        counters: Dict[str, Dict[int, int]] = {k: {} for k in active_keys}
        skipped = {
            "missing_required_column": {col: 0 for col in
                                        (key_to_column[k] for k in active_keys)},
            "malformed_trim_value": {col: 0 for col in
                                     (key_to_column[k] for k in active_keys)},
            "negative_trim_value": {col: 0 for col in
                                    (key_to_column[k] for k in active_keys)},
        }
        dropped_columns: Dict[str, int] = {col: 0 for col in ignored_columns}
        rejected_entries: List[Dict[str, Any]] = []
        warnings: List[str] = []
        ignored_warned = {col: False for col in ignored_columns}
        vj_d_warned = False

        for row_idx, row in enumerate(records):
            # Surface non-zero `v_trim_5` / `j_trim_3` columns: track
            # the count but never consume them. One warning per
            # column across the dataset.
            for col in ignored_columns:
                raw = row.get(col)
                if raw not in (None, "", "0"):
                    try:
                        if int(str(raw).strip()) != 0:
                            dropped_columns[col] += 1
                            if not ignored_warned[col]:
                                warnings.append(
                                    f"{col} column non-zero in input — "
                                    f"contribution dropped (no V_5 / J_3 "
                                    f"trim pass in the engine)"
                                )
                                ignored_warned[col] = True
                    except (TypeError, ValueError):
                        # Malformed in an unused column is uninteresting.
                        pass

            # Warn once per dataset if a VJ cartridge sees populated
            # D-trim columns in the input. Same boundary as
            # `estimate_allele_usage`'s D-call ignore.
            if not chain_has_d:
                for col in ("d_trim_5", "d_trim_3"):
                    raw = row.get(col)
                    if raw not in (None, "", "0") and not vj_d_warned:
                        try:
                            if int(str(raw).strip()) != 0:
                                warnings.append(
                                    "VJ cartridge: d_trim_5 / d_trim_3 "
                                    "columns present in records; "
                                    "contribution ignored"
                                )
                                vj_d_warned = True
                                break
                        except (TypeError, ValueError):
                            pass

            # Field-local validation: each key/column independently.
            for key in active_keys:
                col = key_to_column[key]
                raw = row.get(col)
                if raw is None or str(raw).strip() == "":
                    skipped["missing_required_column"][col] += 1
                    rejected_entries.append(
                        {
                            "stage": "estimate_trim_distributions",
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
                    skipped["malformed_trim_value"][col] += 1
                    rejected_entries.append(
                        {
                            "stage": "estimate_trim_distributions",
                            "row_index": row_idx,
                            "column": col,
                            "key": key,
                            "value": raw,
                            "reason": "malformed_trim_value",
                        }
                    )
                    continue
                if value < 0:
                    skipped["negative_trim_value"][col] += 1
                    rejected_entries.append(
                        {
                            "stage": "estimate_trim_distributions",
                            "row_index": row_idx,
                            "column": col,
                            "key": key,
                            "value": value,
                            "reason": "negative_trim_value",
                        }
                    )
                    continue
                counters[key][value] = counters[key].get(value, 0) + 1

        # min_count + pseudocount + per-key normalisation.
        below_min: Dict[str, int] = {k: 0 for k in active_keys}
        inferred_pairs: Dict[str, List[tuple]] = {k: [] for k in active_keys}

        for key in active_keys:
            raw_counts = counters[key]
            # Drop values strictly below min_count.
            kept = {v: c for v, c in raw_counts.items() if c >= min_count}
            below_min[key] = len(raw_counts) - len(kept)
            # Apply pseudocount to observed values only.
            if pseudocount > 0.0:
                kept = {v: (c + pseudocount) for v, c in kept.items()}
            # Normalise to sum 1.0 if anything survived.
            total = float(sum(kept.values()))
            if total <= 0.0 or not kept:
                inferred_pairs[key] = []
                continue
            normalised = [(v, kept[v] / total) for v in sorted(kept.keys())]
            inferred_pairs[key] = normalised

        if any(below_min[k] for k in active_keys):
            warnings.append(
                f"trim values below min_count={min_count} dropped — "
                + ", ".join(f"{k}={below_min[k]}" for k in active_keys)
            )

        # Build the per-key EmpiricalDistributionSpec instances.
        new_trims: Dict[str, EmpiricalDistributionSpec] = {}
        for key, pairs in inferred_pairs.items():
            if not pairs:
                continue
            spec = EmpiricalDistributionSpec(pairs)
            spec.validate(name=f"trims[{key}]")
            new_trims[key] = spec

        # Attach to the existing reference_models (creating it if
        # absent). Other typed planes are preserved.
        if self._reference_models is None:
            self._reference_models = ReferenceEmpiricalModels()
        self._reference_models = ReferenceEmpiricalModels(
            np_lengths=self._reference_models.np_lengths,
            trims=new_trims,
            np_bases=self._reference_models.np_bases,
            p_nucleotide_lengths=self._reference_models.p_nucleotide_lengths,
            allele_usage=self._reference_models.allele_usage,
        )
        # Validate the full container under the cartridge's chain
        # type so D-on-VJ etc. raise at attach time.
        chain_label = "vdj" if chain_has_d else "vj"
        self._reference_models.validate(chain_type=chain_label)

        self._report.stages.append(
            {
                "stage": "estimate_trim_distributions",
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
                    "below_min_count": {k: below_min.get(k, 0) for k in TRIM_KEYS},
                    "dropped_columns": dropped_columns,
                },
                "warnings": warnings,
            }
        )
        self._report.rejected.extend(rejected_entries)
        return self
