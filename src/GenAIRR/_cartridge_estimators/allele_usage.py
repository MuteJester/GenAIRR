"""Allele-usage estimator sub-mixin (behavior-preserving). The estimator
method is moved verbatim from the original single-file module."""
from __future__ import annotations

from typing import Any, Dict, List

from ._common import (
    _allele_names_for_segment,
    _load_rearrangements,
    _split_tie_set,
)
from ..reference_models import (
    AlleleUsageSpec,
    ReferenceEmpiricalModels,
)


class _AlleleUsageEstimatorMixin:
    __slots__ = ()

    def estimate_allele_usage(
        self,
        rearrangements: Any,
        *,
        min_count: float = 1.0,
        ambiguous: str = "fractional",
        replace: bool = True,
    ) -> "ReferenceCartridgeBuilder":
        """Estimate per-segment allele-usage weights from observed
        AIRR rearrangement records.

        ``rearrangements`` accepts:

        - a list of dicts (each row is one AIRR record),
        - a path-like (filesystem path to an AIRR TSV — parsed via
          ``csv.DictReader`` with tab delimiter), or
        - an open text file handle pointing at AIRR TSV.

        ``ambiguous`` selects the tie-set policy (column values
        like ``"IGHV1*01,IGHV2*01"``):

        - ``"fractional"`` (default): split one row's credit
          ``1.0`` evenly across all known alleles in the tie set.
          Unknown allele names in the tie set are excluded; if
          NO names in the tie set are known to the cartridge, the
          row is recorded as ``unknown_allele`` and skipped.
        - ``"truth_first"``: credit only the first comma-separated
          allele in the tie set. Matches the existing
          :func:`GenAIRR._mcp_summary` convention.
        - ``"reject"``: drop ambiguous (multi-call) rows entirely
          and record them in ``report.rejected``.

        ``min_count`` drops alleles whose final per-segment count
        is strictly below the threshold. Pre-normalisation; the
        per-segment weights remaining after the drop are then
        renormalised to sum to ``1.0``.

        ``replace`` (default ``True``) controls idempotency. When
        ``True``, calling :meth:`estimate_allele_usage` twice
        overwrites the previous spec and writes ``replaced=True``
        on the new stage entry. When ``False``, the second call
        raises :class:`ValueError`.

        Updates ``self._reference_models.allele_usage`` so a
        downstream :meth:`build` carries the estimated spec into
        the cartridge.
        """
        if ambiguous not in ("fractional", "truth_first", "reject"):
            raise ValueError(
                f"ambiguous must be one of 'fractional' / 'truth_first' / "
                f"'reject', got {ambiguous!r}"
            )
        previously_estimated = any(
            entry.get("stage") == "estimate_allele_usage"
            for entry in self._report.stages
        )
        if previously_estimated and not replace:
            raise ValueError(
                "estimate_allele_usage already ran; pass replace=True to "
                "overwrite the previous spec"
            )

        records, source_label = _load_rearrangements(rearrangements)
        chain_has_d = self._chain_type.has_d
        v_pool = _allele_names_for_segment(self._v_alleles)
        d_pool = _allele_names_for_segment(self._d_alleles)
        j_pool = _allele_names_for_segment(self._j_alleles)

        v_counts: Dict[str, float] = {}
        d_counts: Dict[str, float] = {}
        j_counts: Dict[str, float] = {}
        skipped = {
            "missing_required_column": 0,
            "unknown_allele": {"V": 0, "D": 0, "J": 0},
            "missing_d_call_on_vdj": 0,
            "ambiguous_rejected": 0,
        }
        rejected_entries: List[Dict[str, Any]] = []
        warnings: List[str] = []
        d_on_vj_warned = False

        for row_idx, row in enumerate(records):
            v_raw = (row.get("v_call") or "").strip()
            j_raw = (row.get("j_call") or "").strip()
            d_raw = (row.get("d_call") or "").strip()

            # Required column check.
            if not v_raw or not j_raw:
                skipped["missing_required_column"] += 1
                rejected_entries.append(
                    {
                        "stage": "estimate_allele_usage",
                        "row_index": row_idx,
                        "reason": "missing_required_column",
                    }
                )
                continue

            # VDJ requires d_call; VJ ignores any d_call.
            if chain_has_d and not d_raw:
                skipped["missing_d_call_on_vdj"] += 1
                rejected_entries.append(
                    {
                        "stage": "estimate_allele_usage",
                        "row_index": row_idx,
                        "reason": "missing_d_call_on_vdj",
                    }
                )
                continue
            if not chain_has_d and d_raw:
                if not d_on_vj_warned:
                    warnings.append(
                        "d_call column present on a VJ cartridge — D "
                        "contribution ignored for every row"
                    )
                    d_on_vj_warned = True
                d_raw = ""  # ignore D contribution silently after warning

            v_tie = _split_tie_set(v_raw)
            j_tie = _split_tie_set(j_raw)
            d_tie = _split_tie_set(d_raw) if d_raw else []

            if ambiguous == "reject":
                if len(v_tie) > 1 or len(j_tie) > 1 or (chain_has_d and len(d_tie) > 1):
                    skipped["ambiguous_rejected"] += 1
                    rejected_entries.append(
                        {
                            "stage": "estimate_allele_usage",
                            "row_index": row_idx,
                            "reason": "ambiguous_rejected",
                        }
                    )
                    continue

            # truth_first collapses the tie set to its first entry.
            if ambiguous == "truth_first":
                v_tie = v_tie[:1]
                j_tie = j_tie[:1]
                d_tie = d_tie[:1] if d_tie else []

            # Resolve known alleles in each tie set.
            v_known = [n for n in v_tie if n in v_pool]
            j_known = [n for n in j_tie if n in j_pool]
            d_known = (
                [n for n in d_tie if n in d_pool] if d_tie else []
            )

            # Unknown-allele bookkeeping. We record one rejection
            # entry per unknown allele per segment per row so the
            # report names the actual unknown names.
            for n in v_tie:
                if n not in v_pool:
                    skipped["unknown_allele"]["V"] += 1
                    rejected_entries.append(
                        {
                            "stage": "estimate_allele_usage",
                            "row_index": row_idx,
                            "segment": "V",
                            "allele_name": n,
                            "reason": "unknown_allele",
                        }
                    )
            for n in j_tie:
                if n not in j_pool:
                    skipped["unknown_allele"]["J"] += 1
                    rejected_entries.append(
                        {
                            "stage": "estimate_allele_usage",
                            "row_index": row_idx,
                            "segment": "J",
                            "allele_name": n,
                            "reason": "unknown_allele",
                        }
                    )
            if d_tie:
                for n in d_tie:
                    if n not in d_pool:
                        skipped["unknown_allele"]["D"] += 1
                        rejected_entries.append(
                            {
                                "stage": "estimate_allele_usage",
                                "row_index": row_idx,
                                "segment": "D",
                                "allele_name": n,
                                "reason": "unknown_allele",
                            }
                        )

            # Fractional credit (the default and the truth_first
            # collapsed path; both use the same accumulation
            # logic now that truth_first reduced the tie set).
            if v_known:
                share = 1.0 / len(v_known)
                for n in v_known:
                    v_counts[n] = v_counts.get(n, 0.0) + share
            if j_known:
                share = 1.0 / len(j_known)
                for n in j_known:
                    j_counts[n] = j_counts.get(n, 0.0) + share
            if d_known:
                share = 1.0 / len(d_known)
                for n in d_known:
                    d_counts[n] = d_counts.get(n, 0.0) + share

        # min_count filter + normalisation per segment.
        below_min: Dict[str, int] = {"V": 0, "D": 0, "J": 0}

        def _filter_and_normalise(
            counts: Dict[str, float], segment_label: str
        ) -> Dict[str, float]:
            kept = {n: c for n, c in counts.items() if c >= min_count}
            dropped = len(counts) - len(kept)
            below_min[segment_label] = dropped
            total = sum(kept.values())
            if total <= 0.0:
                return {}
            return {n: w / total for n, w in kept.items()}

        v_weights = _filter_and_normalise(v_counts, "V")
        d_weights = _filter_and_normalise(d_counts, "D")
        j_weights = _filter_and_normalise(j_counts, "J")

        if below_min["V"] or below_min["D"] or below_min["J"]:
            warnings.append(
                f"alleles below min_count={min_count} dropped — "
                f"V={below_min['V']}, D={below_min['D']}, J={below_min['J']}"
            )

        spec = AlleleUsageSpec(v=v_weights, d=d_weights, j=j_weights)
        chain_label = "vdj" if chain_has_d else "vj"
        spec.validate(chain_type=chain_label, name="allele_usage")

        # Attach to the existing reference_models (creating it if
        # absent). The cartridge built downstream carries the
        # spec automatically.
        if self._reference_models is None:
            self._reference_models = ReferenceEmpiricalModels()
        self._reference_models = ReferenceEmpiricalModels(
            np_lengths=self._reference_models.np_lengths,
            trims=self._reference_models.trims,
            np_bases=self._reference_models.np_bases,
            p_nucleotide_lengths=self._reference_models.p_nucleotide_lengths,
            allele_usage=spec,
        )

        self._report.stages.append(
            {
                "stage": "estimate_allele_usage",
                "inputs": {
                    "record_count": len(records),
                    "ambiguous": ambiguous,
                    "min_count": float(min_count),
                    "source": source_label,
                    "replaced": previously_estimated,
                },
                "inferred": {
                    "V": v_weights,
                    "D": d_weights,
                    "J": j_weights,
                    "skipped": skipped,
                    "below_min_count": below_min,
                },
                "warnings": warnings,
            }
        )
        self._report.rejected.extend(rejected_entries)
        return self
