"""NP base-model estimator sub-mixin (behavior-preserving). The estimator
method is moved verbatim from the original single-file module."""
from __future__ import annotations

from typing import Any, Dict, List

from ._common import _NP_CANONICAL_BASES, _load_rearrangements
from ..reference_models import (
    NP_KEYS,
    NpBaseModelSpec,
    ReferenceEmpiricalModels,
)


class _NpBaseModelEstimatorMixin:
    __slots__ = ()

    def estimate_np_base_model(
        self,
        rearrangements: Any,
        *,
        kind: str = "markov",
        min_count: int = 1,
        pseudocount: float = 0.0,
        replace: bool = True,
    ) -> "ReferenceCartridgeBuilder":
        """Estimate per-key NP base sampling model from
        observed AIRR rearrangement records.

        Writes ``NpBaseModelSpec`` instances into
        ``self._reference_models.np_bases`` keyed by
        ``"NP1"`` (always) and (on VDJ cartridges)
        ``"NP2"``. VJ cartridges produce only the ``NP1``
        key.

        ``kind`` selects the model family:

        - ``"empirical_first_base"``: estimates a single
          categorical over A/C/G/T from the **full base
          composition** of every observed NP string (every
          position, not just position 0). The engine
          samples every NP position independently from
          this distribution — the name is preserved for
          API stability with existing cartridges; the
          biologically correct estimate is the full-base
          composition.
        - ``"markov"`` (default): estimates a first-base
          row from position 0 of each NP string plus a
          4×4 transition matrix from every observed
          (prev, next) pair.

        AIRR columns consumed (per audit §1.2):

        - ``np1`` (always).
        - ``np2`` (VDJ only — VJ rows with a non-empty
          ``np2`` raise a one-time warning and the
          contribution is skipped).

        Columns deliberately ignored: ``junction`` (audit
        §1.2), ``p_v_3_length`` / ``p_d_5_length`` /
        ``p_d_3_length`` / ``p_j_5_length`` (separate
        P-nucleotide biology), ``np1_length`` /
        ``np2_length`` (length-only — owned by
        :meth:`estimate_np_length_distributions`).

        Per-row validation is **field-local**: a malformed
        or missing value in one NP column drops only that
        key's contribution. Per-row structured entries
        land in ``report.rejected`` with reasons
        ``"missing_required_column"`` (empty / missing
        string) or ``"noncanonical_base"`` (any character
        outside ``{A,C,G,T}`` after uppercasing).

        ``min_count`` (int, default 1) drops first-base
        categories whose observed count is strictly below
        the threshold before normalisation. For ``markov``,
        ``min_count`` applies to the **first-base row
        only**: dropping transition cells could leave a
        from-base row with no positive weights, which the
        :class:`NpBaseModelSpec` validator rejects. v1
        keeps transition rows intact and surfaces the
        first-base drops in ``below_min_count.first_base``.

        ``pseudocount`` (float, default 0.0) adds a uniform
        prior:

        - ``empirical_first_base``: added to every A/C/G/T
          base category before normalisation.
        - ``markov``: added to every A/C/G/T first-base
          category AND to every cell of the 4×4 transition
          matrix.

        ``replace`` (default ``True``) controls idempotency.
        When ``False`` AND a prior typed-plane
        ``np_bases`` is already attached, the call raises
        :class:`ValueError` before consuming any records.
        """
        if kind not in ("empirical_first_base", "markov"):
            raise ValueError(
                f"kind must be one of 'empirical_first_base' / "
                f"'markov', got {kind!r}"
            )
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
                self._reference_models.np_bases
                if self._reference_models is not None
                else None
            )
            if existing:
                raise ValueError(
                    "np_bases model already attached to this cartridge; "
                    "pass replace=True to overwrite"
                )
        previously_estimated = any(
            entry.get("stage") == "estimate_np_base_model"
            for entry in self._report.stages
        )

        records, source_label = _load_rearrangements(rearrangements)
        chain_has_d = self._chain_type.has_d
        active_keys = ("NP1", "NP2") if chain_has_d else ("NP1",)

        key_to_column = {"NP1": "np1", "NP2": "np2"}
        bases = _NP_CANONICAL_BASES  # ("A","C","G","T")

        # Per-key tallies. first_base[key] is a dict over A/C/G/T;
        # transitions[key] is a 4×4 dict-of-dict over A/C/G/T from→to.
        first_base_counts: Dict[str, Dict[str, int]] = {
            key: {b: 0 for b in bases} for key in active_keys
        }
        transition_counts: Dict[str, Dict[str, Dict[str, int]]] = {
            key: {b: {t: 0 for t in bases} for b in bases}
            for key in active_keys
        }
        skipped = {
            "missing_required_column": {
                key_to_column[k]: 0 for k in active_keys
            },
            "noncanonical_base": {
                key_to_column[k]: 0 for k in active_keys
            },
        }
        # VJ chains track non-empty `np2` strings as a dropped
        # column with a single one-time warning.
        dropped_columns: Dict[str, int] = {}
        if not chain_has_d:
            dropped_columns["np2"] = 0
        rejected_entries: List[Dict[str, Any]] = []
        warnings: List[str] = []
        vj_np2_warned = False

        for row_idx, row in enumerate(records):
            # Surface non-empty `np2` on VJ as a dropped column
            # with one warning across the dataset.
            if not chain_has_d:
                raw_np2 = row.get("np2")
                if isinstance(raw_np2, str) and raw_np2.strip():
                    dropped_columns["np2"] += 1
                    if not vj_np2_warned:
                        warnings.append(
                            "VJ cartridge: np2 column non-empty in "
                            "records; contribution ignored (no NP2 "
                            "region on a VJ chain)"
                        )
                        vj_np2_warned = True

            # Field-local validation per active key.
            for key in active_keys:
                col = key_to_column[key]
                raw = row.get(col)
                if raw is None or not isinstance(raw, str) or not raw.strip():
                    skipped["missing_required_column"][col] += 1
                    rejected_entries.append(
                        {
                            "stage": "estimate_np_base_model",
                            "row_index": row_idx,
                            "column": col,
                            "key": key,
                            "reason": "missing_required_column",
                        }
                    )
                    continue
                upper = raw.strip().upper()
                non_canonical = sorted(set(upper) - set(bases))
                if non_canonical:
                    skipped["noncanonical_base"][col] += 1
                    rejected_entries.append(
                        {
                            "stage": "estimate_np_base_model",
                            "row_index": row_idx,
                            "column": col,
                            "key": key,
                            "unknown_chars": non_canonical,
                            "reason": "noncanonical_base",
                        }
                    )
                    continue
                # Tally bases.
                #   - For `empirical_first_base` the engine samples
                #     every NP position independently from the same
                #     distribution, so the biologically correct
                #     estimate is the FULL base composition (every
                #     position).
                #   - For `markov` the first-base row models position
                #     0 only (the seed of the chain); position 1+
                #     are modelled by transitions.
                if kind == "empirical_first_base":
                    for b in upper:
                        first_base_counts[key][b] += 1
                else:
                    first_base_counts[key][upper[0]] += 1
                for i in range(len(upper) - 1):
                    prev_b = upper[i]
                    next_b = upper[i + 1]
                    transition_counts[key][prev_b][next_b] += 1

        # min_count + pseudocount + per-kind spec construction.
        below_min: Dict[str, Dict[str, int]] = {
            key: {"first_base": 0} for key in active_keys
        }
        inferred: Dict[str, Dict[str, Any]] = {
            key: {} for key in active_keys
        }

        new_np_bases: Dict[str, NpBaseModelSpec] = {}
        from_base_errors: List[str] = []

        for key in active_keys:
            fb_counts = first_base_counts[key]
            # Apply min_count filter to first-base row (drop bases
            # whose count is strictly below threshold).
            fb_kept = {b: c for b, c in fb_counts.items() if c >= min_count}
            below_min[key]["first_base"] = len(fb_counts) - len(fb_kept)
            # Apply pseudocount AFTER filter — observed-only support.
            if pseudocount > 0.0:
                for b in bases:
                    if b in fb_kept:
                        fb_kept[b] = fb_kept[b] + pseudocount
                    else:
                        fb_kept[b] = pseudocount
            fb_total = float(sum(fb_kept.values()))
            if fb_total <= 0.0:
                # Nothing observed AND no pseudocount → skip this key.
                inferred[key] = {"first_base": None, "transitions": None}
                continue
            first_base_weights = {
                b: fb_kept[b] / fb_total
                for b in sorted(fb_kept.keys())
                if fb_kept[b] > 0
            }

            spec_payload: Dict[str, Any] = {
                "first_base": first_base_weights,
            }
            if kind == "markov":
                # Apply pseudocount to every cell of every transition row.
                tr = transition_counts[key]
                tr_weights: Dict[str, Dict[str, float]] = {}
                for from_b in bases:
                    row_counts = {t: float(tr[from_b][t]) for t in bases}
                    if pseudocount > 0.0:
                        for t in bases:
                            row_counts[t] += pseudocount
                    row_total = sum(row_counts.values())
                    if row_total <= 0.0:
                        # Pseudocount=0 and from-base never observed: surface a
                        # tagged error so the caller knows to pass pseudocount.
                        from_base_errors.append(
                            f"{key}: transition row for from-base "
                            f"{from_b!r} has no observed transitions "
                            f"and pseudocount=0; pass pseudocount > 0 "
                            f"or supply more data"
                        )
                        continue
                    tr_weights[from_b] = {
                        t: row_counts[t] / row_total
                        for t in bases
                        if row_counts[t] > 0
                    }
                if from_base_errors:
                    # Bail before constructing a partial spec.
                    raise ValueError("; ".join(from_base_errors))
                spec_payload["transitions"] = tr_weights
            else:
                spec_payload["transitions"] = None

            # Construct + validate the spec.
            spec = NpBaseModelSpec(
                kind=kind,
                first_base=spec_payload["first_base"],
                transitions=spec_payload["transitions"],
            )
            spec.validate(name=f"np_bases[{key}]")
            new_np_bases[key] = spec
            inferred[key] = {
                "first_base": dict(spec_payload["first_base"]),
                "transitions": (
                    {k: dict(v) for k, v in spec_payload["transitions"].items()}
                    if spec_payload["transitions"] is not None
                    else None
                ),
            }

        if any(below_min[k]["first_base"] for k in active_keys):
            warnings.append(
                f"first-base values below min_count={min_count} dropped — "
                + ", ".join(
                    f"{k}={below_min[k]['first_base']}"
                    for k in active_keys
                )
            )

        # Attach to existing reference_models (creating if absent);
        # other typed planes preserved.
        if self._reference_models is None:
            self._reference_models = ReferenceEmpiricalModels()
        self._reference_models = ReferenceEmpiricalModels(
            np_lengths=self._reference_models.np_lengths,
            trims=self._reference_models.trims,
            np_bases=new_np_bases,
            p_nucleotide_lengths=self._reference_models.p_nucleotide_lengths,
            allele_usage=self._reference_models.allele_usage,
        )
        chain_label = "vdj" if chain_has_d else "vj"
        self._reference_models.validate(chain_type=chain_label)

        # NP2 entry surfaces in inferred even on VJ (as empty dict)
        # for shape stability — same discipline as the NP-length
        # estimator's empty `D_5` / `D_3` entries.
        np1_inferred = inferred.get("NP1", {"first_base": None, "transitions": None})
        np2_inferred = inferred.get("NP2", {"first_base": None, "transitions": None})

        self._report.stages.append(
            {
                "stage": "estimate_np_base_model",
                "inputs": {
                    "record_count": len(records),
                    "kind": kind,
                    "min_count": int(min_count),
                    "pseudocount": float(pseudocount),
                    "source": source_label,
                    "replaced": previously_estimated,
                },
                "inferred": {
                    "NP1": np1_inferred,
                    "NP2": np2_inferred,
                    "skipped": skipped,
                    "below_min_count": {
                        k: below_min.get(k, {"first_base": 0})
                        for k in NP_KEYS
                    },
                    "dropped_columns": dropped_columns,
                },
                "warnings": warnings,
            }
        )
        self._report.rejected.extend(rejected_entries)
        return self
