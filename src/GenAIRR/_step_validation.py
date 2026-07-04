"""Step-parameter validation helpers for :class:`GenAIRR.experiment.Experiment`.

Extracted verbatim from ``experiment.py`` (behavior-preserving) — pure
leaf validation of ``segment_rates`` / ``v_subregion_rates`` kwargs and the
clonal descendant-phase step classifier, with no dependency on the
``Experiment`` class itself.
"""
from typing import Dict, List, Optional, Tuple


# ──────────────────────────────────────────────────────────────────
# Per-segment SHM rate validation (slice: segment_rates kwarg on
# Experiment.mutate). See ``docs/shm_segment_rate_design.md``.
# ──────────────────────────────────────────────────────────────────

# Canonical bucket order — matches the Rust ``SegmentRateWeights``
# struct field order so positional plumbing through PyO3 stays
# straightforward.
_SEGMENT_RATE_BUCKETS: Tuple[str, ...] = ("V", "D", "J", "NP")

# Default flat-substrate rate vector (1.0 for every bucket). The
# pipeline-IR ``_MutateStep`` carries this verbatim when the user
# omits ``segment_rates``; the Rust passes detect the flat-default
# case and take the existing (pre-slice) fast path so legacy
# pipelines stay byte-identical.
_DEFAULT_SEGMENT_RATES: Tuple[float, float, float, float] = (1.0, 1.0, 1.0, 1.0)


def _validate_segment_rates(
    segment_rates: Optional[Dict[str, float]],
) -> Tuple[float, float, float, float]:
    """Validate the user's ``segment_rates`` dict and return the
    normalised ``(v, d, j, np)`` tuple. Pure helper — no DSL state.

    Validation:

    - ``None`` (or omitted) → flat default ``(1.0, 1.0, 1.0, 1.0)``.
    - Keys must be a subset of ``{"V", "D", "J", "NP"}`` (case-
      sensitive — matches the DSL spec). Other keys raise
      ``ValueError`` naming the offending key.
    - Values must be ``int`` / ``float`` (not bool), finite, and
      ``>= 0``. Negative / NaN / inf raise ``ValueError``.
    - Sparse: omitted keys default to ``1.0``.
    - At least one effective rate must be strictly positive — an
      all-zero (or all-omitted-then-explicitly-zero) configuration
      would make the SHM pass a deterministic no-op, which is
      almost certainly a builder bug. Reject with ``ValueError``.
    """
    if segment_rates is None:
        return _DEFAULT_SEGMENT_RATES

    if not isinstance(segment_rates, dict):
        raise TypeError(
            f"segment_rates must be a dict or None, got "
            f"{type(segment_rates).__name__}"
        )

    # Reject unknown keys first so a typo surfaces with a clear
    # message instead of silently defaulting.
    unknown = sorted(set(segment_rates.keys()) - set(_SEGMENT_RATE_BUCKETS))
    if unknown:
        raise ValueError(
            f"segment_rates: unknown segment key(s) {unknown!r}. "
            f"Allowed: {list(_SEGMENT_RATE_BUCKETS)!r}. "
            "Keys are case-sensitive; 'V' / 'D' / 'J' for the V/D/J "
            "segments and 'NP' for both Np1 and Np2."
        )

    out_list: List[float] = []
    for bucket in _SEGMENT_RATE_BUCKETS:
        if bucket not in segment_rates:
            out_list.append(1.0)
            continue
        value = segment_rates[bucket]
        if isinstance(value, bool) or not isinstance(value, (int, float)):
            raise TypeError(
                f"segment_rates[{bucket!r}]: value must be a finite "
                f"non-negative number, got {type(value).__name__}"
            )
        value_f = float(value)
        # NaN check FIRST — every comparison with NaN is False, so
        # a `value_f < 0` test wouldn't catch it; that's the same
        # rationale as the ``invert_d`` / ``rate`` NaN handling.
        if value_f != value_f:
            raise ValueError(
                f"segment_rates[{bucket!r}]: value must be a finite "
                "number, got NaN"
            )
        if value_f < 0.0 or value_f == float("inf") or value_f == float("-inf"):
            raise ValueError(
                f"segment_rates[{bucket!r}]: value must be a finite "
                f"non-negative number, got {value_f}"
            )
        out_list.append(value_f)

    if sum(out_list) <= 0.0:
        raise ValueError(
            "segment_rates: at least one bucket must have a positive "
            "rate. The supplied configuration zeroes out every "
            "biological segment, which would make the SHM pass a "
            "deterministic no-op."
        )

    return (out_list[0], out_list[1], out_list[2], out_list[3])


# ──────────────────────────────────────────────────────────────────
# Per-V-subregion SHM rate validation (Slice B —
# `v_subregion_rates` kwarg on Experiment.mutate). See
# ``docs/v_subregion_shm_rate_design.md``.
# ──────────────────────────────────────────────────────────────────

# Canonical label order — matches the Rust ``VSubregionRateWeights``
# struct field order (FWR1 / CDR1 / FWR2 / CDR2 / FWR3) so
# positional plumbing through PyO3 stays straightforward.
_V_SUBREGION_RATE_LABELS: Tuple[str, ...] = ("FWR1", "CDR1", "FWR2", "CDR2", "FWR3")

# Two-letter aliases that expand to a group of canonical labels.
# Resolution rule (audit §3): aliases expand first, then explicit
# labels override. So ``{"FWR": 0.5, "FWR2": 2.0}`` resolves to
# ``{FWR1: 0.5, FWR2: 2.0, FWR3: 0.5, CDR1: 1.0, CDR2: 1.0}``.
_V_SUBREGION_RATE_ALIASES: Dict[str, Tuple[str, ...]] = {
    "FWR": ("FWR1", "FWR2", "FWR3"),
    "CDR": ("CDR1", "CDR2"),
}

# Default flat-substrate rate vector (1.0 for every label) — same
# fast-path discipline as ``_DEFAULT_SEGMENT_RATES``.
_DEFAULT_V_SUBREGION_RATES: Tuple[float, float, float, float, float] = (
    1.0,
    1.0,
    1.0,
    1.0,
    1.0,
)


def _validate_v_subregion_rates(
    v_subregion_rates: Optional[Dict[str, float]],
) -> Tuple[float, float, float, float, float]:
    """Validate the user's ``v_subregion_rates`` dict and return
    the normalised ``(FWR1, CDR1, FWR2, CDR2, FWR3)`` tuple.
    Pure helper — no DSL state.

    Validation rules:

    - ``None`` (or omitted) → flat default
      ``(1.0, 1.0, 1.0, 1.0, 1.0)``.
    - Empty dict ``{}`` is treated as ``None`` — same flat default.
    - Keys must be a subset of the five canonical labels
      (``FWR1`` / ``CDR1`` / ``FWR2`` / ``CDR2`` / ``FWR3``) plus
      the two aliases ``FWR`` (expands to FWR1 / FWR2 / FWR3) and
      ``CDR`` (expands to CDR1 / CDR2). Case-sensitive — matches
      the V-subregion annotation surface (Slice 1).
    - Values must be ``int`` / ``float`` (not bool), finite, and
      ``>= 0``. NaN / inf / bool / negative raise ``ValueError``.
    - Sparse: omitted labels default to ``1.0``.
    - Alias expansion happens first, then explicit labels
      override. So ``{"FWR": 0.5, "FWR2": 2.0}`` → ``FWR1=0.5,
      FWR2=2.0, FWR3=0.5, CDR1=1.0, CDR2=1.0``.
    - After expansion, at least one label must be strictly
      positive — an all-zero vector would zero every V site and
      is almost certainly a builder bug. Reject with
      ``ValueError``.
    """
    if v_subregion_rates is None:
        return _DEFAULT_V_SUBREGION_RATES

    if not isinstance(v_subregion_rates, dict):
        raise TypeError(
            f"v_subregion_rates must be a dict or None, got "
            f"{type(v_subregion_rates).__name__}"
        )

    if not v_subregion_rates:
        # Empty dict is equivalent to omitting the kwarg.
        return _DEFAULT_V_SUBREGION_RATES

    accepted_keys = set(_V_SUBREGION_RATE_LABELS) | set(_V_SUBREGION_RATE_ALIASES)
    unknown = sorted(set(v_subregion_rates.keys()) - accepted_keys)
    if unknown:
        raise ValueError(
            f"v_subregion_rates: unknown label(s) {unknown!r}. "
            f"Allowed: {list(_V_SUBREGION_RATE_LABELS)!r} plus the "
            f"aliases {list(_V_SUBREGION_RATE_ALIASES.keys())!r}. "
            "Labels are case-sensitive; CDR3 / FWR4 are out of "
            "scope (CDR3 lives in the junction, FWR4 in the J "
            "segment) — use ``segment_rates`` for those."
        )

    def _coerce(label_for_error: str, value: object) -> float:
        if isinstance(value, bool) or not isinstance(value, (int, float)):
            raise TypeError(
                f"v_subregion_rates[{label_for_error!r}]: value must be a "
                f"finite non-negative number, got {type(value).__name__}"
            )
        v = float(value)
        if v != v:  # NaN
            raise ValueError(
                f"v_subregion_rates[{label_for_error!r}]: value must be a "
                "finite number, got NaN"
            )
        if v < 0.0 or v == float("inf") or v == float("-inf"):
            raise ValueError(
                f"v_subregion_rates[{label_for_error!r}]: value must be a "
                f"finite non-negative number, got {v}"
            )
        return v

    # Phase 1: alias expansion. Aliases populate their labels first;
    # phase 2's explicit labels then override.
    resolved: Dict[str, float] = {}
    for alias, expanded_labels in _V_SUBREGION_RATE_ALIASES.items():
        if alias in v_subregion_rates:
            v = _coerce(alias, v_subregion_rates[alias])
            for lbl in expanded_labels:
                resolved[lbl] = v

    # Phase 2: explicit labels override.
    for lbl in _V_SUBREGION_RATE_LABELS:
        if lbl in v_subregion_rates:
            resolved[lbl] = _coerce(lbl, v_subregion_rates[lbl])

    # Fill any remaining label with the flat default.
    out_list: List[float] = []
    for lbl in _V_SUBREGION_RATE_LABELS:
        out_list.append(resolved.get(lbl, 1.0))

    if sum(out_list) <= 0.0:
        raise ValueError(
            "v_subregion_rates: at least one label must have a "
            "positive rate after alias expansion. The supplied "
            "configuration zeroes every V subregion, which would "
            "drop every V site out of SHM support."
        )

    return (out_list[0], out_list[1], out_list[2], out_list[3], out_list[4])


# ──────────────────────────────────────────────────────────────────
# Clonal ordering-guard table (Slice: descendant-phase guards)
# ──────────────────────────────────────────────────────────────────
#
# Every DSL step in this table is **descendant-phase** — it models
# observation / library-prep / sequencing biology that must be
# sampled independently per clone member. Pre-fork placement either
# silently misreports the AIRR field (Bugs C / E / F: trace-sourced
# fields that don't survive the parent→descendant boundary) or
# collapses descendant diversity (every clone member shares an
# identical effect because the pass ran once on the parent IR).
#
# The clonal fork methods scan the already-appended step list against
# this table; the first match is rejected with a
# message naming the offending DSL method and the canonical fix.
#
# Each entry is ``(predicate, dsl_method_name)`` where the predicate
# inspects one step and returns ``True`` if it came from
# ``dsl_method_name``. Some DSL methods append :class:`_CorruptStep`
# with different ``kind`` discriminators; others append distinct
# step types.
def _descendant_phase_step_classifier(step):
    """Return the DSL method name a descendant-phase ``step`` came
    from, or ``None`` if ``step`` is not a descendant-phase step.

    Single source of truth for the unified guard in the flat clonal
    fork methods. Adding a new descendant-phase DSL method means
    appending a clause here (and adding the
    companion spec test in
    ``tests/test_clonal_descendant_phase_guards.py``).
    """
    from ._pipeline_ir import (
        _CORRUPT_KIND_3PRIME_LOSS,
        _CORRUPT_KIND_5PRIME_LOSS,
        _CORRUPT_KIND_INDEL,
        _CORRUPT_KIND_NS,
        _CORRUPT_KIND_PCR,
        _CORRUPT_KIND_QUALITY,
        _CORRUPT_KIND_REV_COMP,
        _CorruptStep,
        _MutateStep,
        _PairedEndStep,
    )

    if isinstance(step, _MutateStep):
        return "mutate"
    if isinstance(step, _PairedEndStep):
        return "paired_end"
    if isinstance(step, _CorruptStep):
        # Per-kind classification — only the descendant-phase kinds
        # appear in the table. ``contaminant`` is deliberately
        # omitted: the DSL slice that added the descendant-phase
        # ordering guards did not list it, so leave the placement
        # unconstrained until a follow-up explicitly classifies it.
        kind_to_method = {
            _CORRUPT_KIND_PCR: "pcr_amplify",
            _CORRUPT_KIND_QUALITY: "sequencing_errors",
            _CORRUPT_KIND_INDEL: "polymerase_indels",
            _CORRUPT_KIND_5PRIME_LOSS: "end_loss_5prime",
            _CORRUPT_KIND_3PRIME_LOSS: "end_loss_3prime",
            _CORRUPT_KIND_NS: "ambiguous_base_calls",
            _CORRUPT_KIND_REV_COMP: "random_strand_orientation",
        }
        return kind_to_method.get(step.kind)
    return None
