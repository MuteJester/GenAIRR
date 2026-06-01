"""``Experiment.describe()`` rendering helpers.

Each step kind has a `_describe_*_step` function that converts the
typed IR dataclass into a biology-style narrative line. The output
is the diagnostic instrument we use to grade DSL readability: if a
chain reads weird in :meth:`describe`, the chain itself reads weird.

Also contains :func:`_format_declared_contracts` /
:func:`_format_active_contracts` for the contract-bundle line, plus
:func:`_describe_experiment_header` for the one-line refdata-source
header.

This module is private. The output strings are part of the
user-visible :meth:`Experiment.describe()` surface and should remain
stable; the helper signatures here may change without notice.
"""
from __future__ import annotations

from typing import Any, List, Optional, Sequence, Tuple, TYPE_CHECKING

from ._pipeline_ir import (
    _CORRUPT_KIND_3PRIME_LOSS,
    _CORRUPT_KIND_5PRIME_LOSS,
    _CORRUPT_KIND_CONTAMINANT,
    _CORRUPT_KIND_INDEL,
    _CORRUPT_KIND_NS,
    _CORRUPT_KIND_PCR,
    _CORRUPT_KIND_QUALITY,
    _CORRUPT_KIND_REV_COMP,
    _DEFAULT_NP_LENGTHS,
    _ClonalForkStep,
    _CorruptStep,
    _MutateStep,
    _RecombineStep,
)

if TYPE_CHECKING:  # pragma: no cover — type-only imports
    from GenAIRR import _engine
    from .dataconfig import DataConfig


def _describe_distribution(pairs: Tuple[Tuple[int, float], ...]) -> str:
    """Render a ``(value, weight)`` distribution as a short
    biology-friendly string: ``"5"`` for a fixed count, ``"5–15"``
    for a uniform range, ``"0–8 weighted"`` for non-uniform."""
    if not pairs:
        return "<empty>"
    values = [v for v, _ in pairs]
    weights = [w for _, w in pairs]
    lo, hi = min(values), max(values)
    if lo == hi:
        return str(lo)
    uniform = all(abs(w - weights[0]) < 1e-9 for w in weights)
    covers_range = sorted(values) == list(range(lo, hi + 1))
    if uniform and covers_range:
        return f"{lo}–{hi}"
    return f"{lo}–{hi} weighted"


def _describe_recombine_step(step: "_RecombineStep", chain_type: str) -> str:
    segs = "V/J" if chain_type == "vj" else "V/D/J"
    parts: List[str] = [f"sample {segs} alleles"]

    locked_segs = []
    if step.locks_v is not None:
        locked_segs.append(f"V ({len(step.locks_v)})")
    if step.locks_d is not None:
        locked_segs.append(f"D ({len(step.locks_d)})")
    if step.locks_j is not None:
        locked_segs.append(f"J ({len(step.locks_j)})")
    if locked_segs:
        parts.append(f"locked to {', '.join(locked_segs)}")

    weighted_segs = []
    for seg, weights in (("V", step.weights_v), ("D", step.weights_d), ("J", step.weights_j)):
        if weights is not None:
            weighted_segs.append(seg)
    if weighted_segs:
        parts.append(f"custom {'/'.join(weighted_segs)} allele weights")

    trim_bits = []
    if step.trim_v_3:
        trim_bits.append("V3'")
    if step.trim_d_5:
        trim_bits.append("D5'")
    if step.trim_d_3:
        trim_bits.append("D3'")
    if step.trim_j_5:
        trim_bits.append("J5'")
    if trim_bits:
        parts.append(f"empirical exonuclease trim ({', '.join(trim_bits)})")

    np1_desc = _describe_distribution(step.np1_lengths)
    np1_source = _np_source_tag(step.np1_lengths)
    if chain_type == "vdj":
        np2_desc = _describe_distribution(step.np2_lengths)
        np2_source = _np_source_tag(step.np2_lengths)
        parts.append(
            f"insert NP1 ({np1_desc} bases{np1_source}) and "
            f"NP2 ({np2_desc} bases{np2_source})"
        )
    else:
        parts.append(f"insert NP1 ({np1_desc} bases{np1_source})")

    return "V(D)J recombination: " + "; ".join(parts)


def _np_source_tag(np_lengths: Tuple[Tuple[int, float], ...]) -> str:
    """Return ``", empirical"`` or ``", uniform fallback"`` to clarify
    whether the NP-length distribution came from the bound DataConfig
    (empirical) or the synthetic [(0,1.0)..(6,1.0)] fallback used on
    raw RefDataConfig with no empirical data. Empty string when the
    user supplied an explicit distribution we can't classify."""
    if tuple(np_lengths) == tuple(_DEFAULT_NP_LENGTHS):
        return ", uniform fallback"
    # User-supplied or empirical; we can't distinguish them at this
    # point without threading more state, but the "weighted" suffix
    # in the distribution describer already tells the reader it's a
    # non-uniform empirical distribution when the weights aren't equal.
    return ", empirical"


_S5F_KERNEL_LABELS = {
    "hh_s5f": "human heavy-chain (HH_S5F)",
    "hkl_s5f": "human kappa/lambda (HKL_S5F)",
    "mk_rs5nf": "mouse (MK_RS5NF)",
}


def _describe_mutate_step(step: "_MutateStep") -> str:
    if step.rate is not None:
        intensity = f"{step.rate:.1%} of bases (Poisson)"
    elif step.count_pairs is not None:
        intensity = f"{_describe_distribution(step.count_pairs)} mutations"
    else:
        intensity = "<unset>"  # pragma: no cover — guarded at builder time
    if step.model == "s5f":
        kernel_label = _S5F_KERNEL_LABELS.get(step.s5f_model_name, step.s5f_model_name)
        return f"Somatic hypermutation (S5F context model, {kernel_label}): {intensity}/record"
    if step.model == "uniform":
        return f"Somatic hypermutation (uniform, position-independent): {intensity}/record"
    return f"Mutation ({step.model}): {intensity}/record"


_CORRUPT_NARRATIVE_LABELS = {
    _CORRUPT_KIND_PCR: ("PCR substitution errors", "errors/record"),
    _CORRUPT_KIND_QUALITY: ("Sequencing quality errors", "errors/record"),
    _CORRUPT_KIND_5PRIME_LOSS: ("5'-end loss", "bases trimmed"),
    _CORRUPT_KIND_3PRIME_LOSS: ("3'-end loss", "bases trimmed"),
    _CORRUPT_KIND_NS: ("Ambiguous base calls (low-Q positions)", "positions/record"),
}


def _describe_corrupt_step(step: "_CorruptStep") -> str:
    if step.kind in _CORRUPT_NARRATIVE_LABELS:
        label, unit = _CORRUPT_NARRATIVE_LABELS[step.kind]
        intensity = _describe_count_or_rate(step, default_unit=unit)
        return f"{label}: {intensity}"
    if step.kind == _CORRUPT_KIND_INDEL:
        count = _describe_distribution(step.count_pairs or ())
        return (
            f"Polymerase indels (PCR slippage): {count} events/record, "
            f"insertion fraction {step.insertion_prob:.2f}"
        )
    if step.kind == _CORRUPT_KIND_CONTAMINANT:
        return (
            f"Cross-sample contamination: {step.apply_prob:.0%} of records "
            f"replaced with a random allele sequence"
        )
    if step.kind == _CORRUPT_KIND_REV_COMP:
        return f"Random strand orientation: {step.apply_prob:.0%} of records reverse-complemented"
    return f"Corruption ({step.kind})"  # pragma: no cover — guarded at builder time


def _describe_count_or_rate(step: "_CorruptStep", *, default_unit: str) -> str:
    """Render the count / rate / events portion of a corruption step
    line, depending on which mode was selected at builder time."""
    if step.rate is not None:
        return f"{step.rate:.2%} of bases (Poisson)/record"
    return f"{_describe_distribution(step.count_pairs or ())} {default_unit}"


def _describe_clonal_fork_step(step: "_ClonalForkStep") -> str:
    return f"Clonal expansion: {step.n_clones} lineages × {step.size} reads/clone"


def _describe_step(step: Any, chain_type: str) -> str:
    """Dispatch to the right step-describer based on the step type."""
    if isinstance(step, _RecombineStep):
        return _describe_recombine_step(step, chain_type)
    if isinstance(step, _MutateStep):
        return _describe_mutate_step(step)
    if isinstance(step, _CorruptStep):
        return _describe_corrupt_step(step)
    if isinstance(step, _ClonalForkStep):
        return _describe_clonal_fork_step(step)
    return f"Unknown step: {type(step).__name__}"  # pragma: no cover


_CHAIN_TYPE_BIOLOGY_LABELS = {
    "BCR_HEAVY": "heavy-chain BCR",
    "BCR_LIGHT_KAPPA": "kappa light-chain BCR",
    "BCR_LIGHT_LAMBDA": "lambda light-chain BCR",
    "TCR_ALPHA": "alpha-chain TCR",
    "TCR_BETA": "beta-chain TCR",
    "TCR_GAMMA": "gamma-chain TCR",
    "TCR_DELTA": "delta-chain TCR",
}


def _describe_experiment_header(
    refdata: "_engine.RefDataConfig",
    dataconfig: "Optional[DataConfig]",
) -> str:
    """One-line header naming the refdata source and chain type."""
    chain = refdata.chain_type
    if dataconfig is not None and getattr(dataconfig, "metadata", None) is not None:
        meta = dataconfig.metadata
        species_label: Optional[str] = None
        species = getattr(meta, "species", None)
        if species is not None:
            species_label = getattr(species, "value", None) or str(species)
        chain_enum = getattr(meta, "chain_type", None)
        chain_key = getattr(chain_enum, "value", None) if chain_enum is not None else None
        chain_label = _CHAIN_TYPE_BIOLOGY_LABELS.get(chain_key or "", f"{chain.upper()} chain")
        if species_label:
            return f"Experiment on {species_label} {chain_label}"
        return f"Experiment on {chain_label}"
    return f"Experiment on raw RefDataConfig ({chain.upper()})"


def _describe_step_sequence(
    steps: Sequence[Any],
    chain_type: str,
    *,
    indent: str = "  ",
    start_index: int = 1,
) -> List[str]:
    """Render a sequence of steps as numbered narrative lines. The
    clonal-fork step gets a visual divider so the pre/post boundary
    is obvious; subsequent steps continue numbering.

    Returns one string per output line (no trailing newlines).
    """
    out: List[str] = []
    step_no = start_index
    for step in steps:
        if isinstance(step, _ClonalForkStep):
            out.append(f"{indent}── {_describe_clonal_fork_step(step)} ──")
            out.append(f"{indent}    (steps above run once per clone; "
                       f"steps below run once per descendant)")
            continue
        out.append(f"{indent}{step_no}. {_describe_step(step, chain_type)}")
        step_no += 1
    return out


_PRODUCTIVE_BUNDLE_NAMES = frozenset({
    "productive_junction_frame",
    "no_stop_codon_in_junction",
    "anchor_preserved.v",
    "anchor_preserved.j",
})


def _format_declared_contracts(declared: Sequence[str]) -> Optional[str]:
    """Render the contract bundles declared by `.productive_only()`
    etc. on an uncompiled Experiment as a biology-style line.

    Companion to :func:`_format_active_contracts`, which renders the
    compiled simulator's per-pass contract names (the underlying
    primitives the bundles expand into). On the Experiment side we
    show the bundle name directly because that's what the user typed.
    """
    if not declared:
        return None
    labels = []
    for bundle in declared:
        if bundle == "productive":
            labels.append("productive sequences only")
        else:
            labels.append(bundle)
    return "; ".join(labels)


def _format_active_contracts(active: Sequence[str]) -> Optional[str]:
    """Render the compiled simulator's active-contracts tuple as a
    biology-style line, or ``None`` if no contracts are bound. The
    `ga.productive()` bundle's four underlying contracts
    (productive_junction_frame, no_stop_codon_in_junction, and the
    two anchor_preserved contracts) collapse to the single label
    "productive sequences only" so the line reads naturally."""
    if not active:
        return None
    names = list(active)
    extras: List[str] = []
    productive_present = _PRODUCTIVE_BUNDLE_NAMES.issubset(set(names))
    bundle_names = _PRODUCTIVE_BUNDLE_NAMES if productive_present else frozenset()
    for n in names:
        if n in bundle_names:
            continue
        extras.append(n)
    labels: List[str] = []
    if productive_present:
        labels.append("productive sequences only")
    labels.extend(extras)
    if not labels:
        return None
    return "; ".join(labels)
