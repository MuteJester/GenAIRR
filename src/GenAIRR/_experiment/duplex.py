"""DSL checks for finite authored duplex recombination."""
from __future__ import annotations

from .._pipeline_ir import (
    _RecombineStep, _InvertDStep, _MutateStep, _CorruptStep,
    _ClonalForkStep, _LineageForkStep, _RepertoireForkStep,
    _CORRUPT_KIND_CONTAMINANT,
)


def append_duplex_recombine(experiment, model, *, np1_lengths, np2_lengths,
                            v_allele_weights, d_allele_weights, j_allele_weights):
    from .._duplex_model import DuplexRepairModel

    if not isinstance(model, DuplexRepairModel):
        raise TypeError("duplex_repair must be a DuplexRepairModel")
    if experiment.chain_type != "vdj":
        raise ValueError("duplex repair currently requires VDJ reference data")
    if np1_lengths is not None or np2_lengths is not None:
        raise ValueError("duplex_repair and explicit legacy NP lengths are mutually exclusive")
    explicit_weights = any(w is not None for w in (v_allele_weights, d_allele_weights, j_allele_weights))
    if explicit_weights and experiment._genotype is not None:
        raise ValueError("recombine(*_allele_weights=...) and with_genotype() are mutually exclusive")
    # Keep cartridge allele usage as an upstream rule. Junction defaults are
    # deliberately absent from this immutable step: the authored history owns
    # the complete repair product, rather than a second N/P/trim model.
    defaults = experiment._recombine_defaults() or {}
    values = []
    profile_defaults = []
    for segment, weights in (("V", v_allele_weights), ("D", d_allele_weights), ("J", j_allele_weights)):
        if weights is None:
            weights = defaults.get(f"allele_usage_{segment.lower()}")
            if weights is not None and experiment._reference is not None:
                profile_defaults.append(f"allele_usage_{segment}")
        values.append(experiment._resolve_allele_weights(segment, weights))
    defaults_profile = None
    if experiment._reference is not None:
        from .._reference.defaults import PROFILE as defaults_profile
    step = _RecombineStep(
        np1_lengths=(), np2_lengths=(), trim_v_3=None, trim_d_5=None,
        trim_d_3=None, trim_j_5=None, weights_v=values[0], weights_d=values[1],
        weights_j=values[2], duplex_repair_json=model.to_json(),
        defaults_profile=defaults_profile, profile_defaults=tuple(profile_defaults),
    )
    if explicit_weights:
        experiment._user_allele_weights_set = True
    experiment._steps.append(step)
    return experiment


def validate_duplex_steps(steps):
    sources = [s for s in steps if isinstance(s, _RecombineStep) and s.duplex_repair_json is not None]
    if not sources:
        return
    if len(sources) != 1 or sum(isinstance(s, _RecombineStep) for s in steps) != 1:
        raise ValueError("duplex repair requires exactly one recombination source")
    if not steps or steps[0] is not sources[0]:
        raise ValueError("duplex recombination must precede downstream steps")
    source = sources[0]
    if (source.np_lengths_overridden or source.trim_overridden or source.coding_end_mode is not None
            or source.np1_lengths or source.np2_lengths
            or any(getattr(source, field) is not None for field in (
                "trim_v_3", "trim_d_5", "trim_d_3", "trim_j_5", "p_v_3_lengths",
                "p_d_5_lengths", "p_d_3_lengths", "p_j_5_lengths", "np1_base_pairs",
                "np2_base_pairs", "np1_markov_transitions", "np2_markov_transitions"))):
        raise ValueError("duplex repair excludes explicit legacy trim, P and NP processing")
    forks = [step for step in steps if isinstance(
        step, (_ClonalForkStep, _LineageForkStep, _RepertoireForkStep)
    )]
    if len(forks) > 1:
        raise ValueError("duplex repair supports at most one family boundary")
    fork_seen = False
    for step in steps[1:]:
        if isinstance(step, (_ClonalForkStep, _LineageForkStep, _RepertoireForkStep)):
            fork_seen = True
            continue
        if isinstance(step, _InvertDStep) and fork_seen:
            raise ValueError("D inversion must precede the duplex family boundary")
        if isinstance(step, (_InvertDStep, _MutateStep)):
            continue
        if isinstance(step, _CorruptStep) and step.kind != _CORRUPT_KIND_CONTAMINANT:
            continue
        raise ValueError(
            f"{type(step).__name__} is unsupported with duplex repair; "
            "receptor revision, paired-end reads and contamination "
            "require additional provenance adapters"
        )
