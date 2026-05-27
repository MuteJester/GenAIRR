"""Lowering layer between the pure-Python pipeline IR and the Rust
:class:`_engine.PassPlan`.

:mod:`._pipeline_ir` owns the typed step dataclasses as **pure data**
— it has no awareness of the engine ABI. This module is the single
place that knows how to lower each step kind onto an engine-native
:class:`_engine.PassPlan`.

The split lets the pipeline IR be inspected, round-tripped, or
serialized without depending on the compiled Rust extension, which
is the precondition for future work on a `ExperimentSpec` surface,
event/observer pipelines, and durable trace/replay schemas.
"""
from __future__ import annotations

from GenAIRR import _engine  # private Rust extension submodule

from ._pipeline_ir import (
    _CORRUPT_KIND_3PRIME_LOSS,
    _CORRUPT_KIND_5PRIME_LOSS,
    _CORRUPT_KIND_CONTAMINANT,
    _CORRUPT_KIND_INDEL,
    _CORRUPT_KIND_NS,
    _CORRUPT_KIND_PCR,
    _CORRUPT_KIND_QUALITY,
    _CORRUPT_KIND_REV_COMP,
    _ClonalForkStep,
    _CorruptStep,
    _MutateStep,
    _RecombineStep,
)


def lower_step(
    step: object,
    plan: "_engine.PassPlan",
    refdata: "_engine.RefDataConfig",
) -> None:
    """Dispatch on step type and lower into the engine `PassPlan`.

    :class:`_ClonalForkStep` is a control-flow marker and produces no
    engine passes; the experiment's compile loop splits the step
    sequence at the marker and lowers each sub-sequence into its own
    `PassPlan`.
    """
    if isinstance(step, _RecombineStep):
        _lower_recombine(step, plan, refdata)
    elif isinstance(step, _MutateStep):
        _lower_mutate(step, plan)
    elif isinstance(step, _CorruptStep):
        _lower_corrupt(step, plan)
    elif isinstance(step, _ClonalForkStep):
        # Marker only — `Experiment.compile()` partitions the step
        # list at this point. Reaching here means the caller didn't
        # partition, which is a structural bug.
        raise ValueError(
            "_ClonalForkStep reached lower_step(); the compile pipeline "
            "must split the step list at this marker and lower each "
            "partition independently."
        )
    else:  # pragma: no cover — step kinds are guarded at builder time.
        raise TypeError(f"unsupported pipeline step type {type(step).__name__!r}")


def _lower_recombine(
    step: _RecombineStep,
    plan: "_engine.PassPlan",
    refdata: "_engine.RefDataConfig",
) -> None:
    chain = refdata.chain_type
    np1 = list(step.np1_lengths)
    np2 = list(step.np2_lengths)

    v_ids = list(step.locks_v) if step.locks_v is not None else None
    d_ids = list(step.locks_d) if step.locks_d is not None else None
    j_ids = list(step.locks_j) if step.locks_j is not None else None
    v_weights = list(step.weights_v) if step.weights_v is not None else None
    d_weights = list(step.weights_d) if step.weights_d is not None else None
    j_weights = list(step.weights_j) if step.weights_j is not None else None

    if chain == "vj":
        plan.push_sample_allele("V", refdata, allowed_ids=v_ids, weights=v_weights)
        plan.push_sample_allele("J", refdata, allowed_ids=j_ids, weights=j_weights)
        if step.trim_v_3:
            plan.push_trim("V", "3", list(step.trim_v_3))
        if step.trim_j_5:
            plan.push_trim("J", "5", list(step.trim_j_5))
        plan.push_assemble("V")
        plan.push_generate_np("NP1", np1)
        plan.push_assemble("J")
    elif chain == "vdj":
        plan.push_sample_allele("V", refdata, allowed_ids=v_ids, weights=v_weights)
        plan.push_sample_allele("D", refdata, allowed_ids=d_ids, weights=d_weights)
        plan.push_sample_allele("J", refdata, allowed_ids=j_ids, weights=j_weights)
        if step.trim_v_3:
            plan.push_trim("V", "3", list(step.trim_v_3))
        if step.trim_d_5:
            plan.push_trim("D", "5", list(step.trim_d_5))
        if step.trim_d_3:
            plan.push_trim("D", "3", list(step.trim_d_3))
        if step.trim_j_5:
            plan.push_trim("J", "5", list(step.trim_j_5))
        plan.push_assemble("V")
        plan.push_generate_np("NP1", np1)
        plan.push_assemble("D")
        plan.push_generate_np("NP2", np2)
        plan.push_assemble("J")
    else:  # pragma: no cover — RefDataConfig only constructs vj/vdj.
        raise ValueError(f"unsupported chain_type {chain!r}")


def _lower_mutate(step: _MutateStep, plan: "_engine.PassPlan") -> None:
    if step.model == "uniform":
        if step.rate is not None:
            plan.push_mutate_uniform_rate(step.rate)
        else:
            plan.push_mutate_uniform(list(step.count_pairs or ()))
    elif step.model == "s5f":
        from ._s5f_loader import load_builtin_s5f_kernel

        mutability, substitution = load_builtin_s5f_kernel(step.s5f_model_name)
        if step.rate is not None:
            plan.push_mutate_s5f_rate(step.rate, mutability, substitution)
        else:
            plan.push_mutate_s5f(
                list(step.count_pairs or ()), mutability, substitution
            )
    else:  # pragma: no cover — guarded at builder time.
        raise ValueError(f"unsupported mutation model {step.model!r}")


def _lower_corrupt(step: _CorruptStep, plan: "_engine.PassPlan") -> None:
    if step.kind == _CORRUPT_KIND_PCR:
        if step.rate is not None:
            plan.push_corrupt_pcr_rate(step.rate)
        else:
            plan.push_corrupt_pcr(list(step.count_pairs or ()))
    elif step.kind == _CORRUPT_KIND_QUALITY:
        if step.rate is not None:
            plan.push_corrupt_quality_rate(step.rate)
        else:
            plan.push_corrupt_quality(list(step.count_pairs or ()))
    elif step.kind == _CORRUPT_KIND_INDEL:
        plan.push_corrupt_indel(
            list(step.count_pairs or ()), insertion_prob=step.insertion_prob
        )
    elif step.kind == _CORRUPT_KIND_CONTAMINANT:
        plan.push_corrupt_contaminant(step.apply_prob)
    elif step.kind == _CORRUPT_KIND_REV_COMP:
        plan.push_corrupt_rev_comp(step.apply_prob)
    elif step.kind == _CORRUPT_KIND_5PRIME_LOSS:
        plan.push_corrupt_5prime_loss(list(step.count_pairs or ()))
    elif step.kind == _CORRUPT_KIND_3PRIME_LOSS:
        plan.push_corrupt_3prime_loss(list(step.count_pairs or ()))
    elif step.kind == _CORRUPT_KIND_NS:
        plan.push_corrupt_ns(list(step.count_pairs or ()))
    else:  # pragma: no cover — guarded at builder time.
        raise ValueError(f"unsupported corruption kind {step.kind!r}")
