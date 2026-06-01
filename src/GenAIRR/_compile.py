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
    _InvertDStep,
    _MutateStep,
    _PairedEndStep,
    _ReceptorRevisionStep,
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

    :class:`_InvertDStep` is **inlined** into the preceding
    :class:`_RecombineStep` by ``Experiment._build_simulator`` before
    dispatch reaches this function — see
    :func:`_extract_invert_d_prob`. Reaching this branch with a bare
    :class:`_InvertDStep` means the caller forgot to apply the
    pre-pass, so we surface a structural error.
    """
    if isinstance(step, _RecombineStep):
        _lower_recombine(step, plan, refdata)
    elif isinstance(step, _MutateStep):
        _lower_mutate(step, plan)
    elif isinstance(step, _CorruptStep):
        _lower_corrupt(step, plan)
    elif isinstance(step, _InvertDStep):
        raise ValueError(
            "_InvertDStep reached lower_step(); the compile pipeline "
            "must inline the inversion probability into the preceding "
            "_RecombineStep via _extract_invert_d_prob() before lowering."
        )
    elif isinstance(step, _ReceptorRevisionStep):
        raise ValueError(
            "_ReceptorRevisionStep reached lower_step(); the compile "
            "pipeline must inline the receptor-revision probability "
            "into the preceding _RecombineStep via "
            "_extract_receptor_revision_prob() before lowering."
        )
    elif isinstance(step, _PairedEndStep):
        raise ValueError(
            "_PairedEndStep reached lower_step(); the compile pipeline "
            "must extract the step via _extract_paired_end_step() and "
            "push it at the END of the plan via _lower_paired_end()."
        )
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


def _extract_invert_d_prob(steps):
    """Pre-pass over the experiment step list.

    Scans for the (at-most-one) :class:`_InvertDStep`, returns its
    ``prob`` (or ``None`` if absent), and a copy of ``steps`` with
    that step removed. The compile loop drives main lowering off the
    *filtered* list and threads the inversion probability into
    :func:`_lower_recombine`, which then inlines
    ``push_invert_d(prob)`` into the canonical V-NP1-D-NP2-J
    recombination sequence at the position immediately before
    ``push_assemble("D")``.

    Why inline instead of schedule edges: the engine schedule
    analyser uses dependency edges + insertion order as the
    tiebreak. Pushing ``invert_d`` separately and then asking the
    analyser to reorder it via an explicit
    ``before(invert_d, assemble.d)`` edge would also re-promote
    ``invert_d`` past ``generate.np2`` / ``assemble.j``,
    accidentally pushing ``assemble.d`` to the END of the pool — a
    structural change to the V-NP1-D-NP2-J layout. Inlining at
    push time gets the position right by construction without any
    new schedule edges.

    Returns ``(invert_prob, filtered_steps)``.
    """
    invert_prob = None
    filtered = []
    for step in steps:
        if isinstance(step, _InvertDStep):
            if invert_prob is not None:
                raise ValueError(
                    "two _InvertDStep instances reached the compile "
                    "pipeline; Experiment.invert_d() rejects duplicates "
                    "at the DSL boundary, so this is a structural bug."
                )
            invert_prob = step.prob
            continue
        filtered.append(step)
    return invert_prob, filtered


def _extract_receptor_revision_prob(steps):
    """Pre-pass over the experiment step list, mirroring
    :func:`_extract_invert_d_prob`.

    Scans for the (at-most-one) :class:`_ReceptorRevisionStep`,
    returns its ``prob`` (or ``None`` if absent), and a copy of
    ``steps`` with that step removed. The compile loop drives main
    lowering off the *filtered* list and threads the revision
    probability into :func:`_lower_recombine`, which then inlines
    ``push_receptor_revision(prob, refdata)`` immediately after
    ``push_assemble("J")`` in the canonical V-NP1-D-NP2-J
    sequence.

    Why this position: receptor revision rewrites the V slice after
    the initial recombination is complete (V+NP1+D+NP2+J fully
    assembled). Placing the push after ``assemble.j`` means
    subsequent ``_MutateStep`` / ``_CorruptStep`` lowering — which
    happens *after* the recombine step in the user's pipeline —
    produces passes ordered after receptor revision in the plan,
    giving the canonical "recombine → revise → mutate/corrupt"
    order the design doc §2 requires.

    Returns ``(receptor_revision_prob, filtered_steps)``.
    """
    revision_prob = None
    filtered = []
    for step in steps:
        if isinstance(step, _ReceptorRevisionStep):
            if revision_prob is not None:
                raise ValueError(
                    "two _ReceptorRevisionStep instances reached the "
                    "compile pipeline; Experiment.receptor_revision() "
                    "rejects duplicates at the DSL boundary, so this "
                    "is a structural bug."
                )
            revision_prob = step.prob
            continue
        filtered.append(step)
    return revision_prob, filtered


def _lower_recombine(
    step: _RecombineStep,
    plan: "_engine.PassPlan",
    refdata: "_engine.RefDataConfig",
    *,
    invert_d_prob=None,
    receptor_revision_prob=None,
) -> None:
    chain = refdata.chain_type
    np1 = list(step.np1_lengths)
    np2 = list(step.np2_lengths)
    # Typed NP base distributions (Slice — Typed NP base model).
    # ``None`` means UniformBase applies (byte-identical to
    # pre-slice); a non-empty list of ``(base_byte, weight)``
    # pairs lowers to ``CategoricalBase`` at the bridge.
    np1_base_pairs = (
        list(step.np1_base_pairs) if step.np1_base_pairs is not None else None
    )
    np2_base_pairs = (
        list(step.np2_base_pairs) if step.np2_base_pairs is not None else None
    )
    # Slice — Markov NP Base Generator. ``None`` means uniform /
    # empirical_first_base (the existing bridge dispatch). Non-
    # ``None`` is a 4-row matrix in canonical A/C/G/T from-base
    # order; the bridge requires ``base_pairs`` (the first-base
    # row) to be set alongside.
    np1_markov_transitions = (
        [list(row) for row in step.np1_markov_transitions]
        if step.np1_markov_transitions is not None
        else None
    )
    np2_markov_transitions = (
        [list(row) for row in step.np2_markov_transitions]
        if step.np2_markov_transitions is not None
        else None
    )

    v_ids = list(step.locks_v) if step.locks_v is not None else None
    d_ids = list(step.locks_d) if step.locks_d is not None else None
    j_ids = list(step.locks_j) if step.locks_j is not None else None
    v_weights = list(step.weights_v) if step.weights_v is not None else None
    d_weights = list(step.weights_d) if step.weights_d is not None else None
    j_weights = list(step.weights_j) if step.weights_j is not None else None

    # P-nucleotide per-end lengths (Slice — P-nucleotide v1).
    # `None` means no P-extension authored for that end → the
    # corresponding `push_p_addition` call is omitted, byte-
    # identical to the pre-slice baseline. Lists are
    # materialised here so the bridge layer gets a plain
    # `List[Tuple[int, float]]` regardless of the frozen
    # tuple form `_RecombineStep` stores.
    p_v_3 = list(step.p_v_3_lengths) if step.p_v_3_lengths is not None else None
    p_d_5 = list(step.p_d_5_lengths) if step.p_d_5_lengths is not None else None
    p_d_3 = list(step.p_d_3_lengths) if step.p_d_3_lengths is not None else None
    p_j_5 = list(step.p_j_5_lengths) if step.p_j_5_lengths is not None else None

    if chain == "vj":
        if receptor_revision_prob is not None:
            # Defensive: the DSL boundary rejects VJ chains before a
            # _ReceptorRevisionStep is constructed. Reaching here would
            # mean the step was injected past the boundary check.
            raise ValueError(
                "receptor_revision is only valid for VDJ chains; the "
                "DSL boundary should have rejected this earlier."
            )
        plan.push_sample_allele("V", refdata, allowed_ids=v_ids, weights=v_weights)
        plan.push_sample_allele("J", refdata, allowed_ids=j_ids, weights=j_weights)
        if step.trim_v_3:
            plan.push_trim("V", "3", list(step.trim_v_3))
        if step.trim_j_5:
            plan.push_trim("J", "5", list(step.trim_j_5))
        plan.push_assemble("V")
        if p_v_3 is not None:
            plan.push_p_addition("V_3", p_v_3)
        plan.push_generate_np(
            "NP1",
            np1,
            base_pairs=np1_base_pairs,
            markov_transitions=np1_markov_transitions,
        )
        if p_j_5 is not None:
            plan.push_p_addition("J_5", p_j_5)
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
        if p_v_3 is not None:
            plan.push_p_addition("V_3", p_v_3)
        plan.push_generate_np(
            "NP1",
            np1,
            base_pairs=np1_base_pairs,
            markov_transitions=np1_markov_transitions,
        )
        # Inversion decision commits before D's effective_seq
        # is read. invert_d_prob fires the D-inversion pass
        # which updates `AlleleInstance.orientation`. P-addition
        # at D_5 derives its palindrome from D's post-orientation
        # 5' flank, so it MUST run after invert_d — otherwise
        # the inverted-D path would palindrome the wrong end.
        # Order: invert_d → p_addition.d_5 → assemble.D.
        # (Audit's §9.2 diagram had this swap inverted; the
        # implementation slice corrected the order per
        # `docs/p_nucleotide_design.md` §9.3.)
        if invert_d_prob is not None:
            plan.push_invert_d(invert_d_prob)
        if p_d_5 is not None:
            plan.push_p_addition("D_5", p_d_5)
        plan.push_assemble("D")
        if p_d_3 is not None:
            plan.push_p_addition("D_3", p_d_3)
        plan.push_generate_np(
            "NP2",
            np2,
            base_pairs=np2_base_pairs,
            markov_transitions=np2_markov_transitions,
        )
        if p_j_5 is not None:
            plan.push_p_addition("J_5", p_j_5)
        plan.push_assemble("J")
        # Receptor revision rewrites the V slice **after** initial
        # V/D/J/NP assembly is complete. Inlining the push here puts
        # the pass between `assemble.j` and any subsequent
        # `_MutateStep` / `_CorruptStep` lowering, preserving the
        # canonical "recombine → revise → mutate/corrupt" order the
        # design doc §2 specifies. See
        # `_extract_receptor_revision_prob` for the full rationale on
        # inlining instead of a standalone _ReceptorRevisionStep lower
        # path.
        if receptor_revision_prob is not None:
            plan.push_receptor_revision(receptor_revision_prob, refdata)
    else:  # pragma: no cover — RefDataConfig only constructs vj/vdj.
        raise ValueError(f"unsupported chain_type {chain!r}")


def _lower_mutate(step: _MutateStep, plan: "_engine.PassPlan") -> None:
    # Segment-rate vector is always present on the step (default
    # ``(1.0, 1.0, 1.0, 1.0)``) so the engine bindings receive a
    # uniform shape. The Rust passes detect the flat-default case
    # and short-circuit to the existing fast path.
    seg = step.segment_rates
    # V-subregion-rate vector (Slice B). Same discipline: always
    # present on the step (default `(1, 1, 1, 1, 1)`), shape-
    # uniform across the bridge; the engine fast-paths the flat
    # default.
    v_sub = step.v_subregion_rates
    if step.model == "uniform":
        if step.rate is not None:
            plan.push_mutate_uniform_rate(
                step.rate,
                segment_rates=list(seg),
                v_subregion_rates=list(v_sub),
            )
        else:
            plan.push_mutate_uniform(
                list(step.count_pairs or ()),
                segment_rates=list(seg),
                v_subregion_rates=list(v_sub),
            )
    elif step.model == "s5f":
        from ._s5f_loader import load_builtin_s5f_kernel

        mutability, substitution = load_builtin_s5f_kernel(step.s5f_model_name)
        # Thread the kernel's short name through to the engine so
        # the plan signature can carry the kernel identity (Slice A —
        # Pass Parameter Signature). Without this, two runs that
        # share rate + segment_rates but disagree on which S5F
        # kernel was loaded would replay against each other silently.
        if step.rate is not None:
            plan.push_mutate_s5f_rate(
                step.rate,
                mutability,
                substitution,
                segment_rates=list(seg),
                v_subregion_rates=list(v_sub),
                kernel_name=step.s5f_model_name,
            )
        else:
            plan.push_mutate_s5f(
                list(step.count_pairs or ()),
                mutability,
                substitution,
                segment_rates=list(seg),
                v_subregion_rates=list(v_sub),
                kernel_name=step.s5f_model_name,
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


def _extract_paired_end_step(steps):
    """Pre-pass over the experiment step list, mirroring
    :func:`_extract_invert_d_prob` /
    :func:`_extract_receptor_revision_prob`.

    Scans for the (at-most-one) :class:`_PairedEndStep`, returns
    the step itself (or ``None`` if absent) and a copy of
    ``steps`` with that step removed. The compile loop drives main
    lowering off the *filtered* list and pushes the paired-end
    pass at the END of the plan via :func:`_lower_paired_end`.

    Why this position: the paired-end pass is trace-only, so its
    placement doesn't affect the simulation. But it samples three
    Ints that the AIRR builder reads to populate
    :attr:`AirrRecord.r1_sequence` / :attr:`r2_sequence` / etc.
    Recording the choices AFTER every biology / corruption /
    rev-comp step keeps the trace order aligned with the
    biological readout order — design doc §6 specifies the
    pipeline as ``…recombination → mutation → corruption →
    end_loss → paired_end projection`` and the trace must
    reflect that.

    Returns ``(paired_end_step, filtered_steps)``.
    """
    paired_end_step = None
    filtered = []
    for step in steps:
        if isinstance(step, _PairedEndStep):
            if paired_end_step is not None:
                raise ValueError(
                    "two _PairedEndStep instances reached the compile "
                    "pipeline; Experiment.paired_end() rejects "
                    "duplicates at the DSL boundary, so this is a "
                    "structural bug."
                )
            paired_end_step = step
            continue
        filtered.append(step)
    return paired_end_step, filtered


def _lower_paired_end(
    step: _PairedEndStep,
    plan: "_engine.PassPlan",
) -> None:
    """Push the paired-end sampling pass at the **end** of the
    plan. Lowering happens after the main per-step loop in
    :meth:`Experiment._build_simulator` so this push lands after
    every other pass in insertion order."""
    plan.push_paired_end(
        r1_length_pairs=list(step.r1_length),
        r2_length_pairs=list(step.r2_length),
        insert_size_pairs=list(step.insert_size),
    )
