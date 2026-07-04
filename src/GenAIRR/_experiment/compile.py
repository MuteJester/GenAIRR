from __future__ import annotations

import warnings
from typing import Optional

from GenAIRR import _engine

from .._compiled import (
    CompiledClonalExperiment,
    CompiledExperiment,
    CompiledLineageExperiment,
    CompiledRepertoireExperiment,
)
from .._lowering import (
    _extract_invert_d_prob,
    _extract_paired_end_step,
    _extract_receptor_revision_prob,
    _lower_paired_end,
    _lower_recombine,
    lower_step,
)
from .._pipeline_ir import (
    _ClonalForkStep,
    _CorruptStep,
    _LineageForkStep,
    _MutateStep,
    _PairedEndStep,
    _RecombineStep,
    _RepertoireForkStep,
)
from .._refdata_resolver import dataconfig_to_refdata


class _CompileMixin:
    __slots__ = ()

    def compile(self, *, allow_curatable_refdata: Optional[bool] = None):
        """Compile the recorded steps into a reusable
        :class:`CompiledExperiment` (or :class:`CompiledClonalExperiment`
        when the pipeline contains a :meth:`expand_clones`
        fork).

        Idempotent: calling ``compile()`` twice produces two distinct
        compiled instances with structurally-equal simulators.

        Constraints declared via :meth:`productive_only` (or future
        bundle methods) are baked into the compiled simulator at this
        step; they're not runtime knobs. To run without constraints,
        omit the constraint methods from the chain.

        ``allow_curatable_refdata`` selects the refdata validation
        mode. ``None`` (default) inherits the instance flag set by
        :meth:`allow_curatable_refdata`; an explicit ``True`` /
        ``False`` overrides per-call. ``False`` runs the gate in
        strict mode — every issue rejects compile with a
        :class:`ValueError`. ``True`` runs the lenient mode — Fatal
        issues (empty pool, duplicates, invalid byte, anchor out of
        bounds) still reject, but Curatable issues (pseudogene-shape
        anchor anomalies) pass.
        """
        if allow_curatable_refdata is None:
            allow_curatable_refdata = self._allow_curatable_refdata

        # Receptor revision with a phased genotype is now haplotype-aware:
        # the lowering builds a genotype-aware ReceptorRevisionPass that
        # restricts the replacement V to carried alleles on the drawn
        # rearrangement chromosome (see _lower_recombine). No rejection here.

        # Genotype provenance (subject_id / haplotype / result.genotypes)
        # is only threaded through the plain compiled path, not the
        # clonal/lineage/repertoire forked classes. Reject the
        # combination rather than silently dropping provenance (review
        # #9); genotype + clonal cohorts are a planned follow-on.
        if self._genotype is not None and self._has_clonal_fork():
            raise ValueError(
                "with_genotype() is not supported together with expand_clones() / "
                "clonal_lineage() / clonal_repertoire() in this release"
            )
        from dataclasses import replace as _replace

        # On raw RefDataConfig with default-on trim, warn at compile
        # time if the user didn't explicitly call .trim() to either
        # disable or replace the no-op default. By compile() time we
        # know the final shape of the recombine step.
        if self._dataconfig is None:
            for step in self._steps:
                if not isinstance(step, _RecombineStep):
                    continue
                has_trim_data = any(
                    (step.trim_v_3, step.trim_d_5, step.trim_d_3, step.trim_j_5)
                )
                if has_trim_data or step.trim_overridden:
                    # Either trim is set up (DataConfig path or
                    # .trim(v_3=...)) or the user explicitly called
                    # .trim() — both silence the warning.
                    continue
                warnings.warn(
                    "Experiment bound to a raw RefDataConfig has no "
                    "trim distributions; exonuclease trim is a no-op. "
                    "Call .trim(enabled=False) to silence this warning, "
                    "or supply custom distributions via "
                    ".trim(v_3=..., d_5=..., ...).",
                    UserWarning,
                    stacklevel=2,
                )
                break  # one warning per compile() call is enough

        contracts = self._build_contracts()
        any_lock = any(self._locks[seg] is not None for seg in ("V", "D", "J"))

        # if a `_ClonalForkStep` is present, split the step list
        # at it and compile two simulators (pre-fork = per-clone,
        # post-fork = per-descendant).
        fork_idx = next(
            (i for i, s in enumerate(self._steps) if isinstance(s, _ClonalForkStep)),
            None,
        )
        if fork_idx is not None:
            fork_step: _ClonalForkStep = self._steps[fork_idx]
            pre_steps = self._steps[:fork_idx]
            post_steps = self._steps[fork_idx + 1 :]
            pre_simulator = self._build_simulator(
                pre_steps,
                contracts,
                any_lock,
                replace_fn=_replace,
                allow_curatable_refdata=allow_curatable_refdata,
            )
            # The post-fork plan inherits the parent's V/D/J/NP
            # backbone (recombination already happened in the
            # pre-fork half). The recombination-time precondition
            # facts (np.np1.length residues, anchor trim supports)
            # aren't produced here. v2 used to drop the contract
            # bundle entirely on this side to avoid the compile-time
            # precondition failure — but doing so left every
            # post-fork mutation / corruption pass unfiltered, which
            # was the dominant cause of non-productive output under
            # `productive_only()`. We now pass the same bundle
            # through; the engine analyzer (in
            # `compiled/analyze.rs::validate_contract_preconditions`)
            # skips the productive-frame check when no recombination
            # facts are present in the plan, and the runtime
            # admit_with_context / admits_post_event paths still
            # enforce no-stop-codon-in-junction and anchor-preserved
            # for every substitution and indel in the post-fork pipeline.
            post_simulator = self._build_simulator(
                post_steps,
                contracts,
                any_lock=False,
                replace_fn=_replace,
                allow_curatable_refdata=allow_curatable_refdata,
            )
            return CompiledClonalExperiment(
                pre_simulator,
                post_simulator,
                fork_step,
                self._refdata,
                pre_steps=tuple(pre_steps),
                post_steps=tuple(post_steps),
                dataconfig=self._dataconfig,
                metadata=self._metadata,
            )

        # if a `_RepertoireForkStep` is present, split the step list
        # at it and compile the pre-fork (per-clone) simulator plus an
        # optional post-fork (per-read) simulator, mirroring the
        # `_ClonalForkStep` branch. Per-clone sizes are drawn at run
        # time from the heavy-tailed distribution; identical reads are
        # collapsed into `duplicate_count`-carrying records.
        repertoire_idx = next(
            (
                i
                for i, s in enumerate(self._steps)
                if isinstance(s, _RepertoireForkStep)
            ),
            None,
        )
        if repertoire_idx is not None:
            repertoire_step: _RepertoireForkStep = self._steps[repertoire_idx]
            pre_steps = self._steps[:repertoire_idx]
            post_steps = self._steps[repertoire_idx + 1 :]
            pre_simulator = self._build_simulator(
                pre_steps,
                contracts,
                any_lock,
                replace_fn=_replace,
                allow_curatable_refdata=allow_curatable_refdata,
            )
            # post_steps may be empty — that's the pure-copy case
            # (each clone collapses to one record with
            # duplicate_count = size). When present, the post-fork
            # plan inherits the parent's V/D/J/NP backbone, so
            # any_lock=False and no recombination facts on this side
            # (same as the clonal branch).
            post_simulator = None
            if post_steps:
                post_simulator = self._build_simulator(
                    post_steps,
                    contracts,
                    any_lock=False,
                    replace_fn=_replace,
                    allow_curatable_refdata=allow_curatable_refdata,
                )
            return CompiledRepertoireExperiment(
                pre_simulator,
                post_simulator,
                repertoire_step,
                self._refdata,
                dataconfig=self._dataconfig,
                metadata=self._metadata,
            )

        # if a `_LineageForkStep` is present, compile a
        # CompiledLineageExperiment: pre-fork steps (recombine) become
        # the founder simulator; post-fork steps must be empty (the
        # lineage engine handles mutation internally).
        lineage_idx = next(
            (i for i, s in enumerate(self._steps) if isinstance(s, _LineageForkStep)),
            None,
        )
        if lineage_idx is not None:
            lineage_step: _LineageForkStep = self._steps[lineage_idx]
            pre_steps = self._steps[:lineage_idx]
            post_steps = self._steps[lineage_idx + 1:]
            # Steps after clonal_lineage() are per-observed-cell
            # library-prep / sequencing artefact passes (the same
            # post-fork set expand_clones() allows). SHM is internal
            # to the lineage engine, so .mutate() is rejected; the
            # paired-end read layout is not yet wired through the
            # per-cell corruption merge, so reject it for now too.
            for s in post_steps:
                if isinstance(s, _MutateStep):
                    raise ValueError(
                        "SHM is internal to clonal_lineage; do not add "
                        ".mutate() after it. Set the within-lineage SHM rate "
                        "via clonal_lineage(rate=...)."
                    )
                if isinstance(s, _PairedEndStep):
                    raise ValueError(
                        "paired_end not yet supported with clonal_lineage; "
                        "apply per-read library-prep passes (sequencing_errors, "
                        "pcr_amplify, polymerase_indels, end_loss_*, "
                        "ambiguous_base_calls, random_strand_orientation) instead."
                    )
                if not isinstance(s, _CorruptStep):
                    raise ValueError(
                        "Only per-read library-prep / sequencing artefact passes "
                        "may follow clonal_lineage(); got "
                        f"{type(s).__name__}."
                    )
            pre_simulator = self._build_simulator(
                pre_steps,
                contracts,
                any_lock,
                replace_fn=_replace,
                allow_curatable_refdata=allow_curatable_refdata,
            )
            # Build the per-cell corruption simulator from the
            # post-fork steps, mirroring the clonal branch's
            # post_simulator (no recombination facts on this side, so
            # any_lock=False and the analyzer skips the productive
            # precondition check).
            post_simulator = None
            if post_steps:
                post_simulator = self._build_simulator(
                    post_steps,
                    contracts,
                    any_lock=False,
                    replace_fn=_replace,
                    allow_curatable_refdata=allow_curatable_refdata,
                )
            return CompiledLineageExperiment(
                pre_simulator,
                lineage_step,
                self._refdata,
                post_simulator=post_simulator,
                post_steps=tuple(post_steps),
                dataconfig=self._dataconfig,
                metadata=self._metadata,
            )

        # When the attached genotype defines novel/private alleles, compile
        # against an *effective* reference = base catalogue + injected novel
        # alleles, so they become real pool entries the engine samples,
        # assembles, and reports like any allele. No genotype, or a genotype
        # without novel alleles, uses the base refdata unchanged.
        effective_refdata = self._refdata
        if self._genotype is not None and self._genotype.has_novel():
            if self._dataconfig is None:
                raise ValueError(
                    "genotype with novel alleles requires a DataConfig-backed "
                    "experiment (Experiment.on(dataconfig), not a raw RefDataConfig)"
                )
            effective_refdata = dataconfig_to_refdata(
                self._genotype.effective_dataconfig()
            )

        simulator = self._build_simulator(
            self._steps,
            contracts,
            any_lock,
            replace_fn=_replace,
            allow_curatable_refdata=allow_curatable_refdata,
            refdata=effective_refdata,
        )
        return CompiledExperiment(
            simulator,
            effective_refdata,
            steps=tuple(self._steps),
            dataconfig=self._dataconfig,
            metadata=self._metadata,
            genotype=self._genotype,
        )

    def _build_simulator(
        self,
        steps,
        contracts,
        any_lock: bool,
        *,
        replace_fn,
        allow_curatable_refdata: bool = False,
        refdata=None,
    ):
        """Compile a list of steps into a `GenAIRR._engine.CompiledSimulator`.
        Lifted out of `compile()` so the clonal-fork branch can build
        two simulators from sub-step-lists with a shared body.

        ``refdata`` overrides ``self._refdata`` — used when a genotype with
        novel alleles compiles against an *effective* reference (base +
        injected private alleles)."""
        refdata = refdata if refdata is not None else self._refdata
        plan = _engine.PassPlan()
        # Pull the (at-most-one) `_InvertDStep` out of the step
        # sequence and thread its probability into the recombine
        # lowering directly. Inlining the InvertDPass push between
        # `push_generate_np("NP1", ...)` and `push_assemble("D")`
        # (see `_lower_recombine`) keeps the canonical V-NP1-D-NP2-J
        # pool layout intact; a separate `_InvertDStep` lowering
        # combined with `Schedule::before(invert_d, assemble.d)`
        # would re-promote the pass past `assemble.j`, swapping D
        # and J in the pool. See `_extract_invert_d_prob` for the
        # full rationale.
        invert_d_prob, steps = _extract_invert_d_prob(steps)
        # Same inline-into-recombine-lowering pattern as
        # _extract_invert_d_prob: pulling the receptor-revision step
        # out here lets the lowering push it at the exact slot
        # between `assemble.j` and any subsequent mutate/corrupt
        # passes. A standalone lower path would either need a new
        # schedule edge or place the pass at the end of the plan
        # (after corruption), both of which break the design doc §2
        # ordering.
        receptor_revision_prob, receptor_revision_same_haplotype, steps = (
            _extract_receptor_revision_prob(steps)
        )
        # Pull out the (at-most-one) paired-end step too. It must
        # land at the END of the plan, not inline with recombine —
        # see `_extract_paired_end_step` for the rationale on
        # placement vs. trace order.
        paired_end_step, steps = _extract_paired_end_step(steps)
        for step in steps:
            # Inject any allele-locks set via ``.restrict_alleles(...)`` into the
            # recombine step at compile time. Other step types ignore
            # locks.
            if any_lock and isinstance(step, _RecombineStep):
                step = replace_fn(
                    step,
                    locks_v=self._locks["V"],
                    locks_d=self._locks["D"],
                    locks_j=self._locks["J"],
                )
            if isinstance(step, _RecombineStep):
                _lower_recombine(
                    step,
                    plan,
                    refdata,
                    invert_d_prob=invert_d_prob,
                    receptor_revision_prob=receptor_revision_prob,
                    receptor_revision_same_haplotype=receptor_revision_same_haplotype,
                    genotype=self._genotype,
                )
            else:
                lower_step(step, plan, refdata)
        # Paired-end is sequencing-stage / readout-stage: lower
        # it AFTER every biology + corruption pass so the trace
        # records land last. See `_extract_paired_end_step` for
        # the full rationale.
        if paired_end_step is not None:
            _lower_paired_end(paired_end_step, plan)
        return plan.compile(
            refdata=refdata,
            respect=contracts,
            allow_curatable_refdata=allow_curatable_refdata,
        )
