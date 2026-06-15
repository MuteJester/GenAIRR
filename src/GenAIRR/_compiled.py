"""Runtime wrappers around :class:`GenAIRR._engine.CompiledSimulator`.

:class:`CompiledExperiment` is the frozen, executable form of
:class:`GenAIRR.experiment.Experiment`. :class:`CompiledClonalExperiment`
wraps the two-stage parent-then-descendants flow produced by
``expand_clones(...)``.

Both classes are intentionally thin: they hold a reference to a
compiled engine simulator and a few pieces of source context
(``steps``, ``dataconfig``, ``metadata``) so ``describe()`` and AIRR
record export can render a faithful narrative. The execution loops
delegate straight into the engine; the wrappers exist for ergonomics
and to keep the public surface stable across engine refactors.
"""
from __future__ import annotations

from typing import Any, Dict, Iterator, List, Optional, Sequence, Tuple, TYPE_CHECKING

from GenAIRR import _engine  # private Rust extension submodule

from ._describe import (
    _describe_clonal_fork_step,
    _describe_experiment_header,
    _describe_step_sequence,
    _format_active_contracts,
)
from ._pipeline_ir import _ClonalForkStep, _LineageForkStep

if TYPE_CHECKING:
    from .dataconfig import DataConfig
    from .result import SimulationResult


class CompiledExperiment:
    """A frozen ``Experiment`` ready for execution.

    Holds the owning :class:`GenAIRR._engine.CompiledSimulator` and the
    refdata it was built against. Contracts are captured at compile
    time; ``run()`` only accepts execution parameters.
    """

    __slots__ = ("_simulator", "_refdata", "_steps", "_dataconfig", "_metadata")

    def __init__(
        self,
        simulator: "_engine.CompiledSimulator",
        refdata: "_engine.RefDataConfig",
        steps: Sequence[Any] = (),
        dataconfig: Optional["DataConfig"] = None,
        metadata: Optional[Dict[str, Any]] = None,
    ) -> None:
        self._simulator = simulator
        self._refdata = refdata
        # Source steps stashed for `describe()`. The compiled simulator
        # itself can render pass names but not biology — keeping the
        # builder steps around is the cheapest way to give a faithful
        # narrative back.
        self._steps: Tuple[Any, ...] = tuple(steps)
        self._dataconfig = dataconfig
        self._metadata = dict(metadata) if metadata else {}

    @property
    def simulator(self) -> "_engine.CompiledSimulator":
        """The owning Rust compiled simulator."""
        return self._simulator

    @property
    def pass_plan(self) -> Tuple[str, ...]:
        """Read-only pass-name summary of the compiled pipeline."""
        return tuple(self._simulator.pass_names())

    @property
    def pass_names(self) -> Tuple[str, ...]:
        """Stable names of the compiled pass sequence."""
        return self.pass_plan

    @property
    def active_contracts(self) -> Tuple[str, ...]:
        """Stable names of the contract bundle captured at compile time."""
        return tuple(self._simulator.active_contracts())

    @property
    def refdata(self) -> "_engine.RefDataConfig":
        """The :class:`GenAIRR._engine.RefDataConfig` the plan was built against."""
        return self._refdata

    def describe(self) -> str:
        """Render a biology-style narrative of the compiled experiment.

        Equivalent to ``Experiment.describe()`` but additionally
        surfaces compile-time constraints (e.g. ``productive_only``)
        attached via ``compile(respect=...)``. See
        :meth:`Experiment.describe` for the output shape.
        """
        header = _describe_experiment_header(self._refdata, self._dataconfig)
        if not self._steps:
            body = ["  (no steps recorded)"]
        else:
            body = _describe_step_sequence(self._steps, self._refdata.chain_type)
        lines = [header, *body]
        contracts_line = _format_active_contracts(self.active_contracts)
        if contracts_line:
            lines.append(f"  Constraints: {contracts_line}")
        if self._metadata:
            stamps = ", ".join(f"{k}={v!r}" for k, v in self._metadata.items())
            lines.append(f"  Metadata stamped on every record: {stamps}")
        return "\n".join(lines)

    def run(
        self,
        *,
        n: int = 1,
        seed: int = 0,
        strict: bool = False,
    ) -> List["_engine.Outcome"]:
        """Run the compiled simulator ``n`` times — **fresh sampling**.

        Each iteration uses ``seed + i`` as the per-run seed so
        consecutive batches stitch together by offsetting ``seed``.

        ``strict`` controls the failure mode when a pass's
        contract-narrowed candidate set is empty *at sample time*:

        - ``False`` (default, **permissive**) — apply the pass's
          declared ``EmptySupport`` policy. The pass writes a
          documented sentinel value to the trace and continues.
          Common sentinels: indel ``site = -1`` NoOp, NP length ``0``,
          NP base ``N``, trim ``0``; SHM substitution skips the slot
          (no trace record).
        - ``True`` (**strict**) — raise
          :class:`GenAIRR._engine.StrictSamplingError`. The exception's
          ``args`` are a 3-tuple ``(pass_name, address, reason)``;
          ``reason`` is one of ``"support_unavailable"``,
          ``"empty_admissible_support"``, or
          ``"invalid_filtered_support"``.

        **Compile-time precondition failures are separate.** If a
        sampling distribution is *statically* incompatible with the
        active contracts (e.g. every NP1 length in the distribution
        violates frame divisibility), :meth:`Experiment.compile`
        raises :class:`ValueError` *before* this method runs.
        ``ValueError`` and ``StrictSamplingError`` have **no shared
        base class** — catching only one will miss the other. See
        ``docs/productive_failure_mode_audit.md`` §6.1.

        **Strict semantics apply only to fresh sampling.** Trace
        replay via :meth:`CompiledExperiment.replay_from_trace_file`
        consumes the recorded values verbatim and does NOT
        re-evaluate contract admissibility — a permissive sentinel
        trace replays cleanly even under ``strict=True``. See that
        method's docstring.

        Raises ``ValueError`` for ``n < 1``.
        """
        if n < 1:
            raise ValueError(f"n must be at least 1, got {n}")
        return self._simulator.run_batch(n, seed, strict=strict)

    def run_records(
        self,
        *,
        n: int = 1,
        seed: int = 0,
        strict: bool = False,
        expose_provenance: bool = False,
        validate_records: bool = False,
    ) -> "SimulationResult":
        """Run the compiled simulator ``n`` times and return the batch as
        a :class:`SimulationResult` ready for ``.to_csv`` /
        ``.to_fasta`` / ``.to_dataframe`` export.

        Same arguments as :meth:`run`. ``expose_provenance=True``
        appends `truth_v_call/d_call/j_call` columns reflecting the
        originally-sampled allele names.

        ``validate_records=True`` runs
        :meth:`SimulationResult.validate_records` on the freshly
        built batch and raises
        :class:`GenAIRR._validation.RecordValidationFailedError`
        (a :class:`RuntimeError` subclass) when any record fails the
        postcondition validator. Default ``False`` keeps this method
        zero-overhead.
        """
        from .result import SimulationResult

        outcomes = self.run(n=n, seed=seed, strict=strict)
        result = SimulationResult.from_outcomes(
            outcomes, self._refdata, expose_provenance=expose_provenance
        )
        if validate_records:
            from ._validation import _raise_on_validation_failure

            _raise_on_validation_failure(result.validate_records(self._refdata))
        return result

    def stream(
        self,
        *,
        n: Optional[int] = None,
        seed: int = 0,
        strict: bool = False,
    ) -> Iterator["_engine.Outcome"]:
        """Lazily yield :class:`GenAIRR._engine.Outcome` objects one at
        a time, without materialising the full batch in memory.

        Useful for large simulations where holding ``n`` outcomes
        would be wasteful — typical pattern is

        >>> for outcome in compiled.stream(n=1_000_000, seed=0):
        ...     write_to_disk(outcome)

        ``n=None`` (the default) yields outcomes indefinitely; the
        caller is expected to stop with ``itertools.islice``,
        ``break``, or similar. ``n=N`` yields exactly ``N`` outcomes
        with seeds ``seed`` … ``seed + N - 1``.

        ``strict`` behaves as in :meth:`run`.

        Raises ``ValueError`` when ``n`` is set to a value below 1.
        """
        if n is not None and n < 1:
            raise ValueError(f"n must be at least 1, got {n}")
        i = 0
        while n is None or i < n:
            yield self._simulator.run(seed + i, strict=strict)
            i += 1

    def stream_records(
        self,
        *,
        n: Optional[int] = None,
        seed: int = 0,
        strict: bool = False,
        id_prefix: str = "seq",
    ) -> Iterator[Dict[str, Any]]:
        """Lazily yield AIRR-format record dicts (one per outcome).

        Same shape as the records inside a :class:`SimulationResult`,
        but yielded one at a time so callers can write each record to
        disk without retaining the prior ones. Pairs naturally with
        :func:`csv.DictWriter` for streaming TSV/CSV output.

        Each record's ``sequence_id`` is set to
        ``f"{id_prefix}{i}"`` so streamed batches have unique
        AIRR-style identifiers without buffering.
        """
        from ._airr_record import outcome_to_airr_record

        for i, outcome in enumerate(
            self.stream(n=n, seed=seed, strict=strict)
        ):
            yield outcome_to_airr_record(
                outcome, self._refdata, sequence_id=f"{id_prefix}{i}"
            )

    def __repr__(self) -> str:
        return (
            f"<CompiledExperiment plan_len={len(self.pass_plan)} "
            f"chain={self._refdata.chain_type} "
            f"contracts={len(self.active_contracts)}>"
        )


class CompiledClonalExperiment:
    """A compiled experiment with a clonal-fork structure.

    Wraps two :class:`GenAIRR._engine.CompiledSimulator`s — the
    pre-fork plan (run once per clone, typically the recombine
    step) and the post-fork plan (run once per descendant inside
    the clone, typically mutate / corrupt_*).

    :meth:`run_records` orchestrates the parent → descendants loop
    and tags every record with a ``clone_id`` integer so downstream
    clonotype-clustering tools can be benchmarked against the true
    clonal structure.
    """

    __slots__ = (
        "_pre",
        "_post",
        "_fork",
        "_refdata",
        "_pre_steps",
        "_post_steps",
        "_dataconfig",
        "_metadata",
    )

    def __init__(
        self,
        pre_simulator: "_engine.CompiledSimulator",
        post_simulator: "_engine.CompiledSimulator",
        fork: "_ClonalForkStep",
        refdata: "_engine.RefDataConfig",
        pre_steps: Sequence[Any] = (),
        post_steps: Sequence[Any] = (),
        dataconfig: Optional["DataConfig"] = None,
        metadata: Optional[Dict[str, Any]] = None,
    ) -> None:
        self._pre = pre_simulator
        self._post = post_simulator
        self._fork = fork
        self._refdata = refdata
        self._pre_steps: Tuple[Any, ...] = tuple(pre_steps)
        self._post_steps: Tuple[Any, ...] = tuple(post_steps)
        self._dataconfig = dataconfig
        self._metadata = dict(metadata) if metadata else {}

    @property
    def n_clones(self) -> int:
        return self._fork.n_clones

    @property
    def size(self) -> int:
        return self._fork.size

    @property
    def total_records(self) -> int:
        """Number of records produced per :meth:`run_records` call."""
        return self._fork.n_clones * self._fork.size

    @property
    def refdata(self) -> "_engine.RefDataConfig":
        return self._refdata

    def describe(self) -> str:
        """Render a biology-style narrative of the compiled clonal
        experiment, with an explicit divider at the fork. See
        :meth:`Experiment.describe` for the basic shape."""
        header = _describe_experiment_header(self._refdata, self._dataconfig)
        lines = [header]
        # pre-fork section (per-clone)
        pre_lines = _describe_step_sequence(
            self._pre_steps, self._refdata.chain_type, start_index=1
        )
        lines.extend(pre_lines)
        # the fork itself
        lines.append(f"  ── {_describe_clonal_fork_step(self._fork)} ──")
        lines.append(
            "      (steps above run once per clone; "
            "steps below run once per descendant)"
        )
        # post-fork section (per-descendant)
        post_start = sum(
            1 for s in self._pre_steps if not isinstance(s, _ClonalForkStep)
        ) + 1
        post_lines = _describe_step_sequence(
            self._post_steps, self._refdata.chain_type, start_index=post_start
        )
        lines.extend(post_lines)
        if self._metadata:
            stamps = ", ".join(f"{k}={v!r}" for k, v in self._metadata.items())
            lines.append(f"  Metadata stamped on every record: {stamps}")
        return "\n".join(lines)

    def run(
        self,
        *,
        n: Optional[int] = None,
        seed: int = 0,
        strict: bool = False,
    ) -> List["_engine.Outcome"]:
        """Run all clonal descendants and return their outcomes in
        clone-major order (clone 0's descendants 0..size-1, clone 1's
        descendants 0..size-1, …).

        ``n`` is optional: when omitted the runtime expands
        ``n_clones * size`` outcomes. Passing ``n`` is allowed only
        when ``n == n_clones * size`` (otherwise raises).
        """
        total = self.total_records
        if n is not None and n != total:
            raise ValueError(
                f"clonal pipeline produces n_clones * size = "
                f"{self._fork.n_clones} * {self._fork.size} = {total} "
                f"records; passing n={n} is inconsistent. Drop the n "
                f"argument or pass n={total}."
            )

        outcomes: List["_engine.Outcome"] = []
        for clone_idx in range(self._fork.n_clones):
            clone_seed = int(seed) + clone_idx * 1_000_000
            parent = self._pre.run(seed=clone_seed, strict=strict)
            parent_sim = parent.final_simulation()
            for desc_idx in range(self._fork.size):
                desc_seed = clone_seed + 1 + desc_idx
                desc = self._post.run_from(
                    parent_sim, desc_seed, strict=strict
                )
                outcomes.append(desc)
        return outcomes

    def run_records(
        self,
        *,
        n: Optional[int] = None,
        seed: int = 0,
        strict: bool = False,
        expose_provenance: bool = False,
        validate_records: bool = False,
    ) -> "SimulationResult":
        """Same as :meth:`run` but returns a :class:`SimulationResult`
        with each record dict carrying an integer ``clone_id`` field
        in ``[0, n_clones)`` plus a ``parent_id`` integer indexing
        into :attr:`SimulationResult.parents`. ``expose_provenance=True``
        also appends `truth_v_call` / `truth_d_call` / `truth_j_call`
        columns from the originally-sampled allele names.

        The returned :class:`SimulationResult` carries the per-clone
        parent ``Outcome`` objects on its ``.parents`` attribute.
        Each parent holds the pre-fork addressed-choice trace, the
        pre-fork event ledger, and the post-recombination IR — useful
        for replay, lineage analysis, and the upcoming parent-aware
        family validator. The flat ``.outcomes`` list continues to
        carry only the descendant outcomes (one per record);
        parents are exposed separately so the per-record list stays
        the same shape clonal consumers already know.

        ``validate_records=True`` runs
        :meth:`SimulationResult.validate_records` on the freshly
        built batch and raises
        :class:`GenAIRR._validation.RecordValidationFailedError`
        on any postcondition failure. After the per-record gate
        passes, this also runs
        :meth:`SimulationResult.validate_families` and raises the
        sibling
        :class:`GenAIRR._validation.FamilyValidationFailedError`
        if any clonal-family invariant is violated. The two
        gates report separately so users can tell projection bugs
        from family-consistency bugs. Default ``False`` keeps this
        method zero-overhead.
        """
        from ._airr_record import outcome_to_airr_record
        from .result import SimulationResult, _inject_truth_columns

        total = self.total_records
        if n is not None and n != total:
            raise ValueError(
                f"clonal pipeline produces n_clones * size = "
                f"{self._fork.n_clones} * {self._fork.size} = {total} "
                f"records; passing n={n} is inconsistent. Drop the n "
                f"argument or pass n={total}."
            )

        records: List[Dict[str, Any]] = []
        outcomes: List["_engine.Outcome"] = []
        # Slice 2: retain parent outcomes — one per clone. The
        # orchestration loop used to drop them after extracting
        # ``final_simulation()``. We now keep them so the returned
        # :class:`SimulationResult` can expose ``.parents`` for
        # replay / lineage tooling. The parent ``Outcome`` itself
        # is not copied onto each descendant; we hold a single
        # reference per clone in the ``parents`` list.
        parents: List["_engine.Outcome"] = []
        for clone_idx in range(self._fork.n_clones):
            clone_seed = int(seed) + clone_idx * 1_000_000
            parent = self._pre.run(seed=clone_seed, strict=strict)
            parents.append(parent)
            parent_sim = parent.final_simulation()
            for desc_idx in range(self._fork.size):
                desc_seed = clone_seed + 1 + desc_idx
                desc = self._post.run_from(
                    parent_sim, desc_seed, strict=strict
                )
                outcomes.append(desc)
                rec = outcome_to_airr_record(
                    desc,
                    self._refdata,
                    sequence_id=f"clone{clone_idx}_desc{desc_idx}",
                )
                rec["clone_id"] = clone_idx
                # ``parent_id`` is the descendant's index into
                # ``result.parents``. Today clones are dense and
                # zero-based, so ``parent_id == clone_id`` by
                # construction — we stamp both because they carry
                # distinct semantics: ``clone_id`` is the family
                # identity (the existing Slice 0 contract);
                # ``parent_id`` is the addressing scheme into the
                # parent-outcome list (the new Slice 2 contract).
                # Keeping them separate now means a future slice
                # that introduces sparse / non-zero-based family
                # ids doesn't have to retrofit both.
                rec["parent_id"] = clone_idx
                if expose_provenance:
                    _inject_truth_columns(desc, self._refdata, rec)
                records.append(rec)
        result = SimulationResult(records, outcomes=outcomes, parents=parents)
        if validate_records:
            from ._validation import (
                _raise_on_family_validation_failure,
                _raise_on_validation_failure,
            )

            _raise_on_validation_failure(result.validate_records(self._refdata))
            _raise_on_family_validation_failure(result.validate_families())
        return result

    def __repr__(self) -> str:
        return (
            f"<CompiledClonalExperiment n_clones={self._fork.n_clones} "
            f"size={self._fork.size} chain={self._refdata.chain_type}>"
        )


class CompiledLineageExperiment:
    """A compiled experiment that grows BCR affinity-maturation lineage trees.

    Wraps a pre-fork :class:`GenAIRR._engine.CompiledSimulator` (founder
    recombination) and a :class:`~GenAIRR._pipeline_ir._LineageForkStep`
    (lineage parameters). :meth:`run_records` grows one lineage tree per
    clone via the Rust ``simulate_family_outcomes`` kernel and returns a
    :class:`~GenAIRR.result.SimulationResultWithLineages` whose records are
    per-observed-node AIRR dicts with lineage metadata.
    """

    __slots__ = ("_pre", "_step", "_refdata", "_dataconfig", "_metadata")

    def __init__(
        self,
        pre_simulator: "_engine.CompiledSimulator",
        step: "_LineageForkStep",
        refdata: "_engine.RefDataConfig",
        dataconfig: Optional["DataConfig"] = None,
        metadata: Optional[Dict[str, Any]] = None,
    ) -> None:
        self._pre = pre_simulator
        self._step = step
        self._refdata = refdata
        self._dataconfig = dataconfig
        self._metadata = dict(metadata) if metadata else {}

    @property
    def refdata(self) -> "_engine.RefDataConfig":
        return self._refdata

    def run_records(
        self,
        *,
        seed: int = 0,
        strict: bool = False,
    ) -> "SimulationResultWithLineages":
        """Grow lineage trees and return per-observed-node AIRR records.

        Each clone is seeded at ``seed + clone_idx * 1_000_000`` so
        independent clones are reproducible and non-overlapping.
        """
        from GenAIRR import _engine as _eng
        from ._s5f_loader import load_builtin_s5f_kernel
        from .result import SimulationResultWithLineages

        step = self._step
        mutability, substitution = load_builtin_s5f_kernel(step.s5f_model)

        records: List[Dict[str, Any]] = []
        lineage_trees = []

        for clone_idx in range(step.n_clones):
            clone_seed = int(seed) + clone_idx * 1_000_000
            # _pre.run() returns a single Outcome (not a list) — mirror
            # CompiledClonalExperiment which calls self._pre.run(seed=...).
            founder = self._pre.run(seed=clone_seed, strict=strict)
            fam = _eng.simulate_family_outcomes(
                founder,
                mutability,
                substitution,
                step.rate,
                step.lambda_base,
                step.lambda_mut,
                step.max_generations,
                step.n_max,
                step.n_sample,
                clone_seed,
                selection_strength=step.selection_strength,
                beta=step.beta,
                target_aa=step.target_aa,
                mature_substitutions=step.mature_substitutions,
            )
            tree = fam.tree()
            lineage_trees.append(tree)

            # Collect observed nodes in id order (same order as airr_records).
            observed_nodes = [n for n in tree.nodes() if n.observed]
            recs = fam.airr_records(self._refdata)

            for node, rec in zip(observed_nodes, recs):
                rec["clone_id"] = clone_idx
                rec["lineage_node_id"] = node.id
                rec["lineage_parent_id"] = (
                    node.parent_id if node.parent_id is not None else -1
                )
                rec["lineage_generation"] = node.generation
                rec["lineage_abundance"] = node.abundance
                rec["lineage_affinity"] = node.affinity
                rec["sequence_id"] = f"clone{clone_idx}_node{node.id}"
                records.append(rec)

        return SimulationResultWithLineages(records, lineage_trees=lineage_trees)

    def __repr__(self) -> str:
        return (
            f"<CompiledLineageExperiment n_clones={self._step.n_clones} "
            f"max_gen={self._step.max_generations} chain={self._refdata.chain_type}>"
        )
