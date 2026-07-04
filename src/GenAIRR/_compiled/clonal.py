"""Runtime wrapper for the clonal-fork compiled experiment.

:class:`CompiledClonalExperiment` wraps the two-stage
parent-then-descendants flow produced by ``expand_clones(...)``. It
holds two compiled engine simulators (pre-fork and post-fork) plus
source context so ``describe()`` and AIRR record export can render a
faithful narrative. The execution loop delegates straight into the
engine; the wrapper exists for ergonomics and to keep the public
surface stable across engine refactors.
"""
from __future__ import annotations

from typing import Any, Dict, List, Optional, Sequence, Tuple, TYPE_CHECKING

from GenAIRR import _engine  # private Rust extension submodule

from .._describe import (
    _describe_clonal_fork_step,
    _describe_experiment_header,
    _describe_step_sequence,
)
from .._pipeline_ir import _ClonalForkStep

if TYPE_CHECKING:
    from ..dataconfig import DataConfig
    from ..result import SimulationResult


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
        from .._airr_record import outcome_to_airr_record
        from ..result import SimulationResult, _inject_truth_columns

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
            from .._validation import (
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
