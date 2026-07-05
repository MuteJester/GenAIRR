"""Runtime wrapper for the non-tree clonal-repertoire experiment.

:class:`CompiledRepertoireExperiment` draws a heavy-tailed size per clone,
emits that many reads through the post-fork passes, and genotype-collapses
identical reads into AIRR records carrying a standard ``duplicate_count``.
It is intentionally thin: it holds references to the compiled engine
simulators and source context. The execution loop delegates straight into
the engine; the wrapper exists for ergonomics and to keep the public
surface stable across engine refactors.
"""
from __future__ import annotations

from typing import Any, Dict, List, Optional, TYPE_CHECKING

from GenAIRR import _engine  # private Rust extension submodule

from .._pipeline_ir import _RepertoireForkStep

if TYPE_CHECKING:
    from ..dataconfig import DataConfig
    from ..result import SimulationResult


class CompiledRepertoireExperiment:
    """A compiled non-tree clonal-repertoire experiment.

    Wraps a pre-fork :class:`GenAIRR._engine.CompiledSimulator` (the
    founder recombination, run once per clone) and an optional
    post-fork simulator (the per-read library-prep / sequencing passes).
    Per clone a size is drawn from a heavy-tailed distribution via
    ``_engine.sample_clone_sizes``; that many reads are emitted through
    the post-fork passes and identical reads are genotype-collapsed into
    AIRR records carrying a standard ``duplicate_count``.

    When there are no post-fork passes every read of a clone is
    identical, so the clone collapses to a single record whose
    ``duplicate_count`` equals the drawn size (the no-corruption
    shortcut). Every record carries an integer ``clone_id``.
    """

    __slots__ = (
        "_pre",
        "_post",
        "_step",
        "_refdata",
        "_dataconfig",
        "_metadata",
    )

    def __init__(
        self,
        pre_simulator: "_engine.CompiledSimulator",
        post_simulator: Optional["_engine.CompiledSimulator"],
        step: "_RepertoireForkStep",
        refdata: "_engine.RefDataConfig",
        dataconfig: Optional["DataConfig"] = None,
        metadata: Optional[Dict[str, Any]] = None,
    ) -> None:
        self._pre = pre_simulator
        self._post = post_simulator
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
        validate_records: bool = False,
        expose_provenance: bool = False,
    ) -> "SimulationResult":
        """Draw per-clone sizes, emit reads through the post-fork passes,
        collapse identical reads, and return a :class:`SimulationResult`
        whose record dicts carry ``clone_id`` and ``duplicate_count``.
        """
        import GenAIRR._engine as _engine

        from .._airr_record import outcome_to_airr_record
        from ..result import SimulationResult, _inject_truth_columns

        step = self._step
        sizes = _engine.sample_clone_sizes(
            step.n_clones,
            int(seed),
            kind=step.size_distribution,
            exponent=step.exponent,
            mu=step.mu,
            sigma=step.sigma,
            max_size=step.max_size,
            unexpanded_fraction=step.unexpanded_fraction,
        )
        records: List[Dict[str, Any]] = []
        outcomes: List["_engine.Outcome"] = []
        for clone_idx, size in enumerate(sizes):
            clone_seed = int(seed) + clone_idx * 1_000_000
            founder = self._pre.run(seed=clone_seed, strict=strict)
            founder_sim = founder.final_simulation()
            if self._post is None:
                # No post-fork passes: all ``size`` copies are
                # identical -> one record with duplicate_count = size.
                rec = outcome_to_airr_record(
                    founder,
                    self._refdata,
                    sequence_id=f"clone{clone_idx}_read0",
                )
                rec["clone_id"] = clone_idx
                rec["duplicate_count"] = int(size)
                if expose_provenance:
                    _inject_truth_columns(founder, self._refdata, rec)
                records.append(rec)
                outcomes.append(founder)
                continue
            # Post-fork passes present: simulate ``size`` reads and
            # collapse by emitted sequence.
            by_seq: Dict[str, list] = {}
            for read_idx in range(int(size)):
                desc = self._post.run_from(
                    founder_sim, clone_seed + 1 + read_idx, strict=strict
                )
                rec = outcome_to_airr_record(
                    desc,
                    self._refdata,
                    sequence_id=f"clone{clone_idx}_read{read_idx}",
                )
                key = rec["sequence"]
                if key in by_seq:
                    by_seq[key][0] += 1
                else:
                    by_seq[key] = [1, desc, rec]
            for count, desc, rec in by_seq.values():
                rec["clone_id"] = clone_idx
                rec["duplicate_count"] = count
                if expose_provenance:
                    _inject_truth_columns(desc, self._refdata, rec)
                records.append(rec)
                outcomes.append(desc)
        result = SimulationResult(records, outcomes=outcomes)
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
            f"<CompiledRepertoireExperiment n_clones={self._step.n_clones} "
            f"size_distribution={self._step.size_distribution!r} "
            f"chain={self._refdata.chain_type}>"
        )
