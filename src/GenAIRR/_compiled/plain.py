"""Runtime wrapper for the plain (non-forked) compiled experiment.

:class:`CompiledExperiment` is the frozen, executable form of
:class:`GenAIRR.experiment.Experiment`. It is intentionally thin: it
holds a reference to a compiled engine simulator and a few pieces of
source context (``steps``, ``dataconfig``, ``metadata``) so
``describe()`` and AIRR record export can render a faithful narrative.
The execution loops delegate straight into the engine; the wrapper
exists for ergonomics and to keep the public surface stable across
engine refactors.
"""
from __future__ import annotations

from typing import Any, Dict, Iterator, List, Optional, Sequence, Tuple, TYPE_CHECKING

from GenAIRR import _engine  # private Rust extension submodule

from .._describe import (
    _describe_experiment_header,
    _describe_step_sequence,
    _format_active_contracts,
)

if TYPE_CHECKING:
    from ..dataconfig import DataConfig
    from ..result import SimulationResult


class CompiledExperiment:
    """A frozen ``Experiment`` ready for execution.

    Holds the owning :class:`GenAIRR._engine.CompiledSimulator` and the
    refdata it was built against. Contracts are captured at compile
    time; ``run()`` only accepts execution parameters.
    """

    __slots__ = (
        "_simulator",
        "_refdata",
        "_steps",
        "_dataconfig",
        "_metadata",
        "_genotype",
    )

    def __init__(
        self,
        simulator: "_engine.CompiledSimulator",
        refdata: "_engine.RefDataConfig",
        steps: Sequence[Any] = (),
        dataconfig: Optional["DataConfig"] = None,
        metadata: Optional[Dict[str, Any]] = None,
        genotype: Optional[Any] = None,
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
        # Attached single-subject genotype (or None). When set,
        # run_records stamps subject_id + haplotype provenance.
        self._genotype = genotype

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
        from ..result import SimulationResult

        outcomes = self.run(n=n, seed=seed, strict=strict)
        result = SimulationResult.from_outcomes(
            outcomes, self._refdata, expose_provenance=expose_provenance
        )
        if self._genotype is not None:
            self._stamp_genotype_provenance(outcomes, result)
        if validate_records:
            from .._validation import _raise_on_validation_failure

            _raise_on_validation_failure(result.validate_records(self._refdata))
        return result

    def _stamp_genotype_provenance(self, outcomes, result) -> None:
        """Stamp per-record ``subject_id`` + ``haplotype`` (the chromosome
        the rearrangement drew from) and expose the genotype on the
        result. Used only when a genotype is attached."""
        subject = self._genotype.subject_id
        for outcome, rec in zip(outcomes, result._records):
            rec["subject_id"] = subject
            hap = outcome.trace().find("sample_haplotype")
            rec["haplotype"] = hap.value if hap is not None else None
        result._genotypes = [self._genotype]

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
        from .._airr_record import outcome_to_airr_record

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
