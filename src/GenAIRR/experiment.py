"""GenAIRR — fluent DSL for receptor-sequence simulation.

The :class:`Experiment` builder lowers a chain of fluent steps into a
Rust ``PassPlan`` (via :mod:`GenAIRR._engine`), compiles that IR into an
owning ``CompiledSimulator``, and runs it to produce a list of
:class:`GenAIRR._engine.Outcome` objects.

Typical usage::

    import GenAIRR as ga
    outcomes = ga.Experiment.on("human_igh").recombine().run(n=100, seed=42)

``Experiment.on`` accepts:

- a config-name string (``"human_igh"``, ``"mouse_tcrb"``, …) — resolved
  through the builtin :mod:`GenAIRR.data` registry,
- a :class:`GenAIRR.DataConfig` object (already-loaded reference data),
- a :class:`GenAIRR._engine.RefDataConfig` (the engine-native form, used
  primarily by tests and advanced callers).

In all three cases the input is normalised to a
:class:`GenAIRR._engine.RefDataConfig` before any pass is appended.
"""
from __future__ import annotations

import warnings
from typing import Any, Dict, Iterable, Iterator, List, Optional, Tuple, Union

from GenAIRR import _engine  # private Rust extension submodule

from .dataconfig import DataConfig
from ._describe import (
    _describe_experiment_header,
    _describe_step_sequence,
    _format_declared_contracts,
)
from ._normalize import (
    _normalize_count,
    _normalize_lengths,
    _to_immutable_pairs,
)
from ._compile import lower_step
from ._compiled import CompiledClonalExperiment, CompiledExperiment
from ._refdata_resolver import (
    _CONFIG_ALIASES,
    ExperimentInput,
    _coerce_to_refdata_and_dataconfig,
    dataconfig_to_refdata,
)
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


# Sentinel for `using(...)` arguments that the caller did not pass at
# all. We can't use ``None`` for "unchanged" because ``None`` already
# means "clear the lock."
class _Unset:
    __slots__ = ()

    def __repr__(self) -> str:
        return "<UNSET>"


_UNSET: _Unset = _Unset()


# Inputs accepted by :meth:`Experiment.restrict_alleles` per segment. ``None``
# clears any prior lock; a single name is a one-allele lock; an
# iterable of names is a multi-allele subset; ``_UNSET`` (the default)
# leaves the existing lock unchanged.
_LockInput = Union[str, Iterable[str], None, _Unset]


# Refdata-resolver helpers (``_CONFIG_ALIASES``, ``_resolve_config_name``,
# ``dataconfig_to_refdata``, ``_coerce_to_refdata_and_dataconfig``,
# ``ExperimentInput``) live in :mod:`._refdata_resolver`. Re-exported
# below so the public surface (``GenAIRR.dataconfig_to_refdata``) and
# legacy internal callers (e.g. ``mcp_server`` reaching for
# ``_CONFIG_ALIASES``) keep working.


class Experiment:
    """Fluent builder for a simulation pipeline.

    Build an ``Experiment`` from a config name, a :class:`DataConfig`,
    or a :class:`GenAIRR._engine.RefDataConfig` via :meth:`Experiment.on`,
    chain configuration steps (currently just :meth:`recombine`), then
    call :meth:`run` (or :meth:`compile` for explicit two-stage flow).

    The builder is **stateful but not destructive**: each fluent call
    returns ``self`` after appending a step. The same ``Experiment``
    can be ``compile()``-d and ``run()``-d multiple times with
    different seeds.
    """

    __slots__ = (
        "_refdata",
        "_steps",
        "_dataconfig",
        "_locks",
        "_metadata",
        "_contracts",
    )

    def __init__(
        self,
        refdata: "_engine.RefDataConfig",
        dataconfig: Optional[DataConfig] = None,
    ) -> None:
        self._refdata = refdata
        self._dataconfig = dataconfig
        self._steps: List[Union[_RecombineStep, _MutateStep, _CorruptStep]] = []
        # Constraint bundles declared via `.productive_only()` etc.
        # Stored as a list so future bundles can compose. Today only
        # the productive bundle is recognized.
        self._contracts: List[str] = []
        # Per-segment allele-lock subsets set by ``.restrict_alleles(...)``. Each
        # entry is ``None`` (no lock — sample uniformly across the pool)
        # or a tuple of allele IDs to sample uniformly across.
        self._locks: Dict[str, Optional[Tuple[int, ...]]] = {
            "V": None,
            "D": None,
            "J": None,
        }
        # sample-level metadata to inject as columns on every
        # AIRR record (e.g. ``sample_id``, ``donor``,
        # ``repertoire_id``, ``cell_id``). Empty by default.
        self._metadata: Dict[str, Any] = {}

    @classmethod
    def on(cls, source: ExperimentInput) -> "Experiment":
        """Start an experiment against the given reference data.

        ``source`` is one of:
        - a config-name string (e.g. ``"human_igh"``),
        - a :class:`GenAIRR.DataConfig` instance,
        - a :class:`GenAIRR._engine.RefDataConfig`.

        When ``source`` is a config name or a ``DataConfig``, the
        underlying empirical distributions (NP lengths, per-gene
        trims) are kept on the experiment so :meth:`recombine` can
        use them as the default sampling distributions. A bare
        ``RefDataConfig`` has no such backing — :meth:`recombine`
        falls through to its uniform ``[0..6]`` placeholder unless
        the caller passes ``np1_lengths`` / ``np2_lengths``
        explicitly.
        """
        refdata, dataconfig = _coerce_to_refdata_and_dataconfig(source)
        return cls(refdata, dataconfig)

    @property
    def chain_type(self) -> str:
        """Chain type of the attached refdata (``"vj"`` or ``"vdj"``)."""
        return self._refdata.chain_type

    @property
    def step_count(self) -> int:
        """Number of fluent steps recorded on this experiment so far."""
        return len(self._steps)

    @property
    def refdata(self) -> "_engine.RefDataConfig":
        """The engine-native ``RefDataConfig`` this experiment is bound to."""
        return self._refdata

    def pcr_amplify(
        self,
        *,
        count: Optional[Union[int, Tuple[int, int], Iterable[Tuple[int, float]]]] = None,
        rate: Optional[float] = None,
    ) -> "Experiment":
        """Append a PCR-error step modelling substitution errors
        introduced during PCR amplification.

        **Specify intensity with exactly one of ``rate`` or ``count``.**

        ``rate`` is the per-base PCR error probability (e.g.
        ``1e-4`` per base per cycle, depending on polymerase
        fidelity). At execute time the engine draws
        ``count ~ Poisson(rate × pool_len)`` against each record's
        current sequence length — matching how PCR error is
        universally reported in the literature. This is the
        canonical biology-default form.

        ``count`` is the legacy explicit count distribution
        (fixed int, ``(low, high)`` range, or empirical
        ``(count, weight)`` list). Useful for benchmark scripts.

        Each error samples a uniform position in the assembled pool
        and replaces the base with a uniform A/C/G/T draw. Trace
        addresses: ``corrupt.pcr.{count, error_site[i], error_base[i]}``.

        Passing both ``count`` and ``rate`` raises ``ValueError``.
        Passing neither raises ``ValueError``.
        """
        return self._append_count_or_rate_step(
            kind=_CORRUPT_KIND_PCR,
            label="pcr_amplify",
            count=count,
            rate=rate,
        )

    def sequencing_errors(
        self,
        *,
        count: Optional[Union[int, Tuple[int, int], Iterable[Tuple[int, float]]]] = None,
        rate: Optional[float] = None,
    ) -> "Experiment":
        """Append a sequencing-quality-error step modelling base-call
        errors during sequencing readout.

        **Specify intensity with exactly one of ``rate`` or ``count``.**

        ``rate`` is the per-base sequencing error probability (e.g.
        ``1e-3`` for a Q30 base, ``1e-2`` for Q20). Drawn as
        ``count ~ Poisson(rate × pool_len)`` per record — the
        canonical Phred-quality framing in immunoseq.

        ``count`` is the legacy explicit count distribution.

        Same shape as :meth:`pcr_amplify` but each substitution
        writes the destination base in **lowercase** to mark the
        position as corrupted (the sequencing-error convention).
        """
        return self._append_count_or_rate_step(
            kind=_CORRUPT_KIND_QUALITY,
            label="sequencing_errors",
            count=count,
            rate=rate,
        )

    def _append_count_or_rate_step(
        self,
        *,
        kind: str,
        label: str,
        count: Optional[Union[int, Tuple[int, int], Iterable[Tuple[int, float]]]],
        rate: Optional[float],
    ) -> "Experiment":
        """Shared validator + step appender for the corruption methods
        that accept both ``count`` and ``rate`` (PCR errors, sequencing
        errors). See :meth:`mutate` for the same shape on the
        mutation side."""
        if count is not None and rate is not None:
            raise ValueError(
                f"{label}(): pass exactly one of `rate` or `count`, not both. "
                f"`rate` is the canonical biology default (e.g. rate=1e-4 "
                f"for per-base PCR error); `count` is the explicit per-record "
                f"count for benchmark / deterministic-count workflows."
            )
        if count is None and rate is None:
            raise ValueError(
                f"{label}(): pass exactly one of `rate` or `count`."
            )
        if rate is not None:
            if not isinstance(rate, (int, float)) or isinstance(rate, bool):
                raise TypeError(
                    f"rate must be a finite float in [0.0, 1.0], got "
                    f"{type(rate).__name__}"
                )
            if not (0.0 <= float(rate) <= 1.0):
                raise ValueError(
                    f"rate must be in [0.0, 1.0] (got {rate!r}); rate is a "
                    f"per-base error probability, not an absolute count."
                )
            self._steps.append(
                _CorruptStep(
                    kind=kind,
                    rate=float(rate),
                )
            )
            return self
        pairs = _normalize_count(count)
        self._steps.append(
            _CorruptStep(
                kind=kind,
                count_pairs=pairs,
            )
        )
        return self

    def polymerase_indels(
        self,
        *,
        count: Union[int, Tuple[int, int], Iterable[Tuple[int, float]]],
        insertion_prob: float = 0.5,
    ) -> "Experiment":
        """Append a polymerase-slippage indel step (PCR-stage artifact).

        Models the small insertions / deletions that arise from
        polymerase slippage during PCR amplification — single-base
        indels at low rate, typically <1 % per read at standard
        polymerase fidelity. Larger structural indels are rare and
        not modelled here.

        ``count`` is the per-simulation distribution over the total
        number of indel events. Each event independently chooses
        insertion (probability ``insertion_prob``) vs. deletion.
        Insertions sample a uniform A/C/G/T base. Indels change
        pool length and shift downstream region ranges accordingly.

        Raises ``ValueError`` if ``insertion_prob`` is outside
        ``[0.0, 1.0]`` or non-finite.
        """
        if (
            insertion_prob != insertion_prob  # NaN check
            or not (0.0 <= insertion_prob <= 1.0)
        ):
            raise ValueError(
                f"insertion_prob must be a finite number in [0.0, 1.0], "
                f"got {insertion_prob}"
            )
        pairs = _normalize_count(count)
        self._steps.append(
            _CorruptStep(
                kind=_CORRUPT_KIND_INDEL,
                count_pairs=pairs,
                insertion_prob=float(insertion_prob),
                apply_prob=0.0,
            )
        )
        return self

    def primer_trim_5prime(
        self,
        *,
        length: Union[int, Tuple[int, int], Iterable[Tuple[int, float]]],
    ) -> "Experiment":
        """Append a 5'-end loss corruption step.

        Drops bases from the start of the assembled sequence to model
        primer-region trimming. ``length`` accepts the same shapes as
        ``count`` on other corrupt ops:

        - ``length=10`` — strip exactly 10 bases.
        - ``length=(0, 20)`` — strip a uniform integer in ``[0, 20]``.
        - ``length=[(5, 1.0), (10, 1.0), (15, 1.0)]`` — empirical
          (length, weight) distribution.

        The actual loss is clamped to the pool length (so a sample
        larger than the sequence drops the whole pool). The loss is
        permanent for downstream passes — subsequent corruption
        operates on the shorter pool.
        """
        pairs = _normalize_count(length)
        self._steps.append(
            _CorruptStep(
                kind=_CORRUPT_KIND_5PRIME_LOSS,
                count_pairs=pairs,
                insertion_prob=0.0,
                apply_prob=0.0,
            )
        )
        return self

    def primer_trim_3prime(
        self,
        *,
        length: Union[int, Tuple[int, int], Iterable[Tuple[int, float]]],
    ) -> "Experiment":
        """Append a 3'-end loss corruption step.

        Drops bases from the end of the assembled sequence to model
        read-end degradation. Same ``length`` shapes as
        :meth:`primer_trim_5prime`.
        """
        pairs = _normalize_count(length)
        self._steps.append(
            _CorruptStep(
                kind=_CORRUPT_KIND_3PRIME_LOSS,
                count_pairs=pairs,
                insertion_prob=0.0,
                apply_prob=0.0,
            )
        )
        return self

    def ambiguous_base_calls(
        self,
        *,
        count: Union[int, Tuple[int, int], Iterable[Tuple[int, float]]],
    ) -> "Experiment":
        """Append an N-substitution corruption step.

        With the per-simulation count drawn from ``count``, replace
        that many uniform-random pool positions with the ambiguous
        base ``N``. Models the low-quality positions that real
        sequencers emit when the base caller cannot commit.

        ``count`` accepts the same shapes as other corruption ops:
        - ``count=10`` — fixed.
        - ``count=(0, 20)`` — uniform integer range.
        - ``count=[(5, 1.0), (10, 1.0)]`` — empirical (count, weight).

        Sampled sites are with-replacement, so ``count`` is the
        upper bound on resulting `N` count (collisions reduce it).
        """
        pairs = _normalize_count(count)
        self._steps.append(
            _CorruptStep(
                kind=_CORRUPT_KIND_NS,
                count_pairs=pairs,
                insertion_prob=0.0,
                apply_prob=0.0,
            )
        )
        return self

    def random_strand_orientation(self, *, prob: float = 0.5) -> "Experiment":
        """Append a reverse-complement corruption step.

        With probability ``prob`` the AIRR record-builder reverse-
        complements the ``sequence``, ``np1``, ``np2``, and
        ``junction`` fields and flips the corresponding pool-position
        coords (``*_sequence_start/end``, ``junction_start/end``).
        Alignment / germline / CIGAR / identity fields stay in
        forward orientation per the AIRR Rearrangement spec.
        ``prob=0.0`` is a no-op (coin flip recorded but never fires);
        ``prob=1.0`` always flips. The biological model is the ~50%
        of antisense reads in real immune-seq libraries.

        Raises ``ValueError`` if ``prob`` is outside ``[0.0, 1.0]``
        or non-finite.
        """
        if prob != prob or not (0.0 <= prob <= 1.0):
            raise ValueError(
                f"prob must be a finite number in [0.0, 1.0], got {prob}"
            )
        self._steps.append(
            _CorruptStep(
                kind=_CORRUPT_KIND_REV_COMP,
                count_pairs=(),
                insertion_prob=0.0,
                apply_prob=float(prob),
            )
        )
        return self

    def contaminate(self, *, prob: float) -> "Experiment":
        """Append a contaminant-replacement corruption step.

        With probability ``prob`` the entire assembled pool is
        overwritten with uniform A/C/G/T bases. Models primer
        dimers, bacterial DNA, or any non-receptor sequence in the
        library. ``prob=0.0`` is a no-op (coin flip recorded but
        never fires); ``prob=1.0`` always contaminates.

        Raises ``ValueError`` if ``prob`` is outside ``[0.0, 1.0]``
        or non-finite.
        """
        if prob != prob or not (0.0 <= prob <= 1.0):
            raise ValueError(
                f"prob must be a finite number in [0.0, 1.0], got {prob}"
            )
        self._steps.append(
            _CorruptStep(
                kind=_CORRUPT_KIND_CONTAMINANT,
                count_pairs=(),
                insertion_prob=0.0,
                apply_prob=float(prob),
            )
        )
        return self

    def with_metadata(self, **fields: Any) -> "Experiment":
        """Attach sample-level metadata to every AIRR record.

        Standard AIRR Repertoire fields like ``sample_id``,
        ``donor``, ``repertoire_id``, and ``cell_id`` are commonly
        used by multi-sample analysis tools — Change-O, scirpy,
        immcantation pipelines all expect them populated. Custom
        keys are also accepted and pass through unchanged.

        Subsequent ``with_metadata()`` calls **merge** with the
        prior set (per-key replacement). Pass ``key=None`` to clear
        a single key. Values are stored as-is and serialised via
        the standard CSV/TSV writers; non-string values are
        converted with ``str()`` on output.

        Example::

            (Experiment.on("human_igh")
             .with_metadata(sample_id="P1", donor="D001",
                            repertoire_id="R-001-IGH")
             .recombine().mutate(count=8))
        """
        for key, value in fields.items():
            if not isinstance(key, str):
                raise TypeError(
                    f"with_metadata: keys must be strings, got "
                    f"{type(key).__name__}"
                )
            if value is None:
                self._metadata.pop(key, None)
            else:
                self._metadata[key] = value
        return self

    def productive_only(self) -> "Experiment":
        """Require every emitted record to be a productive sequence.

        Attaches the canonical productive-sequence contract bundle to
        the experiment: junction frame in-register, no stop codons in
        the junction, and V/J anchor amino acids preserved. The bundle
        is enforced during recombination, mutation, and corruption
        passes by narrowing each pass's action support before
        sampling. If a constrained support is empty, permissive
        execution uses the pass's explicit no-op / sentinel behavior;
        use ``strict=True`` at run time to raise instead.

        **Failure surfaces** (see
        ``docs/productive_failure_mode_audit.md`` for the full matrix):

        - *Compile-time precondition* — when a sampling distribution
          is statically impossible under the bundle (e.g. every NP1
          length violates frame), :meth:`compile` raises
          ``ValueError`` regardless of the ``strict`` flag.
        - *Runtime fresh strict* (``run(..., strict=True)``) — when
          dynamic state makes a sampler's admissible support empty,
          raises :class:`GenAIRR._engine.StrictSamplingError` with
          structured ``(pass_name, address, reason)`` args.
        - *Runtime fresh permissive* (``strict=False``, default) —
          the pass records its declared sentinel (indel ``site=-1``,
          NP length ``0``, NP base ``N``, trim ``0``) or skips the
          slot; the record continues.
        - *Trace replay* (``replay_from_trace_file``) — consumes
          recorded values verbatim, does not re-evaluate
          admissibility. A permissive-sentinel trace replays cleanly
          even with ``strict=True``.

        **Order-independent.** This method is a constraint declaration,
        not a pipeline step — it can be called anywhere in the chain
        and the result is identical. Convention is to place it last
        (right before ``run_records()``) so the constraint reads as a
        post-hoc requirement on the emitted records.

        TCR refdata accepts the call but raises ``ValueError`` at
        :meth:`compile` time because TCRs don't have somatic
        hypermutation and the productive bundle's anchor checks
        assume BCR semantics. Catch this early at the builder if it
        matters to you.

        Example::

            result = (
                Experiment.on("human_igh")
                .recombine()
                .mutate(count=(5, 15))
                .productive_only()
                .run_records(seed=42)
            )
        """
        if "productive" not in self._contracts:
            self._contracts.append("productive")
        return self

    def expand_clones(
        self,
        *,
        n_clones: int,
        per_clone: int,
    ) -> "Experiment":
        """Expand the pipeline into clonal lineages.

        Marks the per-clone / per-descendant boundary in the chain:
        steps appended *before* this call run **once per clone** —
        typically just :meth:`recombine`, which establishes the
        parent V/D/J + trim + NP + assembled IR for the clonal
        family. Steps appended *after* this call run **once per
        read inside the family** — typically :meth:`mutate` and the
        library-prep / sequencing-stage steps, which introduce
        per-read divergence within the clone.

        Concrete shape::

            exp = (Experiment.on("human_igh")
                   .recombine()
                   .expand_clones(n_clones=10, per_clone=20)
                   .mutate(rate=0.05)
                   .pcr_amplify(count=2))
            result = exp.run_records(seed=0)
            # 10 clones × 20 descendants = 200 records.
            # Each record carries a ``clone_id`` integer in [0, 10).

        ``n`` can be omitted from :meth:`run_records` for a clonal
        experiment — the runtime expands ``n_clones * per_clone``
        records automatically. Passing ``n`` is allowed only when
        ``n == n_clones * per_clone``.

        Constraints:
        - Both ``n_clones`` and ``per_clone`` must be positive ints.
        - At most one expansion per pipeline; calling this method
          twice raises ``ValueError``.

        Implementation note: the runtime forks the parent's IR
        (final ``Simulation`` after the pre-fork plan) into
        descendants by running the post-fork plan from that IR
        with distinct seeds. Within a clone, every descendant
        shares the same recombination provenance (V allele, trim,
        NP bases) and only diverges through the post-fork passes.
        """
        if not isinstance(n_clones, int) or isinstance(n_clones, bool) or n_clones < 1:
            raise ValueError(
                f"n_clones must be a positive int, got {n_clones!r}"
            )
        if not isinstance(per_clone, int) or isinstance(per_clone, bool) or per_clone < 1:
            raise ValueError(f"per_clone must be a positive int, got {per_clone!r}")
        if any(isinstance(s, _ClonalForkStep) for s in self._steps):
            raise ValueError(
                "expand_clones() can only be called once per pipeline"
            )
        self._steps.append(_ClonalForkStep(n_clones=n_clones, size=per_clone))
        return self

    def mutate(
        self,
        *,
        model: str = "s5f",
        count: Optional[Union[int, Tuple[int, int], Iterable[Tuple[int, float]]]] = None,
        rate: Optional[float] = None,
        s5f_model: str = "hh_s5f",
    ) -> "Experiment":
        """Append a somatic-hypermutation step.

        ``model`` selects the mutation kernel:
        - ``"s5f"`` (default) — context-dependent SHM via the bundled
          S5F kernel named in ``s5f_model``. Available kernels:
          ``"hh_s5f"``, ``"hh_s5f_60"``, ``"hh_s5f_opposite"``,
          ``"hkl_s5f"``.
        - ``"uniform"`` — position-independent SHM. Each mutated
          position gets a uniformly drawn A/C/G/T replacement.

        **Specify intensity with exactly one of ``rate`` or ``count``.**

        ``rate`` is the per-base mutation rate (e.g. ``0.03`` for 3 %
        SHM, which is roughly memory B-cell SHM). At execute time the
        engine draws ``count ~ Poisson(rate × pool_len)`` against each
        record's current sequence length — so the realized count
        scales with each record's actual length, matching how
        immunologists report SHM in the literature. This is the
        canonical, biology-default form.

        ``count`` is the legacy explicit count distribution, useful
        for benchmark scripts that want a deterministic count
        independent of record length:
        - ``count=15`` — fixed: every simulation gets exactly 15
          mutations.
        - ``count=(5, 25)`` — uniform integer in ``[5, 25]`` (both
          endpoints inclusive).
        - ``count=[(5, 1.0), (10, 2.0), ...]`` — explicit empirical
          ``(count, weight)`` distribution.

        Passing both ``count`` and ``rate`` raises ``ValueError``.
        Passing neither raises ``ValueError``.

        **TCR guard:** somatic hypermutation is a B-cell
        phenomenon — T-cells do not undergo SHM in the periphery.
        Calling ``.mutate()`` on a TCR-configured experiment raises
        ``ValueError`` to prevent silent biological misuse. Use
        ``pcr_amplify`` / ``sequencing_errors`` for sequencing-error
        realism on TCR data instead.
        """
        if self._is_tcr_refdata():
            raise ValueError(
                "mutate(): somatic hypermutation does not occur in TCR "
                "sequences (T-cells lack AID and the SHM machinery). The "
                "configured refdata is a TCR locus. For sequencing-error "
                "realism on TCR data, use pcr_amplify / sequencing_errors "
                "/ polymerase_indels instead."
            )
        model_lc = model.lower()
        if model_lc not in ("uniform", "s5f"):
            raise ValueError(
                f"model must be 'uniform' or 's5f' (got {model!r})"
            )
        if count is not None and rate is not None:
            raise ValueError(
                "mutate(): pass exactly one of `rate` or `count`, not both. "
                "`rate` is the canonical biology default (e.g. rate=0.03 "
                "for 3% SHM); `count` is the explicit per-record count "
                "for benchmark / deterministic-count workflows."
            )
        if count is None and rate is None:
            raise ValueError(
                "mutate(): pass exactly one of `rate` or `count`. "
                "Suggested default: rate=0.03 (~3% SHM, memory B-cell range)."
            )
        if rate is not None:
            if not isinstance(rate, (int, float)) or isinstance(rate, bool):
                raise TypeError(
                    f"rate must be a finite float in [0.0, 1.0], got "
                    f"{type(rate).__name__}"
                )
            if not (0.0 <= float(rate) <= 1.0):
                raise ValueError(
                    f"rate must be in [0.0, 1.0] (got {rate!r}); rate is a "
                    f"per-base mutation probability, not an absolute count."
                )
            self._steps.append(
                _MutateStep(
                    model=model_lc,
                    s5f_model_name=s5f_model,
                    rate=float(rate),
                )
            )
            return self
        pairs = _normalize_count(count)
        self._steps.append(
            _MutateStep(
                model=model_lc,
                s5f_model_name=s5f_model,
                count_pairs=pairs,
            )
        )
        return self

    def _is_tcr_refdata(self) -> bool:
        """Detect whether the bound refdata is a TCR locus.

        TCR allele names are prefixed with ``TR`` (TRA, TRB, TRG,
        TRD); BCR with ``IG``. Inspecting the first V allele is
        enough — locus is uniform across the pool. Returns ``False``
        when the V pool is empty (defensive — let downstream errors
        surface that condition).
        """
        if self._refdata.v_pool_size() == 0:
            return False
        first_v_name = self._refdata.v_allele(0).name
        return first_v_name.upper().startswith("TR")

    def restrict_alleles(
        self,
        *,
        v: "_LockInput" = _UNSET,
        d: "_LockInput" = _UNSET,
        j: "_LockInput" = _UNSET,
    ) -> "Experiment":
        """Narrow allele sampling to a named subset, per segment.

        Sampling stays uniform — this restricts the *support* of the
        sampling distribution to the listed allele names, not pins a
        single allele. If you supply one name, the result is
        effectively deterministic (uniform over one element); for any
        list of N names, the recombination pass samples uniformly
        over those N alleles.

        Each kwarg accepts:
        - a single allele name string (e.g. ``"IGHV1-2*02"``),
        - a list / tuple of allele names (sample uniformly among them),
        - ``None`` to clear a previously-set restriction for that segment,
        - omitted (default) — the existing restriction for that segment
          is left unchanged.

        The restrictions are applied to the next ``recombine()`` step's
        ``push_sample_allele`` calls at compile time. Calling
        ``restrict_alleles()`` more than once overlays the new values
        onto the previous restrictions (per-segment), so
        ``.restrict_alleles(v="A").restrict_alleles(d="B")`` restricts
        both V and D.

        Raises:
        - ``ValueError`` if an allele name doesn't exist in the
          configured refdata, or if a restriction is set for ``D`` on
          a VJ chain.
        - ``TypeError`` if an unexpected input shape is passed.
        """
        for segment, value in (("V", v), ("D", d), ("J", j)):
            if value is _UNSET:
                continue
            if value is None:
                self._locks[segment] = None
                continue
            ids = self._resolve_lock(segment, value)
            self._locks[segment] = ids
        return self

    def _resolve_lock(
        self,
        segment: str,
        value: Union[str, Iterable[str]],
    ) -> Tuple[int, ...]:
        """Resolve allele-name(s) → tuple of allele IDs against this
        experiment's refdata. Raises on unknown names, on D locks
        for VJ chains, or on non-string inputs."""
        if segment == "D" and self._refdata.chain_type != "vdj":
            raise ValueError(
                f"cannot lock D alleles on a {self._refdata.chain_type!r} chain"
            )

        if isinstance(value, str):
            names: Tuple[str, ...] = (value,)
        else:
            try:
                names = tuple(value)
            except TypeError as exc:
                raise TypeError(
                    f"restrict_alleles(): {segment} lock must be a name string or an "
                    f"iterable of name strings; got {type(value).__name__}"
                ) from exc
            if not all(isinstance(n, str) for n in names):
                raise TypeError(
                    f"restrict_alleles(): {segment} lock entries must all be strings"
                )
            if not names:
                raise ValueError(
                    f"restrict_alleles(): {segment} lock list must be non-empty "
                    "(pass None to clear instead)"
                )

        index = self._allele_name_index(segment)
        ids: List[int] = []
        seen: set = set()
        for name in names:
            if name not in index:
                raise ValueError(
                    f"restrict_alleles(): no {segment} allele named {name!r} in refdata "
                    "(check spelling against `list_alleles` or the loaded config)"
                )
            allele_id = index[name]
            if allele_id in seen:
                raise ValueError(
                    f"restrict_alleles(): duplicate {segment} allele {name!r} in lock list"
                )
            seen.add(allele_id)
            ids.append(allele_id)
        return tuple(ids)

    def _allele_name_index(self, segment: str) -> Dict[str, int]:
        """Build a name → allele_id map for the given segment by
        scanning the refdata pool. Cheap enough to do per-call given
        typical pool sizes (≤ a few hundred) and the rarity of
        ``.restrict_alleles()``."""
        if segment == "V":
            n = self._refdata.v_pool_size()
            getter = self._refdata.v_allele
        elif segment == "D":
            n = self._refdata.d_pool_size()
            getter = self._refdata.d_allele
        elif segment == "J":
            n = self._refdata.j_pool_size()
            getter = self._refdata.j_allele
        else:  # pragma: no cover — guarded above.
            raise ValueError(f"unsupported segment {segment!r}")
        return {getter(i).name: i for i in range(n)}

    def recombine(
        self,
        *,
        np1_lengths: Optional[Iterable[Tuple[int, float]]] = None,
        np2_lengths: Optional[Iterable[Tuple[int, float]]] = None,
        v_allele_weights: Optional[Dict[str, float]] = None,
        d_allele_weights: Optional[Dict[str, float]] = None,
        j_allele_weights: Optional[Dict[str, float]] = None,
    ) -> "Experiment":
        """Append a standard V(D)J recombination step.

        Compiles to:
        - **VJ:** sample V → sample J → (trim V_3, J_5) → assemble V
          → generate NP1 → assemble J.
        - **VDJ:** sample V → sample D → sample J → (trim V_3, D_5,
          D_3, J_5) → assemble V → generate NP1 → assemble D →
          generate NP2 → assemble J.

        ``np1_lengths`` / ``np2_lengths`` default to the species'
        empirical NP-length distributions (from
        ``DataConfig.NP_lengths``) when the experiment is bound to
        a DataConfig. For raw-RefDataConfig experiments where no
        empirical data is available, both fall back to the uniform
        ``[(0, 1.0), ..., (6, 1.0)]`` distribution and emit a
        :class:`UserWarning` so the caller knows the synthetic
        default is being used. Pass an explicit iterable of
        ``(length, weight)`` tuples to override the default.
        Passing ``np2_lengths`` on a VJ chain raises ``ValueError``
        (VJ chains have no NP2 region — there's no D segment to
        bracket).

        **Exonuclease trim** is enabled by default and uses the
        empirical per-segment trim distributions from the bound
        ``DataConfig`` when available. To disable trim, or supply
        custom trim distributions, call :meth:`trim` *after*
        :meth:`recombine` in the chain (before any mutation /
        corruption step). On a raw ``RefDataConfig`` without trim
        data, recombine emits a :class:`UserWarning` and falls back
        to a no-op trim.

        ``v_allele_weights`` / ``d_allele_weights`` /
        ``j_allele_weights`` — optional ``{allele_name:
        weight}`` dicts that bias allele sampling. Listed alleles
        get the supplied positive weight; unlisted alleles default
        to 1.0, so e.g. ``v_allele_weights={"IGHV3-23*01": 100}``
        boosts that allele while keeping every other V allele
        possible at 1/100th the rate. Mutually exclusive with the
        per-segment :meth:`restrict_alleles` restriction. Raises
        ``ValueError`` for unknown allele names or non-positive
        weights.
        """
        # VJ chains have no NP2 region — surface user mistakes loudly
        # instead of silently dropping the argument.
        if np2_lengths is not None and self._refdata.chain_type != "vdj":
            raise ValueError(
                f"np2_lengths is only valid for VDJ chains; the bound "
                f"refdata is {self._refdata.chain_type!r} (no D segment, "
                f"no NP2 region). Drop the np2_lengths kwarg or bind a "
                f"VDJ refdata."
            )

        # Trim is default-on; .trim() may later disable or replace it.
        trim = True
        defaults = self._recombine_defaults() if trim or self._dataconfig else None

        # Raw-RefDataConfig path: there's no DataConfig backing this
        # experiment, so empirical NP lengths and trim distributions
        # don't exist. We fall back to a uniform NP-length default —
        # historically silent — and surface a warning so the synthetic
        # default isn't mistaken for real biology. The trim warning is
        # emitted at compile() time instead, after .trim() has had a
        # chance to disable trim explicitly.
        if self._dataconfig is None:
            if np1_lengths is None or (
                self._refdata.chain_type == "vdj" and np2_lengths is None
            ):
                warnings.warn(
                    "Experiment bound to a raw RefDataConfig with no empirical "
                    "NP-length distribution; falling back to uniform "
                    "[(0, 1.0), ..., (6, 1.0)]. Pass np1_lengths "
                    "(and np2_lengths for VDJ) explicitly to silence this.",
                    UserWarning,
                    stacklevel=2,
                )

        np1 = (
            _normalize_lengths(np1_lengths)
            if np1_lengths is not None
            else self._default_np_lengths(defaults, "np1")
        )
        np2 = (
            _normalize_lengths(np2_lengths)
            if np2_lengths is not None
            else self._default_np_lengths(defaults, "np2")
        )

        if trim and defaults is not None:
            trim_v_3 = _to_immutable_pairs(defaults.get("trim_v_3"))
            trim_d_5 = _to_immutable_pairs(defaults.get("trim_d_5"))
            trim_d_3 = _to_immutable_pairs(defaults.get("trim_d_3"))
            trim_j_5 = _to_immutable_pairs(defaults.get("trim_j_5"))
        else:
            trim_v_3 = trim_d_5 = trim_d_3 = trim_j_5 = None

        weights_v = self._resolve_allele_weights("V", v_allele_weights)
        weights_d = self._resolve_allele_weights("D", d_allele_weights)
        weights_j = self._resolve_allele_weights("J", j_allele_weights)

        step = _RecombineStep(
            np1_lengths=np1,
            np2_lengths=np2,
            trim_v_3=trim_v_3,
            trim_d_5=trim_d_5,
            trim_d_3=trim_d_3,
            trim_j_5=trim_j_5,
            weights_v=weights_v,
            weights_d=weights_d,
            weights_j=weights_j,
        )
        self._steps.append(step)
        return self

    def trim(
        self,
        *,
        enabled: bool = True,
        v_3: Optional[Iterable[Tuple[int, float]]] = None,
        d_5: Optional[Iterable[Tuple[int, float]]] = None,
        d_3: Optional[Iterable[Tuple[int, float]]] = None,
        j_5: Optional[Iterable[Tuple[int, float]]] = None,
    ) -> "Experiment":
        """Configure exonuclease trim on the preceding :meth:`recombine`
        step.

        Trim is a per-segment exonuclease step that lives biologically
        *inside* V(D)J recombination (between allele selection and NP
        insertion). It is on by default with empirical distributions
        sourced from the bound ``DataConfig``. Use this method only
        when you need to override that default:

        - ``trim(enabled=False)`` — disable trim entirely. Equivalent to
          recombining against the raw allele endpoints.
        - ``trim(v_3=..., d_5=..., d_3=..., j_5=...)`` — supply custom
          per-segment trim-length distributions as ``(length, weight)``
          iterables. Omitted segments keep their empirical defaults
          (or no-op on raw RefDataConfig).

        **Position in the chain.** ``trim()`` is configuration applied
        to the most recent ``recombine()`` step. It must appear after
        ``recombine()`` and before any mutation / corruption step.
        Calling it before ``recombine()`` raises ``ValueError``;
        calling it after a mutation step raises ``ValueError``.

        Raises ``ValueError`` if no ``recombine()`` step has been
        appended yet, or if the chain already has steps that would
        biologically follow recombination (mutate, corrupt_*).
        """
        from dataclasses import replace as _replace

        # Find the most recent recombine step.
        rec_idx: Optional[int] = None
        for i in range(len(self._steps) - 1, -1, -1):
            if isinstance(self._steps[i], _RecombineStep):
                rec_idx = i
                break
        if rec_idx is None:
            raise ValueError(
                "trim() must be called after recombine(); no recombine "
                "step is on this Experiment yet."
            )

        # Anything appended after the recombine step is wrong-ordered
        # for trim configuration.
        for j in range(rec_idx + 1, len(self._steps)):
            offending = type(self._steps[j]).__name__
            raise ValueError(
                f"trim() must be called immediately after recombine() — "
                f"before any mutation / corruption / clonal-fork step. "
                f"Found a {offending!r} between the latest recombine() "
                f"and this trim() call."
            )

        prior: _RecombineStep = self._steps[rec_idx]  # type: ignore[assignment]
        if not enabled:
            # Disable all trim slots, keep NP + weights untouched.
            new_step = _replace(
                prior,
                trim_v_3=None,
                trim_d_5=None,
                trim_d_3=None,
                trim_j_5=None,
                trim_overridden=True,
            )
            self._steps[rec_idx] = new_step
            return self

        # Override individual distributions; pass-through on None.
        def _resolve(
            current: Optional[Tuple[Tuple[int, float], ...]],
            override: Optional[Iterable[Tuple[int, float]]],
        ) -> Optional[Tuple[Tuple[int, float], ...]]:
            if override is None:
                return current
            return _normalize_lengths(override)

        new_step = _replace(
            prior,
            trim_v_3=_resolve(prior.trim_v_3, v_3),
            trim_d_5=_resolve(prior.trim_d_5, d_5),
            trim_d_3=_resolve(prior.trim_d_3, d_3),
            trim_j_5=_resolve(prior.trim_j_5, j_5),
            trim_overridden=True,
        )
        self._steps[rec_idx] = new_step
        return self

    def _resolve_allele_weights(
        self,
        segment: str,
        user_weights: Optional[Dict[str, float]],
    ) -> Optional[Tuple[float, ...]]:
        """Build a dense pool-aligned weight vector from the
        user-supplied ``{name: weight}`` dict. Listed alleles get the
        supplied weight; everything else gets ``1.0``. Returns
        ``None`` when no weights were supplied (preserving the
        upstream uniform-default behavior).

        Raises ``ValueError`` for unknown allele names, non-positive
        weights, or D weights on a VJ chain.
        """
        if user_weights is None:
            return None
        if not user_weights:
            raise ValueError(
                f"recombine: {segment.lower()}_allele_weights must contain "
                "at least one (name, weight) entry"
            )
        if segment == "D" and self._refdata.chain_type != "vdj":
            raise ValueError(
                f"recombine: cannot weight D alleles on a "
                f"{self._refdata.chain_type!r} chain"
            )

        index = self._allele_name_index(segment)
        if not index:
            return None  # No alleles for this segment (e.g. VJ + D).

        # Dense vector indexed by allele_id; default 1.0.
        pool_size = max(index.values()) + 1
        weights: List[float] = [1.0] * pool_size
        for name, w in user_weights.items():
            if not isinstance(name, str):
                raise TypeError(
                    f"recombine: {segment.lower()}_allele_weights keys must "
                    f"be allele-name strings, got {type(name).__name__}"
                )
            if isinstance(w, bool) or not isinstance(w, (int, float)):
                raise TypeError(
                    f"recombine: {segment.lower()}_allele_weights['{name}'] "
                    f"must be numeric, got {type(w).__name__}"
                )
            wf = float(w)
            if not (wf > 0.0 and wf == wf and wf != float("inf")):
                raise ValueError(
                    f"recombine: {segment.lower()}_allele_weights['{name}'] "
                    f"must be a finite positive number, got {wf}"
                )
            if name not in index:
                raise ValueError(
                    f"recombine: no {segment} allele named {name!r} in "
                    "refdata (check spelling against `list_alleles` or the "
                    "loaded config)"
                )
            weights[index[name]] = wf
        return tuple(weights)

    def _recombine_defaults(self):
        """Lazy-extract the empirical distributions for this
        experiment's DataConfig. Returns ``None`` for raw-RefDataConfig
        experiments. Cached on first call.
        """
        if self._dataconfig is None:
            return None
        from ._dataconfig_extract import extract_recombine_defaults

        return extract_recombine_defaults(self._dataconfig)

    def _build_contracts(self) -> Optional["_engine.ContractSet"]:
        """Synthesize the engine ``ContractSet`` from declared bundles.

        Today only the productive bundle is recognized. Future bundles
        (e.g. ``.in_frame_only()``) would compose here. Returns
        ``None`` when no constraint methods have been called.
        """
        if not self._contracts:
            return None
        # Composition across bundles isn't supported by the engine yet,
        # so for now exactly one bundle is allowed.
        if len(self._contracts) > 1:
            raise NotImplementedError(
                f"composing multiple constraint bundles is not yet supported; "
                f"got {self._contracts!r}"
            )
        bundle = self._contracts[0]
        if bundle == "productive":
            return _engine.productive()
        raise NotImplementedError(
            f"unknown constraint bundle {bundle!r}"
        )

    @staticmethod
    def _default_np_lengths(
        defaults, key: str
    ) -> Tuple[Tuple[int, float], ...]:
        """Pick the NP-length distribution to use when the user
        didn't pass an explicit one: empirical if available, else
        the uniform placeholder.
        """
        if defaults is not None:
            empirical = defaults.get(key)
            if empirical:
                return tuple(empirical)
        return tuple(_DEFAULT_NP_LENGTHS)

    def compile(self):
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
        """
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
                pre_steps, contracts, any_lock, replace_fn=_replace
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
                post_steps, contracts, any_lock=False, replace_fn=_replace
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

        simulator = self._build_simulator(
            self._steps, contracts, any_lock, replace_fn=_replace
        )
        return CompiledExperiment(
            simulator,
            self._refdata,
            steps=tuple(self._steps),
            dataconfig=self._dataconfig,
            metadata=self._metadata,
        )

    def _build_simulator(
        self,
        steps,
        contracts,
        any_lock: bool,
        *,
        replace_fn,
    ):
        """Compile a list of steps into a `GenAIRR._engine.CompiledSimulator`.
        Lifted out of `compile()` so the clonal-fork branch can build
        two simulators from sub-step-lists with a shared body."""
        plan = _engine.PassPlan()
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
            lower_step(step, plan, self._refdata)
        return plan.compile(refdata=self._refdata, respect=contracts)

    def run_records(
        self,
        *,
        n: Optional[int] = None,
        seed: int = 0,
        strict: bool = False,
        expose_provenance: bool = False,
    ) -> "SimulationResult":
        """Compile and run, then return the batch as a
        :class:`SimulationResult` ready for ``.to_csv`` / ``.to_fasta``
        / ``.to_dataframe`` export.

        For non-clonal experiments ``n`` defaults to 1. For clonal
        experiments (when the pipeline contains :meth:`with_clonal
        _structure`) ``n`` defaults to ``n_clones * size`` and may
        be omitted; passing ``n`` explicitly is allowed only if it
        matches that product.

        ``expose_provenance=True`` appends ``truth_v_call``,
        ``truth_d_call``, ``truth_j_call`` columns containing the
        originally-sampled allele names — distinct from the
        evidence-driven ``v_call`` / ``d_call`` / ``j_call`` fields
        an aligner would produce. Useful for benchmarking aligners
        against ground truth without keeping a side truth file.

        ``strict`` semantics match :meth:`run` — strict-mode applies
        only to **fresh sampling**. Trace replay
        (:meth:`CompiledExperiment.replay_from_trace_file`) consumes
        recorded sentinel values verbatim, so a permissive trace
        replays cleanly even with ``strict=True``. See
        ``docs/productive_failure_mode_audit.md`` §5.

        Returns a :class:`SimulationResult`; clonal records carry
        an integer ``clone_id`` field per row.
        """
        compiled = self.compile()
        if isinstance(compiled, CompiledClonalExperiment):
            result = compiled.run_records(
                n=n,
                seed=seed,
                strict=strict,
                expose_provenance=expose_provenance,
            )
        else:
            result = compiled.run_records(
                n=1 if n is None else n,
                seed=seed,
                strict=strict,
                expose_provenance=expose_provenance,
            )
        if self._metadata:
            for rec in result.records:
                for key, value in self._metadata.items():
                    rec[key] = value
        return result

    def run(
        self,
        *,
        n: Optional[int] = None,
        seed: int = 0,
        strict: bool = False,
    ) -> List["_engine.Outcome"]:
        """Compile and run this experiment ``n`` times.

        Equivalent to
        ``self.compile().run(n=n, seed=seed, strict=strict)``.
        Returns a list of :class:`GenAIRR._engine.Outcome` objects in
        clone-major order for clonal experiments.

        Attach :meth:`productive_only` (or any future constraint
        method) to the chain to require admissible records; the
        runtime filters NP base draws, length samples, and mutation
        / contamination substitutions in real time so the resulting
        sequences satisfy the bundle by construction.

        Statically impossible contract configurations fail during
        ``compile()`` with ``ValueError``. For runtime residue
        — i.e., empty admissible support emerging dynamically at
        sample time — ``strict=False`` (default) lets a pass consume
        the slot as its explicit no-op / sentinel; ``strict=True``
        raises :class:`GenAIRR._engine.StrictSamplingError` instead.

        Note the two error paths use **different exception classes**:
        ``ValueError`` for compile-time preconditions,
        ``StrictSamplingError`` (subclass of ``Exception``, NOT of
        ``ValueError``) for runtime empty-support. A bare
        ``except ValueError:`` will not catch the runtime case. See
        ``docs/productive_failure_mode_audit.md`` §6.1.

        ``strict`` only governs **fresh sampling**. Trace replay
        (:meth:`CompiledExperiment.replay_from_trace_file`) consumes
        recorded values verbatim; a permissive-recorded sentinel
        trace replays cleanly even with ``strict=True``. To re-execute
        a trace under strict-fresh semantics, call
        ``simulator.run(seed=<original_seed>, strict=True)`` instead.
        """
        compiled = self.compile()
        if isinstance(compiled, CompiledClonalExperiment):
            return compiled.run(n=n, seed=seed, strict=strict)
        return compiled.run(n=1 if n is None else n, seed=seed, strict=strict)

    def stream(
        self,
        *,
        n: Optional[int] = None,
        seed: int = 0,
        strict: bool = False,
    ) -> Iterator["_engine.Outcome"]:
        """Compile and lazily yield :class:`GenAIRR._engine.Outcome`
        objects. See :meth:`CompiledExperiment.stream` for full
        semantics."""
        return self.compile().stream(n=n, seed=seed, strict=strict)

    def stream_records(
        self,
        *,
        n: Optional[int] = None,
        seed: int = 0,
        strict: bool = False,
        id_prefix: str = "seq",
    ) -> Iterator[Dict[str, Any]]:
        """Compile and lazily yield AIRR-format record dicts. See
        :meth:`CompiledExperiment.stream_records`."""
        return self.compile().stream_records(
            n=n,
            seed=seed,
            strict=strict,
            id_prefix=id_prefix,
        )

    def describe(self) -> str:
        """Render a biology-style narrative of this experiment.

        The output is one line per step, prefixed with its position
        in the chain. Locks, allele weights, NP-length distributions,
        SHM kernels, and per-corruption rates are all surfaced. Use
        this to sanity-check what a fluent chain actually encodes —
        if it doesn't read like an immunology protocol, the chain is
        too murky.

        Returns a multi-line string ending without a trailing newline.
        Safe to ``print(exp.describe())``.

        Example::

            >>> print(Experiment.on("human_igh").recombine().mutate(count=(5, 15)).describe())
            Experiment on human_igh (vdj, DataConfig)
              1. V(D)J recombination: sample V/D/J alleles; empirical exonuclease trim (V3', D5', D3', J5'); insert NP1 (0–11 weighted bases) and NP2 (0–11 weighted bases)
              2. Somatic hypermutation (S5F context model, human heavy-chain (HH_S5F)): 5–15 mutations/record
        """
        header = _describe_experiment_header(self._refdata, self._dataconfig)
        if not self._steps and not self._contracts:
            return header + "\n  (no steps appended yet)"

        lines = [header]
        resolved = self._steps_with_locks_resolved()
        body_lines = _describe_step_sequence(resolved, self._refdata.chain_type)
        lines.extend(body_lines)
        contracts_line = _format_declared_contracts(self._contracts)
        if contracts_line:
            lines.append(f"  Constraints: {contracts_line}")
        if self._metadata:
            stamps = ", ".join(f"{k}={v!r}" for k, v in self._metadata.items())
            lines.append(f"  Metadata stamped on every record: {stamps}")
        return "\n".join(lines)

    def _steps_with_locks_resolved(self) -> List[Any]:
        """Return a copy of ``self._steps`` with per-segment allele
        locks (from :meth:`restrict_alleles`) injected into the first
        :class:`_RecombineStep`. Compile-time injection lives in
        :meth:`_build_simulator`; ``describe()`` needs the same
        substitution to render lock info correctly without
        side-effecting the live step list."""
        from dataclasses import replace as _replace

        any_lock = any(self._locks[seg] is not None for seg in ("V", "D", "J"))
        if not any_lock:
            return list(self._steps)
        out: List[Any] = []
        injected = False
        for step in self._steps:
            if not injected and isinstance(step, _RecombineStep):
                step = _replace(
                    step,
                    locks_v=self._locks["V"],
                    locks_d=self._locks["D"],
                    locks_j=self._locks["J"],
                )
                injected = True
            out.append(step)
        return out

    def __repr__(self) -> str:
        return f"<Experiment chain={self.chain_type} steps={self.step_count}>"

