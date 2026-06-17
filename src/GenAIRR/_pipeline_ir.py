"""Typed pipeline IR for :class:`GenAIRR.experiment.Experiment`.

Every fluent step on :class:`Experiment` (``recombine``, ``mutate``,
``pcr_amplify``, ``polymerase_indels``, ``expand_clones``, …) appends
one of the dataclasses defined here to ``Experiment._steps``.

This module owns the **step shape** (pure frozen dataclasses,
no engine awareness). Lowering each step onto the engine-native
:class:`_engine.PassPlan` lives in :mod:`._compile`. The split lets
the pipeline IR be inspected, round-tripped, or serialized without
importing the compiled Rust extension.

The module is intentionally private (single leading underscore on
the filename). The dataclasses are themselves prefixed with ``_`` to
discourage external dependency — the typed-IR shape may change in
minor releases.
"""
from __future__ import annotations

from dataclasses import dataclass
from typing import Optional, Sequence, Tuple


# Default NP-length distribution: lengths 0..6 with equal weight.
# Used by raw-RefDataConfig paths with no empirical NP_lengths and by
# the ``describe()`` helpers as the marker for "synthetic fallback."
_DEFAULT_NP_LENGTHS: Sequence[Tuple[int, float]] = tuple((i, 1.0) for i in range(7))


# Names of the corruption passes the Rust kernel exposes. Used both
# as a discriminator on :class:`_CorruptStep` and (indirectly) as
# the Python-facing kwarg name on the :class:`Experiment` builder.
_CORRUPT_KIND_PCR = "pcr"
_CORRUPT_KIND_QUALITY = "quality"
_CORRUPT_KIND_INDEL = "indel"
_CORRUPT_KIND_CONTAMINANT = "contaminant"
_CORRUPT_KIND_REV_COMP = "rev_comp"
_CORRUPT_KIND_5PRIME_LOSS = "5prime_loss"
_CORRUPT_KIND_3PRIME_LOSS = "3prime_loss"
_CORRUPT_KIND_NS = "ns"


@dataclass(frozen=True)
class _RecombineStep:
    """One ``recombine()`` invocation, captured for later compilation.

    ``np1_lengths`` / ``np2_lengths`` are the empirical NP-length
    distributions to use. ``trim_*`` are optional empirical trim
    distributions; when ``None`` the corresponding trim pass is
    omitted (so a synthetic refdata without trim data still works).

    ``locks_*`` are optional allele-ID subsets supplied by
    :meth:`Experiment.restrict_alleles`. When set, the matching
    ``push_sample_allele`` call samples uniformly only from those
    IDs instead of the full pool.
    """

    np1_lengths: Tuple[Tuple[int, float], ...]
    np2_lengths: Tuple[Tuple[int, float], ...]
    trim_v_3: Optional[Tuple[Tuple[int, float], ...]]
    trim_d_5: Optional[Tuple[Tuple[int, float], ...]]
    trim_d_3: Optional[Tuple[Tuple[int, float], ...]]
    trim_j_5: Optional[Tuple[Tuple[int, float], ...]]
    locks_v: Optional[Tuple[int, ...]] = None
    locks_d: Optional[Tuple[int, ...]] = None
    locks_j: Optional[Tuple[int, ...]] = None
    # V(D)J N-addition base sampling pairs, resolved from the
    # cartridge's typed ``ReferenceEmpiricalModels.np_bases`` plane
    # at extract time (Slice — Typed NP base model). Each entry is
    # a tuple of ``(base_byte, weight)`` pairs in the canonical
    # A/C/G/T alphabet that ``push_generate_np(base_pairs=…)``
    # passes through to the Rust ``CategoricalBase`` distribution.
    # ``None`` means no typed base model is configured — the
    # pre-slice ``UniformBase`` (4-way equal) applies, byte-
    # identical to legacy behaviour.
    np1_base_pairs: Optional[Tuple[Tuple[int, float], ...]] = None
    np2_base_pairs: Optional[Tuple[Tuple[int, float], ...]] = None
    # Markov transition matrices (Slice — Markov NP Base
    # Generator). Each is a 4-row tuple in canonical A/C/G/T
    # from-base order; each row is a 4-element tuple of
    # ``(to_base_byte, weight)`` pairs. Set together with
    # ``np{1,2}_base_pairs`` when the cartridge's typed
    # ``NpBaseModelSpec`` has ``kind="markov"`` — the
    # ``base_pairs`` carry the position-0 first-base row, the
    # transitions cover positions 1+. ``None`` for uniform /
    # empirical_first_base.
    np1_markov_transitions: Optional[
        Tuple[Tuple[Tuple[int, float], ...], ...]
    ] = None
    np2_markov_transitions: Optional[
        Tuple[Tuple[Tuple[int, float], ...], ...]
    ] = None
    # Per-end P-nucleotide length distributions (Slice —
    # P-nucleotide v1). Each is a tuple of ``(length, weight)``
    # pairs the engine consumes verbatim. ``None`` means the
    # cartridge didn't author a P-pass for that end, in which
    # case ``_lower_recombine`` omits the corresponding
    # `push_p_addition` call — byte-identical to the pre-slice
    # baseline. VJ chains never carry D-end values (validated
    # at the spec layer).
    p_v_3_lengths: Optional[Tuple[Tuple[int, float], ...]] = None
    p_d_5_lengths: Optional[Tuple[Tuple[int, float], ...]] = None
    p_d_3_lengths: Optional[Tuple[Tuple[int, float], ...]] = None
    p_j_5_lengths: Optional[Tuple[Tuple[int, float], ...]] = None
    weights_v: Optional[Tuple[float, ...]] = None
    weights_d: Optional[Tuple[float, ...]] = None
    weights_j: Optional[Tuple[float, ...]] = None
    # Set to True by ``Experiment.trim(...)``. Drives the
    # raw-RefDataConfig trim-noop warning at compile() time: when the
    # user explicitly called .trim(), silence the warning regardless
    # of whether trim is enabled or disabled.
    trim_overridden: bool = False


@dataclass(frozen=True)
class _MutateStep:
    """One ``mutate()`` invocation, captured for later compilation.

    ``model`` is either ``"uniform"`` (position-independent SHM via
    ``UniformMutationPass``) or ``"s5f"`` (context-dependent SHM via
    ``S5FMutationPass`` with the bundled S5F kernel named in
    ``s5f_model_name``).

    Exactly one of ``count_pairs`` and ``rate`` is set; the other is
    ``None``. ``count_pairs`` is the empirical distribution over the
    number of mutations per simulation, expressed as
    ``((count, weight), ...)``. ``rate`` is a per-base mutation rate
    (e.g. ``0.03`` for 3 % SHM); at execute time the engine draws
    ``count ~ Poisson(rate × pool_len)`` against the current pool
    length so the realized mutation count scales with each record's
    own sequence length — matching how immunologists report SHM in
    the literature.

    ``segment_rates`` is the normalised 4-bucket biological-segment
    rate vector ``(v, d, j, np)`` applied to site selection. Each
    entry is a finite non-negative ``float``; ``(1.0, 1.0, 1.0, 1.0)``
    is the default (flat-substrate) behaviour and produces
    byte-identical output to the pre-slice engine. The DSL
    boundary builds this tuple from the user's optional
    ``{"V"|"D"|"J"|"NP": float}`` dict; values outside the four
    buckets, negatives, NaN/inf, or all-zero are rejected at
    builder time.

    ``v_subregion_rates`` is the normalised 5-label V-subregion
    rate vector ``(fwr1, cdr1, fwr2, cdr2, fwr3)`` applied as an
    additional multiplicative factor to V-site selection (Slice B
    — ``docs/v_subregion_shm_rate_design.md``). Composes
    multiplicatively with ``segment_rates`` on V sites; non-V
    sites are unaffected. V sites on alleles without IMGT
    subregion annotations receive factor ``1.0`` so mixed
    cartridges stay usable. Default ``(1.0, 1.0, 1.0, 1.0, 1.0)``
    short-circuits to the pre-slice fast path.
    """

    model: str
    s5f_model_name: str
    # Exactly one of the two below is non-None.
    count_pairs: Optional[Tuple[Tuple[int, float], ...]] = None
    rate: Optional[float] = None
    # Normalised segment-rate vector (v, d, j, np). Default flat
    # rates so legacy pipelines stay byte-identical.
    segment_rates: Tuple[float, float, float, float] = (1.0, 1.0, 1.0, 1.0)
    # Normalised V-subregion rate vector
    # (fwr1, cdr1, fwr2, cdr2, fwr3). Default flat rates so legacy
    # pipelines stay byte-identical.
    v_subregion_rates: Tuple[float, float, float, float, float] = (
        1.0,
        1.0,
        1.0,
        1.0,
        1.0,
    )


@dataclass(frozen=True)
class _ClonalForkStep:
    """Marks the boundary between "per-clone" passes (run once per
    clonal family — typically `recombine`) and "per-descendant" passes
    (run once per read inside the family — typically `mutate` and the
    library-prep / sequencing-stage steps).

    The compile pipeline splits the experiment's step list at this
    marker, builds two `CompiledExperiment`s, and the runtime
    forks the parent IR into descendants for each clone.
    """

    n_clones: int
    size: int


@dataclass(frozen=True)
class _RepertoireForkStep:
    """Marks a non-tree clonal-repertoire fork.

    Generalizes :class:`_ClonalForkStep`: instead of a fixed
    ``per_clone`` size, each clone draws a size from a heavy-tailed
    distribution (with an unexpanded-singleton fraction). That many
    reads pass through the post-fork library-prep / sequencing passes
    and identical reads collapse into AIRR records carrying a standard
    ``duplicate_count``.

    Causes :meth:`Experiment.compile` to return a
    :class:`~GenAIRR._compiled.CompiledRepertoireExperiment`.
    """

    n_clones: int
    size_distribution: str
    exponent: float
    mu: float
    sigma: float
    max_size: int
    unexpanded_fraction: float


@dataclass(frozen=True)
class _LineageForkStep:
    """Marks an affinity-maturation lineage fork.

    Causes :meth:`Experiment.compile` to return a
    :class:`~GenAIRR._compiled.CompiledLineageExperiment` that grows BCR
    clonal lineage trees via the Rust ``simulate_family_outcomes`` kernel
    and returns per-observed-node AIRR records with lineage metadata.
    """

    n_clones: int
    max_generations: int
    n_max: int
    n_sample: int
    rate: float
    lambda_base: float
    selection_strength: float
    beta: float
    target_aa: Optional[str]
    mature_substitutions: int
    s5f_model: str
    allow_extinction: bool = False


@dataclass(frozen=True)
class _PairedEndStep:
    """One ``paired_end(r1_length=…, r2_length=…, insert_size=…)``
    invocation, captured for later compilation.

    Models the Illumina paired-end read layout (Slice D of the
    paired-end roadmap). The step holds three integer distributions
    — R1 read length, R2 read length, fragment insert size — and
    lowers into one
    :class:`~GenAIRR._engine.PassPlan.push_paired_end` call that
    constructs a trace-only
    :class:`~GenAIRR._engine.PairedEndSamplingPass`. The AIRR
    builder reads the three resulting trace records and applies the
    projection kernel landed in Slice B.

    Pipeline-IR position: the step lowers at the **end** of the
    plan, after every IR-mutating + observation-stage pass
    (recombine, mutation, corruption, end-loss, rev-comp,
    receptor revision, D inversion). Even though the pass is
    trace-only and does not mutate the simulation, lowering it
    last keeps the trace record order aligned with the
    biological/readout order — the canonical "biology → corruption
    → readout" ordering the design doc §6 requires.

    Valid on both VDJ and VJ chains (paired-end is a sequencing-
    stage observable, not a biology mechanism). :meth:`Experiment.paired_end`
    rejects duplicates at the DSL boundary.
    """

    # Each is a tuple of (int, float) pairs (the canonical
    # empirical-distribution shape used elsewhere in the IR).
    r1_length: Tuple[Tuple[int, float], ...]
    r2_length: Tuple[Tuple[int, float], ...]
    insert_size: Tuple[Tuple[int, float], ...]


@dataclass(frozen=True)
class _ReceptorRevisionStep:
    """One ``receptor_revision(prob=...)`` invocation, captured for
    later compilation.

    Models post-recombination V-segment replacement (Slice C of the
    receptor-revision roadmap): with probability ``prob`` the V slot
    is reassigned to a different germline V allele and the V slice
    in the pool is rewritten via
    :class:`~GenAIRR._engine.ReceptorRevisionPass`. The pass records
    a Bool at ``receptor_revision.applied`` for every simulation; on
    ``true`` it additionally records the replacement allele id at
    ``receptor_revision.v_allele`` and the derived 3' trim at
    ``receptor_revision.v_trim_3``.

    Pipeline-IR position: the step lowers into a single
    ``push_receptor_revision(prob, refdata)`` placed immediately
    after ``push_assemble("J")`` in the canonical V-NP1-D-NP2-J
    recombine sequence, so receptor revision sees the fully-
    assembled V/D/J/NP pool. Subsequent ``_MutateStep`` /
    ``_CorruptStep`` lowering pushes the post-recombination passes
    after this one in the plan, preserving the
    "recombine → revise → mutate/corrupt" order the design doc §2
    requires.

    Valid only for VDJ chains; VJ chains have no D pool and the
    biology of receptor revision is heavy-chain-only in v1.
    :meth:`Experiment.receptor_revision` rejects VJ at the DSL
    boundary before this dataclass is constructed.

    ``same_haplotype`` (default ``True``) only matters when a genotype
    is attached: the replacement V is restricted to the carried V
    alleles on the drawn rearrangement chromosome (the cis VH-
    replacement model). ``False`` is a synthetic control that draws
    from both chromosomes. Ignored on the no-genotype path.
    """

    prob: float
    same_haplotype: bool = True


@dataclass(frozen=True)
class _InvertDStep:
    """One ``invert_d(prob=...)`` invocation, captured for later
    compilation.

    Models V(D)J inversion: the D allele is committed in
    ``SegmentOrientation::ReverseComplement`` with probability
    ``prob``, otherwise ``Forward``. The Rust pass
    (``InvertDPass``) records a per-simulation Bool at
    ``sample_allele.d.inverted`` and updates the D
    ``AlleleInstance`` accordingly; Slice B's assembly emission
    consumes the orientation flag.

    Pipeline-IR position: the step must sit *after* the user's
    :class:`_RecombineStep` (lowering depends on the
    ``assemble.d`` pass already existing in the plan so the
    explicit ``before(invert_d, assemble.d)`` schedule edge can
    fire). :class:`Experiment.invert_d` enforces this by
    appending the step to ``self._steps`` directly.

    Valid only for VDJ chains; VJ chains have no D pool and
    rejected at the :meth:`Experiment.invert_d` boundary before
    a step is constructed.
    """

    prob: float


@dataclass(frozen=True)
class _CorruptStep:
    """One library-prep / sequencing-stage step, captured for later compilation.

    ``kind`` selects which Rust corruption pass to construct:
    ``"pcr"``, ``"quality"``, ``"indel"``, or ``"contaminant"``.
    Other fields are pass-specific:

    - PCR / quality each use exactly one of ``count_pairs`` *or*
      ``rate``; when ``rate`` is set the engine draws
      ``count ~ Poisson(rate × pool_len)`` per record.
    - Indel uses ``count_pairs`` + ``insertion_prob``.
    - Contaminant / reverse-complement use ``apply_prob``.
    """

    kind: str
    count_pairs: Optional[Tuple[Tuple[int, float], ...]] = None
    insertion_prob: float = 0.0
    apply_prob: float = 0.0
    rate: Optional[float] = None
