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
    """

    model: str
    s5f_model_name: str
    # Exactly one of the two below is non-None.
    count_pairs: Optional[Tuple[Tuple[int, float], ...]] = None
    rate: Optional[float] = None


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
