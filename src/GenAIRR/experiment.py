"""GenAIRR — fluent DSL for receptor-sequence simulation.

The :class:`Experiment` builder compiles a chain of fluent steps into
a Rust ``PassPlan`` (via :mod:`genairr_engine`) and runs it to
produce a list of :class:`genairr_engine.Outcome` objects.

Typical usage::

    import GenAIRR as ga
    outcomes = ga.Experiment.on("human_igh").recombine().run(n=100, seed=42)

``Experiment.on`` accepts:

- a config-name string (``"human_igh"``, ``"mouse_tcrb"``, …) — resolved
  through the builtin :mod:`GenAIRR.data` registry,
- a :class:`GenAIRR.DataConfig` object (already-loaded reference data),
- a :class:`genairr_engine.RefDataConfig` (the engine-native form, used
  primarily by tests and advanced callers).

In all three cases the input is normalised to a
:class:`genairr_engine.RefDataConfig` before any pass is appended.

This file replaces the V5 ``protocol.py`` + ``experiment.py`` pair.
The V5 fluent surface (``.mutate(...)``, ``.sequence(...)``,
``.observe(...)``) has not yet been ported to V6 — only
:meth:`Experiment.recombine` is wired today. Subsequent phases of
the V6 migration add mutation / corruption / observation steps.
"""
from __future__ import annotations

from dataclasses import dataclass
from typing import Any, Dict, Iterable, Iterator, List, Optional, Sequence, Tuple, Union

import genairr_engine as _engine

from .dataconfig import DataConfig
from .dataconfig.enums import ChainType


# Anything :meth:`Experiment.run` (or :meth:`CompiledExperiment.run`)
# accepts as ``respect=``. The Rust runner ultimately needs a single
# :class:`genairr_engine.ContractSet`; the helpers below normalise
# the convenient list-form ``[productive()]`` (V5 muscle-memory) to
# that single value.
RespectInput = Union["_engine.ContractSet", Sequence["_engine.ContractSet"], None]


# Sentinel for `using(...)` arguments that the caller did not pass at
# all. We can't use ``None`` for "unchanged" because ``None`` already
# means "clear the lock."
class _Unset:
    __slots__ = ()

    def __repr__(self) -> str:
        return "<UNSET>"


_UNSET: _Unset = _Unset()


# Inputs accepted by :meth:`Experiment.using` per segment. ``None``
# clears any prior lock; a single name is a one-allele lock; an
# iterable of names is a multi-allele subset; ``_UNSET`` (the default)
# leaves the existing lock unchanged.
_LockInput = Union[str, Iterable[str], None, _Unset]


# ──────────────────────────────────────────────────────────────────
# Config-name resolver — string → DataConfig
# ──────────────────────────────────────────────────────────────────

# Short aliases (preferred defaults). Heavy / light / TCR for human
# fall back to OGRDB when present, otherwise IMGT V-QUEST.
_CONFIG_ALIASES = {
    # Human (OGRDB preferred for IG, IMGT for TCR)
    "human_igh": "HUMAN_IGH_OGRDB",
    "human_igk": "HUMAN_IGK_OGRDB",
    "human_igl": "HUMAN_IGL_OGRDB",
    "human_tcra": "HUMAN_TCRA_IMGT",
    "human_tcrb": "HUMAN_TCRB_IMGT",
    "human_tcrd": "HUMAN_TCRD_IMGT",
    "human_tcrg": "HUMAN_TCRG_IMGT",
    # Mouse
    "mouse_igh": "MOUSE_IGH_IMGT",
    "mouse_igk": "MOUSE_IGK_IMGT",
    "mouse_igl": "MOUSE_IGL_IMGT",
    "mouse_tcra": "MOUSE_TCRA_IMGT",
    "mouse_tcrb": "MOUSE_TCRB_IMGT",
    "mouse_tcrd": "MOUSE_TCRD_IMGT",
    "mouse_tcrg": "MOUSE_TCRG_IMGT",
    # Rat
    "rat_igh": "RAT_IGH_IMGT",
    "rat_igk": "RAT_IGK_IMGT",
    "rat_igl": "RAT_IGL_IMGT",
    # Rabbit
    "rabbit_igh": "RABBIT_IGH_IMGT",
    "rabbit_igk": "RABBIT_IGK_IMGT",
    "rabbit_igl": "RABBIT_IGL_IMGT",
    "rabbit_tcrb": "RABBIT_TCRB_IMGT",
    # Rhesus macaque
    "rhesus_igh": "RHESUS_IGH_IMGT",
    "rhesus_igk": "RHESUS_IGK_IMGT",
    "rhesus_igl": "RHESUS_IGL_IMGT",
    "rhesus_tcrb": "RHESUS_TCRB_IMGT",
    # Cow
    "cow_igh": "COW_IGH_IMGT",
    "cow_igk": "COW_IGK_IMGT",
    "cow_igl": "COW_IGL_IMGT",
    "cow_tcrb": "COW_TCRB_IMGT",
    # Dog
    "dog_igh": "DOG_IGH_IMGT",
    "dog_igk": "DOG_IGK_IMGT",
    "dog_igl": "DOG_IGL_IMGT",
    "dog_tcrb": "DOG_TCRB_IMGT",
    # Cat
    "cat_igk": "CAT_IGK_IMGT",
    "cat_igl": "CAT_IGL_IMGT",
    "cat_tcrb": "CAT_TCRB_IMGT",
    # Pig
    "pig_igh": "PIG_IGH_IMGT",
    "pig_igk": "PIG_IGK_IMGT",
    "pig_igl": "PIG_IGL_IMGT",
    "pig_tcrb": "PIG_TCRB_IMGT",
}


def _build_full_aliases() -> None:
    """Auto-register every builtin DataConfig under both its full and
    short ``_imgt`` / ``_ogrdb``-stripped lowercase names."""
    from .data import _CONFIG_NAMES

    for name in _CONFIG_NAMES:
        lower = name.lower()
        if lower not in _CONFIG_ALIASES:
            _CONFIG_ALIASES[lower] = name
        for suffix in ("_imgt", "_ogrdb"):
            if lower.endswith(suffix):
                short = lower[: -len(suffix)]
                if short not in _CONFIG_ALIASES:
                    _CONFIG_ALIASES[short] = name


_build_full_aliases()


def _resolve_config_name(name: str) -> DataConfig:
    """Resolve a config-name string to its loaded :class:`DataConfig`.

    Hyphens and case are normalised. Unknown names raise ``ValueError``
    with a hint.
    """
    key = name.lower().replace("-", "_")
    const_name = _CONFIG_ALIASES.get(key)
    if const_name is None:
        from .data import list_configs

        available = list_configs()
        raise ValueError(
            f"Unknown config name {name!r}. "
            f"Available configs ({len(available)} total): "
            f"{', '.join(a.lower() for a in available[:10])}... "
            f"Use GenAIRR.list_configs() for the full list, "
            f"or pass a DataConfig object directly."
        )
    from . import data as _data

    return getattr(_data, const_name)


# ──────────────────────────────────────────────────────────────────
# DataConfig → genairr_engine.RefDataConfig translator
# ──────────────────────────────────────────────────────────────────


def _chain_type_label(chain_type: Optional[ChainType]) -> str:
    """Map a :class:`GenAIRR.dataconfig.enums.ChainType` to the
    lowercase ``"vj"`` / ``"vdj"`` label expected by ``RefDataConfig``.

    Heavy chain (BCR_HEAVY) and TCR β / δ have a D segment → ``"vdj"``.
    Everything else → ``"vj"``. ``None`` (older pickles missing
    metadata) defaults to ``"vdj"`` so the more general builder is
    selected.
    """
    if chain_type is None:
        return "vdj"
    if chain_type in (
        ChainType.BCR_HEAVY,
        ChainType.TCR_BETA,
        ChainType.TCR_DELTA,
    ):
        return "vdj"
    return "vj"


def _push_alleles(alleles_by_gene, add_method) -> int:
    """Walk ``Dict[str, List[Allele]]`` and push every allele through
    ``add_method``. Returns the count of pushed alleles.
    """
    count = 0
    if not alleles_by_gene:
        return count
    for _gene, allele_list in alleles_by_gene.items():
        for allele in allele_list:
            seq_str = allele.ungapped_seq
            seq = seq_str.encode("ascii") if isinstance(seq_str, str) else bytes(seq_str)
            anchor = allele.anchor
            if anchor is None:
                add_method(allele.name, allele.gene, seq)
            else:
                add_method(allele.name, allele.gene, seq, anchor=int(anchor))
            count += 1
    return count


def dataconfig_to_refdata(cfg: DataConfig) -> "_engine.RefDataConfig":
    """Translate a :class:`DataConfig` (V5-style species pickle) into
    an engine-native :class:`genairr_engine.RefDataConfig`.

    All V / D / J alleles are copied; the C segment is dropped (V6 has
    no C-segment passes yet). Anchorless alleles are preserved with
    ``anchor=None`` — they pass through but cannot satisfy the
    ``AnchorPreserved`` contract (matches V5's behaviour).
    """
    chain = _chain_type_label(cfg.metadata.chain_type if cfg.metadata else None)
    refdata = _engine.RefDataConfig(chain)

    _push_alleles(cfg.v_alleles, refdata.add_v_allele)
    if chain == "vdj":
        _push_alleles(cfg.d_alleles, refdata.add_d_allele)
    _push_alleles(cfg.j_alleles, refdata.add_j_allele)
    return refdata


# ──────────────────────────────────────────────────────────────────
# Build steps
# ──────────────────────────────────────────────────────────────────

# Default NP-length distribution: lengths 0..6 with equal weight.
_DEFAULT_NP_LENGTHS: Sequence[Tuple[int, float]] = tuple((i, 1.0) for i in range(7))


@dataclass(frozen=True)
class _RecombineStep:
    """One ``recombine()`` invocation, captured for later compilation.

    ``np1_lengths`` / ``np2_lengths`` are the empirical NP-length
    distributions to use. ``trim_*`` are optional empirical trim
    distributions; when ``None`` the corresponding trim pass is
    omitted (so a synthetic refdata without trim data still works).

    ``locks_*`` are optional allele-ID subsets supplied by
    :meth:`Experiment.using`. When set, the matching
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

    def apply(
        self,
        plan: "_engine.PassPlan",
        refdata: "_engine.RefDataConfig",
    ) -> None:
        chain = refdata.chain_type
        np1 = list(self.np1_lengths)
        np2 = list(self.np2_lengths)

        v_ids = list(self.locks_v) if self.locks_v is not None else None
        d_ids = list(self.locks_d) if self.locks_d is not None else None
        j_ids = list(self.locks_j) if self.locks_j is not None else None

        if chain == "vj":
            plan.push_sample_allele("V", refdata, allowed_ids=v_ids)
            plan.push_sample_allele("J", refdata, allowed_ids=j_ids)
            if self.trim_v_3:
                plan.push_trim("V", "3", list(self.trim_v_3))
            if self.trim_j_5:
                plan.push_trim("J", "5", list(self.trim_j_5))
            plan.push_assemble("V")
            plan.push_generate_np("NP1", np1)
            plan.push_assemble("J")
        elif chain == "vdj":
            plan.push_sample_allele("V", refdata, allowed_ids=v_ids)
            plan.push_sample_allele("D", refdata, allowed_ids=d_ids)
            plan.push_sample_allele("J", refdata, allowed_ids=j_ids)
            if self.trim_v_3:
                plan.push_trim("V", "3", list(self.trim_v_3))
            if self.trim_d_5:
                plan.push_trim("D", "5", list(self.trim_d_5))
            if self.trim_d_3:
                plan.push_trim("D", "3", list(self.trim_d_3))
            if self.trim_j_5:
                plan.push_trim("J", "5", list(self.trim_j_5))
            plan.push_assemble("V")
            plan.push_generate_np("NP1", np1)
            plan.push_assemble("D")
            plan.push_generate_np("NP2", np2)
            plan.push_assemble("J")
        else:  # pragma: no cover — RefDataConfig only constructs vj/vdj.
            raise ValueError(f"unsupported chain_type {chain!r}")


@dataclass(frozen=True)
class _MutateStep:
    """One ``mutate()`` invocation, captured for later compilation.

    ``model`` is either ``"uniform"`` (position-independent SHM via
    ``UniformMutationPass``) or ``"s5f"`` (context-dependent SHM via
    ``S5FMutationPass`` with the bundled S5F kernel named in
    ``s5f_model_name``).

    ``count_pairs`` is the empirical distribution over the number of
    mutations per simulation, expressed as ``((count, weight), ...)``.
    """

    model: str
    count_pairs: Tuple[Tuple[int, float], ...]
    s5f_model_name: str

    def apply(
        self,
        plan: "_engine.PassPlan",
        refdata: "_engine.RefDataConfig",
    ) -> None:
        del refdata  # not needed for mutation passes
        count_list = list(self.count_pairs)
        if self.model == "uniform":
            plan.push_mutate_uniform(count_list)
        elif self.model == "s5f":
            from ._s5f_loader import load_builtin_s5f_kernel

            mutability, substitution = load_builtin_s5f_kernel(self.s5f_model_name)
            plan.push_mutate_s5f(count_list, mutability, substitution)
        else:  # pragma: no cover — guarded at builder time.
            raise ValueError(f"unsupported mutation model {self.model!r}")


# Names of the four corruption passes the Rust kernel exposes. Used
# both as a discriminator on `_CorruptStep` and as the Python-facing
# kwarg name on the `Experiment` builder.
_CORRUPT_KIND_PCR = "pcr"
_CORRUPT_KIND_QUALITY = "quality"
_CORRUPT_KIND_INDEL = "indel"
_CORRUPT_KIND_CONTAMINANT = "contaminant"


@dataclass(frozen=True)
class _CorruptStep:
    """One ``corrupt_*()`` invocation, captured for later compilation.

    ``kind`` selects which Rust corruption pass to construct:
    ``"pcr"``, ``"quality"``, ``"indel"``, or ``"contaminant"``.
    Other fields are pass-specific:

    - PCR / quality / indel use ``count_pairs``.
    - Indel additionally uses ``insertion_prob``.
    - Contaminant uses ``apply_prob`` only (its model is
      probability-driven at the read level, not count-driven).
    """

    kind: str
    count_pairs: Tuple[Tuple[int, float], ...]
    insertion_prob: float
    apply_prob: float

    def apply(
        self,
        plan: "_engine.PassPlan",
        refdata: "_engine.RefDataConfig",
    ) -> None:
        del refdata  # not used by corruption passes
        count_list = list(self.count_pairs)
        if self.kind == _CORRUPT_KIND_PCR:
            plan.push_corrupt_pcr(count_list)
        elif self.kind == _CORRUPT_KIND_QUALITY:
            plan.push_corrupt_quality(count_list)
        elif self.kind == _CORRUPT_KIND_INDEL:
            plan.push_corrupt_indel(count_list, insertion_prob=self.insertion_prob)
        elif self.kind == _CORRUPT_KIND_CONTAMINANT:
            plan.push_corrupt_contaminant(self.apply_prob)
        else:  # pragma: no cover — guarded at builder time.
            raise ValueError(f"unsupported corruption kind {self.kind!r}")


def _normalize_count(
    count: Union[int, Tuple[int, int], Iterable[Tuple[int, float]]],
) -> Tuple[Tuple[int, float], ...]:
    """Normalise the user-friendly ``count=`` argument to the
    list-of-pairs shape ``PyPassPlan.push_mutate_*`` expects.

    Accepts:
    - ``count=15`` — fixed count.
    - ``count=(5, 25)`` — uniform integer in ``[5, 25]``
      (high inclusive, matches the Pythonic ``range``-with-stop
      conventions; both endpoints are valid sample values).
    - ``count=[(5, 1.0), (10, 1.0), (15, 1.0)]`` — explicit
      empirical distribution.
    """
    # Single int: fixed count. Reject bools (which are int subclasses)
    # so a stray True/False raises TypeError instead of silently
    # becoming ``count=1`` / ``count=0``.
    if isinstance(count, bool) or not isinstance(count, (int, tuple, list)):
        raise TypeError(
            f"count: expected int, (low, high) tuple, or list of "
            f"(count, weight) pairs, got {type(count).__name__}"
        )
    if isinstance(count, int):
        if count < 0:
            raise ValueError(f"count must be non-negative, got {count}")
        return ((int(count), 1.0),)
    # Tuple form: must be (low, high) of two ints OR a single
    # (count, weight) pair we promote to a 1-element list.
    if isinstance(count, tuple) and len(count) == 2:
        a, b = count[0], count[1]
        if (
            isinstance(a, int)
            and isinstance(b, int)
            and not isinstance(a, bool)
            and not isinstance(b, bool)
        ):
            low, high = a, b
            if low < 0 or high < low:
                raise ValueError(
                    f"count range must satisfy 0 <= low <= high, got ({low}, {high})"
                )
            return tuple((c, 1.0) for c in range(low, high + 1))
    # Otherwise: empirical (count, weight) list/tuple.
    pairs: List[Tuple[int, float]] = []
    for pair in count:
        if not (isinstance(pair, (tuple, list)) and len(pair) == 2):
            raise TypeError(
                f"count entries must be (count, weight) pairs, got {pair!r}"
            )
        c, w = pair
        if isinstance(c, bool) or not isinstance(c, int) or c < 0:
            raise ValueError(
                f"count entry must have a non-negative int count, got {c!r}"
            )
        pairs.append((int(c), float(w)))
    if not pairs:
        raise ValueError("count distribution must contain at least one entry")
    return tuple(pairs)


def _to_immutable_pairs(
    pairs: Optional[List[Tuple[int, float]]],
) -> Optional[Tuple[Tuple[int, float], ...]]:
    """Convert an optional ``list[(int, float)]`` (the shape returned
    by :mod:`._dataconfig_extract`) into the hashable tuple-of-tuples
    form that :class:`_RecombineStep`'s frozen dataclass requires.
    ``None`` passes through unchanged so the step can omit the trim
    pass entirely when no empirical distribution is available.
    """
    if not pairs:
        return None
    return tuple((int(c), float(w)) for c, w in pairs)


def _normalize_lengths(
    lengths: Optional[Iterable[Tuple[int, float]]],
) -> Tuple[Tuple[int, float], ...]:
    """Convert a user-supplied length iterable into a hashable tuple of
    ``(int, float)`` pairs. ``None`` falls back to the module default;
    empty iterables raise ``ValueError`` here so the error surfaces at
    builder time."""
    if lengths is None:
        return tuple(_DEFAULT_NP_LENGTHS)
    pairs: List[Tuple[int, float]] = []
    for pair in lengths:
        length, weight = pair
        pairs.append((int(length), float(weight)))
    if not pairs:
        raise ValueError(
            "length distribution must contain at least one (length, weight) pair"
        )
    return tuple(pairs)


# ──────────────────────────────────────────────────────────────────
# Public API: Experiment + CompiledExperiment
# ──────────────────────────────────────────────────────────────────

def _coerce_respect(respect: RespectInput) -> Optional["_engine.ContractSet"]:
    """Normalise the ``respect=`` argument to a single ``ContractSet``.

    Accepts:
    - ``None`` → no contracts active.
    - a single ``ContractSet`` → used directly.
    - a length-1 sequence containing one ``ContractSet`` → unwrapped
      (the V5 ``respect=[productive()]`` muscle-memory case).

    Sequences with multiple bundles raise ``NotImplementedError`` —
    contract composition will land in a later phase. Anything else
    raises ``TypeError``.
    """
    if respect is None:
        return None
    if isinstance(respect, _engine.ContractSet):
        return respect
    if isinstance(respect, (list, tuple)):
        if len(respect) == 0:
            return None
        if len(respect) == 1:
            item = respect[0]
            if not isinstance(item, _engine.ContractSet):
                raise TypeError(
                    f"respect[0]: expected genairr_engine.ContractSet, "
                    f"got {type(item).__name__}"
                )
            return item
        raise NotImplementedError(
            f"respect=[...] with {len(respect)} bundles is not yet supported. "
            "Compose contracts inside a single bundle for now."
        )
    raise TypeError(
        f"respect=: expected genairr_engine.ContractSet or sequence of one, "
        f"got {type(respect).__name__}"
    )


# Anything :meth:`Experiment.on` accepts.
ExperimentInput = Union[str, DataConfig, "_engine.RefDataConfig"]


def _coerce_to_refdata_and_dataconfig(
    source: ExperimentInput,
) -> Tuple["_engine.RefDataConfig", Optional[DataConfig]]:
    """Normalise ``Experiment.on`` input to ``(refdata, dataconfig)``.

    The DataConfig — when available — carries the empirical
    distributions (NP lengths, per-gene trims) that ``recombine()``
    uses by default. A bare ``RefDataConfig`` from the user has no
    DataConfig backing, so empirical distributions fall through to
    the uniform ``[0..6]`` placeholder.
    """
    if isinstance(source, _engine.RefDataConfig):
        return source, None
    if isinstance(source, str):
        cfg = _resolve_config_name(source)
        return dataconfig_to_refdata(cfg), cfg
    if isinstance(source, DataConfig):
        return dataconfig_to_refdata(source), source
    raise TypeError(
        f"Experiment.on(...): expected a config-name string, "
        f"GenAIRR.DataConfig, or genairr_engine.RefDataConfig, "
        f"got {type(source).__name__}"
    )


class Experiment:
    """Fluent builder for a simulation pipeline.

    Build an ``Experiment`` from a config name, a :class:`DataConfig`,
    or a :class:`genairr_engine.RefDataConfig` via :meth:`Experiment.on`,
    chain configuration steps (currently just :meth:`recombine`), then
    call :meth:`run` (or :meth:`compile` for explicit two-stage flow).

    The builder is **stateful but not destructive**: each fluent call
    returns ``self`` after appending a step. The same ``Experiment``
    can be ``compile()``-d and ``run()``-d multiple times with
    different seeds.
    """

    __slots__ = ("_refdata", "_steps", "_dataconfig", "_locks")

    def __init__(
        self,
        refdata: "_engine.RefDataConfig",
        dataconfig: Optional[DataConfig] = None,
    ) -> None:
        self._refdata = refdata
        self._dataconfig = dataconfig
        self._steps: List[Union[_RecombineStep, _MutateStep, _CorruptStep]] = []
        # Per-segment allele-lock subsets set by ``.using(...)``. Each
        # entry is ``None`` (no lock — sample uniformly across the pool)
        # or a tuple of allele IDs to sample uniformly across.
        self._locks: Dict[str, Optional[Tuple[int, ...]]] = {
            "V": None,
            "D": None,
            "J": None,
        }

    @classmethod
    def on(cls, source: ExperimentInput) -> "Experiment":
        """Start an experiment against the given reference data.

        ``source`` is one of:
        - a config-name string (e.g. ``"human_igh"``),
        - a :class:`GenAIRR.DataConfig` instance,
        - a :class:`genairr_engine.RefDataConfig`.

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

    def corrupt_pcr(
        self,
        *,
        count: Union[int, Tuple[int, int], Iterable[Tuple[int, float]]],
    ) -> "Experiment":
        """Append a PCR-error corruption step.

        ``count`` is the per-simulation distribution over the number
        of PCR-induced base substitutions. Same input forms as
        :meth:`mutate`'s ``count`` (fixed int, ``(low, high)`` range,
        or empirical ``(count, weight)`` list).

        Each error samples a uniform position in the assembled pool
        and replaces the base with a uniform A/C/G/T draw. The trace
        records the events at ``corrupt.pcr.{count, error_site[i],
        error_base[i]}``.
        """
        pairs = _normalize_count(count)
        self._steps.append(
            _CorruptStep(
                kind=_CORRUPT_KIND_PCR,
                count_pairs=pairs,
                insertion_prob=0.0,
                apply_prob=0.0,
            )
        )
        return self

    def corrupt_quality(
        self,
        *,
        count: Union[int, Tuple[int, int], Iterable[Tuple[int, float]]],
    ) -> "Experiment":
        """Append a sequencing-quality-error corruption step.

        Same shape as :meth:`corrupt_pcr` but each substitution
        writes the destination base in **lowercase** to mark the
        position as corrupted (the V5 sequencing-error convention).
        """
        pairs = _normalize_count(count)
        self._steps.append(
            _CorruptStep(
                kind=_CORRUPT_KIND_QUALITY,
                count_pairs=pairs,
                insertion_prob=0.0,
                apply_prob=0.0,
            )
        )
        return self

    def corrupt_indels(
        self,
        *,
        count: Union[int, Tuple[int, int], Iterable[Tuple[int, float]]],
        insertion_prob: float = 0.5,
    ) -> "Experiment":
        """Append an indel corruption step.

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

    def corrupt_contaminants(self, *, prob: float) -> "Experiment":
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

    def mutate(
        self,
        *,
        model: str = "s5f",
        count: Union[int, Tuple[int, int], Iterable[Tuple[int, float]]],
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

        ``count`` is the per-simulation mutation count distribution.
        Accepts:
        - ``count=15`` — fixed: every simulation gets exactly 15
          mutations.
        - ``count=(5, 25)`` — uniform integer in ``[5, 25]`` (both
          endpoints inclusive).
        - ``count=[(5, 1.0), (10, 2.0), ...]`` — explicit empirical
          ``(count, weight)`` distribution.

        ``count=0`` is a no-op (the pass runs but applies zero
        mutations).
        """
        model_lc = model.lower()
        if model_lc not in ("uniform", "s5f"):
            raise ValueError(
                f"model must be 'uniform' or 's5f' (got {model!r})"
            )
        pairs = _normalize_count(count)
        self._steps.append(
            _MutateStep(
                model=model_lc,
                count_pairs=pairs,
                s5f_model_name=s5f_model,
            )
        )
        return self

    def using(
        self,
        *,
        v: "_LockInput" = _UNSET,
        d: "_LockInput" = _UNSET,
        j: "_LockInput" = _UNSET,
    ) -> "Experiment":
        """Restrict allele sampling to a named subset (allele-locking).

        Each kwarg accepts:
        - a single allele name string (e.g. ``"IGHV1-2*02"``),
        - a list / tuple of allele names (sample uniformly among them),
        - ``None`` to clear a previously-set lock for that segment.
        - omitted (default) — the existing lock for that segment is
          left unchanged.

        The locks are applied to the next ``recombine()`` step's
        ``push_sample_allele`` calls at compile time. Calling
        ``.using(...)`` more than once overlays the new values onto
        the previous locks (per-segment), so
        ``.using(v="A").using(d="B")`` locks both V and D.

        Raises:
        - ``ValueError`` if an allele name doesn't exist in the
          configured refdata, or if a lock is set for ``D`` on a VJ
          chain.
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
                    f"using(): {segment} lock must be a name string or an "
                    f"iterable of name strings; got {type(value).__name__}"
                ) from exc
            if not all(isinstance(n, str) for n in names):
                raise TypeError(
                    f"using(): {segment} lock entries must all be strings"
                )
            if not names:
                raise ValueError(
                    f"using(): {segment} lock list must be non-empty "
                    "(pass None to clear instead)"
                )

        index = self._allele_name_index(segment)
        ids: List[int] = []
        seen: set = set()
        for name in names:
            if name not in index:
                raise ValueError(
                    f"using(): no {segment} allele named {name!r} in refdata "
                    "(check spelling against `list_alleles` or the loaded config)"
                )
            allele_id = index[name]
            if allele_id in seen:
                raise ValueError(
                    f"using(): duplicate {segment} allele {name!r} in lock list"
                )
            seen.add(allele_id)
            ids.append(allele_id)
        return tuple(ids)

    def _allele_name_index(self, segment: str) -> Dict[str, int]:
        """Build a name → allele_id map for the given segment by
        scanning the refdata pool. Cheap enough to do per-call given
        typical pool sizes (≤ a few hundred) and the rarity of
        ``.using()``."""
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
        trim: bool = True,
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
        ``[(0, 1.0), ..., (6, 1.0)]`` distribution. Pass an explicit
        iterable of ``(length, weight)`` tuples to override the
        default. ``np2_lengths`` is silently ignored on VJ chains.

        ``trim=True`` (default) inserts trim passes before assembly
        when the bound DataConfig carries empirical trim
        distributions. Trim distributions are marginalized from the
        per-gene ``trim_dicts`` to a single segment-level
        distribution. Set ``trim=False`` to disable trimming
        entirely (e.g. for tests that expect untrimmed alleles).
        Raw-RefDataConfig experiments always have ``trim`` as a
        no-op since there's no DataConfig to source the
        distributions from.
        """
        defaults = self._recombine_defaults() if trim or self._dataconfig else None

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

        step = _RecombineStep(
            np1_lengths=np1,
            np2_lengths=np2,
            trim_v_3=trim_v_3,
            trim_d_5=trim_d_5,
            trim_d_3=trim_d_3,
            trim_j_5=trim_j_5,
        )
        self._steps.append(step)
        return self

    def _recombine_defaults(self):
        """Lazy-extract the empirical distributions for this
        experiment's DataConfig. Returns ``None`` for raw-RefDataConfig
        experiments. Cached on first call.
        """
        if self._dataconfig is None:
            return None
        from ._dataconfig_extract import extract_recombine_defaults

        return extract_recombine_defaults(self._dataconfig)

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

    def compile(self) -> "CompiledExperiment":
        """Compile the recorded steps into a ``PassPlan`` and return
        a :class:`CompiledExperiment` ready to be ``run()``.

        Idempotent: calling ``compile()`` twice produces two distinct
        ``CompiledExperiment`` instances with structurally-equal plans.
        """
        from dataclasses import replace as _replace

        any_lock = any(self._locks[seg] is not None for seg in ("V", "D", "J"))
        plan = _engine.PassPlan()
        for step in self._steps:
            # Inject any allele-locks set via ``.using(...)`` into the
            # recombine step at compile time. Other step types ignore
            # locks.
            if any_lock and isinstance(step, _RecombineStep):
                step = _replace(
                    step,
                    locks_v=self._locks["V"],
                    locks_d=self._locks["D"],
                    locks_j=self._locks["J"],
                )
            step.apply(plan, self._refdata)
        return CompiledExperiment(plan, self._refdata)

    def run_records(
        self,
        *,
        n: int = 1,
        seed: int = 0,
        respect: RespectInput = None,
        strict: bool = False,
    ) -> "SimulationResult":
        """Compile and run, then return the batch as a
        :class:`SimulationResult` ready for ``.to_csv`` / ``.to_fasta``
        / ``.to_dataframe`` export.

        Same arguments as :meth:`run`. Returns a
        :class:`GenAIRR.SimulationResult` instead of a raw list of
        ``Outcome`` objects.
        """
        from .result import SimulationResult

        outcomes = self.run(n=n, seed=seed, respect=respect, strict=strict)
        return SimulationResult.from_outcomes(outcomes, self._refdata)

    def run(
        self,
        *,
        n: int = 1,
        seed: int = 0,
        respect: RespectInput = None,
        strict: bool = False,
    ) -> List["_engine.Outcome"]:
        """Compile and run this experiment ``n`` times.

        Equivalent to
        ``self.compile().run(n=n, seed=seed, respect=respect, strict=strict)``.
        Returns a list of :class:`genairr_engine.Outcome` objects.

        Pass ``respect=productive()`` (or any other ``ContractSet``) to
        constrain every random draw to admissible values; the runtime
        filters NP base draws, length samples, and mutation /
        contamination substitutions in real time so the resulting
        sequences satisfy the bundle by construction.

        ``strict=False`` (default) lets a pass fall back to
        unconstrained sampling when no admissible candidate exists.
        ``strict=True`` raises
        :class:`genairr_engine.StrictSamplingError` instead, so
        unsatisfiable plans surface as exceptions rather than
        silently-relaxed outputs.
        """
        return self.compile().run(n=n, seed=seed, respect=respect, strict=strict)

    def stream(
        self,
        *,
        n: Optional[int] = None,
        seed: int = 0,
        respect: RespectInput = None,
        strict: bool = False,
    ) -> Iterator["_engine.Outcome"]:
        """Compile and lazily yield :class:`genairr_engine.Outcome`
        objects. See :meth:`CompiledExperiment.stream` for full
        semantics."""
        return self.compile().stream(
            n=n, seed=seed, respect=respect, strict=strict
        )

    def stream_records(
        self,
        *,
        n: Optional[int] = None,
        seed: int = 0,
        respect: RespectInput = None,
        strict: bool = False,
        id_prefix: str = "seq",
    ) -> Iterator[Dict[str, Any]]:
        """Compile and lazily yield AIRR-format record dicts. See
        :meth:`CompiledExperiment.stream_records`."""
        return self.compile().stream_records(
            n=n,
            seed=seed,
            respect=respect,
            strict=strict,
            id_prefix=id_prefix,
        )

    def __repr__(self) -> str:
        return f"<Experiment chain={self.chain_type} steps={self.step_count}>"


class CompiledExperiment:
    """A frozen ``Experiment`` ready for execution.

    Holds the compiled :class:`genairr_engine.PassPlan` and the
    refdata it was built against. ``run()`` is the only public
    behaviour — a ``CompiledExperiment`` is intentionally a leaf
    object with no further configuration.
    """

    __slots__ = ("_plan", "_refdata")

    def __init__(
        self,
        plan: "_engine.PassPlan",
        refdata: "_engine.RefDataConfig",
    ) -> None:
        self._plan = plan
        self._refdata = refdata

    @property
    def plan(self) -> "_engine.PassPlan":
        """The compiled :class:`genairr_engine.PassPlan`."""
        return self._plan

    @property
    def refdata(self) -> "_engine.RefDataConfig":
        """The :class:`genairr_engine.RefDataConfig` the plan was built against."""
        return self._refdata

    def run(
        self,
        *,
        n: int = 1,
        seed: int = 0,
        respect: RespectInput = None,
        strict: bool = False,
    ) -> List["_engine.Outcome"]:
        """Run the compiled plan ``n`` times.

        Each iteration uses ``seed + i`` as the per-run seed so
        consecutive batches stitch together by offsetting ``seed``.

        ``respect`` accepts ``None`` (no contracts), a single
        :class:`genairr_engine.ContractSet`, or a length-1 sequence
        containing one (the V5-style ``[productive()]``). When set,
        every sampling pass filters its candidate distribution
        through the contract bundle.

        ``strict`` controls the failure mode when filtering yields
        no admissible candidate:
        - ``False`` (default, **permissive**) — fall back to
          unconstrained sampling and continue.
        - ``True`` (**strict**) — raise
          :class:`genairr_engine.StrictSamplingError` carrying the
          failing pass name, trace address, and failure reason.

        Raises ``ValueError`` for ``n < 1``.
        """
        if n < 1:
            raise ValueError(f"n must be at least 1, got {n}")
        contracts = _coerce_respect(respect)
        return [
            _engine.run(
                self._plan,
                seed=seed + i,
                refdata=self._refdata,
                respect=contracts,
                strict=strict,
            )
            for i in range(n)
        ]

    def run_records(
        self,
        *,
        n: int = 1,
        seed: int = 0,
        respect: RespectInput = None,
        strict: bool = False,
    ) -> "SimulationResult":
        """Run the compiled plan ``n`` times and return the batch as
        a :class:`SimulationResult` ready for ``.to_csv`` /
        ``.to_fasta`` / ``.to_dataframe`` export.

        Same arguments as :meth:`run`.
        """
        from .result import SimulationResult

        outcomes = self.run(n=n, seed=seed, respect=respect, strict=strict)
        return SimulationResult.from_outcomes(outcomes, self._refdata)

    def stream(
        self,
        *,
        n: Optional[int] = None,
        seed: int = 0,
        respect: RespectInput = None,
        strict: bool = False,
    ) -> Iterator["_engine.Outcome"]:
        """Lazily yield :class:`genairr_engine.Outcome` objects one at
        a time, without materialising the full batch in memory.

        Useful for large simulations where holding ``n`` outcomes
        would be wasteful — typical pattern is

        >>> for outcome in compiled.stream(n=1_000_000, seed=0):
        ...     write_to_disk(outcome)

        ``n=None`` (the default) yields outcomes indefinitely; the
        caller is expected to stop with ``itertools.islice``,
        ``break``, or similar. ``n=N`` yields exactly ``N`` outcomes
        with seeds ``seed`` … ``seed + N - 1``.

        ``respect`` and ``strict`` behave as in :meth:`run`.

        Raises ``ValueError`` when ``n`` is set to a value below 1.
        """
        if n is not None and n < 1:
            raise ValueError(f"n must be at least 1, got {n}")
        contracts = _coerce_respect(respect)
        i = 0
        while n is None or i < n:
            yield _engine.run(
                self._plan,
                seed=seed + i,
                refdata=self._refdata,
                respect=contracts,
                strict=strict,
            )
            i += 1

    def stream_records(
        self,
        *,
        n: Optional[int] = None,
        seed: int = 0,
        respect: RespectInput = None,
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
            self.stream(n=n, seed=seed, respect=respect, strict=strict)
        ):
            yield outcome_to_airr_record(
                outcome, self._refdata, sequence_id=f"{id_prefix}{i}"
            )

    def __repr__(self) -> str:
        return (
            f"<CompiledExperiment plan_len={len(self._plan)} "
            f"chain={self._refdata.chain_type}>"
        )
