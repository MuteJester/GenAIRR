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
from typing import Iterable, List, Optional, Sequence, Tuple, Union

import genairr_engine as _engine

from .dataconfig import DataConfig
from .dataconfig.enums import ChainType


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
    """One ``recombine()`` invocation, captured for later compilation."""

    np1_lengths: Tuple[Tuple[int, float], ...]
    np2_lengths: Tuple[Tuple[int, float], ...]

    def apply(
        self,
        plan: "_engine.PassPlan",
        refdata: "_engine.RefDataConfig",
    ) -> None:
        chain = refdata.chain_type
        np1 = list(self.np1_lengths)
        np2 = list(self.np2_lengths)

        if chain == "vj":
            plan.push_sample_allele("V", refdata)
            plan.push_sample_allele("J", refdata)
            plan.push_assemble("V")
            plan.push_generate_np("NP1", np1)
            plan.push_assemble("J")
        elif chain == "vdj":
            plan.push_sample_allele("V", refdata)
            plan.push_sample_allele("D", refdata)
            plan.push_sample_allele("J", refdata)
            plan.push_assemble("V")
            plan.push_generate_np("NP1", np1)
            plan.push_assemble("D")
            plan.push_generate_np("NP2", np2)
            plan.push_assemble("J")
        else:  # pragma: no cover — RefDataConfig only constructs vj/vdj.
            raise ValueError(f"unsupported chain_type {chain!r}")


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

# Anything :meth:`Experiment.on` accepts.
ExperimentInput = Union[str, DataConfig, "_engine.RefDataConfig"]


def _coerce_to_refdata(source: ExperimentInput) -> "_engine.RefDataConfig":
    """Normalise ``Experiment.on`` input to a Rust ``RefDataConfig``."""
    if isinstance(source, _engine.RefDataConfig):
        return source
    if isinstance(source, str):
        cfg = _resolve_config_name(source)
        return dataconfig_to_refdata(cfg)
    if isinstance(source, DataConfig):
        return dataconfig_to_refdata(source)
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

    __slots__ = ("_refdata", "_steps")

    def __init__(self, refdata: "_engine.RefDataConfig") -> None:
        self._refdata = refdata
        self._steps: List[_RecombineStep] = []

    @classmethod
    def on(cls, source: ExperimentInput) -> "Experiment":
        """Start an experiment against the given reference data.

        ``source`` is one of:
        - a config-name string (e.g. ``"human_igh"``),
        - a :class:`GenAIRR.DataConfig` instance,
        - a :class:`genairr_engine.RefDataConfig`.
        """
        return cls(_coerce_to_refdata(source))

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

    def recombine(
        self,
        *,
        np1_lengths: Optional[Iterable[Tuple[int, float]]] = None,
        np2_lengths: Optional[Iterable[Tuple[int, float]]] = None,
    ) -> "Experiment":
        """Append a standard V(D)J recombination step.

        Compiles to:
        - **VJ:** sample V → sample J → assemble V → generate NP1 →
          assemble J.
        - **VDJ:** sample V → sample D → sample J → assemble V →
          generate NP1 → assemble D → generate NP2 → assemble J.

        ``np1_lengths`` and ``np2_lengths`` default to the uniform
        ``[(0, 1.0), ..., (6, 1.0)]`` distribution. Pass an explicit
        iterable of ``(length, weight)`` tuples to override.
        ``np2_lengths`` is silently ignored on VJ chains.
        """
        step = _RecombineStep(
            np1_lengths=_normalize_lengths(np1_lengths),
            np2_lengths=_normalize_lengths(np2_lengths),
        )
        self._steps.append(step)
        return self

    def compile(self) -> "CompiledExperiment":
        """Compile the recorded steps into a ``PassPlan`` and return
        a :class:`CompiledExperiment` ready to be ``run()``.

        Idempotent: calling ``compile()`` twice produces two distinct
        ``CompiledExperiment`` instances with structurally-equal plans.
        """
        plan = _engine.PassPlan()
        for step in self._steps:
            step.apply(plan, self._refdata)
        return CompiledExperiment(plan, self._refdata)

    def run(self, *, n: int = 1, seed: int = 0) -> List["_engine.Outcome"]:
        """Compile and run this experiment ``n`` times.

        Equivalent to ``self.compile().run(n=n, seed=seed)``. Returns
        a list of :class:`genairr_engine.Outcome` objects.
        """
        return self.compile().run(n=n, seed=seed)

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

    def run(self, *, n: int = 1, seed: int = 0) -> List["_engine.Outcome"]:
        """Run the compiled plan ``n`` times.

        Each iteration uses ``seed + i`` as the per-run seed so
        consecutive batches stitch together by offsetting ``seed``.

        Raises ``ValueError`` for ``n < 1``.
        """
        if n < 1:
            raise ValueError(f"n must be at least 1, got {n}")
        return [
            _engine.run(self._plan, seed=seed + i, refdata=self._refdata)
            for i in range(n)
        ]

    def __repr__(self) -> str:
        return (
            f"<CompiledExperiment plan_len={len(self._plan)} "
            f"chain={self._refdata.chain_type}>"
        )
