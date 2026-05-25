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
from dataclasses import dataclass
from typing import Any, Dict, Iterable, Iterator, List, Optional, Sequence, Tuple, Union

from GenAIRR import _engine  # private Rust extension submodule

from .dataconfig import DataConfig
from .dataconfig.enums import ChainType


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
# DataConfig → GenAIRR._engine.RefDataConfig translator
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
    """Translate a :class:`DataConfig` species pickle into an
    engine-native :class:`GenAIRR._engine.RefDataConfig`.

    All V / D / J alleles are copied; the C segment is dropped (the
    engine has no C-segment passes). Anchorless alleles are preserved
    with ``anchor=None`` — they pass through but cannot satisfy the
    ``AnchorPreserved`` contract.
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
    weights_v: Optional[Tuple[float, ...]] = None
    weights_d: Optional[Tuple[float, ...]] = None
    weights_j: Optional[Tuple[float, ...]] = None
    # Set to True by ``Experiment.trim(...)``. Drives the
    # raw-RefDataConfig trim-noop warning at compile() time: when the
    # user explicitly called .trim(), silence the warning regardless
    # of whether trim is enabled or disabled.
    trim_overridden: bool = False

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
        v_weights = list(self.weights_v) if self.weights_v is not None else None
        d_weights = list(self.weights_d) if self.weights_d is not None else None
        j_weights = list(self.weights_j) if self.weights_j is not None else None

        if chain == "vj":
            plan.push_sample_allele("V", refdata, allowed_ids=v_ids, weights=v_weights)
            plan.push_sample_allele("J", refdata, allowed_ids=j_ids, weights=j_weights)
            if self.trim_v_3:
                plan.push_trim("V", "3", list(self.trim_v_3))
            if self.trim_j_5:
                plan.push_trim("J", "5", list(self.trim_j_5))
            plan.push_assemble("V")
            plan.push_generate_np("NP1", np1)
            plan.push_assemble("J")
        elif chain == "vdj":
            plan.push_sample_allele("V", refdata, allowed_ids=v_ids, weights=v_weights)
            plan.push_sample_allele("D", refdata, allowed_ids=d_ids, weights=d_weights)
            plan.push_sample_allele("J", refdata, allowed_ids=j_ids, weights=j_weights)
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
@dataclass(frozen=True)
class _ClonalForkStep:
    """Marks the boundary between "per-clone" passes (run once per
    clonal family — typically `recombine`) and "per-descendant" passes
    (run once per read inside the family — typically `mutate`,
    `corrupt_*`).

    The compile pipeline splits the experiment's step list at this
    marker, builds two `CompiledExperiment`s, and the runtime
    forks the parent IR into descendants for each clone.
    """

    n_clones: int
    size: int


_CORRUPT_KIND_PCR = "pcr"
_CORRUPT_KIND_QUALITY = "quality"
_CORRUPT_KIND_INDEL = "indel"
_CORRUPT_KIND_CONTAMINANT = "contaminant"
_CORRUPT_KIND_REV_COMP = "rev_comp"
_CORRUPT_KIND_5PRIME_LOSS = "5prime_loss"
_CORRUPT_KIND_3PRIME_LOSS = "3prime_loss"
_CORRUPT_KIND_NS = "ns"


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
        elif self.kind == _CORRUPT_KIND_REV_COMP:
            plan.push_corrupt_rev_comp(self.apply_prob)
        elif self.kind == _CORRUPT_KIND_5PRIME_LOSS:
            plan.push_corrupt_5prime_loss(count_list)
        elif self.kind == _CORRUPT_KIND_3PRIME_LOSS:
            plan.push_corrupt_3prime_loss(count_list)
        elif self.kind == _CORRUPT_KIND_NS:
            plan.push_corrupt_ns(count_list)
        else:  # pragma: no cover — guarded at builder time.
            raise ValueError(f"unsupported corruption kind {self.kind!r}")


# ──────────────────────────────────────────────────────────────────
# describe() helpers — render step configurations as biology-style
# narrative strings. The output is the diagnostic instrument we use
# to grade DSL readability: if a chain reads weird in describe(), the
# chain itself reads weird.
# ──────────────────────────────────────────────────────────────────


def _describe_distribution(pairs: Tuple[Tuple[int, float], ...]) -> str:
    """Render a ``(value, weight)`` distribution as a short
    biology-friendly string: ``"5"`` for a fixed count, ``"5–15"``
    for a uniform range, ``"0–8 weighted"`` for non-uniform."""
    if not pairs:
        return "<empty>"
    values = [v for v, _ in pairs]
    weights = [w for _, w in pairs]
    lo, hi = min(values), max(values)
    if lo == hi:
        return str(lo)
    uniform = all(abs(w - weights[0]) < 1e-9 for w in weights)
    covers_range = sorted(values) == list(range(lo, hi + 1))
    if uniform and covers_range:
        return f"{lo}–{hi}"
    return f"{lo}–{hi} weighted"


def _describe_recombine_step(step: "_RecombineStep", chain_type: str) -> str:
    segs = "V/J" if chain_type == "vj" else "V/D/J"
    parts: List[str] = [f"sample {segs} alleles"]

    locked_segs = []
    if step.locks_v is not None:
        locked_segs.append(f"V ({len(step.locks_v)})")
    if step.locks_d is not None:
        locked_segs.append(f"D ({len(step.locks_d)})")
    if step.locks_j is not None:
        locked_segs.append(f"J ({len(step.locks_j)})")
    if locked_segs:
        parts.append(f"locked to {', '.join(locked_segs)}")

    weighted_segs = []
    for seg, weights in (("V", step.weights_v), ("D", step.weights_d), ("J", step.weights_j)):
        if weights is not None:
            weighted_segs.append(seg)
    if weighted_segs:
        parts.append(f"custom {'/'.join(weighted_segs)} allele weights")

    trim_bits = []
    if step.trim_v_3:
        trim_bits.append("V3'")
    if step.trim_d_5:
        trim_bits.append("D5'")
    if step.trim_d_3:
        trim_bits.append("D3'")
    if step.trim_j_5:
        trim_bits.append("J5'")
    if trim_bits:
        parts.append(f"empirical exonuclease trim ({', '.join(trim_bits)})")

    np1_desc = _describe_distribution(step.np1_lengths)
    if chain_type == "vdj":
        np2_desc = _describe_distribution(step.np2_lengths)
        parts.append(f"insert NP1 ({np1_desc} bases) and NP2 ({np2_desc} bases)")
    else:
        parts.append(f"insert NP1 ({np1_desc} bases)")

    return "V(D)J recombination: " + "; ".join(parts)


_S5F_KERNEL_LABELS = {
    "hh_s5f": "human heavy-chain (HH_S5F)",
    "hkl_s5f": "human kappa/lambda (HKL_S5F)",
    "mk_rs5nf": "mouse (MK_RS5NF)",
}


def _describe_mutate_step(step: "_MutateStep") -> str:
    count = _describe_distribution(step.count_pairs)
    if step.model == "s5f":
        kernel_label = _S5F_KERNEL_LABELS.get(step.s5f_model_name, step.s5f_model_name)
        return f"Somatic hypermutation (S5F context model, {kernel_label}): {count} mutations/record"
    if step.model == "uniform":
        return f"Somatic hypermutation (uniform, position-independent): {count} mutations/record"
    return f"Mutation ({step.model}): {count}/record"


_CORRUPT_NARRATIVE_LABELS = {
    _CORRUPT_KIND_PCR: ("PCR substitution errors", "errors/record"),
    _CORRUPT_KIND_QUALITY: ("Sequencing quality errors", "errors/record"),
    _CORRUPT_KIND_5PRIME_LOSS: ("5'-end loss (primer/adapter trim)", "bases trimmed"),
    _CORRUPT_KIND_3PRIME_LOSS: ("3'-end loss", "bases trimmed"),
    _CORRUPT_KIND_NS: ("Low-quality N-base injection", "bases replaced with N"),
}


def _describe_corrupt_step(step: "_CorruptStep") -> str:
    if step.kind in _CORRUPT_NARRATIVE_LABELS:
        label, unit = _CORRUPT_NARRATIVE_LABELS[step.kind]
        count = _describe_distribution(step.count_pairs)
        return f"{label}: {count} {unit}"
    if step.kind == _CORRUPT_KIND_INDEL:
        count = _describe_distribution(step.count_pairs)
        return (
            f"Structural indels (library prep): {count} events/record, "
            f"insertion fraction {step.insertion_prob:.2f}"
        )
    if step.kind == _CORRUPT_KIND_CONTAMINANT:
        return (
            f"Cross-sample contamination: {step.apply_prob:.0%} of records "
            f"replaced with a random allele sequence"
        )
    if step.kind == _CORRUPT_KIND_REV_COMP:
        return f"Random strand orientation: {step.apply_prob:.0%} of records reverse-complemented"
    return f"Corruption ({step.kind})"  # pragma: no cover — guarded at builder time


def _describe_clonal_fork_step(step: "_ClonalForkStep") -> str:
    return f"Clonal expansion: {step.n_clones} lineages × {step.size} reads/clone"


def _describe_step(step: Any, chain_type: str) -> str:
    """Dispatch to the right step-describer based on the step type."""
    if isinstance(step, _RecombineStep):
        return _describe_recombine_step(step, chain_type)
    if isinstance(step, _MutateStep):
        return _describe_mutate_step(step)
    if isinstance(step, _CorruptStep):
        return _describe_corrupt_step(step)
    if isinstance(step, _ClonalForkStep):
        return _describe_clonal_fork_step(step)
    return f"Unknown step: {type(step).__name__}"  # pragma: no cover


_CHAIN_TYPE_BIOLOGY_LABELS = {
    "BCR_HEAVY": "heavy-chain BCR",
    "BCR_LIGHT_KAPPA": "kappa light-chain BCR",
    "BCR_LIGHT_LAMBDA": "lambda light-chain BCR",
    "TCR_ALPHA": "alpha-chain TCR",
    "TCR_BETA": "beta-chain TCR",
    "TCR_GAMMA": "gamma-chain TCR",
    "TCR_DELTA": "delta-chain TCR",
}


def _describe_experiment_header(
    refdata: "_engine.RefDataConfig",
    dataconfig: Optional[DataConfig],
) -> str:
    """One-line header naming the refdata source and chain type."""
    chain = refdata.chain_type
    if dataconfig is not None and getattr(dataconfig, "metadata", None) is not None:
        meta = dataconfig.metadata
        species_label: Optional[str] = None
        species = getattr(meta, "species", None)
        if species is not None:
            species_label = getattr(species, "value", None) or str(species)
        chain_enum = getattr(meta, "chain_type", None)
        chain_key = getattr(chain_enum, "value", None) if chain_enum is not None else None
        chain_label = _CHAIN_TYPE_BIOLOGY_LABELS.get(chain_key or "", f"{chain.upper()} chain")
        if species_label:
            return f"Experiment on {species_label} {chain_label}"
        return f"Experiment on {chain_label}"
    return f"Experiment on raw RefDataConfig ({chain.upper()})"


def _describe_step_sequence(
    steps: Sequence[Any],
    chain_type: str,
    *,
    indent: str = "  ",
    start_index: int = 1,
) -> List[str]:
    """Render a sequence of steps as numbered narrative lines. The
    clonal-fork step gets a visual divider so the pre/post boundary
    is obvious; subsequent steps continue numbering.

    Returns one string per output line (no trailing newlines).
    """
    out: List[str] = []
    step_no = start_index
    for step in steps:
        if isinstance(step, _ClonalForkStep):
            out.append(f"{indent}── {_describe_clonal_fork_step(step)} ──")
            out.append(f"{indent}    (steps above run once per clone; "
                       f"steps below run once per descendant)")
            continue
        out.append(f"{indent}{step_no}. {_describe_step(step, chain_type)}")
        step_no += 1
    return out


_PRODUCTIVE_BUNDLE_NAMES = frozenset({
    "productive_junction_frame",
    "no_stop_codon_in_junction",
    "anchor_preserved.v",
    "anchor_preserved.j",
})


def _format_declared_contracts(declared: Sequence[str]) -> Optional[str]:
    """Render the contract bundles declared by `.productive_only()`
    etc. on an uncompiled Experiment as a biology-style line.

    Companion to :func:`_format_active_contracts`, which renders the
    compiled simulator's per-pass contract names (the underlying
    primitives the bundles expand into). On the Experiment side we
    show the bundle name directly because that's what the user typed.
    """
    if not declared:
        return None
    labels = []
    for bundle in declared:
        if bundle == "productive":
            labels.append("productive sequences only")
        else:
            labels.append(bundle)
    return "; ".join(labels)


def _format_active_contracts(active: Sequence[str]) -> Optional[str]:
    """Render the compiled simulator's active-contracts tuple as a
    biology-style line, or ``None`` if no contracts are bound. The
    `ga.productive()` bundle's four underlying contracts
    (productive_junction_frame, no_stop_codon_in_junction, and the
    two anchor_preserved contracts) collapse to the single label
    "productive sequences only" so the line reads naturally."""
    if not active:
        return None
    names = list(active)
    extras: List[str] = []
    productive_present = _PRODUCTIVE_BUNDLE_NAMES.issubset(set(names))
    bundle_names = _PRODUCTIVE_BUNDLE_NAMES if productive_present else frozenset()
    for n in names:
        if n in bundle_names:
            continue
        extras.append(n)
    labels: List[str] = []
    if productive_present:
        labels.append("productive sequences only")
    labels.extend(extras)
    if not labels:
        return None
    return "; ".join(labels)


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
        f"GenAIRR.DataConfig, or GenAIRR._engine.RefDataConfig, "
        f"got {type(source).__name__}"
    )


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
        position as corrupted (the sequencing-error convention).
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

    def corrupt_5prime_loss(
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

    def corrupt_3prime_loss(
        self,
        *,
        length: Union[int, Tuple[int, int], Iterable[Tuple[int, float]]],
    ) -> "Experiment":
        """Append a 3'-end loss corruption step.

        Drops bases from the end of the assembled sequence to model
        read-end degradation. Same ``length`` shapes as
        :meth:`corrupt_5prime_loss`.
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

    def corrupt_ns(
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

    def corrupt_reverse_complement(self, *, prob: float = 0.5) -> "Experiment":
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
        is enforced during recombination and SHM passes, so
        non-admissible candidates are filtered before commit (the
        runtime falls back to permissive sampling if no admissible
        candidate exists; use ``strict=True`` at run time to raise
        instead).

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

    def with_clonal_structure(
        self,
        *,
        n_clones: int,
        size: int,
    ) -> "Experiment":
        """Split the pipeline into per-clone (parent) and
        per-descendant phases for clonal-family generation.

        Steps appended *before* this call run **once per clone** —
        typically just :meth:`recombine`, which establishes the
        parent V/D/J + trim + NP + assembled IR for the clonal
        family. Steps appended *after* this call run **once per
        read inside the family** — typically :meth:`mutate` and the
        ``corrupt_*`` ops, which introduce per-read divergence
        within the clone.

        Concrete shape::

            exp = (Experiment.on("human_igh")
                   .recombine()
                   .with_clonal_structure(n_clones=10, size=20)
                   .mutate(count=8)
                   .corrupt_pcr(count=2))
            result = exp.run_records(seed=0)
            # 10 clones × 20 descendants = 200 records.
            # Each record carries a ``clone_id`` integer in [0, 10).

        ``n`` can be omitted from :meth:`run_records` for a clonal
        experiment — the runtime expands ``n_clones * size``
        records automatically. Passing ``n`` is allowed only when
        ``n == n_clones * size``.

        Constraints:
        - Both ``n_clones`` and ``size`` must be positive ints.
        - At most one fork per pipeline; calling this method twice
          raises ``ValueError``.

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
        if not isinstance(size, int) or isinstance(size, bool) or size < 1:
            raise ValueError(f"size must be a positive int, got {size!r}")
        if any(isinstance(s, _ClonalForkStep) for s in self._steps):
            raise ValueError(
                "with_clonal_structure() can only be called once per pipeline"
            )
        self._steps.append(_ClonalForkStep(n_clones=n_clones, size=size))
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

        **TCR guard:** somatic hypermutation is a B-cell
        phenomenon — T-cells do not undergo SHM in the periphery.
        Calling ``.mutate()`` on a TCR-configured experiment raises
        ``ValueError`` to prevent silent biological misuse. Use
        ``corrupt_pcr`` / ``corrupt_quality`` for sequencing-error
        realism on TCR data instead.
        """
        if self._is_tcr_refdata():
            raise ValueError(
                "mutate(): somatic hypermutation does not occur in TCR "
                "sequences (T-cells lack AID and the SHM machinery). The "
                "configured refdata is a TCR locus. For sequencing-error "
                "realism on TCR data, use corrupt_pcr / corrupt_quality / "
                "corrupt_indels instead."
            )
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
        when the pipeline contains a :meth:`with_clonal_structure`
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
            # post-fork plan inherits the parent's V/D/J/NP backbone, so
            # the recombination-time contracts have already been enforced
            # upstream. Re-passing them here would force the compiler to
            # look for support (np.np1.length, anchor trims) in passes
            # the post-fork plan doesn't contain — which would fail any
            # pipeline that combines `with_clonal_structure()` with
            # `.productive_only()`.
            post_simulator = self._build_simulator(
                post_steps, None, any_lock=False, replace_fn=_replace
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
            step.apply(plan, self._refdata)
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
        ``compile()`` with ``ValueError``. For runtime residue,
        ``strict=False`` (default) lets a pass fall back to
        unconstrained sampling when no admissible candidate exists;
        ``strict=True`` raises
        :class:`GenAIRR._engine.StrictSamplingError` instead.
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
        dataconfig: Optional[DataConfig] = None,
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
        """Run the compiled simulator ``n`` times.

        Each iteration uses ``seed + i`` as the per-run seed so
        consecutive batches stitch together by offsetting ``seed``.

        ``strict`` controls the failure mode when filtering yields
        no admissible candidate:
        - ``False`` (default, **permissive**) — fall back to
          unconstrained sampling and continue.
        - ``True`` (**strict**) — raise
          :class:`GenAIRR._engine.StrictSamplingError` carrying the
          failing pass name, trace address, and failure reason.

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
    ) -> "SimulationResult":
        """Run the compiled simulator ``n`` times and return the batch as
        a :class:`SimulationResult` ready for ``.to_csv`` /
        ``.to_fasta`` / ``.to_dataframe`` export.

        Same arguments as :meth:`run`. ``expose_provenance=True``
        appends `truth_v_call/d_call/j_call` columns reflecting the
        originally-sampled allele names.
        """
        from .result import SimulationResult

        outcomes = self.run(n=n, seed=seed, strict=strict)
        return SimulationResult.from_outcomes(
            outcomes, self._refdata, expose_provenance=expose_provenance
        )

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
        dataconfig: Optional[DataConfig] = None,
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
    ) -> "SimulationResult":
        """Same as :meth:`run` but returns a :class:`SimulationResult`
        with each record dict carrying an integer ``clone_id`` field
        in ``[0, n_clones)``. ``expose_provenance=True`` also
        appends `truth_v_call` / `truth_d_call` / `truth_j_call`
        columns from the originally-sampled allele names.
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
                rec = outcome_to_airr_record(
                    desc,
                    self._refdata,
                    sequence_id=f"clone{clone_idx}_desc{desc_idx}",
                )
                rec["clone_id"] = clone_idx
                if expose_provenance:
                    _inject_truth_columns(desc, self._refdata, rec)
                records.append(rec)
        return SimulationResult(records, outcomes=outcomes)

    def __repr__(self) -> str:
        return (
            f"<CompiledClonalExperiment n_clones={self._fork.n_clones} "
            f"size={self._fork.size} chain={self._refdata.chain_type}>"
        )
