"""Resolve ``Experiment.on(...)`` input into an engine-native
:class:`GenAIRR._engine.RefDataConfig`.

Three input shapes flow through this module:

- A **config-name string** (``"human_igh"``, ``"mouse_tcrb"``, …) —
  resolved through the builtin :mod:`GenAIRR.data` registry via
  :data:`_CONFIG_ALIASES`.
- A :class:`GenAIRR.DataConfig` — translated allele-by-allele into a
  :class:`GenAIRR._engine.RefDataConfig`.
- A bare :class:`GenAIRR._engine.RefDataConfig` — already engine-
  native, passed through unchanged.

Output is the typed pair ``(RefDataConfig, Optional[DataConfig])`` —
the DataConfig is retained when available because ``recombine()``
reads its empirical NP-length and per-gene trim distributions for
defaults; a bare engine-side RefDataConfig has no DataConfig backing
and falls through to the synthetic ``[0..6]`` placeholder.
"""
from __future__ import annotations

from typing import Optional, Tuple, Union

from GenAIRR import _engine  # private Rust extension submodule

from .dataconfig import DataConfig
from .dataconfig.enums import ChainType


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
