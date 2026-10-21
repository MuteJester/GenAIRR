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

A portable :class:`GenAIRR.Reference` is dispatched by
:meth:`Experiment.on` itself. A string or path naming a ``.json`` file
is rejected with a pointer to ``GenAIRR.Reference.read(path)``; the
file is not opened.

Output is the typed pair ``(RefDataConfig, Optional[DataConfig])`` —
the DataConfig is retained when available because ``recombine()``
reads its empirical NP-length and per-gene trim distributions for
defaults; a bare engine-side RefDataConfig has no DataConfig backing
and falls through to the synthetic ``[0..6]`` placeholder.
"""
from __future__ import annotations

import os
from typing import TYPE_CHECKING, List, Optional, Tuple, Union

from GenAIRR import _engine  # private Rust extension submodule

from .dataconfig import DataConfig
from .dataconfig.enums import ChainType

if TYPE_CHECKING:
    from ._reference.artifact import Reference


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


def _push_alleles(alleles_by_gene, add_method, *, derive_subregions: bool = False, refdata=None, segment=None) -> int:
    """Walk ``Dict[str, List[Allele]]`` and push every allele through
    ``add_method``. Returns the count of pushed alleles.

    Reads ``allele.functional_status`` (if present) and forwards it
    as the engine-native string. ``None`` and unrecognised values both
    fall through as ``None`` — the bridge never raises here because
    legacy bundled cartridges leave the field unset.

    When ``derive_subregions=True`` (passed by the V-allele branch
    only — D / J alleles have no IMGT V-subregion structure), the
    bridge looks for ``allele.subregions`` first (user-supplied,
    trust-but-validate) and falls back to deriving from an IMGT-gapped
    ``allele.gapped_seq`` via the existing ``imgt_regions`` helper.
    Alleles without IMGT gaps (or where derivation fails) are
    pushed with an empty subregion list — the cartridge stays
    valid; the manifest's coverage statistics surface the gap.
    """
    count = 0
    if not alleles_by_gene:
        return count
    for _gene, allele_list in alleles_by_gene.items():
        for allele in allele_list:
            seq_str = allele.ungapped_seq
            seq = seq_str.encode("ascii") if isinstance(seq_str, str) else bytes(seq_str)
            anchor = allele.anchor
            status = _normalise_functional_status(getattr(allele, "functional_status", None))
            kwargs = {}
            if anchor is not None:
                kwargs["anchor"] = int(anchor)
            if status is not None:
                kwargs["functional_status"] = status
            if derive_subregions:
                sub_entries = _resolve_v_subregions(allele)
                if sub_entries:
                    kwargs["subregions"] = sub_entries
            allele_id = add_method(allele.name, allele.gene, seq, **kwargs)
            evidence = getattr(allele, "reference_evidence", None)
            if evidence is not None and refdata is not None:
                import json
                refdata.set_allele_evidence(segment, allele_id, json.dumps(evidence, sort_keys=True))
            count += 1
    return count


def _resolve_v_subregions(allele) -> "list[tuple[str, int, int]]":
    """Resolve the V-region substructure annotations for a single
    allele as a list of ``(label, start, end)`` tuples ready for the
    Rust bridge.

    Priority:

    1. If the caller explicitly set ``allele.subregions`` (a dict
       ``{label: (start, end)}``), forward it verbatim — the bridge
       will validate ranges + overlap.
    2. Else, if ``gapped_seq`` declares IMGT gaps, derive via
       ``compute_v_region_boundaries`` from
       :mod:`GenAIRR.utilities.imgt_regions`. Ungapped sequences
       (e.g. bare FASTA) have no numbering alignment and stay unknown.
    3. Else, return ``[]`` — the cartridge will carry an empty
       subregion list for this allele.

    Failures in derivation (e.g. unexpected gapped_seq shape) are
    swallowed and surfaced as an empty list rather than aborting
    the cartridge load — a single allele's subregion gap should
    not break the whole bridge.
    """
    explicit = getattr(allele, "subregions", None)
    if explicit:
        # Already-validated user input: a dict ``{label: (start, end)}``.
        # Coerce to the bridge tuple shape; the Rust side validates
        # ranges, overlaps, and labels.
        try:
            return [
                (str(label), int(start), int(end))
                for label, (start, end) in explicit.items()
            ]
        except (TypeError, ValueError, AttributeError):
            return []

    from .utilities.imgt_regions import compute_v_region_boundaries, declares_imgt_gaps

    # Ungapped FASTA is also stored as ``gapped_seq``; only a declared
    # IMGT alignment can carry fixed IMGT subregion positions.
    if not declares_imgt_gaps(getattr(allele, "gapped_seq", None)):
        return []
    try:
        boundaries = compute_v_region_boundaries(allele)
    except Exception:
        return []
    return [
        (label, int(start), int(end))
        for label, (start, end) in boundaries.items()
        if start < end
    ]


# Map every value the Python Allele's ``functional_status`` could
# carry into the lowercase engine-native string. Accepts:
# - ``None`` (cartridge unannotated) → ``None``
# - a string (``"F"``, ``"ORF"``, ``"functional"``, …) — case
#   insensitive; unknown strings collapse to ``None`` so the bridge
#   stays silent on legacy data.
# - a native ``FunctionalStatus`` enum from
#   ``GenAIRR._native._anchor`` — read its ``.name`` attribute.
_KNOWN_STATUSES = {"functional", "orf", "pseudogene", "unknown"}
_STATUS_ALIASES = {
    "f": "functional",
    "p": "pseudogene",
}


def _normalise_functional_status(value):
    if value is None:
        return None
    # Native enum from GenAIRR._native._anchor → use its name.
    raw = getattr(value, "name", value) if not isinstance(value, str) else value
    if not isinstance(raw, str):
        return None
    lowered = raw.lower()
    canonical = _STATUS_ALIASES.get(lowered, lowered)
    if canonical in _KNOWN_STATUSES:
        return canonical
    return None


def dataconfig_to_refdata(cfg: DataConfig) -> "_engine.RefDataConfig":
    """Translate a :class:`DataConfig` species pickle into an
    engine-native :class:`GenAIRR._engine.RefDataConfig`.

    All V / D / J alleles are copied; the C segment is dropped (the
    engine has no C-segment passes). Anchorless alleles are preserved
    with ``anchor=None`` — they pass through but cannot satisfy the
    ``AnchorPreserved`` contract.

    Populates the cartridge's :class:`ReferenceIdentity` from
    ``cfg.metadata`` (species, locus, reference set) and ``cfg.name``.
    The cartridge's declared locus is then the canonical source for
    locus-derived defaults — including J-anchor rule fallback when
    no explicit :class:`ReferenceRulesSpec` is attached.
    """
    chain = _chain_type_label(cfg.metadata.chain_type if cfg.metadata else None)
    refdata = _engine.RefDataConfig(chain)

    _push_alleles(
        cfg.v_alleles, refdata.add_v_allele, derive_subregions=True, refdata=refdata, segment="V"
    )
    if chain == "vdj":
        _push_alleles(cfg.d_alleles, refdata.add_d_allele, refdata=refdata, segment="D")
    _push_alleles(cfg.j_alleles, refdata.add_j_allele, refdata=refdata, segment="J")

    # Stamp identity FIRST so the locus cascade below can consult it.
    _set_identity_from_dataconfig(refdata, cfg)

    rules = getattr(cfg, "reference_rules", None)
    if rules is not None:
        # User-authored ReferenceRulesSpec wins over locus inference.
        # Shape-check before crossing the PyO3 boundary so a bad spec
        # surfaces as a clean Python ValueError instead of being
        # mangled into a Rust panic.
        rules.validate()
        refdata.set_allowed_bases(list(rules.allowed_bases))
        refdata.set_v_anchor_rule(
            expected_aa=list(rules.v_anchor.expected_aa),
            required=rules.v_anchor.required,
            missing_severity=rules.v_anchor.missing_severity,
            mismatch_severity=rules.v_anchor.mismatch_severity,
        )
        refdata.set_j_anchor_rule(
            expected_aa=list(rules.j_anchor.expected_aa),
            required=rules.j_anchor.required,
            missing_severity=rules.j_anchor.missing_severity,
            mismatch_severity=rules.j_anchor.mismatch_severity,
        )
    else:
        # No user spec — derive J anchor expectations from the locus
        # cascade: identity.locus first (cartridge speaks for itself),
        # allele-name inference last (partial/synthetic catalogues).
        # Unknown locus leaves the lenient default in place.
        locus = refdata.identity().get("locus") or _infer_locus_prefix(cfg)
        residues = _j_anchor_residues_for_locus(locus)
        if residues is not None:
            refdata.set_j_anchor_rule(expected_aa=residues)

    return refdata


def _j_anchor_residues_for_locus(locus: Optional[str]) -> Optional[List[str]]:
    """Expected J anchor residues for ``locus``, or ``None`` when unknown."""
    if locus == "IGH":
        return ["W"]
    if locus in ("IGK", "IGL", "TRA", "TRB", "TRG", "TRD"):
        return ["F"]
    return None


# Map ``ChainType.value`` strings (from
# :class:`GenAIRR.dataconfig.enums.ChainType`) to AIRR locus prefixes.
# Locus is data, not behaviour, so this mapping lives at the Python
# boundary rather than being baked into Rust.
_CHAIN_TYPE_TO_LOCUS = {
    "BCR_HEAVY": "IGH",
    "BCR_LIGHT_KAPPA": "IGK",
    "BCR_LIGHT_LAMBDA": "IGL",
    "TCR_ALPHA": "TRA",
    "TCR_BETA": "TRB",
    "TCR_GAMMA": "TRG",
    "TCR_DELTA": "TRD",
}


def _set_identity_from_dataconfig(
    refdata: "_engine.RefDataConfig",
    cfg: DataConfig,
) -> None:
    """Bridge ``cfg.metadata`` + ``cfg.name`` into the Rust
    :class:`ReferenceIdentity`. Missing metadata leaves the
    corresponding slot ``None`` — identity is opt-in and degrades
    gracefully.
    """
    species: Optional[str] = None
    locus: Optional[str] = None
    reference_set: Optional[str] = None

    meta = cfg.metadata
    if meta is not None:
        species_attr = getattr(meta, "species", None)
        species_value = getattr(species_attr, "value", None)
        if isinstance(species_value, str) and species_value:
            species = species_value

        chain_attr = getattr(meta, "chain_type", None)
        chain_value = getattr(chain_attr, "value", None)
        if isinstance(chain_value, str):
            locus = _CHAIN_TYPE_TO_LOCUS.get(chain_value)

        rs = getattr(meta, "reference_set", None)
        if isinstance(rs, str) and rs:
            reference_set = rs

    name = cfg.name if isinstance(cfg.name, str) and cfg.name else None
    refdata.set_identity(
        species=species,
        locus=locus,
        reference_set=reference_set,
        name=name,
        source="DataConfig",
    )


def _infer_locus_prefix(cfg: DataConfig) -> Optional[str]:
    """Return the AIRR locus prefix for ``cfg``, derived from its
    first V allele's name. Returns ``None`` if the catalogue is empty
    or the first V name doesn't match a recognised prefix.

    ``DataConfig.v_alleles`` is a dict-of-lists keyed by gene; this
    helper peeks at the first allele's ``name`` attribute and
    extracts the 3-character locus prefix.
    """
    v_alleles = cfg.v_alleles
    if v_alleles is None:
        return None
    first_name: Optional[str] = None
    if isinstance(v_alleles, dict):
        for _gene, alleles in v_alleles.items():
            if alleles:
                first_name = getattr(alleles[0], "name", None)
                break
    elif isinstance(v_alleles, (list, tuple)) and v_alleles:
        first_name = getattr(v_alleles[0], "name", None)
    if not first_name or len(first_name) < 3:
        return None
    prefix = first_name[:3].upper()
    if prefix in {"IGH", "IGK", "IGL", "TRA", "TRB", "TRG", "TRD"}:
        return prefix
    return None


# Anything :meth:`Experiment.on` accepts. A ``Reference`` is dispatched by
# ``Experiment.on`` itself, because it also supplies the run's defaults.
ExperimentInput = Union[str, DataConfig, "_engine.RefDataConfig", "Reference"]


def _reference_path_error(source) -> TypeError:
    return TypeError(
        f"Experiment.on({str(source)!r}): Experiment.on does not read files. "
        f"Load a reference artifact with GenAIRR.Reference.read(path) and pass "
        f"the Reference to Experiment.on."
    )


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
        # A path to a reference artifact is never opened here; the
        # artifact must be loaded (and checked) by Reference.read.
        if (source.lower().replace("-", "_") not in _CONFIG_ALIASES
                and source.lower().endswith(".json")):
            raise _reference_path_error(source)
        cfg = _resolve_config_name(source)
        return dataconfig_to_refdata(cfg), cfg
    if isinstance(source, DataConfig):
        return dataconfig_to_refdata(source), source
    if isinstance(source, os.PathLike):
        raise _reference_path_error(os.fspath(source))
    raise TypeError(
        f"Experiment.on(...): expected a config-name string, "
        f"GenAIRR.DataConfig, GenAIRR._engine.RefDataConfig or "
        f"GenAIRR.Reference, got {type(source).__name__}"
    )
