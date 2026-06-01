"""Extract empirical distributions from a :class:`DataConfig` into the
``[(value, weight), ...]`` shape that ``GenAIRR._engine.PassPlan`` expects.

A :class:`DataConfig` carries per-species reference distributions.
``recombine()`` reads these through the helpers in this module so
simulations default to the species' empirical NP lengths and trim
amounts instead of the uniform ``[0..6]`` placeholder.

Two layers of granularity:

- **NP lengths** are species-wide ``{length: prob}`` dicts in
  ``cfg.NP_lengths`` keyed by ``"NP1"`` / ``"NP2"``. These map
  directly to the engine's empirical-length-distribution shape.

- **Trim amounts** are nested ``{family: {gene: {trim: prob}}}``
  in ``cfg.trim_dicts`` keyed by ``"V_3" / "D_5" / "D_3" / "J_5"``.
  We **marginalize** across families and genes to a single
  segment-level distribution by averaging the per-gene
  ``{trim: prob}`` dicts. That loses per-gene precision but gives
  us biologically reasonable defaults without needing a per-gene
  dispatch in the engine.
"""
from __future__ import annotations

from collections import defaultdict
from typing import Any, Dict, List, Optional, Tuple


# Per-segment trim-dict keys in `DataConfig.trim_dicts`. Maps the
# Rust trace-address fragment ``(segment, end)`` to the DataConfig
# key. Keep the mapping explicit so readers don't have to infer the
# naming convention.
_TRIM_DICT_KEYS: Dict[Tuple[str, str], str] = {
    ("V", "3"): "V_3",
    ("D", "5"): "D_5",
    ("D", "3"): "D_3",
    ("J", "5"): "J_5",
}


def extract_np_lengths(
    cfg: Any, np_segment: str
) -> Optional[List[Tuple[int, float]]]:
    """Return the empirical length distribution for ``np_segment``
    (``"NP1"`` or ``"NP2"``) as a ``[(length, weight), ...]`` list,
    or ``None`` if the config doesn't carry one.

    The DataConfig stores this as a ``{length: prob}`` dict; the
    engine accepts any positive-weight pairs (it normalises
    internally), so we pass the probabilities through as weights.
    """
    if not getattr(cfg, "NP_lengths", None):
        return None
    inner = cfg.NP_lengths.get(np_segment)
    if not inner:
        return None
    pairs: List[Tuple[int, float]] = []
    for length, prob in inner.items():
        try:
            length_int = int(length)
        except (TypeError, ValueError):
            continue
        weight = float(prob)
        if weight <= 0.0:
            continue
        pairs.append((length_int, weight))
    if not pairs:
        return None
    pairs.sort(key=lambda p: p[0])
    return pairs


def _min_allele_length(alleles_by_gene: Optional[Dict[str, list]]) -> Optional[int]:
    """Smallest allele length across every gene → allele in
    ``alleles_by_gene``. ``None`` when the dict is missing or empty.
    """
    if not alleles_by_gene:
        return None
    smallest: Optional[int] = None
    for _gene, allele_list in alleles_by_gene.items():
        if not allele_list:
            continue
        for a in allele_list:
            length = getattr(a, "ungapped_len", None) or getattr(a, "length", None)
            if length is None:
                continue
            length_int = int(length)
            if smallest is None or length_int < smallest:
                smallest = length_int
    return smallest


def _trim_cap(min_len: Optional[int], *, two_sided: bool) -> Optional[int]:
    """Cap a trim distribution so the assembled segment stays
    long enough to be meaningful even on the shortest allele.

    For one-sided trims (V_3 alone, J_5 alone) the cap leaves a 3 bp
    margin so the conserved Cys / W/F anchor codon survives the
    trim. For two-sided trims (both D_5 and D_3 on the same allele)
    the cap is halved so each side cannot independently consume
    more than half the available room.
    """
    if min_len is None:
        return None
    margin = 3  # protect the anchor codon
    cap = max(0, min_len - margin)
    if two_sided:
        cap //= 2
    return cap


def extract_trim_distribution(
    cfg: Any,
    segment: str,
    end: str,
    *,
    cap: Optional[int] = None,
) -> Optional[List[Tuple[int, float]]]:
    """Return a marginalized trim distribution for ``(segment, end)``
    (e.g. ``("V", "3")``) as a ``[(trim, weight), ...]`` list, or
    ``None`` if not present.

    The DataConfig stores trims as ``{family: {gene: {trim: prob}}}``
    keyed by the segment+end string (``"V_3"`` / ``"D_5"`` / ``"D_3"``
    / ``"J_5"``). We marginalize by averaging the per-gene
    distributions: for each trim amount, weight =
    sum-over-genes(prob[trim]) / number-of-genes.

    ``cap`` (optional) clamps the support so trims past that value
    are filtered out. Used to keep the marginalized distribution
    safe against short alleles in the pool — see :func:`_trim_cap`.
    """
    key = _TRIM_DICT_KEYS.get((segment, end))
    if key is None:
        return None
    if not getattr(cfg, "trim_dicts", None):
        return None
    segment_dict = cfg.trim_dicts.get(key)
    if not segment_dict:
        return None

    total_per_trim: Dict[int, float] = defaultdict(float)
    n_genes = 0
    for _family, genes in segment_dict.items():
        if not genes:
            continue
        for _gene, trim_dist in genes.items():
            if not trim_dist:
                continue
            n_genes += 1
            for trim, prob in trim_dist.items():
                try:
                    trim_int = int(trim)
                except (TypeError, ValueError):
                    continue
                weight = float(prob)
                if weight < 0.0:
                    continue
                total_per_trim[trim_int] += weight
    if n_genes == 0 or not total_per_trim:
        return None
    pairs: List[Tuple[int, float]] = [
        (trim, total / n_genes)
        for trim, total in total_per_trim.items()
        if total > 0.0
    ]
    if cap is not None:
        pairs = [(t, w) for t, w in pairs if t <= cap]
    if not pairs:
        return None
    pairs.sort(key=lambda p: p[0])
    return pairs


def extract_recombine_defaults(cfg: Any) -> Dict[str, Optional[List[Tuple[int, float]]]]:
    """Extract every empirical distribution :meth:`Experiment.recombine`
    might use, in one pass.

    Returns a dict with keys ``np1`` / ``np2`` / ``trim_v_3`` /
    ``trim_d_5`` / ``trim_d_3`` / ``trim_j_5``, each mapping to a
    ``[(value, weight), ...]`` list or ``None``.

    Resolution order per key:

    1. **Explicit ``ReferenceEmpiricalModels``** attached to the
       config (``cfg.reference_models``). Cartridge-authored defaults
       win — they're typed, validated, and live alongside identity /
       rules.
    2. **Legacy nested-dict extraction** from ``cfg.NP_lengths`` /
       ``cfg.trim_dicts`` (the historical path). Trim distributions
       are clamped per-segment based on the pool's shortest allele
       (see ``_trim_cap``).
    3. ``None`` — the caller falls back to the uniform placeholder
       (raw-RefDataConfig path).
    """
    explicit = _explicit_models(cfg)
    if explicit:
        # Validate before consuming. Catches malformed user-authored
        # specs before they propagate into the engine, with a
        # cartridge-tagged error.
        chain_type = _chain_type_label(cfg)
        try:
            explicit.validate(chain_type=chain_type)
        except ValueError as exc:
            raise ValueError(
                f"DataConfig.reference_models failed shape validation: {exc}"
            ) from exc

    cap_v = _trim_cap(_min_allele_length(getattr(cfg, "v_alleles", None)), two_sided=False)
    cap_d = _trim_cap(_min_allele_length(getattr(cfg, "d_alleles", None)), two_sided=True)
    cap_j = _trim_cap(_min_allele_length(getattr(cfg, "j_alleles", None)), two_sided=False)

    return {
        "np1": _np_lengths_from_models(explicit, "NP1") or extract_np_lengths(cfg, "NP1"),
        "np2": _np_lengths_from_models(explicit, "NP2") or extract_np_lengths(cfg, "NP2"),
        # NP base sampling distributions. Slice — Typed NP base
        # model. Cascade is `typed → uniform` (the legacy auto-lift
        # of `NP_first_bases` / `NP_transitions` is deliberately
        # deferred per `docs/junction_n_addition_audit.md` because
        # auto-lifting would silently change output bytes vs the
        # pre-slice baseline). `None` means UniformBase applies.
        "np1_bases": _np_bases_from_models(explicit, "NP1"),
        "np2_bases": _np_bases_from_models(explicit, "NP2"),
        # Markov transition matrices, lowered separately so the
        # bridge can fire `push_generate_np(base_pairs=...,
        # markov_transitions=...)` in lockstep. `None` for
        # uniform / empirical_first_base specs — those flow
        # through `base_pairs` alone, byte-identical to the
        # pre-Markov-slice plan signature.
        "np1_markov_transitions": _np_markov_transitions_from_models(
            explicit, "NP1"
        ),
        "np2_markov_transitions": _np_markov_transitions_from_models(
            explicit, "NP2"
        ),
        # P-nucleotide per-end length distributions (Slice —
        # P-nucleotide v1). `None` means the pipeline omits
        # the P-pass for that end — byte-identical to the
        # pre-slice baseline. Legacy `p_nucleotide_length_probs`
        # is NOT auto-lifted; the orphan boundary holds.
        "p_v_3_lengths": _p_nucleotide_lengths_from_models(
            explicit, "V_3"
        ),
        "p_d_5_lengths": _p_nucleotide_lengths_from_models(
            explicit, "D_5"
        ),
        "p_d_3_lengths": _p_nucleotide_lengths_from_models(
            explicit, "D_3"
        ),
        "p_j_5_lengths": _p_nucleotide_lengths_from_models(
            explicit, "J_5"
        ),
        "trim_v_3": _trim_from_models(explicit, "V_3", cap=cap_v)
            or extract_trim_distribution(cfg, "V", "3", cap=cap_v),
        "trim_d_5": _trim_from_models(explicit, "D_5", cap=cap_d)
            or extract_trim_distribution(cfg, "D", "5", cap=cap_d),
        "trim_d_3": _trim_from_models(explicit, "D_3", cap=cap_d)
            or extract_trim_distribution(cfg, "D", "3", cap=cap_d),
        "trim_j_5": _trim_from_models(explicit, "J_5", cap=cap_j)
            or extract_trim_distribution(cfg, "J", "5", cap=cap_j),
        # Per-segment allele usage weights (Slice — Allele Usage
        # Estimation v1). Each is a ``{allele_name: weight}`` dict
        # OR ``None``. ``None`` means "no cartridge default for
        # this segment — uniform applies unless overridden by the
        # ``Experiment.recombine(*_allele_weights=...)`` kwarg".
        # The empty-segment-dict case (`{}`) is collapsed to
        # ``None`` here so the downstream resolver only sees one
        # "uniform" code path. Legacy ``DataConfig.gene_use_dict``
        # is NOT auto-lifted; the orphan boundary holds.
        "allele_usage_v": _allele_usage_from_models(explicit, "V"),
        "allele_usage_d": _allele_usage_from_models(explicit, "D"),
        "allele_usage_j": _allele_usage_from_models(explicit, "J"),
    }


def _explicit_models(cfg: Any):
    """Return ``cfg.reference_models`` when it's a non-empty
    :class:`ReferenceEmpiricalModels`, else ``None``. The
    ``__getattr__`` shim on legacy DataConfigs returns ``None``;
    new instances default to ``None`` too. An empty container is
    treated as "no explicit model" so legacy fallback still runs.
    """
    rm = getattr(cfg, "reference_models", None)
    if rm is None:
        return None
    # Empty container → no signal, defer to legacy path.
    if (
        not getattr(rm, "np_lengths", None)
        and not getattr(rm, "trims", None)
        and not getattr(rm, "np_bases", None)
        and not getattr(rm, "p_nucleotide_lengths", None)
        and not getattr(rm, "allele_usage", None)
    ):
        return None
    return rm


def _chain_type_label(cfg: Any) -> Optional[str]:
    """Return ``"vj"`` or ``"vdj"`` from ``cfg.metadata.chain_type``,
    or ``None`` when metadata is missing.
    """
    meta = getattr(cfg, "metadata", None)
    if meta is None:
        return None
    has_d_attr = getattr(meta, "has_d", None)
    if isinstance(has_d_attr, bool):
        return "vdj" if has_d_attr else "vj"
    return None


def _np_lengths_from_models(
    models: Any, key: str,
) -> Optional[List[Tuple[int, float]]]:
    """Read ``np_lengths[key]`` off an explicit
    :class:`ReferenceEmpiricalModels` (or return ``None`` if absent).

    The spec is already validated (see ``extract_recombine_defaults``)
    — we just lower it into the engine's ``[(value, weight), ...]``
    shape with a stable sort order.
    """
    if models is None:
        return None
    spec = models.np_lengths.get(key)
    if spec is None:
        return None
    pairs = [(int(v), float(w)) for v, w in spec.values]
    pairs.sort(key=lambda p: p[0])
    return pairs or None


# Canonical A/C/G/T → byte mapping used to lower NP base
# categorical specs into the engine's `(u8, f64)` wire shape.
# Sorted alphabetically so two specs with the same weights but
# different key insertion orders lower to the same pairs (and
# therefore the same plan signature).
_NP_BASE_BYTE: Dict[str, int] = {"A": 65, "C": 67, "G": 71, "T": 84}


def _np_bases_from_models(
    models: Any, key: str,
) -> Optional[List[Tuple[int, float]]]:
    """Read ``np_bases[key]`` off an explicit
    :class:`ReferenceEmpiricalModels` and lower the first-base
    categorical into the engine's ``[(base_byte, weight), ...]``
    shape. Returns ``None`` when the spec is absent or the spec
    is the ``"uniform"`` kind (in both cases the bridge default
    fires and ``UniformNpGenerator`` applies).

    For ``"markov"`` specs this returns the first-base row only;
    the transition matrix is lowered separately via
    :func:`_np_markov_transitions_from_models` so the two
    bridge kwargs (``base_pairs`` + ``markov_transitions``) can
    be filled in lockstep at compile time.
    """
    if models is None:
        return None
    spec = getattr(models, "np_bases", {}).get(key)
    if spec is None:
        return None
    if spec.kind == "uniform":
        return None  # Bridge default — no `base_pairs` kwarg sent.
    # ``empirical_first_base`` AND ``markov`` both carry a
    # ``first_base`` categorical; for Markov this is the
    # position-0 row. The bridge's Markov path enforces full
    # A/C/G/T alphabet coverage, so emit every canonical base
    # (zeros included) for Markov; for empirical_first_base
    # keep the legacy filter that drops zeros so existing plan
    # signatures stay byte-identical.
    include_zeros = spec.kind == "markov"
    out: List[Tuple[int, float]] = []
    for letter in ("A", "C", "G", "T"):
        weight = spec.first_base.get(letter, 0.0)
        if include_zeros or weight > 0.0:
            out.append((_NP_BASE_BYTE[letter], float(weight)))
    return out or None


def _np_markov_transitions_from_models(
    models: Any, key: str,
) -> Optional[List[List[Tuple[int, float]]]]:
    """Read ``np_bases[key]`` off an explicit
    :class:`ReferenceEmpiricalModels` and lower the
    ``transitions`` matrix into the engine's
    ``[ [(to_byte, weight), ...], ... ]`` shape with rows in
    canonical A/C/G/T from-base order. Returns ``None`` when
    the spec is absent or its kind is not ``"markov"``.

    Slice — Markov NP Base Generator. The Python spec layer
    (``NpBaseModelSpec.validate``) already enforces row
    completeness over the full A/C/G/T from-base alphabet;
    here we just deterministically order the rows and the
    per-row pairs so plan signatures are byte-stable across
    runs.
    """
    if models is None:
        return None
    spec = getattr(models, "np_bases", {}).get(key)
    if spec is None:
        return None
    if spec.kind != "markov":
        return None
    rows: List[List[Tuple[int, float]]] = []
    for from_letter in ("A", "C", "G", "T"):
        row_dict = spec.transitions[from_letter]
        row_pairs: List[Tuple[int, float]] = []
        for to_letter in ("A", "C", "G", "T"):
            weight = row_dict.get(to_letter, 0.0)
            # Always emit the full A/C/G/T to-row — zero
            # weights are kept so the bridge's defensive row-
            # completeness check passes. The engine's
            # inverse-CDF sampler ignores zero-weight entries.
            row_pairs.append((_NP_BASE_BYTE[to_letter], float(weight)))
        rows.append(row_pairs)
    return rows


def _p_nucleotide_lengths_from_models(
    models: Any, key: str,
) -> Optional[List[Tuple[int, float]]]:
    """Read ``p_nucleotide_lengths[key]`` off an explicit
    :class:`ReferenceEmpiricalModels` and lower it into the
    engine's ``[(length, weight), ...]`` shape. Returns
    ``None`` when the spec is absent — the lowering then
    omits the `PAdditionPass` for that end, byte-identical
    to the pre-slice baseline.

    Slice — P-nucleotide v1. Legacy
    `DataConfig.p_nucleotide_length_probs` is NOT auto-lifted
    (same boundary the Markov slice respected for legacy NP
    fields); cartridges author the typed plane explicitly.
    """
    if models is None:
        return None
    spec = getattr(models, "p_nucleotide_lengths", {}).get(key)
    if spec is None:
        return None
    pairs = [(int(v), float(w)) for v, w in spec.values]
    pairs.sort(key=lambda p: p[0])
    return pairs or None


def _allele_usage_from_models(
    models: Any, segment: str,
) -> Optional[Dict[str, float]]:
    """Read ``allele_usage`` off an explicit
    :class:`ReferenceEmpiricalModels` and return the per-segment
    ``{allele_name: weight}`` dict, or ``None`` when:

    - the cartridge has no typed plane authored, or
    - the typed plane's per-segment dict is empty.

    The returned dict is the raw author-supplied weights — NOT
    normalised, NOT validated against the refdata pool. The
    bridge layer in :class:`Experiment._resolve_allele_weights`
    does the pool-name lookup + renormalisation when threading
    into the dense vector the engine consumes.

    Slice — Allele Usage Estimation v1. Legacy
    ``DataConfig.gene_use_dict`` is NOT auto-lifted; the orphan
    boundary holds.
    """
    if models is None:
        return None
    spec = getattr(models, "allele_usage", None)
    if spec is None:
        return None
    attr = {"V": "v", "D": "d", "J": "j"}[segment]
    mapping = getattr(spec, attr, None) or {}
    if not mapping:
        return None
    return dict(mapping)


def _trim_from_models(
    models: Any, key: str, *, cap: Optional[int],
) -> Optional[List[Tuple[int, float]]]:
    """Read ``trims[key]`` off an explicit
    :class:`ReferenceEmpiricalModels`. The per-segment cap (driven
    by the pool's shortest allele) is applied here too so an
    explicit model can't sample a trim larger than the pool
    supports — that would silently break recombination on short
    alleles.
    """
    if models is None:
        return None
    spec = models.trims.get(key)
    if spec is None:
        return None
    pairs = [(int(v), float(w)) for v, w in spec.values]
    if cap is not None:
        pairs = [(v, w) for v, w in pairs if v <= cap]
    if not pairs:
        return None
    pairs.sort(key=lambda p: p[0])
    return pairs
