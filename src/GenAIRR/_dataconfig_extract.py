"""Extract empirical distributions from a :class:`DataConfig` into the
``[(value, weight), ...]`` shape that ``genairr_engine.PassPlan`` expects.

A :class:`DataConfig` carries per-species reference distributions
that the V5 simulation used directly. V6's `recombine()` reads these
through the helpers in this module so simulations default to the
species' empirical NP lengths and trim amounts instead of the uniform
``[0..6]`` placeholder.

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
  dispatch in the engine. Per-gene trims land in a follow-up phase.
"""
from __future__ import annotations

from collections import defaultdict
from typing import Any, Dict, List, Optional, Tuple


# Per-segment trim-dict keys in `DataConfig.trim_dicts`. Maps the
# Rust trace-address fragment ``(segment, end)`` to the V5
# DataConfig key. Keep the mapping explicit so readers don't have to
# infer the V5 naming convention.
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

    The DataConfig stores trims as
    ``{family: {gene: {trim: prob}}}`` keyed by the V5 string form
    (``"V_3"`` / ``"D_5"`` / ``"D_3"`` / ``"J_5"``). We marginalize
    by averaging the per-gene distributions: for each trim amount,
    weight = sum-over-genes(prob[trim]) / number-of-genes.

    ``cap`` (optional) clamps the support so trims past that value
    are filtered out. Used to keep the marginalized distribution
    safe against short alleles in the pool (V5's per-gene runtime
    filter is reduced to a compile-time cap here — see
    :func:`_trim_cap`).
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
    ``[(value, weight), ...]`` list or ``None``. Trim distributions
    are clamped per-segment based on the pool's shortest allele so
    sampled trims always leave a viable assembled slice, including
    the two-sided D case where ``D_5 + D_3`` could otherwise
    overshoot the allele.
    """
    cap_v = _trim_cap(_min_allele_length(getattr(cfg, "v_alleles", None)), two_sided=False)
    cap_d = _trim_cap(_min_allele_length(getattr(cfg, "d_alleles", None)), two_sided=True)
    cap_j = _trim_cap(_min_allele_length(getattr(cfg, "j_alleles", None)), two_sided=False)

    return {
        "np1": extract_np_lengths(cfg, "NP1"),
        "np2": extract_np_lengths(cfg, "NP2"),
        "trim_v_3": extract_trim_distribution(cfg, "V", "3", cap=cap_v),
        "trim_d_5": extract_trim_distribution(cfg, "D", "5", cap=cap_d),
        "trim_d_3": extract_trim_distribution(cfg, "D", "3", cap=cap_d),
        "trim_j_5": extract_trim_distribution(cfg, "J", "5", cap=cap_j),
    }
