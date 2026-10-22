"""
IMGT region boundary mapper.

Maps IMGT-defined CDR/FWR boundaries from V allele gapped sequences
to positions in the assembled rearranged sequence.

The V allele's ``gapped_seq`` uses IMGT numbering with dots (``.``) as
gap characters.  IMGT defines fixed gapped nucleotide boundaries so that
the same gapped position always corresponds to the same structural
location, regardless of allele-specific CDR length variation.

This is the same principle used by ``VAllele._find_anchor()`` which reads
``gapped_seq[306:315]`` to locate the conserved Cysteine at IMGT
position 104.
"""

from __future__ import annotations

from typing import Dict, Optional, Tuple

# Standard IMGT gapped nucleotide boundaries (0-indexed, half-open).
# These are invariant across all V alleles and species.
IMGT_GAPPED_BOUNDARIES = {
    'FWR1': (0, 78),     # IMGT aa 1-26
    'CDR1': (78, 114),   # IMGT aa 27-38
    'FWR2': (114, 165),  # IMGT aa 39-55
    'CDR2': (165, 195),  # IMGT aa 56-65
    'FWR3': (195, 312),  # IMGT aa 66-104
}


def gapped_to_ungapped_pos(gapped_seq: str, gapped_pos: int) -> int:
    """Convert an IMGT gapped nucleotide position to the corresponding
    ungapped (actual DNA) position by counting non-gap characters.

    Args:
        gapped_seq: The IMGT-gapped nucleotide sequence (dots for gaps).
        gapped_pos: 0-indexed position in the gapped sequence.

    Returns:
        The number of non-gap characters before ``gapped_pos``.
    """
    return sum(1 for c in gapped_seq[:gapped_pos] if c != '.')


def declares_imgt_gaps(gapped_seq: Optional[str]) -> bool:
    """Whether ``gapped_seq`` declares an IMGT-numbered alignment.

    Follows the ``VAllele._find_anchor()`` convention: IMGT gap dots are
    the declaration. Ungapped FASTA stored as ``gapped_seq`` has no
    numbering alignment, so fixed IMGT positions cannot be mapped onto it.
    """
    return bool(gapped_seq) and '.' in gapped_seq


def imgt_region_boundaries(gapped_seq: str) -> Dict[str, Tuple[int, int]]:
    """Map the fixed IMGT gapped boundaries onto ungapped coordinates.

    Performs no declaration check: the caller must already know that
    ``gapped_seq`` is IMGT-numbered.  Returns ``(ungapped_start, ungapped_end)``
    half-open intervals for each region.
    """
    result = {}
    for region, (g_start, g_end) in IMGT_GAPPED_BOUNDARIES.items():
        # Clamp to actual gapped sequence length
        g_start_safe = min(g_start, len(gapped_seq))
        g_end_safe = min(g_end, len(gapped_seq))
        ug_start = gapped_to_ungapped_pos(gapped_seq, g_start_safe)
        ug_end = gapped_to_ungapped_pos(gapped_seq, g_end_safe)
        result[region] = (ug_start, ug_end)
    return result


def compute_v_region_boundaries(v_allele) -> Dict[str, Tuple[int, int]]:
    """Compute ungapped V-allele-relative boundaries for each IMGT region.

    Uses the allele's gapped sequence to map fixed IMGT gapped positions
    to ungapped (actual DNA) positions.

    Args:
        v_allele: A ``VAllele`` instance with a ``gapped_seq`` attribute.

    Returns:
        Dict mapping region name (``'FWR1'``, ``'CDR1'``, etc.) to
        ``(ungapped_start, ungapped_end)`` half-open intervals within
        the V allele's ungapped sequence.

    Raises:
        ValueError: if ``gapped_seq`` does not declare IMGT gaps.
    """
    gapped = v_allele.gapped_seq
    if not declares_imgt_gaps(gapped):
        raise ValueError(
            f"{v_allele.name}: gapped_seq declares no IMGT gap alignment; "
            "IMGT subregions cannot be derived from an ungapped sequence"
        )
    return imgt_region_boundaries(gapped)


def classify_position(
    pos: int,
    v_seq_start: int,
    v_seq_end: int,
    v_boundaries: Dict[str, Tuple[int, int]],
    junction_start: int,
    junction_end: int,
    j_seq_end: int,
) -> str:
    """Classify an assembled-sequence position into an IMGT region.

    Priority order:
    1. CDR3 (junction boundaries)
    2. FWR4 (J region after junction)
    3. V region sub-regions (FWR1/CDR1/FWR2/CDR2/FWR3)
    4. NP (everything else: NP regions, D region)

    Args:
        pos: 0-indexed position in the assembled sequence.
        v_seq_start: Start of V segment in assembled sequence.
        v_seq_end: End of V segment in assembled sequence.
        v_boundaries: Dict from ``compute_v_region_boundaries()``.
        junction_start: Start of CDR3/junction in assembled sequence.
        junction_end: End of CDR3/junction in assembled sequence.
        j_seq_end: End of J segment in assembled sequence.

    Returns:
        One of ``'FWR1'``, ``'CDR1'``, ``'FWR2'``, ``'CDR2'``,
        ``'FWR3'``, ``'CDR3'``, ``'FWR4'``, or ``'NP'``.
    """
    # CDR3 (junction) takes priority
    if junction_start <= pos < junction_end:
        return 'CDR3'

    # FWR4 (J region after junction end)
    if junction_end <= pos < j_seq_end:
        return 'FWR4'

    # V region sub-regions
    if v_seq_start <= pos < v_seq_end:
        # Map assembled position to ungapped V allele position.
        # v_germline_start is always 0 for V alleles (never 5'-trimmed).
        v_ungapped_pos = pos - v_seq_start
        for region, (ug_start, ug_end) in v_boundaries.items():
            if ug_start <= v_ungapped_pos < ug_end:
                return region

    return 'NP'
