import csv
import logging
import re
from collections import defaultdict
from ..alleles.allele import VAllele, JAllele, DAllele, CAllele
from ..utilities import parse_fasta

logger = logging.getLogger(__name__)


# ─────────────────────────────────────────────────────────────
# Default anchor finders (standard IMGT conventions)
# ─────────────────────────────────────────────────────────────

def default_v_anchor_finder(gapped_seq, ungapped_seq):
    """Return the V anchor position (int), or None if not found.

    Checks IMGT position 104 (gapped nucleotide 309-312) for a Cysteine
    codon (TGT/TGC). This is the conserved Cys at the start of CDR3 in
    all standard mammalian BCR and TCR V genes.
    """
    cys = gapped_seq[309:312]
    if cys not in ("tgt", "tgc"):
        return None
    cys_wider = gapped_seq[306:315]
    pos = ungapped_seq.rfind(cys_wider)
    if pos == -1:
        return None
    return pos + 3


def default_j_anchor_finder(ungapped_seq):
    """Return the J anchor position (int), or None if not found.

    Searches for the conserved Phe/Trp-Gly-X-Gly motif across all three
    reading frames.  This motif marks the end of CDR3 in standard
    mammalian BCR and TCR J genes.
    """
    motif = re.compile(r'(ttt|ttc|tgg)(ggt|ggc|gga|ggg)[a-z]{3}(ggt|ggc|gga|ggg)')
    for frame in range(3):
        match = motif.search(ungapped_seq[frame:])
        if match and match.span()[0] % 3 == 0:
            return match.span()[0] + frame
    return None


# ─────────────────────────────────────────────────────────────
# Segment detection from gene names
# ─────────────────────────────────────────────────────────────

def _detect_segment_from_name(allele_name):
    """Detect segment type (V/D/J/C) from an IMGT-style gene name.

    Scans the family prefix backward from the last non-digit character
    to find the segment indicator letter.

    Examples::

        IGHV1-2*01  -> V
        TRBD1*01    -> D
        TRBJ2-1*01  -> J
        IGKC*01     -> C
        TRAV12-1*01 -> V

    Returns:
        One of ``'V'``, ``'D'``, ``'J'``, ``'C'``, or ``None``.
    """
    if "-" in allele_name:
        family = allele_name.split("-")[0]
    else:
        family = allele_name.split("*")[0]

    # Strip trailing digits to get the alphabetic prefix (e.g., "IGHV1" -> "IGHV")
    prefix = family.upper().rstrip("0123456789")
    if prefix and prefix[-1] in ('V', 'D', 'J', 'C'):
        return prefix[-1]

    return None


# ─────────────────────────────────────────────────────────────
# Main FASTA loading function
# ─────────────────────────────────────────────────────────────

def create_allele_dict(
    fasta,
    *,
    segment_type=None,
    v_anchor_finder=None,
    j_anchor_finder=None,
    keep_anchorless=False,
):
    """
    Load alleles from an IMGT-formatted FASTA reference file.

    Constructs a dictionary mapping gene names to lists of Allele objects.
    Only functional (``F``) and open reading frame (``ORF``) alleles that
    pass anchor validation are included.

    Args:
        fasta: Path to a FASTA file with IMGT-gapped sequences.
        segment_type: Explicit segment type (``'V'``, ``'D'``, ``'J'``,
            or ``'C'``).  When provided, all alleles in the file are
            treated as this segment type instead of inferring from names.
        v_anchor_finder: Callable ``(gapped_seq, ungapped_seq) -> int | None``.
            Returns the V anchor position or ``None`` to reject.  Defaults
            to :func:`default_v_anchor_finder` (IMGT Cys-104).
        j_anchor_finder: Callable ``(ungapped_seq) -> int | None``.
            Returns the J anchor position or ``None`` to reject.  Defaults
            to :func:`default_j_anchor_finder` (Phe/Trp-Gly-X-Gly motif).
        keep_anchorless: If True, V/J alleles whose anchor cannot be found
            are kept with ``anchor=0`` instead of being rejected.  Such
            alleles will never produce productive rearrangements but are
            retained in the reference for completeness.

    Returns:
        dict: ``{gene_name: [Allele, ...]}`` mapping.
    """
    if v_anchor_finder is None:
        v_anchor_finder = default_v_anchor_finder
    if j_anchor_finder is None:
        j_anchor_finder = default_j_anchor_finder

    allele_dict = defaultdict(list)
    allele_class_map = {'V': VAllele, 'D': DAllele, 'J': JAllele, 'C': CAllele}

    # Rejection tracking
    total = 0
    accepted = 0
    rejected_anchor = 0
    rejected_partial = 0
    rejected_functional = 0
    rejected_unknown = 0
    kept_anchorless = 0

    with open(fasta) as f:
        for header, seq in parse_fasta(f):
            total += 1
            header = header.replace(">", "")
            seq = seq.lower()

            # Extract allele name (IMGT pipe-delimited or plain)
            parts = header.split("|")
            allele_name = parts[1] if len(parts) > 1 else parts[0]

            # --- Determine segment type ---
            if segment_type is not None:
                segment = segment_type
            else:
                segment = _detect_segment_from_name(allele_name)

            if segment is None:
                logger.debug("Skipping %s: could not determine segment type", allele_name)
                rejected_unknown += 1
                continue

            # --- Functional status filter ---
            coding = parts[3] if len(parts) > 3 else "F"
            if "partial" in header:
                logger.debug("Skipping partial allele %s", allele_name)
                rejected_partial += 1
                continue
            if coding not in ("F", "ORF"):
                logger.debug("Skipping non-functional allele %s (coding=%s)", allele_name, coding)
                rejected_functional += 1
                continue

            # --- Segment-specific processing ---
            ungapped_seq = seq.replace(".", "")
            ungapped_length = len(ungapped_seq)
            gene = allele_name.split("*")[0]

            if segment == 'D':
                seq = ungapped_seq  # D segments stored ungapped

            if segment == 'V':
                anchor = v_anchor_finder(seq, ungapped_seq)
                if anchor is None:
                    if keep_anchorless:
                        logger.debug(
                            "Keeping anchorless V allele %s (anchor=0)",
                            allele_name,
                        )
                        anchor = 0
                        kept_anchorless += 1
                    else:
                        logger.debug(
                            "Rejecting V allele %s: anchor not found (gapped[309:312]=%r)",
                            allele_name, seq[309:312] if len(seq) > 312 else "<short>",
                        )
                        rejected_anchor += 1
                        continue
                allele_dict[gene].append(
                    VAllele(allele_name, seq, ungapped_length, anchor_override=anchor)
                )
            elif segment == 'J':
                anchor = j_anchor_finder(ungapped_seq)
                if anchor is None:
                    if keep_anchorless:
                        logger.debug(
                            "Keeping anchorless J allele %s (anchor=0)",
                            allele_name,
                        )
                        anchor = 0
                        kept_anchorless += 1
                    else:
                        logger.debug("Rejecting J allele %s: anchor motif not found", allele_name)
                        rejected_anchor += 1
                        continue
                allele_dict[gene].append(
                    JAllele(allele_name, seq, ungapped_length, anchor_override=anchor)
                )
            else:
                allele_dict[gene].append(
                    allele_class_map[segment](allele_name, seq, ungapped_length)
                )

            accepted += 1

    # --- Summary logging ---
    rejected_total = total - accepted
    anchorless_msg = f", {kept_anchorless} anchorless" if kept_anchorless else ""
    logger.info(
        "FASTA %s: loaded %d/%d alleles (rejected: %d anchor, %d partial, "
        "%d non-functional, %d unknown-segment%s)",
        fasta, accepted, total,
        rejected_anchor, rejected_partial, rejected_functional, rejected_unknown,
        anchorless_msg,
    )
    if accepted == 0 and total > 0:
        logger.warning(
            "No alleles accepted from %s (%d total in file). "
            "Check that anchor positions and functional status match your data.",
            fasta, total,
        )

    return dict(allele_dict)
