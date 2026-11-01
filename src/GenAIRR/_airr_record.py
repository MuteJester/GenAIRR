"""Convert a Rust ``Outcome`` into an AIRR-format record dict.

The body of the builder lives in the Rust kernel
(`engine_rs/src/airr_record.rs` → exposed as `Outcome.airr_record()`)
for a >5× speedup at scale. This module is a thin Python wrapper
that delegates to that PyO3 method, plus the column helpers that
`tests/components/experiment` imports directly.

It used to carry a second, complete builder beside the wrapper --
`_legacy_outcome_to_airr_record`, 218 lines, called from nowhere.
That mattered because it answered the same questions differently.
Its productivity rule was

    productive = vj_in_frame and not stop_codon

two of the five predicates the engine actually conjoins, so anyone
reading this file to learn the contract learned a rule the engine
had abandoned. Two estimator pins were also asserting strings
inside it, which made them green on code that could not run.

Field names follow the AIRR Rearrangement schema (MiAIRR).
"""
from __future__ import annotations

from typing import Any, Dict, Optional, Sequence


# Genetic code for junction translation. Matches the engine's
# ``GENETIC_CODE`` constant; replicated here to avoid an extra
# round-trip through the engine for a function that's pure logic.
_GENETIC_CODE = {
    # Phenylalanine
    "TTT": "F", "TTC": "F",
    # Leucine
    "TTA": "L", "TTG": "L", "CTT": "L", "CTC": "L", "CTA": "L", "CTG": "L",
    # Isoleucine
    "ATT": "I", "ATC": "I", "ATA": "I",
    # Methionine / start
    "ATG": "M",
    # Valine
    "GTT": "V", "GTC": "V", "GTA": "V", "GTG": "V",
    # Serine
    "TCT": "S", "TCC": "S", "TCA": "S", "TCG": "S", "AGT": "S", "AGC": "S",
    # Proline
    "CCT": "P", "CCC": "P", "CCA": "P", "CCG": "P",
    # Threonine
    "ACT": "T", "ACC": "T", "ACA": "T", "ACG": "T",
    # Alanine
    "GCT": "A", "GCC": "A", "GCA": "A", "GCG": "A",
    # Tyrosine
    "TAT": "Y", "TAC": "Y",
    # Stop codons
    "TAA": "*", "TAG": "*", "TGA": "*",
    # Histidine
    "CAT": "H", "CAC": "H",
    # Glutamine
    "CAA": "Q", "CAG": "Q",
    # Asparagine
    "AAT": "N", "AAC": "N",
    # Lysine
    "AAA": "K", "AAG": "K",
    # Aspartic acid
    "GAT": "D", "GAC": "D",
    # Glutamic acid
    "GAA": "E", "GAG": "E",
    # Cysteine
    "TGT": "C", "TGC": "C",
    # Tryptophan
    "TGG": "W",
    # Arginine
    "CGT": "R", "CGC": "R", "CGA": "R", "CGG": "R", "AGA": "R", "AGG": "R",
    # Glycine
    "GGT": "G", "GGC": "G", "GGA": "G", "GGG": "G",
}


def _translate_codon(codon: str) -> str:
    """Translate one codon to its single-letter amino acid. Returns
    ``"X"`` for any codon containing a non-A/C/G/T base, ``"*"`` for
    stops. Case-insensitive."""
    upper = codon.upper().replace("U", "T")
    if len(upper) != 3:
        return "X"
    return _GENETIC_CODE.get(upper, "X")


def _alignment_columns(outcome, refdata) -> "list[tuple[str, str, str]]":
    """Walk the simulation IR once and return one tuple per
    alignment column: ``(segment, sa_char, ga_char)``.

    ``segment`` is one of ``"V"``, ``"D"``, ``"J"``, ``"NP1"``,
    ``"NP2"``, or ``""`` for tail-end indel insertions that fell
    outside every region.

    ``sa_char`` is the sequence-alignment character (a base or
    ``"-"`` for a deletion-gap column). ``ga_char`` is the
    germline-alignment character (a base, ``"N"`` for NP, or ``"-"``
    for an insertion-gap column).

    Both alignment-string and CIGAR construction read this single
    column list — no double walks.
    """
    sim = outcome.final_simulation()
    bases_bytes = sim.bases()
    if not bases_bytes:
        return []
    germline_bytes = sim.germline_bases()
    regions = sim.regions()

    def _allele_seq(segment: str) -> Optional[bytes]:
        ids = {
            "V": sim.v_allele_id(),
            "D": sim.d_allele_id(),
            "J": sim.j_allele_id(),
        }
        getter = {
            "V": refdata.v_allele,
            "D": refdata.d_allele,
            "J": refdata.j_allele,
        }
        aid = ids.get(segment)
        if aid is None:
            return None
        try:
            return bytes(getter[segment](aid).seq())
        except (IndexError, AttributeError):
            return None

    # Each coding region's germline positions live in
    # ``[trim_5, allele_len - trim_3)``; ``trim_5`` is the expected
    # starting germline_pos for the walk so we can spot deletion
    # gaps. Our DSL only trims V_3, D_5, D_3, J_5 — V_5 and J_3 are
    # always 0.
    trim_5 = {
        "V": 0,
        "D": _trace_int(outcome, "trim.d_5", 0),
        "J": _trace_int(outcome, "trim.j_5", 0),
    }

    columns: list = []

    # Track which pool positions have been emitted so any tail-end
    # indel insertions that didn't extend a region still get placed in
    # the alignment.
    covered: set = set()

    for region in regions:
        for idx in range(region.start, region.end):
            covered.add(idx)
        seg = region.segment
        if seg in ("V", "D", "J"):
            allele_seq = _allele_seq(seg)
            if allele_seq is None:
                # No source allele assigned — emit raw bases without
                # germline detail to keep lengths aligned.
                for i in range(region.start, region.end):
                    columns.append((seg, chr(bases_bytes[i]), "N"))
                continue
            expected_pos = trim_5[seg]
            # End of the *recombination-surviving* germline range —
            # any position the indel pass later removed inside this
            # range must surface as a ``D`` op.
            end_germ = len(allele_seq) - {
                "V": _trace_int(outcome, "trim.v_3", 0),
                "D": _trace_int(outcome, "trim.d_3", 0),
                "J": 0,
            }[seg]
            for i in range(region.start, region.end):
                germ_pos = sim.germline_position(i)
                base_char = chr(bases_bytes[i])
                if germ_pos is None:
                    # Indel insertion within coding region — gap on
                    # the germline side.
                    columns.append((seg, base_char, "-"))
                else:
                    # Fill any deletion-gap that precedes this base.
                    while expected_pos < germ_pos:
                        columns.append((seg, "-", chr(allele_seq[expected_pos])))
                        expected_pos += 1
                    # ``germline_bytes[i]`` carries the *original*
                    # source base — immune to S5F mutations applied
                    # downstream.
                    columns.append((seg, base_char, chr(germline_bytes[i])))
                    expected_pos = germ_pos + 1
            # Fill any trailing deletion-gap: positions in
            # ``[expected_pos, end_germ)`` were in the assembled segment
            # after recombination but got removed by the indel pass.
            while expected_pos < end_germ:
                columns.append((seg, "-", chr(allele_seq[expected_pos])))
                expected_pos += 1
        elif seg in ("NP1", "NP2"):
            for i in range(region.start, region.end):
                columns.append((seg, chr(bases_bytes[i]), "N"))

    # Tail-end indel insertions (at the very end of the pool) don't
    # always extend a region's end index; they end up outside every
    # region. Emit them as bare insertions.
    for i in range(len(bases_bytes)):
        if i not in covered:
            columns.append(("", chr(bases_bytes[i]), "-"))

    return columns


def _identity_from_columns(
    columns: "list[tuple[str, str, str]]",
) -> "dict[str, Optional[float]]":
    """Per-segment fractional identity, derived from the alignment
    columns we already produced.

    For each V/D/J segment, identity is

        matches / (M + I + D)

    where ``matches`` counts columns in which the case-insensitive
    sequence base equals the germline base. Insertion gaps
    (``ga == "-"``) and deletion gaps (``sa == "-"``) count toward
    the denominator but never toward the numerator. ``N`` germline
    columns shouldn't appear inside V/D/J spans but are skipped from
    matches if they do.

    Returns ``{"V": float | None, "D": float | None, "J": float | None}``.
    A segment is ``None`` when it contributed no alignment columns
    (e.g. ``D`` on a VJ chain).
    """
    matches: dict = {"V": 0, "D": 0, "J": 0}
    totals: dict = {"V": 0, "D": 0, "J": 0}
    for seg, sa_c, ga_c in columns:
        if seg not in ("V", "D", "J"):
            continue
        totals[seg] += 1
        if sa_c == "-" or ga_c == "-" or ga_c == "N":
            continue
        if sa_c.upper() == ga_c.upper():
            matches[seg] += 1
    return {
        seg: (matches[seg] / totals[seg]) if totals[seg] > 0 else None
        for seg in ("V", "D", "J")
    }


def _cigars_from_columns(
    columns: "list[tuple[str, str, str]]",
) -> "tuple[str, str, str]":
    """Build per-segment SAM CIGAR strings from per-column data.

    Each V/D/J alignment column maps to one CIGAR operation:

    - ``M`` — match-or-mismatch (both ``sa`` and ``ga`` carry a base)
    - ``I`` — insertion in query (``ga`` is ``"-"``)
    - ``D`` — deletion from query (``sa`` is ``"-"``)

    Returns ``(v_cigar, d_cigar, j_cigar)``. Empty string when the
    segment isn't assigned. NP-region and orphan columns are not
    represented in any CIGAR — they have no germline reference.
    """
    runs: dict = {"V": [], "D": [], "J": []}
    for seg, sa_c, ga_c in columns:
        if seg not in ("V", "D", "J"):
            continue
        if sa_c == "-":
            op = "D"
        elif ga_c == "-":
            op = "I"
        else:
            op = "M"
        runs[seg].append(op)

    def _runlength(ops: list) -> str:
        if not ops:
            return ""
        out = []
        current = ops[0]
        count = 1
        for op in ops[1:]:
            if op == current:
                count += 1
            else:
                out.append(f"{count}{current}")
                current = op
                count = 1
        out.append(f"{count}{current}")
        return "".join(out)

    return _runlength(runs["V"]), _runlength(runs["D"]), _runlength(runs["J"])


def _aa_slice_for_region(
    region_start: int,
    region_end: int,
    frame_offset: int,
    sequence_aa: str,
) -> str:
    """Return the amino-acid translation of nucleotides ``[region_start,
    region_end)`` interpreted in the global reading frame anchored at
    ``frame_offset`` (V-anchor frame).

    Only complete codons fully contained within the region are
    translated — partial codons at either boundary are dropped, matching
    the AIRR ``np1_aa`` / ``np2_aa`` convention.
    """
    if region_end <= region_start:
        return ""
    # First codon-aligned position at or after `region_start`.
    delta = (frame_offset - region_start) % 3
    codon_start = region_start + delta
    if codon_start >= region_end:
        return ""
    n_codons = (region_end - codon_start) // 3
    if n_codons <= 0:
        return ""
    aa_index = (codon_start - frame_offset) // 3
    return sequence_aa[aa_index : aa_index + n_codons]


# AIRR ``locus`` is one of the seven canonical 3-letter codes. We
# derive it from the V allele name's prefix (IGH..., TRBV..., etc.)
# since every config we ship uses IMGT naming, falling back to J then
# D and finally an empty string when no segment is assigned.
_AIRR_LOCI = ("IGH", "IGK", "IGL", "TRA", "TRB", "TRG", "TRD")


def _derive_locus(v_call: str, j_call: str, d_call: str) -> str:
    """Pick the AIRR ``locus`` value from the assigned allele names.

    Returns one of :data:`_AIRR_LOCI` when any segment name starts
    with a canonical prefix, else an empty string.
    """
    for name in (v_call, j_call, d_call):
        if not name:
            continue
        prefix = name[:3].upper()
        if prefix in _AIRR_LOCI:
            return prefix
    return ""


def _trace_int(outcome, address: str, default: int = 0) -> int:
    """Read a single integer trace value, or ``default`` when the
    address isn't recorded."""
    rec = outcome.trace().find(address)
    if rec is None:
        return default
    val = rec.value
    return int(val) if isinstance(val, int) else default


def outcome_to_airr_record(
    outcome,
    refdata,
    *,
    sequence_id: Optional[str] = None,
) -> Dict[str, Any]:
    """Build an AIRR Rearrangement record dict via the Rust kernel.

    Thin wrapper around ``Outcome.airr_record(refdata, sequence_id=...)``
    — the Rust path produces the canonical 69-field record in one pass
    over the IR. See ``engine_rs/src/airr_record.rs`` for the
    implementation. The function signature is preserved so callers
    that imported this directly keep working.
    """
    return outcome.airr_record(
        refdata, sequence_id=sequence_id if sequence_id is not None else ""
    )

