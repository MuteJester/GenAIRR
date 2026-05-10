"""Convert a Rust ``Outcome`` into an AIRR-format record dict.

The body of the builder lives in the Rust kernel
(`engine_rs/src/airr_record.rs` → exposed as `Outcome.airr_record()`)
for a >5× speedup at scale. This module is a thin Python wrapper
that delegates to that PyO3 method, plus a few helpers kept around
for tests and any direct importers.

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


def _translate_seq(seq: str) -> str:
    """Translate a nucleotide string codon-by-codon, ignoring any
    trailing partial codon."""
    return "".join(
        _translate_codon(seq[i : i + 3]) for i in range(0, len(seq) - 2, 3)
    )


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


def _strings_from_columns(
    columns: "list[tuple[str, str, str]]",
) -> "tuple[str, str, str]":
    """Compose the AIRR alignment-string triplet from per-column data.

    Returns ``(sequence_alignment, germline_alignment,
    germline_alignment_d_mask)`` — uppercase, gap-aware. D-mask
    replaces D-region germline bases with ``"N"`` and preserves
    ``"-"`` gaps.
    """
    sa = []
    ga = []
    dm = []
    for seg, sa_c, ga_c in columns:
        sa.append(sa_c)
        ga.append(ga_c)
        if seg == "D" and ga_c != "-":
            dm.append("N")
        else:
            dm.append(ga_c)
    return "".join(sa).upper(), "".join(ga).upper(), "".join(dm).upper()


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


def _coord_pairs_from_columns(
    columns: "list[tuple[str, str, str]]",
    outcome,
    refdata,
) -> "dict[str, tuple[int, int, int, int]]":
    """Compute per-segment ``(alignment_start, alignment_end,
    germline_start, germline_end)`` tuples.

    All coordinates are **0-based half-open** for consistency with
    the existing ``v_sequence_start`` / ``v_sequence_end`` fields.
    An ``airr_strict`` export flag (in :mod:`result`) converts to
    the 1-based-inclusive convention the AIRR spec requires.

    - ``alignment_start`` / ``alignment_end`` are positions in the
      gap-aware alignment string (``sa`` / ``ga``). With deletion
      gaps, these can exceed the matching ``*_sequence_end``.
    - ``germline_start`` / ``germline_end`` are positions in the
      source allele sequence. Derived from trim values: V uses
      ``[0, allele_len - v_trim_3)``; D uses ``[d_trim_5,
      allele_len - d_trim_3)``; J uses ``[j_trim_5, allele_len)``.

    Segments with no assigned allele are omitted from the result.
    """
    sim = outcome.final_simulation()

    # Walk columns once to find each V/D/J segment's alignment-string
    # span. Columns are emitted in order, and each segment is one
    # contiguous run, so first-and-last suffice.
    align_ranges: dict = {}
    for i, (seg, _, _) in enumerate(columns):
        if seg not in ("V", "D", "J"):
            continue
        if seg not in align_ranges:
            align_ranges[seg] = [i, i + 1]
        else:
            align_ranges[seg][1] = i + 1

    trim_5 = {
        "V": 0,
        "D": _trace_int(outcome, "trim.d_5", 0),
        "J": _trace_int(outcome, "trim.j_5", 0),
    }
    trim_3 = {
        "V": _trace_int(outcome, "trim.v_3", 0),
        "D": _trace_int(outcome, "trim.d_3", 0),
        "J": 0,
    }

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

    result: dict = {}
    for seg, (a_start, a_end) in align_ranges.items():
        aid = ids[seg]
        if aid is None:
            continue
        try:
            allele_len = len(bytes(getter[seg](aid).seq()))
        except (IndexError, AttributeError):
            continue
        g_start = trim_5[seg]
        g_end = allele_len - trim_3[seg]
        result[seg] = (a_start, a_end, g_start, g_end)
    return result


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


def _allele_name(refdata, segment: str, allele_id: Optional[int]) -> str:
    """Look up the IMGT-style allele name for an ID, or return an
    empty string when the segment isn't assigned (e.g. ``d_call`` on
    a VJ chain)."""
    if allele_id is None:
        return ""
    if segment == "V":
        return refdata.v_allele(allele_id).name
    if segment == "D":
        return refdata.d_allele(allele_id).name
    if segment == "J":
        return refdata.j_allele(allele_id).name
    return ""


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


def _allele_anchor(refdata, segment: str, allele_id: Optional[int]) -> Optional[int]:
    """Anchor position in the source allele for the assigned V/J
    allele. ``None`` if the allele is anchorless or not assigned."""
    if allele_id is None:
        return None
    try:
        if segment == "V":
            allele = refdata.v_allele(allele_id)
        elif segment == "J":
            allele = refdata.j_allele(allele_id)
        else:
            return None
    except IndexError:
        return None
    return allele.anchor


def _trace_int(outcome, address: str, default: int = 0) -> int:
    """Read a single integer trace value, or ``default`` when the
    address isn't recorded."""
    rec = outcome.trace().find(address)
    if rec is None:
        return default
    val = rec.value
    return int(val) if isinstance(val, int) else default


def _trace_bool(outcome, address: str, default: bool = False) -> bool:
    """Read a single boolean trace value, or ``default`` when not
    recorded."""
    rec = outcome.trace().find(address)
    if rec is None:
        return default
    return bool(rec.value)


def _region_for_segment(regions: Sequence, segment: str):
    """Return the first region with ``segment == segment``, or
    ``None``. Used for V / D / J / NP1 / NP2 lookup."""
    for r in regions:
        if r.segment == segment:
            return r
    return None


def _has_stop(seq: str) -> bool:
    """``True`` when a codon-aligned walk of ``seq`` hits a stop."""
    for i in range(0, len(seq) - 2, 3):
        codon = seq[i : i + 3].upper().replace("U", "T")
        if codon in ("TAA", "TAG", "TGA"):
            return True
    return False


def _mutation_count(outcome) -> int:
    """Total mutations from S5F / Uniform passes (whichever ran)."""
    s5f = outcome.trace().find("mutate.s5f.count")
    if s5f is not None and isinstance(s5f.value, int):
        return int(s5f.value)
    uniform = outcome.trace().find("mutate.uniform.count")
    if uniform is not None and isinstance(uniform.value, int):
        return int(uniform.value)
    return 0


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


def _legacy_outcome_to_airr_record(
    outcome,
    refdata,
    *,
    sequence_id: Optional[str] = None,
) -> Dict[str, Any]:
    """Build an AIRR Rearrangement record dict from a Rust ``Outcome``.

    ``refdata`` is the engine's ``RefDataConfig`` — used to look up
    allele names + anchor positions for the junction window.

    ``sequence_id`` is the AIRR per-record identifier; pass an
    explicit string to override the default (the empty string).
    """
    sim = outcome.final_simulation()
    bases_bytes = sim.bases()
    sequence = bases_bytes.decode("ascii", errors="replace") if bases_bytes else ""
    regions = sim.regions()

    v_id = sim.v_allele_id()
    d_id = sim.d_allele_id()
    j_id = sim.j_allele_id()

    v_region = _region_for_segment(regions, "V")
    d_region = _region_for_segment(regions, "D")
    j_region = _region_for_segment(regions, "J")
    np1_region = _region_for_segment(regions, "NP1")
    np2_region = _region_for_segment(regions, "NP2")

    # Trim values default to 0 — recombine(trim=False) plans don't
    # record trim.* addresses in the trace.
    v_trim_3 = _trace_int(outcome, "trim.v_3", 0)
    d_trim_5 = _trace_int(outcome, "trim.d_5", 0)
    d_trim_3 = _trace_int(outcome, "trim.d_3", 0)
    j_trim_5 = _trace_int(outcome, "trim.j_5", 0)

    # Junction window: V anchor → J anchor + 3 (inclusive of the
    # J W/F codon). Computable only when both V and J are assigned
    # *and* both have anchors *and* both regions are in the pool.
    junction_start: Optional[int] = None
    junction_end: Optional[int] = None
    junction_nt = ""
    junction_aa = ""
    junction_length: Optional[int] = None
    vj_in_frame: Optional[bool] = None
    stop_codon: Optional[bool] = None

    v_anchor = _allele_anchor(refdata, "V", v_id)
    j_anchor = _allele_anchor(refdata, "J", j_id)
    if (
        v_region is not None
        and j_region is not None
        and v_anchor is not None
        and j_anchor is not None
    ):
        v_trim_5 = 0  # current Experiment.recombine doesn't trim V_5
        j_anchor_in_pool = j_region.start + (j_anchor - j_trim_5)
        v_anchor_in_pool = v_region.start + (v_anchor - v_trim_5)
        if j_anchor_in_pool + 3 > v_anchor_in_pool:
            junction_start = v_anchor_in_pool
            junction_end = j_anchor_in_pool + 3
            junction_nt = sequence[junction_start:junction_end]
            junction_length = len(junction_nt)
            vj_in_frame = junction_length % 3 == 0
            stop_codon = _has_stop(junction_nt) if vj_in_frame else False
            junction_aa = _translate_seq(junction_nt) if vj_in_frame else ""

    productive: Optional[bool] = None
    if vj_in_frame is not None and stop_codon is not None:
        productive = vj_in_frame and not stop_codon

    v_call = _allele_name(refdata, "V", v_id)
    d_call = _allele_name(refdata, "D", d_id)
    j_call = _allele_name(refdata, "J", j_id)
    locus = _derive_locus(v_call, j_call, d_call)

    # Reading frame for ``sequence_aa`` / ``np*_aa`` is anchored at
    # the V-anchor — same frame the junction is computed in. When the
    # junction can't be derived (no V/J anchor) we leave the AA fields
    # empty rather than guessing a frame.
    sequence_aa = ""
    np1_aa = ""
    np2_aa = ""
    if junction_start is not None:
        frame_offset = junction_start % 3
        sequence_aa = _translate_seq(sequence[frame_offset:])
        if np1_region is not None:
            np1_aa = _aa_slice_for_region(
                np1_region.start, np1_region.end, frame_offset, sequence_aa
            )
        if np2_region is not None:
            np2_aa = _aa_slice_for_region(
                np2_region.start, np2_region.end, frame_offset, sequence_aa
            )

    # Ground-truth alignment columns — derived from the engine's
    # per-nucleotide germline provenance instead of re-aligning the
    # output. The per-column list is consumed twice: once to compose
    # the AIRR alignment strings, once to compose the per-segment
    # CIGARs.
    columns = _alignment_columns(outcome, refdata)
    (
        sequence_alignment,
        germline_alignment,
        germline_alignment_d_mask,
    ) = _strings_from_columns(columns)
    v_cigar, d_cigar, j_cigar = _cigars_from_columns(columns)
    coord_pairs = _coord_pairs_from_columns(columns, outcome, refdata)
    v_align = coord_pairs.get("V", (None, None, None, None))
    d_align = coord_pairs.get("D", (None, None, None, None))
    j_align = coord_pairs.get("J", (None, None, None, None))

    # Identity: ground truth, computed from the column walk. Score
    # and support are aligner-specific and don't apply to a
    # simulator → stubbed as ``None``. The columns still appear in
    # output so AIRR-tooling parsers see the canonical schema shape.
    identities = _identity_from_columns(columns)
    v_identity = identities["V"]
    d_identity = identities["D"]
    j_identity = identities["J"]

    record: Dict[str, Any] = {
        # AIRR metadata
        "sequence_id": sequence_id if sequence_id is not None else "",
        "sequence": sequence,
        "sequence_aa": sequence_aa,
        "sequence_alignment": sequence_alignment,
        "germline_alignment": germline_alignment,
        "germline_alignment_d_mask": germline_alignment_d_mask,
        "sequence_length": len(sequence),
        "rev_comp": False,
        "locus": locus,
        # Calls
        "v_call": v_call,
        "v_cigar": v_cigar,
        "v_score": None,
        "v_identity": v_identity,
        "v_support": None,
        "v_sequence_start": v_region.start if v_region is not None else None,
        "v_sequence_end": v_region.end if v_region is not None else None,
        "v_alignment_start": v_align[0],
        "v_alignment_end": v_align[1],
        "v_germline_start": v_align[2],
        "v_germline_end": v_align[3],
        "v_trim_5": 0,
        "v_trim_3": v_trim_3,
        "d_call": d_call,
        "d_cigar": d_cigar,
        "d_score": None,
        "d_identity": d_identity,
        "d_support": None,
        "d_sequence_start": d_region.start if d_region is not None else None,
        "d_sequence_end": d_region.end if d_region is not None else None,
        "d_alignment_start": d_align[0],
        "d_alignment_end": d_align[1],
        "d_germline_start": d_align[2],
        "d_germline_end": d_align[3],
        "d_trim_5": d_trim_5,
        "d_trim_3": d_trim_3,
        "j_call": j_call,
        "j_cigar": j_cigar,
        "j_score": None,
        "j_identity": j_identity,
        "j_support": None,
        "j_sequence_start": j_region.start if j_region is not None else None,
        "j_sequence_end": j_region.end if j_region is not None else None,
        "j_alignment_start": j_align[0],
        "j_alignment_end": j_align[1],
        "j_germline_start": j_align[2],
        "j_germline_end": j_align[3],
        "j_trim_5": j_trim_5,
        "j_trim_3": 0,
        "c_call": "",
        # Junction
        "junction": junction_nt,
        "junction_aa": junction_aa,
        "junction_start": junction_start,
        "junction_end": junction_end,
        "junction_length": junction_length,
        # NP regions
        "np1": (
            sequence[np1_region.start : np1_region.end] if np1_region is not None else ""
        ),
        "np1_aa": np1_aa,
        "np1_length": np1_region.end - np1_region.start if np1_region is not None else 0,
        "np2": (
            sequence[np2_region.start : np2_region.end] if np2_region is not None else ""
        ),
        "np2_aa": np2_aa,
        "np2_length": np2_region.end - np2_region.start if np2_region is not None else 0,
        # Functionality
        "productive": productive,
        "vj_in_frame": vj_in_frame,
        "stop_codon": stop_codon,
        # SHM
        "n_mutations": _mutation_count(outcome),
        # Corruption annotations
        "n_pcr_errors": _trace_int(outcome, "corrupt.pcr.count", 0),
        "n_quality_errors": _trace_int(outcome, "corrupt.quality.count", 0),
        "n_indels": _trace_int(outcome, "corrupt.indel.count", 0),
        "is_contaminant": _trace_bool(outcome, "corrupt.contaminant.applied", False),
    }

    seq_len = record["sequence_length"]
    record["mutation_rate"] = (
        record["n_mutations"] / seq_len if seq_len > 0 else 0.0
    )
    return record
