"""Check the projected alignment against the reference it claims to come from.

The projection fields -- the cigars and the two alignment strings -- are what a
downstream annotator is scored against when GenAIRR is used as the oracle. The
existing tests establish their shape: that only ``M``/``I``/``D`` appear, that
run-lengths collapse, that op counts match the alignment columns, that a VJ
chain has an empty ``d_cigar``. None of them reads the reference, and the
strongest one re-derives the cigar by calling the production helpers and says
in its own comment that "equivalence is guaranteed by construction".

So this module crosses to the reference. Seven identities, each measured before
being written down:

1. ``len(sequence_alignment) == len(germline_alignment)``
2. ``len(sequence_alignment) == len(sequence) + ΣD`` over all three cigars
3. per segment, ``Σ(M, I) == sequence_end - sequence_start``
4. per segment, ``Σ(M, D) == germline_end - germline_start``
5. per segment, walking the cigar from ``alignment_start``: every ``M``
   column's ``germline_alignment`` base is the reference base at the germline
   coordinate the walk has reached
6. ``sequence_alignment`` with gaps removed is the read, case-insensitively --
   or the read's reverse complement when ``rev_comp`` is set, because the
   alignment is written in germline orientation
7. ``germline_alignment_d_mask`` is ``germline_alignment`` with the bases in
   the D span replaced by ``N``, and the gaps in that span left alone

Identities 3 and 4 tie the cigar to the coordinate fields, which are a
different producer. Identity 5 is the one that matters: it crosses the cigar,
the alignment coordinates, the germline alignment and the reference table, so a
wrong allele, a shifted alignment, or a cigar whose op lengths are right and
whose placement is wrong each break it.

**Nothing is skipped.** An inverted D is checked rather than passed over: the
germline alignment carries the reverse complement for an inverted segment, as
it must to align with the read, so identity 5 applies with the expectation
reverse-complemented within that segment's own germline span. Measured, that
matches 320 of 320 columns where the plain reference matches 62. The
base-level oracle next door cannot do this, because the IR holds forward
coordinates with complemented bases.

Validated over 330 records across IGH, IGK and TCRB and the plain, indel,
mutated, reverse-complemented and inverted-D stacks.

Imports nothing from GenAIRR: the reference sequences and the set of inverted
segments are handed in by the caller.
"""
from __future__ import annotations

import re
from dataclasses import dataclass

SEGMENTS = ("v", "d", "j")
GAP = "-"
COMPLEMENT = {"A": "T", "T": "A", "C": "G", "G": "C", "N": "N", GAP: GAP}

_CIGAR_OP = re.compile(r"(\d+)([A-Z=])")


@dataclass(frozen=True)
class Disagreement:
    """One identity that does not hold, and where."""

    sequence_id: str
    segment: str
    identity: str
    recorded: object
    expected: object

    def __str__(self) -> str:
        return (f"{self.sequence_id}: {self.segment} {self.identity}: "
                f"recorded {self.recorded!r}, expected {self.expected!r}")


def cigar_ops(cigar) -> tuple[tuple[int, str], ...]:
    """``((length, op), ...)`` for a cigar string, empty for ``None`` or ``""``."""
    if not cigar:
        return ()
    return tuple((int(length), op) for length, op in _CIGAR_OP.findall(cigar))


def reverse_complement(sequence: str) -> str:
    return "".join(COMPLEMENT.get(base.upper(), "N") for base in reversed(sequence))


def _count(record, field: str) -> int:
    value = record.get(field)
    return 0 if value is None else int(value)


def _expected_germline_base(reference: str, position: int, *, segment_start: int,
                            segment_end: int, inverted: bool) -> str | None:
    """The base the germline alignment should carry at ``position``.

    For a forward segment that is simply the reference's base. For an inverted
    one the germline alignment carries the span reverse-complemented, so the
    base at the i-th position of the span is the complement of the reference's
    base at the span's (len - 1 - i)-th position.
    """
    if inverted:
        position = segment_end - 1 - (position - segment_start)
    if not 0 <= position < len(reference):
        return None
    base = reference[position]
    return COMPLEMENT.get(base, "N") if inverted else base


def _walk(record, segment: str):
    """``(column, germline_position, op)`` per cigar op, from ``alignment_start``.

    A ``D`` advances both the alignment column and the germline coordinate --
    the germline base is present and the read is not. An ``I`` advances only
    the column. Anything else is reported by the caller as an unknown op.
    """
    column = record[f"{segment}_alignment_start"]
    position = record[f"{segment}_germline_start"]
    for length, op in cigar_ops(record.get(f"{segment}_cigar")):
        yield column, position, length, op
        if op in ("M", "=", "X", "D"):
            column += length
            position += length
        elif op == "I":
            column += length
        else:
            column += length


def disagreements(record, allele_sequences, *, inverted=(), segments=SEGMENTS
                  ) -> tuple[Disagreement, ...]:
    """Identities in this record's projection that do not hold.

    ``allele_sequences`` maps a segment letter to the reference sequence of the
    allele this record used, upper case. ``inverted`` names the segments the
    simulation assembled in reverse orientation. Passing both in rather than
    looking them up keeps this module free of any import from the code it
    checks.

    A segment absent from ``allele_sequences``, or with no alignment
    coordinates, is skipped -- that is how a VJ locus's missing D is handled.
    """
    out: list[Disagreement] = []
    identifier = record.get("sequence_id", "<unidentified>")
    inverted = set(inverted)
    sequence = record.get("sequence") or ""
    aligned = record.get("sequence_alignment")
    germline = record.get("germline_alignment")

    if aligned is None or germline is None:
        return ()

    def report(segment, identity, recorded, expected):
        out.append(Disagreement(sequence_id=identifier, segment=segment,
                                identity=identity, recorded=recorded,
                                expected=expected))

    # 1. The two strings are one alignment, so they have one length.
    if len(aligned) != len(germline):
        report("-", "len(sequence_alignment) == len(germline_alignment)",
               len(aligned), len(germline))

    deletions = 0
    for segment in segments:
        reference = allele_sequences.get(segment)
        if reference is None or record.get(f"{segment}_alignment_start") is None:
            continue
        ops = cigar_ops(record.get(f"{segment}_cigar"))
        unknown = sorted({op for _, op in ops} - {"M", "I", "D", "=", "X"})
        if unknown:
            report(segment, "the cigar uses only M, I and D", unknown, "no others")
        read_span = sum(n for n, op in ops if op in ("M", "I", "=", "X"))
        germline_span = sum(n for n, op in ops if op in ("M", "D", "=", "X"))
        deletions += sum(n for n, op in ops if op == "D")

        # 3 and 4. The cigar accounts for exactly the spans the coordinates
        # report. These cross from the cigar to the coordinate fields.
        recorded_read = record[f"{segment}_sequence_end"] - record[f"{segment}_sequence_start"]
        if read_span != recorded_read:
            report(segment, "cigar M+I == read span", read_span, recorded_read)
        recorded_germline = (record[f"{segment}_germline_end"]
                             - record[f"{segment}_germline_start"])
        if germline_span != recorded_germline:
            report(segment, "cigar M+D == germline span", germline_span, recorded_germline)

        # 5. The content check: every M column against the reference.
        span_start = record[f"{segment}_germline_start"]
        span_end = record[f"{segment}_germline_end"]
        is_inverted = segment in inverted
        for column, position, length, op in _walk(record, segment):
            if op not in ("M", "=", "X"):
                continue
            for offset in range(length):
                if column + offset >= len(germline):
                    break
                expected = _expected_germline_base(
                    reference, position + offset, segment_start=span_start,
                    segment_end=span_end, inverted=is_inverted)
                if expected is None:
                    report(segment,
                           "an M column's germline coordinate is inside the allele",
                           position + offset, f"0..{len(reference) - 1}")
                    continue
                found = germline[column + offset].upper()
                if found != expected:
                    report(segment,
                           f"M column {column + offset} carries the reference base "
                           f"at germline {position + offset}", found, expected)

    # 2. The alignment is longer than the read by exactly the deletions.
    if len(aligned) != len(sequence) + deletions:
        report("-", "len(sequence_alignment) == len(sequence) + deletions",
               len(aligned), len(sequence) + deletions)

    # 6. Remove the gaps and the alignment is the read -- reverse complemented
    # when the read is reported on the other strand, because the alignment is
    # written in germline orientation.
    bare = aligned.replace(GAP, "").upper()
    on_other_strand = bool(record.get("rev_comp"))
    expected_read = (reverse_complement(sequence) if on_other_strand
                     else sequence.upper())
    if bare != expected_read.upper():
        report("-",
               "sequence_alignment without gaps is the read"
               + (" reverse complemented" if on_other_strand else ""),
               bare[:40], expected_read.upper()[:40])

    # 7. The D mask replaces the D span's bases and leaves its gaps.
    mask = record.get("germline_alignment_d_mask")
    if mask is not None:
        start = record.get("d_alignment_start")
        end = record.get("d_alignment_end")
        if start is None or end is None:
            start = end = -1
        expected_mask = "".join(
            "N" if (start <= i < end and germline[i] != GAP) else germline[i]
            for i in range(len(germline)))
        if mask != expected_mask:
            report("d", "the D mask is the germline alignment with the D span's "
                        "bases replaced by N", mask[:40], expected_mask[:40])

    return tuple(out)


def columns_compared(record, allele_sequences, *, inverted=(), segments=SEGMENTS
                     ) -> dict[str, int]:
    """Per segment, how many ``M`` columns identity 5 actually compared.

    An empty disagreement tuple means "every column agreed" and also "no column
    was looked at" -- a segment skipped for a missing reference, or an inverted
    D passed over, reports clean either way. Callers assert on this so a
    verdict cannot rest on nothing having been examined.
    """
    counts = {segment: 0 for segment in segments}
    germline = record.get("germline_alignment")
    if germline is None:
        return counts
    for segment in segments:
        reference = allele_sequences.get(segment)
        if reference is None or record.get(f"{segment}_alignment_start") is None:
            continue
        span_start = record[f"{segment}_germline_start"]
        span_end = record[f"{segment}_germline_end"]
        for column, position, length, op in _walk(record, segment):
            if op not in ("M", "=", "X"):
                continue
            for offset in range(length):
                if column + offset >= len(germline):
                    break
                if _expected_germline_base(
                        reference, position + offset, segment_start=span_start,
                        segment_end=span_end, inverted=segment in inverted) is not None:
                    counts[segment] += 1
    return counts


# ── translation ──────────────────────────────────────────────────────────────
#
# The frame is not a property of `sequence_aa` alone, so it cannot be checked
# by inspecting the field. Read from the engine (productivity.rs::coding_span,
# builder.rs:397) and then confirmed against records:
#
#     native_phase = v_allele.anchor % 3          -- the allele's coding phase
#     removed      = v_trim_5 + end_loss_5_length -- bases gone from V's 5' end
#     phase        = (v_sequence_start + native_phase + 3 - removed % 3) % 3
#     subject      = reverse_complement(sequence) if rev_comp else sequence
#     sequence_aa == translate(subject[phase:])
#
# Measured: 150 of 150 records across the plain, indel, mutated and end-loss
# stacks, and 30 of 30 reverse-complemented once the subject is taken in
# germline orientation -- which is where the rest of the projection lives too.
#
# `native_phase` is passed in rather than derived here, because a curated
# reference may carry explicit `coding_phase` evidence that overrides the
# anchor, and this module must not reimplement that precedence.

#: The standard genetic code. Written out rather than generated, because it is
#: the thing that has to be right: a generated table shares any error in the
#: ordering it was generated from, and the engine's table is the one under test.
_CODON_TABLE = {
    "TTT": "F", "TTC": "F", "TTA": "L", "TTG": "L",
    "CTT": "L", "CTC": "L", "CTA": "L", "CTG": "L",
    "ATT": "I", "ATC": "I", "ATA": "I", "ATG": "M",
    "GTT": "V", "GTC": "V", "GTA": "V", "GTG": "V",
    "TCT": "S", "TCC": "S", "TCA": "S", "TCG": "S",
    "CCT": "P", "CCC": "P", "CCA": "P", "CCG": "P",
    "ACT": "T", "ACC": "T", "ACA": "T", "ACG": "T",
    "GCT": "A", "GCC": "A", "GCA": "A", "GCG": "A",
    "TAT": "Y", "TAC": "Y", "TAA": "*", "TAG": "*",
    "CAT": "H", "CAC": "H", "CAA": "Q", "CAG": "Q",
    "AAT": "N", "AAC": "N", "AAA": "K", "AAG": "K",
    "GAT": "D", "GAC": "D", "GAA": "E", "GAG": "E",
    "TGT": "C", "TGC": "C", "TGA": "*", "TGG": "W",
    "CGT": "R", "CGC": "R", "CGA": "R", "CGG": "R",
    "AGT": "S", "AGC": "S", "AGA": "R", "AGG": "R",
    "GGT": "G", "GGC": "G", "GGA": "G", "GGG": "G",
}


def codon_table() -> dict[str, str]:
    """The table, so a test can assert it is complete before relying on it."""
    return dict(_CODON_TABLE)


def translate(sequence: str) -> str:
    """Amino acids for the whole codons of ``sequence``, from position 0.

    A trailing partial codon is dropped, as the engine drops it. An unknown
    codon -- anything with an N or an ambiguity code in it -- becomes ``X``,
    which is distinguishable from a stop.
    """
    upper = sequence.upper()
    return "".join(_CODON_TABLE.get(upper[i:i + 3], "X")
                   for i in range(0, len(upper) - 2, 3))


def translation_phase(record, *, native_phase: int) -> int:
    """The reading frame ``sequence_aa`` is written in."""
    removed = _count(record, "v_trim_5") + _count(record, "end_loss_5_length")
    return (record["v_sequence_start"] + native_phase % 3 + 3 - removed % 3) % 3


def translation_disagreement(record, *, native_phase: int) -> tuple[str, str] | None:
    """``(recorded, expected)`` if ``sequence_aa`` is not the translation.

    ``None`` means it is. Raises when there is nothing to compare -- an empty
    ``sequence_aa`` means the engine found no frame, which is a different
    statement and the caller must decide what it expects.
    """
    recorded = record.get("sequence_aa")
    if not recorded:
        raise ValueError(
            "sequence_aa is empty: the engine found no coding frame for this "
            "record, so there is no translation to check. Gate on it.")
    sequence = record.get("sequence") or ""
    subject = reverse_complement(sequence) if record.get("rev_comp") else sequence
    phase = translation_phase(record, native_phase=native_phase)
    expected = translate(subject[phase:])
    return None if recorded == expected else (recorded, expected)
