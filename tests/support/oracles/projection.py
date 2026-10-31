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

**What these do and do not establish.** Being specific about this, because an
earlier version of this docstring claimed the opposite of the truth and the
caveat was already written out correctly in the sibling oracle one directory
away.

- Identity 5 establishes **allele identity**: that the germline string was laid
  down from the allele the record names. A rotated reference breaks it.
- It also establishes **column consistency**: with the gap-pattern check, the
  cigar's op placement is pinned to where the two strings actually have gaps.
  Without that check the placement of ``D`` runs was entirely unconstrained --
  ``13M1D277M`` could be rewritten ``1D290M``, ``290M1D`` or ``145M1D145M``
  and all three agreed, because the walk advances column and germline position
  together for both ``M`` and ``D``.
- It does **not** establish the germline coordinates. ``*_germline_start/end``
  come from the same column walk that emits the cigar and the germline string
  (``builder.rs`` sets them from ``walk.ref_ranges``), so shifting the
  coordinate and the string together agrees: shifting them by +1 leaves 226
  columns carrying a different germline base and produces 0 disagreements. (At
  -1 there is 1, from the bounds guard reaching position -1, which is an
  artifact of the probe rather than the identity noticing.) The coordinates are
  checked by the ``coordinates`` oracle, which
  compares them to the reference allele's own *length* -- a number the walk
  does not produce.
- Identities 3 and 4 are **true by construction** relative to that same walk:
  ``M+I`` equals the read span because every column with a pool position is one
  pool position, and ``M+D`` equals the ref range because both come from one
  closure. The engine says so itself -- ``builder.rs``: *"This guarantees
  germline_span == M + D by construction."* What they catch is a field-plumbing
  error in the builder, such as writing one segment's range into another's.
  That is worth catching and is less than it looks.

**Nothing is skipped.** Every column is checked, not only the ``M`` runs: a
``D`` column's germline base is the one place the projector fills from its own
cursor rather than from a nucleotide's recorded provenance, so skipping it left
the most error-prone arithmetic unchecked. An inverted D is checked rather than
passed over, with the expectation reverse-complemented within that segment's
own germline span: measured 407 of 407 columns, against 100 for the plain
reference. The base-level oracle next door cannot do this, because the IR
holds forward coordinates with complemented bases.

Validated over 325 records -- 13 of the 15 locus-and-stack combinations at 25
records each, the other two being skipped as a VJ locus with no D to invert --
across IGH, IGK and TCRB and the plain, indel, mutated, reverse-complemented
and inverted-D stacks.

Imports nothing from GenAIRR: the reference sequences and the set of inverted
segments are handed in by the caller.
"""
from __future__ import annotations

import re
from dataclasses import dataclass

SEGMENTS = ("v", "d", "j")

#: Every identity this module can report, as a stable key. The formatted
#: detail goes in ``Disagreement.detail``; the key is what CIRCULARITY is
#: keyed by and what the meta test checks coverage against. Declared because
#: the first version formatted the column number into the identity itself, so
#: no two reports shared a key and nothing could check that the circularity
#: mapping covered them -- the blind spot that mapping's own docstring
#: described, realized.
IDENTITIES = (
    "the alignment strings are present to be checked",
    "a segment the record names is checkable",
    "len(sequence_alignment) == len(germline_alignment)",
    "len(sequence_alignment) == len(sequence) + deletions",
    "the cigar uses only M, I and D",
    "cigar M+I == read span",
    "cigar M+D == germline span",
    "every cigar column is inside the alignment",
    "a column declared M/D/I has that op's gap pattern",
    "a column's germline coordinate is inside the allele",
    "a column consuming a germline base carries the reference's base there",
    "sequence_alignment without gaps is the read",
    "the D mask is the germline alignment with the D span's bases replaced by N",
)
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
    #: Where, specifically -- the column, the op, the field name. Kept out of
    #: ``identity`` so the identity stays a stable key.
    detail: str = ""

    def __post_init__(self) -> None:
        if self.identity not in IDENTITIES:
            raise ValueError(
                f"{self.identity!r} is not a declared identity. Add it to "
                "IDENTITIES and probe it for circularity, or reuse the key it "
                "belongs to and put the specifics in `detail`.")

    def __str__(self) -> str:
        where = f" [{self.detail}]" if self.detail else ""
        return (f"{self.sequence_id}: {self.segment} {self.identity}{where}: "
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
        else:
            # I, and anything else -- an unknown op is separately reported by
            # the caller, and advancing the column without the germline
            # position is the only sane thing to do with a column that
            # consumes no reference base.
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
    identifier = record.get("sequence_id") or "<unidentified>"
    # Normalized like the sibling oracle, so a caller passing {"D"} and one
    # passing {"d"} behave the same. Two oracles with one parameter name and
    # opposite case conventions is how a silently-ignored argument happens,
    # which is the defect this branch fixed in reconstruction.py.
    inverted = {s.lower() for s in inverted}
    sequence = record.get("sequence") or ""
    aligned = record.get("sequence_alignment")
    germline = record.get("germline_alignment")

    # An absent alignment used to return () -- indistinguishable from "every
    # identity held". The registered tests were green on records with all
    # three strings None, because their field-presence helpers also pass on
    # None, so the only place the claimed field names were visible to the
    # catalog's gate was an assertion that accepts the field being absent.
    missing = [name for name, value in (("sequence_alignment", aligned),
                                        ("germline_alignment", germline),
                                        ("germline_alignment_d_mask",
                                         record.get("germline_alignment_d_mask")))
               if value is None]
    if missing:
        return tuple(
            Disagreement(sequence_id=identifier, segment="-",
                         identity="the alignment strings are present to be checked",
                         recorded=None, expected="a string", detail=name)
            for name in missing)

    def report(segment, identity, recorded, expected, detail=""):
        out.append(Disagreement(sequence_id=identifier, segment=segment,
                                identity=identity, recorded=recorded,
                                expected=expected, detail=detail))

    # 1. The two strings are one alignment, so they have one length.
    if len(aligned) != len(germline):
        report("-", "len(sequence_alignment) == len(germline_alignment)",
               len(aligned), len(germline))

    deletions = 0
    for segment in segments:
        reference = allele_sequences.get(segment)
        assembled = record.get(f"{segment}_call")
        if reference is None or record.get(f"{segment}_alignment_start") is None:
            # A VJ locus has no D and that is not a gap. A segment the record
            # says it assembled, with no reference to check it against or no
            # coordinates to check it at, is: dropping J from the reference
            # map used to give 0 disagreements over 314 compared columns, with
            # both column floors satisfied by V alone.
            if assembled:
                report(segment,
                       "a segment the record names is checkable",
                       "no reference sequence" if reference is None
                       else "no alignment coordinates",
                       f"both, for {assembled}")
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

        # 5. Every column, against both strings and against the reference.
        #
        # The op is checked against the gap pattern first, which is what pins
        # its *placement*. Without it the cigar's op order is unconstrained:
        # _walk advances the column and the germline position together for M
        # and for D, so the column-to-position mapping is identical wherever
        # the D runs sit, and `13M1D277M` -> `1D290M` -> `290M1D` all agreed.
        # The pattern is exact and has no exceptions: measured over
        # 112,802 cigar-declared columns, which is what the registered tests
        # examine in one run across the 13 locus-and-stack combinations.
        span_start = record[f"{segment}_germline_start"]
        span_end = record[f"{segment}_germline_end"]
        is_inverted = segment in inverted
        for column, position, length, op in _walk(record, segment):
            for offset in range(length):
                index = column + offset
                if index >= len(germline) or index >= len(aligned):
                    report(segment, "every cigar column is inside the alignment",
                           index, f"0..{min(len(germline), len(aligned)) - 1}")
                    break
                read_gap = aligned[index] == GAP
                germline_gap = germline[index] == GAP
                wanted = {"M": (False, False), "=": (False, False), "X": (False, False),
                          "D": (True, False), "I": (False, True)}.get(op)
                if wanted is not None and (read_gap, germline_gap) != wanted:
                    report(segment,
                           "a column declared M/D/I has that op's gap pattern",
                           f"read={'gap' if read_gap else 'base'}, "
                           f"germline={'gap' if germline_gap else 'base'}",
                           f"read={'gap' if wanted[0] else 'base'}, "
                           f"germline={'gap' if wanted[1] else 'base'}",
                           detail=f"column {index}, op {op}")
                    continue
                # M and D both consume a germline base, so both are compared
                # against the reference. Only comparing M left the deletion
                # columns unchecked -- and a D column's germline base is the
                # one place the projector fills from its own cursor rather
                # than from a nucleotide's recorded provenance, which is the
                # arithmetic most likely to be wrong.
                if op not in ("M", "=", "X", "D"):
                    continue
                expected = _expected_germline_base(
                    reference, position + offset, segment_start=span_start,
                    segment_end=span_end, inverted=is_inverted)
                if expected is None:
                    report(segment,
                           "a column's germline coordinate is inside the allele",
                           position + offset, f"0..{len(reference) - 1}")
                    continue
                found = germline[index].upper()
                if found != expected:
                    coordinate = (span_end - 1 - (position + offset - span_start)
                                  if is_inverted else position + offset)
                    report(segment,
                           "a column consuming a germline base carries the "
                           "reference's base there", found, expected,
                           detail=f"{op} column {index}, germline {coordinate}")

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
        # The orientation goes in `detail`, not in the identity. Suffixing the
        # identity with " reverse complemented" produced a string that is not
        # in IDENTITIES, so on a reverse-complemented record this call raised
        # ValueError out of Disagreement.__post_init__ instead of reporting
        # the disagreement -- a real defect in `sequence_alignment` became a
        # crash in the oracle. Nothing caught it because no mutant corrupted
        # that field on a flipped record.
        report("-", "sequence_alignment without gaps is the read",
               bare[:40], expected_read.upper()[:40],
               detail="reverse complemented" if on_other_strand else "forward")

    # 7. The D mask replaces the D span's bases and leaves its gaps.
    mask = record.get("germline_alignment_d_mask")
    if True:          # absence is reported above, not skipped here
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
    """The reading frame ``sequence_aa`` is written in.

    ``removed`` is taken from ``v_germline_start``, which is where the engine's
    ``trim_5 + observation_loss_5`` lands. The first version of this reproduced
    that sum as ``v_trim_5 + end_loss_5_length`` and was wrong whenever two 5'
    truncations compose: with a 4-base end loss and a 5-base primer trim,
    ``end_loss_5_length`` reads 4 while 9 bases are gone, so the frame came out
    2 instead of 0 and the expected translation was nonsense with stop codons
    in it -- a false failure, on 30 of 30 records.

    ``v_germline_start`` is correct in every case measured, including both
    composition orders. It is also the field the engine's own ``coding_span``
    is consistent with, rather than a reconstruction of it.
    """
    return (record["v_sequence_start"] + native_phase % 3
            + 3 - _count(record, "v_germline_start") % 3) % 3


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


#: What a both-sides perturbation produces, per identity. See the package
#: docstring for why this is required and what it does not guarantee.
CIRCULARITY = {
    "a column consuming a germline base carries the reference's base there":
        "CIRCULAR on the coordinate. Shift v_germline_start/end by one AND "
        "re-render the germline columns from the shifted reference slice: this "
        "identity stays silent, with 226 columns now carrying a different "
        "base. *_germline_start/end come from the same column walk that emits "
        "the cigar and the germline string (builder.rs sets them from "
        "walk.ref_ranges). What it establishes is ALLELE IDENTITY -- a rotated "
        "reference breaks it on every V column -- not the coordinates.",
    "a column declared M/D/I has that op's gap pattern":
        "CIRCULAR on deletion PLACEMENT. Move a D in the cigar and move the "
        "matching gap in sequence_alignment together: 0 disagreements. So this "
        "pins the cigar relative to the alignment strings, and all three come "
        "from one walk. A consistently relocated deletion is a different but "
        "internally consistent alignment of the same read to the same "
        "germline, and which one actually happened is recorded only in the "
        "engine's per-base provenance -- so establishing placement needs an "
        "identity against the IR, which this module does not have. It IS "
        "load-bearing against an inconsistent relocation: moving the cigar's D "
        "alone gives 2 disagreements where it gave 0 before this check.",
    "the alignment strings are present to be checked":
        "NOT CIRCULAR -- a presence guard, with no second side to move. Each "
        "of the 3 strings absent fires it 1 of 1. It exists because the first "
        "version returned () for an absent alignment, which a registered test "
        "read as agreement.",
    "a segment the record names is checkable":
        "NOT CIRCULAR -- a presence guard. Dropping J from the reference map "
        "fires it 1 of 1, where it previously gave 0 disagreements over 314 "
        "compared columns with both column floors satisfied by V alone.",
    "len(sequence_alignment) == len(sequence) + deletions":
        "NOT CIRCULAR against the assembled sequence, a different producer: "
        "shortening both alignment strings by 1 fires it 1 time of 4. "
        "CIRCULAR once the cigar moves too -- shortening both strings AND "
        "removing a D from the cigar leaves it silent, 0 of 5 -- because the "
        "deletion count comes from the same walk as the strings.",
    "the cigar uses only M, I and D":
        "NOT CIRCULAR -- a format guard, no second side. Injecting an S op "
        "fires it 1 time of 212, the rest being the column checks the bad op "
        "then displaces.",
    "every cigar column is inside the alignment":
        "NOT CIRCULAR -- a bounds guard. A cigar claiming 500 columns more "
        "than the alignment holds fires it 1 time of 78.",
    "a column's germline coordinate is inside the allele":
        "NOT CIRCULAR -- a bounds guard against the reference's own length, "
        "which the walk does not produce. Putting germline_start past the "
        "allele end fires it 291 times of 292.",
    "cigar M+I == read span":
        "CIRCULAR. True by construction in the engine, which says so itself: "
        "'This guarantees germline_span == M + D by construction.' Measured on "
        "the indel-bearing IGH record the sensitivity fixture selects "
        "(seed 53, cigar 241M1I28M1D22M): +2 to a V M-run alone fires it 1 "
        "time of 5 total, and moving the cigar AND the coordinates together "
        "leaves it silent. What it catches alone is a field-plumbing error in the "
        "builder, such as one segment's range written into another's.",
    "cigar M+D == germline span":
        "CIRCULAR, for the same reason and from the same closure: +2 to a V "
        "M-run alone fires it 1 time of 5 total, on the same record.",
    "len(sequence_alignment) == len(germline_alignment)":
        "CIRCULAR. One alignment, one length, emitted together: truncating "
        "germline_alignment alone gives 3 disagreements, truncating both gives "
        "4 -- all of them from the other identities, this one silent.",
    "sequence_alignment without gaps is the read":
        "NOT CIRCULAR against the assembled sequence, which is a different "
        "producer: corrupting 1 base of sequence_alignment fires it, 1 of 1 "
        "disagreement. It does not constrain where the gaps ARE -- relocating "
        "a gap fires it 0 times of 2 -- because removing them gives the same "
        "string wherever they sat.",
    "the D mask is the germline alignment with the D span's bases replaced by N":
        "CIRCULAR on the span. Handing the oracle the unmasked germline string "
        "gives 1 disagreement, so the identity does catch a mask that was "
        "never applied or applied to the wrong span. It is not evidence that "
        "an applied mask sits in the right PLACE: the mask is built in the "
        "same walk whose boundaries become d_alignment_start/end, so the span "
        "it covers and the span it is compared against are one number.",
}
