"""Check the recorded coordinates against arithmetic the engine must obey.

The reconstruction oracle next door compares each recorded germline *base* to
the reference base at the coordinate the record claims. That catches an
assembly copied from the wrong allele. It cannot catch a wrong coordinate,
because production lays the byte and the coordinate down from the same two
numbers -- ``Nucleotide::germline(slice[i], slice_start + i)`` -- so a trim
that slices one base too many moves both together and the comparison still
agrees.

This module never compares a base. It checks that four numbers produced by
different parts of the system agree:

- the trim amounts the trimming pass sampled (``*_trim_5``, ``*_trim_3``),
- the end-loss lengths the assay recorded (``end_loss_*_length``),
- the coordinates the record builder exported from the nucleotide pool
  (``*_germline_*``, ``*_sequence_*``, ``*_alignment_*``),
- the allele's own length, read by the caller from the reference table.

The identities, per segment present in the record:

1. ``germline_start == trim_5 + lead``
2. ``germline_end   == allele_length - trim_3 - tail``
3. ``germline_end - germline_start == sequence_end - sequence_start``
4. ``(alignment_start, alignment_end) == (sequence_start, sequence_end)``

and three more over the record as a whole, which pin absolute position rather
than differences -- see ``_chain_disagreements``:

5. the first segment starts at read position 0
6. the gap between consecutive segments is the untemplated length recorded there
7. the read ends where the last segment ends

``lead`` is ``end_loss_5_length`` for V and zero elsewhere; ``tail`` is
``end_loss_3_length`` for J and zero elsewhere, because read truncation can
only shorten the rearrangement at its two outer ends.

Identities 1-3 are exact when the segment carries no indel, and may otherwise
be off by at most ``n_<segment>_indels``: an insertion or deletion is the only
recorded event that displaces a coordinate relative to the germline, and it
cannot displace it by more than once per edit. Identity 4 is exact when the
whole record carries no indel -- an indel anywhere shifts the alignment frame
of every segment downstream of it -- and otherwise requires the alignment
coordinates to sit at or after the read coordinates and span at least as far.

Reverse complement is not an exception to any of them: the coordinate fields
are reported in germline orientation, which these identities assert by holding
unchanged when the read is flipped.

Measured while writing this, over records built from four loci (IGH, IGK, IGL,
TCRB) and the pass stacks the committed tests exercise -- 5,920 segment
observations for identities 1-4 and 1,380 records for 5-7. The committed tests
cover three of those loci (``LOCI`` in the module that registers them) plus
drawn seeds under Hypothesis; IGL was in the exploratory measurement only. Treat the
numbers as the author's evidence for the identities, not as a description of
what runs on every commit.
"""
from __future__ import annotations

from dataclasses import dataclass

SEGMENTS = ("v", "d", "j")


@dataclass(frozen=True)
class Disagreement:
    """One identity that does not hold, and by how much."""

    sequence_id: str
    segment: str
    identity: str
    recorded: object
    expected: object
    tolerated: int = 0

    def __str__(self) -> str:
        return (f"{self.sequence_id}: {self.segment.upper()} {self.identity}: "
                f"recorded {self.recorded}, expected {self.expected}"
                + (f" (±{self.tolerated})" if self.tolerated else ""))


def _count(record, field: str) -> int:
    """A record's integer field, treating an absent or null count as zero."""
    value = record.get(field)
    return 0 if value is None else int(value)


def disagreements(record, allele_lengths, *, segments=SEGMENTS) -> tuple[Disagreement, ...]:
    """Identities in this record's coordinates that do not hold.

    ``allele_lengths`` maps a segment letter to the length of the allele this
    record used, read from the reference table by the caller. Passing it in
    rather than looking it up keeps this module free of any import from the
    code it checks.

    A segment absent from ``allele_lengths`` is skipped, which is how a VJ
    locus's missing D is handled.
    """
    out: list[Disagreement] = []
    identifier = record.get("sequence_id", "<unidentified>")
    record_indels = _count(record, "n_indels")
    lead_loss = _count(record, "end_loss_5_length")
    tail_loss = _count(record, "end_loss_3_length")

    for segment in segments:
        length = allele_lengths.get(segment)
        if length is None or record.get(f"{segment}_sequence_start") is None:
            continue
        germline_start = record[f"{segment}_germline_start"]
        germline_end = record[f"{segment}_germline_end"]
        sequence_start = record[f"{segment}_sequence_start"]
        sequence_end = record[f"{segment}_sequence_end"]
        alignment_start = record[f"{segment}_alignment_start"]
        alignment_end = record[f"{segment}_alignment_end"]
        tolerated = _count(record, f"n_{segment}_indels")
        lead = lead_loss if segment == "v" else 0
        tail = tail_loss if segment == "j" else 0

        # The indel tolerance is not uniform, and granting it uniformly was
        # dead slack. Measured over 720 records across three loci and three
        # indel-bearing stacks: the germline bounds deviate only at an end
        # where read truncation also applies -- an indel inside the lost
        # window moves the first or last retained germline base -- and are
        # exact everywhere else however many indels the segment carries. So a
        # one-base error in d_germline_end is a finding even on a record with
        # three indels in D, which under a uniform tolerance it was not.
        for identity, recorded, expected, slack in (
            ("germline_start == trim_5 + end_loss_5_length",
             germline_start, _count(record, f"{segment}_trim_5") + lead,
             tolerated if lead else 0),
            ("germline_end == allele_length - trim_3 - end_loss_3_length",
             germline_end, length - _count(record, f"{segment}_trim_3") - tail,
             tolerated if tail else 0),
            # The read span is the one an indel genuinely displaces: an
            # insertion or deletion inside the segment changes how many read
            # bases its germline span occupies.
            ("germline span == read span",
             germline_end - germline_start, sequence_end - sequence_start,
             tolerated),
        ):
            if abs(recorded - expected) > slack:
                out.append(Disagreement(sequence_id=identifier, segment=segment,
                                        identity=identity, recorded=recorded,
                                        expected=expected, tolerated=slack))

        if record_indels == 0:
            if (alignment_start, alignment_end) != (sequence_start, sequence_end):
                out.append(Disagreement(
                    sequence_id=identifier, segment=segment,
                    identity="alignment coordinates == read coordinates",
                    recorded=(alignment_start, alignment_end),
                    expected=(sequence_start, sequence_end)))
        elif not (alignment_start >= sequence_start and alignment_end >= sequence_end
                  and alignment_end - alignment_start >= sequence_end - sequence_start):
            out.append(Disagreement(
                sequence_id=identifier, segment=segment,
                identity="alignment coordinates are at or after the read coordinates",
                recorded=(alignment_start, alignment_end),
                expected=f"at or after ({sequence_start}, {sequence_end})"))

    out.extend(_chain_disagreements(record, identifier, allele_lengths, segments))
    return tuple(out)


def _chain_disagreements(record, identifier, allele_lengths, segments):
    """Identities that pin each segment's position in the read absolutely.

    Everything above constrains differences -- spans against spans, one
    coordinate pair against another. That leaves a whole class of defect
    invisible: shift a segment's read and alignment coordinates together by the
    same amount and every difference is preserved, so the oracle sees nothing
    while the coordinates an aligner is scored on are all wrong.

    These three tie the absolute positions to the read itself:

    5. the first segment starts at 0 -- the read begins where the
       rearrangement begins,
    6. the gap between consecutive segments is exactly the untemplated length
       recorded between them, so no read position is unaccounted for,
    7. the read ends where the last segment ends.

    Taken together they fix every segment's absolute position: there is no
    offset you can add to one segment's coordinates that satisfies all three.
    No indel tolerance applies -- measured over 1,380 records across four loci
    and six stacks including indels, end loss and reverse complement, these
    hold exactly.
    """
    present = [s for s in segments
               if allele_lengths.get(s) is not None
               and record.get(f"{s}_sequence_start") is not None]
    if not present:
        return []

    out = []
    first, last = present[0], present[-1]
    if record[f"{first}_sequence_start"] != 0:
        out.append(Disagreement(
            sequence_id=identifier, segment=first,
            identity="the first segment starts at read position 0",
            recorded=record[f"{first}_sequence_start"], expected=0))

    # Three segments means two junctions, np1 then np2. Two segments means one
    # junction, which carries whatever both np fields hold between them.
    np_lengths = [_count(record, "np1_length"), _count(record, "np2_length")]
    gaps = np_lengths if len(present) == 3 else [sum(np_lengths)]
    for index, segment in enumerate(present[1:]):
        previous = present[index]
        gap = record[f"{segment}_sequence_start"] - record[f"{previous}_sequence_end"]
        expected = gaps[index] if index < len(gaps) else 0
        if gap != expected:
            out.append(Disagreement(
                sequence_id=identifier, segment=f"{previous}->{segment}",
                identity="the gap between segments is the untemplated length recorded there",
                recorded=gap, expected=expected))

    trailing = _count(record, "sequence_length") - record[f"{last}_sequence_end"]
    if trailing != 0:
        out.append(Disagreement(
            sequence_id=identifier, segment=last,
            identity="the read ends where the last segment ends",
            recorded=record[f"{last}_sequence_end"],
            expected=_count(record, "sequence_length")))
    return out


#: What a both-sides perturbation produces, per identity. See the package
#: docstring for why this is required and what it does not guarantee.
#:
#: The result here is the interesting one: EVERY identity is individually
#: circular -- perturb both of its own sides and it stays silent -- and the
#: SET is not, because perturbing one pair breaks a different identity. The
#: strength is in the over-determination, not in any single comparison, which
#: is the opposite of how the first draft of this module described itself.
CIRCULARITY = {
    "germline_start == trim_5 + end_loss_5_length":
        "CIRCULAR alone: +3 to both germline_start and trim_5 leaves this "
        "identity silent. Caught by the set -- 1 other identity fires.",
    "germline_end == allele_length - trim_3 - end_loss_3_length":
        "CIRCULAR alone: -3 to germline_end with +3 to trim_3 leaves this "
        "identity silent. Caught by the set -- 1 other identity fires. This "
        "one does additionally cross the reference's own LENGTH: stretching "
        "the allele by 3 fires it, which no other identity here does.",
    "germline span == read span":
        "CIRCULAR alone: -2 to both spans leaves this identity silent. Caught "
        "by the set -- 3 other identities fire.",
    "alignment coordinates == read coordinates":
        "CIRCULAR alone: +2 or +5 to all four of a segment's alignment and "
        "read coordinates leaves this identity silent, because both pairs "
        "come from the walk. Caught by the set -- 2 other identities fire, "
        "the leading anchor and the inter-segment gap.",
    "the first segment starts at read position 0":
        "NOT CIRCULAR -- an absolute anchor, with no second side to move. "
        "Displacing the whole rearrangement by 5 fires it 1 time of 2 total "
        "disagreements, and every difference-based identity above preserves "
        "that displacement, which is why the anchor exists.",
    "the gap between segments is the untemplated length recorded there":
        "NOT CIRCULAR against np1/np2, which are a different producer: "
        "np1_length - 1 fires it 1 time of 1 total disagreement. Circular if "
        "the coordinates and the np length move together.",
    "the read ends where the last segment ends":
        "NOT CIRCULAR -- the second absolute anchor. The same +5 displacement "
        "fires it 1 time of 2, against sequence_length, which the segment "
        "coordinates do not produce.",
}
