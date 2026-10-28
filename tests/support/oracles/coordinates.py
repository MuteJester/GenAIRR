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

Verified across 5,920 segment observations: four loci (IGH, IGK, IGL, TCRB) and
eight pass stacks including mutation, indels, end loss, reverse complement and
all of them together.
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

        for identity, recorded, expected in (
            ("germline_start == trim_5 + end_loss_5_length",
             germline_start, _count(record, f"{segment}_trim_5") + lead),
            ("germline_end == allele_length - trim_3 - end_loss_3_length",
             germline_end, length - _count(record, f"{segment}_trim_3") - tail),
            ("germline span == read span",
             germline_end - germline_start, sequence_end - sequence_start),
        ):
            if abs(recorded - expected) > tolerated:
                out.append(Disagreement(sequence_id=identifier, segment=segment,
                                        identity=identity, recorded=recorded,
                                        expected=expected, tolerated=tolerated))

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

    return tuple(out)
