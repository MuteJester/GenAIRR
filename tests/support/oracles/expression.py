"""Whether the rearrangement could express, re-derived rather than re-read.

The three productivity flags are what a benchmark filters on. The registered
test that covered them asserted an implication between them -- that a
`productive` record is in frame and free of stops -- which is a comparison
between three exported views of one assessment and holds whatever the
assessment was. This module derives each of them instead, from the assembled
molecule and this package's own codon table.

Three identities, each measured before being written down:

E1. ``vj_in_frame == (junction_length % 3 == 0)``. Measured 299/299 on plain
    records, 299/299 under s5f at ten mutations, 100/100 under
    ``productive_only``. No phase term is needed; one was in the design note
    and the measurement did not want it.
E2. ``stop_codon == ('*' in translation of the coding span)``, the span being
    ``v_sequence_start`` rounded up to the record's frame through
    ``j_sequence_end``. Measured 199/199 plain, 199/199 mutated, 200/200 on
    IGK.
E3. ``productive == vj_in_frame and not stop_codon and the read's V anchor
    codon translates to C and its J anchor codon to W or F``. Measured
    299/299 on all three fixtures above, the mutated one included, which is
    where the anchors actually get hit.

**What E2 settles.** ``stop_codon`` is NOT junction-local. Read as "a stop in
``junction_aa``" it disagrees on 2 of 299 plain records and 46 of 299 mutated
ones -- records whose junction translates cleanly while the flag is set and
the engine's own reason is ``coding_stop``. The flag follows the whole coding
sequence, which is the biology and what the reason name says.
"""
from __future__ import annotations

from dataclasses import dataclass

from .projection import translate, translation_phase
from .reconstruction import reverse_complement

IDENTITIES = (
    "vj_in_frame == the junction length being a whole number of codons",
    "stop_codon == a stop in the translated coding span",
    "productive == in frame, no stop, and both anchors intact",
)


@dataclass(frozen=True)
class Disagreement:
    """One identity that does not hold, and where."""

    sequence_id: str
    identity: str
    recorded: object
    expected: object
    detail: str = ""

    def __post_init__(self) -> None:
        if self.identity not in IDENTITIES:
            raise ValueError(
                f"{self.identity!r} is not a declared identity. Add it to "
                "IDENTITIES and probe it for circularity, or reuse the key it "
                "belongs to and put the specifics in `detail`.")


def subject(record) -> str:
    """The strand the coordinates are with respect to."""
    if record.get("rev_comp"):
        return reverse_complement(record["sequence"])
    return record["sequence"]


def coding_translation(record, *, native_phase: int = 0) -> str:
    """The coding span translated with this package's own codon table.

    The span starts at ``v_sequence_start`` advanced to the next position in
    the record's reading frame and runs to ``j_sequence_end``. The frame comes
    from the translation oracle so the precedence rules over curated coding
    phase are not reimplemented here.
    """
    phase = translation_phase(record, native_phase=native_phase)
    start = record["v_sequence_start"] + ((phase - record["v_sequence_start"]) % 3)
    return translate(subject(record)[start:record["j_sequence_end"]])


def anchors_intact(record) -> tuple[bool, bool]:
    """``(V anchor reads as C, J anchor reads as W or F)`` in the molecule.

    Read off the junction window rather than the reference, because the
    question is whether the anchor SURVIVED into this record -- a mutation can
    destroy a cysteine the reference declares.
    """
    body = subject(record)
    start, end = record["junction_start"], record["junction_end"]
    return (translate(body[start:start + 3]) == "C",
            translate(body[end - 3:end]) in ("W", "F"))


def disagreements(record, *, native_phase: int = 0) -> tuple[Disagreement, ...]:
    """Every identity above that does not hold for this record.

    A record whose anchor was trimmed away has no junction to reason from and
    no resolved frame; the junction oracle owns that contract, and this
    returns nothing rather than inventing an expectation for it.
    """
    if record.get("junction_length") is None:
        return ()
    out: list[Disagreement] = []
    sequence_id = record.get("sequence_id", "?")

    in_frame = record["junction_length"] % 3 == 0
    if bool(record["vj_in_frame"]) != in_frame:
        out.append(Disagreement(
            sequence_id,
            "vj_in_frame == the junction length being a whole number of codons",
            record["vj_in_frame"], in_frame,
            detail=f"junction_length={record['junction_length']}"))

    has_stop = "*" in coding_translation(record, native_phase=native_phase)
    if bool(record["stop_codon"]) != has_stop:
        out.append(Disagreement(
            sequence_id, "stop_codon == a stop in the translated coding span",
            record["stop_codon"], has_stop,
            detail="coding span, not the junction"))

    v_intact, j_intact = anchors_intact(record)
    expected = in_frame and not has_stop and v_intact and j_intact
    if bool(record["productive"]) != expected:
        out.append(Disagreement(
            sequence_id,
            "productive == in frame, no stop, and both anchors intact",
            record["productive"], expected,
            detail=f"in_frame={in_frame} stop={has_stop} "
                   f"v_anchor_C={v_intact} j_anchor_WF={j_intact}"))
    return tuple(out)


#: What was perturbed, and what this module said. Measured on 199 mutated
#: ``human_igh`` records (s5f, ten mutations), which produce zero
#: disagreements unperturbed.
CIRCULARITY = {
    "vj_in_frame == the junction length being a whole number of codons":
        "NOT CIRCULAR against the flag, CIRCULAR against the junction "
        "fields. Flipping vj_in_frame alone fires it 199 of 199. Moving "
        "junction_length and vj_in_frame together leaves it silent -- both "
        "are fields of one record written by one producer -- and what fires "
        "instead is the productive identity, 54 of 199, because the "
        "conjunction no longer agrees. The anchor underneath this is in the "
        "junction oracle: its two bound identities are what establish the "
        "length this one reads, and they fire 60 of 60 under a consistent "
        "re-slice.",
    "stop_codon == a stop in the translated coding span":
        "NOT CIRCULAR -- two producers, and the probe that shows it is the "
        "one that stays silent. Flipping stop_codon alone fires 199 of 199. "
        "Rewriting a codon of the molecule to TAA and setting the flag to "
        "match leaves this identity silent, which is the identity HOLDING "
        "rather than a blind spot: the flag followed the sequence, and the "
        "sequence is written by the assembly passes while the flag is "
        "written by the productivity assessment. The productive identity "
        "fires 27 of 199 on that same probe, because a record with a stop "
        "cannot express.\n\n"
        "        This identity also settled what the flag means. Read as 'a "
        "stop in junction_aa', the documented definition, it disagrees on 2 "
        "of 299 plain records and 46 of 299 mutated ones -- junctions that "
        "translate cleanly while the flag is set and the engine's own reason "
        "is `coding_stop`. The flag is whole-coding-span.",
    "productive == in frame, no stop, and both anchors intact":
        "NOT CIRCULAR against the flag, and the identity that catches what "
        "the other two miss. Flipping productive alone fires 199 of 199, and "
        "it is the one that fires on both joint probes above, 54 and 27 of "
        "199. The anchor residues are read from the molecule rather than the "
        "reference on purpose: the question is whether the anchor survived "
        "into this record, and a mutation can destroy a cysteine the "
        "reference declares. Measured 299/299 on plain records, 299/299 "
        "under s5f at ten mutations and 100/100 under productive_only -- the "
        "mutated fixture being the one where anchors actually get hit.",
}
