"""Check the reported junction against the anchors the reference declares.

The junction is what a CDR3 benchmark reads. Clustering thresholds, repertoire
overlap statistics and every CDR3-length distribution are computed from
``junction`` / ``junction_aa`` / ``junction_length``, and the two bound fields
``junction_start`` / ``junction_end`` are what a tool is told the window is.
Nothing registered checked any of them.

Five identities, each measured before being written down. The anchor
coordinate and the anchor residue come from the reference configuration rather
than from the run, which is what makes this more than a comparison between two
exported views of one number:

1. ``junction_start`` is the read index of the first surviving base of the V
   anchor codon.
2. ``junction_end`` is one past the read index of the last surviving base of
   the J anchor codon.
3. ``junction`` is ``subject[junction_start:junction_end]``, the subject being
   the read reverse-complemented when ``rev_comp`` is set, because the
   junction is reported in germline orientation.
4. ``junction_length == junction_end - junction_start``.
5. ``junction_aa`` is the oracle's own translation of ``junction``.

**What the reference contributes, and the one thing it cannot see.** The V
anchor codon is a cysteine and the J anchor codon a tryptophan or a
phenylalanine in every bundled BCR allele -- measured 198/198, 168/168 and
181/181 V alleles across human IGH, IGK and IGL, and 7/7 IGH J codons TGG with
8/8 IGK and 10/10 IGL J codons TTT/TTC. So a junction bound displaced off the
anchor is caught: of 198 IGH V alleles only 2 carry another in-frame cysteine
within nine nucleotides of the anchor, and the figures for IGK and IGL are
1/168 and 11/181.

A *large* displacement is a different matter. 197 of 198 IGH V alleles carry
another in-frame cysteine somewhere -- FR1 Cys23, around seventy nucleotides
upstream -- so a bound moved onto that codon satisfies the residue check and
this oracle cannot distinguish it. The claim these identities support is
therefore "the bound names an anchor codon the reference declares", not "the
bound names the CDR3 anchor specifically".
"""
from __future__ import annotations

from dataclasses import dataclass

from .projection import translate

#: Every identity this module can report. A `Disagreement` may not name one
#: that is not here, because an identity nobody wrote down is an identity
#: nobody probed for circularity.
IDENTITIES = (
    "junction_start == first surviving base of the V anchor codon",
    "junction_end == one past the last surviving base of the J anchor codon",
    "junction == the subject sliced at the recorded bounds",
    "junction_length == junction_end - junction_start",
    "junction_aa == translation of junction",
    "the reference anchor codon carries the conserved residue",
    "a record with no junction reports it empty and says the frame is unresolved",
)

CYSTEINE = frozenset({"TGT", "TGC"})
#: The J anchor residue is W on heavy chains and F on light chains; the
#: junction is defined by whichever the locus carries.
TRYPTOPHAN_OR_PHENYLALANINE = frozenset({"TGG", "TTT", "TTC"})


@dataclass(frozen=True)
class Disagreement:
    """One identity that does not hold, and where."""

    sequence_id: str
    identity: str
    recorded: object
    expected: object
    #: The specifics -- which bound, which allele -- kept out of ``identity``
    #: so the identity stays a stable key.
    detail: str = ""

    def __post_init__(self) -> None:
        if self.identity not in IDENTITIES:
            raise ValueError(
                f"{self.identity!r} is not a declared identity. Add it to "
                "IDENTITIES and probe it for circularity, or reuse the key it "
                "belongs to and put the specifics in `detail`.")


def anchor_codon_positions(simulation, segment_of, segment: str,
                           anchor: int) -> tuple[int, ...]:
    """Read indices of the anchor codon's bases that survived, in order.

    Trimming can remove part of a conserved codon, so this is the surviving
    subset rather than always three positions. ``segment_of`` maps a read
    index to its segment, so a germline coordinate is only matched inside the
    segment that declares the anchor -- without that, a D or J base sharing
    the coordinate would be picked up.
    """
    wanted = {anchor, anchor + 1, anchor + 2}
    found = []
    for index in range(len(simulation.bases())):
        if segment_of.get(index) != segment:
            continue
        position = simulation.germline_position(index)
        if position is not None and position in wanted:
            found.append(index)
    return tuple(found)


def reference_anchor_codon(allele_sequence: bytes, anchor: int) -> str | None:
    """The reference's anchor codon, or ``None`` if the allele is too short."""
    if anchor is None or anchor + 3 > len(allele_sequence):
        return None
    return allele_sequence[anchor:anchor + 3].decode().upper()


def subject(record) -> str:
    """The strand the junction is reported against.

    ``rev_comp`` means the record's coordinates are with respect to the
    reverse complement of ``sequence``, so the junction slice is taken there.
    Probe P4 is what establishes this rather than assuming it.
    """
    from .reconstruction import reverse_complement
    if record.get("rev_comp"):
        return reverse_complement(record["sequence"])
    return record["sequence"]


def disagreements(record, simulation, segment_of, *,
                  v_anchor: int, j_anchor: int,
                  v_sequence: bytes, j_sequence: bytes) -> tuple[Disagreement, ...]:
    """Every identity above that does not hold for this record."""
    out: list[Disagreement] = []
    sequence_id = record.get("sequence_id", "?")
    start = record["junction_start"]
    end = record["junction_end"]

    if start is None or end is None:
        # Trimming can remove a conserved codon outright, and then there is no
        # window to report. The bound identities do not apply; what does is
        # that the whole group goes null together and the frame flag says
        # UNRESOLVED rather than out-of-frame. Measured on 60 of 60 records
        # under trim(v_3=12..25): junction and junction_aa empty,
        # junction_length and junction_end null, vj_in_frame None, productive
        # False, and the reasons naming anchor_removed and
        # vj_frame_unresolved. A consumer reading vj_in_frame as a boolean
        # gets False for "we could not tell", which is why this is an
        # identity and not a skip.
        expected = ("", "", None, None)
        recorded = (record["junction"], record["junction_aa"],
                    record["junction_length"], record["vj_in_frame"])
        if recorded != expected:
            out.append(Disagreement(
                sequence_id,
                "a record with no junction reports it empty and says the "
                "frame is unresolved",
                recorded, expected,
                detail="(junction, junction_aa, junction_length, vj_in_frame)"))
        out.extend(_reference_residue_disagreements(
            sequence_id, v_sequence, v_anchor, j_sequence, j_anchor))
        return tuple(out)

    v_codon = anchor_codon_positions(simulation, segment_of, "V", v_anchor)
    if v_codon and start != v_codon[0]:
        out.append(Disagreement(
            sequence_id,
            "junction_start == first surviving base of the V anchor codon",
            start, v_codon[0],
            detail=f"surviving V anchor bases at {list(v_codon)}"))

    j_codon = anchor_codon_positions(simulation, segment_of, "J", j_anchor)
    if j_codon and end != j_codon[-1] + 1:
        out.append(Disagreement(
            sequence_id,
            "junction_end == one past the last surviving base of the J anchor codon",
            end, j_codon[-1] + 1,
            detail=f"surviving J anchor bases at {list(j_codon)}"))

    sliced = subject(record)[start:end]
    if record["junction"].upper() != sliced.upper():
        out.append(Disagreement(
            sequence_id, "junction == the subject sliced at the recorded bounds",
            record["junction"], sliced,
            detail="rev_comp" if record.get("rev_comp") else "forward"))

    if record["junction_length"] != end - start:
        out.append(Disagreement(
            sequence_id, "junction_length == junction_end - junction_start",
            record["junction_length"], end - start))

    expected_aa = translate(record["junction"])
    if record["junction_aa"].upper() != expected_aa.upper():
        out.append(Disagreement(
            sequence_id, "junction_aa == translation of junction",
            record["junction_aa"], expected_aa))

    out.extend(_reference_residue_disagreements(
        sequence_id, v_sequence, v_anchor, j_sequence, j_anchor))
    return tuple(out)


def _reference_residue_disagreements(sequence_id, v_sequence, v_anchor,
                                     j_sequence, j_anchor):
    """The reference-only identity, which applies whether or not the record
    kept its junction."""
    out = []
    for label, sequence, anchor, allowed in (
            ("V", v_sequence, v_anchor, CYSTEINE),
            ("J", j_sequence, j_anchor, TRYPTOPHAN_OR_PHENYLALANINE)):
        codon = reference_anchor_codon(sequence, anchor)
        if codon is not None and codon not in allowed:
            out.append(Disagreement(
                sequence_id,
                "the reference anchor codon carries the conserved residue",
                codon, sorted(allowed),
                detail=f"{label} anchor at {anchor}"))
    return out


#: What was perturbed, and what this oracle said. Measured on 60 unperturbed
#: ``human_igh`` records, which produce zero disagreements before any probe.
CIRCULARITY = {
    "junction_start == first surviving base of the V anchor codon":
        "NOT CIRCULAR against the reference residue, CIRCULAR against the "
        "slice. Shifting the V anchor by +3 and junction_start by +3 together "
        "leaves this identity silent -- it is an independent implementation "
        "over the same two inputs the engine reads, so a displaced anchor "
        "moves both sides. What fires instead is the residue check, 60 of 60, "
        "because the codon three bases on is not a cysteine. Of 198 IGH V "
        "alleles only 2 carry another in-frame cysteine within nine "
        "nucleotides of the anchor (1 of 168 on IGK, 11 of 181 on IGL), so a "
        "small displacement is caught for all but those. A displacement onto "
        "FR1 Cys23, around seventy nucleotides upstream, is NOT caught: 197 "
        "of 198 IGH V alleles carry such a codon. The claim is therefore that "
        "the bound names an anchor codon the reference declares, not that it "
        "names the CDR3 anchor specifically.",
    "junction_end == one past the last surviving base of the J anchor codon":
        "NOT CIRCULAR against a consistent re-slice. Moving both bounds by +3 "
        "and re-deriving junction, junction_length and junction_aa from the "
        "sequence silences identities 3, 4 and 5 -- they share one producer -- "
        "and this identity fires 60 of 60, as does the V bound above. These "
        "two are the absolute anchors of the set, the same shape as anchors "
        "5-7 in the coordinates oracle.",
    "junction == the subject sliced at the recorded bounds":
        "CIRCULAR once the bounds move with it, NOT CIRCULAR against the "
        "sequence. A consistent re-slice at shifted bounds leaves it silent "
        "(probe above); changing one base of `junction` and re-translating "
        "`junction_aa` fires it on 56 of 56 records where the substitution "
        "landed. Four of 60 probe mutations changed nothing because the "
        "junction base was lowercase and the replacement uppercase -- the "
        "same case trap the mutation counters carried, here in the probe "
        "rather than in the engine.\n\n"
        "        Orientation is part of this identity and was probed "
        "separately: flipping rev_comp and reverse-complementing `sequence` "
        "together produces 0 disagreements over 60 records. That is what "
        "establishes that the junction is reported in germline orientation "
        "rather than assuming it, and why `subject()` exists.",
    "junction_length == junction_end - junction_start":
        "CIRCULAR against a consistent change, NOT CIRCULAR against an "
        "inconsistent one. Shifting junction_start by +3 without re-deriving "
        "the length fires it 60 of 60; shifting both bounds and re-deriving "
        "the length from them leaves it silent, because it reads those two "
        "bounds and nothing else. Kept for the first case: it catches a "
        "length that disagrees with its own bounds, which is what the field "
        "would be if it were derived from the junction string instead.",
    "junction_aa == translation of junction":
        "CIRCULAR against the junction string. Changing one base and "
        "re-translating leaves this silent while identity 3 fires. It "
        "establishes that the translation is the oracle's own table applied "
        "to the recorded nucleotides, and nothing about whether those "
        "nucleotides are the right ones.",
    "a record with no junction reports it empty and says the frame is unresolved":
        "NOT CIRCULAR -- it compares four fields that no single producer "
        "writes together against the null contract, so there is no second "
        "side to move. Probed on 60 of 60 records under trim(v_3=12..25), "
        "every one of which has its V anchor codon trimmed away: turning "
        "vj_in_frame from None into False fires it 60 of 60, populating "
        "junction while the bounds stay null fires 60 of 60, and turning "
        "junction_length from None into 0 fires 60 of 60. Unperturbed, 0 of "
        "60. The first of those three is the mis-read this identity exists "
        "for: None means the frame could not be resolved, and a consumer "
        "reading the field as a boolean gets 'out of frame' for 'we could "
        "not tell'. This path also found a fault in this oracle -- before it "
        "existed the module raised TypeError on such a record, and the "
        "fixtures of the test that shipped it never produced one.",
    "the reference anchor codon carries the conserved residue":
        "NOT CIRCULAR -- it reads the reference only, with no side from the "
        "run to move. This is the identity that makes the two bound checks "
        "mean something, and the one that would catch a mis-curated anchor. "
        "Measured across every bundled BCR locus: 198/198, 168/168 and "
        "181/181 V anchor codons are TGT or TGC on human IGH, IGK and IGL; "
        "7/7 IGH J codons are TGG and 8/8 IGK with 10/10 IGL are TTT or TTC.",
}
