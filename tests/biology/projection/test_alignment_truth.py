"""The projection a downstream tool is scored against must match the reference.

GenAIRR is used to benchmark annotators, and the fields here are what those
annotators are compared to. That makes a wrong cigar the most expensive silent
error in the record: every tool measured against it is wrong in the same
direction, which reads as a property of the tools rather than of the oracle.

The existing cigar tests are structurally sound -- op types, run-length
collapse, canonical ops, op count against alignment columns -- and none of them
reads the reference. Their strongest member says so itself: it re-derives the
cigar by calling the production helpers and notes that "equivalence is
guaranteed by construction". So a cigar whose op lengths are right and whose
placement is wrong passes all of them.

What this module adds is the crossing none of them makes: from the projector's
output to the reference table, at scale, with indels, reverse complement and
inverted D present.
"""
from __future__ import annotations

import pytest

from tests.biology.catalog import verifies
from tests.support.builders import experiment as experiment_builder
from tests.support.oracles import projection

pytestmark = pytest.mark.behavior

LOCI = ("human_igh", "human_igk", "human_tcrb")

#: Each stack reaches a shape the projector handles differently: gaps in both
#: alignment strings, the read reported on the opposite strand, and a D whose
#: germline is carried reverse-complemented.
STACKS = {
    "plain": lambda exp: exp,
    "indels": lambda exp: exp.polymerase_indels(count=3),
    "mutated": lambda exp: exp.mutate(model="s5f", count=6),
    "reverse-complemented": lambda exp: exp.random_strand_orientation(prob=1.0),
    "inverted-d": lambda exp: exp.invert_d(prob=1.0),
}

#: Stacks a TCR locus cannot run: SHM is a B-cell phenomenon and the engine
#: refuses it rather than simulating biology it does not claim.
B_CELL_ONLY = frozenset({"mutated"})

#: Stacks that need a D segment to mean anything. A VJ locus has none, so
#: inverting is a no-op there and the stack would pass by having nothing to do.
NEEDS_D = frozenset({"inverted-d"})

#: Loci with a D segment.
VDJ_LOCI = frozenset({"human_igh", "human_tcrb"})


def _observed(exp, *, n: int, seed: int):
    """``(record, allele_sequences, inverted)`` per outcome.

    The reference sequences are read here and handed to the oracle, which
    imports nothing from GenAIRR. ``inverted`` names the segments assembled in
    reverse orientation, which the oracle needs because the germline alignment
    carries the reverse complement for those.
    """
    refdata = exp.refdata
    getters = {"v": refdata.v_allele, "d": refdata.d_allele, "j": refdata.j_allele}
    for outcome in exp.run(n=n, seed=seed):
        simulation = outcome.final_simulation()
        identifiers = {"v": simulation.v_allele_id(), "d": simulation.d_allele_id(),
                       "j": simulation.j_allele_id()}
        sequences = {}
        for segment, identifier in identifiers.items():
            if identifier is None:
                continue
            allele = getters[segment](identifier)
            raw = allele.seq() if callable(allele.seq) else allele.seq
            sequences[segment] = raw.decode().upper()
        record = outcome.airr_record(refdata)
        # Read from the record, not the simulation: a Simulation exposes no
        # orientation at all. The previous oracle's `skip_inverted` read
        # `getattr(simulation, "d_inverted", None)`, always found None, and so
        # skipped nothing while claiming inverted D was skipped by default.
        inverted = {"d"} if record.get("d_inverted") else set()
        yield record, sequences, inverted


def _cigar_fields_present(record):
    """The three cigars, named literally so the catalog's claimed-field check
    can see which fields this test reads."""
    for field in ("v_cigar", "d_cigar", "j_cigar"):
        value = record.get(field)
        assert value is None or isinstance(value, str), f"{field} is not a string"


def _alignment_fields_present(record):
    """The three alignment strings, named literally. See above."""
    for field in ("sequence_alignment", "germline_alignment",
                  "germline_alignment_d_mask"):
        value = record.get(field)
        assert value is None or isinstance(value, str), f"{field} is not a string"


@verifies("projection.cigar",
          claim="every column a cigar declares has the gap pattern that op "
                "requires in the two alignment strings -- M base/base, D "
                "gap/base, I base/gap -- so an op cannot be declared in a "
                "place the strings do not agree with, and every column "
                "consuming a germline base carries the reference's base there. "
                "The op lengths matching the recorded spans is true by "
                "construction in the engine and catches a field-plumbing "
                "error, not a disagreement between producers; the germline "
                "COORDINATES are not established here, because they come from "
                "the same column walk as the cigar",
          fields=("v_cigar", "d_cigar", "j_cigar"))
@pytest.mark.parametrize("stack", sorted(STACKS))
@pytest.mark.parametrize("locus", LOCI)
def test_every_cigar_describes_the_span_it_claims_against_the_reference(locus, stack):
    if "tcr" in locus and stack in B_CELL_ONLY:
        pytest.skip("somatic hypermutation does not occur in TCR sequences")
    if stack in NEEDS_D and locus not in VDJ_LOCI:
        pytest.skip(f"{locus} is a VJ locus, so there is no D to invert")
    exp = STACKS[stack](experiment_builder.recombined(locus=locus))
    checked = compared = 0
    for record, sequences, inverted in _observed(exp, n=25, seed=37):
        found = projection.disagreements(record, sequences, inverted=inverted)
        assert found == (), (
            f"{locus}/{stack}: {len(found)} projection disagreement(s); "
            f"first: {found[0]}")
        _cigar_fields_present(record)
        # Without this the test is green on records where not one column was
        # compared: the oracle returned () for an absent alignment and the
        # presence helper passes on None, so the whole body agreed about
        # nothing. Both are fixed, and this is the floor that keeps them so.
        columns = projection.columns_compared(record, sequences, inverted=inverted)
        assert sum(columns.values()) > 100, (
            f"{locus}/{stack}: only {sum(columns.values())} columns compared")
        compared += sum(columns.values())
        checked += 1
    assert checked == 25, f"{locus}/{stack}: only {checked} records were checked"
    assert compared > 5000, (
        f"{locus}/{stack}: {compared} columns compared across 25 records, which "
        "is too few for the claim to rest on")


@verifies("projection.germline_alignment",
          claim="the two alignment strings are the same length, are longer than "
                "the read by exactly the deletions the cigars record, reduce to "
                "the read itself once gaps are removed -- to its reverse "
                "complement when the read is on the other strand, because the "
                "alignment is in germline orientation -- and every column that "
                "consumes a germline base carries the reference's base there, "
                "deletion columns included. The D mask is checked only for "
                "covering exactly the D span: that its bases become N is how "
                "the engine builds it in the same walk, so asserting it "
                "restates the construction rather than establishing it",
          fields=("sequence_alignment", "germline_alignment",
                  "germline_alignment_d_mask"))
@pytest.mark.parametrize("stack", sorted(STACKS))
@pytest.mark.parametrize("locus", LOCI)
def test_the_alignment_strings_reduce_to_the_read_and_the_reference(locus, stack):
    if "tcr" in locus and stack in B_CELL_ONLY:
        pytest.skip("somatic hypermutation does not occur in TCR sequences")
    if stack in NEEDS_D and locus not in VDJ_LOCI:
        pytest.skip(f"{locus} is a VJ locus, so there is no D to invert")
    exp = STACKS[stack](experiment_builder.recombined(locus=locus))
    saw_gap = False
    for record, sequences, inverted in _observed(exp, n=25, seed=41):
        found = projection.disagreements(record, sequences, inverted=inverted)
        assert found == (), (
            f"{locus}/{stack}: {len(found)} projection disagreement(s); "
            f"first: {found[0]}")
        _alignment_fields_present(record)
        columns = projection.columns_compared(record, sequences, inverted=inverted)
        assert sum(columns.values()) > 100, (
            f"{locus}/{stack}: only {sum(columns.values())} columns compared")
        saw_gap |= "-" in record["germline_alignment"] or "-" in record["sequence_alignment"]
    if stack == "indels":
        assert saw_gap, (
            f"{locus}: no gap appeared in either alignment string, so the "
            "gap-aware half of the claim was never exercised")


def test_an_inverted_d_is_checked_rather_than_skipped():
    """The base-level oracle skips an inverted D: its coordinates are
    forward-reference while the assembled bases are the reverse complement, so
    a direct lookup does not apply.

    The projection can be checked, because the germline alignment carries the
    reverse complement too -- it has to, to align with the read. Measured with
    this fixture: the germline alignment matches the reference
    reverse-complemented within the D's own span on 407 of 407 columns, against
    100 for the plain reference. This asserts the oracle uses that and is not
    quietly passing because it skipped the segment.
    """
    exp = experiment_builder.recombined().invert_d(prob=1.0)
    inverted_seen = 0
    for record, sequences, inverted in _observed(exp, n=20, seed=43):
        assert "d" in inverted, "the fixture did not invert D"
        found = projection.disagreements(record, sequences, inverted=inverted)
        assert found == (), f"inverted D disagreement: {found[0] if found else ''}"
        columns = projection.columns_compared(record, sequences, inverted=inverted)
        assert columns["d"] > 0, (
            "no D column was compared, so the inverted case passed by being "
            "skipped rather than by agreeing")
        inverted_seen += 1
    assert inverted_seen == 20
