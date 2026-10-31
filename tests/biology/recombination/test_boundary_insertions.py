"""Where a polymerase insertion goes when it lands at a segment's outer edge.

An insertion inside the J shows up the way a consumer expects: the J's read
span covers it and the J cigar carries an `I`. An insertion at the J's 3'
boundary does not. The base is placed outside every region -- the provenance
and the record agree that it belongs to no segment -- so:

- the segment spans stop tiling the read: V, NP1, D, NP2 and J together end
  one base before `sequence_length`;
- `n_j_indels` reports 1 while the J cigar reports none, so the per-segment
  indel counter and the per-segment alignment disagree about where that base
  is.

Measured under `polymerase_indels(count=1, insertion_prob=1.0)`: on IGH, 26 of
200 records put the insertion inside the J and report it in the cigar, and 2 of
200 put it at the boundary and do not; on IGK, 13 and 1.

**This module records the behaviour; it does not endorse it.** Whether the
engine should extend the J span and emit the `I`, or attribute the base to no
segment and leave `n_j_indels` at zero, is a question for the owner. What is
not in question is that the two exported views currently disagree, and a
consumer reconstructing a molecule from segment spans silently loses the base.

Found by the pipeline-space sweep, which drew
`polymerase_indels(insertion_prob=1.0)` with `productive_only()` -- a
combination no fixture in this suite had run.
"""
from __future__ import annotations

import pytest

from tests.support.builders import experiment as experiment_builder
from tests.support.oracles import reconstruction

pytestmark = pytest.mark.behavior

RECORDS = 200
SEED = 0


def _insertion_records(locus):
    experiment = (experiment_builder.recombined(locus=locus)
                  .polymerase_indels(count=1, insertion_prob=1.0))
    for outcome in experiment.run(n=RECORDS, seed=SEED):
        yield outcome, outcome.airr_record(experiment.refdata)


@pytest.mark.parametrize("locus,inside,boundary", [("human_igh", 20, 1),
                                                   ("human_igk", 10, 1)])
def test_an_insertion_inside_a_segment_is_reported_in_its_cigar(
        locus, inside, boundary):
    """The two populations, counted, so a change in the split is visible.

    The floors are below the measured counts (26 and 2 on IGH, 13 and 1 on
    IGK) because the point is the existence of both populations, not the exact
    draw. A run where the boundary population vanished would mean the engine
    changed where it places such an insertion -- which is the change this
    module exists to notice.
    """
    with_cigar = at_boundary = 0
    for _outcome, record in _insertion_records(locus):
        if not record.get("n_j_indels"):
            continue
        trailing = record["sequence_length"] - record["j_sequence_end"]
        if "I" in (record["j_cigar"] or ""):
            with_cigar += 1
            assert trailing == 0, (
                f"{locus}: an insertion reported in the J cigar also left "
                f"{trailing} base(s) outside the J span")
        else:
            at_boundary += 1
            assert trailing == 1, (
                f"{locus}: an insertion attributed to the J appears in "
                f"neither its cigar nor outside its span: trailing={trailing}, "
                f"cigar={record['j_cigar']}")
    assert with_cigar >= inside, (
        f"{locus}: only {with_cigar} insertions landed inside the J")
    assert at_boundary >= boundary, (
        f"{locus}: no insertion landed at the J boundary, so the case this "
        "module records is no longer reachable from this fixture")


def test_a_boundary_insertion_belongs_to_no_segment_in_the_provenance_either():
    """The record is not lying: the engine's own per-base record agrees.

    This is what makes the behaviour a design question rather than a defect in
    the projection. The last base carries no germline coordinate and no
    region, and `j_sequence_end` reports exactly where the J region ends.
    """
    for outcome, record in _insertion_records("human_igh"):
        trailing = record["sequence_length"] - record["j_sequence_end"]
        if not record.get("n_j_indels") or trailing != 1:
            continue
        simulation = outcome.final_simulation()
        segments = reconstruction.segments_by_position(simulation)
        last = len(simulation.bases()) - 1
        assert segments.get(last) is None, (
            f"the trailing base is in segment {segments.get(last)!r}, so the "
            "record's J span is under-reporting it rather than the base being "
            "outside every region")
        assert simulation.germline_position(last) is None, (
            "the trailing base carries a germline coordinate, so it is not an "
            "insertion at all")
        assert chr(simulation.bases()[last]).isupper(), (
            "the trailing base is lowercase, so the case channel calls it "
            "templated while the provenance calls it untemplated")
        return
    pytest.fail("no boundary insertion appeared in 200 records; see the counts "
                "in the test above")


def test_the_segment_spans_tile_the_read_except_for_a_boundary_insertion():
    """The consumer-visible consequence, stated as the arithmetic it breaks.

    A consumer that rebuilds a molecule from the segment spans and the np
    lengths gets `sequence_length` on every record except these, where it is
    one base short. Asserted in both directions so the exception stays
    exactly as wide as it is.
    """
    short = exact = 0
    for _outcome, record in _insertion_records("human_igh"):
        spans = sum(record[f"{segment}_sequence_end"]
                    - record[f"{segment}_sequence_start"]
                    for segment in ("v", "d", "j")
                    if record.get(f"{segment}_sequence_start") is not None)
        gaps = (record["np1_length"] or 0) + (record["np2_length"] or 0)
        rebuilt = spans + gaps
        trailing = record["sequence_length"] - record["j_sequence_end"]
        if trailing:
            short += 1
            assert rebuilt == record["sequence_length"] - trailing, (
                f"a boundary insertion left {record['sequence_length'] - rebuilt} "
                f"bases unaccounted for, not {trailing}")
        else:
            exact += 1
            assert rebuilt == record["sequence_length"], (
                f"the spans and np lengths rebuild {rebuilt} of "
                f"{record['sequence_length']} bases on a record with no "
                "boundary insertion")
    assert short >= 1 and exact >= 150, (
        f"{short} short and {exact} exact; both populations are needed for "
        "this to say anything")
