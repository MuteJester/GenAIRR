"""What `vj_in_frame` means, and where it stops meaning what consumers assume.

Every tool that reads an AIRR record computes in-frame-ness the same way:
`junction_length % 3 == 0`. On records without an indel that is exactly what
GenAIRR's flag says -- 200 of 200 on human_igk -- so a benchmark author has
every reason to believe the two are interchangeable.

They are not. `vj_in_frame` is the global V-J reading frame, and a deletion in
the V body shifts it while leaving the junction window, and therefore
`junction_length`, a whole number of codons. On 54 of 200 records carrying one
indel the two disagree, always with the flag saying out-of-frame and the
length saying in-frame.

Nothing here is wrong in the engine. The record is self-consistent and says
why -- `j_tail_out_of_frame` and `coding_stop` appear in
`productivity_reasons` -- and the frame genuinely is broken, which is the more
useful fact. What was missing is any statement that the relationship breaks,
so a consumer scoring an aligner's `vj_in_frame` against ours on an
indel-bearing repertoire would count 27% false errors and have no way to know
why.

That gap is closed: `coding_frame_start` is exported, and with it the second
clause of `vj_in_frame` -- the junction beginning on a codon boundary of the
frame -- is something a consumer can compute. This module keeps the
measurement, because the naive formula every tool uses is still wrong under
indels and a reader should be able to find that stated;
test_productivity_is_reproducible.py is where the exported form is checked.

Found by the pipeline-space sweep: `polymerase_indels(count=1,
insertion_prob=0.0)` on human_igk, which no fixture in this suite had run
against the productivity flags.
"""
from __future__ import annotations

import pytest

from tests.biology.catalog import verifies
from tests.support.builders import experiment as experiment_builder

pytestmark = pytest.mark.behavior

#: Locus and seed the measurements in this module were taken on.
LOCUS = "human_igk"
SEED = 0
RECORDS = 200


def _flags(experiment, n=RECORDS, seed=SEED):
    """``(flag, length_says_in_frame)`` per record that has a junction."""
    out = []
    for outcome in experiment.run(n=n, seed=seed):
        record = outcome.airr_record(experiment.refdata)
        if record["junction_length"] is None:
            continue
        out.append((record["vj_in_frame"], record["junction_length"] % 3 == 0,
                    record))
    return out


@verifies("expression.productivity",
          claim="vj_in_frame is the global V-J reading frame and not a "
                "restatement of junction_length. Without an indel the two "
                "agree on every record, so a consumer may use either; with "
                "one they disagree on about a quarter, always with the flag "
                "reporting out-of-frame where the junction length alone would "
                "say in-frame, because an indel upstream of the junction "
                "shifts the frame without changing the window. The record "
                "says why in productivity_reasons",
          fields=("vj_in_frame",),
          basis="reference")
def test_the_frame_flag_is_the_reading_frame_and_not_the_junction_length():
    """Both halves, because only the pair is informative.

    The agreement half is what makes the disagreement half a hazard: a
    consumer who checks the relationship on a clean run finds it exact.
    """
    clean = _flags(experiment_builder.recombined(locus=LOCUS))
    assert len(clean) > 150, f"only {len(clean)} clean records had a junction"
    mismatched = [record for flag, by_length, record in clean
                  if bool(flag) != by_length]
    assert mismatched == [], (
        f"{len(mismatched)} of {len(clean)} records disagree with no indel "
        "pass at all, which would make the relationship below unremarkable")

    with_indels = _flags(experiment_builder.recombined(locus=LOCUS)
                         .polymerase_indels(count=1, insertion_prob=0.0))
    disagreeing = [record for flag, by_length, record in with_indels
                   if bool(flag) != by_length]
    assert len(disagreeing) > 20, (
        f"only {len(disagreeing)} of {len(with_indels)} indel-bearing records "
        "disagree; the fixture is no longer reaching the case this documents")

    for record in disagreeing:
        assert record["vj_in_frame"] is False and record["junction_length"] % 3 == 0, (
            "the disagreement reversed direction: the flag now says in-frame "
            "where the junction length says otherwise, which would mean the "
            "flag is reading the length after all")
        assert record["productivity_reasons"], (
            "a record whose frame flag contradicts its own junction length "
            "says nothing about why")
        assert any(reason.endswith("out_of_frame")
                   for reason in record["productivity_reasons"]), (
            f"no reason explains the broken frame: {record['productivity_reasons']}")


def test_an_indel_upstream_of_the_junction_is_what_decouples_them():
    """The mechanism, not just the statistic.

    A deletion inside the V body shifts everything downstream of it; the
    junction window is downstream, so it keeps its length and loses its frame.
    Asserted by finding such a record and reading its V cigar: the deletion is
    before junction_start, which is the whole explanation.
    """
    experiment = (experiment_builder.recombined(locus=LOCUS)
                  .polymerase_indels(count=1, insertion_prob=0.0))
    for outcome in experiment.run(n=RECORDS, seed=SEED):
        record = outcome.airr_record(experiment.refdata)
        if record["junction_length"] is None:
            continue
        if record["vj_in_frame"] is not False or record["junction_length"] % 3:
            continue
        assert record["n_v_indels"] == 1, (
            f"the decoupled record's indel is not in the V: {record['v_cigar']}")
        assert "D" in record["v_cigar"] or "I" in record["v_cigar"], (
            f"v_cigar {record['v_cigar']} records no indel at all")
        consumed = 0
        import re

        for count, op in re.findall(r"(\d+)([MID])", record["v_cigar"]):
            if op == "D":
                break
            consumed += int(count)
        assert consumed < record["junction_start"], (
            f"the indel at read offset {consumed} is not upstream of the "
            f"junction at {record['junction_start']}, so the explanation in "
            "this module's docstring does not hold for this record")
        return
    pytest.fail("no record decoupled the flag from the length; see the "
                "measurements in this module's docstring")


@verifies("expression.productivity",
          claim="the documented conditions for productive -- in frame, no "
                "junction stop, V cysteine preserved, J tryptophan or "
                "phenylalanine preserved -- are necessary but NOT sufficient. "
                "The engine also requires the J tail to be in frame and "
                "reports j_tail_out_of_frame when it is not: 5 of 200 IGH "
                "records with one indel satisfy all four documented "
                "conditions and are unproductive for that reason alone, "
                "against 0 of 200 without indels. A record the engine calls "
                "productive does satisfy all four, which is the direction "
                "that holds",
          fields=("productive", "productivity_reasons"),
          basis="reference")
@pytest.mark.parametrize("locus,expected_at_least", [("human_igh", 3),
                                                     ("human_igk", 2)])
def test_the_documented_conditions_for_productive_are_not_the_whole_rule(
        locus, expected_at_least):
    """The documentation defect, as a measurement.

    The getting-started page calls these four the "full triad", which makes it
    a complete specification, and `j_tail_out_of_frame` appears nowhere in the
    documentation at all. A consumer implementing the documented rule gets a
    different answer from ours on a few percent of indel-bearing records, with
    nothing to tell them why.

    The J tail frame is germline-phase arithmetic: `(j_sequence_end -
    junction_end) % 3` explains neither direction -- 17 records one way and 55
    the other on the IGH fixture -- so this test asserts the engine's reason
    rather than recomputing the condition.
    """
    from tests.support.oracles import expression as expression_oracle

    experiment = (experiment_builder.recombined(locus=locus)
                  .polymerase_indels(count=1, insertion_prob=0.0))
    satisfied_but_unproductive = []
    productive_without_all_four = []
    for outcome in experiment.run(n=RECORDS, seed=SEED):
        record = outcome.airr_record(experiment.refdata)
        if record["junction_length"] is None:
            continue
        v_intact, j_intact = expression_oracle.anchors_intact(record)
        documented = (bool(record["vj_in_frame"]) and not record["stop_codon"]
                      and v_intact and j_intact)
        if documented and not record["productive"]:
            satisfied_but_unproductive.append(record)
        if record["productive"] and not documented:
            productive_without_all_four.append(record)

    assert productive_without_all_four == [], (
        f"{locus}: {len(productive_without_all_four)} records are productive "
        "without satisfying all four documented conditions, so the direction "
        "that was supposed to hold does not")
    assert len(satisfied_but_unproductive) >= expected_at_least, (
        f"{locus}: only {len(satisfied_but_unproductive)} records satisfied "
        "the documented conditions and were unproductive anyway; the fixture "
        "no longer reaches the case this documents")
    for record in satisfied_but_unproductive:
        assert "j_tail_out_of_frame" in record["productivity_reasons"], (
            f"{locus}: a record satisfying all four documented conditions is "
            f"unproductive for an undocumented reason: "
            f"{record['productivity_reasons']}")
