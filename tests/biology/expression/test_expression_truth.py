"""Whether the rearrangement could produce a receptor, and what the record says."""
from __future__ import annotations

import pytest

from tests.biology.catalog import verifies
from tests.support.builders import experiment as experiment_builder

pytestmark = pytest.mark.behavior


@verifies("expression.productivity",
          claim="the productive flag agrees with the frame and stop-codon flags "
                "it is derived from, so a benchmark reading any one of them is "
                "reading a consistent assessment",
          fields=("productive", "stop_codon", "vj_in_frame"),
          basis="consistent")
def test_the_productive_flag_agrees_with_the_flags_it_derives_from():
    """A productive rearrangement must be in frame and free of a stop codon in
    the coding region. The converse does not hold -- an intact anchor is also
    required -- so this checks the implication, not an equivalence."""
    exp = experiment_builder.recombined()
    saw_productive = saw_unproductive = False
    for outcome in exp.run(n=200, seed=8):
        record = outcome.airr_record(exp.refdata)
        if record["productive"]:
            saw_productive = True
            assert record["vj_in_frame"] is True, (
                "a productive rearrangement that is out of frame cannot express")
            assert record["stop_codon"] is False, (
                "a productive rearrangement with a stop codon cannot express")
        else:
            saw_unproductive = True
    assert saw_productive and saw_unproductive, (
        "the sample did not contain both outcomes, so the implication is untested")


@verifies("expression.productivity",
          claim="productive_only() yields only productive records, so a user "
                "asking for a productive repertoire gets one",
          fields=("productive",),
          basis="exercised")
def test_productive_only_yields_only_productive_records():
    exp = experiment_builder.recombined().productive_only()
    records = [o.airr_record(exp.refdata) for o in exp.run(n=60, seed=8)]
    assert records, "no records produced"
    unproductive = [r for r in records if not r["productive"]]
    assert unproductive == [], f"{len(unproductive)} unproductive records survived"


@verifies("expression.receptor_revision",
          claim="when a rearrangement is revised the record says so and keeps the "
                "original call recoverable, so the pre-revision truth is not lost",
          fields=("receptor_revision_applied", "original_v_call"),
          basis="exercised")
def test_a_revised_receptor_keeps_its_original_call():
    exp = experiment_builder.recombined().receptor_revision(prob=1.0)
    revised = 0
    for outcome in exp.run(n=40, seed=9):
        record = outcome.airr_record(exp.refdata)
        if not record["receptor_revision_applied"]:
            continue
        revised += 1
        assert record["original_v_call"], (
            "a revised receptor records no original V call, so the truth it "
            "replaced is unrecoverable")
    assert revised, "no record was revised; the fixture establishes nothing"
