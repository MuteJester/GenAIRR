"""Whether the rearrangement could produce a receptor, and what the record says."""
from __future__ import annotations

import pytest

from tests.biology.catalog import verifies
from tests.support.builders import experiment as experiment_builder
from tests.support.oracles import expression as expression_oracle

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


EXPRESSION_FIXTURES = (
    ("plain", lambda: experiment_builder.recombined()),
    ("mutated", lambda: experiment_builder.recombined().mutate(model="s5f", count=10)),
    ("light_chain", lambda: experiment_builder.recombined(locus="human_igk")),
)


@verifies("expression.productivity",
          claim="each flag is re-derivable rather than merely consistent with "
                "the others: vj_in_frame from the junction length, stop_codon "
                "from translating the coding span with an independent codon "
                "table, and productive as the conjunction of those two with "
                "both anchors surviving in the molecule. Checked under "
                "hypermutation, which is where an anchor actually gets "
                "destroyed, and on a light chain, where the J anchor residue "
                "is phenylalanine rather than tryptophan",
          fields=("productive", "stop_codon", "vj_in_frame"),
          basis="reference")
@pytest.mark.parametrize("label,make", EXPRESSION_FIXTURES,
                         ids=[name for name, _ in EXPRESSION_FIXTURES])
def test_the_productivity_flags_are_each_derivable_from_the_molecule(label, make):
    """The existing test on these fields asserts an implication between them,
    which holds whatever the assessment was. This derives each one.

    `stop_codon` is the one worth watching: it is not junction-local. Read as
    a stop in `junction_aa` -- which is what the getting-started page said
    until this commit -- it disagrees on 46 of 299 mutated records whose
    junctions translate cleanly, because the stop is elsewhere in the coding
    span and the engine's own reason for it is `coding_stop`.
    """
    exp = make()
    found, checked, skipped = [], 0, 0
    for outcome in exp.run(n=120, seed=8):
        record = outcome.airr_record(exp.refdata)
        if record["junction_length"] is None:
            # No anchor, no junction, no resolved frame. The junction oracle
            # owns that contract; see the no-junction test below.
            skipped += 1
            continue
        found.extend(expression_oracle.disagreements(record))
        for field in ("productive", "stop_codon", "vj_in_frame"):
            assert record[field] is not None, (
                f"{field} is null on a record that kept its junction")
        checked += 1
    assert found == [], f"{label}: {found[:3]}"
    assert checked > 100, f"{label}: only {checked} records were checked"


@verifies("expression.productivity",
          claim="a record whose anchor codon was trimmed away reports "
                "vj_in_frame as null rather than false -- the frame could not "
                "be resolved, which is not the same as being out of frame -- "
                "with the junction fields empty and productive false, and the "
                "reasons naming anchor_removed and vj_frame_unresolved",
          fields=("vj_in_frame", "productivity_reasons"),
          basis="provenance")
def test_an_unresolvable_frame_is_null_rather_than_false():
    """The three-valued case, which a consumer treating the flag as a boolean
    gets wrong: `None` means "could not tell" and reads as "out of frame".

    Reachable on one record in three hundred by chance, so the fixture trims
    the V 3' end hard enough to remove the anchor codon on every record. That
    is also how the junction oracle's own fault here was found -- it raised
    TypeError on such a record, and the fixtures of the test that shipped it
    never produced one.
    """
    exp = experiment_builder.recombined().trim(
        v_3=[(length, 1.0) for length in range(12, 26)])
    checked = 0
    for outcome in exp.run(n=40, seed=3):
        record = outcome.airr_record(exp.refdata)
        assert record["junction_length"] is None, (
            "the fixture kept a junction; it is meant to remove the anchor")
        assert record["vj_in_frame"] is None, (
            f"an unresolvable frame reported {record['vj_in_frame']!r} rather "
            f"than null, which a consumer reads as out-of-frame")
        assert record["productive"] is False
        reasons = set(record["productivity_reasons"] or ())
        assert {"anchor_removed", "vj_frame_unresolved"} <= reasons, (
            f"the record does not say why the frame is unresolved: {sorted(reasons)}")
        checked += 1
    assert checked == 40


@verifies("expression.receptor_revision",
          claim="the flag marks a revision that actually happened and "
                "original_v_call is the allele that would have been used "
                "without it. Against the same seed at prob=0.0: the flag is "
                "set on 40 of 40 revised records and none of the unrevised, "
                "v_call differs on 40 of 40, and original_v_call equals the "
                "unrevised run's v_call on 40 of 40",
          fields=("receptor_revision_applied", "original_v_call"),
          basis="cross-run")
def test_the_revision_flag_marks_a_real_replacement():
    """The existing test asserts `original_v_call` is truthy, which a field
    echoing the current call would satisfy.

    The expected value here comes from a second run rather than from the
    record: the same seed with the revision switched off produces the
    rearrangement that would have been, so `original_v_call` can be checked
    against the allele actually used rather than merely against being
    non-empty. Switching the pass off does not disturb the draws before it,
    which is what makes the pairing valid.
    """
    unrevised = experiment_builder.recombined().receptor_revision(prob=0.0)
    revised = experiment_builder.recombined().receptor_revision(prob=1.0)
    before = [outcome.airr_record(unrevised.refdata)
              for outcome in unrevised.run(n=40, seed=5)]
    after = [outcome.airr_record(revised.refdata)
             for outcome in revised.run(n=40, seed=5)]
    assert len(before) == len(after) == 40
    for index, (plain, applied) in enumerate(zip(before, after)):
        assert not plain["receptor_revision_applied"], (
            f"record {index} is flagged as revised with the pass switched off")
        assert not plain["original_v_call"], (
            f"record {index} carries a superseded call with nothing superseded")
        assert applied["receptor_revision_applied"], (
            f"record {index} ran the revision at prob=1.0 and is not flagged")
        assert applied["original_v_call"] == plain["v_call"], (
            f"record {index}: original_v_call is "
            f"{applied['original_v_call']!r} where the unrevised run used "
            f"{plain['v_call']!r}")
        assert applied["v_call"] != plain["v_call"], (
            f"record {index}: the revision is flagged but v_call is unchanged, "
            f"so nothing was replaced")
