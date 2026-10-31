"""Does a pass change only the truth it owns?

Every oracle in this suite reads one record. None of them can ask this
question, because the answer lives in the difference between two runs: a
record whose germline coordinates were quietly perturbed by a sequencing-error
pass is perfectly self-consistent, and every per-record check passes it.

The property is worth stating because a benchmark depends on it. Scoring an
aligner's germline assignment on a corrupted read only means something if the
corruption did not move the germline truth it is scored against. Nothing in
this suite said so before.

What is compared: the same seed, with and without the pass, field by field,
record by record. What is asserted:

- every field of `RECOMBINATION_FACTS` the pass moved is one the registry
  declares it may move, with a reason;
- every field the registry declares it may move actually moved, so a
  declaration cannot rot into a blanket allowance;
- a pass with an off switch, set to off, leaves every exported field
  identical to a run without the pass -- unless the registry declares
  otherwise with a measurement;
- a pass moves no counter belonging to a mechanism it does not exercise,
  which is a property only expressible because the registry knows each pass's
  mechanisms and the catalog knows each mechanism's fields.

That last one found the exception: `invert_d(prob=0.0)` is not the same as no
`invert_d`. See `zero_differs` in the registry for the measurement and the
mechanism.

**Why the comparison is valid.** The engine addresses its choices, so adding a
pass does not disturb the draws of the passes around it -- which is what makes
"same seed, pass toggled" a controlled experiment rather than two unrelated
repertoires. The suite already leans on this in the cross-run basis tests.
Where it does not hold, that is itself the finding, and `invert_d` is the one
place it does not.
"""
from __future__ import annotations

import pytest

from tests.support.pipelines import Pipeline, for_locus
from tests.support.pipelines.registry import RECOMBINATION_FACTS

pytestmark = pytest.mark.behavior

LOCUS = "human_igh"
RECORDS = 8
SEED = 17

#: A parameter value that makes each pass actually do something. The registry's
#: strategies are for generating variety; here the point is to make the pass
#: fire as hard as it can, so a field it perturbs cannot escape by the draw
#: coming up zero.
LOUD = {"count": 4, "length": 4, "prob": 1.0, "insertion_prob": 0.5}

#: Passes this comparison cannot be run on, with the reason. `productive_only`
#: changes WHICH rearrangements are drawn rather than what happens to one, so
#: the two runs hold different rearrangements and a field-by-field diff
#: compares unrelated records: measured 6 of the protected fields moving, all
#: of them because the repertoire is different. Its own property -- that every
#: record satisfies the constraint -- is tested in the expression suite.
NOT_COMPARABLE = {
    "productive_only": "constrains which rearrangements are drawn, so the two "
                       "runs are different repertoires rather than the same "
                       "one with and without a pass",
}


def _records(steps):
    experiment = Pipeline(LOCUS, steps).build()
    return [outcome.airr_record(experiment.refdata)
            for outcome in experiment.run(n=RECORDS, seed=SEED)]


def _loud_arguments(one):
    """The pass's parameters, turned up so the pass certainly fires.

    The registry's own `loud` declaration wins where it has one: a coding end
    emits a P nucleotide only on a negative outcome, which turning a `count`
    up will never reach. The first version of this guessed from keyword names
    alone and so declared four P-length permissions it never exercised --
    which the non-vacuity assertion caught.
    """
    from hypothesis import find

    if one.loud:
        return dict(one.loud)
    kwargs = {key: find(strategy, lambda _value: True)
              for key, strategy in one.parameters.items()}
    for key, value in LOUD.items():
        if key in kwargs:
            kwargs[key] = value
    return kwargs


def _moved(before, after):
    """Every field that differs on any record of the pair."""
    return {field for one, other in zip(before, after)
            for field in one if one[field] != other[field]}


COMPARABLE = [one for one in for_locus(LOCUS) if one.name not in NOT_COMPARABLE]


@pytest.mark.parametrize("one", COMPARABLE, ids=[one.name for one in COMPARABLE])
def test_a_pass_moves_no_recombination_fact_it_does_not_declare(one):
    """The isolation guarantee, per pass.

    Eight of these passes declare nothing and move nothing: a library-prep
    pass may rewrite bases and everything derived from them, and may not touch
    the rearrangement. The four recombination-stage passes and the two
    end-loss passes declare exactly what they move and why.
    """
    baseline = _records(())
    after = _records(((one.name, _loud_arguments(one)),))
    moved = _moved(baseline, after)
    assert moved, (
        f"{one.name} at {_loud_arguments(one)} changed nothing at all, so "
        "this comparison establishes nothing about what it isolates")

    undeclared = sorted((moved & RECOMBINATION_FACTS) - set(one.may_move))
    assert undeclared == [], (
        f"{one.name} moved {undeclared}, which describe the rearrangement "
        "rather than what the library did to it. Either the pass should not "
        "be touching them, or the registry should declare each one with its "
        "reason.")

    unexercised = sorted(set(one.may_move) - moved)
    assert unexercised == [], (
        f"{one.name} declares it may move {unexercised} and did not, so the "
        "declaration is wider than the pass. A declaration nobody exercises "
        "is a permission nobody checked.")


def test_the_library_prep_passes_leave_the_rearrangement_alone():
    """The guarantee stated once, over the whole class, rather than per pass.

    This is the sentence a benchmark author needs: corruption does not move
    the germline truth it is scored against. Asserted as a group so the
    *number* of such passes is visible -- if it fell to one, the parametrized
    test above would still pass while the guarantee had quietly narrowed.
    """
    library_prep = [one for one in COMPARABLE if not one.may_move]
    assert len(library_prep) >= 7, (
        f"only {len(library_prep)} passes claim to leave the rearrangement "
        f"alone: {[one.name for one in library_prep]}")
    baseline = _records(())
    for one in library_prep:
        after = _records(((one.name, _loud_arguments(one)),))
        intact = {field for field in RECOMBINATION_FACTS
                  if all(a[field] == b[field] for a, b in zip(baseline, after))}
        assert intact == set(RECOMBINATION_FACTS), (
            f"{one.name} moved {sorted(set(RECOMBINATION_FACTS) - intact)}")


ZEROABLE = [one for one in for_locus(LOCUS)
            if one.zero_setting and one.name not in NOT_COMPARABLE]


@pytest.mark.parametrize("one", ZEROABLE, ids=[one.name for one in ZEROABLE])
def test_a_pass_configured_off_is_the_same_as_no_pass(one):
    """A pass set to zero must be indistinguishable from its absence.

    This is what makes a seeded ablation possible -- the same repertoire with
    and without one pass -- which is how anyone studies the effect of one
    pass, and how several tests in this suite build their controls.

    `invert_d` is the one pass where it does not hold, and the registry says
    so with the measurement. The assertion is inverted for it rather than
    skipped: if it ever starts behaving like the others, this fails and the
    declaration comes out.
    """
    keyword, off = one.zero_setting
    kwargs = _loud_arguments(one)
    kwargs[keyword] = off
    baseline = _records(())
    after = _records(((one.name, kwargs),))
    moved = sorted(_moved(baseline, after))
    if one.zero_differs:
        assert moved, (
            f"{one.name} at {keyword}={off} is now indistinguishable from its "
            "absence. That is the better behaviour -- delete the "
            "`zero_differs` declaration and this branch with it.")
        return
    assert moved == [], (
        f"{one.name} at {keyword}={off} still changed {moved}. A pass "
        "configured off should be indistinguishable from a pass that is not "
        "there, or no seeded comparison against it means anything.")


def test_the_one_pass_whose_off_switch_is_not_absence_is_pinned_exactly():
    """`invert_d(prob=0.0)` re-draws NP2, measured.

    The inversion coin flip is inlined into the V-NP1-D-NP2-J sequence
    immediately before `assemble(D)`, so its draw falls between the two NP
    generations and NP2 is re-drawn whether or not anything inverts. The V,
    the D and NP1 are untouched, which is what identifies the cause.

    Pinned at this precision because the consequence is narrow and specific: a
    seeded ablation over inversion is comparing different rearrangements, and
    nothing else about either record is wrong.
    """
    baseline = _records(())
    off = _records((("invert_d", {"prob": 0.0}),))

    assert all(record["d_inverted"] is False for record in off), (
        "prob=0.0 inverted a D, which is a different and worse problem")
    untouched = ["v_call", "d_call", "j_call", "np1_length",
                 "d_sequence_start", "d_sequence_end",
                 "d_germline_start", "d_germline_end", "d_trim_5", "d_trim_3"]
    for field in untouched:
        assert all(a[field] == b[field] for a, b in zip(baseline, off)), (
            f"{field} moved too, so the cause is not confined to the NP2 draw")
    moved = _moved(baseline, off)
    assert "np2_length" in moved, (
        "np2_length no longer moves, so the finding this pins is fixed -- "
        "delete this test and invert_d's `zero_differs` declaration")
    assert moved & RECOMBINATION_FACTS == {"np2_length"}, (
        f"the re-draw now reaches further than NP2: "
        f"{sorted(moved & RECOMBINATION_FACTS)}")


def _counter_owners():
    """``counter field -> the mechanisms answerable for it``.

    Read from the catalog rather than listed here, so a mechanism that gains a
    counter is covered by this property the moment it is declared.
    """
    from tests.biology import catalog

    owners: dict[str, set[str]] = {}
    for mechanism in catalog.all_mechanisms():
        for field in mechanism.truth_fields:
            if field.startswith("n_"):
                owners.setdefault(field, set()).add(mechanism.id)
    return owners


@pytest.mark.parametrize("one", COMPARABLE, ids=[one.name for one in COMPARABLE])
def test_a_pass_moves_no_counter_belonging_to_another_mechanism(one):
    """A pass's counters are its own.

    This property is only expressible because two structures exist: the
    registry knows which mechanisms each pass exercises, and the catalog knows
    which fields each mechanism is answerable for. Neither could state it
    alone.

    What it rules out is a counter drifting under a pass that has no business
    touching it -- a sequencing-error pass nudging `n_mutations`, say, which
    would make every hypermutation benchmark on a corrupted read wrong in a
    way no single record reveals. Measured: every pass's counter movement is
    confined to its own mechanisms, with `mutate` moving sixteen,
    `polymerase_indels` two, `pcr_amplify` and `sequencing_errors` one each,
    and the rest none.
    """
    owners = _counter_owners()
    baseline = _records(())
    after = _records(((one.name, _loud_arguments(one)),))
    moved = _moved(baseline, after) & set(owners)
    foreign = sorted(field for field in moved
                     if not (owners[field] & set(one.mechanisms)))
    assert foreign == [], (
        f"{one.name} moved {foreign}, which belong to mechanisms it does not "
        f"exercise ({sorted(one.mechanisms)}). Either the pass is reaching "
        "into another mechanism's truth, or the registry understates what it "
        "does.")


def test_the_pass_with_no_counter_of_its_own_moves_none():
    """`ambiguous_base_calls` masks bases and is reported on no counter.

    Its count lives in the trace (`corrupt.ns.count`) and nowhere in the
    record, which is a declared truth gap on `assay.ambiguous_calls` rather
    than an oversight -- an N is not a mutation and the suite is explicit that
    counting it as one would be wrong. Asserted here because the property
    above is satisfied trivially by a pass that moves no counter at all, and
    the reason it moves none should be visible rather than inferred.
    """
    from tests.biology import catalog

    mechanism = catalog.by_id("assay.ambiguous_calls")
    assert mechanism.truth_gap, (
        "assay.ambiguous_calls no longer declares a truth gap, so its count "
        "may now be exported -- in which case this test should check it moves")
    owners = _counter_owners()
    baseline = _records(())
    after = _records((("ambiguous_base_calls", {"count": 4}),))
    moved = _moved(baseline, after)
    assert moved, "masking four bases changed nothing at all"
    assert moved & set(owners) == set(), (
        f"a counter moved for a pass that has none: "
        f"{sorted(moved & set(owners))}")
