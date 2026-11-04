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
  identical to a run without the pass -- with a drawing pass behind it as
  well as without one;
- a pass moves no counter belonging to a mechanism it does not exercise,
  which is a property only expressible because the registry knows each pass's
  mechanisms and the catalog knows each mechanism's fields.

**What makes the comparison valid.** Every pass draws from its own random
substream, keyed by the pass's identity rather than its position in the
pipeline, so toggling one pass does not move what any other pass draws.
"Same seed, pass toggled" is therefore a controlled experiment about the
molecule and not only about the rearrangement.

That took two corrections to get right, both worth keeping here, because the
module asserted the wrong thing twice in the same place. The first version
claimed the engine already addressed its choices per pass, with `invert_d`
as the single exception -- wrong in the direction that flatters the engine,
because `receptor_revision(prob=0.0)` looked clean only while nothing drew
after it. The second version pinned the true state: one shared stream, so
any pass's presence shifted every later draw, and with a `mutate` behind
either pass the sequence differed on 200 of 200 records. It said in its own
docstring that per-pass substreams would make the sequence half start
passing and that the test should then be replaced by the stronger claim.
That is what happened, and what the last test in this module now asserts.

**What this comparison is not sensitive to.** It runs one seed. A pass whose
effect depends on where its draws land is measured only in the cases that seed
produces, and `polymerase_indels` is where that shows: it moves `np1_length`
and `np2_length` -- documented behaviour, since those fields are the interval
between the segment bounds and an insertion joins it -- and seed 17 put no
indel in an NP interval, so nothing here noticed until per-pass substreams
moved where the indels land. 9 of the first 20 seeds move them.

What that exposes is in `RECOMBINATION_FACTS`, not in the engine: the set
holds the two np lengths, which the product documents as observed intervals,
and holds no field for the amount TdT added at generation time. There is none
to hold -- `np1_n_tdt` counts the TdT bases still present, so a deletion moves
it too, on 30 of 400 records. See `polymerase_indels` in the registry.
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

    The number is at its floor, and it got there the honest way: the group
    held eight until `polymerase_indels` was measured moving `np1_length` and
    `np2_length` -- documented behaviour for an interval an insertion can join
    -- and had to declare them. So the next pass that declares anything fails
    this assertion, which is the point of writing the count down rather than
    deriving it.
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

    This is the SINGLE-PASS form: nothing draws after the pass under test.
    The two-pass form below is the one a real ablation takes, and it is the
    one that used to fail.

    There is no longer an exception. `invert_d(prob=0.0)` was one: with a
    single random stream shared by the plan, its coin flip fell between the
    two NP generations and NP2 was re-drawn whether or not anything inverted.
    Per-pass substreams removed it, so the assertion is now unconditional --
    which is why the registry's `zero_differs` field and the branch that read
    it are gone, along with the test that pinned the exception exactly.
    """
    keyword, off = one.zero_setting
    kwargs = _loud_arguments(one)
    kwargs[keyword] = off
    baseline = _records(())
    after = _records(((one.name, kwargs),))
    moved = sorted(_moved(baseline, after))
    assert moved == [], (
        f"{one.name} at {keyword}={off} still changed {moved}. A pass "
        "configured off should be indistinguishable from a pass that is not "
        "there, or no seeded comparison against it means anything.")


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


def test_the_masking_pass_moves_its_own_counter_and_no_other():
    """`ambiguous_base_calls` masks bases and now reports how many.

    This test used to assert the opposite. `assay.ambiguous_calls` declared a
    truth gap -- its count lived in the trace (`corrupt.ns.count`) and nowhere
    in the record -- and the assertion was that no counter moved at all, with
    a message saying that if the gap ever closed this test should check the
    count moves instead. The gap has closed: `n_ambiguous_bases` is exported,
    so this is that test.

    The counter is *realized*, not requested, which is the part worth pinning.
    A site drawn twice writes `N` over an `N` and masks nothing new, so a
    record can report fewer than the count asked for -- and asking for `N`
    where the reference already holds one is not a masking either. Counting
    the `N`s in the sequence is a different number again, because an allele
    may carry an `N` of its own.
    """
    from tests.biology import catalog

    mechanism = catalog.by_id("assay.ambiguous_calls")
    assert not mechanism.truth_gap, (
        "assay.ambiguous_calls declares a truth gap again, so its count is no "
        "longer exported and this test is asserting against a field that has "
        "gone away")
    assert "n_ambiguous_bases" in mechanism.truth_fields, (
        "the mechanism no longer owns n_ambiguous_bases, so this test is no "
        "longer checking the counter that belongs to it")

    owners = _counter_owners()
    baseline = _records(())
    after = _records((("ambiguous_base_calls", {"count": 4}),))
    moved = _moved(baseline, after)
    assert moved, "masking four bases changed nothing at all"

    counters_moved = moved & set(owners)
    assert counters_moved == {"n_ambiguous_bases"}, (
        f"masking moved {sorted(counters_moved)}; it owns n_ambiguous_bases "
        "and nothing else, so anything more means the pass is reaching into "
        "another mechanism's truth and anything less means its own counter "
        "stopped reporting")

    # Realized, not requested: every record's count is the number of
    # positions that actually differ from the unmasked run, which is at most
    # what was asked for and is frequently less.
    for before, after_record in zip(baseline, after):
        realized = sum(1 for old, new in zip(before["sequence"], after_record["sequence"])
                       if old != new)
        assert after_record["n_ambiguous_bases"] == realized, (
            f"n_ambiguous_bases reports {after_record['n_ambiguous_bases']} "
            f"where {realized} positions changed")
        assert after_record["n_ambiguous_bases"] <= 4, (
            "the realized count exceeded the requested count, which no "
            "sequence of draws can produce")


#: A pass that draws per base, placed after the pass under test, so anything
#: the earlier pass consumed from the stream shows up as a different molecule.
#: `uniform` rather than `s5f` because its draws are the simplest to reason
#: about: one site and one replacement base per mutation.
DOWNSTREAM = ("mutate", {"model": "uniform", "count": 4})


@pytest.mark.parametrize("one", ZEROABLE, ids=[one.name for one in ZEROABLE])
def test_an_off_pass_is_indistinguishable_from_its_absence_with_a_pass_behind_it(one):
    """The two-pass form, which is the form an ablation actually takes.

    The single-pass form above can pass for an uninteresting reason: with
    nothing drawing after the pass under test, a pass that fires nothing has
    nothing to disturb. This form puts a `mutate` behind it, so anything the
    earlier pass did to the shared state of the run shows up as a different
    molecule.

    This test used to assert the opposite, and said so: the engine held one
    random stream, every pass's presence consumed words from it, and the
    sequence therefore differed on every record. Measured then over 200
    human_igh records, for both passes with an off switch at the recombination
    stage: rearrangement identical 200/200, final sequence identical **0/200**.
    The docstring ended by saying that if per-pass substreams were ever added
    the sequence half would start passing and the test should be replaced by
    the stronger claim. They were, and this is that claim.

    Both halves are asserted, because the pair is the useful statement:

    - the REARRANGEMENT survives, as it always did;
    - the MOLECULE survives too, on every record. A seeded ablation now
      compares two runs that differ by the pass under test and by nothing
      else, which is what makes the comparison mean anything.

    Measured after the change on the same configuration: both passes identical
    on 40/40 records, every exported field included.
    """
    keyword, off = one.zero_setting
    if one.stage not in ("post_recombine", "recombination"):
        pytest.skip("a pass after `mutate` cannot shift `mutate`'s own draws; "
                    "this is about a pass that precedes the drawing one")
    kwargs = _loud_arguments(one)
    kwargs[keyword] = off

    def records(steps):
        from tests.support.builders import experiment as experiment_builder

        experiment = experiment_builder.recombined(locus=LOCUS)
        for name, arguments in steps:
            experiment = getattr(experiment, name)(**arguments)
        return [outcome.airr_record(experiment.refdata)
                for outcome in experiment.run(n=40, seed=SEED)]

    without = records((DOWNSTREAM,))
    with_off = records(((one.name, kwargs), DOWNSTREAM))

    facts = sorted(RECOMBINATION_FACTS - set(one.may_move))
    for field in facts:
        moved = [index for index, (a, b) in enumerate(zip(without, with_off))
                 if a[field] != b[field]]
        assert moved == [], (
            f"{one.name} at {keyword}={off} moved {field} on records "
            f"{moved[:5]} once a drawing pass followed it.")

    differing = [index for index, (a, b) in enumerate(zip(without, with_off))
                 if a["sequence"] != b["sequence"]]
    assert differing == [], (
        f"{one.name} at {keyword}={off} changed the molecule on records "
        f"{differing[:5]} of {len(without)} with a drawing pass behind it. A "
        "pass configured off has to be indistinguishable from a pass that is "
        "not there, molecule included, or a seeded ablation over it is "
        "comparing two different draws rather than two different pipelines.")

    moved = sorted(_moved(without, with_off) - set(one.may_move))
    assert moved == [], (
        f"{one.name} at {keyword}={off} moved {moved}, which the sequence "
        "comparison above did not catch -- so the pass is reaching something "
        "that is not derived from the bases.")
