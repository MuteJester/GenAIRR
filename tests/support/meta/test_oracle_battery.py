"""Does the battery run what it claims, and find what it should?

The battery is about to become the judge of every generated pipeline, so its
own failure modes matter more than most. Three of them:

- a check that is wired in but never fires is a check the sweep cannot use;
- a check that fires on a correct record poisons every campaign that runs it;
- a check that quietly sits out, or runs and compares nothing, looks exactly
  like a check that passed.

Each gets a test here. The corruption cases work on the exported record as a
dict, which is what the oracles read, so they establish that the path from a
wrong value to a reported finding is connected end to end.
"""
from __future__ import annotations

import pytest

from tests.support.builders import experiment as experiment_builder
from tests.support.oracles import battery
from tests.support.pipelines import Pipeline

pytestmark = pytest.mark.behavior

#: A pipeline with no steps: every check applies, which is what a hand-built
#: fixture wants and what makes the corruption cases below unambiguous.
PLAIN = Pipeline("human_igh", ())


def _case(locus="human_igh", **passes):
    """One outcome, its record and its reference, from a plain recombination."""
    experiment = experiment_builder.recombined(locus=locus)
    for name, kwargs in passes.items():
        experiment = getattr(experiment, name)(**kwargs)
    outcome = experiment.run(n=1, seed=5)[0]
    return outcome, outcome.airr_record(experiment.refdata), experiment.refdata


def test_a_clean_record_produces_no_findings_and_runs_every_check():
    outcome, record, refdata = _case()
    result = battery.run(outcome, record, refdata, PLAIN)
    assert result.findings == (), result.report()
    assert set(result.ran) == {check.name for check in battery.CHECKS}, (
        f"only {sorted(result.ran)} ran on a plain pipeline")
    assert result.skipped == {}


def test_a_check_that_compared_nothing_says_so():
    """`p_tip` has nothing to compare without a P addition, and the result
    reports that rather than counting it as a pass. Without this the sweep
    would read ten green checks where nine ran."""
    outcome, record, refdata = _case()
    result = battery.run(outcome, record, refdata, PLAIN)
    assert "p_tip" in result.vacuous, result.report()
    assert "compared nothing" in result.report()


def test_a_p_run_is_compared_when_the_pipeline_makes_one():
    """And the same check stops being vacuous once there is a run to check,
    which is what makes the vacuity report above a measurement rather than a
    permanent excuse."""
    experiment = experiment_builder.with_p_nucleotides()
    for outcome in experiment.run(n=20, seed=1):
        record = outcome.airr_record(experiment.refdata)
        if not record.get("p_v_3_length"):
            continue
        result = battery.run(outcome, record, experiment.refdata, PLAIN)
        assert result.findings == (), result.report()
        assert result.compared["p_tip"] == 1
        return
    pytest.fail("no record carried a V 3' P addition; the fixture is wrong")


@pytest.mark.parametrize("field,value,check,identity", [
    ("v_identity", 0.5, "identity_exact", "each segment's identity is its matching fraction"),
    ("np1", "TTTT", "np_fields", "the np fields are the untemplated runs' bases and lengths"),
    ("n_observable_mutations", 7, "observable_mutations",
     "observable_substitutions == the engine's n_observable_mutations"),
    ("junction_start", 3, "junction",
     "junction_start == first surviving base of the V anchor codon"),
    ("sequence_length", 11, "coordinates", "the read ends where the last segment ends"),
])
def test_each_check_fires_on_a_record_corrupted_where_it_looks(
        field, value, check, identity):
    """One field rewritten, one check expected to notice.

    Parametrised over the five checks whose identity a single wrong field
    reaches. A check that cannot be made to fire from the record is a check
    the sweep would run for nothing, and this is what would say so.
    """
    outcome, record, refdata = _case()
    assert battery.run(outcome, record, refdata, PLAIN).findings == (), (
        "the control record already has findings")
    corrupted = {**record, field: value}
    result = battery.run(outcome, corrupted, refdata, PLAIN)
    identities = {finding.identity for finding in result.findings}
    assert identity in identities, (
        f"rewriting {field} to {value!r} did not produce {identity!r}; "
        f"got {sorted(identities)}")


def test_an_identity_above_the_substitution_fraction_is_caught_even_under_indels():
    """The bound check is the one that still applies when the exact one cannot.

    An indel pipeline skips `identity_exact` because the engine's denominator
    counts columns the provenance has no position for. If the inequality were
    skipped too, the indel case would have no identity coverage at all -- so it
    is asserted here on the pipeline that skips its sibling.
    """
    # Mutations as well as indels: without a mutation pass the
    # substitution-only fraction is 1.0, so "identity reported as 1.0" is not
    # a corruption at all and the first version of this test asserted nothing.
    outcome, record, refdata = _case(mutate={"model": "s5f", "count": 10},
                                     polymerase_indels={"count": 3})
    indel_pipeline = Pipeline("human_igh", (
        ("mutate", {"model": "s5f", "count": 10}),
        ("polymerase_indels", {"count": 3})))
    clean = battery.run(outcome, record, refdata, indel_pipeline)
    assert "identity_exact" in clean.skipped
    assert "identity_bound" in clean.ran
    assert clean.findings == (), clean.report()

    corrupted = {**record, "v_identity": 1.0, "d_identity": 1.0, "j_identity": 1.0}
    result = battery.run(outcome, corrupted, refdata, indel_pipeline)
    assert any(finding.identity.startswith("an indel lowers the identity")
               for finding in result.findings), result.report()


@pytest.mark.parametrize("name,kwargs,skipped", [
    ("contaminate", {"prob": 1.0},
     {"junction", "expression", "assembly", "np_fields",
      "observable_mutations", "p_tip"}),
    ("polymerase_indels", {"count": 2},
     {"np_fields", "observable_mutations", "identity_exact", "p_tip"}),
    ("sequencing_errors", {"count": 2}, {"observable_mutations"}),
])
def test_a_pass_that_voids_a_check_makes_the_battery_skip_it_with_a_reason(
        name, kwargs, skipped):
    """Applicability, asserted rather than trusted.

    Each of these passes makes some comparison invalid -- a contaminated read
    is deliberately not its germline, an insertion is untemplated without
    being an NP region, a sequencing miscall is a visible difference that is
    not hypermutation. The battery must skip exactly those and keep the rest,
    because a skip that is too broad loses coverage silently and a skip that is
    too narrow reports noise as findings.
    """
    pipeline = Pipeline("human_igh", ((name, kwargs),))
    outcome, record, refdata = _case(**{name: kwargs})
    result = battery.run(outcome, record, refdata, pipeline)
    assert set(result.skipped) == skipped, (
        f"{name} skipped {sorted(result.skipped)}, expected {sorted(skipped)}")
    for reason in result.skipped.values():
        assert len(reason.split()) >= 8, f"the reason is too thin: {reason!r}"
    assert result.ran, f"{name} left the battery with nothing to run"


def test_every_reporting_oracle_is_in_the_battery_or_declared_out_of_it():
    """A new oracle must be wired in or excluded on purpose.

    The battery is the only place that runs everything, so an oracle nobody
    added to it is an oracle the sweep does not use -- which looks like
    coverage from the outside, since the module exists and has its own tests.
    """
    from tests.support import oracles

    reporting = set(oracles.reporting_modules())
    wired = {check.oracle for check in battery.CHECKS}
    unaccounted = sorted(reporting - wired - set(battery.NOT_IN_BATTERY))
    assert unaccounted == [], (
        f"oracles neither in the battery nor declared out of it: {unaccounted}")
    for name, reason in battery.NOT_IN_BATTERY.items():
        assert len(reason.split()) >= 8, f"{name}: the reason is too thin"
    stale = sorted(set(battery.NOT_IN_BATTERY) & wired)
    assert stale == [], f"declared out of the battery but wired in: {stale}"
