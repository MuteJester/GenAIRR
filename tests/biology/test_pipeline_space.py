"""Run pipelines nobody wrote by hand, and check everything about the result.

Every other module here asserts that recorded truth is correct on a fixture
somebody chose. This one goes looking. It builds pass stacks across the space
the DSL accepts, runs them, and puts every applicable truth check over every
record -- so a defect in an interaction nobody thought to test shows up as a
failure with the pipeline that produced it.

The two defects this project has found the hard way were both interactions.
The junction end was wrong by the indel length when an indel landed inside the
J anchor codon; `end_loss_5_length` counted only the first truncation when
both ends were lost. Each was found by a person comparing things by hand,
once. Neither had a test, because the suite's fixtures ran one pass at a time.

Two techniques, deliberately both:

- **Exhaustive pairwise coverage.** Every legal pass and every legal pair of
  passes, on every locus they run on: 313 pipelines, deterministic, named so a
  failure says which pair. Complete 2-wise coverage of the pass space is a
  stronger statement than any number of random draws, and it is what makes
  "the sweep found nothing" mean something.
- **Property-based depth.** Hypothesis draws stacks of up to five passes with
  drawn parameters, which reaches the three-way interactions and the parameter
  values the matrix holds fixed. Shrinking turns a failure into the smallest
  pipeline that still shows it.

**Reproducibility.** The campaign runs deep on every pull request, so a
coin-flip failure would block unrelated work. The seed is therefore passed in
by CI, derived from the commit: every commit explores new ground, and every
failure reproduces exactly from the seed in the log. Locally the profile's own
seeding applies and `--hypothesis-seed` takes an explicit one.

The drawn test carries the `campaign` marker and is **excluded from the matrix
run**, which runs it once in its own step instead. Two reasons, both learned
the hard way: the matrix's profile keeps Hypothesis's example database, so a
full-suite run replayed a stored failure from before a fix and failed a clean
tree; and running it in both places means two different example sets per
commit, where the point is one deliberate campaign. The exhaustive matrix
above stays in the suite -- it is deterministic and takes 40 seconds.
"""
from __future__ import annotations

import itertools

import pytest
from hypothesis import HealthCheck, event, given, settings
from hypothesis import strategies as st

from tests.support.oracles import battery
from tests.support.pipelines import Pipeline, for_locus
from tests.support.strategies import pipelines as pipeline_strategies

pytestmark = pytest.mark.behavior

#: Records per pipeline. Three rather than one because a per-record draw --
#: whether a coin-flip pass fired, whether an indel landed in the anchor codon
#: -- is what makes a pipeline's interaction visible, and the build cost
#: dominates anyway: 0.02s to build, 0.1s per record.
RECORDS = 3

#: The seed the deterministic matrix runs at. Fixed, because the matrix's job
#: is complete pair coverage rather than variety; the drawn campaign below is
#: where variety comes from.
MATRIX_SEED = 11


def _pairs():
    """Every legal pass and pair of passes, per locus, as pytest parameters."""
    cases = []
    for locus in ("human_igh", "human_igk", "human_tcrb"):
        passes = for_locus(locus)
        for one in passes:
            cases.append(pytest.param(locus, (one.name,),
                                      id=f"{locus}-{one.name}"))
        for first, second in itertools.combinations(passes, 2):
            if second.name in first.excludes or first.name in second.excludes:
                continue
            cases.append(pytest.param(locus, (first.name, second.name),
                                      id=f"{locus}-{first.name}+{second.name}"))
    return cases


def _fixed_arguments(one):
    """A representative value per parameter, so the matrix is deterministic.

    Chosen rather than drawn: the matrix is about which passes run together,
    and a drawn parameter would make a failure's cause ambiguous between the
    pair and the value. The drawn campaign covers the values.
    """
    from hypothesis import find

    return {key: find(strategy, lambda _value: True)
            for key, strategy in one.parameters.items()}


def _check(pipeline, seed):
    """Run a pipeline and battery every record. Returns the failure, or None.

    An explicit refusal is a legal outcome, not a finding: a constrained
    sampler that cannot satisfy its constraint says so rather than emitting
    records that violate it, which is the behaviour anyone would want. Four
    deletions under `productive_only()` is the case the campaign found --
    deletions can only shorten the coding sequence, so the frame cannot
    survive them. Every OTHER exception is still a failure, because a crash
    is exactly what a campaign is for.
    """
    import GenAIRR as ga

    experiment = pipeline.build()
    try:
        outcomes = list(experiment.run(n=RECORDS, seed=seed))
    except ga.StrictSamplingError as refusal:
        event(f"refused:{type(refusal).__name__}")
        return None
    except Exception as error:                      # noqa: BLE001 - reported
        return (f"the pipeline did not run: {type(error).__name__}: {error}")
    for index, outcome in enumerate(outcomes):
        record = outcome.airr_record(experiment.refdata)
        result = battery.run(outcome, record, experiment.refdata, pipeline)
        if result.findings:
            return f"record {index} of {RECORDS}:\n{result.report()}"
    return None


@pytest.mark.parametrize("locus,names", _pairs())
def test_every_legal_pass_and_pair_produces_records_the_oracles_agree_with(
        locus, names):
    """Complete 1-wise and 2-wise coverage of the pass space.

    Deterministic, exhaustive and named: a failure identifies the pair, which
    is the hard part of an interaction defect. The parameters are fixed
    representative values, so a failure here is about the combination and not
    about a drawn number.
    """
    offered = {one.name: one for one in for_locus(locus)}
    steps = tuple((name, _fixed_arguments(offered[name])) for name in names)
    pipeline = Pipeline(locus, tuple(
        sorted(steps, key=lambda step: pipeline_stage_index(offered, step[0]))))
    failure = _check(pipeline, MATRIX_SEED)
    assert failure is None, (
        f"{failure}\n\nreproduce with seed {MATRIX_SEED}:\n{pipeline.render()}")


def pipeline_stage_index(offered, name):
    from tests.support.pipelines.registry import STAGES

    return STAGES.index(offered[name].stage)


@pytest.mark.campaign
@settings(deadline=None,
          suppress_health_check=[HealthCheck.too_slow,
                                 HealthCheck.filter_too_much,
                                 HealthCheck.data_too_large])
@given(pipeline_strategies.pipelines(min_passes=2, max_passes=5),
       st.integers(min_value=0, max_value=2**31 - 1))
def test_a_drawn_pipeline_produces_records_the_oracles_agree_with(pipeline, seed):
    """Depth: three-way and four-way interactions, with drawn parameters.

    Everything the matrix cannot reach. The example count comes from the
    Hypothesis profile, so CI runs this deep and a local run does not.

    Hypothesis events record which passes and which pairs an actual campaign
    reached, visible with --hypothesis-show-statistics. A campaign that reports
    no findings is only meaningful beside what it covered.
    """
    event(f"locus:{pipeline.locus}")
    for name in pipeline.names:
        event(f"pass:{name}")
    for first, second in itertools.combinations(sorted(pipeline.names), 2):
        event(f"pair:{first}+{second}")
    event(f"depth:{len(pipeline.steps)}")

    failure = _check(pipeline, seed)
    assert failure is None, (
        f"{failure}\n\nreproduce with seed {seed}:\n{pipeline.render()}")
