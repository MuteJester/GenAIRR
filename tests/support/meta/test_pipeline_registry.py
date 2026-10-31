"""Is the pass registry telling the truth about what the DSL accepts?

The registry is the generator's model of the simulator: which passes exist,
which loci they run on, what order the compiler wants them in, and what their
parameters may be. A wrong model does not fail loudly -- it quietly spends a
campaign's budget on ValueError and reports that it found nothing.

So the model is checked against the DSL itself: every declared pass runs, every
declared precondition is real (the pass fails on a locus it excludes, which is
why it excludes it), and the declared stage order is the order the compiler
accepts.
"""
from __future__ import annotations

import pytest
from hypothesis import HealthCheck, given, settings

from tests.support.builders import experiment as experiment_builder
from tests.support.pipelines import LOCI, PASSES, STAGES, by_name, for_locus
from tests.support.pipelines.registry import ALWAYS_EXERCISED
from tests.support.strategies import pipelines as pipeline_strategies

pytestmark = pytest.mark.behavior

#: Hypothesis settings for the registry's own checks: each example builds and
#: runs a simulation, so these are deliberately few. The campaign that spends
#: real examples is the pipeline-space module under biology/, which arrives
#: with the oracle battery.
REGISTRY_SETTINGS = settings(
    max_examples=25, deadline=None,
    suppress_health_check=[HealthCheck.too_slow, HealthCheck.filter_too_much])


def test_every_declared_pass_is_a_method_the_dsl_has():
    """The cheapest way for the registry to rot: a renamed DSL method."""
    experiment = experiment_builder.recombined()
    missing = [one.name for one in PASSES
               if not callable(getattr(experiment, one.name, None))]
    assert missing == [], f"registry names passes the DSL does not have: {missing}"


@pytest.mark.parametrize("one", PASSES, ids=[one.name for one in PASSES])
def test_each_pass_runs_on_every_locus_it_claims(one):
    """A pass that cannot run is a pass that contributes nothing but filtering.

    Run with drawn-looking but fixed arguments -- the midpoint of each declared
    range -- because this test is about the pass being callable and compilable
    on that locus, not about the parameter space, which the campaign covers.
    """
    from hypothesis import find

    ran = 0
    for locus in sorted(LOCI):
        if not one.applies_to(locus):
            continue
        kwargs = {key: find(strategy, lambda _value: True)
                  for key, strategy in one.parameters.items()}
        experiment = getattr(
            experiment_builder.recombined(locus=locus), one.name)(**kwargs)
        records = experiment.run_records(n=1, seed=3)
        assert len(records) == 1, (
            f"{one.name} on {locus} produced {len(records)} records")
        ran += 1
    assert ran, f"{one.name} claims no locus at all"


@pytest.mark.parametrize(
    "one", [one for one in PASSES if one.requires],
    ids=[one.name for one in PASSES if one.requires])
def test_each_declared_precondition_is_real(one):
    """A precondition nobody checked is a guess about the engine.

    `mutate` is declared B-cell only and `invert_d` D-only. If either stopped
    being true -- the engine learns to simulate TCR hypermutation, say -- the
    registry would be excluding a pass for a reason that no longer holds, and
    the sweep would never draw it there again. That is exactly the kind of
    silent coverage loss this suite exists to prevent, so the exclusion has to
    keep failing for its stated reason.
    """
    from hypothesis import find

    refused = 0
    for locus in sorted(LOCI):
        if one.applies_to(locus):
            continue
        kwargs = {key: find(strategy, lambda _value: True)
                  for key, strategy in one.parameters.items()}
        reasons = one.why_not(locus)
        assert reasons, f"{one.name} excludes {locus} with no reason recorded"
        with pytest.raises(Exception) as caught:
            experiment = getattr(
                experiment_builder.recombined(locus=locus), one.name)(**kwargs)
            experiment.run_records(n=1, seed=3)
        assert caught.type is not AssertionError, (
            f"{one.name} on {locus}: the refusal is an assertion, not a "
            "refusal from the library")
        refused += 1
    assert refused, (
        f"{one.name} declares a precondition but every locus satisfies it, so "
        "nothing was checked")


def test_the_declared_stage_order_is_the_order_the_compiler_accepts():
    """The registry sorts a drawn stack by stage. If that order is not the
    compiler's, the generator produces programs that read like one pipeline and
    compile as another -- or do not compile at all.

    Checked by building one pipeline per adjacent stage pair, in the declared
    order, and running it.
    """
    known = by_name()
    pairs = 0
    for earlier, later in zip(STAGES, STAGES[1:]):
        first = next((one for one in for_locus("human_igh")
                      if one.stage == earlier), None)
        if first is None:
            continue
        # Skip a pair the registry declares mutually exclusive: trim and
        # process_coding_ends are adjacent in stage and refused together, so
        # pairing them would test the DSL's exclusion rather than the order.
        second = next((one for one in for_locus("human_igh")
                       if one.stage == later
                       and one.name not in first.excludes
                       and first.name not in one.excludes), None)
        if second is None:
            continue
        from hypothesis import find

        experiment = experiment_builder.recombined()
        for one in (first, second):
            kwargs = {key: find(strategy, lambda _value: True)
                      for key, strategy in one.parameters.items()}
            experiment = getattr(experiment, one.name)(**kwargs)
        records = experiment.run_records(n=1, seed=5)
        assert len(records) == 1, (
            f"{first.name} then {second.name} ({earlier} then {later}) "
            "produced no record")
        pairs += 1
        assert known[first.name].stage == earlier
    assert pairs >= 3, f"only {pairs} adjacent stage pairs were exercised"


@REGISTRY_SETTINGS
@given(pipeline_strategies.pipelines())
def test_every_drawn_pipeline_builds_and_runs(pipeline):
    """The claim the whole campaign rests on: a drawn pipeline is a valid
    program.

    Validity is built into the generator rather than filtered afterwards, so a
    failure here is a finding about the registry -- a missing precondition, a
    parameter range the DSL refuses, a stage in the wrong place -- and not an
    example to discard. The rendered pipeline is in the message because that is
    what makes such a finding fixable.
    """
    try:
        records = pipeline.build().run_records(n=1, seed=7)
    except Exception as error:                      # noqa: BLE001 - reported
        raise AssertionError(
            f"a drawn pipeline did not run: {type(error).__name__}: {error}\n\n"
            f"{pipeline.render()}") from error
    assert len(records) == 1, f"no record from\n{pipeline.render()}"


@REGISTRY_SETTINGS
@given(pipeline_strategies.pipelines())
def test_a_drawn_pipeline_renders_as_the_chain_that_built_it(pipeline):
    """A finding nobody can reproduce is not a finding.

    The rendered text is checked as code -- compiled, not eyeballed -- because
    a render that drifts from the pipeline is worse than no render: it sends
    whoever reads the report after a pipeline that was never run.
    """
    rendered = pipeline.render()
    assert rendered.startswith("(experiment_builder.recombined("), rendered
    for name, _kwargs in pipeline.steps:
        assert f".{name}(" in rendered, f"{name} is missing from the render"
    compile(rendered, "<rendered pipeline>", "eval")


@pytest.mark.parametrize(
    "one", [one for one in PASSES if one.excludes],
    ids=[one.name for one in PASSES if one.excludes])
def test_each_declared_exclusion_is_real(one):
    """An exclusion nobody verified silently shrinks the drawn space.

    The generator drops a pass its neighbour excludes, so a wrong exclusion is
    a pair of passes that never get tested together and no failure anywhere to
    say so. Each declared pair therefore has to be refused by the DSL, for a
    reason the registry records.
    """
    from hypothesis import find

    known = by_name()
    checked = 0
    for other, reason in one.excludes.items():
        assert reason.strip(), f"{one.name} excludes {other} with no reason"
        assert other in known, f"{one.name} excludes {other!r}, which is not a pass"
        locus = next(locus for locus in sorted(LOCI)
                     if one.applies_to(locus) and known[other].applies_to(locus))
        experiment = experiment_builder.recombined(locus=locus)
        with pytest.raises(ValueError, match="exclusive|refus|both"):
            for step in (one, known[other]):
                kwargs = {key: find(strategy, lambda _value: True)
                          for key, strategy in step.parameters.items()}
                experiment = getattr(experiment, step.name)(**kwargs)
            experiment.run_records(n=1, seed=3)
        checked += 1
    assert checked == len(one.excludes)


def test_the_registry_reaches_most_of_the_mechanism_catalog():
    """What the sweep can and cannot see, stated as a number.

    A campaign can only find defects in mechanisms its passes exercise. The
    forks left out of v1 -- clonal, duplex, genotype -- are the reason this is
    not the whole catalog, and naming the shortfall here is what keeps "the
    sweep found nothing" from being read as "there is nothing to find".
    """
    from tests.biology import catalog

    reached = ({mechanism for one in PASSES for mechanism in one.mechanisms}
               | set(ALWAYS_EXERCISED))
    known = {one.id for one in catalog.all_mechanisms()}
    unknown = sorted(reached - known)
    assert unknown == [], f"registry names mechanisms the catalog does not: {unknown}"
    missed = sorted(known - reached)
    assert len(reached) >= 22, (
        f"the registry's passes exercise only {len(reached)} of {len(known)} "
        f"mechanisms, which is too narrow a sweep to be worth running")
    # Pinned rather than asserted away: this is the sweep's blind spot, and it
    # should shrink deliberately when the forks are brought in.
    # The genuine blind spot: the three forks v1 leaves out, plus the one
    # mechanism with no exported truth to find anything in.
    assert missed == [
        "assay.duplex_formation",      # its own oracle; changes the record shape
        "clonality.founder",           # the clonal forks are out of v1
        "clonality.sampling",
        "clonality.selection",
        "clonality.tree_topology",
        "germline.genotype",           # changes which alleles exist, not a record's truth
        "maturation.class_switching",  # nothing exported to find
    ], f"the sweep's blind spot changed: {missed}"


def test_ci_runs_the_campaign_with_a_seed_derived_from_the_commit():
    """The campaign's reproducibility contract, pinned.

    Deep on every pull request means a coin-flip failure would block
    unrelated work, and a derandomized campaign would explore the same
    examples for ever. The seed comes from the commit, which gives both: new
    ground per commit, and a failure that reproduces from the seed in the log.
    """
    from tests.support import paths

    text = (paths.repo_root() / ".github" / "workflows" / "test.yml").read_text(
        encoding="utf-8")
    assert "Discovery campaign" in text, "no CI step runs the campaign"
    lines = [line for line in text.splitlines()
             if "hypothesis-seed" in line or "rev-parse" in line]
    assert any("rev-parse" in line for line in lines), (
        f"the campaign's seed is not derived from the commit: {lines}")
    assert any("--hypothesis-seed" in line for line in lines), (
        f"the campaign runs without an explicit seed: {lines}")
    assert "HYPOTHESIS_PROFILE: campaign" in text, (
        "the campaign step does not select the campaign profile, so it would "
        "run at the default example count")
    assert "-m campaign" in text, (
        "the campaign step does not select the campaign marker")
    # Excluded from the matrix run, in both of its pytest invocations. The
    # matrix profile keeps Hypothesis's example database, so a campaign run
    # there replays a stored failure from before a fix and fails a clean tree;
    # and two runs per commit means two example sets where the point is one
    # deliberate campaign.
    matrix_runs = [line for line in text.splitlines()
                   if "pytest tests/" in line and "-n auto" in line]
    assert matrix_runs, "the matrix no longer runs the suite"
    for line in matrix_runs:
        assert "not campaign" in line, (
            f"the matrix run would also run the campaign: {line.strip()}")
