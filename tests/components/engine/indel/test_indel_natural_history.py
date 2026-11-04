"""Fixed-count natural histories, honest resource errors, and strict replay."""

import json

import pytest

import GenAIRR as ga
from GenAIRR import _engine as eng

pytestmark = pytest.mark.behavior


def experiment(*, count=2, insertion_prob=0.5, productive=True, **options):
    cfg = eng.RefDataConfig.vj()
    cfg.add_v_allele("V*01", "V", b"TGTGGG", anchor=0)
    cfg.add_j_allele("J*01", "J", b"TGGAAA", anchor=0)
    exp = (
        ga.Experiment.on(cfg)
        .recombine(np1_lengths=[(0, 1)])
        .trim(enabled=False)
        .polymerase_indels(count=count, insertion_prob=insertion_prob, **options)
    )
    if productive:
        exp.productive_only()
    return exp


@pytest.mark.parametrize("strict", [False, True])
@pytest.mark.parametrize("count", [1, 2, 3])
def test_unresolved_budget_is_never_an_empty_support_or_permissive_noop(strict, count):
    compiled = experiment(count=count, rejection_attempts=0, exact_node_budget=0).compile()
    with pytest.raises(eng.StrictSamplingError) as caught:
        compiled.simulator.run(seed=17, strict=strict)
    error = caught.value
    assert error.args == ("corrupt.indel", "corrupt.indel.count", "sampling_budget_exceeded")
    assert error.requested_count == count
    assert error.explored_states == error.state_budget == error.attempts == 0
    assert error.strategy == "natural_sequential_exact_search"


@pytest.mark.parametrize("attempts", [0, 256])
def test_exact_fallback_and_rejection_have_auditable_strategy_and_replay(attempts):
    compiled = experiment(rejection_attempts=attempts).compile()
    out = compiled.simulator.run(seed=17, strict=True)
    report = out.trace().diagnostics()[0]
    assert report["model"] == "natural_sequential_v1"
    assert report["requested_count"] == report["applied_count"] == 2
    assert report["status"] == "completed"
    assert report["choice_offset"] == len(out.trace())
    assert report["strategy"] == (
        "exact_short_tuple" if attempts == 0 else "independent_tuple_rejection"
    )
    artifact = compiled.simulator.trace_file_from(out, seed=17)
    loaded = eng.TraceFile.from_json(artifact.to_json())
    assert loaded.trace.diagnostics() == out.trace().diagnostics()
    replayed = compiled.simulator.replay_from_trace_file(loaded, strict=True)
    assert replayed.airr_record(compiled.refdata) == out.airr_record(compiled.refdata)
    assert replayed.trace().diagnostics()[0]["strategy"] == "replay"
    assert replayed.trace().diagnostics()[0]["attempts"] == 0


def edited_trace(compiled, changes):
    out = compiled.simulator.run(seed=17)
    raw = json.loads(compiled.simulator.trace_file_from(out, seed=17).to_json())
    records = raw["trace"]["choices"]
    for address, value in changes.items():
        next(record for record in records if record["address"] == address)["value"] = value
    return eng.TraceFile.from_json(json.dumps(raw))


@pytest.mark.parametrize("insertion_prob", [0.0, 1.0])
def test_replay_cannot_wrap_a_site_through_u32(insertion_prob):
    compiled = experiment(count=1, insertion_prob=insertion_prob, productive=False).compile()
    artifact = edited_trace(compiled, {"corrupt.indel.site[0]": {"Int": 2**32 + 5}})
    with pytest.raises(eng.StrictSamplingError, match="invalid_plan_state"):
        compiled.simulator.replay_from_trace_file(artifact, strict=True)


def test_replay_cannot_change_count_to_an_unsupported_positive_value():
    compiled = experiment(count=1, productive=False).compile()
    artifact = edited_trace(compiled, {"corrupt.indel.count": {"Int": 2}})
    with pytest.raises(eng.StrictSamplingError, match="replay_count_outside_positive_support"):
        compiled.simulator.replay_from_trace_file(artifact, strict=True)


def test_new_all_noop_replay_requires_a_complete_empty_support_proof():
    compiled = experiment(count=1, insertion_prob=0.0).compile()
    out = compiled.simulator.run(seed=17, strict=False)
    assert out.trace().diagnostics()[0]["status"] == "proven_empty_support"
    artifact = compiled.simulator.trace_file_from(out, seed=17)
    replayed = compiled.simulator.replay_from_trace_file(artifact, strict=True)
    assert replayed.trace().diagnostics()[0]["status"] == "proven_empty_support"
    assert replayed.trace().diagnostics()[0]["applied_count"] == 0
    assert replayed.airr_record(compiled.refdata)["productive"] is True


@pytest.mark.parametrize("option,value", [
    ("rejection_attempts", -1), ("rejection_attempts", True),
    ("exact_node_budget", -1), ("exact_node_budget", 1.5),
    ("legacy_replay_sentinels", 1),
])
def test_budget_and_compatibility_options_are_checked(option, value):
    with pytest.raises((ValueError, TypeError)):
        experiment(**{option: value})


#: How many examples of the hard case the test insists on seeing, and how far
#: it looks for them. About one seed in 900 draws a five-edit count whose
#: independent-tuple rejection phase runs out and hands the problem to the
#: exact search -- 8 in the first 7,000 -- so the case has to be hunted. The
#: hunt costs about two seconds; the seeds used to be hunted by hand and
#: written in, and then per-pass substreams moved every one of them.
FIVE_EDIT_CASES = 4
FIVE_EDIT_SCAN = 10_000


def test_five_edit_support_survives_the_initial_rejection_phase():
    """The exact search finds a five-edit history the rejection phase could not.

    The premise is narrow: the count drawn has to be 5, and the cheap
    independent-tuple rejection phase has to fail on it, so that the exact
    search is the thing being tested. Searched rather than pinned, because a
    pinned seed states the premise in a number that stops being true the next
    time the draw order moves -- which is exactly what happened to the four
    seeds this test used to name.

    Every example found is checked in full, including the strict replay
    round trip, and the test fails if the scan cannot find enough of them.
    """
    from tests.support.fixtures.budget_experiments import _vdj_recombine

    compiled = (
        _vdj_recombine(canonical=True)
        .polymerase_indels(count=(0, 5), insertion_prob=0.5)
        .productive_only()
        .compile()
    )
    found = []
    for seed in range(FIVE_EDIT_SCAN):
        outcome = compiled.simulator.run(seed=seed, strict=True)
        report = outcome.trace().diagnostics()[0]
        if not (report["requested_count"] == 5
                and report["strategy"] == "exact_search_independent_probe"):
            continue
        found.append(seed)
        row = outcome.airr_record(compiled.refdata)
        assert row["productive"] is True
        assert row["n_indels"] == 5
        assert report["applied_count"] == 5
        assert report["attempts"] > 256, (
            f"seed {seed}: the exact search ran without the rejection phase "
            "having been exhausted first, so this is not the case the test is "
            "about")
        assert 0 < report["explored_states"] < 20_000_000
        artifact = compiled.simulator.trace_file_from(outcome, seed=seed)
        replayed = compiled.simulator.replay_from_trace_file(artifact, strict=True)
        assert replayed.airr_record(compiled.refdata) == row
        assert replayed.trace().diagnostics()[0]["strategy"] == "replay"
        assert replayed.trace().diagnostics()[0]["attempts"] == 0
        if len(found) == FIVE_EDIT_CASES:
            break
    assert len(found) == FIVE_EDIT_CASES, (
        f"only {len(found)} of {FIVE_EDIT_CASES} five-edit histories outlived "
        f"the rejection phase in {FIVE_EDIT_SCAN} seeds (found {found}). "
        "Either the case got rarer or the rejection phase now solves it, and "
        "either way the exact search is no longer being exercised here.")
