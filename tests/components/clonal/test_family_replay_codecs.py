"""Closed native family codecs preserve execution, identity and bounded replay."""
import json

import pytest

from GenAIRR import Experiment, _engine as ge
from tests.support.cases.duplex import history, model, references


def source(refs=None):
    refs = refs or references()
    simulator = (Experiment.on(refs)
                 .recombine(duplex_repair=model(refs, [history("plain"), history("restored", restore=True)]))
                 .compile().simulator)
    return refs, simulator


def choices(outcome):
    return [(choice.address, choice.value) for choice in outcome.trace().choices()]


@pytest.mark.parametrize("strict", [False, True])
def test_native_source_spec_reconstructs_without_external_reference_or_model(strict):
    refs, original = source()
    spec = original.family_replay_spec()
    restored = ge.CompiledSimulator.from_family_replay_spec(spec)
    assert restored.family_replay_spec() == spec
    assert restored.family_replay_refdata().content_hash() == refs.content_hash()
    for seed in [0, 31, 1000003]:
        first = original.run(seed, strict=strict)
        second = restored.run(seed, strict=strict)
        assert choices(first) == choices(second)
        assert first.airr_record(refs) == second.airr_record(restored.family_replay_refdata())
        receipt = original.family_trace_receipt(first, strict)
        replay = restored.replay_family_trace(None, receipt)
        assert choices(replay) == choices(first)
        assert replay.airr_record(refs) == first.airr_record(refs)


@pytest.mark.parametrize("strict", [False, True])
def test_actual_execution_policy_and_prepared_context_are_required(strict):
    refs, simulator = source()
    old = simulator.run(3, strict=strict)
    simulator.family_replay_spec()
    with pytest.raises(ValueError, match="captured execution"):
        simulator.family_trace_receipt(old, strict)
    actual = simulator.run(3, strict=strict)
    with pytest.raises(ValueError, match="binding/policy"):
        simulator.family_trace_receipt(actual, not strict)
    receipt = simulator.family_trace_receipt(actual, strict)
    raw = json.loads(receipt)
    assert raw["strict"] is strict
    assert raw["productivity_policy"] == "duplex_observed_vdj_v1"
    assert raw["source_trace"]["schema_version"] == 4


@pytest.mark.parametrize("strict", [False, True])
def test_postread_receipt_reconstructs_exact_current_query_and_founder_truth(strict):
    refs, founder = source()
    founder.family_replay_spec()
    initial = founder.run(2, strict=strict).final_simulation()
    plan = ge.PassPlan()
    plan.push_corrupt_pcr([(2, 1.0)])
    plan.push_corrupt_quality_rate(0.05)
    plan.push_corrupt_ns([(1, 1.0)])
    plan.push_corrupt_3prime_loss([(1, 1.0)])
    plan.push_corrupt_rev_comp(1.0)
    original = plan.compile(refdata=refs)
    spec = original.family_replay_spec()
    restored = ge.CompiledSimulator.from_family_replay_spec(spec)
    budget = ge.FamilyReplayBudget()
    first = original.run_for_family_capture(31, strict, budget, initial=initial)
    receipt = original.family_trace_receipt(first, strict, budget=budget)
    replay = restored.replay_family_trace(initial, receipt, budget=budget)
    assert first.airr_record(refs) == replay.airr_record(refs)
    assert choices(first) == choices(replay)
    assert replay.validate_record(refs) == []
    wire = json.loads(receipt)
    assert wire["source_trace"] is None and wire["source"] is False
    with pytest.raises(ValueError, match="initial/source"):
        restored.replay_family_trace(None, receipt)


@pytest.mark.parametrize("mutation", [
    lambda spec: spec.update(unknown_order=[]),
    lambda spec: spec["passes"][-1].update(unknown=True),
    lambda spec: spec["reference"].update(unknown=True),
    lambda spec: spec.update(productivity_policy="canonical_observed_vdj_v1"),
    lambda spec: spec["passes"][0]["distribution"].update(cumulative=["7ff0000000000000"]),
])
def test_closed_specs_reject_unknown_fields_and_malformed_context(mutation):
    _, simulator = source()
    raw = json.loads(simulator.family_replay_spec())
    mutation(raw)
    with pytest.raises(ValueError):
        ge.CompiledSimulator.from_family_replay_spec(json.dumps(raw))


def test_unencoded_explicit_execution_order_is_rejected_before_capture():
    refs = references()
    plan = ge.PassPlan()
    plan.push_corrupt_pcr([(1, 1.0)])
    plan.push_corrupt_quality([(1, 1.0)])
    plan.add_edge(1, 0)
    simulator = plan.compile(refdata=refs)
    assert simulator.pass_names() == ["corrupt.quality", "corrupt.pcr"]
    with pytest.raises(ValueError, match="schedule"):
        simulator.family_replay_spec()


def test_caller_budget_bounds_native_reconstruction_and_capture():
    _, simulator = source()
    spec = simulator.family_replay_spec()
    with pytest.raises(ValueError, match="resource limit"):
        ge.CompiledSimulator.from_family_replay_spec(spec, budget=ge.FamilyReplayBudget('{"max_bases":1}'))
    with pytest.raises(ValueError, match="resource limit"):
        simulator.run_for_family_capture(0, False, ge.FamilyReplayBudget('{"max_work":1}'))
    outcome = simulator.run(0)
    with pytest.raises(ValueError, match="resource limit"):
        simulator.family_trace_receipt(outcome, False, budget=ge.FamilyReplayBudget('{"max_choices":0}'))


def test_legacy_reference_delimiter_alias_cannot_cross_typed_stage_binding():
    first = references()
    second = references()
    first.set_identity(species="a\nidentity.locus=b", locus="c")
    second.set_identity(species="a", locus="b\nidentity.locus=c")
    assert first.content_hash() == second.content_hash()
    _, a = source(first)
    _, b = source(second)
    assert json.loads(a.family_replay_spec())["reference"] != json.loads(b.family_replay_spec())["reference"]
    outcome = a.run(7)
    receipt = a.family_trace_receipt(outcome, False)
    with pytest.raises(ValueError, match="context mismatch"):
        b.replay_family_trace(None, receipt)
    with pytest.raises(ValueError, match="binding/policy"):
        b.family_trace_receipt(outcome, False)


@pytest.mark.parametrize("target", ["trace", "choice"])
def test_family_attempt_rejects_unknown_fields_nested_in_legacy_trace_dtos(target):
    refs = references()
    result = (Experiment.on(refs).recombine(duplex_repair=model(refs))
              .clonal_lineage(n_clones=1, max_generations=1, n_max=4, n_sample=2,
                             rate=0.05, lambda_base=3)
              .run_records(seed=7, capture_family_replay=True))
    receipt = json.loads(result.families[0].capture_receipt())
    trace = receipt["edges"][0]["trace"]
    target_object = trace if target == "trace" else trace["choices"][0]
    target_object["unrecognized"] = "must not be silently dropped"
    with pytest.raises(ValueError, match="unknown/noncanonical"):
        ge.replay_family_outcomes(result.parents[0], refs, json.dumps(receipt))


@pytest.mark.parametrize("kind", ["pcr", "ns", "indel"])
@pytest.mark.parametrize("value", [-1, 2**32])
def test_count_descriptors_reject_positive_out_of_range_support_before_execution(kind, value):
    refs = references()
    plan = ge.PassPlan()
    if kind == "indel":
        plan.push_corrupt_indel([(1, 1.0)], insertion_prob=0.5)
    else:
        getattr(plan, "push_corrupt_" + kind)([(1, 1.0)])
    raw = json.loads(plan.compile(refdata=refs).family_replay_spec())
    raw["passes"][0]["count"]["distribution"]["values"] = [value]
    with pytest.raises(ValueError, match="count support"):
        ge.CompiledSimulator.from_family_replay_spec(json.dumps(raw))


def _helper_family():
    refs = references()
    result = (Experiment.on(refs).recombine(duplex_repair=model(refs))
              .clonal_lineage(n_clones=1, max_generations=2, n_max=4, n_sample=3,
                             rate=0.05, lambda_base=3)
              .run_records(seed=7, capture_family_replay=True))
    return refs, result.families[0]


def test_bounded_output_getters_copy_the_same_actual_outcomes_and_share_counters():
    refs, family = _helper_family()
    budget = ge.FamilyReplayBudget()
    assert ge._family_replay_has_observed(family, budget)
    observed = ge._family_replay_observed_outcomes(family, budget)
    expected = family.observed_outcomes()
    assert [o.airr_record(refs) for o in observed] == [o.airr_record(refs) for o in expected]
    copied = json.loads(budget.snapshot())["history_entries"]
    assert copied > 0
    limited = ge.FamilyReplayBudget(json.dumps({"max_history_entries": copied}))
    ge._family_replay_observed_outcomes(family, limited)
    with pytest.raises(ValueError, match="resource limit"):
        ge._family_replay_observed_outcomes(family, limited)
    nodes = ge._family_replay_node_outcomes(family, budget)
    assert len(nodes) == len(family.node_outcomes())
    assert [o is None for o in nodes] == [o is None for o in family.node_outcomes()]
    for name in ("_family_replay_observed_outcomes", "_family_replay_node_outcomes"):
        with pytest.raises(ValueError, match="resource limit"):
            getattr(ge, name)(family, ge.FamilyReplayBudget('{"max_history_entries":0}'))
    with pytest.raises(ValueError, match="resource limit"):
        ge._family_replay_has_observed(family, ge.FamilyReplayBudget('{"max_work":0}'))


def test_bounded_projection_validation_and_initial_state_helpers_preserve_truth():
    refs, family = _helper_family()
    outcome = family.observed_outcomes()[0]
    budget = ge.FamilyReplayBudget()
    assert ge._family_replay_airr_record(outcome, refs, "node", budget) == outcome.airr_record(refs, sequence_id="node")
    assert ge._family_replay_validate_record(outcome, refs, "node", budget) == []
    assert ge._family_replay_final_simulation(outcome, budget).bases() == outcome.final_simulation().bases()
    for name in ("_family_replay_airr_record", "_family_replay_validate_record"):
        with pytest.raises(ValueError, match="resource limit"):
            getattr(ge, name)(outcome, refs, "node", ge.FamilyReplayBudget('{"max_bases":0}'))
    with pytest.raises(ValueError, match="resource limit"):
        ge._family_replay_final_simulation(outcome, ge.FamilyReplayBudget('{"max_history_entries":0}'))


def test_bounded_merge_charges_actual_trace_and_event_copies():
    refs, family = _helper_family()
    base = family.observed_outcomes()[0]
    plan = ge.PassPlan()
    plan.push_corrupt_quality([(1, 1.0)])
    read = plan.compile(refdata=refs).run_from(base.final_simulation(), 31)
    budget = ge.FamilyReplayBudget()
    expected = ge.merge_lineage_corruption(base, read)
    merged = ge._family_replay_merge_outcomes(base, read, budget)
    assert merged.airr_record(refs) == expected.airr_record(refs)
    assert choices(merged) == choices(expected)
    assert json.loads(budget.snapshot())["history_entries"] > base.revision_count()
    with pytest.raises(ValueError, match="resource limit"):
        ge._family_replay_merge_outcomes(base, read, ge.FamilyReplayBudget('{"max_history_entries":0}'))
