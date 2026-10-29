"""Durable finite-source traces bind the actual natural law and provenance."""
import json

import pytest

from GenAIRR import Experiment, _engine as ge
from tests.support.cases.duplex import history, model, references

pytestmark = pytest.mark.composition


def _compile(refs, rows=None, *, weights=None):
    return (
        Experiment.on(refs)
        .recombine(duplex_repair=model(refs, rows), v_allele_weights=weights)
        .mutate(model="uniform", count=1)
        .productive_only()
        .compile(allow_curatable_refdata=True)
    )


@pytest.mark.parametrize("strict", [False, True])
def test_generated_records_stream_and_durable_replay_keep_exact_selected_history(tmp_path, strict):
    refs = references()
    compiled = _compile(refs, [history("blunt"), history("restored", restore=True, weight=3)])
    result = compiled.run_records(n=4, seed=31, strict=strict, validate_records=True)
    assert result.records == list(compiled.stream_records(n=4, seed=31, strict=strict))
    simulator = compiled.simulator
    for i, outcome in enumerate(result.outcomes):
        record = result[i]
        selection = record["duplex_repair"]["selection"]
        assert selection["history_id"] in {"blunt", "restored"}
        assert record["productive"] is True and record["p_v_3_length"] is None
        trace = simulator.trace_file_from(outcome, 31 + i)
        assert trace.productivity_policy == "duplex_observed_vdj_v1"
        target = tmp_path / f"trace-{i}.json"
        trace.write_to(str(target))
        loaded = ge.TraceFile.read_from(str(target))
        for replay in (simulator.replay_from_trace_file, simulator.rerun_from_trace_file):
            reproduced = replay(loaded, strict=strict)
            assert reproduced.airr_record(refs) == outcome.airr_record(refs)
            assert reproduced.validate_record(refs) == []


@pytest.mark.parametrize("field,value", [
    ("productivity_policy", "canonical_observed_vdj_v1"),
    ("refdata_content_hash", None),
    ("address_schema_version", 0),
    ("active_contracts", []),
])
def test_duplex_trace_context_tampering_rejects_replay_and_seed_rerun(field, value):
    refs = references()
    simulator = _compile(refs).simulator
    original = simulator.trace_file_from(simulator.run(3), 3)
    raw = json.loads(original.to_json())
    raw[field] = value
    changed = ge.TraceFile.from_json(json.dumps(raw))
    for replay in (simulator.replay_from_trace_file, simulator.rerun_from_trace_file):
        with pytest.raises(ValueError):
            replay(changed)


@pytest.mark.parametrize("schema", [1, 2, 3])
def test_legacy_trace_compatibility_cannot_erase_duplex_law_binding(schema):
    refs = references()
    simulator = _compile(refs).simulator
    raw = json.loads(simulator.trace_file_from(simulator.run(3), 3).to_json())
    raw.update(schema_version=schema, productivity_policy=None, refdata_content_hash=None,
               active_contracts=None, address_schema_version=0)
    if schema < 3:
        raw["pass_plan_signature"] = "|".join(simulator.pass_names())
    old = ge.TraceFile.from_json(json.dumps(raw))
    for replay in (simulator.replay_from_trace_file, simulator.rerun_from_trace_file):
        with pytest.raises(ValueError, match="policy"):
            replay(old)


def test_true_choice_replay_is_independent_of_the_file_seed():
    refs = references()
    simulator = _compile(refs, [history("first"), history("second", restore=True)]).simulator
    original = simulator.run(31)
    raw = json.loads(simulator.trace_file_from(original, 31).to_json())
    raw["seed"] = 1_000_003
    changed = ge.TraceFile.from_json(json.dumps(raw))
    assert simulator.replay_from_trace_file(changed).airr_record(refs) == original.airr_record(refs)
    # Seed reruns deliberately draw a new history rather than consume choices.
    assert simulator.rerun_from_trace_file(changed).airr_record(refs) != original.airr_record(refs)


def test_duplex_trace_writer_rejects_an_unrelated_ordinary_outcome():
    refs = references()
    simulator = _compile(refs).simulator
    ordinary = ge.PassPlan().compile(refdata=refs).run(3)
    assert ordinary.final_simulation().bases() == b""
    with pytest.raises(ValueError, match="source binding mismatch"):
        simulator.trace_file_from(ordinary, 3)


@pytest.mark.parametrize("changed_law", ["allele_weights", "history_weights", "history_order"])
def test_duplex_trace_binds_upstream_probabilities_and_history_law(changed_law):
    refs = references(v_count=2)
    rows = [
        history("a", v_call="V1*01"),
        history("b", v_call="V2*01", restore=True),
    ]
    first = _compile(refs, rows, weights={"V1*01": 1, "V2*01": 2}).simulator
    outcome = first.run(12)
    trace = first.trace_file_from(outcome, 12)
    if changed_law == "allele_weights":
        second = _compile(refs, rows, weights={"V1*01": 2, "V2*01": 1}).simulator
    elif changed_law == "history_weights":
        second = _compile(refs, [
            history("a", v_call="V1*01", weight=2),
            history("b", v_call="V2*01", restore=True),
        ], weights={"V1*01": 1, "V2*01": 2}).simulator
    else:
        second = _compile(refs, list(reversed(rows)), weights={"V1*01": 1, "V2*01": 2}).simulator
    with pytest.raises(ValueError, match="binding"):
        second.trace_file_from(outcome, 12)
    for replay in (second.replay_from_trace_file, second.rerun_from_trace_file):
        with pytest.raises(ValueError, match="signature"):
            replay(trace)


def test_equal_sequence_from_distinct_founders_cannot_be_merged_as_one_history():
    refs = references()
    # Use empty postplans, so the two physical histories have identical query
    # bytes and mutation counters at the attempted continuation boundary.
    def unmutated(name):
        return (Experiment.on(refs).recombine(duplex_repair=model(refs, [history(name)]))
                .compile(allow_curatable_refdata=True).simulator.run(0))
    a, b = unmutated("first"), unmutated("second")
    assert a.final_simulation().bases() == b.final_simulation().bases()
    empty = ge.PassPlan().compile(refdata=refs, allow_curatable_refdata=True)
    continuation = empty.run_from(b.final_simulation(), 0)
    with pytest.raises(ValueError, match="does not match"):
        ge.merge_outcomes(a, continuation)


def test_genotype_gene_usage_law_is_bound_at_writing_and_both_replay_entrypoints():
    refs = references(v_count=2)
    source = model(refs, [history("v1"), history("v2", v_call="V2*01")])
    def compile_usage(weights):
        plan = ge.PassPlan()
        plan.push_genotype_recombine(
            refs, (0.5, 0.5), "donor", refs.content_hash(),
            v=[(h, allele, 1, 1.0) for h in (0, 1) for allele in (0, 1)],
            d=[(h, 0, 1, 1.0) for h in (0, 1)],
            j=[(h, 0, 1, 1.0) for h in (0, 1)],
            d_required=True, v_weights=weights,
        )
        plan.push_duplex_repair_source(source.to_json())
        return plan.compile(refdata=refs, respect=ge.productive(), allow_curatable_refdata=True)
    first, second = compile_usage([1.0, 2.0]), compile_usage([2.0, 1.0])
    outcome = first.run(19)
    trace = first.trace_file_from(outcome, 19)
    with pytest.raises(ValueError, match="binding"):
        second.trace_file_from(outcome, 19)
    for replay in (second.replay_from_trace_file, second.rerun_from_trace_file):
        with pytest.raises(ValueError, match="signature"):
            replay(trace)
