import json
import re

import pytest

import GenAIRR as ga
from GenAIRR import _engine

pytestmark = pytest.mark.behavior


def experiment(v=b"AAACCCTGT"):
    ref = ga.RefDataConfig.vj()
    ref.add_v_allele("V*01", "V", v, anchor=6)
    ref.add_j_allele("J*01", "J", b"TTTGGG", anchor=0)
    return ga.Experiment.on(ref).recombine(np1_lengths=[(0, 1)]).trim(enabled=False)


def trace(compiled):
    out = compiled.simulator.run(seed=17, strict=True)
    return compiled.simulator.trace_file_from(out, seed=17)


def test_v4_records_policy_and_actual_contract_context():
    compiled = experiment().productive_only().compile()
    artifact = trace(compiled)
    assert artifact.schema_version == 4
    assert artifact.productivity_policy == "canonical_observed_vdj_v1"
    assert artifact.active_contracts == compiled.simulator.active_contracts()
    for method in (compiled.simulator.replay_from_trace_file, compiled.simulator.rerun_from_trace_file):
        row = method(artifact, strict=True).airr_record(compiled.refdata)
        assert row["productive"] is True


@pytest.mark.parametrize("field,value", [("productivity_policy", "future"),
                                          ("productivity_policy", None),
                                          ("active_contracts", None)])
def test_unknown_or_missing_v4_semantics_are_rejected(field, value):
    raw = json.loads(trace(experiment().compile()).to_json())
    raw[field] = value
    with pytest.raises(ValueError, match="(productivity policy|active contracts)"):
        _engine.TraceFile.from_json(json.dumps(raw))


def test_trace_cannot_silently_change_active_contracts():
    artifact = trace(experiment().compile())
    target = experiment().productive_only().compile()
    for method in (target.simulator.replay_from_trace_file, target.simulator.rerun_from_trace_file):
        with pytest.raises(ValueError, match="active contracts mismatch"):
            method(artifact)


def test_seed_rerun_and_choice_replay_both_check_reference_hash_and_addresses():
    artifact = trace(experiment().compile())
    target = experiment(v=b"AAAAGCTGT").compile()
    for method in (target.simulator.replay_from_trace_file, target.simulator.rerun_from_trace_file):
        with pytest.raises(ValueError, match="content.hash mismatch"):
            method(artifact)
    target = experiment().compile()
    raw = json.loads(artifact.to_json())
    raw["address_schema_version"] = 999
    changed = _engine.TraceFile.from_json(json.dumps(raw))
    for method in (target.simulator.replay_from_trace_file, target.simulator.rerun_from_trace_file):
        with pytest.raises(ValueError, match="address schema version mismatch"):
            method(changed)


def test_old_indel_signature_replay_requires_explicit_compatibility():
    compiled = experiment().polymerase_indels(count=2, legacy_replay_sentinels=True).compile()
    artifact = trace(compiled)
    raw = json.loads(artifact.to_json())
    raw["schema_version"] = 3
    raw.pop("productivity_policy")
    raw.pop("active_contracts")
    raw["pass_plan_signature"] = re.sub(
        r"sampler=natural_sequential_v1,(?:search=bounded_interleaved_v1,)?rejection_attempts=\d+,exact_node_budget=\d+,legacy_replay_sentinels=true,",
        "", raw["pass_plan_signature"])
    old = _engine.TraceFile.from_json(json.dumps(raw))
    replayed = compiled.simulator.replay_from_trace_file(old, strict=True)
    assert replayed.airr_record(compiled.refdata) == compiled.simulator.replay_from_trace_file(artifact, strict=True).airr_record(compiled.refdata)
    with pytest.raises(ValueError, match="pass plan signature mismatch"):
        compiled.simulator.rerun_from_trace_file(old)
    target = experiment().polymerase_indels(count=2).compile()
    with pytest.raises(ValueError, match="pass plan signature mismatch"):
        target.simulator.replay_from_trace_file(old)
