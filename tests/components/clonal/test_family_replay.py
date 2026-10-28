"""Whole-family replay boundaries, independent output checks and hostile inputs."""
from copy import deepcopy
from dataclasses import replace
import json
import subprocess
import sys

import pytest

from GenAIRR import Experiment, SimulationResult, _engine
from GenAIRR.experimental import FamilyReplayBundle, FamilyReplayError, FamilyReplayLimits
from GenAIRR._family_replay import _bundle_from_result, _content_decode, _digest
from tests.support.cases.duplex import references, model, history
from tests.support import paths

pytestmark = [pytest.mark.behavior, pytest.mark.slow]


def experiment(*, rows=None, n_clones=2, depth=2, rate=.05, inverted=False, **family):
    refs = references()
    exp = Experiment.on(refs).recombine(duplex_repair=model(refs, rows or [history(restore=True)]))
    if inverted:
        exp.invert_d(prob=1)
    exp.clonal_lineage(n_clones=n_clones, max_generations=depth, n_max=8,
                       n_sample=6, rate=rate, lambda_base=family.pop("lambda_base", 3), **family)
    return refs, exp


@pytest.fixture
def captured():
    refs, exp = experiment()
    result = exp.productive_only().with_metadata(
        assay={"bool": False, "int": 0, "float": 0., "negative_zero": -0.,
               "hex": "8000000000000000", "nested": [None, "é", {"t": "f64"}]},
    ).run_records(seed=7, strict=True, expose_provenance=True, capture_family_replay=True)
    return refs, result, _bundle_from_result(result)


def packed(value):
    return _engine._family_replay_pack(value)


def alter(result, edit):
    payload = json.loads(result._family_replay_seal)
    edit(payload)
    return FamilyReplayBundle._from_payload(payload, FamilyReplayLimits())


@pytest.mark.parametrize("strict", [False, True])
@pytest.mark.parametrize("read", ["identity", "substitutions", "ancestry"])
def test_capture_changes_no_draws_and_replays_all_observation_stages(strict, read):
    refs, exp = experiment()
    if read == "substitutions":
        exp.pcr_amplify(count=1).sequencing_errors(count=1)
    elif read == "ancestry":
        exp.polymerase_indels(count=2).ambiguous_base_calls(count=1).end_loss_5prime(length=1).end_loss_3prime(length=1).random_strand_orientation(prob=1)
    else:
        exp.productive_only()
    normal = exp.run_records(seed=19, strict=strict, validate_records=True)
    result = exp.run_records(seed=19, strict=strict, validate_records=True, capture_family_replay=True)
    assert packed(normal.records) == packed(result.records)
    bundle = _bundle_from_result(result)
    with pytest.raises(FamilyReplayError, match="not_verified"):
        _ = bundle.refdata
    restored = bundle.replay_records()
    assert packed(restored.records) == packed(result.records)
    assert restored.refdata.content_hash() == refs.content_hash()
    assert restored.validate_records(bundle.refdata).ok
    assert restored.validate_families_with_parents(bundle.refdata).ok
    assert len(restored.parents) == len(restored.families) == len(restored.lineage_trees) == 2
    # This literal control also checks actual inherited event counts, rather
    # than deriving mutations by subtracting the final nucleotide strings.
    for record, outcome in zip(restored.records, restored.outcomes):
        nodes = restored.lineage_trees[record["clone_id"]].nodes()
        node = nodes[record["lineage_node_id"]]
        cumulative = 0
        while node.parent_id is not None:
            cumulative += node.mutations_from_parent
            node = nodes[node.parent_id]
        assert record["n_mutations"] == cumulative
        assert record["duplicate_count"] == record["lineage_abundance"]
        if read == "identity":
            assert record["productive"] is True
            assert record["p_v_3_length"] is None
        elif read == "ancestry":
            assert record["n_indels"] == 2
            assert record["end_loss_5_length"] == record["end_loss_3_length"] == 1


def test_bundle_is_self_contained_in_fresh_process(captured, tmp_path):
    _, result, _ = captured
    target = tmp_path / "family.json"
    result.save_family_replay(target)
    expected = tmp_path / "expected.json"
    expected.write_text(packed(result.records))
    code = """
import sys
from GenAIRR import Experiment, _engine
from GenAIRR.experimental import FamilyReplayBundle
import GenAIRR._s5f_loader as loader
def forbidden(*args, **kwargs):
    raise AssertionError('original experiment or named kernel loader used')
Experiment.on = forbidden
loader.load_builtin_s5f_kernel = forbidden
bundle = FamilyReplayBundle.read_from(sys.argv[1])
result = bundle.replay_records()
assert _engine._family_replay_pack(result.records) == open(sys.argv[2]).read()
assert result.validate_records(bundle.refdata).ok
assert result.validate_families_with_parents(bundle.refdata).ok
assert bundle.replay_node(clone_id=0,node_id=0).airr_record(bundle.refdata)['sequence'] == 'AAATGTCCCTTTAAA'
print('fresh-process exact records and native evidence passed')
"""
    completed = subprocess.run([sys.executable, "-c", code, str(target), str(expected)],
                               cwd=paths.repo_root(), text=True,
                               capture_output=True, check=True)
    assert "exact records" in completed.stdout


@pytest.mark.parametrize("inverted", [False, True])
@pytest.mark.parametrize("order", [("vd", "dj"), ("dj", "vd")])
def test_source_row_and_strand_identity_survive(inverted, order):
    refs, exp = experiment(rows=[history("rare", weight=1, inverted=inverted, order=order),
                                 history("common", weight=3, inverted=inverted, order=order, restore=True)],
                           n_clones=1, depth=1, rate=0, inverted=inverted)
    result = exp.run_records(seed=3, capture_family_replay=True)
    restored = _bundle_from_result(result).replay_records()
    assert packed(restored.records) == packed(result.records)
    assert restored.parents[0].airr_record(restored.refdata)["sequence"] == ("AAATGTGGGTTTAAA" if inverted else "AAATGTCCCTTTAAA")
    assert restored.validate_families_with_parents(refs).ok


@pytest.mark.parametrize("target", [None, "KCPFK"])
def test_selected_targets_and_unsampled_cells_are_exact(target):
    _, exp = experiment(n_clones=1, depth=2, selection_strength=1, mature_substitutions=8, target_aa=target)
    result = exp.productive_only().run_records(seed=17, capture_family_replay=True)
    bundle = _bundle_from_result(result)
    restored = bundle.replay_records()
    assert packed(restored.records) == packed(result.records)
    for node in restored.lineage_trees[0].nodes():
        expected = result.families[0].replay_node(node.id, result.refdata)
        cell = bundle.replay_node(clone_id=0, node_id=node.id)
        assert packed(cell.airr_record(bundle.refdata)) == packed(expected.airr_record(result.refdata))


def test_all_failed_attempts_recorded_and_founder_sampled_once():
    _, exp = experiment(n_clones=1, depth=1, rate=0, lambda_base=.4)
    compiled = exp.compile()
    class CountSource:
        calls = 0
        def __init__(self, wrapped):
            self.wrapped = wrapped
        def __getattr__(self, name):
            return getattr(self.wrapped, name)
        def run_for_family_capture(self, *args):
            self.calls += 1
            return self.wrapped.run_for_family_capture(*args)
    proxy = CountSource(compiled._pre)
    compiled._pre = proxy
    # This seed's first draw is extinct; the exact seed is a scheduling label,
    # while each attempt's semantic receipt is independently verified below.
    result = compiled.run_records(seed=0, capture_family_replay=True)
    assert proxy.calls == 1
    payload = json.loads(result._family_replay_seal)
    attempts = payload["clones"][0]["attempts"]
    assert len(attempts) > 1
    assert payload["clones"][0]["accepted_attempt"] == len(attempts) - 1
    assert all(not receipt["samples"] for receipt in attempts[:-1])
    assert attempts[-1]["samples"]
    assert packed(_bundle_from_result(result).replay_records().records) == packed(result.records)


def test_accepted_extinct_clones_retain_parents_trees_and_roots():
    _, exp = experiment(depth=2, rate=0, lambda_base=0, allow_extinction=True)
    result = exp.run_records(seed=0, capture_family_replay=True)
    bundle = _bundle_from_result(result)
    restored = bundle.replay_records()
    assert restored.records == []
    assert len(restored.parents) == len(restored.families) == len(restored.lineage_trees) == 2
    assert all(len(tree.nodes()) == 1 for tree in restored.lineage_trees)
    assert bundle.replay_node(clone_id=1, node_id=0).airr_record(bundle.refdata)["sequence"] == "AAATGTCCCTTTAAA"


def test_experiment_annotations_seal_once_and_direct_compiled_ignores_unused_annotations():
    _, exp = experiment(n_clones=1, depth=1, rate=0)
    exp.with_metadata(unused=object())
    direct = exp.compile().run_records(seed=7, capture_family_replay=True)
    assert all("unused" not in row for row in direct)
    assert packed(_bundle_from_result(direct).replay_records().records) == packed(direct.records)
    with pytest.raises(FamilyReplayError, match="annotations"):
        exp.run_records(seed=7, capture_family_replay=True)


def test_dynamic_metadata_lossless_and_snapshot_does_not_alias_declared_object(captured):
    _, result, bundle = captured
    restored = bundle.replay_records()
    assert packed(restored.records) == packed(result.records)
    metadata = restored[0]["assay"]
    assert type(metadata["bool"]) is bool and type(metadata["int"]) is int and type(metadata["float"]) is float
    assert packed(metadata["negative_zero"]) != packed(metadata["float"])
    assert metadata["hex"] == "8000000000000000"
    _, exp = experiment(n_clones=1, depth=1, rate=0)
    declared = {"nested": [False]}
    run = exp.with_metadata(annotation=declared).run_records(seed=3, capture_family_replay=True)
    declared["nested"].append(1)
    assert run[0]["annotation"] == {"nested": [False]}
    _bundle_from_result(run)


@pytest.mark.parametrize("edit", [
    lambda r: r.records.pop(),
    lambda r: r.records.reverse(),
    lambda r: r.records[0].__setitem__("clone_id", False),
    lambda r: r.records[0]["assay"]["nested"].append("edited"),
    lambda r: r.outcomes.reverse(),
    lambda r: r.parents.reverse(),
    lambda r: r.families.reverse(),
    lambda r: r.lineage_trees.reverse(),
    lambda r: setattr(r, "_family_replay_seal", "{}"),
])
def test_changed_mutable_result_never_becomes_a_new_generating_history(captured, edit, tmp_path):
    _, result, _ = captured
    edit(result)
    destination = tmp_path / "untouched.json"
    destination.write_text("previous artifact")
    with pytest.raises(FamilyReplayError):
        result.save_family_replay(destination, replace=True)
    assert destination.read_text() == "previous artifact"
    assert list(tmp_path.iterdir()) == [destination]


def test_returned_records_are_not_a_future_replay_cache(captured):
    _, _, bundle = captured
    first = bundle.replay_records()
    first[0]["assay"]["nested"].clear()
    first.parents.clear()
    assert bundle.replay_records()[0]["assay"]["nested"] == [None, "é", {"t": "f64"}]


def test_atomic_no_clobber_replace_and_public_export(captured, tmp_path):
    _, result, bundle = captured
    target = tmp_path / "family.json"
    target.write_text("existing")
    with pytest.raises(FamilyReplayError, match="io"):
        result.save_family_replay(target)
    assert target.read_text() == "existing"
    result.save_family_replay(target, replace=True)
    # Artifacts are UTF-8 bytes; the platform default (cp1252 on Windows) is not.
    assert target.read_text(encoding="utf-8") == bundle.to_json()
    restored = FamilyReplayBundle.read_from(target).replay_records()
    restored.to_tsv(tmp_path / "out.tsv")
    assert (tmp_path / "out.tsv").read_text(encoding="utf-8").splitlines()[0].startswith("sequence_id\t")
    assert not list(tmp_path.glob(".*.tmp"))


def test_uncaptured_and_records_only_results_cannot_export_histories(tmp_path):
    _, exp = experiment(n_clones=1, depth=1)
    result = exp.run_records(seed=7)
    for uncaptured in (result, SimulationResult(result.records)):
        with pytest.raises(FamilyReplayError, match="missing_capture"):
            uncaptured.save_family_replay(tmp_path / "no.json")
    assert not list(tmp_path.iterdir())


@pytest.mark.parametrize("kind", ["fixed", "repertoire", "legacy"])
def test_unsupported_capture_rejects_before_sampling(kind, monkeypatch):
    refs = references()
    exp = Experiment.on(refs)
    exp.recombine(**({} if kind == "legacy" else {"duplex_repair": model(refs)}))
    if kind == "fixed":
        with pytest.warns(DeprecationWarning):
            exp.expand_clones(n_clones=1, per_clone=2)
    elif kind == "repertoire":
        exp.clonal_repertoire(n_clones=1)
    else:
        exp.clonal_lineage(n_clones=1)
    import GenAIRR._s5f_loader as loader
    monkeypatch.setattr(loader, "load_builtin_s5f_kernel", lambda *a, **k: pytest.fail("sample preparation reached"))
    with pytest.raises(FamilyReplayError, match="unsupported"):
        exp.run_records(seed=7, capture_family_replay=True)


def test_capture_growth_bound_is_only_a_computational_guard():
    _, exp = experiment(n_clones=1, depth=1, rate=0, lambda_base=65)
    assert exp.run_records(seed=7)
    with pytest.raises(FamilyReplayError, match="offspring rate bound"):
        exp.run_records(seed=7, capture_family_replay=True)


@pytest.mark.parametrize("edit,category", [
    (lambda p: p["clones"][0].__setitem__("clone_id", False), "format"),
    (lambda p: p["clones"][0]["observations"].reverse(), "invalid_scope"),
    (lambda p: p["clones"][0]["attempts"][0]["samples"][0].__setitem__("node_id", 0), "invalid_evidence"),
    (lambda p: p["clones"][0]["attempts"][0]["growth"][0].__setitem__("parent_id", 1), "invalid_evidence"),
    (lambda p: p["clones"][0]["source"].__setitem__("strict", False), "plan_mismatch"),
    (lambda p: p["execution"]["family_config"].__setitem__("rate", "0000000000000000"), "plan_mismatch"),
    (lambda p: p["execution"]["family_config"].__setitem__("max_generations", True), "plan_mismatch"),
    (lambda p: p["clones"][0]["observations"].pop(), "invalid_order"),
    (lambda p: p["clones"][0]["attempts"].append(deepcopy(p["clones"][0]["attempts"][0])), "invalid_order"),
    (lambda p: p["numeric_profile"].__setitem__("unqualified", True), "numerical_incompatibility"),
])
def test_semantically_invalid_evidence_rejects_after_all_digests_are_recomputed(captured, edit, category):
    _, result, _ = captured
    forged = alter(result, edit)
    with pytest.raises(FamilyReplayError) as error:
        forged.replay_records()
    assert error.value.category == category
    assert error.value.path


def test_informative_provenance_does_not_disqualify_compatible_execution(captured):
    _, result, _ = captured
    def edit(payload):
        payload["provenance"].update(package_version="other-compatible-build", engine_version="other", extension_sha256=None)
    assert packed(alter(result, edit).replay_records().records) == packed(result.records)


@pytest.mark.parametrize("field,value", [("schema_version", True), ("schema_version", 999),
                                         ("capabilities", ["future_required_semantics"])])
def test_unknown_bundle_semantics_reject_before_native_execution(captured, field, value):
    wire = json.loads(captured[2].to_json())
    wire[field] = value
    wire["commitment"] = _digest({k: v for k, v in wire.items() if k != "commitment"}, "genairr_family_bundle_v1")
    with pytest.raises(FamilyReplayError, match="unsupported"):
        FamilyReplayBundle.from_json(json.dumps(wire))


def test_duplicate_json_keys_and_extra_components_reject(captured):
    with pytest.raises(FamilyReplayError, match="duplicate"):
        FamilyReplayBundle.from_json('{"artifact":1,"artifact":2}')
    wire = json.loads(captured[2].to_json())
    node = {"kind": "array", "items": [{"literal": "unreachable"}]}
    wire["content"].append({"id": _digest(node, "genairr_family_content_v1"), "value": node})
    wire["commitment"] = _digest({k: v for k, v in wire.items() if k != "commitment"}, "genairr_family_bundle_v1")
    with pytest.raises(FamilyReplayError, match="unreachable"):
        FamilyReplayBundle.from_json(json.dumps(wire)).replay_records()


@pytest.mark.parametrize("limit,value", [("max_file_bytes", 100), ("max_work", 1),
                                        ("max_nodes", 1), ("max_choices", 1), ("max_bases", 1),
                                        ("max_history_entries", 1), ("max_attempts", 1),
                                        ("max_retained_bytes", 100), ("max_metadata_bytes", 1),
                                        ("max_components", 1)])
def test_replay_caller_limits_are_enforced_and_failed_verification_clears_cache(captured, limit, value):
    bundle = captured[2]
    bundle.replay_records()
    with pytest.raises(FamilyReplayError) as error:
        bundle.replay_records(limits=replace(FamilyReplayLimits(), **{limit: value}))
    assert error.value.category == "resource"
    with pytest.raises(FamilyReplayError, match="not_verified"):
        _ = bundle.refdata


def test_capture_resource_failure_has_no_complete_result():
    _, exp = experiment(n_clones=1, depth=1)
    with pytest.raises(FamilyReplayError) as error:
        exp.run_records(seed=7, capture_family_replay=True,
                        family_replay_limits=replace(FamilyReplayLimits(), max_work=1))
    assert error.value.category == "resource"


def test_repeated_annotations_are_charged_before_retaining_cumulative_output(captured):
    limits = replace(FamilyReplayLimits(), max_retained_bytes=1_200_000)
    _, exp = experiment()
    exp.productive_only().with_metadata(annotation="x" * 300_000)
    with pytest.raises(FamilyReplayError, match="cumulative output byte") as error:
        exp.run_records(seed=7, strict=True, capture_family_replay=True, family_replay_limits=limits)
    assert error.value.path == "records"
    def edit(payload):
        annotations = _engine._family_replay_unpack(json.dumps(payload["annotations"]))
        annotations["annotation"] = "x" * 300_000
        payload["annotations"] = json.loads(packed(annotations))
    # The malicious expected rows stay small. Actual regenerated rows must
    # still obey their cumulative output ceiling before the final mismatch.
    with pytest.raises(FamilyReplayError, match="cumulative output byte") as error:
        alter(captured[1], edit).replay_records(limits=limits)
    assert error.value.category == "resource" and error.value.path == "records"


def test_reconstructed_source_cannot_bypass_public_bcr_only_workflow(captured):
    from GenAIRR import RefDataConfig
    from GenAIRR.experimental import DuplexHistory
    from tests.support.cases.duplex import script

    refs = RefDataConfig.vdj()
    refs.add_v_allele("TRBV*01", "TRBV", b"AAATGT", anchor=3)
    refs.add_d_allele("TRBD*01", "TRBD", b"CCC")
    refs.add_j_allele("TRBJ*01", "TRBJ", b"TTTAAA", anchor=0)
    native_source = Experiment.on(refs).recombine(duplex_repair=model(refs, [
        DuplexHistory("tcr", script(), "TRBV*01", "TRBD*01", "TRBJ*01"),
    ])).compile()._simulator
    def edit(payload):
        payload["source_plan"] = json.loads(native_source.family_replay_spec())
    with pytest.raises(FamilyReplayError, match="BCR VDJ"):
        alter(captured[1], edit).replay_records()
    _, exp = experiment(n_clones=1)
    wrapper = exp.compile()
    wrapper._pre, wrapper._refdata = native_source, refs
    with pytest.raises(FamilyReplayError, match="BCR VDJ"):
        wrapper.run_records(capture_family_replay=True)


def test_generic_native_stage_codec_does_not_expand_family_phase_scope():
    _, exp = experiment(n_clones=1)
    exp.pcr_amplify(count=1)
    compiled = exp.compile()
    # A native post-read simulator is executable, but it cannot be installed
    # as a family source or supplemented with an SHM observation stage.
    from GenAIRR._family_replay import _stage_shapes
    with pytest.raises(FamilyReplayError, match="source.plan"):
        _stage_shapes(compiled._post, None)
    refs = references()
    mutation = Experiment.on(refs).mutate(model="uniform", count=1).compile()
    with pytest.raises(FamilyReplayError, match="observations.plan"):
        _stage_shapes(compiled._pre, mutation._simulator)


@pytest.mark.parametrize("value", [0, -1, True, 1.5, 2**63])
def test_resource_limits_require_positive_checked_integers(value):
    with pytest.raises(FamilyReplayError, match="resource"):
        FamilyReplayLimits(max_work=value)


def test_requested_cell_reconstruction_uses_the_same_aggregate_replay_budget(captured):
    bundle = captured[2]
    result = bundle.replay_records()
    consumed = json.loads(bundle._verified[2].snapshot())["work"]
    limits = replace(FamilyReplayLimits(), max_work=consumed)
    # Full verification alone fits exactly; a subsequent requested cell is
    # additional work and cannot reset the original caller's ceiling.
    assert bundle.replay_records(limits=limits).records == result.records
    with pytest.raises(FamilyReplayError) as error:
        bundle.replay_node(clone_id=0, node_id=result[0]["lineage_node_id"], limits=limits)
    assert error.value.category == "resource"


def test_each_actual_growth_attempt_is_charged_once(captured):
    _, result, bundle = captured
    payload = json.loads(result._family_replay_seal)
    actual = sum(len(clone["attempts"]) for clone in payload["clones"])
    assert json.loads(result._family_replay_capture.budget.snapshot())["attempts"] == actual
    restored = bundle.replay_records(limits=replace(FamilyReplayLimits(), max_attempts=actual))
    assert packed(restored.records) == packed(result.records)
    assert json.loads(bundle._verified[2].snapshot())["attempts"] == actual
    _, exp = experiment()
    same = exp.productive_only().with_metadata(assay=result[0]["assay"]).run_records(
        seed=7, strict=True, expose_provenance=True, capture_family_replay=True,
        family_replay_limits=replace(FamilyReplayLimits(), max_attempts=actual),
    )
    assert packed(same.records) == packed(result.records)


def test_survival_does_not_copy_outcomes_but_exports_charge_actual_history(captured):
    family = captured[1].families[0]
    limits = replace(FamilyReplayLimits(), max_history_entries=1)
    budget = _engine.FamilyReplayBudget(limits._json())
    assert _engine._family_replay_has_observed(family, budget=budget)
    assert json.loads(budget.snapshot())["history_entries"] == 0
    with pytest.raises(ValueError, match="resource limit"):
        _engine._family_replay_observed_outcomes(family, budget=budget)
    budget = _engine.FamilyReplayBudget()
    copied = _engine._family_replay_node_outcomes(family, budget=budget)
    assert any(outcome is not None for outcome in copied)
    assert json.loads(budget.snapshot())["history_entries"] > 0


def test_final_validation_work_cannot_escape_cumulative_history_limit(captured):
    bundle = captured[2]
    expected = bundle.replay_records()
    consumed = json.loads(bundle._verified[2].snapshot())["history_entries"]
    assert consumed > 1
    exact = replace(FamilyReplayLimits(), max_history_entries=consumed)
    assert packed(bundle.replay_records(limits=exact).records) == packed(expected.records)
    with pytest.raises(FamilyReplayError) as error:
        bundle.replay_records(limits=replace(exact, max_history_entries=consumed - 1))
    assert error.value.category == "resource"
    assert ".projection" in error.value.path or ".validation" in error.value.path


def test_deep_content_chain_is_iterative_and_expansion_is_preflighted():
    limits = replace(FamilyReplayLimits(), max_depth=2000)
    token, entries = {"literal": None}, []
    for _ in range(1500):
        node = {"kind": "array", "items": [token]}
        identity = _digest(node, "genairr_family_content_v1")
        entries.append({"id": identity, "value": node})
        token = {"ref": identity}
    value = _content_decode(token, entries, limits, _engine.FamilyReplayBudget(limits._json()))
    depth = 0
    while isinstance(value, list):
        depth += 1
        value = value[0]
    assert depth == 1500 and value is None
    # A small DAG can describe exponentially many expanded values. Reject it
    # from checked subtree costs before allocating the expanded structure.
    token, entries = {"literal": None}, []
    for _ in range(30):
        node = {"kind": "array", "items": [token, token]}
        identity = _digest(node, "genairr_family_content_v1")
        entries.append({"id": identity, "value": node})
        token = {"ref": identity}
    with pytest.raises(FamilyReplayError, match="expanded value"):
        _content_decode(token, entries, limits, _engine.FamilyReplayBudget(limits._json()))


@pytest.mark.parametrize("node", [
    {"kind": "object", "items": [["same", {"literal": 0}], ["same", {"literal": 1}]]},
    {"kind": "array", "items": [{"ref": "0" * 64}]},
    {"kind": "array", "items": [{"literal": 0, "ref": "0" * 64}]},
])
def test_content_tokens_are_closed_and_object_keys_unique(node):
    identity = _digest(node, "genairr_family_content_v1")
    limits = FamilyReplayLimits()
    with pytest.raises(FamilyReplayError, match="format"):
        _content_decode({"ref": identity}, [{"id": identity, "value": node}], limits,
                        _engine.FamilyReplayBudget(limits._json()))
