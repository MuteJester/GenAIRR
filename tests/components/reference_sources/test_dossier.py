"""The dossier manifest and the digest graph binding its members together.
"""
from __future__ import annotations
import copy
import inspect
import pytest
from GenAIRR._artifact_io import canonical_json, typed_digest
from GenAIRR._reference_sources import model as m
from tests.support.fixtures.reference_source_models import _raises, lock_dict, manifest_dict, receipt_dict

pytestmark = pytest.mark.behavior

def test_manifest_round_trips_and_binds_the_report(report):
    d = manifest_dict(report)
    manifest = m.BuildManifest.from_dict(d)
    assert manifest.to_dict() == d
    assert canonical_json(m.BuildManifest.from_json(canonical_json(d)).to_dict()) == canonical_json(d)
    assert m.manifest_hash(manifest) == typed_digest(d, "genairr_reference_source_build_manifest_v1")
    m.check_manifest(manifest, report.source_lock, report)


@pytest.mark.parametrize("name", ["", "/report.json", "../report.json", "inputs/../report.json", "inputs\\v.fasta",
                                  "C:report.json", "c:/report.json", "inputs//v.fasta", "./report.json",
                                  "inputs/v.fasta/", "receipts/2026-10-22T09:58:15Z.json",
                                  "manifest.json", "notes.txt", "cache/blobs/sha256/" + "A" * 64,
                                  "cache/blobs/sha256/abc", "inputs/v.fasta\x00", None])
def test_manifest_member_names_are_closed_and_relative(report, name):
    d = manifest_dict(report)
    d["members"][0]["name"] = name
    with _raises("malformed", "/members/0/name"):
        m.BuildManifest.from_dict(d)


@pytest.mark.parametrize("edit,kind,path", [
    (lambda d: d["members"].reverse(), "malformed", "/members"),
    (lambda d: d["members"].insert(0, copy.deepcopy(d["members"][0])), "malformed", "/members"),
    (lambda d: d["members"].pop(_member(d, "report.json")), "malformed", "/members"),
    (lambda d: d["members"][_member(d, "report.json")].update(bytes=-1), "malformed", "/members/"),
    (lambda d: d["members"][_member(d, "report.json")].update(bytes=True), "malformed", "/members/"),
    (lambda d: d["members"][_member(d, "report.json")].update(sha256="sha256:x"), "malformed", "/members/"),
    (lambda d: d["members"][_blob(d)].update(sha256="sha256:" + "0" * 64), "digest_mismatch", "/members/"),
    (lambda d: d["versions"].update(report=2), "unsupported_schema", "/versions/report"),
    (lambda d: d["versions"].update(reference=True), "malformed", "/versions/reference"),
    (lambda d: d["projection"].update(version=9), "unsupported_profile", "/projection"),
    (lambda d: d["builder"].update(version=""), "malformed", "/builder/version"),
    (lambda d: d.update(commitment="sha256:" + "0" * 64), "malformed", ""),
    (lambda d: d["hashes"].pop("report"), "malformed", "/hashes"),
])
def test_manifest_fields_are_strict(report, edit, kind, path):
    d = manifest_dict(report)
    edit(d)
    with _raises(kind) as caught:
        m.BuildManifest.from_dict(d)
    assert caught.error.path.startswith(path), caught.error


def _member(d, name):
    return [i for i, x in enumerate(d["members"]) if x["name"] == name][0]


def _blob(d):
    return [i for i, x in enumerate(d["members"]) if x["name"].startswith("cache/")][0]


def test_check_manifest_against_lock_and_report(report, report_d, inventory, inv):
    lock = report.source_lock
    d = manifest_dict(report)
    d["members"].pop(_blob(d))
    with _raises("missing_asset", "/members"):
        m.check_manifest(m.BuildManifest.from_dict(d), lock, report)
    d = manifest_dict(report)
    d["members"].append({"name": "cache/blobs/sha256/" + "e" * 64, "sha256": "sha256:" + "e" * 64, "bytes": 5})
    d["members"].sort(key=lambda x: x["name"])
    with _raises("identity_mismatch"):
        m.check_manifest(m.BuildManifest.from_dict(d), lock, report)
    d = manifest_dict(report)
    d["members"][_blob(d)]["bytes"] += 1
    with _raises("digest_mismatch"):
        m.check_manifest(m.BuildManifest.from_dict(d), lock, report)
    d = manifest_dict(report, d_fasta=False)
    with _raises("identity_mismatch", "/members"):
        m.check_manifest(m.BuildManifest.from_dict(d), lock, report)
    d = manifest_dict(report)
    d["hashes"]["source_lock"] = "sha256:" + "0" * 64
    with _raises("digest_mismatch", "/hashes/source_lock"):
        m.check_manifest(m.BuildManifest.from_dict(d), lock, report)
    d = manifest_dict(report)
    d["builder"]["version"] = "2.4.0"
    with _raises("identity_mismatch", "/builder"):
        m.check_manifest(m.BuildManifest.from_dict(d), lock, report)


def test_tampered_report_runtime_body_and_report_hash_are_both_caught(report, report_d):
    """Rewrite the stored report's runtime hashes and rehash it: the manifest still binds the
    original report hash; rehashing the manifest's report hash too still leaves its runtime."""
    lock = report.source_lock
    manifest = m.BuildManifest.from_dict(manifest_dict(report))
    tampered_d = copy.deepcopy(report_d)
    tampered_d["runtime"]["execution_hash"] = "sha256:" + "5" * 64
    tampered = m.BuildReport.from_dict(tampered_d)
    with _raises("digest_mismatch", "/hashes/report"):
        m.check_manifest(manifest, lock, tampered)
    rehashed = manifest_dict(report)
    rehashed["hashes"]["report"] = m.report_hash(tampered)
    with _raises("digest_mismatch", "/runtime"):
        m.check_manifest(m.BuildManifest.from_dict(rehashed), lock, tampered)


def test_digest_domains_are_the_six_distinct_graph_domains():
    assert m.DOMAINS == ("genairr_reference_source_lock_v1", "genairr_reference_source_inventory_v1",
                         "genairr_reference_source_recipe_v1", "genairr_reference_source_build_inputs_v1",
                         "genairr_reference_source_build_report_v1", "genairr_reference_source_build_manifest_v1")
    assert len(set(m.DOMAINS)) == 6
    assert not {m.CLAIM_ID_DOMAIN, m.FINDING_ID_DOMAIN, m.MAPPING_ID_DOMAIN} & set(m.DOMAINS)


def test_each_digest_takes_only_its_upstream_models():
    params = {f.__name__: list(inspect.signature(f).parameters)
              for f in (m.lock_hash, m.inventory_hash, m.recipe_hash, m.build_input_hash, m.report_hash,
                        m.manifest_hash)}
    assert params == {"lock_hash": ["lock"], "inventory_hash": ["inventory"], "recipe_hash": ["recipe"],
                      "build_input_hash": ["lock", "inventory", "recipe"], "report_hash": ["report"],
                      "manifest_hash": ["manifest"]}


def test_build_input_identity_is_upstream_of_report_and_manifest(report, report_d, inventory):
    lock, recipe = report.source_lock, report.recipe
    before = m.build_input_hash(lock, inventory, recipe)
    other_d = copy.deepcopy(report_d)
    other_d["runtime"]["artifact_hash"] = "sha256:" + "6" * 64
    other_d["generator"]["version"] = "9.9.9"
    other = m.BuildReport.from_dict(other_d)
    assert m.report_hash(other) != m.report_hash(report)
    assert m.build_input_hash(other.source_lock, other.inventory, other.recipe) == before
    assert report.hashes.build_inputs == before
    inputs = m.build_inputs(lock, inventory, recipe)
    assert set(inputs) == {"source_lock", "inventory", "recipe", "adapter", "projection"}
    assert typed_digest(inputs, "genairr_reference_source_build_inputs_v1") == before
    text = canonical_json(inputs)
    assert m.report_hash(report) not in text and "runtime" not in text and "generator" not in text


def test_lock_inventory_recipe_and_build_input_chain(report, inventory):
    lock, recipe = report.source_lock, report.recipe
    assert inventory.source_lock_hash == m.lock_hash(lock) == recipe.source_lock_hash
    inputs = m.build_inputs(lock, inventory, recipe)
    assert (inputs["source_lock"], inputs["inventory"], inputs["recipe"]) == (
        m.lock_hash(lock), m.inventory_hash(inventory), m.recipe_hash(recipe))


def test_report_and_manifest_hashes_are_downstream(report, report_d):
    manifest = m.BuildManifest.from_dict(manifest_dict(report))
    assert manifest.hashes.report == m.report_hash(report)
    assert m.report_hash(report) not in canonical_json(report.to_dict())
    assert m.manifest_hash(manifest) not in canonical_json(manifest.to_dict())
    changed = manifest_dict(report)
    changed["runtime"]["models_hash"] = "sha256:" + "7" * 64
    assert m.manifest_hash(m.BuildManifest.from_dict(changed)) != m.manifest_hash(manifest)


def test_receipts_enter_no_digest(report, inventory):
    receipt = m.AcquisitionReceipt.from_dict(receipt_dict())
    for function in (m.lock_hash, m.inventory_hash, m.recipe_hash, m.report_hash, m.manifest_hash):
        with pytest.raises(TypeError):
            function(receipt)
    with pytest.raises(TypeError):
        m.build_input_hash(report.source_lock, inventory, receipt)
    for document in (report.source_lock, report.recipe, inventory, report):
        text = canonical_json(document.to_dict())
        assert "event_time_utc" not in text and "genairr_reference_source_receipt" not in text


def test_build_input_digest_refuses_inconsistent_inputs(report, inventory):
    other = m.SourceLock.from_dict(lock_dict(origin="local_file"))
    with _raises("digest_mismatch", "/source_lock_hash"):
        m.build_input_hash(other, inventory, report.recipe)


def test_recipe_or_provenance_changes_move_the_build_input_identity_without_dna_changes(report, inventory):
    """Runtime notes may cite the build-input hash, so a changed recipe can legitimately
    change runtime identity even when every selected sequence is unchanged."""
    d = report.recipe.to_dict()
    d["selection"][0]["reason"] = "Same supplied sequence, reviewed again"
    changed = m.Recipe.from_dict(d)
    assert [e.sequence_sha256 for e in changed.selection] == [e.sequence_sha256 for e in report.recipe.selection]
    assert m.build_input_hash(report.source_lock, inventory, changed) != \
        m.build_input_hash(report.source_lock, inventory, report.recipe)
