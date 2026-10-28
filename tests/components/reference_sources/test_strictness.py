"""Strict parsing across every owned format, and the resource budgets it enforces.
"""
from __future__ import annotations
import copy
import json
import pytest
from GenAIRR._artifact_io import canonical_json, strict_loads, typed_digest
from GenAIRR._reference_sources.errors import KINDS, ReferenceSourceError
from GenAIRR._reference_sources import model as m
from tests.support.fixtures.reference_source_models import ADAPTER, CLAIM_ID, FINDING_ID, FIXED, GENERATOR, IGH, INV_RECORDS, Inv, L, ROLE, ROOT, RUNTIME, SELECTED, SOURCES, _RECORDS, _blocks, _capabilities, _claim, _fasta, _find, _finding, _gapped_index, _inventory_dict, _load_json, _nested, _pointer, _raises, _recipe, _records, _reseal, _sorted, _source_set, decision, inv, inventory, lock_dict, manifest_dict, receipt_dict, recipe_dict, report, report_d, report_decisions, report_dict, selection_entry

pytestmark = pytest.mark.behavior

_WRONG = (None, True, -1, 0, 1.5, "", "x", [], {})


def _leaves(value, path=()):
    """Every (path, value) inside a JSON document, containers included."""
    yield path, value
    if isinstance(value, dict):
        for key, item in value.items():
            yield from _leaves(item, path + (key,))
    elif isinstance(value, list):
        for i, item in enumerate(value):
            yield from _leaves(item, path + (i,))


def _replaced(document, path, value):
    d = copy.deepcopy(document)
    if not path:
        return value
    target = d
    for token in path[:-1]:
        target = target[token]
    target[path[-1]] = value
    return d


def _sweep(reader, document, *, stride=1):
    """Replace leaves with wrong-typed values: the reader either accepts the document or
    raises ReferenceSourceError with a JSON-pointer-like path - never another exception."""
    rejected = 0
    for n, (path, original) in enumerate(_leaves(document)):
        if n % stride:
            continue
        for wrong in _WRONG:
            if type(wrong) is type(original) and wrong == original:
                continue
            try:
                reader(_replaced(document, path, wrong))
            except ReferenceSourceError as exc:
                assert exc.kind in KINDS and (exc.path == "" or exc.path.startswith("/")), exc
                rejected += 1
    assert rejected


def test_every_wrong_value_in_a_lock_receipt_or_recipe_is_a_typed_error():
    _sweep(m.SourceLock.from_dict, lock_dict())
    _sweep(m.AcquisitionReceipt.from_dict, receipt_dict())
    d = recipe_dict()
    d["decisions"] = [decision(), decision("withhold_claim"), decision("accept_claim", claims=["claim:" + "c" * 64])]
    _sweep(m.Recipe.from_dict, d)


def test_every_wrong_value_in_an_inventory_report_or_manifest_is_a_typed_error(inventory, report, report_d):
    _sweep(m.SourceInventory.from_dict, inventory.to_dict(), stride=7)
    _sweep(m.BuildManifest.from_dict, manifest_dict(report))
    small = copy.deepcopy(report_d)
    _sweep(m.BuildReport.from_dict, {k: v for k, v in small.items()}, stride=97)


@pytest.mark.parametrize("cls,make", [
    ("SourceLock", lambda r, i: lock_dict()),
    ("AcquisitionReceipt", lambda r, i: receipt_dict()),
    ("Recipe", lambda r, i: recipe_dict()),
    ("SourceInventory", lambda r, i: i.to_dict()),
    ("BuildReport", lambda r, i: r.to_dict()),
    ("BuildManifest", lambda r, i: manifest_dict(r)),
])
def test_every_owned_json_text_is_strict(report, inventory, tmp_path, cls, make):
    reader = getattr(m, cls)
    d = make(report, inventory)
    text = canonical_json(d)
    assert reader.from_json(text) == reader.from_dict(d)
    first_key = '"artifact":"' + d["artifact"] + '"'
    with _raises("malformed"):
        reader.from_json(text.replace(first_key, first_key + "," + first_key, 1))
    with _raises("malformed"):
        reader.from_json(text[:-1] + ',"x":NaN}')
    with _raises("malformed"):
        reader.from_json(text + "[]")
    with _raises("resource_limit"):
        reader.from_json(text[:-1] + ',"x":' + "[" * 70 + "]" * 70 + "}")
    with _raises("malformed"):
        reader.from_json(None)
    path = tmp_path / "document.json"
    path.write_bytes(text.encode("utf-8"))
    assert reader.read(path) == reader.from_dict(d)
    with _raises("io"):
        reader.read(tmp_path / "missing.json")


def test_owned_documents_bound_their_byte_size(monkeypatch):
    from GenAIRR._reference_sources import recipe as recipe_module
    text = canonical_json(recipe_dict())
    monkeypatch.setattr(recipe_module, "MAX_RECIPE_BYTES", len(text.encode("utf-8")) - 1)
    with _raises("resource_limit"):
        m.Recipe.from_json(text)


def test_containers_are_bounded_by_named_limits(monkeypatch):
    from GenAIRR._reference_sources import recipe as recipe_module
    monkeypatch.setattr(recipe_module, "MAX_SELECTION", 4)
    with _raises("resource_limit", "/selection"):
        m.Recipe.from_dict(recipe_dict())
    monkeypatch.setattr(recipe_module, "MAX_REFS", 1)
    d = recipe_dict(selected=SELECTED[:4])
    d["decisions"] = [decision(findings=[FINDING_ID, "finding:" + "d" * 64])]
    with _raises("resource_limit", "/decisions/0/findings"):
        m.Recipe.from_dict(d)
    d = lock_dict()
    d["assets"] = d["assets"] * 5
    with _raises("resource_limit", "/assets"):
        m.SourceLock.from_dict(d)


def test_strings_are_bounded_by_named_limits():
    d = recipe_dict()
    d["reference"]["name"] = "n" * (m.MAX_NAME + 1)
    with _raises("resource_limit", "/reference/name"):
        m.Recipe.from_dict(d)
    d = recipe_dict()
    d["selection"][0]["label"] = "IGHV" + "1" * 300
    with _raises("malformed", "/selection/0/label"):
        m.Recipe.from_dict(d)
    with _raises("resource_limit", "/value"):
        m.Observation(locator=L("airr", "/Info/title"), record=None, presence="present",
                      value="x" * (m.MAX_VALUE_STRING + 1), interpretation="preserved")


def test_null_and_missing_are_distinct_in_generated_formats(inventory):
    d = _inventory_dict(inventory)
    d["records"][0]["provider_ref"] = None
    assert m.SourceInventory.from_dict(d).records[0].provider_ref is None
    del d["records"][0]["provider_ref"]
    with _raises("malformed", "/records/0"):
        m.SourceInventory.from_dict(d)
    d = _inventory_dict(inventory)
    (index,) = [i for i, o in enumerate(d["observations"]) if o["record"] is None and o["locator"]["asset"] == "airr"]
    del d["observations"][index]["record"]
    with _raises("malformed", f"/observations/{index}"):
        m.SourceInventory.from_dict(d)


def test_independent_builds_are_byte_identical(inventory):
    again = Inv().build()
    assert canonical_json(again.to_dict()) == canonical_json(inventory.to_dict())
    assert m.inventory_hash(again) == m.inventory_hash(inventory)


def test_report_and_inventory_views_are_copies(report, inventory):
    d = report.to_dict()
    d["inventory"]["records"].clear()
    d["records"].clear()
    assert len(report.to_dict()["records"]) == len(inventory.records)
    (extension,) = [c for c in inventory.claims if c.predicate == "extension"]
    copied = extension.value.get()
    copied["length"] = 99
    assert extension.value.get()["length"] == 1
    d["inventory"]["claims"][0]["status"] = "candidate"
    assert report.to_dict()["inventory"]["claims"][0]["status"] == inventory.claims[0].status
    with pytest.raises(AttributeError):
        report.records = ()
    with pytest.raises(AttributeError):
        inventory.claims[0].status = "accepted"


def test_resource_budget_is_shared_by_every_format():
    assert (m.MAX_ASSET_BYTES, m.MAX_SOURCE_PAYLOAD_BYTES, m.MAX_ASSETS) == (16 * 1024 * 1024, 128 * 1024 * 1024, 16)
    assert (m.MAX_SOURCE_RECORDS, m.MAX_JSON_DEPTH, m.MAX_SEQUENCE) == (100_000, 64, 65_535)
    from GenAIRR._reference_sources import inventory, lock, receipt
    assert lock.MAX_ASSET_BYTES is m.MAX_ASSET_BYTES and receipt.MAX_ASSET_BYTES is m.MAX_ASSET_BYTES
    assert inventory.MAX_RECORDS == m.MAX_SOURCE_RECORDS
    loc = L("airr", "/Info/title")
    m.Observation(locator=loc, record=None, presence="present", value=_nested(m.MAX_EMBEDDED_VALUE_DEPTH),
                  interpretation="preserved")


def test_asset_and_aggregate_payload_sizes_are_bounded(monkeypatch):
    d = lock_dict()
    d["assets"][0]["bytes"] = m.MAX_ASSET_BYTES
    m.SourceLock.from_dict(_reseal(d))
    d["assets"][0]["bytes"] = m.MAX_ASSET_BYTES + 1
    with _raises("resource_limit", "/assets/0/bytes") as caught:
        m.SourceLock.from_dict(_reseal(d))
    assert str(m.MAX_ASSET_BYTES) in caught.error.message and str(m.MAX_ASSET_BYTES + 1) in caught.error.message
    from GenAIRR._reference_sources import lock as lock_module
    d = lock_dict()
    total = sum(a["bytes"] for a in d["assets"])
    monkeypatch.setattr(lock_module, "MAX_SOURCE_PAYLOAD_BYTES", total - 1)
    with _raises("resource_limit", "/assets") as caught:
        m.SourceLock.from_dict(d)
    assert str(total) in caught.error.message and str(total - 1) in caught.error.message
    receipt = receipt_dict()
    receipt["bytes"] = m.MAX_ASSET_BYTES + 1
    with _raises("resource_limit", "/bytes"):
        m.AcquisitionReceipt.from_dict(receipt)


def test_profile_scope_lists_the_validated_sets():
    (scope,) = m.PROFILE_SCOPES.values()
    assert list(m.PROFILE_SCOPES) == [("ogrdb_api_v2_human_ig_observed_v1", 1)]
    assert [(t, s, r, tuple(a[0] for a in assets)) for t, s, r, assets in scope] == [
        (9606, "9606.IGH_VDJ", 10, ("airr", "airr_ex", "gapped", "gapped_ex")),
        (9606, "9606.IGKappa_VJ", 5, ("airr", "airr_ex", "gapped", "gapped_ex")),
        (9606, "9606.IGLambda_VJ", 4, ("airr", "airr_ex", "gapped", "gapped_ex"))]
    for set_id in (IGH, "9606.IGKappa_VJ", "9606.IGLambda_VJ"):
        m.check_profile_scope(m.SourceLock.from_dict(lock_dict(set_id)))


def test_profile_scope_pins_exactly_the_pinned_fixture_bytes():
    """The registry literals and the provenance manifest of the pinned exports cannot diverge."""
    (scope,) = m.PROFILE_SCOPES.values()
    expected = []
    for set_id in (IGH, "9606.IGKappa_VJ", "9606.IGLambda_VJ"):
        source = _source_set(set_id)
        by_format = {a["format"]: a for a in source["assets"]}
        expected.append((int(source["species"]["id"].split(":")[1]), set_id, source["release"],
                         tuple((f, "sha256:" + by_format[f]["sha256"], by_format[f]["bytes"])
                               for f in ("airr", "airr_ex", "gapped", "gapped_ex"))))
    assert list(scope) == expected


@pytest.mark.parametrize("edit", [lambda a: a.update(sha256="sha256:" + "e" * 64),
                                  lambda a: a.update(bytes=a["bytes"] + 1)])
def test_drifted_bytes_for_a_validated_release_are_out_of_scope(edit):
    d = lock_dict()
    edit([a for a in d["assets"] if a["format"] == "gapped_ex"][0])
    lock = m.SourceLock.from_dict(_reseal(d))
    with _raises("unsupported_profile", "") as caught:
        m.check_profile_scope(lock)
    assert "gapped_ex" in caught.error.message
    built = Inv(lock)
    with _raises("unsupported_profile", ""):
        m.check_build_inputs(lock, built.build(claims=[], mappings=[]), m.Recipe.from_dict(recipe_dict(lock.to_dict())))
    with _raises("unsupported_profile") as caught:
        m.check_inventory(built.build(), lock)
    assert caught.error.path.endswith("/status")


def test_a_lock_missing_some_validated_exports_stays_in_scope():
    m.check_profile_scope(m.SourceLock.from_dict(lock_dict(formats=("airr", "gapped"))))


def _out_of_scope_lock(**edits):
    d = lock_dict()
    d.update(edits)
    return m.SourceLock.from_dict(_reseal(d))


@pytest.mark.parametrize("edits", [{"release": 9}, {"set_id": "9606.IGH_VDJ_test"},
                                   {"identity": {"taxon": {"scientific_name": "Mus musculus", "taxonomy_id": 10090},
                                                 "locus": "IGH", "chain_type": "vdj"}}])
def test_locks_outside_the_profile_scope_are_valid_but_cannot_build(edits):
    lock = _out_of_scope_lock(**edits)
    with _raises("unsupported_profile", ""):
        m.check_profile_scope(lock)
    inventory = Inv(lock).build(claims=[], mappings=[])
    d = recipe_dict(lock.to_dict())
    d["reference"]["taxon"] = lock.to_dict()["identity"]["taxon"]
    with _raises("unsupported_profile", ""):
        m.check_build_inputs(lock, inventory, m.Recipe.from_dict(d))


def test_accepted_coordinate_claims_need_a_lock_inside_the_profile_scope():
    lock = _out_of_scope_lock(release=9)
    built = Inv(lock)
    with _raises("unsupported_profile") as caught:
        m.check_inventory(built.build(), lock)
    assert caught.error.path.startswith("/claims/") and caught.error.path.endswith("/status")
    coordinate = ("anchor", "coding_phase", "phase_independence", "v_subregions", "coding_interval")
    kept = [c for c in built.claims if c.status != "accepted" or c.predicate not in coordinate]
    m.check_inventory(built.build(claims=kept, mappings=[]), lock)
    assert m.COORDINATE_PREDICATES == coordinate


def test_out_of_scope_coordinate_claims_are_only_candidates_or_unresolved():
    lock = _out_of_scope_lock(release=9)
    built = Inv(lock)
    vcand = built.cand["airr", "IGHV1-58*03"]
    cdr3 = built.obs_at(built.rec["airr", "IGHV1-58*03"], "v_gene_delineations/0/cdr3_start")
    others = [c for c in built.claims if c.predicate not in m.COORDINATE_PREDICATES]
    conflicting = [built.claim(vcand, "anchor", built.anchor, [cdr3], system="coding_sequence", status="conflicting"),
                   built.claim(vcand, "anchor", built.anchor + 3, [L("airr", vcand.pointer)], system="coding_sequence",
                               status="conflicting")]
    with _raises("unsupported_profile") as caught:
        m.check_inventory(built.build(claims=others + conflicting, mappings=[]), lock)
    assert caught.error.path.endswith("/status") and "conflicting" in caught.error.message
    allowed = [built.claim(vcand, "anchor", built.anchor, [cdr3], system="coding_sequence", status="candidate",
                           support="none"),
               built.claim(vcand, "coding_phase", None, [cdr3], system="coding_sequence", status="unresolved",
                           support="none")]
    m.check_inventory(built.build(claims=others + allowed, mappings=[]), lock)


def test_candidate_mention_and_mapping_caps_follow_the_record_budget():
    from GenAIRR._reference_sources import inventory as inventory_module
    assert inventory_module.MAX_RECORDS == m.MAX_SOURCE_RECORDS == 100_000
    assert inventory_module.MAX_CANDIDATES == inventory_module.MAX_MENTIONED == m.MAX_SOURCE_RECORDS
    assert inventory_module.MAX_MAPPINGS == 2 * inventory_module.MAX_CANDIDATES


def test_an_inventory_beyond_65536_records_with_one_candidate_each_is_accepted():
    lock = m.SourceLock.from_dict(lock_dict())
    count = 65_537
    records, candidates = [], []
    for i in range(count):
        loc = L("airr", _pointer(i))
        records.append(m.SourceRecord(locator=loc, status="parsed", problem=None, label=f"IGHV9-{i}*01",
                                      provider_id=None, provider_ref=None, segment="V", taxonomy_id=9606,
                                      locus="IGH"))
        sequence = format(i, "b").replace("0", "A").replace("1", "C") + "GT"
        candidates.append(m.SequenceCandidate(locator=L("airr", loc.pointer + "/coding_sequence"), record=loc,
                                              sequence=sequence, sha256=m.sequence_digest(sequence)))
    inventory = m.SourceInventory(source_lock_hash=lock.commitment, adapter=lock.adapter, records=tuple(records),
                                  candidates=tuple(candidates), mentioned_identities=(), observations=(),
                                  support=(), claims=(), mappings=(), findings=())
    assert len(inventory.candidates) == count
    assert inventory.export_counts("airr")["sequence_candidates"]["V"] == count
