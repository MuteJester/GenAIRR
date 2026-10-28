"""Typed errors, the source lock, and the receipt an acquisition records.
"""
from __future__ import annotations
import copy
import json
import pytest
from GenAIRR._artifact_io import canonical_json, strict_loads, typed_digest
from GenAIRR._reference_sources.errors import KINDS, ReferenceSourceError
from GenAIRR._reference_sources import model as m
from tests.support.fixtures.reference_source_models import ADAPTER, IGH, ROLE, ROOT, SOURCES, _load_json, _raises, _reseal, _source_set, lock_dict, receipt_dict

def test_error_kinds_are_the_stable_vocabulary():
    assert KINDS == ("io", "network", "resource_limit", "malformed", "unsupported_schema",
                     "unsupported_profile", "identity_mismatch", "digest_mismatch", "missing_asset",
                     "incomplete_dossier", "unresolved_selection", "stale_decision", "projection",
                     "capability")


def test_error_carries_kind_path_and_message():
    error = ReferenceSourceError("malformed", "/assets/0/sha256", "expected sha256:<64 hex>")
    assert (error.kind, error.path, error.message) == ("malformed", "/assets/0/sha256",
                                                       "expected sha256:<64 hex>")
    assert str(error) == "/assets/0/sha256: expected sha256:<64 hex>"
    assert isinstance(error, ValueError)
    with pytest.raises(ValueError, match="unknown"):
        ReferenceSourceError("integrity", "", "not a stable kind")


@pytest.mark.parametrize("set_id", [IGH, "9606.IGKappa_VJ", "9606.IGLambda_VJ"])
def test_lock_round_trips_for_every_pinned_set(set_id):
    d = lock_dict(set_id)
    lock = m.SourceLock.from_dict(d)
    assert lock.to_dict() == d
    assert m.SourceLock.from_dict(lock.to_dict()) == lock
    text = canonical_json(lock.to_dict())
    assert m.SourceLock.from_json(text) == lock
    assert canonical_json(m.SourceLock.from_json(text).to_dict()) == text
    assert lock.identity.chain_type == ("vdj" if set_id == IGH else "vj")


def test_lock_commitment_is_the_typed_digest_of_every_preceding_field():
    d = lock_dict()
    lock = m.SourceLock.from_dict(d)
    body = {k: v for k, v in d.items() if k != "commitment"}
    assert lock.commitment == typed_digest(body, "genairr_reference_source_lock_v1")
    assert m.lock_hash(lock) == lock.commitment
    assert lock.body() == body


def test_lock_seal_computes_the_commitment():
    d = lock_dict()
    lock = m.SourceLock.from_dict(d)
    sealed = m.SourceLock.seal(provider=lock.provider, set_id=lock.set_id, release=lock.release,
                               identity=lock.identity, adapter=lock.adapter, assets=lock.assets,
                               missing_assets=lock.missing_assets, terms=lock.terms)
    assert sealed == lock


def test_lock_with_edited_body_and_stale_commitment_is_a_digest_mismatch():
    d = lock_dict()
    d["release"] = 9
    with _raises("digest_mismatch", "/commitment"):
        m.SourceLock.from_dict(d)


def test_lock_asset_digests_and_sizes_come_from_the_pinned_fixtures():
    lock = m.SourceLock.from_dict(lock_dict())
    by_format = {a["format"]: a for a in _source_set(IGH)["assets"]}
    for asset in lock.assets:
        assert asset.sha256 == "sha256:" + by_format[asset.format]["sha256"]
        assert asset.bytes == by_format[asset.format]["bytes"]
    assert [a.asset_id for a in lock.assets] == ["airr", "airr_ex", "gapped", "gapped_ex"]


@pytest.mark.parametrize("key,value", [("retrieved_at_utc", "2026-10-22T09:58:15Z"), ("http_status", 200),
                                       ("final_url", "https://ogrdb.airr-community.org/x"),
                                       ("content_disposition", "attachment; filename=Homo sapiens_IGH_VDJ_rev_10.json"),
                                       ("file_name", "airr.json"), ("path", "/home/user/airr.json")])
def test_lock_assets_cannot_carry_acquisition_or_filesystem_details(key, value):
    d = lock_dict()
    d["assets"][0][key] = value
    with _raises("malformed", "/assets/0") as caught:
        m.SourceLock.from_dict(_reseal(d))
    assert key in caught.error.message


def test_moving_or_refetching_identical_bytes_keeps_the_lock_identity():
    first = m.SourceLock.from_dict(lock_dict(origin="local_file"))
    again = m.SourceLock.seal(provider="ogrdb", set_id=first.set_id, release=first.release, identity=first.identity,
                              adapter=first.adapter, assets=tuple(m.Asset(**a.__dict__) for a in first.assets),
                              missing_assets=first.missing_assets, terms=first.terms)
    assert again.commitment == first.commitment


def test_local_import_and_provider_download_have_distinct_lock_identities():
    provider = m.SourceLock.from_dict(lock_dict(origin="provider_download"))
    local = m.SourceLock.from_dict(lock_dict(origin="local_file"))
    assert [a.sha256 for a in provider.assets] == [a.sha256 for a in local.assets]
    assert provider.commitment != local.commitment
    assert all(a.declared_url is None for a in local.assets)


@pytest.mark.parametrize("key", ["artifact", "schema_version", "provider", "set_id", "release", "identity",
                                 "adapter", "assets", "missing_assets", "terms", "commitment"])
def test_lock_rejects_every_missing_top_level_key(key):
    d = lock_dict()
    del d[key]
    with _raises("malformed"):
        m.SourceLock.from_dict(d)


def test_lock_rejects_unknown_keys_even_when_rehashed():
    d = lock_dict()
    d["retrieved_at"] = "2026-10-22T09:58:15Z"
    with _raises("malformed", "") as caught:
        m.SourceLock.from_dict(_reseal(d))
    assert "retrieved_at" in caught.error.message


def test_lock_schema_discriminator_and_version():
    d = lock_dict()
    d["artifact"] = "genairr_reference_recipe"
    with _raises("malformed", "/artifact"):
        m.SourceLock.from_dict(d)
    for bad_version, kind in ((2, "unsupported_schema"), (0, "unsupported_schema"),
                              (True, "malformed"), ("1", "malformed"), (1.0, "malformed"),
                              (None, "malformed")):
        d = lock_dict()
        d["schema_version"] = bad_version
        with _raises(kind, "/schema_version"):
            m.SourceLock.from_dict(d)


def test_lock_version_check_precedes_key_closure():
    d = lock_dict()
    d["schema_version"] = 2
    d["future_field"] = {}
    with _raises("unsupported_schema"):
        m.SourceLock.from_dict(d)


@pytest.mark.parametrize("release", [0, -1, True, False, "10", 10.0, None, 2 ** 31])
def test_lock_release_is_a_positive_integer(release):
    d = lock_dict()
    d["release"] = release
    with _raises("malformed", "/release"):
        m.SourceLock.from_dict(_reseal(d))


@pytest.mark.parametrize("edit,kind,path", [
    (lambda d: d.update(provider="imgt"), "unsupported_profile", "/provider"),
    (lambda d: d["adapter"].update(id="ogrdb_other"), "unsupported_profile", "/adapter"),
    (lambda d: d["adapter"].update(version=2), "unsupported_profile", "/adapter"),
    (lambda d: d["adapter"]["coordinate_profile"].update(version=2), "malformed", "/adapter/coordinate_profile/id"),
    (lambda d: d["adapter"]["coordinate_profile"].update(id="imgt"), "malformed", "/adapter/coordinate_profile/id"),
    (lambda d: d["adapter"]["coordinate_profile"].update(id="ogrdb_api_v2"), "malformed",
     "/adapter/coordinate_profile/id"),
    (lambda d: d["adapter"]["coordinate_profile"].update(id="imgt_v1"), "unsupported_profile", "/adapter"),
    (lambda d: d["adapter"]["coordinate_profile"].update(id="ogrdb_api_v2_human_ig_observed_v2", version=2),
     "unsupported_profile", "/adapter"),
    (lambda d: d["adapter"].update(version=True), "malformed", "/adapter/version"),
    (lambda d: d["adapter"]["coordinate_profile"].update(version=0), "malformed",
     "/adapter/coordinate_profile/version"),
])
def test_lock_rejects_unknown_semantic_profiles(edit, kind, path):
    d = lock_dict()
    edit(d)
    with _raises(kind, path):
        m.SourceLock.from_dict(_reseal(d))


@pytest.mark.parametrize("edit,path", [
    (lambda d: d.update(set_id=""), "/set_id"),
    (lambda d: d.update(set_id="../9606.IGH_VDJ"), "/set_id"),
    (lambda d: d.update(set_id="9606 IGH"), "/set_id"),
    (lambda d: d["identity"].update(locus="IGX"), "/identity/locus"),
    (lambda d: d["identity"].update(chain_type="vj"), "/identity/chain_type"),
    (lambda d: d["identity"]["taxon"].update(taxonomy_id=None), "/identity/taxon/taxonomy_id"),
    (lambda d: d["identity"]["taxon"].update(taxonomy_id=0), "/identity/taxon/taxonomy_id"),
    (lambda d: d["identity"]["taxon"].update(taxonomy_id=True), "/identity/taxon/taxonomy_id"),
    (lambda d: d["identity"]["taxon"].update(scientific_name="  "), "/identity/taxon/scientific_name"),
    (lambda d: d["identity"]["taxon"].pop("taxonomy_id"), "/identity/taxon"),
])
def test_lock_identity_is_validated(edit, path):
    d = lock_dict()
    edit(d)
    with _raises("malformed", path):
        m.SourceLock.from_dict(_reseal(d))


@pytest.mark.parametrize("asset_id", ["", "/airr", "../airr", "air/r", "a\\b", "C:airr", "c:",
                                      "Airr", "airr ", "air\x00r", 7, None])
def test_lock_rejects_path_forms_in_asset_names(asset_id):
    d = lock_dict()
    d["assets"][0]["asset_id"] = asset_id
    with _raises("malformed", "/assets/0/asset_id"):
        m.SourceLock.from_dict(_reseal(d))


def test_lock_asset_ids_are_generated_from_the_format():
    d = lock_dict()
    d["assets"][0]["asset_id"] = "homo_sapiens_igh"
    with _raises("malformed", "/assets/0/asset_id"):
        m.SourceLock.from_dict(_reseal(d))


@pytest.mark.parametrize("edit,kind,path", [
    (lambda a: a.update(sha256="6889fea1e2da71ff2bde50f8b5e055613319408977e913631c8e6d04f53d44b1"),
     "malformed", "/assets/0/sha256"),
    (lambda a: a.update(sha256="sha256:" + "A" * 64), "malformed", "/assets/0/sha256"),
    (lambda a: a.update(bytes=0), "malformed", "/assets/0/bytes"),
    (lambda a: a.update(bytes=-5), "malformed", "/assets/0/bytes"),
    (lambda a: a.update(bytes=True), "malformed", "/assets/0/bytes"),
    (lambda a: a.update(bytes=1272259.0), "malformed", "/assets/0/bytes"),
    (lambda a: a.update(bytes=m.MAX_ASSET_BYTES + 1), "resource_limit", "/assets/0/bytes"),
    (lambda a: a.update(role="gapped_alignment"), "malformed", "/assets/0/role"),
    (lambda a: a.update(origin_kind="cache"), "malformed", "/assets/0/origin_kind"),
    (lambda a: a.update(declared_url=None), "malformed", "/assets/0/declared_url"),
    (lambda a: a.update(declared_url="http://ogrdb.airr-community.org/x"), "malformed", "/assets/0/declared_url"),
    (lambda a: a.update(declared_url="https://user:pw@ogrdb.airr-community.org/x"), "malformed",
     "/assets/0/declared_url"),
    (lambda a: a.update(declared_url="https://ogrdb.airr-community.org/a b"), "malformed",
     "/assets/0/declared_url"),
    (lambda a: a.update(declared_url="https://x/" + "a" * m.MAX_URL), "resource_limit", "/assets/0/declared_url"),
    (lambda a: a.pop("declared_url"), "malformed", "/assets/0"),
])
def test_lock_asset_fields_are_strict(edit, kind, path):
    d = lock_dict()
    edit(d["assets"][0])
    with _raises(kind, path):
        m.SourceLock.from_dict(_reseal(d))


def test_lock_local_file_may_declare_a_url_but_a_download_must():
    d = lock_dict(origin="local_file")
    d["assets"][0]["declared_url"] = _source_set(IGH)["assets"][0]["url"]
    lock = m.SourceLock.from_dict(_reseal(d))
    assert lock.assets[0].declared_url.startswith("https://ogrdb.airr-community.org/")


def test_lock_assets_are_in_canonical_format_order_without_duplicates():
    d = lock_dict()
    d["assets"].reverse()
    with _raises("malformed", "/assets"):
        m.SourceLock.from_dict(_reseal(d))
    d = lock_dict()
    d["assets"].insert(1, copy.deepcopy(d["assets"][0]))
    with _raises("malformed", "/assets"):
        m.SourceLock.from_dict(_reseal(d))


def test_lock_needs_a_record_export_and_accounts_for_every_evidence_format():
    d = lock_dict(formats=("gapped", "gapped_ex"))
    with _raises("malformed", "/assets"):
        m.SourceLock.from_dict(d)
    d = lock_dict(formats=("airr", "gapped", "gapped_ex"))
    lock = m.SourceLock.from_dict(d)
    assert [x.format for x in lock.missing_assets] == ["airr_ex"]
    d["missing_assets"] = []
    with _raises("malformed", "/missing_assets"):
        m.SourceLock.from_dict(_reseal(d))


def test_lock_missing_assets_are_explicit_and_never_also_present():
    d = lock_dict(formats=("airr", "gapped"))
    lock = m.SourceLock.from_dict(d)
    assert [(x.format, x.role, x.reason) for x in lock.missing_assets] == [
        ("airr_ex", "record_export", "not_requested"), ("gapped_ex", "gapped_alignment", "not_requested")]
    d["missing_assets"].append({"format": "airr", "role": "record_export", "reason": "not_supplied",
                                "detail": None})
    with _raises("malformed", "/missing_assets"):
        m.SourceLock.from_dict(_reseal(d))
    d = lock_dict(formats=("airr", "gapped"))
    d["missing_assets"][0]["reason"] = "download_failed"
    with _raises("malformed", "/missing_assets/0/reason"):
        m.SourceLock.from_dict(_reseal(d))
    d = lock_dict(formats=("airr", "gapped"))
    d["missing_assets"].reverse()
    with _raises("malformed", "/missing_assets"):
        m.SourceLock.from_dict(_reseal(d))


def test_lock_terms_are_recorded_from_an_asset_or_explicitly_unresolved():
    d = lock_dict()
    d["assets"].append({"asset_id": "licence_statement", "role": "data_terms", "format": "licence_statement",
                        "sha256": "sha256:" + SOURCES["data_licence"]["statement_page_sha256"], "bytes": 4096,
                        "declared_url": SOURCES["data_licence"]["statement_url"],
                        "origin_kind": "provider_download"})
    d["terms"] = {"status": "recorded", "assets": ["licence_statement"], "licence": "CC0 1.0 Universal",
                  "citation": SOURCES["citation"]["requested_by_provider"], "reason": None}
    lock = m.SourceLock.from_dict(_reseal(d))
    assert lock.terms.licence == "CC0 1.0 Universal"
    for edit in (lambda t: t.update(assets=[]), lambda t: t.update(licence=None),
                 lambda t: t.update(reason="not captured"), lambda t: t.update(assets=["airr"]),
                 lambda t: t.update(status="unresolved")):
        bad = copy.deepcopy(d)
        edit(bad["terms"])
        with _raises("malformed"):
            m.SourceLock.from_dict(_reseal(bad))
    bad = copy.deepcopy(d)
    bad["terms"] = lock_dict()["terms"]
    with _raises("malformed", "/terms/assets"):
        m.SourceLock.from_dict(_reseal(bad))


def test_lock_json_text_is_strict():
    text = canonical_json(lock_dict())
    duplicated = text.replace('"provider":"ogrdb"', '"provider":"ogrdb","provider":"ogrdb"')
    with _raises("malformed"):
        m.SourceLock.from_json(duplicated)
    with _raises("malformed"):
        m.SourceLock.from_json(text.replace('"release":10', '"release":NaN'))
    with _raises("malformed"):
        m.SourceLock.from_json(text + "x")
    with _raises("resource_limit"):
        m.SourceLock.from_json(text[:-1] + "," + '"pad":"' + "x" * m.MAX_LOCK_BYTES + '"}')
    with _raises("malformed"):
        m.SourceLock.from_json(b"{}")


def test_lock_from_json_accepts_a_byte_order_mark_and_reads_files(tmp_path):
    text = canonical_json(lock_dict())
    assert m.SourceLock.from_json("﻿" + text).to_dict() == lock_dict()
    path = tmp_path / "source.lock.json"
    path.write_bytes(text.encode("utf-8"))
    assert m.SourceLock.read(path).to_dict() == lock_dict()
    with _raises("io"):
        m.SourceLock.read(tmp_path / "absent.json")
    path.write_bytes(b"\xff\xfe{}")
    with _raises("malformed"):
        m.SourceLock.read(path)


def test_lock_objects_are_immutable():
    lock = m.SourceLock.from_dict(lock_dict())
    with pytest.raises(AttributeError):
        lock.release = 9
    assert isinstance(lock.assets, tuple)
    returned = lock.to_dict()
    returned["assets"].clear()
    assert len(lock.to_dict()["assets"]) == 4


def test_lock_direct_construction_is_validated():
    lock = m.SourceLock.from_dict(lock_dict())
    with _raises("malformed", "/assets"):
        m.SourceLock.seal(provider="ogrdb", set_id=lock.set_id, release=lock.release, identity=lock.identity,
                          adapter=lock.adapter, assets=list(lock.assets), missing_assets=lock.missing_assets,
                          terms=lock.terms)
    with _raises("malformed", "/bytes"):
        m.Asset(asset_id="airr", role="record_export", format="airr", sha256="sha256:" + "0" * 64,
                bytes=True, declared_url=None, origin_kind="local_file")


@pytest.mark.parametrize("origin", ["provider_download", "local_file", "verified_cache"])
def test_receipt_round_trips(origin):
    d = receipt_dict(origin=origin)
    receipt = m.AcquisitionReceipt.from_dict(d)
    assert receipt.to_dict() == d
    assert m.AcquisitionReceipt.from_json(canonical_json(d)) == receipt
    assert receipt.origin == origin


@pytest.mark.parametrize("edit,path", [
    (lambda d: d["headers"].update(authorization="Bearer secret"), "/headers"),
    (lambda d: d["headers"].update(cookie="session=1"), "/headers"),
    (lambda d: d["headers"].update({"Content-Type": "text/html"}), "/headers"),
    (lambda d: d.update(requested_url="https://user:pw@ogrdb.airr-community.org/x"), "/requested_url"),
    (lambda d: d.update(http_status=404), "/http_status"),
    (lambda d: d.update(http_status=True), "/http_status"),
    (lambda d: d.update(final_url=None), "/final_url"),
    (lambda d: d.update(event_time_utc="2026-10-22 09:58:15"), "/event_time_utc"),
    (lambda d: d.update(event_time_utc="2026-10-22T09:58:15+02:00"), "/event_time_utc"),
    (lambda d: d.update(redirects=["https://ogrdb.airr-community.org/elsewhere"]), "/redirects"),
    (lambda d: d.update(bytes=0), "/bytes"),
    (lambda d: d.update(origin="cache"), "/origin"),
])
def test_receipt_rejects_credentials_and_inconsistent_events(edit, path):
    d = receipt_dict()
    edit(d)
    with _raises("malformed", path):
        m.AcquisitionReceipt.from_dict(d)


def test_receipt_records_actual_redirects_ending_at_the_final_url():
    d = receipt_dict()
    d["requested_url"] = "https://ogrdb.airr-community.org/api_v2/germline/set/9606.IGH_VDJ/latest/airr"
    d["redirects"] = [d["final_url"]]
    assert m.AcquisitionReceipt.from_dict(d).redirects == (d["final_url"],)
    d["redirects"] = []
    with _raises("malformed", "/redirects"):
        m.AcquisitionReceipt.from_dict(d)


def test_local_and_cached_receipts_carry_no_network_fields():
    for origin in ("local_file", "verified_cache"):
        for key, value in (("requested_url", "https://ogrdb.airr-community.org/x"), ("http_status", 200),
                           ("headers", {"content-type": "application/json"})):
            d = receipt_dict(origin=origin)
            d[key] = value
            with _raises("malformed", "/" + key):
                m.AcquisitionReceipt.from_dict(d)


def test_receipts_with_different_times_describe_the_same_lock():
    first, second = receipt_dict(), receipt_dict()
    second["event_time_utc"] = "2026-10-22T10:10:32.986529+00:00"
    assert m.AcquisitionReceipt.from_dict(first) != m.AcquisitionReceipt.from_dict(second)
    lock = m.SourceLock.from_dict(lock_dict())
    (asset,) = [a for a in lock.assets if a.format == "airr"]
    for receipt in (first, second):
        assert (receipt["sha256"], receipt["bytes"]) == (asset.sha256, asset.bytes)
    assert "event_time_utc" not in canonical_json(lock.to_dict())
