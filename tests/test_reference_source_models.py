"""Closed data contracts for importing published germline sources.

Covers the source lock, acquisition receipt, build recipe, normalized source inventory
(records, sequence candidates, mentioned identities, observations, claims, mappings and
findings), build report and dossier manifest, and the typed digest graph that links
them.  Realistic examples use the pinned OGRDB fixtures in ``tests/data/reference_sources``;
provider payloads are parsed only with the strict JSON loader.
"""
from __future__ import annotations

import copy
import inspect
import json
from pathlib import Path

import pytest

from GenAIRR._artifact_io import canonical_json, strict_loads, typed_digest
from GenAIRR._reference_sources.errors import KINDS, ReferenceSourceError
from GenAIRR._reference_sources import model as m

ROOT = Path(__file__).resolve().parent / "data" / "reference_sources"
IGH = "9606.IGH_VDJ"
ROLE = {"airr": "record_export", "airr_ex": "record_export", "gapped": "gapped_alignment",
        "gapped_ex": "gapped_alignment"}
ADAPTER = {"id": "ogrdb_germline_set", "version": 1,
           "coordinate_profile": {"id": "ogrdb_api_v2_human_ig_observed_v1", "version": 1}}


def _load_json(path: Path, max_bytes: int = 1 << 23):
    return strict_loads(path.read_bytes().decode("utf-8"), max_bytes=max_bytes)


SOURCES = _load_json(ROOT / "sources.json")


def _source_set(set_id):
    (found,) = [s for s in SOURCES["sets"] if s["set_id"] == set_id]
    return found


def _raises(kind, path=None):
    """Context manager: a ReferenceSourceError of ``kind`` (and ``path`` if given)."""
    class _Check:
        def __enter__(self):
            return self

        def __exit__(self, exc_type, exc, tb):
            assert exc_type is not None, f"expected ReferenceSourceError({kind!r})"
            assert issubclass(exc_type, ReferenceSourceError), f"unexpected {exc_type.__name__}: {exc}"
            assert exc.kind == kind, f"expected kind {kind!r}, got {exc.kind!r}: {exc}"
            if path is not None:
                assert exc.path == path, f"expected path {path!r}, got {exc.path!r}: {exc}"
            self.error = exc
            return True
    return _Check()


def _reseal(lock: dict) -> dict:
    """Recompute a lock dict's commitment after an edit (a rehashed document)."""
    body = {k: v for k, v in lock.items() if k != "commitment"}
    lock["commitment"] = typed_digest(body, "genairr_reference_source_lock_v1")
    return lock


def lock_dict(set_id=IGH, formats=("airr", "airr_ex", "gapped", "gapped_ex"),
              origin="provider_download") -> dict:
    source = _source_set(set_id)
    assets = []
    for asset in source["assets"]:
        if asset["format"] in formats:
            fmt = asset["format"]
            assets.append({"asset_id": fmt, "role": ROLE[fmt], "format": fmt,
                           "sha256": "sha256:" + asset["sha256"], "bytes": asset["bytes"],
                           "declared_url": asset["url"] if origin == "provider_download" else None,
                           "origin_kind": origin})
    missing = [{"format": f, "role": ROLE[f], "reason": "not_requested", "detail": None}
               for f in ("airr", "airr_ex", "gapped", "gapped_ex") if f not in formats]
    locus = source["locus"]
    lock = {
        "artifact": "genairr_reference_source_lock", "schema_version": 1, "provider": "ogrdb",
        "set_id": source["set_id"], "release": source["release"],
        "identity": {"taxon": {"scientific_name": source["species"]["label"],
                               "taxonomy_id": int(source["species"]["id"].split(":")[1])},
                     "locus": locus, "chain_type": "vdj" if locus == "IGH" else "vj"},
        "adapter": copy.deepcopy(ADAPTER), "assets": assets, "missing_assets": missing,
        "terms": {"status": "unresolved", "assets": [], "licence": None, "citation": None,
                  "reason": "data terms statement not captured with these exports"},
    }
    return _reseal(lock)


# ── errors ───────────────────────────────────────────────────────────


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


# ── source lock ──────────────────────────────────────────────────────


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


def test_lock_identity_excludes_receipt_time_paths_and_headers():
    d = lock_dict()
    text = canonical_json(d)
    for leaked in ("retrieved_at", "content_disposition", "Homo sapiens_IGH_VDJ_rev_10", "http_status",
                   "final_url", "fixture"):
        assert leaked not in text


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


@pytest.mark.parametrize("asset_id", ["", "/airr", "../airr", "air/r", "a\\b", "C:airr", "c:", "a" * 65,
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


# ── acquisition receipt ──────────────────────────────────────────────


def receipt_dict(fmt="airr", origin="provider_download") -> dict:
    (asset,) = [a for a in _source_set(IGH)["assets"] if a["format"] == fmt]
    provider = origin == "provider_download"
    return {
        "artifact": "genairr_reference_source_receipt", "schema_version": 1,
        "event_time_utc": asset["retrieved_at_utc"], "origin": origin, "format": fmt,
        "sha256": "sha256:" + asset["sha256"], "bytes": asset["bytes"],
        "requested_url": asset["url"] if provider else None,
        "final_url": asset["final_url"] if provider else None,
        "redirects": [],
        "http_status": asset["http_status"] if provider else None,
        "headers": ({"content-type": asset["content_type"], "content-disposition": asset["content_disposition"]}
                    if provider else {}),
    }


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


# ── recipe ───────────────────────────────────────────────────────────

_RECORDS = {}


def _records(fmt, set_id=IGH):
    key = (set_id, fmt)
    if key not in _RECORDS:
        source = _source_set(set_id)
        document = _load_json(ROOT / "ogrdb" / set_id / str(source["release"]) / f"{fmt}.json")
        (germline_set,) = document["GermlineSet"]
        _RECORDS[key] = germline_set["allele_descriptions"]
    return _RECORDS[key]


def _pointer(index):
    return f"/GermlineSet/0/allele_descriptions/{index}"


def _find(label, fmt="airr", set_id=IGH):
    (found,) = [(i, r) for i, r in enumerate(_records(fmt, set_id)) if r["label"] == label]
    return found


def selection_entry(label, fmt="airr", reason="Use this supplied observed sequence for the declared workflow"):
    index, record = _find(label, fmt)
    return {"segment": record["sequence_type"], "label": label,
            "record": {"asset": fmt, "pointer": _pointer(index)},
            "sequence_sha256": m.sequence_digest(record["coding_sequence"]), "reason": reason}


SELECTED = ("IGHV1-58*03", "IGHV1-69*01", "IGHV1-69D*01", "IGHD1-1*01", "IGHJ4*02")


def recipe_dict(lock=None, selected=SELECTED) -> dict:
    lock = lock or lock_dict()
    return {
        "artifact": "genairr_reference_recipe", "schema_version": 1,
        "source_lock_hash": lock["commitment"],
        "adapter": copy.deepcopy(ADAPTER), "projection": {"id": "observed_span", "version": 1},
        "reference": {"name": "human_igh_ogrdb_r10", "version": "1", "taxon": copy.deepcopy(lock["identity"]["taxon"]),
                      "locus": lock["identity"]["locus"], "chain_type": lock["identity"]["chain_type"],
                      "build_date": "2026-10-22"},
        "selection": [selection_entry(label) for label in selected],
        "decisions": [],
        "requires": [{"capability": "observed_span_simulation", "evidence": "assertion"}],
        "defaults_profile": "genairr.unfitted_defaults.v1",
    }


def test_recipe_keys_are_frozen():
    recipe = m.Recipe.from_dict(recipe_dict())
    assert tuple(sorted(recipe.to_dict())) == tuple(sorted((
        "artifact", "schema_version", "source_lock_hash", "adapter", "projection", "reference", "selection",
        "decisions", "requires", "defaults_profile")))
    assert tuple(recipe.to_dict()["adapter"]) == ("id", "version", "coordinate_profile")
    assert tuple(recipe.to_dict()["projection"]) == ("id", "version")
    assert sorted(recipe.to_dict()["selection"][0]) == ["label", "reason", "record", "segment", "sequence_sha256"]


def test_recipe_round_trips_and_is_byte_stable():
    d = recipe_dict()
    recipe = m.Recipe.from_dict(d)
    assert recipe.to_dict() == d
    text = canonical_json(d)
    assert canonical_json(m.Recipe.from_json(text).to_dict()) == text
    assert m.Recipe.from_dict(recipe.to_dict()) == recipe


def test_recipe_reference_matches_the_annotation_sidecar_reference_shape():
    from GenAIRR._reference.model import ReferenceHeader
    recipe = m.Recipe.from_dict(recipe_dict())
    header = ReferenceHeader.from_dict(recipe.reference.to_dict())
    assert (header.locus, header.chain_type, header.taxon.taxonomy_id) == ("IGH", "vdj", 9606)


def test_recipe_digest_is_computed_over_the_complete_recipe_and_never_stored_in_it():
    d = recipe_dict()
    recipe = m.Recipe.from_dict(d)
    assert m.recipe_hash(recipe) == typed_digest(d, "genairr_reference_source_recipe_v1")
    assert m.recipe_hash(recipe) not in canonical_json(recipe.to_dict())
    for key in ("commitment", "recipe_hash", "digest"):
        bad = copy.deepcopy(d)
        bad[key] = m.recipe_hash(recipe)
        with _raises("malformed", ""):
            m.Recipe.from_dict(bad)


def test_recipe_selection_order_is_preserved_as_authored():
    forward = m.Recipe.from_dict(recipe_dict())
    backward = recipe_dict()
    backward["selection"].reverse()
    backward = m.Recipe.from_dict(backward)
    assert [e.label for e in backward.selection] == list(reversed(SELECTED))
    assert m.recipe_hash(forward) != m.recipe_hash(backward)


def test_recipe_selects_distinct_identities_with_identical_dna():
    recipe = m.Recipe.from_dict(recipe_dict())
    by_label = {e.label: e for e in recipe.selection}
    assert by_label["IGHV1-69*01"].sequence_sha256 == by_label["IGHV1-69D*01"].sequence_sha256
    assert by_label["IGHV1-69*01"].record != by_label["IGHV1-69D*01"].record


def test_recipe_rejects_two_alternatives_for_one_executable_label():
    d = recipe_dict()
    d["selection"].append(selection_entry("IGHV1-58*03", fmt="airr_ex"))
    with _raises("unresolved_selection", "/selection"):
        m.Recipe.from_dict(d)
    d = recipe_dict()
    d["selection"].append(dict(d["selection"][0], label="IGHV1-58*99"))
    with _raises("unresolved_selection", "/selection"):
        m.Recipe.from_dict(d)


@pytest.mark.parametrize("edit,kind,path", [
    (lambda e: e.update(segment="C"), "malformed", "/selection/0/segment"),
    (lambda e: e.update(label="IGHV1-58 *03"), "malformed", "/selection/0/label"),
    (lambda e: e.update(label=""), "malformed", "/selection/0/label"),
    (lambda e: e["record"].update(pointer="GermlineSet/0"), "malformed", "/selection/0/record/pointer"),
    (lambda e: e["record"].update(pointer="/GermlineSet/~2"), "malformed", "/selection/0/record/pointer"),
    (lambda e: e["record"].update(pointer="/" + "a" * m.MAX_POINTER), "resource_limit", "/selection/0/record/pointer"),
    (lambda e: e["record"].update(asset="../airr"), "malformed", "/selection/0/record/asset"),
    (lambda e: e["record"].update(asset="IGH_VDJ.json"), "malformed", "/selection/0/record/asset"),
    (lambda e: e.update(sequence_sha256="sha256:abc"), "malformed", "/selection/0/sequence_sha256"),
    (lambda e: e.update(reason=" "), "malformed", "/selection/0/reason"),
    (lambda e: e.update(reason=None), "malformed", "/selection/0/reason"),
    (lambda e: e.update(reason="x" * (m.MAX_TEXT + 1)), "resource_limit", "/selection/0/reason"),
    (lambda e: e.update(predicate="all_functional"), "malformed", "/selection/0"),
    (lambda e: e.pop("reason"), "malformed", "/selection/0"),
])
def test_recipe_selection_entries_are_strict(edit, kind, path):
    d = recipe_dict()
    edit(d["selection"][0])
    with _raises(kind, path):
        m.Recipe.from_dict(d)


def test_recipe_selection_segments_must_fit_the_chain_type():
    lock = lock_dict("9606.IGKappa_VJ")
    d = recipe_dict(lock, selected=())
    d["selection"] = [dict(selection_entry("IGHD1-1*01"), record={"asset": "airr", "pointer": _pointer(0)})]
    with _raises("identity_mismatch", "/selection/0/segment"):
        m.Recipe.from_dict(d)


def test_recipe_null_and_missing_values_are_distinct():
    d = recipe_dict()
    d["reference"]["build_date"] = None
    d["reference"]["version"] = None
    recipe = m.Recipe.from_dict(d)
    assert recipe.reference.build_date is None and recipe.reference.version is None
    for key in ("build_date", "version"):
        d = recipe_dict()
        del d["reference"][key]
        with _raises("malformed", "/reference"):
            m.Recipe.from_dict(d)
    for bad_date in ("today", "2026-02-30", "2026-10-22T00:00:00", 20261022):
        d = recipe_dict()
        d["reference"]["build_date"] = bad_date
        with _raises("malformed", "/reference/build_date"):
            m.Recipe.from_dict(d)
    d = recipe_dict()
    d["reference"]["name"] = None
    with _raises("malformed", "/reference/name"):
        m.Recipe.from_dict(d)


@pytest.mark.parametrize("edit,kind,path", [
    (lambda d: d.update(schema_version=2), "unsupported_schema", "/schema_version"),
    (lambda d: d.update(schema_version=True), "malformed", "/schema_version"),
    (lambda d: d.update(artifact="genairr_reference_source_lock"), "malformed", "/artifact"),
    (lambda d: d["adapter"].update(version=2), "unsupported_profile", "/adapter"),
    (lambda d: d["projection"].update(version=2), "unsupported_profile", "/projection"),
    (lambda d: d["projection"].update(id="full_span"), "unsupported_profile", "/projection"),
    (lambda d: d.update(defaults_profile="genairr.unfitted_defaults.v2"), "unsupported_profile", "/defaults_profile"),
    (lambda d: d.update(source_lock_hash=None), "malformed", "/source_lock_hash"),
    (lambda d: d["reference"].update(locus="IGH", chain_type="vj"), "malformed", "/reference/chain_type"),
    (lambda d: d["reference"]["taxon"].update(taxonomy_id="9606"), "malformed", "/reference/taxon/taxonomy_id"),
])
def test_recipe_rejects_unknown_schema_and_profiles(edit, kind, path):
    d = recipe_dict()
    edit(d)
    with _raises(kind, path):
        m.Recipe.from_dict(d)


@pytest.mark.parametrize("edit,kind,path", [
    (lambda r: r.append({"capability": "all_functional", "evidence": "source"}), "malformed", "/requires/1/capability"),
    (lambda r: r.append({"capability": "canonical_productive", "evidence": "independent"}), "malformed",
     "/requires/1/evidence"),
    (lambda r: r.append({"capability": "observed_span_simulation", "evidence": "source"}), "malformed", "/requires"),
    (lambda r: r[0].update(level="source"), "malformed", "/requires/0"),
])
def test_recipe_requirements_use_the_runtime_vocabulary(edit, kind, path):
    d = recipe_dict()
    edit(d["requires"])
    with _raises(kind, path):
        m.Recipe.from_dict(d)


def test_recipe_requirements_keep_authored_order():
    d = recipe_dict()
    d["requires"] = [{"capability": "source_annotated_junction", "evidence": "source"},
                     {"capability": "observed_span_simulation", "evidence": "assertion"}]
    recipe = m.Recipe.from_dict(d)
    assert [r.capability for r in recipe.requires] == ["source_annotated_junction", "observed_span_simulation"]


CLAIM_ID = "claim:" + "a" * 64
FINDING_ID = "finding:" + "b" * 64


def decision(action="select_alternative", **overrides):
    entry = selection_entry("IGHV1-58*03")
    d = {"action": action, "findings": [FINDING_ID] if action == "select_alternative" else [],
         "claims": [] if action == "select_alternative" else [CLAIM_ID],
         "sequence_sha256": entry["sequence_sha256"],
         "reason": "The airr record ends at gene_end like every other selected V; the extension is not adopted",
         "evidence": [{"asset": "airr", "pointer": entry["record"]["pointer"] + "/coding_sequence"}]}
    d.update(overrides)
    return d


def test_recipe_decisions_support_only_the_three_actions():
    d = recipe_dict()
    d["decisions"] = [decision(), decision("withhold_claim"),
                      decision("accept_claim", claims=["claim:" + "c" * 64])]
    recipe = m.Recipe.from_dict(d)
    assert [x.action for x in recipe.decisions] == ["select_alternative", "withhold_claim", "accept_claim"]
    for action in ("edit_bases", "extend_terminal", "ignore_warnings", "guess_coordinates", "exec"):
        d["decisions"] = [decision(action=action)]
        with _raises("malformed", "/decisions/0/action"):
            m.Recipe.from_dict(d)


@pytest.mark.parametrize("entry,path", [
    (decision(findings=[]), "/decisions/0/findings"),
    (decision("withhold_claim", claims=[]), "/decisions/0/claims"),
    (decision("accept_claim", evidence=[]), "/decisions/0/evidence"),
    (decision("select_alternative", claims=[CLAIM_ID]), "/decisions/0/claims"),
    (decision(findings=["finding:" + "B" * 64]), "/decisions/0/findings/0"),
    (decision(findings=[FINDING_ID, FINDING_ID]), "/decisions/0/findings"),
    (decision("withhold_claim", claims=["IGHV1-58*03/anchor"]), "/decisions/0/claims/0"),
    (decision(sequence_sha256=None), "/decisions/0/sequence_sha256"),
    (decision(reason=""), "/decisions/0/reason"),
    (decision(evidence=[{"asset": "airr", "pointer": "/x"}, {"asset": "airr", "pointer": "/x"}]),
     "/decisions/0/evidence"),
    (dict(decision(), ignore_warnings=True), "/decisions/0"),
])
def test_recipe_decision_entries_are_strict(entry, path):
    d = recipe_dict()
    d["decisions"] = [entry]
    with _raises("malformed", path):
        m.Recipe.from_dict(d)


def test_recipe_one_claim_is_decided_once():
    d = recipe_dict()
    d["decisions"] = [decision("withhold_claim"), decision("accept_claim")]
    with _raises("malformed", "/decisions"):
        m.Recipe.from_dict(d)


def test_check_recipe_accepts_its_own_lock():
    lock = m.SourceLock.from_dict(lock_dict())
    m.check_recipe(m.Recipe.from_dict(recipe_dict()), lock)


def test_rehashed_recipe_for_another_lock_is_a_digest_mismatch():
    lock = m.SourceLock.from_dict(lock_dict())
    other = lock_dict(origin="local_file")
    with _raises("digest_mismatch", "/source_lock_hash"):
        m.check_recipe(m.Recipe.from_dict(recipe_dict(other)), lock)


@pytest.mark.parametrize("edit,path", [
    (lambda r: r["taxon"].update(scientific_name="Mus musculus"), "/reference/taxon"),
    (lambda r: r["taxon"].update(taxonomy_id=10090), "/reference/taxon"),
    (lambda r: r.update(locus="TRB"), "/reference/locus"),
    (lambda r: r.update(locus="IGK", chain_type="vj"), "/reference/locus"),
])
def test_recipe_cannot_relabel_the_source_identity(edit, path):
    lock = m.SourceLock.from_dict(lock_dict())
    d = recipe_dict()
    edit(d["reference"])
    if d["reference"]["chain_type"] == "vj":
        d["selection"] = [e for e in d["selection"] if e["segment"] != "D"]
    with _raises("identity_mismatch", path):
        m.check_recipe(m.Recipe.from_dict(d), lock)


def test_recipe_adapter_must_equal_the_lock_adapter(monkeypatch):
    from GenAIRR._reference_sources import common
    monkeypatch.setitem(common.SUPPORTED_ADAPTERS, ("ogrdb_germline_set", 2),
                        ("ogrdb", (("ogrdb_api_v2_human_ig_observed_v1", 1),)))
    lock = m.SourceLock.from_dict(lock_dict())
    d = recipe_dict()
    d["adapter"]["version"] = 2
    with _raises("identity_mismatch", "/adapter"):
        m.check_recipe(m.Recipe.from_dict(d), lock)


def test_recipe_records_must_name_locked_record_exports():
    lock = m.SourceLock.from_dict(lock_dict(formats=("airr", "gapped", "gapped_ex")))
    d = recipe_dict(lock.to_dict(), selected=SELECTED[:-1])
    d["selection"].append(selection_entry("IGHJ4*02", fmt="airr_ex"))
    with _raises("identity_mismatch", "/selection/4/record/asset"):
        m.check_recipe(m.Recipe.from_dict(d), lock)
    d = recipe_dict(lock.to_dict())
    d["selection"][0]["record"]["asset"] = "gapped"
    with _raises("identity_mismatch", "/selection/0/record/asset"):
        m.check_recipe(m.Recipe.from_dict(d), lock)
    d = recipe_dict(lock.to_dict())
    d["decisions"] = [decision(evidence=[{"asset": "airr_ex", "pointer": "/x"}])]
    with _raises("identity_mismatch", "/decisions/0/evidence/0/asset"):
        m.check_recipe(m.Recipe.from_dict(d), lock)


def test_changed_recipe_text_changes_the_recipe_identity_even_with_identical_dna():
    first = m.Recipe.from_dict(recipe_dict())
    d = recipe_dict()
    d["selection"][0]["reason"] = "Reviewed again; same supplied sequence"
    second = m.Recipe.from_dict(d)
    assert [e.sequence_sha256 for e in first.selection] == [e.sequence_sha256 for e in second.selection]
    assert m.recipe_hash(first) != m.recipe_hash(second)


# ── source inventory ─────────────────────────────────────────────────

L = m.Locator


def _fasta(fmt, set_id=IGH):
    source = _source_set(set_id)
    text = (ROOT / "ogrdb" / set_id / str(source["release"]) / f"{fmt}.fasta").read_bytes().decode("ascii")
    records, name = [], None
    for line in text.split("\n"):
        if line.startswith(">"):
            records.append([line[1:], ""])
        elif line:
            records[-1][1] += line
    return records


def _gapped_index(label, fmt):
    (found,) = [(i, seq) for i, (name, seq) in enumerate(_fasta(fmt)) if name == label]
    return found


def _blocks(gapped):
    """Exact interval map from ungapped positions to gapped columns."""
    blocks, position = [], 0
    for column, char in enumerate(gapped):
        if char == ".":
            continue
        if blocks and blocks[-1][1][1] == column:
            blocks[-1] = ((blocks[-1][0][0], position + 1), (blocks[-1][1][0], column + 1))
        else:
            blocks.append(((position, position + 1), (column, column + 1)))
        position += 1
    return tuple(m.Block(input=i, output=o) for i, o in blocks)


def _sorted(locators):
    return tuple(sorted(locators, key=L.sort_key))


INV_RECORDS = (("airr", "IGHV1-58*03"), ("airr", "IGHV1-69*01"), ("airr", "IGHV1-69D*01"), ("airr", "IGHD1-1*01"),
               ("airr", "IGHJ4*02"), ("airr_ex", "IGHV1-58*03"), ("airr_ex", "IGHV1-69*01"),
               ("airr_ex", "IGHD1-1*01"), ("airr_ex", "IGHJ4*02"))


class Inv:
    """A small realistic inventory over pinned IGH release 10 records, built field by field."""

    def __init__(self, lock=None):
        self.lock = lock or m.SourceLock.from_dict(lock_dict())
        self.records, self.candidates, self.observations, self.support = [], [], [], []
        self.mentioned, self.mappings, self.claims, self.findings = [], [], [], []
        self.rec, self.cand, self.obs = {}, {}, {}
        for fmt, label in INV_RECORDS:
            index, raw = _find(label, fmt)
            loc = L(fmt, _pointer(index))
            self.rec[fmt, label] = loc
            self.records.append(m.SourceRecord(
                locator=loc, status="parsed", problem=None, label=raw["label"],
                provider_id=raw["allele_description_id"], provider_ref=raw["allele_description_ref"],
                segment=raw["sequence_type"], taxonomy_id=int(raw["species"]["id"].split(":")[1]),
                locus=raw["locus"]))
            fields = ["label", "coding_sequence", "functional", "sequence_type"]
            if raw["sequence_type"] == "V":
                fields.append("v_gene_delineations/0/cdr3_start")
            if raw["sequence_type"] == "J":
                fields += ["j_codon_frame", "j_cdr3_end"]
            if raw["paralogs"]:
                fields.append("paralogs/0")
            for field in fields:
                value = raw
                for token in field.split("/"):
                    value = value[int(token)] if isinstance(value, list) else value[token]
                self._observe(L(fmt, f"{loc.pointer}/{field}"), loc, value)
            cand = L(fmt, loc.pointer + "/coding_sequence")
            self.cand[fmt, label] = cand
            self.candidates.append(m.SequenceCandidate(locator=cand, record=loc, sequence=raw["coding_sequence"],
                                                       sha256=m.sequence_digest(raw["coding_sequence"])))
            self.claims.append(self.claim(loc, "functionality", "functional", [self.obs_at(loc, "functional")]))
        # Paralog-only identity in the extended export, mentioned by IGHV1-69*01.
        ex69 = self.rec["airr_ex", "IGHV1-69*01"]
        self.mentioned.append(m.MentionedIdentity(asset="airr_ex", segment="V", label="IGHV1-69D*01", mentions=(
            m.Mention(record=ex69, locator=L("airr_ex", ex69.pointer + "/paralogs/0"), relationship="paralog"),)))
        # Embedded, unretrieved repository support of one record.
        rec5803 = self.rec["airr", "IGHV1-58*03"]
        _, raw5803 = _find("IGHV1-58*03")
        entry = raw5803["unrearranged_support"][0]
        self.support.append(m.SupportReference(
            locator=L("airr", rec5803.pointer + "/unrearranged_support/0"), record=rec5803, kind="unrearranged",
            repository=entry["repository_name"], reference=entry["repository_ref"], status="embedded_unretrieved"))
        # Set-level observation: the API version is not the schema version or the data release.
        info = L("airr", "/Info/version")
        self._observe(info, None, _load_json(ROOT / "ogrdb" / IGH / "10" / "airr.json")["Info"]["version"],
                      interpretation="preserved")
        # Variant-specific gapped alignments, each mapping one record occurrence to its own variant.
        for fmt, gapped_fmt in (("airr", "gapped"), ("airr_ex", "gapped_ex")):
            gindex, gapped = _gapped_index("IGHV1-58*03", gapped_fmt)
            gloc = L(gapped_fmt, f"/{gindex}")
            self._observe(gloc, None, {"header": "IGHV1-58*03", "sequence": gapped})
            cand = self.cand[fmt, "IGHV1-58*03"]
            sequence = gapped.replace(".", "")
            mapping = m.Mapping(source_lock_hash=self.lock.commitment, kind="gapped_alignment", input=cand,
                                output=gloc,
                                input_sha256=m.sequence_digest(sequence),
                                output_sha256=m.text_digest(gapped), input_space="coding_sequence",
                                output_space="imgt_gapped", blocks=_blocks(gapped), transform="remove_imgt_gaps",
                                validation="exact")
            self.mappings.append(mapping)
        # V junction anchor: source-stated, and re-derived from the variant's own gapped record.
        cand, rec = self.cand["airr", "IGHV1-58*03"], self.rec["airr", "IGHV1-58*03"]
        cdr3 = self.obs_at(rec, "v_gene_delineations/0/cdr3_start")
        anchor = _find("IGHV1-58*03")[1]["v_gene_delineations"][0]["cdr3_start"] - 4
        gapped = _gapped_index("IGHV1-58*03", "gapped")[1]
        assert len(gapped[:312].replace(".", "")) - 3 == anchor
        self.anchor = anchor
        self.claims.append(self.claim(cand, "anchor", anchor, [cdr3], system="coding_sequence"))
        self.claims.append(self.claim(cand, "anchor", anchor, [L("gapped", f"/{_gapped_index('IGHV1-58*03', 'gapped')[0]}")],
                                      system="coding_sequence", derivation="imgt_alignment",
                                      mappings=[self.mappings[0].id]))
        self.claims.append(self.claim(cand, "coding_phase", anchor % 3, [cdr3], system="coding_sequence",
                                      derivation="anchor_derived"))
        self.claims.append(self.claim(cand, "phase_independence", "anchor_derived", [cdr3],
                                      derivation="anchor_derived"))
        # J: the reported phase value is accepted; whether it is independent of the landmark is unresolved.
        jcand, jrec = self.cand["airr", "IGHJ4*02"], self.rec["airr", "IGHJ4*02"]
        _, jraw = _find("IGHJ4*02")
        frame = self.obs_at(jrec, "j_codon_frame")
        self.claims.append(self.claim(jcand, "anchor", jraw["j_cdr3_end"] - 1, [self.obs_at(jrec, "j_cdr3_end")],
                                      system="coding_sequence"))
        self.claims.append(self.claim(jcand, "coding_phase", jraw["j_codon_frame"] - 1, [frame],
                                      system="coding_sequence"))
        self.claims.append(self.claim(jcand, "phase_independence", None, [frame], status="unresolved",
                                      support="none"))
        # Extended-export alternative: provider 3' extension beyond the observed span.
        ex = self.cand["airr_ex", "IGHV1-58*03"]
        self.claims.append(self.claim(ex, "extension", {"side": "3", "length": 1, "matches_record_context": True},
                                      [L("airr_ex", ex.pointer), L("airr", cand.pointer)],
                                      derivation="variant_comparison"))
        # D: no universal coding phase; its physical 3' end stays unresolved.
        dcand = self.cand["airr", "IGHD1-1*01"]
        self.claims.append(self.claim(dcand, "boundary_3", None, [L("airr", dcand.pointer)], status="unresolved",
                                      support="none"))
        self.findings.append(m.Finding(
            source_lock_hash=self.lock.commitment, code="variant_coding_sequence_differs", severity="warning", category="contradictory_evidence",
            subjects=_sorted([cand, ex]), locators=_sorted([cand, ex]),
            observations={"label": "IGHV1-58*03", "extension_3": "A"},
            capabilities=("observed_span_simulation", "canonical_productive")))
        self.findings.append(m.Finding(
            source_lock_hash=self.lock.commitment, code="paralog_only_identity", severity="info", category="informational",
            subjects=(self.rec["airr_ex", "IGHV1-69*01"],), locators=(L("airr_ex", ex69.pointer + "/paralogs/0"),),
            observations={"label": "IGHV1-69D*01"}, capabilities=()))
        self.findings.append(m.Finding(
            source_lock_hash=self.lock.commitment, code="api_version_is_not_release", severity="info", category="unsupported_interpretation",
            subjects=(), locators=(info,), observations={"value": "2.0.0"}, capabilities=()))

    def _observe(self, locator, record, value, interpretation="interpreted"):
        self.obs[locator.id] = locator
        self.observations.append(m.Observation(locator=locator, record=record, presence="present", value=value,
                                               interpretation=interpretation))

    def obs_at(self, record, field):
        return self.obs[L(record.asset, f"{record.pointer}/{field}").id]

    def claim(self, subject, predicate, value, evidence, *, system=None, derivation="stated", status="accepted",
              support="source", mappings=()):
        return m.Claim(source_lock_hash=self.lock.commitment, subject=subject, predicate=predicate, value=value, coordinate_system=system, support=support,
                       derivation=derivation, status=status, evidence=_sorted(evidence), mappings=tuple(mappings))

    def build(self, **overrides):
        fields = dict(source_lock_hash=self.lock.commitment, adapter=self.lock.adapter, records=self.records,
                      candidates=self.candidates, mentioned_identities=self.mentioned,
                      observations=self.observations, support=self.support, claims=self.claims,
                      mappings=self.mappings, findings=self.findings)
        fields.update(overrides)
        return m.SourceInventory.assemble(**fields)


@pytest.fixture(scope="module")
def inv():
    return Inv()


@pytest.fixture(scope="module")
def inventory(inv):
    return inv.build()


def _inventory_dict(inventory):
    return copy.deepcopy(inventory.to_dict())


def test_inventory_round_trips_and_is_byte_stable(inventory):
    d = inventory.to_dict()
    again = m.SourceInventory.from_dict(d)
    assert again == inventory
    text = canonical_json(d)
    assert canonical_json(m.SourceInventory.from_json(text).to_dict()) == text
    assert d["artifact"] == "genairr_reference_source_inventory" and d["schema_version"] == 1
    assert m.inventory_hash(inventory) == typed_digest(d, "genairr_reference_source_inventory_v1")


def test_inventory_assemble_puts_every_list_in_canonical_order(inv, inventory):
    shuffled = inv.build(records=list(reversed(inv.records)), claims=list(reversed(inv.claims)),
                         observations=list(reversed(inv.observations)), findings=list(reversed(inv.findings)))
    assert shuffled == inventory
    with _raises("malformed", "/records"):
        m.SourceInventory(**{**_fields_of(inventory), "records": tuple(reversed(inventory.records))})


def _fields_of(inventory):
    return {k: getattr(inventory, k) for k in ("source_lock_hash", "adapter", "records", "candidates",
                                                "mentioned_identities", "observations", "support", "claims",
                                                "mappings", "findings")}


def test_inventory_records_are_occurrences_keyed_by_asset_and_pointer(inventory):
    ids = [r.id for r in inventory.records]
    assert len(ids) == len(set(ids)) == len(INV_RECORDS)
    by_provider = {}
    for record in inventory.records:
        by_provider.setdefault(record.provider_id, []).append(record.locator.asset)
    # The same provider id in both exports stays two record occurrences.
    assert sorted(by_provider["OGRDB:A00014"]) == ["airr", "airr_ex"]


def test_inventory_counts_keep_occurrences_candidates_and_mentions_distinct(inventory):
    counts = inventory.export_counts("airr_ex")
    assert counts == {"record_occurrences": {"V": 2, "D": 1, "J": 1, "unassigned": 0},
                      "sequence_candidates": {"V": 2, "D": 1, "J": 1, "unassigned": 0},
                      "mentioned_identities": {"V": 1, "D": 0, "J": 0, "unassigned": 0}}
    assert inventory.export_counts("airr")["mentioned_identities"]["V"] == 0
    assert inventory.export_counts("airr")["record_occurrences"]["V"] == 3


def test_identical_dna_identities_cannot_cite_each_others_evidence(inventory, inv):
    own, twin = inv.cand["airr", "IGHV1-69D*01"], inv.cand["airr", "IGHV1-69*01"]
    assert inventory.candidate(own).sha256 == inventory.candidate(twin).sha256
    assert inventory.candidate(own).record != inventory.candidate(twin).record
    own_cdr3 = inv.obs_at(inv.rec["airr", "IGHV1-69D*01"], "v_gene_delineations/0/cdr3_start")
    twin_cdr3 = inv.obs_at(inv.rec["airr", "IGHV1-69*01"], "v_gene_delineations/0/cdr3_start")
    inv.build(claims=inv.claims + [inv.claim(own, "anchor", 285, [own_cdr3], system="coding_sequence")])
    with _raises("malformed") as caught:
        inv.build(claims=inv.claims + [inv.claim(own, "anchor", 285, [twin_cdr3], system="coding_sequence")])
    assert caught.error.path.startswith("/claims/") and caught.error.path.endswith("/evidence/0")
    assert "IGHV1-69" in caught.error.message or inv.rec["airr", "IGHV1-69D*01"].id in caught.error.message


def test_evidence_roots_belong_to_the_claim_subject(inv):
    target = inv.cand["airr", "IGHV1-69*01"]
    foreign = inv.obs_at(inv.rec["airr", "IGHV1-58*03"], "v_gene_delineations/0/cdr3_start")
    with _raises("malformed") as caught:
        inv.build(claims=inv.claims + [inv.claim(target, "anchor", inv.anchor, [foreign], system="coding_sequence")])
    assert caught.error.path.endswith("/evidence/0")
    # A set-level field is not evidence about one allele.
    with _raises("malformed"):
        inv.build(claims=inv.claims + [inv.claim(target, "anchor", inv.anchor, [L("airr", "/Info/version")],
                                                 system="coding_sequence")])
    # The same-label record of the other export is citable only by a variant comparison.
    ex, cand = inv.cand["airr_ex", "IGHV1-58*03"], inv.cand["airr", "IGHV1-58*03"]
    with _raises("malformed"):
        inv.build(claims=inv.claims + [inv.claim(ex, "boundary_5", "complete", [L("airr", cand.pointer)])])
    inv.build(claims=inv.claims + [inv.claim(ex, "boundary_3", "partial", [L("airr", cand.pointer)],
                                             derivation="variant_comparison")])
    other_label = inv.cand["airr", "IGHV1-69*01"]
    with _raises("malformed"):
        inv.build(claims=inv.claims + [inv.claim(ex, "boundary_3", "partial", [L("airr", other_label.pointer)],
                                                 derivation="variant_comparison")])


def test_claims_cite_only_mappings_of_their_own_candidate(inv):
    cand = inv.cand["airr", "IGHV1-58*03"]
    ex_mapping = [x for x in inv.mappings if x.input.asset == "airr_ex"][0]
    claims = [c for c in inv.claims if not (c.subject == cand and c.derivation == "imgt_alignment")]
    crossed = inv.claim(cand, "anchor", inv.anchor, [ex_mapping.output], system="coding_sequence",
                        derivation="imgt_alignment", mappings=[ex_mapping.id])
    with _raises("malformed") as caught:
        inv.build(claims=claims + [crossed])
    assert caught.error.path.startswith("/claims/") and caught.error.path.endswith("/mappings/0")
    record = inv.rec["airr", "IGHV1-58*03"]
    functional = inv.claim(record, "functionality", "functional", [inv.obs_at(record, "functional")],
                           mappings=[inv.mappings[0].id])
    with _raises("malformed") as caught:
        inv.build(claims=[c for c in inv.claims if c.subject != record] + [functional])
    assert caught.error.path.endswith("/mappings/0")


def test_accept_claim_evidence_is_bound_to_the_claim_subject(inventory, inv):
    jcand, jrec = inv.cand["airr", "IGHJ4*02"], inv.rec["airr", "IGHJ4*02"]
    phase = _claim(inventory, "coding_phase", jcand)
    sha = inventory.candidate(jcand).sha256
    for evidence in ({"asset": "airr", "pointer": "/Info/version"},
                     inv.obs_at(inv.rec["airr", "IGHV1-58*03"], "functional").to_dict()):
        entry = decision("accept_claim", claims=[phase.id], sequence_sha256=sha, evidence=[evidence])
        with _raises("identity_mismatch", "/decisions/0/evidence/0"):
            m.check_decisions(_recipe(inventory, decisions=[entry]), inventory)
    entry = decision("accept_claim", claims=[phase.id], sequence_sha256=sha,
                     evidence=[inv.obs_at(jrec, "j_codon_frame").to_dict()])
    m.check_decisions(_recipe(inventory, decisions=[entry]), inventory)


def test_phase_value_and_phase_independence_are_separate_claims(inventory, inv):
    jcand = inv.cand["airr", "IGHJ4*02"]
    phase = [c for c in inventory.claims if c.subject == jcand and c.predicate == "coding_phase"]
    independence = [c for c in inventory.claims if c.subject == jcand and c.predicate == "phase_independence"]
    assert [(c.value.get(), c.status, c.derivation) for c in phase] == [(2, "accepted", "stated")]
    assert [(c.value.get(), c.status, c.support) for c in independence] == [(None, "unresolved", "none")]
    vcand = inv.cand["airr", "IGHV1-58*03"]
    (v_independence,) = [c for c in inventory.claims if c.subject == vcand and c.predicate == "phase_independence"]
    assert (v_independence.value.get(), v_independence.status) == ("anchor_derived", "accepted")


def test_claim_and_finding_ids_are_content_derived_and_verified(inventory):
    d = _inventory_dict(inventory)
    claim = d["claims"][0]
    body = {k: v for k, v in claim.items() if k != "id"}
    assert claim["id"] == "claim:" + typed_digest(body, "genairr_reference_source_claim_id_v1").split(":")[1]
    d["claims"][0]["id"] = "claim:" + "0" * 64
    with _raises("digest_mismatch", "/claims/0/id"):
        m.SourceInventory.from_dict(d)
    d = _inventory_dict(inventory)
    (index,) = [i for i, f in enumerate(d["findings"]) if f["code"] == "paralog_only_identity"]
    d["findings"][index]["observations"]["label"] = "IGHV3-23D*01"
    with _raises("digest_mismatch", f"/findings/{index}/id"):
        m.SourceInventory.from_dict(d)
    d = _inventory_dict(inventory)
    d["mappings"][0]["blocks"][0]["output"][0] += 0
    d["mappings"][0]["transform"] = "remove_imgt_gaps"
    assert m.SourceInventory.from_dict(d) == inventory


@pytest.mark.parametrize("edit,kind,path_prefix", [
    (lambda d: d["records"].append(copy.deepcopy(d["records"][0])), "malformed", "/records"),
    (lambda d: d["records"][0].update(label=None), "malformed", "/records/0/label"),
    (lambda d: d["records"][0].update(segment=None), "malformed", "/records/0/segment"),
    (lambda d: d["records"][0].update(status="malformed"), "malformed", "/records/0/problem"),
    (lambda d: d["records"][0].update(problem="unreadable"), "malformed", "/records/0/problem"),
    (lambda d: d["records"][0].update(taxonomy_id=True), "malformed", "/records/0/taxonomy_id"),
    (lambda d: d["records"][0]["locator"].update(asset="gapped"), "malformed", "/records/0/locator"),
    (lambda d: d["candidates"][0].update(sha256="sha256:" + "0" * 64), "digest_mismatch", "/candidates/0/sha256"),
    (lambda d: d["candidates"][0].update(sequence=d["candidates"][0]["sequence"].lower(),
                                         sha256=m.sequence_digest(d["candidates"][0]["sequence"].lower())),
     "malformed", "/candidates/0/sequence"),
    (lambda d: d["candidates"][0]["record"].update(pointer="/GermlineSet/0/allele_descriptions/1"), "malformed",
     "/candidates/0"),
    (lambda d: d["observations"][0].update(presence="absent"), "malformed", "/observations/0/value"),
    (lambda d: d["observations"][0].update(interpretation="trusted"), "malformed",
     "/observations/0/interpretation"),
    (lambda d: d["observations"][0]["record"].update(asset="airr_ex"), "malformed", "/observations/0"),
    (lambda d: d["support"][0].update(status="verified"), "malformed", "/support/0/status"),
    (lambda d: d["mentioned_identities"][0].update(label="IGHV1-69*01"), "malformed", "/mentioned_identities/0"),
    (lambda d: d["mentioned_identities"][0].update(segment="D"), "malformed", "/mentioned_identities/0"),
    (lambda d: d["mentioned_identities"][0].update(mentions=[]), "malformed", "/mentioned_identities/0/mentions"),
    (lambda d: d["mentioned_identities"][0]["mentions"][0].update(relationship="ortholog"), "malformed",
     "/mentioned_identities/0/mentions/0/relationship"),
])
def test_inventory_items_are_strict(inventory, edit, kind, path_prefix):
    d = _inventory_dict(inventory)
    edit(d)
    with _raises(kind) as caught:
        m.SourceInventory.from_dict(d)
    assert caught.error.path.startswith(path_prefix), caught.error


def test_inventory_keeps_malformed_and_unsupported_records(inv):
    loc = L("airr", _pointer(0))
    odd = m.SourceRecord(locator=loc, status="unsupported", problem="sequence_type 'C' is not a V, D or J segment",
                         label="IGHC*01 (constant)", provider_id=None, provider_ref=None, segment=None,
                         taxonomy_id=None, locus="IGH")
    inventory = inv.build(records=inv.records + [odd])
    assert inventory.record(loc).status == "unsupported"
    assert inventory.export_counts("airr")["record_occurrences"]["unassigned"] == 1
    with _raises("malformed", "/label"):
        m.SourceRecord(locator=loc, status="parsed", problem=None, label="IGHC*01 (constant)", provider_id=None,
                       provider_ref=None, segment="V", taxonomy_id=None, locus=None)


def _rebuild_claim(inv, index, **changes):
    claims = list(inv.claims)
    fields = {k: getattr(claims[index], k) for k in ("source_lock_hash", "subject", "predicate", "value", "coordinate_system", "support",
                                                      "derivation", "status", "evidence", "mappings")}
    fields["value"] = fields["value"].get()
    fields.update(changes)
    claims[index] = m.Claim(**fields)
    return inv.build(claims=claims)


def _claim_index(inv, predicate, label="IGHV1-58*03", fmt="airr", **match):
    subjects = {inv.cand[fmt, label], inv.rec[fmt, label]}
    for i, claim in enumerate(inv.claims):
        if claim.predicate == predicate and claim.subject in subjects and all(
                getattr(claim, k) == v for k, v in match.items()):
            return i
    raise AssertionError(predicate)


@pytest.mark.parametrize("predicate,label,changes,path", [
    ("phase_independence", "IGHJ4*02", {"value": "independent"}, "/value"),
    ("coding_phase", "IGHJ4*02", {"value": None}, "/value"),
    ("coding_phase", "IGHJ4*02", {"value": 3}, "/value"),
    ("coding_phase", "IGHJ4*02", {"value": True}, "/value"),
    ("coding_phase", "IGHJ4*02", {"support": "none"}, "/support"),
    ("coding_phase", "IGHJ4*02", {"coordinate_system": None}, "/coordinate_system"),
    ("coding_phase", "IGHJ4*02", {"coordinate_system": "imgt_gapped"}, "/coordinate_system"),
    ("coding_phase", "IGHJ4*02", {"status": "candidate"}, "/support"),
    ("coding_phase", "IGHJ4*02", {"evidence": ()}, "/evidence"),
    ("coding_phase", "IGHJ4*02", {"derivation": "guessed"}, "/derivation"),
    ("anchor", "IGHJ4*02", {"value": -1}, "/value"),
    ("anchor", "IGHJ4*02", {"derivation": "anchor_derived"}, "/derivation"),
    ("functionality", "IGHJ4*02", {"value": "unknown"}, "/value"),
    ("functionality", "IGHJ4*02", {"value": False}, "/value"),
    ("functionality", "IGHJ4*02", {"coordinate_system": "coding_sequence"}, "/coordinate_system"),
    ("extension", "IGHV1-58*03", {"value": {"side": "3", "length": 0, "matches_record_context": True}},
     "/value/length"),
    ("extension", "IGHV1-58*03", {"value": {"side": "3", "length": 1}}, "/value"),
])
def test_claim_values_and_vocabularies_are_strict(inv, predicate, label, changes, path):
    fmt = "airr_ex" if predicate == "extension" else "airr"
    index = _claim_index(inv, predicate, label, fmt)
    with _raises("malformed") as caught:
        _rebuild_claim(inv, index, **changes)
    assert caught.error.path.endswith(path), caught.error


def test_claims_are_checked_against_their_subjects_and_evidence(inv):
    jcand = inv.cand["airr", "IGHJ4*02"]
    dcand = inv.cand["airr", "IGHD1-1*01"]
    cdr3 = inv.obs_at(inv.rec["airr", "IGHV1-58*03"], "v_gene_delineations/0/cdr3_start")
    cases = [
        # No universal D coding phase, and no anchors on D.
        (inv.claim(dcand, "coding_phase", 0, [L("airr", dcand.pointer)], system="coding_sequence"), "malformed"),
        (inv.claim(dcand, "anchor", 3, [L("airr", dcand.pointer)], system="coding_sequence"), "malformed"),
        # v_subregions only on V; functionality only on a record.
        (inv.claim(jcand, "v_subregions", {"FWR1": [0, 10]}, [L("airr", jcand.pointer)], system="coding_sequence",
                   derivation="imgt_alignment", mappings=[inv.mappings[0].id]), "malformed"),
        (inv.claim(jcand, "functionality", "functional", [L("airr", jcand.pointer)]), "malformed"),
        # Anchor beyond the candidate's coding sequence.
        (inv.claim(jcand, "anchor", 46, [L("airr", jcand.pointer)], system="coding_sequence"), "malformed"),
        # Evidence roots are observed source fields.
        (inv.claim(jcand, "anchor", 14, [L("airr", jcand.pointer + "/nowhere")], system="coding_sequence"),
         "malformed"),
        # Mappings must exist; alignment-derived claims need one.
        (inv.claim(jcand, "anchor", 14, [L("airr", jcand.pointer)], system="coding_sequence",
                   mappings=["mapping:" + "0" * 64]), "malformed"),
        # Unknown subject.
        (inv.claim(L("airr", "/GermlineSet/0/allele_descriptions/0/coding_sequence"), "anchor", 1, [cdr3],
                   system="coding_sequence"), "malformed"),
    ]
    for claim, kind in cases:
        with _raises(kind):
            inv.build(claims=inv.claims + [claim])
    with _raises("malformed", "/mappings"):
        inv.claim(inv.cand["airr", "IGHV1-58*03"], "anchor", inv.anchor, [cdr3], system="coding_sequence",
                  derivation="imgt_alignment")


def test_claim_statuses_express_agreement_and_conflict(inv):
    vcand = inv.cand["airr", "IGHV1-58*03"]
    cdr3 = inv.obs_at(inv.rec["airr", "IGHV1-58*03"], "v_gene_delineations/0/cdr3_start")
    disagreeing = inv.claim(vcand, "anchor", inv.anchor + 3, [L("airr", vcand.pointer)], system="coding_sequence",
                            derivation="unique_motif")
    with _raises("malformed"):
        inv.build(claims=inv.claims + [disagreeing])
    lone_conflict = inv.claim(vcand, "boundary_5", "complete", [cdr3], status="conflicting")
    with _raises("malformed"):
        inv.build(claims=inv.claims + [lone_conflict])
    pair = [inv.claim(vcand, "boundary_5", "complete", [cdr3], status="conflicting"),
            inv.claim(vcand, "boundary_5", "partial", [L("airr", vcand.pointer)], status="conflicting")]
    assert sum(c.status == "conflicting" for c in inv.build(claims=inv.claims + pair).claims) == 2
    motif = inv.claim(vcand, "anchor", inv.anchor + 3, [L("airr", vcand.pointer)], system="coding_sequence",
                      derivation="unique_motif", status="candidate", support="none")
    inv.build(claims=inv.claims + [motif])
    independent = inv.claim(vcand, "phase_independence", "independent", [cdr3])
    with _raises("malformed"):
        inv.build(claims=[c for c in inv.claims if not (c.subject == vcand and c.predicate == "phase_independence")]
                  + [independent])


def test_mappings_bind_one_record_occurrence_to_one_gapped_record_of_its_own_variant(inv):
    good = inv.mappings[0]
    fields = {k: getattr(good, k) for k in ("source_lock_hash", "kind", "input", "output", "input_sha256", "output_sha256", "input_space",
                                             "output_space", "blocks", "transform", "validation")}
    ex_gapped = inv.mappings[1].output
    with _raises("malformed", "/output"):
        m.Mapping(**{**fields, "output": ex_gapped})
    with _raises("malformed", "/output"):
        m.Mapping(**{**fields, "output": L("gapped", "/GermlineSet/0")})
    with _raises("malformed", "/blocks"):
        m.Mapping(**{**fields, "blocks": ()})
    with _raises("malformed", "/blocks"):
        m.Mapping(**{**fields, "blocks": (m.Block((0, 10), (0, 10)), m.Block((5, 20), (30, 45)))})
    with _raises("malformed", "/output"):
        m.Block((0, 10), (0, 11))
    with _raises("malformed", "/input"):
        m.Block((5, 5), (0, 0))
    with _raises("malformed", "/input"):
        m.Block((0, True), (0, 1))
    with _raises("malformed", "/transform"):
        m.Mapping(**{**fields, "transform": "align_anything"})
    mismatch = m.Mapping(**{**fields, "blocks": (), "validation": "mismatch"})
    assert mismatch.validation == "mismatch"
    with _raises("digest_mismatch"):
        inv.build(mappings=[m.Mapping(**{**fields, "input_sha256": "sha256:" + "1" * 64}), inv.mappings[1]])
    with _raises("malformed"):
        inv.build(mappings=inv.mappings + [mismatch])
    other = inv.cand["airr", "IGHV1-69*01"]
    with _raises("malformed"):
        inv.build(mappings=inv.mappings + [m.Mapping(**{**fields, "input": other, "blocks": (), "validation": "mismatch",
                                                         "input_sha256": inv.build().candidate(other).sha256})])
    short = fields["blocks"][:-1]
    with _raises("malformed"):
        inv.build(mappings=[m.Mapping(**{**fields, "blocks": short}), inv.mappings[1]],
                  claims=[c for c in inv.claims if not c.mappings])


def test_coding_interval_mapping_converts_record_sequence_coordinates(inv):
    index, raw = _find("IGHJ4*02")
    rec, cand = inv.rec["airr", "IGHJ4*02"], inv.cand["airr", "IGHJ4*02"]
    start, end = raw["gene_start"] - 1, raw["gene_end"]
    assert raw["sequence"][start:end] == raw["coding_sequence"]
    seq_loc = L("airr", rec.pointer + "/sequence")
    observation = m.Observation(locator=seq_loc, record=rec, presence="present", value=raw["sequence"],
                                interpretation="interpreted")
    mapping = m.Mapping(source_lock_hash=inv.lock.commitment, kind="coding_interval", input=seq_loc, output=cand,
                        input_sha256=m.text_digest(raw["sequence"]),
                        output_sha256=m.sequence_digest(raw["coding_sequence"]), input_space="record_sequence",
                        output_space="coding_sequence", blocks=(m.Block((start, end), (0, end - start)),),
                        transform="extract_interval", validation="exact")
    inventory = inv.build(observations=inv.observations + [observation], mappings=inv.mappings + [mapping])
    assert mapping in inventory.mappings
    interval = inv.claim(cand, "coding_interval", [start, end], [seq_loc], system="record_sequence",
                         mappings=[mapping.id])
    inv.build(observations=inv.observations + [observation], mappings=inv.mappings + [mapping],
              claims=inv.claims + [interval])
    with _raises("malformed"):
        m.Mapping(source_lock_hash=inv.lock.commitment, kind="coding_interval", input=seq_loc, output=cand,
                  input_sha256=mapping.input_sha256,
                  output_sha256=mapping.output_sha256, input_space="record_sequence", output_space="coding_sequence",
                  blocks=(m.Block((start, end), (1, end - start + 1)),), transform="extract_interval",
                  validation="exact")


def test_findings_are_strict(inv):
    base = inv.findings[0]
    fields = {k: getattr(base, k) for k in ("source_lock_hash", "code", "severity", "category", "subjects", "locators", "observations",
                                             "capabilities")}
    fields["observations"] = fields["observations"].get()
    for change, path in (({"code": "Variant differs"}, "/code"), ({"severity": "fatal"}, "/severity"),
                         ({"category": "biology"}, "/category"), ({"observations": ["x"]}, "/observations"),
                         ({"capabilities": ("canonical_productive", "observed_span_simulation")}, "/capabilities"),
                         ({"capabilities": ("all_functional",)}, "/capabilities/0"),
                         ({"subjects": (), "locators": ()}, "/subjects"),
                         ({"subjects": tuple(reversed(base.subjects))}, "/subjects")):
        with _raises("malformed", path):
            m.Finding(**{**fields, **change})
    unknown = m.Finding(**{**fields, "subjects": (L("airr", _pointer(0)),)})
    with _raises("malformed"):
        inv.build(findings=inv.findings + [unknown])


def test_preserved_values_are_finite_bounded_json(inv):
    loc = L("airr", "/Info/title")
    for value, kind in ((float("nan"), "malformed"), (float("inf"), "malformed"), ({1: "x"}, "malformed"),
                        ((1, 2), "malformed"), ({"\ud800": 1}, "malformed"),
                        (_nested(m.MAX_JSON_DEPTH + 1), "resource_limit"),
                        (["x"] * (m.MAX_VALUE_NODES + 1), "resource_limit")):
        with _raises(kind):
            m.Observation(locator=loc, record=None, presence="present", value=value, interpretation="preserved")
    a = m.Observation(locator=loc, record=None, presence="present", value={"b": 1, "a": [1.5, None]},
                      interpretation="preserved")
    b = m.Observation(locator=loc, record=None, presence="present", value={"a": [1.5, None], "b": 1},
                      interpretation="preserved")
    assert a == b and a.value.get() == {"a": [1.5, None], "b": 1}
    copy_ = a.value.get()
    copy_["a"].append(2)
    assert a.value.get() == {"a": [1.5, None], "b": 1}
    absent = m.Observation(locator=L("airr", "/Info/missing"), record=None, presence="absent", value=None,
                           interpretation="preserved")
    null = m.Observation(locator=L("airr", "/Info/missing"), record=None, presence="present", value=None,
                         interpretation="preserved")
    assert absent != null and absent.to_dict()["presence"] == "absent" and null.to_dict()["presence"] == "present"


def _nested(depth):
    value = 0
    for _ in range(depth):
        value = [value]
    return value


def test_check_inventory_against_its_lock(inventory):
    m.check_inventory(inventory, m.SourceLock.from_dict(lock_dict()))
    with _raises("digest_mismatch", "/source_lock_hash"):
        m.check_inventory(inventory, m.SourceLock.from_dict(lock_dict(origin="local_file")))
    narrow = Inv(m.SourceLock.from_dict(lock_dict(formats=("airr", "gapped", "gapped_ex")))).build()
    with _raises("identity_mismatch") as caught:
        m.check_inventory(narrow, m.SourceLock.from_dict(lock_dict(formats=("airr", "gapped", "gapped_ex"))))
    assert "airr_ex" in caught.error.message


# ── selection and decisions against the inventory ────────────────────


def _recipe(inventory, *, selected=SELECTED, decisions=(), edit=None):
    d = recipe_dict(selected=selected)
    d["decisions"] = list(decisions)
    if edit:
        edit(d)
    return m.Recipe.from_dict(d)


def test_check_selection_binds_exact_record_and_digest(inventory):
    m.check_selection(_recipe(inventory), inventory)


@pytest.mark.parametrize("edit,kind,path", [
    (lambda d: d["selection"][0].update(label="IGHV1-58*04"), "identity_mismatch", "/selection/0/label"),
    (lambda d: d["selection"][3].update(segment="J"), "identity_mismatch", "/selection/3/segment"),
    (lambda d: d["selection"][0]["record"].update(pointer=_pointer(0)), "unresolved_selection",
     "/selection/0/record"),
    (lambda d: d["selection"][0].update(sequence_sha256=selection_entry("IGHV1-58*03", "airr_ex")["sequence_sha256"]),
     "unresolved_selection", "/selection/0/sequence_sha256"),
])
def test_rehashed_selection_that_disagrees_with_its_record_is_rejected(inventory, edit, kind, path):
    with _raises(kind, path):
        m.check_selection(_recipe(inventory, edit=edit), inventory)


def test_selected_records_must_agree_with_the_source_identity(inv):
    records = list(inv.records)
    index = [i for i, r in enumerate(records) if r.locator == inv.rec["airr", "IGHJ4*02"]][0]
    fields = {k: getattr(records[index], k) for k in ("locator", "status", "problem", "label", "provider_id",
                                                       "provider_ref", "segment", "taxonomy_id", "locus")}
    for change in ({"taxonomy_id": 10090}, {"locus": "IGK"}):
        records[index] = m.SourceRecord(**{**fields, **change})
        with _raises("identity_mismatch", "/selection/4/record"):
            m.check_selection(_recipe(None), inv.build(records=records))
    records[index] = m.SourceRecord(**{**fields, "status": "malformed", "problem": "coding_sequence is not DNA"})
    with _raises("unresolved_selection", "/selection/4/record"):
        m.check_selection(_recipe(None), inv.build(records=records))


def test_check_selection_needs_the_same_lock(inventory):
    other = Inv(m.SourceLock.from_dict(lock_dict(origin="local_file"))).build()
    with _raises("digest_mismatch", "/source_lock_hash"):
        m.check_selection(_recipe(inventory), other)


def _finding(inventory, code):
    (found,) = [f for f in inventory.findings if f.code == code]
    return found


def _claim(inventory, predicate, subject):
    (found,) = [c for c in inventory.claims if c.predicate == predicate and c.subject == subject]
    return found


def test_check_decisions_accepts_current_ids_and_digests(inventory, inv):
    finding = _finding(inventory, "variant_coding_sequence_differs")
    jcand = inv.cand["airr", "IGHJ4*02"]
    j_sha = inventory.candidate(jcand).sha256
    independence = _claim(inventory, "phase_independence", jcand)
    phase = _claim(inventory, "coding_phase", jcand)
    decisions = [decision(findings=[finding.id]),
                 decision("withhold_claim", claims=[independence.id], sequence_sha256=j_sha),
                 decision("accept_claim", claims=[phase.id], sequence_sha256=j_sha,
                          evidence=[inv.obs_at(inv.rec["airr", "IGHJ4*02"], "j_codon_frame").to_dict()])]
    m.check_decisions(_recipe(inventory, decisions=decisions), inventory)


def test_stale_decisions_are_rejected(inventory, inv):
    finding = _finding(inventory, "variant_coding_sequence_differs")
    jcand = inv.cand["airr", "IGHJ4*02"]
    j_sha = inventory.candidate(jcand).sha256
    phase = _claim(inventory, "coding_phase", jcand)
    cases = [
        (decision(findings=["finding:" + "f" * 64]), "stale_decision", "/decisions/0/findings/0"),
        (decision("withhold_claim", claims=["claim:" + "c" * 64], sequence_sha256=j_sha), "stale_decision",
         "/decisions/0/claims/0"),
        (decision("withhold_claim", claims=[phase.id]), "stale_decision", "/decisions/0/sequence_sha256"),
        (decision(findings=[finding.id], evidence=[{"asset": "airr", "pointer": "/GermlineSet/0/nowhere"}]),
         "stale_decision", "/decisions/0/evidence/0"),
        (decision(findings=[_finding(inventory, "paralog_only_identity").id]), "stale_decision",
         "/decisions/0/sequence_sha256"),
        (decision(findings=[finding.id], sequence_sha256=selection_entry("IGHV1-58*03", "airr_ex")["sequence_sha256"]),
         "stale_decision", "/decisions/0/sequence_sha256"),
    ]
    for entry, kind, path in cases:
        with _raises(kind, path):
            m.check_decisions(_recipe(inventory, decisions=[entry]), inventory)


def test_accepted_claims_must_be_representable_by_the_runtime_reference(inventory, inv):
    ex = inv.cand["airr_ex", "IGHV1-58*03"]
    extension = _claim(inventory, "extension", ex)
    jcand = inv.cand["airr", "IGHJ4*02"]
    independence = _claim(inventory, "phase_independence", jcand)
    for claim, sha in ((extension, inventory.candidate(ex).sha256), (independence, inventory.candidate(jcand).sha256)):
        entry = decision("accept_claim", claims=[claim.id], sequence_sha256=sha,
                         evidence=[claim.evidence[0].to_dict()])
        with _raises("projection", "/decisions/0/claims/0"):
            m.check_decisions(_recipe(inventory, decisions=[entry]), inventory)


def test_decisions_on_unselected_sequences_are_stale(inventory, inv):
    ex = inv.cand["airr_ex", "IGHJ4*02"]
    record = inv.rec["airr_ex", "IGHJ4*02"]
    functionality = _claim(inventory, "functionality", record)
    entry = decision("withhold_claim", claims=[functionality.id], sequence_sha256=inventory.candidate(ex).sha256)
    selected = tuple(label for label in SELECTED if label != "IGHJ4*02")
    with _raises("stale_decision", "/decisions/0/claims/0"):
        m.check_decisions(_recipe(inventory, selected=selected, decisions=[entry]), inventory)


# ── build report ─────────────────────────────────────────────────────

RUNTIME = {"catalogue_hash": "sha256:" + "1" * 64, "models_hash": "sha256:" + "2" * 64,
           "execution_hash": "sha256:" + "3" * 64, "artifact_hash": "sha256:" + "4" * 64}
GENERATOR = {"name": "GenAIRR", "version": "2.3.1"}


def _capabilities():
    def status(name, available, evidence=None, unmet=(), note=None):
        return {"name": name, "available": available, "evidence": evidence,
                "unmet": [{"allele": a, "requirement": r, "hint": h} for a, r, h in unmet], "note": note}
    blocker = ("V:IGHV1-69*01", 'accepted boundary_3 == "complete"', 'add "boundary_3": {"value": "complete", ...}')
    return {"statuses": [
        status("catalogue_inspection", True), status("observed_span_simulation", True),
        status("source_annotated_junction", True, "source"),
        status("canonical_productive", True, "assertion", note="phase support: anchor_derived"),
        status("original_coding_ends", False, unmet=[blocker]),
        status("complete_reference_span", False, unmet=[blocker]),
        status("genotype_constrained", False, note="not yet supported for portable references"),
        status("duplex_family_replay", False, note="not yet supported for portable references")],
        "phase_support": "anchor_derived"}


def report_decisions(inventory, inv):
    jcand, dcand = inv.cand["airr", "IGHJ4*02"], inv.cand["airr", "IGHD1-1*01"]
    return [decision(findings=[_finding(inventory, "variant_coding_sequence_differs").id]),
            decision("withhold_claim", claims=[_claim(inventory, "boundary_3", dcand).id],
                     sequence_sha256=inventory.candidate(dcand).sha256,
                     reason="The physical 3' end of this D is not established by the source"),
            decision("accept_claim", claims=[_claim(inventory, "coding_phase", jcand).id],
                     sequence_sha256=inventory.candidate(jcand).sha256,
                     evidence=[inv.obs_at(inv.rec["airr", "IGHJ4*02"], "j_codon_frame").to_dict()],
                     reason="j_codon_frame - 1 is the reported phase under the pinned producer profile")]


def report_dict(inventory, inv, *, decisions=None) -> dict:
    lock = m.SourceLock.from_dict(lock_dict())
    recipe = _recipe(inventory, decisions=report_decisions(inventory, inv) if decisions is None else decisions)
    selected = {e.record.id: e for e in recipe.selection}
    chosen = {}
    for e in recipe.selection:
        (c,) = [c for c in inventory.candidates_of(e.record) if c.sha256 == e.sequence_sha256]
        chosen[e.record.id] = c
    by_label = {(e.segment, e.label): e for e in recipe.selection}
    records = []
    for r in inventory.records:
        if r.id in selected:
            records.append({"record": r.locator.to_dict(), "outcome": "selected", "allele": selected[r.id].allele_id,
                            "exclusion": None, "related": None, "detail": None})
            continue
        if r.status != "parsed":
            records.append({"record": r.locator.to_dict(), "outcome": "rejected", "allele": None,
                            "exclusion": None, "related": None, "detail": r.problem})
            continue
        entry = by_label.get((r.segment, r.label))
        if entry is None:
            records.append({"record": r.locator.to_dict(), "outcome": "excluded", "allele": None,
                            "exclusion": "outside_explicit_selection", "related": None, "detail": None})
            continue
        same = any(c.sha256 == entry.sequence_sha256 for c in inventory.candidates_of(r.locator))
        records.append({"record": r.locator.to_dict(), "outcome": "excluded", "allele": None,
                        "exclusion": "duplicate_occurrence" if same else "alternative_for_selected_identity",
                        "related": entry.record.to_dict(), "detail": None})
    decided = {cid: i for i, d in enumerate(recipe.decisions) for cid in d.claims}
    selected_subjects = set(selected) | {c.id for c in chosen.values()}
    claims = []
    for c in inventory.claims:
        if c.subject.id not in selected_subjects:
            outcome, reason = "not_selected", None
        elif c.predicate not in m.RUNTIME_PREDICATES:
            outcome, reason = "not_representable", "the runtime reference has no such claim"
        elif c.id in decided and recipe.decisions[decided[c.id]].action == "withhold_claim":
            outcome, reason = "withheld", recipe.decisions[decided[c.id]].reason
        elif c.status == "unresolved":
            outcome, reason = "withheld", "unresolved"
        else:
            outcome, reason = "projected", None
        claims.append({"claim": c.id, "outcome": outcome, "decision": decided.get(c.id), "reason": reason})
    exports = [{"asset": a, **inventory.export_counts(a)} for a in ("airr", "airr_ex")]
    seg = {"V": 0, "D": 0, "J": 0, "unassigned": 0}
    for e in recipe.selection:
        seg[e.segment] += 1
    outcomes = {k: sum(r["outcome"] == k for r in records) for k in ("selected", "excluded", "rejected")}
    return {
        "artifact": "genairr_reference_source_build_report", "schema_version": 1,
        "generator": dict(GENERATOR), "source_lock": lock.to_dict(), "recipe": recipe.to_dict(),
        "inventory": inventory.to_dict(),
        "hashes": {"source_lock": m.lock_hash(lock), "inventory": m.inventory_hash(inventory),
                   "recipe": m.recipe_hash(recipe), "build_inputs": m.build_input_hash(lock, inventory, recipe)},
        "runtime": dict(RUNTIME),
        "comparison": {"complete": True, "exports": ["airr", "airr_ex"], "missing": [],
                       "findings": [_finding(inventory, "variant_coding_sequence_differs").id]},
        "counts": {"exports": exports, "selected_alleles": seg, "record_outcomes": outcomes},
        "selection": [{"allele": e.allele_id, "candidate": chosen[e.record.id].locator.to_dict()}
                      for e in recipe.selection],
        "records": records, "claims": claims, "capabilities": _capabilities(),
        "requirements": [{"capability": "observed_span_simulation", "evidence": "assertion", "satisfied": True}],
    }


@pytest.fixture(scope="module")
def report_d(inventory, inv):
    return report_dict(inventory, inv)


@pytest.fixture(scope="module")
def report(report_d):
    return m.BuildReport.from_dict(report_d)


def test_report_round_trips_and_is_byte_stable(report, report_d):
    assert report.to_dict() == report_d
    text = canonical_json(report_d)
    assert canonical_json(m.BuildReport.from_json(text).to_dict()) == text
    assert m.report_hash(report) == typed_digest(report_d, "genairr_reference_source_build_report_v1")


def test_report_exposes_every_record_claim_and_distinct_count(report):
    d = report.to_dict()
    assert len(d["records"]) == len(d["inventory"]["records"])
    assert len(d["claims"]) == len(d["inventory"]["claims"])
    airr_ex = d["counts"]["exports"][1]
    assert airr_ex["record_occurrences"]["V"] == 2 and airr_ex["mentioned_identities"]["V"] == 1
    assert d["counts"]["selected_alleles"] == {"V": 3, "D": 1, "J": 1, "unassigned": 0}
    assert d["counts"]["record_outcomes"] == {"selected": 5, "excluded": 4, "rejected": 0}
    outcomes = {c["outcome"] for c in d["claims"]}
    assert outcomes == {"projected", "withheld", "not_representable", "not_selected"}
    assert {r["exclusion"] for r in d["records"] if r["outcome"] == "excluded"} == {
        "alternative_for_selected_identity", "duplicate_occurrence"}


def _edited_report(report_d, edit):
    d = copy.deepcopy(report_d)
    edit(d)
    return d


def _record_index(d, label, asset):
    for i, r in enumerate(d["records"]):
        rec = [x for x in d["inventory"]["records"] if x["locator"] == r["record"]][0]
        if rec["label"] == label and r["record"]["asset"] == asset:
            return i
    raise AssertionError(label)


def _claim_outcome_index(d, outcome):
    return [i for i, c in enumerate(d["claims"]) if c["outcome"] == outcome][0]


@pytest.mark.parametrize("edit,kind,path", [
    (lambda d: d["hashes"].update(build_inputs="sha256:" + "9" * 64), "digest_mismatch", "/hashes/build_inputs"),
    (lambda d: d["hashes"].update(inventory="sha256:" + "9" * 64), "digest_mismatch", "/hashes/inventory"),
    (lambda d: d["recipe"].update(source_lock_hash="sha256:" + "9" * 64), "digest_mismatch", "/recipe/source_lock_hash"),
    (lambda d: d["recipe"]["reference"].update(locus="TRB"), "identity_mismatch", "/recipe/reference/locus"),
    (lambda d: d["runtime"].update(execution_hash="sha256:XYZ"), "malformed", "/runtime/execution_hash"),
    (lambda d: d["generator"].update(name="NotGenAIRR"), "malformed", "/generator/name"),
    (lambda d: d.update(schema_version=2), "unsupported_schema", "/schema_version"),
    (lambda d: d.update(report_hash="sha256:" + "0" * 64), "malformed", ""),
    (lambda d: d["comparison"].update(complete=False), "malformed", "/comparison/complete"),
    (lambda d: d["comparison"].update(findings=["finding:" + "0" * 64]), "malformed", "/comparison/findings/0"),
    (lambda d: d["counts"]["exports"][0]["record_occurrences"].update(V=209), "malformed", "/counts/exports"),
    (lambda d: d["counts"]["selected_alleles"].update(V=190), "malformed", "/counts/selected_alleles"),
    (lambda d: d["counts"]["record_outcomes"].update(rejected=1), "malformed", "/counts/record_outcomes"),
    (lambda d: d["counts"]["selected_alleles"].update(V=True), "malformed", "/counts/selected_alleles/V"),
    (lambda d: d["counts"]["selected_alleles"].update(V=-1), "malformed", "/counts/selected_alleles/V"),
    (lambda d: d["selection"].pop(), "malformed", "/selection"),
    (lambda d: d["selection"][0].update(allele="V:IGHV1-58*04"), "malformed", "/selection/0/allele"),
    (lambda d: d["selection"][0].update(candidate=d["selection"][1]["candidate"]), "malformed",
     "/selection/0/candidate"),
    (lambda d: d["records"].reverse(), "malformed", "/records/0/record"),
    (lambda d: d["records"][_record_index(d, "IGHV1-58*03", "airr")].update(
        outcome="excluded", allele=None, exclusion="outside_explicit_selection"), "malformed", "/records/"),
    (lambda d: d["records"][_record_index(d, "IGHV1-58*03", "airr_ex")].update(exclusion="duplicate_occurrence"),
     "malformed", "/records/"),
    (lambda d: d["records"][_record_index(d, "IGHV1-69*01", "airr_ex")].update(
        exclusion="alternative_for_selected_identity"), "malformed", "/records/"),
    (lambda d: d["records"][_record_index(d, "IGHV1-69*01", "airr_ex")].update(related=None), "malformed",
     "/records/"),
    (lambda d: d["records"][_record_index(d, "IGHV1-69*01", "airr_ex")].update(exclusion="looked_odd"), "malformed",
     "/records/"),
    (lambda d: d["records"][_record_index(d, "IGHV1-69*01", "airr_ex")].update(outcome="rejected", exclusion=None,
                                                                              related=None), "malformed", "/records/"),
    (lambda d: d["claims"][_claim_outcome_index(d, "not_representable")].update(outcome="projected"), "malformed",
     "/claims/"),
    (lambda d: d["claims"][_claim_outcome_index(d, "not_selected")].update(outcome="projected"), "malformed",
     "/claims/"),
    (lambda d: d["claims"][_claim_outcome_index(d, "projected")].update(outcome="not_selected"), "malformed",
     "/claims/"),
    (lambda d: d["claims"][_claim_outcome_index(d, "withheld")].update(reason=None), "malformed", "/claims/"),
    (lambda d: d["claims"][_claim_outcome_index(d, "withheld")].update(outcome="projected", reason=None), "malformed",
     "/claims/"),
    (lambda d: d["claims"][_claim_outcome_index(d, "projected")].update(decision=0), "malformed", "/claims/"),
    (lambda d: d["claims"][_claim_outcome_index(d, "projected")].update(decision=7), "malformed", "/claims/"),
    (lambda d: d["claims"][_claim_outcome_index(d, "projected")].update(decision=True), "malformed", "/claims/"),
    (lambda d: d["claims"].reverse(), "malformed", "/claims/0/claim"),
    (lambda d: d["capabilities"]["statuses"].reverse(), "malformed", "/capabilities/statuses"),
    (lambda d: d["capabilities"]["statuses"][4].update(available=True), "malformed", "/capabilities/statuses/4"),
    (lambda d: d["capabilities"]["statuses"][6].update(evidence="source"), "malformed",
     "/capabilities/statuses/6/evidence"),
    (lambda d: d["capabilities"].update(phase_support="guessed"), "malformed", "/capabilities/phase_support"),
    (lambda d: d["requirements"][0].update(satisfied=False), "malformed", "/requirements/0/satisfied"),
    (lambda d: d["requirements"][0].update(evidence="source"), "malformed", "/requirements/0"),
    (lambda d: d["requirements"].append({"capability": "original_coding_ends", "evidence": "assertion",
                                         "satisfied": False}), "malformed", "/requirements"),
])
def test_report_rejects_inconsistent_content(report_d, edit, kind, path):
    with _raises(kind) as caught:
        m.BuildReport.from_dict(_edited_report(report_d, edit))
    assert caught.error.path.startswith(path), caught.error


def test_report_withheld_and_accepted_claims_follow_the_decision_ledger(report_d):
    d = report_d
    decided = {c["claim"]: c for c in d["claims"] if c["decision"] is not None}
    actions = [x["action"] for x in d["recipe"]["decisions"]]
    assert sorted((actions[c["decision"]], c["outcome"]) for c in decided.values()) == [
        ("accept_claim", "projected"), ("withhold_claim", "withheld")]
    bad = copy.deepcopy(d)
    (index,) = [i for i, c in enumerate(bad["claims"]) if c["decision"] == 1]
    bad["claims"][index].update(decision=None, reason="unresolved")
    with _raises("malformed", "/claims"):
        m.BuildReport.from_dict(bad)


def test_report_records_rejected_source_records_apart_from_exclusions(inventory, inv):
    loc = L("airr", _pointer(0))
    odd = m.SourceRecord(locator=loc, status="malformed", problem="coding_sequence is not a string", label=None,
                         provider_id=None, provider_ref=None, segment=None, taxonomy_id=None, locus=None)
    _, raw = _find("IGHV1-2*02")
    other = m.SourceRecord(locator=L("airr", _pointer(_find("IGHV1-2*02")[0])), status="parsed", problem=None,
                           label=raw["label"], provider_id=raw["allele_description_id"],
                           provider_ref=raw["allele_description_ref"], segment="V", taxonomy_id=9606, locus="IGH")
    d = report_dict(inv.build(records=inv.records + [odd, other]), inv, decisions=[])
    report = m.BuildReport.from_dict(d)
    by_record = {r.record: r for r in report.records}
    assert by_record[loc].outcome == "rejected" and by_record[loc].detail == "coding_sequence is not a string"
    assert by_record[other.locator].exclusion == "outside_explicit_selection"
    assert report.counts.record_outcomes.to_dict() == {"selected": 5, "excluded": 5, "rejected": 1}
    (index,) = [i for i, r in enumerate(d["records"]) if r["record"] == loc.to_dict()]
    d["records"][index].update(outcome="excluded", exclusion="outside_explicit_selection")
    d["counts"]["record_outcomes"] = {"selected": 5, "excluded": 6, "rejected": 0}
    with _raises("malformed", f"/records/{index}/outcome"):
        m.BuildReport.from_dict(d)


def test_report_incomplete_comparison_lists_the_missing_export():
    lock = m.SourceLock.from_dict(lock_dict(formats=("airr", "gapped")))
    comparison = m.Comparison(complete=False, exports=("airr",), missing=("airr_ex",), findings=())
    comparison.check(lock, None)
    with _raises("malformed", "/complete"):
        m.Comparison(complete=True, exports=("airr",), missing=("airr_ex",), findings=())


# ── dossier manifest ─────────────────────────────────────────────────

FIXED = ("inputs/annotations.json", "inputs/d.fasta", "inputs/j.fasta", "inputs/v.fasta", "recipe.json",
         "reference.genairr-ref.json", "report.json", "source.lock.json")


def manifest_dict(report, *, d_fasta=True) -> dict:
    members = [{"name": name, "sha256": "sha256:" + format(i, "064x"), "bytes": 100 + i}
               for i, name in enumerate(FIXED) if d_fasta or name != "inputs/d.fasta"]
    for asset in report.source_lock.assets:
        members.append({"name": "cache/blobs/sha256/" + asset.blob_name, "sha256": asset.sha256, "bytes": asset.bytes})
    members.sort(key=lambda x: x["name"])
    hashes = report.to_dict()["hashes"]
    return {"artifact": "genairr_reference_source_build_manifest", "schema_version": 1,
            "builder": dict(GENERATOR),
            "versions": {"source_lock": 1, "recipe": 1, "report": 1, "reference": 1, "annotations": 1},
            "adapter": copy.deepcopy(ADAPTER), "projection": {"id": "observed_span", "version": 1},
            "hashes": {**hashes, "report": m.report_hash(report)}, "runtime": dict(RUNTIME), "members": members}


def test_manifest_round_trips_and_binds_the_report(report):
    d = manifest_dict(report)
    manifest = m.BuildManifest.from_dict(d)
    assert manifest.to_dict() == d
    assert canonical_json(m.BuildManifest.from_json(canonical_json(d)).to_dict()) == canonical_json(d)
    assert m.manifest_hash(manifest) == typed_digest(d, "genairr_reference_source_build_manifest_v1")
    m.check_manifest(manifest, report.source_lock, report)


@pytest.mark.parametrize("name", ["", "/report.json", "../report.json", "inputs/../report.json", "inputs\\v.fasta",
                                  "C:report.json", "c:/report.json", "inputs//v.fasta", "./report.json",
                                  "inputs/v.fasta/", "x" * 256, "receipts/2026-10-22T09:58:15Z.json",
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


# ── the digest graph ─────────────────────────────────────────────────


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


# ── strictness across every owned format ─────────────────────────────

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
    value = inventory.observations[0].value.get()
    if isinstance(value, str):
        assert inventory.observations[0].value.get() == value
    with pytest.raises(AttributeError):
        report.records = ()
    with pytest.raises(AttributeError):
        inventory.claims[0].status = "accepted"




# ── shared resource budget and the interpretation-profile scope ──────


def test_resource_budget_is_shared_by_every_format():
    assert (m.MAX_ASSET_BYTES, m.MAX_SOURCE_PAYLOAD_BYTES, m.MAX_ASSETS) == (16 * 1024 * 1024, 128 * 1024 * 1024, 16)
    assert (m.MAX_SOURCE_RECORDS, m.MAX_JSON_DEPTH, m.MAX_SEQUENCE) == (100_000, 64, 65_535)
    from GenAIRR._reference_sources import inventory, lock, receipt
    assert lock.MAX_ASSET_BYTES is m.MAX_ASSET_BYTES and receipt.MAX_ASSET_BYTES is m.MAX_ASSET_BYTES
    assert inventory.MAX_RECORDS == m.MAX_SOURCE_RECORDS
    loc = L("airr", "/Info/title")
    m.Observation(locator=loc, record=None, presence="present", value=_nested(m.MAX_JSON_DEPTH),
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
    assert m.PROFILE_SCOPES == {("ogrdb_api_v2_human_ig_observed_v1", 1): (
        (9606, "9606.IGH_VDJ", 10, ("airr", "airr_ex", "gapped", "gapped_ex")),
        (9606, "9606.IGKappa_VJ", 5, ("airr", "airr_ex", "gapped", "gapped_ex")),
        (9606, "9606.IGLambda_VJ", 4, ("airr", "airr_ex", "gapped", "gapped_ex")))}
    for set_id in (IGH, "9606.IGKappa_VJ", "9606.IGLambda_VJ"):
        m.check_profile_scope(m.SourceLock.from_dict(lock_dict(set_id)))


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



# ── content-derived ids are bound to the lock ────────────────────────


def test_content_ids_change_with_the_lock_even_for_identical_content(inventory):
    other_lock = m.SourceLock.from_dict(lock_dict(origin="local_file"))
    other = Inv(other_lock).build()
    assert [r.locator for r in other.records] == [r.locator for r in inventory.records]
    for name in ("claims", "mappings", "findings"):
        mine, theirs = getattr(inventory, name), getattr(other, name)
        assert len(mine) == len(theirs)
        assert not {x.id for x in mine} & {x.id for x in theirs}, name
        assert {x.source_lock_hash for x in theirs} == {other_lock.commitment}
    body = {k: v for k, v in inventory.claims[0].to_dict().items() if k != "id"}
    assert body["source_lock_hash"] == inventory.source_lock_hash
    assert inventory.claims[0].id == "claim:" + typed_digest(body, "genairr_reference_source_claim_id_v1")[7:]


def test_items_from_another_lock_cannot_join_an_inventory(inv, inventory):
    other = Inv(m.SourceLock.from_dict(lock_dict(origin="local_file")))
    for name in ("claims", "findings"):
        with _raises("digest_mismatch") as caught:
            inv.build(**{name: list(getattr(inv, name))[1:] + [getattr(other, name)[0]]})
        assert caught.error.path.startswith(f"/{name}/") and caught.error.path.endswith("/source_lock_hash")


def test_decisions_go_stale_when_the_recipe_is_pointed_at_another_lock(inventory, inv):
    decisions = report_decisions(inventory, inv)
    lock = m.SourceLock.from_dict(lock_dict())
    m.check_build_inputs(lock, inventory, _recipe(inventory, decisions=decisions))
    other_lock = m.SourceLock.from_dict(lock_dict(origin="local_file"))
    other_inventory = Inv(other_lock).build()
    repointed = _recipe(inventory, decisions=decisions,
                        edit=lambda d: d.update(source_lock_hash=other_lock.commitment))
    with _raises("stale_decision", "/decisions/0/findings/0"):
        m.check_build_inputs(other_lock, other_inventory, repointed)



# ── mappings are verified against the observed values ────────────────


def _replace_mapping(inv, index, **changes):
    """Rebuild the inventory with mapping ``index`` changed (claims follow its new id)."""
    import dataclasses
    old = inv.mappings[index]
    new = dataclasses.replace(old, **changes)
    claims = [dataclasses.replace(c, mappings=tuple(sorted(new.id if x == old.id else x for x in c.mappings)))
              for c in inv.claims]
    mappings = list(inv.mappings)
    mappings[index] = new
    return inv.build(mappings=mappings, claims=claims)


def test_alignment_blocks_are_the_exact_runs_of_ungapped_positions():
    gapped = _gapped_index("IGHV1-58*03", "gapped")[1]
    assert m.alignment_blocks(gapped) == _blocks(gapped)
    assert m.alignment_blocks("AC..G.T") == (m.Block((0, 2), (0, 2)), m.Block((2, 3), (4, 5)), m.Block((3, 4), (6, 7)))
    assert m.IMGT_GAP == "."


def test_fasta_record_observations_have_a_frozen_value_shape():
    loc = L("gapped", "/3")
    m.Observation(locator=loc, record=None, presence="present", value={"header": "IGHV1-2*02", "sequence": "CAG..GTG"},
                  interpretation="interpreted")
    for value in ("CAG..GTG", {"header": "IGHV1-2*02"}, {"header": "IGHV1-2*02", "sequence": 7},
                  {"header": None, "sequence": "CAG"}, {"header": "x", "sequence": "CAG", "description": ""}):
        with _raises("malformed", "/value"):
            m.Observation(locator=loc, record=None, presence="present", value=value, interpretation="interpreted")
    with _raises("malformed", "/presence"):
        m.Observation(locator=loc, record=None, presence="absent", value=None, interpretation="interpreted")


def test_gapped_mapping_output_digest_must_match_the_observed_gapped_record(inv):
    with _raises("digest_mismatch") as caught:
        _replace_mapping(inv, 0, output_sha256="sha256:" + "a" * 64)
    assert caught.error.path.startswith("/mappings/") and caught.error.path.endswith("/output_sha256")


def test_exact_gapped_blocks_must_reproduce_the_candidate(inv):
    n = inv.mappings[0].blocks[-1].input[1]
    with _raises("malformed") as caught:
        _replace_mapping(inv, 0, blocks=(m.Block(input=(0, n), output=(1, n + 1)),))
    assert caught.error.path.endswith("/blocks")
    with _raises("malformed") as caught:
        _replace_mapping(inv, 0, blocks=(), validation="mismatch")
    assert caught.error.path.endswith("/validation")


def test_coding_interval_mapping_is_verified_against_the_record_field(inv):
    index, raw = _find("IGHJ4*02")
    rec, cand = inv.rec["airr", "IGHJ4*02"], inv.cand["airr", "IGHJ4*02"]
    start, end = raw["gene_start"] - 1, raw["gene_end"]
    seq_loc = L("airr", rec.pointer + "/sequence")
    observation = m.Observation(locator=seq_loc, record=rec, presence="present", value=raw["sequence"],
                                interpretation="interpreted")

    def mapping(**changes):
        fields = dict(source_lock_hash=inv.lock.commitment, kind="coding_interval", input=seq_loc, output=cand,
                      input_sha256=m.text_digest(raw["sequence"]), output_sha256=m.sequence_digest(raw["coding_sequence"]),
                      input_space="record_sequence", output_space="coding_sequence",
                      blocks=(m.Block((start, end), (0, end - start)),), transform="extract_interval",
                      validation="exact")
        fields.update(changes)
        return inv.build(observations=inv.observations + [observation], mappings=inv.mappings + [m.Mapping(**fields)])

    mapping()
    with _raises("digest_mismatch") as caught:
        mapping(input_sha256="sha256:" + "b" * 64)
    assert caught.error.path.endswith("/input_sha256")
    with _raises("malformed") as caught:
        mapping(blocks=(m.Block((start + 1, end + 1), (0, end - start)),))
    assert caught.error.path.endswith("/blocks")


def test_record_fields_are_observed_with_their_record(inv):
    rec = inv.rec["airr", "IGHV1-69*01"]
    orphan = m.Observation(locator=L("airr", rec.pointer + "/curation"), record=None, presence="present", value="x",
                           interpretation="preserved")
    with _raises("malformed") as caught:
        inv.build(observations=inv.observations + [orphan])
    assert caught.error.path.startswith("/observations/") and caught.error.path.endswith("/record")


def test_mention_locators_are_observed_fields(inv):
    ex69 = inv.rec["airr_ex", "IGHV1-69*01"]
    unobserved = m.MentionedIdentity(asset="airr_ex", segment="V", label="IGHV1-69D*01", mentions=(
        m.Mention(record=ex69, locator=L("airr_ex", ex69.pointer + "/paralogs/1"), relationship="paralog"),))
    with _raises("malformed", "/mentioned_identities/0/mentions/0/locator"):
        inv.build(mentioned_identities=[unobserved])


def test_alignment_derived_claims_rest_on_exact_mappings(inv):
    cand, rec = inv.cand["airr", "IGHV1-69*01"], inv.rec["airr", "IGHV1-69*01"]
    gindex, gapped = _gapped_index("IGHV1-2*02", "gapped")
    gloc = L("gapped", f"/{gindex}")
    observation = m.Observation(locator=gloc, record=None, presence="present",
                                value={"header": "IGHV1-2*02", "sequence": gapped}, interpretation="interpreted")
    mismatch = m.Mapping(source_lock_hash=inv.lock.commitment, kind="gapped_alignment", input=cand, output=gloc,
                         input_sha256=m.sequence_digest(_find("IGHV1-69*01")[1]["coding_sequence"]),
                         output_sha256=m.text_digest(gapped), input_space="coding_sequence",
                         output_space="imgt_gapped", blocks=(), transform="remove_imgt_gaps", validation="mismatch")
    fields = dict(observations=inv.observations + [observation], mappings=inv.mappings + [mismatch])
    inv.build(**fields)
    cdr3 = inv.obs_at(rec, "v_gene_delineations/0/cdr3_start")
    claim = inv.claim(cand, "anchor", 285, [cdr3, gloc], system="coding_sequence", derivation="imgt_alignment",
                      mappings=[mismatch.id])
    with _raises("malformed") as caught:
        inv.build(claims=inv.claims + [claim], **fields)
    assert caught.error.path.endswith("/mappings/0") and "mismatch" in caught.error.message


# ── identical locked bytes ───────────────────────────────────────────


def _shared_bytes_lock(bytes_delta=0):
    d = lock_dict()
    (gapped,) = [a for a in d["assets"] if a["format"] == "gapped"]
    for asset in d["assets"]:
        if asset["format"] == "gapped_ex":
            asset["sha256"], asset["bytes"] = gapped["sha256"], gapped["bytes"] + bytes_delta
    return _reseal(d)


def test_one_digest_has_one_size_in_a_lock():
    m.SourceLock.from_dict(_shared_bytes_lock())
    with _raises("malformed", "/assets") as caught:
        m.SourceLock.from_dict(_shared_bytes_lock(bytes_delta=1))
    assert "gapped" in caught.error.message


def test_manifest_holds_one_blob_for_assets_with_identical_bytes(report):
    lock = m.SourceLock.from_dict(_shared_bytes_lock())
    d = manifest_dict(report)
    d["hashes"]["source_lock"] = lock.commitment
    names = {"cache/blobs/sha256/" + a.blob_name: a for a in lock.assets}
    d["members"] = [x for x in d["members"] if not x["name"].startswith("cache/")]
    d["members"] += [{"name": n, "sha256": a.sha256, "bytes": a.bytes} for n, a in names.items()]
    d["members"].sort(key=lambda x: x["name"])
    assert len(names) == 3
    m.check_manifest(m.BuildManifest.from_dict(d), lock)
    shared = "cache/blobs/sha256/" + lock.asset("gapped").blob_name
    d["members"] = [x for x in d["members"] if x["name"] != shared]
    with _raises("missing_asset", "/members") as caught:
        m.check_manifest(m.BuildManifest.from_dict(d), lock)
    assert "gapped" in caught.error.message
