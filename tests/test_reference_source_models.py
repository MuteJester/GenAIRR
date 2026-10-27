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
           "coordinate_profile": {"id": "ogrdb_api_v2_coordinates", "version": 1}}


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
    (lambda d: d["adapter"]["coordinate_profile"].update(version=2), "unsupported_profile", "/adapter"),
    (lambda d: d["adapter"]["coordinate_profile"].update(id="imgt"), "unsupported_profile", "/adapter"),
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
         "evidence": [{"asset": "airr", "pointer": entry["record"]["pointer"] + "/gene_end"}]}
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
                        ("ogrdb", (("ogrdb_api_v2_coordinates", 1),)))
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
