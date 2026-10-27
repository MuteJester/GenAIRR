"""OGRDB germline-set adapter: set identity, normalized inventory, accounting and comparison.

The adapter reads the locked bytes of OGRDB API v2 exports.  It validates the set
identity from the payload envelope, keeps every record occurrence (malformed and
unsupported ones included), every supplied coding sequence, every relationship-only
name and every provider field with its locator, and compares the two record exports
label by label without preferring either.  It interprets no coordinates.

Expected values come from the independently authored fixture files
(``expectations.json``, ``human_igh_release10_candidate_selection.json``,
``sources.json``), from literals, or from counting raw payload fields with the standard
library JSON reader; never from the adapter under test.  Provider payloads handed to
the adapter are parsed only by its strict loader.
"""
from __future__ import annotations

import copy
import hashlib
import json
from collections import Counter
from pathlib import Path

import pytest

from GenAIRR._reference_sources import model as m
from GenAIRR._reference_sources import ogrdb
from GenAIRR._reference_sources.acquisition import ingest_files, read_locked
from GenAIRR._reference_sources.errors import ReferenceSourceError

ROOT = Path(__file__).resolve().parent / "data" / "reference_sources"
SOURCES = json.loads((ROOT / "sources.json").read_bytes().decode("utf-8"))
EXPECT = json.loads((ROOT / "expectations.json").read_bytes().decode("utf-8"))
SELECTION = json.loads((ROOT / "human_igh_release10_candidate_selection.json").read_bytes().decode("utf-8"))

IGH, IGK, IGL = ("9606.IGH_VDJ", 10), ("9606.IGKappa_VJ", 5), ("9606.IGLambda_VJ", 4)
SETS = (IGH, IGK, IGL)
EXPORTS = ("airr", "airr_ex")
ALL = ("airr", "airr_ex", "gapped", "gapped_ex")
SUFFIX = {"airr": "json", "airr_ex": "json", "gapped": "fasta", "gapped_ex": "fasta"}
CHAIN = {"IGH": "vdj", "IGK": "vj", "IGL": "vj"}
RECORDS = "/GermlineSet/0/allele_descriptions"
ADAPTER = m.AdapterRef("ogrdb_germline_set", 1, m.CoordinateProfileRef("ogrdb_api_v2_human_ig_observed_v1", 1))


# ── helpers ──────────────────────────────────────────────────────────


def _key(source) -> str:
    return f"{source[0]}/{source[1]}"


def _entry(source) -> dict:
    (entry,) = [s for s in SOURCES["sets"] if (s["set_id"], s["release"]) == source]
    return entry


def raw(source, fmt) -> bytes:
    set_id, release = source
    return (ROOT / "ogrdb" / set_id / str(release) / f"{fmt}.{SUFFIX[fmt]}").read_bytes()


def payload(source, fmt) -> dict:
    """The raw export, read independently with the standard library."""
    return json.loads(raw(source, fmt).decode("utf-8"))


def raw_records(source, fmt) -> list:
    return payload(source, fmt)["GermlineSet"][0]["allele_descriptions"]


def encode(document) -> bytes:
    return json.dumps(document, ensure_ascii=False).encode("utf-8")


def _sha(data: bytes) -> str:
    return "sha256:" + hashlib.sha256(data).hexdigest()


def expected_identity(source) -> m.SourceIdentity:
    entry = _entry(source)
    taxonomy_id = int(entry["species"]["id"].split(":")[1])
    return m.SourceIdentity(m.Taxon(entry["species"]["label"], taxonomy_id), entry["locus"], CHAIN[entry["locus"]])


def seal(source, blobs: dict, *, identity=None, set_id=None, release=None) -> m.SourceLock:
    """A local-file lock over exactly ``blobs`` (format -> bytes)."""
    assets = tuple(m.Asset(asset_id=f, role=m.FORMAT_ROLES[f], format=f, sha256=_sha(blobs[f]),
                           bytes=len(blobs[f]), declared_url=None, origin_kind="local_file")
                   for f in m.FORMATS if f in blobs)
    missing = tuple(m.MissingAsset(format=f, role=m.FORMAT_ROLES[f], reason="not_supplied", detail=None)
                    for f in m.EVIDENCE_FORMATS if f not in blobs)
    terms = m.Terms(status="unresolved", assets=(), licence=None, citation=None,
                    reason="no data-terms statement was supplied")
    return m.SourceLock.seal(provider="ogrdb", set_id=set_id or source[0], release=release or source[1],
                             identity=identity or expected_identity(source), adapter=ADAPTER, assets=assets,
                             missing_assets=missing, terms=terms)


def pinned(source, formats=ALL) -> dict:
    return {f: raw(source, f) for f in formats}


_INVENTORIES: dict = {}


def inventory_of(source, formats=ALL):
    """(lock, inventory) over the pinned bytes, built once per module run."""
    key = (source, tuple(formats))
    if key not in _INVENTORIES:
        blobs = pinned(source, formats)
        lock = seal(source, blobs)
        _INVENTORIES[key] = (lock, ogrdb.build_inventory(lock, blobs))
    return _INVENTORIES[key]


def build(source, blobs):
    lock = seal(source, blobs)
    return lock, ogrdb.build_inventory(lock, blobs)


def error(kind, call, *args, **kwargs) -> ReferenceSourceError:
    with pytest.raises(ReferenceSourceError) as info:
        call(*args, **kwargs)
    assert info.value.kind == kind, f"expected {kind!r}, got {info.value.kind!r}: {info.value}"
    return info.value


def identify(source, blobs, **overrides):
    request = {"set_id": source[0], "release": source[1], "assets": blobs}
    request.update(overrides)
    return ogrdb.identify(**request)


def edited(source, fmt, edit) -> bytes:
    document = payload(source, fmt)
    edit(document)
    return encode(document)


def rec(index) -> str:
    return f"{RECORDS}/{index}"


def records_of(inventory, asset) -> list:
    return [r for r in inventory.records if r.locator.asset == asset]


def by_label(inventory, asset) -> dict:
    found: dict = {}
    for record in records_of(inventory, asset):
        found.setdefault(record.label, []).append(record)
    return found


def findings(inventory, code) -> list:
    return [f for f in inventory.findings if f.code == code]


def observed(inventory, asset, pointer):
    return inventory.observation(m.Locator(asset, pointer))


# ── set identity (usable as the acquisition identity reader) ────────


@pytest.mark.parametrize("source", SETS, ids=_key)
def test_identity_is_read_from_the_payload_envelope(source):
    assert identify(source, pinned(source)) == expected_identity(source)
    for fmt in EXPORTS:
        assert identify(source, {fmt: raw(source, fmt)}) == expected_identity(source)


def test_identity_literals_of_the_pinned_sets():
    identity = identify(IGL, pinned(IGL))
    assert identity.to_dict() == {"taxon": {"scientific_name": "Homo sapiens", "taxonomy_id": 9606},
                                  "locus": "IGL", "chain_type": "vj"}
    assert identify(IGH, pinned(IGH)).chain_type == "vdj"


def test_identity_reader_plugs_into_ingestion_and_offline_reading(tmp_path):
    files = {}
    for fmt in ALL:
        path = tmp_path / f"{fmt}.{SUFFIX[fmt]}"
        path.write_bytes(raw(IGL, fmt))
        files[fmt] = path
    result = ingest_files(set_id=IGL[0], release=IGL[1], profile="ogrdb_api_v2_human_ig_observed_v1",
                          files=files, cache=tmp_path / "cache", identify=ogrdb.identify)
    assert result.lock.identity == expected_identity(IGL)
    again = read_locked(result.lock, cache=tmp_path / "cache", identify=ogrdb.identify)
    assert again.lock == result.lock


@pytest.mark.parametrize("override", [
    {"set_id": "9606.IGKappa_VJ"},          # another set's name
    {"set_id": "10090.IGLambda_VJ"},        # another taxon
    {"set_id": "IGLambda_VJ"},              # no taxonomy prefix to agree with
    {"release": 5},                         # another release
], ids=["set", "taxon", "unprefixed", "release"])
def test_identity_is_checked_against_the_requested_set_and_release(override):
    error_ = error("identity_mismatch", identify, IGL, pinned(IGL), **override)
    assert error_.path.startswith("/assets/airr")


def test_identity_does_not_trust_the_requested_set():
    # IGH bytes requested under the IGL set id and release: the payload decides.
    error("identity_mismatch", identify, IGL, pinned(IGH))


@pytest.mark.parametrize("edit", [
    lambda d: d["GermlineSet"][0].update(locus="IGK"),
    lambda d: d["GermlineSet"][0]["species"].update(id="NCBITAXON:10090"),
    lambda d: d["GermlineSet"][0]["species"].update(label="Mus musculus"),
    lambda d: d["GermlineSet"][0].update(germline_set_id="OGRDB:G99999"),
], ids=["locus", "taxonomy", "scientific-name", "set-id"])
def test_every_record_export_must_agree(edit):
    blobs = {"airr": raw(IGL, "airr"), "airr_ex": edited(IGL, "airr_ex", edit)}
    found = error("identity_mismatch", identify, IGL, blobs)
    assert "airr_ex" in found.path


def test_exports_of_two_different_sets_disagree():
    blobs = {"airr": raw(IGK, "airr"), "airr_ex": raw(IGL, "airr_ex")}
    error("identity_mismatch", identify, IGL, blobs)


@pytest.mark.parametrize("edit,ok", [
    (lambda d: d["GermlineSet"][0].update(germline_set_ref=None), True),
    (lambda d: d["GermlineSet"][0].update(release_version=None), True),
    (lambda d: d["GermlineSet"][0].update(germline_set_name=None), True),
    (lambda d: d["GermlineSet"][0].pop("germline_set_ref"), True),
    (lambda d: d["GermlineSet"][0].update(release_version=4), True),
    (lambda d: d["GermlineSet"][0].update(germline_set_ref="OGRDB:IGLambda_VJ.5"), False),
    (lambda d: d["GermlineSet"][0].update(germline_set_ref="OGRDB:IGKappa_VJ.4"), False),
    (lambda d: d["GermlineSet"][0].update(release_version=4.5), False),
    (lambda d: d["GermlineSet"][0].update(release_version=5.0), False),
    (lambda d: d["GermlineSet"][0].update(germline_set_name="IGKappa_VJ"), False),
    (lambda d: d["GermlineSet"][0].update(germline_set_name="IGLambda_VJ_extended"), False),
], ids=["ref-null", "release-null", "name-null", "ref-absent", "release-int", "ref-release", "ref-set",
        "release-fraction", "release-other", "name-other", "name-extended-in-airr"])
def test_stated_set_identifiers_must_match_the_lock(edit, ok):
    blobs = {"airr": edited(IGL, "airr", edit)}
    if ok:
        assert identify(IGL, blobs) == expected_identity(IGL)
    else:
        error("identity_mismatch", identify, IGL, blobs)


def test_the_extended_export_names_its_set_with_the_extended_suffix():
    blobs = {"airr_ex": edited(IGL, "airr_ex", lambda d: d["GermlineSet"][0].update(germline_set_name="IGLambda_VJ"))}
    assert identify(IGL, blobs) == expected_identity(IGL)
    blobs = {"airr_ex": edited(IGL, "airr_ex", lambda d: d["GermlineSet"][0].update(germline_set_name="Other"))}
    error("identity_mismatch", identify, IGL, blobs)


def _set(d):
    return d["GermlineSet"][0]


@pytest.mark.parametrize("data,kind,path", [
    (b"\xff\xfe{}", "malformed", "/assets/airr"),
    (b"{\"Info\": {}", "malformed", "/assets/airr"),
    (b"[]", "malformed", "/assets/airr"),
    (b"{\"GermlineSet\": [], \"GermlineSet\": []}", "malformed", "/assets/airr"),
    (b"{\"GermlineSet\": [{\"locus\": NaN}]}", "malformed", "/assets/airr"),
    (b"[" * 65 + b"1" + b"]" * 65, "resource_limit", "/assets/airr"),
    (encode({"Info": {}}), "malformed", "/assets/airr/GermlineSet"),
    (encode({"GermlineSet": {}}), "malformed", "/assets/airr/GermlineSet"),
    (encode({"GermlineSet": []}), "malformed", "/assets/airr/GermlineSet"),
    (encode({"GermlineSet": [{}, {}]}), "malformed", "/assets/airr/GermlineSet"),
    (encode({"GermlineSet": ["set"]}), "malformed", "/assets/airr/GermlineSet/0"),
], ids=["not-utf8", "truncated", "root-array", "duplicate-key", "nan", "too-deep", "no-set", "set-object",
        "no-sets", "two-sets", "set-not-object"])
def test_a_malformed_global_envelope_is_a_typed_error(data, kind, path):
    found = error(kind, identify, IGL, {"airr": data})
    assert found.path == path


@pytest.mark.parametrize("edit,path", [
    (lambda d: _set(d).pop("species"), "/species"),
    (lambda d: _set(d).update(species=None), "/species"),
    (lambda d: _set(d).update(species={"id": "9606", "label": "Homo sapiens"}), "/species/id"),
    (lambda d: _set(d).update(species={"id": "NCBITAXON:0", "label": "Homo sapiens"}), "/species/id"),
    (lambda d: _set(d).update(species={"id": "NCBITAXON:9606", "label": ""}), "/species/label"),
    (lambda d: _set(d).update(species={"id": "NCBITAXON:9606"}), "/species/label"),
    (lambda d: _set(d).pop("locus"), "/locus"),
    (lambda d: _set(d).update(locus="IGX"), "/locus"),
    (lambda d: _set(d).update(locus=None), "/locus"),
    (lambda d: _set(d).pop("allele_descriptions"), "/allele_descriptions"),
    (lambda d: _set(d).update(allele_descriptions={}), "/allele_descriptions"),
    (lambda d: _set(d).update(release_version="4"), "/release_version"),
    (lambda d: _set(d).update(release_version=True), "/release_version"),
    (lambda d: _set(d).update(germline_set_ref=4), "/germline_set_ref"),
    (lambda d: _set(d).update(germline_set_name=["IGLambda_VJ"]), "/germline_set_name"),
    (lambda d: _set(d).update(germline_set_id=7), "/germline_set_id"),
], ids=["no-species", "null-species", "taxon-id-form", "taxon-id-zero", "blank-name", "no-name", "no-locus",
        "unknown-locus", "null-locus", "no-records", "records-object", "release-text", "release-bool",
        "ref-number", "name-list", "set-id-number"])
def test_a_malformed_set_identity_is_a_typed_error(edit, path):
    found = error("malformed", identify, IGL, {"airr": edited(IGL, "airr", edit)})
    assert found.path == "/assets/airr/GermlineSet/0" + path


def test_identity_needs_a_record_export():
    found = error("missing_asset", identify, IGL, {"gapped": raw(IGL, "gapped")})
    assert found.path == "/assets"


def test_identity_reader_rejects_arguments_of_the_wrong_type():
    with pytest.raises(TypeError):
        identify(IGL, {"airr": raw(IGL, "airr").decode("utf-8")})
