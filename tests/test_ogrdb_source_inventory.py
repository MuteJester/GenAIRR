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


# ── records, candidates and observations ─────────────────────────────

# Record fields the inventory itself normalizes; every other known field is kept as is.
INTERPRETED_FIELDS = {"label", "sequence_type", "coding_sequence", "allele_description_id",
                      "allele_description_ref", "locus", "species", "aliases", "paralogs",
                      "unrearranged_support", "rearranged_support"}


@pytest.mark.parametrize("fmt", EXPORTS)
@pytest.mark.parametrize("source", SETS, ids=_key)
def test_every_record_occurrence_is_a_source_record(source, fmt):
    _, inventory = inventory_of(source)
    raws = raw_records(source, fmt)
    records = records_of(inventory, fmt)
    expected = EXPECT["sets"][_key(source)][fmt]
    assert len(records) == len(raws) == expected["records"]
    assert [r.locator.pointer for r in records] == [rec(i) for i in range(len(raws))]
    for record, item in zip(records, raws):
        assert (record.status, record.problem) == ("parsed", None)
        assert record.label == item["label"]
        assert record.provider_id == item["allele_description_id"]
        assert record.provider_ref == item["allele_description_ref"]
        assert record.segment == item["sequence_type"]
        assert (record.taxonomy_id, record.locus) == (9606, _entry(source)["locus"])
    counts = Counter(r.segment for r in records)
    assert {s: counts.get(s, 0) for s in "VDJ"} == {s: expected[s] for s in "VDJ"}


@pytest.mark.parametrize("fmt", EXPORTS)
@pytest.mark.parametrize("source", SETS, ids=_key)
def test_each_supplied_coding_sequence_is_a_candidate_of_its_own_record(source, fmt):
    _, inventory = inventory_of(source)
    raws = raw_records(source, fmt)
    candidates = [c for c in inventory.candidates if c.locator.asset == fmt]
    assert len(candidates) == len(raws)
    for i, (candidate, item) in enumerate(zip(candidates, raws)):
        assert candidate.record == m.Locator(fmt, rec(i))
        assert candidate.locator == m.Locator(fmt, rec(i) + "/coding_sequence")
        assert candidate.sequence == item["coding_sequence"]
        assert candidate.sha256 == _sha(item["coding_sequence"].encode("ascii"))
        assert inventory.candidates_of(candidate.record) == (candidate,)


@pytest.mark.parametrize("fmt", EXPORTS)
@pytest.mark.parametrize("source", SETS, ids=_key)
def test_every_record_field_is_observed_with_its_original_value(source, fmt):
    _, inventory = inventory_of(source)
    for i, item in enumerate(raw_records(source, fmt)):
        for key, value in item.items():
            found = observed(inventory, fmt, f"{rec(i)}/{key}")
            assert (found.presence, found.value.get(), found.record) == ("present", value, m.Locator(fmt, rec(i)))
            assert found.interpretation == ("interpreted" if key in INTERPRETED_FIELDS else "preserved"), key
        for field in ("aliases", "paralogs"):
            for j, name in enumerate(item[field] or ()):
                found = observed(inventory, fmt, f"{rec(i)}/{field}/{j}")
                assert (found.value.get(), found.interpretation) == (name, "interpreted")
        for j, delineation in enumerate(item["v_gene_delineations"] or ()):
            for key, value in delineation.items():
                found = observed(inventory, fmt, f"{rec(i)}/v_gene_delineations/{j}/{key}")
                assert (found.value.get(), found.interpretation) == (value, "preserved")


def _record_index(source, fmt, label) -> int:
    (index,) = [i for i, item in enumerate(raw_records(source, fmt)) if item["label"] == label]
    return index


def test_coordinate_fields_are_kept_uninterpreted_at_their_locators():
    sets = {_key(s): s for s in SETS}
    for example in EXPECT["j_examples"]:
        source, fmt = sets[example["set"]], example["export"]
        _, inventory = inventory_of(source)
        base = rec(_record_index(source, fmt, example["label"]))
        for field in ("gene_start", "gene_end", "j_cdr3_end", "j_codon_frame"):
            found = observed(inventory, fmt, f"{base}/{field}")
            assert (found.value.get(), found.interpretation) == (example[field], "preserved")
    for example in EXPECT["v_examples"]:
        source, fmt = sets[example["set"]], example["export"]
        _, inventory = inventory_of(source)
        base = rec(_record_index(source, fmt, example["label"]))
        found = observed(inventory, fmt, f"{base}/v_gene_delineations/0/cdr3_start")
        assert (found.value.get(), found.interpretation) == (example["cdr3_start"], "preserved")
        assert observed(inventory, fmt, f"{base}/functional").value.get() is example["functional"]
    _, inventory = inventory_of(IGL)
    for example in EXPECT["igl4_airr_ex_coding_end_observations"]:
        base = rec(_record_index(IGL, "airr_ex", example["label"]))
        for field in ("gene_end", "v_rs_start"):
            found = observed(inventory, "airr_ex", f"{base}/{field}")
            assert (found.value.get(), found.interpretation) == (example[field], "preserved")


@pytest.mark.parametrize("source", SETS, ids=_key)
def test_set_level_fields_are_observed_without_a_record(source):
    _, inventory = inventory_of(source)
    for fmt in EXPORTS:
        document = payload(source, fmt)
        assert observed(inventory, fmt, "/Info").value.get() == document["Info"]
        for key, value in document["GermlineSet"][0].items():
            if key == "allele_descriptions":
                with pytest.raises(KeyError):
                    observed(inventory, fmt, f"/GermlineSet/0/{key}")
                continue
            found = observed(inventory, fmt, f"/GermlineSet/0/{key}")
            assert (found.record, found.presence, found.value.get()) == (None, "present", value)
    assert observed(inventory, "airr", "/GermlineSet/0/germline_set_ref").value.get() == _entry(source)["germline_set_ref"]


@pytest.mark.parametrize("fmt", EXPORTS)
@pytest.mark.parametrize("source", SETS, ids=_key)
def test_support_entries_are_embedded_unretrieved_references_of_their_own_record(source, fmt):
    _, inventory = inventory_of(source)
    support = [s for s in inventory.support if s.locator.asset == fmt]
    expected = []
    for i, item in enumerate(raw_records(source, fmt)):
        for kind in ("unrearranged", "rearranged"):
            for j, entry in enumerate(item[f"{kind}_support"]):
                pointer = f"{rec(i)}/{kind}_support/{j}"
                expected.append((pointer, rec(i), kind, entry["repository_name"], entry["repository_ref"]))
                assert observed(inventory, fmt, pointer).value.get() == entry
    found = [(s.locator.pointer, s.record.pointer, s.kind, s.repository, s.reference) for s in support]
    assert sorted(found) == sorted(expected)
    assert {s.status for s in support} <= {"embedded_unretrieved"}


@pytest.mark.parametrize("fmt", ["gapped", "gapped_ex"])
@pytest.mark.parametrize("source", SETS, ids=_key)
def test_fasta_records_are_header_and_sequence_observations(source, fmt):
    _, inventory = inventory_of(source)
    blocks = raw(source, fmt).decode("utf-8").split(">")[1:]
    expected = [{"header": b.split("\n", 1)[0], "sequence": b.split("\n", 1)[1].replace("\n", "")} for b in blocks]
    found = [o for o in inventory.observations if o.locator.asset == fmt]
    assert len(found) == EXPECT["sets"][_key(source)]["gapped_records"][fmt]
    assert [o.value.get() for o in found] == expected
    assert [o.locator.pointer for o in found] == [f"/{n}" for n in range(len(expected))]
    assert {(o.record, o.presence, o.interpretation) for o in found} == {(None, "present", "preserved")}


def test_fasta_line_terminators_are_removed_and_everything_else_is_kept():
    data = "﻿>IGLJ1*01 note é\r\nTTA.TG\r\n\r\nccg \n>IGLJ2*01\n\n".encode("utf-8")
    blobs = pinned(IGL, ("airr",))
    blobs["gapped"] = data
    _, inventory = build(IGL, blobs)
    found = [o.value.get() for o in inventory.observations if o.locator.asset == "gapped"]
    assert found == [{"header": "IGLJ1*01 note é", "sequence": "TTA.TGccg "},
                     {"header": "IGLJ2*01", "sequence": ""}]


@pytest.mark.parametrize("data", [b"TTAT\n>IGLJ1*01\nTTAT\n", b">IGLJ1*01\n\xff\n"], ids=["headless", "not-utf8"])
def test_a_malformed_fasta_export_is_a_typed_error(data):
    blobs = pinned(IGL, ("airr",))
    blobs["gapped"] = data
    lock = seal(IGL, blobs)
    assert error("malformed", ogrdb.build_inventory, lock, blobs).path == "/assets/gapped"


def test_the_inventory_is_bound_to_its_lock_and_round_trips():
    lock, inventory = inventory_of(IGL)
    assert inventory.source_lock_hash == m.lock_hash(lock)
    assert inventory.adapter == lock.adapter
    m.check_inventory(inventory, lock)
    assert m.SourceInventory.from_json(json.dumps(inventory.to_dict())) == inventory
    again = ogrdb.build_inventory(lock, pinned(IGL))
    assert m.inventory_hash(again) == m.inventory_hash(inventory)


def test_the_inventory_reads_exactly_the_locked_bytes():
    blobs = pinned(IGL)
    lock = seal(IGL, blobs)
    found = error("digest_mismatch", ogrdb.build_inventory, lock, dict(blobs, gapped=blobs["gapped"] + b"\n"))
    assert found.path == "/assets/gapped"
    found = error("digest_mismatch", ogrdb.build_inventory, lock, dict(blobs, gapped=blobs["gapped_ex"]))
    assert found.path == "/assets/gapped"
    found = error("missing_asset", ogrdb.build_inventory, lock, {k: v for k, v in blobs.items() if k != "airr_ex"})
    assert found.path == "/assets/airr_ex"
    found = error("identity_mismatch", ogrdb.build_inventory, lock, dict(blobs, ungapped=b">x\nACGT\n"))
    assert found.path == "/assets/ungapped"
    with pytest.raises(TypeError):
        ogrdb.build_inventory(lock.to_dict(), blobs)
    with pytest.raises(TypeError):
        ogrdb.build_inventory(lock, list(blobs.items()))


def test_the_lock_identity_must_be_the_payload_identity():
    blobs = pinned(IGL)
    lock = seal(IGL, blobs, identity=expected_identity(IGK))
    assert error("identity_mismatch", ogrdb.build_inventory, lock, blobs).path == "/identity"
    lock = seal(IGL, blobs, set_id="9606.IGKappa_VJ")
    assert error("identity_mismatch", ogrdb.build_inventory, lock, blobs).path.startswith("/assets/airr")


def test_inventories_of_locks_outside_the_validated_profile_scope_are_still_built():
    blobs = {"airr": edited(IGL, "airr", lambda d: d["GermlineSet"][0].update(curation="edited"))}
    lock, inventory = build(IGL, blobs)
    assert m.profile_scope_problem(lock) is not None
    assert len(inventory.records) == EXPECT["sets"][_key(IGL)]["airr"]["records"]


def test_the_identity_reader_and_the_inventory_share_the_ingested_bytes(tmp_path):
    files = {}
    for fmt in ALL:
        path = tmp_path / f"{fmt}.{SUFFIX[fmt]}"
        path.write_bytes(raw(IGL, fmt))
        files[fmt] = path
    result = ingest_files(set_id=IGL[0], release=IGL[1], profile="ogrdb_api_v2_human_ig_observed_v1",
                          files=files, cache=tmp_path / "cache", identify=ogrdb.identify)
    again = read_locked(result.lock, cache=tmp_path / "cache", identify=ogrdb.identify)
    inventory = ogrdb.build_inventory(again.lock, again.assets)
    assert inventory.source_lock_hash == m.lock_hash(result.lock)
    assert m.inventory_hash(inventory) != m.inventory_hash(inventory_of(IGL)[1])  # another lock, other ids


def test_relationship_only_names_are_mentioned_identities_without_sequence():
    _, inventory = inventory_of(IGL)
    j2 = rec(_record_index(IGL, "airr_ex", "IGLJ2*01"))
    (identity,) = inventory.mentioned_identities
    assert (identity.asset, identity.segment, identity.label) == ("airr_ex", "J", "IGLJ3*01")
    assert [(x.record.pointer, x.locator.pointer, x.relationship) for x in identity.mentions] == [
        (j2, j2 + "/paralogs/0", "paralog")]
    assert "IGLJ3*01" not in by_label(inventory, "airr_ex")
    assert [r.label for r in records_of(inventory, "airr")].count("IGLJ3*01") == 1


# ── malformed, unsupported and unknown content ──────────────────────


def with_records(source, fmt, edit) -> bytes:
    return edited(source, fmt, lambda d: edit(d["GermlineSet"][0]["allele_descriptions"]))


def _broken_records(records):
    j1 = copy.deepcopy(records[0])
    absent = copy.deepcopy(j1)
    absent.pop("label")
    records.extend([
        "not a record",
        dict(j1, label=None),
        absent,
        dict(j1, label="IGLJ1 *01"),
        dict(j1, label="IGLJX*01", sequence_type="C"),
        dict(j1, label="IGLDX*01", sequence_type="D"),
        dict(j1, label="IGLJX*02", sequence_type=7),
        dict(j1, label="IGLJX*03", coding_sequence="ACGN"),
        dict(j1, label="IGLJX*04", coding_sequence=None),
        dict(j1, label="IGLJX*05", allele_description_id=17, locus="IG\nL", species={"id": "taxon"}),
    ])


def test_malformed_and_unsupported_records_are_kept_with_diagnostics():
    lock, inventory = build(IGL, {"airr": with_records(IGL, "airr", _broken_records)})
    m.check_inventory(inventory, lock)
    n = len(raw_records(IGL, "airr"))
    records = records_of(inventory, "airr")
    assert len(records) == n + 10
    assert [(r.status, r.label, r.segment) for r in records[n:]] == [
        ("malformed", None, None),
        ("malformed", None, "J"),
        ("malformed", None, "J"),
        ("malformed", "IGLJ1 *01", "J"),
        ("unsupported", "IGLJX*01", None),
        ("unsupported", "IGLDX*01", "D"),
        ("malformed", "IGLJX*02", None),
        ("parsed", "IGLJX*03", "J"),
        ("parsed", "IGLJX*04", "J"),
        ("parsed", "IGLJX*05", "J"),
    ]
    for record in records[n:n + 7]:
        assert record.problem
        (finding,) = [f for f in inventory.findings if record.locator in f.subjects
                      and f.code == f"{record.status}_record"]
        assert finding.observations.get() == {"problem": record.problem}
    assert records[n].problem == "the record is a JSON string, not an object"
    assert records[n + 2].problem == "the record has no label"
    assert "D record cannot belong to a IGL (vj) set" in records[n + 5].problem
    without = {f.subjects[0].pointer: f.observations.get()["reason"]
               for f in findings(inventory, "record_without_sequence_candidate")}
    assert without == {rec(n + 7): "not a sequence of A, C, G and T", rec(n + 8): "null"}
    candidates = {c.record.pointer for c in inventory.candidates if c.locator.asset == "airr"}
    assert candidates == {rec(i) for i in range(n + 10)} - {rec(n), rec(n + 7), rec(n + 8)}
    assert observed(inventory, "airr", rec(n + 1) + "/label").interpretation == "invalid"
    absent = observed(inventory, "airr", rec(n + 2) + "/label")
    assert (absent.presence, absent.value.get(), absent.interpretation) == ("absent", None, "invalid")
    assert observed(inventory, "airr", rec(n + 7) + "/coding_sequence").value.get() == "ACGN"
    last = records[n + 9]
    assert (last.provider_id, last.locus, last.taxonomy_id) == (None, None, None)
    invalid = {f.observations.get()["field"] for f in findings(inventory, "invalid_field_value")
               if f.subjects == (last.locator,)}
    assert invalid == {"allele_description_id", "locus", "species"}
    with pytest.raises(KeyError):
        observed(inventory, "airr", rec(n) + "/label")


def test_relationship_and_support_values_of_the_wrong_form_are_kept_and_flagged():
    def edit(records):
        records[0]["aliases"] = "IGLJ1*01x"
        records[1]["paralogs"] = ["ok*01", "not ok", 5]
        records[2]["unrearranged_support"] = ["entry", {"repository_name": "", "repository_ref": "R1"}]
    lock, inventory = build(IGL, {"airr": with_records(IGL, "airr", edit)})
    assert observed(inventory, "airr", rec(0) + "/aliases").interpretation == "invalid"
    assert [observed(inventory, "airr", rec(1) + f"/paralogs/{j}").interpretation for j in range(3)] == [
        "interpreted", "invalid", "invalid"]
    (identity,) = [i for i in inventory.mentioned_identities if i.asset == "airr"]
    assert identity.label == "ok*01"
    support = [s for s in inventory.support if s.record == m.Locator("airr", rec(2))]
    assert [(s.locator.pointer, s.repository, s.reference) for s in support] == [
        (rec(2) + "/unrearranged_support/1", None, "R1")]
    assert observed(inventory, "airr", rec(2) + "/unrearranged_support/0").interpretation == "invalid"
    flagged = {(f.subjects[0].pointer, f.observations.get()["field"]) for f in findings(inventory, "invalid_field_value")}
    assert flagged == {(rec(0), "aliases"), (rec(1), "paralogs"), (rec(2), "unrearranged_support"),
                       (rec(2), "unrearranged_support/repository_name")}


def test_a_record_whose_identity_differs_from_its_set_stays_a_record():
    def edit(records):
        records[0]["locus"] = "IGK"
        records[1]["species"] = {"id": "NCBITAXON:10090", "label": "Mus musculus"}
    _, inventory = build(IGL, {"airr": with_records(IGL, "airr", edit)})
    found = {f.subjects[0].pointer: f.observations.get() for f in findings(inventory, "record_identity_differs_from_set")}
    assert found == {rec(0): {"locus": {"record": "IGK", "set": "IGL"}},
                     rec(1): {"taxonomy_id": {"record": 10090, "set": 9606}}}
    assert [r.status for r in records_of(inventory, "airr")[:2]] == ["parsed", "parsed"]
    assert records_of(inventory, "airr")[1].taxonomy_id == 10090


def test_unknown_provider_fields_are_kept_as_unsupported_observations_and_recorded():
    def edit(d):
        d["Extra"] = {"note": 1}
        d["GermlineSet"][0]["new_set_field"] = [1, 2]
        d["GermlineSet"][0]["species"]["authority"] = "x"
        records = d["GermlineSet"][0]["allele_descriptions"]
        for i in range(3):
            records[i]["new_record_field"] = {"v": i}
        records[3]["species"]["rank"] = "species"
        v = next(r for r in records if r["v_gene_delineations"])
        v["v_gene_delineations"][0]["new_delineation_field"] = 5
        supported = next(r for r in records if r["unrearranged_support"])
        supported["unrearranged_support"][0]["new_support_field"] = "s"
    lock, inventory = build(IGL, {"airr": edited(IGL, "airr", edit)})
    m.check_inventory(inventory, lock)
    extra = observed(inventory, "airr", "/Extra")
    assert (extra.record, extra.value.get(), extra.interpretation) == (None, {"note": 1}, "unsupported")
    assert observed(inventory, "airr", "/GermlineSet/0/new_set_field").interpretation == "unsupported"
    for i in range(3):
        found = observed(inventory, "airr", rec(i) + "/new_record_field")
        assert (found.record, found.value.get(), found.interpretation) == (m.Locator("airr", rec(i)), {"v": i},
                                                                           "unsupported")
    recorded = {f.observations.get()["field"]: f for f in findings(inventory, "unknown_provider_field")}
    assert set(recorded) == {"/Extra", "/GermlineSet/0/new_set_field", "/GermlineSet/0/species/authority",
                             RECORDS + "/*/new_record_field", RECORDS + "/*/species/rank",
                             RECORDS + "/*/v_gene_delineations/*/new_delineation_field",
                             RECORDS + "/*/unrearranged_support/*/new_support_field"}
    many = recorded[RECORDS + "/*/new_record_field"]
    assert many.observations.get()["occurrences"] == 3
    assert [x.pointer for x in many.locators] == [rec(i) + "/new_record_field" for i in range(3)]
    assert {(f.severity, f.category) for f in recorded.values()} == {("info", "informational")}
    assert [r.status for r in records_of(inventory, "airr")] == ["parsed"] * len(raw_records(IGL, "airr"))


def test_unknown_fields_no_pointer_can_name_are_recorded_at_their_record():
    long_key, control_key = "k" * 2000, "bad\u0001key"
    def edit(records):
        records[0][long_key] = 1
        records[0][control_key] = 2
    _, inventory = build(IGL, {"airr": with_records(IGL, "airr", edit)})
    recorded = {f.observations.get()["field"]: f for f in findings(inventory, "unknown_provider_field")}
    for key in (long_key, control_key):
        finding = recorded[RECORDS + "/*/" + key]
        assert finding.locators == (m.Locator("airr", rec(0)),)


def test_the_aggregate_source_record_budget_spans_both_record_exports():
    from GenAIRR._reference_sources import _fields
    assert ogrdb.MAX_SOURCE_RECORDS == _fields.MAX_SOURCE_RECORDS == 100_000

    def sized(fmt, count):
        return edited(IGL, fmt, lambda d: d["GermlineSet"][0].update(allele_descriptions=[{}] * count))
    blobs = {"airr": sized("airr", 50_000), "airr_ex": sized("airr_ex", 50_001)}
    found = error("resource_limit", ogrdb.build_inventory, seal(IGL, blobs), blobs)
    assert found.path == "/assets" and "100001" in found.message


def test_the_record_budget_boundary_is_inclusive(monkeypatch):
    monkeypatch.setattr(ogrdb, "MAX_SOURCE_RECORDS", 4)
    blobs = {fmt: edited(IGL, fmt, lambda d: d["GermlineSet"][0].update(
        allele_descriptions=d["GermlineSet"][0]["allele_descriptions"][:2])) for fmt in EXPORTS}
    _, inventory = build(IGL, blobs)
    assert len(inventory.records) == 4
    blobs["airr"] = edited(IGL, "airr", lambda d: d["GermlineSet"][0].update(
        allele_descriptions=d["GermlineSet"][0]["allele_descriptions"][:3]))
    error("resource_limit", ogrdb.build_inventory, seal(IGL, blobs), blobs)


def test_more_observations_than_an_inventory_holds_is_a_resource_limit(monkeypatch):
    monkeypatch.setattr(ogrdb, "MAX_OBSERVATIONS", 100)
    blobs = pinned(IGL, ("airr",))
    assert error("resource_limit", ogrdb.build_inventory, seal(IGL, blobs), blobs).path == "/assets"


def _nested(depth):
    value = 1
    for _ in range(depth):
        value = [value]
    return value


def test_provider_values_are_kept_to_the_parser_depth_limit():
    # A record field sits five containers deep; 59 more reach the 64-level parser limit.
    blobs = {"airr": with_records(IGL, "airr", lambda records: records[0].update(curation=_nested(59)))}
    _, inventory = build(IGL, blobs)
    assert observed(inventory, "airr", rec(0) + "/curation").value.get() == _nested(59)
    blobs = {"airr": with_records(IGL, "airr", lambda records: records[0].update(curation=_nested(60)))}
    found = error("resource_limit", ogrdb.build_inventory, seal(IGL, blobs), blobs)
    assert found.path == "/assets/airr"


def test_an_envelope_value_nested_beyond_a_preserved_value_is_a_resource_limit():
    blobs = {"airr": edited(IGL, "airr", lambda d: d.update(Info=_nested(62)))}
    found = error("resource_limit", ogrdb.build_inventory, seal(IGL, blobs), blobs)
    assert found.path.startswith("/assets/airr/Info")


@pytest.mark.parametrize("text", [
    '{"GermlineSet": [{"allele_descriptions": [{"label": "a", "label": "b"}]}]}',
    '{"GermlineSet": [{"allele_descriptions": [{"gene_start": Infinity}]}]}',
], ids=["duplicate-key-in-a-record", "non-finite-number"])
def test_provider_json_is_parsed_strictly(text):
    blobs = {"airr": text.encode("utf-8")}
    lock = seal(IGL, blobs)
    assert error("malformed", ogrdb.build_inventory, lock, blobs).path == "/assets/airr"


# ── accounting: every count reconciles to the raw payloads ──────────


def account_of(source, fmt, formats=ALL):
    lock, inventory = inventory_of(source, formats)
    return ogrdb.account(lock, inventory).export(fmt)


def _items(item, field):
    value = item.get(field)
    return value if isinstance(value, list) else []


@pytest.mark.parametrize("fmt", EXPORTS)
@pytest.mark.parametrize("source", SETS, ids=_key)
def test_every_conservation_equation_reconciles_to_the_raw_payload(source, fmt):
    lock, inventory = inventory_of(source)
    account = ogrdb.account(lock, inventory).export(fmt)
    raws = raw_records(source, fmt)
    expected = EXPECT["sets"][_key(source)][fmt]
    # Record occurrences: every raw record, by status and by segment.
    assert account.record_occurrences == len(raws) == expected["records"]
    assert account.record_occurrences == (account.parsed_records + account.malformed_records
                                          + account.unsupported_records)
    assert account.segments["record_occurrences"] == {"V": expected["V"], "D": expected["D"], "J": expected["J"],
                                                     "unassigned": 0}
    # Candidates: at most one per record.
    assert account.record_occurrences == account.sequence_candidates + account.records_without_candidate
    assert account.sequence_candidates == sum(1 for item in raws if item["coding_sequence"])
    # Parsed records: first supply of each identity's sequence, repeats, or no candidate.
    assert account.parsed_records == (account.identity_sequences + account.duplicate_occurrences
                                      + account.parsed_without_candidate)
    assert account.identity_sequences == (account.record_identities - account.identities_without_candidate
                                          + account.alternative_sequences)
    assert account.record_identities == len({(item["sequence_type"], item["label"]) for item in raws})
    # Named identities: own records plus relationship-only names; aliases add none.
    paralog_only = expected["paralog_only_identities"]
    assert account.mentioned_only_identities == len(paralog_only)
    assert account.named_identities == account.record_identities + account.mentioned_only_identities
    assert account.alias_names == sum(len(_items(item, "aliases")) for item in raws)
    assert account.paralog_names == sum(len(_items(item, "paralogs")) for item in raws)
    assert account.paralog_names == (account.paralog_names_of_records + account.paralog_names_of_mentioned_only
                                     + account.paralog_names_unresolved)
    mentions = sum(len(i.mentions) for i in inventory.mentioned_identities if i.asset == fmt)
    assert account.paralog_names_of_mentioned_only == mentions
    labels = {item["label"] for item in raws}
    assert account.paralog_names_of_records == sum(name in labels for item in raws for name in _items(item, "paralogs"))
    assert account.paralog_names_unresolved == 0
    # Support entries and provider refs.
    assert account.support_references == sum(len(item["unrearranged_support"]) + len(item["rearranged_support"])
                                              for item in raws)
    assert account.provider_refs_with_ext_suffix == expected["refs_ending_ext"]
    assert account.provider_refs_with_ext_suffix == sum(item["allele_description_ref"].endswith("_ext")
                                                        for item in raws)


@pytest.mark.parametrize("source", SETS, ids=_key)
def test_identity_counts_of_both_exports_close_over_relationship_only_names(source):
    airr, airr_ex = account_of(source, "airr"), account_of(source, "airr_ex")
    # The extended export names, as paralogs only, exactly the identities it has no record for.
    assert airr_ex.named_identities == airr.named_identities == airr.record_occurrences
    assert airr.mentioned_only_identities == 0


def test_igl_release4_explicit_records_and_mentioned_identities():
    airr, airr_ex = account_of(IGL, "airr"), account_of(IGL, "airr_ex")
    assert (airr_ex.record_occurrences, airr_ex.sequence_candidates, airr_ex.named_identities,
            airr_ex.mentioned_only_identities) == (97, 97, 98, 1)
    assert (airr.record_occurrences, airr.sequence_candidates, airr.named_identities,
            airr.mentioned_only_identities) == (98, 98, 98, 0)
    # Executable alleles are chosen later, by an explicit selection; nothing here counts them.
    assert not hasattr(airr, "selected_alleles")


@pytest.mark.parametrize("fmt", EXPORTS)
@pytest.mark.parametrize("source", SETS, ids=_key)
def test_paralog_only_names_are_the_expected_mentioned_identities(source, fmt):
    _, inventory = inventory_of(source)
    found = sorted(i.label for i in inventory.mentioned_identities if i.asset == fmt)
    assert found == sorted(EXPECT["sets"][_key(source)][fmt]["paralog_only_identities"])
    raws = raw_records(source, fmt)
    labels = {item["label"] for item in raws}
    for identity in (i for i in inventory.mentioned_identities if i.asset == fmt):
        assert identity.label not in labels
        assert inventory.candidates_of(m.Locator(fmt, "/" + identity.label)) == ()
        for mention in identity.mentions:
            item = raws[int(mention.record.pointer.rsplit("/", 1)[1])]
            assert identity.label in item["paralogs"]
            assert (mention.relationship, identity.segment) == ("paralog", item["sequence_type"])


@pytest.mark.parametrize("source", SETS, ids=_key)
def test_identical_dna_under_distinct_labels_stays_distinct_identities(source):
    lock, inventory = inventory_of(source)
    expected = sorted(tuple(sorted(group)) for group in EXPECT["sets"][_key(source)]["airr_identical_coding_groups"])
    assert ogrdb.account(lock, inventory).export("airr").shared_sequences == tuple(expected)
    assert ogrdb.account(lock, inventory).export("airr_ex").shared_sequences == ()
    found = findings(inventory, "identical_sequence_distinct_identities")
    assert sorted(tuple(i["label"] for i in f.observations.get()["identities"]) for f in found) == expected
    labels = {r.locator: r.label for r in inventory.records}
    for finding in found:
        records = [s for s in finding.subjects if s in labels]
        candidates = [s for s in finding.subjects if s not in labels]
        assert sorted(labels[r] for r in records) == sorted(i["label"] for i in finding.observations.get()["identities"])
        assert sorted(c.pointer for c in candidates) == sorted(r.pointer + "/coding_sequence" for r in records)
        assert len({inventory.candidate(c).sha256 for c in candidates}) == 1
        assert (finding.severity, finding.category) == ("info", "informational")


def test_j2_and_j3_keep_their_own_records_ids_and_support():
    lock, inventory = inventory_of(IGL)
    j2, j3 = by_label(inventory, "airr")["IGLJ2*01"], by_label(inventory, "airr")["IGLJ3*01"]
    (j2,), (j3,) = j2, j3
    (j2_candidate,), (j3_candidate,) = inventory.candidates_of(j2.locator), inventory.candidates_of(j3.locator)
    assert j2_candidate.sha256 == j3_candidate.sha256 and j2_candidate.locator != j3_candidate.locator
    assert (j2.provider_id, j3.provider_id) == ("OGRDB:A02870", "OGRDB:A02871")
    own = {}
    for support in inventory.support:
        own.setdefault(support.record, []).append(support.reference)
    assert own[j2.locator] == ["D87023", "Gibson:HG01106", "GRCh38:NC_000022.11"]
    assert own[j3.locator] == ["Gibson:HG02572"]
    # In airr_ex, IGLJ3*01 is only a paralog name: no record, candidate or support of its own.
    (ex_j2,) = by_label(inventory, "airr_ex")["IGLJ2*01"]
    assert [s.reference for s in inventory.support if s.record == ex_j2.locator] == own[j2.locator]
    j2_entries = [inventory.observation(s.locator).value for s in inventory.support if s.record == j2.locator]
    holders = {s.record for s in inventory.support if inventory.observation(s.locator).value in j2_entries}
    assert holders == {j2.locator, ex_j2.locator}
    (j3_mentioned,) = [i for i in inventory.mentioned_identities if i.label == "IGLJ3*01"]
    assert j3_mentioned.asset == "airr_ex"
    # No J2 field can be evidence about J3, in either export.
    for target in (j3.locator, j3_candidate.locator):
        claim = m.Claim(source_lock_hash=inventory.source_lock_hash, subject=target,
                        predicate="functionality" if target == j3.locator else "boundary_3",
                        value=None, coordinate_system=None, support="none", derivation="stated",
                        status="unresolved", evidence=(m.Locator("airr", j3.locator.pointer + "/functional"),),
                        mappings=())
        for borrowed in (j2.locator.pointer + "/unrearranged_support/0", j2.locator.pointer + "/functional"):
            for asset in EXPORTS:
                assert inventory.evidence_binding_problem(claim, m.Locator(asset, borrowed)) is not None
        assert inventory.evidence_binding_problem(claim, m.Locator("airr", j3.locator.pointer + "/functional")) is None


def test_declared_aliases_are_names_of_their_own_record_not_identities():
    _, inventory = inventory_of(IGL)
    for fmt in EXPORTS:
        raws = raw_records(IGL, fmt)
        aliases = [(i, j, name) for i, item in enumerate(raws) for j, name in enumerate(_items(item, "aliases"))]
        assert len(aliases) == 6
        mentioned = {x.label for x in inventory.mentioned_identities}
        for i, j, name in aliases:
            assert observed(inventory, fmt, f"{rec(i)}/aliases/{j}").value.get() == name
            assert name not in mentioned
        assert account_of(IGL, fmt).alias_names == 6
    assert findings(inventory, "alias_names_explicit_record") == []


def test_an_alias_naming_another_record_is_flagged_and_merges_nothing():
    def edit(records):
        records[0]["aliases"] = [records[1]["label"]]
    lock, inventory = build(IGL, {"airr": with_records(IGL, "airr", edit)})
    (finding,) = findings(inventory, "alias_names_explicit_record")
    labels = raw_records(IGL, "airr")
    assert finding.observations.get() == {"label": labels[0]["label"], "alias": labels[1]["label"]}
    assert finding.subjects == (m.Locator("airr", rec(0)), m.Locator("airr", rec(1)))
    assert ogrdb.account(lock, inventory).export("airr").record_identities == len(labels)


def test_duplicate_occurrences_and_same_label_alternatives_are_counted_apart():
    base = raw_records(IGL, "airr")
    first = base[0]
    def edit(records):
        records.append(copy.deepcopy(first))
        records.append(copy.deepcopy(first))
        records.append(dict(copy.deepcopy(first), coding_sequence=first["coding_sequence"] + "A"))
    lock, inventory = build(IGL, {"airr": with_records(IGL, "airr", edit)})
    account = ogrdb.account(lock, inventory).export("airr")
    n = len(base)
    assert (account.record_occurrences, account.sequence_candidates, account.record_identities) == (n + 3, n + 3, n)
    assert (account.identity_sequences, account.duplicate_occurrences) == (n + 1, 2)
    assert (account.alternative_identities, account.alternative_sequences) == (1, 1)
    assert account.parsed_records == account.identity_sequences + account.duplicate_occurrences
    (duplicate,) = findings(inventory, "duplicate_record_occurrence")
    assert duplicate.observations.get()["occurrences"] == 3
    assert [s.pointer for s in duplicate.subjects] == [rec(0), rec(0) + "/coding_sequence", rec(n),
                                                       rec(n) + "/coding_sequence", rec(n + 1),
                                                       rec(n + 1) + "/coding_sequence"]
    (alternatives,) = findings(inventory, "same_label_sequence_alternatives")
    observations = alternatives.observations.get()
    assert (observations["label"], observations["segment"]) == (first["label"], first["sequence_type"])
    assert sorted((a["length"], a["records"]) for a in observations["alternatives"]) == [
        (len(first["coding_sequence"]), 3), (len(first["coding_sequence"]) + 1, 1)]
    assert len(alternatives.subjects) == 8
    assert (alternatives.severity, alternatives.category) == ("warning", "contradictory_evidence")
    (shared,) = findings(inventory, "duplicate_provider_id")
    assert shared.observations.get() == {"provider_id": first["allele_description_id"], "records": 4}


def test_duplicate_or_missing_provider_ids_never_merge_records():
    def edit(records):
        records[0]["allele_description_id"] = None
        records[1].pop("allele_description_id")
        records[3]["allele_description_id"] = records[2]["allele_description_id"]
    lock, inventory = build(IGL, {"airr": with_records(IGL, "airr", edit)})
    records = records_of(inventory, "airr")
    assert len(records) == len(raw_records(IGL, "airr"))
    assert [r.provider_id for r in records[:2]] == [None, None]
    assert len({r.id for r in records}) == len(records)
    (finding,) = findings(inventory, "duplicate_provider_id")
    assert [s.pointer for s in finding.subjects] == [rec(2), rec(3)]
    absent = observed(inventory, "airr", rec(1) + "/allele_description_id")
    assert (absent.presence, absent.interpretation) == ("absent", "interpreted")


@pytest.mark.parametrize("source", SETS, ids=_key)
def test_fasta_records_are_counted_per_export(source):
    lock, inventory = inventory_of(source)
    assert dict(ogrdb.account(lock, inventory).fasta) == EXPECT["sets"][_key(source)]["gapped_records"]


def test_accounts_belong_to_the_inventory_of_their_lock():
    lock, inventory = inventory_of(IGL)
    other, _ = inventory_of(IGK)
    error("digest_mismatch", ogrdb.account, other, inventory)
    with pytest.raises(TypeError):
        ogrdb.account(lock, inventory.to_dict())


def test_paralog_names_that_resolve_to_nothing_are_counted_as_unresolved():
    raws = raw_records(IGL, "airr")
    assert [raws[i]["sequence_type"] for i in range(4)] == ["J", "V", "V", "V"]

    def edit(records):
        records[0]["paralogs"] = ["bad name", "IGLX*01", "IGLZ*01"]   # invalid; mentioned only; named by J ...
        records[1]["paralogs"] = ["IGLZ*01"]                           # ... and by V: no single segment
        records[3]["label"] = "bad label"                              # a malformed record naming a paralog
        records[3]["paralogs"] = ["IGLY*01"]
    lock, inventory = build(IGL, {"airr": with_records(IGL, "airr", edit)})
    account = ogrdb.account(lock, inventory).export("airr")
    assert (account.paralog_names, account.paralog_names_of_records, account.paralog_names_of_mentioned_only,
            account.paralog_names_unresolved) == (7, 2, 1, 4)
    assert [(i.label, i.segment) for i in inventory.mentioned_identities] == [("IGLX*01", "J")]
    (conflict,) = findings(inventory, "mentioned_identity_segment_conflict")
    assert conflict.observations.get() == {"label": "IGLZ*01", "segments": ["J", "V"]}
    assert [s.pointer for s in conflict.subjects] == [rec(0), rec(1)]


# ── endpoint comparison: airr and airr_ex, label by label ───────────

PAIR_FINDINGS = ("export_record_sequence_differs", "export_support_differs", "export_provider_ref_differs",
                 "export_annotation_differs")


def _pairs(source):
    """Raw records of both exports by (segment, label), read with the standard library."""
    return {fmt: {(item["sequence_type"], item["label"]): item for item in raw_records(source, fmt)}
            for fmt in EXPORTS}


def coding_differences(inventory) -> dict:
    return {f.observations.get()["label"]: f for f in findings(inventory, "export_coding_sequence_differs")}


@pytest.mark.parametrize("source", SETS, ids=_key)
def test_coding_sequence_differences_are_the_expected_extensions(source):
    _, inventory = inventory_of(source)
    expected = EXPECT["sets"][_key(source)]["airr_ex_coding_suffix_beyond_airr"]
    found = coding_differences(inventory)
    assert sorted(found) == sorted(expected)
    pairs = _pairs(source)
    for label, finding in found.items():
        observations = finding.observations.get()
        assert (observations["relation"], observations["extra"]) == ("airr_ex_extends_3prime", expected[label])
        segment = observations["segment"]
        airr, airr_ex = pairs["airr"][(segment, label)], pairs["airr_ex"][(segment, label)]
        assert observations["airr"]["sequences"] == [{"sha256": _sha(airr["coding_sequence"].encode("ascii")),
                                                      "length": len(airr["coding_sequence"])}]
        assert observations["airr_ex"]["sequences"][0]["length"] == len(airr_ex["coding_sequence"])
        records = {r.locator.asset: r for r in inventory.records if r.label == label}
        assert set(finding.subjects) == {records["airr"].locator, records["airr_ex"].locator,
                                         m.Locator("airr", records["airr"].locator.pointer + "/coding_sequence"),
                                         m.Locator("airr_ex", records["airr_ex"].locator.pointer + "/coding_sequence")}
        assert finding.capabilities == ("original_coding_ends",)
        assert (finding.severity, finding.category) == ("warning", "contradictory_evidence")


def test_igl_release4_six_coding_differences_match_the_coding_end_table():
    _, inventory = inventory_of(IGL)
    found = coding_differences(inventory)
    table = {row["label"]: row["coding_extension"] for row in EXPECT["igl4_airr_ex_coding_end_observations"]}
    assert len(found) == 6
    assert {label: f.observations.get()["extra"] for label, f in found.items()} == {**table, "IGLV3-9*01": "CACCC"}


def test_igh_release10_extensions_match_the_candidate_selection():
    _, inventory = inventory_of(IGH)
    found = coding_differences(inventory)
    assert len(found) == 23 and all(f.observations.get()["segment"] == "V" for f in found.values())
    chosen = {c["label"]: c for c in SELECTION["record_choices"]}
    unchosen = {o["label"] for o in SELECTION["omitted"]
                if any(r["rule"] == "explicit_record_choice" for r in o["failed_rules"])}
    assert set(found) == set(chosen) | unchosen and not set(chosen) & unchosen
    for label, choice in chosen.items():
        observations = found[label].observations.get()
        assert observations["extra"] == choice["alternative"]["extension_3prime"]
        assert observations["airr"]["sequences"][0]["length"] == choice["chosen_coding_length"]
        assert observations["airr_ex"]["sequences"][0]["length"] == choice["alternative"]["coding_length"]


@pytest.mark.parametrize("source", SETS, ids=_key)
def test_pinned_pairs_differ_only_in_provider_ref_and_coding_sequence(source):
    pairs = _pairs(source)
    for key in pairs["airr"].keys() & pairs["airr_ex"].keys():
        a, b = pairs["airr"][key], pairs["airr_ex"][key]
        assert {k for k in set(a) | set(b) if a.get(k) != b.get(k)} <= {"allele_description_ref", "coding_sequence"}
    _, inventory = inventory_of(source)
    for code in PAIR_FINDINGS:
        assert findings(inventory, code) == [], code


@pytest.mark.parametrize("source", SETS, ids=_key)
def test_the_ext_suffix_is_recorded_once_for_the_comparison(source):
    _, inventory = inventory_of(source)
    expected = EXPECT["sets"][_key(source)]
    pairs = len(_pairs(source)["airr"].keys() & _pairs(source)["airr_ex"].keys())
    (finding,) = findings(inventory, "export_provider_ref_suffix")
    assert finding.observations.get() == {
        "suffix": "_ext", "pairs": pairs, "suffix_only": pairs, "identical": 0, "other": 0,
        "airr": {"records": expected["airr"]["records"], "refs_with_suffix": expected["airr"]["refs_ending_ext"]},
        "airr_ex": {"records": expected["airr_ex"]["records"],
                    "refs_with_suffix": expected["airr_ex"]["refs_ending_ext"]}}
    assert finding.locators == (m.Locator("airr", "/GermlineSet/0/germline_set_ref"),
                                m.Locator("airr_ex", "/GermlineSet/0/germline_set_ref"))


def test_igl_release4_ext_suffix_on_every_extended_ref_and_no_plain_ref():
    _, inventory = inventory_of(IGL)
    (finding,) = findings(inventory, "export_provider_ref_suffix")
    observations = finding.observations.get()
    assert (observations["airr_ex"], observations["airr"]) == ({"records": 97, "refs_with_suffix": 97},
                                                               {"records": 98, "refs_with_suffix": 0})


def absent_labels(inventory) -> dict:
    return {f.observations.get()["label"]: f for f in findings(inventory, "export_label_absent")}


@pytest.mark.parametrize("source", SETS, ids=_key)
def test_labels_present_in_one_export_only_are_named_as_paralogs_in_the_other(source):
    _, inventory = inventory_of(source)
    expected = EXPECT["sets"][_key(source)]
    found = absent_labels(inventory)
    assert sorted(found) == sorted(expected["labels_only_in_airr"])
    assert expected["labels_only_in_airr_ex"] == []
    pairs = _pairs(source)
    for label, finding in found.items():
        observations = finding.observations.get()
        assert (observations["present_in"], observations["absent_from"]) == ("airr", "airr_ex")
        (named,) = observations["named_in_absent_export"]
        assert named["relationship"] == "paralog"
        counterpart = pairs["airr_ex"][(observations["segment"], named["by"])]
        assert label in counterpart["paralogs"]
        # Only the identity's own record and candidate are subjects; never its counterpart's.
        (own,) = [r for r in inventory.records if r.locator.asset == "airr" and r.label == label]
        assert set(finding.subjects) == {own.locator, m.Locator("airr", own.locator.pointer + "/coding_sequence")}
        (locator,) = finding.locators
        assert observed(inventory, "airr_ex", locator.pointer).value.get() == label


def test_igh_release10_airr_only_identities_and_their_counterparts():
    _, inventory = inventory_of(IGH)
    found = absent_labels(inventory)
    expected = {o["label"]: o["identical_coding_sequence_to"] for o in SELECTION["selected_record_observations"]
                if o["observation"] == "identical_coding_sequence_identity_without_airr_ex_record"}
    assert len(found) == len(expected) == 9
    assert {label: f.observations.get()["named_in_absent_export"][0]["by"] for label, f in found.items()} == expected


def test_igl_release4_j3_is_absent_from_the_extended_export_and_named_there_by_j2():
    _, inventory = inventory_of(IGL)
    finding = absent_labels(inventory)["IGLJ3*01"]
    assert finding.observations.get() == {
        "label": "IGLJ3*01", "segment": "J", "present_in": "airr", "absent_from": "airr_ex",
        "named_in_absent_export": [{"relationship": "paralog", "by": "IGLJ2*01"}]}
    (j3,) = by_label(inventory, "airr")["IGLJ3*01"]
    assert set(finding.subjects) == {j3.locator, m.Locator("airr", j3.locator.pointer + "/coding_sequence")}


def _label_index(records, label):
    (index,) = [i for i, item in enumerate(records) if item["label"] == label]
    return index


def test_each_kind_of_per_label_difference_is_recorded_without_preferring_an_export():
    airr_records, ex_records = raw_records(IGL, "airr"), raw_records(IGL, "airr_ex")
    shared = [item["label"] for item in ex_records if item["label"] != "IGLJ7*02" and item["coding_sequence"] in
              {a["coding_sequence"] for a in airr_records if a["label"] == item["label"]}]
    one, two, three, four, five = shared[:5]

    def edit_ex(records):
        records[_label_index(records, one)]["functional"] = not records[_label_index(records, one)]["functional"]
        records[_label_index(records, one)]["new_field"] = 1
        records[_label_index(records, two)]["sequence"] += "A"
        records[_label_index(records, three)]["unrearranged_support"] = []
        records[_label_index(records, four)]["allele_description_ref"] = "OGRDB:other"

    def edit_airr(records):
        records[_label_index(records, five)]["coding_sequence"] += "GG"
        del records[_label_index(records, "IGLJ7*02")]

    blobs = {"airr": with_records(IGL, "airr", edit_airr), "airr_ex": with_records(IGL, "airr_ex", edit_ex)}
    _, inventory = build(IGL, blobs)
    annotation = {f.observations.get()["label"]: f.observations.get()
                  for f in findings(inventory, "export_annotation_differs")}
    assert annotation == {one: {"label": one, "segment": "V" if one.startswith("IGLV") else "J",
                                "fields": ["functional", "new_field"]}}
    (sequence,) = findings(inventory, "export_record_sequence_differs")
    assert sequence.observations.get()["label"] == two
    (support,) = findings(inventory, "export_support_differs")
    observations = support.observations.get()
    assert observations["label"] == three
    assert observations["unrearranged_support"]["airr_ex"] == 0
    assert observations["unrearranged_support"]["equal"] is False and observations["rearranged_support"]["equal"]
    (ref,) = findings(inventory, "export_provider_ref_differs")
    assert ref.observations.get()["label"] == four and ref.observations.get()["airr_ex"] == "OGRDB:other"
    longer = coding_differences(inventory)[five].observations.get()
    assert (longer["relation"], longer["extra"]) == ("airr_extends_3prime", "GG")
    absent = absent_labels(inventory)["IGLJ7*02"].observations.get()
    assert (absent["present_in"], absent["absent_from"], absent["named_in_absent_export"]) == ("airr_ex", "airr", [])
    (suffix,) = findings(inventory, "export_provider_ref_suffix")
    assert suffix.observations.get()["other"] == 1


@pytest.mark.parametrize("source", SETS, ids=_key)
def test_a_complete_comparison_lists_its_findings(source):
    lock, inventory = inventory_of(source)
    comparison = ogrdb.comparison(lock, inventory)
    assert (comparison.complete, comparison.exports, comparison.missing) == (True, EXPORTS, ())
    codes = set(ogrdb.COMPARISON_CODES)
    assert comparison.findings == tuple(sorted(f.id for f in inventory.findings if f.code in codes))
    comparison.check(lock, inventory)
    summary = ogrdb.account(lock, inventory).comparison
    airr, airr_ex = (ogrdb.account(lock, inventory).export(f) for f in EXPORTS)
    assert airr.record_identities == summary.identities_in_both + len(summary.only_in["airr"])
    assert airr_ex.record_identities == summary.identities_in_both + len(summary.only_in["airr_ex"])
    assert sorted(summary.only_in["airr"]) == sorted(EXPECT["sets"][_key(source)]["labels_only_in_airr"])


@pytest.mark.parametrize("present", EXPORTS)
def test_single_export_inspection_succeeds_and_marks_the_comparison_incomplete(present):
    gapped = {"airr": "gapped", "airr_ex": "gapped_ex"}[present]
    lock, inventory = inventory_of(IGL, (present, gapped))
    comparison = ogrdb.comparison(lock, inventory)
    missing = "airr_ex" if present == "airr" else "airr"
    assert (comparison.complete, comparison.exports, comparison.missing) == (False, (present,), (missing,))
    comparison.check(lock, inventory)
    (incomplete,) = findings(inventory, "endpoint_comparison_incomplete")
    assert incomplete.observations.get() == {"exports": [present],
                                             "missing": [{"format": missing, "reason": "not_supplied"}]}
    assert incomplete.locators == (m.Locator(present, "/GermlineSet/0/germline_set_ref"),)
    assert (incomplete.severity, incomplete.category) == ("warning", "unavailable_evidence")
    assert incomplete.id in comparison.findings
    for code in ("export_coding_sequence_differs", "export_label_absent", "export_provider_ref_suffix"):
        assert findings(inventory, code) == []
    without = findings(inventory, "mentioned_identity_without_record")
    if present == "airr_ex":
        (finding,) = without
        assert finding.observations.get() == {"asset": "airr_ex", "label": "IGLJ3*01", "segment": "J"}
        (j2,) = by_label(inventory, "airr_ex")["IGLJ2*01"]
        assert finding.subjects == (j2.locator,)
    else:
        assert without == []
    summary = ogrdb.account(lock, inventory).comparison
    assert (summary.complete, summary.identities_in_both, summary.only_in) == (False, 0, {})


def test_comparison_needs_the_inventory_of_its_lock():
    lock, inventory = inventory_of(IGL)
    error("digest_mismatch", ogrdb.comparison, inventory_of(IGK)[0], inventory)


def test_an_exceeded_inventory_bound_is_located_in_the_inventory(monkeypatch):
    def over(**fields):
        raise ReferenceSourceError("resource_limit", "/findings", "65537 items exceed the limit of 65536")
    monkeypatch.setattr(ogrdb.SourceInventory, "assemble", staticmethod(over))
    blobs = pinned(IGL, ("airr",))
    found = error("resource_limit", ogrdb.build_inventory, seal(IGL, blobs), blobs)
    assert found.path == "/inventory/findings"


def test_findings_naming_more_than_a_finding_holds_are_split_into_parts(monkeypatch):
    monkeypatch.setattr(ogrdb, "MAX_REFS", 3)
    def edit(records):
        for i in range(7):
            records[i]["new_record_field"] = i
    _, inventory = build(IGL, {"airr": with_records(IGL, "airr", edit)})
    parts = [f for f in findings(inventory, "unknown_provider_field")]
    assert sorted((f.observations.get()["part"], f.observations.get()["parts"], len(f.locators)) for f in parts) == [
        (1, 3, 3), (2, 3, 3), (3, 3, 1)]
    assert {f.observations.get()["occurrences"] for f in parts} == {7}
    located = sorted((x for f in parts for x in f.locators), key=m.Locator.sort_key)
    assert [x.pointer for x in located] == [rec(i) + "/new_record_field" for i in range(7)]


@pytest.mark.parametrize("value", [None, "F", "ORF", 0, {"code": "P"}], ids=["null", "text", "orf-text", "number",
                                                                           "object"])
def test_provider_functionality_is_kept_raw_and_never_classified(value):
    blobs = {"airr": with_records(IGL, "airr", lambda records: records[0].update(functional=value))}
    _, inventory = build(IGL, blobs)
    found = observed(inventory, "airr", rec(0) + "/functional")
    assert (found.value.get(), found.interpretation) == (value, "preserved")
    assert records_of(inventory, "airr")[0].status == "parsed"
    assert not [f for f in inventory.findings if m.Locator("airr", rec(0)) in f.subjects]
    assert inventory.claims == () and inventory.mappings == ()
