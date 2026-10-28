"""The public published-source API: ``GenAIRR.reference_sources``.

``OgrdbSource`` acquires (``fetch``, the only network path), records (``from_files``) or reads
(``read``) one locked OGRDB germline-set release, writes its lock, inspects it and builds an
explicit recipe into a ``ReferenceBuild``.  Builds are checked against the projected reference
and their build report.

Expected values come from the independently authored fixtures under
``tests/data/reference_sources/`` (``sources.json``, ``expectations.json``, the selection
fixture), from literals, or from reading raw payload fields with the standard library; never
from the code under test.  No test touches the network: ``fetch`` runs over a fake transport
injected through the private transport factory, and every offline operation runs with sockets
and the transport factory disabled.
"""
from __future__ import annotations

import datetime
import hashlib
import http.client
import json
import pathlib
import re
import shutil
import socket

import pytest

import GenAIRR
from GenAIRR import reference_sources as public
from GenAIRR._reference_sources import api
from GenAIRR._reference_sources import model as m
from GenAIRR._reference_sources.acquisition import HttpsTransport
from GenAIRR._reference_sources.errors import KINDS, ReferenceSourceError
from GenAIRR.reference_sources import OgrdbSource, ReferenceBuild
from tests.test_ogrdb_source_coordinates import in_scope, small
from tests.test_reference_source_acquisition import FakeResponse, FakeTransport, serve
from tests.test_reference_source_projection import MINI, entry, recipe

DATA = pathlib.Path(__file__).resolve().parent / "data" / "reference_sources"
SOURCES = json.loads((DATA / "sources.json").read_bytes().decode("utf-8"))
SUFFIX = {"airr": "json", "airr_ex": "json", "gapped": "fasta", "gapped_ex": "fasta"}
PROFILE = "ogrdb_api_v2_human_ig_observed_v1"
ALL = ("airr", "airr_ex", "gapped", "gapped_ex")
IGH, IGK, IGL = ("9606.IGH_VDJ", 10), ("9606.IGKappa_VJ", 5), ("9606.IGLambda_VJ", 4)
T0 = datetime.datetime(2026, 10, 27, 9, 0, 0, tzinfo=datetime.timezone.utc)
T1 = datetime.datetime(2026, 10, 28, 17, 30, 0, tzinfo=datetime.timezone.utc)


# ── helpers ──────────────────────────────────────────────────────────


def fixture_path(source, fmt) -> pathlib.Path:
    return DATA / "ogrdb" / source[0] / str(source[1]) / f"{fmt}.{SUFFIX[fmt]}"


def fixture_bytes(source, fmt) -> bytes:
    return fixture_path(source, fmt).read_bytes()


def pinned_asset(source, fmt) -> dict:
    (entry_,) = [s for s in SOURCES["sets"] if (s["set_id"], s["release"]) == source]
    (asset,) = [a for a in entry_["assets"] if a["format"] == fmt]
    return asset


def sha(data: bytes) -> str:
    return "sha256:" + hashlib.sha256(data).hexdigest()


def from_files(tmp_path, source=IGL, formats=ALL, *, cache="cache", files=None) -> OgrdbSource:
    exports = files if files is not None else {fmt: fixture_path(source, fmt) for fmt in formats}
    return OgrdbSource.from_files(set_id=source[0], release=source[1], profile=PROFILE, exports=exports,
                                  cache=tmp_path / cache)


def written(tmp_path, blobs: dict, folder="snapshot") -> dict:
    """``blobs`` (format -> bytes) written as local files; returns ``{format: path}``."""
    directory = tmp_path / folder
    directory.mkdir(parents=True, exist_ok=True)
    paths = {}
    for fmt, data in blobs.items():
        paths[fmt] = directory / f"{fmt}.{SUFFIX[fmt]}"
        paths[fmt].write_bytes(data)
    return paths


def mini_source(tmp_path, monkeypatch, labels=MINI, **edits):
    """A small synthetic IGL source over ``labels``, registered as validated for this test.

    Returns ``(source, blobs)``.  The lock registered in scope is the one ``from_files`` seals
    (local files: the scope compares digests and sizes only)."""
    from tests.test_ogrdb_source_inventory import seal as seal_lock
    blobs = small(IGL, set(labels), **edits)
    in_scope(monkeypatch, seal_lock(IGL, blobs))
    return from_files(tmp_path, files=written(tmp_path, blobs)), blobs


def full_recipe(tmp_path, source: OgrdbSource, labels=MINI, **kwargs) -> dict:
    """A recipe selecting ``labels`` from the pinned IGL release 4 ``airr`` export."""
    blobs = {fmt: fixture_bytes(IGL, fmt) for fmt in ALL}
    return recipe(source.lock, [entry(blobs, "airr", label) for label in labels], **kwargs)


def error(kind, call, *args, **kwargs) -> ReferenceSourceError:
    with pytest.raises(ReferenceSourceError) as info:
        call(*args, **kwargs)
    assert info.value.kind == kind, f"expected {kind}, got {info.value.kind}: {info.value}"
    return info.value


def tree(directory: pathlib.Path) -> dict:
    """Every file below ``directory`` (relative POSIX name -> bytes)."""
    return {p.relative_to(directory).as_posix(): p.read_bytes() for p in sorted(directory.rglob("*"))
            if p.is_file()}


@pytest.fixture
def offline(monkeypatch):
    """Any attempt to open a connection, or to create the fetch transport, fails the test."""
    def refuse(*args, **kwargs):
        pytest.fail("an offline operation tried to reach the network")

    monkeypatch.setattr(socket, "create_connection", refuse)
    monkeypatch.setattr(socket, "getaddrinfo", refuse)
    monkeypatch.setattr(socket.socket, "connect", refuse)
    monkeypatch.setattr(http.client.HTTPSConnection, "connect", refuse)
    monkeypatch.setattr(http.client.HTTPConnection, "connect", refuse)
    monkeypatch.setattr(HttpsTransport, "get", refuse)
    monkeypatch.setattr(api, "_transport", refuse)


@pytest.fixture(scope="module")
def igl_build(tmp_path_factory):
    """A two-allele build over the complete pinned IGL release 4 exports (no scope patch)."""
    root = tmp_path_factory.mktemp("igl")
    source = from_files(root)
    return source, source.build(full_recipe(root, source))


# ── the public module ────────────────────────────────────────────────


def test_the_public_module_exports_exactly_the_source_api():
    assert public.__all__ == ("OgrdbSource", "ReferenceBuild", "ReferenceSourceError")
    assert public.ReferenceSourceError is ReferenceSourceError
    assert public.OgrdbSource is api.OgrdbSource
    assert issubclass(ReferenceSourceError, ValueError)
    assert KINDS == ("io", "network", "resource_limit", "malformed", "unsupported_schema", "unsupported_profile",
                     "identity_mismatch", "digest_mismatch", "missing_asset", "incomplete_dossier",
                     "unresolved_selection", "stale_decision", "projection", "capability")


def test_sources_and_builds_come_only_from_their_operations():
    with pytest.raises(TypeError):
        OgrdbSource()
    with pytest.raises(TypeError):
        ReferenceBuild()


# ── local snapshots ──────────────────────────────────────────────────


def test_from_files_records_the_supplied_bytes_without_claiming_a_download(tmp_path, offline):
    source = from_files(tmp_path)
    assert (source.set_id, source.release, source.profile) == ("9606.IGLambda_VJ", 4, PROFILE)
    lock = source.lock
    assert lock.asset_ids == ALL
    for asset in lock.assets:
        pinned = pinned_asset(IGL, asset.format)
        assert (asset.sha256, asset.bytes) == ("sha256:" + pinned["sha256"], pinned["bytes"])
        assert (asset.origin_kind, asset.declared_url) == ("local_file", None)
    assert isinstance(source.receipts, tuple)
    assert [(r.format, r.origin, r.sha256) for r in source.receipts] == [
        (a.format, "local_file", a.sha256) for a in lock.assets]
    for fmt in ALL:                       # stored exactly, content-addressed
        assert (tmp_path / "cache" / "blobs" / "sha256" / pinned_asset(IGL, fmt)["sha256"]).read_bytes() \
            == fixture_bytes(IGL, fmt)


def test_a_source_is_immutable(tmp_path):
    source = from_files(tmp_path)
    for name in ("set_id", "release", "profile", "lock", "receipts"):
        with pytest.raises(AttributeError):
            setattr(source, name, None)
    with pytest.raises(AttributeError):
        source._assets = {}
    assert not hasattr(source, "assets")


def test_a_moved_and_renamed_snapshot_has_the_same_lock(tmp_path):
    first = from_files(tmp_path)
    moved = {fmt: tmp_path / "elsewhere" / f"renamed-{i}.bin" for i, fmt in enumerate(ALL)}
    moved[ALL[0]].parent.mkdir()
    for fmt, path in moved.items():
        shutil.copyfile(fixture_path(IGL, fmt), path)
    second = from_files(tmp_path, cache="other-cache", files=moved)
    assert second.lock == first.lock
    first.write_lock(tmp_path / "a.lock.json")
    second.write_lock(tmp_path / "b.lock.json")
    text = (tmp_path / "a.lock.json").read_bytes()
    assert text == (tmp_path / "b.lock.json").read_bytes()
    for fragment in (str(tmp_path).encode("utf-8"), b"renamed-", b"airr.json", b"elsewhere"):
        assert fragment not in text


def test_write_lock_writes_canonical_lf_text_and_never_replaces_silently(tmp_path, offline):
    source = from_files(tmp_path)
    path = tmp_path / "source.lock.json"
    source.write_lock(path)
    data = path.read_bytes()
    assert data.endswith(b"}\n") and data.count(b"\n") == 1 and b"\r" not in data
    document = json.loads(data.decode("utf-8"))
    assert data.decode("utf-8") == json.dumps(document, ensure_ascii=False, sort_keys=True,
                                              separators=(",", ":")) + "\n"
    assert document["artifact"] == "genairr_reference_source_lock"
    assert document["commitment"] == source.lock.commitment
    error("io", source.write_lock, path)
    assert path.read_bytes() == data
    with pytest.raises(TypeError):
        source.write_lock(path, replace="yes")
    source.write_lock(path, replace=True)
    assert path.read_bytes() == data


def test_read_checks_the_lock_and_every_blob_offline(tmp_path, offline):
    source = from_files(tmp_path)
    source.write_lock(tmp_path / "source.lock.json")
    again = OgrdbSource.read(tmp_path / "source.lock.json", cache=tmp_path / "cache")
    assert again.lock == source.lock
    assert [r.origin for r in again.receipts] == ["verified_cache"] * len(ALL)


def test_an_offline_cache_miss_or_corrupt_blob_is_a_typed_error(tmp_path, offline):
    source = from_files(tmp_path)
    source.write_lock(tmp_path / "source.lock.json")
    exc = error("missing_asset", OgrdbSource.read, tmp_path / "source.lock.json", cache=tmp_path / "empty")
    assert "reading never downloads" in exc.message
    blob = tmp_path / "cache" / "blobs" / "sha256" / pinned_asset(IGL, "gapped")["sha256"]
    blob.write_bytes(blob.read_bytes()[:-1] + b"X")
    exc = error("digest_mismatch", OgrdbSource.read, tmp_path / "source.lock.json", cache=tmp_path / "cache")
    assert exc.path == "/assets/2"


def test_a_lock_of_another_schema_version_is_unsupported(tmp_path):
    source = from_files(tmp_path)
    document = source.lock.to_dict()
    document["schema_version"] = 2
    (tmp_path / "future.lock.json").write_bytes(json.dumps(document).encode("utf-8"))
    error("unsupported_schema", OgrdbSource.read, tmp_path / "future.lock.json", cache=tmp_path / "cache")


# ── acquisition (the only network path, over a fake transport) ───────


def provider_url(source, fmt) -> str:
    return f"https://ogrdb.airr-community.org/api_v2/germline/set/{source[0]}/{source[1]}/{fmt}"


def fake(monkeypatch, routes) -> FakeTransport:
    transport = FakeTransport(routes)
    monkeypatch.setattr(api, "_transport", lambda: transport)
    return transport


def fetch(tmp_path, **kwargs) -> OgrdbSource:
    kwargs.setdefault("formats", ALL)
    return OgrdbSource.fetch(set_id=IGL[0], release=IGL[1], profile=PROFILE, cache=tmp_path / "cache", **kwargs)


def test_the_fetch_transport_is_the_standard_library_client():
    assert type(api._transport()) is HttpsTransport
    doc = OgrdbSource.fetch.__doc__
    assert "proxy" in doc and "from_files" in doc


def test_fetch_downloads_exactly_the_requested_formats(tmp_path, monkeypatch):
    transport = fake(monkeypatch, serve(IGL))
    source = fetch(tmp_path)
    assert transport.urls == [provider_url(IGL, fmt) for fmt in ALL]
    for asset in source.lock.assets:
        pinned = pinned_asset(IGL, asset.format)
        assert (asset.sha256, asset.bytes, asset.origin_kind, asset.declared_url) == (
            "sha256:" + pinned["sha256"], pinned["bytes"], "provider_download", provider_url(IGL, asset.format))
    assert [r.origin for r in source.receipts] == ["provider_download"] * len(ALL)


def test_fetch_defaults_to_every_validated_export(tmp_path, monkeypatch):
    """Without the gapped exports no V junction anchor can be accepted, so the default asks for all four."""
    transport = fake(monkeypatch, serve(IGL))
    source = OgrdbSource.fetch(set_id=IGL[0], release=IGL[1], profile=PROFILE, cache=tmp_path / "cache")
    assert transport.urls == [provider_url(IGL, fmt) for fmt in ("airr", "airr_ex", "gapped", "gapped_ex")]
    assert source.lock.asset_ids == ("airr", "airr_ex", "gapped", "gapped_ex")
    assert source.lock.missing_assets == ()


@pytest.mark.parametrize("override, kind, path", [
    ({"release": "latest"}, "malformed", "/release"),
    ({"release": True}, "malformed", "/release"),
    ({"release": 0}, "malformed", "/release"),
    ({"release": -4}, "malformed", "/release"),
    ({"profile": "ogrdb_latest"}, "unsupported_profile", "/profile"),
    ({"formats": ("ungapped",)}, "unsupported_profile", "/formats/0"),
])
def test_fetch_rejects_implicit_or_unsupported_requests_before_any_transport_call(tmp_path, monkeypatch,
                                                                                   override, kind, path):
    transport = fake(monkeypatch, serve(IGL))
    arguments = dict(set_id=IGL[0], release=IGL[1], profile=PROFILE, cache=tmp_path / "cache", formats=ALL)
    arguments.update(override)
    exc = error(kind, OgrdbSource.fetch, **arguments)
    assert exc.path == path
    assert transport.calls == []
    assert not (tmp_path / "cache").exists()


def test_fetch_against_an_expected_lock_checks_the_reacquired_bytes(tmp_path, monkeypatch):
    fake(monkeypatch, serve(IGL))
    first = fetch(tmp_path)
    first.write_lock(tmp_path / "pin.lock.json")
    transport = fake(monkeypatch, serve(IGL))
    again = fetch(tmp_path, expected_lock=tmp_path / "pin.lock.json")
    assert again.lock == first.lock
    assert transport.urls == [provider_url(IGL, fmt) for fmt in ALL]      # reacquired, never trusted
    again = fetch(tmp_path, expected_lock=first.lock)
    assert again.lock == first.lock


def test_changed_bytes_under_the_same_release_are_drift_and_change_nothing(tmp_path, monkeypatch):
    fake(monkeypatch, serve(IGL))
    first = fetch(tmp_path)
    first.write_lock(tmp_path / "pin.lock.json")
    before = (tmp_path / "pin.lock.json").read_bytes(), tree(tmp_path / "cache")
    routes = serve(IGL)
    changed = fixture_bytes(IGL, "airr").replace(b"IGLV1-36*01", b"IGLV1-36*09", 1)
    routes[provider_url(IGL, "airr")] = lambda: FakeResponse(200, changed)
    fake(monkeypatch, routes)
    exc = error("digest_mismatch", fetch, tmp_path, expected_lock=tmp_path / "pin.lock.json")
    assert "drift" in exc.message and sha(changed) in exc.message
    assert ((tmp_path / "pin.lock.json").read_bytes(), tree(tmp_path / "cache")) == before


def test_an_expected_lock_of_the_wrong_type_is_a_programming_error(tmp_path, monkeypatch):
    transport = fake(monkeypatch, serve(IGL))
    with pytest.raises(TypeError):
        fetch(tmp_path, expected_lock=42)
    assert transport.calls == []


# ── building ─────────────────────────────────────────────────────────


def test_build_takes_a_recipe_document_or_a_recipe_file(tmp_path, offline):
    source = from_files(tmp_path)
    document = full_recipe(tmp_path, source)
    path = tmp_path / "igl.recipe.json"
    path.write_bytes(json.dumps(document, indent=2).encode("utf-8"))
    from_document, from_file = source.build(document), source.build(path)
    assert from_document.reference.to_json() == from_file.reference.to_json()
    assert from_document.recipe == from_file.recipe == m.Recipe.from_dict(document)
    assert from_document.lock == source.lock
    labels = [a.label for a in from_document.reference.alleles]
    assert labels == list(MINI)
    with pytest.raises(TypeError):
        source.build(42)


def test_a_build_exposes_read_only_results(igl_build):
    source, build = igl_build
    assert build.report.source_lock == build.lock == source.lock
    assert build.report.recipe == build.recipe
    runtime = build.report.runtime
    assert (runtime.catalogue_hash, runtime.models_hash, runtime.execution_hash, runtime.artifact_hash) == (
        build.reference.catalogue_hash, build.reference.models_hash, build.reference.execution_hash,
        build.reference.artifact_hash)
    assert build.report.generator.version == GenAIRR.__version__
    assert build.manifest.runtime == runtime
    assert build.manifest.builder == build.report.generator
    assert build.checked_with is None
    for name in ("reference", "report", "lock", "recipe", "manifest", "checked_with"):
        with pytest.raises(AttributeError):
            setattr(build, name, None)


def test_a_build_manifest_lists_every_deterministic_member(igl_build):
    source, build = igl_build
    names = [member.name for member in build.manifest.members]
    blobs = sorted("cache/blobs/sha256/" + pinned_asset(IGL, fmt)["sha256"] for fmt in ALL)
    assert names == blobs + ["inputs/annotations.json", "inputs/j.fasta", "inputs/v.fasta", "recipe.json",
                             "reference.genairr-ref.json", "report.json", "source.lock.json"]
    by_name = {member.name: member for member in build.manifest.members}
    for fmt in ALL:
        pinned = pinned_asset(IGL, fmt)
        assert (by_name["cache/blobs/sha256/" + pinned["sha256"]].sha256,
                by_name["cache/blobs/sha256/" + pinned["sha256"]].bytes) == ("sha256:" + pinned["sha256"],
                                                                              pinned["bytes"])
    reference = build.reference.to_json().encode("utf-8")
    assert (by_name["reference.genairr-ref.json"].sha256, by_name["reference.genairr-ref.json"].bytes) == (
        sha(reference), len(reference))


def test_build_errors_keep_recipe_paths_and_name_the_recipe_file(tmp_path, offline):
    source = from_files(tmp_path)
    document = full_recipe(tmp_path, source)
    document["selection"][1]["sequence_sha256"] = "sha256:" + "0" * 64
    exc = error("unresolved_selection", source.build, document)
    assert exc.path.startswith("/selection/1")
    path = tmp_path / "wrong.recipe.json"
    path.write_bytes(json.dumps(document).encode("utf-8"))
    exc = error("unresolved_selection", source.build, path)
    assert exc.path.startswith("/selection/1") and "wrong.recipe.json" in exc.message
    path.write_bytes(b'{"artifact": "genairr_reference_recipe", "artifact": 1}')
    exc = error("malformed", source.build, path)
    assert "wrong.recipe.json" in exc.message
    document = full_recipe(tmp_path, source)
    document["projection"]["version"] = 2
    path.write_bytes(json.dumps(document).encode("utf-8"))
    exc = error("unsupported_profile", source.build, path)
    assert exc.path == "/projection" and "wrong.recipe.json" in exc.message
    exc = error("io", source.build, tmp_path / "absent.recipe.json")
    assert "absent.recipe.json" in exc.message


def test_a_build_compares_both_record_exports(tmp_path, offline):
    source = from_files(tmp_path, formats=("airr", "gapped"))
    exc = error("missing_asset", source.build, full_recipe(tmp_path, source))
    assert exc.path == "/source_lock/assets" and "airr_ex" in exc.message


def test_a_failed_capability_requirement_is_a_typed_build_failure(tmp_path, offline):
    source = from_files(tmp_path)
    exc = error("capability", source.build,
                full_recipe(tmp_path, source, requires=(("complete_reference_span", "assertion"),)))
    assert exc.path == "/requires/0" and "complete_reference_span" in exc.message


# ── inspection ───────────────────────────────────────────────────────

EXPECT = json.loads((DATA / "expectations.json").read_bytes().decode("utf-8"))


def raw_records(source, fmt) -> list:
    return json.loads(fixture_bytes(source, fmt).decode("utf-8"))["GermlineSet"][0]["allele_descriptions"]


def test_inspection_reports_the_complete_source(tmp_path, offline):
    source = from_files(tmp_path)
    inspection = source.inspect()
    report = inspection.to_dict()
    assert (report["provider"], report["set_id"], report["release"]) == ("ogrdb", "9606.IGLambda_VJ", 4)
    assert report["identity"] == {"taxon": {"scientific_name": "Homo sapiens", "taxonomy_id": 9606}, "locus": "IGL",
                                  "chain_type": "vj"}
    assert report["adapter"]["coordinate_profile"] == {"id": PROFILE, "version": 1}
    assert report["lock_hash"] == source.lock.commitment
    assert report["scope"] == {"in_scope": True, "problem": None}
    assert [(a["asset_id"], a["sha256"], a["bytes"]) for a in report["assets"]] == [
        (fmt, "sha256:" + pinned_asset(IGL, fmt)["sha256"], pinned_asset(IGL, fmt)["bytes"]) for fmt in ALL]
    assert report["missing_assets"] == []
    assert report["comparison"]["complete"] is True
    assert (report["comparison"]["exports"], report["comparison"]["missing"]) == (["airr", "airr_ex"], [])
    expected = EXPECT["sets"]["9606.IGLambda_VJ/4"]
    for counts in report["exports"]:
        stated = expected[counts["asset"]]
        assert counts["record_occurrences"] == {"V": stated["V"], "D": stated["D"], "J": stated["J"], "unassigned": 0}
        assert sum(counts["mentioned_identities"].values()) == len(stated["paralog_only_identities"])
    records = [r for fmt in ("airr", "airr_ex") for r in raw_records(IGL, fmt)]
    assert len(report["records"]) == len(records) == 98 + 97
    assert [r["label"] for r in report["records"]] == [r["label"] for r in records]
    assert [r["segment"] for r in report["records"]] == [r["sequence_type"] for r in records]
    assert {r["status"] for r in report["records"]} == {"parsed"}
    assert report["accepted_coordinate_claims"] > 0
    assert {(c["predicate"], c["status"]) for c in report["claims"]} >= {("anchor", "accepted"),
                                                                         ("functionality", "unresolved")}
    assert all(f["count"] > 0 and f["severity"] in ("info", "warning", "error") for f in report["findings"])
    assert inspection.comparison_complete is True and inspection.in_scope is True


def test_inspection_is_an_immutable_copy_returning_report(tmp_path):
    source = from_files(tmp_path)
    inspection = source.inspect()
    first = inspection.to_dict()
    first["records"].clear()
    first["scope"]["in_scope"] = False
    assert inspection.to_dict() != first
    assert inspection.to_dict() == source.inspect().to_dict()
    json.dumps(inspection.to_dict(), allow_nan=False)
    with pytest.raises(dataclasses_frozen_error()):
        inspection.set_id = "other"
    assert str(inspection) == str(source.inspect())


def dataclasses_frozen_error():
    import dataclasses
    return dataclasses.FrozenInstanceError


def test_the_inspection_text_covers_identity_assets_counts_and_every_record(tmp_path):
    source = from_files(tmp_path)
    text = str(source.inspect())
    assert "9606.IGLambda_VJ release 4" in text and source.lock.commitment in text
    assert PROFILE in text and "validated" in text
    for fmt in ALL:
        assert "sha256:" + pinned_asset(IGL, fmt)["sha256"] in text
    assert "comparison: complete" in text.lower()
    for record in raw_records(IGL, "airr") + raw_records(IGL, "airr_ex"):
        assert record["label"] in text
    lines = text.splitlines()
    assert sum(1 for line in lines if "#/GermlineSet/0/allele_descriptions/" in line) == 98 + 97
    assert all(len(line) <= 200 for line in lines)


def test_a_single_export_can_be_inspected_but_its_comparison_is_incomplete(tmp_path):
    source = from_files(tmp_path, formats=("airr",))
    inspection = source.inspect()
    report = inspection.to_dict()
    assert report["comparison"]["complete"] is False
    assert (report["comparison"]["exports"], report["comparison"]["missing"]) == (["airr"], ["airr_ex"])
    assert [(m["format"], m["reason"]) for m in report["missing_assets"]] == [
        ("airr_ex", "not_supplied"), ("gapped", "not_supplied"), ("gapped_ex", "not_supplied")]
    assert len(report["records"]) == 98
    assert inspection.comparison_complete is False
    text = str(inspection).lower()
    assert "comparison: incomplete" in text and "airr_ex" in text


def test_an_out_of_scope_release_is_inspected_and_reported_as_such(tmp_path):
    blobs = small(IGL, set(MINI))                         # bytes the profile was never validated on
    source = from_files(tmp_path, files=written(tmp_path, blobs))
    inspection = source.inspect()
    report = inspection.to_dict()
    assert report["scope"]["in_scope"] is False and "not the validated" in report["scope"]["problem"]
    assert report["accepted_coordinate_claims"] == 0
    assert not [c for c in report["claims"] if c["predicate"] in ("anchor", "coding_phase") and
                c["status"] == "accepted"]
    assert inspection.in_scope is False
    text = str(inspection)
    assert "outside the validated scope" in text and "no coordinate claim is accepted" in text.lower()
    error("unsupported_profile", source.build, full_recipe(tmp_path, source))


# ── the readable build report ────────────────────────────────────────

SELECTION = json.loads((DATA / "human_igh_release10_candidate_selection.json").read_bytes().decode("utf-8"))
LOSS_CODES = ("v_subregions", "gapped_sequences", "affirmed_spans", "genomic_context")
SECTIONS = ("Source", "Hashes", "Endpoint comparison", "Counts", "Selection", "Source records",
            "Mentioned-only identities", "Claims", "Claim provenance", "Findings", "Runtime capabilities",
            "Required capabilities", "Not in the runtime reference")


@pytest.fixture(scope="module")
def igh_starter(tmp_path_factory):
    """The IGH release 10 starter: the selection fixture as a recipe over the pinned exports."""
    from tests.test_reference_source_projection import interpreted, starter_recipe
    root = tmp_path_factory.mktemp("igh")
    source = from_files(root, IGH)
    blobs = {fmt: fixture_bytes(IGH, fmt) for fmt in ALL}
    build = source.build(starter_recipe(source.lock, interpreted(blobs, source.lock)))
    return build, str(build.report)


def section(text: str, title: str) -> str:
    """The lines of one titled section (up to the next section title)."""
    lines = text.splitlines()
    start = next(i for i, line in enumerate(lines) if line.startswith(title))
    end = next((i for i in range(start + 1, len(lines)) if any(lines[i].startswith(t) for t in SECTIONS)), len(lines))
    return "\n".join(lines[start:end])


def test_the_report_text_has_its_sections_in_order(igh_starter):
    _, text = igh_starter
    lines = text.splitlines()
    starts = [next(i for i, line in enumerate(lines) if line.startswith(title)) for title in SECTIONS]
    assert starts == sorted(starts)


def test_the_report_text_names_the_source_versions_hashes_and_generator(igh_starter):
    build, text = igh_starter
    head = section(text, "Source")
    assert "OGRDB 9606.IGH_VDJ release 10" in head and "Homo sapiens" in head and "IGH" in head
    assert f"GenAIRR {GenAIRR.__version__}" in text
    for fragment in ("ogrdb_germline_set v1", f"{PROFILE} v1", "observed_span v1", "genairr.unfitted_defaults.v1"):
        assert fragment in head
    for fmt in ALL:
        pinned = pinned_asset(IGH, fmt)
        assert f"sha256:{pinned['sha256']}" in head and str(pinned["bytes"]) in head
    hashes = section(text, "Hashes")
    report = build.report
    for value in (*report.hashes.to_dict().values(), *report.runtime.to_dict().values()):
        assert value in hashes
    assert "Reference.report()" in text and "runtime catalogue" in text
    assert "functional allele" not in text.lower()
    generator = next(line for line in text.splitlines() if line.startswith("Generator: "))
    assert "a dossier check verifies content, not this statement" in generator
    assert "the lock writer's statements" in head


def test_the_report_text_accounts_for_every_source_record(igh_starter):
    build, text = igh_starter
    expected = EXPECT["sets"]["9606.IGH_VDJ/10"]
    total = expected["airr"]["records"] + expected["airr_ex"]["records"]
    selected = sum(len(labels) for labels in SELECTION["selected"].values())
    records = section(text, "Source records")
    assert f"Source records: {total} (selected {selected}, excluded {total - selected}, rejected 0)" in records
    for fmt in ("airr", "airr_ex"):
        for index in range(expected[fmt]["records"]):
            locator = f"{fmt}#/GermlineSet/0/allele_descriptions/{index}"
            assert f"{locator} " in text or f"{locator}/" in text or text.endswith(locator), locator
    for outcome in build.report.records:
        if outcome.exclusion:
            assert outcome.exclusion in records
        if outcome.detail:
            assert outcome.detail in records, outcome.record.id


def test_the_report_text_lists_the_ordered_selection(igh_starter):
    _, text = igh_starter
    lines = [line.split()[1] for line in section(text, "Selection").splitlines()[1:]]
    assert lines == [f"{segment}:{label}" for segment in ("V", "D", "J") for label in SELECTION["selected"][segment]]


def test_the_report_text_counts_every_claim_and_points_to_the_provenance_of_unselected_ones(igh_starter):
    build, text = igh_starter
    claims = section(text, "Claims")
    tallies = {}
    for outcome, claim in zip(build.report.claims, build.report.inventory.claims):
        key = (outcome.outcome, claim.predicate, claim.status)
        tallies[key] = tallies.get(key, 0) + 1
    titles = {"projected": "projected", "withheld": "withheld", "not_representable": "not representable at runtime",
              "not_selected": "not selected"}
    for outcome, title in titles.items():
        expected = sorted(f"{p} {s} {n}" for (o, p, s), n in tallies.items() if o == outcome)
        (line,) = [line.strip() for line in claims.splitlines() if line.strip().startswith(title + ":")]
        assert sorted(line[len(title) + 1:].strip().split(", ")) == expected, outcome
    flat = " ".join(claims.split())
    assert "full provenance" in flat and "to_dict()" in flat and "report.json" in flat


def raw_gapped_index(fmt: str, label: str) -> int:
    headers = [chunk.partition("\n")[0].rstrip("\r")
               for chunk in fixture_bytes(IGH, fmt).decode("utf-8").split(">")[1:]]
    (index,) = [i for i, header in enumerate(headers) if header == label]
    return index


RECORD_FIELD = re.compile(r"(/GermlineSet/0/allele_descriptions/[0-9]+)/(.+)")


def expected_evidence(locators) -> str:
    """Evidence as the report states it: ``asset#<record pointer>: field, field; asset#pointer``, with
    consecutive fields of one record export record grouped under it."""
    groups = []
    for locator in locators:
        match = RECORD_FIELD.fullmatch(locator.pointer) if locator.asset in ("airr", "airr_ex") else None
        if match is None:
            groups.append((None, f"{locator.asset}#{locator.pointer}"))
        elif groups and groups[-1][0] == f"{locator.asset}#{match.group(1)}":
            groups[-1][1].append(match.group(2))
        else:
            groups.append((f"{locator.asset}#{match.group(1)}", [match.group(2)]))
    return "; ".join(text if head is None else f"{head}: {', '.join(text)}" for head, text in groups)


def provenance_lines(text: str) -> dict:
    """``{claim id: line}`` of the claim provenance section; every claim id appears on exactly one line."""
    lines = {}
    for line in section(text, "Claim provenance").splitlines():
        if line.startswith("    claim:"):
            claim_id = line.split()[0]
            assert claim_id not in lines, claim_id
            lines[claim_id] = line
    return lines


def test_every_claim_of_a_selected_allele_has_its_provenance_line(igh_starter):
    build, text = igh_starter
    report = build.report
    lines = provenance_lines(text)
    shown = [(o, c) for o, c in zip(report.claims, report.inventory.claims) if o.outcome != "not_selected"]
    assert set(lines) == {c.id for _, c in shown} and len(shown) == 394 + 690 + 427
    for outcome, claim in shown:
        line = lines[claim.id]
        value = json.dumps(claim.value.get(), ensure_ascii=False, sort_keys=True, separators=(",", ":"))
        assert f" {claim.predicate} = {value} " in line
        assert f"({claim.coordinate_system or 'no coordinate system'})" in line
        assert f"[{claim.status}, support {claim.support}, {claim.derivation}]" in line
        assert f"-> {outcome.outcome}" in line
        if outcome.decision is not None:
            assert f"recipe decision {outcome.decision}" in line
        if outcome.reason:
            assert line.endswith(f"reason: {outcome.reason}")
        assert f"; evidence {expected_evidence(claim.evidence)}; mappings " in line
        for mapping_id in claim.mappings:
            mapping = report.inventory.mapping(mapping_id)
            assert f"{mapping_id} ({mapping.kind}, {mapping.validation})" in line
    assert "functional allele" not in text.lower()


def test_a_selected_v_anchor_line_cites_its_raw_source_fields(igh_starter):
    """IGHV1-2*02: the anchor is cdr3_start - 4 of its own airr record, bound to its own gapped record."""
    _, text = igh_starter
    records = raw_records(IGH, "airr")
    (index,) = [i for i, r in enumerate(records) if r["label"] == "IGHV1-2*02"]
    (delineation,) = [j for j, d in enumerate(records[index]["v_gene_delineations"]) if d["delineation_scheme"] == "IMGT"]
    cdr3_start = records[index]["v_gene_delineations"][delineation]["cdr3_start"]
    block = section(text, "Claim provenance").split("  V:IGHV1-2*02 ", 1)[1].split("\n  V:", 1)[0]
    (anchor,) = [line for line in block.splitlines() if " anchor = " in line]
    assert f" anchor = {cdr3_start - 4} (coding_sequence) [accepted, support source, stated] -> projected" in anchor
    record = f"airr#/GermlineSet/0/allele_descriptions/{index}"
    assert (f"evidence {record}: v_gene_delineations/{delineation}/cdr3_start, "
            f"v_gene_delineations/{delineation}/delineation_scheme; gapped#/{raw_gapped_index('gapped', 'IGHV1-2*02')}"
            in anchor)
    assert "(gapped_alignment, exact)" in anchor
    assert {line.split()[1] for line in block.splitlines() if line.startswith("    claim:")} == {
        "functionality", "anchor", "coding_phase", "phase_independence", "boundary_5", "boundary_3", "coding_interval"}


def test_mentioned_only_identities_are_listed_with_the_records_that_name_them(igh_starter):
    _, text = igh_starter
    mentioned = section(text, "Mentioned-only identities")
    for fmt in ("airr", "airr_ex"):
        names = EXPECT["sets"]["9606.IGH_VDJ/10"][fmt]["paralog_only_identities"]
        listed = [line for line in mentioned.splitlines() if line.startswith(f"  {fmt} ")]
        assert len(listed) == len(names)
        for name in names:
            (line,) = [line for line in listed if f":{name}: " in line]
            assert f"paralog named by {fmt}#/GermlineSet/0/allele_descriptions/" in line


def test_the_report_text_shows_capabilities_with_blockers_and_hints(igh_starter):
    build, text = igh_starter
    capabilities = section(text, "Runtime capabilities")
    selected = {segment: len(labels) for segment, labels in SELECTION["selected"].items()}
    for status in build.report.capabilities.statuses:
        assert status.name in capabilities
        for unmet in status.unmet:
            assert unmet.hint in capabilities and unmet.requirement in capabilities
    blockers = selected["V"] + selected["J"] + 2 * selected["D"]          # V 3', J 5', both D ends
    assert f"original_coding_ends: unavailable; {blockers} blockers" in capabilities
    assert "phase support anchor_derived" in capabilities
    required = section(text, "Required capabilities")
    for capability in ("observed_span_simulation", "source_annotated_junction", "canonical_productive"):
        assert f"{capability} at evidence" in required
    assert required.count("satisfied") == 3 and "not satisfied" not in required


def test_the_report_text_names_every_projection_loss(igh_starter):
    from GenAIRR._reference_sources.projection import PROJECTION_LOSSES
    _, text = igh_starter
    losses = section(text, "Not in the runtime reference")
    assert [code for code, _ in PROJECTION_LOSSES] == list(LOSS_CODES)
    for code in LOSS_CODES:
        assert f"{code}:" in losses
    assert "n_v_unannotated_mutations" in " ".join(losses.split())


def test_the_report_text_is_deterministic_and_depends_only_on_the_report(igh_starter):
    from GenAIRR._reference_sources.report import BuildReport
    build, text = igh_starter
    again = BuildReport.from_json(json.dumps(build.report.to_dict()))
    assert str(again) == text
    assert all(line == line.rstrip() for line in text.splitlines())


def grouped(text: str) -> list:
    """``(declared count, listed items)`` for every group header (four-space indent, ``(n)``) or
    ``n records`` line followed by a wrapped list (six-space indent)."""
    import re
    groups, lines = [], text.splitlines()
    for i, line in enumerate(lines):
        if not line.startswith("    ") or line.startswith("      "):
            continue
        found = re.search(r"\] \((\d+)\)|^    (\d+) records|\((\d+)\); would be resolved", line)
        if not found:
            continue
        items = []
        for follow in lines[i + 1:]:
            if not follow.startswith("      "):
                break
            items += [item for item in follow.strip().rstrip(",").split(", ") if item]
        groups.append((int(next(g for g in found.groups() if g)), items))
    return groups


def test_grouped_lists_name_every_record_and_blocker(igh_starter):
    build, text = igh_starter
    report = build.report
    for title in ("Source records", "Runtime capabilities"):
        groups = grouped(section(text, title))
        assert groups and all(declared == len(items) for declared, items in groups), title
    blockers = grouped(section(text, "Runtime capabilities"))
    assert sum(declared for declared, _ in blockers) == sum(len(s.unmet) for s in report.capabilities.statuses)


# ── one inventory build per projection ───────────────────────────────


def test_a_build_reads_the_source_bytes_into_an_inventory_once(tmp_path, monkeypatch):
    from GenAIRR._reference_sources import evidence, ogrdb
    source = from_files(tmp_path)
    document = full_recipe(tmp_path, source)
    calls = []
    real = ogrdb.build_inventory

    def counted(lock, assets):
        calls.append(lock.commitment)
        return real(lock, assets)

    monkeypatch.setattr(ogrdb, "build_inventory", counted)
    monkeypatch.setattr(evidence, "build_inventory", counted)
    source.build(document)
    assert calls == [source.lock.commitment]


def test_interpretation_of_a_supplied_inventory_still_rebuilds_and_compares_it(tmp_path):
    import dataclasses
    from GenAIRR._reference_sources import evidence, ogrdb
    source = from_files(tmp_path)
    blobs = {fmt: fixture_bytes(IGL, fmt) for fmt in ALL}
    inventory = ogrdb.build_inventory(source.lock, blobs)
    assert inventory.findings
    edited = dataclasses.replace(inventory, findings=inventory.findings[:-1])
    exc = error("identity_mismatch", evidence.interpret, source.lock, blobs, edited)
    assert exc.path == "/inventory"


def test_inspection_reads_the_source_bytes_into_an_inventory_once(tmp_path, monkeypatch):
    from GenAIRR._reference_sources import evidence, ogrdb
    source = from_files(tmp_path)
    calls = []
    real = ogrdb.build_inventory

    def counted(lock, assets):
        calls.append(lock.commitment)
        return real(lock, assets)

    monkeypatch.setattr(ogrdb, "build_inventory", counted)
    monkeypatch.setattr(evidence, "build_inventory", counted)
    source.inspect()
    assert calls == [source.lock.commitment]
