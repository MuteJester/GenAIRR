"""Portable build dossiers: writing, checking and rebuilding fully offline.

A dossier is written into a new directory, member by member with the atomic file writer,
``manifest.json`` last.  ``ReferenceBuild.read_dossier`` refuses a directory without a manifest,
verifies every member, regenerates the inventory, inputs, runtime reference and report from the
locked bytes and the recipe, and compares them with what is stored.

Expected values come from the fixtures under ``tests/data/reference_sources/``, from literals,
or from bytes the test itself wrote; never from the code under test.  Tampered dossiers are
rehashed with the published digest functions, as an editor hiding a change would do.
"""
from __future__ import annotations

import datetime
import hashlib
import json
import pathlib

import pytest

import GenAIRR
from GenAIRR import _artifact_io
from GenAIRR._reference import artifact as artifact_module
from GenAIRR._reference_sources import build as build_module
from GenAIRR.reference_sources import ReferenceBuild
from tests.test_reference_source_api import (
    ALL, IGL, MINI, error, fixture_bytes, from_files, full_recipe, offline, pinned_asset, sha, tree,  # noqa: F401
)

T0 = datetime.datetime(2026, 10, 27, 9, 0, 0, tzinfo=datetime.timezone.utc)


# ── helpers ──────────────────────────────────────────────────────────


@pytest.fixture(scope="module")
def igl(tmp_path_factory):
    """A two-allele build over the complete pinned IGL release 4 exports."""
    root = tmp_path_factory.mktemp("igl")
    source = from_files(root)
    return source, source.build(full_recipe(root, source))


def recording_writer(monkeypatch, fail_at=None) -> list:
    """Record every atomic write (the dossier writer and the reference artifact writer share one
    implementation); with ``fail_at``, the write with that index raises instead of writing."""
    calls = []
    real = _artifact_io.write_atomic

    def write(path, data, *, replace):
        calls.append(pathlib.Path(path))
        if fail_at is not None and len(calls) - 1 == fail_at:
            raise _artifact_io.ArtifactIOError("io", f"simulated failure writing {pathlib.Path(path).name}")
        real(path, data, replace=replace)

    monkeypatch.setattr(build_module, "write_atomic", write)
    monkeypatch.setattr(artifact_module, "write_atomic", write)
    return calls


def member_names(directory: pathlib.Path, calls) -> list:
    return [p.relative_to(directory).as_posix() for p in calls]


# ── writing ──────────────────────────────────────────────────────────


def test_a_dossier_holds_the_manifest_members_and_receipts(igl, tmp_path):
    source, build = igl
    build.write_dossier(tmp_path / "dossier")
    files = tree(tmp_path / "dossier")
    blobs = {"cache/blobs/sha256/" + pinned_asset(IGL, fmt)["sha256"]: fixture_bytes(IGL, fmt) for fmt in ALL}
    for name, data in blobs.items():
        assert files[name] == data
    assert files["reference.genairr-ref.json"] == build.reference.to_json().encode("utf-8")
    assert json.loads(files["recipe.json"].decode("utf-8")) == full_recipe(tmp_path, source)
    assert json.loads(files["source.lock.json"].decode("utf-8"))["commitment"] == source.lock.commitment
    assert [name for name in files if name.startswith("inputs/")] == [
        "inputs/annotations.json", "inputs/j.fasta", "inputs/v.fasta"]
    assert files["inputs/v.fasta"].startswith(b">IGLV1-36*01\n")
    receipts = sorted(name for name in files if name.startswith("receipts/"))
    assert receipts == [f"receipts/{i}.json" for i in range(len(ALL))]
    assert [json.loads(files[name].decode("utf-8"))["format"] for name in receipts] == list(ALL)
    manifest = json.loads(files["manifest.json"].decode("utf-8"))
    listed = {member["name"]: member for member in manifest["members"]}
    assert set(files) == set(listed) | {"manifest.json"} | set(receipts)
    for name, member in listed.items():
        assert (member["sha256"], member["bytes"]) == (sha(files[name]), len(files[name]))
    for name in ("manifest.json", "report.json", "source.lock.json", "recipe.json", "inputs/annotations.json",
                 *receipts):
        text = files[name].decode("utf-8")
        assert text == json.dumps(json.loads(text), ensure_ascii=False, sort_keys=True, separators=(",", ":")) + "\n"
    for data in files.values():
        assert b"\r" not in data or data in blobs.values()


def test_members_are_published_with_the_atomic_writer_and_the_manifest_last(igl, tmp_path, monkeypatch):
    _, build = igl
    calls = recording_writer(monkeypatch)
    build.write_dossier(tmp_path / "dossier")
    names = member_names(tmp_path / "dossier", calls)
    blobs = sorted("cache/blobs/sha256/" + pinned_asset(IGL, fmt)["sha256"] for fmt in ALL)
    assert sorted(names[:4]) == blobs
    assert names[4:] == ["source.lock.json", "recipe.json", "inputs/v.fasta", "inputs/j.fasta",
                         "inputs/annotations.json", "reference.genairr-ref.json", "report.json",
                         "receipts/0.json", "receipts/1.json", "receipts/2.json", "receipts/3.json", "manifest.json"]


def test_a_dossier_is_written_into_a_new_directory_only(igl, tmp_path):
    _, build = igl
    existing = tmp_path / "existing"
    existing.mkdir()
    (existing / "notes.txt").write_bytes(b"kept")
    exc = error("io", build.write_dossier, existing)
    assert "new directory" in exc.message
    assert tree(existing) == {"notes.txt": b"kept"}
    error("io", build.write_dossier, tmp_path / "missing-parent" / "dossier")
    assert not (tmp_path / "missing-parent").exists()
    (tmp_path / "a-file").write_bytes(b"kept")
    error("io", build.write_dossier, tmp_path / "a-file")
    assert (tmp_path / "a-file").read_bytes() == b"kept"
    with pytest.raises(TypeError):
        build.write_dossier(42)


def test_a_second_writer_fails_before_writing_anything(igl, tmp_path, monkeypatch):
    _, build = igl
    build.write_dossier(tmp_path / "dossier")
    before = tree(tmp_path / "dossier")
    calls = recording_writer(monkeypatch)
    error("io", build.write_dossier, tmp_path / "dossier")
    assert calls == []
    assert tree(tmp_path / "dossier") == before


def test_a_failed_write_leaves_an_incomplete_dossier_and_removes_nothing(igl, tmp_path, monkeypatch):
    _, build = igl
    complete = tmp_path / "complete"
    build.write_dossier(complete)
    expected = tree(complete)
    steps = len(expected)
    assert steps == 4 + 7 + 4 + 1          # blobs, lock to report, receipts, manifest
    (tmp_path / "foreign.txt").write_bytes(b"another writer's file")
    for step in range(steps):
        directory = tmp_path / f"failed-{step}"
        calls = recording_writer(monkeypatch, fail_at=step)
        exc = error("io", build.write_dossier, directory)
        failed = member_names(directory, calls)[step]
        assert exc.path == "/" + failed
        assert "incomplete" in exc.message and "manifest.json" in exc.message
        written = tree(directory)
        assert "manifest.json" not in written
        assert set(written) == set(member_names(directory, calls)[:step])
        for name, data in written.items():
            assert data == expected[name]
        assert (tmp_path / "foreign.txt").read_bytes() == b"another writer's file"
