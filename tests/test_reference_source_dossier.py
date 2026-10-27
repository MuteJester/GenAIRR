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


def recording_writer(monkeypatch, fail_at=None, *, after=False) -> list:
    """Record every atomic write (the dossier writer and the reference artifact writer share one
    implementation); with ``fail_at``, the write with that index raises instead of writing, or with
    ``after``, raises after publishing its file (as when the directory cannot be synced)."""
    calls = []
    real = _artifact_io.write_atomic

    def write(path, data, *, replace):
        calls.append(pathlib.Path(path))
        failing = fail_at is not None and len(calls) - 1 == fail_at
        if failing and not after:
            raise _artifact_io.ArtifactIOError("io", f"simulated failure writing {pathlib.Path(path).name}")
        real(path, data, replace=replace)
        if failing:
            raise _artifact_io.ArtifactIOError("io", f"simulated failure after publishing {pathlib.Path(path).name}")

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


def test_a_failure_after_the_manifest_is_published_says_the_dossier_is_complete(igl, tmp_path, monkeypatch):
    _, build = igl
    recording_writer(monkeypatch, fail_at=15, after=True)              # the manifest is the 16th write
    exc = error("io", build.write_dossier, tmp_path / "dossier")
    monkeypatch.undo()
    assert exc.path == "/manifest.json"
    assert "has no manifest.json" not in exc.message and "complete" in exc.message
    assert "simulated failure after publishing manifest.json" in exc.message
    assert ReferenceBuild.read_dossier(tmp_path / "dossier").manifest == build.manifest


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
        error("incomplete_dossier", ReferenceBuild.read_dossier, directory)
        assert tree(directory) == written
        assert (tmp_path / "foreign.txt").read_bytes() == b"another writer's file"


# ── checking: round trip and incomplete or unsafe dossiers ───────────


@pytest.mark.usefixtures("offline")
def test_a_written_dossier_reads_back_as_the_stored_build(igl, tmp_path):
    _, build = igl
    build.write_dossier(tmp_path / "dossier")
    before = tree(tmp_path / "dossier")
    checked = ReferenceBuild.read_dossier(tmp_path / "dossier")
    assert tree(tmp_path / "dossier") == before                  # checking writes nothing
    assert checked.reference.to_json() == build.reference.to_json()
    assert checked.report.to_dict() == build.report.to_dict()
    assert checked.manifest == build.manifest
    assert (checked.lock, checked.recipe) == (build.lock, build.recipe)
    assert checked.checked_with == GenAIRR.__version__
    checked.write_dossier(tmp_path / "copy")                     # a checked build writes the same dossier
    assert tree(tmp_path / "copy") == before


def test_a_directory_without_a_manifest_is_an_incomplete_dossier(igl, tmp_path):
    _, build = igl
    build.write_dossier(tmp_path / "dossier")
    (tmp_path / "dossier" / "manifest.json").unlink()
    exc = error("incomplete_dossier", ReferenceBuild.read_dossier, tmp_path / "dossier")
    assert "manifest.json" in exc.message
    (tmp_path / "empty").mkdir()
    error("incomplete_dossier", ReferenceBuild.read_dossier, tmp_path / "empty")


def test_only_an_existing_directory_is_a_dossier(tmp_path):
    error("io", ReferenceBuild.read_dossier, tmp_path / "absent")
    (tmp_path / "file").write_bytes(b"{}")
    error("io", ReferenceBuild.read_dossier, tmp_path / "file")
    with pytest.raises(TypeError):
        ReferenceBuild.read_dossier(None)


def test_a_listed_member_that_is_absent_is_named(igl, tmp_path):
    _, build = igl
    build.write_dossier(tmp_path / "dossier")
    (tmp_path / "dossier" / "inputs" / "j.fasta").unlink()
    exc = error("missing_asset", ReferenceBuild.read_dossier, tmp_path / "dossier")
    assert exc.path == "/inputs/j.fasta"


def test_a_member_that_differs_from_the_manifest_is_a_digest_mismatch(igl, tmp_path):
    _, build = igl
    build.write_dossier(tmp_path / "dossier")
    path = tmp_path / "dossier" / "inputs" / "v.fasta"
    data = path.read_bytes()
    path.write_bytes(data.replace(b"IGLV1-36*01", b"IGLV1-36*02"))           # same size, other digest
    exc = error("digest_mismatch", ReferenceBuild.read_dossier, tmp_path / "dossier")
    assert exc.path == "/inputs/v.fasta" and sha(data) in exc.message
    path.write_bytes(data + b"\n")
    exc = error("digest_mismatch", ReferenceBuild.read_dossier, tmp_path / "dossier")
    assert exc.path == "/inputs/v.fasta" and str(len(data)) in exc.message


def test_unlisted_files_and_directories_are_refused(igl, tmp_path):
    _, build = igl
    build.write_dossier(tmp_path / "dossier")
    for extra in ("notes.txt", "inputs/extra.fasta", "cache/blobs/sha256/" + "0" * 64, "receipts/x.txt",
                  ".report.json.0123456789abcdef.tmp"):
        path = tmp_path / "dossier" / extra
        path.write_bytes(b"unlisted")
        exc = error("malformed", ReferenceBuild.read_dossier, tmp_path / "dossier")
        assert exc.path == "/" + extra
        path.unlink()
    (tmp_path / "dossier" / "inputs" / "nested").mkdir()
    exc = error("malformed", ReferenceBuild.read_dossier, tmp_path / "dossier")
    assert exc.path == "/inputs/nested"


def symlink_or_skip(target: pathlib.Path, link: pathlib.Path, *, directory=False) -> None:
    try:
        link.symlink_to(target, target_is_directory=directory)
    except (OSError, NotImplementedError) as exc:
        pytest.skip(f"symbolic links are unavailable here: {exc}")


def test_a_symbolic_link_anywhere_in_a_dossier_is_refused(mini, tmp_path):
    _, dossier = mini
    member = dossier / "inputs" / "v.fasta"
    outside = tmp_path / "outside.fasta"
    outside.write_bytes(member.read_bytes())
    member.unlink()
    symlink_or_skip(outside, member)
    exc = error("io", ReferenceBuild.read_dossier, dossier)
    assert exc.path == "/inputs/v.fasta" and "symbolic link" in exc.message
    member.unlink()
    member.write_bytes(outside.read_bytes())
    ReferenceBuild.read_dossier(dossier)
    inputs = dossier / "inputs"
    inputs.rename(tmp_path / "real-inputs")
    symlink_or_skip(tmp_path / "real-inputs", inputs, directory=True)
    exc = error("io", ReferenceBuild.read_dossier, dossier)
    assert exc.path == "/inputs"
    inputs.unlink()
    (tmp_path / "real-inputs").rename(inputs)
    manifest = dossier / "manifest.json"
    manifest.rename(tmp_path / "manifest.json")
    symlink_or_skip(tmp_path / "manifest.json", manifest)
    exc = error("io", ReferenceBuild.read_dossier, dossier)
    assert exc.path == "/manifest.json"


def test_declared_sizes_are_bounded_before_any_member_is_read(igl, tmp_path, monkeypatch):
    _, build = igl
    build.write_dossier(tmp_path / "dossier")
    manifest = load_manifest(tmp_path / "dossier")
    for member in manifest["members"]:
        if member["name"] == "report.json":
            member["bytes"] = 256 * 1024 * 1024 + 1
    save_manifest(tmp_path / "dossier", manifest)
    reads = []
    real_open = build_module.os.open
    monkeypatch.setattr(build_module.os, "open", lambda path, *a, **k: reads.append(path) or real_open(path, *a, **k))
    exc = error("resource_limit", ReferenceBuild.read_dossier, tmp_path / "dossier")
    assert exc.path == "/report.json" and "268435457" in exc.message
    assert [pathlib.Path(p).name for p in reads] == ["manifest.json"]


def test_generated_members_must_be_in_canonical_form(mini, tmp_path):
    _, dossier = mini
    report = (dossier / "report.json").read_bytes()
    rewrite(dossier, "report.json", json.dumps(json.loads(report), indent=1).encode("utf-8"))
    exc = error("malformed", ReferenceBuild.read_dossier, dossier)
    assert exc.path == "/report.json" and "canonical" in exc.message
    rewrite(dossier, "report.json", report)
    ReferenceBuild.read_dossier(dossier)
    rewrite(dossier, "reference.genairr-ref.json",
            b"\xef\xbb\xbf" + (dossier / "reference.genairr-ref.json").read_bytes())
    exc = error("malformed", ReferenceBuild.read_dossier, dossier)
    assert exc.path == "/reference.genairr-ref.json"


def test_receipts_are_checked_but_stay_outside_the_manifest(mini, tmp_path):
    build, dossier = mini
    receipt = dossier / "receipts" / "0.json"
    document = json.loads(receipt.read_bytes())
    receipt.write_bytes(b'{"artifact": "genairr_reference_source_receipt"}')
    exc = error("malformed", ReferenceBuild.read_dossier, dossier)
    assert exc.path.startswith("/receipts/0.json")
    document["sha256"] = "sha256:" + "1" * 64
    receipt.write_bytes(canonical(document))
    exc = error("identity_mismatch", ReferenceBuild.read_dossier, dossier)
    assert exc.path == "/receipts/0.json"
    receipt.unlink()
    assert ReferenceBuild.read_dossier(dossier).manifest == build.manifest


# ── checking: rehashed tampering is caught by regeneration ───────────


def canonical(document) -> bytes:
    return (json.dumps(document, ensure_ascii=False, sort_keys=True, separators=(",", ":")) + "\n").encode("utf-8")


def load_manifest(dossier: pathlib.Path) -> dict:
    return json.loads((dossier / "manifest.json").read_bytes().decode("utf-8"))


def save_manifest(dossier: pathlib.Path, manifest: dict) -> None:
    (dossier / "manifest.json").write_bytes(canonical(manifest))


def rewrite(dossier: pathlib.Path, name: str, data: bytes, *, hashes=None, runtime=None) -> None:
    """Replace member ``name`` and rehash the manifest: the member's digest and size, and any given
    manifest ``hashes`` and ``runtime`` fields, as an editor hiding a change would."""
    (dossier / name).write_bytes(data)
    manifest = load_manifest(dossier)
    (member,) = [x for x in manifest["members"] if x["name"] == name]
    member["sha256"], member["bytes"] = sha(data), len(data)
    manifest["hashes"].update(hashes or {})
    manifest["runtime"].update(runtime or {})
    save_manifest(dossier, manifest)


def rewrite_report(dossier: pathlib.Path, change) -> None:
    """Replace ``report.json`` by ``change(report)`` (a new, self-consistent report) and rehash."""
    from GenAIRR._reference_sources.digests import report_hash
    from GenAIRR._reference_sources.report import BuildReport
    report = change(BuildReport.from_json((dossier / "report.json").read_bytes().decode("utf-8")))
    rewrite(dossier, "report.json", canonical(report.to_dict()),
            hashes={**report.hashes.to_dict(), "report": report_hash(report)}, runtime=report.runtime.to_dict())


def with_inventory(report, *, claim=lambda c: c, keep_mapping=lambda x: True):
    """The report over an edited inventory: ``claim(c)`` replaces each claim, mappings failing
    ``keep_mapping`` are dropped; claim outcomes follow their claims and every input hash is recomputed."""
    import dataclasses
    from GenAIRR._reference_sources import model
    from GenAIRR._reference_sources._fields import digest
    from GenAIRR._reference_sources.digests import BUILD_INPUTS_DOMAIN, unchecked_build_inputs
    from GenAIRR._reference_sources.report import BuildHashes
    inventory, renamed, claims = report.inventory, {}, []
    for old in inventory.claims:
        new = claim(old)
        renamed[old.id] = new.id
        claims.append(new)
    fields = {name: getattr(inventory, name) for name in ("source_lock_hash", "adapter", "records", "candidates",
                                                          "mentioned_identities", "observations", "support",
                                                          "findings")}
    edited = model.SourceInventory.assemble(claims=tuple(claims), **fields,
                                            mappings=tuple(x for x in inventory.mappings if keep_mapping(x)))
    outcomes = {renamed[o.claim]: dataclasses.replace(o, claim=renamed[o.claim]) for o in report.claims}
    inputs = unchecked_build_inputs(report.source_lock, edited, report.recipe)
    hashes = BuildHashes(inputs["source_lock"], inputs["inventory"], inputs["recipe"],
                         digest(inputs, BUILD_INPUTS_DOMAIN))
    return dataclasses.replace(report, inventory=edited, claims=tuple(outcomes[c.id] for c in edited.claims),
                               hashes=hashes)


def rebuilt_reference(dossier: pathlib.Path, tmp_path: pathlib.Path, edit=lambda files: files):
    """A runtime reference built from the dossier's inputs after ``edit({name: bytes})``."""
    from GenAIRR._reference.artifact import Reference
    files = edit({p.name: p.read_bytes() for p in (dossier / "inputs").iterdir()})
    folder = tmp_path / "edited-inputs"
    folder.mkdir(exist_ok=True)
    for name, data in files.items():
        (folder / name).write_bytes(data)
    return Reference.from_fasta(v=folder / "v.fasta", j=folder / "j.fasta", annotations=folder / "annotations.json",
                                d=folder / "d.fasta" if "d.fasta" in files else None)


def runtime_of(reference) -> dict:
    return {"catalogue_hash": reference.catalogue_hash, "models_hash": reference.models_hash,
            "execution_hash": reference.execution_hash, "artifact_hash": reference.artifact_hash}


def with_runtime(report, reference):
    import dataclasses
    from GenAIRR._reference_sources.common import RuntimeHashes
    return dataclasses.replace(report, runtime=RuntimeHashes(**runtime_of(reference)))


J7 = "IGLJ7*02"
V_RECORD, J_RECORD, OTHER = ("airr#/GermlineSet/0/allele_descriptions/1", "airr#/GermlineSet/0/allele_descriptions/2",
                             "airr#/GermlineSet/0/allele_descriptions/0")


@pytest.fixture
def mini(tmp_path, monkeypatch):
    """A synthetic IGL source over V, J1 and J7 (registered as validated for the test); V and J1 are
    selected, J7 is outside the selection.  Returns ``(build, dossier directory)``."""
    from tests.test_reference_source_api import mini_source
    from tests.test_reference_source_projection import entry, recipe
    source, blobs = mini_source(tmp_path, monkeypatch, labels=MINI + (J7,))
    build = source.build(recipe(source.lock, [entry(blobs, "airr", label) for label in MINI]))
    build.write_dossier(tmp_path / "dossier")
    return build, tmp_path / "dossier"


def rejected(dossier, kind, member) -> ReferenceBuild:
    exc = error(kind, ReferenceBuild.read_dossier, dossier)
    assert exc.path == "/" + member or exc.path.startswith("/" + member + "/"), exc
    return exc


def test_the_untampered_synthetic_dossier_reads(mini):
    build, dossier = mini
    assert [(r.record.id, r.outcome, r.exclusion) for r in build.report.records][:3] == [
        (OTHER, "excluded", "outside_explicit_selection"), (V_RECORD, "selected", None), (J_RECORD, "selected", None)]
    assert ReferenceBuild.read_dossier(dossier).manifest == build.manifest


def test_an_altered_mapping_is_rejected(mini):
    _, dossier = mini

    def unbind(report):
        dropped = {x.id for x in report.inventory.mappings
                   if x.kind == "gapped_alignment" and x.input.id == V_RECORD + "/coding_sequence"}
        assert len(dropped) == 1
        return with_inventory(report, keep_mapping=lambda x: x.id not in dropped, claim=lambda c: __import__(
            "dataclasses").replace(c, mappings=tuple(i for i in c.mappings if i not in dropped),
                                   evidence=tuple(e for e in c.evidence if e.asset != "gapped"))
            if dropped & set(c.mappings) else c)

    rewrite_report(dossier, unbind)
    exc = rejected(dossier, "digest_mismatch", "report.json")
    assert "regenerat" in exc.message


def test_an_altered_claim_support_is_rejected(mini):
    _, dossier = mini
    import dataclasses

    def weaker(report):
        return with_inventory(report, claim=lambda c: dataclasses.replace(c, support="assertion")
                              if (c.subject.id, c.predicate) == (J_RECORD + "/coding_sequence", "anchor") else c)

    rewrite_report(dossier, weaker)
    exc = rejected(dossier, "digest_mismatch", "report.json")
    assert "regenerat" in exc.message


def test_a_hidden_exclusion_is_rejected(mini):
    _, dossier = mini
    import dataclasses

    def hide(report):
        records = tuple(dataclasses.replace(r, exclusion="structural_incompatibility",
                                            detail="the record cannot be projected")
                        if r.record.id == OTHER else r for r in report.records)
        return dataclasses.replace(report, records=records)

    rewrite_report(dossier, hide)
    exc = rejected(dossier, "digest_mismatch", "report.json")
    assert exc.path == "/report.json/records/0/detail" or exc.path == "/report.json/records/0/exclusion"


def test_a_selected_label_changed_everywhere_but_the_source_is_rejected(mini, tmp_path):
    _, dossier = mini
    relabel = {"j.fasta": lambda d: d.replace(b">IGLJ1*01", b">IGLJ1*99"),
               "annotations.json": lambda d: d.replace(b'"label":"IGLJ1*01"', b'"label":"IGLJ1*99"')}

    def edit(files):
        return {name: relabel.get(name, lambda d: d)(data) for name, data in files.items()}

    reference = rebuilt_reference(dossier, tmp_path, edit)
    assert [a.label for a in reference.alleles] == ["IGLV1-36*01", "IGLJ1*99"]
    for name in ("j.fasta", "annotations.json"):
        rewrite(dossier, "inputs/" + name, edit({name: (dossier / "inputs" / name).read_bytes()})[name])
    rewrite(dossier, "reference.genairr-ref.json", reference.to_json().encode("utf-8"))
    rewrite_report(dossier, lambda report: with_runtime(report, reference))
    exc = rejected(dossier, "digest_mismatch", "inputs/j.fasta")             # the first regenerated input
    assert "regenerat" in exc.message


def test_altered_native_output_is_rejected(mini, tmp_path):
    _, dossier = mini

    def mutate(files):
        header, _, rest = files["v.fasta"].partition(b"\n")
        base = rest[:1]
        files["v.fasta"] = header + b"\n" + (b"A" if base != b"A" else b"C") + rest[1:]
        old = "sha256:" + hashlib.sha256(rest.replace(b"\n", b"")).hexdigest()
        new = "sha256:" + hashlib.sha256(files["v.fasta"].partition(b"\n")[2].replace(b"\n", b"")).hexdigest()
        files["annotations.json"] = files["annotations.json"].replace(old.encode(), new.encode())
        return files

    reference = rebuilt_reference(dossier, tmp_path, mutate)
    rewrite(dossier, "reference.genairr-ref.json", reference.to_json().encode("utf-8"))
    rewrite_report(dossier, lambda report: with_runtime(report, reference))
    exc = rejected(dossier, "digest_mismatch", "reference.genairr-ref.json")
    assert "regenerat" in exc.message


def test_a_report_claiming_another_runtime_identity_is_rejected(mini):
    _, dossier = mini
    import dataclasses
    from GenAIRR._reference_sources.common import RuntimeHashes

    def other(report):
        runtime = dict(report.runtime.to_dict(), artifact_hash="sha256:" + "ab" * 32)
        return dataclasses.replace(report, runtime=RuntimeHashes(**runtime))

    rewrite_report(dossier, other)
    exc = rejected(dossier, "digest_mismatch", "report.json")
    assert exc.path == "/report.json/runtime/artifact_hash"


def test_a_runtime_body_resealed_under_another_generator_is_rejected(mini, tmp_path, monkeypatch):
    _, dossier = mini
    with monkeypatch.context() as patch:
        patch.setattr(GenAIRR, "__version__", "0.0.1+forged")
        reference = rebuilt_reference(dossier, tmp_path)
    assert json.loads(reference.to_json())["generator"]["version"] == "0.0.1+forged"
    rewrite(dossier, "reference.genairr-ref.json", reference.to_json().encode("utf-8"))
    exc = rejected(dossier, "digest_mismatch", "report.json")             # the report names the stored artifact
    assert exc.path == "/report.json/runtime/artifact_hash"
    rewrite_report(dossier, lambda report: with_runtime(report, reference))
    exc = rejected(dossier, "identity_mismatch", "reference.genairr-ref.json")   # and its generator is not the builder
    assert "0.0.1+forged" in exc.message


def test_a_recipe_that_is_not_the_reports_recipe_is_rejected(mini):
    _, dossier = mini
    recipe = json.loads((dossier / "recipe.json").read_bytes())
    recipe["selection"][0]["reason"] = "a reason edited after the build"
    rewrite(dossier, "recipe.json", canonical(recipe))
    exc = rejected(dossier, "digest_mismatch", "recipe.json")
    assert "manifest.json commits to" in exc.message


# ── acquisition receipts never enter an identity ─────────────────────


def test_builds_with_different_receipts_have_identical_identities(tmp_path, monkeypatch):
    from GenAIRR._reference_sources import api
    builds = []
    for moment in (T0, T0 + datetime.timedelta(days=3, microseconds=7)):
        monkeypatch.setattr(api, "_now", lambda moment=moment: moment)
        source = from_files(tmp_path, cache=f"cache-{len(builds)}")
        builds.append(source.build(full_recipe(tmp_path, source)))
        builds[-1].write_dossier(tmp_path / f"dossier-{len(builds)}")
    first, second = builds
    assert first.manifest == second.manifest
    assert first.report.to_dict() == second.report.to_dict()
    assert runtime_of(first.reference) == runtime_of(second.reference)
    one, two = tree(tmp_path / "dossier-1"), tree(tmp_path / "dossier-2")
    receipts = [name for name in one if name.startswith("receipts/")]
    assert receipts and all(one[name] != two[name] for name in receipts)
    assert {k: v for k, v in one.items() if k not in receipts} == {k: v for k, v in two.items() if k not in receipts}
    assert b"2026-10-27T09:00:00" in one["receipts/0.json"] and b"2026-10-30T09:00:00" in two["receipts/0.json"]


# ── offline rebuilding and generator changes ─────────────────────────


@pytest.mark.usefixtures("offline")
def test_a_clean_cache_build_and_a_moved_offline_dossier_agree(tmp_path):
    """The public offline-rebuild composition: read the dossier, then rebuild from its own lock,
    cache and recipe after moving it; the deterministic files are byte-identical."""
    from GenAIRR.reference_sources import OgrdbSource
    source = from_files(tmp_path, cache="clean-cache")
    build = source.build(full_recipe(tmp_path, source))
    build.write_dossier(tmp_path / "igl-build")
    (tmp_path / "igl-build").rename(tmp_path / "moved")
    for path in sorted((tmp_path / "clean-cache").rglob("*"), reverse=True):
        path.unlink() if path.is_file() else path.rmdir()
    checked = ReferenceBuild.read_dossier(str(tmp_path / "moved"))
    assert checked.reference.execution_hash == build.reference.execution_hash
    dossier = pathlib.Path(str(tmp_path / "moved"))
    locked = OgrdbSource.read(dossier / "source.lock.json", cache=dossier / "cache")
    rebuilt = locked.build(dossier / "recipe.json")
    rebuilt.write_dossier(str(tmp_path / "igl-rebuilt"))
    original, again = tree(dossier), tree(tmp_path / "igl-rebuilt")
    listed = {member["name"] for member in load_manifest(dossier)["members"]} | {"manifest.json"}
    assert {k: original[k] for k in listed} == {k: again[k] for k in listed}
    assert all(b"verified_cache" in again[name] for name in again if name.startswith("receipts/"))


def test_an_older_generators_dossier_checks_and_rebuilds_under_the_installed_one(igl, tmp_path, monkeypatch):
    from GenAIRR._reference_sources.digests import manifest_hash, report_hash
    from GenAIRR.reference_sources import OgrdbSource
    old = "1.9.0"
    current = GenAIRR.__version__
    assert current != old
    source, _ = igl
    monkeypatch.setattr(GenAIRR, "__version__", old)
    build = source.build(full_recipe(tmp_path, source))
    build.write_dossier(tmp_path / "old")
    monkeypatch.undo()
    assert GenAIRR.__version__ == current
    stored = tree(tmp_path / "old")
    for name in ("reference.genairr-ref.json", "report.json", "manifest.json"):
        assert f'"version":"{old}"'.encode() in stored[name]

    checked = ReferenceBuild.read_dossier(tmp_path / "old")
    assert checked.checked_with == current
    assert json.loads(checked.reference.to_json())["generator"]["version"] == old       # the stored artifact
    assert (checked.report.generator.version, checked.manifest.builder.version) == (old, old)
    assert tree(tmp_path / "old") == stored

    dossier = tmp_path / "old"
    rebuilt = OgrdbSource.read(dossier / "source.lock.json", cache=dossier / "cache").build(dossier / "recipe.json")
    rebuilt.write_dossier(tmp_path / "rebuilt")
    assert (rebuilt.report.generator.version, rebuilt.manifest.builder.version) == (current, current)
    assert json.loads(rebuilt.reference.to_json())["generator"]["version"] == current
    assert rebuilt.reference.execution_hash == checked.reference.execution_hash
    assert (rebuilt.reference.catalogue_hash, rebuilt.reference.models_hash) == (
        checked.reference.catalogue_hash, checked.reference.models_hash)
    assert rebuilt.reference.artifact_hash != checked.reference.artifact_hash
    assert report_hash(rebuilt.report) != report_hash(checked.report)
    assert manifest_hash(rebuilt.manifest) != manifest_hash(checked.manifest)
    for name, data in tree(tmp_path / "rebuilt").items():
        assert f'"version":"{old}"'.encode() not in data, name
    assert ReferenceBuild.read_dossier(tmp_path / "rebuilt").manifest == rebuilt.manifest


# ── a moved dossier in a fresh offline process ───────────────────────

# Runs in a fresh isolated interpreter: every network attempt is counted and refused before
# GenAIRR is imported, and every warning is an error.
OFFLINE_SCRIPT = r'''
import json
import socket
import sys
import warnings
from pathlib import Path

ATTEMPTS = []


def _refuse(*args, **kwargs):
    ATTEMPTS.append(repr(args)[:200])          # counted even if a caller swallows the error
    raise OSError("network access is disabled in this offline check")


class _NoSocket(socket.socket):
    def __init__(self, *args, **kwargs):
        _refuse(*args)


socket.socket = _NoSocket
socket.create_connection = _refuse
socket.getaddrinfo = _refuse
warnings.simplefilter("error")

from GenAIRR import Experiment
from GenAIRR.reference_sources import OgrdbSource, ReferenceBuild

dossier, out = Path(sys.argv[1]), Path(sys.argv[2])
checked = ReferenceBuild.read_dossier(dossier)
locked = OgrdbSource.read(dossier / "source.lock.json", cache=dossier / "cache")
rebuilt = locked.build(dossier / "recipe.json")
rebuilt.write_dossier(out)
reference = ReferenceBuild.read_dossier(out).reference
result = Experiment.on(reference).recombine().run_records(n=25, seed=7, validate_records=True)
print(json.dumps({
    "checked_execution_hash": checked.reference.execution_hash,
    "rebuilt_execution_hash": rebuilt.reference.execution_hash,
    "reference_content_hash": result.reference_content_hash,
    "execution_hash": reference.execution_hash,
    "records": len(result),
    "alleles": sorted({record["v_call"] for record in result} | {record["j_call"] for record in result}),
    "network_attempts": ATTEMPTS,
}))
'''


def test_a_moved_dossier_checks_rebuilds_and_simulates_in_a_fresh_offline_process(igl, tmp_path):
    import shutil
    import subprocess
    import sys
    _, build = igl
    build.write_dossier(tmp_path / "written")
    moved = tmp_path / "elsewhere" / "igl-build"
    shutil.copytree(tmp_path / "written", moved)
    shutil.rmtree(tmp_path / "written")
    script = tmp_path / "offline_check.py"
    script.write_bytes(OFFLINE_SCRIPT.encode("utf-8"))
    process = subprocess.run([sys.executable, "-I", str(script), str(moved), str(tmp_path / "rebuilt")],
                             capture_output=True, timeout=600, cwd=tmp_path)
    assert process.returncode == 0, process.stderr.decode("utf-8", "replace")
    outcome = json.loads(process.stdout.decode("utf-8"))
    assert outcome["network_attempts"] == []
    assert outcome["checked_execution_hash"] == outcome["rebuilt_execution_hash"] == build.reference.execution_hash
    assert outcome["reference_content_hash"] == outcome["execution_hash"] == build.reference.execution_hash
    assert outcome["records"] == 25
    assert set(outcome["alleles"]) <= set(MINI)
    original, rebuilt = tree(moved), tree(tmp_path / "rebuilt")
    listed = {member["name"] for member in load_manifest(moved)["members"]} | {"manifest.json"}
    assert {k: original[k] for k in listed} == {k: rebuilt[k] for k in listed}


# ── packaged dossiers ────────────────────────────────────────────────

STARTER_ID = "human_igh_ogrdb_10_observed_v1"


def test_the_bundled_registry_is_closed_and_exact():
    from GenAIRR._reference_sources import bundled
    assert dict(bundled.REGISTRY) == {}
    with pytest.raises(TypeError):
        bundled.REGISTRY["x"] = None
    exc = error("missing_asset", ReferenceBuild.bundled, STARTER_ID)
    assert exc.path == "" and STARTER_ID in exc.message and "available: none" in exc.message
    with pytest.raises(TypeError):
        ReferenceBuild.bundled(None)


def packaged(tmp_path, monkeypatch, dossier: pathlib.Path, *, name="genairr_test_pkg", zipped=False, extra=()):
    """``dossier`` packaged as resource ``builds/starter`` of a throwaway importable package; ``extra``
    adds ``(entry name below the package, bytes)`` entries to a zipped package."""
    import shutil
    import sys
    import warnings
    import zipfile
    root = tmp_path / "site"
    root.mkdir()
    if zipped:
        archive = root / "packages.zip"
        with zipfile.ZipFile(archive, "w") as bundle, warnings.catch_warnings():
            warnings.simplefilter("ignore")            # duplicate entry names are written on purpose
            bundle.writestr(f"{name}/__init__.py", "")
            for path in sorted(dossier.rglob("*")):
                if path.is_file():
                    bundle.write(path, f"{name}/builds/starter/{path.relative_to(dossier).as_posix()}")
            for entry_name, data in extra:
                bundle.writestr(f"{name}/{entry_name}", data)
        monkeypatch.syspath_prepend(str(archive))
    else:
        (root / name).mkdir()
        (root / name / "__init__.py").write_bytes(b"")
        shutil.copytree(dossier, root / name / "builds" / "starter")
        monkeypatch.syspath_prepend(str(root))
    monkeypatch.delitem(sys.modules, name, raising=False)
    return root


def register(monkeypatch, entry):
    """Register ``entry`` as the only packaged dossier, with an empty cache of checked builds."""
    from GenAIRR._reference_sources import bundled
    monkeypatch.setattr(bundled, "REGISTRY", __import__("types").MappingProxyType({STARTER_ID: entry}))
    monkeypatch.setattr(build_module, "_BUNDLED", {})


def test_a_registered_dossier_loads_through_the_dossier_reader(mini, tmp_path, monkeypatch):
    from GenAIRR._reference_sources import bundled
    from GenAIRR._reference_sources.digests import manifest_hash
    build, dossier = mini
    root = packaged(tmp_path, monkeypatch, dossier)
    register(monkeypatch, bundled.BundledDossier("genairr_test_pkg", "builds/starter", manifest_hash(build.manifest)))
    before = tree(root)
    reads = []
    real = ReferenceBuild.read_dossier.__func__
    monkeypatch.setattr(ReferenceBuild, "read_dossier",
                        classmethod(lambda cls, directory: reads.append(directory) or real(cls, directory)))
    if hasattr(__import__("os"), "geteuid"):
        for path in sorted(root.rglob("*"), reverse=True):
            path.chmod(0o555 if path.is_dir() else 0o444)
    try:
        loaded = ReferenceBuild.bundled(STARTER_ID)
        assert ReferenceBuild.bundled(STARTER_ID) is loaded               # checked once, never reopened
    finally:
        for path in sorted(root.rglob("*")):
            path.chmod(0o755 if path.is_dir() else 0o644)
    assert len(reads) == 1 and pathlib.Path(reads[0]) == root / "genairr_test_pkg" / "builds" / "starter"
    assert loaded.manifest == build.manifest and loaded.reference.to_json() == build.reference.to_json()
    assert tree(root) == before
    exc = error("missing_asset", ReferenceBuild.bundled, STARTER_ID.upper())
    assert STARTER_ID in exc.message


def test_a_dossier_packaged_in_an_archive_is_read_from_a_private_copy(mini, tmp_path, monkeypatch):
    from GenAIRR._reference_sources import bundled
    from GenAIRR._reference_sources.digests import manifest_hash
    build, dossier = mini
    packaged(tmp_path, monkeypatch, dossier, name="genairr_zipped_pkg", zipped=True)
    register(monkeypatch, bundled.BundledDossier("genairr_zipped_pkg", "builds/starter", manifest_hash(build.manifest)))
    copies = []
    real = ReferenceBuild.read_dossier.__func__
    monkeypatch.setattr(ReferenceBuild, "read_dossier",
                        classmethod(lambda cls, directory: copies.append(pathlib.Path(directory)) or real(cls, directory)))
    loaded = ReferenceBuild.bundled(STARTER_ID)
    assert loaded.manifest == build.manifest
    (copy,) = copies
    assert not copy.exists()                                     # the private copy is removed afterwards
    assert tmp_path not in copy.parents


def test_a_packaged_dossier_must_have_its_registered_manifest_before_it_is_checked(mini, tmp_path, monkeypatch):
    from GenAIRR._reference_sources import bundled
    _, dossier = mini
    packaged(tmp_path, monkeypatch, dossier, name="genairr_other_pkg")
    register(monkeypatch, bundled.BundledDossier("genairr_other_pkg", "builds/starter", "sha256:" + "0" * 64))
    checks = []
    real = ReferenceBuild.read_dossier.__func__
    monkeypatch.setattr(ReferenceBuild, "read_dossier",
                        classmethod(lambda cls, directory: checks.append(directory) or real(cls, directory)))
    exc = error("digest_mismatch", ReferenceBuild.bundled, STARTER_ID)
    assert "sha256:" + "0" * 64 in exc.message and exc.path == "/manifest.json"
    assert checks == []                                          # refused before the full check


def test_concurrent_first_loads_check_a_packaged_dossier_once(mini, tmp_path, monkeypatch):
    import threading
    import time
    from GenAIRR._reference_sources import bundled
    from GenAIRR._reference_sources.digests import manifest_hash
    build, dossier = mini
    packaged(tmp_path, monkeypatch, dossier, name="genairr_threaded_pkg")
    register(monkeypatch, bundled.BundledDossier("genairr_threaded_pkg", "builds/starter",
                                                 manifest_hash(build.manifest)))
    checks = []
    real = ReferenceBuild.read_dossier.__func__

    def slow(cls, directory):
        checks.append(directory)
        time.sleep(0.3)
        return real(cls, directory)

    monkeypatch.setattr(ReferenceBuild, "read_dossier", classmethod(slow))
    start = threading.Barrier(3)
    results = []

    def load():
        start.wait()
        results.append(ReferenceBuild.bundled(STARTER_ID))

    threads = [threading.Thread(target=load) for _ in range(3)]
    for thread in threads:
        thread.start()
    for thread in threads:
        thread.join(timeout=60)
    assert len(results) == 3 and len(checks) == 1
    assert results[0] is results[1] is results[2]


@pytest.mark.parametrize("resource", ["/abs", "../up", "a\\b", "", "C:/x", "a/../b"])
def test_a_registry_entry_names_a_resource_below_its_package(resource):
    from GenAIRR._reference_sources import bundled
    with pytest.raises(Exception) as info:
        bundled.BundledDossier("GenAIRR", resource, "sha256:" + "0" * 64)
    assert getattr(info.value, "kind", None) in ("malformed", "resource_limit")


# ── manifests that cannot be trusted ────────────────────────────────


@pytest.mark.parametrize("name", ["/etc/passwd", "../outside.json", "inputs\\v.fasta", "C:/report.json",
                                  "inputs/../report.json", "receipts/0.json", "manifest.json"])
def test_a_manifest_naming_an_unsafe_member_is_refused(mini, name):
    _, dossier = mini
    manifest = load_manifest(dossier)
    manifest["members"].insert(0, {"name": name, "sha256": "sha256:" + "0" * 64, "bytes": 1})
    save_manifest(dossier, manifest)
    exc = error("malformed", ReferenceBuild.read_dossier, dossier)
    assert exc.path == "/manifest.json/members/0/name"


@pytest.mark.parametrize("field, value, kind, path", [
    ("versions", {"report": 2}, "unsupported_schema", "/manifest.json/versions/report"),
    ("projection", {"version": 2}, "unsupported_profile", "/manifest.json/projection"),
    ("adapter", {"version": 2}, "unsupported_profile", "/manifest.json/adapter"),
    ("builder", {"version": "0.0.1"}, "digest_mismatch", "/report.json/generator/version"),
    ("runtime", {"execution_hash": "sha256:" + "1" * 64}, "digest_mismatch", "/manifest.json/runtime"),
    ("hashes", {"report": "sha256:" + "2" * 64}, "digest_mismatch", "/manifest.json/hashes/report"),
])
def test_a_manifest_that_does_not_commit_to_its_members_is_refused(mini, field, value, kind, path):
    _, dossier = mini
    manifest = load_manifest(dossier)
    manifest[field].update(value)
    save_manifest(dossier, manifest)
    exc = error(kind, ReferenceBuild.read_dossier, dossier)
    assert exc.path == path


def test_a_blob_absent_from_the_dossier_and_its_manifest_is_a_missing_asset(mini):
    build, dossier = mini
    gapped = build.lock.asset("gapped")
    name = "cache/blobs/sha256/" + gapped.blob_name
    (dossier / name).unlink()
    manifest = load_manifest(dossier)
    manifest["members"] = [m for m in manifest["members"] if m["name"] != name]
    save_manifest(dossier, manifest)
    exc = error("missing_asset", ReferenceBuild.read_dossier, dossier)
    assert exc.path == "/manifest.json/members" and gapped.sha256 in exc.message


def test_a_dossier_outside_this_versions_validated_scope_is_unsupported(mini, monkeypatch):
    from GenAIRR._reference_sources import common, consistency
    _, dossier = mini
    monkeypatch.setattr(consistency, "PROFILE_SCOPES", common.PROFILE_SCOPES)
    exc = error("unsupported_profile", ReferenceBuild.read_dossier, dossier)
    assert exc.path == "/source.lock.json" and "not been validated" in exc.message


@pytest.mark.parametrize("tag, link", [
    (0xA000000C, True),                 # IO_REPARSE_TAG_SYMLINK
    (0xA0000003, True),                 # IO_REPARSE_TAG_MOUNT_POINT (directory junction)
    (0xA000001D, True),                 # IO_REPARSE_TAG_LX_SYMLINK
    (0x9000601A, False),                # a cloud-file placeholder: an ordinary file, read without following
    (0x80000013, False),                # IO_REPARSE_TAG_DEDUP
    (None, True),                       # a reparse point whose tag is unavailable: refused
])
def test_only_windows_name_surrogate_reparse_points_count_as_links(tag, link):
    import stat
    import types
    attributes = {"st_file_attributes": stat.FILE_ATTRIBUTE_REPARSE_POINT, "st_mode": stat.S_IFREG | 0o644}
    if tag is not None:
        attributes["st_reparse_tag"] = tag
    entry = types.SimpleNamespace(is_symlink=lambda: False,
                                  stat=lambda follow_symlinks=True: types.SimpleNamespace(**attributes))
    plain = types.SimpleNamespace(is_symlink=lambda: False,
                                  stat=lambda follow_symlinks=True: types.SimpleNamespace(st_mode=stat.S_IFREG))
    assert build_module._is_link(entry) is link
    assert build_module._is_link(plain) is False



# ── the cost of a check ──────────────────────────────────────────────


def count_report_work(monkeypatch) -> dict:
    """Count stored-report parses, report-hash computations and report serialisations."""
    from GenAIRR._reference_sources import _fields, digests
    from GenAIRR._reference_sources.report import ARTIFACT, BuildReport
    counts = {"parsed": 0, "hashed": 0, "serialised": 0}
    real_from_dict = BuildReport.from_dict.__func__

    def from_dict(cls, d, path=""):
        counts["parsed"] += 1
        return real_from_dict(cls, d, path)

    real_digest = digests.digest

    def digest(value, domain):
        counts["hashed"] += domain == digests.REPORT_DOMAIN
        return real_digest(value, domain)

    real_bytes = build_module.document_bytes

    def document_bytes(document):
        counts["serialised"] += document.get("artifact") == ARTIFACT
        return real_bytes(document)

    monkeypatch.setattr(BuildReport, "from_dict", classmethod(from_dict))
    monkeypatch.setattr(digests, "digest", digest)
    monkeypatch.setattr(build_module, "document_bytes", document_bytes)
    assert _fields.digest is real_digest
    return counts


def test_a_successful_check_compares_report_bytes_and_hashes_the_report_once(igl, tmp_path, monkeypatch):
    _, build = igl
    build.write_dossier(tmp_path / "dossier")
    counts = count_report_work(monkeypatch)
    checked = ReferenceBuild.read_dossier(tmp_path / "dossier")
    assert counts == {"parsed": 0, "hashed": 1, "serialised": 1}
    assert checked.report.to_dict() == build.report.to_dict() and checked.manifest == build.manifest


def test_a_build_hashes_its_report_once(igl, tmp_path, monkeypatch):
    source, _ = igl
    document = full_recipe(tmp_path, source)
    counts = count_report_work(monkeypatch)
    source.build(document)
    assert (counts["hashed"], counts["serialised"]) == (1, 1)


# ── platform metadata and leftovers ──────────────────────────────────

PLATFORM_FILES = (".DS_Store", "inputs/Thumbs.db", "cache/blobs/desktop.ini", "receipts/._0.json", "._manifest.json",
                  "cache/blobs/sha256/.DS_Store")


def test_platform_metadata_files_are_ignored_and_never_read(mini, monkeypatch):
    build, dossier = mini
    for name in PLATFORM_FILES:
        (dossier / name).write_bytes(b"\x00 written by a file manager \xff")
    opened = []
    real_open = build_module.os.open
    monkeypatch.setattr(build_module.os, "open", lambda path, *a, **k: opened.append(pathlib.Path(path).name)
                        or real_open(path, *a, **k))
    checked = ReferenceBuild.read_dossier(dossier)
    assert checked.manifest == build.manifest
    assert not {pathlib.Path(name).name for name in PLATFORM_FILES} & set(opened)
    assert all(pathlib.Path(name).name in opened for name in ("manifest.json", "report.json", "v.fasta"))


@pytest.mark.parametrize("name", [".DS_Store2", "inputs/thumbs.db", "Desktop.ini", "_.foo", "receipts/desktop.ini.bak"])
def test_names_near_platform_metadata_are_still_refused(mini, name):
    _, dossier = mini
    (dossier / name).write_bytes(b"x")
    exc = error("malformed", ReferenceBuild.read_dossier, dossier)
    assert exc.path == "/" + name


def test_platform_metadata_names_must_be_regular_files(mini, tmp_path):
    _, dossier = mini
    (dossier / "inputs" / ".DS_Store").mkdir()
    exc = error("malformed", ReferenceBuild.read_dossier, dossier)
    assert exc.path == "/inputs/.DS_Store"
    (dossier / "inputs" / ".DS_Store").rmdir()
    (tmp_path / "elsewhere").write_bytes(b"x")
    symlink_or_skip(tmp_path / "elsewhere", dossier / "._report.json")
    exc = error("io", ReferenceBuild.read_dossier, dossier)
    assert exc.path == "/._report.json"


def test_a_leftover_temporary_file_is_named_as_such(mini):
    _, dossier = mini
    (dossier / ".manifest.json.0123456789abcdef.tmp").write_bytes(b"{}")
    exc = error("malformed", ReferenceBuild.read_dossier, dossier)
    assert exc.path == "/.manifest.json.0123456789abcdef.tmp"
    assert "leftover temporary file from an interrupted write" in exc.message
    assert "removed once no writer is active" in exc.message
    (dossier / ".manifest.json.0123456789abcdef.tmp").unlink()
    (dossier / ".manifest.json.not-a-temporary.tmp").write_bytes(b"{}")
    exc = error("malformed", ReferenceBuild.read_dossier, dossier)
    assert "leftover temporary" not in exc.message


@pytest.mark.skipif(not hasattr(__import__("os"), "mkfifo"), reason="FIFOs are unavailable here")
def test_a_fifo_swapped_in_after_the_scan_never_blocks_the_reader(mini, monkeypatch):
    import os
    import threading
    from GenAIRR.reference_sources import ReferenceSourceError
    _, dossier = mini
    member = dossier / "inputs" / "v.fasta"
    real_scan = build_module._scan

    def swapping(root):
        found = real_scan(root)
        member.unlink()
        os.mkfifo(member)
        return found

    monkeypatch.setattr(build_module, "_scan", swapping)
    outcome = []

    def check():
        try:
            ReferenceBuild.read_dossier(dossier)
        except ReferenceSourceError as exc:
            outcome.append(exc)

    worker = threading.Thread(target=check, daemon=True)
    worker.start()
    worker.join(timeout=10)
    if worker.is_alive():
        with open(member, "wb"):              # releases the blocked reader
            pass
        worker.join(timeout=10)
        pytest.fail("opening a dossier member blocked on a FIFO")
    (exc,) = outcome
    assert (exc.kind, exc.path) == ("io", "/inputs/v.fasta") and "regular file" in exc.message


def test_a_packaged_resource_path_is_joined_one_part_at_a_time(mini, tmp_path, monkeypatch):
    """Resource directories join one name per call, as importlib resources do on Python 3.9."""
    import importlib.resources
    from GenAIRR._reference_sources import bundled
    from GenAIRR._reference_sources.digests import manifest_hash
    build, dossier = mini
    packaged(tmp_path, monkeypatch, dossier, name="genairr_joined_pkg", zipped=True)
    register(monkeypatch, bundled.BundledDossier("genairr_joined_pkg", "builds/starter", manifest_hash(build.manifest)))
    real_files = importlib.resources.files

    class OneName:
        def __init__(self, inner):
            self._inner = inner

        def joinpath(self, child):
            return OneName(self._inner.joinpath(child))

        def iterdir(self):
            return (OneName(x) for x in self._inner.iterdir())

        def __getattr__(self, name):
            return getattr(self._inner, name)

    monkeypatch.setattr(bundled.importlib.resources, "files", lambda package: OneName(real_files(package)))
    assert ReferenceBuild.bundled(STARTER_ID).manifest == build.manifest


def test_a_packaged_copy_never_reads_beyond_its_byte_budget(mini, tmp_path, monkeypatch):
    import importlib.resources
    from GenAIRR._reference_sources import bundled
    from GenAIRR._reference_sources.digests import manifest_hash
    build, dossier = mini
    packaged(tmp_path, monkeypatch, dossier, name="genairr_budget_pkg", zipped=True)
    register(monkeypatch, bundled.BundledDossier("genairr_budget_pkg", "builds/starter", manifest_hash(build.manifest)))
    monkeypatch.setattr(bundled, "_MAX_BYTES", 4096)
    real_files = importlib.resources.files
    delivered = []

    class Counted:
        def __init__(self, inner):
            self._inner = inner

        def joinpath(self, child):
            return Counted(self._inner.joinpath(child))

        def iterdir(self):
            return (Counted(x) for x in self._inner.iterdir())

        def read_bytes(self):
            data = self._inner.read_bytes()
            delivered.append(len(data))
            return data

        def open(self, mode="r", *args, **kwargs):
            handle = self._inner.open(mode, *args, **kwargs)
            real_read = handle.read

            def read(size=-1):
                data = real_read(size)
                delivered.append(len(data))
                return data

            handle.read = read
            return handle

        def __getattr__(self, name):
            return getattr(self._inner, name)

    monkeypatch.setattr(bundled.importlib.resources, "files", lambda package: Counted(real_files(package)))
    error("resource_limit", ReferenceBuild.bundled, STARTER_ID)
    assert delivered and max(delivered) <= 4096 + 1


@pytest.mark.parametrize("extra, fragment", [
    ((("builds/starter/../escaped.txt", b"x"),), "'..'"),
    ((("builds/starter/inputs/v.fasta", b">again\nACGT\n"),), "twice"),
])
def test_a_packaged_entry_with_an_unsafe_or_repeated_name_is_a_typed_error(mini, tmp_path, monkeypatch, extra,
                                                                             fragment):
    from GenAIRR._reference_sources import bundled
    from GenAIRR._reference_sources.digests import manifest_hash
    build, dossier = mini
    packaged(tmp_path, monkeypatch, dossier, name="genairr_unsafe_pkg", zipped=True, extra=extra)
    register(monkeypatch, bundled.BundledDossier("genairr_unsafe_pkg", "builds/starter", manifest_hash(build.manifest)))
    exc = error("malformed", ReferenceBuild.bundled, STARTER_ID)
    assert fragment in exc.message
    assert not (tmp_path / "escaped.txt").exists()


@pytest.mark.parametrize("zipped", [False, True])
def test_an_absent_packaged_dossier_is_a_missing_asset(mini, tmp_path, monkeypatch, zipped):
    from GenAIRR._reference_sources import bundled
    build, dossier = mini
    name = "genairr_absent_zip" if zipped else "genairr_absent_pkg"
    packaged(tmp_path, monkeypatch, dossier, name=name, zipped=zipped)
    register(monkeypatch, bundled.BundledDossier(name, "builds/elsewhere", "sha256:" + "0" * 64))
    exc = error("missing_asset", ReferenceBuild.bundled, STARTER_ID)
    assert "builds/elsewhere" in exc.message
