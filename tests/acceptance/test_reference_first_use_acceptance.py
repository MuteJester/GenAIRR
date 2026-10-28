"""First-use acceptance for portable local references.

A small human IGH reference (V, D and J FASTA plus an annotation sidecar) is built once,
written as a ``*.genairr-ref.json`` artifact, moved, and run in a fresh offline process; a
gapped variant derives the same anchors from a declared IMGT alignment; a VJ light-chain
reference with local allele names exercises record accounting, missing-evidence
explanations and relationship-only names; an unlisted species runs without any registry.

The fixtures under ``tests/data/reference_first_use`` and the golden artifact are written
by ``tests/data/reference_first_use/make_fixtures.py``.  They are software fixtures copied
from GenAIRR's bundled cartridges: every annotation is an assertion, never source evidence.
"""
import json
import os
import shutil
import subprocess
import sys
from pathlib import Path, PurePosixPath, PureWindowsPath

import pytest

from GenAIRR import Experiment, Reference, ReferenceCapabilityError, Species
from GenAIRR._artifact_io import canonical_json, typed_digest
from tests.support import paths

TESTS = paths.tests_root()
FIXTURES = TESTS / "data" / "reference_first_use"
IGH = FIXTURES / "human_igh"
IGH_GAPPED = FIXTURES / "human_igh_gapped"
IGK = FIXTURES / "igk_local"
GOLDEN = TESTS / "golden" / "reference" / "human_igh.genairr-ref.json"
ARTIFACT_DOMAIN = "genairr_reference_artifact_v1"

# Pasted from the golden artifact (see the module docstring for how it is regenerated).
IGH_CATALOGUE_HASH = "sha256:08362156352554255239a4882a6738e430ac1c26dd38cbb8ce9c6238a1d1873c"
IGH_EXECUTION_HASH = "sha256:843f9736d07894e36a9da363ebf1ccc9c4390741f4794870f086d9e4f07eaad9"

IGH_V = ("IGHV1-2*02", "IGHV3-23*01", "IGHV4-34*01")
IGH_D = ("IGHD3-10*01", "IGHD2-2*01")
IGH_J = ("IGHJ4*02", "IGHJ6*02")
IGH_METHOD = "software fixture copied from GenAIRR HUMAN_IGH_IMGT; not independent evidence"
SUBREGIONS = ["FWR1", "CDR1", "FWR2", "CDR2", "FWR3"]

# Runs in a fresh interpreter: no network (every attempt is counted), every warning is an error.
SCRIPT = r'''
import json
import socket
import sys
import warnings


NETWORK_CALLS = []


def _offline(*args, **kwargs):
    NETWORK_CALLS.append(args)   # counted even if a caller swallows the error below
    raise OSError("network access is disabled in this acceptance test")


socket.socket = _offline
socket.create_connection = _offline
warnings.simplefilter("error")

from GenAIRR import Experiment, Reference

ref = Reference.read(sys.argv[1])
result = Experiment.on(ref).recombine().run_records(n=200, seed=7, validate_records=True)
print(json.dumps({
    "catalogue_hash": ref.catalogue_hash,
    "execution_hash": ref.execution_hash,
    "reference_content_hash": result.reference_content_hash,
    "n": len(result),
    "productive_unknown": sum(1 for record in result if record["productive"] is None),
    "network_calls": len(NETWORK_CALLS),
}))
'''


def build_igh(directory=IGH):
    return Reference.from_fasta(v=directory / "v.fasta", d=directory / "d.fasta", j=directory / "j.fasta",
                                annotations=directory / "annotations.json")


def build_igh_gapped():
    return Reference.from_fasta(v=IGH_GAPPED / "v.fasta", d=IGH / "d.fasta", j=IGH / "j.fasta",
                                annotations=IGH_GAPPED / "annotations.json")


def build_igk(directory=IGK):
    return Reference.from_fasta(v=directory / "v.fasta", j=directory / "j.fasta",
                                annotations=directory / "annotations.json")


@pytest.fixture(scope="module")
def igh():
    return build_igh()


@pytest.fixture(scope="module")
def igh_gapped():
    return build_igh_gapped()


@pytest.fixture(scope="module")
def igk():
    return build_igk()


def _segment(ref, segment):
    return [a for a in ref.alleles if a.segment == segment]


def _names(call):
    return [name for name in call.split(",") if name]


def _strings(value):
    """Every string in a JSON document, keys included."""
    if isinstance(value, str):
        yield value
    elif isinstance(value, dict):
        for key, item in value.items():
            yield key
            yield from _strings(item)
    elif isinstance(value, list):
        for item in value:
            yield from _strings(item)


# ── Human IGH: built once, moved, run offline in a fresh process ─────


def test_igh_first_use_offline_fresh_process(igh, tmp_path):
    built = tmp_path / "build" / "human_igh.genairr-ref.json"
    built.parent.mkdir()
    igh.write(built)
    moved = tmp_path / "elsewhere" / built.name
    moved.parent.mkdir()
    shutil.move(str(built), str(moved))
    assert not built.exists()

    try:
        proc = subprocess.run([sys.executable, "-c", SCRIPT, str(moved)], capture_output=True, text=True,
                              encoding="utf-8", check=True, cwd=moved.parent)
    except subprocess.CalledProcessError as exc:
        pytest.fail(f"fresh offline process failed (exit {exc.returncode}):\n{exc.stderr}")
    assert json.loads(proc.stdout) == {
        "catalogue_hash": igh.catalogue_hash,
        "execution_hash": igh.execution_hash,
        "reference_content_hash": igh.execution_hash,
        "n": 200,
        "productive_unknown": 0,
        "network_calls": 0,
    }


def test_igh_identity_and_annotations(igh):
    assert (igh.name, igh.chain_type, igh.locus) == ("human_igh_fixture", "vdj", "IGH")
    assert (igh.taxon.scientific_name, igh.taxon.taxonomy_id) == ("Homo sapiens", 9606)
    assert [a.label for a in igh.alleles] == [*IGH_V, *IGH_D, *IGH_J]
    assert [a.anchor.value for a in _segment(igh, "V")] == [285, 285, 282]
    assert [a.anchor.value for a in _segment(igh, "J")] == [14, 29]
    assert [a.coding_phase.value for a in _segment(igh, "V")] == [0, 0, 0]
    assert [a.coding_phase.value for a in _segment(igh, "J")] == [2, 2]
    expected_boundaries = {"V": ("unknown", "complete"), "D": ("complete", "complete"),
                           "J": ("complete", "unknown")}
    for a in igh.alleles:
        assert a.functionality.value == "functional"
        boundaries = tuple("unknown" if b is None else b.value for b in (a.boundary_5, a.boundary_3))
        assert boundaries == expected_boundaries[a.segment], a.label
        for name in ("anchor", "functionality", "coding_phase", "boundary_5", "boundary_3"):
            summary = getattr(a, name)
            if summary is None:
                continue
            # Software fixture: every annotation is the author's assertion, never source evidence.
            assert (summary.support, summary.source) == ("assertion", None), (a.label, name)
            expected_derivation = "anchor_derived" if name == "coding_phase" else "stated"
            assert summary.derivation == expected_derivation, (a.label, name)
    side = json.loads((IGH / "annotations.json").read_text(encoding="utf-8"))
    methods = {entry[name]["method"] for entry in side["alleles"]
               for name in ("anchor", "functionality", "coding_phase", "boundary_5", "boundary_3")
               if entry[name] is not None}
    assert methods == {IGH_METHOD}
    assert all(r.status == "accepted" for r in igh.records)
    assert igh.entities == () and igh.findings == ()


def test_igh_capabilities_and_evidence_levels(igh):
    caps = igh.capabilities()
    assert caps.available() == ("catalogue_inspection", "observed_span_simulation", "canonical_productive",
                                "original_coding_ends")
    for name in caps.available():
        status = caps.status(name)
        if name in ("catalogue_inspection", "observed_span_simulation"):
            assert status.evidence is None, name
        else:
            assert status.evidence == "assertion", name
    assert caps.available(evidence="source") == ("catalogue_inspection", "observed_span_simulation")
    assert caps.phase_support == "anchor_derived"

    span = caps.status("complete_reference_span")
    assert not span.available and span.evidence is None
    assert [(u.allele, u.requirement) for u in span.unmet] == (
        [(f"V:{label}", 'accepted boundary_5 == "complete"') for label in IGH_V]
        + [(f"J:{label}", 'accepted boundary_3 == "complete"') for label in IGH_J])

    # Every anchor is an assertion, so the junction is not source-annotated.
    junction = caps.status("source_annotated_junction")
    assert not junction.available and junction.evidence is None
    assert [(u.allele, u.requirement) for u in junction.unmet] == (
        [(f"V:{label}", "anchor") for label in IGH_V] + [(f"J:{label}", "anchor") for label in IGH_J])

    igh.require("observed_span_simulation", "canonical_productive", "original_coding_ends")
    igh.require("observed_span_simulation", evidence="source")
    with pytest.raises(ReferenceCapabilityError):
        igh.require("canonical_productive", evidence="source")
    with pytest.raises(ReferenceCapabilityError):
        igh.require("source_annotated_junction")
    with pytest.raises(ReferenceCapabilityError):
        igh.require("complete_reference_span")


# ── Declared IMGT alignment: the same anchors, derived rather than stated ──


def test_declared_alignment_fixture_derives_same_anchors(igh_gapped, igh):
    side = json.loads((IGH_GAPPED / "annotations.json").read_text(encoding="utf-8"))
    assert side["conventions"]["v_alignment"] == "imgt_gapped"
    assert all(entry["anchor"] is None for entry in side["alleles"] if entry["segment"] == "V")

    v = _segment(igh_gapped, "V")
    assert [a.label for a in v] == list(IGH_V)
    assert [a.anchor.value for a in v] == [285, 285, 282]
    assert {(a.anchor.support, a.anchor.derivation, a.anchor.source) for a in v} == {
        ("assertion", "imgt_alignment", None)}
    assert [a.sequence_sha256 for a in v] == [a.sequence_sha256 for a in _segment(igh, "V")]

    doc = json.loads(igh_gapped.to_json())
    subregion_claims = [c for c in doc["catalogue"]["claims"] if c["predicate"] == "v_subregions"]
    assert [c["subject"] for c in subregion_claims] == [f"V:{label}" for label in IGH_V]
    assert {(c["support"], c["derivation"], c["status"]) for c in subregion_claims} == {
        ("assertion", "imgt_alignment", "accepted")}
    native = igh_gapped.native_refdata()
    for index in range(native.v_pool_size()):
        assert [label for label, _, _ in native.v_allele(index).subregions] == SUBREGIONS
    assert all(not igh.native_refdata().v_allele(i).subregions for i in range(3))

    assert igh_gapped.capabilities().available() == igh.capabilities().available()
    result = Experiment.on(igh_gapped).recombine().run_records(n=50, seed=5, validate_records=True)
    assert result.reference_content_hash == igh_gapped.execution_hash


# ── VJ light chain with local allele names ───────────────────────────


def test_igk_local_vj_accounting_names_and_explanations(igk):
    v_records = [r for r in igk.records if r.id.startswith("rec:V:")]
    assert [(r.header, r.label, r.status, r.reason, r.duplicate_of) for r in v_records] == [
        ("kappa_v_a", "kappa_v_a", "accepted", None, None),
        ("kappa_v_b", "kappa_v_b", "accepted", None, None),
        ("kappa_v_a", "kappa_v_a", "duplicate", None, "rec:V:0"),
        ("kappa_v_bad", "kappa_v_bad", "rejected", "unsupported_base: 'N' at position 100", None),
        ("kappa,v", None, "rejected", "unrepresentable_label: contains ',' at position 5", None),
    ]
    assert [(r.label, r.status) for r in igk.records if r.id.startswith("rec:J:")] == [
        ("kappa_j_1", "accepted"), ("kappa_j_2", "accepted")]
    assert [a.label for a in igk.alleles] == ["kappa_v_a", "kappa_v_b", "kappa_j_1", "kappa_j_2"]
    assert (igk.chain_type, igk.locus) == ("vj", "IGK")

    caps = igk.capabilities()
    assert caps.available() == ("catalogue_inspection", "observed_span_simulation")
    explanation = caps.explain()
    lines = explanation.splitlines()
    assert any(line.startswith("- V:kappa_v_b: accepted anchor") for line in lines)
    (j2_line,) = {line for line in lines if line.startswith("- J:kappa_j_2: accepted anchor")}
    assert "candidate at 8" in j2_line

    result = Experiment.on(igk).recombine().run_records(n=200, seed=11, validate_records=True)
    assert result.reference_content_hash == igk.execution_hash
    called = {name for record in result for name in _names(record["v_call"])}
    assert called == {"kappa_v_a", "kappa_v_b"}   # only the local view labels, and both of them
    using_b = [record for record in result if _names(record["v_call"]) == ["kappa_v_b"]]
    assert using_b
    for record in using_b:
        assert record["productive"] is None
        assert "coding_phase_unresolved" in record["productivity_reasons"]


def test_igk_relationship_only_entity_is_visible_not_executable(igk):
    assert [(e.id, e.name, e.executable) for e in igk.entities] == [
        ("entity:lab_alias_kva", "lab_alias_kva", False)]
    report = igk.report()
    assert [r.to_dict() for r in report.relationships] == [
        {"kind": "alias", "allele": "V:kappa_v_a", "entity": "entity:lab_alias_kva", "source": None}]
    assert "  lab_alias_kva: alias of V:kappa_v_a" in str(report).splitlines()
    assert "lab_alias_kva" not in {a.label for a in igk.alleles}

    native = igk.native_refdata()
    assert native.v_pool_size() == 2
    assert [native.v_allele(i).name for i in range(native.v_pool_size())] == ["kappa_v_a", "kappa_v_b"]
    assert json.loads(igk.to_json())["models"]["allele_usage"]["V"] == ["kappa_v_a", "kappa_v_b"]


# ── An unlisted species needs no registry ────────────────────────────


def test_unlisted_species_reference_runs(tmp_path, igk):
    # A software fixture with a species GenAIRR does not list; not dolphin biology.
    side = json.loads((IGK / "annotations.json").read_text(encoding="utf-8"))
    side["reference"]["name"] = "synthetic_unlisted_species_fixture"
    side["reference"]["taxon"] = {"scientific_name": "Tursiops truncatus", "taxonomy_id": 9739}
    directory = tmp_path / "unlisted"
    directory.mkdir()
    for name in ("v.fasta", "j.fasta"):
        shutil.copyfile(IGK / name, directory / name)
    (directory / "annotations.json").write_bytes((json.dumps(side, indent=2) + "\n").encode("utf-8"))
    ref = build_igk(directory)

    assert ref.name == "synthetic_unlisted_species_fixture"
    assert (ref.taxon.scientific_name, ref.taxon.taxonomy_id) == ("Tursiops truncatus", 9739)
    species = ref.native_refdata().identity()["species"]
    assert species == "Tursiops truncatus"
    assert species not in {member.value for member in Species}
    assert ref.catalogue_hash != igk.catalogue_hash

    result = Experiment.on(ref).recombine().run_records(n=20, seed=3, validate_records=True)
    assert len(result) == 20
    assert result.reference_content_hash == ref.execution_hash
    assert result.validate_records(ref).ok


# ── Save and reload ──────────────────────────────────────────────────


@pytest.mark.parametrize("fixture", ["igh", "igh_gapped", "igk"])
def test_save_reload_preserves_identity_order_annotations(fixture, request, tmp_path):
    ref = request.getfixturevalue(fixture)
    path = tmp_path / f"{fixture}.genairr-ref.json"
    ref.write(path)
    assert path.read_bytes() == ref.to_json().encode("utf-8")
    again = Reference.read(path)

    assert again.to_json() == ref.to_json()
    assert Reference.from_json(ref.to_json()).to_json() == ref.to_json()
    for name in ("name", "chain_type", "locus", "taxon", "catalogue_hash", "models_hash", "execution_hash",
                 "artifact_hash", "records", "entities", "findings", "admission"):
        assert getattr(again, name) == getattr(ref, name), name
    assert [a.to_dict() for a in again.alleles] == [a.to_dict() for a in ref.alleles]
    assert [(a.segment, a.label) for a in again.alleles] == [(a.segment, a.label) for a in ref.alleles]
    assert again.capabilities() == ref.capabilities()
    assert again.report().to_dict() == ref.report().to_dict()
    assert again.recombine_defaults() == ref.recombine_defaults()

    native, original = again.native_refdata(), ref.native_refdata()
    assert native.content_hash() == original.content_hash() == ref.execution_hash
    for segment in ("v", "d", "j"):
        size = getattr(native, f"{segment}_pool_size")()
        assert size == getattr(original, f"{segment}_pool_size")()
        names = [getattr(native, f"{segment}_allele")(i).name for i in range(size)]
        assert names == [a.label for a in ref.alleles if a.segment == segment.upper()]


# ── Golden artifact and portability ──────────────────────────────────


def test_golden_igh_artifact(igh):
    golden = GOLDEN.read_bytes()
    doc = json.loads(igh.to_json())
    doc["generator"]["version"] = json.loads(golden.decode("utf-8"))["generator"]["version"]
    doc["commitment"] = typed_digest({k: v for k, v in doc.items() if k != "commitment"}, ARTIFACT_DOMAIN)
    assert (canonical_json(doc) + "\n").encode("utf-8") == golden
    assert doc["catalogue_hash"] == IGH_CATALOGUE_HASH
    assert doc["native"]["execution_hash"] == IGH_EXECUTION_HASH
    assert (igh.catalogue_hash, igh.execution_hash) == (IGH_CATALOGUE_HASH, IGH_EXECUTION_HASH)

    stored = Reference.read(GOLDEN)
    assert (stored.catalogue_hash, stored.execution_hash) == (IGH_CATALOGUE_HASH, IGH_EXECUTION_HASH)
    assert stored.artifact_hash == doc["commitment"]


def test_artifact_contains_no_local_paths(igh, tmp_path):
    directory = tmp_path / "private_lab_share" / "igh_build_area"
    directory.mkdir(parents=True)
    for name in ("v.fasta", "d.fasta", "j.fasta", "annotations.json"):
        shutil.copyfile(IGH / name, directory / name)
    ref = build_igh(directory)
    text = ref.to_json()

    doc = json.loads(text)
    assert [s["file_name"] for s in doc["catalogue"]["sources"]] == ["v.fasta", "d.fasta", "j.fasta"]
    for needle in (str(tmp_path), "private_lab_share", "igh_build_area", str(TESTS.parent),
                   str(Path.home()), "reference_first_use"):
        assert needle not in text, needle
    # No string is an absolute POSIX, Windows-drive or UNC path ("V:label" ids are not: no root).
    assert [s for s in _strings(doc) if PurePosixPath(s).is_absolute() or PureWindowsPath(s).is_absolute()] == []
    # The artifact does not depend on where its inputs were read from.
    assert text == igh.to_json()
    # Nor on the directory it is written to or read from.
    elsewhere = tmp_path / "other" / "place.genairr-ref.json"
    elsewhere.parent.mkdir()
    ref.write(elsewhere)
    assert Reference.read(elsewhere).to_json() == text
    assert os.path.basename(elsewhere) not in text
