import copy
import dataclasses
import datetime
import errno
import hashlib
import json
import os
import pathlib
import re
import stat

import pytest

import GenAIRR
from GenAIRR import Reference, ReferenceArtifactError, ReferenceCapabilityError
from GenAIRR._artifact_io import canonical_json, typed_digest
from GenAIRR._reference import artifact as artifact_module
from GenAIRR._reference.artifact import ClaimSummary, ReferenceAllele, ReferenceReport
from GenAIRR._reference.defaults import models_hash
from GenAIRR._reference.lowering import PROFILE_ID, lower, seal
from GenAIRR._reference.model import Catalogue
from tests.support import paths

# A human IGHV1-2*02-like V (TGT Cys codon at 285), two IGHD-like Ds and two IGHJ-like Js
# (TGG-GGC-CAG/A-GG motif at 14 and 29).  Software fixtures, not evidence.
V1 = ("CAGGTGCAGCTGGTGCAGTCTGGGGCTGAGGTGAAGAAGCCTGGGGCCTCAGTGAAGGTCTCCTGCAAGGCTTCTGGATACACCTTCACC"
      "GGCTACTATATGCACTGGGTGCGACAGGCCCCTGGACAAGGGCTTGAGTGGATGGGATGGATCAACCCTAACAGTGGTGGCACAAACTAT"
      "GCACAGAAGTTTCAGGGCAGGGTCACCATGACCAGGGACACGTCCATCAGCACAGCCTACATGGAGCTGAGCAGGCTGAGATCTGACGAC"
      "ACGGCCGTGTATTACTGTGCGAGAGA")
V2 = "G" + V1[1:]
V_BAD = V1[:100] + "N" + V1[101:]
D1 = "GTATTACTATGGTTCGGGGAGTTATTATAAC"
D2 = "AGGATATTGTAGTAGTACCAGCTGCTATGCC"
J1 = "ACTACTTTGACTACTGGGGCCAGGGAACCCTGGTCACCGTCTCCTCAG"
J2 = "ATTACTACTACTACTACGGTATGGACGTCTGGGGCCAAGGGACCACGGTCACCGTCTCCTCAG"
GAPPED_V = "AAA" + "..." + "AAA" * 101 + "TGT" + "AAA"   # IMGT Cys104 lands at ungapped 306

V_FASTA = (f">IGHV1-2*02 lab clone\n{V1}\n>IGHV1-2*02 second copy\n{V1}\n>IGHV1-3*01\n{V2}\n"
           f">IGHV_bad\n{V_BAD}\n>IGHV1-2_copy\n{V1}\n")
D_FASTA = f">IGHD3-10*01\n{D1}\n>IGHD2-2*01\n{D2}\n"
J_FASTA = f">IGHJ4*02\n{J1}\n>IGHJ6*02\n{J2}\n"
SOURCE = "doi:10.1000/lab-reference"
ARTIFACT_DOMAIN = "genairr_reference_artifact_v1"
GOLDEN = paths.golden("reference", "human_igh.genairr-ref.json")


def digest(seq):
    return "sha256:" + hashlib.sha256(seq.encode("ascii")).hexdigest()


def ann(value, support="assertion", source=None, method="test fixture"):
    return {"value": value, "support": support, "source": source, "method": method}


def phase(value, derivation="anchor_derived"):
    return {**ann(value), "derivation": derivation}


def entry(segment, label, seq, **annotations):
    base = {"segment": segment, "label": label, "sequence_sha256": digest(seq), "gene": None,
            "functionality": None, "anchor": None, "coding_phase": None,
            "boundary_5": None, "boundary_3": None, "notes": []}
    base.update(annotations)
    return base


def header(chain_type="vdj", locus="IGH", build_date=None):
    return {"name": "lab_igh", "version": "1", "locus": locus, "chain_type": chain_type,
            "taxon": {"scientific_name": "Homo sapiens", "taxonomy_id": 9606}, "build_date": build_date}


def default_alleles(v_anchor_support="assertion"):
    v_src = SOURCE if v_anchor_support == "source" else None
    return [
        entry("V", "IGHV1-2*02", V1, anchor=ann(285, v_anchor_support, v_src),
              functionality=ann("functional", "source", SOURCE), coding_phase=phase(0),
              boundary_3=ann("complete")),
        entry("V", "IGHV1-3*01", V2, anchor=ann(285, v_anchor_support, v_src)),
        entry("V", "IGHV1-2_copy", V1, anchor=ann(285, v_anchor_support, v_src)),
        entry("J", "IGHJ4*02", J1, anchor=ann(14, "source", SOURCE),
              boundary_5=ann("complete", "source", SOURCE)),
        entry("J", "IGHJ6*02", J2, anchor=ann(29, v_anchor_support, v_src)),
    ]


def sidecar(alleles, *, chain_type="vdj", locus="IGH", v_alignment="none", build_date=None, relationships=()):
    return {"schema": "genairr_reference_annotations", "schema_version": 1,
            "reference": header(chain_type, locus, build_date),
            "conventions": {"coordinates": "zero_based_half_open", "v_alignment": v_alignment},
            "alleles": list(alleles), "relationships": list(relationships), "suggestions": []}


def fasta_files(directory, *, v=V_FASTA, d=D_FASTA, j=J_FASTA):
    directory.mkdir(parents=True, exist_ok=True)
    paths = {}
    for name, text in (("v", v), ("d", d), ("j", j)):
        if text is not None:
            paths[name] = directory / f"{name}.fasta"
            paths[name].write_bytes(text.encode("utf-8"))
    return paths


def build(directory, *, alleles=None, v=V_FASTA, d=D_FASTA, j=J_FASTA, **side_kw):
    paths = fasta_files(directory, v=v, d=d, j=j)
    side = sidecar(default_alleles() if alleles is None else alleles, **side_kw)
    annotations = directory / "annotations.json"
    annotations.write_bytes(json.dumps(side).encode("utf-8"))
    return Reference.from_fasta(v=paths["v"], j=paths["j"], d=paths.get("d"), annotations=annotations)


@pytest.fixture
def ref(tmp_path):
    return build(tmp_path / "ref")


def doc_of(reference):
    return json.loads(reference.to_json())


def _commit(doc):
    body = {k: v for k, v in doc.items() if k != "commitment"}
    doc["commitment"] = typed_digest(body, ARTIFACT_DOMAIN)
    return doc


def _recommit(doc):
    """Recompute catalogue_hash, models_hash and the commitment; leaves native.context stale."""
    doc["catalogue_hash"] = typed_digest(doc["catalogue"], "genairr_reference_catalogue_v1")
    doc["models_hash"] = typed_digest(doc["models"], "genairr_reference_models_v1")
    return _commit(doc)


def _full_recommit(doc):
    """Recompute every hash plus admission, native projection, context and execution hash.

    This is what a careful forger would do, so only the regeneration of derived catalogue
    content (and the catalogue's own validation) can reject such an artifact.
    """
    validate = Catalogue.validate
    Catalogue.validate = lambda self: None
    try:
        catalogue = Catalogue.from_dict(doc["catalogue"])
    finally:
        Catalogue.validate = validate
    doc["catalogue_hash"] = catalogue.digest()
    doc["models_hash"] = models_hash(doc["models"])
    lowered = lower(catalogue, doc["view"]["alleles"])
    doc["view"]["admission"] = lowered.admission
    doc["native"]["reference_spec"] = lowered.spec
    context = {"schema_version": 1, "catalogue_hash": doc["catalogue_hash"], "view": "default",
               "profile": PROFILE_ID, "models_hash": doc["models_hash"]}
    doc["native"]["context"] = context
    doc["native"]["execution_hash"] = seal(lowered, context).content_hash()
    return _commit(doc)


def claim_of(doc, subject, predicate, derivation):
    (claim,) = [c for c in doc["catalogue"]["claims"] if (c["subject"], c["predicate"], c["derivation"])
                == (subject, predicate, derivation)]
    return claim


def load(doc):
    return Reference.from_json(canonical_json(doc) + "\n")


def rejected(doc, kind, match=None):
    with pytest.raises(ReferenceArtifactError) as e:
        load(doc)
    assert e.value.kind == kind, (e.value.kind, str(e.value))
    if match is not None:
        assert re.search(match, str(e.value)), str(e.value)
    return e.value


# ── Serialization ────────────────────────────────────────────────────


def test_write_read_round_trip_is_byte_identical(ref, tmp_path):
    path = tmp_path / "lab.genairr-ref.json"
    ref.write(path)
    data = path.read_bytes()
    assert data == ref.to_json().encode("utf-8")
    assert data.endswith(b"\n") and data.count(b"\n") == 1
    again = Reference.read(path)
    assert again.to_json() == ref.to_json()
    for name in ("name", "chain_type", "locus", "taxon", "catalogue_hash", "models_hash",
                 "execution_hash", "artifact_hash", "alleles", "records", "entities", "findings",
                 "admission"):
        assert getattr(again, name) == getattr(ref, name), name
    assert again.native_refdata().content_hash() == ref.execution_hash
    again.write(tmp_path / "copy.genairr-ref.json")
    assert (tmp_path / "copy.genairr-ref.json").read_bytes() == data
    assert Reference.from_json(data.decode("utf-8")).artifact_hash == ref.artifact_hash


def test_artifact_layout(ref):
    doc = doc_of(ref)
    assert list(doc) == sorted(doc)
    assert set(doc) == {"artifact", "schema_version", "generator", "catalogue", "catalogue_hash",
                        "view", "models", "models_hash", "native", "commitment"}
    assert doc["artifact"] == "genairr_reference" and doc["schema_version"] == 1
    assert doc["generator"] == {"name": "GenAIRR", "version": GenAIRR.__version__}
    assert set(doc["view"]) == {"id", "profile", "alleles", "admission"}
    assert doc["view"]["id"] == "default" and doc["view"]["profile"] == PROFILE_ID
    assert doc["view"]["alleles"] == {"V": ["V:IGHV1-2*02", "V:IGHV1-3*01", "V:IGHV1-2_copy"],
                                      "D": ["D:IGHD3-10*01", "D:IGHD2-2*01"],
                                      "J": ["J:IGHJ4*02", "J:IGHJ6*02"]}
    assert set(doc["native"]) == {"reference_spec", "context", "execution_hash"}
    assert doc["native"]["context"] == {"schema_version": 1, "catalogue_hash": ref.catalogue_hash,
                                        "view": "default", "profile": PROFILE_ID,
                                        "models_hash": ref.models_hash}
    assert doc["catalogue_hash"] == ref.catalogue_hash
    assert doc["models_hash"] == ref.models_hash
    assert doc["native"]["execution_hash"] == ref.execution_hash
    assert doc["commitment"] == ref.artifact_hash
    body = {k: v for k, v in doc.items() if k != "commitment"}
    assert typed_digest(body, ARTIFACT_DOMAIN) == ref.artifact_hash
    assert doc["models"]["profile"] == "genairr.unfitted_defaults.v1"
    assert doc["models"]["allele_usage"]["V"] == ["IGHV1-2*02", "IGHV1-3*01", "IGHV1-2_copy"]


def test_build_twice_is_deterministic(tmp_path):
    first = build(tmp_path / "a")
    second = build(tmp_path / "b" / "elsewhere")
    assert first.to_json() == second.to_json()
    assert first.artifact_hash == second.artifact_hash
    assert first.execution_hash == second.execution_hash
    assert str(tmp_path) not in first.to_json()


def test_no_wall_clock_date(tmp_path):
    text = build(tmp_path / "undated").to_json()
    assert datetime.date.today().isoformat() not in text
    assert not re.search(r"\d{4}-\d{2}-\d{2}", text)
    dated = build(tmp_path / "dated", build_date="2001-02-03").to_json()
    assert re.findall(r"\d{4}-\d{2}-\d{2}", dated) == ["2001-02-03"]


def test_vj_reference_round_trip(tmp_path):
    j_kappa = "TGGACGTTCGGCCAAGGGACCAAGGTGGAAATCAAAC"
    alleles = [entry("V", "kv1", V1, anchor=ann(285)), entry("J", "kj1", j_kappa, anchor=ann(7))]
    ref = build(tmp_path, alleles=alleles, v=f">kv1\n{V1}\n", d=None, j=f">kj1\n{j_kappa}\n",
                chain_type="vj", locus="IGK")
    assert ref.chain_type == "vj" and ref.locus == "IGK"
    assert [a.segment for a in ref.alleles] == ["V", "J"]
    doc = doc_of(ref)
    assert doc["view"]["alleles"] == {"V": ["V:kv1"], "J": ["J:kj1"]}
    assert "D" not in doc["models"]["allele_usage"]
    assert Reference.from_json(ref.to_json()).to_json() == ref.to_json()


# ── Integrity and meaning checks ─────────────────────────────────────


def test_commitment_tamper_detected(ref):
    doc = doc_of(ref)
    doc["catalogue"]["header"]["name"] = "someone_else"
    rejected(doc, "integrity", "commitment")
    doc = doc_of(ref)
    doc["commitment"] = "sha256:" + "0" * 64
    rejected(doc, "integrity", "commitment")
    doc = doc_of(ref)
    doc["generator"]["version"] = "999.0"
    rejected(doc, "integrity", "commitment")
    assert load(_commit(doc)).to_json() != ref.to_json()  # a different generator version still loads


def test_recommitted_np_support_0_to_100_rejected(ref):
    doc = doc_of(ref)
    doc["models"]["np_lengths"]["NP1"]["support"] = list(range(101))
    rejected(_recommit(doc), "integrity", "models differ from profile genairr.unfitted_defaults.v1")


def test_recommitted_models_hash_mismatch_rejected(ref):
    doc = doc_of(ref)
    doc["models_hash"] = "sha256:" + "1" * 64
    err = rejected(_commit(doc), "integrity")
    assert err.path == "models_hash"


def test_catalogue_hash_mismatch_rejected(ref):
    doc = doc_of(ref)
    doc["catalogue_hash"] = "sha256:" + "2" * 64
    err = rejected(_commit(doc), "integrity")
    assert err.path == "catalogue_hash"


def test_recommitted_admission_edit_rejected(ref):
    doc = doc_of(ref)
    assert doc["view"]["admission"] == {"profile": PROFILE_ID, "accepted_curatable_findings": []}
    doc["view"]["admission"]["accepted_curatable_findings"].append(
        {"kind": "MissingAnchor", "segment": "V", "allele": "V:IGHV1-3*01"})
    err = rejected(_recommit(doc), "integrity", "admission")
    assert err.path == "view.admission"


def test_recommitted_profile_id_edit_rejected(ref):
    doc = doc_of(ref)
    doc["view"]["profile"] = "canonical_productive"
    err = rejected(_recommit(doc), "integrity")
    assert err.path == "view.profile"
    doc = doc_of(ref)
    doc["view"]["id"] = "custom"
    assert rejected(_recommit(doc), "integrity").path == "view.id"


def test_recommitted_view_reorder_rejected(ref):
    doc = doc_of(ref)
    doc["view"]["alleles"]["V"].reverse()
    err = rejected(_recommit(doc), "integrity")
    assert err.path == "view.alleles"
    doc = doc_of(ref)
    doc["view"]["alleles"]["J"].pop()
    assert rejected(_recommit(doc), "integrity").path == "view.alleles"


def test_recommitted_claim_to_missing_allele_rejected(ref):
    doc = doc_of(ref)
    doc["catalogue"]["claims"].append({
        "subject": "V:not_in_catalogue", "predicate": "anchor", "value": 3, "support": "assertion",
        "derivation": "stated", "status": "accepted", "source": None, "method": "m"})
    err = rejected(_recommit(doc), "integrity", "unknown allele")
    assert err.path.startswith("catalogue.claims[")


def test_recommitted_sequence_edit_without_digest_rejected(ref):
    doc = doc_of(ref)
    allele = doc["catalogue"]["alleles"][0]
    allele["sequence"] = ("A" if allele["sequence"][0] != "A" else "C") + allele["sequence"][1:]
    err = rejected(_recommit(doc), "integrity", "digest does not match sequence")
    assert err.path == "catalogue.alleles[0].sequence_sha256"


def test_recommitted_evidence_upgrade_rejected(ref):
    # A consistent-looking catalogue edit that upgrades evidence must not survive the native projection.
    doc = doc_of(ref)
    for claim in doc["catalogue"]["claims"]:
        if claim["subject"] == "V:IGHV1-3*01" and claim["predicate"] == "anchor":
            claim.update(support="source", source=SOURCE)
    rejected(_recommit(doc), "integrity", "native projection differs from catalogue lowering")


def test_recommitted_spec_edit_rejected(ref):
    doc = doc_of(ref)
    doc["native"]["reference_spec"]["v"][1]["anchor"] = 282
    err = rejected(_recommit(doc), "integrity", "native projection differs from catalogue lowering")
    assert err.path == "native.reference_spec"
    doc = doc_of(ref)
    doc["native"]["reference_spec"]["v"][0]["evidence"]["anchor_status"] = "verified"
    rejected(_recommit(doc), "integrity", "native projection differs from catalogue lowering")


def test_recommitted_context_edit_rejected(ref):
    doc = doc_of(ref)
    doc["native"]["context"]["view"] = "custom"
    err = rejected(_recommit(doc), "integrity")
    assert err.path == "native.context"
    doc = doc_of(ref)
    doc["native"]["context"]["schema_version"] = True
    assert rejected(_recommit(doc), "integrity").path == "native.context"


def test_recommitted_execution_hash_edit_rejected(ref):
    doc = doc_of(ref)
    doc["native"]["execution_hash"] = "sha256:" + "3" * 64
    err = rejected(_recommit(doc), "integrity")
    assert err.path == "native.execution_hash"


def test_duplicate_key_and_unknown_nested_key_rejected(ref):
    text = ref.to_json()
    duplicated = text.replace('"artifact":"genairr_reference"',
                              '"artifact":"genairr_reference","artifact":"genairr_reference"', 1)
    assert duplicated != text
    with pytest.raises(ReferenceArtifactError) as e:
        Reference.from_json(duplicated)
    assert e.value.kind == "malformed"
    for section in ("generator", "view", "native"):
        doc = doc_of(ref)
        doc[section]["extra"] = 1
        assert rejected(_recommit(doc), "malformed", "unknown keys").path == section
    doc = doc_of(ref)
    doc["extra"] = 1
    rejected(_recommit(doc), "malformed", "unknown keys")
    doc = doc_of(ref)
    del doc["view"]["admission"]
    rejected(_recommit(doc), "malformed", "missing keys")
    doc = doc_of(ref)
    doc["catalogue"]["alleles"][0]["extra"] = 1
    err = rejected(_recommit(doc), "malformed", "unknown keys")
    assert err.path == "catalogue.alleles[0]"
    for bad in ("[]", "not json", '{"artifact": NaN}'):
        with pytest.raises(ReferenceArtifactError) as e:
            Reference.from_json(bad)
        assert e.value.kind == "malformed"


def test_schema_version_2_unsupported(ref):
    doc = doc_of(ref)
    doc["schema_version"] = 2
    rejected(_commit(doc), "unsupported_schema", "schema version 2")
    for bad in ("1", True, 1.0, None):
        doc = doc_of(ref)
        doc["schema_version"] = bad
        rejected(_commit(doc), "unsupported_schema")
    doc = doc_of(ref)
    del doc["schema_version"]
    rejected(doc, "unsupported_schema")
    doc = doc_of(ref)
    doc["artifact"] = "genairr_family_replay"
    assert rejected(_commit(doc), "malformed", "not a GenAIRR reference artifact").path == "artifact"
    doc = doc_of(ref)
    del doc["artifact"]
    rejected(doc, "malformed", "not a GenAIRR reference artifact")
    doc = doc_of(ref)
    doc["generator"]["name"] = "Other"
    assert rejected(_commit(doc), "malformed").path == "generator.name"


def test_v2_shaped_artifact_with_extra_key_is_unsupported_schema(ref):
    doc = doc_of(ref)
    doc["schema_version"] = 2
    doc["lineage"] = {"parents": []}
    del doc["models_hash"]
    err = rejected(doc, "unsupported_schema", "schema version 2")
    assert err.path == "schema_version"


def test_size_bound(ref, tmp_path, monkeypatch):
    path = tmp_path / "big.genairr-ref.json"
    ref.write(path)
    monkeypatch.setattr(artifact_module, "MAX_BYTES", 100)
    with pytest.raises(ReferenceArtifactError) as e:
        Reference.from_json(ref.to_json())
    assert e.value.kind == "io" and "limit" in str(e.value)
    with pytest.raises(ReferenceArtifactError) as e:
        Reference.read(path)
    assert e.value.kind == "io"
    monkeypatch.undo()
    assert artifact_module.MAX_BYTES == 64 * 1024 * 1024


def test_read_missing_file_is_io_error(tmp_path):
    with pytest.raises(ReferenceArtifactError) as e:
        Reference.read(tmp_path / "missing.genairr-ref.json")
    assert e.value.kind == "io"


def test_write_refuses_overwrite_unless_replace(ref, tmp_path):
    path = tmp_path / "lab.genairr-ref.json"
    path.write_bytes(b"previous")
    with pytest.raises(ReferenceArtifactError) as e:
        ref.write(path)
    assert e.value.kind == "io"
    assert path.read_text(encoding="utf-8") == "previous"
    ref.write(path, replace=True)
    assert path.read_bytes() == ref.to_json().encode("utf-8")
    with pytest.raises(TypeError):
        ref.write(path, replace="yes")
    with pytest.raises(ReferenceArtifactError) as e:
        ref.write(tmp_path / "no_such_dir" / "x.genairr-ref.json")
    assert e.value.kind == "io"


def test_write_without_hard_links_is_an_io_error_unless_replacing(ref, tmp_path, monkeypatch):
    # A filesystem without hard links (FAT/exFAT, some SMB shares) offers no atomic
    # no-overwrite publication outside Windows.
    from GenAIRR import _artifact_io

    def link(source, destination, *args, **kwargs):
        raise PermissionError(errno.EPERM, "Operation not permitted")

    monkeypatch.setattr(_artifact_io, "_WINDOWS", False, raising=False)
    monkeypatch.setattr(os, "link", link)
    out = tmp_path / "out"
    out.mkdir()
    path = out / "lab.genairr-ref.json"
    with pytest.raises(ReferenceArtifactError) as e:
        ref.write(path)
    assert e.value.kind == "io" and e.value.path == path.name
    assert "does not support atomic no-overwrite publication" in e.value.message
    assert list(out.iterdir()) == []
    ref.write(path, replace=True)
    assert path.read_bytes() == ref.to_json().encode("utf-8")


@pytest.mark.skipif(os.name == "nt", reason="POSIX permission bits")
@pytest.mark.parametrize("mask,mode", [(0o022, 0o644), (0o077, 0o600)])
def test_written_artifact_mode_follows_umask(ref, tmp_path, mask, mode):
    # A reference written on a shared host is as readable as any file the user writes.
    path = tmp_path / "lab.genairr-ref.json"
    previous = os.umask(mask)
    try:
        ref.write(path)
        ref.write(path, replace=True)
    finally:
        os.umask(previous)
    assert stat.S_IMODE(path.stat().st_mode) == mode
    assert Reference.read(path).to_json() == ref.to_json()


# ── Object contract ──────────────────────────────────────────────────


def test_reference_is_immutable_and_native_sealed(ref):
    with pytest.raises(AttributeError):
        ref.name = "other"
    with pytest.raises(AttributeError):
        ref._catalogue = None
    with pytest.raises(AttributeError):
        del ref._catalogue
    with pytest.raises(AttributeError):
        ref.anything = 1
    assert not hasattr(ref, "__dict__")
    assert copy.copy(ref) is ref and copy.deepcopy(ref) is ref
    with pytest.raises(TypeError):
        Reference()
    for name in ("alleles", "records", "entities", "findings"):
        assert isinstance(getattr(ref, name), tuple), name
    with pytest.raises(dataclasses.FrozenInstanceError):
        ref.alleles[0].label = "x"
    admission = ref.admission
    admission["profile"] = "changed"
    admission["accepted_curatable_findings"].append({})
    assert ref.admission == {"profile": PROFILE_ID, "accepted_curatable_findings": []}
    defaults = ref.recombine_defaults()
    defaults["np1"].clear()
    defaults["allele_usage_v"]["IGHV1-2*02"] = 99.0
    assert ref.recombine_defaults()["np1"] == [(v, 1.0) for v in range(7)]
    assert ref.recombine_defaults()["allele_usage_v"] == {
        "IGHV1-2*02": 1.0, "IGHV1-3*01": 1.0, "IGHV1-2_copy": 1.0}
    assert ref.recombine_defaults()["trim_v_3"] is None
    native = ref.native_refdata()
    assert native is ref.native_refdata()
    assert native.is_sealed()
    assert native.content_hash() == ref.execution_hash
    assert json.loads(native.reference_context_json())["catalogue_hash"] == ref.catalogue_hash
    with pytest.raises(ValueError, match="sealed"):
        native.add_v_allele("Vx", "Vx", V1.encode("ascii"), anchor=285)
    with pytest.raises(ValueError, match="sealed"):
        native.set_identity(name="x")
    assert [native.v_allele(i).name for i in range(native.v_pool_size())] == \
        ["IGHV1-2*02", "IGHV1-3*01", "IGHV1-2_copy"]


def test_identity_properties(ref):
    assert (ref.name, ref.chain_type, ref.locus) == ("lab_igh", "vdj", "IGH")
    assert ref.taxon.scientific_name == "Homo sapiens" and ref.taxon.taxonomy_id == 9606
    for h in (ref.catalogue_hash, ref.models_hash, ref.execution_hash, ref.artifact_hash):
        assert re.fullmatch(r"sha256:[0-9a-f]{64}", h)
    assert len({ref.catalogue_hash, ref.models_hash, ref.execution_hash, ref.artifact_hash}) == 4
    assert [r.status for r in ref.records if r.id.startswith("rec:V:")] == \
        ["accepted", "duplicate", "accepted", "rejected", "accepted"]
    assert [f.kind for f in ref.findings] == ["identical_sequence_distinct_labels"]
    assert ref.entities == ()


def test_alleles_in_view_order_with_recorded_evidence(ref):
    assert [(a.segment, a.label) for a in ref.alleles] == [
        ("V", "IGHV1-2*02"), ("V", "IGHV1-3*01"), ("V", "IGHV1-2_copy"),
        ("D", "IGHD3-10*01"), ("D", "IGHD2-2*01"), ("J", "IGHJ4*02"), ("J", "IGHJ6*02")]
    v = ref.alleles[0]
    assert isinstance(v, ReferenceAllele)
    assert (v.gene, v.length, v.sequence_sha256) == ("IGHV1-2", len(V1), digest(V1))
    assert v.anchor == ClaimSummary(285, "assertion", "stated", None)
    assert v.functionality == ClaimSummary("functional", "source", "stated", SOURCE)
    assert v.coding_phase == ClaimSummary(0, "assertion", "anchor_derived", None)
    assert v.boundary_3 == ClaimSummary("complete", "assertion", "stated", None)
    assert v.boundary_5 is None
    d = ref.alleles[3]
    assert d.anchor is d.functionality is d.coding_phase is d.boundary_5 is d.boundary_3 is None
    j = ref.alleles[5]
    assert j.anchor == ClaimSummary(14, "source", "stated", SOURCE)


def test_allele_summary_never_overstates_mixed_or_conflicting_claims(tmp_path):
    ungapped = GAPPED_V.replace(".", "")
    alleles = [entry("V", "gv", ungapped, anchor=ann(306, "source", SOURCE)),
               entry("V", "gw", ungapped[:-3] + "TTT", anchor=ann(303, "source", SOURCE))] + \
        default_alleles()[3:]
    v = f">gv\n{GAPPED_V}\n>gw\n{GAPPED_V[:-3]}TTT\n"
    ref = build(tmp_path, alleles=alleles, v=v, v_alignment="imgt_gapped")
    agreeing, conflicting = ref.alleles[0], ref.alleles[1]
    # A source-stated anchor agreeing with an alignment assertion is only as strong as the assertion.
    assert agreeing.anchor == ClaimSummary(306, "assertion", "imgt_alignment", None)
    assert conflicting.anchor is None
    text = str(ref.report())
    assert "conflict" in text
    assert ref.report().conflicts == (("V:gw", "anchor"),)


def test_require_evidence_levels(tmp_path):
    ref = build(tmp_path / "mixed")
    ref.require()
    ref.require("catalogue_inspection", "observed_span_simulation", evidence="source")
    ref.require("canonical_productive")
    assert ref.capabilities().status("canonical_productive").evidence == "assertion"
    # Most anchors are the author's assertions, so the junction is not source-annotated.
    assert not ref.capabilities().status("source_annotated_junction").available
    for evidence in ("assertion", "source"):
        with pytest.raises(ReferenceCapabilityError) as e:
            ref.require("catalogue_inspection", "source_annotated_junction", evidence=evidence)
        assert isinstance(e.value, ReferenceArtifactError)
        assert e.value.kind == "capability" and e.value.path == "source_annotated_junction"
        assert e.value.message == ref.capabilities().explain()
    with pytest.raises(ReferenceCapabilityError) as e:
        ref.require("original_coding_ends")
    assert e.value.path == "original_coding_ends"
    with pytest.raises(ReferenceCapabilityError):
        ref.require("duplex_family_replay")
    with pytest.raises(ValueError, match="unknown capability"):
        ref.require("orignal_coding_ends")
    with pytest.raises(ValueError):
        ref.require("catalogue_inspection", evidence="strong")
    sourced = build(tmp_path / "sourced", alleles=default_alleles(v_anchor_support="source"))
    sourced.require("source_annotated_junction", evidence="source")
    assert sourced.capabilities().status("source_annotated_junction").evidence == "source"


def test_report_lists_accounting_support_and_capabilities(ref):
    report = ref.report()
    assert isinstance(report, ReferenceReport)
    with pytest.raises(dataclasses.FrozenInstanceError):
        report.name = "x"
    text = str(report)
    for value in ("lab_igh", "Homo sapiens", "IGH", "vdj", ref.catalogue_hash, ref.models_hash,
                  ref.execution_hash, ref.artifact_hash):
        assert value in text, value
    assert "V: 5 records: 3 accepted, 1 duplicate, 1 rejected" in text
    assert "D: 2 records: 2 accepted, 0 duplicate, 0 rejected" in text
    assert "J: 2 records: 2 accepted, 0 duplicate, 0 rejected" in text
    assert "duplicate of rec:V:0" in text
    assert "unsupported_base: 'N' at position 100" in text
    assert "285 assertion/stated" in text
    assert "14 source/stated" in text
    assert "functional source/stated" in text
    assert "0 assertion/anchor_derived" in text
    assert "identical_sequence_distinct_labels" in text and "V:IGHV1-2_copy" in text
    assert "genairr.unfitted_defaults.v1" in text
    assert "uniform over 0..6" in text and "disabled" in text
    assert ref.capabilities().explain() in text
    assert "verified" not in text  # no evidence is promoted beyond its recorded support
    data = report.to_dict()
    json.dumps(data)
    assert data["hashes"] == {"catalogue_hash": ref.catalogue_hash, "models_hash": ref.models_hash,
                              "execution_hash": ref.execution_hash, "artifact_hash": ref.artifact_hash}
    assert data["records"]["counts"]["V"] == {"accepted": 3, "duplicate": 1, "rejected": 1}
    assert [r["id"] for r in data["records"]["not_accepted"]] == ["rec:V:1", "rec:V:3"]
    assert data["alleles"][0]["anchor"] == {"value": 285, "support": "assertion",
                                            "derivation": "stated", "source": None}
    assert data["alleles"][0]["boundary_5"] is None
    assert data["defaults"]["profile"] == "genairr.unfitted_defaults.v1"
    assert data["capabilities"] == ref.capabilities().to_dict()
    assert data["findings"][0]["kind"] == "identical_sequence_distinct_labels"


def test_template_load_build_round_trip(tmp_path):
    paths = fasta_files(tmp_path)
    template = Reference.annotation_template(v=paths["v"], d=paths["d"], j=paths["j"], reference=header())
    assert [(a["segment"], a["label"]) for a in template["alleles"]] == [
        ("V", "IGHV1-2*02"), ("V", "IGHV1-3*01"), ("V", "IGHV1-2_copy"),
        ("D", "IGHD3-10*01"), ("D", "IGHD2-2*01"), ("J", "IGHJ4*02"), ("J", "IGHJ6*02")]
    assert all(a[k] is None for a in template["alleles"]
               for k in ("functionality", "anchor", "coding_phase", "boundary_5", "boundary_3"))
    assert {(s["allele"], s["value"]) for s in template["suggestions"]} == {("J:IGHJ4*02", 14),
                                                                            ("J:IGHJ6*02", 29)}
    assert template["conventions"] == {"coordinates": "zero_based_half_open", "v_alignment": "none"}
    annotations = tmp_path / "annotations.json"
    annotations.write_bytes(canonical_json(template).encode("utf-8"))
    ref = Reference.from_fasta(v=paths["v"], d=paths["d"], j=paths["j"], annotations=annotations)
    caps = ref.capabilities()
    assert caps.available() == ("catalogue_inspection", "observed_span_simulation")
    explain = caps.explain()
    assert "- J:IGHJ4*02: accepted anchor" in explain and "candidate at 14 (unique J motif)" in explain
    assert "- J:IGHJ6*02: accepted anchor" in explain and "candidate at 29 (unique J motif)" in explain
    assert {f["kind"] for f in ref.admission["accepted_curatable_findings"]} == {
        "MissingAnchor", "ReferenceEvidenceProblem"}
    assert all(a.anchor is None for a in ref.alleles)
    assert Reference.from_json(ref.to_json()).to_json() == ref.to_json()


def test_template_validates_reference_and_segments(tmp_path):
    paths = fasta_files(tmp_path)
    with pytest.raises(ReferenceArtifactError) as e:
        Reference.annotation_template(v=paths["v"], j=paths["j"], reference=header())
    assert e.value.kind == "malformed" and e.value.path == "d"
    with pytest.raises(ReferenceArtifactError) as e:
        Reference.annotation_template(v=paths["v"], d=paths["d"], j=paths["j"],
                                      reference=header(chain_type="vj"))
    assert e.value.kind == "malformed" and e.value.path == "d"
    bad = header()
    del bad["locus"]
    with pytest.raises(ReferenceArtifactError) as e:
        Reference.annotation_template(v=paths["v"], d=paths["d"], j=paths["j"], reference=bad)
    assert e.value.kind == "malformed" and e.value.path == "reference"
    with pytest.raises(ReferenceArtifactError) as e:
        Reference.annotation_template(v=paths["v"], d=paths["d"], j=paths["j"], reference=header(),
                                      v_alignment="gapped")
    assert e.value.kind == "malformed"
    gapped = Reference.annotation_template(v=paths["v"], d=paths["d"], j=paths["j"], reference=header(),
                                           v_alignment="imgt_gapped")
    assert gapped["conventions"]["v_alignment"] == "imgt_gapped"


def test_public_exports():
    for name in ("Reference", "ReferenceArtifactError", "ReferenceCapabilityError"):
        assert name in GenAIRR.__all__
        assert hasattr(GenAIRR, name)
    assert issubclass(ReferenceCapabilityError, ReferenceArtifactError)


# ── Regeneration of derived catalogue content ────────────────────────


def test_full_recommit_helper_is_faithful(ref):
    assert load(_full_recommit(doc_of(ref))).to_json() == ref.to_json()


def test_recommitted_finding_deleted_or_forged_rejected(ref):
    doc = doc_of(ref)
    doc["catalogue"]["findings"] = []
    err = rejected(_full_recommit(doc), "integrity", "identical_sequence_distinct_labels")
    assert err.path.startswith("catalogue.findings")
    doc = doc_of(ref)
    doc["catalogue"]["findings"].append({
        "kind": "identical_sequence_distinct_labels", "subjects": ["V:IGHV1-3*01", "J:IGHJ4*02"],
        "detail": "2 alleles share one sequence under different labels"})
    err = rejected(_full_recommit(doc), "integrity", "finding")
    assert err.path == "catalogue.findings[1]"


def test_recommitted_anchor_derived_phase_disagreeing_with_anchor_rejected(ref):
    doc = doc_of(ref)
    claim_of(doc, "V:IGHV1-2*02", "coding_phase", "anchor_derived")["value"] = 1
    err = rejected(_full_recommit(doc), "integrity", "anchor-derived coding phase 1")
    assert "285" in str(err) and "V:IGHV1-2*02" in str(err)


def _golden_with_stated_source_phases(disagreeing):
    """The golden artifact with each anchor-derived phase restated as a source-supported
    independent phase, plus (if ``disagreeing``) a second such claim of another value."""
    doc = json.loads(GOLDEN.read_text(encoding="utf-8"))
    claims = []
    for c in doc["catalogue"]["claims"]:
        if c["predicate"] != "coding_phase":
            claims.append(c)
            continue
        stated = {**c, "derivation": "stated", "support": "source", "source": SOURCE}
        claims.append(stated)
        if disagreeing:
            claims.append({**stated, "value": (c["value"] + 1) % 3})
    doc["catalogue"]["claims"] = claims
    return doc


def test_recommitted_conflicting_phases_are_not_source_phase_support():
    # Every V and J carries two stated, source-supported phases that disagree: the report
    # lists five phase conflicts, so no allele has a usable phase.
    reference = load(_full_recommit(_golden_with_stated_source_phases(disagreeing=True)))
    v_and_j = [f"{a.segment}:{a.label}" for a in reference.alleles if a.segment != "D"]
    assert len(v_and_j) == 5
    assert reference.report().conflicts == tuple((a, "coding_phase") for a in v_and_j)
    assert all(a.coding_phase is None for a in reference.alleles)
    caps = reference.capabilities()
    assert caps.phase_support == "none"
    assert caps.status("canonical_productive").note == "phase support: none"
    assert "phase support: none" in str(reference.report())
    # The same phases without the disagreeing copies are independent source support.
    agreeing = load(_full_recommit(_golden_with_stated_source_phases(disagreeing=False)))
    assert agreeing.report().conflicts == ()
    assert agreeing.capabilities().phase_support == "independent_source"
    # Its anchors are assertions, so canonical productivity stays assertion-level.
    assert agreeing.capabilities().status("canonical_productive").evidence == "assertion"


def test_recommitted_motif_candidate_moved_rejected(ref):
    doc = doc_of(ref)
    claim_of(doc, "J:IGHJ4*02", "anchor", "unique_motif")["value"] = 15
    err = rejected(_full_recommit(doc), "integrity", r"J:IGHJ4\*02")
    assert err.path.startswith("catalogue.claims[")
    doc = doc_of(ref)
    doc["catalogue"]["claims"] = [c for c in doc["catalogue"]["claims"] if c["derivation"] != "unique_motif"]
    rejected(_full_recommit(doc), "integrity", "claim")


def _gapped(tmp_path, relationships=()):
    ungapped = GAPPED_V.replace(".", "")
    alleles = [entry("V", "gv", ungapped, anchor=ann(306, "source", SOURCE)),
               entry("V", "gw", ungapped[:-3] + "TTT", anchor=ann(303, "source", SOURCE))] + \
        default_alleles()[3:]
    v = f">gv\n{GAPPED_V}\n>gw\n{GAPPED_V[:-3]}TTT\n"
    return build(tmp_path, alleles=alleles, v=v, v_alignment="imgt_gapped", relationships=relationships)


def test_recommitted_alignment_claims_and_conflict_statuses_regenerated(tmp_path):
    ref = _gapped(tmp_path)
    doc = doc_of(ref)
    for c in doc["catalogue"]["claims"]:
        if c["subject"] == "V:gw" and c["predicate"] == "anchor":
            assert c["status"] == "conflict"
            c["status"] = "accepted"
    rejected(_full_recommit(doc), "integrity", "V:gw")
    doc = doc_of(ref)
    doc["catalogue"]["claims"].remove(claim_of(doc, "V:gv", "v_subregions", "imgt_alignment"))
    rejected(_full_recommit(doc), "integrity", "claim")
    doc = doc_of(ref)
    claim_of(doc, "V:gv", "v_subregions", "imgt_alignment")["value"]["FWR1"] = [0, 60]
    rejected(_full_recommit(doc), "integrity", "V:gv")


def test_recommitted_entities_relationships_and_genes_regenerated(tmp_path):
    rel = {"kind": "alias", "allele": "V:IGHV1-2*02", "entity": "lab_alias", "source": None}
    ref = build(tmp_path, relationships=[rel])
    assert [(e.id, e.name) for e in ref.entities] == [("entity:lab_alias", "lab_alias")]
    doc = doc_of(ref)
    doc["catalogue"]["entities"].append({"id": "entity:orphan", "name": "orphan", "executable": False})
    assert rejected(_full_recommit(doc), "integrity", "orphan").path.startswith("catalogue.entities")
    doc = doc_of(ref)
    doc["catalogue"]["entities"][0]["name"] = "renamed"
    rejected(_full_recommit(doc), "integrity", "entit")
    doc = doc_of(ref)
    doc["catalogue"]["entities"][0]["id"] = "ent:lab_alias"
    doc["catalogue"]["relationships"][0]["entity"] = "ent:lab_alias"
    rejected(_full_recommit(doc), "integrity", "relationship")
    doc = doc_of(ref)
    doc["catalogue"]["alleles"][0]["gene"] = "IGHV9"
    err = rejected(_full_recommit(doc), "integrity", "gene")
    assert err.path == "catalogue.alleles[0].gene"



def test_paralog_between_executable_alleles_links_alleles_without_an_entity(tmp_path):
    from GenAIRR._reference.importer import regenerate

    rel = {"kind": "paralog", "allele": "V:IGHV1-2*02", "entity": "IGHV1-3*01", "source": None}
    ref = build(tmp_path, relationships=[rel])
    assert ref.entities == ()
    assert {"V:IGHV1-2*02", "V:IGHV1-3*01"} <= {f"{a.segment}:{a.label}" for a in ref.alleles}
    report = ref.report()
    assert [r.to_dict() for r in report.relationships] == [
        {"kind": "paralog", "allele": "V:IGHV1-2*02", "entity": "V:IGHV1-3*01", "source": None}]
    lines = str(report).splitlines()
    assert "  V:IGHV1-2*02 paralog of V:IGHV1-3*01" in lines
    entities_at = lines.index("Entities (named, not executable)")
    assert lines[entities_at + 1] == "  (none)"
    assert not any(line.startswith("  IGHV1-3*01:") for line in lines)
    again = Reference.from_json(ref.to_json())
    assert again.to_json() == ref.to_json()
    derived = regenerate(ref._catalogue)
    assert (derived.relationships, derived.entities) == (ref._catalogue.relationships, ())
    # The earlier representation (a phantom entity named like the allele) does not load.
    doc = doc_of(ref)
    doc["catalogue"]["entities"] = [{"id": "entity:IGHV1-3*01", "name": "IGHV1-3*01", "executable": False}]
    doc["catalogue"]["relationships"][0]["entity"] = "entity:IGHV1-3*01"
    err = rejected(_full_recommit(doc), "integrity", "relationship")
    assert err.path.startswith("catalogue.relationships")


def test_recommitted_allele_with_two_accepted_records_rejected(ref):
    doc = doc_of(ref)
    records = doc["catalogue"]["records"]
    assert (records[1]["status"], records[1]["duplicate_of"]) == ("duplicate", "rec:V:0")
    records[1].update(status="accepted", duplicate_of=None)
    err = rejected(_full_recommit(doc), "integrity", "exactly one accepted record")
    assert err.path == "catalogue.alleles[0].records"


def test_recommitted_accepted_record_with_another_label_rejected():
    # An accepted record names its executable allele: it carries that allele's label.
    doc = json.loads(GOLDEN.read_text(encoding="utf-8"))
    record = doc["catalogue"]["records"][0]
    assert (record["status"], record["allele"], record["label"]) == ("accepted", "V:IGHV1-2*02", "IGHV1-2*02")
    record["label"] = "invented_alias"
    err = rejected(_full_recommit(doc), "integrity", "invented_alias")
    assert err.path == "catalogue.records[0].label"
    assert "IGHV1-2*02" in err.message


def test_recommitted_duplicate_with_another_label_rejected(ref):
    doc = doc_of(ref)
    doc["catalogue"]["records"][1]["label"] = "IGHV1-2_other"
    err = rejected(_full_recommit(doc), "integrity", "label")
    assert err.path == "catalogue.records[1].label"


@pytest.mark.parametrize("file_name", ["lab/v.fasta", "..\\v.fasta", "C:v.fasta", "/abs/v.fasta",
                                       "\\\\server\\share\\v.fasta", "..", ".", ""])
def test_recommitted_source_path_rejected(ref, file_name):
    doc = doc_of(ref)
    doc["catalogue"]["sources"][0]["file_name"] = file_name
    err = rejected(_full_recommit(doc), "integrity", "bare file name")
    assert err.path == "catalogue.sources[0].file_name"


def test_out_of_bounds_sidecar_anchor_is_located(tmp_path):
    alleles = default_alleles()
    alleles[1]["anchor"] = ann(900)
    with pytest.raises(ReferenceArtifactError) as e:
        build(tmp_path, alleles=alleles)
    assert e.value.kind == "conflict"
    assert e.value.path == "alleles[1].anchor"
    assert all(text in e.value.message for text in ("V:IGHV1-3*01", "900", f"length {len(V2)}")), \
        e.value.message
    alleles[1]["anchor"] = ann(len(V2) - 2)     # the codon would run past the end
    with pytest.raises(ReferenceArtifactError) as e:
        build(tmp_path / "edge", alleles=alleles)
    assert (e.value.kind, e.value.path) == ("conflict", "alleles[1].anchor")
    alleles[1]["anchor"] = ann(len(V2) - 3)     # the last codon fits
    assert build(tmp_path / "fits", alleles=alleles).alleles[1].anchor.value == len(V2) - 3


def test_allele_notes_round_trip_and_reach_the_report(tmp_path):
    alleles = default_alleles()
    alleles[0]["notes"] = ["second 5' UTR variant seen in lab clone 17", "re-sequenced 2026"]
    alleles[4]["notes"] = ["low coverage"]
    ref = build(tmp_path / "build", alleles=alleles)
    plain = build(tmp_path / "plain")
    assert ref.catalogue_hash != plain.catalogue_hash      # notes are part of the catalogue
    stored = {a["id"]: a["notes"] for a in doc_of(ref)["catalogue"]["alleles"]}
    assert stored["V:IGHV1-2*02"] == alleles[0]["notes"]
    assert stored["J:IGHJ6*02"] == ["low coverage"]
    assert stored["V:IGHV1-3*01"] == []
    path = tmp_path / "notes.genairr-ref.json"
    ref.write(path)
    again = Reference.read(path)
    assert again.to_json() == ref.to_json()
    assert [a.notes for a in again.alleles] == [tuple(alleles[0]["notes"]), (), (), (), (), (), ("low coverage",)]
    text = str(again.report())
    lines = text.splitlines()
    notes_at = lines.index("Notes")
    assert lines[notes_at + 1:notes_at + 4] == [
        "  V:IGHV1-2*02: second 5' UTR variant seen in lab clone 17",
        "  V:IGHV1-2*02: re-sequenced 2026",
        "  J:IGHJ6*02: low coverage"]
    assert again.report().to_dict()["alleles"][0]["notes"] == alleles[0]["notes"]
    assert "Notes" not in str(plain.report()).splitlines()
    # Notes are base data: a recommitted edit changes the catalogue rather than failing regeneration.
    doc = doc_of(ref)
    doc["catalogue"]["alleles"][0]["notes"] = ["edited"]
    assert load(_full_recommit(doc)).alleles[0].notes == ("edited",)


def test_recommitted_malformed_notes_rejected(ref):
    for bad, match in ((None, "expected an array"), ([1], "expected a string"), (["x" * 4097], "4096")):
        doc = doc_of(ref)
        doc["catalogue"]["alleles"][0]["notes"] = bad
        err = rejected(_recommit(doc), "malformed", match)
        assert err.path.startswith("catalogue.alleles[0].notes")
    doc = doc_of(ref)
    del doc["catalogue"]["alleles"][0]["notes"]
    rejected(_recommit(doc), "malformed", "missing keys")


# ── Empty selections ─────────────────────────────────────────────────

EMPTY = {"V": f">v_n\n{V_BAD}\n>v_e\n\n>v,comma\n{V1}\n>v_n2\nACGN\n",
         "D": ">d_n\nGGNTT\n",
         "J": f">j_e\n\n>j_n\n{J1[:-1]}N\n"}


@pytest.mark.parametrize("segment, counts, reasons", [
    ("V", "0 accepted, 0 duplicate, 4 rejected",
     ["unsupported_base (2)", "empty_sequence (1)", "unrepresentable_label (1)"]),
    ("D", "0 accepted, 0 duplicate, 1 rejected", ["unsupported_base (1)"]),
    ("J", "0 accepted, 0 duplicate, 2 rejected", ["empty_sequence (1)", "unsupported_base (1)"]),
])
def test_empty_selection_names_segment_and_accounting(tmp_path, monkeypatch, segment, counts, reasons):
    fastas = {"v": V_FASTA, "d": D_FASTA, "j": J_FASTA, segment.lower(): EMPTY[segment]}
    alleles = [a for a in default_alleles() if a["segment"] != segment]

    def no_lowering(*args, **kwargs):
        raise AssertionError("native lowering ran before the selection check")

    monkeypatch.setattr(artifact_module, "lower", no_lowering)
    with pytest.raises(ReferenceArtifactError) as e:
        build(tmp_path, alleles=alleles, **fastas)
    assert e.value.kind == "empty_selection" and e.value.path == segment
    message = e.value.message
    assert f"no executable {segment} alleles" in message and counts in message
    assert all(r in message for r in reasons), message
    positions = [message.index(r) for r in reasons]
    assert positions == sorted(positions)


def test_empty_selection_counts_same_label_conflicts(tmp_path):
    fastas = {"v": f">IGHV1-2*02\n{V_BAD}\n>IGHV1-2*02\n{V1}\n", "d": D_FASTA, "j": J_FASTA}
    alleles = [a for a in default_alleles() if a["segment"] != "V"]
    with pytest.raises(ReferenceArtifactError) as e:
        build(tmp_path, alleles=alleles, **fastas)
    assert e.value.kind == "empty_selection" and e.value.path == "V"
    assert "0 accepted, 0 duplicate, 2 rejected" in e.value.message
    assert "rejection reasons: conflicting_sequence_for_label (2)" in e.value.message


# ── Untampered round trips with a declared IMGT alignment ───────────


def _subregion_order(reference):
    native = reference.native_refdata()
    return [[s[0] for s in native.v_allele(i).subregions] for i in range(native.v_pool_size())]


def test_gapped_reference_with_subregions_round_trips_byte_identical(tmp_path):
    ref = _gapped(tmp_path / "build")
    full = ["FWR1", "CDR1", "FWR2", "CDR2", "FWR3"]
    assert _subregion_order(ref) == [full, full]
    again = Reference.from_json(ref.to_json())
    assert again.to_json() == ref.to_json()
    assert _subregion_order(again) == [full, full]
    path = tmp_path / "gapped.genairr-ref.json"
    ref.write(path)
    read = Reference.read(path)
    for name in ("catalogue_hash", "models_hash", "execution_hash", "artifact_hash"):
        assert getattr(read, name) == getattr(ref, name), name
    read.write(tmp_path / "copy.genairr-ref.json")
    assert (tmp_path / "copy.genairr-ref.json").read_bytes() == path.read_bytes() == \
        ref.to_json().encode("utf-8")
    assert _subregion_order(read) == [full, full]
