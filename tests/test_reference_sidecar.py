import copy
import hashlib
import json

import pytest

from GenAIRR._artifact_io import canonical_json
from GenAIRR._reference.errors import ReferenceArtifactError, ReferenceCapabilityError
from GenAIRR._reference.model import Catalogue, check_label
from GenAIRR._reference.sidecar import annotation_template, load_sidecar

V_SEQ = "ACGT" * 8
J_SEQ = "TGGGGCCAAGGGACC"


def digest(seq):
    return "sha256:" + hashlib.sha256(seq.encode("ascii")).hexdigest()


def minimal_sidecar():
    return {
        "schema": "genairr_reference_annotations",
        "schema_version": 1,
        "reference": {
            "name": "lab_human_igh", "version": "1", "locus": "IGH", "chain_type": "vdj",
            "taxon": {"scientific_name": "Homo sapiens", "taxonomy_id": 9606},
            "build_date": "2026-10-21",
        },
        "conventions": {"coordinates": "zero_based_half_open", "v_alignment": "none"},
        "alleles": [{
            "segment": "V", "label": "IGHV1-2*02", "sequence_sha256": digest(V_SEQ), "gene": None,
            "functionality": {"value": "functional", "support": "assertion", "source": None, "method": "m"},
            "anchor": {"value": 5, "support": "source", "source": "doi:x", "method": "m"},
            "coding_phase": {"value": 0, "derivation": "anchor_derived", "support": "assertion",
                             "source": None, "method": "m"},
            "boundary_5": None,
            "boundary_3": {"value": "complete", "support": "assertion", "source": None, "method": "m"},
            "notes": [],
        }],
        "relationships": [{"kind": "alias", "allele": "V:IGHV1-2*02", "entity": "lab_v1", "source": None}],
        "suggestions": [{"allele": "J:IGHJ4*02", "predicate": "anchor", "value": 14,
                         "derivation": "unique_motif", "note": "candidate"}],
    }


def load(d):
    return load_sidecar(json.dumps(d))


def test_minimal_sidecar_loads():
    out = load(minimal_sidecar())
    assert out["alleles"][0]["label"] == "IGHV1-2*02"
    d = minimal_sidecar()
    del d["relationships"], d["suggestions"]
    out = load(d)
    assert out["relationships"] == [] and out["suggestions"] == []


def test_sidecar_with_bom_loads():
    assert load_sidecar("﻿" + json.dumps(minimal_sidecar()))["schema_version"] == 1


def test_unknown_key_rejected_at_every_level():
    def at(path):
        d = minimal_sidecar()
        node = d
        for p in path:
            node = node[p]
        node["extra"] = 1
        return d
    for path in [(), ("reference",), ("reference", "taxon"), ("conventions",), ("alleles", 0),
                 ("alleles", 0, "anchor"), ("relationships", 0), ("suggestions", 0)]:
        with pytest.raises(ReferenceArtifactError) as e:
            load(at(path))
        assert e.value.kind == "malformed", path


def test_source_support_requires_source():
    for blank in (None, "", "  "):
        d = minimal_sidecar()
        d["alleles"][0]["anchor"]["source"] = blank
        with pytest.raises(ReferenceArtifactError) as e:
            load(d)
        assert e.value.kind == "malformed"


def test_bare_string_boundary_rejected():
    d = minimal_sidecar()
    d["alleles"][0]["boundary_3"] = "complete"
    with pytest.raises(ReferenceArtifactError) as e:
        load(d)
    assert e.value.kind == "malformed"


def test_unsupported_coordinate_convention():
    d = minimal_sidecar()
    d["conventions"]["coordinates"] = "one_based_closed"
    with pytest.raises(ReferenceArtifactError) as e:
        load(d)
    assert e.value.kind == "unsupported_convention"


def test_schema_version_2_unsupported():
    d = minimal_sidecar()
    d["schema_version"] = 2
    with pytest.raises(ReferenceArtifactError) as e:
        load(d)
    assert e.value.kind == "unsupported_schema"


def test_duplicate_allele_entry_rejected():
    d = minimal_sidecar()
    d["alleles"].append(copy.deepcopy(d["alleles"][0]))
    with pytest.raises(ReferenceArtifactError) as e:
        load(d)
    assert e.value.kind == "malformed"


def test_malformed_json_and_values():
    with pytest.raises(ReferenceArtifactError) as e:
        load_sidecar('{"a":1,"a":2}')
    assert e.value.kind == "malformed"
    for mutate in (
        lambda d: d["alleles"][0]["anchor"].__setitem__("value", -1),
        lambda d: d["alleles"][0]["coding_phase"].__setitem__("value", 3),
        lambda d: d["alleles"][0]["functionality"].__setitem__("value", "maybe"),
        lambda d: d["alleles"][0]["boundary_3"].__setitem__("value", "full"),
        lambda d: d["alleles"][0]["anchor"].__setitem__("value", True),
        lambda d: d["alleles"][0].__setitem__("label", "a,b"),
        lambda d: d["alleles"][0].__setitem__("sequence_sha256", "abc"),
        lambda d: d["reference"].__setitem__("chain_type", "vdj2"),
        lambda d: d["reference"].__setitem__("build_date", "2026-13-40"),
        lambda d: d["conventions"].__setitem__("v_alignment", "other"),
        lambda d: d["alleles"][0]["anchor"].__setitem__("derivation", "independent"),
        lambda d: d["alleles"][0]["coding_phase"].__setitem__("derivation", "stated"),
    ):
        d = minimal_sidecar()
        mutate(d)
        with pytest.raises(ReferenceArtifactError) as e:
            load(d)
        assert e.value.kind == "malformed"


def test_suggestions_shape_validated():
    for mutate in (
        lambda s: s.pop("note"),
        lambda s: s.__setitem__("predicate", "nope"),
        lambda s: s.__setitem__("derivation", "independent"),
        lambda s: s.__setitem__("allele", "X:foo"),
        lambda s: s.__setitem__("value", "complete"),
    ):
        d = minimal_sidecar()
        mutate(d["suggestions"][0])
        with pytest.raises(ReferenceArtifactError) as e:
            load(d)
        assert e.value.kind == "malformed"
    d = minimal_sidecar()
    d["suggestions"] = {}
    with pytest.raises(ReferenceArtifactError):
        load(d)


def test_template_loads_and_prefills_nothing():
    s = minimal_sidecar()
    template = annotation_template(
        records_by_segment={
            "V": [{"label": "IGHV1-2*02", "sequence_sha256": digest(V_SEQ)}],
            "D": [],
            "J": [{"label": "IGHJ4*02", "sequence_sha256": digest(J_SEQ)}],
        },
        header=s["reference"], conventions=s["conventions"], suggestions=s["suggestions"],
    )
    loaded = load_sidecar(canonical_json(template))
    for a in loaded["alleles"]:
        for k in ("functionality", "anchor", "coding_phase", "boundary_5", "boundary_3", "gene"):
            assert a[k] is None
        assert a["notes"] == []
    stripped = dict(template, suggestions=[])
    text = canonical_json(stripped)
    for word in ('"source"', '"complete"', '"independent"'):
        assert word not in text
    assert [a["segment"] for a in template["alleles"]] == ["V", "J"]


def test_check_label():
    assert check_label("kv_a") is None
    assert check_label("IGHV1-2*02") is None
    assert check_label("a,b") == "contains ',' at position 1"
    assert "whitespace" in check_label("a b")
    assert "whitespace" in check_label("a\tb")
    assert "format" in check_label("a​b") and "position 1" in check_label("a​b")
    assert "control" in check_label("a\x00b")
    assert check_label("") is not None
    assert check_label("x" * 256) is None
    assert check_label("x" * 257) is not None


def test_error_types():
    e = ReferenceCapabilityError("capability", "p", "m")
    assert isinstance(e, ReferenceArtifactError) and isinstance(e, ValueError)
    assert (e.kind, e.path) == ("capability", "p")


def catalogue_dict():
    return {
        "header": {"name": "r", "version": "1",
                   "taxon": {"scientific_name": "Homo sapiens", "taxonomy_id": 9606},
                   "locus": "IGK", "chain_type": "vj", "build_date": None},
        "conventions": {"coordinates": "zero_based_half_open", "v_alignment": "imgt_gapped"},
        "sources": [{"id": "src:V", "role": "V", "file_name": "v.fasta", "sha256": digest("x"), "bytes": 10},
                    {"id": "src:J", "role": "J", "file_name": "j.fasta", "sha256": digest("y"), "bytes": 10}],
        "records": [
            {"id": "rec:V:0", "source": "src:V", "index": 0, "header": "kv_a", "label": "kv_a",
             "status": "accepted", "reason": None, "allele": "V:kv_a", "duplicate_of": None},
            {"id": "rec:V:1", "source": "src:V", "index": 1, "header": "kv_a", "label": "kv_a",
             "status": "duplicate", "reason": None, "allele": None, "duplicate_of": "rec:V:0"},
            {"id": "rec:V:2", "source": "src:V", "index": 2, "header": "bad", "label": None,
             "status": "rejected", "reason": "bad char", "allele": None, "duplicate_of": None},
            {"id": "rec:J:0", "source": "src:J", "index": 0, "header": "kj_a", "label": "kj_a",
             "status": "accepted", "reason": None, "allele": "J:kj_a", "duplicate_of": None},
        ],
        "alleles": [
            {"id": "V:kv_a", "segment": "V", "label": "kv_a", "gene": "kv", "gene_basis": "label_prefix",
             "sequence": V_SEQ, "sequence_sha256": digest(V_SEQ),
             "gapped_sequence": "ACGT..ACGT" + V_SEQ[8:],
             "transforms": [], "records": ["rec:V:0"]},
            {"id": "J:kj_a", "segment": "J", "label": "kj_a", "gene": "kj", "gene_basis": "label_prefix",
             "sequence": J_SEQ, "sequence_sha256": digest(J_SEQ), "gapped_sequence": None,
             "transforms": [{"kind": "uppercase", "input_sha256": digest("acgt"),
                             "output_sha256": digest(J_SEQ)}],
             "records": ["rec:J:0"]},
        ],
        "claims": [
            {"subject": "V:kv_a", "predicate": "anchor", "value": 6, "support": "source",
             "derivation": "stated", "status": "accepted", "source": "doi:x", "method": "m"},
            {"subject": "V:kv_a", "predicate": "coding_phase", "value": 0, "support": "assertion",
             "derivation": "anchor_derived", "status": "accepted", "source": None, "method": None},
            {"subject": "V:kv_a", "predicate": "v_subregions",
             "value": {"FWR1": [0, 10], "CDR1": [10, 15]}, "support": "assertion",
             "derivation": "imgt_alignment", "status": "accepted", "source": None, "method": None},
            {"subject": "J:kj_a", "predicate": "anchor", "value": 3, "support": "none",
             "derivation": "unique_motif", "status": "candidate", "source": None, "method": None},
        ],
        "entities": [{"id": "ent:lab_v1", "name": "lab_v1", "executable": False}],
        "relationships": [{"kind": "alias", "allele": "V:kv_a", "entity": "ent:lab_v1", "source": None}],
        "findings": [{"kind": "identical_sequence_distinct_labels", "subjects": ["V:kv_a"], "detail": "d"}],
    }


def test_catalogue_round_trip_and_digest():
    d = catalogue_dict()
    cat = Catalogue.from_dict(d)
    assert cat.to_dict() == json.loads(json.dumps(d))
    assert Catalogue.from_dict(cat.to_dict()).digest() == cat.digest()
    assert cat.digest().startswith("sha256:")
    assert cat.allele("V:kv_a").label == "kv_a"
    assert len(cat.claims_for("V:kv_a", "anchor")) == 1
    assert cat.claims_for("V:kv_a", "boundary_3") == ()
    d2 = catalogue_dict()
    d2["claims"][0]["value"] = 7
    assert Catalogue.from_dict(d2).digest() != cat.digest()


def _set(path, value):
    def f(d):
        node = d
        for p in path[:-1]:
            node = node[p]
        node[path[-1]] = value
    return f


def _pop(path):
    def f(d):
        node = d
        for p in path[:-1]:
            node = node[p]
        del node[path[-1]]
    return f


CASES = {
    "claim_missing_allele": _set(("claims", 0, "subject"), "V:nope"),
    "anchor_out_of_bounds": _set(("claims", 0, "value"), len(V_SEQ) - 2),
    "anchor_negative": _set(("claims", 0, "value"), -1),
    "source_without_source": _set(("claims", 0, "source"), " "),
    "digest_mismatch": _set(("alleles", 0, "sequence_sha256"), digest("ACGT")),
    "alias_targets_executable": lambda d: (
        d["alleles"].append(dict(d["alleles"][0], id="V:kv_b", label="kv_b", records=[])),
        d["entities"].__setitem__(0, {"id": "ent:lab_v1", "name": "kv_a", "executable": False}),
        d["relationships"].__setitem__(0, {"kind": "alias", "allele": "V:kv_b",
                                           "entity": "ent:lab_v1", "source": None})),
    "d_in_vj": lambda d: (
        d["alleles"].append(dict(d["alleles"][0], id="D:kd", segment="D", label="kd", records=[]))),
    "anchor_derived_on_boundary": lambda d: (
        d["claims"].append({"subject": "V:kv_a", "predicate": "boundary_3", "value": "complete",
                            "support": "assertion", "derivation": "anchor_derived", "status": "accepted",
                            "source": None, "method": None})),
    "candidate_with_support": _set(("claims", 3, "support"), "assertion"),
    "none_support_accepted": _set(("claims", 3, "status"), "accepted"),
    "bad_phase": _set(("claims", 1, "value"), 3),
    "overlapping_subregions": _set(("claims", 2, "value"), {"FWR1": [0, 10], "CDR1": [5, 15]}),
    "subregion_out_of_bounds": _set(("claims", 2, "value"), {"FWR1": [0, 99]}),
    "bad_subregion_label": _set(("claims", 2, "value"), {"XYZ": [0, 3]}),
    "bad_sequence_chars": lambda d: d["alleles"][0].update(
        sequence="ACGN" * 8, sequence_sha256=digest("ACGN" * 8)),
    "id_mismatch": _set(("alleles", 0, "id"), "V:other"),
    "duplicate_allele_id": lambda d: d["alleles"].append(copy.deepcopy(d["alleles"][0])),
    "bad_label": lambda d: d["alleles"][0].update(label="a,b", id="V:a,b"),
    "transform_output_mismatch": _set(("alleles", 1, "transforms", 0, "output_sha256"), digest("zzz")),
    "gapped_without_imgt": _set(("conventions", "v_alignment"), "none"),
    "accepted_record_missing_allele": _set(("records", 0, "allele"), "V:none"),
    "duplicate_of_not_accepted": _set(("records", 1, "duplicate_of"), "rec:V:2"),
    "rejected_without_reason": _set(("records", 2, "reason"), None),
    "rejected_with_allele": _set(("records", 2, "allele"), "V:kv_a"),
    "executable_entity": _set(("entities", 0, "executable"), True),
    "relationship_missing_allele": _set(("relationships", 0, "allele"), "V:none"),
    "float_phase": _set(("claims", 1, "value"), 1.0),
    "float_anchor": _set(("claims", 0, "value"), 6.0),
    "float_subregion": _set(("claims", 2, "value"), {"FWR1": [0.0, 10]}),
    "float_taxonomy_id": _set(("header", "taxon", "taxonomy_id"), 9606.0),
    "float_record_index": _set(("records", 0, "index"), 0.0),
    "float_source_bytes": _set(("sources", 0, "bytes"), 10.0),
    "allele_does_not_list_record": _set(("alleles", 0, "records"), []),
    "allele_lists_foreign_record": _set(("alleles", 0, "records"), ["rec:V:0", "rec:J:0"]),
    "record_segment_differs_from_allele": lambda d: (
        d["records"][0].__setitem__("id", "rec:J:0"), d["records"][3].__setitem__("id", "rec:V:3"),
        d["records"][3].__setitem__("index", 3), d["records"][0].__setitem__("index", 0),
        d["alleles"][0].__setitem__("records", ["rec:J:0"]),
        d["alleles"][1].__setitem__("records", ["rec:V:3"])),
    "finding_unknown_subject": _set(("findings", 0, "subjects"), ["V:nope"]),
    "unknown_key": _set(("alleles", 0, "extra"), 1),
    "missing_key": _pop(("claims", 0, "method")),
    "bad_chain_type": _set(("header", "chain_type"), "xx"),
    "bad_functionality": lambda d: d["claims"].append({
        "subject": "V:kv_a", "predicate": "functionality", "value": "meh", "support": "assertion",
        "derivation": "stated", "status": "accepted", "source": None, "method": None}),
}


@pytest.mark.parametrize("name", sorted(CASES))
def test_catalogue_integrity_rejections(name):
    d = catalogue_dict()
    CASES[name](d)
    with pytest.raises(ReferenceArtifactError) as e:
        Catalogue.from_dict(d)
    assert e.value.kind in ("integrity", "malformed")


def test_gapped_sequence_accepted_under_imgt_gapped():
    d = catalogue_dict()
    d["conventions"]["v_alignment"] = "imgt_gapped"
    d["alleles"][0]["gapped_sequence"] = "ACGT..ACGT" + V_SEQ[8:]
    # dots removed gives ACGTACGT + rest == V_SEQ
    Catalogue.from_dict(d)
    d["alleles"][0]["gapped_sequence"] = "ACGT..ACGA" + V_SEQ[8:]
    with pytest.raises(ReferenceArtifactError):
        Catalogue.from_dict(d)


def test_float_integers_rejected_in_sidecar():
    for name, value in (("coding_phase", 1.0), ("coding_phase", 2.0), ("anchor", 285.0)):
        d = minimal_sidecar()
        d["alleles"][0][name]["value"] = value
        with pytest.raises(ReferenceArtifactError) as e:
            load(d)
        assert e.value.kind == "malformed"
    d = minimal_sidecar()
    d["reference"]["taxon"]["taxonomy_id"] = 9606.0
    with pytest.raises(ReferenceArtifactError):
        load(d)
    d = minimal_sidecar()
    d["schema_version"] = 1.0
    with pytest.raises(ReferenceArtifactError):
        load(d)
