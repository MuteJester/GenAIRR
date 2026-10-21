import hashlib
import json

import pytest

from GenAIRR import data
from GenAIRR._refdata_resolver import dataconfig_to_refdata
from GenAIRR._reference import lowering
from GenAIRR._reference.errors import ReferenceArtifactError
from GenAIRR._reference.lowering import PROFILE_ID, default_view, lower, seal
from GenAIRR._reference.model import (
    Allele, Catalogue, Claim, RecordEntry, ReferenceHeader, SourceFile, Taxon,
)

V_SEQ = "AAATGTAAAGGG"   # TGT at 3
J_SEQ = "TGGGGCCAAGGG"   # TGG at 0
D_SEQ = "GGGTTTAAA"
CTX = {"schema_version": 1, "catalogue_hash": "sha256:" + "a" * 64, "view": "default",
       "profile": PROFILE_ID, "models_hash": "sha256:" + "b" * 64}


def _allele(seg, label, seq):
    return Allele(f"{seg}:{label}", seg, label, label, "declared", seq,
                  "sha256:" + hashlib.sha256(seq.encode("ascii")).hexdigest(), None, (), ())


def claim(subject, predicate, value, support="source", derivation="stated", status="accepted"):
    src = "doi:x" if support == "source" else None
    return Claim(subject, predicate, value, support, derivation, status, src, "m")


def build(claims=(), locus="IGH", chain="vdj", version="1", labels=("V1", "J1"), alignment="none",
          extra_v=()):
    specs = [("V", labels[0], V_SEQ)] + [("V", x, V_SEQ[:-1] + "C") for x in extra_v]
    if chain == "vdj":
        specs.append(("D", "D1", D_SEQ))
    specs.append(("J", labels[1], J_SEQ))
    counts = {}
    alleles, records = [], []
    for seg, label, seq in specs:
        i = counts[seg] = counts.get(seg, -1) + 1
        rid = f"rec:{seg}:{i}"
        alleles.append(Allele(f"{seg}:{label}", seg, label, label, "declared", seq,
                              "sha256:" + hashlib.sha256(seq.encode("ascii")).hexdigest(), None, (), (rid,)))
        records.append(RecordEntry(rid, "src", i, ">" + label, label, "accepted", None, f"{seg}:{label}", None))
    cat = Catalogue(
        ReferenceHeader("t", version, Taxon("Mus musculus", 10090), locus, chain, None),
        {"coordinates": "zero_based_half_open", "v_alignment": alignment},
        (SourceFile("src", "v", "v.fasta", "sha256:" + "0" * 64, 1),), tuple(records), tuple(alleles),
        tuple(claims), (), (), ())
    return Catalogue.from_dict(cat.to_dict())


def anchors(support="source"):
    return [claim("V:V1", "anchor", 3, support), claim("J:J1", "anchor", 0, support)]


def evidence(low, seg, i=0):
    return json.loads(low.refdata.allele_evidence_json(seg, i))


def test_lowered_pool_order_and_names_match_view():
    cat = build(labels=("VB", "VA"))
    view = default_view(cat)
    assert view == {"V": ["V:VB"], "D": ["D:D1"], "J": ["J:VA"]}
    low = lower(cat, view)
    spec = low.spec
    assert [a["name"] for a in spec["v"]] == ["VB"]
    assert [a["name"] for a in spec["j"]] == ["VA"]
    assert [a["name"] for a in spec["d"]] == ["D1"]
    assert not low.refdata.is_sealed()
    assert spec["rules"]["v_anchor"]["expected_amino_acids"] == ["C"]
    assert spec["rules"]["j_anchor"]["expected_amino_acids"] == ["W"]


def test_default_view_vj_has_no_d_and_keeps_record_order():
    cat = build(chain="vj", locus="IGK", extra_v=("A0",))
    assert default_view(cat) == {"V": ["V:V1", "V:A0"], "J": ["J:J1"]}


def test_identity_uses_free_form_taxon():
    ident = lower(build(), default_view(build())).spec["identity"]
    assert ident == {"species": "Mus musculus", "locus": "IGH", "reference_set": "t@1",
                     "name": "t", "source": "genairr_reference_v1"}
    cat = build(version=None)
    assert lower(cat, default_view(cat)).spec["identity"]["reference_set"] == "t"


def test_unknown_claims_stay_unknown_natively():
    low = lower(build(), default_view(build()))
    for seg in ("V", "J"):
        ev = evidence(low, seg)
        assert ev["anchor_status"] == "unknown" and ev["coding_phase"] is None
        assert ev["boundary_5"] == ev["boundary_3"] == "unknown"
        assert ev["source"] is None and ev["method"] == "genairr_reference_v1:none"
        assert ev["notes"] == []
    assert low.spec["v"][0]["anchor"] is None


def test_assertion_anchor_is_inferred_and_source_anchor_verified():
    cat = build([claim("V:V1", "anchor", 3, "source"), claim("J:J1", "anchor", 0, "assertion")])
    low = lower(cat, default_view(cat))
    v, j = evidence(low, "V"), evidence(low, "J")
    assert v["anchor_status"] == "verified" and v["source"] == "doi:x"
    assert v["method"] == "genairr_reference_v1:source:stated"
    assert j["anchor_status"] == "inferred" and j["source"] is None
    assert j["method"] == "genairr_reference_v1:assertion:stated"
    assert low.spec["v"][0]["anchor"] == 3


def test_anchor_derived_phase_not_projected():
    cat = build(anchors() + [claim("V:V1", "coding_phase", 0, "assertion", "anchor_derived")])
    ev = evidence(lower(cat, default_view(cat)), "V")
    assert ev["coding_phase"] is None
    assert ev["notes"] == ["coding_phase=assertion:anchor_derived"]


def test_stated_phase_boundaries_functionality_projected():
    cat = build(anchors() + [
        claim("V:V1", "coding_phase", 0), claim("V:V1", "boundary_5", "complete"),
        claim("V:V1", "boundary_3", "partial", "assertion"),
        claim("V:V1", "functionality", "orf")])
    low = lower(cat, default_view(cat))
    ev = evidence(low, "V")
    assert ev["coding_phase"] == 0
    assert (ev["boundary_5"], ev["boundary_3"]) == ("complete", "partial")
    assert ev["notes"] == sorted(["boundary_3=assertion:stated", "boundary_5=source:stated",
                                  "coding_phase=source:stated", "functionality=source:stated"])
    assert low.spec["v"][0]["functionality"] == "orf"


def test_subregions_only_from_declared_alignment():
    assert lower(build(), default_view(build())).spec["v"][0]["subregions"] == []
    cat = build(anchors() + [claim("V:V1", "v_subregions", {"FWR1": [0, 3], "CDR1": [3, 6]})])
    sub = lower(cat, default_view(cat)).spec["v"][0]["subregions"]
    assert len(sub) == 2


def test_phase_anchor_disagreement_is_conflict():
    cat = build(anchors() + [claim("V:V1", "coding_phase", 1)])
    with pytest.raises(ReferenceArtifactError) as e:
        lower(cat, default_view(cat))
    assert e.value.kind == "conflict"


def test_admission_records_missing_anchor_findings():
    cat = build()
    low = lower(cat, default_view(cat))
    assert low.admission == {"profile": PROFILE_ID, "accepted_curatable_findings": [
        {"kind": "MissingAnchor", "segment": "J", "allele": "J:J1"},
        {"kind": "MissingAnchor", "segment": "V", "allele": "V:V1"},
        {"kind": "ReferenceEvidenceProblem", "segment": "J", "allele": "J:J1"},
        {"kind": "ReferenceEvidenceProblem", "segment": "V", "allele": "V:V1"}]}
    cat = build(anchors())
    assert lower(cat, default_view(cat)).admission["accepted_curatable_findings"] == []


def test_unadmitted_curatable_kind_is_capability_error(monkeypatch):
    monkeypatch.setattr(lowering, "ADMITTED_CURATABLE_KINDS", frozenset())
    cat = build()
    with pytest.raises(ReferenceArtifactError) as e:
        lower(cat, default_view(cat))
    assert e.value.kind == "capability"


def test_seal_makes_hash_include_context():
    cat = build(anchors())
    low = lower(cat, default_view(cat))
    before = low.refdata.content_hash()
    sealed = seal(low, CTX)
    assert sealed.is_sealed() and sealed.content_hash() != before


def test_resolver_rule_extraction_preserves_bundled_hash():
    assert dataconfig_to_refdata(data.HUMAN_IGK_OGRDB).content_hash() == \
        "sha256:6d126831c55da4f7a7f74550123709866a2835c427e16ae1d1b30c1068157b5d"
