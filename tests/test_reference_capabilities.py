import hashlib

from GenAIRR._reference.capabilities import CAPABILITIES, evaluate
from GenAIRR._reference.model import (
    Allele, Catalogue, Claim, ReferenceHeader, Taxon,
)

V_SEQ = "AAATGTAAAGGG"   # TGT (C) at 3
J_SEQ = "TGGGGCCAAGGG"   # TGG (W) at 0
D_SEQ = "GGGTTTAAA"
VIEW = {"V": ["V:V1"], "D": ["D:D1"], "J": ["J:J1"]}


def _allele(seg, label, seq):
    return Allele(f"{seg}:{label}", seg, label, label, "declared", seq,
                  "sha256:" + hashlib.sha256(seq.encode("ascii")).hexdigest(), None, (), ())


def claim(subject, predicate, value, support="source", derivation="stated", status="accepted"):
    src = "doi:x" if support == "source" else None
    return Claim(subject, predicate, value, support, derivation, status, src, "m")


def build(claims=(), locus="IGH", chain="vdj"):
    alleles = [_allele("V", "V1", V_SEQ), _allele("J", "J1", J_SEQ)]
    if chain == "vdj":
        alleles.insert(1, _allele("D", "D1", D_SEQ))
    cat = Catalogue(
        ReferenceHeader("t", "1", Taxon("Homo sapiens", 9606), locus, chain, None),
        {"coordinates": "zero_based_half_open", "v_alignment": "none"},
        (), (), tuple(alleles), tuple(claims), (), (), ())
    return Catalogue.from_dict(cat.to_dict())


def anchors(support="source", v=3, j=0):
    return [claim("V:V1", "anchor", v, support), claim("J:J1", "anchor", j, support)]


def test_bare_catalogue_supports_only_inspection_and_observed_span():
    r = evaluate(build(), VIEW)
    assert r.available() == ("catalogue_inspection", "observed_span_simulation")
    assert r.phase_support == "none"
    assert tuple(s.name for s in r.statuses) == CAPABILITIES


def test_source_anchors_give_source_evidence():
    r = evaluate(build(anchors()), VIEW)
    assert r.status("source_annotated_junction").evidence == "source"
    assert r.status("canonical_productive").available
    assert r.status("canonical_productive").evidence == "source"
    assert "canonical_productive" in r.available(evidence="source")


def test_assertion_anchor_downgrades_evidence():
    r = evaluate(build([claim("V:V1", "anchor", 3, "source"),
                        claim("J:J1", "anchor", 0, "assertion")]), VIEW)
    st = r.status("canonical_productive")
    assert st.available and st.evidence == "assertion"
    assert "canonical_productive" not in r.available(evidence="source")
    assert r.available(evidence="source") == ("catalogue_inspection", "observed_span_simulation")


def test_unsupported_complete_boundary_is_assertion_level_only():
    cl = [claim("V:V1", "boundary_3", "complete", "assertion"),
          claim("J:J1", "boundary_5", "complete", "source"),
          claim("D:D1", "boundary_5", "complete", "source"),
          claim("D:D1", "boundary_3", "complete", "source")]
    st = evaluate(build(cl), VIEW).status("original_coding_ends")
    assert st.available and st.evidence == "assertion"


def test_anchor_derived_phase_is_not_independent():
    cl = anchors() + [claim("V:V1", "coding_phase", 0, "assertion", "anchor_derived"),
                      claim("J:J1", "coding_phase", 0, "assertion", "anchor_derived")]
    r = evaluate(build(cl), VIEW)
    assert r.phase_support == "anchor_derived"
    assert "phase support: anchor_derived" in r.explain()
    assert "anchor_derived" in r.status("canonical_productive").note


def test_phase_claims_do_not_change_capabilities():
    base = evaluate(build(), VIEW)
    cl = [claim("V:V1", "coding_phase", 0), claim("J:J1", "coding_phase", 1)]
    r = evaluate(build(cl), VIEW)
    assert r.phase_support == "independent_source"
    assert r.available() == base.available()
    assert not r.status("canonical_productive").available
    cl = [claim("V:V1", "coding_phase", 0, "assertion"), claim("J:J1", "coding_phase", 1, "source")]
    assert evaluate(build(cl), VIEW).phase_support == "independent_assertion"
    r2 = evaluate(build(anchors() + cl), VIEW)
    assert r2.status("canonical_productive").evidence == "source"


def test_conflict_blocks_junction():
    cl = anchors() + [claim("J:J1", "anchor", 1, "source", status="conflict")]
    r = evaluate(build(cl), VIEW)
    assert not r.status("source_annotated_junction").available
    assert r.status("source_annotated_junction").unmet[0].allele == "J:J1"


def test_noncanonical_v_anchor_blocks_canonical_only():
    r = evaluate(build(anchors(v=0)), VIEW)
    assert r.status("source_annotated_junction").available
    assert not r.status("canonical_productive").available
    r = evaluate(build(anchors(j=3)), VIEW)  # GGC is G, not W
    assert not r.status("canonical_productive").available
    # IGK expects F, so a W codon is not canonical there
    r = evaluate(build(anchors(), locus="IGK", chain="vj"), {"V": ["V:V1"], "J": ["J:J1"]})
    assert r.status("source_annotated_junction").available
    assert not r.status("canonical_productive").available


def test_motif_candidate_hint():
    cl = [claim("V:V1", "anchor", 3), Claim("J:J1", "anchor", 0, "none", "unique_motif",
                                            "candidate", None, "motif")]
    r = evaluate(build(cl), VIEW)
    (u,) = r.status("source_annotated_junction").unmet
    assert u.hint == ('add an "anchor" annotation (support "source" with a citation, or "assertion")'
                      '; candidate at 0 (unique J motif) — copy it only if your evidence confirms it')
    assert f"- J:J1: accepted anchor — {u.hint}" in r.explain()


def test_original_coding_ends_requires_d_both_ends():
    cl = [claim("V:V1", "boundary_3", "complete"), claim("J:J1", "boundary_5", "complete"),
          claim("D:D1", "boundary_5", "complete")]
    st = evaluate(build(cl), VIEW).status("original_coding_ends")
    assert not st.available
    assert [(u.allele, u.requirement) for u in st.unmet] == [("D:D1", 'accepted boundary_3 == "complete"')]
    assert st.unmet[0].hint == 'add "boundary_3": {"value": "complete", ...} with its support'


def test_fatal_native_findings_block_observed_span():
    r = evaluate(build(), VIEW, native_fatal=[{"kind": "x", "detail": "bad"}])
    assert not r.status("observed_span_simulation").available
    assert r.available() == ("catalogue_inspection",)


def test_c4_capabilities_unavailable_with_note():
    r = evaluate(build(anchors()), VIEW)
    for n in ("genotype_constrained", "duplex_family_replay"):
        st = r.status(n)
        assert not st.available and st.note == "not yet supported for portable references"
    assert "unavailable genotype_constrained: not yet supported for portable references" in r.explain()
    assert r.explain().splitlines()[-1] == "phase support: none"
    assert r.to_dict()["phase_support"] == "none"
