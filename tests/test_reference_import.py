import copy
import hashlib
import json

import pytest

from GenAIRR._reference.errors import ReferenceArtifactError
from GenAIRR._reference.importer import import_local, scan_fasta
from GenAIRR._reference.model import Catalogue

V_SEQ = "ACGT" * 10
J_SEQ = "TGGGGCCAAGGGACC"
J_KAPPA = "TGGTTCGGCGGAGGGACCAAG"
GAPPED_V = "AAA" + "..." + "AAA" * 101 + "TGT" + "AAA"


def digest(seq):
    return "sha256:" + hashlib.sha256(seq.encode("ascii")).hexdigest()


def entry(segment, label, seq, **extra):
    base = {"segment": segment, "label": label, "sequence_sha256": digest(seq), "gene": None,
            "functionality": None, "anchor": None, "coding_phase": None,
            "boundary_5": None, "boundary_3": None, "notes": []}
    base.update(extra)
    return base


def stated(value, support="assertion", source=None, method="m"):
    return {"value": value, "support": support, "source": source, "method": method}


def sidecar(alleles=(), *, chain_type="vj", v_alignment="none", locus="IGH", relationships=()):
    return {
        "schema": "genairr_reference_annotations", "schema_version": 1,
        "reference": {"name": "lab", "version": "1", "locus": locus, "chain_type": chain_type,
                      "taxon": {"scientific_name": "Homo sapiens", "taxonomy_id": 9606},
                      "build_date": "2026-10-21"},
        "conventions": {"coordinates": "zero_based_half_open", "v_alignment": v_alignment},
        "alleles": list(alleles), "relationships": list(relationships), "suggestions": [],
    }


def write(tmp_path, name, text, binary=False):
    p = tmp_path / name
    if binary:
        p.write_bytes(text)
    else:
        p.write_text(text, encoding="utf-8")
    return p


def run(tmp_path, v="", j=f">J1\n{J_SEQ}\n", d=None, side=None):
    side = side if side is not None else sidecar(chain_type="vj" if d is None else "vdj")
    return import_local(
        v=write(tmp_path, "v.fasta", v), j=write(tmp_path, "j.fasta", j),
        d=None if d is None else write(tmp_path, "d.fasta", d),
        annotations=write(tmp_path, "ann.json", json.dumps(side)))


def test_every_record_is_accounted(tmp_path):
    v = f">a\n{V_SEQ}\n>a\n{V_SEQ}\n>n\nACGN\n>e\n\n>bad,label\n{V_SEQ}\n"
    cat = run(tmp_path, v)
    vrecs = [r for r in cat.records if r.id.startswith("rec:V:")]
    assert [r.status for r in vrecs] == ["accepted", "duplicate", "rejected", "rejected", "rejected"]
    assert vrecs[1].duplicate_of == vrecs[0].id
    assert "unsupported_base" in vrecs[2].reason and "'N'" in vrecs[2].reason and "3" in vrecs[2].reason
    assert vrecs[3].reason == "empty_sequence"
    assert vrecs[4].reason.startswith("unrepresentable_label: ")
    assert len(cat.alleles) == 2  # V:a and J:J1


def test_conflicting_label_rejects_all_copies(tmp_path):
    cat = run(tmp_path, f">a\n{V_SEQ}\n>a\n{V_SEQ[:-4]}\n>b\n{V_SEQ}\n")
    vrecs = [r for r in cat.records if r.id.startswith("rec:V:")]
    assert [r.status for r in vrecs] == ["rejected", "rejected", "accepted"]
    assert all(r.reason == "conflicting_sequence_for_label" for r in vrecs[:2])
    assert [a.label for a in cat.alleles if a.segment == "V"] == ["b"]


def test_crlf_fasta_matches_lf(tmp_path):
    lf = run(tmp_path, f">a desc\n{V_SEQ}\n")
    crlf_dir = tmp_path / "crlf"
    crlf_dir.mkdir()
    crlf = run(crlf_dir, f">a desc\n{V_SEQ}\n".replace("\n", "\r\n"), j=f">J1\r\n{J_SEQ}\r\n")
    assert [a.to_dict() for a in lf.alleles] == [a.to_dict() for a in crlf.alleles]
    assert [r.label for r in lf.records] == [r.label for r in crlf.records] == ["a", "J1"]
    assert all("\r" not in r.header for r in crlf.records)


def test_imgt_pipe_header_label(tmp_path):
    header = ">X1|IGHV1-2*02|Homo sapiens|F|V-REGION|1..296|296 nt|1| | | | |296+0=296| | |"
    cat = run(tmp_path, f"{header}\n{V_SEQ.lower()}\n")
    allele = cat.alleles[0]
    assert allele.label == "IGHV1-2*02" and allele.id == "V:IGHV1-2*02"
    assert allele.sequence == V_SEQ
    assert [t.kind for t in allele.transforms] == ["uppercase"]
    other = run(tmp_path / ".", f">plain tail\n{V_SEQ}\n")
    assert other.alleles[0].transforms == ()
    assert other.alleles[0].label == "plain"


def test_dots_without_declared_alignment_rejected(tmp_path):
    cat = run(tmp_path, f">a\n{GAPPED_V}\n")
    rec = cat.records[0]
    assert rec.status == "rejected" and rec.reason == "gap_without_declared_alignment"


def declared(tmp_path, anchor=None):
    extra = {"anchor": anchor} if anchor else {}
    side = sidecar([entry("V", "v1", GAPPED_V.replace(".", ""), **extra)], v_alignment="imgt_gapped")
    return run(tmp_path, f">v1\n{GAPPED_V}\n", side=side)


def test_declared_alignment_gives_cys104_anchor_and_subregions(tmp_path):
    cat = declared(tmp_path)
    allele = cat.allele("V:v1")
    assert allele.sequence == GAPPED_V.replace(".", "")
    assert allele.gapped_sequence == GAPPED_V
    assert [t.kind for t in allele.transforms] == ["remove_imgt_gaps"]
    (anchor,) = cat.claims_for("V:v1", "anchor")
    assert (anchor.value, anchor.support, anchor.derivation, anchor.status) == \
        (306, "assertion", "imgt_alignment", "accepted")
    (sub,) = cat.claims_for("V:v1", "v_subregions")
    assert sub.value == {"FWR1": [0, 75], "CDR1": [75, 111], "FWR2": [111, 162],
                         "CDR2": [162, 192], "FWR3": [192, 309]}
    assert (sub.support, sub.derivation, sub.status) == ("assertion", "imgt_alignment", "accepted")


def test_stated_vs_alignment_anchor_conflict(tmp_path):
    cat = declared(tmp_path, anchor=stated(303))
    claims = cat.claims_for("V:v1", "anchor")
    assert sorted((c.derivation, c.value, c.status) for c in claims) == \
        [("imgt_alignment", 306, "conflict"), ("stated", 303, "conflict")]
    agree = declared(tmp_path / ".", anchor=stated(306))
    assert {c.status for c in agree.claims_for("V:v1", "anchor")} == {"accepted"}


def test_stale_sidecar_digest_error_names_expected_form(tmp_path):
    side = sidecar([entry("V", "a", V_SEQ[:-4])])
    with pytest.raises(ReferenceArtifactError) as err:
        run(tmp_path, f">a\n{V_SEQ}\n", side=side)
    assert err.value.kind == "conflict"
    msg = err.value.message
    assert "V:a" in msg and digest(V_SEQ) in msg and digest(V_SEQ[:-4]) in msg
    assert "digest is over the uppercase ungapped A/C/G/T sequence" in msg


def test_motif_candidate_not_accepted(tmp_path):
    side = sidecar([entry("J", "k1", J_KAPPA)], chain_type="vj", locus="IGK")
    cat = run(tmp_path, f">v\n{V_SEQ}\n", j=f">k1\n{J_KAPPA}\n", side=side)
    (claim,) = cat.claims_for("J:k1", "anchor")
    assert (claim.support, claim.derivation, claim.status) == ("none", "unique_motif", "candidate")
    scan = scan_fasta(v=tmp_path / "v.fasta", j=tmp_path / "j.fasta", d=None,
                      v_alignment="none", locus="IGK")
    records, suggestions = scan
    assert records["J"] == [{"label": "k1", "sequence_sha256": digest(J_KAPPA)}]
    assert [s["allele"] for s in suggestions] == ["J:k1"]
    assert suggestions[0]["predicate"] == "anchor" and suggestions[0]["derivation"] == "unique_motif"
    assert suggestions[0]["value"] == claim.value and isinstance(suggestions[0]["note"], str)


def test_alias_between_executable_records_rejected(tmp_path):
    side = sidecar([entry("V", "a", V_SEQ), entry("V", "b", V_SEQ + "AC")],
                   relationships=[{"kind": "alias", "allele": "V:a", "entity": "b", "source": None}])
    with pytest.raises(ReferenceArtifactError) as err:
        run(tmp_path, f">a\n{V_SEQ}\n>b\n{V_SEQ}AC\n", side=side)
    assert err.value.kind == "conflict"
    assert "alias names an executable record; keep one record and declare the other name as an alias entity" \
        in err.value.message


def test_paralogs_with_identical_sequences_stay_separate(tmp_path):
    side = sidecar(relationships=[{"kind": "paralog", "allele": "V:a", "entity": "b", "source": None}])
    cat = run(tmp_path, f">a\n{V_SEQ}\n>b\n{V_SEQ}\n", side=side)
    assert [a.id for a in cat.alleles if a.segment == "V"] == ["V:a", "V:b"]
    (finding,) = cat.findings
    assert finding.kind == "identical_sequence_distinct_labels"
    assert set(finding.subjects) == {"V:a", "V:b"}
    assert [r.kind for r in cat.relationships] == ["paralog"]


def test_relationship_only_entity_not_executable(tmp_path):
    side = sidecar(relationships=[{"kind": "alias", "allele": "V:a", "entity": "lab_name", "source": "s"}])
    cat = run(tmp_path, f">a\n{V_SEQ}\n", side=side)
    assert [(e.name, e.executable) for e in cat.entities] == [("lab_name", False)]
    assert "V:lab_name" not in {a.id for a in cat.alleles}
    assert cat.relationships[0].source == "s"


def test_labels_never_renamed_and_gene_basis(tmp_path):
    side = sidecar([entry("V", "lab_x", V_SEQ, gene="GENE-X")])
    cat = run(tmp_path, f">lab_x\n{V_SEQ}\n>IGHV3-1*01\n{V_SEQ}AC\n>plain\n{V_SEQ}GG\n", side=side)
    by_label = {a.label: a for a in cat.alleles if a.segment == "V"}
    assert set(by_label) == {"lab_x", "IGHV3-1*01", "plain"}
    assert (by_label["lab_x"].gene, by_label["lab_x"].gene_basis) == ("GENE-X", "declared")
    assert (by_label["IGHV3-1*01"].gene, by_label["IGHV3-1*01"].gene_basis) == ("IGHV3-1", "label_prefix")
    assert (by_label["plain"].gene, by_label["plain"].gene_basis) == ("plain", "label")


def test_topology_and_sources(tmp_path):
    with pytest.raises(ReferenceArtifactError) as err:
        run(tmp_path, f">a\n{V_SEQ}\n", d=f">d1\n{V_SEQ}\n", side=sidecar(chain_type="vj"))
    assert err.value.kind == "malformed"
    with pytest.raises(ReferenceArtifactError) as err:
        run(tmp_path, f">a\n{V_SEQ}\n", side=sidecar(chain_type="vdj"))
    assert err.value.kind == "malformed"
    cat = run(tmp_path, f">a\n{V_SEQ}\n")
    assert {s.file_name for s in cat.sources} == {"v.fasta", "j.fasta"}
    assert all("/" not in s.file_name for s in cat.sources)


def test_sidecar_entry_without_allele_is_conflict(tmp_path):
    with pytest.raises(ReferenceArtifactError) as err:
        run(tmp_path, f">a\n{V_SEQ}\n", side=sidecar([entry("V", "ghost", V_SEQ)]))
    assert err.value.kind == "conflict"


def test_stated_claims_copied(tmp_path):
    side = sidecar([entry("V", "a", V_SEQ, functionality=stated("functional", "source", "doi:x"),
                          coding_phase={**stated(1), "derivation": "anchor_derived"})])
    cat = run(tmp_path, f">a\n{V_SEQ}\n", side=side)
    (f,) = cat.claims_for("V:a", "functionality")
    assert (f.derivation, f.support, f.source, f.method, f.status) == \
        ("stated", "source", "doi:x", "m", "accepted")
    (p,) = cat.claims_for("V:a", "coding_phase")
    assert p.derivation == "anchor_derived"


def test_round_trip_validates(tmp_path):
    cat = declared(tmp_path, anchor=stated(306))
    again = Catalogue.from_dict(copy.deepcopy(cat.to_dict()))
    assert again.digest() == cat.digest()


def test_vdj_import_includes_d_segment(tmp_path):
    cat = run(tmp_path, f">a\n{V_SEQ}\n", d=f">d1\n{V_SEQ[:12]}\n", side=sidecar(chain_type="vdj"))
    assert [a.id for a in cat.alleles] == ["V:a", "D:d1", "J:J1"]
    assert {s.id for s in cat.sources} == {"src:V", "src:D", "src:J"}
