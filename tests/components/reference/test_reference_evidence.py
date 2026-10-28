"""Sequence-bound evidence and native authoring regressions."""
import copy
import hashlib
import json
import pickle

import pytest

from GenAIRR import _engine, data
from GenAIRR.alleles import VAllele, JAllele
from GenAIRR.cartridge_builder import ReferenceCartridgeBuilder
from GenAIRR._refdata_resolver import dataconfig_to_refdata


def _builder():
    # IMGT nucleotide coordinate 309 is conserved residue104, with 3 gaps.
    v = "AAA" + "..." + "AAA" * 101 + "TGT" + "AAA"
    assert v[309:312] == "TGT"
    return ReferenceCartridgeBuilder.from_fasta(
        v_fasta=">ACC1|TRAV1*01|Mus musculus|P|V-REGION|\n" + v + "\n",
        j_fasta=">ACC2|TRAJ1*01|Mus musculus|F|J-REGION|\nTTTGGGAAAGGG\n",
        chain_type="TCR_ALPHA",
    )


def test_native_authoring_retains_header_and_unknown_cartridge_identity():
    cfg = _builder().build()
    v = next(iter(cfg.v_alleles.values()))[0]
    j = next(iter(cfg.j_alleles.values()))[0]
    assert v.anchor == 306
    assert v.functional_status == "pseudogene" and j.functional_status == "functional"
    assert v.species == "Mus musculus" and v.source == "ACC1"
    assert cfg.metadata.species is None  # no explicit cartridge-wide declaration
    assert "Unknown" in repr(cfg.metadata)
    native = dataconfig_to_refdata(cfg)
    evidence = json.loads(native.allele_evidence_json("V", 0))
    assert evidence == v.reference_evidence
    assert evidence["coding_phase"] is None
    assert evidence["anchor_status"] == "inferred"
    assert evidence["boundary_3"] == "unknown"
    assert "source_header=" in evidence["notes"][-1]
    assert cfg.cartridge_manifest()["reference_evidence"]["legacy_unannotated_count"] == 0
    assert cfg.checksum_scheme == "pickle4_roundtrip_v2"
    restored = pickle.loads(pickle.dumps(cfg))
    restored.verify_integrity()
    pickle.loads(pickle.dumps(restored)).verify_integrity()
    assert dataconfig_to_refdata(restored).allele_evidence_json("V", 0) == native.allele_evidence_json("V", 0)


def test_bare_v_does_not_guess_terminal_cysteine_and_duplicate_j_is_unknown():
    v = VAllele("TRAV1*01", "AAATGTAAA", 9)
    assert v.anchor is None and v.reference_evidence["anchor_status"] == "unknown"
    j = JAllele("TRAJ1*01", "TTTGGGAAAGGGTTTGGGAAAGGG", 24)
    assert j.anchor is None and j.anchor_meta["candidates"] == [0, 12]


def test_stale_evidence_fails_at_authoring_bridge_and_cannot_be_silently_ignored():
    cfg = _builder().build()
    v = next(iter(cfg.v_alleles.values()))[0]
    v.ungapped_seq = "C" + v.ungapped_seq[1:]
    with pytest.raises(ValueError, match="checksum"):
        dataconfig_to_refdata(cfg)


def test_explicit_phase_is_independent_and_consistency_checked():
    builder = _builder()
    v = next(iter(builder._v_alleles.values()))[0]
    e = copy.deepcopy(v.reference_evidence)
    e.update(coding_phase=0, anchor_status="verified", source="synthetic source fixture", method="explicit_fixture_annotation", boundary_3="complete")
    cfg = builder.with_reference_evidence({v.name: e}).build()
    assert json.loads(dataconfig_to_refdata(cfg).allele_evidence_json("V", 0))["coding_phase"] == 0
    e["coding_phase"] = 1
    with pytest.raises(ValueError, match="coding_phase"):
        _builder().with_reference_evidence({v.name: e}).build()


def test_every_archived_cartridge_has_checksum_bound_status_without_rewriting():
    statuses = [data.audit_config(n) for n in data.list_configs()]
    assert len(statuses) == 106
    assert all(e["source_evidence"] == "authoritative_source_snapshot_not_bundled" for e in statuses)
    bad = data.audit_config("MOUSE_TCRA_IMGT")
    assert bad["status"] == "quarantined" and bad["issue_counts"]["VAnchorNotCys"] == 197
    cfg = data.MOUSE_TCRA_IMGT
    cfg.verify_integrity()
    assert cfg.cartridge_manifest()["reference_audit"]["status"] == "quarantined"


def test_legacy_reference_hash_is_not_upgraded_by_missing_evidence():
    cfg = data.HUMAN_IGH_OGRDB
    checksum = cfg.compute_checksum()
    native = dataconfig_to_refdata(cfg)
    assert native.allele_evidence_json("V", 0) is None
    assert cfg.compute_checksum() == checksum == cfg.schema_sha256


def test_known_partial_coding_tip_cannot_be_overridden_by_boolean_assertion():
    import GenAIRR as ga
    cfg = _builder().build()
    v = next(iter(cfg.v_alleles.values()))[0]
    v.reference_evidence["boundary_3"] = "partial"
    native = dataconfig_to_refdata(cfg)
    # Partial reference remains valid for ordinary observed-span simulations.
    native.validate_strict()
    with pytest.raises(ValueError, match="partial coding boundaries"):
        (ga.Experiment.on(native).recombine()
         .process_coding_ends(reference_ends_are_coding_ends=True).compile())


def test_new_checksum_survives_repeated_roundtrips_and_fresh_hash_seed(tmp_path):
    import os
    import subprocess
    import sys
    cfg = _builder().build()
    expected = cfg.schema_sha256
    for _ in range(4):
        cfg = pickle.loads(pickle.dumps(cfg, protocol=4))
        cfg.verify_integrity()
        assert cfg.compute_checksum() == expected
    path = tmp_path / "evidence.pkl"
    path.write_bytes(pickle.dumps(cfg, protocol=4))
    script = "import pickle,sys; c=pickle.load(open(sys.argv[1],'rb')); c.verify_integrity(); print(c.compute_checksum())"
    for seed in ("1", "71", "random"):
        env = dict(os.environ, PYTHONHASHSEED=seed)
        actual = subprocess.check_output([sys.executable, "-c", script, str(path)], env=env, text=True).strip()
        assert actual == expected
    next(iter(cfg.v_alleles.values()))[0].reference_evidence["boundary_3"] = "partial"
    with pytest.raises(Exception, match="checksum mismatch"):
        cfg.verify_integrity()


def _v_only(v_header, v_seq):
    cfg = ReferenceCartridgeBuilder.from_fasta(
        v_fasta=f"{v_header}\n{v_seq}\n",
        j_fasta=">ACC2|TRAJ1*01|Mus musculus|F|J-REGION|\nTTTGGGAAAGGG\n",
        chain_type="TCR_ALPHA",
    ).build()
    return next(iter(cfg.v_alleles.values()))[0]


def test_headerless_imgt_gapped_v_keeps_its_validated_cys104_anchor():
    # e.g. an OGRDB "gapped" FASTA: IMGT gaps, allele name only in the header.
    v = _v_only(">TRAV1*01", "AAA" + "..." + "AAA" * 101 + "TGT" + "AAA")
    assert v.anchor == 306
    assert v.reference_evidence["method"] == "imgt_gapped_104"


def test_imgt_header_does_not_turn_ungapped_v_into_an_imgt_alignment():
    # Ungapped, length >= 312, with TGT at 309: not IMGT codon 104.
    seq = "GCT" * 103 + "TGT" + "GCT"
    assert seq[309:312] == "TGT" and "." not in seq
    v = _v_only(">ACC1|TRAV1*01|Mus musculus|F|V-REGION|", seq)
    assert v.anchor is None
    assert v.reference_evidence["anchor_status"] == "unknown"
