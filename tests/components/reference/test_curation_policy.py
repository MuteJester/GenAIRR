"""Curation policy v1 — explicit catalogue filtering.

The cartridge model now separates three concerns:
  - Validation describes what's in the catalogue.
  - Curation selects which subset participates in simulation.
  - Simulation runs against the curated cartridge.

This module pins curation behaviour:
  - ``functional_anchors_only`` drops V/J alleles with anchor
    problems (missing, out-of-bounds, codon AA mismatch).
  - D pools pass through unchanged.
  - Structural corruption (duplicates, invalid bytes) is NOT
    fixed by curation — it still surfaces from the validator.
  - Curation tags ``identity.source`` and changes the content
    hash so trace files distinguish raw from curated artefacts.
"""
from __future__ import annotations

import pytest

import GenAIRR as ga
from GenAIRR import _engine as ge


# ──────────────────────────────────────────────────────────────────
# Direct API on engine-native RefDataConfig
# ──────────────────────────────────────────────────────────────────


def _mixed_vj_cfg() -> "ge.RefDataConfig":
    """V pool: one valid (Cys) + one non-Cys.
    J pool: one valid (W) + one anchorless.
    """
    cfg = ge.RefDataConfig.vj()
    cfg.add_v_allele("MYV-good*01", "MYV-good", b"TGTAAACCC", anchor=0)
    cfg.add_v_allele("MYV-gly*01", "MYV-gly", b"GGGAAACCC", anchor=0)
    cfg.add_j_allele("MYJ-good*01", "MYJ-good", b"TGGAAACCC", anchor=0)
    cfg.add_j_allele("MYJ-orphan*01", "MYJ-orphan", b"GGG", anchor=None)
    return cfg


def test_curated_functional_anchors_only_drops_non_cys_v() -> None:
    cfg = _mixed_vj_cfg()
    curated = cfg.curated("functional_anchors_only")
    assert curated.v_pool_size() == 1
    assert curated.v_allele(0).name == "MYV-good*01"


def test_curated_functional_anchors_only_drops_missing_anchor_j() -> None:
    cfg = _mixed_vj_cfg()
    curated = cfg.curated("functional_anchors_only")
    assert curated.j_pool_size() == 1
    assert curated.j_allele(0).name == "MYJ-good*01"


def test_curated_raw_is_identity() -> None:
    cfg = _mixed_vj_cfg()
    curated = cfg.curated("raw")
    assert curated.v_pool_size() == 2
    assert curated.j_pool_size() == 2


def test_curate_in_place_is_equivalent_to_curated() -> None:
    a = _mixed_vj_cfg()
    a.curate("functional_anchors_only")

    b = _mixed_vj_cfg().curated("functional_anchors_only")

    assert a.v_pool_size() == b.v_pool_size()
    assert a.j_pool_size() == b.j_pool_size()
    assert a.content_hash() == b.content_hash()


def test_curation_tags_identity_source() -> None:
    cfg = _mixed_vj_cfg()
    curated = cfg.curated("functional_anchors_only")
    assert curated.identity()["source"] == "curated:functional_anchors_only"


def test_curation_extends_existing_identity_source() -> None:
    cfg = _mixed_vj_cfg()
    cfg.set_identity(source="DataConfig")
    curated = cfg.curated("functional_anchors_only")
    assert curated.identity()["source"] == "DataConfig|curated:functional_anchors_only"


def test_curation_preserves_species_locus_reference_set() -> None:
    cfg = _mixed_vj_cfg()
    cfg.set_identity(species="Human", locus="IGK", reference_set="OGRDB")
    curated = cfg.curated("functional_anchors_only")
    id_ = curated.identity()
    assert id_["species"] == "Human"
    assert id_["locus"] == "IGK"
    assert id_["reference_set"] == "OGRDB"


def test_curation_changes_content_hash() -> None:
    cfg = _mixed_vj_cfg()
    raw_hash = cfg.content_hash()
    curated_hash = cfg.curated("functional_anchors_only").content_hash()
    assert raw_hash != curated_hash


def test_unknown_policy_raises_value_error() -> None:
    cfg = _mixed_vj_cfg()
    with pytest.raises(ValueError, match="unknown curation policy"):
        cfg.curated("bogus")


# ──────────────────────────────────────────────────────────────────
# Curation does NOT silently fix structural corruption
# ──────────────────────────────────────────────────────────────────


def test_curation_does_not_remove_invalid_byte_alleles() -> None:
    cfg = _mixed_vj_cfg()
    cfg.add_v_allele("MYV-bad*01", "MYV-bad", b"TGT.AAACC", anchor=0)
    curated = cfg.curated("functional_anchors_only")
    issues = curated.validate()
    assert any(i["kind"] == "InvalidAlleleByte" for i in issues), (
        "invalid byte must remain Fatal post-curation; curation never "
        f"silently hides corruption. issues={issues}"
    )


def test_curation_does_not_dedupe_allele_names() -> None:
    cfg = _mixed_vj_cfg()
    # Both Cys-anchored — they pass curation, then DuplicateAlleleName
    # surfaces.
    cfg.add_v_allele("MYV-good*01", "MYV-good", b"TGTAAACCC", anchor=0)
    curated = cfg.curated("functional_anchors_only")
    issues = curated.validate()
    assert any(
        i["kind"] == "DuplicateAlleleName" and i["name"] == "MYV-good*01"
        for i in issues
    )


def test_curation_does_not_silence_locus_chain_mismatch() -> None:
    cfg = _mixed_vj_cfg()
    cfg.set_identity(locus="IGH")  # VJ + IGH → mismatch
    curated = cfg.curated("functional_anchors_only")
    issues = curated.validate()
    assert any(i["kind"] == "LocusChainTypeMismatch" for i in issues)


# ──────────────────────────────────────────────────────────────────
# Compile-time behaviour: curation empties a required pool
# ──────────────────────────────────────────────────────────────────


def test_curation_emptying_v_pool_fails_compile_with_empty_required() -> None:
    cfg = ge.RefDataConfig.vj()
    # All V alleles are non-Cys → curation drops everything.
    cfg.add_v_allele("MYV-gly*01", "MYV-gly", b"GGGAAACCC", anchor=0)
    cfg.add_j_allele("MYJ*01", "MYJ", b"TGGAAACCC", anchor=0)
    curated = cfg.curated("functional_anchors_only")
    assert curated.v_pool_size() == 0
    issues = curated.validate()
    empty_v = [
        i for i in issues
        if i["kind"] == "EmptyRequiredPool" and i["segment"] == "V"
    ]
    assert empty_v, f"expected EmptyRequiredPool[V] after curation, got {issues}"


# ──────────────────────────────────────────────────────────────────
# Experiment.curate_refdata — end-to-end with bundled mouse_igh
# ──────────────────────────────────────────────────────────────────


def test_bundled_mouse_igh_fails_strict_compile_raw() -> None:
    """Bundled mouse_igh contains pseudogene-shape V alleles whose
    anchors don't translate to Cys. Strict compile must fail with
    the curatable remediation hint pointing the user at either
    allow_curatable_refdata or filter_functional_alleles."""
    with pytest.raises(ValueError, match="allow_curatable_refdata"):
        ga.Experiment.on("mouse_igh").recombine().compile()


def test_bundled_mouse_igh_passes_strict_after_curation() -> None:
    compiled = (
        ga.Experiment.on("mouse_igh")
        .curate_refdata("functional_anchors_only")
        .recombine()
        .compile()
    )
    assert compiled is not None


def test_curated_mouse_igh_identity_records_curation() -> None:
    rd = (
        ga.Experiment.on("mouse_igh")
        .curate_refdata("functional_anchors_only")
        .refdata
    )
    assert rd.identity()["source"] == "DataConfig|curated:functional_anchors_only"


def test_curated_mouse_igh_drops_some_v_alleles() -> None:
    raw = ga.Experiment.on("mouse_igh").refdata
    curated = ga.Experiment.on("mouse_igh").curate_refdata(
        "functional_anchors_only"
    ).refdata
    assert curated.v_pool_size() < raw.v_pool_size()
    assert curated.v_pool_size() > 0


def test_curate_refdata_returns_self_for_chaining() -> None:
    exp = ga.Experiment.on("mouse_igh")
    chained = exp.curate_refdata("functional_anchors_only")
    assert chained is exp


# ──────────────────────────────────────────────────────────────────
# Trace replay — refdata signature distinguishes raw vs curated
# ──────────────────────────────────────────────────────────────────


def test_trace_replay_rejects_when_cartridge_curation_changes() -> None:
    """A trace file built against a curated cartridge must not
    replay against the raw one (or vice versa) — the content hash
    differs, so the trace-file integrity check trips."""
    compiled_curated = (
        ga.Experiment.on("human_igh")
        .curate_refdata("functional_anchors_only")
        .recombine()
        .compile()
    )
    outcome = compiled_curated.simulator.run(seed=42)
    tf = compiled_curated.simulator.trace_file_from(outcome, seed=42)

    # Now try to replay against the RAW human_igh cartridge.
    compiled_raw = ga.Experiment.on("human_igh").recombine().compile()
    with pytest.raises(ValueError):
        compiled_raw.simulator.replay_from_trace_file(tf)
