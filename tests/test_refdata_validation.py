"""RefData validator — sanity gates plus the bundled-preset suite.

Mirrors the cargo lib tests on the Python boundary so consumers can
hit `refdata.validate()` / `validate_strict()` directly. The bundled
human IGH / IGK / IGL presets MUST validate clean — any regression
here means the validator has drifted from the data the engine ships
with.
"""
from __future__ import annotations

import pytest

import GenAIRR as ga


# ──────────────────────────────────────────────────────────────────
# Bundled-preset gate. If this ever fires, the validator and the
# bundled refdata are out of sync — fix one or the other.
# ──────────────────────────────────────────────────────────────────


@pytest.mark.parametrize("preset", ["human_igh", "human_igk", "human_igl"])
def test_bundled_preset_passes_validation(preset: str) -> None:
    refdata = ga.Experiment.on(preset).refdata
    issues = refdata.validate()
    assert issues == [], (
        f"{preset} validation failed: {len(issues)} issue(s) - "
        f"first={issues[0] if issues else None}"
    )
    refdata.validate_strict()  # must not raise


# ──────────────────────────────────────────────────────────────────
# Synthetic fixtures — small, deterministic, exercise each rule.
# ──────────────────────────────────────────────────────────────────


def _make_minimal_vdj():
    """Minimal valid VDJ refdata: builder helper for failure tests."""
    from GenAIRR._engine import RefDataConfig

    cfg = RefDataConfig.vdj()
    cfg.add_v_allele("IGHV1-1*01", "IGHV1-1", b"TGTAAACCC", anchor=0)
    cfg.add_d_allele("IGHD1-1*01", "IGHD1-1", b"GGGCCCAAA")
    cfg.add_j_allele("IGHJ1*01", "IGHJ1", b"TGGAAACCC", anchor=0)
    return cfg


def test_minimal_vdj_passes() -> None:
    cfg = _make_minimal_vdj()
    assert cfg.validate() == []


def test_empty_v_pool_flagged() -> None:
    from GenAIRR._engine import RefDataConfig

    cfg = RefDataConfig.vdj()
    cfg.add_d_allele("IGHD1*01", "IGHD1", b"GGGCCCAAA")
    cfg.add_j_allele("IGHJ1*01", "IGHJ1", b"TGGAAACCC", anchor=0)
    kinds = [i["kind"] for i in cfg.validate()]
    assert "EmptyRequiredPool" in kinds


def test_empty_d_pool_rejected_on_vdj() -> None:
    from GenAIRR._engine import RefDataConfig

    cfg = RefDataConfig.vdj()
    cfg.add_v_allele("IGHV1*01", "IGHV1", b"TGTAAACCC", anchor=0)
    cfg.add_j_allele("IGHJ1*01", "IGHJ1", b"TGGAAACCC", anchor=0)
    kinds = [(i["kind"], i.get("segment")) for i in cfg.validate()]
    assert ("EmptyRequiredPool", "D") in kinds


def test_empty_d_pool_allowed_on_vj() -> None:
    from GenAIRR._engine import RefDataConfig

    cfg = RefDataConfig.vj()
    cfg.add_v_allele("IGKV1*01", "IGKV1", b"TGTAAACCC", anchor=0)
    cfg.add_j_allele("IGKJ1*01", "IGKJ1", b"TTCAAACCC", anchor=0)
    issues = cfg.validate()
    # VJ chains have empty D pools by design — must NOT be flagged.
    assert not any(
        i["kind"] == "EmptyRequiredPool" and i.get("segment") == "D"
        for i in issues
    )


def test_invalid_byte_reports_exact_position() -> None:
    cfg = _make_minimal_vdj()
    cfg.add_v_allele("IGHV2*01", "IGHV2", b"TGT.AAA", anchor=0)
    issues = cfg.validate()
    bad = [i for i in issues if i["kind"] == "InvalidAlleleByte"]
    assert bad, f"expected an InvalidAlleleByte issue, got: {issues}"
    only = bad[0]
    assert only["pos"] == 3
    assert only["byte"] == ord(".")


def test_v_anchor_non_cys_flagged() -> None:
    cfg = _make_minimal_vdj()
    # GGG (Gly) at anchor.
    cfg.add_v_allele("IGHV-bad*01", "IGHV-bad", b"GGGAAACCC", anchor=0)
    issues = cfg.validate()
    assert any(i["kind"] == "VAnchorNotCys" and i["aa"] == "G" for i in issues)


def test_j_anchor_unexpected_aa_for_igh_flagged() -> None:
    cfg = _make_minimal_vdj()
    # Under the new ReferenceRules layer, J anchor expectations are
    # config-level (no per-allele name inference). To exercise the
    # IGH convention here, configure the rule explicitly.
    cfg.set_j_anchor_rule(expected_aa=["W"])
    # IGH J expects W; TTC = F → unexpected.
    cfg.add_j_allele("IGHJ-bad*01", "IGHJ-bad", b"TTCAAACCC", anchor=0)
    issues = cfg.validate()
    j_unexpected = [i for i in issues if i["kind"] == "JAnchorUnexpectedAa"]
    assert j_unexpected, f"expected J anchor mismatch, got: {issues}"
    assert j_unexpected[0]["aa"] == "F"
    assert "W" in j_unexpected[0]["expected"]


def test_duplicate_allele_name_flagged() -> None:
    cfg = _make_minimal_vdj()
    cfg.add_v_allele("IGHV1-1*01", "IGHV1-1", b"TGTAAA", anchor=0)
    issues = cfg.validate()
    assert any(
        i["kind"] == "DuplicateAlleleName" and i["name"] == "IGHV1-1*01"
        for i in issues
    )


def test_validate_strict_raises_on_issue() -> None:
    cfg = _make_minimal_vdj()
    cfg.add_v_allele("IGHV1-1*01", "IGHV1-1", b"TGTAAA", anchor=0)
    with pytest.raises(ValueError, match=r"validation issue"):
        cfg.validate_strict()
