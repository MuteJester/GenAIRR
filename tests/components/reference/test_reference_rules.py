"""Reference Rules v1 — Python surface tests.

The Rust ``RefDataConfig`` now carries a typed :class:`ReferenceRules`
slice (anchor expectations + allowed alphabet). The validator and the
compile gate consult it instead of inferring expectations from allele
names. Python exposes:

  - ``set_v_anchor_rule`` / ``set_j_anchor_rule`` / ``set_allowed_bases``
  - ``v_anchor_rule`` / ``j_anchor_rule`` / ``allowed_bases`` readers

Bundled refdata is stamped with locus-appropriate rules at
``dataconfig_to_refdata`` time so the default ``Experiment.compile()``
preserves prior behaviour without per-allele name inference.
"""
from __future__ import annotations

import pytest

import GenAIRR as ga
from GenAIRR import _engine as ge


# ──────────────────────────────────────────────────────────────────
# Default rules
# ──────────────────────────────────────────────────────────────────


def test_default_rules_match_legacy_expectations() -> None:
    cfg = ge.RefDataConfig.vj()
    v = cfg.v_anchor_rule()
    j = cfg.j_anchor_rule()
    assert v["expected_aa"] == ["C"]
    assert v["required"] is True
    assert v["missing_severity"] == "curatable"
    assert v["mismatch_severity"] == "curatable"
    assert j["expected_aa"] == ["W", "F"]
    assert j["required"] is True
    assert cfg.allowed_bases() == ["A", "C", "G", "T", "N"]


# ──────────────────────────────────────────────────────────────────
# Bundled presets still get the right rule at load time
# ──────────────────────────────────────────────────────────────────


@pytest.mark.parametrize(
    "preset,expected_j",
    [
        ("human_igh", ["W"]),
        ("human_igk", ["F"]),
        ("human_igl", ["F"]),
        ("human_tcrb", ["F"]),
        ("mouse_igh", ["W"]),
    ],
)
def test_bundled_preset_gets_locus_rule(preset: str, expected_j: list[str]) -> None:
    refdata = ga.Experiment.on(preset).refdata
    assert refdata.j_anchor_rule()["expected_aa"] == expected_j
    assert refdata.v_anchor_rule()["expected_aa"] == ["C"]


def test_default_igh_still_accepts_j_w() -> None:
    """Spec point: default IGH still accepts J `W`. Bundled IGH J
    alleles all carry TGG (W); strict compile must not flag them."""
    ga.Experiment.on("human_igh").recombine().compile()


def test_default_igk_still_accepts_j_f() -> None:
    """Spec point: default IGK still accepts J `F`. Bundled IGK J
    alleles all carry TTC/TTT (F); strict compile must accept them."""
    ga.Experiment.on("human_igk").recombine().compile()


# ──────────────────────────────────────────────────────────────────
# Custom rules — direct API
# ──────────────────────────────────────────────────────────────────


def _custom_vj_minimal() -> "ge.RefDataConfig":
    """Returns a VJ refdata with a Cys-anchored V and an empty J pool
    — caller adds J alleles to exercise specific J rules."""
    cfg = ge.RefDataConfig.vj()
    cfg.add_v_allele("MYV1*01", "MYV1", b"TGTAAACCC", anchor=0)
    return cfg


def test_custom_j_rule_y_accepts_tat_and_tac() -> None:
    """Spec point: a `['Y']` J anchor rule accepts `TAT`/`TAC`."""
    cfg = _custom_vj_minimal()
    cfg.set_j_anchor_rule(expected_aa=["Y"])
    cfg.add_j_allele("MYJ-tat*01", "MYJ-tat", b"TATAAACCC", anchor=0)
    cfg.add_j_allele("MYJ-tac*01", "MYJ-tac", b"TACAAACCC", anchor=0)
    issues = cfg.validate()
    assert issues == [], f"TAT/TAC J anchors should pass under Y rule: {issues}"


def test_custom_j_rule_y_rejects_tgg() -> None:
    """Spec point: a `['Y']` J anchor rule rejects `TGG` (W)."""
    cfg = _custom_vj_minimal()
    cfg.set_j_anchor_rule(expected_aa=["Y"])
    cfg.add_j_allele("MYJ-bad*01", "MYJ-bad", b"TGGAAACCC", anchor=0)
    issues = cfg.validate()
    assert any(i["kind"] == "JAnchorUnexpectedAa" and i["aa"] == "W" for i in issues)


def test_custom_v_rule_extends_expected_aa() -> None:
    """A custom V rule accepting both `C` and `H` lets non-Cys but
    still-canonical anchors through."""
    cfg = ge.RefDataConfig.vj()
    cfg.set_v_anchor_rule(expected_aa=["C", "H"])
    # CAT → H.
    cfg.add_v_allele("MYV1*01", "MYV1", b"CATAAACCC", anchor=0)
    cfg.add_j_allele("MYJ1*01", "MYJ1", b"TGGAAA", anchor=0)
    assert cfg.validate() == []


def test_anchor_rule_required_false_suppresses_missing_anchor() -> None:
    """`required=False` lets anchorless alleles pass without
    emitting `MissingAnchor`."""
    cfg = _custom_vj_minimal()
    cfg.set_j_anchor_rule(expected_aa=["W", "F"], required=False)
    cfg.add_j_allele("MYJ-orphan*01", "MYJ-orphan", b"GGG", anchor=None)
    assert cfg.validate() == []


def test_missing_severity_follows_rule() -> None:
    """Spec point: missing anchor severity follows the rule. With
    Fatal severity, AllowCuratable mode cannot opt the issue out."""
    cfg = _custom_vj_minimal()
    cfg.set_j_anchor_rule(expected_aa=["W"], missing_severity="fatal")
    cfg.add_j_allele("MYJ-orphan*01", "MYJ-orphan", b"GGG", anchor=None)
    issues = cfg.validate()
    miss = next((i for i in issues if i["kind"] == "MissingAnchor"), None)
    assert miss is not None
    assert miss["severity"] == "fatal"
    # AllowCuratable cannot rescue a Fatal-tagged issue.
    with pytest.raises(ValueError):
        cfg.validate_with_mode(mode="allow_curatable")


def test_mismatch_severity_fatal_blocks_allow_curatable() -> None:
    cfg = ge.RefDataConfig.vj()
    cfg.set_v_anchor_rule(expected_aa=["C"], mismatch_severity="fatal")
    cfg.add_v_allele("MYV-bad*01", "MYV-bad", b"GGGAAACCC", anchor=0)
    cfg.add_j_allele("MYJ*01", "MYJ", b"TGGAAA", anchor=0)
    with pytest.raises(ValueError):
        cfg.validate_with_mode(mode="allow_curatable")


def test_invalid_byte_remains_fatal_regardless_of_rule() -> None:
    """Spec point: invalid byte stays Fatal regardless of rule
    permissiveness. Structural issues are not rule-controlled."""
    cfg = ge.RefDataConfig.vj()
    # Extremely permissive anchor expectations.
    cfg.set_v_anchor_rule(expected_aa=["A", "C", "D", "E", "F", "G", "H", "I"])
    cfg.set_j_anchor_rule(expected_aa=["A", "C", "D", "E", "F", "G", "H", "I"])
    cfg.add_v_allele("MYV-bad*01", "MYV-bad", b"TGT.AAACC", anchor=0)
    cfg.add_j_allele("MYJ*01", "MYJ", b"TGGAAA", anchor=0)
    issues = cfg.validate()
    bad_byte = next((i for i in issues if i["kind"] == "InvalidAlleleByte"), None)
    assert bad_byte is not None
    assert bad_byte["severity"] == "fatal"


def test_set_allowed_bases_extends_alphabet() -> None:
    cfg = ge.RefDataConfig.vj()
    cfg.set_allowed_bases(["A", "C", "G", "T", "N", "R"])
    assert "R" in cfg.allowed_bases()
    cfg.add_v_allele("MYV*01", "MYV", b"TGTRAA", anchor=0)
    cfg.add_j_allele("MYJ*01", "MYJ", b"TGGAAA", anchor=0)
    issues = cfg.validate()
    assert not any(i["kind"] == "InvalidAlleleByte" for i in issues)


# ──────────────────────────────────────────────────────────────────
# Content hash includes rules
# ──────────────────────────────────────────────────────────────────


def test_refdata_content_hash_changes_when_rules_change() -> None:
    """Spec point: `refdata_content_hash()` changes when rules change."""
    cfg_a = ge.RefDataConfig.vj()
    cfg_a.add_v_allele("MYV*01", "MYV", b"TGTAAA", anchor=0)
    cfg_a.add_j_allele("MYJ*01", "MYJ", b"TGGAAACCC", anchor=0)
    hash_a = cfg_a.content_hash()

    cfg_b = ge.RefDataConfig.vj()
    cfg_b.add_v_allele("MYV*01", "MYV", b"TGTAAA", anchor=0)
    cfg_b.add_j_allele("MYJ*01", "MYJ", b"TGGAAACCC", anchor=0)
    # Same catalogue, different J anchor rule.
    cfg_b.set_j_anchor_rule(expected_aa=["W"])
    hash_b = cfg_b.content_hash()

    assert hash_a != hash_b, "rule change must change the content hash"


def test_refdata_content_hash_changes_with_alphabet_change() -> None:
    cfg_a = ge.RefDataConfig.vj()
    cfg_a.add_v_allele("MYV*01", "MYV", b"TGTAAA", anchor=0)
    cfg_a.add_j_allele("MYJ*01", "MYJ", b"TGGAAA", anchor=0)
    hash_a = cfg_a.content_hash()

    cfg_b = ge.RefDataConfig.vj()
    cfg_b.add_v_allele("MYV*01", "MYV", b"TGTAAA", anchor=0)
    cfg_b.add_j_allele("MYJ*01", "MYJ", b"TGGAAA", anchor=0)
    cfg_b.set_allowed_bases(["A", "C", "G", "T", "N", "R"])
    hash_b = cfg_b.content_hash()

    assert hash_a != hash_b


def test_refdata_content_hash_stable_for_identical_rules() -> None:
    """Two configs with the same catalogue + same rules produce
    identical hashes. Without this, the change-detection above
    would be meaningless."""
    def make() -> "ge.RefDataConfig":
        cfg = ge.RefDataConfig.vj()
        cfg.add_v_allele("MYV*01", "MYV", b"TGTAAA", anchor=0)
        cfg.add_j_allele("MYJ*01", "MYJ", b"TGGAAA", anchor=0)
        cfg.set_j_anchor_rule(expected_aa=["W"])
        return cfg
    assert make().content_hash() == make().content_hash()
