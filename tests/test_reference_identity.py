"""Reference Identity tests — cartridge self-description.

The cartridge now carries species/locus/reference_set/name/source
identity in Rust ``RefDataConfig.identity``. The Python loader stamps
it from ``DataConfig.metadata`` + ``DataConfig.name``; user code can
also set it directly via ``cfg.set_identity(...)``.

Identity drives:
  - the cartridge content hash (two configs with identical catalogues
    but different declared identity are different cartridges),
  - the locus-fallback cascade for J-anchor expectations when no
    explicit ``ReferenceRulesSpec`` is attached (cartridge.identity
    speaks first, then allele-name inference),
  - a Fatal ``LocusChainTypeMismatch`` validation issue when the
    declared locus disagrees with the chain topology.
"""
from __future__ import annotations

import pytest

import GenAIRR as ga
from GenAIRR import _engine as ge


# ──────────────────────────────────────────────────────────────────
# Direct API — set / read identity on engine-native RefDataConfig
# ──────────────────────────────────────────────────────────────────


def test_default_identity_is_all_none() -> None:
    cfg = ge.RefDataConfig.vj()
    id_ = cfg.identity()
    assert id_ == {
        "species": None,
        "locus": None,
        "reference_set": None,
        "name": None,
        "source": None,
    }


def test_set_identity_populates_one_field() -> None:
    cfg = ge.RefDataConfig.vj()
    cfg.set_identity(locus="IGL")
    assert cfg.identity()["locus"] == "IGL"
    # Other fields stay None.
    assert cfg.identity()["species"] is None


def test_set_identity_leaves_unchanged_fields_alone() -> None:
    cfg = ge.RefDataConfig.vj()
    cfg.set_identity(species="Llama", locus="IGL", reference_set="custom")
    cfg.set_identity(name="my_llama_kappa")  # only updates name
    id_ = cfg.identity()
    assert id_["species"] == "Llama"
    assert id_["locus"] == "IGL"
    assert id_["reference_set"] == "custom"
    assert id_["name"] == "my_llama_kappa"


def test_set_identity_empty_string_clears_field() -> None:
    cfg = ge.RefDataConfig.vj()
    cfg.set_identity(locus="IGK")
    cfg.set_identity(locus="")  # empty clears
    assert cfg.identity()["locus"] is None


# ──────────────────────────────────────────────────────────────────
# Bridge — DataConfig metadata → ReferenceIdentity
# ──────────────────────────────────────────────────────────────────


@pytest.mark.parametrize(
    "preset,expected_locus,expected_species",
    [
        ("human_igh", "IGH", "Human"),
        ("human_igk", "IGK", "Human"),
        ("human_igl", "IGL", "Human"),
        ("human_tcrb", "TRB", "Human"),
        ("mouse_igh", "IGH", "Mouse"),
    ],
)
def test_bundled_preset_carries_identity(
    preset: str, expected_locus: str, expected_species: str,
) -> None:
    rd = ga.Experiment.on(preset).refdata
    id_ = rd.identity()
    assert id_["locus"] == expected_locus
    assert id_["species"] == expected_species
    assert id_["reference_set"]  # present, non-empty
    assert id_["source"] == "DataConfig"


def test_bundled_preset_identity_carries_reference_set_when_metadata_has_it() -> None:
    """Bundled IGH/IGK/IGL all have a populated reference_set. The
    bridge must transfer it verbatim."""
    for preset in ("human_igh", "human_igk", "human_igl"):
        rs = ga.Experiment.on(preset).refdata.identity()["reference_set"]
        assert isinstance(rs, str) and rs


# ──────────────────────────────────────────────────────────────────
# Content hash — identity is part of cartridge identity
# ──────────────────────────────────────────────────────────────────


def _minimal_vj_cartridge() -> "ge.RefDataConfig":
    cfg = ge.RefDataConfig.vj()
    cfg.add_v_allele("MYV1*01", "MYV1", b"TGTAAACCC", anchor=0)
    cfg.add_j_allele("MYJ1*01", "MYJ1", b"TGGAAACCC", anchor=0)
    return cfg


def test_content_hash_changes_when_identity_changes() -> None:
    cfg_a = _minimal_vj_cartridge()
    cfg_b = _minimal_vj_cartridge()
    cfg_b.set_identity(species="Human", locus="IGK", reference_set="OGRDB")
    assert cfg_a.content_hash() != cfg_b.content_hash()


def test_content_hash_stable_for_same_identity() -> None:
    cfg_a = _minimal_vj_cartridge()
    cfg_b = _minimal_vj_cartridge()
    for cfg in (cfg_a, cfg_b):
        cfg.set_identity(species="Human", locus="IGK", reference_set="OGRDB")
    assert cfg_a.content_hash() == cfg_b.content_hash()


def test_content_hash_sensitive_to_locus_only() -> None:
    """Hash differs when ONLY the locus identity field differs."""
    cfg_a = _minimal_vj_cartridge()
    cfg_b = _minimal_vj_cartridge()
    cfg_a.set_identity(locus="IGK")
    cfg_b.set_identity(locus="IGL")
    assert cfg_a.content_hash() != cfg_b.content_hash()


# ──────────────────────────────────────────────────────────────────
# Locus cascade — identity wins over allele-name inference
# ──────────────────────────────────────────────────────────────────


def test_identity_locus_drives_rule_fallback_when_no_spec() -> None:
    """``dataconfig_to_refdata`` consults ``identity.locus`` BEFORE
    allele-name inference. Loading human_igk via the bridge sets
    identity.locus=IGK, so the J rule narrows to ['F'] even though
    the test catalogue's allele names don't necessarily round-trip
    through the inference heuristic.
    """
    rd = ga.Experiment.on("human_igk").refdata
    assert rd.identity()["locus"] == "IGK"
    assert rd.j_anchor_rule()["expected_aa"] == ["F"]


def test_identity_locus_present_even_if_no_v_alleles_inferred() -> None:
    """When the loader can stamp identity from metadata, the
    cartridge should not need allele-name heuristics to pick a rule."""
    rd = ga.Experiment.on("human_igh").refdata
    assert rd.identity()["locus"] == "IGH"
    assert rd.j_anchor_rule()["expected_aa"] == ["W"]


# ──────────────────────────────────────────────────────────────────
# Validation — LocusChainTypeMismatch
# ──────────────────────────────────────────────────────────────────


def test_vj_chain_with_igh_locus_fails_validation() -> None:
    """``locus="IGH"`` on ``RefDataConfig.vj()`` is structurally
    inconsistent and must fail validation."""
    cfg = _minimal_vj_cartridge()
    cfg.set_identity(locus="IGH")
    issues = cfg.validate()
    mismatches = [i for i in issues if i["kind"] == "LocusChainTypeMismatch"]
    assert mismatches, f"expected LocusChainTypeMismatch, got: {issues}"
    assert mismatches[0]["locus"] == "IGH"
    assert mismatches[0]["chain_type"] == "vj"
    assert mismatches[0]["severity"] == "fatal"


def test_vdj_chain_with_igk_locus_fails_validation() -> None:
    cfg = ge.RefDataConfig.vdj()
    cfg.add_v_allele("MYV1*01", "MYV1", b"TGTAAACCC", anchor=0)
    cfg.add_d_allele("MYD1*01", "MYD1", b"GGGCCC")
    cfg.add_j_allele("MYJ1*01", "MYJ1", b"TGGAAACCC", anchor=0)
    cfg.set_identity(locus="IGK")
    issues = cfg.validate()
    assert any(i["kind"] == "LocusChainTypeMismatch" for i in issues)


def test_unknown_locus_does_not_fail_validation() -> None:
    cfg = _minimal_vj_cartridge()
    cfg.set_identity(locus="XYZ")  # unknown — no expectation enforced
    assert not any(
        i["kind"] == "LocusChainTypeMismatch" for i in cfg.validate()
    )


def test_no_locus_does_not_fail_validation() -> None:
    cfg = _minimal_vj_cartridge()
    assert not any(
        i["kind"] == "LocusChainTypeMismatch" for i in cfg.validate()
    )


def test_locus_chain_mismatch_blocks_strict_compile() -> None:
    """End-to-end: a cartridge with a declared-locus / chain-topology
    mismatch must fail ``Experiment.compile()`` with a clear message
    that the LocusChainTypeMismatch is fatal."""
    cfg = _minimal_vj_cartridge()
    cfg.set_identity(locus="IGH")
    with pytest.raises(ValueError) as exc:
        ga.Experiment.on(cfg).recombine().compile()
    msg = str(exc.value)
    assert "[fatal]" in msg
    assert "IGH" in msg
    assert "chain_type" in msg or "Vj" in msg
