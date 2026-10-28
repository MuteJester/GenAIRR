"""DataConfig ↔ ReferenceRules bridge tests.

The slice extends Python ``DataConfig`` with an optional
:class:`ReferenceRulesSpec`. ``dataconfig_to_refdata`` transfers it
verbatim into the Rust ``RefDataConfig.rules`` slice; the spec
shape-checks before crossing the PyO3 boundary, and the absence
(``None``) preserves the bundled-locus inference path so legacy
pickles continue to load with their old checksum.
"""
from __future__ import annotations

import pytest

import GenAIRR as ga
from GenAIRR import _engine as ge
from GenAIRR._refdata_resolver import dataconfig_to_refdata
from GenAIRR.dataconfig.data_config import DataConfig
from GenAIRR.reference_rules import AnchorRuleSpec, ReferenceRulesSpec


# ──────────────────────────────────────────────────────────────────
# Fixture helpers
# ──────────────────────────────────────────────────────────────────


def _bundled_igh_with_rules(spec=None) -> DataConfig:
    """Load the bundled IGH DataConfig and attach ``spec`` (which may
    be ``None``). Builtin DataConfigs are full, validated catalogues
    with real V/D/J alleles — perfect to exercise the bridge end-to-end
    without hand-rolling Allele instances.

    We deepcopy so each test gets a fresh, mutable cfg. Mutation
    persists across calls without the copy because
    ``_resolve_config_name`` returns the cached module-level instance.
    """
    from GenAIRR._refdata_resolver import _resolve_config_name
    cfg = _resolve_config_name("human_igh").copy()
    cfg.reference_rules = spec
    return cfg


# ──────────────────────────────────────────────────────────────────
# Spec validation (Python-side shape checks)
# ──────────────────────────────────────────────────────────────────


def test_anchor_rule_spec_defaults_validate() -> None:
    AnchorRuleSpec(["C"]).validate("V")
    AnchorRuleSpec(["W", "F"]).validate("J")


def test_anchor_rule_spec_rejects_multi_char_aa() -> None:
    with pytest.raises(ValueError, match="one-character"):
        AnchorRuleSpec(["XYZ"]).validate("V")


def test_anchor_rule_spec_rejects_unknown_aa_letter() -> None:
    with pytest.raises(ValueError, match="recognised amino-acid"):
        AnchorRuleSpec(["Z"]).validate("V")


def test_anchor_rule_spec_rejects_bad_severity() -> None:
    with pytest.raises(ValueError, match="missing_severity"):
        AnchorRuleSpec(["C"], missing_severity="warning").validate("V")
    with pytest.raises(ValueError, match="mismatch_severity"):
        AnchorRuleSpec(["C"], mismatch_severity="info").validate("V")


def test_anchor_rule_spec_required_true_needs_nonempty_aa() -> None:
    with pytest.raises(ValueError, match="non-empty"):
        AnchorRuleSpec([], required=True).validate("J")


def test_anchor_rule_spec_required_false_allows_empty_aa() -> None:
    AnchorRuleSpec([], required=False).validate("J")


def test_reference_rules_spec_defaults_validate() -> None:
    ReferenceRulesSpec().validate()


def test_reference_rules_spec_rejects_missing_canonical_base() -> None:
    with pytest.raises(ValueError, match="canonical DNA bases"):
        ReferenceRulesSpec(allowed_bases=["A", "C", "G", "N"]).validate()


def test_reference_rules_spec_rejects_non_letter_alphabet() -> None:
    with pytest.raises(ValueError, match="ASCII letters"):
        ReferenceRulesSpec(allowed_bases=["A", "C", "G", "T", "N", "."]).validate()


def test_reference_rules_spec_rejects_multi_char_base() -> None:
    with pytest.raises(ValueError, match="single-character"):
        ReferenceRulesSpec(allowed_bases=["A", "C", "G", "T", "NN"]).validate()


# ──────────────────────────────────────────────────────────────────
# DataConfig integration — field + backcompat shim + checksum
# ──────────────────────────────────────────────────────────────────


def test_dataconfig_legacy_object_returns_none_for_reference_rules() -> None:
    """Pickled DataConfigs predating this field don't have the key in
    ``__dict__``. ``__getattr__`` shim must return None so older
    consumers stay happy."""
    cfg = DataConfig(name="legacy")
    # Simulate legacy pickle: key absent from instance dict entirely.
    cfg.__dict__.pop("reference_rules", None)
    assert cfg.reference_rules is None


def test_dataconfig_default_reference_rules_is_none() -> None:
    cfg = DataConfig(name="modern")
    assert cfg.reference_rules is None


def test_default_none_preserves_legacy_checksum() -> None:
    """A DataConfig with default ``reference_rules=None`` must produce
    the same checksum it would have before the field was added.
    Simulate the legacy by popping the key from ``__dict__``."""
    cfg_new = DataConfig(name="dummy")  # has `reference_rules=None` in __dict__
    cfg_legacy = DataConfig(name="dummy")
    cfg_legacy.__dict__.pop("reference_rules", None)
    assert cfg_new.compute_checksum() == cfg_legacy.compute_checksum()


def test_setting_reference_rules_changes_dataconfig_checksum() -> None:
    """Non-None ``reference_rules`` is part of the cartridge identity
    and must change the checksum."""
    base = DataConfig(name="dummy")
    base_hash = base.compute_checksum()
    base.reference_rules = ReferenceRulesSpec(
        j_anchor=AnchorRuleSpec(["Y"]),
    )
    assert base.compute_checksum() != base_hash


def test_compute_checksum_remains_deterministic_with_rules() -> None:
    spec = ReferenceRulesSpec(j_anchor=AnchorRuleSpec(["W"]))
    cfg_a = DataConfig(name="dup", reference_rules=spec)
    cfg_b = DataConfig(name="dup", reference_rules=ReferenceRulesSpec(
        j_anchor=AnchorRuleSpec(["W"])
    ))
    assert cfg_a.compute_checksum() == cfg_b.compute_checksum()


# ──────────────────────────────────────────────────────────────────
# Bridge — dataconfig_to_refdata applies the spec
# ──────────────────────────────────────────────────────────────────


def test_bridge_transfers_custom_j_y_rule() -> None:
    spec = ReferenceRulesSpec(j_anchor=AnchorRuleSpec(["Y"]))
    cfg = _bundled_igh_with_rules(spec=spec)
    rd = dataconfig_to_refdata(cfg)
    assert rd.j_anchor_rule()["expected_aa"] == ["Y"]
    assert rd.v_anchor_rule()["expected_aa"] == ["C"]


def test_bridge_transfers_custom_alphabet() -> None:
    spec = ReferenceRulesSpec(allowed_bases=["A", "C", "G", "T", "N", "R"])
    cfg = _bundled_igh_with_rules(spec=spec)
    rd = dataconfig_to_refdata(cfg)
    assert "R" in rd.allowed_bases()


def test_bridge_transfers_severity_settings() -> None:
    spec = ReferenceRulesSpec(
        j_anchor=AnchorRuleSpec(
            ["W"],
            missing_severity="fatal",
            mismatch_severity="fatal",
        ),
    )
    cfg = _bundled_igh_with_rules(spec=spec)
    rd = dataconfig_to_refdata(cfg)
    j = rd.j_anchor_rule()
    assert j["missing_severity"] == "fatal"
    assert j["mismatch_severity"] == "fatal"


def test_bridge_alphabet_via_extended_set_compiles_through_compile_gate() -> None:
    """End-to-end: cfg with an extended-alphabet spec compiles cleanly
    through ``Experiment.compile()`` — bundled bases are all A/C/G/T/N
    so the added ``R`` doesn't matter to the catalogue, but the
    rule transfer must not block compile."""
    spec = ReferenceRulesSpec(allowed_bases=["A", "C", "G", "T", "N", "R"])
    cfg = _bundled_igh_with_rules(spec=spec)
    ga.Experiment.on(cfg).recombine().compile()


def test_bridge_without_spec_keeps_locus_default() -> None:
    """Spec point: no ``reference_rules`` preserves the current
    bundled-locus default. For the IGH-named V allele below, the
    locus inference narrows J to ['W']."""
    cfg = _bundled_igh_with_rules(spec=None)
    rd = dataconfig_to_refdata(cfg)
    assert rd.j_anchor_rule()["expected_aa"] == ["W"]
    assert rd.v_anchor_rule()["expected_aa"] == ["C"]


def test_bridge_spec_overrides_locus_default() -> None:
    """When a user attaches a spec, it must win over locus inference
    even on a bundled-named catalogue."""
    spec = ReferenceRulesSpec(j_anchor=AnchorRuleSpec(["F"]))  # IGH would default to ['W']
    cfg = _bundled_igh_with_rules(spec=spec)
    rd = dataconfig_to_refdata(cfg)
    assert rd.j_anchor_rule()["expected_aa"] == ["F"]


def test_bridge_rejects_invalid_spec_before_pyo3_crossing() -> None:
    """Spec point: invalid rule shape raises a clean Python
    ``ValueError`` — not a Rust panic or PyO3 conversion error."""
    spec = ReferenceRulesSpec(j_anchor=AnchorRuleSpec(["Z"]))  # Z is not a valid AA
    cfg = _bundled_igh_with_rules(spec=spec)
    with pytest.raises(ValueError, match="recognised amino-acid"):
        dataconfig_to_refdata(cfg)


def test_bridge_rejects_invalid_severity_with_clear_message() -> None:
    """Spec point: invalid rule severity raises a clear ValueError."""
    spec = ReferenceRulesSpec(
        j_anchor=AnchorRuleSpec(["W"], missing_severity="bogus"),
    )
    cfg = _bundled_igh_with_rules(spec=spec)
    with pytest.raises(ValueError, match="missing_severity"):
        dataconfig_to_refdata(cfg)


# ──────────────────────────────────────────────────────────────────
# Cartridge identity — hash changes when rules differ
# ──────────────────────────────────────────────────────────────────


def test_bridge_refdata_hash_differs_when_rules_differ() -> None:
    """Two cartridges with the same catalogue but different spec
    produce different refdata content hashes."""
    cfg_a = _bundled_igh_with_rules(
        spec=ReferenceRulesSpec(j_anchor=AnchorRuleSpec(["W"])),
    )
    cfg_b = _bundled_igh_with_rules(
        spec=ReferenceRulesSpec(j_anchor=AnchorRuleSpec(["F"])),
    )
    assert dataconfig_to_refdata(cfg_a).content_hash() != dataconfig_to_refdata(cfg_b).content_hash()


def test_bridge_refdata_hash_stable_when_rules_match() -> None:
    spec_a = ReferenceRulesSpec(j_anchor=AnchorRuleSpec(["W"]))
    spec_b = ReferenceRulesSpec(j_anchor=AnchorRuleSpec(["W"]))
    cfg_a = _bundled_igh_with_rules(spec=spec_a)
    cfg_b = _bundled_igh_with_rules(spec=spec_b)
    assert dataconfig_to_refdata(cfg_a).content_hash() == dataconfig_to_refdata(cfg_b).content_hash()
