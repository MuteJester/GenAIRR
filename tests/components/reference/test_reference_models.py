"""Reference Empirical Models v1 — Python cartridge defaults plane.

DataConfig now carries an optional :class:`ReferenceEmpiricalModels`
bundle (NP-length + trim distributions). The extractor consumes it
first; falls back to the legacy ``NP_lengths`` / ``trim_dicts``
nested-dict path; finally falls through to the uniform placeholder
when neither is present.

Shape-only validation lives Python-side. The Rust ``RefDataConfig``
does NOT carry empirical models in v1 — the engine consumes
already-lowered pass distributions. A later slice may surface a
normalized summary into the Rust cartridge for trace identity.
"""
from __future__ import annotations

import copy

import pytest

from GenAIRR._dataconfig_extract import extract_recombine_defaults
from GenAIRR._refdata_resolver import _resolve_config_name
from GenAIRR.dataconfig.data_config import DataConfig
from GenAIRR.reference_models import (
    EmpiricalDistributionSpec,
    ReferenceEmpiricalModels,
)


# ──────────────────────────────────────────────────────────────────
# Spec shape validation (pure Python)
# ──────────────────────────────────────────────────────────────────


def test_empirical_distribution_spec_defaults_validate() -> None:
    EmpiricalDistributionSpec([(0, 1.0), (1, 2.0), (2, 0.5)]).validate("NP1")


def test_empirical_distribution_spec_rejects_empty() -> None:
    with pytest.raises(ValueError, match="non-empty"):
        EmpiricalDistributionSpec([]).validate("NP1")


def test_empirical_distribution_spec_rejects_negative_value() -> None:
    with pytest.raises(ValueError, match="non-negative"):
        EmpiricalDistributionSpec([(-1, 1.0)]).validate("V_3")


def test_empirical_distribution_spec_rejects_zero_weight() -> None:
    with pytest.raises(ValueError, match="> 0"):
        EmpiricalDistributionSpec([(0, 0.0)]).validate("V_3")


def test_empirical_distribution_spec_rejects_negative_weight() -> None:
    with pytest.raises(ValueError, match="> 0"):
        EmpiricalDistributionSpec([(0, -1.0)]).validate("V_3")


def test_empirical_distribution_spec_rejects_nan_weight() -> None:
    with pytest.raises(ValueError, match="finite"):
        EmpiricalDistributionSpec([(0, float("nan"))]).validate("NP1")


def test_empirical_distribution_spec_rejects_inf_weight() -> None:
    with pytest.raises(ValueError, match="finite"):
        EmpiricalDistributionSpec([(0, float("inf"))]).validate("NP1")


def test_empirical_distribution_spec_rejects_non_int_value() -> None:
    with pytest.raises(ValueError, match="non-negative int"):
        EmpiricalDistributionSpec([(1.5, 1.0)]).validate("NP1")  # type: ignore[arg-type]


def test_empirical_distribution_spec_rejects_bool_value() -> None:
    # Booleans are technically int subclasses in Python; rejecting
    # them keeps the surface from accepting True / False as "values".
    with pytest.raises(ValueError, match="non-negative int"):
        EmpiricalDistributionSpec([(True, 1.0)]).validate("NP1")  # type: ignore[list-item]


def test_reference_empirical_models_defaults_validate() -> None:
    ReferenceEmpiricalModels().validate(chain_type="vdj")
    ReferenceEmpiricalModels().validate(chain_type="vj")


def test_reference_empirical_models_rejects_unknown_np_key() -> None:
    rm = ReferenceEmpiricalModels(
        np_lengths={"NP3": EmpiricalDistributionSpec([(0, 1.0)])},
    )
    with pytest.raises(ValueError, match="np_lengths key"):
        rm.validate(chain_type="vdj")


def test_reference_empirical_models_rejects_unknown_trim_key() -> None:
    rm = ReferenceEmpiricalModels(
        trims={"V_5": EmpiricalDistributionSpec([(0, 1.0)])},
    )
    with pytest.raises(ValueError, match="trims key"):
        rm.validate(chain_type="vdj")


def test_reference_empirical_models_rejects_d_trim_on_vj_chain() -> None:
    rm = ReferenceEmpiricalModels(
        trims={"D_5": EmpiricalDistributionSpec([(0, 1.0)])},
    )
    with pytest.raises(ValueError, match="VJ chain"):
        rm.validate(chain_type="vj")


def test_reference_empirical_models_accepts_d_trim_on_vdj_chain() -> None:
    rm = ReferenceEmpiricalModels(
        trims={
            "V_3": EmpiricalDistributionSpec([(0, 1.0)]),
            "D_5": EmpiricalDistributionSpec([(0, 1.0)]),
            "D_3": EmpiricalDistributionSpec([(0, 1.0)]),
            "J_5": EmpiricalDistributionSpec([(0, 1.0)]),
        },
    )
    rm.validate(chain_type="vdj")


# ──────────────────────────────────────────────────────────────────
# DataConfig integration — field + backcompat shim + checksum
# ──────────────────────────────────────────────────────────────────


def test_legacy_dataconfig_returns_none_for_reference_models() -> None:
    cfg = DataConfig(name="legacy")
    cfg.__dict__.pop("reference_models", None)
    assert cfg.reference_models is None


def test_default_new_dataconfig_has_none_reference_models() -> None:
    cfg = DataConfig(name="modern")
    assert cfg.reference_models is None


def test_default_none_models_preserves_checksum_vs_legacy() -> None:
    """A default ``reference_models=None`` instance must hash the
    same as a legacy pickle without the key (soft-transition)."""
    cfg_new = DataConfig(name="dummy")
    cfg_legacy = DataConfig(name="dummy")
    cfg_legacy.__dict__.pop("reference_models", None)
    assert cfg_new.compute_checksum() == cfg_legacy.compute_checksum()


def test_setting_reference_models_changes_checksum() -> None:
    cfg = DataConfig(name="dummy")
    h_none = cfg.compute_checksum()
    cfg.reference_models = ReferenceEmpiricalModels(
        np_lengths={"NP1": EmpiricalDistributionSpec([(3, 1.0)])},
    )
    assert cfg.compute_checksum() != h_none


def test_two_equal_reference_models_share_checksum() -> None:
    rm = lambda: ReferenceEmpiricalModels(
        np_lengths={"NP1": EmpiricalDistributionSpec([(3, 1.0)])},
    )
    cfg_a = DataConfig(name="dup", reference_models=rm())
    cfg_b = DataConfig(name="dup", reference_models=rm())
    assert cfg_a.compute_checksum() == cfg_b.compute_checksum()


# ──────────────────────────────────────────────────────────────────
# Extractor — reference_models > legacy NP_lengths > uniform
# ──────────────────────────────────────────────────────────────────


def _human_igh() -> DataConfig:
    """Deepcopy of the bundled IGH DataConfig — exposes its legacy
    NP_lengths / trim_dicts so the legacy-path tests still have
    meaningful empirical data to fall back to.
    """
    return _resolve_config_name("human_igh").copy()


def test_explicit_np1_model_overrides_legacy_extraction() -> None:
    cfg = _human_igh()
    cfg.reference_models = ReferenceEmpiricalModels(
        np_lengths={"NP1": EmpiricalDistributionSpec([(7, 1.0)])},
    )
    defaults = extract_recombine_defaults(cfg)
    # Explicit model lowered into the engine's flat pair-list shape;
    # NOT a marginalised approximation of the legacy nested dict.
    assert defaults["np1"] == [(7, 1.0)]


def test_explicit_trim_v3_model_overrides_legacy_extraction() -> None:
    cfg = _human_igh()
    cfg.reference_models = ReferenceEmpiricalModels(
        trims={"V_3": EmpiricalDistributionSpec([(2, 0.7), (5, 0.3)])},
    )
    defaults = extract_recombine_defaults(cfg)
    pairs = defaults["trim_v_3"]
    assert pairs is not None
    # Cap from shortest V allele still applies, but the values are
    # small enough that none are clamped. Order is sorted.
    assert pairs[0][0] == 2
    assert any(v == 5 for v, _ in pairs)


def test_missing_explicit_key_falls_back_to_legacy() -> None:
    """``reference_models`` may be partial — keys it doesn't carry
    fall back to the legacy extraction path. Both can coexist."""
    cfg = _human_igh()
    cfg.reference_models = ReferenceEmpiricalModels(
        np_lengths={"NP1": EmpiricalDistributionSpec([(11, 1.0)])},
    )
    defaults = extract_recombine_defaults(cfg)
    # NP1: explicit.
    assert defaults["np1"] == [(11, 1.0)]
    # NP2 unauthored: legacy IGH NP_lengths extraction still runs.
    assert defaults["np2"] is not None and len(defaults["np2"]) > 0


def test_no_reference_models_keeps_legacy_extraction() -> None:
    cfg = _human_igh()
    assert cfg.reference_models is None
    defaults = extract_recombine_defaults(cfg)
    assert defaults["np1"] is not None
    assert defaults["trim_v_3"] is not None


def test_extractor_rejects_invalid_reference_models() -> None:
    cfg = _human_igh()
    cfg.reference_models = ReferenceEmpiricalModels(
        np_lengths={"NP1": EmpiricalDistributionSpec([(0, -1.0)])},  # bad weight
    )
    with pytest.raises(ValueError, match="reference_models failed shape validation"):
        extract_recombine_defaults(cfg)


def test_extractor_rejects_d_trim_on_vj_metadata() -> None:
    from GenAIRR.dataconfig.config_info import ConfigInfo
    from GenAIRR.dataconfig.enums import ChainType, Species
    from datetime import date

    cfg = DataConfig(
        name="dummy_vj",
        metadata=ConfigInfo(
            species=Species.HUMAN,
            chain_type=ChainType.BCR_LIGHT_KAPPA,
            reference_set="test",
            last_updated=date(2024, 1, 1),
            has_d=False,
        ),
        reference_models=ReferenceEmpiricalModels(
            trims={"D_5": EmpiricalDistributionSpec([(1, 1.0)])},
        ),
    )
    with pytest.raises(ValueError, match="VJ chain"):
        extract_recombine_defaults(cfg)


def test_trim_cap_applies_to_explicit_model() -> None:
    """An explicit model can't sample a trim larger than the pool's
    shortest allele permits — the cap clamps both legacy and explicit
    paths so recombination stays viable."""
    cfg = _human_igh()
    # Bundled IGH V alleles are ~300 bp; cap is ~150 (half is N/A,
    # one-sided V_3). Authoring a 9999-bp trim must be filtered.
    cfg.reference_models = ReferenceEmpiricalModels(
        trims={"V_3": EmpiricalDistributionSpec([(2, 1.0), (9999, 1.0)])},
    )
    defaults = extract_recombine_defaults(cfg)
    pairs = defaults["trim_v_3"]
    assert pairs is not None
    assert all(v != 9999 for v, _ in pairs)


# ──────────────────────────────────────────────────────────────────
# End-to-end — recombine uses the explicit model
# ──────────────────────────────────────────────────────────────────


def test_recombine_uses_explicit_np_model_end_to_end() -> None:
    """Build an Experiment from a DataConfig carrying an explicit
    NP1 model and confirm the recombine step lowers that model into
    the compiled pass plan (no legacy extraction in between)."""
    import GenAIRR as ga

    cfg = _human_igh()
    cfg.reference_models = ReferenceEmpiricalModels(
        np_lengths={"NP1": EmpiricalDistributionSpec([(13, 1.0)])},
    )
    exp = ga.Experiment.on(cfg).recombine()
    # The recombine step stores the lowered np1_lengths on its
    # internal step descriptor; inspecting the experiment confirms
    # the explicit model survived the cascade.
    recombine_steps = [s for s in exp._steps if type(s).__name__ == "_RecombineStep"]
    assert recombine_steps, "recombine() must have appended exactly one _RecombineStep"
    assert recombine_steps[0].np1_lengths == ((13, 1.0),)
