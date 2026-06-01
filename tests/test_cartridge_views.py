"""Cartridge plane views on :class:`DataConfig`.

These tests pin the read-only re-exposure of `DataConfig` fields
under the four cartridge planes (identity, catalogue, rules,
empirical models). No behaviour change in this slice — same pickle,
same checksum, same load path — but docs/tests/users can now talk
about the cartridge model directly on a `DataConfig` instance.
"""
from __future__ import annotations

from dataclasses import FrozenInstanceError

import pytest

from GenAIRR.dataconfig import (
    CartridgeCatalogueView,
    CartridgeIdentityView,
    CartridgeModelsView,
    CartridgeRulesView,
    DataConfig,
)
from GenAIRR._refdata_resolver import _resolve_config_name
from GenAIRR.reference_models import (
    EmpiricalDistributionSpec,
    ReferenceEmpiricalModels,
)
from GenAIRR.reference_rules import (
    AnchorRuleSpec,
    ReferenceRulesSpec,
)


# ──────────────────────────────────────────────────────────────────
# 1. Each view exists on a fresh DataConfig
# ──────────────────────────────────────────────────────────────────


def test_new_dataconfig_exposes_all_four_views() -> None:
    cfg = DataConfig(name="empty")
    assert isinstance(cfg.cartridge_identity, CartridgeIdentityView)
    assert isinstance(cfg.cartridge_catalogue, CartridgeCatalogueView)
    assert isinstance(cfg.cartridge_rules, CartridgeRulesView)
    assert isinstance(cfg.cartridge_models, CartridgeModelsView)


def test_views_are_importable_from_dataconfig_package() -> None:
    """Top-level discoverability — a user can `from GenAIRR.dataconfig
    import CartridgeIdentityView` without reaching into private
    submodules."""
    # Already imported at module load; this assertion just pins the names.
    for cls in (
        CartridgeIdentityView,
        CartridgeCatalogueView,
        CartridgeRulesView,
        CartridgeModelsView,
    ):
        assert cls.__name__.startswith("Cartridge")


# ──────────────────────────────────────────────────────────────────
# 2. Identity / rules / models views are read-only (frozen)
# ──────────────────────────────────────────────────────────────────


@pytest.mark.parametrize(
    "view_factory,attr",
    [
        (lambda c: c.cartridge_identity, "name"),
        (lambda c: c.cartridge_catalogue, "v_alleles"),
        (lambda c: c.cartridge_rules, "reference_rules"),
        (lambda c: c.cartridge_models, "reference_models"),
    ],
)
def test_views_are_frozen(view_factory, attr) -> None:
    cfg = DataConfig(name="frozen_test")
    view = view_factory(cfg)
    with pytest.raises(FrozenInstanceError):
        setattr(view, attr, "tampered")


# ──────────────────────────────────────────────────────────────────
# 3. Identity view surfaces name + metadata
# ──────────────────────────────────────────────────────────────────


def test_identity_view_carries_name_and_metadata() -> None:
    cfg = _resolve_config_name("human_igh")
    view = cfg.cartridge_identity
    assert view.name == cfg.name
    assert view.metadata is cfg.metadata


def test_identity_view_handles_missing_metadata() -> None:
    cfg = DataConfig(name="bare")
    view = cfg.cartridge_identity
    assert view.name == "bare"
    assert view.metadata is None


# ──────────────────────────────────────────────────────────────────
# 4. Catalogue view counts match the legacy `number_of_*_alleles`
# ──────────────────────────────────────────────────────────────────


def test_catalogue_view_counts_match_legacy_properties() -> None:
    cfg = _resolve_config_name("human_igh")
    view = cfg.cartridge_catalogue
    assert view.number_of_v_alleles == cfg.number_of_v_alleles
    assert view.number_of_d_alleles == cfg.number_of_d_alleles
    assert view.number_of_j_alleles == cfg.number_of_j_alleles
    assert view.number_of_c_alleles == cfg.number_of_c_alleles


def test_catalogue_view_references_same_underlying_dicts() -> None:
    """Views are zero-copy snapshots — they reference the same
    storage the DataConfig owns. Mutating the underlying storage
    isn't supported (and shouldn't happen in production), but the
    identity guarantee makes the views genuinely cheap."""
    cfg = _resolve_config_name("human_igh")
    view = cfg.cartridge_catalogue
    assert view.v_alleles is cfg.v_alleles
    assert view.d_alleles is cfg.d_alleles
    assert view.j_alleles is cfg.j_alleles
    assert view.c_alleles is cfg.c_alleles


def test_catalogue_view_returns_zero_counts_for_empty_pools() -> None:
    cfg = DataConfig(name="empty_pools")
    view = cfg.cartridge_catalogue
    assert view.number_of_v_alleles == 0
    assert view.number_of_d_alleles == 0
    assert view.number_of_j_alleles == 0
    assert view.number_of_c_alleles == 0


# ──────────────────────────────────────────────────────────────────
# 5. Rules view — legacy + modern + None paths
# ──────────────────────────────────────────────────────────────────


def test_rules_view_reflects_attached_reference_rules() -> None:
    cfg = DataConfig(name="rules_test")
    spec = ReferenceRulesSpec(j_anchor=AnchorRuleSpec(["Y"]))
    cfg.reference_rules = spec
    assert cfg.cartridge_rules.reference_rules is spec


def test_rules_view_is_none_for_default_dataconfig() -> None:
    cfg = DataConfig(name="default")
    assert cfg.cartridge_rules.reference_rules is None


def test_rules_view_is_none_for_legacy_pickled_object() -> None:
    """Legacy pickles don't have ``reference_rules`` in
    ``__dict__``. The view must return ``None`` (via the
    ``__getattr__`` shim) without raising."""
    cfg = DataConfig(name="legacy")
    cfg.__dict__.pop("reference_rules", None)
    assert cfg.cartridge_rules.reference_rules is None


# ──────────────────────────────────────────────────────────────────
# 6. Models view exposes both new + legacy surfaces
# ──────────────────────────────────────────────────────────────────


def test_models_view_exposes_typed_spec_when_attached() -> None:
    cfg = DataConfig(name="models_test")
    spec = ReferenceEmpiricalModels(
        np_lengths={"NP1": EmpiricalDistributionSpec([(3, 1.0)])},
    )
    cfg.reference_models = spec
    view = cfg.cartridge_models
    assert view.reference_models is spec


def test_models_view_exposes_legacy_dicts() -> None:
    cfg = _resolve_config_name("human_igh")
    view = cfg.cartridge_models
    # Bundled human_igh carries the legacy nested-dict shapes.
    assert isinstance(view.legacy_np_lengths, dict)
    assert isinstance(view.legacy_trim_dicts, dict)
    assert "NP1" in view.legacy_np_lengths
    assert "V_3" in view.legacy_trim_dicts


def test_models_view_legacy_dicts_are_same_storage() -> None:
    cfg = _resolve_config_name("human_igh")
    view = cfg.cartridge_models
    assert view.legacy_np_lengths is cfg.NP_lengths
    assert view.legacy_trim_dicts is cfg.trim_dicts


def test_models_view_reference_models_is_none_on_legacy_object() -> None:
    cfg = DataConfig(name="legacy_models")
    cfg.__dict__.pop("reference_models", None)
    assert cfg.cartridge_models.reference_models is None


def test_models_view_shows_both_surfaces_simultaneously() -> None:
    """A cartridge can carry both — the typed spec wins at
    extraction, but the view shows both shapes side-by-side so a
    user can see exactly what's authored vs what's legacy.

    Copy before mutating: ``_resolve_config_name`` returns a cached
    module-level instance, so attaching a spec to the live object
    leaks state into other tests.
    """
    cfg = _resolve_config_name("human_igh").copy()
    cfg.reference_models = ReferenceEmpiricalModels(
        np_lengths={"NP1": EmpiricalDistributionSpec([(7, 1.0)])},
    )
    view = cfg.cartridge_models
    assert view.reference_models is cfg.reference_models
    assert view.legacy_np_lengths is cfg.NP_lengths


# ──────────────────────────────────────────────────────────────────
# 7. Views do not change pickle / checksum behaviour
# ──────────────────────────────────────────────────────────────────


def test_view_access_does_not_change_checksum() -> None:
    """Constructing views must be a pure read — the underlying
    checksum stays identical."""
    cfg = DataConfig(name="checksum_stability")
    h_before = cfg.compute_checksum()
    _ = cfg.cartridge_identity
    _ = cfg.cartridge_catalogue
    _ = cfg.cartridge_rules
    _ = cfg.cartridge_models
    h_after = cfg.compute_checksum()
    assert h_before == h_after
