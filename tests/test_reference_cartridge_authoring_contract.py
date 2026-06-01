"""Contract pins for the Reference Cartridge Authoring /
Inference API audit.

Companion to
[`docs/reference_cartridge_authoring_audit.md`](../docs/reference_cartridge_authoring_audit.md).

Pin set:

- ``pin_scaffold_*`` — live inference-adjacent helpers the
  new builder will reuse (`parse_fasta`,
  `compute_v_region_boundaries`, `cartridge_manifest`,
  `verify_integrity`, `RefDataConfig.{vj,vdj}` +
  `add_*_allele`, `ReferenceEmpiricalModels`,
  `ReferenceRulesSpec`, `dataconfig_to_refdata` bridge,
  bundled-pickle loader).
- ``pin_present_*`` — dead-reference inventory
  (`DataConfig.build_report` field exists but is unconditionally
  `None`; docstring still names a removed class; private
  build script imports a dead module).
- ``pin_absence_*`` — the gaps the implementation slice
  closes (no public `ReferenceCartridgeBuilder`, no
  `from_fasta` / `from_airr` constructors, no `infer_*` /
  `estimate_*` public step methods, no
  `CartridgeBuildReport` dataclass).

**Pre-flight verdict (audit §3): clean-yes — no hidden live
builder.** The historical `RandomDataConfigBuilder` /
`CustomDataConfigBuilder` are REMOVED from the live source.
Three dead-reference sites remain (docstring + private script
+ build-cache mirror); the implementation slice cleans them
up in lockstep with the new builder.
"""
from __future__ import annotations

import importlib
import inspect
import pickle
import subprocess
from pathlib import Path

import pytest

import GenAIRR as ga
from GenAIRR.dataconfig.data_config import DataConfig
from GenAIRR.reference_models import (
    EmpiricalDistributionSpec,
    NpBaseModelSpec,
    ReferenceEmpiricalModels,
)
from GenAIRR.reference_rules import ReferenceRulesSpec


_REPO_ROOT = Path(__file__).resolve().parent.parent
_AUDIT_DOC = _REPO_ROOT / "docs" / "reference_cartridge_authoring_audit.md"


# ──────────────────────────────────────────────────────────────────
# Section A — pin_scaffold_* (live surfaces the builder will reuse)
# ──────────────────────────────────────────────────────────────────


def test_pin_scaffold_parse_fasta_is_public_and_parses_headers_correctly() -> None:
    """`parse_fasta(file)` is the existing public FASTA parser
    in `utilities/misc.py`, re-exported through
    `GenAIRR.utilities`. The builder's `from_fasta` constructor
    will feed it verbatim."""
    from GenAIRR.utilities import parse_fasta

    assert callable(parse_fasta)
    import io

    records = list(
        parse_fasta(
            io.StringIO(
                ">IGHV1-2*02\nGAGGTG\nCAGCTG\n>IGHV3-23*01\nGAGGTG\n"
            )
        )
    )
    assert records == [
        (">IGHV1-2*02", "GAGGTGCAGCTG"),
        (">IGHV3-23*01", "GAGGTG"),
    ]


def test_pin_scaffold_compute_v_region_boundaries_derives_subregions_from_gapped_seq() -> None:
    """`compute_v_region_boundaries(v_allele)` maps IMGT-
    gapped positions → ungapped (start, end) intervals. The
    builder's `infer_v_subregions()` stage feeds this for
    every V allele with `gapped_seq`."""
    from GenAIRR.utilities.imgt_regions import compute_v_region_boundaries

    # Pick a real bundled V allele with gapped_seq.
    cfg = ga.HUMAN_IGH_OGRDB
    for gene_alleles in cfg.v_alleles.values():
        for allele in gene_alleles:
            if getattr(allele, "gapped_seq", None):
                bounds = compute_v_region_boundaries(allele)
                assert set(bounds) == {"FWR1", "CDR1", "FWR2", "CDR2", "FWR3"}
                for label, (s, e) in bounds.items():
                    assert isinstance(s, int) and isinstance(e, int)
                    assert 0 <= s <= e
                return
    pytest.skip("no bundled V allele with gapped_seq found — refdata drift")


def test_pin_scaffold_cartridge_manifest_is_json_clean() -> None:
    """`DataConfig.cartridge_manifest()` returns a JSON-clean
    dict. The builder's build report includes a
    `manifest_snapshot` block populated from this method."""
    import json

    m = ga.HUMAN_IGH_OGRDB.cartridge_manifest()
    # JSON-clean round-trip works.
    s = json.dumps(m)
    again = json.loads(s)
    assert again["models"]["np_base_models"]["supported_kinds"] == [
        "uniform",
        "empirical_first_base",
        "markov",
    ]


def test_pin_scaffold_verify_integrity_and_compute_checksum_gate_corrupted_pickles() -> None:
    """`DataConfig.verify_integrity` and `compute_checksum`
    are the existing integrity gate. The builder's `.build()`
    calls `verify_integrity()` so a malformed `DataConfig`
    surfaces at build time."""
    cfg = ga.HUMAN_IGH_OGRDB
    cfg.verify_integrity()  # no raise on a bundled cartridge
    assert isinstance(cfg.compute_checksum(), str)
    assert len(cfg.compute_checksum()) == 64  # sha256 hex


def test_pin_scaffold_refdata_config_vj_vdj_constructors_and_add_allele_methods() -> None:
    """`RefDataConfig.{vj,vdj}()` + `add_v_allele` /
    `add_d_allele` / `add_j_allele` are the existing manual
    Rust-backed construction path. The new builder does NOT
    replace this surface; tests and perf workloads continue to
    use it directly."""
    cfg_vj = ga.RefDataConfig.vj()
    cfg_vj.add_v_allele("v1*01", "v1", b"GAGGTG", anchor=0)
    cfg_vj.add_j_allele("j1*01", "j1", b"TGGGGC", anchor=0)
    cfg_vdj = ga.RefDataConfig.vdj()
    cfg_vdj.add_v_allele("v1*01", "v1", b"GAGGTG", anchor=0)
    cfg_vdj.add_d_allele("d1*01", "d1", b"AAACCC")
    cfg_vdj.add_j_allele("j1*01", "j1", b"TGGGGC", anchor=0)


def test_pin_scaffold_reference_empirical_models_accepts_every_typed_plane() -> None:
    """`ReferenceEmpiricalModels` is the typed authoring
    entry point for NP / trim / P-nucleotide / NP-base
    distributions. The new builder constructs this object."""
    rm = ReferenceEmpiricalModels(
        np_lengths={"NP1": EmpiricalDistributionSpec([(0, 1.0)])},
        trims={"V_3": EmpiricalDistributionSpec([(0, 1.0)])},
        np_bases={"NP1": NpBaseModelSpec(kind="uniform")},
        p_nucleotide_lengths={"V_3": EmpiricalDistributionSpec([(0, 1.0)])},
    )
    rm.validate(chain_type="vdj")


def test_pin_scaffold_reference_rules_spec_is_the_rules_authoring_entry_point() -> None:
    """`ReferenceRulesSpec` is the existing programmable
    interpretation-layer authoring surface. The builder's
    `infer_rules()` stage emits this object. Pinned at the
    dataclass + validate-method level so the new builder's
    rules output flows through the existing validation."""
    spec = ReferenceRulesSpec()
    spec.validate()
    # The dataclass exposes the documented field surface.
    sig = inspect.signature(ReferenceRulesSpec)
    for required in ("allowed_bases", "v_anchor", "j_anchor"):
        assert required in sig.parameters, (
            f"ReferenceRulesSpec missing {required!r} kwarg — the audit's "
            f"rules-authoring surface drifted"
        )


def test_pin_scaffold_dataconfig_to_refdata_bridge_validates_malformed_cartridges() -> None:
    """`dataconfig_to_refdata` is the existing compile-time
    bridge that validates the cartridge shape. The new
    builder's output is a plain `DataConfig` that flows
    through this bridge unchanged."""
    from GenAIRR import dataconfig_to_refdata

    rd = dataconfig_to_refdata(ga.HUMAN_IGH_OGRDB)
    assert rd is not None


def test_pin_scaffold_bundled_cartridge_lazy_loader_lists_106_configs() -> None:
    """The lazy `__getattr__` cartridge loader in
    `GenAIRR.data` ships 106 bundled configs. The builder's
    output joins them as path 5 (per audit §4)."""
    names = ga.list_configs()
    # Allow for additions; pin a sane lower bound.
    assert len(names) >= 100, f"expected at least 100 bundled configs, got {len(names)}"
    assert "HUMAN_IGH_OGRDB" in names


# ──────────────────────────────────────────────────────────────────
# Section B — pin_present_* (dead-reference inventory)
# ──────────────────────────────────────────────────────────────────


def test_pin_present_data_config_build_report_field_is_none_on_every_bundled_cartridge() -> None:
    """`DataConfig.build_report` field exists (added at the
    time the historical builder shipped) but is
    unconditionally `None` on every bundled cartridge — there
    is no live producer. The new builder fills this in.
    Bundled configs are accessed via `GenAIRR.data.<NAME>`
    (top-level `getattr` only covers the 5 most-common
    cartridges)."""
    from GenAIRR import data as ga_data

    for name in ga.list_configs()[:10]:  # sample
        cfg = getattr(ga_data, name)
        assert getattr(cfg, "build_report", "<missing>") is None, (
            f"{name}.build_report is now non-None — a builder may have "
            f"started populating it. Update the audit doc + flip this "
            f"pin to a scaffold pin asserting the producer's output shape."
        )


def test_pin_present_build_report_docstring_now_references_new_builder() -> None:
    """Post-slice — the `DataConfig.build_report` field's
    docstring now names
    `GenAIRR.cartridge_builder.ReferenceCartridgeBuilder.build`
    instead of the dead `RandomDataConfigBuilder`. Flipped
    from the prior dead-reference present-pin."""
    src = (
        _REPO_ROOT / "src" / "GenAIRR" / "dataconfig" / "data_config.py"
    ).read_text(encoding="utf-8")
    assert "RandomDataConfigBuilder" not in src, (
        "RandomDataConfigBuilder reference still present in "
        "data_config.py — dead-reference cleanup regressed"
    )
    assert "ReferenceCartridgeBuilder" in src, (
        "data_config.py docstring no longer references the new "
        "ReferenceCartridgeBuilder — verify the build_report docstring"
    )
    # Dead module still does not exist (cleanup didn't accidentally
    # resurrect it).
    with pytest.raises(ModuleNotFoundError):
        importlib.import_module("GenAIRR.dataconfig.make.random")


def test_pin_present_private_build_script_now_raises_explicit_legacy_error() -> None:
    """Post-slice — the `.private/scripts/build_imgt_configs.py`
    script raises an explicit `NotImplementedError` at module-
    load time pointing at the new builder, rather than failing
    with a deep `ModuleNotFoundError` deep inside an obsolete
    import. Flipped from the prior dead-import present-pin."""
    script = _REPO_ROOT / ".private" / "scripts" / "build_imgt_configs.py"
    if not script.exists():
        pytest.skip(".private/scripts/build_imgt_configs.py absent")
    src = script.read_text(encoding="utf-8")
    # The dead import is no longer the load-time failure point.
    # The explicit raise comes first.
    assert "raise NotImplementedError(" in src, (
        "private build script no longer raises explicit "
        "NotImplementedError at module load — verify the dead-"
        "reference cleanup landed and the script's failure mode "
        "is still self-documenting"
    )
    assert "ReferenceCartridgeBuilder" in src, (
        "private build script does not reference the new builder — "
        "the porting hint regressed"
    )
    raise_pos = src.find("raise NotImplementedError(")
    legacy_import_pos = src.find("from GenAIRR.dataconfig.make.random import")
    assert legacy_import_pos == -1 or raise_pos < legacy_import_pos, (
        "the dead import would fire BEFORE the explicit raise — "
        "the legacy guard needs to come first or the dead import "
        "needs to be moved into the unreachable body"
    )


# ──────────────────────────────────────────────────────────────────
# Section C — pin_scaffold_* (historical builder gone from live)
# ──────────────────────────────────────────────────────────────────


def test_pin_scaffold_historical_random_builder_module_is_gone() -> None:
    """`GenAIRR.dataconfig.make.random` does not exist in
    the live source tree. Its historical implementation
    survives only in the build-cache mirror at
    `docs/build/lib.../GenAIRR/dataconfig/make/`, which is
    NOT on the import path. Pinned so a future regression
    that resurrects the dead namespace surfaces here for
    explicit review."""
    with pytest.raises(ModuleNotFoundError):
        importlib.import_module("GenAIRR.dataconfig.make.random")


def test_pin_scaffold_historical_custom_builder_module_is_gone() -> None:
    """Companion to the random-builder pin — the
    `CustomDataConfigBuilder` companion is also gone from
    the live source tree."""
    with pytest.raises(ModuleNotFoundError):
        importlib.import_module("GenAIRR.dataconfig.make.custom")


def test_pin_scaffold_historical_builders_dataconfig_make_namespace_absent() -> None:
    """The whole `GenAIRR.dataconfig.make` namespace is
    absent. Pinned at the namespace level so a regression
    that re-introduces ANY of the auxiliary builders
    (`TrimmingProbabilityGenerator`, `NPMarkovParameterBuilder`,
    etc.) surfaces here."""
    with pytest.raises(ModuleNotFoundError):
        importlib.import_module("GenAIRR.dataconfig.make")


# ──────────────────────────────────────────────────────────────────
# Section D — pin_absence_* (gaps the implementation slice closes)
# ──────────────────────────────────────────────────────────────────


def test_pin_present_reference_cartridge_builder_module_landed_at_cartridge_builder() -> None:
    """Post-slice — `GenAIRR.cartridge_builder` is the new
    builder module. The two sibling location candidates
    (`GenAIRR.cartridge` / `GenAIRR.dataconfig.cartridge_builder`)
    stay absent so a future regression that duplicates the
    namespace surfaces here. Flipped from the prior absence
    pin."""
    mod = importlib.import_module("GenAIRR.cartridge_builder")
    assert hasattr(mod, "ReferenceCartridgeBuilder")
    assert hasattr(mod, "CartridgeBuildReport")
    for forbidden_path in (
        "GenAIRR.cartridge",
        "GenAIRR.dataconfig.cartridge_builder",
    ):
        with pytest.raises(ModuleNotFoundError):
            importlib.import_module(forbidden_path)


def test_pin_present_reference_cartridge_builder_class_at_top_level() -> None:
    """Post-slice — `GenAIRR.ReferenceCartridgeBuilder` is
    the only public builder name at the top level. The
    historical / alternative names stay absent so a future
    regression that introduces parallel surfaces surfaces
    here."""
    assert hasattr(ga, "ReferenceCartridgeBuilder")
    assert hasattr(ga, "CartridgeBuildReport")
    # Historical / alternative names remain absent.
    for forbidden_name in (
        "CartridgeBuilder",
        "DataConfigBuilder",
        "RandomDataConfigBuilder",
        "CustomDataConfigBuilder",
    ):
        assert not hasattr(ga, forbidden_name), (
            f"GenAIRR.{forbidden_name} is now reachable — a parallel "
            f"builder surface may have been introduced. Update the "
            f"audit doc + pins or remove the duplicate."
        )


def test_pin_present_from_fasta_constructor_on_builder_only() -> None:
    """Post-slice — `ReferenceCartridgeBuilder.from_fasta(...)`
    is the new FASTA constructor; `DataConfig` and
    `RefDataConfig` still do NOT carry `from_fasta` (the
    builder is the canonical entry point, not the cartridge
    types themselves)."""
    assert hasattr(ga.ReferenceCartridgeBuilder, "from_fasta")
    assert callable(ga.ReferenceCartridgeBuilder.from_fasta)
    for owner_name in ("DataConfig", "RefDataConfig"):
        owner = getattr(ga, owner_name, None)
        if owner is None:
            continue
        assert not hasattr(owner, "from_fasta"), (
            f"{owner_name}.from_fasta now exists — the audit recommended "
            f"the entry point live on the builder, not the cartridge "
            f"types. Verify."
        )


def test_pin_absence_no_from_airr_constructor_in_v1() -> None:
    """v1 boundary — no `from_airr` / `from_rearrangement_tsv`
    constructor on EITHER the cartridge types or the builder.
    The audit §11 defers AIRR-rearrangement-data inference
    to a follow-up slice. Pinned at all three locations so
    a future slice introducing it surfaces here for
    coordinated review."""
    for owner_name in ("DataConfig", "RefDataConfig", "ReferenceCartridgeBuilder"):
        owner = getattr(ga, owner_name, None)
        if owner is None:
            continue
        for method_name in ("from_airr", "from_rearrangement_tsv"):
            assert not hasattr(owner, method_name), (
                f"{owner_name}.{method_name} now exists — v1 audit boundary "
                f"deferred AIRR-rearrangement inference. Verify."
            )


def test_pin_present_infer_step_methods_on_builder_not_cartridge() -> None:
    """Post-slice — `infer_*` step methods live on
    `ReferenceCartridgeBuilder` (per audit §8 design), NOT on
    `DataConfig` / `RefDataConfig`. v1 ships
    `infer_identity` + `infer_v_subregions`; the audit's
    `infer_anchors_from_imgt_gapped` is deferred to a
    follow-up slice."""
    builder_cls = ga.ReferenceCartridgeBuilder
    assert hasattr(builder_cls, "infer_identity")
    assert hasattr(builder_cls, "infer_v_subregions")
    # The cartridge types stay clean.
    for owner_name in ("DataConfig", "RefDataConfig"):
        owner = getattr(ga, owner_name, None)
        if owner is None:
            continue
        infer_methods = [m for m in dir(owner) if m.startswith("infer_")]
        assert not infer_methods, (
            f"{owner_name} now exposes infer_* methods: {infer_methods}. "
            f"Verify the audit doc's builder-owns-inference boundary "
            f"held."
        )


def test_pin_estimate_step_method_boundary_post_p_nucleotide_length_slice() -> None:
    """Estimator-method boundary, post P-nucleotide-length slice.

    `estimate_allele_usage`, `estimate_trim_distributions`,
    `estimate_np_length_distributions`,
    `estimate_np_base_model`, and
    `estimate_p_nucleotide_lengths` are the five
    statistical estimators in the v1 surface (see the
    respective design docs and contract files). Only
    `estimate_shm_rates` remains deferred from the
    original audit's deferred list. Pinned at all three
    locations so a future slice introducing more
    estimators surfaces here for explicit review."""
    allowed_estimators = {
        "estimate_allele_usage",
        "estimate_trim_distributions",
        "estimate_np_length_distributions",
        "estimate_np_base_model",
        "estimate_p_nucleotide_lengths",
    }
    for owner_name in ("DataConfig", "RefDataConfig", "ReferenceCartridgeBuilder"):
        owner = getattr(ga, owner_name, None)
        if owner is None:
            continue
        estimate_methods = {m for m in dir(owner) if m.startswith("estimate_")}
        # On DataConfig / RefDataConfig: nothing should land.
        if owner_name in ("DataConfig", "RefDataConfig"):
            assert not estimate_methods, (
                f"{owner_name} now exposes estimate_* methods: "
                f"{sorted(estimate_methods)}. Estimators belong on "
                f"the builder, not the config surface — verify."
            )
            continue
        # On the builder: only the two shipped estimators are in scope.
        unexpected = estimate_methods - allowed_estimators
        assert not unexpected, (
            f"{owner_name} exposes unexpected estimate_* methods: "
            f"{sorted(unexpected)}. v1 audit boundary deferred these — "
            f"verify the new slice updated the audit doc + pins in "
            f"lockstep."
        )


def test_pin_present_cartridge_build_report_dataclass_landed() -> None:
    """Post-slice — `CartridgeBuildReport` dataclass is the
    typed container for the audit trail. Re-exported at the
    top level for discoverability and importable from the
    new builder module."""
    assert hasattr(ga, "CartridgeBuildReport")
    from GenAIRR.cartridge_builder import CartridgeBuildReport

    assert CartridgeBuildReport is ga.CartridgeBuildReport
    # Sibling locations remain absent.
    for forbidden_path in (
        "GenAIRR.cartridge",
        "GenAIRR.dataconfig.cartridge_builder",
        "GenAIRR.build_report",
    ):
        with pytest.raises(ModuleNotFoundError):
            importlib.import_module(forbidden_path)


# ──────────────────────────────────────────────────────────────────
# Section E — pin_scaffold_* (integration boundary)
# ──────────────────────────────────────────────────────────────────


def test_pin_scaffold_pickle_round_trip_preserves_build_report() -> None:
    """The new builder writes a typed `CartridgeBuildReport`
    onto `cfg.build_report`. The audit requires that pickle
    round-trip preserves the field (so a `pickle.load`
    of a builder-produced cartridge carries the report)."""
    # Pre-slice: build_report is None and round-trips as None.
    # Pinned so a future slice that adds a producer must also
    # confirm round-trip with the actual report dataclass.
    cfg = ga.HUMAN_IGH_OGRDB
    blob = pickle.dumps(cfg)
    cfg2 = pickle.loads(blob)
    assert cfg.build_report == cfg2.build_report  # both None today


def test_pin_scaffold_manual_dataconfig_construction_remains_supported() -> None:
    """The new builder is the FIFTH cartridge-creation path
    (per audit §4). The existing fourth path (manual
    `DataConfig(...)` construction) MUST keep working — tests
    and lightweight programmatic creation depend on it."""
    cfg = DataConfig(name="audit-fixture")
    assert cfg.name == "audit-fixture"
    assert cfg.build_report is None
    # The empty cartridge survives `cartridge_manifest()`.
    m = cfg.cartridge_manifest()
    assert isinstance(m, dict)
    assert m["identity"]["name"] == "audit-fixture"


# ──────────────────────────────────────────────────────────────────
# Section F — Doc anchor
# ──────────────────────────────────────────────────────────────────


def test_pin_scaffold_audit_doc_exists_and_references_contract() -> None:
    assert _AUDIT_DOC.exists(), "reference_cartridge_authoring_audit.md missing"
    doc = _AUDIT_DOC.read_text(encoding="utf-8")
    assert "test_reference_cartridge_authoring_contract.py" in doc, (
        "audit doc no longer references the contract file"
    )
    for marker in (
        "## 1. Current creation paths",
        "## 3. Dead / historical references",
        "## 8. Q5",
        "## 12. Test surface",
        "## 14. Summary table",
    ):
        assert marker in doc, f"audit doc missing section marker {marker!r}"
