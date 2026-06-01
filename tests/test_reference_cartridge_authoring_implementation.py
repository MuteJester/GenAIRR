"""End-to-end implementation tests for the
**ReferenceCartridgeBuilder v1** slice.

Companion to
[`docs/reference_cartridge_authoring_audit.md`](../docs/reference_cartridge_authoring_audit.md)
and the contract pins in
[`tests/test_reference_cartridge_authoring_contract.py`](test_reference_cartridge_authoring_contract.py).

Covers the 12 behaviours the implementation slice promised:

1. Top-level imports.
2. ``from_fasta(...).build()`` round-trip.
3. Built config compiles through ``Experiment.on(cfg)``.
4. ``infer_identity()`` populates manifest identity.
5. ``infer_v_subregions()`` derives coverage from gapped V FASTA.
6. Missing gapped V → report warning, not crash.
7. ``with_rules()`` / ``with_models()`` attach the specs.
8. ``build_report`` survives pickle round-trip.
9. Report ``.to_dict()`` is JSON-clean.
10. Manual ``DataConfig(...)`` construction remains unchanged.
11. Dead-reference pins flipped (covered by the contract file;
    re-asserted here as an integration probe).
12. Statistical-estimator absence pins remain (contract file
    is the source of truth; integration probe here).
"""
from __future__ import annotations

import io
import json
import pickle

import pytest

import GenAIRR as ga
from GenAIRR.reference_models import (
    EmpiricalDistributionSpec,
    NpBaseModelSpec,
    ReferenceEmpiricalModels,
)
from GenAIRR.reference_rules import ReferenceRulesSpec


# ──────────────────────────────────────────────────────────────────
# Helpers — synthesise FASTAs from a bundled cartridge so the
# builder has realistic IMGT-gapped V sequences to chew on.
# ──────────────────────────────────────────────────────────────────


def _fastas_from_bundled(
    src: "ga.DataConfig",
    *,
    n_v: int = 3,
    n_j: int = 2,
    n_d: int = 0,
    gapped_v: bool = True,
) -> "tuple[str, str, str | None]":
    """Build raw FASTA strings from the first ``n_v`` / ``n_j``
    genes' first allele in a bundled cartridge. The V FASTA
    uses gapped or ungapped sequences per ``gapped_v``."""
    v_lines: list[str] = []
    for gene_alleles in list(src.v_alleles.values())[:n_v]:
        a = gene_alleles[0]
        seq = a.gapped_seq if gapped_v and a.gapped_seq else a.ungapped_seq
        v_lines.append(f">{a.name}\n{seq}\n")
    j_lines: list[str] = []
    for gene_alleles in list(src.j_alleles.values())[:n_j]:
        a = gene_alleles[0]
        j_lines.append(f">{a.name}\n{a.ungapped_seq}\n")
    d_fasta: str | None = None
    if n_d > 0 and src.d_alleles:
        d_lines: list[str] = []
        for gene_alleles in list(src.d_alleles.values())[:n_d]:
            a = gene_alleles[0]
            d_lines.append(f">{a.name}\n{a.ungapped_seq}\n")
        d_fasta = "".join(d_lines)
    return "".join(v_lines), "".join(j_lines), d_fasta


# ──────────────────────────────────────────────────────────────────
# 1. Top-level imports work
# ──────────────────────────────────────────────────────────────────


def test_top_level_imports_resolve() -> None:
    """``from GenAIRR import ReferenceCartridgeBuilder,
    CartridgeBuildReport`` works and resolves to the same
    objects as the dotted-module path."""
    from GenAIRR import CartridgeBuildReport, ReferenceCartridgeBuilder

    assert ReferenceCartridgeBuilder is ga.ReferenceCartridgeBuilder
    assert CartridgeBuildReport is ga.CartridgeBuildReport
    # And module-level access.
    import GenAIRR.cartridge_builder as cb

    assert cb.ReferenceCartridgeBuilder is ReferenceCartridgeBuilder
    assert cb.CartridgeBuildReport is CartridgeBuildReport


# ──────────────────────────────────────────────────────────────────
# 2. from_fasta(...).build() round-trip
# ──────────────────────────────────────────────────────────────────


def test_from_fasta_build_returns_usable_dataconfig() -> None:
    """The minimal happy path — supply V + J FASTA, call
    ``build()``, get back a valid ``DataConfig`` with the
    expected counts. Cartridge survives
    ``verify_integrity()`` (called inside ``build()``)."""
    v_fasta, j_fasta, _ = _fastas_from_bundled(ga.HUMAN_IGK_OGRDB)
    cfg = (
        ga.ReferenceCartridgeBuilder.from_fasta(
            v_fasta=v_fasta,
            j_fasta=j_fasta,
            chain_type="BCR_LIGHT_KAPPA",
        )
        .build()
    )
    assert isinstance(cfg, ga.DataConfig)
    assert len(cfg.v_alleles) == 3
    assert len(cfg.j_alleles) == 2
    assert cfg.build_report is not None
    assert cfg.schema_sha256
    # And the round-trip checksum verifies.
    cfg.verify_integrity()


# ──────────────────────────────────────────────────────────────────
# 3. Built config compiles through Experiment
# ──────────────────────────────────────────────────────────────────


def test_built_config_compiles_through_experiment() -> None:
    """The cartridge produced by ``build()`` is a drop-in for
    ``Experiment.on(cfg)``. Synthetic alleles built without
    the native C anchor resolver lack anchors, so the
    compile step needs ``allow_curatable_refdata()`` — that's
    the expected path for new cartridges per the audit (§7,
    "anchor inference when gapped_seq doesn't follow IMGT
    convention")."""
    v_fasta, j_fasta, _ = _fastas_from_bundled(ga.HUMAN_IGK_OGRDB)
    cfg = (
        ga.ReferenceCartridgeBuilder.from_fasta(
            v_fasta=v_fasta,
            j_fasta=j_fasta,
            chain_type="BCR_LIGHT_KAPPA",
        )
        .infer_identity(species="HUMAN", reference_set="ITER1")
        .build()
    )
    compiled = (
        ga.Experiment.on(cfg)
        .allow_curatable_refdata()
        .recombine()
        .compile()
    )
    assert compiled is not None


# ──────────────────────────────────────────────────────────────────
# 4. infer_identity() populates manifest identity
# ──────────────────────────────────────────────────────────────────


def test_infer_identity_populates_manifest_identity() -> None:
    """``infer_identity()`` writes the species / locus /
    reference_set / name / source onto the cartridge's
    metadata + allele provenance. Manifest identity reflects
    every authored field."""
    v_fasta, j_fasta, _ = _fastas_from_bundled(ga.HUMAN_IGK_OGRDB)
    cfg = (
        ga.ReferenceCartridgeBuilder.from_fasta(
            v_fasta=v_fasta,
            j_fasta=j_fasta,
            chain_type="BCR_LIGHT_KAPPA",
        )
        .infer_identity(
            species="HUMAN",
            locus="IGK",
            reference_set="OGRDB-derived-2026Q1",
            name="MY_HUMAN_IGK",
            source="ReferenceCartridgeBuilder",
        )
        .build()
    )
    manifest = cfg.cartridge_manifest()
    identity = manifest["identity"]
    assert identity["name"] == "MY_HUMAN_IGK"
    assert identity["species"] == "HUMAN"
    assert identity["reference_set"] == "OGRDB-derived-2026Q1"


# ──────────────────────────────────────────────────────────────────
# 5. infer_v_subregions() derives coverage from gapped V FASTA
# ──────────────────────────────────────────────────────────────────


def test_infer_v_subregions_derives_coverage_from_gapped_v_fasta() -> None:
    """V alleles supplied with IMGT-gapped sequences get
    their five canonical subregion intervals (FWR1 / CDR1 /
    FWR2 / CDR2 / FWR3) computed and attached. Report records
    the annotated count."""
    v_fasta, j_fasta, _ = _fastas_from_bundled(ga.HUMAN_IGK_OGRDB, gapped_v=True)
    builder = ga.ReferenceCartridgeBuilder.from_fasta(
        v_fasta=v_fasta,
        j_fasta=j_fasta,
        chain_type="BCR_LIGHT_KAPPA",
    )
    builder.infer_v_subregions()
    cfg = builder.build()
    annotated = 0
    for gene_alleles in cfg.v_alleles.values():
        for allele in gene_alleles:
            if allele.subregions:
                assert set(allele.subregions) == {
                    "FWR1",
                    "CDR1",
                    "FWR2",
                    "CDR2",
                    "FWR3",
                }
                annotated += 1
    assert annotated == 3
    # Report records the count.
    report = builder.report()
    subregion_stage = next(
        s for s in report.stages if s["stage"] == "infer_v_subregions"
    )
    assert subregion_stage["inferred"]["alleles_annotated"] == 3


# ──────────────────────────────────────────────────────────────────
# 6. Missing gapped V emits report warning, not crash
# ──────────────────────────────────────────────────────────────────


def test_missing_gapped_v_emits_warning_not_crash() -> None:
    """When a V allele lacks ``gapped_seq``,
    ``infer_v_subregions`` skips it AND records the skip in
    the per-stage report, rather than raising. The skip
    count is exposed for downstream auditing."""
    # Build a V FASTA from UNGAPPED sequences only (the
    # bundled allele's ungapped_seq contains no dots).
    src = ga.HUMAN_IGK_OGRDB
    v_lines: list[str] = []
    for gene_alleles in list(src.v_alleles.values())[:2]:
        a = gene_alleles[0]
        v_lines.append(f">{a.name}\n{a.ungapped_seq}\n")
    v_fasta = "".join(v_lines)
    _, j_fasta, _ = _fastas_from_bundled(src)

    builder = ga.ReferenceCartridgeBuilder.from_fasta(
        v_fasta=v_fasta,
        j_fasta=j_fasta,
        chain_type="BCR_LIGHT_KAPPA",
    )
    # Should not raise.
    builder.infer_v_subregions()
    report = builder.report()
    subregion_stage = next(
        s for s in report.stages if s["stage"] == "infer_v_subregions"
    )
    # No alleles annotated; skips reported per allele.
    assert subregion_stage["inferred"]["alleles_skipped_no_gapped"] == 2
    # And a stage-level warning surfaces the gap.
    assert any(
        "lack gapped_seq" in w for w in subregion_stage["warnings"]
    )


# ──────────────────────────────────────────────────────────────────
# 7. with_rules() / with_models() attach the specs
# ──────────────────────────────────────────────────────────────────


def test_with_rules_and_with_models_attach_specs() -> None:
    """User-authored ``ReferenceRulesSpec`` /
    ``ReferenceEmpiricalModels`` flow through to the built
    cartridge. Stage entries record the typed plane keys."""
    v_fasta, j_fasta, _ = _fastas_from_bundled(ga.HUMAN_IGK_OGRDB)
    rules = ReferenceRulesSpec()  # default — valid as-is
    models = ReferenceEmpiricalModels(
        np_lengths={"NP1": EmpiricalDistributionSpec([(0, 1.0), (1, 1.0)])},
        np_bases={
            "NP1": NpBaseModelSpec(
                kind="empirical_first_base",
                first_base={"A": 1.0, "C": 1.0, "G": 1.0, "T": 1.0},
            )
        },
        p_nucleotide_lengths={"V_3": EmpiricalDistributionSpec([(0, 1.0), (1, 1.0)])},
    )
    builder = (
        ga.ReferenceCartridgeBuilder.from_fasta(
            v_fasta=v_fasta,
            j_fasta=j_fasta,
            chain_type="BCR_LIGHT_KAPPA",
        )
        .with_rules(rules)
        .with_models(models)
    )
    cfg = builder.build()
    assert cfg.reference_rules is rules
    assert cfg.reference_models is models
    report = builder.report()
    rules_stage = next(s for s in report.stages if s["stage"] == "with_rules")
    models_stage = next(s for s in report.stages if s["stage"] == "with_models")
    assert rules_stage["inferred"]["rules_attached"] is True
    assert "NP1" in models_stage["inferred"]["np_length_keys"]
    assert "V_3" in models_stage["inferred"]["p_nucleotide_length_keys"]


# ──────────────────────────────────────────────────────────────────
# 8. build_report survives pickle round-trip
# ──────────────────────────────────────────────────────────────────


def test_build_report_survives_pickle_round_trip() -> None:
    """A cartridge produced by the builder, pickled and
    unpickled, retains its ``build_report`` with the same
    stages list (deepcopy via pickle preserves the dataclass
    by construction; this test pins the discipline)."""
    v_fasta, j_fasta, _ = _fastas_from_bundled(ga.HUMAN_IGK_OGRDB)
    cfg = (
        ga.ReferenceCartridgeBuilder.from_fasta(
            v_fasta=v_fasta,
            j_fasta=j_fasta,
            chain_type="BCR_LIGHT_KAPPA",
        )
        .infer_identity(species="HUMAN", reference_set="ITER1")
        .infer_v_subregions()
        .build()
    )
    blob = pickle.dumps(cfg)
    cfg2 = pickle.loads(blob)
    assert cfg2.build_report is not None
    # Stage list deep-equal.
    assert [s["stage"] for s in cfg.build_report.stages] == [
        s["stage"] for s in cfg2.build_report.stages
    ]
    # Manifest snapshot equality at the identity layer.
    a_identity = cfg.build_report.manifest_snapshot["identity"]
    b_identity = cfg2.build_report.manifest_snapshot["identity"]
    assert a_identity == b_identity


# ──────────────────────────────────────────────────────────────────
# 9. Report .to_dict() is JSON-clean
# ──────────────────────────────────────────────────────────────────


def test_report_to_dict_is_json_serialisable() -> None:
    """``CartridgeBuildReport.to_dict()`` produces a dict
    that round-trips through ``json.dumps`` / ``json.loads``
    without raising. Required for downstream tooling
    (CI artifacts, audit dashboards)."""
    v_fasta, j_fasta, _ = _fastas_from_bundled(ga.HUMAN_IGK_OGRDB)
    cfg = (
        ga.ReferenceCartridgeBuilder.from_fasta(
            v_fasta=v_fasta,
            j_fasta=j_fasta,
            chain_type="BCR_LIGHT_KAPPA",
        )
        .infer_identity(species="HUMAN", reference_set="ITER1")
        .infer_v_subregions()
        .build()
    )
    payload = cfg.build_report.to_dict()
    s = json.dumps(payload)
    again = json.loads(s)
    assert isinstance(again, dict)
    assert again["stages"][0]["stage"] == "from_fasta"
    # The manifest snapshot is dictionarised by `cartridge_manifest`
    # itself — the report doesn't have to coerce it.
    assert "manifest_snapshot" in again


# ──────────────────────────────────────────────────────────────────
# 10. Manual DataConfig construction remains unchanged
# ──────────────────────────────────────────────────────────────────


def test_manual_dataconfig_construction_remains_supported() -> None:
    """The new builder is the FIFTH cartridge-creation path
    per audit §4. The existing fourth path (manual
    ``DataConfig(name=...)``) MUST keep working — tests and
    lightweight programmatic creation depend on it."""
    cfg = ga.DataConfig(name="iter1-fixture")
    assert cfg.name == "iter1-fixture"
    # Bare manual cartridge has no build_report.
    assert cfg.build_report is None
    # But manifest still works.
    m = cfg.cartridge_manifest()
    assert isinstance(m, dict)
    assert m["identity"]["name"] == "iter1-fixture"


# ──────────────────────────────────────────────────────────────────
# 11. Dead-reference cleanup integration probe
# ──────────────────────────────────────────────────────────────────


def test_dead_reference_cleanup_data_config_docstring_and_private_script() -> None:
    """Integration check that the lockstep cleanup landed:

    - ``DataConfig.build_report`` docstring names the new
      builder.
    - ``.private/scripts/build_imgt_configs.py`` raises
      ``NotImplementedError`` at module-load time rather
      than ``ModuleNotFoundError`` on the dead import.

    The contract file's `pin_present_*` pins are the
    authoritative source; this test is the integration
    probe so a regression that breaks one but not the
    other still surfaces here."""
    from pathlib import Path

    repo = Path(__file__).resolve().parent.parent
    dc_src = (
        repo / "src" / "GenAIRR" / "dataconfig" / "data_config.py"
    ).read_text(encoding="utf-8")
    assert "RandomDataConfigBuilder" not in dc_src
    assert "ReferenceCartridgeBuilder" in dc_src

    script = repo / ".private" / "scripts" / "build_imgt_configs.py"
    if script.exists():
        src = script.read_text(encoding="utf-8")
        assert "raise NotImplementedError(" in src
        assert "ReferenceCartridgeBuilder" in src


# ──────────────────────────────────────────────────────────────────
# 12. Statistical-estimator boundary post trim-distribution slice
# ──────────────────────────────────────────────────────────────────


def test_v1_defers_remaining_statistical_estimators() -> None:
    """Estimator-method boundary, post P-nucleotide-length slice.

    ``estimate_allele_usage``, ``estimate_trim_distributions``,
    ``estimate_np_length_distributions``,
    ``estimate_np_base_model``, and
    ``estimate_p_nucleotide_lengths`` have landed (see the
    respective design docs and implementation test files).
    Only ``estimate_shm_rates`` remains unimplemented from
    the original audit's deferred list. Calling it must
    produce a clean ``AttributeError``, not a partial
    half-implementation."""
    v_fasta, j_fasta, _ = _fastas_from_bundled(ga.HUMAN_IGK_OGRDB)
    builder = ga.ReferenceCartridgeBuilder.from_fasta(
        v_fasta=v_fasta,
        j_fasta=j_fasta,
        chain_type="BCR_LIGHT_KAPPA",
    )
    # The five estimators that DID land.
    assert hasattr(builder, "estimate_allele_usage")
    assert hasattr(builder, "estimate_trim_distributions")
    assert hasattr(builder, "estimate_np_length_distributions")
    assert hasattr(builder, "estimate_np_base_model")
    assert hasattr(builder, "estimate_p_nucleotide_lengths")
    # Only estimate_shm_rates stays deferred.
    assert not hasattr(builder, "estimate_shm_rates"), (
        "builder exposes estimate_shm_rates — v1 audit boundary "
        "deferred this. Verify."
    )
    with pytest.raises(AttributeError):
        getattr(builder, "estimate_shm_rates")


# ──────────────────────────────────────────────────────────────────
# Bonus — error paths that the audit's §10 specified
# ──────────────────────────────────────────────────────────────────


def test_from_fasta_rejects_d_fasta_on_vj_chain() -> None:
    """User-error catch: supplying ``d_fasta`` to a VJ chain
    surfaces a clear ``ValueError`` rather than silently
    producing a malformed cartridge."""
    v_fasta, j_fasta, _ = _fastas_from_bundled(ga.HUMAN_IGK_OGRDB)
    with pytest.raises(ValueError, match="d_fasta supplied"):
        ga.ReferenceCartridgeBuilder.from_fasta(
            v_fasta=v_fasta,
            j_fasta=j_fasta,
            d_fasta=">d1*01\nGGGTAT\n",
            chain_type="BCR_LIGHT_KAPPA",
        )


def test_from_fasta_requires_d_fasta_on_vdj_chain() -> None:
    """Symmetric — a VDJ chain without ``d_fasta`` is a
    cartridge-incomplete authoring error."""
    v_fasta, j_fasta, _ = _fastas_from_bundled(ga.HUMAN_IGH_OGRDB)
    with pytest.raises(ValueError, match="d_fasta is required"):
        ga.ReferenceCartridgeBuilder.from_fasta(
            v_fasta=v_fasta,
            j_fasta=j_fasta,
            chain_type="BCR_HEAVY",
        )


def test_duplicate_allele_names_are_rejected_with_report_entry() -> None:
    """Duplicate allele names in a FASTA produce a structured
    rejection entry on the build report, not a silent
    overwrite or crash."""
    src = ga.HUMAN_IGK_OGRDB
    first_v = list(src.v_alleles.values())[0][0]
    v_fasta = (
        f">{first_v.name}\n{first_v.ungapped_seq}\n"
        f">{first_v.name}\n{first_v.ungapped_seq}\n"
    )
    _, j_fasta, _ = _fastas_from_bundled(src)
    builder = ga.ReferenceCartridgeBuilder.from_fasta(
        v_fasta=v_fasta,
        j_fasta=j_fasta,
        chain_type="BCR_LIGHT_KAPPA",
    )
    report = builder.report()
    dupes = [
        r for r in report.rejected
        if r.get("reason") == "duplicate_name" and r.get("segment") == "V"
    ]
    assert len(dupes) == 1, report.rejected
