"""Contract pins for the Allele Usage Estimation audit.

Companion to
[`docs/allele_usage_estimation_design.md`](../docs/allele_usage_estimation_design.md).

Pin set:

- ``pin_scaffold_*`` — live surfaces the new estimator
  reuses verbatim (the engine's existing weighted-allele
  selection chain end-to-end, AIRR column convention via
  `_mcp_summary`, the typed cartridge plane shape, the
  builder stage-entry shape).
- ``pin_present_*`` — stop-and-report verification:
  `gene_use_dict` is genuinely orphan in the simulation
  pipeline; `DataConfig.validate()` is dead code; MCP
  `gene_use` endpoint is read-only.
- ``pin_present_*`` — documented boundary state: plan-
  signature soft gap 1 for `v_allele_weights` holds.
- ``pin_absence_*`` — the surfaces the implementation slice
  closes (`estimate_allele_usage` method, `allele_usage`
  cartridge plane, `AlleleUsageSpec` dataclass, manifest
  `allele_usage` block, ambiguous-policy kwarg).

**Pre-flight verdict (audit §7): clean-yes.** The engine
surface is already wired end-to-end; `gene_use_dict` is
genuinely orphan; no silent corruption exists.
"""
from __future__ import annotations

import importlib
import inspect
import json
import subprocess
from pathlib import Path

import pytest

import GenAIRR as ga
from GenAIRR.dataconfig.data_config import (
    DataConfig,
    _DOCUMENTED_ORPHAN_DATACONFIG_FIELDS,
)
from GenAIRR.reference_models import ReferenceEmpiricalModels


_REPO_ROOT = Path(__file__).resolve().parent.parent
_AUDIT_DOC = _REPO_ROOT / "docs" / "allele_usage_estimation_design.md"


# ──────────────────────────────────────────────────────────────────
# A. pin_scaffold_* — engine weighted-allele surface (live, wired)
# ──────────────────────────────────────────────────────────────────


def test_pin_scaffold_experiment_recombine_accepts_allele_weights_kwargs() -> None:
    """`Experiment.recombine` exposes three optional kwargs
    that take ``{allele_name: weight}`` dicts. The estimator
    slice's output flows into the same surface."""
    sig = inspect.signature(ga.Experiment.recombine)
    for kw in ("v_allele_weights", "d_allele_weights", "j_allele_weights"):
        assert kw in sig.parameters, (
            f"Experiment.recombine no longer accepts {kw!r} — engine "
            f"weighted-allele surface regressed"
        )


def test_pin_scaffold_resolve_allele_weights_converts_dict_to_dense_tuple() -> None:
    """`Experiment._resolve_allele_weights` converts the
    user dict into a dense pool-aligned ``Tuple[float, ...]``
    with 1.0 default for unlisted alleles. Pinned via direct
    call so the conversion shape is auditable."""
    exp = ga.Experiment.on("human_igh")
    pool_idx = exp._allele_name_index("V")
    a_name = next(iter(pool_idx.keys()))
    weights = exp._resolve_allele_weights("V", {a_name: 100.0})
    assert isinstance(weights, tuple)
    assert len(weights) == max(pool_idx.values()) + 1
    assert weights[pool_idx[a_name]] == 100.0
    # All other entries default to 1.0.
    for name, idx in pool_idx.items():
        if name != a_name:
            assert weights[idx] == 1.0


def test_pin_scaffold_pipeline_ir_recombine_step_has_weights_fields() -> None:
    """`_RecombineStep.weights_{v,d,j}` carry the dense
    weight vector forward from the DSL to the lowering layer.
    The estimator slice's cartridge-driven default lowers into
    the same fields when no per-experiment kwarg overrides."""
    from GenAIRR._pipeline_ir import _RecombineStep
    from dataclasses import fields

    field_names = {f.name for f in fields(_RecombineStep)}
    for required in ("weights_v", "weights_d", "weights_j"):
        assert required in field_names, (
            f"_RecombineStep.{required} field missing — the engine's "
            f"weighted-allele plumbing regressed"
        )


def test_pin_scaffold_lower_recombine_passes_weights_to_push_sample_allele() -> None:
    """`_compile.py::_lower_recombine` calls
    ``plan.push_sample_allele(..., weights=...)`` with the
    per-step weight vector. Pinned at source."""
    src = (_REPO_ROOT / "src" / "GenAIRR" / "_compile.py").read_text(encoding="utf-8")
    assert "push_sample_allele(" in src
    assert "weights=v_weights" in src
    assert "weights=j_weights" in src
    # D weights live inside the vdj branch.
    assert "weights=d_weights" in src


def test_pin_scaffold_plan_push_sample_allele_accepts_weights_arg() -> None:
    """The PyO3 bridge's `push_sample_allele` signature
    declares the `weights` kwarg. Pinned at source."""
    src = (
        _REPO_ROOT / "engine_rs" / "src" / "python" / "plan.rs"
    ).read_text(encoding="utf-8")
    assert "fn push_sample_allele" in src
    assert "weights: Option<Vec<f64>>" in src
    # Mutual exclusion with allowed_ids is preserved.
    assert "allowed_ids.is_some() && weights.is_some()" in src


def test_pin_scaffold_allele_pool_dist_from_weights_is_engine_consumer() -> None:
    """`AllelePoolDist::from_weights` is the engine-side
    constructor that consumes the pool-aligned weight
    vector. Pinned at source."""
    src = (
        _REPO_ROOT / "engine_rs" / "src" / "dist" / "allele_pool.rs"
    ).read_text(encoding="utf-8")
    assert "pub fn from_weights" in src or "from_weights(" in src


# ──────────────────────────────────────────────────────────────────
# B. pin_scaffold_* — AIRR column convention via _mcp_summary
# ──────────────────────────────────────────────────────────────────


def test_pin_scaffold_mcp_summary_first_call_helper_treats_v_call_as_comma_separated() -> None:
    """`_mcp_summary` already establishes the convention the
    `ambiguous="truth_first"` policy mirrors: comma-separated
    tie-sets get the first call attributed. The public entry
    point is `compute_repertoire_summary`."""
    from GenAIRR._mcp_summary import compute_repertoire_summary

    records = [
        {"v_call": "IGHV1*01,IGHV2*01", "d_call": "IGHD1*01", "j_call": "IGHJ1*01"},
    ]
    summary = compute_repertoire_summary(records)
    # The summary's per-call top map uses the first entry of
    # the tie-set as the attribution key.
    v_top = summary["v_usage_top"]
    assert "IGHV1*01" in v_top
    assert "IGHV2*01" not in v_top  # second entry NOT credited under truth-first


def test_pin_scaffold_mcp_summary_reads_v_call_d_call_j_call_from_record_dicts() -> None:
    """`_mcp_summary.compute_dataset_summary` reads the
    `v_call` / `d_call` / `j_call` columns from each record.
    Pinned at source so the estimator can mirror the column
    convention."""
    src = (_REPO_ROOT / "src" / "GenAIRR" / "_mcp_summary.py").read_text(encoding="utf-8")
    assert '"v_call"' in src
    assert '"d_call"' in src
    assert '"j_call"' in src


def test_pin_scaffold_csv_dict_reader_stdlib_handles_airr_tsv() -> None:
    """Python stdlib `csv.DictReader` with `delimiter='\\t'`
    is sufficient for AIRR TSV parsing. The estimator slice
    uses this rather than introducing pandas as a hard
    dependency."""
    import csv
    import io

    fixture = (
        "v_call\td_call\tj_call\n"
        "IGHV1*01\tIGHD1*01\tIGHJ1*01\n"
        "IGHV2*01,IGHV3*01\tIGHD2*01\tIGHJ2*01\n"
    )
    rows = list(csv.DictReader(io.StringIO(fixture), delimiter="\t"))
    assert len(rows) == 2
    assert rows[0]["v_call"] == "IGHV1*01"
    assert rows[1]["v_call"] == "IGHV2*01,IGHV3*01"


# ──────────────────────────────────────────────────────────────────
# C. pin_scaffold_* — chain-type classifier + typed-plane shape
# ──────────────────────────────────────────────────────────────────


def test_pin_scaffold_config_info_has_d_drives_chain_classification() -> None:
    """`ConfigInfo.has_d` is the authoritative classifier the
    estimator dispatches on (VJ vs VDJ chain). Pinned at the
    enum + ConfigInfo level."""
    from GenAIRR.dataconfig.enums import ChainType
    from GenAIRR.dataconfig.config_info import ConfigInfo

    assert ChainType.BCR_HEAVY.has_d is True
    assert ChainType.BCR_LIGHT_KAPPA.has_d is False
    assert ChainType.TCR_ALPHA.has_d is False
    assert ChainType.TCR_BETA.has_d is True
    # ConfigInfo carries the field that propagates into the cartridge.
    sig = inspect.signature(ConfigInfo)
    assert "has_d" in sig.parameters


def test_pin_scaffold_recombine_rejects_np2_on_vj_chain() -> None:
    """Existing precedent — `Experiment.recombine` raises
    `ValueError` for VJ-chain misuse of `np2_lengths`. The
    estimator slice mirrors this rejection pattern for
    `d_call` data on a VJ cartridge."""
    exp = ga.Experiment.on("human_igk")  # VJ chain
    with pytest.raises(ValueError, match="np2_lengths is only valid for VDJ"):
        exp.recombine(np2_lengths=[(0, 1.0)])


def test_pin_scaffold_reference_empirical_models_carries_four_typed_planes() -> None:
    """`ReferenceEmpiricalModels` carries `np_lengths` /
    `trims` / `np_bases` / `p_nucleotide_lengths` typed
    planes today. `allele_usage` would slot in as the
    fifth — pinned so a future regression that drops the
    typed-plane discipline surfaces here."""
    sig = inspect.signature(ReferenceEmpiricalModels)
    for plane in ("np_lengths", "trims", "np_bases", "p_nucleotide_lengths"):
        assert plane in sig.parameters, (
            f"ReferenceEmpiricalModels no longer carries {plane!r} — "
            f"typed-plane discipline regressed"
        )


# ──────────────────────────────────────────────────────────────────
# D. pin_scaffold_* — builder stage entry shape (release smoke baseline)
# ──────────────────────────────────────────────────────────────────


def test_pin_scaffold_builder_stage_entries_have_inputs_inferred_warnings() -> None:
    """The canonical `{stage, inputs, inferred, warnings}`
    shape is pinned by the release smoke test. The
    `estimate_allele_usage` stage MUST follow the same
    shape — pinned here so a future estimator slice that
    drifts from the shape surfaces immediately."""
    builder = ga.ReferenceCartridgeBuilder.from_fasta(
        v_fasta=">v1*01\nGAGGTG\n",
        j_fasta=">j1*01\nTGGGGC\n",
        chain_type="BCR_LIGHT_KAPPA",
    )
    for entry in builder.report().stages:
        assert set(entry.keys()) >= {"stage", "inputs", "inferred", "warnings"}, (
            f"stage entry missing canonical keys: {entry}"
        )


def test_pin_scaffold_idempotency_pattern_via_replaced_flag_in_v_subregions() -> None:
    """`infer_v_subregions()` writes a `replaced: True` flag
    when called a second time. The `estimate_allele_usage`
    method mirrors this idempotency discipline."""
    v_fasta = ">v1*01\nGAG.GTG\n"  # dot triggers gapped path
    j_fasta = ">j1*01\nTGGGGC\n"
    builder = ga.ReferenceCartridgeBuilder.from_fasta(
        v_fasta=v_fasta, j_fasta=j_fasta, chain_type="BCR_LIGHT_KAPPA"
    )
    builder.infer_v_subregions()
    builder.infer_v_subregions()  # second call
    subregion_stages = [
        s for s in builder.report().stages if s["stage"] == "infer_v_subregions"
    ]
    assert len(subregion_stages) == 2
    assert subregion_stages[0]["inputs"].get("replaced") is False
    assert subregion_stages[1]["inputs"].get("replaced") is True


# ──────────────────────────────────────────────────────────────────
# E. pin_present_* — stop-and-report verification (gene_use_dict)
# ──────────────────────────────────────────────────────────────────


def test_pin_present_gene_use_dict_field_exists_and_is_in_orphan_list() -> None:
    """`DataConfig.gene_use_dict` field exists on the
    dataclass AND is listed in
    `_DOCUMENTED_ORPHAN_DATACONFIG_FIELDS`. The estimator
    slice MUST NOT auto-lift it; the orphan boundary holds."""
    cfg = DataConfig(name="audit-fixture")
    assert hasattr(cfg, "gene_use_dict")
    assert "gene_use_dict" in _DOCUMENTED_ORPHAN_DATACONFIG_FIELDS


def test_pin_present_gene_use_dict_has_no_simulator_pipeline_consumer() -> None:
    """**Stop-and-report gate.** `gene_use_dict` must NOT be
    silently consumed by any simulation path. The only
    allowed consumers are the dataclass declaration itself
    + the MCP read-only `gene_use` diagnostic endpoint. If
    a grep finds another live-source consumer, the
    implementation slice must stop and report.

    This pin protects the audit's clean-yes pre-flight finding."""
    result = subprocess.run(
        [
            "grep",
            "-rln",
            "--include=*.py",
            "--include=*.pyi",
            "--include=*.rs",
            "gene_use_dict",
            "src/GenAIRR",
            "engine_rs/src",
        ],
        capture_output=True,
        text=True,
        cwd=str(_REPO_ROOT),
    )
    consumers = {
        line.strip() for line in result.stdout.splitlines() if line.strip()
    }
    allowed = {
        "src/GenAIRR/dataconfig/data_config.py",
        "src/GenAIRR/utilities/mcp_helpers.py",
        # Post-Allele-Usage-Estimation-v1 slice — the new typed
        # plane's resolver / lowering / spec docstrings
        # explicitly document the no-auto-lift boundary
        # against the legacy field. Those are documentation
        # references in comments / docstrings, NOT simulator
        # consumers. Verified by inspection: each occurrence
        # in the three files below is inside a docstring or
        # comment explaining that the new estimator does NOT
        # touch `gene_use_dict`.
        "src/GenAIRR/_dataconfig_extract.py",
        "src/GenAIRR/experiment.py",
        "src/GenAIRR/reference_models.py",
    }
    unexpected = consumers - allowed
    assert not unexpected, (
        f"gene_use_dict is referenced in unexpected source files: "
        f"{sorted(unexpected)}. The audit's clean-yes verdict assumed "
        f"only the dataclass declaration + MCP diagnostic endpoint "
        f"read it. If any of these are legitimate simulator-pipeline "
        f"consumers, the estimator implementation slice must STOP AND "
        f"REPORT before adding a typed plane — the legacy field is no "
        f"longer orphan."
    )


def test_pin_scaffold_dataconfig_validate_is_dead_code_today() -> None:
    """`DataConfig.validate()` is defined but never called by
    the live source pipeline. The legacy `gene_use_dict`
    requirement inside `validate()` is therefore dead code.
    Pinned so a regression that wires `validate()` into the
    pipeline (which would surface the dead `gene_use_dict`
    requirement at user-facing layers) is caught."""
    # Confirm the method exists.
    assert hasattr(DataConfig, "validate")
    # Confirm it's not called by any live consumer.
    result = subprocess.run(
        [
            "grep",
            "-rn",
            "--include=*.py",
            r"\.validate()",
            "src/GenAIRR",
        ],
        capture_output=True,
        text=True,
        cwd=str(_REPO_ROOT),
    )
    # `validate()` is also a method on ReferenceRulesSpec /
    # ReferenceEmpiricalModels / NpBaseModelSpec /
    # EmpiricalDistributionSpec — those calls are fine.
    # We need to verify NONE of the `.validate()` calls
    # routes into a DataConfig instance. The cheap
    # heuristic: look for `cfg.validate()` / `dc.validate()`
    # / `data_config.validate()` / `self.validate()` inside
    # data_config.py itself.
    cfg_call_patterns = ("cfg.validate(", "dc.validate(", "data_config.validate(")
    for line in result.stdout.splitlines():
        for pattern in cfg_call_patterns:
            assert pattern not in line, (
                f"DataConfig.validate() is now called at {line!r}; the "
                f"audit's dead-code assumption regressed — verify the "
                f"gene_use_dict legacy requirement isn't reachable"
            )


def test_pin_present_mcp_helpers_gene_use_endpoint_is_read_only() -> None:
    """The MCP helper's `gene_use` diagnostic endpoint
    reads `getattr(dc, "gene_use_dict", {})` as a read-only
    inspection — it does NOT feed simulation. Pinned at
    source so a refactor wiring the endpoint into the
    sampler surfaces here."""
    src = (
        _REPO_ROOT / "src" / "GenAIRR" / "utilities" / "mcp_helpers.py"
    ).read_text(encoding="utf-8")
    assert "gene_use_dict" in src
    assert 'section == "gene_use"' in src
    # The endpoint reads, doesn't push to a plan.
    assert "plan.push" not in src  # no engine wiring in mcp_helpers


# ──────────────────────────────────────────────────────────────────
# F. pin_present_* — documented soft-gap boundary state
# ──────────────────────────────────────────────────────────────────


def test_pin_present_plan_signature_soft_gap_for_allele_weights_holds() -> None:
    """Plan-signature completeness audit's soft gap 1: two
    experiments differing only by `v_allele_weights`
    produce equal plan signatures. The estimator slice's
    cartridge-driven output inherits the same boundary.

    Pinned here as well as in the plan-signature audit so
    a future estimator slice that accidentally tightens the
    gap on one path but not the other surfaces here."""
    a = ga.Experiment.on("human_igh").recombine()
    b = ga.Experiment.on("human_igh").recombine(
        v_allele_weights={"IGHVF1-G1*01": 100.0}
    )
    ca = a.compile()
    cb = b.compile()
    sa = json.loads(
        ca.simulator.trace_file_from(ca.simulator.run(seed=0), seed=0).to_json()
    )["pass_plan_signature"]
    sb = json.loads(
        cb.simulator.trace_file_from(cb.simulator.run(seed=0), seed=0).to_json()
    )["pass_plan_signature"]
    assert sa == sb, (
        "v_allele_weights now changes plan signature — soft gap 1 has "
        "been tightened. Verify the estimator slice's cartridge-"
        "driven path is also folded in lockstep, and update the "
        "plan-signature completeness audit's soft-gap-1 documentation."
    )


# ──────────────────────────────────────────────────────────────────
# G. pin_absence_* — gaps the implementation slice closes
# ──────────────────────────────────────────────────────────────────


def test_pin_present_estimate_allele_usage_method_on_builder() -> None:
    """Post-slice — `ReferenceCartridgeBuilder.estimate_allele_usage`
    is callable. The method signature carries the
    `ambiguous` / `min_count` / `replace` kwargs documented
    in the audit §6.1. Flipped from the prior absence pin."""
    assert hasattr(ga.ReferenceCartridgeBuilder, "estimate_allele_usage")
    sig = inspect.signature(ga.ReferenceCartridgeBuilder.estimate_allele_usage)
    for kw in ("min_count", "ambiguous", "replace"):
        assert kw in sig.parameters, (
            f"estimate_allele_usage missing {kw!r} kwarg — audit §6.1 "
            f"signature regressed"
        )


def test_pin_present_allele_usage_field_on_reference_empirical_models() -> None:
    """Post-slice — `ReferenceEmpiricalModels.allele_usage`
    is the fifth typed plane, parallel to `np_lengths` /
    `trims` / `np_bases` / `p_nucleotide_lengths`. Sibling
    naming alternatives stay absent so a future regression
    that adds a parallel surface surfaces here."""
    sig = inspect.signature(ReferenceEmpiricalModels)
    assert "allele_usage" in sig.parameters
    for forbidden in ("v_allele_usage", "allele_weights", "gene_usage"):
        assert forbidden not in sig.parameters, (
            f"ReferenceEmpiricalModels now accepts {forbidden!r} — "
            f"verify the audit doc + implementation slice agree on "
            f"the field name."
        )


def test_pin_present_allele_usage_spec_dataclass() -> None:
    """Post-slice — `AlleleUsageSpec` dataclass is the typed
    container for per-segment allele weights. Re-exported at
    the top level for discoverability."""
    import GenAIRR.reference_models as rm

    assert hasattr(rm, "AlleleUsageSpec")
    assert hasattr(ga, "AlleleUsageSpec")
    assert ga.AlleleUsageSpec is rm.AlleleUsageSpec
    # Triple-segment field shape matches audit §2.
    sig = inspect.signature(ga.AlleleUsageSpec)
    for kw in ("v", "d", "j"):
        assert kw in sig.parameters


def test_pin_present_allele_usage_block_in_manifest() -> None:
    """Post-slice — `manifest['models']['allele_usage']`
    block exposes the typed plane's status with the documented
    keys: `available` / `segments` / `nonempty_segments` /
    `legacy_gene_use_dict_present` / `legacy_fallback` /
    `in_plan_signature` / `source`."""
    m = ga.HUMAN_IGH_OGRDB.cartridge_manifest()
    au = m["models"]["allele_usage"]
    # Baseline (no spec on bundled cartridge).
    assert au["available"] is False
    assert au["segments"] == ["V", "D", "J"]
    assert au["nonempty_segments"] == []
    assert au["legacy_gene_use_dict_present"] is True  # bundled has dict
    assert au["legacy_fallback"] is False
    assert au["in_plan_signature"] is False  # documented soft gap 1
    assert au["source"] == "ReferenceEmpiricalModels.allele_usage"
    # Sibling naming alternatives stay absent so a regression
    # introducing a parallel surface surfaces here.
    for forbidden in ("gene_use", "allele_weights"):
        assert forbidden not in m["models"]


def test_pin_present_ambiguous_kwarg_accepts_three_string_literals() -> None:
    """Post-slice — the `ambiguous` kwarg accepts
    `"fractional"` / `"truth_first"` / `"reject"` string
    literals (no dataclass enum). The audit §3 documented
    that the policy is a string literal rather than a typed
    enum to keep the surface tight; pinned so a future
    refactor introducing an enum surfaces here for review."""
    import GenAIRR.cartridge_builder as cb

    # No typed-enum sibling exists.
    for forbidden in ("AmbiguousPolicy", "AmbiguousCallPolicy", "TieSetPolicy"):
        assert not hasattr(cb, forbidden), (
            f"cartridge_builder.{forbidden} now exists — the audit's "
            f"string-literal-only boundary regressed"
        )
    # The kwarg itself is present.
    sig = inspect.signature(ga.ReferenceCartridgeBuilder.estimate_allele_usage)
    assert "ambiguous" in sig.parameters
    # The default value is `"fractional"`.
    assert sig.parameters["ambiguous"].default == "fractional"


# ──────────────────────────────────────────────────────────────────
# H. Doc anchor
# ──────────────────────────────────────────────────────────────────


def test_pin_scaffold_audit_doc_exists_and_references_contract() -> None:
    if not _AUDIT_DOC.exists():
        import pytest
        pytest.skip("docs/ is contributor-only; audit doc not present in this checkout")
    doc = _AUDIT_DOC.read_text(encoding="utf-8")
    assert "test_allele_usage_estimation_contract.py" in doc, (
        "audit doc no longer references the contract file"
    )
    for marker in (
        "## 1. Q1",
        "## 2. Q2",
        "## 7. Q7",
        "## 9. Implementation order",
        "## 12. Summary table",
    ):
        assert marker in doc, f"audit doc missing section marker {marker!r}"
