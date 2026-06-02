"""Contract pins for the NP Length Distribution Estimation
audit.

Companion to
[`docs/np_length_estimation_design.md`](../docs/np_length_estimation_design.md).

Pin set:

- ``pin_scaffold_*`` — live surfaces the new estimator
  reuses verbatim:
  - Rust `AirrRecord` exposing `np1_length` / `np2_length`
    integer fields + `np1` / `np2` sequence fields.
  - Python AIRR projection populating both length fields.
  - `result.py` canonical column order including both.
  - `ReferenceEmpiricalModels.np_lengths` typed plane +
    `NP_KEYS` constants + key-membership validator.
  - `_dataconfig_extract.extract_recombine_defaults` +
    `_np_lengths_from_models` resolver.
  - `Experiment.recombine(np1_lengths=..., np2_lengths=...)`
    kwargs + VJ-NP2 rejection.
  - Engine `GenerateNPPass` parameter signature folding
    via `fmt_int_dist`.
  - Builder stage entry shape + idempotency.
  - `csv.DictReader` AIRR-TSV ingestion.
- ``pin_scaffold_*`` — provenance distinction: NP lengths
  vs P-nucleotide lengths at every layer.
- ``pin_present_*`` — stop-and-report verification:
  bundled cartridges populate `np1_length` / `np2_length`
  at significant rates; VJ chains hard-zero `np2_length`.
- ``pin_present_*`` — documented surface state: NP-length
  distributions fold into the plan signature (NO inherited
  soft gap).
- ``pin_present_*`` — documented asymmetry: the typed-plane
  validator does NOT chain-type-reject `NP2` on VJ at
  attach time (asymmetric with the `trims` plane). The
  estimator must enforce the boundary at estimation time.
- ``pin_absence_*`` — the surfaces the implementation
  slice closes (`estimate_np_length_distributions` method,
  `np_length_models` manifest block, sibling-class/module).

**Pre-flight verdict (audit §7): clean-yes.** The
estimator can land without disturbing any existing surface.
"""
from __future__ import annotations

import copy
import inspect
import json
from pathlib import Path

import pytest

import GenAIRR as ga
from GenAIRR.reference_models import (
    EmpiricalDistributionSpec,
    ReferenceEmpiricalModels,
)


_REPO_ROOT = Path(__file__).resolve().parent.parent
_AUDIT_DOC = _REPO_ROOT / "audit-docs" / "np_length_estimation_design.md"


# ──────────────────────────────────────────────────────────────────
# A. pin_scaffold_* — AIRR NP fields (Rust + Python)
# ──────────────────────────────────────────────────────────────────


def test_pin_scaffold_airr_record_carries_np1_and_np2_length_fields() -> None:
    """The Rust `AirrRecord` struct exposes `np1_length` and
    `np2_length` as `i64` fields. The estimator's input
    columns reduce to these two."""
    src = (
        _REPO_ROOT / "engine_rs" / "src" / "airr_record" / "record.rs"
    ).read_text(encoding="utf-8")
    assert "pub np1_length: i64" in src
    assert "pub np2_length: i64" in src


def test_pin_scaffold_airr_record_carries_np1_and_np2_sequence_fields() -> None:
    """The Rust `AirrRecord` also exposes `np1` / `np1_aa`
    / `np2` / `np2_aa` String fields. The estimator
    deliberately does NOT consume these — direct integer
    length fields are canonical (per audit §1.2)."""
    src = (
        _REPO_ROOT / "engine_rs" / "src" / "airr_record" / "record.rs"
    ).read_text(encoding="utf-8")
    for field in ("np1", "np1_aa", "np2", "np2_aa"):
        assert f"pub {field}: String" in src, (
            f"AirrRecord.{field} missing — NP sequence surface "
            f"regressed"
        )


def test_pin_scaffold_python_airr_projection_emits_np1_length_np2_length() -> None:
    """The Python AIRR projection at `_airr_record.py`
    populates both length fields from region span. Pinned
    at source — the estimator inherits the same field
    names."""
    src = (
        _REPO_ROOT / "src" / "GenAIRR" / "_airr_record.py"
    ).read_text(encoding="utf-8")
    assert '"np1_length"' in src
    assert '"np2_length"' in src


def test_pin_scaffold_result_column_order_includes_np_length_fields() -> None:
    """`result.py`'s canonical column order declares both
    length fields. The estimator's output is consumed by
    `Experiment.on(cfg).recombine().run_records()` which
    will emit AIRR rows carrying these columns."""
    src = (
        _REPO_ROOT / "src" / "GenAIRR" / "result.py"
    ).read_text(encoding="utf-8")
    assert '"np1_length"' in src
    assert '"np2_length"' in src


# ──────────────────────────────────────────────────────────────────
# B. pin_scaffold_* — typed plane + validator
# ──────────────────────────────────────────────────────────────────


def test_pin_scaffold_reference_empirical_models_np_lengths_plane_exists() -> None:
    """`ReferenceEmpiricalModels` carries `np_lengths` — the
    existing typed surface the estimator writes into."""
    sig = inspect.signature(ReferenceEmpiricalModels)
    assert "np_lengths" in sig.parameters


def test_pin_scaffold_np_keys_constant_holds() -> None:
    """`NP_KEYS = ("NP1", "NP2")` is the canonical
    NP-length-key vocabulary. Pinned so a future drift adds
    a key without the matching engine pass."""
    from GenAIRR.reference_models import NP_KEYS

    assert NP_KEYS == ("NP1", "NP2"), (
        f"NP_KEYS drifted: {NP_KEYS!r}; audit §2.1 documented the "
        f"two-key vocabulary"
    )


def test_pin_scaffold_np_lengths_validator_rejects_unknown_keys() -> None:
    """The `np_lengths` plane validator rejects keys outside
    `NP_KEYS`."""
    spec = ReferenceEmpiricalModels(
        np_lengths={"NP3": EmpiricalDistributionSpec([(0, 1.0)])}
    )
    with pytest.raises(ValueError, match=r"np_lengths key .* is not recognised"):
        spec.validate()


def test_pin_present_np_lengths_validator_does_not_reject_np2_on_vj() -> None:
    """**Documented asymmetry vs `trims`.** The current
    `np_lengths` plane validator only checks key membership
    in `NP_KEYS` — it does NOT chain-type-reject `NP2` on
    a VJ cartridge at attach time. The trim validator DOES
    reject `D_5` / `D_3` on VJ. Pinned so a future
    tightening that closes this asymmetry surfaces here
    for explicit review (and updates both surfaces in
    lockstep)."""
    spec = ReferenceEmpiricalModels(
        np_lengths={"NP2": EmpiricalDistributionSpec([(0, 1.0)])}
    )
    # No exception — the spec passes validation on a VJ
    # chain even though the cartridge will refuse to compile.
    spec.validate(chain_type="vj")
    # The estimator must therefore enforce the boundary at
    # estimation time, not rely on the validator.


def test_pin_scaffold_empirical_distribution_spec_rejects_negative_values() -> None:
    """Reused — the per-spec validator rejects negative
    integer values. The estimator's per-row parser catches
    these first, but the spec validator is the final gate."""
    with pytest.raises(ValueError):
        EmpiricalDistributionSpec([(-1, 1.0)]).validate(name="np_lengths[NP1]")


def test_pin_scaffold_empirical_distribution_spec_rejects_non_positive_weights() -> None:
    """Reused — the per-spec validator rejects zero / negative
    weights."""
    with pytest.raises(ValueError):
        EmpiricalDistributionSpec([(0, 0.0)]).validate(name="np_lengths[NP1]")
    with pytest.raises(ValueError):
        EmpiricalDistributionSpec([(0, -1.0)]).validate(name="np_lengths[NP1]")


# ──────────────────────────────────────────────────────────────────
# C. pin_scaffold_* — bridge resolver + lowering precedence
# ──────────────────────────────────────────────────────────────────


def test_pin_scaffold_extract_recombine_defaults_consumes_np_lengths_plane() -> None:
    """`extract_recombine_defaults` returns `np1` / `np2`
    keys with the typed-plane resolver taking precedence
    over the legacy nested-dict path."""
    from GenAIRR._dataconfig_extract import extract_recombine_defaults

    defaults = extract_recombine_defaults(ga.HUMAN_IGH_OGRDB)
    for k in ("np1", "np2"):
        assert k in defaults, (
            f"extract_recombine_defaults no longer returns {k!r} — "
            f"the bridge plumbing the estimator relies on regressed"
        )


def test_pin_scaffold_np_lengths_from_models_resolver_is_typed_plane_adapter() -> None:
    """`_dataconfig_extract._np_lengths_from_models` reads
    the typed plane and returns a flat `[(int, float), ...]`
    list. Pinned at source."""
    src = (
        _REPO_ROOT / "src" / "GenAIRR" / "_dataconfig_extract.py"
    ).read_text(encoding="utf-8")
    assert "def _np_lengths_from_models(" in src, (
        "_np_lengths_from_models helper missing — typed-plane "
        "lowering for np_lengths regressed"
    )
    assert "models.np_lengths.get(key)" in src


def test_pin_scaffold_typed_plane_takes_precedence_over_legacy_np_lengths() -> None:
    """The bridge precedence is typed plane > legacy nested
    `cfg.NP_lengths` > None. Verified by direct comparison:
    a cartridge with both a typed plane AND legacy dict
    produces the typed pairs from
    `extract_recombine_defaults`."""
    from GenAIRR._dataconfig_extract import extract_recombine_defaults

    cfg = copy.deepcopy(ga.HUMAN_IGH_OGRDB)
    # Sanity: legacy NP_lengths populated on the bundled cartridge.
    assert cfg.NP_lengths, "bundled cartridge has no legacy NP_lengths"
    cfg.reference_models = ReferenceEmpiricalModels(
        np_lengths={"NP1": EmpiricalDistributionSpec([(0, 1.0)])}
    )
    defaults = extract_recombine_defaults(cfg)
    typed = defaults["np1"]
    assert typed is not None
    assert typed == [(0, 1.0)], (
        f"np1 didn't honour the typed plane override: {typed!r}; "
        f"the precedence the estimator relies on regressed"
    )


def test_pin_scaffold_experiment_recombine_accepts_np1_lengths_and_np2_lengths_kwargs() -> None:
    """`Experiment.recombine(np1_lengths=..., np2_lengths=...)`
    is the per-experiment override surface — priority 1
    in the precedence chain. The estimator's typed-plane
    output is at priority 2."""
    sig = inspect.signature(ga.Experiment.recombine)
    for kw in ("np1_lengths", "np2_lengths"):
        assert kw in sig.parameters, (
            f"Experiment.recombine no longer accepts {kw!r} — "
            f"the NP-length DSL override surface regressed"
        )


def test_pin_scaffold_recombine_rejects_np2_lengths_on_vj() -> None:
    """`Experiment.recombine(np2_lengths=...)` raises
    `ValueError` on VJ chains — confirmed at source. The
    estimator inherits this safety net: a VJ cartridge
    with an `NP2` plane key would crash at recombine
    compile time."""
    exp = ga.Experiment.on("human_igk")  # VJ chain
    with pytest.raises(ValueError, match="np2_lengths is only valid for VDJ"):
        exp.recombine(np2_lengths=[(0, 1.0)])


# ──────────────────────────────────────────────────────────────────
# D. pin_scaffold_* — engine surface (GenerateNPPass + paramsig)
# ──────────────────────────────────────────────────────────────────


def test_pin_scaffold_engine_generate_np_pass_folds_length_dist_into_signature() -> None:
    """`GenerateNPPass.parameter_signature` folds the
    length distribution via `fmt_int_dist`. Pinned at
    source so a future fold-optimisation that elides
    constant distributions surfaces here."""
    src = (
        _REPO_ROOT / "engine_rs" / "src" / "passes" / "generate_np.rs"
    ).read_text(encoding="utf-8")
    assert "fn parameter_signature" in src
    assert "fmt_int_dist" in src, (
        "GenerateNPPass no longer folds via fmt_int_dist — "
        "NP-length distributions may have a soft-gap inheritance now"
    )


# ──────────────────────────────────────────────────────────────────
# E. pin_scaffold_* — provenance distinction (NP vs P-nucleotide)
# ──────────────────────────────────────────────────────────────────


def test_pin_scaffold_p_nucleotide_length_fields_distinct_from_np_length_fields() -> None:
    """The Rust AirrRecord exposes `p_v_3_length` /
    `p_d_5_length` / `p_d_3_length` / `p_j_5_length` as
    SEPARATE fields from `np1_length` / `np2_length`. The
    estimator MUST NOT consume `p_*_length` columns."""
    src = (
        _REPO_ROOT / "engine_rs" / "src" / "airr_record" / "record.rs"
    ).read_text(encoding="utf-8")
    for p_field in ("p_v_3_length", "p_d_5_length",
                    "p_d_3_length", "p_j_5_length"):
        assert f"pub {p_field}: i64" in src
    # And both NP fields live alongside, not in place of.
    for np_field in ("np1_length", "np2_length"):
        assert f"pub {np_field}: i64" in src


def test_pin_scaffold_p_addition_pass_module_distinct_from_generate_np() -> None:
    """`engine_rs/src/passes/p_addition.rs` is the
    P-nucleotide engine pass, distinct from
    `passes/generate_np.rs`. The two surfaces share the
    EmpiricalDistributionSpec container shape but are
    biologically distinct."""
    p_addition_rs = (
        _REPO_ROOT / "engine_rs" / "src" / "passes" / "p_addition.rs"
    )
    generate_np_rs = (
        _REPO_ROOT / "engine_rs" / "src" / "passes" / "generate_np.rs"
    )
    assert p_addition_rs.exists(), f"{p_addition_rs} missing"
    assert generate_np_rs.exists(), f"{generate_np_rs} missing"


def test_pin_scaffold_p_nucleotide_lengths_plane_distinct_from_np_lengths() -> None:
    """`ReferenceEmpiricalModels` carries
    `p_nucleotide_lengths` as a SEPARATE plane from
    `np_lengths`."""
    sig = inspect.signature(ReferenceEmpiricalModels)
    assert "np_lengths" in sig.parameters
    assert "p_nucleotide_lengths" in sig.parameters
    assert sig.parameters["np_lengths"].name != sig.parameters[
        "p_nucleotide_lengths"
    ].name


# ──────────────────────────────────────────────────────────────────
# F. pin_scaffold_* — chain-type classifier + builder shape
# ──────────────────────────────────────────────────────────────────


def test_pin_scaffold_config_info_has_d_drives_chain_classification() -> None:
    """`ConfigInfo.has_d` — the authoritative classifier."""
    from GenAIRR.dataconfig.enums import ChainType
    from GenAIRR.dataconfig.config_info import ConfigInfo

    assert ChainType.BCR_HEAVY.has_d is True
    assert ChainType.BCR_LIGHT_KAPPA.has_d is False
    assert ChainType.TCR_ALPHA.has_d is False
    assert ChainType.TCR_BETA.has_d is True
    sig = inspect.signature(ConfigInfo)
    assert "has_d" in sig.parameters


def test_pin_scaffold_builder_stage_entries_have_inputs_inferred_warnings() -> None:
    """Reused — canonical `{stage, inputs, inferred,
    warnings}` shape from the release smoke."""
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
    """Reused — `replaced=True` flag for idempotent stages."""
    builder = ga.ReferenceCartridgeBuilder.from_fasta(
        v_fasta=">v1*01\nGAG.GTG\n",
        j_fasta=">j1*01\nTGGGGC\n",
        chain_type="BCR_LIGHT_KAPPA",
    )
    builder.infer_v_subregions()
    builder.infer_v_subregions()
    stages = [s for s in builder.report().stages if s["stage"] == "infer_v_subregions"]
    assert len(stages) == 2
    assert stages[0]["inputs"].get("replaced") is False
    assert stages[1]["inputs"].get("replaced") is True


def test_pin_scaffold_csv_dict_reader_stdlib_handles_airr_tsv() -> None:
    """Reused — `csv.DictReader` AIRR-TSV ingestion."""
    import csv
    import io

    fixture = (
        "v_call\tnp1_length\td_call\tnp2_length\tj_call\n"
        "IGHV1*01\t3\tIGHD1*01\t5\tIGHJ1*01\n"
        "IGHV2*01\t1\tIGHD2*01\t0\tIGHJ2*01\n"
    )
    rows = list(csv.DictReader(io.StringIO(fixture), delimiter="\t"))
    assert len(rows) == 2
    assert rows[0]["np1_length"] == "3"
    assert rows[1]["np2_length"] == "0"


# ──────────────────────────────────────────────────────────────────
# G. pin_present_* — stop-and-report verification
# ──────────────────────────────────────────────────────────────────


def test_pin_present_genairr_populates_vdj_np_length_fields_reliably() -> None:
    """**Stop-and-report gate.** Bundled VDJ cartridge
    `HUMAN_IGH_OGRDB` populates `np1_length` and
    `np2_length` at significant rates over a 100-record
    fixed-seed smoke."""
    result = (
        ga.Experiment.on("human_igh")
        .recombine()
        .run_records(n=100, seed=42)
    )
    np1_nonzero = sum(1 for r in result.records if r["np1_length"] > 0)
    np2_nonzero = sum(1 for r in result.records if r["np2_length"] > 0)
    assert np1_nonzero >= 50, (
        f"VDJ cartridge populates np1_length on only {np1_nonzero}/100 "
        f"records (expected ≥ 50) — STOP AND REPORT before "
        f"implementing the estimator"
    )
    assert np2_nonzero >= 50, (
        f"VDJ cartridge populates np2_length on only {np2_nonzero}/100 "
        f"records (expected ≥ 50) — STOP AND REPORT"
    )


def test_pin_present_genairr_populates_vj_np1_length_field_reliably() -> None:
    """**Stop-and-report gate, VJ side.** Bundled VJ
    cartridge `HUMAN_IGK_OGRDB` populates `np1_length` at
    significant rates; `np2_length` is hard-zero (no NP2
    region on a VJ chain)."""
    result = (
        ga.Experiment.on("human_igk")
        .recombine()
        .run_records(n=100, seed=42)
    )
    np1_nonzero = sum(1 for r in result.records if r["np1_length"] > 0)
    assert np1_nonzero >= 50, (
        f"VJ cartridge populates np1_length on only {np1_nonzero}/100 "
        f"records (expected ≥ 50) — STOP AND REPORT"
    )


def test_pin_present_vj_np2_length_is_hard_zero() -> None:
    """`np2_length` is hard-zero on every record of a VJ
    cartridge — confirms the documented "no NP2 on VJ"
    boundary. The estimator MUST NOT write an `NP2` key
    to a VJ cartridge's plane."""
    result = (
        ga.Experiment.on("human_igk")
        .recombine()
        .run_records(n=100, seed=42)
    )
    for r in result.records:
        assert r["np2_length"] == 0, (
            f"VJ cartridge: np2_length nonzero — engine pipeline "
            f"now runs an NP2 pass on VJ, audit's hard-zero "
            f"boundary regressed"
        )


# ──────────────────────────────────────────────────────────────────
# H. pin_present_* — documented surface state (NO soft-gap inherited)
# ──────────────────────────────────────────────────────────────────


def test_pin_present_np_length_distributions_fold_into_plan_signature() -> None:
    """Two cartridges differing only in `np_lengths["NP1"]`
    produce different plan signatures via
    `GenerateNPPass.parameter_signature`. The estimator's
    output benefits from cross-cartridge replay protection
    automatically — no inherited soft gap (unlike
    allele-usage)."""
    cfg_a = copy.deepcopy(ga.HUMAN_IGH_OGRDB)
    cfg_b = copy.deepcopy(ga.HUMAN_IGH_OGRDB)
    cfg_a.reference_models = ReferenceEmpiricalModels(
        np_lengths={"NP1": EmpiricalDistributionSpec([(0, 1.0)])}
    )
    cfg_b.reference_models = ReferenceEmpiricalModels(
        np_lengths={"NP1": EmpiricalDistributionSpec([(0, 0.5), (1, 0.5)])}
    )
    exp_a = ga.Experiment.on(cfg_a).recombine()
    exp_b = ga.Experiment.on(cfg_b).recombine()
    ca = exp_a.compile()
    cb = exp_b.compile()
    sa = json.loads(
        ca.simulator.trace_file_from(ca.simulator.run(seed=0), seed=0).to_json()
    )["pass_plan_signature"]
    sb = json.loads(
        cb.simulator.trace_file_from(cb.simulator.run(seed=0), seed=0).to_json()
    )["pass_plan_signature"]
    assert sa != sb, (
        "Two cartridges differing only in np_lengths['NP1'] now "
        "produce EQUAL plan signatures — NP-length distributions "
        "no longer fold into the signature; the estimator slice "
        "now inherits a soft gap"
    )


# ──────────────────────────────────────────────────────────────────
# I. pin_scaffold_* — manifest already exposes minimal NP block
# ──────────────────────────────────────────────────────────────────


def test_pin_scaffold_manifest_exposes_np_length_keys_today() -> None:
    """The bundled cartridges already expose minimal
    `np_length_keys` (list) and `legacy_np_lengths_present`
    (bool). The implementation slice extends this into a
    structured `np_length_models` block per audit §8.2 —
    pinned here so the existing minimal entries remain
    backwards compatible."""
    m = ga.HUMAN_IGH_OGRDB.cartridge_manifest()
    assert "np_length_keys" in m["models"]
    assert isinstance(m["models"]["np_length_keys"], list)
    assert "legacy_np_lengths_present" in m["models"]
    assert m["models"]["legacy_np_lengths_present"] is True


# ──────────────────────────────────────────────────────────────────
# J. pin_absence_* — gaps the implementation slice closes
# ──────────────────────────────────────────────────────────────────


def test_pin_present_estimate_np_length_distributions_method_on_builder() -> None:
    """Post-slice —
    `ReferenceCartridgeBuilder.estimate_np_length_distributions`
    is callable. The method signature carries the
    `min_count` / `pseudocount` / `replace` kwargs
    documented in the user brief. Flipped from the prior
    absence pin."""
    assert hasattr(
        ga.ReferenceCartridgeBuilder, "estimate_np_length_distributions"
    )
    sig = inspect.signature(
        ga.ReferenceCartridgeBuilder.estimate_np_length_distributions
    )
    for kw in ("min_count", "pseudocount", "replace"):
        assert kw in sig.parameters, (
            f"estimate_np_length_distributions missing {kw!r} kwarg — "
            f"user brief signature regressed"
        )
    assert sig.parameters["min_count"].default == 1
    assert sig.parameters["pseudocount"].default == 0.0
    assert sig.parameters["replace"].default is True


def test_pin_present_np_length_models_block_in_manifest() -> None:
    """Post-slice —
    `manifest['models']['np_length_models']` block exposes
    the typed plane's status with the documented keys:
    `keys` / `source` / `in_plan_signature` /
    `legacy_np_lengths_present` / `legacy_fallback`.
    Flipped from the prior absence pin."""
    m = ga.HUMAN_IGH_OGRDB.cartridge_manifest()
    assert "np_length_models" in m["models"], (
        "manifest['models']['np_length_models'] block missing — "
        "user brief specified this block name; verify the manifest "
        "extension landed"
    )
    block = m["models"]["np_length_models"]
    # Baseline (no typed-plane spec on bundled cartridge).
    assert block["keys"] == []
    assert block["source"] == "ReferenceEmpiricalModels.np_lengths"
    # NO inherited soft gap — distinguishes from allele-usage.
    assert block["in_plan_signature"] is True
    assert block["legacy_np_lengths_present"] is True
    assert block["legacy_fallback"] is False
    # Minimal top-level entries stay for backwards compatibility.
    assert "np_length_keys" in m["models"]
    assert "legacy_np_lengths_present" in m["models"]


def test_pin_absence_no_np_length_estimator_module_or_class() -> None:
    """No sibling module / class with the estimator's name
    was accidentally introduced. The method lives on
    `ReferenceCartridgeBuilder`, not as a free function."""
    import GenAIRR.cartridge_builder as cb

    for forbidden in ("NPLengthEstimator", "NpLengthEstimator",
                      "estimate_np_length_distributions"):
        assert not hasattr(cb, forbidden), (
            f"cartridge_builder.{forbidden} now exists — verify "
            f"the estimator surface is owned by the builder method, "
            f"not a parallel free-function / dataclass surface"
        )


# ──────────────────────────────────────────────────────────────────
# K. Doc anchor
# ──────────────────────────────────────────────────────────────────


def test_pin_scaffold_audit_doc_exists_and_references_contract() -> None:
    """The audit doc exists and references the contract
    file by name; section structure intact."""
    doc = _AUDIT_DOC.read_text(encoding="utf-8")
    assert "test_np_length_estimation_contract.py" in doc, (
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
