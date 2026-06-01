"""Contract pins for the NP Base Model Estimation audit.

Companion to
[`docs/np_base_model_estimation_design.md`](../docs/np_base_model_estimation_design.md).

Pin set:

- ``pin_scaffold_*`` — live surfaces the new estimator
  reuses verbatim:
  - Rust `AirrRecord` exposing `np1` / `np2` String fields.
  - `projection.rs::unclaimed_np_string` walking the
    structural NP region only.
  - `result.py` canonical column order including both.
  - `ReferenceEmpiricalModels.np_bases` typed plane.
  - `NpBaseModelSpec` supporting `uniform` /
    `empirical_first_base` / `markov`.
  - Spec validator rejecting non-A/C/G/T bases +
    incomplete Markov transition rows.
  - `_dataconfig_extract` bridge resolvers
    (`_np_bases_from_models` +
    `_np_markov_transitions_from_models`).
  - Engine `GenerateNPPass.parameter_signature` folding
    via canonical base weights.
  - Manifest's existing perfectly-shaped
    `np_base_models` block.
  - Builder stage entry shape + idempotency + CSV ingestion.
- ``pin_present_*`` — **the critical stop-and-report
  verification**: `len(np1) == np1_length` and
  `len(np2) == np2_length` even under a max-P plane.
  P-bytes are NOT contaminating the AIRR strings.
- ``pin_present_*`` — Markov wired end-to-end (the
  spec docstring's "deferred" claim is stale).
- ``pin_present_*`` — plan signature folds NP base model
  (no soft gap inherited).
- ``pin_present_*`` — legacy orphan boundary held.
- ``pin_absence_*`` — the surfaces the implementation
  slice closes (`estimate_np_base_model` method, `kind`
  / `pseudocount` kwarg surface, sibling class).

**Pre-flight verdict (audit §7): clean-yes.** Both
anticipated tricky points (P-contamination, Markov
status) resolved cleanly at audit time.
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
    NpBaseModelSpec,
    ReferenceEmpiricalModels,
)


_REPO_ROOT = Path(__file__).resolve().parent.parent
_AUDIT_DOC = _REPO_ROOT / "docs" / "np_base_model_estimation_design.md"


# ──────────────────────────────────────────────────────────────────
# A. pin_scaffold_* — AIRR np1/np2 fields (Rust + Python)
# ──────────────────────────────────────────────────────────────────


def test_pin_scaffold_airr_record_carries_np1_and_np2_string_fields() -> None:
    """The Rust `AirrRecord` struct exposes `np1` and `np2`
    as `String` fields — the estimator's input columns."""
    src = (
        _REPO_ROOT / "engine_rs" / "src" / "airr_record" / "record.rs"
    ).read_text(encoding="utf-8")
    assert "pub np1: String" in src
    assert "pub np2: String" in src


def test_pin_scaffold_unclaimed_np_string_walks_structural_region_only() -> None:
    """`projection.rs::unclaimed_np_string` walks
    `region.start..region.end` over the STRUCTURAL NP
    region. By construction, the range excludes P-byte
    ranges (P bytes don't add to `sim.sequence.regions`).
    Pinned at source so a future refactor that
    accidentally widens the walk surfaces here."""
    src = (
        _REPO_ROOT / "engine_rs" / "src" / "airr_record" / "projection.rs"
    ).read_text(encoding="utf-8")
    assert "pub(super) fn unclaimed_np_string(" in src
    # The walk iterates region.start..region.end of the supplied region.
    assert "let r_start = region.start.index()" in src
    assert "let r_end = region.end.index()" in src
    assert "for i in r_start..r_end" in src


def test_pin_scaffold_result_column_order_includes_np1_and_np2_strings() -> None:
    """`result.py`'s canonical column order declares both
    `np1` and `np2` string fields. The estimator inherits
    the column names directly."""
    src = (
        _REPO_ROOT / "src" / "GenAIRR" / "result.py"
    ).read_text(encoding="utf-8")
    assert '"np1"' in src
    assert '"np2"' in src


# ──────────────────────────────────────────────────────────────────
# B. pin_scaffold_* — typed plane + NpBaseModelSpec
# ──────────────────────────────────────────────────────────────────


def test_pin_scaffold_reference_empirical_models_np_bases_plane_exists() -> None:
    """`ReferenceEmpiricalModels` carries `np_bases`."""
    sig = inspect.signature(ReferenceEmpiricalModels)
    assert "np_bases" in sig.parameters


def test_pin_scaffold_np_base_model_spec_supports_three_kinds() -> None:
    """`NpBaseModelSpec` recognises three kinds — `uniform`,
    `empirical_first_base`, `markov`. The estimator's
    `kind` kwarg accepts `empirical_first_base` and
    `markov` (uniform is not data-derived)."""
    from GenAIRR.reference_models import _NP_BASE_MODEL_KINDS

    assert set(_NP_BASE_MODEL_KINDS) == {
        "uniform", "empirical_first_base", "markov"
    }, (
        f"_NP_BASE_MODEL_KINDS drifted: {_NP_BASE_MODEL_KINDS!r}; "
        f"audit §2.2 documented the three-kind vocabulary"
    )


def test_pin_scaffold_np_base_model_spec_validates_first_base_alphabet() -> None:
    """The spec validator rejects non-A/C/G/T bases in the
    `first_base` distribution. The estimator's
    `noncanonical_base` rejection logic feeds clean A/C/G/T
    counts in before construction."""
    with pytest.raises(ValueError):
        NpBaseModelSpec(
            kind="empirical_first_base",
            first_base={"A": 1.0, "X": 1.0, "G": 1.0, "T": 1.0},  # X invalid
        )


def test_pin_scaffold_np_base_model_spec_validates_markov_row_coverage() -> None:
    """The spec validator rejects a Markov transition matrix
    that lacks any of the four A/C/G/T from-bases. Pinned so
    a future loosening of the row-coverage check (which would
    let a partial estimator output through) surfaces here."""
    with pytest.raises(ValueError):
        NpBaseModelSpec(
            kind="markov",
            first_base={"A": 1.0, "C": 1.0, "G": 1.0, "T": 1.0},
            transitions={
                "A": {"A": 0.25, "C": 0.25, "G": 0.25, "T": 0.25},
                "C": {"A": 0.25, "C": 0.25, "G": 0.25, "T": 0.25},
                # G + T rows missing → spec must reject.
            },
        )


# ──────────────────────────────────────────────────────────────────
# C. pin_scaffold_* — bridge resolver + lowering
# ──────────────────────────────────────────────────────────────────


def test_pin_scaffold_extract_recombine_defaults_consumes_np_bases_plane() -> None:
    """`extract_recombine_defaults` returns `np1_bases` /
    `np2_bases` keys. The estimator's output flows into
    these resolver-produced defaults verbatim."""
    from GenAIRR._dataconfig_extract import extract_recombine_defaults

    defaults = extract_recombine_defaults(ga.HUMAN_IGH_OGRDB)
    for k in ("np1_bases", "np2_bases"):
        assert k in defaults, (
            f"extract_recombine_defaults no longer returns {k!r} — "
            f"the bridge plumbing the estimator relies on regressed"
        )


def test_pin_scaffold_np_bases_from_models_resolver_exists() -> None:
    """`_dataconfig_extract._np_bases_from_models` reads
    the typed plane. Pinned at source."""
    src = (
        _REPO_ROOT / "src" / "GenAIRR" / "_dataconfig_extract.py"
    ).read_text(encoding="utf-8")
    assert "def _np_bases_from_models(" in src
    assert "def _np_markov_transitions_from_models(" in src, (
        "Markov transitions resolver missing — the kind='markov' "
        "estimator output would not reach the engine"
    )


def test_pin_scaffold_engine_generate_np_pass_folds_base_payload_into_signature() -> None:
    """`GenerateNPPass.parameter_signature` folds the base
    generator's payload — both first-base row and (for
    markov) the 4×4 transition matrix. Pinned at source so
    a fold-optimisation that elides constant generators
    surfaces here."""
    src = (
        _REPO_ROOT / "engine_rs" / "src" / "passes" / "generate_np.rs"
    ).read_text(encoding="utf-8")
    assert "fn parameter_signature" in src
    # The base generator's signature contribution is folded in.
    assert "base_generator" in src or "base_pairs" in src, (
        "GenerateNPPass.parameter_signature no longer includes the "
        "base generator payload — replay protection for NP base "
        "models regressed"
    )


# ──────────────────────────────────────────────────────────────────
# D. pin_present_* — STOP-AND-REPORT: P-byte contamination check
# ──────────────────────────────────────────────────────────────────


def test_pin_present_np1_length_matches_len_np1_string_under_max_p_plane() -> None:
    """**The critical P-cross-contamination gate.**

    With a max-P plane authored (3 P-bytes at each of 4
    ends, 12 P-bytes per record), `len(rec["np1"])` MUST
    equal `rec["np1_length"]` and `len(rec["np2"])` MUST
    equal `rec["np2_length"]` on every record. If P-bytes
    leaked into the AIRR strings, the lengths would
    diverge.

    The estimator's plan to use `rec["np1"]` / `rec["np2"]`
    verbatim as P-clean A/C/G/T streams depends on this
    invariant. If it breaks, the implementation slice must
    STOP AND REPORT and add P-aware filtering."""
    cfg = copy.deepcopy(ga.HUMAN_IGH_OGRDB)
    cfg.reference_models = ReferenceEmpiricalModels(
        p_nucleotide_lengths={
            "V_3": EmpiricalDistributionSpec([(3, 1.0)]),
            "D_5": EmpiricalDistributionSpec([(3, 1.0)]),
            "D_3": EmpiricalDistributionSpec([(3, 1.0)]),
            "J_5": EmpiricalDistributionSpec([(3, 1.0)]),
        },
    )
    result = (
        ga.Experiment.on(cfg)
        .recombine()
        .run_records(n=50, seed=4242)
    )
    np1_mismatch = sum(
        1 for r in result.records
        if len(r.get("np1", "")) != r.get("np1_length", 0)
    )
    np2_mismatch = sum(
        1 for r in result.records
        if len(r.get("np2", "")) != r.get("np2_length", 0)
    )
    assert np1_mismatch == 0, (
        f"len(np1) != np1_length on {np1_mismatch}/50 records under a "
        f"max-P plane — P-bytes are NOW contaminating the np1 string. "
        f"STOP AND REPORT before implementing the estimator; the slice "
        f"must add P-aware filtering."
    )
    assert np2_mismatch == 0, (
        f"len(np2) != np2_length on {np2_mismatch}/50 records under a "
        f"max-P plane — P-bytes are NOW contaminating the np2 string. "
        f"STOP AND REPORT."
    )


def test_pin_present_bundled_cartridges_produce_canonical_acgt_in_np_strings() -> None:
    """The bundled cartridges produce only canonical A/C/G/T
    characters in `np1` / `np2`. The estimator's
    `noncanonical_base` rejection logic is exercised only
    by user inputs from external simulators, not by
    GenAIRR's own AIRR records."""
    canonical = set("ACGT")
    for cfg_name in ("human_igh", "human_igk"):
        result = (
            ga.Experiment.on(cfg_name)
            .recombine()
            .run_records(n=50, seed=42)
        )
        for r in result.records:
            for col in ("np1", "np2"):
                value = r.get(col, "") or ""
                non_canonical = set(value) - canonical
                assert not non_canonical, (
                    f"{cfg_name}: record carries non-A/C/G/T characters "
                    f"in {col}: {non_canonical!r}"
                )


# ──────────────────────────────────────────────────────────────────
# E. pin_present_* — Markov wired end-to-end (NOT deferred)
# ──────────────────────────────────────────────────────────────────


def test_pin_present_markov_np_base_model_reaches_engine_and_biases_output() -> None:
    """Despite the stale "Markov deferred" docstring at
    `reference_models.py:145-158`, the engine
    `MarkovBaseGenerator` is wired end-to-end. A forced
    cyclic Markov matrix (A→T, T→C, C→G, G→A) reproduces
    on the majority of observed transitions in NP1 output.

    If this regresses, the estimator's `kind="markov"`
    output would not bias the engine — STOP AND REPORT
    before implementing."""
    cfg = copy.deepcopy(ga.HUMAN_IGH_OGRDB)
    cfg.reference_models = ReferenceEmpiricalModels(
        np_bases={
            "NP1": NpBaseModelSpec(
                kind="markov",
                first_base={"A": 1.0, "C": 1.0, "G": 1.0, "T": 1.0},
                transitions={
                    "A": {"T": 1.0, "A": 0.01, "C": 0.01, "G": 0.01},
                    "C": {"G": 1.0, "A": 0.01, "C": 0.01, "T": 0.01},
                    "G": {"A": 1.0, "C": 0.01, "G": 0.01, "T": 0.01},
                    "T": {"C": 1.0, "A": 0.01, "G": 0.01, "T": 0.01},
                },
            )
        }
    )
    result = (
        ga.Experiment.on(cfg)
        .recombine()
        .run_records(n=50, seed=42)
    )
    transitions: dict[str, int] = {}
    for r in result.records:
        np1 = r.get("np1", "")
        for i in range(len(np1) - 1):
            pair = np1[i:i+2]
            transitions[pair] = transitions.get(pair, 0) + 1
    # Forced pairs should dominate.
    forced = ("AT", "TC", "CG", "GA")
    forced_count = sum(transitions.get(p, 0) for p in forced)
    total = sum(transitions.values()) or 1
    forced_ratio = forced_count / total
    assert forced_ratio >= 0.8, (
        f"Markov model failed to bias engine output: forced "
        f"transitions {forced} captured {forced_ratio:.2%} of "
        f"observed pairs (expected ≥ 80%). STOP AND REPORT — "
        f"the `kind='markov'` estimator output would not reach "
        f"the engine."
    )


# ──────────────────────────────────────────────────────────────────
# F. pin_present_* — plan signature + chain-type
# ──────────────────────────────────────────────────────────────────


def test_pin_present_np_base_distribution_changes_plan_signature() -> None:
    """Two cartridges differing only in `np_bases["NP1"]`
    produce different plan signatures. No soft gap
    inherited."""
    cfg_a = copy.deepcopy(ga.HUMAN_IGH_OGRDB)
    cfg_b = copy.deepcopy(ga.HUMAN_IGH_OGRDB)
    cfg_a.reference_models = ReferenceEmpiricalModels(
        np_bases={
            "NP1": NpBaseModelSpec(
                kind="empirical_first_base",
                first_base={"A": 1.0, "C": 1.0, "G": 1.0, "T": 1.0},
            )
        }
    )
    cfg_b.reference_models = ReferenceEmpiricalModels(
        np_bases={
            "NP1": NpBaseModelSpec(
                kind="empirical_first_base",
                first_base={"A": 1.0, "C": 0.1, "G": 1.0, "T": 1.0},
            )
        }
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
        "Two cartridges differing only in np_bases['NP1'] produce "
        "EQUAL plan signatures — NP base distributions no longer "
        "fold into the signature; the estimator inherits a soft gap"
    )


def test_pin_present_vj_np2_string_is_empty_on_bundled_cartridge() -> None:
    """`np2` is the empty string on every record of a VJ
    cartridge — no NP2 region. The estimator MUST NOT
    write an `NP2` key to a VJ cartridge's plane."""
    result = (
        ga.Experiment.on("human_igk")
        .recombine()
        .run_records(n=50, seed=42)
    )
    for r in result.records:
        assert r["np2"] == "", (
            f"VJ cartridge: np2 non-empty — engine pipeline now runs "
            f"an NP2 generation pass on VJ; the audit's hard-empty "
            f"boundary regressed"
        )


# ──────────────────────────────────────────────────────────────────
# G. pin_present_* — legacy orphan boundary
# ──────────────────────────────────────────────────────────────────


def test_pin_present_legacy_np_first_bases_and_transitions_are_orphan() -> None:
    """`NP_first_bases` and `NP_transitions` are listed in
    `_DOCUMENTED_ORPHAN_DATACONFIG_FIELDS`. The estimator
    MUST NOT auto-lift them; the orphan boundary holds —
    same discipline every prior estimator slice respected."""
    from GenAIRR.dataconfig.data_config import (
        _DOCUMENTED_ORPHAN_DATACONFIG_FIELDS,
    )
    assert "NP_first_bases" in _DOCUMENTED_ORPHAN_DATACONFIG_FIELDS
    assert "NP_transitions" in _DOCUMENTED_ORPHAN_DATACONFIG_FIELDS


def test_pin_present_legacy_np_first_bases_does_not_auto_lift_to_typed_plane() -> None:
    """The bridge resolver does NOT lift the legacy
    `NP_first_bases` / `NP_transitions` dicts into the
    typed plane. Verified by `_np_bases_from_models`
    returning `None` for a bundled cartridge (which carries
    the legacy dicts but no typed spec)."""
    from GenAIRR._dataconfig_extract import _np_bases_from_models

    cfg = ga.HUMAN_IGH_OGRDB
    # Bundled cartridge has the legacy dicts populated...
    assert getattr(cfg, "NP_first_bases", None), (
        "bundled cartridge no longer carries legacy NP_first_bases — "
        "this pin's premise regressed"
    )
    # ...but no typed cartridge_models with np_bases authored.
    rm = getattr(cfg, "reference_models", None)
    # The resolver returns None for NP1 / NP2 — no auto-lift happens.
    assert _np_bases_from_models(rm, "NP1") is None
    assert _np_bases_from_models(rm, "NP2") is None


# ──────────────────────────────────────────────────────────────────
# H. pin_scaffold_* — manifest already shaped (no extension needed)
# ──────────────────────────────────────────────────────────────────


def test_pin_scaffold_manifest_np_base_models_block_already_carries_full_shape() -> None:
    """Unlike the trim / NP-length slices, the manifest's
    `np_base_models` block ALREADY carries the full shape
    the estimator needs:

    - `models: Dict[key, {"kind": str}]`
    - `supported_kinds`, `deferred_kinds`
    - `legacy_fallback`, `legacy_np_*_present` flags
    - `in_plan_signature`, `in_content_hash`

    The estimator slice's output (a populated
    `np_bases` plane) flows through this existing block
    without any manifest helper change."""
    m = ga.HUMAN_IGH_OGRDB.cartridge_manifest()
    assert "np_base_models" in m["models"]
    block = m["models"]["np_base_models"]
    for key in ("models", "supported_kinds", "deferred_kinds",
                "legacy_fallback", "legacy_np_transitions_present",
                "legacy_np_first_bases_present",
                "in_plan_signature", "in_content_hash"):
        assert key in block, (
            f"manifest np_base_models block missing key {key!r}; the "
            f"estimator slice's plan to ride the existing block "
            f"regressed"
        )
    assert set(block["supported_kinds"]) == {
        "uniform", "empirical_first_base", "markov"
    }
    assert block["in_plan_signature"] is True


# ──────────────────────────────────────────────────────────────────
# I. pin_scaffold_* — chain-type classifier + builder shape (reused)
# ──────────────────────────────────────────────────────────────────


def test_pin_scaffold_config_info_has_d_drives_chain_classification() -> None:
    """`ConfigInfo.has_d` — the authoritative classifier."""
    from GenAIRR.dataconfig.enums import ChainType
    from GenAIRR.dataconfig.config_info import ConfigInfo

    assert ChainType.BCR_HEAVY.has_d is True
    assert ChainType.BCR_LIGHT_KAPPA.has_d is False
    sig = inspect.signature(ConfigInfo)
    assert "has_d" in sig.parameters


def test_pin_scaffold_builder_stage_entries_have_inputs_inferred_warnings() -> None:
    """Reused — canonical `{stage, inputs, inferred,
    warnings}` shape."""
    builder = ga.ReferenceCartridgeBuilder.from_fasta(
        v_fasta=">v1*01\nGAGGTG\n",
        j_fasta=">j1*01\nTGGGGC\n",
        chain_type="BCR_LIGHT_KAPPA",
    )
    for entry in builder.report().stages:
        assert set(entry.keys()) >= {"stage", "inputs", "inferred", "warnings"}


def test_pin_scaffold_idempotency_pattern_via_replaced_flag_in_v_subregions() -> None:
    """Reused — `replaced=True` idempotency."""
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
    """Reused — `csv.DictReader` AIRR-TSV ingestion shape
    that includes the `np1` / `np2` columns."""
    import csv
    import io

    fixture = (
        "v_call\tnp1\td_call\tnp2\tj_call\n"
        "IGHV1*01\tACGT\tIGHD1*01\tTTTT\tIGHJ1*01\n"
        "IGHV2*01\t\tIGHD2*01\tGGGA\tIGHJ2*01\n"
    )
    rows = list(csv.DictReader(io.StringIO(fixture), delimiter="\t"))
    assert len(rows) == 2
    assert rows[0]["np1"] == "ACGT"
    assert rows[1]["np2"] == "GGGA"


# ──────────────────────────────────────────────────────────────────
# J. pin_absence_* — gaps the implementation slice closes
# ──────────────────────────────────────────────────────────────────


def test_pin_present_estimate_np_base_model_method_on_builder() -> None:
    """Post-slice —
    `ReferenceCartridgeBuilder.estimate_np_base_model`
    is callable. The method signature carries the
    `kind` / `min_count` / `pseudocount` / `replace`
    kwargs documented in the user brief. Flipped from
    the prior absence pin."""
    assert hasattr(
        ga.ReferenceCartridgeBuilder, "estimate_np_base_model"
    )
    sig = inspect.signature(
        ga.ReferenceCartridgeBuilder.estimate_np_base_model
    )
    for kw in ("kind", "min_count", "pseudocount", "replace"):
        assert kw in sig.parameters, (
            f"estimate_np_base_model missing {kw!r} kwarg — user "
            f"brief signature regressed"
        )
    assert sig.parameters["kind"].default == "markov"
    assert sig.parameters["min_count"].default == 1
    assert sig.parameters["pseudocount"].default == 0.0
    assert sig.parameters["replace"].default is True


def test_pin_present_kind_kwarg_owned_by_np_base_model_estimator() -> None:
    """Post-slice — `kind` is owned by the NP-base-model
    estimator only. Pinned so a future slice that
    smuggles `kind` onto a sibling method (where it would
    have different semantics) surfaces here."""
    method_owners_kind: list[str] = []
    for method_name in dir(ga.ReferenceCartridgeBuilder):
        if method_name.startswith("_"):
            continue
        method = getattr(ga.ReferenceCartridgeBuilder, method_name, None)
        if not callable(method):
            continue
        try:
            sig = inspect.signature(method)
        except (TypeError, ValueError):
            continue
        if "kind" in sig.parameters:
            method_owners_kind.append(method_name)
    assert method_owners_kind == ["estimate_np_base_model"], (
        f"`kind` kwarg owned by methods other than "
        f"estimate_np_base_model: {method_owners_kind}"
    )


def test_pin_absence_no_np_base_model_estimator_module_or_class() -> None:
    """No sibling module / class with the estimator's name
    was accidentally introduced. The method lives on
    `ReferenceCartridgeBuilder`."""
    import GenAIRR.cartridge_builder as cb

    for forbidden in ("NpBaseModelEstimator", "MarkovEstimator",
                      "estimate_np_base_model"):
        assert not hasattr(cb, forbidden), (
            f"cartridge_builder.{forbidden} now exists — verify "
            f"the estimator surface is owned by the builder method"
        )


# ──────────────────────────────────────────────────────────────────
# K. Doc anchor
# ──────────────────────────────────────────────────────────────────


def test_pin_scaffold_audit_doc_exists_and_references_contract() -> None:
    """The audit doc exists and references the contract
    file by name; section structure intact."""
    assert _AUDIT_DOC.exists(), "np_base_model_estimation_design.md missing"
    doc = _AUDIT_DOC.read_text(encoding="utf-8")
    assert "test_np_base_model_estimation_contract.py" in doc, (
        "audit doc no longer references the contract file"
    )
    for marker in (
        "## 1. Q1",
        "## 2. Q2",
        "## 7. Q7",
        "## 10. Implementation order",
        "## 13. Summary table",
    ):
        assert marker in doc, f"audit doc missing section marker {marker!r}"
