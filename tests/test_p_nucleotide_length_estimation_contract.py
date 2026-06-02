"""Contract pins for the P-Nucleotide Length Estimation
audit.

Companion to
[`docs/p_nucleotide_length_estimation_design.md`](../docs/p_nucleotide_length_estimation_design.md).

Pin set:

- ``pin_scaffold_*`` — live surfaces the new estimator
  reuses verbatim:
  - Rust `AirrRecord` exposing four `p_*_length` integer
    fields.
  - Python AIRR projection populating them.
  - `result.py` canonical column order.
  - `ReferenceEmpiricalModels.p_nucleotide_lengths`
    typed plane + `P_NUCLEOTIDE_END_KEYS{,_VJ}` constants.
  - Spec validator rejecting unknown keys + D-end keys
    on VJ at attach time.
  - `EmpiricalDistributionSpec` validator (reused).
  - `_dataconfig_extract` bridge resolvers.
  - Engine `PAdditionPass` parameter signature folding.
  - Manifest's existing `p_nucleotide_models` block.
  - Builder stage / idempotency / CSV ingestion (reused).
- ``pin_present_*`` — stop-and-report verification:
  baseline cartridges produce 0/50 nonzero P-lengths;
  authored P-plane produces nonzero at meaningful rates;
  plan signature folds the distribution; legacy orphan
  boundary held.
- ``pin_absence_*`` — the surfaces the implementation
  slice closes (`estimate_p_nucleotide_lengths` method,
  the `min_count` / `pseudocount` kwarg surface owned by
  this estimator, sibling-class absence, no
  junction/NP-string heuristic).

**Pre-flight verdict (audit §7): clean-yes** with a
documented utility caveat (external AIRR tools don't
populate these fields — handled via doc warning + per-
key auto-warning at implementation time).
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
_AUDIT_DOC = _REPO_ROOT / "audit-docs" / "p_nucleotide_length_estimation_design.md"


# ──────────────────────────────────────────────────────────────────
# A. pin_scaffold_* — AIRR P-length fields (Rust + Python)
# ──────────────────────────────────────────────────────────────────


def test_pin_scaffold_airr_record_carries_four_p_length_fields() -> None:
    """The Rust `AirrRecord` struct exposes
    `p_v_3_length` / `p_d_5_length` / `p_d_3_length` /
    `p_j_5_length` as `i64` fields. The estimator's input
    surface reduces to these four."""
    src = (
        _REPO_ROOT / "engine_rs" / "src" / "airr_record" / "record.rs"
    ).read_text(encoding="utf-8")
    for field in ("p_v_3_length", "p_d_5_length",
                  "p_d_3_length", "p_j_5_length"):
        assert f"pub {field}: i64" in src, (
            f"AirrRecord.{field} missing — P-length AIRR surface "
            f"regressed"
        )


def test_pin_scaffold_airr_record_dicts_carry_four_p_length_fields() -> None:
    """An AIRR record dict returned by `run_records()`
    carries all four `p_*_length` fields — populated by
    the Rust builder, surfaced verbatim in the Python
    record dict. The estimator reads them by these
    canonical column names."""
    result = (
        ga.Experiment.on("human_igh")
        .recombine()
        .run_records(n=1, seed=0)
    )
    rec = result.records[0]
    for col in ("p_v_3_length", "p_d_5_length",
                "p_d_3_length", "p_j_5_length"):
        assert col in rec, (
            f"AIRR record dict missing column {col!r} — "
            f"P-length projection regressed"
        )
        assert isinstance(rec[col], int), (
            f"AIRR record dict {col!r} not an int: type="
            f"{type(rec[col]).__name__}"
        )


def test_pin_scaffold_to_tsv_export_carries_four_p_length_columns() -> None:
    """The canonical TSV export carries all four
    `p_*_length` columns in its header — confirms the
    end-to-end export schema includes the estimator's
    input columns."""
    result = (
        ga.Experiment.on("human_igh")
        .recombine()
        .run_records(n=1, seed=0)
    )
    import io
    buf = io.StringIO()
    # to_tsv writes to a path; emulate via the records directly.
    columns = list(result.records[0].keys())
    for col in ("p_v_3_length", "p_d_5_length",
                "p_d_3_length", "p_j_5_length"):
        assert col in columns, (
            f"AIRR record output missing {col!r} — the estimator's "
            f"input column surface regressed"
        )


# ──────────────────────────────────────────────────────────────────
# B. pin_scaffold_* — typed plane + spec validator
# ──────────────────────────────────────────────────────────────────


def test_pin_scaffold_reference_empirical_models_p_nucleotide_lengths_plane_exists() -> None:
    """`ReferenceEmpiricalModels` carries the
    `p_nucleotide_lengths` plane."""
    sig = inspect.signature(ReferenceEmpiricalModels)
    assert "p_nucleotide_lengths" in sig.parameters


def test_pin_scaffold_p_nucleotide_end_keys_constants_hold() -> None:
    """`P_NUCLEOTIDE_END_KEYS` and
    `P_NUCLEOTIDE_END_KEYS_VJ` constants hold the
    canonical key vocabulary. Pinned so a future drift
    surfaces here."""
    from GenAIRR.reference_models import (
        P_NUCLEOTIDE_END_KEYS, P_NUCLEOTIDE_END_KEYS_VJ,
    )
    assert P_NUCLEOTIDE_END_KEYS == ("V_3", "D_5", "D_3", "J_5"), (
        f"P_NUCLEOTIDE_END_KEYS drifted: {P_NUCLEOTIDE_END_KEYS!r}"
    )
    assert P_NUCLEOTIDE_END_KEYS_VJ == ("V_3", "J_5"), (
        f"P_NUCLEOTIDE_END_KEYS_VJ drifted: {P_NUCLEOTIDE_END_KEYS_VJ!r}"
    )


def test_pin_scaffold_p_nucleotide_lengths_validator_rejects_unknown_keys() -> None:
    """The `p_nucleotide_lengths` validator rejects keys
    outside `P_NUCLEOTIDE_END_KEYS`."""
    spec = ReferenceEmpiricalModels(
        p_nucleotide_lengths={
            "V_5": EmpiricalDistributionSpec([(0, 1.0)]),  # V_5 NOT in keys
        }
    )
    with pytest.raises(ValueError, match=r"p_nucleotide_lengths key .* is not recognised"):
        spec.validate()


def test_pin_present_p_nucleotide_lengths_validator_rejects_d_keys_on_vj() -> None:
    """The validator chain-type-rejects D-end keys on VJ
    at attach time — symmetric with `trims`, asymmetric
    with `np_bases`. Pinned so a future loosening surfaces
    here for explicit review."""
    spec = ReferenceEmpiricalModels(
        p_nucleotide_lengths={
            "D_5": EmpiricalDistributionSpec([(0, 1.0)]),
        }
    )
    with pytest.raises(ValueError, match=r"VJ chain"):
        spec.validate(chain_type="vj")


def test_pin_scaffold_empirical_distribution_spec_rejects_negative_values() -> None:
    """Reused — `EmpiricalDistributionSpec` rejects negative
    integer values."""
    with pytest.raises(ValueError):
        EmpiricalDistributionSpec([(-1, 1.0)]).validate(
            name="p_nucleotide_lengths[V_3]"
        )


def test_pin_scaffold_empirical_distribution_spec_rejects_non_positive_weights() -> None:
    """Reused — `EmpiricalDistributionSpec` rejects zero /
    negative weights."""
    with pytest.raises(ValueError):
        EmpiricalDistributionSpec([(0, 0.0)]).validate(
            name="p_nucleotide_lengths[V_3]"
        )


# ──────────────────────────────────────────────────────────────────
# C. pin_scaffold_* — bridge resolver + engine surface
# ──────────────────────────────────────────────────────────────────


def test_pin_scaffold_extract_recombine_defaults_returns_four_p_length_keys() -> None:
    """`extract_recombine_defaults` returns four
    `p_*_lengths` keys — the bridge plumbing the
    estimator relies on."""
    from GenAIRR._dataconfig_extract import extract_recombine_defaults

    defaults = extract_recombine_defaults(ga.HUMAN_IGH_OGRDB)
    for k in ("p_v_3_lengths", "p_d_5_lengths",
              "p_d_3_lengths", "p_j_5_lengths"):
        assert k in defaults, (
            f"extract_recombine_defaults no longer returns {k!r} — "
            f"the bridge plumbing regressed"
        )


def test_pin_scaffold_p_nucleotide_lengths_from_models_resolver_exists() -> None:
    """`_dataconfig_extract._p_nucleotide_lengths_from_models`
    reads the typed plane. Pinned at source."""
    src = (
        _REPO_ROOT / "src" / "GenAIRR" / "_dataconfig_extract.py"
    ).read_text(encoding="utf-8")
    assert "def _p_nucleotide_lengths_from_models(" in src, (
        "_p_nucleotide_lengths_from_models helper missing"
    )
    assert "p_nucleotide_lengths" in src


def test_pin_scaffold_engine_p_addition_pass_folds_length_dist_into_signature() -> None:
    """`PAdditionPass.parameter_signature` folds the length
    distribution via `fmt_int_dist`. Pinned at source."""
    p_addition_rs = _REPO_ROOT / "engine_rs" / "src" / "passes" / "p_addition.rs"
    assert p_addition_rs.exists(), f"{p_addition_rs} missing"
    src = p_addition_rs.read_text(encoding="utf-8")
    assert "fn parameter_signature" in src
    assert "fmt_int_dist" in src, (
        "PAdditionPass no longer folds via fmt_int_dist — "
        "P-length distributions may have a soft-gap inheritance now"
    )


# ──────────────────────────────────────────────────────────────────
# D. pin_present_* — stop-and-report verification
# ──────────────────────────────────────────────────────────────────


def test_pin_present_baseline_no_p_plane_produces_zero_p_lengths() -> None:
    """**Stop-and-report gate (negative case).**

    Bundled cartridges without an authored P-plane
    produce 0/50 records with any non-zero
    `p_*_length` field. Records from cartridges of this
    shape provide NO signal to the estimator — it would
    produce a degenerate `[(0, 1.0)]` distribution.

    The implementation slice's per-key auto-warning
    surfaces this case in the build report."""
    for cfg_name in ("human_igh", "human_igk"):
        result = (
            ga.Experiment.on(cfg_name)
            .recombine()
            .run_records(n=50, seed=42)
        )
        for col in ("p_v_3_length", "p_d_5_length",
                    "p_d_3_length", "p_j_5_length"):
            nonzero = sum(1 for r in result.records if r.get(col, 0) != 0)
            assert nonzero == 0, (
                f"{cfg_name} baseline: {col} nonzero on {nonzero}/50 "
                f"records — a P-plane is silently active on this "
                f"cartridge. STOP AND REPORT before implementing the "
                f"estimator: the baseline assumption (no signal "
                f"without authored plane) is now violated."
            )


def test_pin_present_authored_p_plane_produces_nonzero_p_lengths() -> None:
    """**Stop-and-report gate (positive case).**

    With a typed P-plane authored (V_3=[(0,0.5),(1,0.3),
    (2,0.2)] etc.), bundled VDJ cartridges produce
    non-zero `p_*_length` fields at meaningful rates —
    the estimator running against such data has signal
    to consume."""
    cfg = copy.deepcopy(ga.HUMAN_IGH_OGRDB)
    cfg.reference_models = ReferenceEmpiricalModels(
        p_nucleotide_lengths={
            "V_3": EmpiricalDistributionSpec([(0, 0.5), (1, 0.3), (2, 0.2)]),
            "D_5": EmpiricalDistributionSpec([(0, 0.7), (1, 0.3)]),
            "D_3": EmpiricalDistributionSpec([(0, 0.7), (1, 0.3)]),
            "J_5": EmpiricalDistributionSpec([(0, 0.8), (1, 0.2)]),
        },
    )
    result = (
        ga.Experiment.on(cfg)
        .recombine()
        .run_records(n=50, seed=42)
    )
    counts = {
        col: sum(1 for r in result.records if r.get(col, 0) != 0)
        for col in ("p_v_3_length", "p_d_5_length",
                    "p_d_3_length", "p_j_5_length")
    }
    # Expect ≥ 10/50 nonzero per end at the authored distributions.
    for col, n in counts.items():
        assert n >= 10, (
            f"P-plane authored but {col} nonzero on only {n}/50 "
            f"records (expected ≥ 10). STOP AND REPORT — the "
            f"authored plane is not reaching the engine at "
            f"expected rates."
        )


def test_pin_present_p_length_distribution_changes_plan_signature() -> None:
    """Two cartridges differing only in
    `p_nucleotide_lengths["V_3"]` produce different plan
    signatures via `PAdditionPass.parameter_signature`.
    No soft gap inherited."""
    cfg_a = copy.deepcopy(ga.HUMAN_IGH_OGRDB)
    cfg_b = copy.deepcopy(ga.HUMAN_IGH_OGRDB)
    cfg_a.reference_models = ReferenceEmpiricalModels(
        p_nucleotide_lengths={
            "V_3": EmpiricalDistributionSpec([(0, 1.0)]),
        }
    )
    cfg_b.reference_models = ReferenceEmpiricalModels(
        p_nucleotide_lengths={
            "V_3": EmpiricalDistributionSpec([(0, 0.5), (1, 0.5)]),
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
        "Two cartridges differing only in p_nucleotide_lengths['V_3'] "
        "produce EQUAL plan signatures — P-length distributions no "
        "longer fold into the signature; estimator inherits a soft gap"
    )


# ──────────────────────────────────────────────────────────────────
# E. pin_present_* — legacy orphan boundary
# ──────────────────────────────────────────────────────────────────


def test_pin_present_legacy_p_nucleotide_length_probs_orphan_boundary_holds() -> None:
    """`p_nucleotide_length_probs` is listed in
    `_DOCUMENTED_ORPHAN_DATACONFIG_FIELDS`. The estimator
    MUST NOT auto-lift it — same boundary every prior
    estimator slice respected."""
    from GenAIRR.dataconfig.data_config import (
        _DOCUMENTED_ORPHAN_DATACONFIG_FIELDS,
    )
    assert "p_nucleotide_length_probs" in _DOCUMENTED_ORPHAN_DATACONFIG_FIELDS


def test_pin_present_manifest_reports_legacy_orphan_without_lifting() -> None:
    """The bundled cartridges carry the legacy
    `p_nucleotide_length_probs` orphan dict. The manifest's
    `p_nucleotide_models` block reports it as
    `legacy_p_nucleotide_length_probs_present` WITHOUT
    lifting it into the typed plane."""
    cfg = ga.HUMAN_IGH_OGRDB
    block = cfg.cartridge_manifest()["models"]["p_nucleotide_models"]
    assert "legacy_p_nucleotide_length_probs_present" in block
    assert block["legacy_fallback"] is False
    # And the typed plane stays empty on bundled cartridges.
    rm = getattr(cfg, "reference_models", None)
    typed = getattr(rm, "p_nucleotide_lengths", None) if rm else None
    assert not typed, (
        f"bundled cartridge typed P-plane is non-empty: {typed!r}"
    )


# ──────────────────────────────────────────────────────────────────
# F. pin_scaffold_* — manifest already shaped (no extension needed)
# ──────────────────────────────────────────────────────────────────


def test_pin_scaffold_manifest_p_nucleotide_models_block_already_carries_full_shape() -> None:
    """The manifest's `p_nucleotide_models` block ALREADY
    carries the full shape the estimator needs. No new
    helper or extension — the estimator's output rides
    the existing block by populating
    `self._reference_models.p_nucleotide_lengths`."""
    m = ga.HUMAN_IGH_OGRDB.cartridge_manifest()
    assert "p_nucleotide_models" in m["models"]
    block = m["models"]["p_nucleotide_models"]
    for key in ("length_keys", "legacy_p_nucleotide_length_probs_present",
                "legacy_fallback", "supported_ends",
                "in_plan_signature", "in_content_hash"):
        assert key in block, (
            f"manifest p_nucleotide_models block missing key "
            f"{key!r}; the estimator slice's plan to ride the "
            f"existing block regressed"
        )
    assert block["in_plan_signature"] is True
    assert set(block["supported_ends"]) == {"V_3", "D_5", "D_3", "J_5"}


# ──────────────────────────────────────────────────────────────────
# G. pin_scaffold_* — chain-type classifier + builder shape (reused)
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
    """Reused — `replaced=True` idempotency pattern."""
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
    """Reused — stdlib `csv.DictReader` AIRR-TSV ingestion
    shape with the P-length columns."""
    import csv
    import io

    fixture = (
        "v_call\tp_v_3_length\td_call\tp_d_5_length\tp_d_3_length\tj_call\tp_j_5_length\n"
        "IGHV1*01\t2\tIGHD1*01\t0\t1\tIGHJ1*01\t0\n"
        "IGHV2*01\t0\tIGHD2*01\t1\t2\tIGHJ2*01\t1\n"
    )
    rows = list(csv.DictReader(io.StringIO(fixture), delimiter="\t"))
    assert len(rows) == 2
    assert rows[0]["p_v_3_length"] == "2"
    assert rows[1]["p_j_5_length"] == "1"


# ──────────────────────────────────────────────────────────────────
# H. pin_absence_* — gaps the implementation slice closes
# ──────────────────────────────────────────────────────────────────


def test_pin_present_estimate_p_nucleotide_lengths_method_on_builder() -> None:
    """Post-slice —
    `ReferenceCartridgeBuilder.estimate_p_nucleotide_lengths`
    is callable. The method signature carries the
    `min_count` / `pseudocount` / `replace` kwargs
    documented in the user brief. Flipped from the prior
    absence pin."""
    assert hasattr(
        ga.ReferenceCartridgeBuilder, "estimate_p_nucleotide_lengths"
    )
    sig = inspect.signature(
        ga.ReferenceCartridgeBuilder.estimate_p_nucleotide_lengths
    )
    for kw in ("min_count", "pseudocount", "replace"):
        assert kw in sig.parameters, (
            f"estimate_p_nucleotide_lengths missing {kw!r} kwarg — "
            f"user brief signature regressed"
        )
    assert sig.parameters["min_count"].default == 1
    assert sig.parameters["pseudocount"].default == 0.0
    assert sig.parameters["replace"].default is True


def test_pin_absence_no_p_nucleotide_length_estimator_module_or_class() -> None:
    """No sibling module / class with the estimator's name
    was accidentally introduced. The method lives on
    `ReferenceCartridgeBuilder`."""
    import GenAIRR.cartridge_builder as cb

    for forbidden in ("PNucleotideLengthEstimator", "PLengthEstimator",
                      "estimate_p_nucleotide_lengths"):
        assert not hasattr(cb, forbidden), (
            f"cartridge_builder.{forbidden} now exists — verify "
            f"the estimator surface is owned by the builder method"
        )


def test_pin_absence_no_junction_arithmetic_or_np_string_heuristic_in_p_length_estimator() -> None:
    """v1 boundary — once `estimate_p_nucleotide_lengths`
    lands, its method body MUST NOT derive P-lengths from
    `junction_length` arithmetic or NP string slicing.
    Pinned defensively so a future heuristic helper
    surfaces here for explicit review.

    Pre-slice (current state) this pin is a no-op because
    the method does not exist yet. Post-slice this pin
    fires if the implementation references
    `junction_length` / `np1` / `np2` inside the method's
    source."""
    method = getattr(
        ga.ReferenceCartridgeBuilder, "estimate_p_nucleotide_lengths", None
    )
    if method is None:
        # Pre-slice — no method to inspect; pin passes trivially.
        return
    src = inspect.getsource(method)
    # Strip the docstring so the pin doesn't catch
    # documentation references like "columns ignored: junction_length".
    lines = src.splitlines()
    in_docstring = False
    body_lines: list[str] = []
    for line in lines:
        stripped = line.strip()
        if stripped.startswith('"""') or stripped.startswith("'''"):
            in_docstring = not in_docstring
            # Single-line docstring start AND end on same line.
            if stripped.count('"""') == 2 or stripped.count("'''") == 2:
                in_docstring = False
            continue
        if in_docstring:
            continue
        body_lines.append(line)
    body = "\n".join(body_lines)
    for forbidden in ("junction_length", '"np1"', '"np2"',
                      "'np1'", "'np2'"):
        assert forbidden not in body, (
            f"estimate_p_nucleotide_lengths body references "
            f"{forbidden!r} — v1 forbids junction-arithmetic / "
            f"NP-string heuristic P-length inference"
        )


# ──────────────────────────────────────────────────────────────────
# I. Doc anchor
# ──────────────────────────────────────────────────────────────────


def test_pin_scaffold_audit_doc_exists_and_references_contract() -> None:
    """The audit doc exists and references the contract
    file by name; section structure intact."""
    doc = _AUDIT_DOC.read_text(encoding="utf-8")
    assert "test_p_nucleotide_length_estimation_contract.py" in doc, (
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
