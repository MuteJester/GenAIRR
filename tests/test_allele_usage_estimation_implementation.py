"""End-to-end implementation tests for the **Allele Usage
Estimation v1** slice.

Companion to
[`docs/allele_usage_estimation_design.md`](../docs/allele_usage_estimation_design.md)
and the contract pins in
[`tests/test_allele_usage_estimation_contract.py`](test_allele_usage_estimation_contract.py).

Covers the 15-test surface from the user brief:

 1. `AlleleUsageSpec` validation accepts valid weights and
    rejects bad names / weights.
 2. Explicit typed plane affects recombination sampling.
 3. Explicit `recombine(v_allele_weights=...)` kwarg overrides
    the typed plane.
 4. `gene_use_dict` remains orphan / no auto-lift.
 5. Builder `ambiguous="fractional"` splits tie-set credit.
 6. Builder `ambiguous="truth_first"` uses first call only.
 7. Builder `ambiguous="reject"` records ambiguous rows in
    `report.rejected`.
 8. Unknown alleles recorded in `report.rejected`.
 9. VJ chain ignores D contribution with a warning.
10. VDJ chain missing D recorded in `report.rejected`.
11. `min_count` drops low-support alleles.
12. Stage entry shape matches audit §6.2.
13. Manifest `models.allele_usage` block reflects estimator
    output.
14. Pickle round-trip preserves estimated plane + report.
15. Contract absence pins flipped; soft-gap pin holds.
"""
from __future__ import annotations

import copy
import json
import math
import pickle

import pytest

import GenAIRR as ga
from GenAIRR.reference_models import (
    AlleleUsageSpec,
    ReferenceEmpiricalModels,
)


# ──────────────────────────────────────────────────────────────────
# Shared FASTA fixtures (synthesised from the bundled HUMAN_IGK_OGRDB
# cartridge so allele anchors / sequences are realistic without
# needing the native C resolver).
# ──────────────────────────────────────────────────────────────────


def _fastas_and_names(src: "ga.DataConfig", n_v: int = 3, n_j: int = 2):
    v_names = [list(src.v_alleles.values())[i][0].name for i in range(n_v)]
    j_names = [list(src.j_alleles.values())[i][0].name for i in range(n_j)]
    v_fasta = "".join(
        f">{a.name}\n{a.gapped_seq}\n"
        for gene_alleles in list(src.v_alleles.values())[:n_v]
        for a in gene_alleles[:1]
    )
    j_fasta = "".join(
        f">{a.name}\n{a.ungapped_seq}\n"
        for gene_alleles in list(src.j_alleles.values())[:n_j]
        for a in gene_alleles[:1]
    )
    return v_fasta, j_fasta, v_names, j_names


def _builder_vj(src=None):
    src = src or ga.HUMAN_IGK_OGRDB
    v_fasta, j_fasta, v_names, j_names = _fastas_and_names(src)
    builder = ga.ReferenceCartridgeBuilder.from_fasta(
        v_fasta=v_fasta, j_fasta=j_fasta, chain_type="BCR_LIGHT_KAPPA"
    )
    return builder, v_names, j_names


# ──────────────────────────────────────────────────────────────────
# 1. AlleleUsageSpec validation accepts valid + rejects bad
# ──────────────────────────────────────────────────────────────────


def test_allele_usage_spec_validates_and_rejects_bad_inputs() -> None:
    """Valid spec passes; bad name shape / non-positive weight /
    D-on-VJ all raise `ValueError` with field-tagged messages."""
    spec = AlleleUsageSpec(v={"IGHV1*01": 0.7, "IGHV2*01": 0.3})
    spec.validate(chain_type="vj")  # OK — only V populated

    # Empty name rejected.
    with pytest.raises(ValueError, match="non-empty strings"):
        AlleleUsageSpec(v={"": 1.0}).validate()
    # Non-positive weight rejected.
    with pytest.raises(ValueError, match="strictly positive"):
        AlleleUsageSpec(v={"X*01": 0.0}).validate()
    # Non-finite weight rejected.
    with pytest.raises(ValueError, match="finite"):
        AlleleUsageSpec(v={"X*01": math.inf}).validate()
    # D-on-VJ rejected with the standard cartridge boundary.
    with pytest.raises(ValueError, match="VJ chain"):
        AlleleUsageSpec(d={"IGHD1*01": 1.0}).validate(chain_type="vj")


# ──────────────────────────────────────────────────────────────────
# 2. Explicit typed plane biases recombination sampling
# ──────────────────────────────────────────────────────────────────


def test_explicit_typed_plane_biases_recombination_sampling() -> None:
    """Attaching `allele_usage` on `reference_models` shifts the
    V-call distribution the recombination pass produces. Strict
    one-allele restriction (weight 100, everyone else 1.0)
    makes the boosted allele dominate the v_call tie sets."""
    cfg = copy.deepcopy(ga.HUMAN_IGH_OGRDB)
    v_names = [a.name for fam in cfg.v_alleles.values() for a in fam]
    boosted = v_names[0]
    cfg.reference_models = ReferenceEmpiricalModels(
        allele_usage=AlleleUsageSpec(v={boosted: 100.0})
    )
    exp = ga.Experiment.on(cfg).recombine()
    result = exp.run_records(n=40, seed=4242)
    boosted_hits = sum(
        1 for r in result.records if boosted in (r.get("v_call") or "")
    )
    assert boosted_hits >= 10, (
        f"boosted V allele {boosted!r} only appeared in {boosted_hits}/40 "
        f"records; the cartridge plane's bias isn't reaching the engine"
    )


# ──────────────────────────────────────────────────────────────────
# 3. Explicit recombine kwarg overrides typed plane
# ──────────────────────────────────────────────────────────────────


def test_explicit_recombine_kwarg_overrides_typed_plane() -> None:
    """Per audit §2 precedence: kwarg > cartridge plane > uniform.
    With BOTH the typed plane biasing allele A AND the kwarg
    biasing allele B, the engine must sample biased toward B."""
    cfg = copy.deepcopy(ga.HUMAN_IGH_OGRDB)
    v_names = [a.name for fam in cfg.v_alleles.values() for a in fam]
    plane_boost = v_names[0]
    kwarg_boost = v_names[5]
    cfg.reference_models = ReferenceEmpiricalModels(
        allele_usage=AlleleUsageSpec(v={plane_boost: 100.0})
    )
    exp = ga.Experiment.on(cfg).recombine(
        v_allele_weights={kwarg_boost: 100.0}
    )
    result = exp.run_records(n=40, seed=4242)
    kwarg_hits = sum(
        1 for r in result.records if kwarg_boost in (r.get("v_call") or "")
    )
    plane_hits = sum(
        1 for r in result.records if plane_boost in (r.get("v_call") or "")
    )
    assert kwarg_hits >= plane_hits, (
        f"kwarg override should dominate plane bias; got kwarg_hits="
        f"{kwarg_hits}, plane_hits={plane_hits}"
    )
    assert kwarg_hits >= 10


# ──────────────────────────────────────────────────────────────────
# 4. gene_use_dict remains orphan / no auto-lift
# ──────────────────────────────────────────────────────────────────


def test_gene_use_dict_remains_orphan_under_estimator_slice() -> None:
    """The bundled cartridge carries a populated `gene_use_dict`.
    Attaching an `allele_usage` spec MUST NOT auto-lift the
    legacy dict; the orphan boundary holds."""
    from GenAIRR.dataconfig.data_config import (
        _DOCUMENTED_ORPHAN_DATACONFIG_FIELDS,
    )
    assert "gene_use_dict" in _DOCUMENTED_ORPHAN_DATACONFIG_FIELDS
    # Even after attaching a spec, the legacy dict is untouched.
    cfg = copy.deepcopy(ga.HUMAN_IGH_OGRDB)
    pre = dict(cfg.gene_use_dict)
    v_names = [a.name for fam in cfg.v_alleles.values() for a in fam]
    cfg.reference_models = ReferenceEmpiricalModels(
        allele_usage=AlleleUsageSpec(v={v_names[0]: 100.0})
    )
    assert cfg.gene_use_dict == pre, "gene_use_dict was silently mutated"


# ──────────────────────────────────────────────────────────────────
# 5. Builder fractional policy splits tie-set credit
# ──────────────────────────────────────────────────────────────────


def test_fractional_policy_splits_tieset_credit_evenly() -> None:
    """A row `v_call="A,B"` under `ambiguous="fractional"`
    credits 0.5 to A and 0.5 to B (when both are known)."""
    builder, v, j = _builder_vj()
    records = [
        {"v_call": f"{v[0]},{v[1]}", "j_call": j[0]},  # split 0.5 / 0.5
        {"v_call": f"{v[0]},{v[1]}", "j_call": j[0]},  # another split
    ]
    builder.estimate_allele_usage(records, ambiguous="fractional", min_count=0.5)
    stage = next(
        s for s in builder.report().stages
        if s["stage"] == "estimate_allele_usage"
    )
    v_weights = stage["inferred"]["V"]
    assert set(v_weights.keys()) == {v[0], v[1]}
    assert abs(v_weights[v[0]] - 0.5) < 1e-9
    assert abs(v_weights[v[1]] - 0.5) < 1e-9


# ──────────────────────────────────────────────────────────────────
# 6. Truth-first uses first call only
# ──────────────────────────────────────────────────────────────────


def test_truth_first_uses_first_call_only() -> None:
    """`ambiguous="truth_first"` collapses each tie set to its
    first allele. The output weights should reflect only the
    first-position calls."""
    builder, v, j = _builder_vj()
    records = [
        {"v_call": f"{v[0]},{v[1]}", "j_call": j[0]},  # → v[0] only
        {"v_call": f"{v[1]},{v[0]}", "j_call": j[0]},  # → v[1] only
        {"v_call": f"{v[1]},{v[0]}", "j_call": j[0]},  # → v[1] only
    ]
    builder.estimate_allele_usage(records, ambiguous="truth_first", min_count=0.5)
    stage = next(
        s for s in builder.report().stages
        if s["stage"] == "estimate_allele_usage"
    )
    v_weights = stage["inferred"]["V"]
    # v[0] appears once as first-call; v[1] appears twice.
    assert abs(v_weights[v[0]] - 1 / 3) < 1e-9
    assert abs(v_weights[v[1]] - 2 / 3) < 1e-9


# ──────────────────────────────────────────────────────────────────
# 7. Reject policy records ambiguous rows in report.rejected
# ──────────────────────────────────────────────────────────────────


def test_reject_policy_records_ambiguous_rows_in_rejected() -> None:
    """`ambiguous="reject"` drops every multi-call row and emits
    a structured `ambiguous_rejected` entry per row in the
    report's rejected list."""
    builder, v, j = _builder_vj()
    records = [
        {"v_call": f"{v[0]},{v[1]}", "j_call": j[0]},  # ambiguous → rejected
        {"v_call": v[0], "j_call": j[0]},              # unambiguous → kept
        {"v_call": v[1], "j_call": f"{j[0]},{j[1]}"},  # ambiguous J → rejected
    ]
    builder.estimate_allele_usage(records, ambiguous="reject", min_count=0.5)
    stage = next(
        s for s in builder.report().stages
        if s["stage"] == "estimate_allele_usage"
    )
    skipped = stage["inferred"]["skipped"]
    assert skipped["ambiguous_rejected"] == 2
    # Only the unambiguous record contributes — v[0] gets all weight.
    assert stage["inferred"]["V"] == {v[0]: 1.0}
    # Rejected entries surface in the report top-level list.
    ambig = [
        r for r in builder.report().rejected
        if r.get("reason") == "ambiguous_rejected"
    ]
    assert len(ambig) == 2


# ──────────────────────────────────────────────────────────────────
# 8. Unknown alleles recorded in report.rejected
# ──────────────────────────────────────────────────────────────────


def test_unknown_alleles_recorded_in_rejected_entries() -> None:
    """An AIRR record naming an allele NOT in the cartridge's
    pool drops a structured rejection entry naming the allele
    + segment + reason."""
    builder, v, j = _builder_vj()
    records = [
        {"v_call": "IGKV-NONEXISTENT*01", "j_call": j[0]},
        {"v_call": v[0], "j_call": "IGKJ-MISSING*99"},
    ]
    builder.estimate_allele_usage(records, min_count=0.5)
    unknown_v = [
        r for r in builder.report().rejected
        if r.get("reason") == "unknown_allele" and r.get("segment") == "V"
    ]
    unknown_j = [
        r for r in builder.report().rejected
        if r.get("reason") == "unknown_allele" and r.get("segment") == "J"
    ]
    assert len(unknown_v) == 1
    assert unknown_v[0]["allele_name"] == "IGKV-NONEXISTENT*01"
    assert len(unknown_j) == 1
    assert unknown_j[0]["allele_name"] == "IGKJ-MISSING*99"


# ──────────────────────────────────────────────────────────────────
# 9. VJ chain ignores D contribution with a warning
# ──────────────────────────────────────────────────────────────────


def test_vj_chain_ignores_d_call_with_warning() -> None:
    """On a VJ cartridge, AIRR records that name a D allele get
    their D contribution silently dropped AFTER a single
    stage-level warning surfaces in the report."""
    builder, v, j = _builder_vj()
    records = [
        {"v_call": v[0], "d_call": "IGHD1*01", "j_call": j[0]},
        {"v_call": v[0], "d_call": "IGHD2*01", "j_call": j[0]},
    ]
    builder.estimate_allele_usage(records, min_count=0.5)
    stage = next(
        s for s in builder.report().stages
        if s["stage"] == "estimate_allele_usage"
    )
    # D contribution is empty under VJ.
    assert stage["inferred"]["D"] == {}
    # And the per-stage warnings carry a single notice.
    assert any(
        "VJ" in w and "d_call" in w for w in stage["warnings"]
    ), f"expected VJ-D warning in {stage['warnings']!r}"
    # The warning fires once even if MANY rows have a D entry.
    d_warnings = [w for w in stage["warnings"] if "d_call" in w]
    assert len(d_warnings) == 1


# ──────────────────────────────────────────────────────────────────
# 10. VDJ chain missing D recorded
# ──────────────────────────────────────────────────────────────────


def test_vdj_chain_missing_d_call_skips_row_with_report_entry() -> None:
    """On a VDJ cartridge, AIRR records that LACK a `d_call`
    are skipped with a structured rejection. Per audit §5 they
    do NOT contribute to V/J counts either (the row is
    treated as malformed at the chain-type-validation gate)."""
    src = ga.HUMAN_IGH_OGRDB
    v_fasta, j_fasta, v_names, j_names = _fastas_and_names(src)
    d_fasta = "".join(
        f">{a.name}\n{a.ungapped_seq}\n"
        for gene_alleles in list(src.d_alleles.values())[:2]
        for a in gene_alleles[:1]
    )
    builder = ga.ReferenceCartridgeBuilder.from_fasta(
        v_fasta=v_fasta, j_fasta=j_fasta, d_fasta=d_fasta, chain_type="BCR_HEAVY"
    )
    d_names = [list(src.d_alleles.values())[i][0].name for i in range(2)]
    records = [
        {"v_call": v_names[0], "d_call": "", "j_call": j_names[0]},
        {"v_call": v_names[0], "d_call": d_names[0], "j_call": j_names[0]},
    ]
    builder.estimate_allele_usage(records, min_count=0.5)
    stage = next(
        s for s in builder.report().stages
        if s["stage"] == "estimate_allele_usage"
    )
    assert stage["inferred"]["skipped"]["missing_d_call_on_vdj"] == 1
    # The malformed row didn't contribute to V/J either.
    v_count_in_kept = 1  # only the well-formed row
    assert stage["inferred"]["V"] == {v_names[0]: 1.0}
    # Rejection entry names the reason.
    missing_d = [
        r for r in builder.report().rejected
        if r.get("reason") == "missing_d_call_on_vdj"
    ]
    assert len(missing_d) == 1


# ──────────────────────────────────────────────────────────────────
# 11. min_count drops low-support alleles
# ──────────────────────────────────────────────────────────────────


def test_min_count_drops_low_support_alleles_before_normalisation() -> None:
    """An allele whose raw count is strictly below `min_count`
    is dropped before normalisation. The kept alleles' weights
    renormalise to sum to 1.0; the dropped count surfaces in
    the report's `below_min_count` block."""
    builder, v, j = _builder_vj()
    # v[0] hit twice, v[1] hit once, v[2] hit 5×.
    records = (
        [{"v_call": v[0], "j_call": j[0]}] * 2
        + [{"v_call": v[1], "j_call": j[0]}] * 1
        + [{"v_call": v[2], "j_call": j[0]}] * 5
    )
    builder.estimate_allele_usage(records, min_count=2.0)
    stage = next(
        s for s in builder.report().stages
        if s["stage"] == "estimate_allele_usage"
    )
    v_weights = stage["inferred"]["V"]
    # v[1] (count=1) drops; v[0] (count=2) and v[2] (count=5) kept.
    assert set(v_weights.keys()) == {v[0], v[2]}
    assert abs(v_weights[v[0]] - 2 / 7) < 1e-9
    assert abs(v_weights[v[2]] - 5 / 7) < 1e-9
    # below_min_count surfaces the drop.
    assert stage["inferred"]["below_min_count"]["V"] == 1


# ──────────────────────────────────────────────────────────────────
# 12. Stage entry shape matches audit §6.2
# ──────────────────────────────────────────────────────────────────


def test_stage_entry_shape_matches_audit() -> None:
    """The stage entry carries the canonical
    `{stage, inputs, inferred, warnings}` keys per audit §6.2.
    `inputs` has `record_count` / `ambiguous` / `min_count` /
    `source` / `replaced`; `inferred` has V / D / J /
    `skipped` / `below_min_count`."""
    builder, v, j = _builder_vj()
    records = [{"v_call": v[0], "j_call": j[0]}]
    builder.estimate_allele_usage(records, min_count=0.5)
    stage = next(
        s for s in builder.report().stages
        if s["stage"] == "estimate_allele_usage"
    )
    assert set(stage.keys()) >= {"stage", "inputs", "inferred", "warnings"}
    for k in ("record_count", "ambiguous", "min_count", "source", "replaced"):
        assert k in stage["inputs"], f"missing inputs key {k!r}"
    for k in ("V", "D", "J", "skipped", "below_min_count"):
        assert k in stage["inferred"], f"missing inferred key {k!r}"
    # source label is structured.
    assert stage["inputs"]["source"].startswith("records:")


# ──────────────────────────────────────────────────────────────────
# 13. Manifest block reflects estimator output
# ──────────────────────────────────────────────────────────────────


def test_manifest_allele_usage_block_reflects_estimator_output() -> None:
    """After estimation, the built cartridge's
    `manifest['models']['allele_usage']` reflects which
    segments are populated; `available=True`,
    `nonempty_segments=["V", "J"]` (D empty on VJ chain)."""
    builder, v, j = _builder_vj()
    records = [{"v_call": v[0], "j_call": j[0]}]
    builder.estimate_allele_usage(records, min_count=0.5)
    cfg = builder.build()
    au = cfg.cartridge_manifest()["models"]["allele_usage"]
    assert au["available"] is True
    assert "V" in au["nonempty_segments"]
    assert "J" in au["nonempty_segments"]
    assert "D" not in au["nonempty_segments"]  # VJ chain, no D
    # Soft-gap pin unchanged.
    assert au["in_plan_signature"] is False


# ──────────────────────────────────────────────────────────────────
# 14. Pickle round-trip preserves estimated plane + report
# ──────────────────────────────────────────────────────────────────


def test_pickle_round_trip_preserves_plane_and_report() -> None:
    """A cartridge produced by `estimate_allele_usage` →
    `build()`, pickled and unpickled, retains both:

    - `cfg.reference_models.allele_usage` (the typed plane).
    - `cfg.build_report` carrying the estimator's stage entry.
    """
    builder, v, j = _builder_vj()
    records = [{"v_call": v[0], "j_call": j[0]}] * 3
    builder.estimate_allele_usage(records, min_count=0.5)
    cfg = builder.build()
    blob = pickle.dumps(cfg)
    cfg2 = pickle.loads(blob)
    # Typed plane survives.
    assert cfg2.reference_models.allele_usage is not None
    assert cfg2.reference_models.allele_usage.v == {v[0]: 1.0}
    # Build report survives + carries the estimator stage.
    assert cfg2.build_report is not None
    stages = [s["stage"] for s in cfg2.build_report.stages]
    assert "estimate_allele_usage" in stages


# ──────────────────────────────────────────────────────────────────
# 15. Contract absence pins flipped; soft-gap pin holds
# ──────────────────────────────────────────────────────────────────


def test_contract_pins_post_slice_state() -> None:
    """Integration probe — the contract file's absence pins
    are flipped to present-state, while the soft-gap pin
    (plan-signature does NOT change on `v_allele_weights`)
    still holds. The contract file is the source of truth;
    this test is the safety net so a regression in either
    direction surfaces in the implementation test file too."""
    # Builder method present.
    assert hasattr(ga.ReferenceCartridgeBuilder, "estimate_allele_usage")
    # Spec dataclass + top-level export present.
    assert hasattr(ga, "AlleleUsageSpec")
    # Manifest block present.
    m = ga.HUMAN_IGH_OGRDB.cartridge_manifest()
    assert "allele_usage" in m["models"]
    # Plan-signature soft gap 1 still holds — empirical reprise.
    exp_a = ga.Experiment.on("human_igh").recombine()
    exp_b = ga.Experiment.on("human_igh").recombine(
        v_allele_weights={"IGHVF1-G1*01": 100.0}
    )
    sa = json.loads(
        exp_a.compile().simulator.trace_file_from(
            exp_a.compile().simulator.run(seed=0), seed=0
        ).to_json()
    )["pass_plan_signature"]
    sb = json.loads(
        exp_b.compile().simulator.trace_file_from(
            exp_b.compile().simulator.run(seed=0), seed=0
        ).to_json()
    )["pass_plan_signature"]
    assert sa == sb, (
        "soft gap 1 closed unexpectedly — verify the cartridge-plane "
        "path is folded into the plan signature in lockstep"
    )
