"""End-to-end implementation tests for the **Trim Distribution
Estimation v1** slice.

Companion to
[`docs/trim_distribution_estimation_design.md`](../docs/trim_distribution_estimation_design.md)
and the contract pins in
[`tests/test_trim_distribution_estimation_contract.py`](test_trim_distribution_estimation_contract.py).

Covers the 12-test surface from the user brief:

 1. VJ estimates `V_3` / `J_5` only.
 2. VDJ estimates all four keys.
 3. Missing / malformed fields recorded but other fields
    still counted (field-local validation).
 4. `min_count` drops low-support values.
 5. `pseudocount` affects weights (additive over observed
    values only).
 6. `replace=False` rejects existing trims.
 7. Built cartridge uses estimated trims in recombination
    defaults.
 8. Explicit `.trim(...)` override still wins.
 9. End-loss fields ignored.
10. Report JSON-clean and stage shape stable.
11. Manifest `models.trim_models` block reflects estimated
    keys.
12. Contract absence pins flip (integration probe).
"""
from __future__ import annotations

import json
import math

import pytest

import GenAIRR as ga
from GenAIRR.reference_models import (
    EmpiricalDistributionSpec,
    ReferenceEmpiricalModels,
)


# ──────────────────────────────────────────────────────────────────
# Inline FASTA fixtures — small enough to keep tests fast.
# Synthetic anchors are absent; recombination compiles under
# `allow_curatable_refdata()` for the few tests that exercise it.
# ──────────────────────────────────────────────────────────────────

_VJ_V_FASTA = (
    ">IGKV1-MOCK*01\n"
    "GACATCCAGATGACCCAGTCTCCATCCTCCCTGTCTGCATCTGTAGGAGACAGAGTCACC"
    "ATCACTTGCCGGGCAAGTCAGAGCATTAGCAGCTATTTAAATTGGTATCAGCAGAAACCA\n"
    ">IGKV2-MOCK*01\n"
    "GAAATTGTGTTGACGCAGTCTCCAGGCACCCTGTCTTTGTCTCCAGGGGAAAGAGCCACC"
    "CTCTCCTGCAGGGCCAGTCAGAGTGTTAGCAGCAGCTACTTAGCCTGGTACCAGCAGAAA\n"
)
_VJ_J_FASTA = (
    ">IGKJ1-MOCK*01\n"
    "GTGGACGTTCGGCCAAGGGACCAAGGTGGAAATCAAAC\n"
    ">IGKJ2-MOCK*01\n"
    "TGTACACTTTTGGCCAGGGGACCAAGCTGGAGATCAAAC\n"
)

_VDJ_V_FASTA = (
    ">IGHV1-MOCK*01\n"
    "GAGGTGCAGCTGGTGGAGTCTGGGGGAGGCTTGGTACAGCCTGGGGGGTCCCTGAGACTC"
    "TCCTGTGCAGCCTCTGGATTCACCTTCAGTAGCTATGCCATGAGCTGGGTCCGCCAGGCT\n"
    ">IGHV2-MOCK*01\n"
    "CAGGTCAACTTAAGGGAGTCTGGTCCTGCGCTGGTGAAACCCACACAGACCCTCACACTG"
    "ACCTGCACCTTCTCTGGGTTCTCACTCAGCACTAGTGGAATGTGTGTGAGCTGGATCCGT\n"
)
_VDJ_D_FASTA = (
    ">IGHD1-MOCK*01\n"
    "GGGTATAGCAGCAGCTGGTAC\n"
    ">IGHD2-MOCK*01\n"
    "AGGATATTGTAGTGGTGGTAGCTGCTACTCC\n"
)
_VDJ_J_FASTA = (
    ">IGHJ1-MOCK*01\n"
    "TACTACTACGGTATGGACGTCTGGGGCCAAGGGACCACGGTCACCGTCTCCTCAG\n"
    ">IGHJ2-MOCK*01\n"
    "ACTACTGGTACTTCGATCTCTGGGGCCGTGGCACCCTGGTCACTGTCTCCTCAG\n"
)


def _builder_vj() -> "ga.ReferenceCartridgeBuilder":
    return ga.ReferenceCartridgeBuilder.from_fasta(
        v_fasta=_VJ_V_FASTA, j_fasta=_VJ_J_FASTA,
        chain_type="BCR_LIGHT_KAPPA",
    )


def _builder_vdj() -> "ga.ReferenceCartridgeBuilder":
    return ga.ReferenceCartridgeBuilder.from_fasta(
        v_fasta=_VDJ_V_FASTA, d_fasta=_VDJ_D_FASTA, j_fasta=_VDJ_J_FASTA,
        chain_type="BCR_HEAVY",
    )


# ──────────────────────────────────────────────────────────────────
# 1. VJ estimates V_3 / J_5 only
# ──────────────────────────────────────────────────────────────────


def test_vj_cartridge_writes_only_v3_and_j5_keys() -> None:
    """On a VJ cartridge, the estimator emits only `V_3` /
    `J_5` typed-plane specs. `D_5` / `D_3` are absent (no D
    segment to trim). The stage entry still carries the
    four-key `inferred` block for shape stability, but the
    D entries are empty lists."""
    builder = _builder_vj()
    records = [
        {"v_trim_3": 1, "j_trim_5": 0},
        {"v_trim_3": 2, "j_trim_5": 1},
        {"v_trim_3": 0, "j_trim_5": 2},
    ]
    builder.estimate_trim_distributions(records, min_count=1)

    # Stage entry shape — V_3 / D_5 / D_3 / J_5 all present
    # as keys; D entries are empty.
    stage = next(
        s for s in builder.report().stages
        if s["stage"] == "estimate_trim_distributions"
    )
    inferred = stage["inferred"]
    for key in ("V_3", "D_5", "D_3", "J_5"):
        assert key in inferred
    assert inferred["V_3"] and inferred["J_5"]
    assert inferred["D_5"] == []
    assert inferred["D_3"] == []

    # Built cartridge: only V_3 / J_5 in the typed plane.
    cfg = builder.build()
    assert set(cfg.reference_models.trims.keys()) == {"V_3", "J_5"}


# ──────────────────────────────────────────────────────────────────
# 2. VDJ estimates all four keys
# ──────────────────────────────────────────────────────────────────


def test_vdj_cartridge_writes_all_four_trim_keys() -> None:
    """On a VDJ cartridge, the estimator produces a typed
    plane carrying all four `V_3` / `D_5` / `D_3` / `J_5`
    keys when the input records cover them."""
    builder = _builder_vdj()
    records = [
        {"v_trim_3": 1, "d_trim_5": 0, "d_trim_3": 2, "j_trim_5": 1},
        {"v_trim_3": 2, "d_trim_5": 1, "d_trim_3": 0, "j_trim_5": 3},
        {"v_trim_3": 0, "d_trim_5": 0, "d_trim_3": 1, "j_trim_5": 2},
    ]
    builder.estimate_trim_distributions(records, min_count=1)
    cfg = builder.build()
    assert set(cfg.reference_models.trims.keys()) == {"V_3", "D_5", "D_3", "J_5"}
    for key in ("V_3", "D_5", "D_3", "J_5"):
        spec = cfg.reference_models.trims[key]
        assert isinstance(spec, EmpiricalDistributionSpec)
        # Weights normalised to sum 1.0.
        total = sum(w for _, w in spec.values)
        assert abs(total - 1.0) < 1e-9, (key, total)


# ──────────────────────────────────────────────────────────────────
# 3. Missing / malformed fields recorded; other fields counted
# ──────────────────────────────────────────────────────────────────


def test_missing_and_malformed_fields_are_field_local() -> None:
    """A row missing OR malformed in ONE trim column drops
    that column's contribution only — the row's other
    well-formed columns still feed their respective
    distributions."""
    builder = _builder_vdj()
    records = [
        # Three well-formed rows.
        {"v_trim_3": 1, "d_trim_5": 0, "d_trim_3": 1, "j_trim_5": 2},
        {"v_trim_3": 1, "d_trim_5": 0, "d_trim_3": 1, "j_trim_5": 2},
        {"v_trim_3": 1, "d_trim_5": 0, "d_trim_3": 1, "j_trim_5": 2},
        # v_trim_3 malformed — only that field's contribution drops.
        {"v_trim_3": "bad", "d_trim_5": 0, "d_trim_3": 1, "j_trim_5": 2},
        # d_trim_5 missing — only that field's contribution drops.
        {"v_trim_3": 1, "d_trim_3": 1, "j_trim_5": 2},
        # d_trim_3 negative — only that field's contribution drops.
        {"v_trim_3": 1, "d_trim_5": 0, "d_trim_3": -1, "j_trim_5": 2},
    ]
    builder.estimate_trim_distributions(records, min_count=1)
    stage = next(
        s for s in builder.report().stages
        if s["stage"] == "estimate_trim_distributions"
    )
    sk = stage["inferred"]["skipped"]
    assert sk["malformed_trim_value"]["v_trim_3"] == 1
    assert sk["missing_required_column"]["d_trim_5"] == 1
    assert sk["negative_trim_value"]["d_trim_3"] == 1

    # Field-local: the 6 records contributed
    # V_3: 5 / 6 (the malformed row dropped only v_trim_3)
    # D_5: 5 / 6 (the missing-d_5 row dropped only d_trim_5)
    # D_3: 5 / 6 (the negative-d_3 row dropped only d_trim_3)
    # J_5: 6 / 6 (every row has a valid j_trim_5)
    cfg = builder.build()
    j5_values = dict(cfg.reference_models.trims["J_5"].values)
    # All six records contribute j_trim_5=2 → distribution is {2: 1.0}.
    assert j5_values == {2: 1.0}

    # V_3 distribution: five rows all v_trim_3=1 → {1: 1.0}.
    v3_values = dict(cfg.reference_models.trims["V_3"].values)
    assert v3_values == {1: 1.0}

    # Structured rejection entries land per (row, column).
    reasons = [r["reason"] for r in builder.report().rejected]
    assert reasons.count("malformed_trim_value") == 1
    assert reasons.count("missing_required_column") == 1
    assert reasons.count("negative_trim_value") == 1


# ──────────────────────────────────────────────────────────────────
# 4. min_count drops low-support trim values
# ──────────────────────────────────────────────────────────────────


def test_min_count_drops_low_support_trim_values() -> None:
    """A trim value whose observed count is strictly below
    `min_count` is dropped before normalisation; surfaces in
    `below_min_count` per key."""
    builder = _builder_vj()
    # v_trim_3=0 observed 5×, v_trim_3=1 observed 1×.
    records = (
        [{"v_trim_3": 0, "j_trim_5": 0}] * 5
        + [{"v_trim_3": 1, "j_trim_5": 0}] * 1
    )
    builder.estimate_trim_distributions(records, min_count=2)
    stage = next(
        s for s in builder.report().stages
        if s["stage"] == "estimate_trim_distributions"
    )
    # v_trim_3=1 (count=1) drops; v_trim_3=0 (count=5) kept.
    v3 = dict(stage["inferred"]["V_3"])
    assert set(v3.keys()) == {0}
    assert abs(v3[0] - 1.0) < 1e-9
    assert stage["inferred"]["below_min_count"]["V_3"] == 1


# ──────────────────────────────────────────────────────────────────
# 5. pseudocount affects weights
# ──────────────────────────────────────────────────────────────────


def test_pseudocount_applies_to_observed_values_only() -> None:
    """`pseudocount` is added to every observed value's count
    before normalisation. v1 does NOT expand the support
    (unobserved values remain unobserved). With raw counts
    {0: 3, 1: 1} and pseudocount=1.0:

    - kept = {0: 4, 1: 2} (3+1 / 1+1)
    - normalised = {0: 4/6, 1: 2/6}
    """
    builder = _builder_vj()
    records = (
        [{"v_trim_3": 0, "j_trim_5": 0}] * 3
        + [{"v_trim_3": 1, "j_trim_5": 0}] * 1
    )
    builder.estimate_trim_distributions(records, min_count=1, pseudocount=1.0)
    stage = next(
        s for s in builder.report().stages
        if s["stage"] == "estimate_trim_distributions"
    )
    v3 = dict(stage["inferred"]["V_3"])
    assert math.isclose(v3[0], 4 / 6, rel_tol=1e-9)
    assert math.isclose(v3[1], 2 / 6, rel_tol=1e-9)
    # Unobserved values (e.g. v_trim_3=2) stay absent.
    assert 2 not in v3


# ──────────────────────────────────────────────────────────────────
# 6. replace=False rejects existing trims
# ──────────────────────────────────────────────────────────────────


def test_replace_false_rejects_existing_trims() -> None:
    """With `replace=False`, a second `estimate_trim_distributions`
    call AND a manual `with_models(...)` attach of a typed
    plane both block re-entry."""
    builder = _builder_vj()
    records = [{"v_trim_3": 0, "j_trim_5": 0}]
    builder.estimate_trim_distributions(records, min_count=1)
    # Second call with replace=False should reject.
    with pytest.raises(ValueError, match="already attached"):
        builder.estimate_trim_distributions(records, min_count=1, replace=False)
    # First call with replace=True (default) overwrites.
    builder.estimate_trim_distributions(records, min_count=1)


def test_replace_false_rejects_manually_attached_trims() -> None:
    """`replace=False` also rejects if `with_models()` already
    attached an explicit typed plane carrying `trims`."""
    builder = _builder_vj()
    builder.with_models(
        ReferenceEmpiricalModels(
            trims={"V_3": EmpiricalDistributionSpec([(0, 1.0)])}
        )
    )
    records = [{"v_trim_3": 0, "j_trim_5": 0}]
    with pytest.raises(ValueError, match="already attached"):
        builder.estimate_trim_distributions(records, min_count=1, replace=False)


# ──────────────────────────────────────────────────────────────────
# 7. Built cartridge uses estimated trims in recombination defaults
# ──────────────────────────────────────────────────────────────────


def test_built_cartridge_uses_estimated_trims_in_recombination_defaults() -> None:
    """The bridge resolver's typed-plane > legacy-dict
    precedence picks up the estimator's output. Confirmed
    via `extract_recombine_defaults`: the returned
    `trim_v_3` distribution matches what the estimator
    wrote."""
    from GenAIRR._dataconfig_extract import extract_recombine_defaults

    builder = _builder_vj()
    records = [{"v_trim_3": 0, "j_trim_5": 0}] * 3 + [
        {"v_trim_3": 1, "j_trim_5": 0}
    ]
    builder.estimate_trim_distributions(records, min_count=1)
    cfg = builder.build()
    defaults = extract_recombine_defaults(cfg)
    # The typed plane's V_3 distribution: {0: 0.75, 1: 0.25}.
    trim_v_3 = defaults["trim_v_3"]
    assert trim_v_3 is not None
    as_dict = dict(trim_v_3)
    assert math.isclose(as_dict[0], 0.75, rel_tol=1e-9)
    assert math.isclose(as_dict[1], 0.25, rel_tol=1e-9)


# ──────────────────────────────────────────────────────────────────
# 8. Explicit .trim(...) override still wins
# ──────────────────────────────────────────────────────────────────


def test_explicit_trim_dsl_override_takes_precedence_over_typed_plane() -> None:
    """`Experiment.trim(v_3=...)` chained after `recombine()`
    overrides the cartridge plane. Verified via the
    compiled pipeline IR's `trim_v_3` field — same surface
    the engine consumes."""
    builder = _builder_vj()
    # Cartridge plane: v_3 always 0.
    builder.estimate_trim_distributions(
        [{"v_trim_3": 0, "j_trim_5": 0}] * 4, min_count=1
    )
    cfg = builder.build()
    # Build an Experiment + .trim() override: v_3 always 5.
    exp = (
        ga.Experiment.on(cfg)
        .allow_curatable_refdata()
        .recombine()
        .trim(v_3=[(5, 1.0)])
    )
    # Inspect the IR — last step is the recombine step (override
    # mutated it via `_replace` in `.trim()`).
    from GenAIRR._pipeline_ir import _RecombineStep
    recombine_step = next(
        s for s in exp._steps if isinstance(s, _RecombineStep)
    )
    assert recombine_step.trim_v_3 is not None
    # The override won — distribution carries (5, 1.0), NOT (0, 1.0).
    pairs = dict(recombine_step.trim_v_3)
    assert pairs == {5: 1.0}, (
        f"explicit .trim(v_3=...) override didn't beat the typed "
        f"plane; pipeline IR carries {pairs!r}"
    )


# ──────────────────────────────────────────────────────────────────
# 9. End-loss fields ignored
# ──────────────────────────────────────────────────────────────────


def test_end_loss_columns_are_ignored_by_the_estimator() -> None:
    """`end_loss_5_length` / `end_loss_3_length` are
    observation-stage corruption fields — completely
    separate from recombination trims. The estimator MUST
    NOT consume them; their presence in input records is a
    silent no-op (no error, no contribution, no warning)."""
    builder = _builder_vj()
    records = [
        {"v_trim_3": 0, "j_trim_5": 0,
         "end_loss_5_length": 7, "end_loss_3_length": 11},
        {"v_trim_3": 1, "j_trim_5": 2,
         "end_loss_5_length": 99, "end_loss_3_length": 99},
    ]
    builder.estimate_trim_distributions(records, min_count=1)
    cfg = builder.build()
    v3 = dict(cfg.reference_models.trims["V_3"].values)
    j5 = dict(cfg.reference_models.trims["J_5"].values)
    # Confirmed: estimator wrote the trim columns, not end_loss.
    assert v3 == {0: 0.5, 1: 0.5}
    assert j5 == {0: 0.5, 2: 0.5}
    # Stage's warnings + dropped_columns do NOT mention end_loss.
    stage = next(
        s for s in builder.report().stages
        if s["stage"] == "estimate_trim_distributions"
    )
    for w in stage["warnings"]:
        assert "end_loss" not in w
    assert "end_loss_5_length" not in stage["inferred"]["dropped_columns"]
    assert "end_loss_3_length" not in stage["inferred"]["dropped_columns"]


def test_v_trim_5_and_j_trim_3_columns_drop_with_one_warning_each() -> None:
    """Non-zero `v_trim_5` / `j_trim_3` columns in input
    surface as `dropped_columns` counters AND a single
    one-time warning per column. The columns are tracked
    because they're recombination-trim-shaped (unlike
    end-loss, which is silently ignored)."""
    builder = _builder_vj()
    records = (
        [{"v_trim_3": 0, "j_trim_5": 0, "v_trim_5": 5, "j_trim_3": 7}] * 3
        + [{"v_trim_3": 0, "j_trim_5": 0}]  # no ignored columns
    )
    builder.estimate_trim_distributions(records, min_count=1)
    stage = next(
        s for s in builder.report().stages
        if s["stage"] == "estimate_trim_distributions"
    )
    assert stage["inferred"]["dropped_columns"]["v_trim_5"] == 3
    assert stage["inferred"]["dropped_columns"]["j_trim_3"] == 3
    # One warning per column, regardless of how many rows.
    assert sum(1 for w in stage["warnings"] if "v_trim_5" in w) == 1
    assert sum(1 for w in stage["warnings"] if "j_trim_3" in w) == 1


# ──────────────────────────────────────────────────────────────────
# 10. Report JSON-clean and stage shape stable
# ──────────────────────────────────────────────────────────────────


def test_stage_entry_shape_matches_audit_and_report_is_json_clean() -> None:
    """The stage entry carries the canonical `{stage,
    inputs, inferred, warnings}` keys per audit §6.2.
    `inputs` has `record_count` / `min_count` /
    `pseudocount` / `source` / `replaced`. `inferred` has
    `V_3` / `D_5` / `D_3` / `J_5` / `skipped` /
    `below_min_count` / `dropped_columns`. The whole
    report round-trips through `json.dumps`."""
    builder = _builder_vdj()
    records = [{"v_trim_3": 0, "d_trim_5": 0, "d_trim_3": 0, "j_trim_5": 0}]
    builder.estimate_trim_distributions(records, min_count=1)
    stage = next(
        s for s in builder.report().stages
        if s["stage"] == "estimate_trim_distributions"
    )
    assert set(stage.keys()) >= {"stage", "inputs", "inferred", "warnings"}
    for k in ("record_count", "min_count", "pseudocount", "source", "replaced"):
        assert k in stage["inputs"], f"missing inputs key {k!r}"
    for k in ("V_3", "D_5", "D_3", "J_5", "skipped",
              "below_min_count", "dropped_columns"):
        assert k in stage["inferred"], f"missing inferred key {k!r}"
    # JSON-clean round-trip.
    blob = json.dumps(builder.report().to_dict())
    again = json.loads(blob)
    assert isinstance(again, dict)
    assert "stages" in again


# ──────────────────────────────────────────────────────────────────
# 11. Manifest block reflects estimated keys
# ──────────────────────────────────────────────────────────────────


def test_manifest_trim_models_block_reflects_estimator_output() -> None:
    """After estimation, the built cartridge's
    `manifest['models']['trim_models']` reflects the typed
    plane's keys; `in_plan_signature=True` (no soft gap
    inherited)."""
    builder = _builder_vdj()
    records = [{"v_trim_3": 0, "d_trim_5": 0, "d_trim_3": 0, "j_trim_5": 0}]
    builder.estimate_trim_distributions(records, min_count=1)
    cfg = builder.build()
    block = cfg.cartridge_manifest()["models"]["trim_models"]
    assert set(block["keys"]) == {"V_3", "D_5", "D_3", "J_5"}
    assert block["source"] == "ReferenceEmpiricalModels.trims"
    assert block["in_plan_signature"] is True  # no soft gap
    assert block["legacy_fallback"] is False
    # Built cartridge has no legacy `trim_dicts`.
    assert block["legacy_trim_dicts_present"] is False


# ──────────────────────────────────────────────────────────────────
# 12. Contract absence pins flip (integration probe)
# ──────────────────────────────────────────────────────────────────


def test_contract_pins_post_slice_state() -> None:
    """Integration probe — the contract file's absence pins
    are flipped to present-state, while the soft-gap pin
    (trim distributions DO fold into plan signature — NO
    inherited soft gap) still holds. The contract file is
    the source of truth; this test is the safety net so a
    regression in either direction surfaces in the
    implementation test file too."""
    # Builder method present.
    assert hasattr(ga.ReferenceCartridgeBuilder, "estimate_trim_distributions")
    # Manifest block present.
    m = ga.HUMAN_IGH_OGRDB.cartridge_manifest()
    assert "trim_models" in m["models"]
    assert m["models"]["trim_models"]["in_plan_signature"] is True
