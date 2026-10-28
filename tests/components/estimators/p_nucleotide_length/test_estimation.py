"""End-to-end implementation tests for the **P-Nucleotide
Length Estimation v1** slice.

Companion to
[`docs/p_nucleotide_length_estimation_design.md`](../docs/p_nucleotide_length_estimation_design.md)
and the contract pins in
[`tests/test_p_nucleotide_length_estimation_contract.py`](test_p_nucleotide_length_estimation_contract.py).

Covers the 14-test surface from the user brief:

 1. VJ estimates V_3 / J_5 only.
 2. VDJ estimates all four keys.
 3. Missing / malformed fields skipped per key
    (field-local).
 4. VJ nonzero D-end fields warn/skip.
 5. `min_count` drops low-support lengths.
 6. `pseudocount` affects weights.
 7. `replace=False` rejects existing model.
 8. Built cartridge uses estimated P lengths in
    recombination defaults.
 9. Direct fixed P fields produce expected distributions.
10. Zero-heavy input triggers ≥ 95% warning per key.
11. Generic P-naive records produce zero-heavy warnings,
    not heuristic inference.
12. Report JSON-clean and stage shape stable.
13. Contract absence pins flip (integration probe).
14. Heuristic-inference absence pin holds — method body
    contains no `junction_length` / `np1` / `np2`
    references.
"""
from __future__ import annotations

import inspect
import json
import math

import pytest

import GenAIRR as ga
from GenAIRR.reference_models import (
    EmpiricalDistributionSpec,
    ReferenceEmpiricalModels,
)


# ──────────────────────────────────────────────────────────────────
# Inline FASTA fixtures — same shape as the prior estimator
# implementation tests.
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
    """On a VJ cartridge the estimator emits only `V_3`
    and `J_5` typed-plane specs. `D_5` / `D_3` are absent.
    The stage entry carries all four keys for shape
    stability — D entries empty."""
    builder = _builder_vj()
    records = [
        {"p_v_3_length": 1, "p_j_5_length": 0},
        {"p_v_3_length": 2, "p_j_5_length": 1},
        {"p_v_3_length": 0, "p_j_5_length": 2},
    ]
    builder.estimate_p_nucleotide_lengths(records, min_count=1)

    stage = next(
        s for s in builder.report().stages
        if s["stage"] == "estimate_p_nucleotide_lengths"
    )
    inferred = stage["inferred"]
    for key in ("V_3", "D_5", "D_3", "J_5"):
        assert key in inferred
    assert inferred["V_3"] and inferred["J_5"]
    assert inferred["D_5"] == []
    assert inferred["D_3"] == []

    cfg = builder.build()
    assert set(cfg.reference_models.p_nucleotide_lengths.keys()) == {
        "V_3", "J_5"
    }


# ──────────────────────────────────────────────────────────────────
# 2. VDJ estimates all four keys
# ──────────────────────────────────────────────────────────────────


def test_vdj_cartridge_writes_all_four_p_length_keys() -> None:
    """On a VDJ cartridge the estimator produces a typed
    plane carrying all four `V_3` / `D_5` / `D_3` / `J_5`
    keys when input records cover them."""
    builder = _builder_vdj()
    records = [
        {"p_v_3_length": 1, "p_d_5_length": 0, "p_d_3_length": 2, "p_j_5_length": 1},
        {"p_v_3_length": 2, "p_d_5_length": 1, "p_d_3_length": 0, "p_j_5_length": 3},
        {"p_v_3_length": 0, "p_d_5_length": 0, "p_d_3_length": 1, "p_j_5_length": 2},
    ]
    builder.estimate_p_nucleotide_lengths(records, min_count=1)
    cfg = builder.build()
    assert set(cfg.reference_models.p_nucleotide_lengths.keys()) == {
        "V_3", "D_5", "D_3", "J_5"
    }
    for key in ("V_3", "D_5", "D_3", "J_5"):
        spec = cfg.reference_models.p_nucleotide_lengths[key]
        assert isinstance(spec, EmpiricalDistributionSpec)
        total = sum(w for _, w in spec.values)
        assert abs(total - 1.0) < 1e-9, (key, total)


# ──────────────────────────────────────────────────────────────────
# 3. Missing / malformed fields skipped per key (field-local)
# ──────────────────────────────────────────────────────────────────


def test_missing_and_malformed_fields_are_field_local() -> None:
    """A row malformed / missing / negative in ONE
    P-length column drops only that key's contribution —
    the row's other columns still feed their
    distributions."""
    builder = _builder_vdj()
    records = [
        # Three well-formed rows.
        {"p_v_3_length": 1, "p_d_5_length": 0, "p_d_3_length": 1, "p_j_5_length": 2},
        {"p_v_3_length": 1, "p_d_5_length": 0, "p_d_3_length": 1, "p_j_5_length": 2},
        {"p_v_3_length": 1, "p_d_5_length": 0, "p_d_3_length": 1, "p_j_5_length": 2},
        # p_v_3_length malformed.
        {"p_v_3_length": "bad", "p_d_5_length": 0, "p_d_3_length": 1, "p_j_5_length": 2},
        # p_d_5_length missing.
        {"p_v_3_length": 1, "p_d_3_length": 1, "p_j_5_length": 2},
        # p_d_3_length negative.
        {"p_v_3_length": 1, "p_d_5_length": 0, "p_d_3_length": -1, "p_j_5_length": 2},
    ]
    builder.estimate_p_nucleotide_lengths(records, min_count=1)
    stage = next(
        s for s in builder.report().stages
        if s["stage"] == "estimate_p_nucleotide_lengths"
    )
    sk = stage["inferred"]["skipped"]
    assert sk["malformed_length_value"]["p_v_3_length"] == 1
    assert sk["missing_required_column"]["p_d_5_length"] == 1
    assert sk["negative_length_value"]["p_d_3_length"] == 1

    cfg = builder.build()
    # J_5 across 6 contributing rows: always 2 → {2: 1.0}.
    j5 = dict(cfg.reference_models.p_nucleotide_lengths["J_5"].values)
    assert j5 == {2: 1.0}
    # V_3 across 5 contributing rows: always 1 → {1: 1.0}.
    v3 = dict(cfg.reference_models.p_nucleotide_lengths["V_3"].values)
    assert v3 == {1: 1.0}

    reasons = [r["reason"] for r in builder.report().rejected]
    assert reasons.count("malformed_length_value") == 1
    assert reasons.count("missing_required_column") == 1
    assert reasons.count("negative_length_value") == 1


# ──────────────────────────────────────────────────────────────────
# 4. VJ nonzero D-end fields warn/skip
# ──────────────────────────────────────────────────────────────────


def test_vj_nonzero_d_end_fields_warn_and_skip() -> None:
    """On a VJ cartridge, records carrying nonzero
    `p_d_5_length` or `p_d_3_length` surface as
    `dropped_columns` counters plus a single one-time
    warning per column."""
    builder = _builder_vj()
    records = (
        [{"p_v_3_length": 0, "p_j_5_length": 0,
          "p_d_5_length": 5, "p_d_3_length": 7}] * 3
        + [{"p_v_3_length": 0, "p_j_5_length": 0}]
    )
    builder.estimate_p_nucleotide_lengths(records, min_count=1)
    stage = next(
        s for s in builder.report().stages
        if s["stage"] == "estimate_p_nucleotide_lengths"
    )
    assert stage["inferred"]["dropped_columns"]["p_d_5_length"] == 3
    assert stage["inferred"]["dropped_columns"]["p_d_3_length"] == 3
    # One warning per column, regardless of row count.
    d5_warnings = [w for w in stage["warnings"] if "p_d_5_length" in w]
    d3_warnings = [w for w in stage["warnings"] if "p_d_3_length" in w]
    assert len(d5_warnings) == 1
    assert len(d3_warnings) == 1


# ──────────────────────────────────────────────────────────────────
# 5. min_count drops low-support lengths
# ──────────────────────────────────────────────────────────────────


def test_min_count_drops_low_support_p_length_values() -> None:
    """A length value whose observed count is strictly
    below `min_count` is dropped before normalisation;
    surfaces in `below_min_count` per key."""
    builder = _builder_vj()
    # p_v_3_length=0 observed 5×, p_v_3_length=4 observed 1×.
    records = (
        [{"p_v_3_length": 0, "p_j_5_length": 0}] * 5
        + [{"p_v_3_length": 4, "p_j_5_length": 0}] * 1
    )
    builder.estimate_p_nucleotide_lengths(records, min_count=2)
    stage = next(
        s for s in builder.report().stages
        if s["stage"] == "estimate_p_nucleotide_lengths"
    )
    v3 = dict(stage["inferred"]["V_3"])
    assert set(v3.keys()) == {0}
    assert math.isclose(v3[0], 1.0, rel_tol=1e-9)
    assert stage["inferred"]["below_min_count"]["V_3"] == 1


# ──────────────────────────────────────────────────────────────────
# 6. pseudocount affects weights
# ──────────────────────────────────────────────────────────────────


def test_pseudocount_applies_to_observed_values_only() -> None:
    """`pseudocount` is added to each observed length
    value's count before normalisation; unobserved values
    stay unobserved."""
    builder = _builder_vj()
    records = (
        [{"p_v_3_length": 0, "p_j_5_length": 0}] * 3
        + [{"p_v_3_length": 1, "p_j_5_length": 0}] * 1
    )
    builder.estimate_p_nucleotide_lengths(
        records, min_count=1, pseudocount=1.0
    )
    stage = next(
        s for s in builder.report().stages
        if s["stage"] == "estimate_p_nucleotide_lengths"
    )
    v3 = dict(stage["inferred"]["V_3"])
    # observed counts {0: 3, 1: 1} + pseudocount 1.0 each:
    # kept = {0: 4, 1: 2}, total = 6 → {0: 4/6, 1: 2/6}.
    assert math.isclose(v3[0], 4 / 6, rel_tol=1e-9)
    assert math.isclose(v3[1], 2 / 6, rel_tol=1e-9)
    assert 5 not in v3  # unobserved → not in support


# ──────────────────────────────────────────────────────────────────
# 7. replace=False rejects existing model
# ──────────────────────────────────────────────────────────────────


def test_replace_false_rejects_existing_p_nucleotide_lengths() -> None:
    """`replace=False` blocks re-entry when the typed
    plane's `p_nucleotide_lengths` already carries a
    spec."""
    builder = _builder_vj()
    records = [{"p_v_3_length": 0, "p_j_5_length": 0}]
    builder.estimate_p_nucleotide_lengths(records, min_count=1)
    with pytest.raises(ValueError, match="already attached"):
        builder.estimate_p_nucleotide_lengths(
            records, min_count=1, replace=False
        )
    # replace=True (default) allows overwrite.
    builder.estimate_p_nucleotide_lengths(records, min_count=1)


def test_replace_false_rejects_manually_attached_p_nucleotide_lengths() -> None:
    """`replace=False` also rejects if `with_models()`
    already attached the typed plane."""
    builder = _builder_vj()
    builder.with_models(
        ReferenceEmpiricalModels(
            p_nucleotide_lengths={
                "V_3": EmpiricalDistributionSpec([(0, 1.0)]),
            }
        )
    )
    with pytest.raises(ValueError, match="already attached"):
        builder.estimate_p_nucleotide_lengths(
            [{"p_v_3_length": 0, "p_j_5_length": 0}],
            min_count=1, replace=False,
        )


# ──────────────────────────────────────────────────────────────────
# 8. Built cartridge uses estimated P lengths in recombination defaults
# ──────────────────────────────────────────────────────────────────


def test_built_cartridge_uses_estimated_p_lengths_in_recombination_defaults() -> None:
    """The bridge resolver's typed-plane precedence picks
    up the estimator's output. Confirmed via
    `extract_recombine_defaults`: the returned
    `p_v_3_lengths` distribution matches what the
    estimator wrote."""
    from GenAIRR._dataconfig_extract import extract_recombine_defaults

    builder = _builder_vj()
    records = (
        [{"p_v_3_length": 0, "p_j_5_length": 0}] * 3
        + [{"p_v_3_length": 1, "p_j_5_length": 0}]
    )
    builder.estimate_p_nucleotide_lengths(records, min_count=1)
    cfg = builder.build()
    defaults = extract_recombine_defaults(cfg)
    p_v_3 = defaults["p_v_3_lengths"]
    assert p_v_3 is not None
    as_dict = dict(p_v_3)
    assert math.isclose(as_dict[0], 0.75, rel_tol=1e-9)
    assert math.isclose(as_dict[1], 0.25, rel_tol=1e-9)


# ──────────────────────────────────────────────────────────────────
# 9. Direct fixed P fields produce expected distributions
# ──────────────────────────────────────────────────────────────────


def test_direct_fixed_p_fields_produce_expected_distributions() -> None:
    """A controlled dataset with known P-length frequencies
    produces the expected normalised distribution per
    key — no surprises in the count → normalisation
    pipeline."""
    builder = _builder_vdj()
    # Hand-designed input: 50% length-0, 30% length-1, 20% length-2 on V_3.
    records = (
        [{"p_v_3_length": 0, "p_d_5_length": 1, "p_d_3_length": 0, "p_j_5_length": 0}] * 5
        + [{"p_v_3_length": 1, "p_d_5_length": 1, "p_d_3_length": 0, "p_j_5_length": 0}] * 3
        + [{"p_v_3_length": 2, "p_d_5_length": 1, "p_d_3_length": 0, "p_j_5_length": 0}] * 2
    )
    builder.estimate_p_nucleotide_lengths(records, min_count=1)
    stage = next(
        s for s in builder.report().stages
        if s["stage"] == "estimate_p_nucleotide_lengths"
    )
    v3 = dict(stage["inferred"]["V_3"])
    assert math.isclose(v3[0], 0.5, rel_tol=1e-9)
    assert math.isclose(v3[1], 0.3, rel_tol=1e-9)
    assert math.isclose(v3[2], 0.2, rel_tol=1e-9)
    # D_5: every row contributes 1 → {1: 1.0}.
    d5 = dict(stage["inferred"]["D_5"])
    assert d5 == {1: 1.0}


# ──────────────────────────────────────────────────────────────────
# 10. Zero-heavy input triggers ≥ 95% warning per key
# ──────────────────────────────────────────────────────────────────


def test_zero_heavy_input_triggers_per_key_provenance_warning() -> None:
    """When a key's contributing rows reported zero in
    ≥ 95% of cases, a stage-level warning surfaces
    naming the key. The threshold is exclusive of the
    fraction itself (≥ 0.95 fires); the message format is
    the canonical `p_nucleotide_lengths[<key>] is >=95%
    zero; ...` shape."""
    builder = _builder_vdj()
    # V_3: 20 rows at 0, 0 nonzero → 100% zero, warning fires.
    # D_5: 19 zero + 1 nonzero → 95% zero, warning fires.
    # D_3: 18 zero + 2 nonzero → 90% zero, NO warning.
    # J_5: 19 zero + 1 nonzero → 95% zero, warning fires.
    records = []
    for _ in range(20):
        records.append({
            "p_v_3_length": 0,
            "p_d_5_length": 0,
            "p_d_3_length": 0,
            "p_j_5_length": 0,
        })
    # Patch one entry of D_5 / J_5 to nonzero (95% zero ratio).
    records[0]["p_d_5_length"] = 1
    records[0]["p_j_5_length"] = 1
    # Patch two D_3 entries to nonzero (90% zero ratio — no warning).
    records[0]["p_d_3_length"] = 1
    records[1]["p_d_3_length"] = 1

    builder.estimate_p_nucleotide_lengths(records, min_count=1)
    stage = next(
        s for s in builder.report().stages
        if s["stage"] == "estimate_p_nucleotide_lengths"
    )
    warnings = stage["warnings"]
    # V_3 + D_5 + J_5 each get a >=95% warning; D_3 does not.
    for key in ("V_3", "D_5", "J_5"):
        matching = [w for w in warnings
                    if f"p_nucleotide_lengths[{key}]" in w
                    and ">=95% zero" in w]
        assert len(matching) == 1, (
            f"expected exactly one >=95% warning for {key}, got "
            f"{matching!r} (all warnings: {warnings!r})"
        )
    d3_matching = [w for w in warnings
                   if "p_nucleotide_lengths[D_3]" in w
                   and ">=95% zero" in w]
    assert len(d3_matching) == 0, (
        f"D_3 zero fraction is 90% — should NOT trigger the >=95% "
        f"warning, got {d3_matching!r}"
    )
    # zero_fraction matches expectation.
    zf = stage["inferred"]["zero_fraction"]
    assert math.isclose(zf["V_3"], 1.0, rel_tol=1e-9)
    assert math.isclose(zf["D_5"], 0.95, rel_tol=1e-9)
    assert math.isclose(zf["D_3"], 0.90, rel_tol=1e-9)
    assert math.isclose(zf["J_5"], 0.95, rel_tol=1e-9)


# ──────────────────────────────────────────────────────────────────
# 11. Generic P-naive records produce zero-heavy warnings
# ──────────────────────────────────────────────────────────────────


def test_generic_p_naive_records_produce_zero_heavy_warnings_not_heuristic_inference() -> None:
    """Records lacking P-length columns OR populating them
    as zero (the external-AIRR-tool case) produce a
    rejection storm + zero-heavy warnings, NOT a heuristic
    inference from junction / NP sibling columns. The
    estimator never infers P-lengths from non-P fields."""
    builder = _builder_vdj()
    # Generic AIRR-like records — junction populated, NP fields
    # populated, but ZERO P-length values. This simulates the
    # external-tool case where P columns are written as 0 per
    # AIRR-C schema's optional-field default.
    records = [
        {
            "junction": "TGCAAGAGGCAGTCAGGGGAGGTGACTACTAC",
            "junction_length": 32,
            "np1": "AAGAG",
            "np2": "GCAG",
            "np1_length": 5,
            "np2_length": 4,
            "p_v_3_length": 0,
            "p_d_5_length": 0,
            "p_d_3_length": 0,
            "p_j_5_length": 0,
        }
        for _ in range(50)
    ]
    builder.estimate_p_nucleotide_lengths(records, min_count=1)
    cfg = builder.build()
    # All four keys: estimator wrote {0: 1.0} (degenerate),
    # NOT a derived distribution from junction / NP arithmetic.
    for key in ("V_3", "D_5", "D_3", "J_5"):
        spec = cfg.reference_models.p_nucleotide_lengths[key]
        pairs = dict(spec.values)
        assert pairs == {0: 1.0}, (
            f"{key} distribution diverged from `[(0, 1.0)]` — the "
            f"estimator inferred P-lengths from a sibling column: "
            f"{pairs!r}"
        )
    # All four keys emit the >= 95% zero warning.
    stage = next(
        s for s in builder.report().stages
        if s["stage"] == "estimate_p_nucleotide_lengths"
    )
    for key in ("V_3", "D_5", "D_3", "J_5"):
        matching = [w for w in stage["warnings"]
                    if f"p_nucleotide_lengths[{key}]" in w
                    and ">=95% zero" in w]
        assert len(matching) == 1, (
            f"P-naive input: {key} should trigger zero-heavy "
            f"warning, got {matching!r}"
        )


# ──────────────────────────────────────────────────────────────────
# 12. Report JSON-clean and stage shape stable
# ──────────────────────────────────────────────────────────────────


def test_stage_entry_shape_matches_brief_and_report_is_json_clean() -> None:
    """Canonical `{stage, inputs, inferred, warnings}`
    shape. `inputs` has `record_count` / `min_count` /
    `pseudocount` / `source` / `replaced`. `inferred` has
    `V_3` / `D_5` / `D_3` / `J_5` / `skipped` /
    `below_min_count` / `zero_fraction` /
    `dropped_columns`. Full report round-trips through
    `json.dumps`."""
    builder = _builder_vdj()
    records = [{"p_v_3_length": 1, "p_d_5_length": 0,
                "p_d_3_length": 1, "p_j_5_length": 2}]
    builder.estimate_p_nucleotide_lengths(records, min_count=1)
    stage = next(
        s for s in builder.report().stages
        if s["stage"] == "estimate_p_nucleotide_lengths"
    )
    assert set(stage.keys()) >= {"stage", "inputs", "inferred", "warnings"}
    for k in ("record_count", "min_count", "pseudocount", "source", "replaced"):
        assert k in stage["inputs"], f"missing inputs key {k!r}"
    for k in ("V_3", "D_5", "D_3", "J_5", "skipped",
              "below_min_count", "zero_fraction", "dropped_columns"):
        assert k in stage["inferred"], f"missing inferred key {k!r}"
    blob = json.dumps(builder.report().to_dict())
    again = json.loads(blob)
    assert isinstance(again, dict)
    assert "stages" in again


# ──────────────────────────────────────────────────────────────────
# 13. Contract absence pins flip (integration probe)
# ──────────────────────────────────────────────────────────────────


def test_contract_pins_post_slice_state() -> None:
    """Integration probe — the contract file's absence pins
    flip; manifest's existing `p_nucleotide_models` block
    reflects the estimator's output without any extension."""
    assert hasattr(
        ga.ReferenceCartridgeBuilder, "estimate_p_nucleotide_lengths"
    )
    m = ga.HUMAN_IGH_OGRDB.cartridge_manifest()
    assert "p_nucleotide_models" in m["models"]


# ──────────────────────────────────────────────────────────────────
# 14. Heuristic-inference absence pin holds
# ──────────────────────────────────────────────────────────────────


def test_estimator_method_body_does_not_reference_junction_or_np_fields() -> None:
    """The estimator method body MUST NOT reference
    `junction_length` / `np1` / `np2` — v1 explicitly
    bans junction-arithmetic / NP-string heuristic
    P-length inference. Docstring mentions are allowed
    (they document the forbidden columns); only code
    references would fail this pin."""
    method = ga.ReferenceCartridgeBuilder.estimate_p_nucleotide_lengths
    src = inspect.getsource(method)
    lines = src.splitlines()
    in_docstring = False
    body_lines: list[str] = []
    for line in lines:
        stripped = line.strip()
        if stripped.startswith('"""') or stripped.startswith("'''"):
            triple = '"""' if '"""' in stripped else "'''"
            if stripped.count(triple) == 2:
                continue  # single-line docstring, skip
            in_docstring = not in_docstring
            continue
        if in_docstring:
            continue
        body_lines.append(line)
    body = "\n".join(body_lines)
    for forbidden in ("junction_length", '"np1"', '"np2"',
                      "'np1'", "'np2'", "np1_length", "np2_length"):
        assert forbidden not in body, (
            f"estimate_p_nucleotide_lengths body references "
            f"{forbidden!r} — v1 forbids junction-arithmetic / "
            f"NP-string heuristic P-length inference"
        )
