"""End-to-end implementation tests for the **NP Length
Distribution Estimation v1** slice.

Companion to
[`docs/np_length_estimation_design.md`](../docs/np_length_estimation_design.md)
and the contract pins in
[`tests/test_np_length_estimation_contract.py`](test_np_length_estimation_contract.py).

Covers the 12-test surface from the user brief:

 1. VJ estimates NP1 only and warns/skips NP2.
 2. VDJ estimates both NP1 and NP2.
 3. Missing / malformed fields skipped per key, not per row.
 4. `min_count` drops low-support lengths.
 5. `pseudocount` affects weights.
 6. `replace=False` rejects existing model.
 7. Built cartridge uses estimated NP lengths in
    recombination defaults.
 8. Explicit `recombine(np1_lengths=...)` override still
    wins.
 9. P-nucleotide / sequence fields ignored.
10. Report JSON-clean and stage shape stable.
11. Manifest `models.np_length_models` block reflects
    estimated keys.
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
# Inline FASTA fixtures — same shape as the trim implementation
# tests, kept small enough to keep the file fast.
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
# 1. VJ estimates NP1 only and warns/skips NP2
# ──────────────────────────────────────────────────────────────────


def test_vj_cartridge_writes_only_np1_key_and_warns_on_np2() -> None:
    """On a VJ cartridge the estimator emits only `NP1`.
    `NP2` is absent from the plane regardless of whether
    the input records carry `np2_length`. Non-zero
    `np2_length` in source rows surfaces as
    `dropped_columns["np2_length"]` with a single warning
    across the dataset."""
    builder = _builder_vj()
    records = [
        {"np1_length": 1, "np2_length": 5},  # non-zero NP2 → dropped, warned once
        {"np1_length": 2, "np2_length": 7},  # non-zero NP2 → dropped (warning re-used)
        {"np1_length": 1, "np2_length": 0},  # zero NP2 → not flagged
        {"np1_length": 3},                   # NP2 absent → not flagged
    ]
    builder.estimate_np_length_distributions(records, min_count=1)

    stage = next(
        s for s in builder.report().stages
        if s["stage"] == "estimate_np_length_distributions"
    )
    inferred = stage["inferred"]

    # Stage carries the canonical NP1 / NP2 keys for shape stability.
    assert "NP1" in inferred and "NP2" in inferred
    assert inferred["NP1"], "NP1 distribution must be populated"
    assert inferred["NP2"] == [], "NP2 must be empty on a VJ cartridge"

    # Non-zero `np2_length` flagged as dropped column + one warning.
    assert stage["inferred"]["dropped_columns"]["np2_length"] == 2
    np2_warnings = [w for w in stage["warnings"] if "np2_length" in w]
    assert len(np2_warnings) == 1, (
        f"expected exactly one np2_length warning across the "
        f"dataset, got {len(np2_warnings)}: {np2_warnings!r}"
    )

    # Built cartridge: only NP1 in the typed plane.
    cfg = builder.build()
    assert set(cfg.reference_models.np_lengths.keys()) == {"NP1"}


# ──────────────────────────────────────────────────────────────────
# 2. VDJ estimates both NP1 and NP2
# ──────────────────────────────────────────────────────────────────


def test_vdj_cartridge_writes_both_np1_and_np2_keys() -> None:
    """On a VDJ cartridge the estimator produces a typed
    plane carrying both `NP1` and `NP2` keys when the input
    records cover them. Each per-key spec normalises to
    sum 1.0."""
    builder = _builder_vdj()
    records = [
        {"np1_length": 1, "np2_length": 5},
        {"np1_length": 2, "np2_length": 5},
        {"np1_length": 3, "np2_length": 4},
    ]
    builder.estimate_np_length_distributions(records, min_count=1)
    cfg = builder.build()
    assert set(cfg.reference_models.np_lengths.keys()) == {"NP1", "NP2"}
    for key in ("NP1", "NP2"):
        spec = cfg.reference_models.np_lengths[key]
        assert isinstance(spec, EmpiricalDistributionSpec)
        total = sum(w for _, w in spec.values)
        assert abs(total - 1.0) < 1e-9, (key, total)


# ──────────────────────────────────────────────────────────────────
# 3. Missing / malformed fields skipped per key, not per row
# ──────────────────────────────────────────────────────────────────


def test_missing_and_malformed_fields_are_field_local() -> None:
    """A row malformed OR missing in ONE NP column drops
    that key's contribution only — the row's other
    well-formed column still feeds its distribution."""
    builder = _builder_vdj()
    records = [
        # Three well-formed rows.
        {"np1_length": 1, "np2_length": 5},
        {"np1_length": 1, "np2_length": 5},
        {"np1_length": 1, "np2_length": 5},
        # np1_length malformed — drops only NP1; NP2 still counts.
        {"np1_length": "bad", "np2_length": 5},
        # np2_length missing — drops only NP2; NP1 still counts.
        {"np1_length": 1},
        # np1_length negative — drops only NP1; NP2 still counts.
        {"np1_length": -1, "np2_length": 5},
    ]
    builder.estimate_np_length_distributions(records, min_count=1)
    stage = next(
        s for s in builder.report().stages
        if s["stage"] == "estimate_np_length_distributions"
    )
    sk = stage["inferred"]["skipped"]
    assert sk["malformed_length_value"]["np1_length"] == 1
    assert sk["missing_required_column"]["np2_length"] == 1
    assert sk["negative_length_value"]["np1_length"] == 1

    # Field-local: NP1 contributed 4/6 (3 well-formed + the
    # missing-np2 row); NP2 contributed 5/6 (3 well-formed +
    # malformed-np1 row + negative-np1 row).
    cfg = builder.build()
    np1 = dict(cfg.reference_models.np_lengths["NP1"].values)
    np2 = dict(cfg.reference_models.np_lengths["NP2"].values)
    # All four NP1-contributing rows had np1_length=1.
    assert np1 == {1: 1.0}
    # All five NP2-contributing rows had np2_length=5.
    assert np2 == {5: 1.0}

    # Structured rejection entries.
    reasons = [r["reason"] for r in builder.report().rejected]
    assert reasons.count("malformed_length_value") == 1
    assert reasons.count("missing_required_column") == 1
    assert reasons.count("negative_length_value") == 1


# ──────────────────────────────────────────────────────────────────
# 4. min_count drops low-support lengths
# ──────────────────────────────────────────────────────────────────


def test_min_count_drops_low_support_length_values() -> None:
    """A length value whose observed count is strictly
    below `min_count` is dropped before normalisation;
    surfaces in `below_min_count` per key."""
    builder = _builder_vj()
    # NP1=0 observed 5×, NP1=7 observed 1×.
    records = (
        [{"np1_length": 0}] * 5
        + [{"np1_length": 7}] * 1
    )
    builder.estimate_np_length_distributions(records, min_count=2)
    stage = next(
        s for s in builder.report().stages
        if s["stage"] == "estimate_np_length_distributions"
    )
    np1 = dict(stage["inferred"]["NP1"])
    # NP1=7 (count=1) drops; NP1=0 (count=5) kept.
    assert set(np1.keys()) == {0}
    assert math.isclose(np1[0], 1.0, rel_tol=1e-9)
    assert stage["inferred"]["below_min_count"]["NP1"] == 1


# ──────────────────────────────────────────────────────────────────
# 5. pseudocount affects weights
# ──────────────────────────────────────────────────────────────────


def test_pseudocount_applies_to_observed_values_only() -> None:
    """`pseudocount` is added to each observed length
    value's count before normalisation. v1 does NOT expand
    the support (unobserved values stay unobserved). With
    raw counts {0: 3, 1: 1} and pseudocount=1.0:

    - kept = {0: 4, 1: 2}
    - normalised = {0: 4/6, 1: 2/6}
    """
    builder = _builder_vj()
    records = (
        [{"np1_length": 0}] * 3
        + [{"np1_length": 1}] * 1
    )
    builder.estimate_np_length_distributions(records, min_count=1, pseudocount=1.0)
    stage = next(
        s for s in builder.report().stages
        if s["stage"] == "estimate_np_length_distributions"
    )
    np1 = dict(stage["inferred"]["NP1"])
    assert math.isclose(np1[0], 4 / 6, rel_tol=1e-9)
    assert math.isclose(np1[1], 2 / 6, rel_tol=1e-9)
    # Unobserved values (e.g. NP1=5) stay absent.
    assert 5 not in np1


# ──────────────────────────────────────────────────────────────────
# 6. replace=False rejects existing model
# ──────────────────────────────────────────────────────────────────


def test_replace_false_rejects_existing_np_lengths() -> None:
    """With `replace=False`, a second
    `estimate_np_length_distributions` call AND a manual
    `with_models(...)` attach of a typed plane both block
    re-entry."""
    builder = _builder_vj()
    records = [{"np1_length": 0}]
    builder.estimate_np_length_distributions(records, min_count=1)
    with pytest.raises(ValueError, match="already attached"):
        builder.estimate_np_length_distributions(
            records, min_count=1, replace=False
        )
    # Default replace=True allows the overwrite.
    builder.estimate_np_length_distributions(records, min_count=1)


def test_replace_false_rejects_manually_attached_np_lengths() -> None:
    """`replace=False` also rejects if `with_models()`
    already attached a typed plane carrying
    `np_lengths`."""
    builder = _builder_vj()
    builder.with_models(
        ReferenceEmpiricalModels(
            np_lengths={"NP1": EmpiricalDistributionSpec([(0, 1.0)])}
        )
    )
    records = [{"np1_length": 0}]
    with pytest.raises(ValueError, match="already attached"):
        builder.estimate_np_length_distributions(
            records, min_count=1, replace=False
        )


# ──────────────────────────────────────────────────────────────────
# 7. Built cartridge uses estimated NP lengths in recombination defaults
# ──────────────────────────────────────────────────────────────────


def test_built_cartridge_uses_estimated_np_lengths_in_recombination_defaults() -> None:
    """The bridge resolver's typed-plane > legacy precedence
    picks up the estimator's output. Confirmed via
    `extract_recombine_defaults`: the returned `np1`
    distribution matches what the estimator wrote."""
    from GenAIRR._dataconfig_extract import extract_recombine_defaults

    builder = _builder_vj()
    records = [{"np1_length": 0}] * 3 + [{"np1_length": 1}]
    builder.estimate_np_length_distributions(records, min_count=1)
    cfg = builder.build()
    defaults = extract_recombine_defaults(cfg)
    np1 = defaults["np1"]
    assert np1 is not None
    as_dict = dict(np1)
    assert math.isclose(as_dict[0], 0.75, rel_tol=1e-9)
    assert math.isclose(as_dict[1], 0.25, rel_tol=1e-9)


# ──────────────────────────────────────────────────────────────────
# 8. Explicit recombine override still wins
# ──────────────────────────────────────────────────────────────────


def test_explicit_recombine_np_lengths_kwarg_takes_precedence() -> None:
    """`Experiment.recombine(np1_lengths=[(8, 1.0)])`
    overrides the cartridge plane. Verified via the
    compiled pipeline IR's `np1_lengths` field — the
    surface the engine consumes."""
    builder = _builder_vj()
    builder.estimate_np_length_distributions(
        [{"np1_length": 0}] * 4, min_count=1
    )
    cfg = builder.build()
    exp = (
        ga.Experiment.on(cfg)
        .allow_curatable_refdata()
        .recombine(np1_lengths=[(8, 1.0)])
    )
    from GenAIRR._pipeline_ir import _RecombineStep
    recombine_step = next(s for s in exp._steps if isinstance(s, _RecombineStep))
    pairs = dict(recombine_step.np1_lengths)
    assert pairs == {8: 1.0}, (
        f"explicit recombine(np1_lengths=...) override didn't beat "
        f"the typed plane; pipeline IR carries {pairs!r}"
    )


# ──────────────────────────────────────────────────────────────────
# 9. P-nucleotide / sequence fields ignored
# ──────────────────────────────────────────────────────────────────


def test_p_nucleotide_and_sequence_fields_are_ignored() -> None:
    """`p_v_3_length` / `p_d_5_length` / `p_d_3_length` /
    `p_j_5_length` / `np1` / `np2` / `junction_length`
    are all silently ignored by the estimator (no error,
    no contribution, no warning). The estimator's input
    surface is exactly the two integer length columns."""
    builder = _builder_vj()
    records = [
        {
            "np1_length": 0,
            # P-nucleotide fields (separate biology).
            "p_v_3_length": 7, "p_d_5_length": 11,
            "p_d_3_length": 13, "p_j_5_length": 17,
            # Sequence fields (post-claim-reabsorption issue).
            "np1": "AAAAAAAAA",  # len=9 — disagrees with np1_length=0
            "np2": "TTTTTT",
            # Junction-length arithmetic.
            "junction_length": 51,
        },
        {"np1_length": 1, "p_v_3_length": 99, "np1": "GG"},
    ]
    builder.estimate_np_length_distributions(records, min_count=1)
    cfg = builder.build()
    np1 = dict(cfg.reference_models.np_lengths["NP1"].values)
    # Estimator wrote np1_length values, not derived from p_* / np1.
    assert np1 == {0: 0.5, 1: 0.5}, (
        f"NP1 distribution diverged from `np1_length` source: {np1!r}; "
        f"a sibling column (p_*, np1 string, junction_length) was "
        f"consumed by mistake"
    )
    # No warnings mention P-nucleotide / sequence / junction fields.
    stage = next(
        s for s in builder.report().stages
        if s["stage"] == "estimate_np_length_distributions"
    )
    for w in stage["warnings"]:
        for forbidden in ("p_v_3", "p_d_5", "p_d_3", "p_j_5",
                          "junction_length", "'np1'", "'np2'"):
            assert forbidden not in w, (
                f"warning {w!r} mentions {forbidden!r} — the "
                f"estimator should silently ignore these columns"
            )
    # dropped_columns is empty on VDJ; on VJ it only tracks np2_length.
    assert "p_v_3_length" not in stage["inferred"]["dropped_columns"]
    assert "junction_length" not in stage["inferred"]["dropped_columns"]
    assert "np1" not in stage["inferred"]["dropped_columns"]


# ──────────────────────────────────────────────────────────────────
# 10. Report JSON-clean and stage shape stable
# ──────────────────────────────────────────────────────────────────


def test_stage_entry_shape_matches_audit_and_report_is_json_clean() -> None:
    """The stage entry carries the canonical `{stage,
    inputs, inferred, warnings}` keys per the user brief.
    `inputs` has `record_count` / `min_count` /
    `pseudocount` / `source` / `replaced`. `inferred` has
    `NP1` / `NP2` / `skipped` / `below_min_count` /
    `dropped_columns`. The whole report round-trips
    through `json.dumps`."""
    builder = _builder_vdj()
    records = [{"np1_length": 0, "np2_length": 0}]
    builder.estimate_np_length_distributions(records, min_count=1)
    stage = next(
        s for s in builder.report().stages
        if s["stage"] == "estimate_np_length_distributions"
    )
    assert set(stage.keys()) >= {"stage", "inputs", "inferred", "warnings"}
    for k in ("record_count", "min_count", "pseudocount", "source", "replaced"):
        assert k in stage["inputs"], f"missing inputs key {k!r}"
    for k in ("NP1", "NP2", "skipped", "below_min_count", "dropped_columns"):
        assert k in stage["inferred"], f"missing inferred key {k!r}"
    # JSON-clean round-trip.
    blob = json.dumps(builder.report().to_dict())
    again = json.loads(blob)
    assert isinstance(again, dict)
    assert "stages" in again


# ──────────────────────────────────────────────────────────────────
# 11. Manifest block reflects estimated keys
# ──────────────────────────────────────────────────────────────────


def test_manifest_np_length_models_block_reflects_estimator_output() -> None:
    """After estimation, the built cartridge's
    `manifest['models']['np_length_models']` reflects the
    typed plane's keys; `in_plan_signature=True` (no soft
    gap inherited)."""
    builder = _builder_vdj()
    records = [{"np1_length": 0, "np2_length": 0}]
    builder.estimate_np_length_distributions(records, min_count=1)
    cfg = builder.build()
    block = cfg.cartridge_manifest()["models"]["np_length_models"]
    assert set(block["keys"]) == {"NP1", "NP2"}
    assert block["source"] == "ReferenceEmpiricalModels.np_lengths"
    assert block["in_plan_signature"] is True  # no soft gap
    assert block["legacy_fallback"] is False
    # Built cartridge has no legacy `NP_lengths`.
    assert block["legacy_np_lengths_present"] is False


# ──────────────────────────────────────────────────────────────────
# 12. Contract absence pins flip (integration probe)
# ──────────────────────────────────────────────────────────────────


def test_contract_pins_post_slice_state() -> None:
    """Integration probe — the contract file's absence pins
    are flipped to present-state, while the soft-gap pin
    (NP-length distributions DO fold into plan signature
    — NO inherited soft gap) still holds."""
    assert hasattr(
        ga.ReferenceCartridgeBuilder, "estimate_np_length_distributions"
    )
    m = ga.HUMAN_IGH_OGRDB.cartridge_manifest()
    assert "np_length_models" in m["models"]
    assert m["models"]["np_length_models"]["in_plan_signature"] is True
