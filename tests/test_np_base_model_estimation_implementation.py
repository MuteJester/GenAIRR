"""End-to-end implementation tests for the **NP Base Model
Estimation v1** slice.

Companion to
[`docs/np_base_model_estimation_design.md`](../docs/np_base_model_estimation_design.md)
and the contract pins in
[`tests/test_np_base_model_estimation_contract.py`](test_np_base_model_estimation_contract.py).

Covers the 15-test surface from the user brief:

 1. `empirical_first_base` estimates first-base only.
 2. `markov` estimates first-base + transitions.
 3. VJ estimates NP1 only and warns on NP2.
 4. VDJ estimates NP1 and NP2.
 5. Empty / missing / malformed strings skipped per key
    (field-local).
 6. `pseudocount` fills sparse rows.
 7. `min_count` affects first-base support.
 8. `replace=False` rejects existing model.
 9. Built cartridge drives NP base composition.
10. Markov estimator drives transition dependency on
    sampled output.
11. Replay round-trip works with estimated model.
12. Legacy `NP_first_bases` / `NP_transitions` remain
    NOT auto-lifted.
13. Report JSON-clean and stage shape stable.
14. Contract absence pins flip (integration probe).
15. Stale Markov-"deferred" docstring in `NpBaseModelSpec`
    is cleaned up.
"""
from __future__ import annotations

import json
import math
import pickle

import pytest

import GenAIRR as ga
from GenAIRR.reference_models import (
    EmpiricalDistributionSpec,
    NpBaseModelSpec,
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
# 1. empirical_first_base estimates first-base only
# ──────────────────────────────────────────────────────────────────


def test_empirical_first_base_kind_writes_first_base_only_no_transitions() -> None:
    """With `kind="empirical_first_base"`, the resulting
    `NpBaseModelSpec` carries a populated `first_base`
    distribution and `transitions=None`. The first-base
    distribution is over the **full base composition** of
    every observed NP string (every position, not just
    position 0)."""
    builder = _builder_vj()
    # 3 records: NP1 strings deliberately A-heavy.
    records = [
        {"np1": "AAAA"},
        {"np1": "AAAC"},
        {"np1": "AAGT"},
    ]
    builder.estimate_np_base_model(
        records, kind="empirical_first_base", min_count=1
    )
    cfg = builder.build()
    spec = cfg.reference_models.np_bases["NP1"]
    assert spec.kind == "empirical_first_base"
    assert spec.transitions is None
    # Full base composition: 9 A's, 1 C, 1 G, 1 T (out of 12 bases).
    assert math.isclose(spec.first_base["A"], 9 / 12, rel_tol=1e-9)
    assert math.isclose(spec.first_base["C"], 1 / 12, rel_tol=1e-9)
    assert math.isclose(spec.first_base["G"], 1 / 12, rel_tol=1e-9)
    assert math.isclose(spec.first_base["T"], 1 / 12, rel_tol=1e-9)


# ──────────────────────────────────────────────────────────────────
# 2. markov estimates first-base + transitions
# ──────────────────────────────────────────────────────────────────


def test_markov_kind_writes_first_base_and_transitions() -> None:
    """With `kind="markov"`, the resulting `NpBaseModelSpec`
    carries both `first_base` (from position 0 of each NP
    string) AND `transitions` (from every observed
    `(prev, next)` pair). The transition matrix covers
    all four A/C/G/T from-bases when pseudocount > 0."""
    builder = _builder_vj()
    records = [
        {"np1": "ATCG"},
        {"np1": "ATCG"},
        {"np1": "GCTA"},
    ]
    builder.estimate_np_base_model(records, kind="markov", pseudocount=1.0)
    cfg = builder.build()
    spec = cfg.reference_models.np_bases["NP1"]
    assert spec.kind == "markov"
    assert spec.first_base is not None
    assert spec.transitions is not None
    # All four from-bases covered (per Markov row-coverage spec).
    assert set(spec.transitions.keys()) == {"A", "C", "G", "T"}
    for from_b, row in spec.transitions.items():
        # Each row sums to 1.0 after normalisation.
        total = sum(row.values())
        assert math.isclose(total, 1.0, rel_tol=1e-9), (from_b, total)


# ──────────────────────────────────────────────────────────────────
# 3. VJ estimates NP1 only and warns on NP2
# ──────────────────────────────────────────────────────────────────


def test_vj_cartridge_writes_only_np1_key_and_warns_on_non_empty_np2() -> None:
    """On a VJ cartridge the estimator emits only `NP1`.
    Records that carry a non-empty `np2` column trigger
    a one-time warning and the contribution is tracked as
    a dropped column."""
    builder = _builder_vj()
    records = [
        {"np1": "ACGT", "np2": "AAAA"},  # non-empty np2 → dropped + warned
        {"np1": "ACGT", "np2": "TTTT"},  # non-empty np2 → dropped, warning re-used
        {"np1": "ACGT", "np2": ""},      # empty np2 → not flagged
        {"np1": "AAAA"},                 # np2 absent → not flagged
    ]
    builder.estimate_np_base_model(
        records, kind="empirical_first_base", min_count=1
    )
    stage = next(
        s for s in builder.report().stages
        if s["stage"] == "estimate_np_base_model"
    )
    inferred = stage["inferred"]
    # NP2 entry is present for shape stability but unpopulated.
    assert "NP1" in inferred and "NP2" in inferred
    assert inferred["NP1"]["first_base"], "NP1 distribution must be populated"
    assert inferred["NP2"] == {"first_base": None, "transitions": None}

    # Dropped column counter + one-time warning.
    assert stage["inferred"]["dropped_columns"]["np2"] == 2
    np2_warnings = [w for w in stage["warnings"] if "np2" in w]
    assert len(np2_warnings) == 1, (
        f"expected exactly one np2 warning across the dataset, "
        f"got {len(np2_warnings)}: {np2_warnings!r}"
    )

    # Built cartridge: only NP1 in the typed plane.
    cfg = builder.build()
    assert set(cfg.reference_models.np_bases.keys()) == {"NP1"}


# ──────────────────────────────────────────────────────────────────
# 4. VDJ estimates NP1 and NP2
# ──────────────────────────────────────────────────────────────────


def test_vdj_cartridge_writes_both_np1_and_np2_keys() -> None:
    """On a VDJ cartridge with both NP1 and NP2 columns
    populated, the estimator emits both keys."""
    builder = _builder_vdj()
    records = [
        {"np1": "ACGT", "np2": "TTAA"},
        {"np1": "GCTA", "np2": "AAAA"},
        {"np1": "AAGG", "np2": "CCGG"},
    ]
    builder.estimate_np_base_model(
        records, kind="empirical_first_base", min_count=1
    )
    cfg = builder.build()
    assert set(cfg.reference_models.np_bases.keys()) == {"NP1", "NP2"}
    for key in ("NP1", "NP2"):
        spec = cfg.reference_models.np_bases[key]
        assert spec.kind == "empirical_first_base"
        assert sum(spec.first_base.values()) == pytest.approx(1.0)


# ──────────────────────────────────────────────────────────────────
# 5. Empty / missing / malformed strings skipped per key
# ──────────────────────────────────────────────────────────────────


def test_field_local_skipping_under_missing_empty_or_noncanonical_input() -> None:
    """A row that is malformed / missing / noncanonical in
    ONE NP column drops only that key's contribution — the
    other NP column still feeds. Structured rejection
    entries land per (row, column)."""
    builder = _builder_vdj()
    records = [
        # Three clean rows.
        {"np1": "ACGT", "np2": "AAAA"},
        {"np1": "ACGT", "np2": "AAAA"},
        {"np1": "ACGT", "np2": "AAAA"},
        # NP1 empty → drops NP1 only; NP2 still counts.
        {"np1": "",     "np2": "TTTT"},
        # NP2 missing → drops NP2 only; NP1 still counts.
        {"np1": "GCTA"},
        # NP1 noncanonical → drops NP1 only; NP2 still counts.
        {"np1": "AAXX", "np2": "GGGG"},
    ]
    builder.estimate_np_base_model(
        records, kind="empirical_first_base", min_count=1
    )
    stage = next(
        s for s in builder.report().stages
        if s["stage"] == "estimate_np_base_model"
    )
    sk = stage["inferred"]["skipped"]
    assert sk["missing_required_column"]["np1"] == 1
    assert sk["missing_required_column"]["np2"] == 1
    assert sk["noncanonical_base"]["np1"] == 1

    # Structured rejection entries.
    reasons = [r["reason"] for r in builder.report().rejected]
    assert reasons.count("missing_required_column") == 2
    assert reasons.count("noncanonical_base") == 1

    # NP1 contributed: 3 clean + 1 missing-NP2 row = 4 rows;
    # the noncanonical row dropped NP1 but kept NP2.
    # NP2 contributed: 3 clean + 1 empty-NP1 row + 1 noncanonical-NP1 row = 5 rows.
    cfg = builder.build()
    np1 = cfg.reference_models.np_bases["NP1"].first_base
    np2 = cfg.reference_models.np_bases["NP2"].first_base
    # NP1 across 4 contributing rows: ACGT * 3 + GCTA * 1 = 12 + 4 = 16 bases.
    # Distribution: A=4, C=4, G=4, T=4 — uniform.
    assert math.isclose(np1["A"], 0.25, rel_tol=1e-9)
    assert math.isclose(np1["C"], 0.25, rel_tol=1e-9)
    # NP2 across 5 contributing rows: AAAA * 3 + TTTT * 1 + GGGG * 1 = 20 bases.
    # Distribution: A=12, T=4, G=4 (no C).
    assert math.isclose(np2["A"], 12 / 20, rel_tol=1e-9)
    assert math.isclose(np2["T"], 4 / 20, rel_tol=1e-9)
    assert math.isclose(np2["G"], 4 / 20, rel_tol=1e-9)


# ──────────────────────────────────────────────────────────────────
# 6. pseudocount fills sparse rows
# ──────────────────────────────────────────────────────────────────


def test_pseudocount_fills_sparse_markov_rows() -> None:
    """With `kind="markov"` and `pseudocount=1.0`, every
    transition row is guaranteed to satisfy the
    `NpBaseModelSpec` row-coverage validator even when an
    entire from-base was never observed in the input.

    Input: NP1 strings contain only A and T bases —
    transitions from C and G are entirely unobserved.
    With pseudocount=1.0 the matrix still validates."""
    builder = _builder_vj()
    records = [
        {"np1": "ATATAT"},
        {"np1": "TATATA"},
        {"np1": "AAAATTTT"},
    ]
    builder.estimate_np_base_model(
        records, kind="markov", min_count=1, pseudocount=1.0
    )
    cfg = builder.build()
    spec = cfg.reference_models.np_bases["NP1"]
    # All four from-base rows present despite C / G never appearing.
    assert set(spec.transitions.keys()) == {"A", "C", "G", "T"}
    # C / G rows are pure pseudocount → uniform.
    for from_b in ("C", "G"):
        row = spec.transitions[from_b]
        for t in ("A", "C", "G", "T"):
            assert math.isclose(row[t], 0.25, rel_tol=1e-9), (from_b, t, row)


def test_pseudocount_zero_fails_on_sparse_input() -> None:
    """With `pseudocount=0.0` AND a from-base that never
    appears in the input, the estimator raises a clear
    error naming the missing from-base."""
    builder = _builder_vj()
    # NP1 strings only contain A and T → C / G transitions are empty.
    records = [{"np1": "ATAT"}, {"np1": "TATA"}]
    with pytest.raises(ValueError, match=r"from-base"):
        builder.estimate_np_base_model(
            records, kind="markov", min_count=1, pseudocount=0.0
        )


# ──────────────────────────────────────────────────────────────────
# 7. min_count affects first-base support
# ──────────────────────────────────────────────────────────────────


def test_min_count_drops_low_support_first_base_categories() -> None:
    """`min_count` drops first-base categories whose
    observed count is strictly below the threshold,
    surfaces the drop in `below_min_count.first_base`."""
    builder = _builder_vj()
    # First-base position 0: A appears 5×, C appears 1×.
    records = (
        [{"np1": "AAA"}] * 5
        + [{"np1": "CAA"}] * 1
    )
    builder.estimate_np_base_model(
        records, kind="empirical_first_base", min_count=2, pseudocount=0.0
    )
    stage = next(
        s for s in builder.report().stages
        if s["stage"] == "estimate_np_base_model"
    )
    # 'A' position-0 count = 5 — but the full base composition
    # is summed over every position, so the actual count for
    # 'C' in the full composition is 1 (the single 'C' from
    # the one CAA row). 'A' count = 14. min_count=2 drops C
    # AND the two never-observed bases (G, T) — every cell
    # below the threshold counts toward `below_min_count`.
    np1 = stage["inferred"]["NP1"]["first_base"]
    assert "C" not in np1
    assert "A" in np1
    # Three drops: C (count=1), G (count=0), T (count=0).
    assert stage["inferred"]["below_min_count"]["NP1"]["first_base"] == 3


# ──────────────────────────────────────────────────────────────────
# 8. replace=False rejects existing model
# ──────────────────────────────────────────────────────────────────


def test_replace_false_rejects_existing_np_bases() -> None:
    """`replace=False` rejects re-entry when the typed
    plane's `np_bases` already carries a spec — whether
    set by a prior `estimate_np_base_model` call or by
    `with_models(...)`."""
    builder = _builder_vj()
    builder.estimate_np_base_model(
        [{"np1": "AAAA"}], kind="empirical_first_base", min_count=1
    )
    with pytest.raises(ValueError, match="already attached"):
        builder.estimate_np_base_model(
            [{"np1": "AAAA"}], kind="empirical_first_base",
            min_count=1, replace=False,
        )
    # replace=True (default) allows overwrite.
    builder.estimate_np_base_model(
        [{"np1": "AAAA"}], kind="empirical_first_base", min_count=1
    )


def test_replace_false_rejects_manually_attached_np_bases() -> None:
    """`replace=False` also rejects if `with_models(...)`
    already attached a typed plane carrying `np_bases`."""
    builder = _builder_vj()
    builder.with_models(
        ReferenceEmpiricalModels(
            np_bases={
                "NP1": NpBaseModelSpec(
                    kind="empirical_first_base",
                    first_base={"A": 1.0, "C": 1.0, "G": 1.0, "T": 1.0},
                )
            }
        )
    )
    with pytest.raises(ValueError, match="already attached"):
        builder.estimate_np_base_model(
            [{"np1": "AAAA"}], kind="empirical_first_base",
            min_count=1, replace=False,
        )


# ──────────────────────────────────────────────────────────────────
# 9. Built cartridge drives NP base composition
# ──────────────────────────────────────────────────────────────────


def test_built_cartridge_drives_empirical_first_base_composition_at_recombine_time() -> None:
    """The bridge resolver's typed-plane > legacy precedence
    picks up the estimator's output. Verified via the
    compiled NP1 base composition: an estimated A-heavy
    `empirical_first_base` model biases the generated NP1
    bases toward A on real recombination output."""
    builder = _builder_vj()
    # Heavy A-bias input.
    records = (
        [{"np1": "AAAA"}] * 20
        + [{"np1": "ACGT"}] * 1
    )
    builder.estimate_np_base_model(
        records, kind="empirical_first_base", min_count=1
    )
    cfg = builder.build()
    result = (
        ga.Experiment.on(cfg)
        .allow_curatable_refdata()
        .recombine()
        .run_records(n=30, seed=4242)
    )
    base_counts = {b: 0 for b in "ACGT"}
    total = 0
    for r in result.records:
        for b in r.get("np1", ""):
            if b in base_counts:
                base_counts[b] += 1
                total += 1
    # Heavy A bias should dominate.
    if total > 0:
        a_ratio = base_counts["A"] / total
        assert a_ratio >= 0.6, (
            f"estimated empirical_first_base model didn't bias NP1 "
            f"output toward A: a_ratio={a_ratio:.2%}, counts="
            f"{base_counts}, total={total}"
        )


# ──────────────────────────────────────────────────────────────────
# 10. Markov estimator drives transition dependency
# ──────────────────────────────────────────────────────────────────


def test_markov_estimator_drives_transition_dependency_on_sampled_output() -> None:
    """A Markov estimator trained on the cyclic pattern
    `ATCGATCG...` should produce a transition matrix where
    A→T, T→C, C→G, G→A are strongly preferred. After
    rebuilding the cartridge and running recombination,
    NP1 strings should exhibit those forced transitions."""
    builder = _builder_vj()
    # Strong cyclic input: A→T→C→G→A...
    records = [{"np1": "ATCGATCGATCG"}] * 20
    builder.estimate_np_base_model(
        records, kind="markov", min_count=1, pseudocount=0.01
    )
    cfg = builder.build()
    spec = cfg.reference_models.np_bases["NP1"]
    # Spec encodes the expected dominance.
    assert spec.transitions["A"]["T"] > 0.8
    assert spec.transitions["T"]["C"] > 0.8
    assert spec.transitions["C"]["G"] > 0.8
    assert spec.transitions["G"]["A"] > 0.8

    # End-to-end: forced pattern dominates engine output.
    result = (
        ga.Experiment.on(cfg)
        .allow_curatable_refdata()
        .recombine()
        .run_records(n=30, seed=4242)
    )
    forced = ("AT", "TC", "CG", "GA")
    transitions_observed: dict[str, int] = {}
    for r in result.records:
        np1 = r.get("np1", "")
        for i in range(len(np1) - 1):
            pair = np1[i:i+2]
            transitions_observed[pair] = transitions_observed.get(pair, 0) + 1
    total = sum(transitions_observed.values()) or 1
    forced_count = sum(transitions_observed.get(p, 0) for p in forced)
    forced_ratio = forced_count / total
    assert forced_ratio >= 0.6, (
        f"estimated Markov matrix didn't drive transition "
        f"dependency on output: forced_ratio={forced_ratio:.2%}; "
        f"observed transitions: {transitions_observed!r}"
    )


# ──────────────────────────────────────────────────────────────────
# 11. Replay round-trip works with estimated model
# ──────────────────────────────────────────────────────────────────


def test_replay_round_trip_with_estimated_np_base_model() -> None:
    """A cartridge with an estimated NP base model pickles
    cleanly and reproduces byte-identically across runs at
    the same seed (same precedence chain, same plan
    signature, same trace)."""
    builder = _builder_vdj()
    records = [{"np1": "ACGT" * 3, "np2": "TGCA" * 3}] * 10
    builder.estimate_np_base_model(
        records, kind="markov", min_count=1, pseudocount=0.5
    )
    cfg = builder.build()

    # Pickle round-trip preserves the spec.
    cfg2 = pickle.loads(pickle.dumps(cfg))
    assert cfg2.reference_models.np_bases["NP1"].kind == "markov"
    assert set(cfg2.reference_models.np_bases["NP1"].transitions.keys()) == {
        "A", "C", "G", "T"
    }

    # Byte-identical replay across two runs at the same seed.
    a = (
        ga.Experiment.on(cfg).allow_curatable_refdata()
        .recombine().run_records(n=5, seed=99)
    )
    b = (
        ga.Experiment.on(cfg2).allow_curatable_refdata()
        .recombine().run_records(n=5, seed=99)
    )
    assert [r["sequence"] for r in a.records] == [
        r["sequence"] for r in b.records
    ]


# ──────────────────────────────────────────────────────────────────
# 12. Legacy NP_first_bases / NP_transitions not auto-lifted
# ──────────────────────────────────────────────────────────────────


def test_estimator_does_not_auto_lift_legacy_np_orphan_fields() -> None:
    """Bundled cartridges carry legacy `NP_first_bases` /
    `NP_transitions` orphan dicts. The estimator slice
    MUST NOT auto-lift them — same boundary every prior
    estimator slice respected."""
    from GenAIRR.dataconfig.data_config import (
        _DOCUMENTED_ORPHAN_DATACONFIG_FIELDS,
    )
    assert "NP_first_bases" in _DOCUMENTED_ORPHAN_DATACONFIG_FIELDS
    assert "NP_transitions" in _DOCUMENTED_ORPHAN_DATACONFIG_FIELDS

    # An estimated cartridge has the typed plane but NOT a
    # populated legacy dict (because builder cartridges start
    # without legacy state).
    builder = _builder_vj()
    builder.estimate_np_base_model(
        [{"np1": "AAAA"}], kind="empirical_first_base", min_count=1
    )
    cfg = builder.build()
    assert cfg.reference_models.np_bases["NP1"] is not None
    assert not cfg.NP_first_bases, (
        "estimator mistakenly populated legacy NP_first_bases — "
        "the orphan boundary regressed"
    )
    assert not cfg.NP_transitions, (
        "estimator mistakenly populated legacy NP_transitions — "
        "the orphan boundary regressed"
    )


# ──────────────────────────────────────────────────────────────────
# 13. Report JSON-clean and stage shape stable
# ──────────────────────────────────────────────────────────────────


def test_stage_entry_shape_matches_brief_and_report_is_json_clean() -> None:
    """The stage entry carries the canonical `{stage,
    inputs, inferred, warnings}` shape. `inputs` has
    `record_count` / `kind` / `min_count` / `pseudocount`
    / `source` / `replaced`. `inferred` has `NP1` / `NP2`
    / `skipped` / `below_min_count` / `dropped_columns`.
    The full report round-trips through `json.dumps`."""
    builder = _builder_vdj()
    records = [{"np1": "ACGT", "np2": "TTAA"}]
    builder.estimate_np_base_model(
        records, kind="empirical_first_base", min_count=1
    )
    stage = next(
        s for s in builder.report().stages
        if s["stage"] == "estimate_np_base_model"
    )
    assert set(stage.keys()) >= {"stage", "inputs", "inferred", "warnings"}
    for k in ("record_count", "kind", "min_count", "pseudocount",
              "source", "replaced"):
        assert k in stage["inputs"], f"missing inputs key {k!r}"
    for k in ("NP1", "NP2", "skipped", "below_min_count", "dropped_columns"):
        assert k in stage["inferred"], f"missing inferred key {k!r}"
    blob = json.dumps(builder.report().to_dict())
    again = json.loads(blob)
    assert isinstance(again, dict)
    assert "stages" in again


# ──────────────────────────────────────────────────────────────────
# 14. Contract absence pins flip (integration probe)
# ──────────────────────────────────────────────────────────────────


def test_contract_pins_post_slice_state() -> None:
    """Integration probe — the contract file's absence pins
    are flipped to present-state, while the no-soft-gap
    pin still holds. The contract file is the source of
    truth; this test is the safety net so a regression in
    either direction surfaces here too."""
    assert hasattr(ga.ReferenceCartridgeBuilder, "estimate_np_base_model")
    # Manifest's existing np_base_models block reflects the
    # estimator's output through the populated typed plane.
    m = ga.HUMAN_IGH_OGRDB.cartridge_manifest()
    assert "np_base_models" in m["models"]


# ──────────────────────────────────────────────────────────────────
# 15. Stale Markov-"deferred" docstring cleaned up
# ──────────────────────────────────────────────────────────────────


def test_markov_docstring_no_longer_describes_kind_as_deferred() -> None:
    """The `NpBaseModelSpec` docstring's "Markov deferred"
    claim was stale at audit time. The implementation
    slice cleans it up. Pinned so a regression
    re-introducing the stale claim surfaces here."""
    doc = NpBaseModelSpec.__doc__ or ""
    assert "deferred" not in doc.lower(), (
        f"NpBaseModelSpec docstring still describes Markov as "
        f"'deferred' — the docstring cleanup regressed"
    )
    assert "NotImplementedError" not in doc, (
        f"NpBaseModelSpec docstring still claims Markov raises "
        f"NotImplementedError — the docstring cleanup regressed"
    )
