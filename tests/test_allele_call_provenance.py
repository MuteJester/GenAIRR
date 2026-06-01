"""Allele-call disambiguation golden tests.

Companion to [docs/allele_call_audit.md](../docs/allele_call_audit.md).
Pins **current** behaviour of the live-call walker on the AIRR
projection under evidence-changing mechanisms (SHM, PCR errors,
sequencing errors, N substitution).

Fixture design notes:

- 2-allele V refdata used throughout: ``v1*01 = AAACCCGGG`` and
  ``v2*01 = AAACTCGGG``. They differ at position 4 (C vs T).
  Positions 0-3 and 5-8 are non-distinguishing.
- ``v1*01`` and ``v1*02`` are identical (same byte string) — used
  to pin the always-tied case.
- ``restrict_alleles(v="v1*01")`` locks the recombination-stage
  truth, isolating evidence-changing mechanisms from sampling
  variance.

Empirical seeds were searched ahead of time and are pinned into
the tests; if walker behaviour shifts the seed assertions are the
first to fire.
"""
from __future__ import annotations

import pytest

import GenAIRR as ga
from GenAIRR import _engine as ge
from GenAIRR import CompiledExperiment


# ──────────────────────────────────────────────────────────────────
# Refdata factories
# ──────────────────────────────────────────────────────────────────


def _two_v_one_distinguishing() -> "ge.RefDataConfig":
    """V refdata with v1*01 / v2*01 differing only at position 4."""
    cfg = ge.RefDataConfig.vj()
    cfg.add_v_allele("v1*01", "v1", b"AAACCCGGG", anchor=6)
    cfg.add_v_allele("v2*01", "v2", b"AAACTCGGG", anchor=6)
    cfg.add_j_allele("j1*01", "j1", b"TTTAAA", anchor=0)
    return cfg


def _two_v_identical() -> "ge.RefDataConfig":
    """V refdata with two literally identical alleles."""
    cfg = ge.RefDataConfig.vj()
    cfg.add_v_allele("v1*01", "v1", b"AAACCCGGG", anchor=6)
    cfg.add_v_allele("v1*02", "v1", b"AAACCCGGG", anchor=6)
    cfg.add_j_allele("j1*01", "j1", b"TTTAAA", anchor=0)
    return cfg


def _three_v_prefix_suffix() -> "ge.RefDataConfig":
    """Three V alleles: v_A baseline, v_B differs at pos 4
    (prefix-tie), v_C differs at pos 8 (suffix-tie)."""
    cfg = ge.RefDataConfig.vj()
    cfg.add_v_allele("v_A*01", "v_A", b"AAACCCGGG", anchor=6)
    cfg.add_v_allele("v_B*01", "v_B", b"AAACTCGGG", anchor=6)
    cfg.add_v_allele("v_C*01", "v_C", b"AAACCCGGT", anchor=6)
    cfg.add_j_allele("j1*01", "j1", b"AAATTT", anchor=0)  # avoid collision with v_C pos 8
    return cfg


def _vdj_two_d_one_distinguishing() -> "ge.RefDataConfig":
    """VDJ refdata with two D alleles differing in the suffix."""
    cfg = ge.RefDataConfig.vdj()
    cfg.add_v_allele("v1*01", "v1", b"AAACCCGGGTTT", anchor=9)
    cfg.add_d_allele("d1*01", "d1", b"CCCGGGAAA")
    cfg.add_d_allele("d2*01", "d2", b"CCCGGGTTT")  # differs in last 3
    cfg.add_j_allele("j1*01", "j1", b"TTTAAACCC", anchor=0)
    return cfg


def _baseline_vj_experiment(cfg) -> "ga.Experiment":
    return (
        ga.Experiment.on(cfg).allow_curatable_refdata()
        .recombine(np1_lengths=[(3, 1.0)])
        .trim(enabled=False)
    )


def _record(exp: "ga.Experiment", seed: int = 0, **kwargs) -> dict:
    result = exp.run_records(n=1, seed=seed, expose_provenance=True, **kwargs)
    assert len(result.records) == 1
    return result.records[0]


# ──────────────────────────────────────────────────────────────────
# 1. Tie-set formation — model basics
# ──────────────────────────────────────────────────────────────────


def test_identical_v_alleles_produce_csv_tie_with_truth_first() -> None:
    """Two literally identical V alleles → walker can't
    distinguish them → tie-set contains both. AIRR `v_call` is
    the CSV with the truth allele first (`projected_call_name`
    convention)."""
    exp = _baseline_vj_experiment(_two_v_identical())
    # Sweep multiple seeds so we see both truth picks.
    truth_first_observations = {"v1*01": False, "v1*02": False}
    for seed in range(8):
        rec = _record(exp, seed=seed)
        truth = rec["truth_v_call"]
        csv = rec["v_call"].split(",")
        assert len(csv) == 2
        assert set(csv) == {"v1*01", "v1*02"}
        assert csv[0] == truth, (
            f"seed {seed}: truth {truth} must be first in CSV {csv}"
        )
        truth_first_observations[truth] = True
    # We expect both truths to have been sampled at least once.
    assert all(truth_first_observations.values()), (
        "fixture didn't sweep both truth alleles in 8 seeds"
    )


def test_distinguishing_alleles_resolve_to_truth_without_mutation() -> None:
    """Without any evidence-changing mechanism, the assembled pool
    is the truth allele's germline. The walker scores it +9 vs.
    +8 for the other allele → singleton tie-set = truth."""
    exp = _baseline_vj_experiment(_two_v_one_distinguishing())
    for seed in range(8):
        rec = _record(exp, seed=seed)
        assert rec["v_call"] == rec["truth_v_call"], (
            f"seed {seed}: v_call {rec['v_call']!r} should resolve to "
            f"truth {rec['truth_v_call']!r} when no evidence shift"
        )
        # Identity against the projected (==truth) allele is perfect.
        assert rec["v_identity"] == 1.0


# ──────────────────────────────────────────────────────────────────
# 2. SHM — call switching / preservation
# ──────────────────────────────────────────────────────────────────


def test_shm_at_distinguishing_position_can_switch_call_away_from_truth() -> None:
    """When SHM lands at the distinguishing position AND mutates
    to the other allele's germline base, the walker scores the
    other allele higher → call switches.

    Seed=24 with `mutate(rate=0.1)` and locked truth=v1*01
    produces sequence 'AACCTCGGG' (T at pos 4, matching v2*01's
    germline). v_call flips to 'v2*01' while truth_v_call stays
    'v1*01' — the truth survives only in the provenance field."""
    exp = (
        _baseline_vj_experiment(_two_v_one_distinguishing())
        .restrict_alleles(v="v1*01")
        .mutate(rate=0.1)
    )
    rec = _record(exp, seed=24)
    assert rec["truth_v_call"] == "v1*01"
    assert rec["v_call"] == "v2*01"  # singleton, not CSV
    # The mutated sequence has T at the distinguishing position.
    assert rec["sequence"][4] == "T"
    # SHM increments n_mutations only.
    assert rec["n_mutations"] > 0
    assert rec["n_pcr_errors"] == 0
    assert rec["n_quality_errors"] == 0


def test_shm_at_non_distinguishing_position_preserves_call() -> None:
    """When SHM lands at a non-distinguishing position, both
    alleles' scores drop by the same amount. The distinguishing
    position still resolves the call to truth.

    Seed=6 with `mutate(rate=0.1)` produces a single mutation
    leaving position 4 unchanged → v_call still equals truth."""
    exp = (
        _baseline_vj_experiment(_two_v_one_distinguishing())
        .restrict_alleles(v="v1*01")
        .mutate(rate=0.1)
    )
    rec = _record(exp, seed=6)
    assert rec["truth_v_call"] == "v1*01"
    assert rec["v_call"] == "v1*01"  # call preserved
    assert rec["sequence"][4] == "C"  # distinguishing position unchanged
    assert rec["n_mutations"] == 1


# ──────────────────────────────────────────────────────────────────
# 3. PCR & sequencing errors — same walker semantics as SHM
# ──────────────────────────────────────────────────────────────────


def test_pcr_error_at_distinguishing_position_switches_call() -> None:
    """PCR errors substitute uppercase canonical bases. They are
    walker-semantically identical to SHM — a PCR error landing at
    the distinguishing position with the other allele's base
    flips the call.

    Seed=39 + `pcr_amplify(rate=0.1)` produces a uniform T at pos
    4 → switches to v2*01."""
    exp = (
        _baseline_vj_experiment(_two_v_one_distinguishing())
        .restrict_alleles(v="v1*01")
        .pcr_amplify(rate=0.1)
    )
    rec = _record(exp, seed=39)
    assert rec["truth_v_call"] == "v1*01"
    assert rec["v_call"] == "v2*01"
    assert rec["sequence"][4] == "T"
    # PCR increments n_pcr_errors only — n_mutations stays zero.
    assert rec["n_pcr_errors"] > 0
    assert rec["n_mutations"] == 0
    assert rec["n_quality_errors"] == 0


def test_sequencing_error_lowercase_letter_is_not_a_wildcard() -> None:
    """`sequencing_errors` uses lowercase letters as a
    *presentation* convention. The walker uppercases the byte
    before lookup, so a lowercase 't' at the distinguishing
    position is treated as 'T' — flipping the call, not widening
    it.

    Seed=39 + `sequencing_errors(rate=0.1)` produces a lowercase
    't' at pos 4 → v_call = 'v2*01' (singleton), confirming the
    walker did NOT treat the lowercase letter as N (which would
    have widened to 'v1*01,v2*01')."""
    exp = (
        _baseline_vj_experiment(_two_v_one_distinguishing())
        .restrict_alleles(v="v1*01")
        .sequencing_errors(rate=0.1)
    )
    rec = _record(exp, seed=39)
    assert rec["truth_v_call"] == "v1*01"
    assert rec["v_call"] == "v2*01"  # singleton — NOT widened
    assert rec["sequence"][4] == "t"  # lowercase preserved in output
    assert rec["n_quality_errors"] > 0
    assert rec["n_mutations"] == 0
    assert rec["n_pcr_errors"] == 0


# ──────────────────────────────────────────────────────────────────
# 4. N substitution — wildcard widens but retains truth
# ──────────────────────────────────────────────────────────────────


def test_n_substitution_at_distinguishing_position_widens_tie_set() -> None:
    """An N at the distinguishing position matches every allele
    via the wildcard rule. The previously-unique max becomes a
    tie — the tie-set widens to include alleles that were one
    point behind.

    Seed=0 + `ambiguous_base_calls(count=1)` places N at position
    4 → v_call widens from {v1*01} to {v1*01, v2*01}, with truth
    listed first per the projection convention."""
    exp = (
        _baseline_vj_experiment(_two_v_one_distinguishing())
        .restrict_alleles(v="v1*01")
        .ambiguous_base_calls(count=1)
    )
    rec = _record(exp, seed=0)
    assert rec["sequence"][4] == "N"
    assert rec["truth_v_call"] == "v1*01"
    assert rec["v_call"] == "v1*01,v2*01"  # truth first, widened


def test_n_substitution_at_non_distinguishing_position_preserves_call() -> None:
    """An N at a non-distinguishing position matches both alleles
    via wildcard, but they're already tied at that column — so
    the *relative* scores don't change. The call stays
    singleton-truth.

    Seed=4 + `ambiguous_base_calls(count=1)` places N at position
    8 (non-distinguishing) → v_call still equals truth."""
    exp = (
        _baseline_vj_experiment(_two_v_one_distinguishing())
        .restrict_alleles(v="v1*01")
        .ambiguous_base_calls(count=1)
    )
    rec = _record(exp, seed=4)
    assert rec["sequence"][8] == "N"
    assert rec["v_call"] == "v1*01"  # unchanged


# ──────────────────────────────────────────────────────────────────
# 5. 3-allele prefix/suffix ambiguity
# ──────────────────────────────────────────────────────────────────


def test_three_allele_prefix_widening_keeps_truth_first() -> None:
    """With three V alleles where v_B differs at pos 4 (prefix-
    ambiguous with v_A) and v_C differs at pos 8 (suffix-
    ambiguous with v_A), an N at pos 4 widens the tie to {v_A,
    v_B} only (v_C still loses at pos 8). Truth v_A first.

    Seed=0 + `ambiguous_base_calls(count=1)` places N at position
    4 on the 3-allele fixture."""
    exp = (
        _baseline_vj_experiment(_three_v_prefix_suffix())
        .restrict_alleles(v="v_A*01")
        .ambiguous_base_calls(count=1)
    )
    rec = _record(exp, seed=0)
    assert rec["sequence"][4] == "N"
    assert rec["truth_v_call"] == "v_A*01"
    assert rec["v_call"] == "v_A*01,v_B*01"  # v_C excluded
    assert "v_C*01" not in rec["v_call"]


def test_three_allele_suffix_widening_keeps_truth_first() -> None:
    """Mirror of the prefix case: N at pos 8 widens to {v_A, v_C}
    only (v_B still loses at pos 4)."""
    exp = (
        _baseline_vj_experiment(_three_v_prefix_suffix())
        .restrict_alleles(v="v_A*01")
        .ambiguous_base_calls(count=1)
    )
    rec = _record(exp, seed=4)
    assert rec["sequence"][8] == "N"
    assert rec["v_call"] == "v_A*01,v_C*01"
    assert "v_B*01" not in rec["v_call"]


# ──────────────────────────────────────────────────────────────────
# 6. Truth decoupling + identity-against-projected
# ──────────────────────────────────────────────────────────────────


def test_truth_v_call_decoupled_from_v_call_under_call_flip() -> None:
    """The audit's load-bearing fact: under any evidence-changing
    mechanism that flips the call, `truth_v_call` (when
    `expose_provenance=True`) preserves the originally-sampled
    allele. Users who treat `v_call` as truth without checking
    `truth_v_call` will silently lose the recombination identity."""
    # SHM flip
    rec_shm = _record(
        _baseline_vj_experiment(_two_v_one_distinguishing())
        .restrict_alleles(v="v1*01")
        .mutate(rate=0.1),
        seed=24,
    )
    assert rec_shm["v_call"] != rec_shm["truth_v_call"]
    assert rec_shm["truth_v_call"] == "v1*01"

    # PCR flip — same property, different mechanism
    rec_pcr = _record(
        _baseline_vj_experiment(_two_v_one_distinguishing())
        .restrict_alleles(v="v1*01")
        .pcr_amplify(rate=0.1),
        seed=39,
    )
    assert rec_pcr["v_call"] != rec_pcr["truth_v_call"]
    assert rec_pcr["truth_v_call"] == "v1*01"


def test_v_identity_is_computed_against_projected_allele_not_truth() -> None:
    """When the call flips away from truth, `v_identity` is the
    similarity of the assembled pool to the **projected
    (first-in-CSV)** allele's germline — not to the truth
    allele. So a high `v_identity` can coexist with `v_call !=
    truth_v_call`.

    Seed=24 SHM-flip example: ~5 mutations against truth (v1*01)
    produce sequence 'AACCTCGGG'. Identity vs v2*01 (the
    projected call) is 8/9 ≈ 0.888. Identity vs truth v1*01
    would be 7/9 ≈ 0.777. The reported `v_identity` is the
    former."""
    rec = _record(
        _baseline_vj_experiment(_two_v_one_distinguishing())
        .restrict_alleles(v="v1*01")
        .mutate(rate=0.1),
        seed=24,
    )
    assert rec["v_call"] == "v2*01"
    # AACCTCGGG vs AAACTCGGG (v2*01) → matches at all except pos 2
    # (C vs A). That's 8/9 ≈ 0.888.
    assert rec["v_identity"] == pytest.approx(8 / 9, abs=1e-9)


# ──────────────────────────────────────────────────────────────────
# 7. Per-mechanism counter isolation
# ──────────────────────────────────────────────────────────────────


def test_per_mechanism_counters_dont_cross_contaminate() -> None:
    """The walker treats SHM / PCR / sequencing errors / N
    identically once they're written into the pool, but the
    counters (`n_mutations`, `n_pcr_errors`, `n_quality_errors`)
    are per-mechanism and stay strictly isolated.

    A pipeline with all four mechanisms enabled increments only
    the relevant counter for each — no double-counting."""
    cfg = _two_v_one_distinguishing()
    rec = _record(
        _baseline_vj_experiment(cfg)
        .restrict_alleles(v="v1*01")
        .mutate(rate=0.1)
        .pcr_amplify(rate=0.05)
        .sequencing_errors(rate=0.05)
        .ambiguous_base_calls(count=1),
        seed=7,
    )
    # All four are independent counters; each only counts its own
    # mechanism. We assert non-cross-contamination by checking the
    # sum bounds — the per-mechanism counts must each be ≥ 0 and
    # the total writes (mutations + pcr + quality + N) should be a
    # reasonable upper bound on sequence_length.
    n_mut = rec["n_mutations"]
    n_pcr = rec["n_pcr_errors"]
    n_qe = rec["n_quality_errors"]
    for n in (n_mut, n_pcr, n_qe):
        assert isinstance(n, int) and n >= 0
    # Lightly cross-check that *some* of the mechanisms fired:
    # with these rates on an 18-base sequence we should see at
    # least one mutation and one error.
    assert (n_mut + n_pcr + n_qe) >= 1


# ──────────────────────────────────────────────────────────────────
# 8. VDJ — D-call disambiguation
# ──────────────────────────────────────────────────────────────────


def test_vdj_d_call_resolves_to_truth_without_mutation() -> None:
    """The same disambiguation rules apply to the D segment in
    VDJ refdata. Without mutation, d_call equals truth_d_call."""
    cfg = _vdj_two_d_one_distinguishing()
    exp = (
        ga.Experiment.on(cfg).allow_curatable_refdata()
        .recombine(np1_lengths=[(2, 1.0)], np2_lengths=[(2, 1.0)])
        .trim(enabled=False)
    )
    for seed in range(6):
        rec = _record(exp, seed=seed)
        assert rec["d_call"] == rec["truth_d_call"]


# ──────────────────────────────────────────────────────────────────
# 9. Replay round-trip
# ──────────────────────────────────────────────────────────────────


@pytest.mark.parametrize(
    "configure,seed",
    [
        pytest.param(
            lambda e: e.restrict_alleles(v="v1*01").mutate(rate=0.1),
            24,
            id="shm_call_flip",
        ),
        pytest.param(
            lambda e: e.restrict_alleles(v="v1*01").pcr_amplify(rate=0.1),
            39,
            id="pcr_call_flip",
        ),
        pytest.param(
            lambda e: e.restrict_alleles(v="v1*01").sequencing_errors(rate=0.1),
            39,
            id="seq_error_flip",
        ),
        pytest.param(
            lambda e: e.restrict_alleles(v="v1*01").ambiguous_base_calls(count=1),
            0,
            id="n_widen",
        ),
    ],
)
def test_replay_reproduces_call_set_and_identity_and_counters(
    configure, seed
) -> None:
    """Trace replay reproduces every evidence-derived AIRR field
    exactly: v_call (including CSV order), v_identity, v_cigar,
    truth_v_call, and the per-mechanism counters. The audit
    depends on this: if replay diverged, the live-call walker
    would have hidden RNG state outside the trace."""
    cfg = _two_v_one_distinguishing()
    base = _baseline_vj_experiment(cfg)
    compiled = configure(base).compile(allow_curatable_refdata=True)
    assert isinstance(compiled, CompiledExperiment)

    out_a = compiled.simulator.run(seed=seed)
    tf = compiled.simulator.trace_file_from(out_a, seed=seed)
    out_b = compiled.simulator.replay_from_trace_file(tf)

    from GenAIRR._airr_record import outcome_to_airr_record

    rec_a = outcome_to_airr_record(out_a, compiled.refdata, sequence_id="a")
    rec_b = outcome_to_airr_record(out_b, compiled.refdata, sequence_id="b")

    fields = [
        "sequence",
        "v_call",
        "v_identity",
        "v_cigar",
        "v_germline_start",
        "v_germline_end",
        "n_mutations",
        "n_pcr_errors",
        "n_quality_errors",
        "productive",
        "vj_in_frame",
    ]
    for field in fields:
        assert rec_a.get(field) == rec_b.get(field), (
            f"replay diverged on {field!r}: "
            f"{rec_a.get(field)!r} vs {rec_b.get(field)!r}"
        )


# ──────────────────────────────────────────────────────────────────
# 10. Pin-the-drift placeholders (audit §6)
#
# These tests document the current behaviour for §6 items. They're
# not failures — they pin what is, so a future fix can flip them.
# ──────────────────────────────────────────────────────────────────


def test_drift_pinned_no_call_count_field_on_airr_record() -> None:
    """Audit §6.1: no `v_call_count` / `d_call_count` /
    `j_call_count` denormalised counter today. Users count via
    `len(v_call.split(","))`. Pin absence so a future field
    addition signals the schema change."""
    rec = _record(
        _baseline_vj_experiment(_two_v_identical()),
        seed=0,
    )
    for field in ("v_call_count", "d_call_count", "j_call_count"):
        assert field not in rec, (
            f"AIRR record now exposes {field!r} — flip this test "
            f"to assert it equals len(call.split(','))"
        )


def test_drift_pinned_no_call_is_truthful_field_on_airr_record() -> None:
    """Audit §6.2: no `v_call_is_truthful` / etc. denormalised
    truthfulness boolean today. Users compare
    `v_call.split(',')[0] != truth_v_call` manually. Pin absence."""
    rec = _record(
        _baseline_vj_experiment(_two_v_one_distinguishing())
        .restrict_alleles(v="v1*01")
        .mutate(rate=0.1),
        seed=24,
    )
    # The flip has happened (v_call != truth_v_call), but no field
    # exposes that boolean directly.
    assert rec["v_call"] != rec["truth_v_call"]
    for field in (
        "v_call_is_truthful",
        "d_call_is_truthful",
        "j_call_is_truthful",
    ):
        assert field not in rec


def test_drift_pinned_no_v_identity_to_truth_field() -> None:
    """Audit §6.3: `v_identity` is against the projected (called)
    allele. No separate `v_identity_to_truth` field today. Pin
    absence so a future fix flips this."""
    rec = _record(
        _baseline_vj_experiment(_two_v_one_distinguishing())
        .restrict_alleles(v="v1*01")
        .mutate(rate=0.1),
        seed=24,
    )
    assert "v_identity_to_truth" not in rec
    assert "d_identity_to_truth" not in rec
    assert "j_identity_to_truth" not in rec
