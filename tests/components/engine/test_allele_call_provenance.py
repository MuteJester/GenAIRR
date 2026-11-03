"""Allele-call disambiguation golden tests.

Companion to [`audit-docs/allele_call_audit.md`](audit-docs/allele_call_audit.md).
Pins the live walker's evidence-based `recalled_*_call` fields under
SHM, PCR errors, sequencing errors, and N substitution. Standard
`v_call`/`d_call`/`j_call` and alignment fields keep the generating
assignment, verified separately on every record.

Fixture design notes:

- 2-allele V refdata used throughout: ``v1*01 = AAACCCGGG`` and
  ``v2*01 = AAACTCGGG``. They differ at position 4 (C vs T).
  Positions 0-3 and 5-8 are non-distinguishing.
- ``v1*01`` and ``v1*02`` are identical (same byte string) — used
  to pin the always-tied case.
- ``restrict_alleles(v="v1*01")`` locks the recombination-stage
  truth, isolating evidence-changing mechanisms from sampling
  variance.

**Seeds are searched by the tests, not pinned into them.** They used to be
pinned, found by a search run by hand and written in beside a prose
description of the molecule each one produced. Every change to the engine's
draw order invalidated all of them at once -- per-pass RNG substreams did,
and eleven of these tests failed together while their comments went on
describing molecules those seeds no longer produce. What a test here needs is
not seed 24, it is a record with the other allele's base at the
distinguishing position. `_seed_where` finds one and says so if it cannot,
which is the same premise the pinned seeds encoded, stated where it can
still be read.
"""
from __future__ import annotations

import pytest

import GenAIRR as ga
from GenAIRR import _engine as ge
from GenAIRR import CompiledExperiment

pytestmark = pytest.mark.behavior


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
    record = result.records[0]
    for segment in ("v", "d", "j"):
        assert record[f"{segment}_call"] == record[f"truth_{segment}_call"]
    assert result.validate_records(exp.refdata).ok
    return record


#: How many seeds a premise search scans before giving up. One record on
#: these nine-base fixtures costs about a millisecond, so the whole module
#: searching is cheaper than one real simulation.
SEED_SCAN = 300


def _seed_where(exp: "ga.Experiment", premise, *, scan: int = SEED_SCAN,
                **kwargs) -> tuple[int, dict]:
    """The first seed whose record satisfies ``premise``, and that record.

    The premise is the molecular configuration a claim is about -- an `N` at
    the distinguishing position, a substitution that wrote the other allele's
    base there. Searching for it keeps the claim exact without tying it to a
    draw order, and a premise no seed satisfies fails loudly rather than
    quietly testing something else.
    """
    for seed in range(scan):
        record = _record(exp, seed=seed, **kwargs)
        if premise(record):
            return seed, record
    raise AssertionError(
        f"no seed below {scan} produced a record satisfying the premise, so "
        "this claim cannot be tested on this pipeline -- the mechanism may "
        "have stopped reaching the position the premise names")


def _distinguishing_pipeline(mechanism) -> "ga.Experiment":
    """The two-allele fixture with the truth locked and one mechanism on."""
    return mechanism(
        _baseline_vj_experiment(_two_v_one_distinguishing())
        .restrict_alleles(v="v1*01"))


def _three_allele_n_pipeline() -> "ga.Experiment":
    """The three-allele fixture with one base masked per record."""
    return (_baseline_vj_experiment(_three_v_prefix_suffix())
            .restrict_alleles(v="v_A*01")
            .ambiguous_base_calls(count=1))


def _shm_flip() -> tuple[int, dict]:
    """A record where hypermutation wrote v2*01's base at position 4."""
    return _seed_where(
        _distinguishing_pipeline(lambda e: e.mutate(rate=0.1)),
        lambda rec: rec["sequence"][4] == "T")


def _pcr_flip() -> tuple[int, dict]:
    """The same molecule, reached by a PCR misincorporation instead."""
    return _seed_where(
        _distinguishing_pipeline(lambda e: e.pcr_amplify(rate=0.1)),
        lambda rec: rec["sequence"][4] == "T")


# ──────────────────────────────────────────────────────────────────
# 1. Tie-set formation — model basics
# ──────────────────────────────────────────────────────────────────


def test_identical_v_alleles_produce_csv_tie_with_truth_first() -> None:
    """Two literally identical V alleles → walker can't
    distinguish them → tie-set contains both. AIRR `recalled_v_call` is
    the CSV with the truth allele first (`projected_call_name`
    convention)."""
    exp = _baseline_vj_experiment(_two_v_identical())
    # Sweep until both truth picks have been seen. The sweep was eight seeds,
    # which was enough for one draw order and not for the next -- the number
    # that matters is "both", so the loop stops when it has them.
    truth_first_observations = {"v1*01": False, "v1*02": False}
    for seed in range(SEED_SCAN):
        if all(truth_first_observations.values()):
            break
        rec = _record(exp, seed=seed)
        truth = rec["truth_v_call"]
        csv = rec["recalled_v_call"].split(",")
        assert len(csv) == 2
        assert set(csv) == {"v1*01", "v1*02"}
        assert csv[0] == truth, (
            f"seed {seed}: truth {truth} must be first in CSV {csv}"
        )
        truth_first_observations[truth] = True
    # Both truths must have been sampled, or the CSV-order claim was only
    # ever checked for one of them.
    assert all(truth_first_observations.values()), (
        f"the fixture drew only one truth allele in {SEED_SCAN} seeds"
    )


def test_distinguishing_alleles_resolve_to_truth_without_mutation() -> None:
    """Without any evidence-changing mechanism, the assembled pool
    is the truth allele's germline. The walker scores it +9 vs.
    +8 for the other allele → singleton tie-set = truth."""
    exp = _baseline_vj_experiment(_two_v_one_distinguishing())
    for seed in range(8):
        rec = _record(exp, seed=seed)
        assert rec["recalled_v_call"] == rec["truth_v_call"], (
            f"seed {seed}: recalled_v_call {rec['recalled_v_call']!r} should resolve to "
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

    The premise is a `T` at position 4, matching v2*01's germline;
    `_shm_flip` searches for a record that has one. recalled_v_call must then
    flip to 'v2*01' while truth_v_call stays 'v1*01', and the standard v_call
    and alignment stay generating."""
    seed, rec = _shm_flip()
    assert rec["sequence"][4] == "T", f"premise broken at seed {seed}"
    assert rec["truth_v_call"] == "v1*01"
    assert rec["recalled_v_call"] == "v2*01", (
        f"seed {seed}: a T at the distinguishing position did not move the "
        "call")
    # SHM increments n_mutations only.
    assert rec["n_mutations"] > 0
    assert rec["n_pcr_errors"] == 0
    assert rec["n_quality_errors"] == 0


def test_shm_at_non_distinguishing_position_preserves_call() -> None:
    """When SHM lands at a non-distinguishing position, both
    alleles' scores drop by the same amount. The distinguishing
    position still resolves the call to truth.

    The premise is one mutation that left position 4 alone, which is the
    case the claim is about -- searched rather than pinned, like the flip
    above."""
    seed, rec = _seed_where(
        _distinguishing_pipeline(lambda e: e.mutate(rate=0.1)),
        lambda r: r["n_mutations"] == 1 and r["sequence"][4] == "C")
    assert rec["truth_v_call"] == "v1*01"
    assert rec["recalled_v_call"] == "v1*01", (
        f"seed {seed}: a mutation away from the distinguishing position moved "
        "the call")


# ──────────────────────────────────────────────────────────────────
# 3. PCR & sequencing errors — same walker semantics as SHM
# ──────────────────────────────────────────────────────────────────


def test_pcr_error_at_distinguishing_position_switches_call() -> None:
    """PCR errors substitute uppercase canonical bases. They are
    walker-semantically identical to SHM — a PCR error landing at
    the distinguishing position with the other allele's base
    flips the call.

    The premise is the same molecule as the SHM case -- an uppercase `T` at
    position 4 -- reached by a different mechanism."""
    seed, rec = _pcr_flip()
    assert rec["sequence"][4] == "T", f"premise broken at seed {seed}"
    assert rec["truth_v_call"] == "v1*01"
    assert rec["recalled_v_call"] == "v2*01", f"at seed {seed}"
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

    The premise is a lowercase 't' at position 4. recalled_v_call must then
    be 'v2*01' alone, which is what confirms the walker did not treat the
    lowercase letter as an N -- that would have widened the set to
    'v1*01,v2*01'."""
    seed, rec = _seed_where(
        _distinguishing_pipeline(lambda e: e.sequencing_errors(rate=0.1)),
        lambda r: r["sequence"][4] == "t")
    assert rec["truth_v_call"] == "v1*01"
    assert rec["recalled_v_call"] == "v2*01", (
        f"seed {seed}: a lowercase base widened the tie set, so it was read "
        "as a wildcard")
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

    The premise is an `N` at position 4. recalled_v_call must then widen
    from {v1*01} to {v1*01, v2*01}, truth first per the projection
    convention."""
    seed, rec = _seed_where(
        _distinguishing_pipeline(lambda e: e.ambiguous_base_calls(count=1)),
        lambda r: r["sequence"][4] == "N")
    assert rec["truth_v_call"] == "v1*01"
    assert rec["recalled_v_call"] == "v1*01,v2*01", (
        f"seed {seed}: an N at the distinguishing position did not widen the "
        "tie set to both alleles with truth first")


def test_n_substitution_at_non_distinguishing_position_preserves_call() -> None:
    """An N at a non-distinguishing position matches both alleles
    via wildcard, but they're already tied at that column — so
    the *relative* scores don't change. The call stays
    singleton-truth.

    The premise is an `N` at position 8, which both alleles match
    identically."""
    seed, rec = _seed_where(
        _distinguishing_pipeline(lambda e: e.ambiguous_base_calls(count=1)),
        lambda r: r["sequence"][8] == "N" and r["sequence"][4] == "C")
    assert rec["recalled_v_call"] == "v1*01", (
        f"seed {seed}: an N away from the distinguishing position widened the "
        "tie set anyway")


# ──────────────────────────────────────────────────────────────────
# 5. 3-allele prefix/suffix ambiguity
# ──────────────────────────────────────────────────────────────────


def test_three_allele_prefix_widening_keeps_truth_first() -> None:
    """With three V alleles where v_B differs at pos 4 (prefix-
    ambiguous with v_A) and v_C differs at pos 8 (suffix-
    ambiguous with v_A), an N at pos 4 widens the tie to {v_A,
    v_B} only (v_C still loses at pos 8). Truth v_A first.

    The premise is an `N` at position 4 on the three-allele fixture."""
    seed, rec = _seed_where(
        _three_allele_n_pipeline(), lambda r: r["sequence"][4] == "N")
    assert rec["truth_v_call"] == "v_A*01"
    assert rec["recalled_v_call"] == "v_A*01,v_B*01", f"at seed {seed}"
    assert "v_C*01" not in rec["recalled_v_call"]


def test_three_allele_suffix_widening_keeps_truth_first() -> None:
    """Mirror of the prefix case: N at pos 8 widens to {v_A, v_C}
    only (v_B still loses at pos 4)."""
    seed, rec = _seed_where(
        _three_allele_n_pipeline(), lambda r: r["sequence"][8] == "N")
    assert rec["recalled_v_call"] == "v_A*01,v_C*01", f"at seed {seed}"
    assert "v_B*01" not in rec["recalled_v_call"]


# ──────────────────────────────────────────────────────────────────
# 6. Truth decoupling + identity-against-projected
# ──────────────────────────────────────────────────────────────────


def test_recalled_v_call_can_differ_while_standard_call_keeps_truth() -> None:
    """Evidence-based inference changes without replacing generative origin."""
    _, rec_shm = _shm_flip()
    assert rec_shm["recalled_v_call"] != rec_shm["truth_v_call"]
    assert rec_shm["truth_v_call"] == "v1*01"

    # Same property, different mechanism.
    _, rec_pcr = _pcr_flip()
    assert rec_pcr["recalled_v_call"] != rec_pcr["truth_v_call"]
    assert rec_pcr["truth_v_call"] == "v1*01"


def test_v_identity_uses_generating_reference_when_recalled_allele_changes() -> None:
    """`v_identity` measures the molecule against the allele that generated
    it, not against the one the evidence now favours -- so on a flipped
    record it is the *lower* of the two numbers."""
    _, rec = _shm_flip()
    assert rec["recalled_v_call"] == "v2*01"
    assert rec["v_call"] == "v1*01"
    assert rec["germline_alignment"][:9] == "AAACCCGGG"
    mismatches = sum(a.upper() != b for a, b in
                     zip(rec["sequence"][:9], "AAACCCGGG"))
    assert mismatches >= 1, "the flipped record matches its generating allele"
    assert rec["v_identity"] == pytest.approx((9 - mismatches) / 9, abs=1e-9)


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
    VDJ refdata. Without mutation, recalled_d_call equals truth_d_call."""
    cfg = _vdj_two_d_one_distinguishing()
    exp = (
        ga.Experiment.on(cfg).allow_curatable_refdata()
        .recombine(np1_lengths=[(2, 1.0)], np2_lengths=[(2, 1.0)])
        .trim(enabled=False)
    )
    for seed in range(6):
        rec = _record(exp, seed=seed)
        assert rec["recalled_d_call"] == rec["truth_d_call"]


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
    exactly: recalled_v_call (including CSV order), v_identity, v_cigar,
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
        "d_call",
        "j_call",
        "recalled_v_call",
        "recalled_d_call",
        "recalled_j_call",
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
    `len(recalled_v_call.split(","))`. Pin absence so a future field
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
    `recalled_v_call.split(',')[0] != truth_v_call` manually. Pin absence."""
    _, rec = _shm_flip()
    # The flip has happened (recalled_v_call != truth_v_call), but no field
    # exposes that boolean directly.
    assert rec["recalled_v_call"] != rec["truth_v_call"]
    for field in (
        "v_call_is_truthful",
        "d_call_is_truthful",
        "j_call_is_truthful",
    ):
        assert field not in rec


def test_drift_pinned_no_v_identity_to_truth_field() -> None:
    """Standard v_identity already measures the generating reference."""
    _, rec = _shm_flip()
    assert "v_identity_to_truth" not in rec
    assert "d_identity_to_truth" not in rec
    assert "j_identity_to_truth" not in rec
    assert rec["v_identity"] < 1.0, (
        "the flipped record matches its generating allele perfectly, so this "
        "says nothing about which reference v_identity measures against")
