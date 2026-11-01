"""Performance budget tests — regression guards.

Companion to
[`audit-docs/performance_baseline.md`](audit-docs/performance_baseline.md).
Each test runs a representative workload for a fixed N records,
asserts the wall time is under the budget set at ~5–10× the
developer-machine baseline. Goal: catch 10× regressions without
flaking on slow CI.

These are NOT microbenchmarks — they're amortised wall-time
bounds. Per-record noise averages out over the workload's N.

Tests are marked ``@pytest.mark.performance`` so a CI pipeline
that wants to skip them can run ``pytest -m "not performance"``.

If a budget fails, **profile first** (py-spy / samply per
``MEMORY.md::performance``), find the regression's root cause,
and only update the budget when the slowdown is intentional.
"""
from __future__ import annotations

import time

import pytest

import GenAIRR as ga
from GenAIRR import _engine as ge
from tests.support.fixtures.budget_experiments import _vdj_recombine


pytestmark = [pytest.mark.performance, pytest.mark.composition]


# ──────────────────────────────────────────────────────────────────
# Reference refdata: human-IGH-ish single allele per segment
# ──────────────────────────────────────────────────────────────────


def _vj() -> "ge.RefDataConfig":
    cfg = ge.RefDataConfig.vj()
    cfg.add_v_allele(
        "v1*01", "v1",
        b"GAGGTGCAGCTGGTGGAGTCTGGGGGAGGCTTGGTACAGCCTGGGGGGTCCCTGAGACTCTCCTGTGCAGCCTCTGGATTCACCTTC",
        anchor=78,
    )
    cfg.add_j_allele(
        "j1*01", "j1",
        b"TACTACTACGGTATGGACGTCTGGGGCCAAGGGACCACGGTCACCGTCTCCTCAG",
        anchor=8,
    )
    return cfg


def _vj_recombine() -> "ga.Experiment":
    return (
        ga.Experiment.on(_vj()).allow_curatable_refdata()
        .recombine(np1_lengths=[(L, 1.0) for L in range(20)])
        .trim(
            v_3=[(L, 1.0) for L in range(5)],
            j_5=[(L, 1.0) for L in range(5)],
        )
    )


def _vdj_full_stack(*, canonical=False) -> "ga.Experiment":
    return (
        _vdj_recombine(canonical=canonical)
        .mutate(rate=0.01)
        .pcr_amplify(rate=0.001)
        .sequencing_errors(rate=0.001)
    )


# ──────────────────────────────────────────────────────────────────
# Common harness
# ──────────────────────────────────────────────────────────────────


WARMUP = 50  # records discarded before timing


def _time_simulation_loop(compiled, n: int, seed_offset: int = 0) -> float:
    """Run `n` simulations with seeds ``[seed_offset, seed_offset + n)``.
    Returns wall-clock seconds."""
    t0 = time.perf_counter()
    for i in range(n):
        compiled.simulator.run(seed=seed_offset + i)
    return time.perf_counter() - t0


def _measure(exp, n: int, label: str) -> float:
    """Compile, warm up, time `n` records. Returns wall seconds."""
    compiled = exp.compile(allow_curatable_refdata=True)
    assert hasattr(compiled, "simulator"), (
        f"{label}: clonal CompiledExperiment not expected here"
    )
    # Warmup absorbs first-call effects.
    _time_simulation_loop(compiled, WARMUP, seed_offset=0)
    return _time_simulation_loop(compiled, n, seed_offset=WARMUP)


# ──────────────────────────────────────────────────────────────────
# W1: VJ recombine only
# ──────────────────────────────────────────────────────────────────


def test_w1_vj_recombine_only_budget() -> None:
    """W1: VJ recombine (no corruption). Cheapest workload.
    Measured 0.025 ms/rec on 2026-10-29; budget 1.0s for 5000 records,
    so 7.9× headroom. The recorded baseline was 0.02 ms/rec (~10×)."""
    wall = _measure(_vj_recombine(), n=5000, label="W1")
    assert wall < 1.0, (
        f"W1 (VJ recombine only) took {wall*1000:.0f}ms for 5000 records, "
        f"budget 1000ms. See docs/performance_baseline.md."
    )


# ──────────────────────────────────────────────────────────────────
# W2: VDJ recombine only
# ──────────────────────────────────────────────────────────────────


def test_w2_vdj_recombine_only_budget() -> None:
    """W2: VDJ recombine. Adds D + NP2 vs W1.
    Measured 0.041 ms/rec on 2026-10-29; budget 1.5s for 5000 records,
    so 7.2× headroom. The recorded baseline was 0.036 ms/rec (~8×)."""
    wall = _measure(_vdj_recombine(), n=5000, label="W2")
    assert wall < 1.5, (
        f"W2 (VDJ recombine only) took {wall*1000:.0f}ms for 5000 records, "
        f"budget 1500ms."
    )


# ──────────────────────────────────────────────────────────────────
# W3: VDJ productive full stack
# ──────────────────────────────────────────────────────────────────


def test_w3_vdj_productive_full_stack_budget() -> None:
    """W3: VDJ + SHM + PCR + sequencing + productive_only.
    The contract-narrowing workload — most expensive per-pass.
    Measured 0.759 ms/rec on 2026-11-01, down from 0.839 on the same box when
    the N-length feasibility walk was quadratic. Still 3.2× the recorded
    0.235 ms/rec, so the drift this docstring has recorded since 2026-10-29 is
    only partly accounted for. The budget value is deliberately unchanged,
    because the drift is the finding.

    **This budget barely sees the cost it is named for.** Its refdata is one
    allele per segment with a narrow N-length support, and the quadratic walk
    that cost packaged TCRB 2.7× showed up here as 10%. A constrained budget
    on a packaged reference is what covers that, and it is W9 below.

    Regressions in NP-base admit mask, S5F per-site filtering,
    or indel DP retry budget would land here first."""
    exp = _vdj_full_stack(canonical=True).productive_only()
    wall = _measure(exp, n=1000, label="W3")
    assert wall < 2.0, (
        f"W3 (VDJ productive full stack) took {wall*1000:.0f}ms for 1000 records, "
        f"budget 2000ms. Profile with py-spy --native to find the regression."
    )


# ──────────────────────────────────────────────────────────────────
# W4: VDJ non-productive full stack
# ──────────────────────────────────────────────────────────────────


def test_w4_vdj_non_productive_full_stack_budget() -> None:
    """W4: same as W3 minus productive_only(). Isolates the
    contract-narrowing cost (W3 ≈ 4× W4 in baseline).
    Measured 0.064 ms/rec on 2026-10-29; budget 1.5s for 3000 records,
    so 7.9× headroom. The recorded baseline was 0.058 ms/rec (~8.6×)."""
    wall = _measure(_vdj_full_stack(), n=3000, label="W4")
    assert wall < 1.5, (
        f"W4 (VDJ non-productive full stack) took {wall*1000:.0f}ms "
        f"for 3000 records, budget 1500ms."
    )


# ──────────────────────────────────────────────────────────────────
# W5: VDJ high SHM (rate=0.10)
# ──────────────────────────────────────────────────────────────────


def test_w5_vdj_high_shm_budget() -> None:
    """W5: SHM at rate=0.10 (~10 mutations per ~160-bp sequence).
    Stresses the S5F candidate enumeration hot path.
    Measured 0.139 ms/rec on 2026-10-29; budget 2.0s for 2000 records,
    so 7.2× headroom. The recorded baseline was 0.122 ms/rec (~8.2×).

    Regression in S5F enumeration (e.g. re-scanning the full pool
    on every mutation instead of using cached mutability) would
    land here first."""
    exp = _vdj_recombine().mutate(rate=0.10)
    wall = _measure(exp, n=2000, label="W5")
    assert wall < 2.0, (
        f"W5 (VDJ high SHM) took {wall*1000:.0f}ms for 2000 records, "
        f"budget 2000ms."
    )


# ──────────────────────────────────────────────────────────────────
# W6: VDJ indel-heavy + productive
# ──────────────────────────────────────────────────────────────────


def test_w6_vdj_indel_heavy_productive_budget() -> None:
    """W6: polymerase_indels(count=(0,5)) + productive_only.
    Exercises complete natural-history proposals and bounded exact search
    with independent restart probes after the initial rejection phase.
    Measured 0.606 ms/rec on 2026-10-29 -- 4.9× the recorded 0.124 ms/rec.
    Against the original 1.0s budget that left 1.6× headroom rather than the
    8.1× recorded here, and a macOS CI runner measured 1207ms and failed it
    while eleven other matrix legs passed. The budget is raised to 2.0s: 3.3×
    headroom locally, ~1.7× on the slowest runner observed, so it still fails
    on a further regression of that size instead of failing on runner
    variance. See the note at the end of this module.

    Regression that makes contract-rejection retries more frequent
    (or each retry more expensive) would land here."""
    exp = (
        _vdj_recombine(canonical=True)
        .polymerase_indels(count=(0, 5), insertion_prob=0.5)
        .productive_only()
    )
    wall = _measure(exp, n=1000, label="W6")
    assert wall < 2.0, (
        f"W6 (VDJ indel + productive) took {wall*1000:.0f}ms for 1000 records, "
        f"budget 2000ms."
    )


# ──────────────────────────────────────────────────────────────────
# W7: Replay-vs-fresh ratio
# ──────────────────────────────────────────────────────────────────


def test_w7_replay_is_substantially_faster_than_fresh() -> None:
    """W7: trace replay must remain substantially faster than
    fresh sampling. Baseline ratio ~0.27 (replay ≈ 3.6× faster).

    Budget: replay wall ≤ 0.8 × fresh wall — generous enough to
    survive CI noise, tight enough to catch a regression that
    makes replay re-run the contract filter on every recorded
    value (which would collapse the ratio toward 1.0).

    Run on the W3 plan (productive full stack) because that's
    where replay's contract-skip win is largest."""
    exp = _vdj_full_stack(canonical=True).productive_only()
    compiled = exp.compile(allow_curatable_refdata=True)

    N = 500  # smaller than W3 to keep total test time low
    # Warmup
    _time_simulation_loop(compiled, WARMUP)

    # Fresh wall.
    fresh_wall = _time_simulation_loop(compiled, N, seed_offset=WARMUP)

    # Build trace files from a fresh batch.
    trace_files = []
    for seed in range(WARMUP, WARMUP + N):
        out = compiled.simulator.run(seed=seed)
        trace_files.append(compiled.simulator.trace_file_from(out, seed=seed))

    # Replay wall.
    t0 = time.perf_counter()
    for tf in trace_files:
        compiled.simulator.replay_from_trace_file(tf)
    replay_wall = time.perf_counter() - t0

    ratio = replay_wall / fresh_wall
    assert ratio < 0.8, (
        f"W7 replay/fresh ratio {ratio:.3f} ≥ 0.8 — replay is no longer "
        f"substantially faster than fresh sampling. "
        f"(fresh={fresh_wall*1000:.0f}ms, replay={replay_wall*1000:.0f}ms, "
        f"baseline ratio ~0.27)"
    )


# ──────────────────────────────────────────────────────────────────
# W8: AIRR projection overhead (sanity)
# ──────────────────────────────────────────────────────────────────


def test_w8_airr_projection_is_negligible_vs_simulation() -> None:
    """W8: AIRR projection (outcome → dict) must remain
    negligible vs simulation cost. Baseline ~0.02 ms/rec
    (~10% of W3 simulation cost). Budget: projection wall ≤
    0.5 × simulation wall on a productive full-stack outcome.

    A regression that bloats the projection (e.g. recomputing
    the junction for every dict access) would land here."""
    from GenAIRR._airr_record import outcome_to_airr_record

    exp = _vdj_full_stack(canonical=True).productive_only()
    compiled = exp.compile(allow_curatable_refdata=True)
    _time_simulation_loop(compiled, WARMUP)

    N = 500
    # Generate outcomes (and time fresh-sim wall as the reference).
    outcomes = []
    t0 = time.perf_counter()
    for seed in range(WARMUP, WARMUP + N):
        outcomes.append(compiled.simulator.run(seed=seed))
    sim_wall = time.perf_counter() - t0

    # Project all outcomes.
    t0 = time.perf_counter()
    for out in outcomes:
        outcome_to_airr_record(out, compiled.refdata, sequence_id="x")
    proj_wall = time.perf_counter() - t0

    ratio = proj_wall / sim_wall
    assert ratio < 0.5, (
        f"AIRR projection wall is {proj_wall*1000:.0f}ms vs simulation "
        f"{sim_wall*1000:.0f}ms (ratio {ratio:.2f}); projection should be "
        f"negligible relative to simulation (baseline ~10%)."
    )


# ──────────────────────────────────────────────────────────────────
# W9: the constraint cost on a packaged reference
# ──────────────────────────────────────────────────────────────────


def test_w9_productive_only_on_a_packaged_reference() -> None:
    """W9: `productive_only()` on packaged human_tcrb, as a ratio.

    Every budget above builds its own refdata: one allele per segment, a
    narrow N-length support. That keeps them fast and stable and it also
    made them blind. A quadratic walk in the N-length feasibility solver
    cost packaged TCRB 3.1x on this very harness, and W3 -- the budget
    named for the contract-narrowing workload -- moved 10%.

    **Why a ratio and not just a wall.** The wall here is 0.36s for 300
    records on the development box, and the quadratic version was 1.13s. A
    limit loose enough not to flake on a slow CI runner is above 1.13s, so
    it would have reported "pass" for both. The ratio against the same
    locus unconstrained on the same box cancels the machine: it was 26x
    before the fix and 7-10x after, which is a gap no runner variance
    closes. That is the guard. The wall is kept beside it as a sanity
    bound, which is all it can honestly be.

    TCRB rather than IGH because it is where the cost showed: IGH sat at
    5-8x throughout. Measured 2026-11-01 on a box at load 10-12, so these
    are upper bounds and the ratio is the claim.
    """
    locus = "human_tcrb"
    plain = _measure(ga.Experiment.on(locus).recombine(), n=1200, label="W9-plain")
    constrained = _measure(
        ga.Experiment.on(locus).recombine().productive_only(), n=300, label="W9"
    )
    assert constrained < 1.5, (
        f"W9 (productive_only on packaged {locus}) took {constrained*1000:.0f}ms "
        f"for 300 records, budget 1500ms."
    )
    ratio = (constrained / 300) / (plain / 1200)
    assert ratio < 18.0, (
        f"W9: productive_only() on packaged {locus} costs {ratio:.1f}x the "
        f"unconstrained run ({constrained/0.3:.2f} vs {plain/1.2:.2f} ms/rec), "
        f"over the 18x limit. Measured 7-10x on 2026-11-01; 26x is what the "
        f"quadratic N-length feasibility walk cost. Profile with "
        f"py-spy --native on a build made per engine_rs/Cargo.toml's recipe -- "
        f"the installed .so is stripped."
    )


# ──────────────────────────────────────────────────────────────────
# What the recorded baselines turned out to say
# ──────────────────────────────────────────────────────────────────
#
# Re-measured on 2026-10-29, after a macOS CI runner failed W6 at 1207ms
# against its 1000ms budget while eleven other matrix legs passed:
#
#     budget  per-rec now  recorded  drift   uses productive_only()
#     W1        0.025ms     0.020ms   1.3x   no
#     W2        0.041ms     0.036ms   1.2x   no
#     W4        0.064ms     0.058ms   1.1x   no
#     W5        0.139ms     0.122ms   1.1x   no
#     W3        0.860ms     0.235ms   3.7x   YES
#     W6        0.606ms     0.124ms   4.9x   YES
#
# Four of six sit within 1.1-1.3x of what was recorded, which is the spread
# you expect from different hardware. The two that drifted are exactly the
# two that narrow for productivity, and nothing else moved -- so this is the
# productive_only() path having become 3.7-4.9x more expensive since these
# baselines were taken, not a machine difference.
#
# WHY it became more expensive is no longer open, and nothing regressed.
# f066b44 (2026-10-20) added the exact productive-completion solver in
# engine_rs/src/feasibility/p_productive.rs -- 669 lines, created whole, no
# prior version -- and added CodingSequenceIntegrity as a fifth contract in
# the productive() bundle. Both are correctness work, and both are reached
# only under productive_only(), which is why these two budgets moved and the
# other four did not. The solver runs even with no mutation pass configured:
# measured at 810 p_lengths calls per record on W6. The contract is reached
# only when a mutation pass runs, so it affects W3 (`.mutate(rate=0.01)`)
# and not W6.
#
# That commit measured its own cost on the day, in
# audit-docs/productivity_performance_2026_10_20.json: W3 0.985 ms/rec and
# W6 0.792 ms/rec. So the increase was known, intentional and recorded when
# it landed -- what was never updated is the baseline these budgets compare
# against, which predates the solver by five months. The 2026-10-29 figures
# above (0.860 and 0.606) are BELOW the October measurement, so the path has
# got faster since the work landed, not slower.
#
# Treat the May numbers as what they are: measurements of an engine that did
# not yet do this. They are not a target to get back to.
#
# These budgets were written to catch a change of this size. They did: the
# numbers above are their output. What they could not do is announce it,
# because the increase never crossed a limit on the runners in use until a
# slower one appeared.
# A budget whose recorded headroom is 8x and whose real headroom is 1.6x
# reports "pass" either way, which is why the figures here are now measured
# values with a date rather than a remembered baseline.
#
# W6's limit is raised so CI stops failing on runner variance. W3's is not:
# it still passes, and leaving it at 2.3x headroom keeps the drift visible
# where raising it would hide it.
#
# ── 2026-11-01: part of the cause found, and the blind spot closed ──
#
# The note above handed the investigation to its own branch. One finding came
# back: `np_exists` in the solver that commit added ran its N-length automaton
# once per candidate length, where the automaton for length L is a strict
# prefix of the one for L+1. With the 51-value support every packaged locus
# declares, that was 1275 position-steps per call instead of 50.
#
# Fixing it bought 2.7x on packaged human_tcrb and 1.2x on human_tcrd. It
# bought W3 10%, which is the point worth keeping: **these budgets could not
# see it.** Every fixture above builds its own refdata with one allele per
# segment and a narrow N-length support, so the term that was quadratic
# barely appeared. The budgets did not fail to announce a 2.7x cost because
# their limits were loose -- they never measured it.
#
# W9 above is the answer: packaged refdata, and a ratio against the same
# locus unconstrained rather than a wall, because a wall loose enough for CI
# sat above the quadratic version too. It reports 26x before the fix and
# 7-10x after. It is also the first budget here that would have caught this
# class at all.
#
# W3's remaining gap to its May baseline (3.2x, down from 3.7x) is still not
# fully accounted for, and the May numbers remain what they were: a
# measurement of an engine that did not do this work.
