"""Performance budget tests — regression guards.

Companion to
[docs/performance_baseline.md](../docs/performance_baseline.md).
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


pytestmark = pytest.mark.performance


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


def _vdj() -> "ge.RefDataConfig":
    cfg = ge.RefDataConfig.vdj()
    cfg.add_v_allele(
        "v1*01", "v1",
        b"GAGGTGCAGCTGGTGGAGTCTGGGGGAGGCTTGGTACAGCCTGGGGGGTCCCTGAGACTCTCCTGTGCAGCCTCTGGATTCACCTTC",
        anchor=78,
    )
    cfg.add_d_allele("d1*01", "d1", b"GGGTATAGCAGCAGCTGGTAC")
    cfg.add_j_allele(
        "j1*01", "j1",
        b"TACTACTACGGTATGGACGTCTGGGGCCAAGGGACCACGGTCACCGTCTCCTCAG",
        anchor=8,
    )
    return cfg


def _vj_recombine() -> "ga.Experiment":
    return (
        ga.Experiment.on(_vj())
        .recombine(np1_lengths=[(L, 1.0) for L in range(20)])
        .trim(
            v_3=[(L, 1.0) for L in range(5)],
            j_5=[(L, 1.0) for L in range(5)],
        )
    )


def _vdj_recombine() -> "ga.Experiment":
    return (
        ga.Experiment.on(_vdj())
        .recombine(
            np1_lengths=[(L, 1.0) for L in range(20)],
            np2_lengths=[(L, 1.0) for L in range(20)],
        )
        .trim(
            v_3=[(L, 1.0) for L in range(5)],
            d_5=[(L, 1.0) for L in range(5)],
            d_3=[(L, 1.0) for L in range(5)],
            j_5=[(L, 1.0) for L in range(5)],
        )
    )


def _vdj_full_stack() -> "ga.Experiment":
    return (
        _vdj_recombine()
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
    compiled = exp.compile()
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
    Baseline ~0.02 ms/rec; budget 1.0s for 5000 records (~10×)."""
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
    Baseline ~0.036 ms/rec; budget 1.5s for 5000 records (~8×)."""
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
    Baseline ~0.235 ms/rec; budget 2.0s for 1000 records (~8.5×).

    Regressions in NP-base admit mask, S5F per-site filtering,
    or indel DP retry budget would land here first."""
    exp = _vdj_full_stack().productive_only()
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
    Baseline ~0.058 ms/rec; budget 1.5s for 3000 records (~8.6×)."""
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
    Baseline ~0.122 ms/rec; budget 2.0s for 2000 records (~8.2×).

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
    Exercises the indel DP tuple sampler + the
    POST_EVENT_RETRY_BUDGET=16 fallback path.
    Baseline ~0.124 ms/rec; budget 1.0s for 1000 records (~8.1×).

    Regression that makes contract-rejection retries more frequent
    (or each retry more expensive) would land here."""
    exp = (
        _vdj_recombine()
        .polymerase_indels(count=(0, 5), insertion_prob=0.5)
        .productive_only()
    )
    wall = _measure(exp, n=1000, label="W6")
    assert wall < 1.0, (
        f"W6 (VDJ indel + productive) took {wall*1000:.0f}ms for 1000 records, "
        f"budget 1000ms."
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
    exp = _vdj_full_stack().productive_only()
    compiled = exp.compile()

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

    exp = _vdj_full_stack().productive_only()
    compiled = exp.compile()
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
