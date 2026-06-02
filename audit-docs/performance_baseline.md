# Performance baseline + budget audit

**Status:** baseline measured, budgets enforced as pytest assertions in
[`tests/test_performance_budgets.py`](../tests/test_performance_budgets.py).

This audit complements the correctness audits
([primer_trim_end_loss](primer_trim_end_loss_audit.md),
[indel_provenance](indel_provenance_audit.md),
[allele_call](allele_call_audit.md),
[productive_failure_mode](productive_failure_mode_audit.md),
[junction_call](junction_call_audit.md),
[distribution_invariant](distribution_invariant_audit.md))
by proving the engine **remains practical at scale** after every
validation, contract-narrowing, event-capture, and replay
mechanism added since the v3.0 architecture.

The goal is **not** a microbenchmark; it is a regression guard.
Budgets are set at ~5–10× the current observed mean so a 10×
regression (e.g. accidentally re-allocating per-record, removing
a fast path, doubling RNG draws) is caught while normal CI noise
is absorbed.

---

## 1. Methodology

- **Single-thread Python orchestrator** calling the compiled
  `_engine.abi3.so` via PyO3.
- **Warmup**: 50 records discarded before timing to absorb the
  first compile/JIT-ish effects of the extension import.
- **Timer**: `time.perf_counter()` around each run loop.
- **Workload N chosen** so the total wall time is ≥ 100ms per
  workload — small enough that the full audit runs in a few
  seconds, large enough that single-record noise averages out.
- **Released-profile binary** (`maturin develop --release`) is the
  reference build. Debug builds will be substantially slower; do
  not compare debug-build measurements to these budgets.

---

## 2. Representative workloads

Each workload exercises a slightly different cross-section of the
engine. Refdata is a single human-IGH-ish V allele (86bp,
anchor=78), single D allele (21bp), single J allele (54bp,
anchor=8) — large enough to be representative but small enough
that allele-pool lookups don't dominate.

| ID | Workload                                | What it stresses                                                          |
|----|------------------------------------------|---------------------------------------------------------------------------|
| W1 | VJ recombine only                       | Plan compile + allele sample + NP length + NP base + trim — minimum cost. |
| W2 | VDJ recombine only                      | W1 + D segment + NP2 length/base.                                          |
| W3 | VDJ productive full stack               | W2 + SHM + PCR + sequencing errors + `productive_only()` narrowing.        |
| W4 | VDJ non-productive full stack           | Same as W3 minus `productive_only()`. Isolates contract-narrowing cost.    |
| W5 | VDJ high SHM (rate=0.10)                | ~10 mutations/record on the 161-base sequence; S5F enumeration hot path.   |
| W6 | VDJ indel-heavy + productive            | `polymerase_indels(count=(0,5))` + productive — exercises the DP tuple sampler. |
| W7 | Replay-vs-fresh ratio (W3 plan)         | Trace-injected replay vs fresh sampling for the same plan + seeds.         |

---

## 3. Baseline measurements (developer machine, 2026-05-27)

Reference machine: Linux 6.17, single-thread Python 3.x, release
build. Numbers below are **mean per-record latency over N records
post-warmup**; CI machines may be 2–3× slower.

| Workload | N records | Wall (ms) | mean ms/rec | rec/s     | trace len (avg) | event count (avg) |
|----------|-----------|-----------|-------------|-----------|-----------------|-------------------|
| W1       | 5,000     | 102       | 0.020       | 48,810    | 14              | 154               |
| W2       | 5,000     | 179       | 0.036       | 28,000    | 28              | 186               |
| W3       | 1,000     | 234       | 0.235       | 4,264     | 35              | 189               |
| W4       | 3,000     | 174       | 0.058       | 17,262    | 35              | 189               |
| W5       | 2,000     | 244       | 0.122       | 8,196     | 63              | 204               |
| W6       | 1,000     | 124       | 0.124       | 8,086     | 35              | 189               |
| W7 fresh | 1,000     | 228       | 0.228       | 4,381     | 35              | 189               |
| W7 replay| 1,000     | 62        | 0.062       | 16,077    | 35              | 189               |

Bonus: **AIRR projection** (`outcome_to_airr_record`) mean
0.02 ms/record. Negligible relative to simulation.

### 3.1 Observations from the baseline

- **W3 vs W4**: productive narrowing adds ~4× latency (0.235 vs
  0.058 ms/rec). Most of the cost is the contract-aware sampler
  paths (NP-base mask, S5F site filtering, indel DP).
- **W6 ≈ W5 ≈ 0.12 ms/rec**: SHM at rate=0.10 and indel-heavy
  cost about the same per record, despite the indel DP being
  more algorithmically complex. The S5F site enumeration is the
  dominant per-mutation cost.
- **W7 replay is ~3.7× faster than fresh** (0.062 vs 0.228 ms).
  Trace replay skips most RNG draws and per-event contract
  filtering — it's reading recorded proposals and validating
  against current support.
- **Trace length scales roughly with sample slot count**: ~14 in
  VJ recombine, ~28 in VDJ recombine, ~35 in VDJ + corruption.
  SHM at rate=0.10 doubles trace length (~63) due to per-mutation
  site/base records.
- **Event count is dominated by `BasePushed`** (one per
  assembled pool byte): ~150 for VJ, ~190 for VDJ. Indels and
  mutations add a few; nothing exotic.

---

## 4. Budget assertions

Budgets are set at **~5–10× the current dev-machine mean**, giving
CI machines (~3× slower) headroom while still catching 10×
regressions.

| Workload | N records | Budget (ms) | Approx multiplier of baseline |
|----------|-----------|-------------|-------------------------------|
| W1       | 5,000     | 1,000       | 9.8×                           |
| W2       | 5,000     | 1,500       | 8.4×                           |
| W3       | 1,000     | 2,000       | 8.5×                           |
| W4       | 3,000     | 1,500       | 8.6×                           |
| W5       | 2,000     | 2,000       | 8.2×                           |
| W6       | 1,000     | 1,000       | 8.1×                           |
| W7 ratio | n/a       | replay ≤ 0.8 × fresh | (1.0 - 0.273 ≈ 3× margin)      |

These are **wall-time bounds, not per-record means**. Aggregating
over the workload's N records absorbs single-record noise.

### 4.1 Anti-patterns the budgets catch

Each budget protects against a specific failure mode:

- **W1/W2 ballooning**: a regression that adds Python-side overhead
  to the per-record path (e.g. dict marshalling in the hot loop)
  would hit W1's tight 0.02 ms/rec baseline first.
- **W3 ballooning**: productive narrowing is the most expensive
  per-pass operation. A regression in per-query contract mask
  caching or NP-base admit-mask precomputation would land here.
- **W5 ballooning**: S5F candidate enumeration is `O(positions × 4)`
  per mutation; a regression that scans the full pool on every
  draw would hit here.
- **W6 ballooning**: the indel DP retry budget (`POST_EVENT_RETRY_BUDGET = 16`)
  catches contract-rejection retries. A regression that makes
  retries more expensive (or more frequent) would land here.
- **W7 replay slowdown**: replay is intentionally much faster
  than fresh; if a future change makes replay re-run the
  contract filter on every recorded value, this ratio would
  collapse.

### 4.2 What the budgets explicitly do NOT cover

- **First-record latency**: extension import + first compile
  takes ~100ms; the warmup absorbs this. Budget is for amortised
  throughput, not cold-start.
- **Memory**: not measured here. The Rust kernel allocates per-pass
  IR revisions, but the production hot path uses `Arc` sharing to
  amortise; if memory regresses, a separate `tracemalloc`-based
  audit is the right venue.
- **Concurrency / GIL effects**: tests run single-threaded; no
  GIL-release benchmarks. Multi-thread scaling is a separate
  concern.

---

## 5. Updating the baseline

When a refactor or feature legitimately changes performance:

1. Re-measure the affected workloads on the same developer
   machine (or document the machine).
2. Update §3's table with the new numbers and the date.
3. Recompute §4's budgets at ~5–10× the new baseline.
4. Update [`tests/test_performance_budgets.py`](../tests/test_performance_budgets.py)'s
   `BUDGETS` dict.
5. Note in the commit message which workload's baseline shifted
   and why (e.g. "added per-event contract verify; W3 budget +20%").

**Do not bump budgets reactively to make CI green.** If a budget
fails, profile first (`py-spy --native record -o flame.svg --
python -m pytest tests/test_performance_budgets.py::test_w3...`)
and find the regression's root cause. Only update the baseline
when the slowdown is intentional and accepted.

---

## 6. Drift identified

### 6.1 No Criterion / cargo-bench harness

The Rust kernel has no `cargo bench` harness — only the
Python-side wall-time tests. A focused regression in a single
hot function (e.g. `sample_filtered_result`'s inverse-CDF loop)
might be diluted by the surrounding Python overhead and slip
past the budget.

**Possible fix:** add `engine_rs/benches/` with Criterion
microbenchmarks for the central helpers (`sample_filtered_result`,
`sample_base_with_admit_mask`, `sample_admissible_tuple`).
Out of scope for this slice — the Python-side budgets are the
first line of defence.

### 6.2 No memory budget

Workload memory isn't measured. Per-record memory for the IR
+ trace + event ledger is non-trivial (Arc-shared but still
allocated per simulation). A regression that triggers per-event
clone could double memory without showing in wall time.

**Possible fix:** `tracemalloc`-based audit in a separate slice
once the wall-time budgets are stable.

### 6.3 Single-machine baseline

The §3 table is from a single developer machine. CI runners may
be 2–3× slower; the budgets are sized for that but no CI
calibration has been done. The first time these tests run on the
actual CI infrastructure, the budgets should be reviewed.

**Possible fix:** capture CI's mean for each workload after one
green run and tighten the multiplier if there's headroom (or
loosen if CI hits the budget).
