# Indel provenance audit

**Status:** audit + golden tests + §6.1 fix + §6.2 fix landed.
§6.3 partially addressed by the new test file.

This audit catalogues how the engine handles polymerase-stage
indels (insertion + deletion events from
[`IndelPass`](../engine_rs/src/passes/corrupt/indel.rs)), what
AIRR fields each event affects, how the event ledger surfaces it,
and where the user-visible projection might miss provenance the
engine has internally. The companion test file
[`tests/test_indel_provenance.py`](../tests/test_indel_provenance.py)
pins current behaviour so a future fix doesn't silently change
observable semantics.

Pattern: mirrors the
[primer-trim / end-loss audit](primer_trim_end_loss_audit.md) —
discovery first, code change later (or never).

---

## 1. What's there today

### 1.1 Engine implementation

[`IndelPass`](../engine_rs/src/passes/corrupt/indel.rs) is the
single canonical indel mechanism. Per-event semantics:

- **Count** sampled from a per-simulation distribution
  (`corrupt.indel.count`).
- **Kind** sampled per event: insertion (`Bool(true)`) or deletion
  (`Bool(false)`), with `insertion_prob` configurable.
- **Site** sampled uniformly within the current pool's range — so
  later events target positions that didn't exist before earlier
  events. The pass walks the *current* pool state each step.
- **Base** for insertions only (`corrupt.indel.base[i]` —
  uniform A/C/G/T by default).

Inserted nucleotides carry `flag::INDEL_INSERTED` and
`germline_pos = NONE`. Deletions remove whatever byte happens to
be at the chosen pool position (could be a real germline byte or
a previously-inserted synthetic byte).

**Edge cases:**

- Empty pool: deletion is a no-op for that step; the kind/site
  records are still written. Insertions at position 0 of an empty
  pool succeed.
- Contract-aware: each candidate event's hypothetical post-state
  is verified against active contracts. Under `productive_only`,
  the pass narrows kind/site/base tuples to those that preserve
  the productive frame.

### 1.2 Trace addresses (D3)

| Address                       | Value kind  | Meaning                                                |
|-------------------------------|-------------|--------------------------------------------------------|
| `corrupt.indel.count`         | `Int`       | Total indel events the pass attempted.                 |
| `corrupt.indel.kind[i]`       | `Bool`      | `true` = insertion, `false` = deletion (per event).    |
| `corrupt.indel.site[i]`       | `Int`       | Pool position chosen for the event.                    |
| `corrupt.indel.base[i]`       | `Base`      | Inserted base (only present for insertion events).     |

### 1.3 Simulation events

Each applied indel emits exactly one
[`SimulationEvent`](../engine_rs/src/ir/sim_event.rs):

- `IndelInserted { at, base, segment, flags }` — for insertions.
- `IndelDeleted { at, removed_base }` — for deletions.

These appear in the pass's `EventRecord.simulation_events`.
The [`LiveCallRefreshHook`](../engine_rs/src/live_call/refresh_hook.rs)
plan-derives `LiveCallRefreshStep::AllStructural` from either
variant — a structural indel triggers a full V/D/J refresh.

---

## 2. AIRR fields affected

What a downstream analyst sees in the projected AIRR record:

| AIRR field                                       | Indel contribution                                                |
|--------------------------------------------------|-------------------------------------------------------------------|
| `sequence`, `sequence_aa`                        | Yes (length changes; insertion bases appear, deletions remove).   |
| `sequence_length`                                | Yes (pool length delta).                                          |
| `n_indels`                                       | **Trace-derived count** of attempted events. See §3.1.            |
| `v_cigar`, `d_cigar`, `j_cigar`                  | **Yes** — walker emits `I` (insertion) and `D` (deletion) ops inside the affected segment. |
| `v_sequence_start/end` (etc.)                    | Indirect — segment-region pool coordinates shift after an indel.  |
| `v_germline_start/end` (etc.)                    | Walker-derived ref ranges. Deletions appear as covered ref positions (`D` ops). Insertions don't consume ref positions. |
| `v_identity`, `j_identity`, etc.                 | Yes — indels increment the denominator (total ops) but not the match numerator. |
| `productive`, `vj_in_frame`, `stop_codon`        | Yes — frame shifts and stop codons follow from net indel offset.  |
| `junction`, `junction_aa`, `junction_length`     | Yes — recomputed after indel.                                     |

### 2.1 Walker handling

[`walk_alignment_columns`](../engine_rs/src/airr_record/walk/mod.rs)
detects indels structurally:

- **Insertion** (`nuc.germline_pos == NONE` inside a V/D/J region)
  → emits `I` CIGAR op for that pool byte; no germline credit;
  identity denominator increments but numerator doesn't.
- **Deletion** (a germline position that doesn't appear in the
  region's bytes when walking sequentially) → fill loop emits `D`
  CIGAR op for each missing germline position; identity
  denominator increments.

CIGAR provenance is therefore **per-segment** and structurally
correct. A reader inspecting `v_cigar` can count `I`/`D` ops to
get the within-V indel count without reading the trace.

---

## 3. AIRR fields NOT affected — the provenance gaps

### 3.1 `n_indels` now reports applied changes — resolved in §6.1 fix

**Historical (pre-fix):** `n_indels` was derived from the trace's
`corrupt.indel.count`, i.e. the *attempted* count. Under deletion-
heavy fixtures where the pool empties mid-pass, the no-op slots
were counted in `n_indels` but emitted no `SimulationEvent`. Same
for `productive_only` no-ops where the contract rejects every
candidate (the trace records the attempt, no event fires).

**Current behaviour:** `n_indels` is the count of `IndelInserted`
+ `IndelDeleted` events emitted by the `corrupt.indel` pass — the
realized structural-change count. See
[`build_airr_record`](../engine_rs/src/airr_record/builder.rs)
near the `n_indels` assignment. End-loss reuses `IndelDeleted` as
its IR primitive (see [`SimulationEvent::BaseDeleted`](../engine_rs/src/ir/sim_event.rs)
which is "Reserved for future emission"), so the count is scoped
to `pass_name == "corrupt.indel"` to exclude end-loss removal.

### 3.2 Per-segment counters now exposed — resolved in §6.2 fix

**Historical (pre-fix):** `n_indels` was the only indel-count
field on the AIRR record. Users wanting "how many insertions
happened inside V" had to parse `v_cigar` for `I` ops.

**Current behaviour:** the AIRR record exposes `n_v_indels`,
`n_d_indels`, `n_j_indels` — applied indel counts attributed by
the event's `segment` field. CIGAR remains authoritative for
*where* inside a segment the indel landed; the counters
denormalize the per-segment frequency. Per-kind splits
(`v_n_insertions` vs `v_n_deletions`) and `n_indels_in_junction`
remain unsplit — CIGAR-parsing is still the pathway for those
finer breakdowns.

NP1 / NP2 events are excluded from per-segment counters by
design (no V/D/J ownership). See §6.2 for the boundary rule.

### 3.3 No insertion/deletion ratio counter

Two parameters of the IndelPass — `count` and `insertion_prob` —
are not directly visible on the record. The user can recover both
from the trace (count records + kind boolean prevalence) but
neither is a standalone AIRR field.

### 3.4 Site/base of each event not surfaced

`corrupt.indel.site[i]` and `corrupt.indel.base[i]` are
trace-only. The CIGAR position implicitly carries the site;
inserted bases appear in `sequence` at the insertion position;
neither is denormalized.

This is consistent with the AIRR convention: per-event provenance
is encoded structurally, not as scalar fields.

---

## 4. Expected biological semantics

| Question                                | Answer                                                              |
|-----------------------------------------|---------------------------------------------------------------------|
| Does `n_indels` count attempts or applied changes? | **Applied changes** (post-§6.1 fix). Equal to `len(EventRecord.simulation_events filter Indel*)` over the `corrupt.indel` pass record. Trace's `corrupt.indel.count` still records attempts. |
| Does CIGAR distinguish insertion/deletion? | **Yes** — `I` vs `D` ops in the affected segment's CIGAR.        |
| Does the productive contract preserve frame across indels? | **Yes,** via per-site classification: each site is labeled `FrameNeutral` (before V_region.start OR ≥ J_region.start) or `FrameDelta(±1)` (inside V/D/NP). The DP picks tuples with net change 0 mod 3, then validates anchor preservation. See §5. |
| Are indel events deterministic under trace replay? | **Yes.** Trace carries kind/site/base; replay reproduces the exact byte changes. |
| Does live-call refresh fire after an indel? | **Yes** — `LiveCallRefreshHook` consumes `IndelInserted`/`IndelDeleted` events and triggers a full V/D/J refresh (`AllStructural` step). |
| Are indel events part of `simulation_events`? | **Yes** — one event per applied change, in pool-position order. |

---

## 5. Productive-contract interaction

Under `productive_only`, the productive bundle narrows indel
candidates **per-site**, not just on net frame shift. The
admissibility classifier in
[`ProductiveJunctionFrame::admissible_indel_class_at`](../engine_rs/src/contract/productive_junction_frame.rs)
labels each site as either `FrameNeutral` or `FrameDelta(±1)`
based on its position relative to V/J region starts; the rest of
the bundle further excludes sites that would corrupt the V or J
anchor codon.

### 5.1 Per-site frame classification

For an indel at site `s`:

- `s < V_region.start` → **FrameNeutral** (both V and J region
  starts shift; junction length unchanged). Typically empty since
  `V_region.start = 0`.
- `V_region.start ≤ s < J_region.start` → **FrameDelta(±1)** (V
  stays, J shifts; junction length ±1). This is V + NP + (D + NP2
  if VDJ).
- `s ≥ J_region.start` → **FrameNeutral** (neither V nor J start
  shifts; junction length unchanged). This is the J region and
  the end-of-pool insertion slot.

The mod-3 DP in the indel pass aggregates this classification
across the candidate tuple — net frame change must be 0 mod 3 —
and then validates the assembled post-state against the rest of
the contract bundle (anchor preservation, no junction stop
codons).

### 5.2 Single-indel admissibility (audit correction)

The original §5 claimed "one indel: net shift ±1 → contract
rejects every candidate." That is **too broad**: it conflates
"the indel always shifts the *junction* frame" with "any indel
landing anywhere shifts the frame." The classifier admits any
FrameNeutral site, and FrameNeutral sites exist whenever the J
region has any pool bases.

Empirically pinned by
[`test_productive_only_count_1_lands_outside_junction_and_anchors`](../tests/test_indel_provenance.py):
on the 18-base baseline (V=9 with V_anchor=6, J=6 with J_anchor=0)
single-indel sites under `productive_only` land in
`[J_anchor_pool + 3, pool.len()] = [15, 18]` — past the J anchor
codon `[12, 15)` and outside the junction window `[6, 15)`.

The admissible-site set for a count=1 indel under `productive_only`
is therefore:

`{s ∈ valid_range : FrameNeutral(s) ∧ s ∉ V_anchor_codon ∧ s ∉ J_anchor_codon}`

This set is **non-empty whenever** the J region extends past its
anchor codon (i.e., J has more than 3 bases) **or** an end-of-pool
insertion is admissible. It is **empty** when (a) the J region is
exactly the J anchor codon (no post-anchor bases) AND (b) the
kind is "delete-only" (no end-of-pool insertion slot to use).

### 5.3 Strict vs. permissive in the unreachable-set case

Pinned by the short-J fixture
([`test_productive_only_count_1_delete_only_short_j_no_ops_permissive`](../tests/test_indel_provenance.py)
and
[`...raises_strict`](../tests/test_indel_provenance.py)):
J=3 bases with anchor=0 → J anchor codon = entire J → no
post-anchor zone. Single-deletion candidates are all inadmissible.

- **Permissive** (default): the slot is recorded as `IndelEvent::NoOp`
  — trace shows `corrupt.indel.kind[0] = false`,
  `corrupt.indel.site[0] = -1`, no `SimulationEvent::IndelDeleted`
  fires. Post-§6.1 fix, `n_indels` reports 0 for these records.
- **Strict**: the indel pass raises `StrictSamplingError` with
  reason `empty_admissible_support`. The error tuple is
  `(pass_name, choice_address, reason)` — here
  `("corrupt.indel", "corrupt.indel.site[0]", "empty_admissible_support")`.

Switching the same fixture to *insertion-only* opens a single
admissible site (`s = pool.len() = 15`, end-of-pool, FrameNeutral,
outside both anchor codons) — pinned by
[`test_productive_only_count_1_insert_only_short_j_applies_at_end_of_pool`](../tests/test_indel_provenance.py).
Both modes succeed.

### 5.4 Count=2 and count=3

- **Count=2**: net shift can be 0 (insertion + deletion paired in
  V/NP/D space, or both indels FrameNeutral). The DP picks a
  net-zero tuple if one exists.
- **Count=3**: net shift ±3 = 0 mod 3 (all insertions or all
  deletions in FrameDelta space) — or 1 FrameDelta + 2
  FrameNeutral, etc. Contract permits.

The productive-contract stress matrix
([`tests/test_productive_stress_matrix.py`](../tests/test_productive_stress_matrix.py))
exercises `indel_count_0/2/3` under both `productive_only` and
unconstrained, and confirms records stay productive under the
contract. ✅

---

## 6. Drift identified

Items below are **not fixed in this slice**. They're catalogued
so the next slice can scope a focused fix (or decide no fix is
needed):

### 6.1 `n_indels` overcounts under pool-empty no-ops — **RESOLVED**

`n_indels` is now derived from the event ledger: it counts
`IndelInserted` + `IndelDeleted` `SimulationEvent`s emitted by
the `corrupt.indel` pass. Pool-empty / contract-rejected slots
record the attempt in the trace but emit no event, so they no
longer count toward `n_indels`.

**Behaviour change:** "records with at least one indel" filters
now reflect *realized* structural changes. The trace's
`corrupt.indel.count` continues to record attempts, so users who
need the attempted count can still read it directly.

**Implementation scope:** in [`build_airr_record`](../engine_rs/src/airr_record/builder.rs)
the `n_indels` assignment iterates `outcome.events()`, filters to
the `corrupt.indel` pass (to exclude `EndLossPass`'s reuse of
`IndelDeleted` as its IR primitive), and counts the matching
`SimulationEvent::IndelInserted` / `IndelDeleted` variants.

**Pinned by:** `test_n_indels_excludes_pool_empty_no_op_attempts`
in [`tests/test_indel_provenance.py`](../tests/test_indel_provenance.py)
— uses `count=30, insertion_prob=0.2, seed=0` to drive a
pool-emptying fixture where trace attempts (30) ≠ events fired
(25), and asserts `n_indels == 25`.

### 6.2 No per-segment indel counter — **RESOLVED**

The AIRR record now exposes three per-segment counters:
`n_v_indels`, `n_d_indels`, `n_j_indels`. Each counts the
`IndelInserted` + `IndelDeleted` events from the `corrupt.indel`
pass whose `segment` field attributes the event to V / D / J.
The CIGAR remains authoritative for site-level provenance; the
counters denormalize the per-segment frequency for ergonomic
SQL / dataframe filtering.

**Attribution convention:**

- **Deletion:** the deleted nucleotide's own segment (read from
  the pre-event pool). To support this without reaching back into
  pre-event simulation state, `SimulationEvent::IndelDeleted` now
  carries a `segment: Segment` field, symmetric with
  `IndelInserted.segment`. The IR builder's `delete_indel`
  populates it from the removed `Nucleotide`'s segment.
- **Insertion:** computed by [`insertion_segment`](../engine_rs/src/passes/corrupt/indel/event.rs)
  — region containment first (start ≤ at < end, strict), then
  fallback to the byte at `at` (if in pool), then to `at-1`. As a
  consequence:
  - Insertion at exactly `region.end` attributes to the *next*
    region (since the in-region check is strict on the end).
  - Insertion at `pool.len()` (end-of-pool) attributes to the
    last segment via the `at-1` fallback.

**Excluded:** events landing in NP1 / NP2 increment `n_indels`
but NOT any per-segment counter. Per-NP counters are deferred
indefinitely until a downstream consumer asks for them — the
scope is pinned by
`test_no_per_np_segment_counter_on_airr_record`.

**Implementation scope:**

- [`engine_rs/src/ir/sim_event.rs`](../engine_rs/src/ir/sim_event.rs):
  `IndelDeleted` gains the `segment` field.
- [`engine_rs/src/ir/builder.rs`](../engine_rs/src/ir/builder.rs):
  `delete_indel` populates `segment` from `removed.segment`.
- [`engine_rs/src/airr_record/record.rs`](../engine_rs/src/airr_record/record.rs):
  `n_v_indels`, `n_d_indels`, `n_j_indels` added to `AirrRecord`.
- [`engine_rs/src/airr_record/builder.rs`](../engine_rs/src/airr_record/builder.rs):
  same event-iteration loop that populates `n_indels` also
  matches on `segment` to bucket per-segment counts.
- [`engine_rs/src/python/outcome.rs`](../engine_rs/src/python/outcome.rs)
  and [`src/GenAIRR/result.py`](../src/GenAIRR/result.py):
  expose through the Python dict and CSV schema.

**Pinned by:**

- `test_per_segment_counters_exposed_on_airr_record` — fields are
  present and non-negative.
- `test_per_segment_counters_sum_to_global_n_indels` — invariant
  `n_v + n_d + n_j ≤ n_indels` across 20 seeds.
- `test_single_deletion_inside_v_increments_n_v_indels` (seed=1).
- `test_single_insertion_inside_j_increments_n_j_indels` (seed=0).
- `test_single_indel_inside_np1_excluded_from_per_segment_counters`
  (seed=16) — NP boundary case: `n_indels == 1` but all per-
  segment counters stay zero.
- `test_no_per_np_segment_counter_on_airr_record` — pins the
  V/D/J-only scope.

### 6.3 Test coverage gap for boundary cases

There are no existing tests pinning:

- An insertion at pool position 0 (immediately before V).
- A deletion at pool position 0 (removing V's first byte).
- An insertion/deletion at the V-NP1 boundary.
- An indel inside the junction window.

The productive-contract stress matrix exercises the existence of
indels, not the boundary correctness. The new test file
catalogued in §7 adds the first round of boundary pins.

---

## 7. Test coverage in this slice

[`tests/test_indel_provenance.py`](../tests/test_indel_provenance.py)
pins the current behaviour with golden tests. Coverage:

- **Single-event invariants:** insertion increments
  `sequence_length` by 1 and emits one `IndelInserted` event;
  deletion decrements by 1 and emits one `IndelDeleted` event.
- **Trace ↔ event consistency:** `n_indels` from trace matches
  the count of `IndelInserted` + `IndelDeleted` events in the
  pass's `EventRecord`, under typical (non-pool-empty) fixtures.
- **CIGAR carries indel provenance:** when an indel falls inside
  V's region, `v_cigar` contains an `I` or `D` op.
- **Productive interaction:** under `productive_only`, a
  count=2 indel pass with `insertion_prob=0.5` produces records
  that satisfy the productive triad (frame preserved by
  insertion-deletion pairing).
- **Replay round-trip:** `trace_file_from` + `replay_from_trace_file`
  reproduces sequence, every key AIRR field, AND the per-pass
  `EventRecord.simulation_events` indel count.
- **Pin-the-drift:** explicit assertion that `n_indels` matches
  the event count (currently true for the deterministic fixtures
  used; if a future fix changes the semantic, this test signals
  the change).

The audit doc and the test file are designed to be read
together: doc explains why, tests show what.
