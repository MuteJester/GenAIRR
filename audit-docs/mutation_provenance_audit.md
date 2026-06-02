# Mutation Provenance Counters — Architecture Contract

**Status: shipped.** The implementation slice landed with the
shape this audit specified: four AIRR `i64` fields
(`n_v_mutations` / `n_d_mutations` / `n_j_mutations` /
`n_np_mutations`) aggregated from `outcome.events()` filtered
to the SHM passes only, with the documented sum invariant
`n_v + n_d + n_j + n_np == n_mutations` enforced via five new
validator issue kinds (`N{V,D,J,Np}MutationsMismatch` +
`MutationCountSumMismatch`). NP1 + NP2 roll together as
specified; corruption-stage `BaseChanged` events are excluded
by the pass-name filter; replay reproduces all four fields
byte-identically.

Release-tier consolidation closed the slice the same way the
targeted-SHM slice closed: full-stack productive IGH +
non-default segment rates + sum-invariant check + corruption
isolation in [`tests/test_release_validation.py`](../tests/test_release_validation.py),
plus a row in [`docs/validation_matrix.md`](validation_matrix.md)
naming the audit doc / contract file / implementation tests /
Rust kernel sites.

---

The original pre-implementation audit follows verbatim below.
Pins today's global mutation-counter shape and audits the
per-segment counter extension — the shape that ultimately
shipped.

This audit is the natural follow-up to the targeted-SHM slice
(`docs/shm_segment_rate_design.md`, shipped): now that users can
restrict SHM to V / D / J / NP via `mutate(segment_rates=...)`,
the next question is *how many mutations actually landed in
each segment* — a question today's `n_mutations` global field
can't answer.

Companion to
[`tests/test_mutation_provenance_contract.py`](../tests/test_mutation_provenance_contract.py)
which freezes today's contract (`pin_scaffold_*`) and the gaps
the future slice closes (`pin_absence_*`).

**Pre-flight check:** the audit confirmed
`SimulationEvent::BaseChanged` already carries `segment:
Segment` + `germline_pos: Option<u16>`. **The implementation
slice is a clean Rust-side aggregation against the existing
event ledger** — no new event field required, no new IR
surface, no replay-format change.

---

## Existing scaffolding the contract pins

| Surface | Where | What it pins |
|---|---|---|
| `Simulation.mutation_count: u32` | [`engine_rs/src/ir/simulation.rs:43`](../engine_rs/src/ir/simulation.rs#L43) | IR-side SHM counter; survives parent→descendant boundary. Sourced from `MutationTransaction::add_to_mutation_count`. |
| `MutationTransaction::add_to_mutation_count` call sites | [`engine_rs/src/passes/mutate/uniform.rs`](../engine_rs/src/passes/mutate/uniform.rs) + [`s5f/execution.rs`](../engine_rs/src/passes/mutate/s5f/execution.rs) | The two — and only the two — passes that bump biological-mutation counts. Source-grepped by the SHM model audit's `pin_scaffold_only_mutate_passes_call_add_to_mutation_count`. |
| `SimulationEvent::BaseChanged { handle, old_base, new_base, segment, germline_pos }` | [`engine_rs/src/ir/sim_event.rs:119`](../engine_rs/src/ir/sim_event.rs#L119) | **Carries segment** + `germline_pos` for free. The per-event provenance the future slice consumes. |
| `EventRecord.pass_name` + `.simulation_events` | [`engine_rs/src/pass/`](../engine_rs/src/pass/) | The pass-name filter pattern indel counters use: `outcome.events()` iterator, filter `pass_name == address::CORRUPT_INDEL`, walk `simulation_events` for typed events. |
| AIRR `n_mutations` (global) | [`engine_rs/src/airr_record/builder.rs:211`](../engine_rs/src/airr_record/builder.rs#L211) | IR-sourced from `sim.mutation_count`. |
| AIRR `n_indels` / `n_v_indels` / `n_d_indels` / `n_j_indels` | [`engine_rs/src/airr_record/builder.rs:225-252`](../engine_rs/src/airr_record/builder.rs#L225-L252) | **The reference implementation pattern** for the proposed per-segment SHM counters. Filters events by `pass_name == address::CORRUPT_INDEL`, switches on `IndelInserted` / `IndelDeleted` segment, increments per-bucket counts. NP1/NP2 currently roll up into `n_indels` total only — same recommendation here. |
| AIRR `n_pcr_errors` / `n_quality_errors` | [`engine_rs/src/airr_record/builder.rs:217-218`](../engine_rs/src/airr_record/builder.rs#L217-L218) | Trace-sourced: `trace_int_choice(trace, ChoiceAddress::CorruptPcrCount)`. Counts the *attempted* count drawn at sample time, not the realised count after constraint filtering. |
| AIRR `end_loss_5_length` / `end_loss_3_length` | [`engine_rs/src/airr_record/builder.rs:258-261`](../engine_rs/src/airr_record/builder.rs#L258-L261) | Trace-sourced similarly. |
| Address constants | [`engine_rs/src/address.rs:88-130`](../engine_rs/src/address.rs#L88-L130) | `MUTATE_UNIFORM = "mutate.uniform"` + `MUTATE_S5F = "mutate.s5f"` — the pass-name filter values the slice needs. |
| Python `EventRecord.simulation_event_count` | [`engine_rs/src/python/event.rs:128-135`](../engine_rs/src/python/event.rs#L128-L135) | Python-side surface exposes the count but NOT the typed events. Python consumers can't compute per-segment SHM independently today; the slice doesn't need to (Rust side aggregates), but the gap is documented. |

---

## 1. Q1 — Current global counters

The AIRR record's counter fields today, with their source and
biological interpretation:

| Field | Source | Counts |
|---|---|---|
| `n_mutations` | IR (`sim.mutation_count`) | Realised biological SHM substitutions (uniform + S5F only). |
| `mutation_rate` | Derived (`n_mutations / sequence_length`) | Per-base SHM rate of the final sequence. |
| `n_pcr_errors` | Trace (`corrupt.pcr.count`) | Attempted PCR substitution count (sampled, not realised). |
| `n_quality_errors` | Trace (`corrupt.quality.count`) | Attempted quality-error count. |
| `n_indels` | Events (`pass_name == "corrupt.indel"`, `IndelInserted/Deleted`) | Realised structural indels — INCLUDES NP1/NP2 events. |
| `n_v_indels` | Events filtered to `Segment::V` | Realised V-region indels. |
| `n_d_indels` | Events filtered to `Segment::D` | Realised D-region indels. |
| `n_j_indels` | Events filtered to `Segment::J` | Realised J-region indels. |
| `end_loss_5_length` | Trace (`corrupt.end_loss.5`) | Attempted 5' end-loss length. |
| `end_loss_3_length` | Trace (`corrupt.end_loss.3`) | Attempted 3' end-loss length. |
| `is_contaminant` | Trace (`corrupt.contaminant.applied`) | Bool flag. |

### Two source kinds: trace vs events

- **Trace-sourced counters** (PCR / quality / end-loss): the
  count drawn at sample time. Constraint filtering or empty
  support can reduce the *realised* count below the recorded
  count. The AIRR field reflects the **attempted** intent, not
  the actual outcome.
- **Event-sourced counters** (indels): the realised count
  aggregated from the per-pass event ledger. Constraint-rejected
  sites don't emit events, so the count reflects what actually
  happened.

`n_mutations` is the third pattern — **IR-sourced**: the pass
itself bumps `Simulation.mutation_count` only for realised
substitutions (via `add_to_mutation_count(applied)` where
`applied` excludes silently-skipped slots). The AIRR field reads
this counter directly. Realised semantics, no event-walk.

### Recommendation for the slice

Use the **event-sourced** pattern for per-segment SHM counters
— identical to `n_v/d/j_indels`. Three reasons:

1. **Realised counts.** Same semantic as the global `n_mutations`
   (which is realised). A user comparing `n_v_mutations +
   n_d_mutations + n_j_mutations + n_np_mutations` against
   `n_mutations` should see the totals add up.
2. **Segment information already on the event.** `BaseChanged`
   carries `segment` — no IR lookup needed.
3. **Architectural consistency.** Indels are the precedent; SHM
   should match.

Pinned by `pin_scaffold_n_mutations_is_realised_count` (already
covered by the SHM model audit).

---

## 2. Q2 — What counts as biological mutation?

The audit's classification matrix — which `BaseChanged` events
contribute to `n_mutations` (and therefore to the proposed
per-segment counters):

| Pass | Emits `BaseChanged`? | Counts as SHM? | Counter source |
|---|---|---|---|
| `UniformMutationPass` (`mutate.uniform`) | Yes | **Yes** | `add_to_mutation_count` + events |
| `S5FMutationPass` (`mutate.s5f`) | Yes | **Yes** | Same |
| `PCRErrorPass` (`corrupt.pcr`) | Yes | No (observation-stage) | Trace `n_pcr_errors` only |
| `QualityErrorPass` (`corrupt.quality`) | Yes | No | Trace `n_quality_errors` only |
| `NCorruptionPass` (`corrupt.ns`) | Yes | No | None today |
| `IndelPass` (`corrupt.indel`) | No (emits `IndelInserted/Deleted`) | No | Events filtered by pass name |
| `EndLossPass` (`corrupt.end_loss.*`) | No (emits `IndelDeleted` as primitive) | No | Trace |
| `RevCompPass` (`corrupt.rev_comp`) | No (projection-only) | No | Trace bool |
| `ContaminantPass` (`corrupt.contaminant`) | Yes (whole pool overwrite) | No | Trace bool |
| `ReceptorRevisionPass` (`receptor_revision`) | No (emits `SegmentReplaced` + `AssignmentChanged`) | No | IR `original_v_call` |
| `InvertDPass` (`invert_d`) | No (mutates `assignments.d.orientation`) | No | IR `d_inverted` |

### The classification rule

**Biological SHM events = `BaseChanged` events emitted by passes
whose `pass_name` is in `{address::MUTATE_UNIFORM,
address::MUTATE_S5F}`.**

Every other `BaseChanged`-emitting pass is observation-stage
corruption and must NOT contribute. This is identical to the
indel pattern (`pass_name == address::CORRUPT_INDEL`) — pin the
filter rule in the contract test.

### Edge cases

- **`add_to_mutation_count` invariant**: only the two SHM passes
  call it (source-greppped pin from the SHM model audit). So a
  future biology pass that *should* count as SHM must call both
  `add_to_mutation_count` (for `n_mutations`) AND fire
  `BaseChanged` events (for per-segment counts). The slice should
  pin both call-site invariants — the source-grep already covers
  the first; the event-emission convention is documented
  per-pass.
- **Indels are NOT SHM** even when they're biologically real (AID
  can produce indels). v1 of GenAIRR doesn't model SHM indels;
  `polymerase_indels` is library-prep. Pin this as a documented
  modeling choice — the per-segment SHM counters count only
  substitutions.

Pinned by `pin_scaffold_only_mutate_passes_emit_shm_basechanged`
— behavioural pin walking the events ledger after a known
configuration and asserting only mutate.* passes contribute.

---

## 3. Q3 — Segment attribution rule

**`BaseChanged` carries `segment` directly.** No post-hoc
coordinate lookup needed.

Looking at the IR event spec
([`sim_event.rs:114-125`](../engine_rs/src/ir/sim_event.rs#L114-L125)):

```rust
/// The base byte at `handle` was replaced. `segment` and
/// `germline_pos` describe the nucleotide at `handle` —
/// they are unchanged by the substitution (only the base byte
/// flips), so the same values apply both before and after.
BaseChanged {
    handle: NucHandle,
    old_base: u8,
    new_base: u8,
    segment: Segment,
    germline_pos: Option<u16>,
}
```

For substitutions, the pre- and post-event segment are the same
(substitution doesn't change which biological segment a position
belongs to). So the per-event `segment` field is unambiguous.

### Mapping to AIRR buckets

Following the four-bucket convention from the segment-rates
slice:

| `Segment` variant | AIRR counter field |
|---|---|
| `V` | `n_v_mutations` |
| `D` | `n_d_mutations` |
| `J` | `n_j_mutations` |
| `Np1` + `Np2` | `n_np_mutations` (rolled together) |

The same rollup pattern the `SegmentRateWeights` config uses —
NP1 and NP2 are tuned together at the cartridge / experiment
level, so reporting them separately at the AIRR level would be
gratuitous granularity. A future per-NP analysis slice can ride
a separate field if demand emerges.

### Indels-vs-end-loss precedent

`n_indels` filters `IndelInserted` / `IndelDeleted` events by
`pass_name == "corrupt.indel"` so the end-loss pass's
`IndelDeleted` events don't leak in. The proposed per-segment
SHM counters should similarly filter by `pass_name in {"mutate.
uniform", "mutate.s5f"}` so PCR / quality / contaminant
`BaseChanged` events don't leak in.

### Under indels / end-loss

Indels shift the pool coordinate space, but `BaseChanged.segment`
is set at event-emission time and is therefore stable —
post-indel coordinate shifts don't invalidate it. The slice is
safe under any combination of SHM + indels + end-loss in any
order.

Pinned by `pin_scaffold_basechanged_carries_segment_field` —
source-level grep on the IR event definition.

---

## 4. Q4 — Do events carry enough data?

**Yes.** The Rust event ledger has everything the slice needs:

- `outcome.events()` returns `Vec<EventRecord>`.
- Each `EventRecord` has `pass_name: String` + `simulation_events: Vec<SimulationEvent>`.
- `SimulationEvent::BaseChanged` carries `segment: Segment` (and `germline_pos: Option<u16>` for free).

The indel counter implementation in
[`airr_record/builder.rs:229-247`](../engine_rs/src/airr_record/builder.rs#L229-L247)
is the canonical pattern:

```rust
for record in outcome.events() {
    if record.pass_name != address::CORRUPT_INDEL {
        continue;
    }
    for ev in &record.simulation_events {
        let segment = match ev {
            SimulationEvent::IndelInserted { segment, .. } => *segment,
            SimulationEvent::IndelDeleted { segment, .. } => *segment,
            _ => continue,
        };
        n_indels += 1;
        match segment {
            Segment::V => n_v_indels += 1,
            Segment::D => n_d_indels += 1,
            Segment::J => n_j_indels += 1,
            Segment::Np1 | Segment::Np2 => {}  // rolls into total only
        }
    }
}
```

The SHM counter slice translates this verbatim:

```rust
let mut n_v_mutations = 0i64;
let mut n_d_mutations = 0i64;
let mut n_j_mutations = 0i64;
let mut n_np_mutations = 0i64;
for record in outcome.events() {
    if record.pass_name != address::MUTATE_UNIFORM
        && record.pass_name != address::MUTATE_S5F
    {
        continue;
    }
    for ev in &record.simulation_events {
        let SimulationEvent::BaseChanged { segment, .. } = ev else {
            continue;
        };
        match segment {
            Segment::V => n_v_mutations += 1,
            Segment::D => n_d_mutations += 1,
            Segment::J => n_j_mutations += 1,
            Segment::Np1 | Segment::Np2 => n_np_mutations += 1,
        }
    }
}
```

### Python observability gap

The Python `PyEventRecord` exposes
`simulation_event_count` but **not** the typed event list. So a
Python consumer can see how many simulation events a pass
emitted, but cannot inspect per-event segment from Python today.

For the proposed slice this is fine — the Rust builder
aggregates and projects four AIRR fields. Python consumers read
the fields directly without needing to walk events.

For a *deeper* "where did each mutation land" surface (e.g. a
list of per-mutation positions), the Python `EventRecord`
surface would need extending — but that's a separate slice. The
v1 counter slice doesn't need it.

Pinned by `pin_absence_python_eventrecord_lacks_typed_events`.

---

## 5. Q5 — Counter shape recommendation

### Recommended: four new AIRR fields

```text
n_v_mutations   : int
n_d_mutations   : int
n_j_mutations   : int
n_np_mutations  : int
```

Sums to `n_mutations` (global) by construction — derive the
total from the same event walk to guarantee:

```text
n_mutations == n_v_mutations + n_d_mutations + n_j_mutations + n_np_mutations
```

The global `n_mutations` STAYS — it's the canonical "how much
SHM did this record see?" field consumers already use. The four
per-segment fields are additive provenance.

### Rejected: per-NP1 / per-NP2 split

The segment-rates DSL collapses NP1+NP2 into a single "NP"
bucket. The counter surface should match that grouping so users
have one mental model. A future "per-NP" slice can add `n_np1_*`
/ `n_np2_*` fields if a clear use case emerges; not in v1.

### Rejected: rate fields per segment

A `v_mutation_rate` / `d_mutation_rate` / etc. surface would
double the field count without adding information (the rate is
trivially `n_<seg>_mutations / <seg>_length`, and the segment
lengths are already on the AIRR record via the V/D/J coordinate
fields). Defer to consumer code.

### Field types

`i64`, matching `n_mutations` and `n_v/d/j_indels`. Default to
`0` for records that never ran SHM (consistent with
`n_indels=0` default on records that didn't run indels).

Pinned by `pin_absence_no_per_segment_mutation_counters` —
behavioural pin that the four field names don't yet exist on
the AIRR record.

---

## 6. Q6 — Validator surface

### Today's validator coverage

The per-record validator
([`engine_rs/src/airr_record/validate.rs`](../engine_rs/src/airr_record/validate.rs))
re-derives every AIRR field from the outcome state and compares
field-by-field. It catches:

- `n_mutations` drift (re-derived from `sim.mutation_count`).
- `n_indels` / `n_v_indels` / `n_d_indels` / `n_j_indels` drift
  (re-derived from the event ledger filtered by
  `corrupt.indel`).
- `n_pcr_errors` / `n_quality_errors` / `end_loss_*_length`
  drift (re-derived from the trace).

### Proposed extension

Add the same per-segment SHM counter re-derivation. Issue kinds:

- `VMutationCountMismatch { reported, expected }`
- `DMutationCountMismatch { reported, expected }`
- `JMutationCountMismatch { reported, expected }`
- `NpMutationCountMismatch { reported, expected }`

Plus an existing-style cross-check:

- `MutationCountSumMismatch { reported, expected }` — fires when
  `n_v_mutations + n_d_mutations + n_j_mutations + n_np_mutations
  != n_mutations`.

The cross-check is the load-bearing pin: if any one of the four
new fields drifts but the global stays consistent (or vice
versa), the sum invariant catches it. Same shape as the
indel-sum invariant the existing validator already runs.

### Source-tag format

Matching the existing `details.source` convention:

- For the per-bucket fields: `"events:mutate.{uniform,s5f}:base_changed:{V,D,J,NP}"`.
- For the sum mismatch: `"derived:n_v_mutations+n_d_mutations+n_j_mutations+n_np_mutations"`.

Pinned by `pin_absence_no_mutation_count_validator_kinds`.

---

## 7. Q7 — Trace / replay

### Counters are event-derived, not trace-derived

The proposed slice consumes `outcome.events()`, not the trace.
This matters for replay semantics:

- **Same trace + same plan**: replay reproduces the same event
  ledger, so the per-segment counters reproduce exactly. Byte-
  identical output.
- **Same trace + different plan** (e.g. different
  `segment_rates`): the plan signature check already catches
  this — replay fails the signature gate before touching the
  counters.
- **Refdata swap**: same shape — refdata signature catches it.

The counters don't introduce any new replay-format change.
**No new trace addresses**, no new IR field.

### Replay determinism pin

Same shape as the segment-rates slice's replay round-trip pin:
record a fresh outcome, build a `TraceFile`, rerun via
`rerun_from_trace_file`, assert the replayed AIRR record's four
new fields match the fresh ones exactly. The implementation
slice's release-tier test should include this.

Pinned by the existing `pin_scaffold_shm_replay_byte_deterministic`
(re-used across SHM audits).

---

## 8. Edge cases the implementation slice must handle

1. **Records with zero SHM.** All four fields default to `0`;
   sum invariant holds trivially.

2. **Heavy SHM with productive constraint.** Permissive mode
   skips constraint-rejected sites silently; the event ledger
   reflects only realised mutations. Per-segment counters
   reflect realised counts. Sum invariant holds.

3. **Segment-rate zero on NP.** Under
   `segment_rates={"NP": 0.0}`, `n_np_mutations` should be `0`
   for every record. The sum invariant holds because no
   `BaseChanged` events fire on NP positions (the pass excludes
   them from support before sampling).

4. **Clonal pipelines.** Each descendant's event ledger is
   independent (descendants don't carry parent events — pinned
   by the clonal-parent audit). Per-segment counters are
   per-descendant, computed from each descendant's own events.

5. **TCR pipelines.** `mutate` is rejected on TCR refdata at
   the DSL boundary, so per-segment SHM counters stay at `0`
   for TCR records by construction.

6. **Repeated mutate calls.** GenAIRR doesn't currently support
   stacking two `mutate(...)` calls; only one is allowed per
   pipeline. The counter aggregation walks all
   `mutate.{uniform,s5f}` pass-name events anyway — forward-
   compatible if duplicate-mutate is ever allowed.

7. **Indels after SHM shifting positions.** Doesn't matter:
   `BaseChanged.segment` is recorded at SHM time. Post-indel
   coordinate shifts don't invalidate the segment classification.

8. **End-loss truncating mutated positions.** Same as #7. The
   event was recorded; if end-loss then chops the position off
   the pool, the event ledger still has it. So the counter
   reflects "mutations that happened" not "mutations still
   present in the final sequence". This matches the indel
   counter's behaviour — pinned as the architectural contract.

---

## 9. Performance

### Cost of the counter walk

`outcome.events()` is iterated once per AIRR record projection
anyway (the existing indel counter walks the same iterator). The
slice adds a second filter (matching `mutate.uniform` /
`mutate.s5f` pass names) and a per-event segment switch.

Per-record cost: O(n_event_records × avg_events_per_pass) —
dominated by the existing indel walk for any realistic
configuration. The slice adds a constant factor (~2× the existing
walk).

For high-SHM batches (`mutate(count=50)`), the per-mutate event
count is ~50, which is O(50 events) per record. Negligible vs.
the cost of the mutations themselves.

### No regression risk

The existing indel counter walk is already in the projection
path. Adding the SHM counter walk doesn't change the algorithmic
complexity. Pin via the existing performance baseline
(`tests/test_performance_budgets.py`).

---

## 10. Manifest integration

The manifest's `models.shm` block today carries `available_models`
+ S5F kernel inventory + `segment_rate_support`. A future
extension surfaces the per-segment counter capability:

```python
manifest["models"]["shm"]["per_segment_counters"] = {
    "available": True,
    "fields": ["n_v_mutations", "n_d_mutations", "n_j_mutations", "n_np_mutations"],
    "sum_invariant": "n_mutations",
}
```

Out of scope for this audit; the slice can decide whether the
manifest extension lands in the same PR.

Pinned by `pin_absence_no_per_segment_counters_in_manifest`.

---

## 11. Implementation order (shipped)

One focused slice landed; followed by docs + release-tier
consolidation.

### Slice 1 — Per-segment SHM counter AIRR fields (shipped)

Scope as built — matches the audit's recommendation exactly:

1. **Shipped** — four `i64` fields on `AirrRecord`:
   `n_v_mutations`, `n_d_mutations`, `n_j_mutations`,
   `n_np_mutations`. See
   [`engine_rs/src/airr_record/record.rs`](../engine_rs/src/airr_record/record.rs).
2. **Shipped** — event-ledger aggregation in
   [`engine_rs/src/airr_record/builder.rs`](../engine_rs/src/airr_record/builder.rs)
   filtered to `mutate.{uniform,s5f}` pass names; verbatim
   translation of the indel counter pattern.
3. **Shipped** — five new validator issue kinds in
   [`engine_rs/src/airr_record/validate.rs`](../engine_rs/src/airr_record/validate.rs):
   `NVMutationsMismatch`, `NDMutationsMismatch`,
   `NJMutationsMismatch`, `NNpMutationsMismatch`,
   `MutationCountSumMismatch`. PyO3 serialisation in
   [`engine_rs/src/python/outcome.rs`](../engine_rs/src/python/outcome.rs)
   carries the documented `details.source` tags.
4. **Shipped** — Python `_DEFAULT_COLUMN_ORDER` extended with
   the four columns so TSV / CSV / DataFrame exports surface
   them deterministically.
5. **Deferred** — manifest extension. The slice intentionally
   skipped the optional `per_segment_counters` block — counter
   availability is a property of the engine + AIRR projection,
   not a per-cartridge property; the manifest's
   `models.shm.segment_rate_support` block already advertises
   the targeting capability, and downstream consumers can
   inspect AIRR records directly for counter availability.

### Slice 2 — Consolidation (shipped)

- Doc updates marking this audit shipped + §11 from
  recommended to implemented.
- README mention naming the four new fields + the
  biological-only filter.
- New row in [`docs/validation_matrix.md`](validation_matrix.md)
  pointing at audit + contract + implementation tests + Rust
  kernel sites.
- One release-tier test
  (`test_per_segment_mutation_counters_full_stack_validates_and_partitions`)
  in [`tests/test_release_validation.py`](../tests/test_release_validation.py)
  asserting validator-clean, sum invariant, and corruption
  isolation under the canonical full IGH stack with non-default
  segment rates.

### Out of scope (still deferred)

- **Per-NP1 / per-NP2 split.** Pinned absent by
  `test_pin_absence_no_per_np_split_counters_in_record`. A
  future biology slice can add `n_np1_mutations` /
  `n_np2_mutations` if demand emerges.
- **Per-segment mutation rates.** Derivable from
  `n_<seg>_mutations / <seg>_length`; pinned absent by
  `test_pin_absence_no_per_segment_mutation_rate_fields`.
- **Python `EventRecord` typed-event exposure.** Not needed
  for v1 — Rust aggregates, Python reads AIRR fields. Pinned
  absent by
  `test_pin_absence_python_eventrecord_lacks_typed_events`.
- **Mutation-ledger validator** (per-site mutation-event vs
  pool-diff cross-check). Out of scope for v1; structural
  protection via the `MutationTransaction` boundary lockdown
  test remains the only static guard.
- **Indel-style "biological" reclassification of polymerase
  indels** (audit §2). SHM indels remain unmodelled; current
  classification (substitutions only) holds.
- **Per-pass cost / counter for any other observation-stage
  pass.** Out of scope.

This audit only proposes the architecture and pins absences;
the implementation slice flips the relevant `pin_absence_*` to
`pin_present_*` in lockstep.

---

## 12. Test surface — what this audit pins

Mirrored in
[`tests/test_mutation_provenance_contract.py`](../tests/test_mutation_provenance_contract.py).

### `pin_scaffold_*` — today's contract

1. `BaseChanged` carries `segment` + `germline_pos` (source-level
   pin on the IR event definition).
2. Only `mutate.uniform` / `mutate.s5f` passes emit `BaseChanged`
   events that should count toward `n_mutations` (behavioural —
   walk the ledger after a known config).
3. Pass-name constants (`MUTATE_UNIFORM`, `MUTATE_S5F`,
   `CORRUPT_INDEL`, `CORRUPT_PCR`, `CORRUPT_QUALITY`) exist with
   their canonical values.
4. Existing global counters (`n_mutations`, `n_indels`,
   `n_v_indels`, `n_d_indels`, `n_j_indels`, `n_pcr_errors`,
   `n_quality_errors`) all present.
5. `n_indels` includes NP1/NP2 events; `n_v/d/j_indels` don't
   (the precedent for the proposed `n_np_mutations` rollup).
6. `n_pcr_errors` / `n_quality_errors` are trace-sourced
   (attempted count semantics, not realised).
7. Receptor revision, D inversion, contaminant DO emit
   IR/trace state changes but do NOT bump `n_mutations`.
8. Segment-rate zero-rate exclusion: with `{"NP": 0.0}`,
   `n_mutations` counts only V/D/J mutations (the future
   `n_np_mutations` would be 0).

### `pin_absence_*` — gaps Slice 1 closes

9. No `n_v_mutations` / `n_d_mutations` / `n_j_mutations` /
   `n_np_mutations` AIRR fields.
10. No `*MutationCountMismatch` / `MutationCountSumMismatch`
    validator issue kinds.
11. No `per_segment_counters` field in
    `cartridge_manifest()["models"]["shm"]`.
12. Python `EventRecord` doesn't expose typed event list (only
    `simulation_event_count`).

### Doc anchor

13. The audit doc continues to exist and references the contract
    file; the 14-section structure stays intact.

---

## 13. Out of scope

Documented here so a future contributor doesn't accidentally
expand the slice.

- **Per-NP1 / per-NP2 split.** v1 collapses both into a single
  `n_np_mutations`. A future slice can add `n_np1_mutations` /
  `n_np2_mutations` if demand emerges.
- **Per-segment mutation rates** (`v_mutation_rate`, etc.). Out
  — derivable from `n_<seg>_mutations / <seg>_length`.
- **Python typed-event exposure.** The slice doesn't need it;
  separate observability decision.
- **Reclassifying polymerase indels as biological.** SHM indels
  are a real biology gap, but `polymerase_indels` is library-
  prep; the audit's recommendation is to leave indels OUT of the
  SHM counter rollup until a dedicated SHM-indel pass exists.
- **CDR/FR sub-segment counters.** Within-V positional
  resolution is a future biology slice.
- **Per-pass-instance counters.** Multiple `mutate()` calls
  aren't supported today; the audit's aggregation is over the
  union of all `mutate.{uniform,s5f}` events.
- **Cross-cartridge / cross-locus counter harmonisation.** TCR
  pipelines have no SHM and counters stay at 0; no special
  handling.

---

## 14. Summary table

| Concern | Answer |
|---|---|
| Current global mutation counter | `n_mutations` (i64) — IR-sourced from `sim.mutation_count`; bumped only by `add_to_mutation_count` calls in `UniformMutationPass` + `S5FMutationPass`. |
| Other current counters | `n_indels` + `n_v/d/j_indels` (event-sourced, indel-pass-filtered); `n_pcr_errors`, `n_quality_errors`, `end_loss_*_length` (trace-sourced, attempted counts). |
| What counts as biological mutation | `BaseChanged` events from passes whose `pass_name ∈ {mutate.uniform, mutate.s5f}`. Receptor revision, D inversion, PCR, quality, N-corruption, contaminant, indels NOT included. |
| Segment attribution rule | `BaseChanged.segment` carried at event-emission time. Unambiguous (segment doesn't change under substitution). Stable across post-event indels / end-loss. |
| Events carry enough data? | **Yes.** `BaseChanged { segment, germline_pos, ... }` already exists. No new IR field needed. |
| Recommended counter shape | Four new AIRR i64 fields: `n_v_mutations`, `n_d_mutations`, `n_j_mutations`, `n_np_mutations` (NP1+NP2 rolled). Sum equals `n_mutations`. |
| Recommended validator additions | Five issue kinds: `V/D/J/NpMutationCountMismatch` + `MutationCountSumMismatch`. Source tag `events:mutate.{uniform,s5f}:base_changed:{V,D,J,NP}`. |
| Trace / replay impact | None. Counters are event-derived; existing plan + refdata signatures cover the slice's reproducibility. |
| Performance | O(n_event_records × events_per_pass), constant factor over the existing indel walk. Negligible. |
| Manifest extension | Optional: `models.shm.per_segment_counters = {available, fields, sum_invariant}`. Slice author decides. |
| Python EventRecord typed-event exposure | NOT needed for v1 — Rust aggregates, Python reads AIRR fields. Documented Python observability gap. |
| Required for slice to land | Just `airr_record/builder.rs` + `validate.rs` + `record.rs` field additions. Roughly 50–80 lines of Rust. |
| Pre-flight bugs found | **None.** Audit clears the slice for implementation. |

The per-segment SHM counter slice is **architecturally trivial**:
the event ledger already carries everything needed, the indel
counter is the canonical reference implementation, and the
recommended four-field shape composes cleanly with both the
segment-rates DSL and the existing global `n_mutations`. The
slice can be a focused ~50-line Rust diff + four Python column-
order entries.

After this audit ships, the implementation slice should close
the "users can target SHM by region but can't measure where it
landed" loop the targeted-SHM slice opened.
