# Paired-End / Read Layout — Design + Pre-Implementation Audit

**Status: shipped.** Slices A–D + the release-consolidation slice
are merged; the user-facing surface is
[`Experiment.paired_end(r1_length=…, r2_length=…, insert_size=…)`](../src/GenAIRR/experiment.py)
and the AIRR record now carries eight paired-end fields gated by
the Slice B projection kernel + Slice C trace integration. Per-slice
deliverables:

| Slice | What landed | Verified by |
|---|---|---|
| A — Record schema + default validator | Eight `AirrRecord` fields (`read_layout` / `r1_sequence` / `r2_sequence` / `r1_start/end` / `r2_start/end` / `insert_size`) with the additive-default precedent. Five `RecordValidationIssue` variants (`PairedEndFieldWithoutLayout` enforced; four `Read*` variants reserved for Slice B). | `engine_rs/src/airr_record/tests/projection.rs` (5 Slice A tests) + [`tests/test_paired_end_schema.py`](../tests/test_paired_end_schema.py) (7 tests) |
| B — Projection kernel + geometry validator | `project_paired_end_layout(rec, r1_length, r2_length, insert_size) -> PairedEndProjection`, validator dispatches on `read_layout` and surfaces the four reserved variants for `"paired_end"` records. `"single_end"` reserved no-op; unknown values → `ReadLayoutMismatch`. | `engine_rs/src/airr_record/sequence.rs` + `engine_rs/src/airr_record/tests/projection.rs::slice_b_*` (17 Slice B tests covering kernel correctness + every geometry tamper) |
| C — Trace addresses + builder integration | Three new `ChoiceAddress` variants (`PairedEndR1Length` / `PairedEndR2Length` / `PairedEndInsertSize`), `PairedEndLayoutSpec` (three int distributions), `PairedEndSamplingPass` (trace-only, no IR mutation), AIRR builder reads trace + applies kernel | `engine_rs/src/passes/paired_end.rs::tests` (12 Slice C tests: fresh, replay, missing/wrong address, bad geometry, baseline) + `airr_record::tests::projection::slice_c_*` (5 builder integration tests) |
| D — DSL | `Experiment.paired_end(*, r1_length, r2_length=None, insert_size)` with three input shapes (int / `(low, high)` / empirical pairs), at-most-once + fixed-value geometry guards, end-of-plan lowering | [`tests/test_paired_end_dsl.py`](../tests/test_paired_end_dsl.py) (20 e2e tests) |
| Consolidation | README row, validation-matrix row, release-tier IGH stack composed with receptor revision + D inversion + mutation + end-loss + rev-comp, full-stack replay round-trip with variable insert sizes, distribution invariant on `paired_end.insert_size` (uniform-int mean ±5σ + endpoint bounds) | [`tests/test_release_validation.py`](../tests/test_release_validation.py), [`tests/test_distribution_invariants.py`](../tests/test_distribution_invariants.py) |

The contract test suite [`tests/test_paired_end_contract.py`](../tests/test_paired_end_contract.py)
closed the arc with every Slice-gated `pin_absence_*` test
flipped to `pin_scaffold_*`. Out-of-scope pins (no FASTQ
exporter, no two-row per-fragment output, no `.single_end()` DSL
method, no per-read quality strings) remain absences per §11 of
this document — and are intentionally left as `pin_absence_*` so
a future contributor who lands one of those features surfaces it
through an explicit lockstep flip.

---

A read-only architecture review of every engine surface that
paired-end read-layout projection will touch, written **before**
any implementation. Companion to
[`tests/test_paired_end_contract.py`](../tests/test_paired_end_contract.py)
which pins today's absence + the existing scaffolding and will
flip to behavior-pinning as implementation slices land.

Paired-end is the dominant Illumina output layout in real
immune-seq libraries: each sequenced fragment yields two reads
(R1 reading forward from the 5' adapter, R2 reading
reverse-complemented from the 3' adapter), and an *insert size*
that is typically larger than 2 × read_length so the two reads do
not overlap. Modelling paired-end output is the natural follow-up
to the [end-loss / primer-trim audit](primer_trim_end_loss_audit.md):
both are observation-layer transforms applied to an already-
assembled biological molecule.

Unlike receptor revision (which mutated an already-assembled
biological truth), paired-end is **projection-only**: it
introduces no new IR mutation, no live-call invalidation, and no
contract narrowing. The whole mechanism lives in the AIRR
projection layer + a thin trace surface for the per-record
geometry choices.

---

## Existing scaffolding already in the engine

Surfaces today's pipeline already provides, which the design
relies on. Each is pinned by a `pin_scaffold_*` test in the
companion contract file.

| Existing surface | Where | What it gives the paired-end design |
|---|---|---|
| `AirrRecord.sequence` is the **single canonical post-pipeline molecule string** | `engine_rs/src/airr_record/record.rs:11` | The substrate every paired-end window indexes into. Paired-end never reshapes the molecule; it only carves windows out of `sequence`. |
| Watson-Crick `reverse_complement` helper | `engine_rs/src/airr_record/sequence.rs:59` (private) and `complement_base` byte-level kernel at `engine_rs/src/ir/nucleotide.rs:26` | The deterministic transform R2 needs. Both helpers are already exercised by the D-inversion (`InvertDPass`) and rev-comp projection paths. |
| `apply_rev_comp_projection` (post-build flip of `sequence` / coords / AA) | `engine_rs/src/airr_record/sequence.rs:14` | Pins the contract that all "strand-orientation" transforms apply **post-build, after every other field is finalised**. Paired-end window selection sits in this same post-build phase. |
| `end_loss_5_length` / `end_loss_3_length` fields | `engine_rs/src/airr_record/record.rs:116` | Demonstrates the additive-AIRR-field pattern (defaults to 0; populated only when the matching pass ran). Paired-end fields will follow the same pattern. |
| `random_strand_orientation()` DSL surface (alias for `corrupt.rev_comp.applied`) | `src/GenAIRR/experiment.py:488` | The user-facing knob for whole-record orientation. Paired-end must compose with it cleanly (see §7). |
| `CorruptEndLoss(PrimeEnd::Five)` / `(Three)` trace addresses + `EndLossPass` | `engine_rs/src/passes/corrupt/end_loss.rs` | Establishes the observation-stage pattern: trace records per-record per-side amounts; projection respects them. Paired-end will mirror this pattern for `paired_end.insert_size` / `paired_end.r1_length` / `paired_end.r2_length`. |
| AIRR validator's structured `details.source` surface | `engine_rs/src/python/outcome.rs::issue_to_pydict` | Lets a future `PairedEndWindowMismatch` issue land with `details.source: "trace:paired_end.insert_size"` etc., mirroring `OriginalVCallMismatch`'s shape. |

---

## 1. Current read model

Today's user-facing model is **one molecule, one AIRR row, one
read**:

- One `Outcome` → one record in `SimulationResult.records`.
- `AirrRecord.sequence` is the full pool projection of the
  post-pipeline `Simulation` (after recombine, mutation,
  corruption, end-loss, and rev-comp projection).
- `random_strand_orientation(prob=…)` flips the projection
  surface — `sequence`, `np1`, `np2`, `junction`, and coords —
  via `apply_rev_comp_projection`. It does **not** touch the IR
  pool.
- `end_loss_5prime` / `end_loss_3prime` delete pool bytes
  *before* projection. The molecule the projection observes is
  already shortened; `end_loss_5_length` / `end_loss_3_length`
  record what was removed.
- No separate read objects. No fragment / insert concept. No
  R1 / R2 windows.

The "molecule" the AIRR record reports is the substrate downstream
aligners would consume as if it were a single Sanger-style read
covering the entire receptor.

---

## 2. Output shape

**v1 keeps one AIRR row per simulated molecule** and adds
read-layout fields to that row:

| Field | Type | Default when no `.paired_end()` step ran |
|---|---|---|
| `r1_sequence` | `String` | `""` |
| `r2_sequence` | `String` | `""` |
| `r1_start` | `Option<i64>` | `None` |
| `r1_end` | `Option<i64>` | `None` |
| `r2_start` | `Option<i64>` | `None` |
| `r2_end` | `Option<i64>` | `None` |
| `insert_size` | `Option<i64>` | `None` |
| `read_layout` | `String` | `""` (literal enum: `""` / `"single_end"` / `"paired_end"`) |

### 2.1 Rejected: two AIRR rows per molecule

Real paired-end output is two FASTQ entries (`_R1.fq`, `_R2.fq`)
sharing a fragment id. The naïve translation is "emit two AIRR
rows per simulation, joined on a shared `sequence_id`."

That approach was considered and rejected for v1:

- **Trace grouping breaks.** Every existing AIRR row maps
  one-to-one to a `Trace` / `EventRecord`. Splitting one
  simulation across two rows requires a parent-child trace ID
  scheme that doesn't exist today; every downstream consumer
  (replay, validator, clonal-fork) assumes the one-to-one map.
- **Clone IDs and clonal_fork compose unpredictably.** A
  per-clone simulation already produces N descendant records;
  paired-end would double that, and the AIRR spec doesn't have
  a clean "two reads, one descendant" shape without inventing
  a new `read_index` column.
- **AIRR-spec semantics.** AIRR Rearrangement v1 is fundamentally
  one row per recombination event. Two-row paired-end is more
  natural in MiAIRR's `RepairableRead` extension, which GenAIRR
  doesn't model.
- **Validation explosion.** Every existing per-record
  invariant (`v_call`, `junction`, contracts) would need to
  decide whether R1's row and R2's row should both pass, or one
  is a "phantom" entry. The single-row shape sidesteps the
  question.

The one-row shape preserves the existing
{outcome → trace → AIRR record → validator} chain. A future
exporter (`SimulationResult.to_paired_fastq(...)`) can split the
single row into two FASTQs at export time without complicating the
internal model. That's a downstream concern, not a v1 audit
question.

### 2.2 `read_layout` as an enum string

Rather than a separate `is_paired_end: bool`, the design uses one
enum-style string. Values:

- `""` — no paired-end step ran (the existing default; baseline
  records carry the empty string).
- `"single_end"` — a future `.single_end(read_length=…)` step
  could land an explicit single-read window without paired
  geometry. Reserved.
- `"paired_end"` — both R1 and R2 windows are populated.

Single field + enum keeps the projection layer's compatibility
straightforward: a future `to_paired_fastq` reads
`read_layout == "paired_end"` to decide whether the row produces
one or two FASTQ entries.

---

## 3. DSL

**Proposed v1 signature:**

```python
Experiment.paired_end(
    read_length: int = 150,
    insert_size: Union[int, Tuple[int, int], Sequence[Tuple[int, float]]] = (250, 450),
)
```

Rules (mirror `invert_d` / `receptor_revision`):

- **VDJ + VJ both supported.** Paired-end is a sequencing-stage
  observable, not a biology mechanism; it makes sense on every
  chain.
- **At most once per pipeline.** A second call raises
  `ValueError` (no last-wins).
- **Must be the last appended step.** Paired-end is a projection
  option; appending it before a mutation/corruption step is a
  user error (the windows would be drawn against the
  pre-mutation molecule). v1 raises `ValueError` if the user
  tries to chain another step after `paired_end`.
  - The lowering enforces this by inlining the projection
    metadata into a per-record post-build phase, *not* a regular
    `PassPlan` push (see §5).

`insert_size` shapes:

- **Single integer** (e.g. `insert_size=300`) → fixed insert size.
- **Tuple of two ints** (e.g. `(250, 450)`) → uniform draw from
  the closed interval `[min, max]`. This is the recommended
  shape for v1; matches the existing `EndLossPass` length-pair
  convention.
- **Sequence of (insert_size, weight) pairs** → empirical
  distribution. Optional; matches the `length_pairs` shape.

### 3.1 Future ergonomic extension (out of v1 scope)

```python
.paired_end(r1_length=150, r2_length=140, insert_size=(...))
```

Asymmetric R1/R2 lengths are common in real sequencing runs (R2
quality drops faster, so libraries are often `2 × 150` with R2
hard-trimmed to 140). v1 ships symmetric `read_length` only; the
constructor signature reserves `r1_length` / `r2_length` for the
follow-up slice. Mixing `read_length` with either asymmetric
keyword raises `ValueError`.

### 3.2 Rejected: `.paired_end_from_real_data(distribution_file=…)`

Loading empirical insert distributions from a JSON file is a
distribution-provider concern, not a DSL concern. The empirical-
pairs `insert_size` form covers the use case; a downstream helper
can load the file into pairs.

---

## 4. Trace addresses (minimal replayable choice set)

Three new addresses, under the `paired_end.*` namespace:

| Address | Type | Records |
|---|---|---|
| `paired_end.insert_size` | `Int` | The sampled insert size for this simulation. |
| `paired_end.r1_length` | `Int` | The R1 read length. |
| `paired_end.r2_length` | `Int` | The R2 read length. |

### 4.1 Why record fixed-constant lengths too

If the user passes `read_length=150` (a constant, not a
distribution), the recorded `r1_length` and `r2_length` are both
150 on every record. Recording them anyway costs ~16 bytes per
trace and gives one important guarantee: **a trace fully explains
the read layout it produced**. A downstream consumer reading the
trace alone (no Experiment object, no original `Experiment`
config) can reconstruct the exact R1/R2 substrings the AIRR
record carries. This matches the existing convention for
constant `apply_prob` values (`corrupt.contaminant.applied` is
always recorded, even when `apply_prob=0.0`).

### 4.2 Record order

For each `Outcome`, the three records are emitted in this order:

1. `paired_end.r1_length`
2. `paired_end.r2_length`
3. `paired_end.insert_size`

The lengths land first because in the asymmetric-future case they
constrain the insert-size support (`insert_size >= r1 + r2`); the
ordering keeps replay validation linear.

### 4.3 Address-schema-version policy

Additive (per the
[`address.rs` top-of-file docs](../engine_rs/src/address.rs)). Old
trace files don't reference `paired_end.*`; new traces parse on
engines that postdate this slice. No `ADDRESS_SCHEMA_VERSION`
bump.

---

## 5. IR or projection? Architectural choice.

**Paired-end is projection / read-layout only, not IR mutation.**

Concrete consequences of this commitment:

| Behaviour | What this design forbids |
|---|---|
| Paired-end **must not** call `builder.delete_indel` / `with_indel_*` / `replace_segment` / any persistent IR mutator. |
| Paired-end **must not** push a `Pass` that lowers into the engine schedule. The projection metadata is threaded through `PassContext.event_log_sink` at most (and even that may not be needed — see §5.2). |
| Paired-end **must not** invalidate live-call evidence. `SegmentLiveCall`s computed against the molecule remain valid; the R1/R2 windows are substrings the validator re-derives from `rec.sequence` after the live-call layer has settled. |
| Paired-end **must not** narrow contracts. Productive / junction-frame / anchor-preserved all run on the molecule before paired-end window selection; the molecule's productivity is independent of where the sequencer reads it. |

### 5.1 Why this is the right call

- **Observation, not biology.** A real Illumina sequencer
  observes a fragment; it doesn't edit it. Modelling paired-end
  as projection matches the biological reality.
- **No live-call invalidation.** Receptor revision's Slice B
  taught us that an IR-mutating projection-like mechanism pulls
  in `LiveCallRefreshHook` complexity. Avoiding that for
  paired-end keeps the change surface narrow.
- **Compatibility.** Users who don't call `.paired_end()` see
  zero behaviour change; the AIRR record gains empty/default
  fields per §10.
- **Validator simplicity.** All paired-end checks reduce to
  substring equality + arithmetic on `rec.sequence` — no
  separate truth oracle, no event-stream re-derivation.

### 5.2 Where the projection metadata lives

Two possible homes for the post-build R1/R2 window write:

1. **`build_airr_record`** reads the trace, draws the
   per-record geometry, and writes the new fields directly.
   No new pass at all.
2. **A thin `PairedEndPass`** that runs at seal time, samples
   the geometry into the trace, and exposes a public hook the
   builder reads.

The audit recommends **option 1** for v1: the geometry is so
simple it's clearer as a builder helper than as a one-line pass.
The trace records are written through a small `PairedEndPlan`
struct the Experiment lowering attaches to the `Outcome` (the same
shape `Outcome.trace_file_from` already pulls from). Slice C
revisits this if a future pass-based shape becomes natural.

---

## 6. Relationship to end-loss

End-loss and paired-end are **two separate observation-stage
transforms applied in a fixed order**:

```text
…recombination → mutation → corruption → end_loss → paired_end projection
```

| Stage | What happens to the molecule |
|---|---|
| `end_loss_5/3prime` | **Pool bytes deleted.** The molecule the AIRR record sees is shorter than the assembled receptor. `end_loss_5_length` / `end_loss_3_length` record the lost amounts. |
| `paired_end` | **Pool bytes preserved.** The molecule is unchanged; R1 / R2 are *windows* into the existing `rec.sequence`. `insert_size` is the position-from-5' of the R2-window end (the "fragment 3' end"); the unobserved middle region between R1 and R2 is the inner mate-pair gap. |

### 6.1 Insert-size geometry on a post-end-loss molecule

`insert_size` is measured **against the final, post-end-loss
`rec.sequence_length`**. A 200-bp end-loss on a 600-bp molecule
followed by `paired_end(insert_size=350)` reads:

- `sequence_length = 400`
- R1 window: `sequence[0:read_length]`
- R2 window: `sequence[insert_size - read_length : insert_size]`
- Insert ends at position 350 (within the 400-bp molecule).

Insert sizes **larger than the post-end-loss molecule length**
raise `ValueError` at projection time. The user is asking the
sequencer to read past the end of the fragment; v1 fails loudly.

A future slice could clamp + record a `read_layout_truncated:
bool` flag instead. v1 keeps the strict behaviour for
correctness-first.

### 6.2 R1/R2 lengths longer than the insert

If `r1_length + r2_length > insert_size`, the two reads overlap
in the middle. v1 supports this — the two windows are still drawn
independently, and downstream consumers can decide whether to
collapse the overlap. The AIRR record carries both reads
verbatim.

Only the degenerate case `r1_length > insert_size` (R1 alone runs
past the fragment 3' end) raises `ValueError`.

---

## 7. Reverse complement: R2 representation + composition with `random_strand_orientation`

Standard Illumina paired-end layout: R1 reads forward from the
fragment 5'; R2 reads reverse-complemented from the fragment 3'.

```
Fragment (5' → 3'):   AAAACCCCGGGGTTTT...

R1 (forward, 5'):     AAAACCCCGGGG
                      ^
                      sequence[0:r1_length]

R2 (reverse-complement, 3'):       TTTT...
                                   ^                 ^
                                   sequence[i-r2_len : i]
                                   then reverse_complement applied
```

**Decision: paired-end window selection runs AFTER any
`random_strand_orientation` projection has been applied to the
molecule.**

Concrete consequences:

- `rec.sequence` is the final user-visible molecule. If
  `corrupt.rev_comp.applied = Bool(true)`, `rec.sequence` is
  already the reverse-complement of the assembled pool.
- R1 is `rec.sequence[r1_start:r1_end]` — forward in the
  *projected* orientation.
- R2 is `reverse_complement(rec.sequence[r2_start:r2_end])` —
  reverse-complemented from the *projected* 3' end.
- `rec.rev_comp` is unchanged in meaning: the molecule itself
  was emitted antisense. R1 / R2 are *windows into that
  antisense molecule*, not orientation-flipped twice.

### 7.1 Why "after rev-comp projection"

The alternative — apply paired-end on the pre-rev-comp
molecule, then flip R1/R2 + swap them — was considered and
rejected:

- **Downstream consumer expectation.** Aligners receive R1/R2 as
  pairs against the reference; they don't track whether the
  source molecule was rev-comped. If R1/R2 are flipped before
  rev-comp, the aligner sees the same bytes either way; the
  field semantics become harder to explain.
- **AIRR field consistency.** `rec.sequence` is already
  rev-comped post-build; `r1_sequence` being a substring of it
  composes obviously.
- **Reverse-complement composition.** Three rev-comp transforms
  in sequence (build → rev-comp molecule → carve windows
  → rev-comp R2) is easier to audit than three transforms in a
  different order.

The validator pins this composition via §8.

---

## 8. Validator

Two new `RecordValidationIssue` variants (mirroring the
`receptor_revision_applied` / `original_v_call` shape from Slice E
of receptor revision):

```rust
PairedEndWindowMismatch {
    side: PairedEndRead,          // R1 or R2
    reported: String,             // record.r1_sequence / r2_sequence
    expected: String,             // recomputed from rec.sequence[start..end]
}

PairedEndGeometryMismatch {
    field: PairedEndGeometryField, // R1Start / R1End / R2Start / R2End / InsertSize / ReadLayout
    reported: String,              // stringified for serialization
    expected: String,
}
```

Where:

```rust
pub enum PairedEndRead { R1, R2 }

pub enum PairedEndGeometryField {
    R1Start,
    R1End,
    R2Start,
    R2End,
    InsertSize,
    ReadLayout,
}
```

### 8.1 Checks the validator runs

For every record (regardless of whether a `.paired_end()` step
ran):

| Check | When `read_layout == ""` | When `read_layout == "paired_end"` |
|---|---|---|
| All paired-end fields default | All-or-nothing: every field must be empty / `None`. Else surface as `PairedEndGeometryMismatch { field: ReadLayout, ... }`. | (skip) |
| `r1_start`, `r1_end`, `r2_start`, `r2_end` in `[0, sequence_length]` and `r1_start <= r1_end`, `r2_start <= r2_end` | (skip — fields are `None`) | Out-of-bounds or inverted → `PairedEndGeometryMismatch`. |
| `r1_sequence == rec.sequence[r1_start:r1_end]` byte-for-byte | (skip) | Mismatch → `PairedEndWindowMismatch { side: R1, ... }`. |
| `r2_sequence == reverse_complement(rec.sequence[r2_start:r2_end])` | (skip) | Mismatch → `PairedEndWindowMismatch { side: R2, ... }`. |
| `insert_size == r2_end` (insert-size is the position-from-5' of the R2 window end) | (skip) | Mismatch → `PairedEndGeometryMismatch { field: InsertSize, ... }`. |
| `paired_end.insert_size` / `r1_length` / `r2_length` trace records match the record's geometry | (skip — addresses absent in trace) | Mismatch → `PairedEndGeometryMismatch { source: trace:paired_end.* }`. |

### 8.2 Structured `details.source` strings

Following the `OriginalVCallMismatch` precedent from receptor
revision Slice E. Strings pinned in the contract file so MCP /
dashboard consumers can match by prefix:

| Issue variant | `details.source` |
|---|---|
| `PairedEndWindowMismatch { side: R1 }` | `"projection:sequence[r1_start:r1_end]"` |
| `PairedEndWindowMismatch { side: R2 }` | `"projection:reverse_complement(sequence[r2_start:r2_end])"` |
| `PairedEndGeometryMismatch { field: InsertSize }` | `"trace:paired_end.insert_size"` |
| `PairedEndGeometryMismatch { field: R1Start/R1End }` | `"trace:paired_end.r1_length"` |
| `PairedEndGeometryMismatch { field: R2Start/R2End }` | `"trace:paired_end.r2_length"` |
| `PairedEndGeometryMismatch { field: ReadLayout }` | `"projection:read_layout"` |

---

## 9. Replay

Replay must reproduce R1 / R2 bytes exactly. The mechanism:

1. The trace carries three `Int` records:
   `paired_end.r1_length`, `paired_end.r2_length`,
   `paired_end.insert_size`.
2. Replay consumes them in order via `expect_int`.
3. The post-replay AIRR builder reads the same records and
   re-derives the windows from `rec.sequence` — which is
   itself the rev-comp-stable projection of the pool the trace
   replay reconstructed.
4. Therefore R1 / R2 round-trip bit-for-bit.

Replay errors that can fire:

- Trace records out-of-bounds for the post-replay
  `sequence_length` (refdata swap shrunk the molecule):
  `PassError::InvalidPlanState`.
- Trace records present but `read_layout = ""` in the replay
  config (a hand-crafted trace that contradicts the
  Experiment shape): `PassError::InvalidPlanState`.

The DSL boundary already enforces "paired_end is the last step,"
so a replay cursor that sees the three `paired_end.*` records can
expect them at a known position in the consume sequence.

---

## 10. Compatibility

**Existing users see zero behaviour change.** Concrete pins
(every one tested by `pin_scaffold_*` in the contract file):

- Baseline runs without `.paired_end()` keep emitting one AIRR
  row per outcome (unchanged).
- The new fields default to `""` / `None` / `0` / `""` for
  records that didn't go through paired-end projection.
- The validator's existing checks (sequence content, junction,
  v_call, …) are unchanged for records that don't carry paired-
  end fields.
- The trace's existing addresses are unchanged; `paired_end.*`
  appears only when the user calls the DSL method.
- `random_strand_orientation` still flips `rec.sequence` etc.
  for the whole molecule; paired-end composes with it per §7.
- `end_loss_5/3prime` still shorten the molecule before paired-
  end window selection (per §6).
- `content_hash` stability: adding `.paired_end()` to a plan
  whose record path doesn't actually carry paired-end fields
  (a defensive future Experiment that lowers conditionally) is
  not a v1 concern; v1 always populates the fields when the DSL
  method is called.

**Trace-file schema:** additive. The three new addresses don't
collide with existing prefixes; `ADDRESS_SCHEMA_VERSION` stays
at the current value.

**`AirrRecord` struct layout:** the new fields are appended after
`original_v_call` (the last Slice E field). Downstream Rust
consumers that pattern-match on the struct destructively will
need to add the new fields; consumers that read fields by name
are unaffected.

**Python column list:** `r1_sequence`, `r2_sequence`,
`r1_start`, `r1_end`, `r2_start`, `r2_end`, `insert_size`,
`read_layout` appended to `result._DEFAULT_COLUMN_ORDER` in the
same order as the Rust struct. Dataframe / CSV exports include
the new columns automatically.

---

## 11. Out of scope

Documented here so a future contributor doesn't accidentally
expand the slice:

- **FASTQ / two-row export.** A downstream
  `SimulationResult.to_paired_fastq()` exporter is a useful
  follow-up but is not part of the v1 paired-end audit. It can
  be added without changing the AIRR record schema.
- **Quality scores per base.** Real FASTQ carries per-base
  quality. v1 leaves quality to the existing
  `corrupt.quality` pass on the molecule; paired-end windows
  inherit those quality positions when an exporter later
  materialises them.
- **Adapter trimming on R1 / R2.** Real libraries can leave
  adapter bases at the end of short fragments. v1 reads the
  exact `read_length` bytes regardless; adapter sim is a
  separate observation-layer concern.
- **Variable per-read length within one library.** v1 records
  fixed r1/r2 lengths per pipeline. A future slice could draw
  from a distribution per record.
- **Multi-read layouts** (e.g. 10x Genomics single-cell:
  R1=barcode+UMI, R2=cDNA). Out of scope; would need a
  separate `read_layout = "10x"` arm.
- **Paired-end + `clonal_fork`** interaction. v1 supports
  paired-end on every descendant of a clonal fork (the
  projection runs per-record after fork expansion). No
  per-clone shared insert sampling; each descendant samples
  independently. Documented here so a future "shared per-clone
  layout" request can scope it cleanly.
- **Strand bias** (R1/R2 forward orientation correlated with
  pool's `rec.rev_comp` flag). v1 doesn't model strand-
  specific library preparation; downstream consumers can
  filter on `rec.rev_comp` if needed.

---

## 12. Implementation order

Recommended slice sequence; each independently revertible and
gated on a green test sweep before the next:

1. **Slice A — Record/schema fields + validator defaults.**
   Add the eight new fields to `AirrRecord`. Default to empty /
   `None`. Wire PyO3 dict + Python `result.py` column list.
   Validator default behaviour: when `read_layout == ""`, the
   eight fields must be at their defaults. No DSL yet, no
   projection logic yet. Audit contract file's `pin_absence_*`
   tests for AIRR fields flip to scaffolding.
2. **Slice B — Read-layout projection helper.** Add the
   `paired_end_window_projection(rec, r1_length, r2_length,
   insert_size)` helper in `airr_record/sequence.rs`. Unit
   tests cover §6 / §7 / §8 cases against synthetic
   `AirrRecord`s. Helper not wired to any production caller
   yet.
3. **Slice C — Trace addresses + replay path.** Three new
   `ChoiceAddress` variants + frozen-spelling tests. A new
   `PairedEndPlan` struct on `Outcome` (sized to hold the
   three `Int`s). `paired_end_window_projection` is invoked
   from the builder when the plan is `Some`. Replay consumes
   the three records.
4. **Slice D — DSL.** `Experiment.paired_end(read_length=...,
   insert_size=...)`. Last-step guard, at-most-once guard,
   prob-validation chain. The lowering attaches a
   `PairedEndPlan` to the compiled simulator rather than
   pushing a Pass (per §5.2).
5. **Slice E — Release consolidation + distribution invariant.**
   README row, validation-matrix row, release-tier test, full-
   stack replay round-trip, Bernoulli-equivalent (insert-size
   distribution within 5σ) test. Same shape as the receptor-
   revision consolidation.

Each slice ≤ one PR-sized chunk. The receptor-revision arc took
five slices to land in the same shape.

---

## 13. Test surface — what the implementation slices must pin

A non-exhaustive list, mirrored as TODO markers in the companion
[`tests/test_paired_end_contract.py`](../tests/test_paired_end_contract.py):

1. **Field defaults.** Baseline records (no `.paired_end()`)
   carry the eight fields at their defaults. Pinned per-field.
2. **Window byte-for-byte equality.**
   `rec.r1_sequence == rec.sequence[rec.r1_start:rec.r1_end]`;
   `rec.r2_sequence == reverse_complement(rec.sequence[rec.r2_start:rec.r2_end])`.
3. **Insert-size geometry.** `rec.insert_size == rec.r2_end`.
4. **Composition with `random_strand_orientation`.** When the
   molecule is rev-comped, R1 reads from the projected 5' end
   (not the pre-rev-comp 5' end). Trace replay reproduces the
   composition.
5. **Composition with `end_loss`.** R1/R2 windows index into
   the post-end-loss `rec.sequence`, never into the pre-end-loss
   pool.
6. **Trace replay.** A trace file emitted with
   `paired_end(read_length=150, insert_size=(250, 450))`
   replays bit-for-bit, including the three `Int` choices and
   the R1/R2 substrings.
7. **VJ + VDJ both supported.** The DSL accepts both chain
   types; the windows index into the molecule regardless of
   chain shape.
8. **At-most-once + last-step DSL guards.** Calling
   `.paired_end()` twice or chaining a step after it raises
   `ValueError`.
9. **Validator tampering surfaces structured issues.**
   Each of the six `PairedEndGeometryField` arms is reachable
   by manually editing one field after the canonical builder
   ran.
10. **Distribution invariant.** With
    `insert_size=(250, 450)` over N=4000 seeds, the empirical
    insert-size distribution is uniform within 5σ per bucket
    (matches the existing `EndLossPass` length-pair
    invariant pattern in
    `tests/test_distribution_invariants.py`).
11. **Content-hash stability.** Adding `.paired_end()` to a
    plan that doesn't read paired-end fields doesn't change
    existing test fixtures' AIRR output.

---

## 14. Out of scope (cross-reference)

See §11 for the full out-of-scope list. The cross-reference
exists so a future audit can find the "what we *didn't* model"
catalogue in the same section number as receptor revision.

---

## Summary table

| Question | Decision |
| --- | --- |
| Scope | Read-layout projection only; no IR mutation, no live-call invalidation, no contract narrowing. |
| Output shape | One AIRR row per simulated molecule; eight new fields default empty / `None`. |
| DSL | `Experiment.paired_end(read_length=150, insert_size=(250, 450))`. Both VJ + VDJ. At most once. Must be the last appended step. |
| Pipeline position | Final projection option; runs after every IR-mutating + observation-stage pass (recombine, mutation, corruption, end-loss). |
| Trace | `paired_end.r1_length` (Int), `.r2_length` (Int), `.insert_size` (Int). Schema-version unchanged. |
| IR carrier | None. Projection-layer only — no `Pass`, no `SimulationEvent`. |
| Event | None. Paired-end fires no `SimulationEvent`s. |
| Live-call refresh | Unaffected. R1/R2 windows are substrings of `rec.sequence` after live calls have finalised. |
| AIRR | Add `r1_sequence`, `r2_sequence`, `r1_start/end`, `r2_start/end`, `insert_size`, `read_layout`. |
| Contracts | Unaffected. Paired-end runs after every contract has settled. |
| Replay | Three-record consume via `expect_int`; absent records → `read_layout = ""` (back-compat). |
| Validator | New `PairedEndWindowMismatch` + `PairedEndGeometryMismatch` issues with structured `details.source` strings. |
| `random_strand_orientation` interaction | Paired-end windows are drawn AFTER `apply_rev_comp_projection`. R1 reads forward in the projected orientation; R2 is reverse-complement of the projected 3' end. |
| `end_loss` interaction | Paired-end windows index into the post-end-loss `rec.sequence` (the molecule the user actually observed). |
| Backwards compat | No content-hash change for baseline runs, no trace-schema bump, no cartridge change, no AIRR-field-default flip. |

Paired-end is the second observation-stage audit (after primer-
trim / end-loss). The existing scaffolding —
`reverse_complement` helper, `apply_rev_comp_projection`
composition pattern, additive AIRR-field defaults, additive
trace-schema policy — anticipated this; the audit confirms the
gap is mostly wiring + a single architectural commitment:
"paired-end is projection, not biology."

Greenlight Slice A when ready.
