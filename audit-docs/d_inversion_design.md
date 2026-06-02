# D Inversion — Design + Pre-Implementation Audit

**Status: shipped + live-call cleanup landed.** The original
Slice A → E arc shipped the IR primitives, assembly emission,
pass + trace + replay, the public DSL, and the `d_inverted`
AIRR field. A follow-up "live-call / allele-call cleanup" slice
then removed the remaining release-tier exception
(`AlleleCallTieSetMismatch{segment: D}` filtered out under
inversion) by teaching the walker and the validator's oracle to
score inverted D bytes against the original allele coordinates
via a new orientation-aware scoring primitive
(`matches_observed_with_orientation`) and per-orientation
monotonic-direction checks. Per-byte semantics:

- D germline coordinates stay in original allele orientation
  (`d_germline_start/end` unchanged).
- `d_inverted: bool` continues to record the provenance.
- Live-call scoring compares observed inverted D bytes against
  the **reverse-complement interpretation** of the candidate
  allele via the orientation-aware primitive.
- v1 of the cleanup disables NP-region extensions under RC
  orientation (extension geometry's flipped-trim-cap semantics
  needs its own audit slice); structural-region scoring alone
  is what the validator oracle compares against, so the walker
  and the oracle stay in lockstep.

After the cleanup, the IGH + `invert_d()` full-stack release tests
pass `validate_records(refdata)` without issue stripping —
matching the closure standard of receptor revision and paired-end.

---

A read-only architecture review of every engine surface that
D-segment inversion will touch, written **before** any
implementation. Companion to
[`tests/test_d_inversion_contract.py`](../tests/test_d_inversion_contract.py)
which pins today's absence and will flip to behavior-pinning when
the slice lands.

D inversion is V(D)J recombination's RSS-symmetric inversion event:
during V→D→J joining, the RSS heptamers around D can pair head-to-head
instead of head-to-tail, the D segment flips, and the assembled
junction contains the **reverse complement** of D's reference
sequence. Prevalence is low (~1–5% per literature) but biologically
real. Existing detection tools (IgBLAST, MiXCR) emit an inverted-D
call when they see this signature; GenAIRR cannot generate such
sequences today, so any pipeline trying to *test* an inverted-D
detector against GenAIRR-generated truth has nothing to point at.

---

## Existing scaffolding already in the engine

Three places where prior contributors anticipated this feature
and left a slot for it:

| Surface | Today | Slot |
| --- | --- | --- |
| [`NucFlags::INVERTED`](../engine_rs/src/ir/nucleotide.rs) | Defined as `1 << 3`; never set or read. | One-line per-base "this came from inverted D" marker. |
| [`ChoiceValue::Bool`](../engine_rs/src/trace.rs) | Docstring lists "D inversion: yes/no, receptor revision: yes/no" as canonical use cases. | Trace value type ready for the per-simulation orientation choice. |
| [`AnchorPreserved::new(Segment::D)`](../engine_rs/src/contract/anchor_preserved.rs) | Compiles, has a stable name (`"anchor_preserved.d"`). NOT in the default `productive()` bundle (which is V + J only). | Contract surface that may need an inversion-aware override later, but does not gate v1. |

`InvertDPass`, the `SegmentOrientation` enum, and the `d_inverted`
AIRR field do **not** exist today. Adding them is exactly the
implementation work this audit precedes.

---

## 1. Where in the pipeline?

D inversion is a **recombination-stage** decision: it shapes
which junction bytes downstream assembly emits. It is not a
corruption pass (corruption mutates an assembled record;
inversion is upstream of assembly).

**Recommended DSL surface:**

```python
Experiment.on("human_igh").recombine().invert_d(prob=0.05)
```

Explicit `.invert_d(prob=...)` step (a new method on
`Experiment`), not a kwarg on `recombine(...)`. Reasons:

- Matches the established corruption-DSL shape (`pcr_amplify`,
  `mutate`, `end_loss_5prime`, …): one step → one biology.
  `recombine()` already takes ~8 kwargs; piling on another
  raises the cognitive cost of reading a chain.
- Future SHM / orientation extensions (receptor revision,
  inverted V/J catalogues) can mirror the same pattern without
  re-shaping `recombine()`.
- Lets the user position the step explicitly: a future contract
  bundle that wants D forward-only can declare its requirement
  against the named step.

**Pipeline position:** `InvertDPass` runs **after**
`SampleAllelePass(Segment::D)` and **before**
`AssembleSegmentPass(Segment::D)`. The trim passes for D run
between sample and assemble, and inversion is independent of
trim choice (trim coordinates remain in original allele space —
see §5), so the relative ordering with `TrimPass(Segment::D, …)`
is irrelevant to correctness. Topological constraints in
`Schedule`:

- `requires` D allele assignment (so `SampleAllelePass(D)` runs first).
- `effects` D orientation assignment (so
  `AssembleSegmentPass(D)` consumes it).

**Chain-type guard.** `InvertDPass.precondition` rejects VJ
catalogues at compile time (no D pool → nothing to invert).
DSL-level guard surfaces this as a `ValueError` before run.

---

## 2. Trace addresses

One new entry in [`ChoiceAddress`](../engine_rs/src/address.rs):

```rust
SampleAlleleDInverted,  // displays as "sample_allele.d.inverted"
```

Why `sample_allele.d.inverted` and not `recombine.d.inverted`:

- There is no `recombine.*` namespace today. Every existing
  trace address sits under a concrete pass name (`sample_allele`,
  `trim`, `assemble`, `generate_np`, `mutate`, `corrupt`).
  Inventing a `recombine.*` namespace just for one address
  forces every downstream tool to learn a new prefix.
- The orientation decision is attached to D's sampled allele —
  it's a *property of the D assignment*, conceptually
  parallel to `trim.d_5` / `trim.d_3` (also attributes of the
  D assignment, also under their own pass namespace).
- `sample_allele.d.inverted` extends the existing
  `sample_allele.d` namespace by one level; the pattern
  matches `np.np1.bases[N]` (an addressable choice that
  belongs to a parent pass).

**Value:** `ChoiceValue::Bool(true)` when D is inverted,
`Bool(false)` when forward. Both ends recorded explicitly so the
trace reader never has to infer absence-means-false.

**Address-schema-version bump.** Adding a new `ChoiceAddress`
variant whose `Display` string does not collide with existing
addresses is **not** a schema bump per the docs at the top of
[`address.rs`](../engine_rs/src/address.rs). Old traces have no
`sample_allele.d.inverted` records; new traces are parseable by
future engines. Keep `ADDRESS_SCHEMA_VERSION = 1`.

---

## 3. SimulationEvent

Add one new variant — **the general form**, per user preference:

```rust
SimulationEvent::OrientationChanged {
    segment: Segment,
    old: SegmentOrientation,
    new: SegmentOrientation,
}
```

Marked `**Reserved for future emission.**` initially, like
[`AssignmentChanged`](../engine_rs/src/ir/sim_event.rs#L153)
and the other `with_*` events. `InvertDPass` becomes the first
caller; future inverted-V / inverted-J / receptor-revision
slices reuse the same variant by passing a different `segment`.

Sinks that only care about D pattern-match on
`segment == Segment::D`. The walker observer, event-log
observer, and dirty-window tracker all have a precedent here:
they pattern-match on the event's discriminant and ignore the
rest.

Emission point: inside `AlleleAssignments::with_orientation(…)`
(parallel to the existing `with_trim`), so any builder that
calls into the persistent setter gets the event for free.

---

## 4. IR representation

Per user spec: do **not** mutate the reference allele byte slice
or branch the refdata pool. Orientation is per-simulation state,
not per-cartridge state.

Add a new enum in [`engine_rs/src/assignment.rs`](../engine_rs/src/assignment.rs):

```rust
#[derive(Copy, Clone, Eq, PartialEq, Debug, Default)]
pub enum SegmentOrientation {
    #[default]
    Forward,
    ReverseComplement,
}
```

Extend [`AlleleInstance`](../engine_rs/src/assignment.rs#L52):

```rust
pub struct AlleleInstance {
    pub allele_id: AlleleId,
    pub trim_5: u16,
    pub trim_3: u16,
    pub orientation: SegmentOrientation,
}
```

Default constructor sets `orientation: SegmentOrientation::Forward`.
A new persistent setter `with_orientation(o)` matches the existing
`with_trim_5` / `with_trim_3` shape.

**Why on `AlleleInstance` and not on the assignment dict.** D
orientation is a property of the *sampled assignment* — same
allele can be sampled forward in one simulation and reversed in
the next. Carrying it on the instance keeps the (allele_id,
trim_5, trim_3, orientation) tuple as the complete per-simulation
state of that segment, and the persistent setter pattern (`with_*`
returning a new instance) extends without ceremony.

**Trim semantics stay in original allele coordinates.** `trim_5`
removes N bases from position 0 of the reference allele; `trim_3`
removes N bases from the high end. After inversion the
*emitted* byte order is reversed and complemented, but the
trim numbers still describe the original allele coordinates.
This keeps:
- `cfg.reference_models.trims[("D", "5")]` empirical
  distributions valid as-is — they sampled in original coords.
- the validator's anchor-window math identical between forward
  and reversed instances (the contract reads the same
  `(trim_5, trim_3)` bytes in either case).
- DSL ergonomics (`.trim(d_5=2)`) consistent.

---

## 5. Assembly behavior

[`AssembleSegmentPass`](../engine_rs/src/passes/assemble_segment/execution.rs#L20)
reads the assignment and emits the post-trim slice. The
inversion-aware change is localized to the emit loop
(lines 120–123 in execution.rs):

```rust
// Forward (today's only path):
for (i, &base) in slice.iter().enumerate() {
    let allele_pos = slice_start + i as u32;
    builder.push_nucleotide(Nucleotide::germline(base, allele_pos as u16, seg));
}

// With orientation support:
match inst.orientation {
    SegmentOrientation::Forward => {
        for (i, &base) in slice.iter().enumerate() {
            let allele_pos = slice_start + i as u32;
            builder.push_nucleotide(Nucleotide::germline(base, allele_pos as u16, seg));
        }
    }
    SegmentOrientation::ReverseComplement => {
        for (i, &base) in slice.iter().rev().enumerate() {
            // First emitted byte comes from `slice_end - 1`, complemented.
            // Next from `slice_end - 2`, etc. germline_pos is the
            // *source* position in the original allele coordinate.
            let allele_pos = (slice_end - 1) as u32 - i as u32;
            let complement = wc_complement(base);
            builder.push_nucleotide(
                Nucleotide::germline(complement, allele_pos as u16, seg)
                    .with_flags(NucFlags::INVERTED),
            );
        }
    }
}
```

Key invariants under reversal:

| Property | Forward | Reverse-complement |
| --- | --- | --- |
| `region.len()` | `slice_len` | `slice_len` (unchanged) |
| `region.start..end` in pool | sequential | sequential (regions are pool-coordinate) |
| `germline_pos` of i-th emitted base | `slice_start + i` | `slice_end - 1 - i` |
| Base byte | `allele.seq[allele_pos]` | `wc_complement(allele.seq[allele_pos])` |
| `NucFlags::INVERTED` | unset | **set** |
| Walker / live-call lookup | uses original allele bytes via `ref_index` | **must complement at lookup** — see §6 |

**Helper:** add a `wc_complement(u8) -> u8` function in
`ir/nucleotide.rs` next to the existing flag definitions
(`A↔T`, `C↔G`, `N→N`; panic on unrecognised byte to match the
strict-alphabet discipline already in `ReferenceAlphabet`).

**Frame phase math is unchanged.** Reversal doesn't change the
emitted region's length, so the cumulative-length-mod-3
calculation that computes `frame_phase` is byte-for-byte
identical to the forward path.

---

## 6. AIRR projection

The hardest design question, per user. Three sub-questions:

### 6.1 `d_germline_alignment` orientation

**Decision:** keep `d_germline_*` coordinates in **original
allele orientation**. The post-inversion *emitted* bytes carry
`germline_pos` pointing into the original allele coords (see
§5), so the walker can derive `d_germline_start` as the
*smallest* `germline_pos` it saw on D, and `d_germline_end` as
the *largest + 1*. Under inversion that yields the same numeric
range it would have produced for the forward case — `[trim_5,
allele_len - trim_3)` — which is what an external aligner like
IgBLAST emits when it reports `d_call="IGHD…"` with
`d_orientation="inverted"`.

### 6.2 `d_cigar` representation

**Decision:** `d_cigar` continues to describe alignment in
original allele coordinates, traversed in pool-byte-emission
order. For an inverted D the walker sees emitted bytes
`[c̄_{end-1}, c̄_{end-2}, …, c̄_{start}]` and produces a CIGAR
identical in shape to the forward case (length = D pool span,
ops = M/D/I) because each byte still matches its original-coord
position. The cigar reader doesn't need to know about
inversion at all.

This makes life easy for downstream tools that interpret
CIGAR independently of orientation — they read the same
3M-D-5M cigar whether the D was forward or reversed.

### 6.3 New non-AIRR field: `d_inverted: bool`

Add a GenAIRR-specific field at the end of the AIRR record:

```rust
pub d_inverted: bool,  // false by default
```

Lives next to `is_contaminant` in the "GenAIRR additions"
section of [`AirrRecord`](../engine_rs/src/airr_record/record.rs).
Projection rule: `d_inverted = sim.assignments.get(Segment::D)
.map(|i| matches!(i.orientation, ReverseComplement))
.unwrap_or(false)`.

VJ chains (no D) report `d_inverted = false`; same shape as
existing `d_*` field defaults.

**Python surface:** the field flows through
[`engine_rs/src/python/outcome.rs`](../engine_rs/src/python/outcome.rs)
as `dict.set_item("d_inverted", rec.d_inverted)?`, and
appears as a column in the `SimulationResult` DataFrame /
CSV / FASTA exports (parallel to `is_contaminant`).

---

## 7. Contracts

The default `productive()` bundle is **untouched** by D
inversion. Reasoning per contract:

- `ProductiveJunctionFrame` — junction length is a pool-byte
  count; reversal preserves length; frame mod-3 unchanged.
- `NoStopCodonInJunction` — scans pool bytes for in-frame
  TAA/TAG/TGA. The pool already contains the reverse-complemented
  D bytes by the time the contract checks. Same scan, different
  bytes — exactly the behavior the slice is supposed to model.
- `AnchorPreserved::V`, `AnchorPreserved::J` — read V / J
  bytes only, never touch D.
- `AnchorPreserved::D` — exists but **not** in `productive()`.
  Real D alleles have no canonical conserved codon
  (`anchor: None` for every bundled D), so the contract is a
  no-op for the bundled cartridges regardless of inversion.

**No new contracts in v1.** A future slice may want
`DOrientationAllowed(set: HashSet<SegmentOrientation>)` to let
users force forward-only or compare inverted-D detector
outputs — but that's a follow-up, not a blocker.

**Productive interaction worth pinning.** A D inversion can
turn a sequence productive (by shifting which junction codons
become stops) or unproductive (same mechanism). Under
`productive_only`, the engine's existing constraint-narrow
machinery will reject inverted-D samples that produce
stop-codons-in-junction at the InvertDPass commit time, just as
it rejects bad trim/NP combinations today. The constraint
machinery is generic; no D-inversion-specific hook needed.

---

## 8. Replay

[`replay_from_trace_file`](../engine_rs/src/replay.rs) consumes
the trace's choices verbatim and re-runs the pipeline. The new
`sample_allele.d.inverted` record must replay deterministically.
Three validation duties:

### 8.1 Pass-level validation

`InvertDPass::execute_one_replay(cursor)` calls
`cursor.expect_bool(ChoiceAddress::SampleAlleleDInverted)` (the
same idiom [`CorruptRevCompPass`](../engine_rs/src/passes/corrupt/rev_comp.rs#L60)
uses for `corrupt.rev_comp.applied`). Misshaped values
(`Int`, `Bases`, etc.) surface as `PassError::trace_mismatch` —
the canonical replay-time mismatch error.

### 8.2 Chain-type validation

The plan-compile step must reject `InvertDPass` on a VJ
catalogue. There is no D pool, no D assignment, nothing to
invert. The check matches today's chain-type guards on
`AssembleSegmentPass(D)` / `SampleAllelePass(D)`.

### 8.3 Schedule ordering

`InvertDPass` declares `requires: AlleleAssignment(D)` and
`effects: OrientationAssignment(D)`.
`AssembleSegmentPass(D)` declares `requires: OrientationAssignment(D)`.
The dep graph then auto-derives the edge
`InvertDPass → AssembleSegmentPass(D)` and rejects
plans that omit `InvertDPass` after `SampleAllelePass(D)`
*if* the user attached `invert_d` to the chain (the new
effect is opt-in; defaulting to absent means today's plans
continue to compile and execute identically — see §9).

---

## 9. Backwards compatibility

Every existing simulation, every existing AIRR record, every
existing trace file pre-dates this slice and must keep working.

- Default `AlleleInstance::orientation = Forward` means
  `AlleleInstance::new(id)` is byte-equivalent to today.
- `AirrRecord::d_inverted = false` by default; existing
  comparisons / fixtures continue to pass.
- Existing trace files have no `sample_allele.d.inverted`
  record; replay treats absence as "no `InvertDPass` ran" and
  the assembled D is forward.
- Bundled cartridges' D pools do not gain any new metadata;
  catalogues are agnostic to whether inversion is enabled.
- The content hash format adds **no** new bytes — orientation
  is per-instance state, not per-cartridge state.
- `refdata_signature` (the structural cartridge fingerprint)
  is unaffected; v2 trace replay continues without bump.

---

## 10. Test surface — what the implementation slice must pin

A non-exhaustive list, mirrored as TODO markers in the companion
[`tests/test_d_inversion_contract.py`](../tests/test_d_inversion_contract.py):

1. **Engine round-trip.** Inverted D nucleotides carry
   `NucFlags::INVERTED`; `germline_pos` of i-th emitted base
   equals `slice_end - 1 - i`.
2. **Walker lookup.** Live-call walker on inverted D resolves
   the same `d_call` set the forward path would have, when the
   D allele is the truth.
3. **AIRR field.** `rec.d_inverted == True` for every inverted
   sample; `False` for every forward and every VJ sample.
4. **Junction content.** With a fixed seed and D=`AAATTTGGG`
   anchor=None, forward D contributes `AAATTTGGG` to the
   junction; reverse-complemented D contributes `CCCAAATTT`.
5. **Productive interaction.** Under `productive_only`,
   inverted-D samples that produce in-frame stops in junction
   are rejected at sample time.
6. **Replay determinism.** A trace file emitted with
   `invert_d(prob=0.5)` replays bit-for-bit, including the
   per-record orientation choice.
7. **Content hash stability.** Adding `invert_d(prob=0.0)` to
   a plan does not change AIRR output records (orientation
   stays Forward, so behavior is identical).
8. **Trace address vocabulary.** A round-trip through
   `ChoiceAddress::parse` / `Display` preserves
   `sample_allele.d.inverted`.

---

## 11. Implementation order

Five-slice sequence shipped. Each slice was independently revertible
and gated on a green test sweep before the next started.

1. **Slice A — IR primitives** ✅ shipped.
   [`SegmentOrientation`](../engine_rs/src/assignment.rs) enum,
   [`AlleleInstance.orientation`](../engine_rs/src/assignment.rs)
   field + `with_orientation` setter, `complement_base` helper in
   [`engine_rs/src/ir/nucleotide.rs`](../engine_rs/src/ir/nucleotide.rs),
   `NucFlags::INVERTED` round-trip tests. Default values kept every
   existing test passing.
2. **Slice B — Assembly emission** ✅ shipped.
   [`AssembleSegmentPass`](../engine_rs/src/passes/assemble_segment/execution.rs)
   honors `orientation` for D only (V/J `ReverseComplement` silently
   no-ops, by design). Rust unit tests at the pass level pin
   reverse-complemented bytes + descending `germline_pos` + the
   `NucFlags::INVERTED` flag.
3. **Slice C — Trace + replay** ✅ shipped.
   [`ChoiceAddress::SampleAlleleDInverted`](../engine_rs/src/address.rs),
   [`InvertDPass`](../engine_rs/src/passes/invert_d.rs) with fresh +
   replay paths, `SimulationEvent::OrientationChanged` generalised
   over all segments, `SimulationBuilder::update_allele_orientation`.
   No DSL yet.
4. **Slice D — DSL** ✅ shipped.
   [`Experiment.invert_d(prob=…)`](../src/GenAIRR/experiment.py)
   with VJ-chain + duplicate-call guards.
   `_InvertDStep` dataclass + inline lowering inside
   [`_lower_recombine`](../src/GenAIRR/_compile.py) (via
   `_extract_invert_d_prob` pre-pass) — keeps the canonical
   V-NP1-D-NP2-J pool layout intact without explicit schedule
   edges. End-to-end pytest in
   [`tests/test_invert_d_dsl.py`](../tests/test_invert_d_dsl.py).
   First two pin-the-absence audit tests flipped.
5. **Slice E — AIRR projection** ✅ shipped.
   `AirrRecord.d_inverted: bool` field sourced from the final IR
   (not the trace), `SimulationResult` column wired in
   [`result.py`](../src/GenAIRR/result.py),
   `RecordValidationIssue::DInvertedMismatch` validator hook.
   Lockstep updates to this §11, §10 invariants below, and the
   last pin-the-absence audit test flipped.

The five slices added 0 new schedule edges, 0 new content-hash
bytes, and 0 trace-schema bumps — every byte produced before
Slice A still round-trips through replay today.

---

## 12. Out of scope

Documented here so a future contributor doesn't accidentally
expand the slice:

- **Inverted V or J.** Biologically attested but rare; the
  `SegmentOrientation` enum is general enough that adding
  `invert_v` / `invert_j` later is mechanical, but v1 is
  D-only.
- **Receptor revision** (re-rearrangement of an already-assembled
  receptor). Different mechanism; different orchestration; not
  in this slice.
- **Orientation in the cartridge catalogue.** Some catalogues
  ship pre-inverted D alleles as separate entries. v1 leaves
  catalogue choice unopinionated — users who want both
  orientations in their pool can keep adding both today.
- **Detection-tool comparison harness.** Comparing GenAIRR
  inverted-D output against IgBLAST / MiXCR detector calls
  is a downstream evaluation slice, not part of generation.

---

## Summary table

| Question | Decision |
| --- | --- |
| Pipeline position | New `InvertDPass` after `SampleAllelePass(D)`, before `AssembleSegmentPass(D)`. |
| DSL surface | `Experiment.invert_d(prob=…)` as an explicit step. |
| Trace address | `sample_allele.d.inverted`, `ChoiceValue::Bool`. |
| SimulationEvent | New general-form `OrientationChanged { segment, old, new }`. |
| IR carrier | `AlleleInstance.orientation: SegmentOrientation` (default Forward). |
| Assembly | Reverse iteration + per-base WC complement under `ReverseComplement`; `germline_pos` keeps original-allele coords; emit with `NucFlags::INVERTED`. |
| AIRR | Keep D coords in original orientation; add GenAIRR-specific `d_inverted: bool`. |
| Productive bundle | Untouched. New `DOrientationAllowed` contract deferred. |
| Replay | `expect_bool` at the new address; chain-type guard; auto-derived schedule edge to assembly. |
| Backwards compat | Defaults preserve every existing simulation, trace, and AIRR record. |

D inversion has been silently expected by the engine since the
`NucFlags::INVERTED` flag and `ChoiceValue::Bool` docstring
landed. Slices A → E converted that expectation into a working
feature without breaking a single existing test, trace, or AIRR
record. Future slices (V/J inversion, receptor revision,
catalogue-side inverted alleles) can reuse the
`SegmentOrientation` enum, the `OrientationChanged` event, and
the `d_inverted` projection pattern without re-shaping any of
the architecture this five-slice arc put in place.
