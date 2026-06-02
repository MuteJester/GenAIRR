# Receptor Revision — Design + Pre-Implementation Audit

**Status: shipped.** Slices A–E + the release-consolidation slice
are merged; the user-facing surface is
[`Experiment.receptor_revision(prob=…)`](../src/GenAIRR/experiment.py)
and the validator gates the two new AIRR provenance fields
(`receptor_revision_applied`, `original_v_call`). Per-slice
deliverables:

| Slice | What landed | Verified by |
|---|---|---|
| A — IR primitives | `SimulationEvent::SegmentReplaced`, `SimulationBuilder::replace_segment`, `Simulation::with_segment_replaced` setter, `PoolRange::after_segment_replacement` kernel | `engine_rs/src/passes/receptor_revision.rs` + `airr_record/tests/projection.rs` (Slice A unit tests) |
| B — Refresh plan | `LiveCallRefreshStep::SegmentReplaced(Segment)` + per-segment dedup + AllStructural-equivalent hook execution | `engine_rs/src/live_call/refresh_plan.rs::tests`, `engine_rs/src/live_call/refresh_hook.rs::tests` |
| C — Pass + trace + replay | `ReceptorRevisionPass`, three `receptor_revision.*` trace addresses, replay path validating retained length | `engine_rs/src/passes/receptor_revision.rs::tests` (17 unit tests) |
| D — DSL | `Experiment.receptor_revision(prob=…)` (VDJ-only, at-most-once), inline-lowering after `assemble.j` | [`tests/test_receptor_revision_dsl.py`](../tests/test_receptor_revision_dsl.py) (16 e2e tests) |
| E — AIRR provenance | `receptor_revision_applied: bool` + `original_v_call: String` fields, `ReceptorRevisionAppliedMismatch` + `OriginalVCallMismatch` validator issues, structured PyO3 dict / source strings | `airr_record/tests/projection.rs` (8 unit tests) + DSL Slice E tests |
| Consolidation | README row, validation-matrix row, release-tier IGH stack test + full-stack replay round-trip, Bernoulli draw distribution test (prob=0.25, N=4000, ±5σ) | [`tests/test_release_validation.py`](../tests/test_release_validation.py), [`tests/test_distribution_invariants.py`](../tests/test_distribution_invariants.py) |

The 19-test [`tests/test_receptor_revision_contract.py`](../tests/test_receptor_revision_contract.py)
finished the arc as a fully-green scaffold suite — every pin-the-
absence test flipped to pin-the-presence in lockstep with its
implementation slice.

---

A read-only architecture review of every engine surface that
V-replacement receptor revision will touch, written **before**
any implementation. Companion to
[`tests/test_receptor_revision_contract.py`](../tests/test_receptor_revision_contract.py)
which pins today's absence + the existing scaffolding and will
flip to behavior-pinning as implementation slices land.

Receptor revision is V(D)J biology's "second chance" mechanism:
after a B-cell completes initial V(D)J recombination, the
existing V allele can be replaced by a more-upstream V via a
cryptic RSS heptamer in the existing V's 3' end. The downstream
D, NP1, NP2, J, and junction context is preserved (often
modified by additional N-nucleotide addition); only the V part
changes. Functionally it's a way for the immune system to edit
auto-reactive receptors before commitment. Prevalence varies by
tissue (single-digit % in peripheral B cells; up to ~25% in
germinal-center clones at some loci).

Unlike D inversion (Slice A → E), receptor revision changes an
**already-assembled biological truth**: the V allele the
simulation committed to is silently swapped post-recombination.
That makes it the right audit to run before any future "edit the
final truth" mechanism (gene conversion, somatic deletion,
isotype switching).

---

## Existing scaffolding already in the engine

Five places where prior contributors anticipated this feature
and left a slot for it — more than D inversion started with.

| Surface | Today | Slot |
| --- | --- | --- |
| [`SimulationEvent::RegionReplaced`](../engine_rs/src/ir/sim_event.rs#L178) | Variant defined with explicit `old: Region, new: Region` payload. Marked "Reserved for future emission" in the parent module doc. | The structural-edit channel a revision pass would broadcast over. |
| [`SimulationBuilder::replace_region`](../engine_rs/src/ir/builder.rs#L334) | Persistent setter wired up: captures `old`, emits `RegionReplaced`, delegates to `Simulation::with_region_replaced_for_segment`. Today only fires from internal tests. | The IR commit path a revision pass would call. |
| [`Simulation::with_region_replaced_for_segment`](../engine_rs/src/ir/simulation.rs#L178) | Matches on `(segment, start)`; returns `self.clone()` unchanged when no match. | The persistent-IR primitive `replace_region` delegates to. |
| [`SimulationBuilder::assign_allele`](../engine_rs/src/ir/builder.rs#L279) | Emits `AssignmentChanged { old, new }` with `old: Option<AlleleInstance>`. `Some(prior)` fires on replacement. Today only `SampleAllelePass` calls it (always with `old=None`). | The assignment-replacement path that propagates the new V instance into the IR. |
| [`RecordValidationIssue::MultipleRegionsForSegment`](../engine_rs/src/airr_record/validate.rs#L138) | Postcondition fires when `sim.sequence.regions` contains > 1 region for the same V / D / J segment. | The validator-side enforcement of the one-region-per-segment invariant we want to preserve. |

What's **NOT** present today:
- `Experiment.receptor_revision` DSL method.
- `ReceptorRevisionPass` (or any pass that calls `replace_region` in production).
- Any `receptor_revision.*` trace address.
- AIRR record fields like `receptor_revision_applied` / `original_v_call`.
- Bulk pool-bytes replacement primitive (`SegmentReplaced` event or a `splice_region_bytes` builder method).
- `LiveCallRefreshPlan` reaction to `RegionReplaced` — see §6.

The scaffolding gap is much smaller than D inversion's. The new
work is mostly **wiring**, with one genuine design decision
(length-changing pool surgery) and two contract surfaces that need
to be made revision-aware (refresh plan + admit-mask observer).

---

## 1. Biological scope

**Recommended v1**: V-replacement only, same locus, post-assembly /
pre-mutation.

- The committed V allele (and its assembled bytes) is replaced
  with a different V from the **same** V pool.
- D, J, NP1, NP2, and the existing junction window are preserved
  by construction (no resampling).
- The pipeline position sits between the end of `recombine()`
  and the start of `mutate()` / corruption — so SHM operates on
  the *revised* sequence.
- VDJ chains only initially. Light-chain V replacement is
  biologically attested but the v1 surface stays on heavy-chain
  IGH where the literature data is densest.

**Explicitly out of scope for v1** (documented for future
contributors so the slice doesn't accidentally expand):

- **Full secondary rearrangement** with a new J / new junction
  / fresh NP. Different mechanism, different orchestration —
  matches the D-inversion §12 "out of scope" discipline.
- **Light-chain receptor editing** (κ→λ chain switch). Requires
  a paired-chain cartridge surface that doesn't exist yet.
- **Iterative revision** — v1 is at most one revision per
  simulation. Wire the "at most one" guard in the DSL just like
  `Experiment.invert_d`.

Heavy-chain V-only is the smallest slice that captures the
distinctive feature: an already-assembled biological truth
mutates *post-recombination*.

---

## 2. Pipeline position + DSL surface

**Recommended DSL:**

```python
Experiment.on("human_igh").recombine().receptor_revision(prob=0.05)
```

Named `receptor_revision`, not `revise_v` — the biological term is
standard in the literature and a future light-chain extension
wouldn't force a rename. The docstring must document the v1
heavy-chain-V-only scope so users with light-chain catalogues see
a `ValueError` immediately, not a silent forward.

Pipeline position is after the `recombine()` step and the
implicit `trim()`, but before `mutate()` / corruption / clonal
fork. The lowering inserts a `ReceptorRevisionPass(prob)` between
the last `assemble.*` pass of recombine and the first `mutate.*`
pass. Insertion-order tiebreak + the auto-derived
`AlleleAssignment(V) → RegionReplaced(V)` ordering (once the
new `PassCompileEffect::ReplaceRegion(Segment)` lands — see §5)
keeps the position right by construction.

DSL-level guards mirror `invert_d`'s exactly:

- VJ chains → `ValueError("receptor_revision is only valid for
  VDJ chains")`.
- At most one revision step per pipeline → `ValueError("receptor
  revision already configured")` on duplicate call.
- `prob` finite, `[0.0, 1.0]`, NaN-rejected before range check.

---

## 3. Trace addresses (minimal replayable choice set)

Replay must reproduce every decision the fresh pass made. For
receptor revision the **minimal** set is:

| Address | Type | Required for replay? |
| --- | --- | --- |
| `receptor_revision.applied` | `Bool(bool)` | Yes — gates the entire pass. |
| `receptor_revision.v_allele` | `AlleleId(u32)` | Yes when applied=true. The replacement V identity. |
| `receptor_revision.v_trim_3` | `Int(i64)` | Yes when applied=true. The new V's 3' trim length. |

What we **don't** record:
- `v_trim_5` — biology says V5' is locked by the first
  recombination event; revision doesn't move it. Pin
  `trim_5 = 0` (or whatever the original V carried) by
  construction in the pass; no trace record needed.
- New NP1 length / NP1 bases. Junction is preserved by
  construction; NP1 stays.
- New D, J, NP2, junction bytes — none of these change.

That's a **three-record** payload: one Bool + one AlleleId + one
Int per simulation. Mirrors `InvertDPass`'s ratio (one Bool +
zero per-allele records) plus the structural pair the new V
allele demands. Cheap to emit, cheap to replay, no schema bump
needed (additive vocabulary, just like D inversion's
`sample_allele.d.inverted`).

**Address namespacing:** `receptor_revision.*` rather than nesting
under `sample_allele.v.*`. The decision isn't a sub-choice of the
original V sampling; it's a separate biological event that
*replaces* the V. The trace-reader vocabulary should reflect that.

`ADDRESS_SCHEMA_VERSION` stays at 1 — additive variants don't
bump per the policy at the top of [`address.rs`](../engine_rs/src/address.rs).

---

## 4. IR representation

Per user spec: **replace the existing V region in place**,
emitting `AssignmentChanged` + `RegionReplaced` events. The
one-region-per-segment invariant is preserved by construction
(the validator's `MultipleRegionsForSegment` check then keeps the
post-condition green for free — see [`validate.rs:797`](../engine_rs/src/airr_record/validate.rs#L797)).

**The hard sub-question: what happens when the new V has a
different length than the old V?**

- New V is **shorter**: pool bytes get removed from the V slice,
  every downstream `Region.start/end` shifts left. NP1, D, NP2,
  J coordinates all need updating.
- New V is **longer**: pool bytes get inserted into the V slice,
  every downstream region shifts right.
- New V is **same length**: pool bytes get replaced 1:1, no
  downstream shifting.

Three approaches, ordered by complexity:

### 4.1 Recommendation: bulk-replacement primitive

Add a new IR primitive **`Simulation::with_segment_replaced`**
that takes:
- the target `Segment` to replace,
- the new region's start coordinate (= old region's start),
- the new bytes to install in the pool,
- the new region's allele coordinates (`germline_pos` per byte).

It does the pool surgery (delete old bytes, insert new bytes)
and the downstream-region shifting in one atomic persistent
revision. Emits a new event:

```rust
SimulationEvent::SegmentReplaced {
    segment: Segment,
    old_region: Region,
    new_region: Region,
    bytes_delta: i32,  // new_len - old_len
}
```

**Why a new event variant.** `RegionReplaced` payload already has
`(old, new): Region`, but the existing event was designed for
*coordinate-only* replacement (same-length, e.g., reframe a
region after a downstream indel). A length-changing replacement
crosses into the structural-event category (alongside
`IndelInserted` / `IndelDeleted`) and the refresh path needs to
treat it as such. Reusing `RegionReplaced` for both would conflate
"metadata edit" with "pool surgery."

### 4.2 Rejected: emit a sequence of IndelInserted/IndelDeleted

Could implement V-replacement as N `IndelDeleted` events
followed by M `IndelInserted` events. Pros: reuses existing
infrastructure (refresh plan already reacts to indels). Cons:
the trace bloats by ~50× (one event per byte), the events
attribute each byte to the V segment when they're really part
of a single revision act, and the AIRR validator's `n_indels`
counter would over-count.

### 4.3 Rejected: forbid length-changing replacement

Require revision to pick a new V with the same trimmed length
as the old V. Pros: no pool surgery. Cons: catalogues have V
alleles of varying lengths (≈285 bp range across human IGHV);
restricting to same-length matches collapses the candidate set
to single digits and the revision rate effectively becomes
zero on real catalogues.

---

## 5. Events

**Expected event stream from a single revision firing:**

1. `AssignmentChanged { segment: V, old: Some(prior_v), new: revised_v }`
   — emitted by `builder.assign_allele(Segment::V, revised_instance)`.
2. `SegmentReplaced { segment: V, old_region, new_region, bytes_delta }`
   — emitted by the new `builder.replace_segment` (Slice C of
   the implementation roadmap).
3. Zero or more `BasePushed` (if the new V is longer than the
   old V) or `BaseDeleted` (if shorter) events for the changed
   bytes themselves — emitted by the bulk-replacement primitive
   inside its single atomic revision.
4. `TrimChanged` if the new `v_trim_3` differs from the old.

Three **new** infrastructure pieces this requires:

| Piece | Where | Effort |
| --- | --- | --- |
| `SimulationEvent::SegmentReplaced` variant | [sim_event.rs](../engine_rs/src/ir/sim_event.rs) | Mirror of `RegionReplaced`; 1-line enum addition + observer arms. |
| `SimulationBuilder::replace_segment` | [builder.rs](../engine_rs/src/ir/builder.rs) | New method; ~30 lines, follows the `assign_allele` / `replace_region` template. |
| `LiveCallRefreshStep::SegmentReplaced(Segment)` | [refresh_plan.rs](../engine_rs/src/live_call/refresh_plan.rs) | New step kind; the plan must invalidate the V live call AND any downstream segment whose coordinates shifted. Largest single piece of new work in this audit. |

`PassCompileEffect::ReplaceRegion(Segment)` belongs in
[`metadata.rs`](../engine_rs/src/pass/metadata.rs) — gives the
schedule analyser something to reason about so a later pass that
consumes `AlleleAssignment(V)` (e.g., a future `mutate.*`
variant that pre-builds an admit mask from V's anchor) knows it
must run after revision.

---

## 6. Live-call refresh

Today's [`LiveCallRefreshPlan`](../engine_rs/src/live_call/refresh_plan.rs#L137)
silently ignores `RegionReplaced`. That has been correct so far
because no production pass emits the event; replacing V in v1
implementation breaks this assumption.

**Required updates:**

1. **The V live call must re-walk** from scratch against the
   revised V region's bytes. This is a `LiveCallRefreshStep::SingleSegment(V)`
   trigger — equivalent to what `RegionAdded { segment: V }` does
   today during initial assembly.
2. **Downstream segments whose coordinates shifted** (NP1, D,
   NP2, J under §4.1's length-changing replacement) **don't** need
   their walker scores recomputed — the bytes themselves didn't
   change, only the `Region.start/end` coordinates did. The
   walker stores its hypotheses in terms of pool handles, which
   were already shifted by the persistent IR revision. A
   refresh-plan no-op for D / J under `SegmentReplaced(V)` is
   correct.
3. **Edge case**: V's anchor codon position changes (revision
   moves the V Cys). The junction window's V-anchor boundary
   must be re-derived. The existing junction-truth re-derivation
   in [`junction.rs`](../engine_rs/src/junction.rs#L114) reads
   `sim.assignments.get(V).anchor` against the V region — both
   updated by the revision act — so the projection is correct
   the next time AIRR is built. No new infrastructure here.

The refresh-plan update is the **largest concrete piece of new
engine work** the audit identifies. Worth a dedicated slice
(Slice B of the implementation roadmap, after the event variant
+ builder method).

---

## 7. AIRR provenance

Two GenAIRR-specific fields, sitting next to `d_inverted` /
`is_contaminant` in the "GenAIRR additions" section of
[`AirrRecord`](../engine_rs/src/airr_record/record.rs):

```rust
pub receptor_revision_applied: bool,  // default false
pub original_v_call: String,          // empty when no revision fired
```

Population rules:

- `receptor_revision_applied`: sourced from the **trace**
  (`receptor_revision.applied` Bool), not the final IR — the IR
  carries only the post-revision V assignment, so the only
  surviving provenance after the revision commits is the trace
  record itself. Symmetric with `is_contaminant`'s sourcing.
- `original_v_call`: when applied=true, the V-allele name from
  the *initial* `sample_allele.v` trace record (mapped through
  `refdata.v_pool` for the name lookup). When applied=false:
  empty string.

`v_call` itself continues to carry the **post-revision** allele
identity — i.e., the V the simulation actually committed to and
that downstream code (mutate, projection, validation) acted on.
This matches the existing `d_inverted` precedent (`d_call` is
the truth-allele name; `d_inverted` is the provenance bit).

**Rejected alternative**: a single `original_v_call` field with
empty-string semantics for "no revision." Less ergonomic — every
consumer would have to check `original_v_call != ""` instead of
`receptor_revision_applied is True`. The Bool + name pair mirrors
the established `d_inverted: bool` / `d_call: str` shape.

---

## 8. Contracts

The default `productive()` bundle interacts with revision in
three ways:

### 8.1 Anchor preservation

`AnchorPreserved::V` already checks that V's Cys codon survives
trims and SHM. Under revision, the V assignment changes. The
contract's existing `admits_sample_allele_candidate` path
already reads the candidate's anchor against `refdata`, so
replaying the contract over the new assignment is correct as
long as the revision pass routes the candidate through the
contract's `admits_typed` check before committing.

### 8.2 Productive junction frame

`ProductiveJunctionFrame` reads the junction length post-hoc.
Under revision the V tail length can change (if the new V has
different `trim_3`), which shifts the junction window. The
contract sees the post-revision state when junction is computed
at projection time — correct by construction.

### 8.3 No-stop-codon-in-junction

`NoStopCodonInJunction` scans pool bytes inside the junction
window. The window changes; the scan happens at projection time
against the revised pool. No new work.

**The hard contract**: under `productive_only()`, the revision
pass needs to constrain its candidate V set to alleles whose
anchor codon would *also* survive the existing trim choices.
Without this, the constraint-aware narrow that recombination
respected gets silently violated by revision. This is
**`ReceptorRevisionPass::admit_mask_observer` work** — the v1
implementation must register a contract-narrowing observer on
the V-candidate distribution, matching how `SampleAllelePass`
does today.

### 8.4 Admit-mask observer reaction to RegionReplaced

[`ProductiveAdmitMaskObserver`](../engine_rs/src/contract/admit_mask_observer.rs#L268)
currently silently ignores `RegionReplaced`. Once revision fires
in production, the observer must re-derive the admit mask for
any pass *downstream* of revision (e.g., a SHM pass that pre-
builds a mask from V's bytes). For v1, since SHM doesn't
pre-build per-V masks (it reads bytes at sample time), this is
a TODO but not a blocker for the first revision slice.

---

## 9. Replay

`InvertDPass`'s replay path is the template. For
`ReceptorRevisionPass`:

```rust
let applied = cursor.expect_bool(ChoiceAddress::ReceptorRevisionApplied)?;
if applied {
    let v_id = cursor.expect_allele_id(ChoiceAddress::ReceptorRevisionVAllele)?;
    let v_trim_3 = cursor.expect_int(ChoiceAddress::ReceptorRevisionVTrim3)?;
    // commit through builder.assign_allele + builder.replace_segment
    // (Slice C of the implementation roadmap)
}
ctx.trace.record_choice(...);
```

Replay validation duties:

1. **Schema-version compatibility**: pre-revision trace files
   have no `receptor_revision.*` records; the replay path treats
   absence as `applied = false` and silently no-ops. Bundled
   golden traces continue to round-trip without re-generation.
2. **Chain-type guard**: the pass declares
   `PassRequirement::AlleleAssignment(V)` and panics-with-
   structured-error (`PassError::InvalidPlanState`) on VJ chains.
   Matches `InvertDPass`'s D-only narrowing.
3. **Schedule edges**: `ReceptorRevisionPass.requirements =
   [AlleleAssignment(V)]` + `effects = [ReplaceRegion(V)]`.
   The latter is the new `PassCompileEffect` variant (§5); it
   lets downstream passes declare a requirement *on* the revised
   V if they need to. Today no pass does, but the slot stays
   open.

---

## 10. Validator

`validate_airr_record` already exercises:

- `MultipleRegionsForSegment` — pinned today; revision preserves
  the invariant by construction (replaces in place).
- `AlleleCallTieSetMismatch` — V tie-set oracle. Re-walks the
  V live call from scratch against the post-revision IR. Same
  expected behaviour as for D inversion's filter: the v1 walker
  may or may not perfectly score the revised V depending on
  trim choices; the audit flags this as a Slice F "known-
  incomplete" carry forward, not a blocker for Slice A.

**New validator check** (Slice E of the implementation roadmap):

```rust
RecordValidationIssue::ReceptorRevisionAppliedMismatch {
    reported: bool,
    expected: bool,  // = trace.find("receptor_revision.applied").value
},
RecordValidationIssue::OriginalVCallMismatch {
    reported: String,
    expected: String,  // = refdata.v_pool[trace.receptor_revision.v_allele].name
},
```

Same shape as `DInvertedMismatch`. Catches a downstream consumer
that manually flips `receptor_revision_applied` post-build, or a
fork of the AIRR builder that forgets to populate the field.

---

## 11. Backwards compatibility

Same five-slice arc D inversion demonstrated, same compat
guarantees:

- `AlleleAssignments::with_assigned` already allows replacement
  (the `slots.set(...)` path overwrites without panic) — no IR
  schema change, no content-hash change.
- `AirrRecord::receptor_revision_applied = false` and
  `original_v_call = ""` by default — every existing test
  fixture continues to compare equal.
- Existing trace files have no `receptor_revision.*` records;
  replay treats absence as `applied = false`. v1 / v2 trace
  files continue to round-trip without bumps.
- The new `SegmentReplaced` event variant is additive; existing
  observers' fall-through arms ignore it until the refresh-plan
  / admit-mask observer updates land in their dedicated slices.
- Bundled cartridges: no changes — receptor revision is a
  per-simulation behaviour, not a per-cartridge property.

---

## 12. Implementation order

Recommended slice sequence; each independently revertible and
gated on a green test sweep before the next:

1. **Slice A — IR primitives.** `SegmentReplaced` event variant
   + `SimulationBuilder::replace_segment` + the
   `Simulation::with_segment_replaced` persistent setter +
   `PassCompileEffect::ReplaceRegion(Segment)`. Default values
   keep every existing test passing.
2. **Slice B — Refresh plan + admit mask.** Wire `SegmentReplaced(V)`
   into `LiveCallRefreshPlan` (re-walk V; no-op for D/J).
   Update `ProductiveAdmitMaskObserver` to react when a future
   downstream pass needs admit-mask invalidation. Pure
   plumbing; no DSL yet.
3. **Slice C — `ReceptorRevisionPass` + trace addresses.** Add
   `ChoiceAddress::ReceptorRevisionApplied` / `ReceptorRevisionVAllele`
   / `ReceptorRevisionVTrim3`. Pass commits through
   `builder.assign_allele` + `builder.replace_segment`. Replay
   safe. No DSL yet.
4. **Slice D — DSL.** `Experiment.receptor_revision(prob=…)`
   with VJ + duplicate guards. End-to-end Python tests.
   Flip the relevant pin-the-absence audit tests.
5. **Slice E — AIRR provenance.** `AirrRecord.receptor_revision_applied`
   + `original_v_call` + `ReceptorRevisionAppliedMismatch`
   validator hook + `SimulationResult` columns. Lockstep docs.

Each slice ≤ one PR-sized chunk; the D inversion arc took five
slices to land in the same shape.

---

## 13. Test surface — what the implementation slices must pin

A non-exhaustive list, mirrored as TODO markers in the companion
[`tests/test_receptor_revision_contract.py`](../tests/test_receptor_revision_contract.py):

1. **Engine round-trip.** `SegmentReplaced(V)` event payload
   carries the correct `(old_region, new_region, bytes_delta)`
   tuple; the persistent IR revision shifts D / NP / J handles
   atomically.
2. **One-region-per-V invariant.** `validate_records()` continues
   to pass the `MultipleRegionsForSegment` check before and
   after revision.
3. **AIRR field.** `receptor_revision_applied = True` ↔ trace
   carries `receptor_revision.applied = Bool(true)`. `v_call`
   equals the *new* V; `original_v_call` equals the V from the
   first `sample_allele.v` trace record.
4. **Trace replay.** A trace file emitted with
   `receptor_revision(prob=0.5)` replays bit-for-bit, including
   the per-record V-allele choice and trim choice.
5. **Schedule ordering.** Compile report places
   `receptor_revision` between the last `assemble.*` of recombine
   and the first `mutate.*` of SHM. Trace order pins this.
6. **Content-hash stability.** Adding `receptor_revision(prob=0.0)`
   to a plan does not change AIRR output — the new V is identical
   to the original V because the revision never fires.
7. **Productive interaction.** Under `productive_only()`,
   revised V's that would break junction frame / introduce a
   stop codon are rejected at sample time by the contract's
   admit-mask narrowing.
8. **Distribution invariant.** With `prob=0.25` over N=4000 seeds,
   the empirical revision rate matches 0.25 within 5σ (same
   pattern as `invert_d` D9 from the consolidation slice).

---

## 14. Out of scope

Documented here so a future contributor doesn't accidentally
expand the slice:

- **Light-chain receptor editing.** Requires a paired-chain
  cartridge surface; out of scope until that lands.
- **Iterative revision** (multiple successive V replacements in
  one simulation). v1 caps at one revision per pipeline.
- **Revision in the cartridge catalogue.** Some catalogues
  annotate cryptic RSS sites; v1 leaves catalogue choice
  unopinionated. Future: a `revisable=true` flag on `Allele`
  could narrow the candidate V set.
- **N-nucleotide addition during revision.** Real receptor
  editing can add new junctional N nucleotides at the cryptic
  RSS site. v1 preserves the existing NP1 untouched; future
  work could model the added N's as a separate
  `revision_added_n_nucleotides: i64` field.
- **Per-tissue revision rate priors.** Different B-cell
  compartments have different revision rates; v1 ships a single
  `prob` knob. A future presets slice could expose
  `germinal_center_default_prob = 0.20` etc.
- **Detection-tool comparison harness.** Comparing GenAIRR
  revised-V output against IgBLAST / MiXCR detector calls is a
  downstream evaluation slice, not part of generation.

---

## Summary table

| Question | Decision |
| --- | --- |
| Scope | Heavy-chain V replacement only, same locus, post-recombine / pre-SHM. |
| DSL | `Experiment.receptor_revision(prob=…)` (NOT `revise_v`). VDJ-only, at most once per pipeline. |
| Pipeline position | Between the last `assemble.*` of recombine and the first `mutate.*` of SHM. |
| Trace | `receptor_revision.applied` (Bool), `.v_allele` (AlleleId), `.v_trim_3` (Int). Schema-version unchanged. |
| IR carrier | Replace existing V region via new `Simulation::with_segment_replaced`; assignment goes through `builder.assign_allele`. |
| Event | New `SegmentReplaced { segment, old_region, new_region, bytes_delta }` variant; `AssignmentChanged` also fires. `RegionReplaced` stays for metadata-only edits. |
| Live-call refresh | New `LiveCallRefreshStep::SegmentReplaced(V)` triggers V re-walk; D/J no-op (their bytes didn't change). |
| AIRR | Add `receptor_revision_applied: bool` and `original_v_call: String`. `v_call` carries post-revision identity. |
| Contracts | `AnchorPreserved::V` already correct. New: revision pass registers an admit-mask observer to constrain candidate V's under `productive_only`. |
| Replay | Three-record consume via `expect_bool` / `expect_allele_id` / `expect_int`; chain-type guard; absent records → `applied = false` (back-compat). |
| Validator | New `ReceptorRevisionAppliedMismatch` + `OriginalVCallMismatch` issues; existing `MultipleRegionsForSegment` continues to defend the invariant. |
| Backwards compat | No content-hash change, no trace-schema bump, no cartridge change, no AIRR-field-default flip. |

Receptor revision crosses one architectural threshold D inversion
didn't: it **mutates an already-assembled biological truth**. The
existing scaffolding (RegionReplaced event, replace_region
builder, AssignmentChanged with `old: Option<...>`,
MultipleRegionsForSegment validator) anticipated this — the
audit confirms the gap is mostly wiring, plus one genuinely-new
piece (`SegmentReplaced` event + refresh-plan reaction).
Greenlight Slice A when ready.
