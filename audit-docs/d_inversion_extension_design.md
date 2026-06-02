# D Inversion — NP-Extension Scoring Audit

**Status: shipped.** The implementation slice landed with the
shape this audit specified: `walk_left_extension` and
`walk_right_extension` are orientation-aware (no RC early-return),
the trim caps swap by allele direction under RC per §3, the
`extension_narrows_tie_set` gate is reused unchanged, and the
validator's `score_alleles_with_extensions` mirrors the same
logic so walker and oracle stay in lockstep. The v1 Rust pin
`rc_extensions_are_disabled_in_v1` flipped to
`rc_extensions_narrow_truth_allele`; a companion regression test
`rc_extension_trim_cap_swap_pool_left_consumes_trim_3` pins the
audit's §3 trim-cap mapping at the kernel level. The Python
contract suite (`tests/test_d_inversion_extension_contract.py`)
ran the lockstep flip from `pin_v1_boundary_*` to
`pin_narrowing_*` in a single commit with the implementation.

After this slice, inverted D's adjacent NP evidence narrows the
allele-call tie set when it strictly distinguishes alleles. The
release-tier IGH stack with `invert_d(prob=1.0)` continues to run
`report.ok` clean with zero exception filtering. Cache parity
under `invert_d(prob=1.0)` stays green; walker and validator
agree on every record.

---

A read-only architecture review of the one remaining
documented-boundary in the orientation-aware live-call walker:
NP-region extensions under `ReverseComplement` D orientation.
Companion to
[`tests/test_d_inversion_extension_contract.py`](../tests/test_d_inversion_extension_contract.py)
which pins today's behaviour and the synthetic-fixture cases the
implementation slice satisfied.

The D-inversion arc closed at the
[live-call / allele-call cleanup slice](d_inversion_design.md):
walker and validator now route per-byte scoring through
`matches_observed_with_orientation` so the **structural** D region
contributes correct evidence under either orientation. The cleanup
intentionally left the **NP-region extension walks disabled under
`ReverseComplement`** as a temporary boundary. This audit defined
the per-orientation semantics that closed that boundary; the
implementation slice landed exactly the four-file change the
audit proposed.

---

## Existing scaffolding the audit relies on

Surfaces today's engine provides, which the proposed
implementation will compose against. Each is pinned by a
`pin_scaffold_*` test in the companion contract file.

| Existing surface | Where | What it gives the extension audit |
|---|---|---|
| `matches_observed_with_orientation` primitive | `engine_rs/src/live_call/scoring.rs:108` | Single-byte comparison rule that handles both `Forward` and `ReverseComplement` — extensions inherit it for free. |
| `extension_narrows_tie_set` gate | `engine_rs/src/live_call/scoring.rs:229` | Orientation-agnostic: it inspects the score vector and the per-byte match mask. The proposed RC extensions reuse it unchanged. |
| `alleles_compatible_at(alleles, ref_pos, observed, orientation)` | `engine_rs/src/live_call/scoring.rs:201` | The slow per-allele scan the validator's oracle uses. Already orientation-aware as of the previous slice; extensions just need to point it at the right `ref_pos`. |
| `SegmentRefIndex::compatible_alleles_at_oriented` | `engine_rs/src/live_call/reference_index.rs:159` | The walker's inverted-index lookup. Same orientation contract — extensions only need to compute the correct `ref_pos`. |
| `ExtensionWalkState.orientation` field | `engine_rs/src/live_call/walker/extensions.rs:55` | Already plumbed through. The proposed implementation reads it directly. |
| `find_left_extension` / `find_right_extension` | `engine_rs/src/live_call/scoring.rs:524` (and matching pool helpers) | Locate the NP regions adjacent to the structural D region. The proposed implementation reuses these without modification. |
| AssembleSegmentPass(D) under invert | `engine_rs/src/passes/assemble_segment/execution.rs:138` | Documents the emission convention: pool-left D byte has the highest `germline_pos`, pool-right D byte has the lowest. The audit's pool-vs-allele direction reasoning hinges on this convention. |

---

## 1. Current state — what v1 ships

After the live-call / allele-call cleanup slice, the engine
processes inverted D as follows:

1. **Assembly.** `AssembleSegmentPass(D)` under `inst.orientation
   == ReverseComplement` emits the retained slice in **reverse
   allele order**, complementing each byte:
   ```text
   pool[seq_start + i].base = complement(allele[slice_end - 1 - i])
   pool[seq_start + i].germline_pos = slice_end - 1 - i
   ```
   So the **pool-left** D byte has `germline_pos = slice_end - 1`
   (the highest allele coord retained), and the **pool-right** D
   byte has `germline_pos = slice_start` (the lowest).
2. **Structural scoring.** Both `WalkerObserverState::on_base_pushed`
   and `walker::call_from_region` route per-byte lookups through
   `SegmentRefIndex::compatible_alleles_at_oriented`, which
   pre-complements the observed byte under `ReverseComplement`.
   The result equals `allele[germline_pos]`, so the inverted-index
   bitset at the original allele coordinate is the right place to
   look. **Score evidence on the structural region is correct.**
3. **Monotonic check.** Under `ReverseComplement` the monotonic
   direction is flipped: ref_pos is expected to decrease
   monotonically (gap-down allowed for deletions). `ref_start`
   tracks `min(seen)`, `ref_end` tracks `max(seen) + 1`. Same
   half-open `[ref_start, ref_end)` semantics as Forward.
4. **Extensions.** `walk_left_extension` and `walk_right_extension`
   in `walker/extensions.rs` early-return at entry when
   `state.orientation == ReverseComplement`. The validator's
   `score_alleles_with_extensions` mirrors the skip with a
   `skip_extensions = matches!(orientation, ReverseComplement)`
   guard around both extension blocks. So **under inverted D,
   walker and oracle both produce structural-only scores**, which
   means they agree on the tie set — releasable behaviour, but
   one bit of evidence is being silently dropped.

The structural-only path is **what makes the validator and walker
agree** without exception filtering. That agreement is the
prerequisite for closing this audit's proposed implementation
without re-introducing release-tier flake.

---

## 2. Architectural question 1: pool direction ↔ allele direction under RC

**Q.** For inverted D, does pool-left extension correspond to
allele-right extension, and pool-right extension correspond to
allele-left extension?

**A. Yes.** The emission convention pins this:

```
Forward D:
  pool position [s, s+1, s+2, ..., s+L-1]
  germline_pos  [trim_5, trim_5+1, ..., slice_end-1]
                ↑ low ref ──────────────→ high ref

ReverseComplement D:
  pool position [s, s+1, s+2, ..., s+L-1]
  germline_pos  [slice_end-1, slice_end-2, ..., trim_5]
                ↑ high ref ──────────────→ low ref
```

So extending **left in the pool** (one position before the
structural D's `seq_start`) corresponds to:
- **Forward**: the immediate predecessor in allele order, i.e.
  `ref_pos = trim_5 - 1`. This is what `walk_left_extension`
  already implements.
- **ReverseComplement**: the immediate successor in allele order,
  i.e. `ref_pos = slice_end - 1 + 1 = slice_end`. The next allele
  position past the highest one the structural walk scored.

And extending **right in the pool** corresponds to:
- **Forward**: `ref_pos = slice_end` (the next allele position
  past the highest structural one).
- **ReverseComplement**: `ref_pos = trim_5 - 1` (the immediate
  predecessor in allele order — the next-lower coordinate the
  walker hasn't seen yet).

The natural geometry: **the pool-direction extension walks are
swapped in allele-coordinate space under `ReverseComplement`.**

The proposed implementation **does not** swap pool iteration
direction. It swaps the **candidate ref_pos formula** inside
`walk_left_extension` / `walk_right_extension`:

| Walk | Forward candidate ref_pos | ReverseComplement candidate ref_pos |
|---|---|---|
| `walk_left_extension` (extending into NP region whose end touches the structural region's `seq_start`) | `*state.ref_start - 1` (and decrement on accepted step) | `*state.ref_end` (and increment on accepted step) |
| `walk_right_extension` (extending into NP region whose start touches the structural region's `seq_end`) | `*state.ref_end` (and increment) | `*state.ref_start - 1` (and decrement) |

Both walks still iterate pool bytes in pool order; the only
orientation-dependent piece is which ref boundary moves and in
which direction.

### 2.1 Why not "iterate pool in the opposite direction under RC"?

The pool's NP-region geometry isn't symmetric — `find_left_extension`
returns the NP region whose `.end.index() == state.seq_start`, which
is the one to the *pool-left* of the structural region. The
"opposite-pool-direction" alternative would require
reinterpreting which NP region each helper returns. That ripples
through `score_alleles_with_extensions` and any future NP-aware
consumer; the cleaner factoring keeps the pool helpers
orientation-agnostic and concentrates the orientation logic in
the ref-coord computation.

---

## 3. Architectural question 2: trim-cap swap under RC

**Q.** How do `trim_5` and `trim_3` caps swap under
reverse-complement orientation?

**A. They swap by allele direction, not by pool direction.** The
caps cap **how many allele-coord positions the extension can
walk past the structural boundary in either direction**:

- `trim_5` caps the number of allele positions BELOW the lowest
  scored allele coord (`ref_start`).
- `trim_3` caps the number of allele positions ABOVE the highest
  scored allele coord (`ref_end`).

Under **Forward** the pool walks left into NP1 to reach allele
positions below `ref_start` (so the left walk consumes the
`trim_5` budget), and pool walks right into NP2 to reach allele
positions above `ref_end` (so the right walk consumes the
`trim_3` budget). This is the existing implementation.

Under **ReverseComplement** the pool-direction-to-allele-direction
mapping flips (per §2):
- Pool walks left → allele positions ABOVE `ref_end` → consumes
  the `trim_3` budget.
- Pool walks right → allele positions BELOW `ref_start` →
  consumes the `trim_5` budget.

So the call sites in `call_from_region` and
`score_alleles_with_extensions` need to **swap which cap is passed
to which walk** under RC:

```rust
// Forward
walk_left_extension(sim, idx, np_region, trim_cap_5, &mut state);
walk_right_extension(sim, idx, np_region, trim_cap_3, &mut state);

// ReverseComplement
walk_left_extension(sim, idx, np_region, trim_cap_3, &mut state);
walk_right_extension(sim, idx, np_region, trim_cap_5, &mut state);
```

The walks themselves only need a single cap, and the audit's
recommendation is to keep the **single-cap signature** rather than
introducing a "both caps + orientation" composite. The orientation
lives at the call site; the walk function stays cap-direction-
agnostic.

---

## 4. Architectural question 3: narrowing gate under RC

**Q.** Should extension bytes still follow the "only consume if
the byte narrows the current max-score tie set" rule?

**A. Yes, unchanged.** `extension_narrows_tie_set(scores, matched)`
inspects:
1. `pre_max > 0` — the structural walk has seeded scores.
2. At least one max-tied allele matches the byte (some current
   max-allele's reference at `candidate_ref_pos` equals the
   orientation-transformed observed byte).
3. At least one max-tied allele doesn't.

None of these depend on orientation. The orientation transforms
the byte at the comparison point — once `matched` is computed,
the gate is direction-agnostic. The proposed implementation reuses
`extension_narrows_tie_set` verbatim.

Specifically: under `ReverseComplement`, "matches" means
`allele.seq[candidate_ref_pos] == complement(observed)`. The
matched mask is then the input to the unchanged narrowing gate.

---

## 5. Architectural question 4: extensions affect tie sets, not coordinates

**Q.** Can extension evidence change only the D tie set, without
changing D coordinates / CIGAR?

**A. Yes.** Today the extension walks update both `state.scores`
(tie-set evidence) AND `state.ref_start` / `state.ref_end` /
`state.seq_start` / `state.seq_end` (boundaries). The boundary
updates feed:
- `germline_alignment` columns
- `d_germline_start` / `d_germline_end`
- `d_cigar` op lengths
- `BOUNDARY_ELASTIC` / `OVERLAPS_OTHER_SEGMENT` flags

For inverted D, the boundary updates are non-trivial: extending
into NP corresponds to ALLELE-coord movement in the opposite
direction from Forward. The audit's question is whether boundary
updates under RC should be **bidirectional** or **scoring-only**.

### 5.1 Recommendation: bidirectional boundary updates

D inversion design §6 already pinned that `d_germline_start/end`
stay in original allele orientation. Under RC, an accepted
left-pool extension step:
- Increments `state.ref_end` (matches allele orientation: extending
  past the highest scored allele coord)
- Decrements `state.seq_start` (matches pool orientation: the new
  pool byte is to the left)

And an accepted right-pool extension step:
- Decrements `state.ref_start` (matches allele orientation:
  extending past the lowest scored allele coord)
- Increments `state.seq_end` (matches pool orientation)

The half-open invariants stay:
- `[seq_start, seq_end)` is the pool-coord range the scored region
  covers (contiguous regardless of orientation).
- `[ref_start, ref_end)` is the allele-coord range scored.

This is the simplest path and preserves the d_germline_start/end
contract for inverted D — the bytes the validator's
`germline_alignment` walker materialises against the (forward)
allele reference still come from the correct allele positions.

### 5.2 Alternative considered: scoring-only updates

Updating only `state.scores` (not boundaries) under RC would
"silently" tighten the tie set without touching CIGAR or
coordinates. But it would diverge the walker's behavior from the
Forward path in a non-orientation-related way: under Forward the
boundaries update when extensions narrow the tie set, and the
non-update under RC would create an asymmetry that future readers
would have to remember. The audit recommends consistency.

---

## 6. Architectural question 5: missed-evidence fixtures

**Q.** Are there fixtures where disabling RC extensions
under-calls ambiguity or over-reports a tie?

**A. Yes — and the contract file pins one explicitly.** The
canonical case: two D alleles tie at every retained structural
byte, but differ at the trimmed-off byte that an adjacent NP byte
would have matched. Under Forward, the extension narrows to the
truth allele; under RC v1, the extension is disabled and the tie
remains. The companion test
`test_pin_v1_boundary_rc_extensions_disabled_drops_narrowing_evidence`
constructs the fixture, asserts the structural-only tie set IS
the current behaviour, and labels the alternative-tie-set the
proposed implementation must produce. When the implementation
slice lands, this test flips from pin-the-boundary to
pin-the-narrowing-evidence.

A non-fixture observation: the **over-reports-a-tie** case is
strictly less likely under RC v1 than under the alternative
"silently include extensions without orientation awareness." The
extension narrowing rule is conservative — it only consumes
bytes that strictly shrink the tie set. So the v1 disable is
**evidence-faithful but loose**: it never includes spurious
narrowing evidence; it just omits real narrowing evidence. That
asymmetry is what makes v1 acceptable as a temporary boundary.

### 6.1 Why this matters

The release-tier IGH stack ran clean under v1 because:
1. NP regions are short (typically 0–6 bases per side).
2. The structural region has enough bytes to narrow the tie set
   for most allele pools.
3. The truth-allele projection fallback handles the residual
   "tied with truth" case correctly.

But the engine's documented promise is **evidence-faithful allele
ambiguity** — the call set reflects what the evidence supports,
no more, no less. The v1 disable trades evidence faithfulness for
implementation simplicity, and the audit's purpose is to close
that trade.

---

## 7. Proposed implementation

A single-PR-sized slice, scoped to the extension walks only. No
DSL changes, no new trace addresses, no API surface beyond the
existing `state.orientation`.

### 7.1 Changes to `walker/extensions.rs`

`walk_left_extension` and `walk_right_extension` both:

1. Drop the early-return at entry when
   `state.orientation == ReverseComplement`.
2. Compute `candidate_ref_pos` per orientation (see §2).
3. Compute the boundary update direction per orientation (see §5).
4. Reuse `compatible_alleles_at_oriented` for the per-byte lookup
   (unchanged).
5. Reuse `extension_narrows_tie_set` for the gate (unchanged).

The pool-direction iteration stays identical to v1. The
`crossed_into_overlap` flag stays meaningful: it still records
whether the extension reached past the NP region into the
preceding/following structural region.

### 7.2 Changes to `walker/mod.rs::call_from_region`

The two `walk_left_extension` / `walk_right_extension` call sites
need orientation-aware cap routing (see §3): the cap passed to
`walk_left_extension` becomes `trim_cap_3` under RC (vs.
`trim_cap_5` under Forward), and symmetrically for
`walk_right_extension`.

### 7.3 Changes to `scoring.rs::score_alleles_with_extensions`

Mirror the validator-side equivalent of §7.1 + §7.2:

1. Drop the `skip_extensions` guard.
2. Compute `candidate_ref_pos` per orientation in the inlined
   left and right extension blocks.
3. Update `current_ref_start` / `current_ref_end` per the §5
   boundary rule.
4. Swap which cap drives the left vs. right block (see §3).

The validator oracle and the walker remain in lockstep by sharing
the orientation rule — no `_strip_*` filter or release-tier
exception can be reintroduced.

### 7.4 Changes to tests

1. `rc_extensions_are_disabled_in_v1` (currently in
   `engine_rs/src/live_call/tests/scoring.rs`) flips from
   "pin the boundary" to "pin the narrowing-evidence
   contract under RC."
2. New `engine_rs/src/live_call/tests/scoring.rs::rc_extensions_narrow_truth_allele`
   constructs the synthetic fixture from §6 directly against the
   scoring kernel and asserts the post-implementation tie set is
   `{truth}` rather than `{truth, sibling}`.
3. `tests/test_d_inversion_extension_contract.py::test_pin_v1_boundary_rc_extensions_disabled_drops_narrowing_evidence`
   flips from pin-the-boundary to pin-the-narrowing-evidence and
   adds an inverse assertion (the Forward equivalent of the same
   fixture already narrows).

### 7.5 Out of scope for this implementation slice

- Extension behaviour under `Forward` orientation is unchanged.
- Walker-vs-oracle lockstep continues to hold (both routes through
  the shared kernel).
- No changes to assembly emission, AIRR fields, validator issues,
  or the public DSL.

---

## 8. Edge cases the implementation must handle

| Edge case | Forward behaviour | RC expected behaviour |
|---|---|---|
| Empty NP region (zero NP bytes) | Walk loop iterates zero times; no extension. | Same. |
| `trim_3 = 0` under RC + non-trivial NP1 | Right-pool walk caps at 0 → zero extension steps. | Right-pool walk (which corresponds to allele-left under RC) caps at `trim_5`. So a Forward-trim_3=0 fixture **does not** disable RC extensions; the *Forward-trim_5* and *Forward-trim_3* of the assigned `AlleleInstance` are the relevant caps under RC. |
| Extension byte fails canonical-base check | Walk halts (current behaviour). | Same — the per-byte check happens before the orientation transform. |
| `OVERLAPS_OTHER_SEGMENT` flag | Set when extension reaches past NP into the preceding/following structural V or J region. | Same condition; flag interpretation is pool-direction-only and unchanged by orientation. |
| Walker rebuilds via `from_existing_region` after an indel-edit pass | Re-walks the structural region with the same orientation; no extension changes. | Same — the rebuild path doesn't enter the extension walks. |

---

## 9. Replay determinism

The extension walks are deterministic functions of `(sim,
segment_index, np_region, trim_cap, state)`. None of those inputs
are trace-recorded; the trace records the high-level pass choices
(allele id, NP length, NP bytes). Replay reproduces the same
`sim`, which reproduces the same extension behaviour — regardless
of orientation.

The implementation slice changes only the scoring kernel +
walker; no new trace addresses are needed, and the existing
trace-file schema (version unchanged) parses old + new traces
identically.

---

## 10. Validator integration

The validator's `oracle_check_segment` calls
`score_alleles_with_extensions` and compares the resulting tie
set against `record.d_call`. As long as the walker and the kernel
share the same orientation rule (which §7 ensures), the oracle
agrees with the walker by construction. No new validator issue
variant is needed; no new structured `details.source` string;
no new column on the JSON dict.

This is the same lockstep guarantee the receptor revision and
paired-end audits established. The implementation slice neither
expands nor narrows the validator's issue surface.

---

## 11. Backwards compatibility

The implementation slice's only behaviour change is:
- Under `ReverseComplement` D orientation, the live-call tie set
  may be **narrower** than today (extensions can now narrow it).

Concretely:
- `d_call` strings may shrink (one allele instead of N), making
  the projection more specific.
- The validator continues to agree with the walker. No issue
  variants surface that didn't surface before.
- `d_inverted` semantics unchanged.
- `d_germline_start/end` semantics unchanged.
- `d_cigar` lengths may grow by up to `trim_5 + trim_3` bytes
  under RC if extensions reach past the structural boundary —
  same shape as the Forward extension behaviour.

No content-hash change is **expected** for records that don't
trigger the v1 disable in the first place (i.e., Forward
records); every existing Forward fixture continues to produce
byte-identical output.

For records that DID trigger the v1 disable (inverted D with
adjacent NP narrowing evidence), the post-implementation tie set
will be narrower. Existing tests that snapshot inverted-D
records' `d_call` would surface the diff. The audit recommends:
- Replace any inverted-D snapshot tests with structural
  invariants (validator-clean, walker ↔ oracle parity, NP byte
  identity), or
- Re-snapshot after the implementation lands and add a comment
  explaining why the value changed.

---

## 12. Implementation order

This audit's full implementation is a single slice, intentionally
small:

1. **Implementation slice — RC extensions.** Implements §7's
   four-step change to `walk_left_extension`,
   `walk_right_extension`, `call_from_region`, and
   `score_alleles_with_extensions`. Flips the contract pins. No
   doc churn except marking this audit as "shipped."

The audit deliberately rules out staged sub-slices (e.g.,
"enable but boundary-update-only" or "boundary-update but no
scoring narrowing") because the walker ↔ oracle parity
requirement makes any half-step unstable.

---

## 13. Test surface — what the implementation slice must pin

A non-exhaustive list, mirrored as TODO markers in the companion
[`tests/test_d_inversion_extension_contract.py`](../tests/test_d_inversion_extension_contract.py):

1. **Synthetic narrowing under RC.** A 2-allele D pool tied on
   the structural region but distinguished by an NP byte. v1
   leaves the tie at 2 alleles; the implementation slice narrows
   to 1.
2. **Forward control fixture stays unchanged.** The same fixture
   under Forward already narrows to 1; the implementation slice
   must not change this.
3. **Walker ↔ oracle lockstep under RC.** Both routes through
   the shared kernel produce the same tie set on every fixture.
4. **Trim-cap swap respected.** A fixture where the Forward
   `trim_3` cap exhausts before the Forward `trim_5` cap; the
   inverted-D counterpart's extension is bounded by `trim_3`
   (which corresponds to the pool-LEFT walk under RC).
5. **`d_germline_start/end` and `d_cigar` reflect the elastic
   boundary correctly under RC.** Pool-direction boundary updates
   are the same as Forward; allele-direction boundary updates are
   orientation-aware.
6. **Release-tier IGH + invert_d still validates clean.** No
   new `AlleleCallTieSetMismatch{D}` after the slice.
7. **Cache parity under `invert_d(prob=1.0)` stays clean.** The
   incremental walker's `from_existing_region` rebuild path
   handles the orientation correctly (already true today; pin it
   here so a refactor that adds RC extensions to the walker but
   not to the rebuild path surfaces).
8. **`rc_extensions_are_disabled_in_v1` flips.** The Rust unit
   test changes from "pin the boundary" to "pin the narrowing
   evidence under RC."

---

## 14. Out of scope

Documented here so a future contributor doesn't accidentally
expand the slice:

- **V / J inversion.** This audit covers D only; V and J
  inversion semantics aren't modelled and the assembly pass's
  orientation check is D-specific.
- **Extension narrowing in C (constant region) live calls.**
  Constant region live-call evidence isn't surfaced as a tie set
  today; whether it should is out of scope.
- **Non-canonical observed bases.** The orientation transform
  preserves canonical / wildcard / invalid classification; the
  audit doesn't change the per-byte semantics.
- **Inverted-D AIRR fields beyond `d_inverted`.** The audit
  doesn't propose new AIRR columns. Downstream consumers that
  want a more granular signal can read the
  `assignments.d.orientation` from the trace.

---

## Summary table

| Question | Decision |
| --- | --- |
| Pool direction ↔ allele direction under RC | Swapped. Pool-left extension corresponds to allele-right; pool-right corresponds to allele-left. |
| `trim_5` / `trim_3` cap mapping under RC | Swapped. Pool-left walk consumes `trim_3` budget; pool-right walk consumes `trim_5` budget. |
| Narrowing gate under RC | Unchanged. `extension_narrows_tie_set` is orientation-agnostic. |
| Tie set vs. coordinates | Both update under accepted extension steps. `d_germline_start/end` stay in original allele orientation; pool-direction boundary updates are the same as Forward. |
| Missed-evidence fixtures | A 2-allele D pool tied structurally but distinguished by an NP byte — Forward narrows, RC v1 doesn't. The contract test pins this case. |
| Implementation surface | Four files: `walker/extensions.rs`, `walker/mod.rs`, `scoring.rs`, and `airr_record/validate.rs` (no validator-side change beyond inheriting from the kernel update). |
| Walker ↔ oracle parity | Maintained by sharing the single scoring kernel. No new exception filter. |
| Backwards compatibility | Forward records byte-identical. Inverted-D records may have narrower tie sets; existing `d_call` snapshots may need re-snapshotting. |
| Trace schema | Unchanged. No new addresses. |
| Validator issues | Unchanged. No new variants. |

D inversion's allele-call promise — evidence-faithful ambiguity
that reflects every byte the engine has scored — is one
implementation slice away from being fully integrated under
either orientation. Greenlight when ready.
