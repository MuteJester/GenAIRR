# AIRR-record postcondition validator

**Status:** validator + golden tests landed. Read-only — no generation
behavior changes.

A single function that answers the question:

> Given a final `(Simulation, Outcome, RefDataConfig)` and the
> `AirrRecord` produced from it, is every reported field internally
> consistent and biologically derivable?

The validator independently re-derives each field from upstream
sources (events, trace, sim, refdata) and reports any divergence as
a structured `RecordValidationIssue`. It complements the existing
audit-doc-driven golden tests by giving a *single source of truth*
for "is this record valid" — instead of dozens of scattered
assertions per audit slice, every record produced by the engine can
be checked with one call.

---

## 1. API

```rust
// engine_rs/src/airr_record/validate.rs
pub fn validate_airr_record(
    record: &AirrRecord,
    outcome: &Outcome,
    refdata: &RefDataConfig,
) -> Vec<RecordValidationIssue>;
```

Empty vector = record passes every check. A non-empty list names the
specific divergences with structured payloads.

**Python wrapper** on `Outcome`:

```python
outcome = compiled.simulator.run(seed=42)
issues = outcome.validate_record(refdata, sequence_id="r1")
# issues: List[str] — each entry is the Debug-formatted issue variant.
assert not issues, issues
```

The validator builds the AIRR record internally and returns issues as
strings. A future iteration can expose the structured enum to Python;
the string form is sufficient for CI assertions today.

---

## 2. Check catalogue

Five categories, each tied to one of the existing audit docs:

### C1: Structural record invariants

Internal consistency of the record's own fields (no outcome
required beyond pool comparison):

- `sequence_length` matches `len(record.sequence)`.
- `record.sequence` matches the simulation pool (case-insensitive,
  since sequencing errors lowercase).
- V/D/J `*_sequence_start/end` and `*_germline_start/end` pairs are
  well-ordered (`start ≥ 0`, `end ≥ start`).
- CIGAR strings parse cleanly (M/I/D/S/N/P/X/= ops only).
- CIGAR query span (M + I op lengths) matches the segment's
  reported sequence span.

### C2: Counter provenance

Each counter is re-derived from its canonical source:

| Counter                 | Source                                                                   |
|-------------------------|--------------------------------------------------------------------------|
| `n_mutations`           | `sim.mutation_count` (set by S5F / Uniform at seal).                     |
| `n_pcr_errors`          | Trace at `ChoiceAddress::CorruptPcrCount`.                               |
| `n_quality_errors`      | Trace at `ChoiceAddress::CorruptQualityCount`.                           |
| `n_indels`              | Count of `IndelInserted` + `IndelDeleted` events from the `corrupt.indel` pass (per indel-audit §6.1). |
| `n_v_indels` / `n_d_indels` / `n_j_indels` | Same events, filtered by `event.segment` (per indel-audit §6.2). |
| `end_loss_5_length` / `end_loss_3_length` | Trace at `ChoiceAddress::CorruptEndLoss(Five/Three)`. |

### C3: Junction truth

Re-derives the junction window using the **builder's
germline-position scan** (`anchor_pool_position`), not the contract-
side `crate::junction::compute_junction`. The two diverge under
indels and end-loss — the validator must match the projection's
choice so it tests the projection, not an alternative.

Checks:
- `junction` content matches the recomputed pool slice.
- `junction_length` matches `len(junction)`.
- `vj_in_frame` matches `length % 3 == 0`.
- `stop_codon` matches `junction_has_stop(junction)` (only when in-
  frame; vacuously False otherwise).
- `productive` matches the full triad: in-frame ∧ no junction stop
  ∧ V anchor amino acid preserved ∧ J anchor amino acid preserved.
- On `productive=False`, the validator's enum payload records WHICH
  predicate fired (`OutOfFrame` / `JunctionStopCodon` /
  `VAnchorAaChanged` / `JAnchorAaChanged`).

### C4: Allele-call oracle

**Independent reimplementation** of the walker's max-match-count
selection, including the walker's NP-region extension scoring. For
each segment (V/D/J):

1. Walk the segment's structural region in `sim.sequence.regions`.
   For each pool byte with a non-None `germline_pos`, count which
   alleles in the refdata's pool match at that reference position.
   Match semantics: `A/C/G/T` canonical (case-insensitive); `N`
   wildcard matches any canonical reference base; everything else
   is no-match (per allele-call audit §1.3).
2. **NP-extension scoring** (mirrors
   `walker/extensions.rs::walk_left_extension` /
   `walk_right_extension`): walk into the adjacent NP region(s),
   capped by the assigned allele's `trim_5` / `trim_3`. Extension
   bytes are scored only when they **strictly narrow** the current
   max-score tie set (see `scoring::extension_narrows_tie_set`) —
   matching the walker's conservative "discriminate but never
   widen" policy.
3. Build the tie-set at max score (sorted by allele id).
4. Apply the projection's ordering convention: truth allele first
   when in the tie-set, otherwise ascending allele id.
5. Compare to the record's `v_call` / `d_call` / `j_call` CSV.

Both tie-set membership (order-insensitive equality) and the first-
in-CSV (the projected allele) are checked. Mismatches report the
reported and expected lists side-by-side.

**Shared scoring kernel.** Match semantics and per-byte scoring
live in [`live_call/scoring.rs`](../engine_rs/src/live_call/scoring.rs)
— `classify_base`, `matches_observed`,
`score_alleles_with_extensions`, `extension_narrows_tie_set`. Both
the walker (via `reference_index.rs::observed_base_kind` and
`walker/extensions.rs`) and the validator's C4 oracle route through
the kernel; pinned by `live_call::tests::scoring`.

**Skipped for rev-comp records.** Rev-comp post-projection flips the
sequence after the live call was computed; re-walking the
post-flip pool would compare against the wrong reference bytes.
The allele-call audit covers rev-comp with dedicated tests.

**Validator-driven runtime fix.** The validator's IGK J residual
that first surfaced when the extension-aware oracle landed was
traced to a strict-inequality bug in
`ir/builder.rs::segment_region_overlaps_dirty`: an `IndelDeleted`
at pool position `at` shrinks the region's end to `at` (when the
deletion is at the boundary), but the overlap check used
`w.start < region_end` — false at `at == region_end` — so the
walker's freshly-rebuilt sealed state was discarded and the stale
pre-deletion live call leaked through. Fixed by making the upper
bound inclusive (`w.start <= region_end`). End-loss 3' delete at
the J boundary now correctly triggers the refresh.

Without this fix the validator would have stayed quietly wrong;
with it, IGK / IGL / IGH all validate at 100% across 500 seeds.
The validator caught a real engine bug the existing audit suite
missed — that's exactly what it's for.

### C5: Region / live-call structural invariants

Two silent invariants the engine has relied on without explicit
checks until now. The validator surfaces violations:

#### §5A: At most one structural region per V/D/J segment

- Live-call (`engine_rs/src/live_call/call.rs::latest_region_for_segment`)
  picks the LATEST region for a segment via `.rev().find()`.
- AIRR projection
  (`engine_rs/src/airr_record/projection.rs`) picks the FIRST via
  `.find()`.
- The two agree by construction because no API path adds a second
  region for the same segment to `sim.sequence.regions`. The
  available methods are `with_region_added` (append once) and
  `with_region_replaced` (replace at index). Indel handling adjusts
  bounds but doesn't duplicate.

The validator pins this with `MultipleRegionsForSegment`. If a
future refactor breaks the invariant (e.g. introduces a "region
history" feature without aligning the two selectors), the validator
flags it before users see silent divergence.

#### §5B: At most one PlacementHypothesis per live call

- `SegmentLiveCall.hypotheses: Vec<PlacementHypothesis>` is a vector,
  but every production code path produces `vec![single]` or
  `Vec::new()`. (Tests can construct multi-hypothesis, but no
  pipeline does.)
- AIRR projection silently uses `hypotheses.first()`. A multi-
  hypothesis call would lose hypotheses ≥ 1 without warning.

The validator pins this with `MultipleHypothesesInLiveCall`. If
multi-hypothesis becomes real (e.g. an aligner-style scoring
extension), the validator flags every record carrying lost
information so projection can be updated deliberately.

---

## 3. Empirical sweep (160 records, 0 issues)

Across 8 representative configurations × 20 seeds:

| Config                          | Issues |
|---------------------------------|--------|
| VJ recombine only               | 0      |
| VDJ recombine only              | 0      |
| VJ + SHM                        | 0      |
| VJ + PCR + sequencing errors    | 0      |
| VJ + indels                     | 0      |
| VJ + N-corruption               | 0      |
| VJ + end-loss                   | 0      |
| VJ productive full stack        | 0      |

The validator surfaced one real divergence during development: the
initial implementation used `crate::junction::compute_junction`
(offset arithmetic) for C3, which diverges from the projection's
`anchor_pool_position` (germline_pos scan) under indels and
end-loss. Fixing the validator to match the projection's approach
was the right call — the projection's behaviour is correct
biology; the contract-side function is correct for a different
purpose (compile-time precondition checks).

This is the kind of bug the validator is designed to catch: a
subtle drift between two co-evolved code paths. It surfaced
immediately on the empirical sweep, not in a hand-curated
adversarial fixture.

---

## 4. What the validator does NOT cover

Deliberate scope limits:

- **Generation behavior**: read-only. The validator doesn't change
  what records the engine produces; it just checks them.
- **Distribution invariants**: those are statistical bounds
  ([distribution_invariant_audit.md](distribution_invariant_audit.md))
  per-record validity, not per-batch distribution.
- **Performance budgets**: separate
  ([performance_baseline.md](performance_baseline.md)).
- **Trace replay reproducibility**: each audit's replay tests cover
  that; the validator runs on a single outcome.
- **MCP-server schema checks**: the existing Python `_mcp_validators.py`
  validates the record's schema independently of the engine
  (required fields, CIGAR op alphabet, locus consistency). The two
  validators coexist: schema validator answers "is this AIRR
  compliant?", postcondition validator answers "is this record
  internally consistent with its origin?".

---

## 5. Drift items

### 5.1 Scoring model centralised — **RESOLVED**

The walker's match semantics (canonical / wildcard / lowercase),
per-byte scoring, and NP-extension scoring now live in one shared
kernel at [`engine_rs/src/live_call/scoring.rs`](../engine_rs/src/live_call/scoring.rs).
Both the walker (via `reference_index.rs::observed_base_kind` and
`walker/extensions.rs`) and the validator's C4 oracle route their
match decisions through the kernel — if the scoring model changes,
both consumers follow.

The oracle uses the per-byte-scan path (slower but independent);
the walker uses the inverted-index path (faster). Same rules, two
implementations of the rules — that independence is the validation
method.

### 5.2 IGK J residual under full mutation/indel pipeline — **RESOLVED**

The residual was traced to a boundary-condition bug in
`ir/builder.rs::segment_region_overlaps_dirty`. The function
decides whether to commit a freshly-sealed walker state by
checking if any dirty window overlaps the segment's current
region. The check used `w.start < region_end`, which failed at
exactly `w.start == region_end` — the case where an `IndelDeleted`
at the boundary shrunk the region's end to the dirty site's
position. Result: the walker's correctly-rebuilt state was
discarded and the stale pre-deletion live call stayed committed.

Fix: change the upper bound to inclusive (`w.start <= region_end`).
The lower bound stays strict because a deletion just BEFORE the
region shifts the region without changing its content.

Pinned by `test_end_loss_three_prime_does_not_strand_stale_live_call`
in `tests/test_release_validation.py`. All four release-tier
configurations (IGH VDJ, IGK / IGL VJ, IGH non-productive) now
validate at 100% across 500 seeds.

### 5.2 Validator runs over one record at a time, not the whole batch

For batch-level invariants (e.g. "every productive record's
junction_aa is non-empty"), a downstream user wraps
`validate_record` in a loop. A `validate_records(outcomes,
refdata)` convenience function returning a per-record issue list +
a batch-level summary could land later.

### 5.3 Issue payloads expose Rust Debug-format strings to Python

The Python binding returns `List[str]` where each string is the
Debug-formatted enum variant. A typed Python class (or even a
dict per issue) would be friendlier for downstream tooling
(dashboards, alerting). Deferred until a consumer asks.

---

## 6. Test coverage in this slice

[`tests/test_airr_record_validator.py`](../tests/test_airr_record_validator.py)
(20 tests):

- **Clean records pass**: parametrised over 9 configurations
  (recombine only, productive, SHM, PCR, sequencing errors,
  N-corruption, indels, end-loss, full stack); each runs 20 seeds.
- **VDJ clean**: same shape with V/D/J refdata.
- **Rev-comp skips C3/C4 by design** and validates clean.
- **§5A invariant pinned**: single region per V/D/J across 40
  records × 2 configs.
- **§5B invariant pinned**: single hypothesis in live call across
  20 indel-heavy records.
- **Counter provenance** parametrised over 5 mechanisms (mutate,
  pcr, sequencing, indels, end-loss); each runs 15 seeds and
  asserts no counter-class issues.
- **Productive triad** under high SHM (rate=0.5, 40 seeds) — the
  walker's anchor-preservation, frame, and stop-codon predicates
  all agree with the validator's re-derivation.
- **Allele oracle** under SHM rate=0.3 (40 seeds) — the
  independent rescoring agrees with the projection's `v_call`
  CSV.

All 20 pass; 0 issues across 200+ records sampled in this slice.
