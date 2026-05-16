# Eliminate `Nucleotide::NO_GERMLINE_POS` sentinel via `GermlinePos` newtype

**Date:** 2026-05-16
**Status:** Design approved, ready for implementation plan
**Scope:** Rust engine — `engine_rs/src/`

## Problem

`Nucleotide.germline_pos: u16` uses `u16::MAX` as a magic sentinel for "no
germline provenance" (NP, P-nuc, contaminant, indel-inserted bases). Every
consumer that reads `germline_pos` must remember to check
`== Nucleotide::NO_GERMLINE_POS` before using the value. Nothing in the type
system enforces this — a forgotten check produces a silent wrong-answer bug
rather than a compile error.

This is a real risk, not a hypothetical one. Bug 4 in the project memory
(S5F mutating indel-inserted nucleotides and writing a null byte into the
annotation string) was caused by exactly this shape — a code path that
forgot the sentinel check.

## Goal

Replace the `u16` field with a newtype that makes the absence state
distinguishable in the type system. After the change, every read site must
explicitly handle the `None` case; the compiler enforces it.

## Non-goals

- No performance change (zero overhead expected; the migration is purely
  representational).
- No external API change (Python's `germline_position(idx) -> Optional[int]`
  surface is preserved).
- No related cleanups (we don't touch the `germline: u8` field, the
  `flags: NucFlags` discipline, or anything else). Stay focused.

## Design

### The new type

In `src/ir.rs`, alongside the existing `Nucleotide` definition:

```rust
/// Position in the source allele, or `NONE` for synthetic bases
/// (NP, P-nuc, contaminant, indel-inserted) with no germline
/// provenance. Replaces the previous `germline_pos: u16` field
/// where `u16::MAX` was a magic sentinel.
#[derive(Copy, Clone, Eq, PartialEq, Hash, Debug)]
#[repr(transparent)]
pub struct GermlinePos(u16);

impl GermlinePos {
    /// Sentinel value meaning "no germline provenance".
    pub const NONE: Self = Self(u16::MAX);

    /// Construct a real allele position. Panics if `v == u16::MAX`
    /// — that value is reserved for `NONE`.
    pub const fn pos(v: u16) -> Self {
        assert!(v != u16::MAX, "GermlinePos::pos: u16::MAX is reserved for NONE");
        Self(v)
    }

    /// Project to `Option<u16>`. `Some(v)` for a real position,
    /// `None` for `NONE`.
    pub const fn get(self) -> Option<u16> {
        if self.0 == u16::MAX { None } else { Some(self.0) }
    }

    pub const fn is_none(self) -> bool { self.0 == u16::MAX }
    pub const fn is_some(self) -> bool { !self.is_none() }
}
```

`#[repr(transparent)]` keeps the layout identical to a raw `u16`, so
`Nucleotide` stays at 8 bytes. The `pub const NONE` and `pub const fn pos`
are `const` so they compose into existing `const fn` constructors
(`Nucleotide::germline`, `Nucleotide::synthetic`).

### `Nucleotide` struct change

```rust
// before
pub struct Nucleotide {
    pub base: u8,
    pub germline: u8,
    pub germline_pos: u16,    // ← magic u16::MAX
    pub segment: Segment,
    pub flags: NucFlags,
}

impl Nucleotide {
    pub const NO_GERMLINE_POS: u16 = u16::MAX;   // ← deleted
    // ...
}

// after
pub struct Nucleotide {
    pub base: u8,
    pub germline: u8,
    pub germline_pos: GermlinePos,    // ← typed
    pub segment: Segment,
    pub flags: NucFlags,
}
// Nucleotide::NO_GERMLINE_POS constant deleted — sentinel now lives on the type.
```

Field stays `pub` (consistent with siblings `base`, `germline`, `segment`,
`flags`); the type itself enforces the discipline.

### Constructor signatures

`Nucleotide::germline(base, germline_pos: u16, segment)` keeps its `u16`
signature. Internally it calls `GermlinePos::pos(germline_pos)`, which
panics on `u16::MAX`. This way:

- Callers that already pass a known-real u16 position (`pos as u16`, `0`,
  `42`, etc.) don't need to wrap.
- The wrapping happens at the type boundary in one place.
- A misuse (`Nucleotide::germline(b'A', u16::MAX, Segment::V)`) panics at
  the boundary, not silently flowing into the pool.

`Nucleotide::synthetic` already uses `Self::NO_GERMLINE_POS` internally;
this becomes `GermlinePos::NONE` and the public signature is unchanged.

This keeps the migration churn minimal — callers of these constructors are
unchanged — while still enforcing the discipline at every read site.

### Migration patterns

59 grep hits in 18 files, but most are doc comments, function names, or
test-fixture parameter names that contain "germline_pos" as a string but
don't need code changes. Actual code-changes count to ~25 sites across
8 distinct patterns:

| Pattern | ~Sites | Old | New |
|---|---|---|---|
| 1. Sentinel check (boolean) | 4 | `n.germline_pos == Nucleotide::NO_GERMLINE_POS` | `n.germline_pos.is_none()` |
| 2. Read + extract (combined gate) | 4 | `let p = nuc.germline_pos as usize` (under prior `!= MAX` check) | `let Some(p) = nuc.germline_pos.get() else { return …; }; let p = p as usize;` |
| 3. Constructor body update | 2 | `Self::NO_GERMLINE_POS` inside `Nucleotide::synthetic`; raw u16 inside `Nucleotide::germline` | `GermlinePos::NONE`; `GermlinePos::pos(germline_pos)` |
| 4. Compare against computed u16 | 2 | `nuc.germline_pos != expected_germline_pos` (where `expected: u16`) | `nuc.germline_pos != GermlinePos::pos(expected_germline_pos)` |
| 5a. `assert_eq!` against `NO_GERMLINE_POS` | 4 | `assert_eq!(n.germline_pos, Nucleotide::NO_GERMLINE_POS)` | `assert!(n.germline_pos.is_none())` |
| 5b. `assert_eq!` against literal u16 | 6 | `assert_eq!(n.germline_pos, 12)` | `assert_eq!(n.germline_pos, GermlinePos::pos(12))` |
| 6. Python projection (PyO3) | 1 | manual `if == MAX { None } else { Some(_) }` | `n.germline_pos.get()` |
| 7. Doc comments referencing `u16::MAX` | 3 | prose mentions `u16::MAX` / `Nucleotide::NO_GERMLINE_POS` | prose mentions `GermlinePos::NONE` |

Patterns 1, 2, and 4 are the **load-bearing** changes — the migration
exists for these. Pattern 5b is the largest by count but is mechanical
churn in test assertions; the wrapper just makes literal-position
expectations explicit.

Pattern 2 delivers the most safety win: the previous "gate-then-extract"
form (`if != MAX { use as usize }`) becomes a single
`let Some(…) = ….get() else { … }`, making it structurally impossible to
extract without first handling absence.

Items NOT touched (despite matching grep):
- `passes/echo.rs` has its own `germline_pos: u16` field on a test-helper
  type that flows into `Nucleotide::germline(base, germline_pos, segment)`.
  Since the constructor keeps its `u16` signature, this helper is
  unchanged.
- `pass.rs` test-helper builders (multiple `germline_pos: 0`, `: 1`, etc.
  in test fixtures) — same reason.
- Function names and identifiers like `productive_uses_germline_pos_anchor_after_indels`
  in test code — they describe the *concept*, which is unchanged.

### Files touched

**Production code (7 files):**
- `src/ir.rs` — type definition, constant removal, constructor body
  updates (patterns 3, 7)
- `src/live_call/walker.rs` — sentinel checks + extracts (patterns 1, 2, 7)
- `src/airr_record/walk.rs` — sentinel check + extract (patterns 1, 2)
- `src/airr_record/mod.rs` — sentinel check (pattern 1)
- `src/contract/anchor_preserved.rs` — compare against computed value
  (pattern 4)
- `src/compiled/execute.rs` — doc comment only (pattern 7)
- `src/python/simulation.rs` — PyO3 projection (pattern 6)

**Test code (9 files):**
- `src/ir_tests.rs` — pattern 5a, 5b, plus four new tests for the newtype
  itself (see below)
- `src/compiled_tests.rs` — pattern 7 (doc comments / fn name reference)
- `src/live_call_tests.rs` — pattern 7 (doc comment)
- `src/passes/corrupt/indel.rs` inline test — pattern 5a, 7
- `src/passes/generate_np.rs` inline test — pattern 5a
- `src/passes/sample_base.rs` inline test — pattern 5a
- `src/passes/assemble_segment.rs` inline test — pattern 5b ×2, 7
- `tests/c10_recombination_integration.rs` — pattern 5a ×2, 5b ×3
- `tests/e8_property_invariants.rs` — pattern 5a, 7

### New tests for the newtype

Four small `#[test]` functions in `src/ir_tests.rs`:

1. `germline_pos_pos_rejects_max` — `GermlinePos::pos(u16::MAX)` panics (uses
   `#[should_panic]`).
2. `germline_pos_none_projects_to_option_none` — `GermlinePos::NONE.get()`
   returns `None`.
3. `germline_pos_pos_round_trips` — `GermlinePos::pos(42).get()` returns
   `Some(42)`.
4. `nucleotide_size_unchanged` — `std::mem::size_of::<Nucleotide>() == 8`.
   Pins the layout guarantee from `#[repr(transparent)]` so a future
   refactor can't silently regress memory cost.

## Risk & rollout

- **Single commit, all-at-once migration.** The blast radius (25 sites) is
  small enough that an incremental rollout (keeping `NO_GERMLINE_POS` as a
  `#[deprecated]` alias during transition) would only add risk — review of
  half-migrated code is harder, and a forgotten site would silently use the
  old style.

- **Type-driven correctness.** Every misuse is a compile error: you can't
  index, compare to a `u16`, or pattern-match without going through the
  newtype's API. The full 607 cargo + 547 pytest existing suite serves as
  the regression net; no behaviour change should land here.

- **PyO3 layer.** Python's `germline_position(idx) -> Optional[int]`
  projection is preserved. Internally the implementation simplifies from a
  manual `if == MAX { None } else { Some(_) }` to `n.germline_pos.get()`.

- **`#[repr(transparent)]` guarantee.** Ensures `Nucleotide` stays 8 bytes.
  Verified by `nucleotide_size_unchanged` test.

- **Estimate:** ~100 lines added, ~80 lines removed. Reviewable in under 15
  minutes. Implementable in a single ~1 hour session.

## What this design explicitly does NOT do

- Does not touch the `germline: u8` field. Synthetic bases set
  `germline = base` rather than using a sentinel, which is a different
  pattern and out of scope.
- Does not touch the `flags: NucFlags` discipline (e.g.,
  `NUC_FLAG_INDEL_INS` is the canonical "this is an inserted base" signal).
- Does not introduce new tests for the existing walker / pass behaviour.
  That's a separate phase (property tests for walker, deferred).
- Does not change PyO3 method signatures or Python-facing behaviour.
- Does not address other audit items (visibility tightening, panic
  reduction in hot loops). One audit item per phase keeps reviews tight.

## Success criteria

1. `cargo build --tests --lib` clean.
2. All 607 existing cargo tests pass.
3. All 547 pytest tests pass after `maturin develop --release`.
4. `Nucleotide::NO_GERMLINE_POS` constant is gone (verified by
   `grep -rn 'NO_GERMLINE_POS' engine_rs/` returning zero hits).
5. No remaining direct comparisons of `germline_pos` against `u16::MAX`
   outside the `GermlinePos` impl itself (verified by
   `grep -rE 'germline_pos.*u16::MAX|u16::MAX.*germline_pos' engine_rs/`
   returning zero hits — the `Self(u16::MAX)` inside `pub const NONE`
   doesn't match this pattern since it doesn't reference `germline_pos`).
6. `nucleotide_size_unchanged` test passes.
