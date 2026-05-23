# Metadata-truth test suite

GenAIRR's flagship guarantee is that every simulated sequence carries
metadata (allele calls, region bounds, junction coordinates, productive
flags, mutation counts, CIGAR strings, ...) that exactly reflects what a
perfect aligner / biologist / oracle would extract from the sequence
given only its bases.

This suite exists to catch any drift between the recorded metadata and
what the post-pass IR actually supports. It's the canonical place to
add new tests whenever you find a metadata-correctness edge case.

## Suite layout

| File                  | Scope                                                              |
|-----------------------|--------------------------------------------------------------------|
| `common.rs`           | Shared fixtures: curated VJ / VDJ refdata with built-in ambiguity. |
| `trim_ambiguity.rs`   | Trim widens the allele call set to all indistinguishable alleles.  |
| `np_extension.rs`     | NP random bases recreate trimmed germline → call narrows back.     |
| `cross_segment.rs`    | V↔D and D↔J overlap via matching prefixes/suffixes.                |
| `mutation_effects.rs` | SHM / PCR / quality / ncorrupt widening, narrowing, flipping calls.|
| `indel_effects.rs`    | Insertions / deletions; coordinate shifts; truth-allele retention. |
| `junction_anchor.rs`  | Junction bounds, frame, stop codon, productive flag invariants.    |
| `corruption_stack.rs` | Heavy combined pipeline; self-consistency of every AIRR field.     |
| `airr_projection.rs`  | Live-call → AIRR record surface; fallback / empty / format rules.  |

Each category file's top-of-module docstring lists the specific
invariants it should cover and points to adjacent in-tree tests that
exercise overlapping ground (those tests aren't duplicated here; they
keep running where they live).

## Why this lives in `engine_rs/tests/` and not inline under `src/`

In-tree tests under `engine_rs/src/*/tests/` are bound to one module's
private fixtures (`pub(super)` helpers, internal enum variants, etc.).
The metadata-truth suite needs to be self-contained at the public crate
boundary so a new contributor can read one folder and understand the
full set of guarantees the engine offers.

## Existing in-tree coverage (do not duplicate)

The following in-tree tests already cover substantial metadata-truth
ground. New tests in this suite should ADD coverage, not duplicate.

### V-end extension into NP
- [src/compiled/tests/v_extensions.rs](../../src/compiled/tests/v_extensions.rs)
  - `v_call_shrinks_when_np1_recreates_trimmed_suffix`
  - `v_call_stays_widened_when_np1_does_not_match_any_allele`
  - `v_call_partially_extends_when_np1_matches_only_a_prefix`
  - `v_call_extension_no_op_when_no_trim`
  - `append_region_np1_bumps_live_call_version`

### D-end extension into NP (both sides)
- [src/compiled/tests/d_extensions.rs](../../src/compiled/tests/d_extensions.rs)
  - `d_call_shrinks_when_np1_recreates_trimmed_prefix`
  - `d_call_shrinks_when_np2_recreates_trimmed_suffix`
  - `d_call_shrinks_via_both_sides_simultaneously`
  - `d_call_stays_widened_when_neither_np_matches`

### J-end extension into NP
- [src/compiled/tests/j_extensions.rs](../../src/compiled/tests/j_extensions.rs)
  - `j_call_shrinks_when_np1_recreates_trimmed_prefix_vj`
  - `j_call_stays_widened_when_np1_does_not_recreate_prefix_vj`
  - `j_call_partially_extends_when_np1_matches_only_a_suffix_of_prefix`
  - `j_left_extension_works_for_vdj_chain_via_np2`

### Cross-segment overlap
- [src/compiled/tests/overlaps.rs](../../src/compiled/tests/overlaps.rs)
  - `v_right_overlaps_d_when_d_starts_with_v_suffix`
  - `v_right_does_not_overlap_when_d_does_not_match_v_suffix`
  - `d_left_overlaps_v_when_v_ends_with_d_prefix`
  - `overlap_walker_halts_at_pool_end`

### Indel call updates
- [src/compiled/tests/indels.rs](../../src/compiled/tests/indels.rs)
  - `deletion_inside_v_widens_live_call_when_distinguishing_base_removed`
  - `deletion_preserves_call_when_distinguishing_base_remains`
  - `insertion_inside_v_does_not_fail_the_call`
  - `combined_insertion_and_deletion_recomputes_correctly`

### End-to-end curated metadata (30+ tests, GOLD STANDARD)
- [src/compiled/tests/curated.rs](../../src/compiled/tests/curated.rs)
  - Trim → widen call → NP recreate → narrow back (V/D/J each)
  - Mutation switches V call to a different allele
  - Indel widens / preserves calls
  - Germline span / sequence end / alignment end / CIGAR / NP-string
    cross-validation
  - Junction locates anchors via germline_pos under indel
  - `curated_full_corruption_stack_keeps_metadata_self_consistent`

### Walker observer ↔ oracle equivalence
- [src/live_call/tests/walker_observer.rs](../../src/live_call/tests/walker_observer.rs)
  - 10 tests pairing observer with `call_from_region` for trim-ambiguity,
    mutation-flip, N-wildcard, NP-extension, cross-segment, indel-skip.
- [src/live_call/tests/walker_observer_change_base.rs](../../src/live_call/tests/walker_observer_change_base.rs)
- [src/live_call/tests/walker_observer_from_existing_region.rs](../../src/live_call/tests/walker_observer_from_existing_region.rs)

### Caller invariants
- [src/live_call/tests/caller.rs](../../src/live_call/tests/caller.rs)
  - 8 tests including `assembled_segment_live_call_keeps_trim_induced_ambiguity`
  - `y6_truth_allele_retained_under_single_position_mutation`
  - `y6_v_call_diverges_from_truth_when_mutations_flip_toward_another_allele`
  - `y6_trim_ambiguity_narrows_via_np_extension`

### AIRR projection
- [src/airr_record/tests/projection.rs](../../src/airr_record/tests/projection.rs)
- [src/airr_record/tests/anchors.rs](../../src/airr_record/tests/anchors.rs)

### Property invariants over many seeds
- [tests/e8_property/](../e8_property/) — codon_rail, regions, indel,
  productive, persistent, trace, determinism, distribution_support.

### Golden regressions
- [tests/golden/](../../../tests/golden/) — TSV snapshots of small
  fixed-seed pipelines for `human_igh`, `human_igk`, `human_tcrb`, and a
  20-record mixed-feature `human_igh` set.

## Running the suite

```bash
cd engine_rs && cargo test --release --test metadata_truth_invariants
```

## Adding new tests

1. Pick the category file that fits the invariant you're testing.
2. Read its top-of-module docstring for scope and adjacent coverage.
3. Use the fixtures in `common.rs` when they fit; build a one-off
   `RefDataConfig` locally when the curated ambiguity layout doesn't
   match your case.
4. Name the test for the invariant it asserts, not the mechanism:
   - Good: `truth_allele_retained_after_three_silent_mutations`
   - Bad:  `test_mutation_3`
5. Prefer asserting on AIRR-level field values (the user-facing
   contract) over engine internals. Internals can shift; AIRR fields
   are the published surface.
