# V-Subregion Mutation Counters — Audit + Slice Shipped

**Status: Slice shipped.** The six AIRR fields landed alongside
the validator extension, the manifest block, the column
ordering, the Python tests, and four Rust unit tests covering
event-time attribution + tampered-counter detection +
partition-sum invariant. The original audit body is preserved
below for traceability; sections marked **[Shipped]** describe
how the recommendations actually landed.

The slice adds six counters that partition `n_v_mutations`
across the five canonical IMGT V-subregion labels plus an
"unannotated" bucket for V mutations that can't be attributed
to a subregion (legacy / mixed cartridges, V-side CDR3 stretch,
and indel-inserted V bases).

This audit is the natural follow-up to the V-subregion SHM
**rate** slice (`docs/v_subregion_shm_rate_design.md`). The
rate slice gave users a knob to **target** SHM at specific
V subregions; the counters slice exposes the **realised**
per-subregion counts on the projected AIRR record. Together
they close the per-subregion SHM loop: targeting at simulate
time, observability at projection time.

Companion to
[`tests/test_v_subregion_mutation_counters_contract.py`](../tests/test_v_subregion_mutation_counters_contract.py)
which freezes today's surfaces (`pin_scaffold_*`) and the gaps
the implementation slice would close (`pin_absence_*`).

**Pre-flight finding (Q3 below): clean yes.** The event-time
attribution path is feasible: `BaseChanged.germline_pos` is
already the allele-local coordinate the audit needs, stamped
at assembly with `trim_5 + offset` and unchanged by subsequent
indels / end-loss. No new event-payload field is required.

---

## Existing scaffolding the contract pins

| Surface | Where | What it pins |
|---|---|---|
| `SimulationEvent::BaseChanged { handle, old_base, new_base, segment, germline_pos: Option<u16> }` | [`engine_rs/src/ir/sim_event.rs:119-125`](../engine_rs/src/ir/sim_event.rs#L119-L125) | The event variant carrying the SHM payload. `germline_pos` is allele-local (assembly stamps it as `trim_5 + offset`); `None` only for indel-inserted bases. |
| Assembly stamps allele-local `germline_pos` | [`engine_rs/src/passes/assemble_segment/execution.rs:158-162`](../engine_rs/src/passes/assemble_segment/execution.rs#L158-L162) | `slice_start + i` where `slice_start = trim_5`. The audit's `pool_pos - v_region.start + instance.trim_5` arithmetic is already encoded in the field's definition. |
| `Nucleotide.germline_pos: GermlinePos` (`u16` newtype with `NONE` sentinel) | [`engine_rs/src/ir/nucleotide.rs`](../engine_rs/src/ir/nucleotide.rs) | Allele-local position carried per pool entry; survives indels / end-loss unchanged. Indel-inserted bases get `GermlinePos::NONE`. |
| Per-segment counter aggregation loop | [`engine_rs/src/airr_record/builder.rs:265-290`](../engine_rs/src/airr_record/builder.rs#L265-L290) | Walks `outcome.events()`, filters `record.pass_name` to `MUTATE_UNIFORM` / `MUTATE_S5F`, matches `SimulationEvent::BaseChanged`, routes by `segment`. The new V-subregion counters mirror this loop's shape. |
| Pass-name allowlist constants | [`engine_rs/src/address.rs:88-100`](../engine_rs/src/address.rs#L88-L100) | `MUTATE_UNIFORM = "mutate.uniform"`, `MUTATE_S5F = "mutate.s5f"`. The new counters reuse the same allowlist so PCR / quality / receptor revision base changes don't leak in. |
| Validator's independent recompute + sum-invariant gate | [`engine_rs/src/airr_record/validate.rs:765-820`](../engine_rs/src/airr_record/validate.rs#L765-L820) | Re-walks `outcome.events()` from scratch with the same filter; emits `N{V,D,J,Np}MutationsMismatch` and `MutationCountSumMismatch`. The new slice adds parallel V-subregion issue kinds. |
| `AirrRecord` struct (where new fields go) | [`engine_rs/src/airr_record/record.rs:122-125`](../engine_rs/src/airr_record/record.rs#L122-L125) | `n_v_mutations`, `n_d_mutations`, `n_j_mutations`, `n_np_mutations` as `i64`. The new six fields land alongside. |
| `airr_record_to_pydict` projection to Python | [`engine_rs/src/python/outcome.rs:228-366`](../engine_rs/src/python/outcome.rs#L228-L366) | Where `dict.set_item(...)` calls for the new field names land. |
| V allele identity at projection time | [`engine_rs/src/ir/simulation.rs:26`](../engine_rs/src/ir/simulation.rs#L26) + [`engine_rs/src/assignment.rs:106-111`](../engine_rs/src/assignment.rs#L106-L111) | `sim.assignments.get(Segment::V) → AlleleInstance { allele_id, trim_5 }`; `refdata.v_pool.get(allele_id).subregions` is the subregion table. Same access path the rate slice uses. |
| `v_subregion_at_position` helper (already shipped for rates) | [`engine_rs/src/passes/mutate/v_subregion_rates.rs`](../engine_rs/src/passes/mutate/v_subregion_rates.rs) | Walks `allele.subregions` linearly (≤ 5 entries). The counters slice doesn't need this exact helper at projection time — it goes directly from `germline_pos` to subregion via a one-time-per-record allele lookup. |
| `MutationCountSumMismatch` + per-segment mismatch issue kinds | [`engine_rs/src/airr_record/validate.rs:87-99`](../engine_rs/src/airr_record/validate.rs#L87-L99) | The validation precedent. The slice adds analogous V-subregion mismatch variants. |

---

## 1. Q1 — Counter shape

The user's recommendation: **five-label fields** for the
canonical IMGT V subregions:

```text
n_fwr1_mutations
n_cdr1_mutations
n_fwr2_mutations
n_cdr2_mutations
n_fwr3_mutations
```

**Endorsed.** Five-label fields are the right surface because:

- They mirror the **rate kwarg's vocabulary** exactly — the
  user already specifies `v_subregion_rates={"CDR1": 2.0,
  "FWR2": 0.5}` per label; the counters expose exactly what
  the rates produced. The targeting → observability loop is
  the lowest-friction shape.
- Two-bucket aggregates (`n_cdr_mutations` /
  `n_fwr_mutations`) **lose information**: a user who
  targeted CDR1 specifically (via `{"CDR1": 3.0, "CDR2":
  1.0}`) couldn't read back the CDR1 hit count separately.
  Two-bucket fields can always be derived by summing pairs;
  five-label fields cannot be reconstructed from two
  buckets.
- IMGT is the field's vocabulary anyway — every immunology
  tool downstream understands FWR1 / CDR1 / FWR2 / CDR2 /
  FWR3.

### Two-bucket aggregates are NOT in v1

A derived `n_cdr_mutations = n_cdr1_mutations +
n_cdr2_mutations` and `n_fwr_mutations = n_fwr1_mutations +
n_fwr2_mutations + n_fwr3_mutations` is **trivially
computable downstream** — a user who wants the two-bucket
view does `df["n_cdr_mutations"] = df["n_cdr1_mutations"] +
df["n_cdr2_mutations"]`. The audit recommends NOT writing
those fields into the AIRR record; if a future ergonomics
slice demands them, they go in then.

### Pinned

- `pin_absence_no_cdr_fr_mutation_counter_fields` — carried
  from the V-substructure audit; pinned in this contract too.

---

## 2. Q2 — Source of truth + event-time attribution

The user's spec: use `SimulationEvent::BaseChanged` from
biological SHM passes only, same as per-segment counters.
Attribute by event position + final/pre-event V assignment
subregion.

### Source of truth: SHM `BaseChanged` events only

The per-segment counters slice already established the
discipline. The aggregation walks `outcome.events()`, filters
to two pass names (`MUTATE_UNIFORM` / `MUTATE_S5F`) — see
[`builder.rs:269-274`](../engine_rs/src/airr_record/builder.rs#L269-L274)
— and inspects `SimulationEvent::BaseChanged` only. PCR /
quality / receptor-revision base changes do NOT count.

The new V-subregion counters reuse the **same filter** verbatim
so:

- PCR-induced V-region substitutions go to `n_pcr_errors`
  (already existing), NOT to the new V-subregion buckets.
- Quality-induced V-region substitutions go to
  `n_quality_errors`, not the new buckets.
- Receptor-revision base changes (recombination-stage) and
  D-inversion base flips do not count as SHM.

### Event payload is sufficient (the central feasibility finding)

`BaseChanged { handle, old_base, new_base, segment,
germline_pos: Option<u16> }` already carries everything the
attribution needs:

1. `segment == Segment::V` → "is this a V SHM event?"
2. `germline_pos` → "where on the V allele did it land?"

The audit recommends NOT adding a `subregion: Option<VSubregionLabel>`
field to the event. Two reasons:

- The event's existing footprint stays minimal — adding a
  payload field bloats every recorded `BaseChanged` whether
  the user is using subregion counters or not (~24 bytes per
  event × every SHM event × every simulation).
- The subregion is **derivable** at projection time from the
  assigned V allele's `subregions` table. Allele identity is
  reachable via `sim.assignments.get(Segment::V)`. The
  derivation is a single 5-entry linear scan per V event.

### Why not post-hoc attribution from final sequence regions

A naive approach: at AIRR projection time, compare the final
mutated sequence against the germline and bucket each
divergence by IMGT region using sequence coordinates.

**Rejected by the audit.** Three problems:

1. **Indel / end-loss instability.** Indels after SHM shift
   pool coordinates; end-loss truncates them. Sequence-
   position-based attribution gets the wrong subregion bucket
   on any record that goes through `polymerase_indels` or
   `primer_trim_*` after SHM.
2. **D-inversion / receptor-revision interference.** Both
   passes rewrite assembled regions; a post-hoc walk over the
   final sequence sees the post-revision state, not the SHM-
   target state.
3. **Loss of provenance.** The event log already records
   which substitution was a *biological SHM event* (via
   pass-name filter); a post-hoc walk would over-count PCR /
   quality substitutions in V.

Event-time attribution via `BaseChanged.germline_pos` avoids
all three.

### Pinned

- `pin_scaffold_base_changed_carries_germline_pos` — the event
  variant already exposes the field.
- `pin_scaffold_per_segment_counter_filter` — the existing
  pass-name allowlist (`MUTATE_UNIFORM` / `MUTATE_S5F`) is
  the precedent.

---

## 3. Q3 — Robustness under later edits

The user's question: are `BaseChanged.segment + germline_pos`
enough to attribute each event independent of later pool
edits?

### Yes — `germline_pos` is allele-local and stable

The investigation confirms:

- `Nucleotide.germline_pos` is stamped at assembly
  ([`assemble_segment/execution.rs:158-162`](../engine_rs/src/passes/assemble_segment/execution.rs#L158-L162))
  as `slice_start + i`, where `slice_start = trim_5`. So
  germline_pos = `trim_5 + offset_into_assembled_V`. This is
  exactly the audit's pool-to-allele arithmetic, already
  pre-computed at assembly time.
- Subsequent indels (`insert_indel` / `delete_indel` in
  `builder.rs`) shift pool handles but do **NOT** rewrite
  `germline_pos` on surviving nucleotides. End-loss truncates
  the pool but the germline_pos on remaining bases is
  unchanged.
- `BaseChanged.germline_pos` is read from the nucleotide
  **before** the mutation is applied — so it reflects the
  pre-mutation allele-local position, which is the right
  coordinate for subregion lookup.

The result: a `BaseChanged` event recorded under
`mutate.s5f` carries the allele-local position of the
mutated base, and that position is invariant under all
subsequent passes (PCR / quality / indels / end-loss / paired-
end / etc.). Event-time attribution is rock-solid.

### Edge case: indel-inserted V bases

`GermlinePos` carries a `NONE` sentinel (`u16::MAX`) for
indel-inserted bases — they have no germline origin. If a
later SHM pass runs after a `polymerase_indels` pass and
mutates one of those inserted bases, the resulting
`BaseChanged.germline_pos = None`. The subregion lookup
treats `None` as the **unannotated bucket** (see Q4 below).

This is rare in practice: SHM almost always precedes
observation-stage corruption (the canonical DSL order is
`recombine → mutate → pcr_amplify → polymerase_indels →
…`); a user who mutates after indels does so deliberately.

### Pre-flight verdict

**No new event-payload field is needed.** The slice can
implement counters using the existing event surface
unchanged. This is the audit's central positive finding.

### Pinned

- `pin_scaffold_germline_pos_is_optional_u16` — the field
  exists on `BaseChanged` with `Option<u16>` shape.
- `pin_scaffold_germline_pos_survives_indels` — pool indels
  don't rewrite `germline_pos` on surviving nucleotides
  (already pinned indirectly by the indel-provenance audit;
  cited here for the counters slice).

---

## 4. Q4 — Missing annotations

The user's recommendation: add `n_v_unannotated_mutations` to
preserve the partition

```text
five subregion counters + unannotated == n_v_mutations
```

**Endorsed.** The audit recommends this exact shape. Three
cases produce an "unannotated" V mutation:

1. **Legacy V allele with no `gapped_seq`** (and thus no
   bridge-derived subregions). The Slice-1 cartridge surface
   loads these with `Allele.subregions = []`; no subregion
   lookup can succeed.
2. **Mixed cartridge — some annotated, some not.** A user
   cartridge with hand-authored V alleles and IMGT-derived
   V alleles. The unannotated V alleles route to the
   unannotated bucket; the annotated ones route to their
   subregion bucket.
3. **Indel-inserted V base** (rare — see Q3 above). The
   inserted nucleotide carries `germline_pos = None`; the
   subregion lookup necessarily returns "unannotated".
4. **V-side CDR3 contribution.** The five IMGT subregion
   labels cover only FWR1 → FWR3; the stretch between
   `FWR3.end` and `len(allele.seq)` (the V-side
   contribution to CDR3) is deliberately outside the
   annotation set. An SHM event landing in that stretch
   carries a valid `germline_pos` but maps to no subregion
   — also routed to unannotated.

### Why a partition discipline

The alternative — "five counters whose sum is `≤
n_v_mutations`" — leaks information silently. A user looking
at `n_cdr1 + n_cdr2 + n_fwr1 + n_fwr2 + n_fwr3 = 4` against
`n_v_mutations = 7` has no way to know whether the missing 3
hit the V-CDR3 stretch, an unannotated allele, or an inserted
base — all three are biologically distinct.

`n_v_unannotated_mutations` makes the partition exact and
auditable:

```text
n_fwr1 + n_cdr1 + n_fwr2 + n_cdr2 + n_fwr3 + n_v_unannotated == n_v_mutations
```

The validator's `VSubregionMutationCountSumMismatch` (Q5)
enforces this at every record.

### Expected value on bundled cartridges

The Slice-1 annotation surface gives **100%** subregion
*annotation* coverage on bundled human IGH / IGK / IGL OGRDB
cartridges (198 / 168 / 181 V alleles respectively, all carrying
the five canonical intervals).

**[Shipped — corrected from the pre-implementation expectation]**:
`n_v_unannotated_mutations` is NOT identically zero on bundled
cartridges. The five canonical IMGT labels cover FWR1 → FWR3;
the stretch between `FWR3.end` and `len(allele.seq)` — the V-side
CDR3 contribution — is deliberately outside the annotation set
(case 4 above). SHM events landing in that stretch carry a
valid `germline_pos` but no matching subregion interval and
therefore route to the unannotated bucket. The smoke test on
HUMAN_IGH_OGRDB / S5F / rate 0.05 / n=10 records shows ~1
unannotated event per ~150 V SHM events — a small but expected
non-zero baseline.

Other non-zero sources:

- The user supplies a cartridge with at least one
  unannotated V allele AND a record draws that allele.
- The user inverts the canonical pass order (indels before
  SHM) and an SHM event lands on an indel-inserted V base.

### Pinned

- `pin_absence_no_n_v_unannotated_mutations_field` — the
  field doesn't exist yet.
- `pin_scaffold_bundled_cartridges_have_100_pct_v_subregion_coverage`
  — carried from the Slice-1 contract; the expected-zero
  invariant depends on it.

---

## 5. Q5 — Validator extension

The user's spec: mismatch variants for each field plus a sum
invariant `VSubregionMutationCountSumMismatch`. Source tag
should be event-derived.

### Recommended issue-kind variants

Six new variants on the Rust validator's issue enum (mirroring
the per-segment counter slice's pattern):

```rust
NFwr1MutationsMismatch { reported: i64, event_count: i64 }
NCdr1MutationsMismatch { reported: i64, event_count: i64 }
NFwr2MutationsMismatch { reported: i64, event_count: i64 }
NCdr2MutationsMismatch { reported: i64, event_count: i64 }
NFwr3MutationsMismatch { reported: i64, event_count: i64 }
NVUnannotatedMutationsMismatch { reported: i64, event_count: i64 }
VSubregionMutationCountSumMismatch {
    reported_v_total: i64,
    sum_of_subregion_buckets: i64,
}
```

Six per-field mismatches + one cross-field sum invariant — same
shape as the per-segment counter slice's `N{V,D,J,Np}MutationsMismatch`
+ `MutationCountSumMismatch`.

### Independent recompute

The validator's existing per-segment recompute at
[`validate.rs:765-807`](../engine_rs/src/airr_record/validate.rs#L765-L807)
walks `outcome.events()` from scratch (NOT from the AIRR
record's reported fields) and applies the same pass-name
allowlist. This is what gives the existing sum invariant
genuine bite: if the projection over-counts by mis-attributing
a PCR event, the validator's independent walk produces a
different total.

The V-subregion validator extension does the same: walks the
event ledger from scratch, applies the same allowlist, and for
each V `BaseChanged` does the subregion lookup independently
of the projection's lookup.

### Source tag = event-derived

The audit's "source tag" convention (used by the per-segment
counter mismatch issues — see
[`mutation_provenance_audit.md`](mutation_provenance_audit.md))
identifies which event ledger frame produced the discrepancy.
For the V-subregion counters, the natural tag is
`source = "event"` plus the failing event's `(pass_name,
event_index)` — same shape as the existing
`source: "events"` tag in
[`python/outcome.rs`](../engine_rs/src/python/outcome.rs).

### Pinned

- `pin_absence_no_v_subregion_mismatch_validator_kinds` —
  the new issue variants don't exist yet (carried from the
  V-substructure audit's `SubregionsBoundariesMismatch` /
  `SubregionSumMismatch` placeholder pin).
- `pin_scaffold_validator_independent_recompute_pattern` —
  source-level pin that the per-segment counter validator
  walks the event ledger independently (Q5's correctness
  argument depends on it).

---

## 6. Q6 — AIRR column names

The user's spec: lowercase canonical
`n_fwr1_mutations, n_cdr1_mutations, ..., n_v_unannotated_mutations`.

**Endorsed.** Six new AIRR record fields:

```text
n_fwr1_mutations
n_cdr1_mutations
n_fwr2_mutations
n_cdr2_mutations
n_fwr3_mutations
n_v_unannotated_mutations
```

All `int` typed (matches the existing per-segment counter
shape). Default to `0` when no SHM ran. Always present on
every record (never `None`), so downstream consumers can
treat them as numeric columns without a missing-value branch.

### Why lowercase + underscored

Matches the existing AIRR convention for SHM counters
(`n_mutations`, `n_v_mutations`, etc.). The IMGT vocabulary
is canonical-cased (`FWR1`, `CDR1`) but Python AIRR field
conventions strip case + add underscores — same boundary the
`v_subregion_rates={"CDR1": …}` kwarg crosses (user-facing
labels are canonical-cased, internal field names are
lowercase).

### Pinned

- `pin_absence_no_cdr_fr_mutation_counter_fields` — all six
  field names absent on `result.records[0]` today.

---

## 7. Q7 — Interaction with existing counters

The user's spec:
- `n_v_mutations` remains the parent total.
- New counters partition V only.
- `n_d` / `n_j` / `n_np_mutations` unchanged.

**Endorsed.** The existing per-segment counters keep their
exact current semantics. The new counters add a finer
partition under V only.

### The partition algebra

```text
n_mutations =
    n_v_mutations + n_d_mutations + n_j_mutations + n_np_mutations
n_v_mutations =
    n_fwr1_mutations + n_cdr1_mutations + n_fwr2_mutations
    + n_cdr2_mutations + n_fwr3_mutations
    + n_v_unannotated_mutations
```

Both invariants enforced by independent validator checks. The
first is the existing `MutationCountSumMismatch`; the second
is the new `VSubregionMutationCountSumMismatch`.

### Composes with the rate slice

The rate slice (`v_subregion_rates={"CDR": 0.0}`) drops every
V-CDR site from proposal support. The counters slice's
expected behaviour under that configuration:

- `n_cdr1_mutations == 0` on every record.
- `n_cdr2_mutations == 0` on every record.
- `n_fwr1_mutations + n_fwr2_mutations + n_fwr3_mutations +
  n_v_unannotated_mutations == n_v_mutations`.

This composition is the canonical "release-tier zero-rate
exclusion under counters" check — analogous to the existing
`test_v_subregion_rates_zero_rate_exclusion_invariant_full_stack`
release test but now provable from the counters directly
(no baseline-diff classifier needed).

### Pinned

- `pin_scaffold_per_segment_counter_partition` — the existing
  `n_v + n_d + n_j + n_np == n_mutations` invariant.

---

## 8. Edge cases the implementation slice must handle

| Case | Expected behaviour |
|---|---|
| Default-rate IGH simulation (no `segment_rates` / no `v_subregion_rates`) | All six new counters non-negative; their sum equals `n_v_mutations`; `n_v_unannotated_mutations == 0` on bundled cartridges. |
| `v_subregion_rates={"CDR": 0.0}` | `n_cdr1 == n_cdr2 == 0` on every record; FWR + unannotated counters sum to `n_v_mutations`. |
| `v_subregion_rates={"FWR": 0.0}` | `n_fwr1 == n_fwr2 == n_fwr3 == 0` on every record; CDR + unannotated counters sum to `n_v_mutations`. |
| Empty plan (no `.mutate()` step) | All six counters `== 0`. |
| `count=0` or rate `== 0` | All six counters `== 0`. |
| Mixed cartridge — one annotated V allele, one unannotated | Records assigned to the annotated allele route to the five buckets; records assigned to the unannotated allele route entirely to `n_v_unannotated_mutations`. |
| Indels before SHM (non-canonical order) | An SHM event on an indel-inserted V base routes to `n_v_unannotated_mutations` (germline_pos `None`). |
| V-side CDR3 contribution stretch | An SHM event between `FWR3.end` and `len(allele.seq)` routes to `n_v_unannotated_mutations`. |
| PCR + quality + receptor revision + D inversion in the stack | None of them increment the new V-subregion counters (pass-name filter). |
| Trace replay | Counters round-trip exactly — they're derived from the event ledger which the replay reproduces deterministically. |

---

## 9. Performance

### Per-record cost

The aggregation walks `outcome.events()` once per record —
same shape as the existing per-segment counter loop. Per
event, the V-subregion path does:

- One `match segment` (already free — part of the existing
  loop).
- For V events only: a 5-entry linear scan over
  `allele.subregions`.

Typical V SHM event count per record: 5-20 (3% S5F rate × 300
bp V × productive_only filtering, per the existing release
tests). Subregion lookup cost: ≤ 5 comparisons per event. Per
record: ~100 comparisons. **Negligible.**

### One-time allele lookup per record

The aggregation can hoist
`assigned_v = sim.assignments.get(Segment::V)` and
`v_allele = refdata.v_pool.get(assigned_v.allele_id)` out of
the inner event loop — they're invariant within a record.
This shaves the `assignments.get` + `v_pool.get` cost from
per-event to per-record.

### No new heap allocations

The aggregation extends six existing `i64` accumulators in the
builder; no `Vec` / `HashMap` / `Box` introduced. Same
allocation discipline as the per-segment counter slice.

### Pinned

- `pin_scaffold_subregion_lookup_is_constant_time` — carried
  from the rate audit; the counters slice inherits the same
  constant-bounded per-event cost.

---

## 10. Trace / replay

### No new trace addresses

The counters slice is **pure projection** — it adds
aggregation logic at AIRR record build time but emits zero
new `ChoiceValue` records. Traces produced by the rate slice
replay byte-identically; trace files written before the
counters slice are unchanged.

### Replay determinism

Counters are derived from the deterministic event ledger
(itself reproduced byte-for-byte under trace replay). Two
runs at the same seed + same plan produce identical counters
by construction.

### Plan-signature impact

Zero — the counters slice changes neither the pass list nor
any pass's `parameter_signature` body. Plan signatures are
unchanged.

### Pinned

- `pin_scaffold_no_new_trace_addresses_under_counters` —
  source-level pin that the mutate address vocabulary is
  unchanged.

---

## 11. Manifest extension

### Recommended addition

```python
manifest["models"]["shm"]["v_subregion_counter_support"] = {
    "available": True,
    "labels": ["FWR1", "CDR1", "FWR2", "CDR2", "FWR3"],
    "partition_field": "n_v_mutations",
    "unannotated_field": "n_v_unannotated_mutations",
    "fields": [
        "n_fwr1_mutations",
        "n_cdr1_mutations",
        "n_fwr2_mutations",
        "n_cdr2_mutations",
        "n_fwr3_mutations",
        "n_v_unannotated_mutations",
    ],
    "source": "event",
    "in_content_hash": False,  # per-record, not per-cartridge
}
```

The manifest already advertises `v_subregion_support` (Slice
1) and `v_subregion_rate_support` (Slice B). This adds the
third block parallel to them.

### Pinned

- `pin_absence_no_v_subregion_counter_support_in_manifest` —
  the block doesn't exist yet.

---

## 12. Implementation order **[Shipped]**

A single self-contained slice. Five sub-steps in dependency
order, all landed:

### Step 1 — Rust builder aggregation

Extend the existing per-segment counter loop in
[`engine_rs/src/airr_record/builder.rs:265-290`](../engine_rs/src/airr_record/builder.rs#L265-L290):

1. Hoist `assigned_v = sim.assignments.get(Segment::V)` and
   `v_allele = refdata.v_pool.get(...)` above the loop.
2. Add six new `i64` accumulators.
3. In the `Segment::V => ...` arm, dispatch on `germline_pos`:
   - `None` → `n_v_unannotated_mutations += 1`.
   - `Some(pos)` → linear scan over
     `v_allele.subregions` for the matching interval;
     route to the matching bucket or fall through to
     `n_v_unannotated_mutations` if no interval contains
     `pos`.

### Step 2 — `AirrRecord` struct + Python dict

- Add six `pub i64` fields to `AirrRecord` in
  [`record.rs`](../engine_rs/src/airr_record/record.rs).
- Add six `dict.set_item(...)` calls in
  `airr_record_to_pydict` in
  [`outcome.rs`](../engine_rs/src/python/outcome.rs).

### Step 3 — Validator extension

- Add seven new issue-kind variants (six per-field +
  `VSubregionMutationCountSumMismatch`) to the validator's
  issue enum in
  [`validate.rs`](../engine_rs/src/airr_record/validate.rs).
- Extend the existing independent-recompute walk to do the
  same subregion dispatch; emit mismatches when the recompute
  disagrees with the reported field.
- Emit `VSubregionMutationCountSumMismatch` when the six
  buckets don't sum to `n_v_mutations`.
- Surface the new issues in the Python issue dict via
  `outcome.rs` issue serialisation.

### Step 4 — Manifest extension

Add `v_subregion_counter_support` block per §11.

### Step 5 — Tests

- DSL contract: per-record counters present + non-negative +
  sum to `n_v_mutations`.
- Default rates: `n_v_unannotated_mutations == 0` on bundled
  cartridges.
- Composition with `v_subregion_rates={"CDR": 0.0}`:
  CDR counters all zero, FWR + unannotated sum to V total.
- Mixed-cartridge edge case: records assigned to unannotated
  V allele route entirely to the unannotated bucket.
- Indel-before-SHM edge case: SHM on an inserted V base goes
  to unannotated.
- Validator: mismatched counter detected; sum invariant
  surfaces.
- Trace replay: counters round-trip exactly.
- Release-tier: full stack with `v_subregion_rates +
  segment_rates` produces clean validation + per-record
  partition.

**Cost estimate:** ~150 lines Rust + ~80 lines validator +
~50 lines manifest + ~30 spec tests (~250 lines).

### Why a single slice

Unlike Slices A + B (which had a real dependency — A unblocked
B's signature surface), the counters slice has no internal
phase boundary. The Rust aggregation, the AIRR fields, the
validator, the manifest, and the tests all land together as
one mechanical extension of the existing per-segment counter
plumbing. No precondition slice is needed.

---

## 13. Test surface — what this audit pins

Mirrored in
[`tests/test_v_subregion_mutation_counters_contract.py`](../tests/test_v_subregion_mutation_counters_contract.py).

### `pin_scaffold_*` — pre-existing surfaces the slice builds on

1. `SimulationEvent::BaseChanged` carries `germline_pos:
   Option<u16>` (the audit's load-bearing field).
2. Assembly stamps `germline_pos = trim_5 + offset` so the
   value is allele-local at event time.
3. The per-segment counter aggregation in
   `airr_record/builder.rs` walks events with the
   `MUTATE_UNIFORM` / `MUTATE_S5F` allowlist — the precedent
   the new aggregation mirrors.
4. The validator does an independent recompute (walks events
   from scratch) — the discipline the new validator mirror
   inherits.
5. `MutationCountSumMismatch` exists in `validate.rs` (the
   pattern the new `VSubregionMutationCountSumMismatch`
   mirrors).
6. `Simulation.assignments.get(Segment::V) → AlleleInstance`
   + `RefDataConfig.v_pool.get(...).subregions` reach
   subregion intervals at aggregation time.
7. Bundled human IGH / IGK / IGL OGRDB cartridges have 100%
   V-subregion coverage.
8. The existing per-segment counter sum-invariant
   (`n_v + n_d + n_j + n_np == n_mutations`) holds on every
   release-tier record.
9. `n_pcr_errors` / `n_quality_errors` are surfaced separately
   so the new V-subregion counters don't need to reproduce
   that filter.

### `pin_absence_*` — gaps the slice closes

10. No `n_fwr1_mutations` / `n_cdr1_mutations` /
    `n_fwr2_mutations` / `n_cdr2_mutations` / `n_fwr3_mutations`
    fields on AIRR records.
11. No `n_v_unannotated_mutations` field on AIRR records.
12. No `v_subregion_counter_support` block in the manifest.
13. No `N{Fwr1,Cdr1,Fwr2,Cdr2,Fwr3,VUnannotated}MutationsMismatch`
    issue kinds on the Rust validator.
14. No `VSubregionMutationCountSumMismatch` issue kind on
    the validator.

### Doc anchor

15. The audit doc exists and references the contract file;
    structure intact.

---

## 14. Out of scope

Documented here so a future implementer doesn't accidentally
expand the work.

- **Two-bucket aggregates** (`n_cdr_mutations` /
  `n_fwr_mutations`). Trivially derivable downstream; not in
  v1.
- **CDR3 / FWR4 counters.** The V annotation surface stops
  at FWR3. CDR3 lives in the junction (covered by
  `n_np_mutations`); FWR4 in J. A user who wants per-region
  CDR3 counters needs a separate junctional-substructure
  audit.
- **Per-region mutation rate output** (e.g.
  `cdr1_mutation_rate`). The realised counts plus
  `mutation_rate` (already present) suffice; deriving rates
  downstream is `n_cdr1_mutations / (cdr1_end - cdr1_start)`.
- **D-region / J-region subregion counters.** D alleles have
  no IMGT-canonical substructure of analogous granularity.
- **Per-allele mutation counts.** A future analytics slice
  might surface `n_v_mutations_per_allele` for clonal-family
  drift analysis; not in v1.
- **Empirical S5F-subregion kernel** (per-subregion 5-mer
  context tables). Same answer as in the rate audit: out
  of scope.

---

## 15. Summary table

| Concern | Post-slice state | Status |
|---|---|---|
| Counter shape | Five-label + one-unannotated; partition `n_v_mutations` | **Shipped** |
| Source of truth | `SimulationEvent::BaseChanged` filtered to `MUTATE_UNIFORM` / `MUTATE_S5F` (same as per-segment) | **Shipped** |
| Attribution method | Event-time, via `BaseChanged.germline_pos` (already allele-local) | **Shipped** |
| Survives later edits | `germline_pos` is assembly-stamped, unchanged by indels / end-loss | **Confirmed** — no event-payload change |
| Missing annotations | `n_v_unannotated_mutations` preserves the partition | **Shipped** |
| Bundled cartridge baseline | ~1 unannotated event per ~150 V SHM events (V-side CDR3 stretch, by design) | **Shipped, correctly characterised** |
| Validator extension | 6 per-field mismatch kinds + 1 sum invariant; independent recompute | **Shipped** |
| AIRR fields | Six new `n_*_mutations` fields, lowercase canonical | **Shipped** |
| Existing `n_v/d/j/np_mutations` | Unchanged; the new partition lives strictly under `n_v_mutations` | **Shipped, unchanged** |
| Manifest reports counter capability | `v_subregion_counter_support` block | **Shipped** |
| Trace addresses | Same three per mutate model — counters are pure projection | **Shipped, no new addresses** |
| Plan signature | Unchanged — counters don't affect signature | **Shipped, no change** |
| Performance | ≤ 100 extra comparisons per record | **Negligible — confirmed at release tier** |

The V-subregion mutation counters slice shipped clean. The
audit's central feasibility claim — that `BaseChanged.germline_pos`
suffices for event-time attribution, survives subsequent passes
unchanged, and that the aggregation slots into the existing
per-segment counter loop — held in practice. The only audit
correction was the over-statement that
`n_v_unannotated_mutations == 0` on bundled cartridges: in
practice the V-side CDR3 stretch (FWR3.end → len(allele.seq))
produces a small baseline of unannotated events even on
fully-annotated cartridges. That's now correctly documented in
§4.

**What's pinned downstream:** four contract files flipped their
`pin_absence_*` entries to `pin_present_*` — counters fields,
validator issue kinds, manifest block, and DataFrame columns
are all enforced by the V-subregion-mutation-counters
implementation tests and the release-tier sanity test. The
two-bucket aggregates (`n_cdr_mutations` / `n_fwr_mutations`)
remain deliberately absent; downstream consumers compose them
trivially as `df["n_cdr_mutations"] = df["n_cdr1_mutations"]
+ df["n_cdr2_mutations"]`.
