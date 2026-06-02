# V-Subregion SHM Rate — Audit + Slices A & B Shipped

**Status: Slices A and B shipped.** Slice A (Pass Parameter
Signature) extended `pass_plan_signature` to fold each pass's
compile-time parameters and bumped the on-disk trace schema to
v3. Slice B (the `v_subregion_rates` kwarg on
`Experiment.mutate`) layers per-V-subregion SHM rate scalars on
top of the existing per-segment rates: a V site's final weight
is `base × segment_rate(V) × v_subregion_rate(label)`. Non-V
sites and unannotated V alleles receive identity factor `1.0`.
Pinned end-to-end by
[`tests/test_pass_parameter_signature.py`](../tests/test_pass_parameter_signature.py)
and
[`tests/test_v_subregion_rates_implementation.py`](../tests/test_v_subregion_rates_implementation.py).
Per-region SHM **counters** (`n_cdr1_mutations` etc.) remain
deferred and are pinned absent in the contract file.

The body below preserves the original pre-implementation
framing for traceability; sections marked **[Shipped]**
describe how the recommendation actually landed.

This is the natural follow-up to two prior slices: the
**Targeted SHM** slice ([`shm_segment_rate_design.md`](shm_segment_rate_design.md))
added per-segment rate scalars on the V / D / J / NP buckets,
and the **V-Subregion Cartridge Annotation Surface** slice
([`v_region_substructure_audit.md`](v_region_substructure_audit.md))
made per-V-allele IMGT FWR1 / CDR1 / FWR2 / CDR2 / FWR3
intervals a first-class cartridge property. Both shipped. The
remaining gap is connecting the two: letting users target SHM
to specific V subregions, matching how immunologists think
about somatic hypermutation (CDRs evolve fast, FWRs stay
conserved).

Companion to
[`tests/test_v_subregion_shm_rate_contract.py`](../tests/test_v_subregion_shm_rate_contract.py)
which freezes today's surfaces (`pin_scaffold_*`) and the gaps
the implementation slice would close (`pin_absence_*`).

**Pre-flight finding (Q3 below):** the mapping from a pool
position to the assigned V allele's subregion coordinates is
implementable as a clean sibling of `segment_at_position` with
all required inputs already in scope at the SHM weight-walk
call site. No architectural friction. The audit recommends
proceeding to implementation.

---

## Existing scaffolding the contract pins

| Surface | Where | What it pins |
|---|---|---|
| `Experiment.mutate(segment_rates=...)` | [`src/GenAIRR/experiment.py`](../src/GenAIRR/experiment.py) | Per-segment SHM rate kwarg. Validated by `_validate_segment_rates`. Buckets: `V` / `D` / `J` / `NP`. The architectural precedent for `v_subregion_rates`. |
| `SegmentRateWeights` | [`engine_rs/src/passes/mutate/segment_rates.rs`](../engine_rs/src/passes/mutate/segment_rates.rs#L32-L92) | Four-field rate vector. `is_default()` short-circuits to the pre-slice fast path. `rate_for(Segment)` does the lookup, folding NP1+NP2 into a single `np` bucket. |
| `segment_at_position(sequence, pos)` | [`engine_rs/src/passes/mutate/segment_rates.rs`](../engine_rs/src/passes/mutate/segment_rates.rs#L106-L116) | Linear scan over `sequence.regions` returning the `Segment` containing `pos`. Takes only `&Sequence` + `pos` — receives no allele identity. The natural sibling site for `v_subregion_at_position`. |
| `S5FMutationPass::build_profile` | [`engine_rs/src/passes/mutate/s5f/sampling.rs`](../engine_rs/src/passes/mutate/s5f/sampling.rs#L33-L67) | The single per-position weight walk. Calls `segment_at_position` once per non-zero-mutability position when `segment_rates_active`. Multiplies the segment factor into the S5F mutability before adding to the profile. Zero-weighted positions never enter the profile (constrain-before-propose). |
| `UniformMutationPass::with_segment_rates` | [`engine_rs/src/passes/mutate/uniform.rs`](../engine_rs/src/passes/mutate/uniform.rs) | The uniform-model equivalent: passes `Some(&self.segment_rates)` into `MutationTransaction::substitute_position_constrained`. Same constrain-before-propose ordering. |
| `Simulation.assignments: AlleleAssignments` | [`engine_rs/src/ir/simulation.rs:26`](../engine_rs/src/ir/simulation.rs#L26) | Carries the sampled `AlleleInstance` per segment (V/D/J/C). Set before the mutate pass runs. **The V allele identity is reachable** at every point where `build_profile` and `segment_at_position` execute. |
| `AlleleInstance { allele_id, trim_5, trim_3 }` | [`engine_rs/src/assignment.rs:106-111`](../engine_rs/src/assignment.rs#L106-L111) | `allele_id` indexes `refdata.v_pool`; `trim_5` is the number of leading bases trimmed off before assembly. Both are needed to map a pool position to the V allele's ungapped coordinate. |
| `Allele.subregions: Vec<VSubregion>` | [`engine_rs/src/refdata.rs:198-261`](../engine_rs/src/refdata.rs#L198-L261) | The Slice-1 substructure annotation. Coordinates are ungapped, allele-relative, half-open. Reached via `refdata.v_pool.get(instance.allele_id).subregions`. |
| `cartridge_manifest()["models"]["shm"]` | [`src/GenAIRR/dataconfig/data_config.py`](../src/GenAIRR/dataconfig/data_config.py) | Carries `segment_rate_support` and `v_subregion_support`. **No `v_subregion_rate_support` key.** |
| Trace addresses for mutate | [`engine_rs/src/address.rs:345-350`](../engine_rs/src/address.rs#L345-L350) | `MutateUniformCount`, `MutateUniformSite(u32)`, `MutateUniformBase(u32)`, `MutateS5fCount`, `MutateS5fSite(u32)`, `MutateS5fBase(u32)`. The segment-rates slice introduced no new addresses; this slice must not either. |
| AIRR mutation counter fields | [`engine_rs/src/airr_record/record.rs`](../engine_rs/src/airr_record/record.rs) | `n_v_mutations` / `n_d_mutations` / `n_j_mutations` / `n_np_mutations`. **No `n_cdr*_mutations` / `n_fwr*_mutations`.** That is the Slice-3 counters slice from the V-Substructure audit; out of scope here. |

---

## 1. Q1 — Rate model shape

The user's recommendation:

```text
final_weight(site) =
  base_model_weight(site)
  * segment_rate(segment_of(site))
  * v_subregion_rate(subregion_of(site))   # iff segment == V and annotation exists
```

**Endorsed.** Three independent multiplicative factors compose
the per-site weight; each is a no-op when its rate vector is at
the flat default. Sites outside V never look up a subregion
rate. V sites belonging to an allele with no subregion
annotation (Slice-1 says these load with empty
`subregions: Vec<VSubregion>`) receive **subregion factor
`1.0`** — same as today, no penalty.

### Why multiplicative and not additive

The segment-rate slice's audit §2 (`shm_segment_rate_design.md`)
established the multiplicative discipline: rate scalars compose
with the base model rather than overriding it, so a user who
sets `segment_rates={"V": 2.0}` doubles V-targeted SHM relative
to the base model's V vs non-V distribution without losing the
base model's within-V structure. The same argument applies one
level finer: `v_subregion_rates={"CDR1": 3.0}` triples CDR1
weight relative to the base model's within-V site distribution,
preserving the S5F context-dependence within CDR1.

### Default behaviour

When `v_subregion_rates` is omitted, the rate vector is the
flat default `{FWR1: 1.0, CDR1: 1.0, FWR2: 1.0, CDR2: 1.0,
FWR3: 1.0}` and the engine takes the pre-slice fast path: no
per-position subregion lookup, byte-identical sampling
distribution to today's `mutate()` call. Same pattern as
`SegmentRateWeights::is_default`.

### Zero-rate semantics

`v_subregion_rates={"CDR2": 0.0}` zeroes every position inside
the assigned V allele's CDR2 interval. The position drops out
of the proposal support before contract admissibility runs —
same constrain-before-propose ordering as the segment-rates
slice (audit §4). Positions in CDR1 / FWR1 / FWR2 / FWR3 are
unaffected.

### Pinned

- `pin_absence_no_v_subregion_rates_kwarg_on_mutate` — the kwarg
  doesn't exist yet. Already pinned in the V-Substructure
  contract; carried forward here.

---

## 2. Q2 — Feasibility: pool position → V-allele subregion

**The central question.** Can the mutate pass, while picking
SHM site weights, cleanly map a pool position to the assigned
V allele's subregion coordinates?

### Verdict: **Yes, clean.**

All four inputs needed for a `v_subregion_at_position` helper
are already in scope at every call site where segment rates
apply.

| Input | Source | Cite |
|---|---|---|
| `sequence: &Sequence` | Already passed to `build_profile`. | [`s5f/sampling.rs:36`](../engine_rs/src/passes/mutate/s5f/sampling.rs#L36) |
| `assignments: &AlleleAssignments` | On `Simulation`. `tx.peek().assignments` mirrors the same snapshot as `tx.peek().sequence`. | [`ir/simulation.rs:26`](../engine_rs/src/ir/simulation.rs#L26) |
| `refdata: &RefDataConfig` | On `PassContext`. Already split-borrowed by the constrained-mutation execution path. | [`passes/mutate/s5f/execution.rs:236`](../engine_rs/src/passes/mutate/s5f/execution.rs#L236) |
| `Allele.subregions: Vec<VSubregion>` | Slice-1 cartridge surface; reached via `refdata.v_pool.get(instance.allele_id).subregions`. | [`refdata.rs:205`](../engine_rs/src/refdata.rs#L205) |

### The pool-position → allele-local arithmetic

For a pool position `pos` inside the V region of an assembled
sequence:

```rust
let v_region = sequence.regions.iter().find(|r| r.segment == Segment::V)?;
let instance = assignments.get(Segment::V)?;
let allele_local: u32 = pos - v_region.start.index() + instance.trim_5 as u32;
let allele = refdata.v_pool.get(instance.allele_id)?;
allele.subregions.iter().find(|s| {
    (s.start as u32) <= allele_local && allele_local < (s.end as u32)
})
```

Linear over `subregions.len()` (always ≤ 5) — same cost shape
as `segment_at_position`'s linear scan over `regions` (always
≤ 5).

### Edge cases the helper must handle

| Case | Return | Why |
|---|---|---|
| `pos` is outside the V region | `None` | Site is in D / J / NP — caller skips subregion factor. |
| V region exists but `instance.trim_5 + (pos - v_start)` falls outside `[0, len(allele.seq))` | `None` | Defence-in-depth; shouldn't happen if assembly was correct. |
| `allele.subregions` is empty (legacy V allele) | Caller's policy — see Q4 below. | Annotation coverage gap. |
| `allele_local` falls between two subregions (e.g. CDR3 region, post-FWR3) | `None` | The five labels cover only FWR1 → FWR3; the FWR3-end → V-end stretch is unannotated by design (CDR3 lives in the junction, not the V allele). |

### What this slice does NOT need

- **No allele-local position cache.** The arithmetic is two
  subtractions and one add per call. No data structure.
- **No per-pool-position subregion index.** A linear scan over
  ≤ 5 intervals beats a hashmap at this scale.
- **No new event types.** Subregion classification is
  in-pass weight computation, not state.

### Pinned

- `pin_scaffold_segment_at_position_takes_only_sequence_and_pos` —
  today's classifier signature, so the slice's helper sits
  beside it as a clean parallel.
- `pin_scaffold_simulation_assignments_carries_allele_identity` —
  `Simulation.assignments` exposes per-segment `AlleleInstance`
  including `allele_id` + `trim_5`.
- `pin_scaffold_v_subregion_lookup_is_constant_time_per_call` —
  the V allele's `subregions` is at most 5 entries; linear scan
  is the right shape.

---

## 3. Q3 — DSL shape

The user's recommendation:

```python
.mutate(
    model="s5f",
    rate=0.03,
    segment_rates={"V": 1.0, "NP": 0.0},
    v_subregion_rates={"CDR1": 2.0, "CDR2": 3.0, "FWR": 1.0},
)
```

**Endorsed.** Mirrors `segment_rates`: a dict keyed by canonical
label, omitted keys default to `1.0`, the kwarg itself defaults
to `None` (flat). Validated at the DSL boundary the same way
(`_validate_v_subregion_rates`).

### Canonical labels

The dict accepts the **five canonical V-subregion labels**
(`FWR1`, `CDR1`, `FWR2`, `CDR2`, `FWR3`) — the same vocabulary
the cartridge annotation surface uses
([`engine_rs/src/refdata.rs:208-249`](../engine_rs/src/refdata.rs#L208-L249)).
Case-sensitive (parallel to the bridge's existing
case-sensitivity for `VSubregionLabel::from_str`).

### `"FWR"` and `"CDR"` aliases

The user's recommendation: support `"FWR"` as an alias for
FWR1/FWR2/FWR3 because users think in CDR vs framework, not in
the five-label subdivision.

**Endorsed for v1**, with the following discipline:

- `"FWR"` is shorthand for `{FWR1: x, FWR2: x, FWR3: x}` —
  applies the same multiplier to all three framework regions.
- `"CDR"` is shorthand for `{CDR1: x, CDR2: x}` — symmetry
  argument; the V-region CDRs are CDR1 + CDR2 (CDR3 lives in
  the junction and is **not** in the V annotation surface).
- A dict mixing aliases with specific labels — e.g.
  `{"FWR": 0.5, "FWR2": 2.0}` — **expands the alias first, then
  applies the specific label as an override**. So that example
  resolves to `{FWR1: 0.5, FWR2: 2.0, FWR3: 0.5, CDR1: 1.0,
  CDR2: 1.0}`. Override semantics rather than rejecting; lets
  the user say "all framework regions slow, except FWR2".
- An alias keyed twice (`{"FWR": 0.5, "FWR": 1.0}`) is a
  Python-side duplicate-key error before reaching validation.

### Why not require explicit labels only

User ergonomics. The literature talks about
"hypermutation in the CDRs vs the framework regions"; forcing
users to type out `{CDR1: 2.0, CDR2: 2.0, FWR1: 0.5, FWR2: 0.5,
FWR3: 0.5}` to express that intent is friction without
correctness benefit. The five-label expansion is unambiguous;
the validator pins the expansion rule with a clear example in
its error messages.

### Why not allow `"CDR3"` / `"FWR4"`

Out of scope. The V annotation surface stops at FWR3. CDR3 is a
junctional concept (V + NP1 + D + NP2 + J anchor → W/F),
covered by the existing `NP` segment rate bucket. FWR4 lives in
the J segment. A user who wants to target CDR3 SHM should use
`segment_rates={"NP": 5.0}` (or zero out NP1/NP2 trim
distributions upstream). The validator rejects `"CDR3"` /
`"FWR4"` with a message that points to the audit doc.

### Pinned

- `pin_absence_no_v_subregion_rates_kwarg_on_mutate` —
  the kwarg doesn't exist.
- `pin_absence_no_v_subregion_rate_validator_in_experiment_py` —
  `_validate_v_subregion_rates` doesn't exist.

---

## 4. Q4 — Missing-annotation policy

When `v_subregion_rates` is provided but the assigned V allele
has no subregions (empty `Vec<VSubregion>`):

### Compile-time check: cartridge has zero annotated V alleles

If the user passes a non-default `v_subregion_rates` against a
cartridge where **no** V allele has subregion annotations,
**reject at refdata-resolve / compile time** with a clear error
pointing at the cartridge manifest's
`v_subregion_support.annotated_v_count`. The user's intent is
unsatisfiable: no V site will ever see a subregion factor.

The check lives at the same DSL boundary as the existing zero-
weight-on-all-segments check (segment_rates audit §2 "all-zero
rejects with `StrictSamplingError`").

### Runtime fallback: mixed cartridges, annotated + unannotated V alleles

If **some** V alleles have subregions and **some** don't
(possible for user-authored cartridges; never for bundled IGH /
IGK / IGL today), an unannotated V allele's sites receive
**subregion factor 1.0** by default — same as a site outside V
entirely. No error, no warning.

Rationale: the cartridge advertises partial coverage; the user
asked for subregion targeting on the subset that supports it.
Silently falling back is the same composition the segment-rates
slice uses when a segment has zero candidate alleles (drops out
of proposal support but doesn't error).

### Future strict mode

A future ergonomics slice could add `strict_subregions=True` on
`mutate()` that escalates "assigned V has no annotations" from
silent fallback to a `StrictSamplingError` per record. Not in
v1. Pinned absent.

### The unannotated-stretch caveat

Even on an **annotated** V allele the five subregions cover
only the FWR1 → FWR3 stretch (positions 0 → ~ungapped pos 315
or wherever the allele's `FWR3.end` lands). Positions between
`FWR3.end` and the allele's actual `len()` — i.e. the CDR3
contribution from the V allele tail — are outside the
annotation set. The lookup returns `None` there and the caller
applies factor `1.0` (default fallback). The validator should
NOT reject this case — it's not a coverage bug, it's the IMGT
definition of where the V annotation ends.

### Pinned

- `pin_absence_no_v_subregion_rates_kwarg_on_mutate` (carried).
- `pin_absence_no_compile_time_rejection_for_unsatisfiable_rates` —
  no validator currently checks that subregion rates are
  satisfiable; the slice's responsibility.

---

## 5. Q5 — Productive-only interaction

### Constrain-before-propose ordering preserved

The site's final weight is

```
base_model_weight × segment_rate × v_subregion_rate
```

Zero in any factor zeros the site's proposal probability;
zero-weighted sites drop out of the support before contract
admissibility runs. Same v3.0 constrain-before-propose
ordering established by the segment-rates slice (audit §4) and
extended here with one more multiplicative factor.

### Empty-support semantics

If `v_subregion_rates` zeros out all positions in every
assigned V allele's annotated subregions, AND segment rates do
not compensate by routing SHM elsewhere, AND `productive_only`
filters reject everything else, the pass may end up with
**empty admissible support** for a record. Behaviour follows
the existing mutation empty-support rules
([`engine_rs/src/passes/mutate/s5f/execution.rs:152`](../engine_rs/src/passes/mutate/s5f/execution.rs#L152)):

- **Strict mode** → `StrictSamplingError` at the record.
- **Permissive mode** → skip the mutation slot for that record,
  produce a record with fewer SHM than the count distribution
  asked for.

The slice introduces **no new** empty-support error variants.
`FilteredSampleError::EmptyAdmissibleSupport` already covers
the case.

### Composition with the existing rules

| Combination | Behaviour |
|---|---|
| `v_subregion_rates={"CDR1": 0}` + `productive_only()` | CDR1 sites drop out of SHM proposal. Other V sites still mutate. AIRR record validates clean. |
| `v_subregion_rates={"CDR1": 0, "CDR2": 0, "FWR1": 0, "FWR2": 0, "FWR3": 0}` + `segment_rates={"V": 0}` redundant case | V is fully zeroed by both factors; D/J/NP can still mutate. No empty-support error unless those are also zeroed. |
| `v_subregion_rates={"CDR1": 0, …, "FWR3": 0}` + `segment_rates={"D": 0, "J": 0, "NP": 0}` | Every segment zeroed → empty support per record → strict raises / permissive skips. |
| Annotated V + unannotated V in the same cartridge | Annotated V applies the factor; unannotated V is factor 1.0 (Q4 above). Mixed populations stable. |

### Pinned

- `pin_scaffold_productive_only_preserves_triad_under_shm` —
  same scaffold pin as the segment-rates slice; carried
  forward.
- `pin_scaffold_filtered_sample_error_variants_unchanged` —
  the slice must not add new empty-support error kinds.

---

## 6. Q6 — Trace / replay

### No new trace addresses

The slice's site-weighting machinery composes multiplicatively
into the existing `MutateS5fSite(u32)` /
`MutateUniformSite(u32)` addresses. The same constraint that
held for the segment-rates slice: rates change *which sites
have non-zero weight*, not *what address space is recorded*.

The existing replay path consumes the recorded site value
verbatim; under a flat-default `v_subregion_rates` the replay
is byte-identical.

### Plan-signature situation **[Closed by Slice A]**

The pre-implementation audit surfaced a finding:
`pass_plan_signature` joined pass *names* only, despite the
segment-rates audit's promise that compile-time pass
parameters would enter the signature. Today, replaying a trace
recorded under `segment_rates={"V": 2.0}` against a fresh plan
with `segment_rates={"V": 1.0}` would silently sample different
bases at the same recorded sites — no signature mismatch
fires.

**Status: shipped.** Slice A (Pass Parameter Signature) landed
the fix. The v_subregion_rates slice (Slice B) inherits the
fixed surface and does not need to touch the signature
mechanism.

What landed:

- **`Pass::parameter_signature(&self) -> String`** is a new
  trait method ([`engine_rs/src/pass/traits.rs`](../engine_rs/src/pass/traits.rs)),
  default `""`.
- **`pass_plan_signature`** now emits `name(params)` per pass.
  ([`engine_rs/src/trace_file.rs`](../engine_rs/src/trace_file.rs))
- Trace-file **schema bumped to v3**.
- Every parameterized pass implements `parameter_signature`
  via shared formatter helpers in
  [`engine_rs/src/passes/paramsig.rs`](../engine_rs/src/passes/paramsig.rs).
  Default rate vectors / probabilities collapse to the empty
  contribution so behaviourally-equivalent inputs collide.
- **S5F kernel name** threads through from
  `_MutateStep.s5f_model_name` → `push_mutate_s5f(kernel_name=…)`
  → `S5FMutationPass.kernel_name: Option<String>` →
  `parameter_signature`, so two passes that share rate +
  segment_rates but disagree on the loaded kernel produce
  different signatures.
- v1/v2 fixtures stay loadable + replayable via the legacy
  `pass_plan_signature_names_only` comparator inside the
  `validate_against` / `replay_from_trace_file` /
  `rerun_from_trace_file` paths (schema-version-aware
  comparator).

Pinned by
[`tests/test_pass_parameter_signature.py`](../tests/test_pass_parameter_signature.py)
(end-to-end: rate-vector divergence flips signature; mismatched
replay raises; default ≡ all-ones; kernel name participates;
schema is v3; distribution-bearing passes vary their signature
when distributions change).

**Implication for Slice B (v_subregion_rates):** the slice
inherits the fixed surface. Its `VSubregionRateWeights`
adds another factor inside `parameter_signature` for both
mutate passes; the existing `paramsig::fmt_segment_rates`
pattern is the template (default → empty short-circuit; otherwise
canonical key=value rendering). No new mechanism needed.

### Replay determinism under default rates

Flat-default `v_subregion_rates` (or `None` kwarg) reproduces
today's sampling distribution exactly. Existing traces replay
unchanged. **Backwards-compatible.**

### Pinned

- `pin_scaffold_pass_plan_signature_omits_compile_time_params_today` —
  freezes today's actual signature shape so the implementation
  slice flips this pin in lockstep with the fix.
- `pin_scaffold_shm_replay_byte_deterministic_today` — carry-
  over from the segment-rates audit; default-rate path must
  reproduce existing traces byte-identically.

---

## 7. Q7 — Manifest extension

### Today's state

[`cartridge_manifest()["models"]["shm"]`](../src/GenAIRR/dataconfig/data_config.py)
carries `segment_rate_support` (the per-segment capability
block) and `v_subregion_support` (the Slice-1 annotation
coverage block). It does NOT carry a `v_subregion_rate_support`
key.

### What the slice should add

```python
"v_subregion_rate_support": {
    "available": True,
    "labels": ["FWR1", "CDR1", "FWR2", "CDR2", "FWR3"],
    "aliases": {
        "FWR": ["FWR1", "FWR2", "FWR3"],
        "CDR": ["CDR1", "CDR2"],
    },
    "default": {"FWR1": 1.0, "CDR1": 1.0, "FWR2": 1.0,
                "CDR2": 1.0, "FWR3": 1.0},
    "annotated_v_count": <int>,       # mirror of v_subregion_support
    "requires_annotated_v": True,     # documents that the kwarg
                                       # needs the Slice-1 surface
    "in_content_hash": False,         # per-experiment, not cartridge
},
```

**`in_content_hash=False`** because the rate vector itself is a
per-experiment parameter (passed to `Experiment.mutate(...)`),
not a cartridge property. This mirrors how `segment_rate_support`
declares its rates outside the content hash.

### Pinned

- `pin_present_v_subregion_support_in_manifest` (carried from
  Slice 1 — the annotation surface).
- `pin_absence_no_v_subregion_rate_support_in_manifest` —
  the rate-support key doesn't exist yet.

---

## 8. Q8 — Counters (deferred)

The brief defers CDR/FWR mutation counters to a future slice:

> No CDR/FWR counters in v1. Future fields could be
> `n_fwr_mutations`, `n_cdr_mutations` or full five-label
> counters; audit should recommend later.

**Endorsed.** The rate slice should NOT introduce new AIRR
counter fields. The user-visible signal of subregion targeting
is the **rate vector itself** (inspectable, replayable, in the
manifest) plus the resulting biological distribution of
mutations across regions (which downstream consumers can compute
from the existing `n_v_mutations` + the `airr_record` + the
cartridge's V subregions).

### Recommendation for a later "subregion counters" slice

When users start asking for per-region SHM counts in AIRR
records, the design choice is between:

- **Two-bucket counters** — `n_cdr_mutations` (CDR1 + CDR2) and
  `n_fwr_mutations` (FWR1 + FWR2 + FWR3). Matches biology
  vocabulary; minimal field count.
- **Five-label counters** — `n_cdr1_mutations`,
  `n_cdr2_mutations`, `n_fwr1_mutations`, `n_fwr2_mutations`,
  `n_fwr3_mutations`. Matches the rate vector shape; lets
  users see CDR1 vs CDR2 separately (which the rate slice
  enables targeting separately).

**The audit recommends five-label counters** for symmetry with
the rate kwarg (one counter per rate bucket the user can
target), with a derived `n_cdr_mutations` / `n_fwr_mutations`
property on the AIRR record dict if the two-bucket vocabulary
turns out to be more frequently consumed in practice.

Either way, that's a separate slice with its own contract
file and audit — `v_region_substructure_audit.md` already pins
their absence as `pin_absence_no_cdr_fr_mutation_counter_fields`.
This audit carries that pin forward unchanged.

### Pinned

- `pin_absence_no_cdr_fr_mutation_counter_fields` — carried
  from the V-Substructure audit.

---

## 9. Edge cases the implementation slice must handle

| Case | Expected behaviour |
|---|---|
| `v_subregion_rates=None` | Flat-default fast path; byte-identical to pre-slice. |
| `v_subregion_rates={}` | Empty dict is equivalent to flat default. Validator accepts. |
| `v_subregion_rates={"FWR": 0.5}` | Expands to `{FWR1: 0.5, FWR2: 0.5, FWR3: 0.5, CDR1: 1.0, CDR2: 1.0}`. |
| `v_subregion_rates={"FWR": 0.5, "FWR2": 2.0}` | Alias expands; specific label overrides. Result: `{FWR1: 0.5, FWR2: 2.0, FWR3: 0.5, CDR1: 1.0, CDR2: 1.0}`. |
| `v_subregion_rates={"CDR3": 1.0}` | `ValueError` at the DSL boundary — CDR3 is out of scope for V annotation. |
| All five labels set to `0.0` + V is the only segment in the cartridge (VJ chain, V dominates) | Empty support per record → strict raises, permissive skips. |
| Cartridge with zero annotated V alleles + non-default `v_subregion_rates` | DSL-side compile-time `ValueError` — unsatisfiable. |
| Cartridge with mixed annotated/unannotated V alleles + non-default rates | Unannotated V → factor 1.0 silently. |
| A V site between `FWR3.end` and `len(allele.seq)` (V-side CDR3 contribution) | Factor 1.0 — outside annotation, not an error. |
| Trace replay with mismatched `v_subregion_rates` | See §6 above — depends on which Option (A / B) the slice adopts for plan-signature folding. |
| User authors a cartridge with `Allele.subregions` overlapping the assembled CDR3 region | The Slice-1 bridge validator (`parse_subregions`) rejects overlapping intervals; can't reach the SHM rate slice. |

---

## 10. Performance

### Per-position lookup cost

`build_profile` already calls `segment_at_position` once per
non-zero-mutability position when `segment_rates_active`. The
slice adds one more call to `v_subregion_at_position` for
positions where `segment == V`. Both are linear scans over ≤ 5
intervals — constant-bounded per call. Per-pool overhead is
`O(n_positions) × O(5)` = `O(n_positions)`, same shape as
existing segment_at_position cost.

### Fast-path preservation

Default `v_subregion_rates` (`is_default()` returns true) takes
the existing fast path — no per-position subregion lookup,
byte-identical to today. Confirmed by the segment-rates slice's
behavior with `SegmentRateWeights::default()`.

### No new allocations

The lookup returns an `Option<&VSubregion>` or
`Option<VSubregionLabel>` — borrows into the `Allele.subregions`
vec. No heap allocation per pool position.

### Pinned

- `pin_scaffold_subregion_lookup_is_constant_time` — see §2.

---

## 11. Implementation order

Two slices, in order of dependency. Slice A has shipped; Slice
B remains the next implementable step.

### Slice A — Pass Parameter Signature **[Shipped]**

The §6 finding about `pass_plan_signature` not folding
compile-time parameters has been closed. The implementation:

- Added `Pass::parameter_signature(&self) -> String` to the
  trait with default `""`
  ([`engine_rs/src/pass/traits.rs`](../engine_rs/src/pass/traits.rs)).
- Implemented on every parameterized pass: `UniformMutationPass`,
  `S5FMutationPass` (folds segment_rates + count_source +
  kernel name), `TrimPass`, `GenerateNPPass`,
  `PairedEndSamplingPass`, `PCRErrorPass`, `QualityErrorPass`,
  `IndelPass`, `NCorruptionPass`, `EndLossPass`,
  `ContaminantPass`, `RevCompPass`, `InvertDPass`,
  `ReceptorRevisionPass`.
- Threaded S5F kernel name through PyO3
  (`push_mutate_s5f(kernel_name=…)`,
  `push_mutate_s5f_rate(kernel_name=…)`) and the Python
  compiler (`_lower_mutate` passes `step.s5f_model_name`).
- Bumped `TRACE_FILE_SCHEMA_VERSION` from 2 to 3.
- Added `pass_plan_signature_names_only` as the legacy
  v1/v2 comparator; schema-version-aware comparator in
  `validate_against`, `replay_from_trace_file`, and
  `rerun_from_trace_file`.
- v1 golden fixtures preserved as backwards-compat anchors;
  byte-equality tests gated on `loaded.schema_version ==
  live.schema_version` to handle cross-schema replay
  gracefully.
- Pinned by
  [`tests/test_pass_parameter_signature.py`](../tests/test_pass_parameter_signature.py)
  (8 spec tests).

### Slice B — V-subregion SHM rate kwarg **[Shipped]**

What landed:

1. **Python DSL surface** ([src/GenAIRR/experiment.py](../src/GenAIRR/experiment.py)):
   - `Experiment.mutate(..., v_subregion_rates: Optional[Dict[str, float]] = None)`.
   - `_validate_v_subregion_rates(...)` enforces the five
     canonical labels (`FWR1` / `CDR1` / `FWR2` / `CDR2` /
     `FWR3`) plus the two aliases (`FWR` → FWR1/FWR2/FWR3,
     `CDR` → CDR1/CDR2). Alias-then-override expansion:
     `{"FWR": 0.5, "FWR2": 2.0}` resolves to `FWR1=0.5,
     FWR2=2.0, FWR3=0.5, CDR1=1.0, CDR2=1.0`. Rejects
     unknown keys, `bool`, `NaN`, `inf`, negative, and
     all-zero-after-expansion at the DSL boundary.
   - `Experiment._check_v_subregion_rates_satisfiable(...)`
     fires at `mutate()` call time and rejects a non-default
     rate vector against a cartridge with zero annotated V
     alleles. Default rates always pass (no-op).
   - `_MutateStep.v_subregion_rates` carries the
     5-tuple in canonical order (FWR1, CDR1, FWR2, CDR2,
     FWR3) for stable PyO3 plumbing.
   - `_lower_mutate` in [_compile.py](../src/GenAIRR/_compile.py)
     forwards the tuple through all four `push_mutate_*` paths
     (uniform + S5F × rate + count).

2. **Rust** ([engine_rs/src/passes/mutate/](../engine_rs/src/passes/mutate/)):
   - `VSubregionRateWeights { fwr1, cdr1, fwr2, cdr2, fwr3 }`
     in
     [v_subregion_rates.rs](../engine_rs/src/passes/mutate/v_subregion_rates.rs)
     with `default()`, `from_tuple(...)`, `is_default()`, and
     `rate_for(VSubregionLabel)` — mirror of
     `SegmentRateWeights`.
   - `v_subregion_at_position(sequence, assignments, refdata,
     pos) -> Option<VSubregionLabel>` sibling of
     `segment_at_position`. Implements the audit §2 arithmetic
     `allele_local = pos - v_region.start + instance.trim_5`
     and walks `Allele.subregions`.
   - `S5FMutationPass` and `UniformMutationPass` each carry
     `v_subregion_rates: VSubregionRateWeights` plus a
     `with_v_subregion_rates(...)` setter.
   - A unified `combined_site_factor(...)` helper in
     [mutation_transaction/substitution.rs](../engine_rs/src/passes/mutation_transaction/substitution.rs)
     computes `segment_rate × v_subregion_rate` per position;
     short-circuits to `1.0` when both vectors are flat
     default. Both mutation passes call it from a single site
     in `build_profile` / `substitute_position_constrained`
     instead of separate per-vector code paths.
   - Replay validation paths (both S5F and uniform) reject
     recorded sites that the current `combined_site_factor` →
     `0.0` would have excluded from support.

3. **Plan signature** (Slice A surface):
   - `paramsig::fmt_v_subregion_rates(...)` formats the rate
     vector as `v_subregion_rates=(fwr1:…,cdr1:…,fwr2:…,cdr2:…,fwr3:…)`.
   - Default rates short-circuit to the empty string, so
     `v_subregion_rates=None` and an explicit all-ones dict
     collide on equal signatures (and equal sampling
     distributions).
   - Both `S5FMutationPass::parameter_signature` and
     `UniformMutationPass::parameter_signature` fold it in
     alongside `segment_rates` / kernel name / count source.
   - Trace replay with mismatched `v_subregion_rates` fails
     at the plan-signature gate before consuming any choice.

4. **PyO3 wire format** ([engine_rs/src/python/plan.rs](../engine_rs/src/python/plan.rs)):
   - All four `push_mutate_*` (uniform + S5F × rate + count)
     accept `v_subregion_rates: Option<Vec<f64>>`.
   - `build_v_subregion_rate_weights(...)` is the
     defence-in-depth re-validator (mirror of
     `build_segment_rate_weights`) for callers that bypass
     the DSL.

5. **Manifest** ([src/GenAIRR/dataconfig/data_config.py](../src/GenAIRR/dataconfig/data_config.py)):
   - `cartridge_manifest()["models"]["shm"]["v_subregion_rate_support"]`
     carries `available`, `labels`, `aliases`, `default`,
     `in_plan_signature=True`, `in_content_hash=False`.

6. **Tests** (~50 specs total):
   - [tests/test_v_subregion_rates_implementation.py](../tests/test_v_subregion_rates_implementation.py)
     (27 spec tests): DSL validation (label set, aliases,
     expansion order, bool/NaN/inf/negative/all-zero,
     unsatisfiable-rates rejection on an
     annotation-stripped cartridge), alias expansion + override,
     FWR-only / CDR-only configurations validate clean,
     zero-rate CDR drops CDR mutations (with baseline-diff
     classifier), zero-rate FWR drops FWR mutations, non-V
     sites unaffected, both uniform and S5F models work,
     productive_only triad invariants preserved, matching-rate
     replay succeeds, mismatched-rate replay fails at the
     plan-signature gate, default ≡ explicit all-ones,
     manifest advertises support, CDR/FR counter fields stay
     absent.
   - [tests/test_v_subregion_shm_rate_contract.py](../tests/test_v_subregion_shm_rate_contract.py)
     pins flipped: kwarg present, validator present, Rust
     struct present, helper present, manifest support
     present, compile-time rejection present.

What this slice deliberately did NOT do (still deferred):

- **CDR/FR mutation counter AIRR fields** —
  `n_cdr1_mutations` / `n_cdr2_mutations` / `n_fwr1_mutations`
  / `n_fwr2_mutations` / `n_fwr3_mutations` (or any
  two-bucket aggregate like `n_cdr_mutations` /
  `n_fwr_mutations`). These are a separate future counters
  slice; the audit's §8 recommendation stands.

### Why this order

Slice A was independent and closed an existing gap (Slice A's
fix makes the segment_rates AND v_subregion_rates kwargs both
fail-fast on replay-vs-record mismatch). Slice B builds the
user-facing biology surface on top.

---

## 12. Test surface — what this audit pins

Mirrored in
[`tests/test_v_subregion_shm_rate_contract.py`](../tests/test_v_subregion_shm_rate_contract.py).

### `pin_scaffold_*` — pre-existing surfaces Slice B builds on

1. `segment_at_position` takes `(sequence, pos)` only — no
   allele identity.
2. `Simulation.assignments` carries `AlleleInstance` per
   segment, including `allele_id` + `trim_5`.
3. The bundled human IGH / IGK / IGL cartridges have 100% V-
   subregion coverage (Slice 1 prerequisite for the rate
   slice).
4. Each V allele in bundled cartridges has the five canonical
   labels (FWR1 / CDR1 / FWR2 / CDR2 / FWR3).
5. Existing `segment_rates` still works against subregion-
   annotated cartridges (no regression from Slice 1).
6. `Experiment.mutate(segment_rates={"V": 2.0})` against
   `HUMAN_IGH_OGRDB` produces a clean simulation that
   validates.
7. The mutate trace addresses today are
   `MutateS5fSite(u32)` / `MutateS5fBase(u32)` /
   `MutateS5fCount` and the uniform equivalents — Slice B
   does not introduce any new addresses.

### `pin_present_*` — Slice A / Slice B landed surfaces

8. **Pass-plan signature folds compile-time pass parameters**
   (Slice A — `Pass::parameter_signature`, `name(params)`
   envelope, schema v3 with v1/v2 fallback).
9. **`S5F build_profile` applies the unified
   `combined_site_factor`** (segment × V-subregion) — single
   call site for both rate vectors.
10. **`Experiment.mutate(v_subregion_rates=...)` kwarg
    exists** with default `None`.
11. **`_validate_v_subregion_rates`** exists in
    `experiment.py` with the canonical labels, `FWR` / `CDR`
    alias constants, and the alias-expansion + override
    rule.
12. **`VSubregionRateWeights` struct exists** in the Rust
    crate (sibling of `SegmentRateWeights`).
13. **`v_subregion_at_position` helper exists** in the Rust
    crate (sibling of `segment_at_position`).
14. **`v_subregion_rate_support` block exists** in the
    cartridge manifest with the six documented keys.
15. **Compile-time rejection fires** when non-default rates
    are supplied against a cartridge with zero annotated V
    alleles.

### `pin_absence_*` — gap reserved for a future counters slice

16. No CDR/FR mutation counter AIRR fields
    (`n_cdr1_mutations`, `n_cdr2_mutations`,
    `n_fwr1_mutations`, `n_fwr2_mutations`,
    `n_fwr3_mutations`, `n_cdr_mutations`, `n_fwr_mutations`).
    Pinned absent so the next biology slice (subregion
    counters) flips this in lockstep.

### Doc anchor

17. The audit doc references the contract file; structure
    intact.

---

## 13. Out of scope

Documented here so a future implementer doesn't accidentally
expand the work.

- **CDR/FWR mutation counter AIRR fields.** Deferred per §8.
- **CDR3 / FWR4 targeting.** The V annotation surface stops at
  FWR3. CDR3 lives in the junction (covered by `segment_rates`
  bucket `NP`); FWR4 is in the J segment.
- **D-region subregions.** D alleles have no IMGT-canonical
  substructure of analogous granularity. Out of scope.
- **C-region subregions.** The engine doesn't model C-region
  SHM passes. Out of scope.
- **Per-allele rate overrides.** A future slice could let the
  user say "double SHM on this specific allele" via
  `allele_rates={"IGHV1-2*02": 2.0}`. Composes
  multiplicatively with subregion rates; not in this slice.
- **Empirical S5F-subregion kernel.** S5F's 5-mer kernel is
  region-agnostic; a region-aware kernel (e.g. fitting separate
  CDR vs FWR 5-mer tables) would be a separate biology slice
  well past this audit.
- **AID hotspot motif targeting.** Same answer as in the
  V-Substructure audit: out of scope.
- **TCR support.** SHM doesn't occur in TCRs (already enforced
  by `Experiment.mutate`'s TCR guard); subregion rates inherit
  the same enforcement.
- **Replay with `strict_subregions=True`.** That's the future
  strict-mode ergonomics slice from §4; not in v1.

---

## 14. Summary table

| Concern | Post-Slice-B state | Status |
|---|---|---|
| Rate model shape | `final = base × segment_rate × v_subregion_rate` (V only); non-V sites unaffected | **Shipped** |
| DSL kwarg | `v_subregion_rates={...}` parallel to `segment_rates`; sparse keys default to 1.0 | **Shipped** |
| Canonical labels | Five IMGT labels + `"FWR"` / `"CDR"` aliases with alias-then-override semantics | **Shipped** |
| Default behaviour | `None`, `{}`, and explicit all-ones all short-circuit to the pre-slice fast path | **Shipped** |
| Zero-rate semantics | Sites drop out of proposal support before contract admissibility | **Shipped** |
| Productive-only interaction | Empty-support → strict raises, permissive skips; no new error variants | **Shipped** |
| Unsatisfiable rates against unannotated cartridge | DSL-side `ValueError` at `mutate()` time | **Shipped** |
| Unannotated V allele in mixed cartridge | Silent fallback to factor `1.0` | **Shipped** |
| Trace addresses | Three per model; no new variants | **Shipped** (no change) |
| Plan signature includes rates | `pass_plan_signature` folds `v_subregion_rates` (and `segment_rates`) via Slice A | **Shipped** |
| Manifest reports rate capability | `v_subregion_rate_support` block alongside `segment_rate_support` | **Shipped** |
| Pool → V-subregion mapping helper | `v_subregion_at_position` sibling of `segment_at_position`; `combined_site_factor` unifies the two factors | **Shipped** |
| Mutation counters per region | None | Deferred — separate counters slice |

The v_subregion_rates feature shipped clean. The plan-signature
gap surfaced by the audit (`pass_plan_signature` joining only
pass names) was closed first by Slice A, and Slice B inherits
the fixed surface — `v_subregion_rates` divergence between a
recorded trace and a replayed plan now fails at the signature
gate before any choice is consumed.

The natural next biology slice is **subregion mutation
counters** — five `n_<region>_mutations` AIRR fields per the
audit's §8 recommendation (or the two-bucket
`n_cdr_mutations` / `n_fwr_mutations` aggregate). That slice
will land on top of the already-pinned annotation surface
(Slice 1) + rate-targeting surface (Slice B) and continues to
be pinned absent in the contract file.
