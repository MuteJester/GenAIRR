# Per-Segment SHM Rate Scalars — Architecture Contract

**Status: shipped.** The implementation slice landed with the
shape this audit specified: `Experiment.mutate(...,
segment_rates: Optional[Dict[str, float]] = None)` accepting the
canonical four-bucket V / D / J / NP dict, threaded through
`_MutateStep` → PyO3 → both Rust `UniformMutationPass` /
`S5FMutationPass` passes via a `SegmentRateWeights` value type.
The constrain-before-propose ordering held; the cartridge
manifest's `models.shm.segment_rate_support` block advertises
the capability with `in_content_hash=False` documenting the v1
boundary.

Release-tier consolidation closed the slice the same way D
inversion / receptor revision / paired-end closed: full-stack
IGH + targeted-SHM + replay round-trip + zero-rate exclusion
invariant in [`test_release_validation.py`](../tests/test_release_validation.py),
plus a row in [`docs/validation_matrix.md`](validation_matrix.md)
naming the audit doc / contract file / implementation file /
Rust kernel surfaces.

---

The original pre-implementation audit follows verbatim below.
Pins today's full-pool SHM targeting (pre-slice baseline), audits
the proposed per-segment-rate biology extension, and recommends
the slice's exact shape — the shape that ultimately shipped.

This audit closes the modeling gap surfaced by the SHM Model
Audit (`docs/shm_model_audit.md` §4 Q4): uniform and S5F SHM
currently target the full assembled pool, treating V / NP / D / J
as a flat substrate. Real B-cell SHM concentrates on V (with
internal rate spikes at CDR1/CDR2) and some on J; NP regions
should be largely SHM-quiet. This audit specifies the
*architecture contract* for adding segment-level rate scalars
without breaking replay determinism, productive-only semantics,
or the existing trace format.

Companion to
[`tests/test_shm_segment_rate_contract.py`](../tests/test_shm_segment_rate_contract.py)
which freezes every claim — today's full-pool targeting,
absence of `segment_rates` DSL surface, the no-per-segment-AIRR-
counter status, replay determinism — as either `pin_scaffold_*`
or `pin_absence_*`.

The pre-flight check confirms the implementation is feasible:
the site→segment lookup is cleanly determinable from
`sim.sequence.regions` (a `Vec<Region>` carrying
`{segment, start, end}`). No bug to report.

---

## Existing scaffolding the contract pins

| Surface | Where | What it pins |
|---|---|---|
| `Region` struct | [`engine_rs/src/ir/region.rs:18`](../engine_rs/src/ir/region.rs#L18) | Per-region IR record carrying `segment: Segment`, `start: NucHandle`, `end: NucHandle`. Site→segment lookup keys off this. |
| `Sequence.regions: Vec<Region>` | [`engine_rs/src/ir/sequence.rs:8`](../engine_rs/src/ir/sequence.rs#L8) | Biological assembly order: V → NP1 → D → NP2 → J for VDJ; V → NP1 → J for VJ. |
| `Segment` enum | [`engine_rs/src/ir/segment.rs:13`](../engine_rs/src/ir/segment.rs#L13) | Five variants: `V` / `Np1` / `D` / `Np2` / `J`. Discriminant indices stable so `PerSegment` indexing is O(1). |
| `UniformMutationPass` site sampling | [`engine_rs/src/passes/mutate/uniform.rs:115-128`](../engine_rs/src/passes/mutate/uniform.rs#L115-L128) | Uniform across `[0, pool_len)`. No segment awareness today. |
| `S5FMutationPass` site sampling | [`engine_rs/src/passes/mutate/s5f/execution.rs:282-313`](../engine_rs/src/passes/mutate/s5f/execution.rs#L282-L313) | Weighted by per-position mutability ratio. No segment awareness today. |
| `MutationTransaction::substitute_position_constrained` | [`engine_rs/src/passes/mutation_transaction/substitution.rs`](../engine_rs/src/passes/mutation_transaction/substitution.rs) | Constrain-before-propose primitive. Filters candidate bases against `ContractSet` at each site; the proposed slice plugs into the SAME primitive with a segment-rate mask. |
| Trace addresses | [`engine_rs/src/address.rs`](../engine_rs/src/address.rs) | `mutate.{uniform,s5f}.{count,site[i],base[i]}`. Segment rates do NOT introduce new addresses (compile-time parameter, not a sampled value). |
| `n_mutations` AIRR field | [`engine_rs/src/airr_record/builder.rs:211`](../engine_rs/src/airr_record/builder.rs#L211) | IR-sourced; remains global under the proposed slice. |
| `DataConfig.cartridge_manifest()["models"]["shm"]` | [`src/GenAIRR/dataconfig/data_config.py`](../src/GenAIRR/dataconfig/data_config.py) | Today carries `available_models`, `s5f_kernels_available`, `default_s5f_kernel`, `s5f_kernel_digest`, `in_content_hash`. A future field for `segment_rate_support` would extend this. |

---

## 1. Q1 — Current targeting

### Confirmed: SHM targets the full assembled pool

`UniformMutationPass.execute_with_sampling_mode`
([`uniform.rs:81-139`](../engine_rs/src/passes/mutate/uniform.rs#L81-L139)):

```rust
let pool_len = sim.pool.len() as u32;
let count = sample_validated_count(&self.count_source, ctx, pool_len, ...);
...
for i in 0..count {
    let wrote = tx.substitute_position_constrained(
        self.base_dist.as_ref(),
        address::ChoiceAddress::MutateUniformSite(i),
        address::ChoiceAddress::MutateUniformBase(i),
        None,
    )?;
    ...
}
```

The site range is `[0, pool_len)` — the entire assembled pool,
spanning V + NP1 + D + NP2 + J indiscriminately.

`S5FMutationPass` builds a per-position mutability profile of
length `pool_len` and samples by cumulative weight. Same range,
just weighted instead of flat.

### Behavioural confirmation: mutations land in every region

Audit probe (canonical IGH `recombine().mutate(count=50)`, 3 records):

| Record | V | NP1 | D | NP2 | J |
|---|---|---|---|---|---|
| 0 | 32 | 1 | 3 | 2 | 8 |
| 1 | 34 | 0 | 5 | 0 | 8 |
| 2 | 36 | 1 | 4 | 1 | 3 |

V dominates by length (~78% of the pool) but NP1 / D / NP2 / J
each receive proportional hits. This pins today's "flat substrate"
behaviour as the audit baseline.

### Pre-flight: site→segment lookup is cleanly determinable

Given a `NucHandle` position `pos` and `sim.sequence.regions`:

```rust
sim.sequence.regions.iter()
    .find(|r| r.start.index() <= pos && pos < r.end.index())
    .map(|r| r.segment)
```

Linear scan over ~5 regions per VDJ pool; effectively O(1).
Alternative: precompute `Vec<Segment>` of length `pool_len` once
per pass execution (O(pool_len), ~400 bytes per record). Both
viable; choice is an implementation detail for the slice.

**No architectural blocker.** The audit greenlight is unblocked
on the lookup question.

Pinned by:
- `pin_scaffold_shm_targets_full_pool_today` — behavioural
  confirmation of per-region distribution under high count.
- `pin_scaffold_site_segment_lookup_is_cleanly_determinable` —
  source-level pin on `Region.segment` field's presence.

---

## 2. Q2 — Desired model: per-segment rates

### Recommendation: segment-family rates, four buckets

The user's suggested shape — **endorsed**:

```python
segment_rates = {"V": 1.0, "D": 1.0, "J": 1.0, "NP": 1.0}
```

with `0.0` allowed to disable a region class.

### Why not per-region (NP1 / NP2 separate)?

The audit explored three alternatives:

| Shape | Buckets | Pro | Con |
|---|---|---|---|
| **Segment-family** (recommended) | 4: V / D / J / NP | Matches biological literature (SHM rates are usually quoted per V/D/J + "non-templated"). Compact DSL. | Can't independently tune NP1 vs NP2. |
| Per-segment | 5: V / Np1 / D / Np2 / J | Matches the engine's `Segment` enum 1:1. Lowest impedance. | Two NP rates are redundant (the literature reports them as a class). |
| Per-region | Variable | Lets future CDR1/CDR2/FR1/FR2 SHM hotspots ride this surface. | Massive DSL surface; needs per-cartridge region naming convention. |

**Recommendation: segment-family.** The bucket count matches
biology (4 is the canonical SHM-rate breakdown in immunology
papers); the NP1/NP2 split is a cartridge-internal detail that
biologists don't tune separately. A future "intra-V hotspot
rates" slice can ride a separate `cdr_rates` parameter without
touching `segment_rates`.

### Internal mapping

`segment_rates` keys map to `Segment` variants as:

| DSL key | `Segment` variants |
|---|---|
| `"V"` | `Segment::V` |
| `"D"` | `Segment::D` |
| `"J"` | `Segment::J` |
| `"NP"` | `Segment::Np1`, `Segment::Np2` |

### Default behaviour

Omitting `segment_rates` from the DSL call produces today's
behaviour: all segments mutate at rate 1.0 (the "flat substrate"
baseline). The slice is purely additive — no existing pipeline
changes meaning.

### Zero-rate semantics

`segment_rates={"NP": 0.0}` excludes NP1 and NP2 positions from
the support entirely. The realised mutation count reflects only
the sites the engine could pick from V / D / J. This is the
biologically-correct behaviour for "SHM doesn't hit NP" — the
pass doesn't silently re-target the rejected mass to other
segments.

### Pinned

- `pin_scaffold_no_segment_rates_kwarg_yet` — current `mutate()`
  signature lacks `segment_rates`.
- `pin_scaffold_segment_enum_has_five_variants` — `V / Np1 / D /
  Np2 / J`, the targets the future surface maps onto.

---

## 3. Q3 — S5F interaction

The proposed weighting rules:

### S5F + segment rates

```
final_weight(site) = S5F_context_mutability(site) * segment_rate(segment_at(site))
```

The S5F kernel's per-position mutability is **multiplied** by the
segment rate. Zero-rate sites drop to zero weight → excluded from
support before the cumulative-weight sampling step.

This is the natural extension of the existing S5F formula:
today's `final_weight = S5F_context_mutability` becomes
`final_weight = S5F_context_mutability * 1.0` under default
rates (no behaviour change), and any non-default rate scales the
per-segment mass.

### Uniform + segment rates

Today's uniform model picks a position uniformly across the pool
(`final_weight(site) = 1`). Under segment rates:

```
final_weight(site) = segment_rate(segment_at(site))
```

— a piecewise-constant distribution: positions within a segment
share the same weight; segment-level mass scales with the rate.

Default rates of 1.0 across the four buckets reproduce today's
uniform behaviour exactly.

### Zero-rate exclusion

If a segment's rate is `0.0`, sites in that segment are
**dropped from the support** before normalisation. The
cumulative-weight sampler never visits them. This is the same
shape as `substitute_position_constrained` handles
contract-rejected sites today — the support shrinks and the
sampler operates on the surviving mass.

### Edge case: all-zero rates

If `segment_rates = {"V": 0, "D": 0, "J": 0, "NP": 0}`, the
support is empty across every site. This is the same shape as a
universally contract-rejected proposal: permissive mode silently
skips (realised count = 0); strict mode raises a structured
`ConstraintSampling` error. The slice does NOT introduce a new
empty-support sentinel.

Pinned by `pin_scaffold_segment_rate_zero_excludes_sites` —
specification-level pin on the proposed semantics (the slice
itself will land the behaviour).

---

## 4. Q4 — Productive-only interaction

The slice must compose with the existing
constrain-before-propose primitive. The user's spec dictates:

> segment-rate filtering happens **before** contract
> admissibility.

This is the correct ordering. The site-weight computation
produces a mass distribution over positions; the constraint
filter then evaluates per-base admissibility at the chosen
site. Reversing the order would compute contract-admissibility
masks at zero-rate sites only to discard them — wasted work.

### Concrete flow

For each mutation slot `i`:

1. Build site-weight vector: `w[site] = base_weight(site) *
   segment_rate(segment_at(site))` for `site in [0, pool_len)`.
2. Drop sites where `w[site] == 0`.
3. Sample a candidate site via cumulative-weight draw.
4. At the chosen site, ask the contract bundle for the
   admissible base mask.
5. If the mask is non-empty, sample a base from it and write.
6. If the mask is empty, permissive mode silently skips this
   slot; strict mode raises `ConstraintSampling`.

### All-rejected-sites edge case

If every nonzero-rate site is contract-rejected (e.g.
`productive_only()` + `segment_rates={"V": 0, "J": 0}` on a
cartridge where D/NP can't introduce a productive substitution
without breaking the frame), the behaviour matches today's
empty-support rules: permissive skips the slot, strict raises.
No new sentinel.

### Pinned

`pin_scaffold_productive_only_preserves_triad_under_segment_rates` —
spec-level pin that the slice MUST preserve the productive
triad under any non-degenerate rate combination. The
implementation tests will exercise concrete rate vectors.

---

## 5. Q5 — Trace / replay

The user's spec:

> no new trace addresses if rates are compile-time pass
> parameters.

**Endorsed.** Rates flow into the pass at compile time as part
of the `_MutateStep` dataclass; the pass reads them once at
construction and uses them for site weighting. No per-site
trace record is needed — the recorded `mutate.{model}.site[i]`
addresses still carry the chosen position; the rate vector
shapes the distribution at sample time, not the trace.

### Replay semantics

Existing replay path: `replay_from_trace_file` consumes the
recorded site + base values verbatim. Under the slice:

- Replay against a plan built with the **same** segment rates:
  the recorded sites land within their original support; the
  pass writes the same bases. **Byte-identical to the original
  run.**
- Replay against a plan built with **different** segment rates:
  the recorded sites might NOT all lie within the new support
  (e.g. a site recorded in NP under `np_rate=1.0` would be
  outside the support under `np_rate=0.0`). The validator
  should treat this as a structural mismatch — a
  `Replay::SiteOutsideSupport` error or equivalent — rather
  than silently accepting.

### Plan-signature implication

**Status: closed by Slice A (Pass Parameter Signature) —
[`v_subregion_shm_rate_design.md`](v_subregion_shm_rate_design.md)
§6 / §11.** The segment-rates audit originally promised that
`segment_rates` would enter the plan signature so a replay
against a mismatched rate vector fails the signature check
before consuming trace records. That part of the slice did not
actually ship at the time — `pass_plan_signature` continued to
join pass *names* only. The v_subregion_rates pre-implementation
audit re-discovered the gap and Slice A closed it.

What landed:

- `Pass::parameter_signature(&self) -> String` is a new trait
  method (default `""`); each parameterized pass implements it
  via the shared formatter helpers in
  `engine_rs/src/passes/paramsig.rs`.
- `pass_plan_signature` now emits `name(params)` per pass and
  joins with `|`. Trace-file schema bumped to v3.
- v1/v2 fixtures stay loadable + replayable via the legacy
  `pass_plan_signature_names_only` comparator inside
  `validate_against` / `replay_from_trace_file` /
  `rerun_from_trace_file`.
- `UniformMutationPass::parameter_signature` and
  `S5FMutationPass::parameter_signature` include the
  segment-rate vector via
  `paramsig::fmt_segment_rates(&self.segment_rates)`. Default
  rates short-circuit to the empty string so
  `segment_rates=None`, `segment_rates={}`, and explicit
  all-ones all collide on the same signature — behavioural
  equivalence preserved at the signature layer.

Pinned at the surface by
[`tests/test_pass_parameter_signature.py`](../tests/test_pass_parameter_signature.py)
and inside the v_subregion_shm_rate contract
([`tests/test_v_subregion_shm_rate_contract.py`](../tests/test_v_subregion_shm_rate_contract.py)
`test_pin_present_pass_plan_signature_folds_compile_time_params`).

### Replay determinism under default rates

Default rates `{V: 1, D: 1, J: 1, NP: 1}` reproduce today's
uniform sampling distribution exactly. Existing trace files
recorded under the current (pre-slice) `mutate()` API replay
correctly against the slice's `mutate()` with default rates.
**Backwards-compatible.**

### Pinned

- `pin_scaffold_shm_replay_byte_deterministic_today` — current
  same-seed replay produces byte-identical sequences. Slice
  must preserve under default rates.

---

## 6. Q6 — DSL shape

The user's preferred shape:

```python
.mutate(model="s5f", rate=0.03, segment_rates={"V": 1.0, "D": 0.2, "J": 0.5, "NP": 0.0})
```

The audit explored three alternatives:

### Option A — Dict-keyed `segment_rates` (recommended)

```python
.mutate(model="s5f", rate=0.03, segment_rates={"V": 1.0, "D": 0.2, "J": 0.5, "NP": 0.0})
```

**Pros:**
- Matches the cartridge-manifest's JSON-clean shape (dict of
  string → float).
- Easy to round-trip through configuration files / CLI args.
- Sparse: `segment_rates={"NP": 0.0}` defaults the other three
  to 1.0 implicitly.
- Single kwarg — doesn't bloat the `mutate()` signature.

**Cons:**
- Spelling errors in the key names (`"v"` vs `"V"`, `"np"` vs
  `"NP"`) need careful validation at the DSL boundary.

### Option B — Explicit kwargs

```python
.mutate(model="s5f", rate=0.03, v_rate=1.0, d_rate=0.2, j_rate=0.5, np_rate=0.0)
```

**Pros:**
- Pythonic — kwargs auto-validate on misspelling.
- IDE autocomplete shows the names.

**Cons:**
- Bloats the `mutate()` signature with four new kwargs.
- Doesn't survive JSON round-trip cleanly (would need to
  re-pack into a dict for config files).

### Option C — Hybrid (dict + per-key kwarg shortcut)

```python
.mutate(model="s5f", segment_rates={"V": 1.0, "NP": 0.0})  # dict form
.mutate(model="s5f", v_rate=1.0, np_rate=0.0)              # kwarg form
```

**Pros:**
- Accommodates both styles.

**Cons:**
- Two surfaces to maintain; mutual-exclusion logic at the
  boundary.

### Recommendation: Option A — dict-keyed `segment_rates`

The dict form composes cleanly with config-driven workflows and
matches the cartridge-manifest convention. Spelling errors
should be caught at the DSL boundary with a clear error:

```python
ValueError(
    "segment_rates: unknown segment key 'NP1'. "
    "Allowed: 'V', 'D', 'J', 'NP'."
)
```

Also at the boundary:
- All values must be `float | int`, non-negative, finite.
- At least one key must have a positive value (otherwise
  `mutate()` is a deterministic no-op; reject as
  `ValueError("segment_rates: all rates are zero; "
  "mutation pass would be a no-op")`).

### Pinned

`pin_scaffold_mutate_signature_lacks_segment_rates_today` —
exact source-level signature check on `Experiment.mutate(...)`
showing no `segment_rates` kwarg.

---

## 7. Q7 — AIRR / counters

The user's preferred shape:

> `n_mutations` stays global. No per-segment counters in v1.

**Endorsed.** Reasoning:

- The slice's biology change is at the *targeting* level
  (which sites get picked), not at the *counter* level. The
  realised count is still a single number per record.
- A future "per-segment counter" surface would be useful for
  comparing two cartridges' SHM profiles, but it's
  observational not enforcement — better to add it as a
  separate diagnostic slice once the targeting biology has
  stabilised.
- Adding per-segment counters in v1 doubles the AIRR record
  surface footprint without a clear consumer. Defer.

### Future shape (NOT in v1)

```
record["n_v_mutations"]
record["n_d_mutations"]
record["n_j_mutations"]
record["n_np_mutations"]
```

These would be derived at projection time by walking
`outcome.events()` and bucketing `BaseChanged` events by their
position's region. No new IR field needed; the events ledger
already carries the per-pass site indices.

### Pinned

- `pin_absence_no_per_segment_mutation_counters_today` — AIRR
  records carry only `n_mutations`.

---

## 8. Edge cases the implementation slice must handle

1. **VJ chains (no D / NP2).** Segment-rate dict carrying `"D"`
   on a VJ cartridge: the slice should ignore `"D"` (the pool
   has no D region; segment-at-site never returns D). No error.
   Same for `"NP"` keys mapping to Np1 only.

2. **Empty rate dict.** `segment_rates={}` treated identically
   to omitting the kwarg: all rates default to 1.0.

3. **Negative or NaN rates.** Reject at DSL boundary with a
   clear `ValueError`.

4. **Per-cartridge segment unavailability.** A future
   constant-region (C) addition to the pool would need the
   segment_rates dict to either silently ignore unknown keys or
   require explicit opt-in. Recommendation: silently ignore
   unknown keys (forward-compat with future segment additions).
   Pin this decision in the slice's contract.

5. **Receptor revision + segment rates.** Receptor revision
   replaces the V slice. SHM after receptor revision still
   operates on the post-revision pool. Segment rates are
   evaluated against the post-revision region layout —
   correct by construction (the pass reads `sim.sequence.regions`
   each execution).

6. **Indels + segment rates.** End-loss / polymerase indels
   shift region boundaries. SHM after indels still gets the
   correct segment-at-site lookup from the updated regions.
   Correct by construction.

7. **Clonal pipelines.** The slice composes with `expand_clones`
   the same way today's `mutate` does — SHM is descendant-phase
   (audit pinned), runs per descendant with its own segment-rate
   sampling. Parent IR is unaffected. No special handling.

---

## 9. Performance

### Site-weight computation cost

Today's S5F builds a per-position mutability profile in
O(pool_len). The slice adds a per-position segment lookup +
multiply:

- Linear scan over regions: O(pool_len × log n_regions) or
  O(pool_len × n_regions) — with ~5 regions, effectively
  O(pool_len).
- Precomputed `Vec<Segment>`: O(pool_len) build + O(1) per
  lookup.

Either approach is a constant-factor overhead to existing S5F
cost. Negligible.

### Default-rate fast path

When all rates are 1.0 (the default), the multiplier is
constant 1.0 and the profile is byte-identical to today's. A
constant-time check at pass construction can skip the per-
position multiply entirely. Recommendation: implement the fast
path so default-rate pipelines see zero overhead.

### Pinned

`pin_scaffold_performance_baseline_exists` — the release-tier
budget file is still present (no slice-specific perf pin
proposed by this audit; the implementation slice should add a
"default rates ≤ today's cost" pin).

---

## 10. Manifest integration

The cartridge manifest's `models.shm` block surfaces what the
engine supports, not what an experiment used. The proposed
slice would extend it with:

```python
manifest["models"]["shm"]["segment_rate_support"] = {
    "available": True,
    "buckets": ["V", "D", "J", "NP"],
}
```

(Or `False` before the slice ships.)

This documents the slice's capability without coupling the
manifest to any specific rate vector — the rate vector is
per-experiment, not per-cartridge.

Pinned by `pin_absence_no_segment_rate_support_in_manifest`.

---

## 11. Implementation order (recommended)

After the audit, one focused biology slice:

### Slice 1 — Per-segment SHM rate scalars

Scope:

1. Add `segment_rates: Optional[Dict[str, float]] = None` kwarg
   to `Experiment.mutate(...)`.
2. DSL-boundary validation: keys in `{"V", "D", "J", "NP"}`,
   values non-negative + finite, at least one positive value.
3. Plumb the rate dict through `_MutateStep` to the Rust pass.
4. In `UniformMutationPass.execute_with_sampling_mode` and
   `S5FMutationPass.execute_with_sampling_mode`: build the
   site-weight vector with the segment-rate multiplier; drop
   zero-weight sites from support; sample as before.
5. Default-rate fast path: skip the per-position multiply when
   all rates are 1.0.
6. Plan signature includes `segment_rates`.
7. Manifest `models.shm.segment_rate_support = {available: True,
   buckets: [...]}`.

Out of scope for v1:
- Per-region (CDR/FR) rates.
- Per-allele SHM hotspots.
- Per-segment AIRR counter fields.
- AID/UNG biology beyond what S5F already approximates.

This audit only proposes the architecture and pins absences;
the implementation slice flips the relevant `pin_absence_*` to
`pin_present_*` in lockstep.

---

## 12. Test surface — what this audit pins

Mirrored in
[`tests/test_shm_segment_rate_contract.py`](../tests/test_shm_segment_rate_contract.py).

### `pin_scaffold_*` — today's contract

1. **SHM mutates the full pool today** — behavioural: under
   `count=50`, mutations land in V, NP1, D, NP2, J in
   proportion to region length.
2. **Site→segment lookup is structurally available** — source-
   level pin on `Region.segment` field's presence.
3. **`Segment` enum has the five variants** the future rate
   surface maps onto.
4. **No `segment_rates` kwarg today** — source-level pin on
   `Experiment.mutate(...)` signature.
5. **Replay byte-deterministic for SHM** — re-pinned for cross-
   doc traceability with the slice's "preserve under default
   rates" requirement.
6. **Productive-only preserves triad under SHM today** — re-
   pinned baseline; slice must preserve under any non-degenerate
   rate vector.
7. **`n_mutations` is the global biological counter today** —
   re-pinned from the SHM model audit.

### `pin_absence_*` — gaps Slice 1 closes

8. **No per-segment mutation counters** — AIRR records carry
   `n_mutations` only.
9. **No `segment_rate_support` field in cartridge manifest** —
   the manifest's `models.shm` block doesn't advertise the
   capability.
10. **No `_segment_rates` field on `_MutateStep`** — source-
    level pin on the pipeline-IR dataclass.

### Doc anchor

11. The audit doc references the contract file; the 14-section
    structure stays intact.

---

## 13. Out of scope

Documented here so a future contributor doesn't accidentally
expand the slice.

- **Per-region (CDR / FR) rates** within V. Real SHM has rate
  spikes at CDR1 / CDR2; v1 keeps segment-level resolution.
- **Per-allele SHM weighting.** Some allele families are known
  to mutate at different rates; v1 doesn't capture this.
- **AID/UNG biology.** The S5F context kernel already
  approximates the empirical outcome; adding AID directly is a
  separate biology slice.
- **Affinity-driven SHM ramping.** Selection coefficients
  shape SHM in germinal centres; v1 doesn't model selection.
- **Per-segment AIRR counters.** Deferred to a separate
  diagnostic slice (audit §7).
- **C-segment SHM.** The engine has no C-region passes; this
  audit doesn't introduce them.
- **Indel SHM.** Real SHM produces a small fraction of indel
  events; v1 only produces substitutions (matching today's
  passes).

---

## 14. Summary table

| Concern | Recommendation |
|---|---|
| Current targeting | Full assembled pool; confirmed mutations land in V, NP1, D, NP2, J. |
| Site→segment lookup | Cleanly determinable via `sim.sequence.regions`; no architectural blocker. |
| Per-segment rate granularity | **Segment-family** (V / D / J / NP); 4 buckets. NP1 / NP2 combined under "NP". |
| S5F + segment rates weighting | `final = S5F_mutability(site) * segment_rate(segment_at(site))`. |
| Uniform + segment rates | `final = segment_rate(segment_at(site))`; piecewise-constant. |
| Zero-rate sites | Dropped from support before normalisation. |
| Productive-only ordering | Segment-rate filtering happens BEFORE contract admissibility. |
| All-rejected-sites edge | Permissive: silent skip. Strict: `ConstraintSampling` error. Same as today's empty-support rules. |
| Trace addresses | No new addresses; rates are compile-time pass parameters. |
| Plan signature | Includes `segment_rates` so cross-rate replay fails the signature check. |
| Replay determinism | Byte-identical under same rates; signature mismatch on rate drift. |
| DSL shape | **Dict-keyed**: `mutate(segment_rates={"V": 1.0, "D": 0.2, "J": 0.5, "NP": 0.0})`. Sparse; unknown keys raise. |
| Default behaviour | Omitting the kwarg = today's flat-substrate behaviour (all 1.0). Backwards-compatible. |
| AIRR counters | `n_mutations` stays global. Per-segment counters deferred to a separate slice. |
| Manifest | Add `models.shm.segment_rate_support = {available, buckets}` when the slice lands. |
| Performance | Default rates: fast path skips the multiply. Non-default: O(pool_len) overhead. |
| Pre-flight bugs found | **None.** Audit clears the slice for implementation. |

The per-segment SHM rate slice is **architecturally
straightforward**: the IR already carries segment info per
position, the constrain-before-propose primitive plugs in
cleanly, and replay/productive-only/clonal compositions all
hold by construction. The biology decision is the rate-bucket
granularity (segment-family vs per-region), and the audit
recommends the simpler four-bucket form.

After this audit ships, the implementation slice should be a
clean, focused diff: one DSL kwarg, two pass-execution
modifications, one plan-signature update, one manifest field.
No new trace addresses, no new validator surface, no AIRR
record changes. The biology lands without architectural
disturbance.
