# Markov NP Base Generator — Audit + Implementation Shipped

**Status: shipped.** The audit's recommended slice landed
end-to-end. The engine now carries an `NpBaseGenerator` trait
with per-position support keyed on the previously emitted
base; `kind="markov"` cartridges lower through the bridge's
new `markov_transitions` kwarg into `MarkovBaseGenerator`;
plan-signature folding, replay determinism, and the
productive-only contract triad all hold under Markov. Legacy
`NP_transitions` / `NP_first_bases` auto-lift remains
deferred — see §7 below.

Companion artefacts:

- [`tests/test_np_markov_base_generator_contract.py`](../tests/test_np_markov_base_generator_contract.py)
  — 20 pins (flipped to post-implementation state).
- [`tests/test_np_markov_base_generator_implementation.py`](../tests/test_np_markov_base_generator_implementation.py)
  — 13 behaviour tests covering the user's greenlight surface
  (lowering, dependency, replay round-trip, signature-gate
  mismatch, productive-only, byte-identical legacy
  signatures, manifest flip).
- [`tests/test_np_markov_release.py`](../tests/test_np_markov_release.py)
  — release-tier composition tests (productive IGH full
  stack + replay round-trip + deterministic walk).

## Implementation surface

| Layer | Where | What it carries |
|---|---|---|
| **Trait + generators** | [`engine_rs/src/passes/generate_np/np_base_generator.rs`](../engine_rs/src/passes/generate_np/np_base_generator.rs) | `pub trait NpBaseGenerator { fn support(position, previous); fn signature(); }`. Concretes: `UniformNpGenerator`, `CategoricalNpGenerator`, `MarkovBaseGenerator { first_base: [f64;4], transitions: [[f64;4];4] }`. Private `DistributionNpGenerator` adapter keeps the legacy `GenerateNPPass::new(_, _, Box<dyn Distribution<Output=u8>>)` constructor working. |
| **Pass field + signature fold** | [`engine_rs/src/passes/generate_np.rs`](../engine_rs/src/passes/generate_np.rs) | `base_generator: Box<dyn NpBaseGenerator>` replaces `base_dist`. `parameter_signature` folds `base_generator.signature()`; legacy wrappers return byte-identical `fmt_byte_dist`-shaped strings, Markov flattens 5 rows. |
| **Per-position loop** | [`engine_rs/src/passes/generate_np/execution.rs`](../engine_rs/src/passes/generate_np/execution.rs) | `let mut previous: Option<u8> = None;` ahead of the loop; updated to `Some(base)` only after the trace record is committed. |
| **Sample + replay validator** | [`engine_rs/src/passes/generate_np/sampling.rs`](../engine_rs/src/passes/generate_np/sampling.rs) | `sample_base` and `validate_replayed_np_base` take `previous: Option<u8>`. Per-position support is materialised once via `base_generator.support(index, previous)` and fed through the unchanged `sample_base_with_admit_mask` / `sample_filtered_with_policy` helpers via a private `SupportPairsDist` adapter — zero churn on the generic helpers. |
| **PyO3 bridge** | [`engine_rs/src/python/plan.rs`](../engine_rs/src/python/plan.rs) | `push_generate_np(..., markov_transitions: Option<Vec<Vec<(u8,f64)>>>=None)` kwarg. When present, requires `base_pairs` (first-base row) and builds `MarkovBaseGenerator::from_first_and_rows`. New `parse_canonical_base_weights` helper enforces row completeness. |
| **Lowering** | [`src/GenAIRR/_dataconfig_extract.py`](../src/GenAIRR/_dataconfig_extract.py) | `_np_bases_from_models` returns the first-base row (all 4 entries for Markov; positive-only for `empirical_first_base`); new `_np_markov_transitions_from_models` lowers the 4×4 matrix in canonical A/C/G/T order. The prior `NotImplementedError` stop-and-report is gone. |
| **Pipeline IR thread** | [`src/GenAIRR/_pipeline_ir.py`](../src/GenAIRR/_pipeline_ir.py), [`src/GenAIRR/experiment.py`](../src/GenAIRR/experiment.py), [`src/GenAIRR/_compile.py`](../src/GenAIRR/_compile.py) | `_RecombineStep` carries `np{1,2}_markov_transitions`; `Experiment.recombine` reads them off the cartridge defaults; `_lower_recombine` passes them to `push_generate_np`. |
| **Normalisation helper** | [`src/GenAIRR/_normalize.py`](../src/GenAIRR/_normalize.py) | `_to_immutable_byte_pair_matrix` freezes the matrix into the hashable tuple-of-tuples-of-tuples shape `_RecombineStep` requires. |
| **Manifest** | [`src/GenAIRR/dataconfig/data_config.py`](../src/GenAIRR/dataconfig/data_config.py) | `models.np_base_models.supported_kinds = ["uniform", "empirical_first_base", "markov"]`; `deferred_kinds = []`. `legacy_fallback` stays `False`. |

**Pre-flight finding (Q9 below): clean yes — no stop-and-report
condition.** The full call graph from
`GenerateNPPass::execute_with_sampling_mode` through the
admit-mask / replay / fast / slow paths consumes the NP base
distribution's `Vec<(u8, f64)>` raw support pairs at every
non-trivial site, NEVER as a `&dyn Distribution<Output = u8>`
trait object passed deep into the narrowing logic. The
generic helpers `sample_base_with_admit_mask`
([`engine_rs/src/dist/filtered.rs:220-228`](../engine_rs/src/dist/filtered.rs#L220-L228))
and `sample_filtered_with_policy`
([`engine_rs/src/dist/filtered.rs:153-165`](../engine_rs/src/dist/filtered.rs#L153-L165))
are themselves generic over `D: Distribution + ?Sized` but
immediately materialise `Vec<(u8, f64)>` via `dist.support()`
and work entirely on the raw pairs from there. The admit-mask
observer and `JunctionStopState` never see a `Distribution`
reference. The refactor is one self-contained engine slice.

---

## Existing scaffolding the contract pins

| Surface | Where | What it pins |
|---|---|---|
| `GenerateNPPass::execute_with_sampling_mode` single entry point | [`engine_rs/src/passes/generate_np/execution.rs:9`](../engine_rs/src/passes/generate_np/execution.rs#L9) | The only execution entry. Serves replay / fast / slow / unconstrained paths via guarded branches inside `sample_base` — no separate fast-RNG shortcut function. The new slice modifies this one entry. |
| Per-position loop body | [`engine_rs/src/passes/generate_np/execution.rs:98-122`](../engine_rs/src/passes/generate_np/execution.rs#L98-L122) | `for i in 0..length` walks positions sequentially, draws via `sample_base(...)`, records via `record_choice(ChoiceValue::Base(base))`, pushes via `push_nucleotide(...)`. **No per-iteration state** beyond `i` survives across iterations — the previous emitted base is read from the pool / local variable, not from a generator's internal field. |
| `sample_base` (the per-position draw) | [`engine_rs/src/passes/generate_np/sampling.rs:97`](../engine_rs/src/passes/generate_np/sampling.rs#L97) | Routes to replay (`replay_cursor.is_some()`), fast-path (`admit_mask.is_some()`), slow path (contracts but no mask), or unconstrained `self.base_dist.sample(ctx.rng)`. The replay path's signature is the audit's **central wiring change** — needs `previous: Option<u8>` threaded in. |
| `sample_base_with_admit_mask` is generic over `D: Distribution + ?Sized` | [`engine_rs/src/dist/filtered.rs:220-228`](../engine_rs/src/dist/filtered.rs#L220-L228) | Reads `dist.support()` once, then filters / weights / inverse-CDF over `Vec<(u8, f64)>`. **Never calls `dist.sample(rng)`.** A Markov generator that returns position-conditional support drops in without touching this helper. |
| `sample_filtered_with_policy` is generic over `D: Distribution + ?Sized` | [`engine_rs/src/dist/filtered.rs:153-165`](../engine_rs/src/dist/filtered.rs#L153-L165) | Same pattern — calls `dist.support()` once, predicate-filters the pairs, samples weighted from survivors. Markov drops in cleanly. |
| `JunctionStopState::build` — single-shot precompute | [`engine_rs/src/passes/generate_np/execution.rs:65`](../engine_rs/src/passes/generate_np/execution.rs#L65) | Called ONCE per `execute_with_sampling_mode` invocation. Computes the junction frame / anchor / codon layout for the chosen NP length. Returns a state object stashed for per-position queries. |
| Per-position admit-mask via observer | [`engine_rs/src/passes/generate_np/execution.rs:108`](../engine_rs/src/passes/generate_np/execution.rs#L108) | `builder.current_admit_mask()` returns the per-position `u8` mask lazily — depends on the committed pool prefix, so the mask differs per position even within a single execute call. Composes naturally with per-position Markov state. |
| Trace recording per position | [`engine_rs/src/passes/generate_np/execution.rs:119`](../engine_rs/src/passes/generate_np/execution.rs#L119) | `ctx.trace.record_choice(base_choice_address, ChoiceValue::Base(base))` emits one record per position. Addresses: `np.np1.bases[i]` / `np.np2.bases[i]` (Slice 1 surface). **No new addresses needed.** Replay walks the recorded sequence; `previous_base` at position `i` is the recorded value at position `i-1`. |
| `validate_replayed_np_base` | [`engine_rs/src/passes/generate_np/sampling.rs:221`](../engine_rs/src/passes/generate_np/sampling.rs#L221) | Replay-side counterpart to `sample_base`. Mirrors the three branches (mask / contracts / no-contracts). The **central gap**: receives `index: u32` but no `previous: Option<u8>`. The new slice adds the param so the Markov generator's per-position support can be reconstructed deterministically. |
| `parameter_signature` folds `base_dist` via `fmt_byte_dist` | [`engine_rs/src/passes/generate_np.rs:113-122`](../engine_rs/src/passes/generate_np.rs#L113-L122) + [`engine_rs/src/passes/paramsig.rs:fmt_byte_dist`](../engine_rs/src/passes/paramsig.rs) | Slice A discipline. `fmt_byte_dist` calls `support()` and renders `[(v:w),(v:w),…]`. A Markov generator's signature MUST flatten 5 rows (first_base + 4 transition rows) into one deterministic string — `fmt_byte_dist` as-is would only see one of them via a single `support()` call. A new `fmt_np_base_generator` helper is the cleanest extension. |
| Existing `kind="markov"` rejection at lowering | [`src/GenAIRR/_dataconfig_extract.py::_np_bases_from_models`](../src/GenAIRR/_dataconfig_extract.py) | Raises `NotImplementedError` with the audit-doc-pointing message. The implementation slice REMOVES this raise and routes through a new bridge function that wires the typed spec into `MarkovBaseGenerator`. |
| Python `NpBaseModelSpec(kind="markov")` validation | [`src/GenAIRR/reference_models.py:NpBaseModelSpec.validate`](../src/GenAIRR/reference_models.py) | Already complete: validates per-base alphabet, weight finiteness, partial-matrix rejection. **The implementation slice changes nothing in the spec layer.** |

---

## 1. Q1 — Generator abstraction

The user's recommendation:

```rust
trait NpBaseGenerator {
    fn support(&self, position: usize, previous: Option<u8>) -> Vec<u8>;
    fn weight(&self, position: usize, previous: Option<u8>, base: u8) -> f64;
    fn sample(...);
}
```

### Recommended trait shape (refined)

```rust
pub trait NpBaseGenerator {
    /// Enumerate the discrete `(base, weight)` support at this
    /// position given the previously emitted base (`None` only at
    /// position 0). Same shape as `Distribution::support()` — the
    /// existing admit-mask / `sample_filtered_with_policy` helpers
    /// consume this verbatim.
    fn support(&self, position: usize, previous: Option<u8>) -> Vec<(u8, f64)>;

    /// Sample one base at this position. Used by the unconstrained
    /// fast path (no contracts, no admit mask). Implementors that
    /// derive `sample` from `support` are encouraged to use the
    /// shared `inverse_cdf_from_support` helper.
    fn sample(&self, rng: &mut Rng, position: usize, previous: Option<u8>) -> u8;

    /// Deterministic per-cartridge identity string for plan-signature
    /// folding. Markov generators flatten all rows in canonical
    /// A → C → G → T order with row separators; uniform / independent
    /// generators return the existing single-support shape so legacy
    /// trace signatures stay byte-identical.
    fn signature(&self) -> String;
}
```

### Why a NEW trait and not extend `Distribution`

Three reasons:

1. **`Distribution<Output = u8>` is shared with PCR / quality /
   contaminant passes** ([`engine_rs/src/passes/corrupt/*.rs`](../engine_rs/src/passes/corrupt/)).
   Those are stateless by design — adding a `previous` parameter
   would force every consumer to thread an `Option<u8>` that is
   meaningless to them. Even with a default of `None`, the trait
   contract becomes incoherent.
2. **The Markov generator is intrinsically a per-pass concept.**
   The trait name documents intent.
3. **Migration is symmetric:** `UniformBase` and `CategoricalBase`
   stay as `Distribution` implementations for their existing
   consumers; new wrapper types `UniformNpGenerator` and
   `CategoricalNpGenerator` implement `NpBaseGenerator` for NP
   use. The wrappers ignore `previous` and `position` (they
   delegate to the inner `Distribution::support()`).

### Field migration on `GenerateNPPass`

Today:

```rust
pub struct GenerateNPPass {
    np_segment: Segment,
    length_dist: Box<dyn Distribution<Output = i64>>,
    base_dist: Box<dyn Distribution<Output = u8>>,
}
```

After the slice:

```rust
pub struct GenerateNPPass {
    np_segment: Segment,
    length_dist: Box<dyn Distribution<Output = i64>>,
    base_generator: Box<dyn NpBaseGenerator>,
}
```

### Pinned

- `pin_absence_no_np_base_generator_trait` — the trait doesn't
  exist yet.
- `pin_absence_no_markov_base_generator_type` — the
  `MarkovBaseGenerator` concrete type doesn't exist yet.

---

## 2. Q2 — Replay

### No new trace addresses

The existing per-position recording at
[`execution.rs:119`](../engine_rs/src/passes/generate_np/execution.rs#L119)
captures every emitted base in order at sequential
`np.np1.bases[i]` / `np.np2.bases[i]` addresses. Replay walks
the sequence; the previous-base state at position `i` is just
`trace[i-1]` — a sequential walk reconstructs it without any
new trace addresses or payloads.

### Replay-validator signature change

The central wiring change in the slice. Today:

```rust
fn validate_replayed_np_base(
    base_dist: &dyn Distribution<Output = u8>,
    index: u32,
    admit_mask: Option<u8>,
    junction_stop_state: Option<&JunctionStopState>,
    recorded_byte: u8,
    /* ... */
) -> Result<(), ReplayError>
```

After the slice:

```rust
fn validate_replayed_np_base(
    base_generator: &dyn NpBaseGenerator,
    index: u32,
    previous: Option<u8>,                       // NEW
    admit_mask: Option<u8>,
    junction_stop_state: Option<&JunctionStopState>,
    recorded_byte: u8,
    /* ... */
) -> Result<(), ReplayError>
```

The replay loop in `sample_base` maintains a `let mut previous:
Option<u8> = None` across iterations and forwards it.

### Replay determinism guarantee

Same-seed replay reproduces all NP bases byte-for-byte under
all three kinds (uniform / empirical_first_base / Markov)
because the trace records the emitted base, not the
conditioning state. The generator's `sample()` at position
`i` consults `previous = Some(trace[i-1])` for the Markov
case; the recorded value re-emerges deterministically.

### Pinned

- `pin_scaffold_np_trace_addresses_unchanged_under_markov_proposal`
- `pin_absence_no_previous_base_param_on_validate_replayed_np_base`

---

## 3. Q3 — Existing generators (concrete types)

After the slice, three implementations of `NpBaseGenerator`:

### `UniformNpGenerator`

Wraps `UniformBase`. Ignores `previous` and `position`. Returns
the canonical 4-way support `[(b'A', 1.0), (b'C', 1.0), (b'G',
1.0), (b'T', 1.0)]`. **Plan signature byte-identical to today's
`UniformBase`** — legacy trace files replay unchanged.

### `CategoricalNpGenerator` (the existing first-base case)

Wraps `CategoricalBase`. Ignores `previous` and `position`.
Returns the cartridge-owned categorical's support. **Plan
signature byte-identical to today's `CategoricalBase`** —
existing cartridges with `kind="empirical_first_base"` replay
unchanged.

### `MarkovBaseGenerator`

New. Owns:

- `first_base: [f64; 4]` — A/C/G/T weights for position 0.
- `transitions: [[f64; 4]; 4]` — row-major `from_base × to_base`
  weights for positions 1+. Indexed in canonical alphabetical
  order: `0 → A`, `1 → C`, `2 → G`, `3 → T`.

`support(position, previous)`:
- `position == 0 || previous.is_none()` → returns `first_base`
  as `Vec<(u8, f64)>`.
- `position > 0 && previous == Some(b)` → returns
  `transitions[base_to_idx(b)]` as `Vec<(u8, f64)>`.

### Existing `UniformBase` / `CategoricalBase` stay

These distributions remain in
[`engine_rs/src/dist/`](../engine_rs/src/dist/) for use by PCR
/ quality / contaminant / sample_base / etc. — they're not
NP-specific. The new trait introduces wrapper types; the
underlying `Distribution` implementations are unchanged.

### Pinned

- `pin_absence_no_markov_base_generator_type`
- `pin_absence_no_uniform_np_generator_wrapper`
- `pin_absence_no_categorical_np_generator_wrapper`

---

## 4. Q4 — Productive-only constraints

### Generator support × admit mask intersection

Today:

```rust
// Slow path (contracts active, no mask)
sample_filtered_with_policy(rng, base_dist.as_ref(), |b| contracts.admits_typed(b), ...)
```

`sample_filtered_with_policy` materialises `base_dist.support()`
once, predicate-filters by `contracts.admits_typed`, then
weighted-samples from survivors. **The intersection happens
before sampling**, not after — empty intersection routes to
the `NP_BASE_EMPTY_SUPPORT` sentinel (`b'N'`).

The Markov generator's per-position support replaces
`base_dist.support()` in this call. The intersection then
becomes `markov_generator.support(i, previous) ×
admit_mask(position i)`. Same composition, same sentinel
behaviour.

### Empty-support behaviour preserved

`EmptySupport::Sentinel(b'N')` ([`sampling.rs:20`](../engine_rs/src/passes/generate_np/sampling.rs#L20))
remains the sentinel. Strict mode raises
`PassError::ConstraintSampling`; permissive emits `b'N'`. No
change.

### Edge case: Markov with admit mask zeroing the categorical's support

If the Markov row from a given `previous` base zeroes every
admit-mask survivor (e.g. `transitions[A] = {A: 1.0, C: 1.0,
G: 1.0, T: 1.0}` filtered by an admit mask that admits only
`A`, but the row's `A` weight is 0 because the cartridge
authored a degenerate matrix), the same empty-support sentinel
fires. The generator's `support()` returning a zero-weight
pair is treated identically to the categorical filtering
yielding nothing — predicate-filter sees no survivor with
positive weight, sentinel.

### Pinned

- `pin_scaffold_admit_mask_intersection_runs_before_sampling`
- `pin_scaffold_np_base_empty_support_sentinel_unchanged`

---

## 5. Q5 — Plan signature

### Markov transition matrix in `parameter_signature`

Markov has TWO supports to fold:

1. The first-base categorical (4 weights, used at position 0).
2. Four transition rows (4 × 4 weights, used at positions 1+).

The audit recommends a flat deterministic encoding:

```text
base=markov:first=[(65:w_A),(67:w_C),(71:w_G),(84:w_T)]
  |from=A:[(65:w_AA),(67:w_AC),(71:w_AG),(84:w_AT)]
  |from=C:[...]
  |from=G:[...]
  |from=T:[...]
```

Row ordering canonical (A → C → G → T) inside both the
top-level "from" iteration AND inside each row's `(to:weight)`
list — matches the existing `fmt_byte_dist` convention.

The implementation slice adds a helper:

```rust
pub fn fmt_np_base_generator(g: &dyn NpBaseGenerator) -> String {
    g.signature()
}
```

Each implementor controls its own signature shape; uniform /
categorical implementors return the existing
`fmt_byte_dist`-compatible string for byte-identical legacy
signatures.

### Backwards compatibility (the key constraint)

For the slice to ship without breaking existing trace
replays, `UniformNpGenerator::signature()` and
`CategoricalNpGenerator::signature()` MUST return the exact
strings `fmt_byte_dist` produces today:

- Uniform: `[(65:1.0),(67:1.0),(71:1.0),(84:1.0)]`
- Categorical: the existing `(b, w)` shape from
  `CategoricalBase::support()`.

This is a hard release-gate requirement; the contract pins
freeze the expected substrings.

### Pinned

- `pin_scaffold_uniform_plan_signature_substring_is_canonical`
- `pin_absence_no_markov_signature_format_today`

---

## 6. Q6 — DSL / cartridge

### Spec-layer change: NONE

`NpBaseModelSpec(kind="markov", first_base=…, transitions=…)`
already validates per the prior slice. The implementation
slice's only Python-side change is in
`_dataconfig_extract._np_bases_from_models`:

```python
# Before (today):
if spec.kind == "markov":
    raise NotImplementedError(...)

# After (the slice):
if spec.kind == "markov":
    return _markov_payload(spec)  # → (first_base, transitions) tuple
```

The resolver's return shape extends from `Optional[List[Tuple[int,
float]]]` to a discriminated union that carries either the
flat categorical pairs OR the markov payload. The bridge
dispatches.

### Bridge-layer change

A new PyO3 entry point or an extended `push_generate_np`
signature. The audit recommends the latter — extend
`base_pairs` with an optional `markov_transitions` kwarg:

```python
plan.push_generate_np(
    "NP1",
    length_pairs,
    *,
    base_pairs=None,           # first-base / uniform → existing
    markov_transitions=None,   # NEW — None | List[(from_base, [(to_base, weight)])]
)
```

When `markov_transitions is not None`, the bridge requires
`base_pairs is not None` (the first-base row) and instantiates
`MarkovBaseGenerator`. When both are `None`, `UniformBase` /
`UniformNpGenerator` applies (existing fast path). When only
`base_pairs is not None`, the existing `CategoricalBase` /
`CategoricalNpGenerator` applies.

### Pinned

- `pin_scaffold_np_base_model_spec_already_validates_markov`
- `pin_absence_no_markov_transitions_kwarg_on_push_generate_np`

---

## 7. Q7 — Compatibility

### Uniform / no-spec byte-identical

Confirmed cleanly:

- The bridge defaults `base_pairs=None` to `UniformBase`
  (today) → `UniformNpGenerator` (after).
- `UniformNpGenerator::sample()` delegates to
  `UniformBase::sample()` — same RNG consumption shape.
- `UniformNpGenerator::signature()` returns the exact existing
  4-way support string.

### `empirical_first_base` byte-identical

Same pattern. `CategoricalNpGenerator` is a `Distribution`
adapter that ignores `previous` and `position`, delegating to
the inner `CategoricalBase`. Replay byte-identical.

### Legacy `NP_transitions` / `NP_first_bases` still NOT auto-lifted

The prior slice's stop-and-report on auto-lift remains in
force. This audit's Markov implementation slice does NOT
change the resolver cascade in
`_dataconfig_extract.py` — the cascade stays **typed →
uniform**. The bundled cartridges' orphan `NP_transitions`
fields stay unconsumed. A separate slice may later add
auto-lift under an explicit opt-in flag.

### Pinned

- `pin_scaffold_legacy_np_transitions_still_orphan_after_markov_slice`

---

## 8. Q8 — Edge cases the implementation slice must handle

| Case | Expected behaviour |
|---|---|
| `length == 0` | Loop body doesn't run; no `previous` state needed; no trace records emitted. Markov resolver doesn't even instantiate the generator's per-position lookup. |
| `length == 1` | Single position; `previous = None`; uses `first_base` row. Symmetric with `kind="empirical_first_base"` for that single position. |
| `length > 1` | Position 0 uses `first_base`; positions 1+ use `transitions[previous]`. |
| Pre-existing trace recorded under uniform NP base | Replay against a Markov cartridge fails the plan-signature gate before consuming choices (Slice A discipline). |
| Pre-existing trace recorded under Markov | Replay against the same Markov cartridge passes byte-identically. |
| Productive-only with Markov | Per-position admit mask intersects per-position Markov support. Empty intersection at any position triggers the existing `b'N'` sentinel. |
| Markov generator with one row producing zero total weight (e.g. `transitions[A] = {C: 0.0, G: 0.0, T: 0.0}`) | Spec-layer validation rejects this at `NpBaseModelSpec.__post_init__` ("at least one weight must be strictly positive"). |
| Markov with `b'N'` recorded mid-stream from a prior empty-support iteration | The validator's `previous = Some(b'N')` is outside the canonical alphabet. The generator must either (a) treat unknown previous as `None` (degrade to `first_base`) or (b) raise. The audit recommends (a): mid-stream `b'N'` was already a sentinel; reverting to `first_base` for the next position is the cleanest fallback that doesn't propagate failure. |
| Replay with a recorded base that the Markov support at that position assigns zero weight | Treated as a recorded value out of support → replay raises (same as today's `validate_replayed_np_base` no-contracts path). |
| Mismatched plan signature (matrices differ) | Replay fails the signature gate. |
| Two Markov cartridges with identical first_base but different transitions | Different signatures → reject; one cartridge's trace cannot replay against the other. |

---

## 9. Q9 — Performance

### Per-position cost

Markov generator's `support(i, previous)` does one row lookup
(`transitions[base_to_idx(previous)]`) — constant cost. The
returned `Vec<(u8, f64)>` is 4 entries; same shape as today's
`UniformBase` / `CategoricalBase` support.

The admit-mask × support intersection cost is unchanged — the
filter / weight / sample logic in `sample_base_with_admit_mask`
and `sample_filtered_with_policy` walks a 4-entry vec.

### Allocation cost

`support()` allocates a 4-entry `Vec` per call. Per
simulation, NP1 + NP2 average ~6 bases × 2 = 12 allocations.
The audit recommends accepting this allocation cost for v1;
a follow-up optimisation could return `&[(u8, f64)]` via an
inline-array if profiling shows it's hot.

### Plan-signature cost

Markov signature is ~120-150 bytes of ASCII (5 rows × ~30
bytes each). Added once per plan compile. Negligible.

### Pinned

- `pin_scaffold_existing_sample_helpers_take_raw_pairs`
- `pin_scaffold_per_position_admit_mask_recomputed_each_iteration`

---

## 10. Trace / replay impact

### Zero new addresses

The existing per-position `np.np1.bases[i]` / `np.np2.bases[i]`
addresses carry the full Markov state. Replay reconstructs
previous-base from the prior record's `ChoiceValue::Base`.
Pre-existing traces replay unchanged.

### Replay-determinism guarantee

Same-seed + same-cartridge replay reproduces all NP bases
byte-for-byte under Markov, just as it does under uniform /
empirical_first_base today. The trace records the emitted
base; the generator's `support()` consults `previous = trace[i-1]`
when re-replaying. Deterministic by construction.

### Plan-signature folding

`MarkovBaseGenerator::signature()` returns a deterministic
flat string covering all 5 rows. A same-cartridge replay
against a different cartridge's matrix fails the signature
gate before any choice is consumed.

### Pinned

- `pin_scaffold_np_trace_addresses_unchanged_under_markov_proposal`

---

## 11. Manifest extension

The `models.np_base_models` block updates two fields:

```python
manifest["models"]["np_base_models"] = {
    "models": {
        "NP1": {"kind": "markov", "first_base_keys": [...], "row_keys": [...]},
        ...
    },
    "legacy_fallback": False,  # unchanged
    "legacy_np_transitions_present": True,  # unchanged
    "legacy_np_first_bases_present": True,  # unchanged
    "supported_kinds": ["uniform", "empirical_first_base", "markov"],  # CHANGED — markov added
    "deferred_kinds": [],  # CHANGED — empty now
    "in_plan_signature": True,
    "in_content_hash": False,
}
```

Per-region entries gain optional `first_base_keys` and
`row_keys` for Markov authors who want to see at a glance
which transitions are populated. Backwards compatible — the
keys are additive.

### Pinned

- `pin_absence_no_markov_in_supported_kinds_today`

---

## 12. Implementation order (recommended)

A single self-contained engine slice. Six sub-steps:

1. **New trait + wrappers** —
   `engine_rs/src/passes/generate_np/np_base_generator.rs`:
   - `pub trait NpBaseGenerator { fn support / fn sample / fn signature; }`.
   - `UniformNpGenerator` wrapper (byte-identical signature).
   - `CategoricalNpGenerator` wrapper (byte-identical signature).
   - `MarkovBaseGenerator` with `first_base: [f64; 4]` +
     `transitions: [[f64; 4]; 4]` + `support()` + `sample()` +
     `signature()`.

2. **`GenerateNPPass` field migration** —
   change `base_dist: Box<dyn Distribution<Output = u8>>` to
   `base_generator: Box<dyn NpBaseGenerator>`. Update
   constructors. The struct shape changes; tests that
   construct `GenerateNPPass` directly need updating.

3. **Loop body wiring** —
   `execution.rs::execute_with_sampling_mode` adds
   `let mut previous: Option<u8> = None;` before the loop;
   sets `previous = Some(base);` at the bottom of each
   iteration after the trace record.

4. **`sample_base` + `validate_replayed_np_base` signature
   change** — both take `previous: Option<u8>`. Forward to
   `base_generator.support(i, previous)` instead of
   `base_dist.support()`.

5. **PyO3 bridge** — extend `push_generate_np` with an
   optional `markov_transitions` kwarg. `_compile.py` reads
   the typed spec's `kind` and dispatches.

6. **Resolver update** — remove the `NotImplementedError` in
   `_np_bases_from_models` for `kind="markov"`; return the
   Markov payload.

7. **Manifest update** — flip `supported_kinds` / `deferred_kinds`.

8. **Tests** — implementation tests for Markov: per-region
   ordering, conditional sampling produces base-dependent
   distributions (e.g. matrix `A → 100% C` produces all `C`
   bases except position 0; matrix `* → uniform first_base`
   recovers `empirical_first_base` behaviour byte-equivalent);
   replay round-trip; plan signature divergence between
   distinct matrices; productive-only triad preservation;
   manifest exposure.

Cost estimate: ~150 lines Rust (trait + wrappers + Markov
generator + tests) + ~80 lines Python (bridge + resolver + manifest)
+ ~50 lines new tests. Single self-contained slice.

---

## 13. Test surface — what this audit pins

Mirrored in
[`tests/test_np_markov_base_generator_contract.py`](../tests/test_np_markov_base_generator_contract.py).

### `pin_scaffold_*` — pre-existing surfaces the slice builds on

1. `GenerateNPPass::execute_with_sampling_mode` is the single
   execution entry point.
2. Per-position loop carries no per-iteration state — the
   previous base is a candidate for a clean refactor-local
   `let mut`.
3. `sample_base_with_admit_mask` is generic over
   `D: Distribution + ?Sized` and consumes raw support pairs.
4. `sample_filtered_with_policy` follows the same pattern.
5. NP trace addresses are unchanged from prior slices.
6. `JunctionStopState::build` runs once per execute; admit
   masks recompute per position via the observer.
7. The `EmptySupport::Sentinel(b'N')` policy is unchanged.
8. Spec-layer `NpBaseModelSpec(kind="markov")` already
   validates.
9. Legacy `NP_transitions` / `NP_first_bases` remain orphans.
10. Plan signature folds `base_dist.support()` via `fmt_byte_dist`.

### `pin_present_*` — current stop-and-report behaviour

11. `kind="markov"` currently raises `NotImplementedError` at
    lowering. The implementation slice flips this in lockstep.

### `pin_absence_*` — gaps the slice closes

12. No `NpBaseGenerator` trait exists.
13. No `MarkovBaseGenerator` type exists.
14. No `UniformNpGenerator` / `CategoricalNpGenerator`
    wrapper types exist.
15. No `previous: Option<u8>` parameter on
    `validate_replayed_np_base` or `sample_base`.
16. No `markov_transitions` kwarg on `push_generate_np`.
17. `models.np_base_models.supported_kinds` does NOT include
    `"markov"` today.

### Doc anchor

18. The audit doc exists and references the contract file;
    structure intact.

---

## 14. Out of scope

Documented here so a future implementer doesn't accidentally
expand the work.

- **P-nucleotides / palindromic additions.** Separate audit;
  needs new event payload / region semantics.
- **Multi-step Markov (2-step or 3-step lookback).** v1
  Markov is 1-step. A 2-step model would need
  `previous_pair: (Option<u8>, Option<u8>)` and 4² = 16 rows;
  separate slice.
- **VJ-chain-specific or locus-specific Markov.** Cartridge-
  level only — author per-cartridge.
- **Length-conditional Markov.** Out of scope.
- **Legacy `NP_transitions` auto-lift.** Stays deferred —
  this audit's Markov slice does NOT enable auto-lift. A
  separate slice may later land it under an explicit opt-in.
- **AIRR `np1_markov_realised` provenance field.** Some
  future analytics slice could surface "this NP1's
  conditional probability density"; not in scope here.

---

## 15. Summary table

| Concern | Today | Recommendation |
|---|---|---|
| NP base distribution trait | `Distribution<Output = u8>` (stateless) | Add `NpBaseGenerator` trait with `previous: Option<u8>` context |
| Replace or supplement `Distribution`? | n/a | **Supplement** — `Distribution` stays for non-NP passes |
| Generator concrete types | `UniformBase` / `CategoricalBase` | + `MarkovBaseGenerator`; existing two wrap as `UniformNpGenerator` / `CategoricalNpGenerator` |
| `GenerateNPPass.base_dist` field | `Box<dyn Distribution<Output = u8>>` | Renamed to `base_generator: Box<dyn NpBaseGenerator>` |
| Productive-only admit-mask coupling to trait | None — helpers take raw `(u8, f64)` pairs and `dist.support()` call only | Unchanged |
| Replay validator signature | `validate_replayed_np_base(base_dist, index, mask, ...)` | Adds `previous: Option<u8>` param |
| Trace addresses | `np.npN.length` + `np.npN.bases[i]` | Unchanged — no new addresses |
| Plan-signature folding | `fmt_byte_dist` folds `support()` | `MarkovBaseGenerator::signature()` flattens 5 rows; uniform / categorical wrappers return existing byte-identical strings |
| `NpBaseModelSpec(kind="markov")` Python spec | Validates; lowering raises `NotImplementedError` | Lowering routes through new bridge function — no spec change |
| Bridge `push_generate_np` | `base_pairs: Option<Vec<(u8, f64)>>` only | Add `markov_transitions: Option<...>` kwarg |
| Manifest `supported_kinds` | `["uniform", "empirical_first_base"]` | `["uniform", "empirical_first_base", "markov"]` |
| Manifest `deferred_kinds` | `["markov"]` | `[]` |
| Legacy `NP_transitions` auto-lift | Deferred (orphan field) | Still deferred — separate slice |
| Backwards compatibility (uniform / empirical replay) | Existing traces replay byte-identical | **Hard constraint:** preserved via byte-identical signatures from the wrapper types |
| Performance | n/a | 4-entry Vec allocation per position; ~12 per simulation; negligible |
| Pre-flight bugs found | **None.** The coupling check is clean — admit-mask / `sample_filtered_with_policy` / `sample_base_with_admit_mask` all consume raw support pairs. The Markov slice is one focused engine refactor with no architectural surprises. | — |

The Markov NP base generator is **architecturally tractable**:
the engine's existing `Vec<(u8, f64)>` consumption discipline
means the new trait can compute its per-position support and
hand it to the existing admit-mask / sampler machinery
verbatim. The slice's only non-trivial wiring change is
threading `previous: Option<u8>` from the per-position loop
into `sample_base` / `validate_replayed_np_base` and updating
the field type on `GenerateNPPass`. Backwards compatibility
is preserved via the wrapper types — existing traces replay
byte-identically as long as their plan signatures are
byte-equal.

The recommended next step is the single Markov implementation
slice per §12. P-nucleotides, multi-step Markov, and legacy
auto-lift remain explicitly out of scope per §14.
