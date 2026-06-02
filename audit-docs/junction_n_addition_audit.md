# Junction / N-Addition Modeling — Audit + Typed NP Base Model + Markov Generator Shipped

**Status: Slices 1 + 2 shipped — all three NP base kinds run
end-to-end (`uniform`, `empirical_first_base`, `markov`).
Legacy `NP_transitions` / `NP_first_bases` auto-lift remains
explicitly deferred.**

- **Slice 1 — Typed NP base model** landed end-to-end
  alongside 27 implementation tests
  ([`tests/test_np_base_model_implementation.py`](../tests/test_np_base_model_implementation.py))
  and pin flips in the audit contract. The
  `empirical_first_base` kind lets cartridges author per-NP-
  region weighted-categorical sampling that replaces the
  hardcoded `UniformBase` byte-identically when no spec is
  configured.
- **Slice 2 — Markov NP base generator** landed via a new
  `NpBaseGenerator` trait, `MarkovBaseGenerator` concrete type,
  and byte-identical-signature `UniformNpGenerator` /
  `CategoricalNpGenerator` wrappers. Per-position support is
  conditioned on the previously emitted base; replay
  reconstructs `previous` from the prior recorded
  `np.np{1,2}.bases[i-1]` value — no new trace addresses.
  Pinned by
  [`tests/test_np_markov_base_generator_contract.py`](../tests/test_np_markov_base_generator_contract.py)
  + 13 behaviour tests in
  [`tests/test_np_markov_base_generator_implementation.py`](../tests/test_np_markov_base_generator_implementation.py).
  See companion design
  [`docs/np_markov_base_generator_design.md`](np_markov_base_generator_design.md).

**Supported `NpBaseModelSpec.kind` values today:** `"uniform"`,
`"empirical_first_base"`, `"markov"`.

**Legacy auto-lift remains deferred.** The bundled cartridges'
`DataConfig.NP_transitions` / `DataConfig.NP_first_bases`
dicts are still listed in
`_DOCUMENTED_ORPHAN_DATACONFIG_FIELDS` and the manifest still
reports `np_base_models.legacy_fallback=False`. Auto-lift
would silently change output bytes vs the pre-slice baseline
and reroute every legacy cartridge through Markov sampling
without an explicit opt-in; that decision is a separate
cartridge-migration slice, not this consolidation. Authors who
want Markov today populate `ReferenceEmpiricalModels.np_bases`
explicitly with a typed `NpBaseModelSpec(kind="markov", ...)`.

The audit body below is preserved for traceability; sections
marked **[Shipped]** describe how the recommendations actually
landed.

The audit is upstream-biology, not export-layer — unlike the
recently-shipped FASTQ writer slice, this work *will*
eventually touch engine internals (the `GenerateNPPass` base
distribution wiring, possibly new cartridge spec types). v1
keeps the slice tight to cartridge-owned distributional
defaults; P-nucleotides / palindromic additions / per-context
Markov chains are out of scope per §13.

Companion to
[`tests/test_junction_n_addition_contract.py`](../tests/test_junction_n_addition_contract.py)
which freezes today's surfaces (`pin_scaffold_*`) and the gaps
the implementation slice would close (`pin_absence_*`).

**Pre-flight finding (Q3 below): clean yes — no stop-and-report
condition.** The legacy `DataConfig.NP_transitions` and
`DataConfig.NP_first_bases` fields exist on bundled cartridges
(populated dicts on HUMAN_IGH_OGRDB) but are **never read by
the simulation pipeline**. They are explicitly listed in
`_DOCUMENTED_ORPHAN_DATACONFIG_FIELDS` and surface only as
documented completeness gaps in the cartridge manifest. The
current NP base distribution at every NP1/NP2 position is a
hardcoded uniform `UniformBase` (4-way equal A/C/G/T)
constructed at the PyO3 bridge. No replay-safety hole exists.

---

## Existing scaffolding the contract pins

| Surface | Where | What it pins |
|---|---|---|
| `GenerateNPPass` struct | [`engine_rs/src/passes/generate_np.rs:51-55`](../engine_rs/src/passes/generate_np.rs#L51-L55) | Three fields: `np_segment: Segment`, `length_dist: Box<dyn Distribution<Output = i64>>`, `base_dist: Box<dyn Distribution<Output = u8>>`. The new slice extends `base_dist` to a cartridge-owned typed distribution without changing the struct shape. |
| PyO3 bridge `push_generate_np` | [`engine_rs/src/python/plan.rs:373-392`](../engine_rs/src/python/plan.rs#L373-L392) | Hardcodes `Box::new(UniformBase)` as the `base_dist`. No Python-side parameter for base distribution today; the new slice adds one. |
| `UniformBase` distribution | [`engine_rs/src/dist/uniform.rs:19-26`](../engine_rs/src/dist/uniform.rs#L19-L26) | Uniform-weighted A/C/G/T; `support()` returns `Some([(b'A',1.0),(b'C',1.0),(b'G',1.0),(b'T',1.0)])`. The new typed model returns the same shape with non-uniform weights. |
| NP length resolution order | [`src/GenAIRR/_dataconfig_extract.py:175-224`](../src/GenAIRR/_dataconfig_extract.py#L175-L224) | `reference_models.np_lengths["NP1"]` (typed) → `cfg.NP_lengths["NP1"]` (legacy dict) → uniform `[(0,1.0),...,(6,1.0)]` placeholder. The new typed base model's resolution mirrors this fallback chain. |
| `EmpiricalDistributionSpec` + `ReferenceEmpiricalModels.np_lengths` | [`src/GenAIRR/reference_models.py`](../src/GenAIRR/reference_models.py) | The typed plane for empirical distributions; the new slice adds a parallel typed plane for NP base / first-base distributions. |
| Trace addresses for NP | [`engine_rs/src/address.rs:80-86`](../engine_rs/src/address.rs#L80-L86) | `NP1_LENGTH = "np.np1.length"`, `NP2_LENGTH = "np.np2.length"`, `NP1_BASES_INDEX_PREFIX = "np.np1.bases["`, `NP2_BASES_INDEX_PREFIX = "np.np2.bases["`. The new slice does NOT introduce new addresses. |
| NP `parameter_signature` | [`engine_rs/src/passes/generate_np.rs:113-122`](../engine_rs/src/passes/generate_np.rs#L113-L122) | Folds both `length_dist` and `base_dist` via `fmt_int_dist` + `fmt_byte_dist`. Slice A's discipline holds — a new base distribution's `support()` produces a distinct signature automatically, so replay safety is enforced without any new code in the signature path. |
| `JunctionStopState` admit mask | [`engine_rs/src/passes/generate_np/execution.rs:63-87`](../engine_rs/src/passes/generate_np/execution.rs#L63-L87) | When `productive_only()` is active, the pass attaches an admit-mask observer that narrows NP base sampling to bases that preserve junction triad (no stop, in-frame, anchor preserved). Empty-support sentinel `EmptySupport::Sentinel(b'N')` ([`generate_np/sampling.rs:20`](../engine_rs/src/passes/generate_np/sampling.rs#L20)); strict mode raises, permissive emits `N`. |
| AIRR `np1`/`np2`/`np1_length`/`np2_length` fields | [`engine_rs/src/airr_record/record.rs:77-82`](../engine_rs/src/airr_record/record.rs#L77-L82) | Already first-class — `np1_length` and `np2_length` carry the realised post-trim NP length (unclaimed bases after V/D/J live-call extensions reabsorb adjacent NP positions). No projection work needed for the next slice. |
| `_DOCUMENTED_ORPHAN_DATACONFIG_FIELDS` | [`src/GenAIRR/dataconfig/data_config.py:43-51`](../src/GenAIRR/dataconfig/data_config.py#L43-L51) | Explicit list of `DataConfig` fields that exist on bundled pickles but are NOT consumed by any view / the bridge / the engine. Today's entries include `gene_use_dict`, **`NP_transitions`**, **`NP_first_bases`**, `correction_maps`, `asc_tables`, `p_nucleotide_length_probs`, `dj_pairing_map`. The new slice consumes (and removes) `NP_transitions` / `NP_first_bases` from this list. |
| Manifest empirical-models block | [`src/GenAIRR/dataconfig/data_config.py:530`](../src/GenAIRR/dataconfig/data_config.py#L530) (`np_length_keys`) | Manifest reports `np_length_keys` from `reference_models.np_lengths` only. There is no `np_base_keys` / `np_first_base_keys` block today — the new slice adds one parallel to the lengths surface. |

---

## 1. Q1 — Current N-addition model

### Length distribution

- **Typed source (winning path):** `cfg.reference_models.np_lengths["NP1"]`
  and `["NP2"]` as `EmpiricalDistributionSpec` objects.
- **Legacy fallback:** `cfg.NP_lengths["NP1"]` / `["NP2"]` as
  `{length: prob}` dicts. Consumed by
  `extract_recombine_defaults` when `reference_models` is
  `None` or empty.
- **Placeholder:** uniform `[(0, 1.0), …, (6, 1.0)]` when
  both are absent (legacy `DataConfig` with no NP defaults).
- **Per-experiment override:** the DSL `recombine(np1_lengths=…,
  np2_lengths=…)` kwargs override the cartridge default
  entirely.

### Base distribution

- **Hardcoded:** `UniformBase` (4-way equal A/C/G/T sampling),
  constructed unconditionally at
  [`engine_rs/src/python/plan.rs:387-391`](../engine_rs/src/python/plan.rs#L387-L391).
  Identical for NP1 and NP2.
- No cartridge-owned base distribution today.
- No DSL kwarg to override base distribution today.
- No locus / chain-type variation today.

### Chain-type behaviour

VJ chains run `generate_np.np1` only (no NP2 because no D
segment). VDJ chains run both `np1` and `np2`. Same uniform
base distribution applies to both segments in both chain
types.

### Pinned

- `pin_scaffold_np_length_resolution_typed_wins_then_legacy_falls_back`
- `pin_scaffold_np_base_dist_is_hardcoded_uniform_base`
- `pin_scaffold_chain_type_vj_runs_np1_only_vdj_runs_both`

---

## 2. Q2 — Current missing biology

### P-nucleotides / palindromic additions

**Not modelled.** No DSL surface, no pass, no trace address,
no AIRR field. The legacy `DataConfig.p_nucleotide_length_probs`
field exists on the dataclass with `default_factory=dict` and
is also listed in `_DOCUMENTED_ORPHAN_DATACONFIG_FIELDS`. P-
nucleotide additions would require new biology: a separate
sampling step before NP generation, a new event payload to
distinguish P-nucleotides from N-nucleotides at projection
time, possibly new region semantics (P-nuc bases inherit the
adjacent V/D/J segment's germline lineage).

This is the heaviest of the deferred features and the audit
recommends NOT taking it in v1.

### TdT context dependency / per-context Markov chains

**[Shipped]** The Markov NP base generator ships in
[`docs/np_markov_base_generator_design.md`](np_markov_base_generator_design.md).
Cartridges author per-NP-region Markov sampling via
`NpBaseModelSpec(kind="markov", first_base=…, transitions=…)`;
the engine's new `NpBaseGenerator` trait carries the per-
position support and the `GenerateNPPass` loop threads a
`previous: Option<u8>` through `sample_base` and the replay
validator. Pre-existing trace addresses
(`np.np{1,2}.bases[i]`) are unchanged — `previous` is
reconstructed from the prior recorded base, so legacy traces
keep replaying byte-identically.

**Legacy `DataConfig.NP_transitions` / `NP_first_bases`
remain orphan.** The bundled cartridges' dicts are still
dead code; auto-lift requires a separate cartridge-
migration slice. Authors who want Markov today populate
`ReferenceEmpiricalModels.np_bases` explicitly.

### Chain / locus-specific N distributions

**Not present.** The bundled `HUMAN_IGH_OGRDB` cartridge has
`NP_transitions["NP1"]` and `NP_transitions["NP2"]` distinct
from `HUMAN_IGK_OGRDB`'s, and presumably from MOUSE_IGH_IMGT,
COW_IGH_IMGT, etc. But because the bridge ignores
`NP_transitions` entirely, every cartridge effectively emits
4-way-uniform NP bases today, losing the locus-specific
biology that the bundled pickles encode.

The audit recommends the new typed model expose this
per-cartridge variation.

### Base-transition Markov behaviour status today

`NP_transitions` is present on the Python `DataConfig` (line
120 of `data_config.py`, `default_factory=dict`). It is
populated on bundled cartridges (HUMAN_IGH_OGRDB carries
non-empty `dict[str, dict[str, float]]` for both `NP1` and
`NP2`). It is **never read by the simulation pipeline** —
the `_dataconfig_extract.py` extraction never references it,
`_compile.py` never reads it, the PyO3 bridge has no kwarg
for it. The MCP helper (`mcp_helpers.py:676-678`) reads it
for display purposes only.

`NP_first_bases` follows the same pattern: present, populated,
ignored.

### Pinned

- `pin_scaffold_np_transitions_field_populated_on_bundled_cartridges`
- `pin_scaffold_np_transitions_listed_as_orphan_field`
- `pin_absence_no_p_nucleotide_pass_or_address`
- `pin_absence_no_locus_specific_base_distribution_today`

---

## 3. Q3 — Cartridge ownership (the central audit question)

### What's typed in `ReferenceEmpiricalModels` today

- `np_lengths: Dict[str, EmpiricalDistributionSpec]` — typed,
  validated, surfaces in manifest as `np_length_keys`.
- `trims: Dict[str, EmpiricalDistributionSpec]` — same shape
  for V/D/J trim distributions.

### What's still legacy-dict-only

- `NP_transitions: Dict[str, Any]` (line 120 of `data_config.py`)
  — orphan; not consumed.
- `NP_first_bases: Dict[str, Any]` (line 121) — orphan; not
  consumed.
- `NP_lengths: Dict[str, Any]` — legacy fallback for
  `reference_models.np_lengths` (still consumed in absence
  of typed model).
- `trim_dicts: Dict[str, Any]` — legacy fallback for
  `reference_models.trims` (still consumed).
- The four other orphan fields: `gene_use_dict`,
  `correction_maps`, `asc_tables`,
  `p_nucleotide_length_probs`, `dj_pairing_map`.

### The stop-condition check

The audit needed to confirm that `NP_transitions` and
`NP_first_bases` do not affect simulation output without
being in the plan signature. **Confirmed: they don't affect
output at all today**, because the bridge never reads them.
No replay-safety hole. The audit may proceed.

### Recommended typed shape

For the implementation slice (deferred), the audit recommends:

```python
@dataclass(frozen=True)
class NpBaseModelSpec:
    """Cartridge-owned NP base sampling model. v1 supports two
    shapes:

    - ``uniform`` (default) — the current behaviour; emits 4-way
      equal A/C/G/T.
    - ``empirical_first_base`` — per-NP-region first-base
      categorical distribution (consumes `NP_first_bases` from
      legacy cartridges).
    - ``markov`` — per-NP-region 4×4 transition matrix on the
      previous emitted base (consumes `NP_transitions`).
    """
    kind: str  # "uniform" | "empirical_first_base" | "markov"
    first_bases: Optional[Dict[str, float]] = None  # per NP region
    transitions: Optional[Dict[str, Dict[str, float]]] = None
```

Stored on `ReferenceEmpiricalModels` as:

```python
np_bases: Dict[str, NpBaseModelSpec]  # keys: "NP1" / "NP2"
```

The legacy `DataConfig.NP_transitions` and
`DataConfig.NP_first_bases` fields are consumed by the
resolver to populate the typed spec when the cartridge
doesn't already provide one (mirrors the
`np_lengths` legacy-fallback pattern). After the slice
ships, both fields are removed from
`_DOCUMENTED_ORPHAN_DATACONFIG_FIELDS`.

### Pinned

- `pin_scaffold_reference_models_carries_typed_np_lengths`
- `pin_absence_no_typed_np_base_model_spec`
- `pin_absence_no_np_bases_field_on_reference_models`
- `pin_scaffold_documented_orphan_list_carries_np_transitions_and_first_bases`

---

## 4. Q4 — Replay and trace

### Current addresses

- `np.np1.length` (single `Int` choice per simulation).
- `np.np1.bases[i]` for each `i in 0..length` (one `Base` choice
  per emitted base).
- `np.np2.length` + `np.np2.bases[i]` for VDJ chains.

Length is recorded at `execution.rs:57-58`; each base is
recorded at `execution.rs:119-120` inside the per-base loop.

### Replay exactness

Same-seed replay reproduces lengths AND bases byte-for-byte.
The trace cursor consumes the length value, then walks the
per-base loop and consumes each recorded base in order. Under
`UniformBase`, the `parameter_signature` is constant across
all runs and the signature gate passes silently. The
existing test `generate_np_pass_is_deterministic_under_same_seed`
([`engine_rs/src/passes/generate_np.rs:296-319`](../engine_rs/src/passes/generate_np.rs#L296-L319))
pins this property.

### Whether Markov dependencies would be represented in trace

The current per-base recording (`np.np1.bases[i]`) records
the *emitted* base, not the conditioning state — that's
correct for any model the audit considers, because the
conditioning state is recoverable from the prior emitted base
in the same `bases[…]` sequence. A future Markov base model
needs no new trace address; the previous-base lookup happens
at sampling time and the resulting weighted-categorical draw
records the same `Base` choice value.

### Plan-signature folding (the replay-safety guarantee)

`GenerateNPPass::parameter_signature` folds both
`length_dist` and `base_dist` (Slice A discipline). A typed
NP base distribution would have a distinct `support()` shape
from `UniformBase`'s 4-way uniform — the plan signature would
fold those weights in automatically, and a same-cartridge
replay against a different cartridge's NP base model would
fail the plan-signature gate before consuming any choices.
This is exactly the behaviour the audit needs.

### Pinned

- `pin_scaffold_np_trace_addresses_exist_and_replay_byte_identical`
- `pin_scaffold_generate_np_parameter_signature_folds_base_dist`

---

## 5. Q5 — Productive-only constraints

### NP length narrowing under `productive_only()`

The length distribution is NOT narrowed today — `productive_only`
operates on per-position base admissibility, not on the
length draw itself. A length that produces a junction longer
than necessary will still be drawn; the per-base narrowing
afterward avoids stop codons.

### NP base narrowing

`JunctionStopState::build` is called at the start of the pass
when contracts are active. It computes the admit mask per
position based on:

- Distance to the junction anchor (V Cys position + 3 for VDJ
  bands).
- Frame requirements (the junction must remain multiple of 3
  to preserve the productive frame).
- No-stop requirement (no AA codon containing TAA / TAG / TGA
  at any frame position).
- Anchor-preservation requirement (the V Cys + J W/F codons
  stay in their respective frame slots).

The mask is consulted at every per-base sample call via the
fast path `sample_base_with_admit_mask`
([`sampling.rs:166-173`](../engine_rs/src/passes/generate_np/sampling.rs#L166-L173)).
If the canonical 4-way support narrows to zero (i.e. every
A/C/G/T at this position would create a stop / break frame),
the slow path `sample_filtered_with_policy` triggers the
empty-support sentinel.

### Empty-support behaviour

`NP_BASE_EMPTY_SUPPORT = EmptySupport::Sentinel(b'N')`
([`sampling.rs:20`](../engine_rs/src/passes/generate_np/sampling.rs#L20)).

- **Strict mode:** raises `PassError::ConstraintSampling`.
- **Permissive mode:** emits `b'N'` as the sentinel base.

Replay validation accepts `b'N'` only when the admit-mask was
empty (`validate_replayed_np_base` at
[`sampling.rs:221-313`](../engine_rs/src/passes/generate_np/sampling.rs#L221-L313)).
A recorded `b'N'` against a non-empty admit-mask would be a
replay validator failure — same shape as the segment-rate
slice's replay validation.

### Composition with a typed NP base model

A non-uniform base distribution composes cleanly with the
admit-mask path: the mask filters the categorical
distribution's support to admit-only bases, then samples
weighted from the survivors. If the mask × distribution
intersection is empty, the same `EmptySupport::Sentinel(b'N')`
fires. No new edge case introduced.

### Pinned

- `pin_scaffold_productive_only_preserves_triad_under_np_generation`
- `pin_scaffold_np_base_empty_support_sentinel_is_capital_n`

---

## 6. Q6 — AIRR projection

### Today's NP-related AIRR fields

[`engine_rs/src/airr_record/record.rs:77-82`](../engine_rs/src/airr_record/record.rs#L77-L82):

```rust
pub np1: String,        // unclaimed NP1 bases
pub np1_aa: String,
pub np1_length: i64,    // unclaimed NP1 base count
pub np2: String,
pub np2_aa: String,
pub np2_length: i64,
```

`np1_length` and `np2_length` are first-class on every AIRR
record, populated at projection time from `unclaimed_np_string()`
([`builder.rs:198-205`](../engine_rs/src/airr_record/builder.rs#L198-L205)).
**"Unclaimed" semantics:** the realised NP length after V/D/J
live-call extensions reabsorb adjacent NP positions. This is
the right semantics for an immunology-aware consumer — what
matters is how many NP-attributed positions are visible in
the final AIRR record, not how many bases the
`GenerateNPPass` originally emitted.

A downstream consumer doesn't need to derive NP lengths from
coordinate differences (`d_sequence_start - v_sequence_end`
etc.). The dedicated fields are correct and authoritative.

### NP region surface in junction

`junction` (the AIRR field) carries the V-end + NP1 + D + NP2 +
J-start concatenation; the NP regions are interpolated into
the junction by position, not labelled separately. A consumer
can recover NP boundaries from `junction[v_anchor_in_junction
+ 1 : d_anchor_in_junction]` and similar coord arithmetic, but
the canonical fields are `np1` / `np2` / `np1_length` /
`np2_length`.

### No new AIRR fields needed for the recommended slice

The typed NP base model doesn't change the projected NP fields.
The new biology surfaces in *which bases* appear in `np1` /
`np2`, not in any new field.

### Pinned

- `pin_scaffold_np1_length_and_np2_length_present_on_airr_records`
- `pin_scaffold_np1_and_np2_string_fields_present_on_airr_records`

---

## 7. Q7 — Recommended next implementation slice

### Slice 1 (recommended) — Typed NP base model in `ReferenceEmpiricalModels`

Scope:

1. **Python typed spec** — add `NpBaseModelSpec` to
   `reference_models.py` with three kinds: `uniform` (default),
   `empirical_first_base`, `markov`. Add `np_bases:
   Dict[str, NpBaseModelSpec]` to `ReferenceEmpiricalModels`.

2. **Legacy bridge** — `_dataconfig_extract.py` reads
   `cfg.NP_first_bases` / `cfg.NP_transitions` as the legacy
   fallback when `reference_models.np_bases` is absent —
   mirror of the `np_lengths` resolution chain.

3. **PyO3 bridge** — extend `push_generate_np` with an
   optional `base_distribution: Option<Vec<(u8, f64)>>` for
   the first-base case OR with a separate
   `push_generate_np_with_markov(transitions: …)` factory for
   the Markov case. The audit recommends an enum tag pattern
   to keep the API surface narrow.

4. **Rust distribution types** — add `CategoricalBase` (a
   weighted 4-way categorical over A/C/G/T) for the
   first-base model. Add `MarkovBase` (a transition matrix
   over the previous emitted base) for the Markov model.

5. **`GenerateNPPass`** — accept a `Box<dyn Distribution>` for
   `base_dist`; the existing constructor signature handles
   this verbatim. No struct shape change.

6. **Plan signature** — `fmt_byte_dist` already folds the
   distribution's `support()` weights; both new distributions
   produce distinct signatures from `UniformBase`. No
   `paramsig.rs` changes needed.

7. **Manifest** — extend the empirical-models block with
   `np_base_keys: List[str]` paralleling `np_length_keys`.
   Remove `NP_transitions` / `NP_first_bases` from
   `_DOCUMENTED_ORPHAN_DATACONFIG_FIELDS`.

8. **Productive-only composition** — no Rust changes; the
   admit-mask path already composes with any
   `Box<dyn Distribution>` because it filters by `support()`
   before sampling.

9. **Tests** — per-base distribution selection composes with
   `productive_only()`; the typed spec validates at load
   time; replay reproduces non-uniform NP bases exactly; the
   manifest reports the new capability; legacy cartridges
   round-trip through the new bridge unchanged (UniformBase
   stays the default for cartridges without typed NP base
   models).

**Cost estimate:** ~150 lines Python (spec + extract + manifest)
+ ~80 lines Rust (CategoricalBase / MarkovBase + bridge
hookup) + ~30 lines tests.

### Why first-base + Markov, not P-nucleotides

P-nucleotide modeling (Slice 2) is much larger:

- Requires a new event payload (or address namespace) to
  distinguish P-additions from N-additions at projection
  time.
- The biology spec is more contested — different sources
  disagree on whether the P-nucleotide adjacent to V is the
  RC of the V-end base or the WC of the V-end base, on
  whether VJ chains have V-side P-nucleotides at all, etc.
- The AIRR consumer experience needs design work — does
  `np1` carry both N and P bases (current shape) or split
  into a new `p_nucleotides` field?

The audit defers P-nucleotides to a separate slice that gets
its own audit doc.

### Pinned

- `pin_absence_no_typed_np_base_model_spec`
- `pin_absence_no_np_bases_field_on_reference_models`
- `pin_absence_no_np_base_keys_in_manifest`

---

## 8. Edge cases the slice must handle

| Case | Expected behaviour |
|---|---|
| Cartridge with no typed `np_bases` and no legacy `NP_transitions` | `UniformBase` continues to apply — byte-identical to today's behaviour. |
| Cartridge with legacy `NP_transitions` only | Resolver lifts to typed `NpBaseModelSpec(kind="markov", transitions=…)`. |
| Cartridge with legacy `NP_first_bases` only | Resolver lifts to typed `NpBaseModelSpec(kind="empirical_first_base", first_bases=…)`. |
| Cartridge with BOTH `NP_transitions` AND `NP_first_bases` | Resolver prefers `transitions` (the Markov model is more biology-rich). Surface a manifest note documenting the choice. |
| Cartridge with typed `np_bases` AND legacy fields | Typed wins; legacy ignored (mirrors `np_lengths` precedence). |
| User passes a `np_base_model=` kwarg on `.recombine(...)` | Per-experiment override; mirrors the existing `np1_lengths` kwarg pattern. |
| `productive_only()` admit mask zeroes the categorical's support | Same empty-support sentinel as today (`b'N'` permissive / raise strict). |
| Replay against a cartridge with different `np_bases` | Plan-signature gate fires the mismatch BEFORE consuming any choice. |
| Markov model produces a zero-weight row (no valid `to_base` from a given `from_base`) | Treat as malformed at load time — raise `ValidationError` from `NpBaseModelSpec` constructor. |
| First-base model with weights that sum to zero | Same — load-time validation rejects. |

---

## 9. Performance

The categorical / Markov base sampling adds at most one
4×4 weight lookup per emitted NP base (the Markov case looks
up the row for the previous base, then samples). Per
simulation: typically 0–6 NP1 bases + 0–6 NP2 bases — under
20 weight lookups total. Negligible.

`UniformBase` stays the fast path for cartridges without typed
NP base models. The `parameter_signature` short-circuit
discipline from Slice A means `UniformBase` continues to
serialise to the same canonical 4-way support string — no
regression in plan-signature generation.

---

## 10. Trace / replay impact

### No new trace addresses

The slice reuses `np.np1.bases[i]` / `np.np2.bases[i]`
verbatim. The Markov / categorical sampling happens at the
existing per-base call site; the recorded base is whatever
the distribution emitted.

### Replay determinism

Same-seed + same cartridge replay reproduces byte-identical
NP bases — the trace records the emitted base, not the
conditioning state. The next-base draw consults the
deterministic prior base (already in the recorded trace), so
the Markov conditioning state is recoverable without recording
it separately.

### Plan-signature folding

Already in place. A new `CategoricalBase([(b'A', 0.3),
(b'C', 0.2), (b'G', 0.2), (b'T', 0.3)])` produces a different
`support()` string from `UniformBase`, so the plan signature
folds the per-cartridge variation automatically (Slice A
discipline).

### Pinned

- `pin_scaffold_np_trace_addresses_exist_and_replay_byte_identical`
  (carry).

---

## 11. Manifest extension

The slice extends the empirical-models block:

```python
manifest["models"]["empirical"] = {
    "has_reference_models": True,
    "np_length_keys": ["NP1", "NP2"],
    "np_base_keys": ["NP1", "NP2"],   # NEW
    "np_base_model_kind": {           # NEW
        "NP1": "markov",
        "NP2": "markov",
    },
    "trim_keys": [...],
    "legacy_np_lengths_present": False,
    "legacy_trim_dicts_present": False,
    "legacy_np_transitions_present": False,   # NEW
    "legacy_np_first_bases_present": False,   # NEW
}
```

The new keys advertise the per-NP-region model kind so a
downstream consumer can detect which biology the cartridge
encodes (uniform vs first-base vs Markov). Two new
`legacy_…_present` flags follow the existing convention so
a cartridge built from legacy dicts is still detectable in
the manifest.

### Pinned

- `pin_absence_no_np_base_keys_in_manifest`
- `pin_absence_no_legacy_np_transitions_present_flag_in_manifest`

---

## 12. Implementation order **[Shipped]**

A single self-contained slice. Six sub-steps, all landed except
where noted:

1. **Typed spec** — `NpBaseModelSpec` + `np_bases:
   Dict[str, NpBaseModelSpec]` on `ReferenceEmpiricalModels`.
2. **Resolver fallback** — `_dataconfig_extract.py` lifts
   legacy `NP_transitions` / `NP_first_bases` when the typed
   field is absent.
3. **PyO3 bridge** — extend `push_generate_np` to thread the
   typed base-distribution data through.
4. **Rust distributions** — `CategoricalBase` and `MarkovBase`
   types; both implement `Distribution<Output = u8>`.
5. **Manifest extension** — new keys per §11.
6. **Tests** — per-cartridge byte-distribution check
   (KL-divergence vs uniform on a high-N batch), replay
   round-trip, plan-signature gate, productive-only
   composition, manifest exposure.

Cost estimate: ~250 lines total, all Python + Rust (no MCP /
docs work in v1 — the slice's docs change is the audit
file's "Slice 1 shipped" flip).

### Why single slice

Same argument as the V-subregion counters slice: no internal
phase boundary. The typed spec, resolver, bridge, Rust
distributions, manifest, and tests all land together. No
precondition slice needed.

---

## 13. Test surface — what this audit pins

Mirrored in
[`tests/test_junction_n_addition_contract.py`](../tests/test_junction_n_addition_contract.py).

### `pin_scaffold_*` — pre-existing surfaces the slice builds on

1. `GenerateNPPass` has the documented `(np_segment,
   length_dist, base_dist)` struct shape.
2. `UniformBase` is the hardcoded base distribution at the
   PyO3 bridge today.
3. `ReferenceEmpiricalModels.np_lengths` exists as typed
   plane.
4. Legacy fallback chain `typed → legacy dict → uniform
   placeholder` exists in `_dataconfig_extract`.
5. Trace addresses `np.np1.length` / `np.np1.bases[i]` (and
   NP2 equivalents) exist and replay byte-identical.
6. `GenerateNPPass::parameter_signature` folds both
   `length_dist` and `base_dist` (Slice A discipline).
7. Productive-only narrowing via `JunctionStopState` works
   today and the `EmptySupport::Sentinel(b'N')` policy is
   pinned.
8. AIRR fields `np1`, `np1_length`, `np2`, `np2_length` are
   present on every record.
9. `_DOCUMENTED_ORPHAN_DATACONFIG_FIELDS` lists both
   `NP_transitions` and `NP_first_bases`.
10. Bundled `HUMAN_IGH_OGRDB` carries non-empty
    `NP_transitions` and `NP_first_bases` (the data is
    present and recoverable; the slice consumes it).
11. VJ chains run `np1` only; VDJ chains run both `np1` and
    `np2`.

### `pin_absence_*` — gaps the slice closes

12. No `NpBaseModelSpec` typed dataclass.
13. No `np_bases` field on `ReferenceEmpiricalModels`.
14. No `np_base_keys` / `np_base_model_kind` in the
    cartridge manifest.
15. No DSL `recombine(np_base_model=…)` kwarg.
16. No P-nucleotide DSL / pass / address (deferred entirely
    per §13.1 — not closed by this slice).

### Doc anchor

17. The audit doc exists and references the contract file;
    structure intact.

---

## §Markov-deferred — stop-and-report explanation

**Why `kind="markov"` raises at lowering time, not at validation.**

The user's brief on the implementation slice contained a
gating check:

> Important check: if `GenerateNPPass` samples each base
> independently through a stateless distribution, true Markov
> dependency cannot be represented with only
> `Distribution<Output = u8>`. If that is the case, stop and
> report before faking Markov.

The check fires positively:

- The Rust `Distribution` trait
  ([`engine_rs/src/dist/mod.rs:37-69`](../engine_rs/src/dist/mod.rs#L37-L69))
  exposes `fn sample(&self, rng: &mut Rng) -> Self::Output` —
  **no context parameter.** There is no place to thread a
  "previous base" into the sample call without changing the
  trait, which is invasive (every distribution in the engine
  would need to opt in or be wrapped).
- The brief's alternative — "Markov needs a pass-level change,
  the correct design is `NpBaseGenerator::sample_base(rng,
  previous_base, …)`" — is a real but bigger change. It would
  need: a new trait, a `GenerateNPPass.base_generator` field
  alternative, integration with the existing
  `JunctionStopState` admit-mask narrowing
  ([`engine_rs/src/passes/generate_np/sampling.rs:166-313`](../engine_rs/src/passes/generate_np/sampling.rs#L166-L313)),
  and integration with the replay-time admit-mask validator.

The v1 slice ships `uniform` + `empirical_first_base` only.
For `kind="markov"`:

1. The Python `NpBaseModelSpec.__post_init__` accepts and
   validates the spec (cartridges can author it; partial
   matrices are still rejected; weights are still validated).
   This preserves authoring fidelity for future engines that
   support the model.
2. The resolver
   ([`src/GenAIRR/_dataconfig_extract.py::_np_bases_from_models`](../src/GenAIRR/_dataconfig_extract.py))
   raises `NotImplementedError` with a message that:
   - Names the exact `np_bases[key]` that failed.
   - Explains why the trait doesn't support it.
   - Points at this doc section.
   - Suggests switching to `empirical_first_base` for marginal
     categorical sampling without faking Markov.

The audit's "fake Markov" failure mode the brief named — silent
degradation to row-marginal sampling — is the explicit
antipattern the stop-and-report avoids. A user who configures
`kind="markov"` either fixes their spec to the supported kind
or waits for the follow-up slice; the engine never silently
produces output that doesn't match the configured model.

The follow-up Markov slice will land:

- A new `NpBaseGenerator` trait with `sample_base(rng,
  previous_base: Option<u8>, …)` — `previous_base = None`
  means "first position, sample from `first_base`".
- A second field on `GenerateNPPass` (or a constructor
  overload) accepting the generator.
- Integration with the admit-mask narrowing (the mask filters
  the post-conditioning categorical at each position).
- Replay-time admit-mask validation extended to consult the
  same previous-base state.
- Plan signature folds the transition matrix's row supports
  alongside the first-base support (Slice A discipline
  composes).

That work is its own contract + audit.

---

## 14. Out of scope

Documented here so a future implementer doesn't accidentally
expand the work.

- **P-nucleotides / palindromic additions.** Separate audit;
  needs new event payload / region semantics / AIRR field
  design.
- **TdT context dependency beyond 1-step Markov.** A 2- or
  3-step lookback Markov model would require a different
  storage shape and more careful empirical calibration; v1
  caps at 1-step.
- **VJ-chain-specific N base distributions** (different
  empirical models for VJ vs VDJ NP1). Out of scope; the
  recommended slice uses the cartridge-level `np_bases` field
  uniformly.
- **Length-conditional base distributions** (e.g. "if length
  is 0 the model emits nothing; if length is 5 the first
  base is biased toward G"). Out of scope.
- **DJ pairing** (the existing `dj_pairing_map` orphan
  field). Separate biology slice.
- **Gene-use weights** (`gene_use_dict` orphan field).
  Separate biology slice.
- **AIRR `np1_length` semantics change.** The "unclaimed
  bases" semantics stays; the slice doesn't redefine the
  field.

---

## 15. Summary table

| Concern | Post-slice state | Status |
|---|---|---|
| NP length distribution | Typed `reference_models.np_lengths` (winning) → legacy `NP_lengths` dict (fallback) → uniform placeholder | **Unchanged** |
| NP base distribution | `NpBaseModelSpec` accepted with three kinds; `uniform` + `empirical_first_base` lower to `UniformBase` / `CategoricalBase` at the bridge | **Shipped** (uniform + empirical_first_base); `markov` deferred |
| Cartridge ownership of base model | `ReferenceEmpiricalModels.np_bases: Dict[str, NpBaseModelSpec]` | **Shipped** |
| Legacy `NP_transitions` consumption | Still dead code; legacy auto-lift deferred (would silently change output bytes; brief's stop-and-report) | **Deliberately not auto-lifted** |
| Legacy `NP_first_bases` consumption | Still dead code; same rationale | **Deliberately not auto-lifted** |
| Chain / locus-specific NP base biology | Surfaces per-cartridge via the typed spec when authored | **Shipped (opt-in)** |
| Trace addresses for NP | Unchanged — `np.np1.length` / `np.np1.bases[i]` / NP2 mirror | **Unchanged** |
| Replay determinism | Byte-identical same-seed reproduction across all kinds (uniform + empirical_first_base) | **Shipped** |
| Plan-signature folding | `CategoricalBase.support()` automatically participates via `fmt_byte_dist` (Slice A discipline) | **Shipped** |
| Productive-only NP base narrowing | `JunctionStopState` admit mask composes with `CategoricalBase` unchanged | **Shipped, no admit-path changes** |
| AIRR projection | `np1`, `np1_length`, `np2`, `np2_length` first-class | **Unchanged** |
| Manifest block | `np_base_models` block under `models` with per-region kinds + `supported_kinds` / `deferred_kinds` / `legacy_fallback` flag | **Shipped** |
| P-nucleotides | Not modelled | **Deferred — separate audit** |
| Markov sampling | Spec validates; lowering raises `NotImplementedError` pointing at this doc | **Stop-and-report — deferred** |
| Pre-flight bugs found | **None.** | — |

The Slice 1 N-addition surface shipped clean. The two
deliberate scope cuts (legacy auto-lift, Markov) are each
documented as stop-and-report decisions: the former because
auto-lifting would silently change output bytes vs the
pre-slice baseline, the latter because the engine's
`Distribution<Output = u8>` trait cannot express previous-base
conditioning without a pass-level architectural change. Both
are tracked as named follow-up slices.

The natural next step is the **Markov NP base model slice** —
ship the `NpBaseGenerator` trait pass-level extension so
`kind="markov"` lowers to real previous-base-conditional
sampling, removing the stop-and-report on
`NotImplementedError`. P-nucleotide modelling stays out of
scope across this whole track and remains in its own future
audit.
