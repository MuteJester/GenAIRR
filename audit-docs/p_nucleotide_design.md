# P-Nucleotide / Palindromic Addition Modeling — Audit + v1 Implementation Shipped

**Status: v1 implementation shipped** — closes the remaining
material junction-biology gap after the N-addition story.
The engine ships `PAdditionPass`, `PEnd`, `PRegionAdded` event,
four `p.*.length` trace addresses, four `p_*_length` AIRR
fields, and a `PLengthMismatch` validator issue kind. The
cartridge gains a typed `ReferenceEmpiricalModels.p_nucleotide_lengths`
plane and a `p_nucleotide_models` manifest block. Legacy
`p_nucleotide_length_probs` auto-lift remains explicitly
deferred (audit's stop-and-report boundary preserved).

Companion artefacts:

- [`tests/test_p_nucleotide_contract.py`](../tests/test_p_nucleotide_contract.py)
  — 21 pins, all flipped to post-implementation state.
- [`tests/test_p_nucleotide_implementation.py`](../tests/test_p_nucleotide_implementation.py)
  — 10 behaviour tests covering default-byte-identical, VJ +
  VDJ P emission, `P_NUC` flag, replay round-trip, productive-
  only triad, tampered-field validation, manifest exposure,
  legacy-orphan boundary, signature folding.

**Audit body preserved below for traceability.** The
v1-vs-out-of-scope split documented at §15 holds verbatim;
pre-trim mode, per-base AIRR strings, aggregate
`n_p_nucleotides`, and legacy auto-lift all remain deferred.

P-nucleotides are **not just another NP base distribution.**
They are templated palindromic extensions of the coding ends
created when the V(D)J hairpin is opened off-centre during
recombination. That difference drives every recommendation in
this audit:

- Lengths are sampled per end (`V_3`, `D_5`, `D_3`, `J_5`); bases
  are deterministic from the underlying allele's post-trim
  coding end (reverse-complement of the last *N* coding bases).
- Provenance is biologically distinct from N-bases, so events,
  trace, and AIRR projection should expose P-bases on their
  own surface (separate length fields per end) even though the
  sequence placement is adjacent to NP1 / NP2.
- The plan signature folds P-length distributions per end so
  same-cartridge replay is byte-deterministic; the actual
  P-bases need not be in the trace because they replay
  deterministically from `(allele, trim, length)`.

Companion to
[`tests/test_p_nucleotide_contract.py`](../tests/test_p_nucleotide_contract.py)
which freezes today's surfaces (`pin_scaffold_*`), the existing
absence of any P-nucleotide pass / event / address / field
(`pin_absence_*`), and the documented-orphan status of the
legacy `DataConfig.p_nucleotide_length_probs` (`pin_present_*`).

**Pre-flight finding (Q7 below): clean yes — no stop-and-report
condition.** The legacy `DataConfig.p_nucleotide_length_probs`
field exists on the dataclass with `default_factory=lambda:
dict(DEFAULT_P_NUCLEOTIDE_LENGTH_PROBS)` and is explicitly
listed in `_DOCUMENTED_ORPHAN_DATACONFIG_FIELDS`. The only
consumer in the entire repo is the MCP helper's read-only
diagnostic endpoint
([`src/GenAIRR/utilities/mcp_helpers.py:682-685`](../src/GenAIRR/utilities/mcp_helpers.py#L682-L685))
which surfaces the dict in a `p_nucleotides` inspection
section — it never feeds the simulation pipeline. No pass
reads it, no bridge kwarg threads it, no AIRR field depends on
it. No replay-safety hole exists.

---

## Existing scaffolding the contract pins

| Surface | Where | What it pins |
|---|---|---|
| `Segment` enum + `Segment::COUNT = 5` invariant | [`engine_rs/src/ir/segment.rs:13-25`](../engine_rs/src/ir/segment.rs#L13-L25) | Today: `V / Np1 / D / Np2 / J`. The audit recommends NOT extending this enum with `Pv3`/`Pd5`/`Pd3`/`Pj5` variants — that would break `PerSegment` storage and every existing `match Segment` site. Use a new typed `PEnd` enum + a new `PRegionAdded` event variant instead (§4). |
| `ChoiceAddress` enum | [`engine_rs/src/address.rs:340-392`](../engine_rs/src/address.rs#L340-L392) | Carries `NpLength` / `NpBase` + every other built-in address. The slice adds one new variant per P-end (`PLength { end: PEnd }`); P-bases are deterministic and need no per-base address (§3). |
| `SimulationEvent` enum + `RegionAdded { region }` precedent | [`engine_rs/src/ir/sim_event.rs:105-173`](../engine_rs/src/ir/sim_event.rs#L105-L173) | Existing `RegionAdded { region }` is what NP / V / D / J emission already uses. The slice adds a new `PRegionAdded { end: PEnd, region }` variant so AIRR projection + counters can attribute P-bases without overloading the existing `RegionAdded` surface. The "Reserved for future emission" `RegionReplaced` variant precedent shows the codebase's discipline of distinct variants per biological act. |
| `complement_base(b: u8) -> u8` | [`engine_rs/src/ir/nucleotide.rs:26`](../engine_rs/src/ir/nucleotide.rs#L26) | The palindrome computation primitive. P-extension at the V 3' end = `complement_base` applied to the V's last `length` coding bases in reverse order. Reuse verbatim. |
| `flag::P_NUC` already reserved | [`engine_rs/src/ir/nucleotide.rs:85-87`](../engine_rs/src/ir/nucleotide.rs#L85-L87) | The `P_NUC` `NucFlags` bit (1 << 0) is **already declared** with the exact biological semantics the slice needs: *"templated palindromic complement of the segment edge during V(D)J recombination"*. No production pass sets it today (`grep` confirms only the flag-arithmetic unit tests + a docstring reference in `sample_base.rs` mention it). The slice's job at the flag layer reduces to actually flipping this pre-reserved bit on every emitted P-base — no new flag declaration needed. Same disposition as the `Nucleotide::GermlinePos::NONE` sentinel that already documents "synthetic bases (NP, P-nuc, contaminant, indel-inserted) with no germline provenance". |
| Existing pipeline lowering order | [`src/GenAIRR/_compile.py:220-276`](../src/GenAIRR/_compile.py#L220-L276) | VJ: `sample → trim → assemble.V → np.NP1 → assemble.J`. VDJ: `sample → trim → assemble.V → np.NP1 → [invert_d] → assemble.D → np.NP2 → assemble.J`. P-passes interleave at the assembly boundaries (§1.4). |
| `GenerateNPPass.execute_with_sampling_mode` per-position discipline | [`engine_rs/src/passes/generate_np/execution.rs`](../engine_rs/src/passes/generate_np/execution.rs) | Precedent for per-position emission with admit-mask + JunctionStopState composition. P-addition's per-end length filter follows the same pattern. |
| `JunctionStopState::build` — single-shot precompute | [`engine_rs/src/contract/junction_stop_state.rs`](../engine_rs/src/contract/junction_stop_state.rs) | Used today by NP base sampling to filter candidates against the productive triad. P-addition needs the same filter applied to the LENGTH distribution (not per-base): a P-length that would produce a junction-disrupting frame or stop codon is rejected before sampling. |
| `_DOCUMENTED_ORPHAN_DATACONFIG_FIELDS` | [`src/GenAIRR/dataconfig/data_config.py:43-51`](../src/GenAIRR/dataconfig/data_config.py#L43-L51) | Explicit list of orphan legacy fields. `p_nucleotide_length_probs` is in this list. The slice does NOT auto-lift this orphan; cartridges author the typed P-length plane explicitly (§7). |
| MCP `p_nucleotides` diagnostic endpoint | [`src/GenAIRR/utilities/mcp_helpers.py:682-688`](../src/GenAIRR/utilities/mcp_helpers.py#L682-L688) | Read-only surface for `getattr(dc, "p_nucleotide_length_probs", {})`. Does NOT feed the simulator. Stays untouched by the slice — it surfaces the legacy field for cartridge inspection, not for sampling. |
| Existing absence pin in junction audit | [`tests/test_junction_n_addition_contract.py::test_pin_absence_no_p_nucleotide_surface`](../tests/test_junction_n_addition_contract.py) | Already pins that no `p_nucleotides` kwarg / method / AIRR field exists. The new contract file `test_p_nucleotide_contract.py` is a per-slice elaboration of this pin set. |

---

## 1. Q1 — Biological placement

### 1.1 Hairpin vs trim — what the biology says

In real V(D)J recombination, RAG1/2 cuts the coding flank to
form a hairpin loop; Artemis opens the hairpin asymmetrically,
which can extend the coding end by 0–4 palindromic
(P-)nucleotides. **Exonucleases then trim the extended end**,
potentially removing both P-bases and original coding bases.
TdT-mediated N-addition fills the resulting gap.

True biological order:

```text
sample allele  →  hairpin open (adds 0–N P-bases)  →  trim (removes 0–M bases)  →  N-addition
```

### 1.2 The trade-off

Honoring true biology in the engine means the trim distribution
operates on a **post-P** coding end. But the cartridges' trim
distributions in [`ReferenceEmpiricalModels.trims`](../src/GenAIRR/reference_models.py)
were empirically calibrated against observed coding-end loss in
sequencing data — which is the *net* loss after hairpin → trim
in the underlying biology. The cartridges' trim distributions
are therefore already a convolution of (`P_length`, `trim_back`)
realisations. Inserting P-addition *before* trim with those same
distributions would double-count the P-extension and produce
systematically more 3'-extended coding ends than the underlying
biology produced.

### 1.3 Recommendation: P-addition AFTER trim

**v1 places P-addition AFTER trim**, not before. Trade-offs
documented openly:

| Consideration | Pre-trim P (true biology) | Post-trim P (audit recommendation) |
|---|---|---|
| Biological faithfulness | Higher — matches RAG/Artemis sequence | Lower — operational simplification |
| Reuses existing trim distributions byte-identically | **No** — distributions would need recalibration | **Yes** — trim distributions stay valid |
| Marginal coding-end-loss distribution | Convolution of P-length × trim distributions | Trim distribution alone, then additive P |
| Replay determinism | Same | Same |
| Plan-signature stability for existing cartridges (uniform / empirical_first_base NP, no P) | Identical (no P → no extra factor) | Identical |
| User intuition | "Hairpin first, then trim" | "Trimmed coding end, then extend with P, then N" |

The "post-trim P" model is operationally clean and keeps the
existing trim distributions valid. The mathematical relationship
to true biology is: the engine's modelled P-length distribution
parametrises a *residual* hairpin extension that survives any
trimming — which is exactly the observable P-nucleotide
contribution in sequencing data.

A future enrichment could add a pre-trim P mode under an opt-in
flag with recalibrated trim distributions. v1 ships post-trim.

### 1.4 Ends that receive P-bases

| Chain | Ends | Pipeline position relative to existing passes |
|---|---|---|
| VJ | V 3', J 5' | After `assemble.V` + before `generate_np.NP1` (V 3'); after NP1 + before `assemble.J` (J 5') |
| VDJ | V 3', D 5', D 3', J 5' | After `assemble.V` + before NP1 (V 3'); after NP1 + before `assemble.D` (D 5'); after `assemble.D` + before NP2 (D 3'); after NP2 + before `assemble.J` (J 5') |

The recommended VDJ lowering order (one pass per end):

```text
sample_allele.{V,D,J}
trim.{v_3, d_5, d_3, j_5}
assemble.V
p_addition.v_3          ← NEW
generate_np.NP1
[invert_d]              ← commits D orientation BEFORE p_addition.d_5 reads it
p_addition.d_5          ← NEW (reads post-inversion D effective_seq)
assemble.D
p_addition.d_3          ← NEW
generate_np.NP2
p_addition.j_5          ← NEW
assemble.J
```

**Critical ordering note (corrected from this doc's earlier
draft):** `p_addition.d_5` MUST run AFTER `invert_d` so that
D's `effective_seq` is read under the post-inversion
orientation. An earlier revision of this section showed
`p_addition.d_5 → invert_d → assemble.D`; that ordering would
palindrome the wrong end of D whenever the inversion fires.
The implementation slice corrected the order per §9.3 — see
the
[contract pin](../tests/test_p_nucleotide_contract.py)
`test_pin_present_pipeline_order_has_p_addition_at_audited_positions`
which freezes the corrected order at the lowering source.

For VJ:

```text
sample_allele.{V,J}
trim.{v_3, j_5}
assemble.V
p_addition.v_3          ← NEW
generate_np.NP1
p_addition.j_5          ← NEW
assemble.J
```

### 1.5 Pinned

- `pin_absence_no_p_addition_pass_in_engine`
- `pin_absence_no_p_addition_in_pipeline_lowering`
- `pin_scaffold_pipeline_order_today_has_no_p_addition`

---

## 2. Q2 — Relationship to NP regions

### 2.1 Should P-bases live inside NP1 / NP2 regions, or in their own regions?

**Recommendation: separate `Region(segment=Pv3 | Pd5 | Pd3 |
Pj5)` regions, NOT extensions of NP1 / NP2.**

Two separate concerns:

1. **Sequence placement** — adjacent to NP1 / NP2. The
   P-extension at V 3' sits at pool positions
   `[v.end, v.end + p_v_3.length)`; NP1 then occupies
   `[v.end + p_v_3.length, v.end + p_v_3.length + np1.length)`;
   the D 5' P-extension sits at
   `[v.end + p_v_3.length + np1.length,
     v.end + p_v_3.length + np1.length + p_d_5.length)`;
   then assembled D coding bases; etc.

2. **Provenance** — P-bases and N-bases are biologically
   distinct (templated vs randomly added). Counters,
   plan-signature folding, AIRR projection, and downstream
   tooling want to distinguish them. Merging into NP1 / NP2
   would lose the distinction.

The audit recommends NEW region segments via a typed `PEnd`
enum:

```rust
pub enum PEnd { V3, D5, D3, J5 }
```

…carried in a new `PRegionAdded { end: PEnd, region: Region }`
event variant. This sidesteps the `Segment` enum extension
trap (§Q4).

### 2.2 Why not just stick P-bases inside NP1 / NP2 with a flag?

Considered and rejected:

- **`Nucleotide.flags` could carry `P_NUC` alongside `N_NUC`**
  — but the per-region `Region.segment` would still say `Np1`,
  forcing every consumer (AIRR projection, counters, plan
  signature folder) to walk into the pool and check per-base
  flags to attribute correctly. Cross-cutting.
- **Existing AIRR `np1` / `np2` fields would silently include
  P-bases** — breaking the v1 boundary that `np1.length` is
  the count of N-additions. Backwards-incompatible projection
  change.

Separate regions keep the existing `Region.segment` /
`Nucleotide.segment` discrimination clean.

### 2.3 Sequence placement diagram

For a VDJ record:

```text
pool indices:  [─── V ───]  [P_V3]  [── NP1 ──]  [P_D5]  [─── D ───]  [P_D3]  [── NP2 ──]  [P_J5]  [─── J ───]
biology:        coding       templ.    N-add      templ.    coding      templ.    N-add      templ.    coding
                 (V)         (V tail)  (TdT)     (D head)    (D)        (D tail)  (TdT)     (J head)    (J)
```

Each `P_*` region carries P-bases that are the reverse-
complement of the adjacent coding flank:

- `P_V3`: complement of the LAST `length` bases of V's trimmed
  coding sequence (in reverse). So if V's trimmed end is
  `...XYZ` and `length=3`, the P-extension reads `Z'Y'X'`
  (where `X' = complement_base(X)`).
- `P_D5`: complement of the FIRST `length` bases of D's
  trimmed coding sequence (in reverse). If D's trimmed 5' is
  `ABC...`, the P-extension reads `C'B'A'` and sits just
  upstream of D in the pool.
- `P_D3`: same as `P_V3` but at the D's 3' end.
- `P_J5`: same as `P_D5` but at the J's 5' end.

The palindrome property: each P-base is the complement of the
mirror-position coding base, so reading the boundary forward
gives `...XYZ Z'Y'X'...` which is a (short) inverted repeat —
the molecular signature of the hairpin opening.

### 2.4 Pinned

- `pin_absence_no_p_region_segment_variants_today`
- `pin_absence_no_p_region_added_event_variant_today`

---

## 3. Q3 — Trace model

### 3.1 Sampled values

**Only lengths are sampled.** P-bases are deterministic from
`(allele, trim, length)` — the post-trim allele coding bases
provide the source for `complement_base` + reversal.

### 3.2 Trace addresses

Recommended four new `ChoiceAddress` variants:

```rust
pub enum PEnd { V3, D5, D3, J5 }

pub enum ChoiceAddress {
    // ...existing...
    PLength { end: PEnd },
}
```

On-disk spellings (parallel to the existing `np.np1.length` /
`trim.v_3` convention):

```text
p.v_3.length
p.d_5.length
p.d_3.length
p.j_5.length
```

### 3.3 No per-base addresses

The Markov NP slice introduced per-base addresses
(`np.np1.bases[i]`) because each base was independently
sampled. P-bases are deterministic from the trace's recorded
length + the assembled allele's post-trim coding bases — both
of which are already replayable. **No per-base address is
needed.**

The replay-side validation is symmetric: the validator
recomputes the expected P-bases from `(allele, trim,
recorded_length)` and asserts they match the pool's
`P_*` region byte-for-byte.

### 3.4 Empty-support / length-0 semantics

A P-length distribution that returns 0 means "no P-extension at
this end". The pass executes a no-op (no region added, no
event emitted, no pool mutation). Symmetric with the existing
`NP_LENGTH_EMPTY_SUPPORT = EmptySupport::Sentinel(0)` policy
for NP-length sampling.

### 3.5 Plan-signature folding

Each P-end's length distribution folds into
`PAdditionPass::parameter_signature` via the existing
`fmt_int_dist` helper:

```text
p_addition.v_3(length=[(0:0.5),(1:0.25),(2:0.15),(3:0.07),(4:0.03)])
```

…parallel to the existing `generate_np.np1(length=...)` shape.
Two cartridges with the same NP biology but different P-length
distributions produce different plan signatures, so replay
against the wrong cartridge fails the signature gate.

### 3.6 Pinned

- `pin_absence_no_p_length_choice_address_today`
- `pin_absence_no_p_length_address_strings_today`

---

## 4. Q4 — Events

### 4.1 Recommendation: new event variant, NOT new Segment variants

**Adding `Segment::Pv3` / `Pd5` / `Pd3` / `Pj5` would break
`Segment::COUNT = 5` invariant** which every `PerSegment`
storage site, every `match Segment` site, and the
[`engine_rs/src/ir/per_segment.rs`](../engine_rs/src/ir/per_segment.rs)
slot array depends on. The audit explicitly recommends NOT
taking that route.

Instead: add a typed `PEnd` enum + a new `SimulationEvent`
variant that carries it:

```rust
#[derive(Copy, Clone, Eq, PartialEq, Hash, Debug)]
pub enum PEnd { V3, D5, D3, J5 }

pub enum SimulationEvent {
    // ...existing...
    /// A new P-nucleotide region was appended at the named end.
    PRegionAdded { end: PEnd, region: Region },
    // ...existing...
}
```

`Region.segment` for P-regions: open question. Two options:

| Option | Trade-off |
|---|---|
| **A.** Reuse `Segment::Np1` / `Segment::Np2` for the region's `segment` field; rely on `PRegionAdded.end` for provenance | Pool walks see one continuous "junction" segment range; minimal change to existing AIRR projection code. P-base discrimination requires checking the event ledger or per-base flags. |
| **B.** Add a single `Segment::P` variant + `PEnd` to discriminate at the event layer | Cleanest discrimination at the pool level; breaks `Segment::COUNT = 5` and every match site. |
| **C.** Reuse the adjacent segment (V for `P_V3`, D for `P_D5` + `P_D3`, J for `P_J5`) | Region segment matches biology (P-bases descend from the named segment's coding end); pool walks can treat P + coding as one "extended segment". |

**Audit recommendation: Option C.** P-bases are templated
extensions of the named segment's coding end; carrying the
adjacent segment in `Region.segment` aligns with biology
("these bases are V-derived") and avoids the `Segment` enum
extension. The pre-reserved `Nucleotide::flags::P_NUC` bit
already exists for exactly this — the slice flips it on
every emitted P-base; per-base discrimination is then a
single bit test (no event-ledger walk).

### 4.2 Cross-checks against existing observers

The `DirtySignalObserver` watches `PassEffect::AppendRegion(seg)`.
Reusing the adjacent `Segment` for P-regions means the dirty-
signal accounting is implicit. A new `PassEffect::AppendPRegion {
end }` variant is OPTIONAL — the audit recommends NOT adding it
in v1 to keep the dirty-signal surface unchanged; live-call
parity falls out of the existing AppendRegion-keyed bookkeeping.

### 4.3 Pinned

- `pin_absence_no_p_region_added_event_variant_today`
- `pin_absence_no_p_end_enum_today`
- `pin_present_p_nuc_flag_already_reserved_no_pass_emits_it`

---

## 5. Q5 — AIRR fields

### 5.1 Recommendation: four new int fields

Length-only projection at the AIRR layer. Per-end:

```text
p_v_3_length: int    (VJ + VDJ; 0 when no P extension)
p_d_5_length: int    (VDJ only; absent or 0 in VJ — recommend 0 for uniform record shape)
p_d_3_length: int    (VDJ only)
p_j_5_length: int    (VJ + VDJ)
```

Each field is the count of templated P-bases inserted at the
named end. Sum is the total P-base contribution to the
junction; useful for downstream filtering.

### 5.2 NOT recommended for v1

| Considered | Reason rejected for v1 |
|---|---|
| `p_v_3: str` (the actual P-base sequence) | Deterministic from `(v_call, v_trim_3, p_v_3_length)`; downstream tools can recompute. Adds projection cost + bloat without new information. Additive in v2 if requested. |
| `n_p_nucleotides: int` (aggregate count) | Easy to derive as `p_v_3_length + p_d_5_length + p_d_3_length + p_j_5_length`. Don't pre-aggregate at projection time. |
| Modifying `np1.length` to include P-bases | Backwards-incompatible: existing tools assume `np1.length` is N-additions only. Hard breakage. |
| `junction_length` excluding P-bases | The junction by definition spans the full inter-anchor region; P-bases are part of the junction sequence. Don't carve them out. |

### 5.3 Junction inclusion

The `junction` / `junction_aa` / `junction_length` fields
**include** P-bases (P-bases are part of the realised junction
sequence). The new `p_*_length` fields are provenance overlays
on top of the existing junction projection.

### 5.4 Validator

The AIRR validator gains four invariants (one per end):

```text
PLengthMismatch { end: PEnd, recorded: i64, recomputed: i64 }
```

The recomputation walks the event ledger filtered to
`PRegionAdded { end, region }` and asserts each recorded
`p_*_length` matches `sum(region.len())` for matching ends.
Symmetric with the existing per-segment-mutation /
NP-mutation counter validator surfaces.

### 5.5 Pinned

- `pin_absence_no_p_length_airr_fields_today`
- `pin_absence_no_p_length_mismatch_issue_kinds_today`

---

## 6. Q6 — Productive-only

### 6.1 What P-bases do to the junction

- They change junction LENGTH (each P adds 1 to the inter-
  anchor span).
- They change junction CONTENT (deterministically — the
  reverse-complement of the adjacent allele's coding flank).

Both factors are visible to the existing productive-only
contracts:

| Contract | How it composes |
|---|---|
| `ProductiveJunctionFrame` | Junction length must be `% 3 == 0`. P-bases add to junction length, so the productive admit-set over `(p_v_3, p_d_5, p_d_3, p_j_5, np1, np2)` lengths is the joint set where `sum(all P + NP lengths) + anchor_codon_contribution ≡ 0 (mod 3)`. |
| `NoStopCodonInJunction` | The junction codon sequence must not contain TAA / TAG / TGA. P-bases are deterministic from `(allele, trim)`, so for a given `(allele, trim, P-length)` the contributed codons are determined; the contract filters P-lengths whose contribution would create a stop. |
| `AnchorPreserved` | P-bases sit between V/J anchor codons and the NP regions; they don't displace anchors. Composes trivially. |

### 6.2 Admit-mask analogue

The Markov NP slice's admit-mask machinery filters per-position
bases. P-addition's analogue filters per-length: the
`JunctionStopState`-driven `ProductiveAdmitMaskObserver` can
extend its precompute to cover P-length admissibility before
the P-length distribution is sampled.

Symmetric with how `ProductiveJunctionFrame` already filters NP1
length: the contract precomputes the admissible-length set
once; the sampler intersects the cartridge's L-distribution
with this set, then inverse-CDF-samples.

For P-length sampling the precompute extends to: at each
sampling site (`p.v_3.length`, ..., `p.j_5.length`), the
admit-set is `{ L : the junction stays productive under the
already-committed P + NP lengths and the to-be-deterministic
P-bases for this L }`.

### 6.3 Coupling with NP-length sampling

NP-length sampling today already runs under
`JunctionStopState`. The composition is straightforward: each
P-length sample commits a new effective coding extension, which
the JunctionStopState observer registers, and subsequent NP-
length samples (NP1, NP2) intersect with the updated state.

The sampling ORDER documented in §1.4 (`p.v_3 → np.np1 → p.d_5
→ ... → p.j_5 → assemble.J`) means each contract narrows in
sequence, so no per-pass joint optimisation is needed — the
existing per-pass admit-mask infrastructure carries the
discipline forward.

### 6.4 Edge case: empty admit set

If a particular P-length distribution is entirely incompatible
with the productive triad under the committed state (e.g.
every candidate length would create a stop codon), the
existing `NP_LENGTH_EMPTY_SUPPORT = EmptySupport::Sentinel(0)`
policy generalises to P-length: emit length 0 (no P-extension)
under permissive; surface
`PassError::ConstraintSampling` under strict. Symmetric with
NP-length empty-support policy.

### 6.5 Pinned

- `pin_scaffold_junction_stop_state_already_precomputes_length_admit_set`
- `pin_scaffold_existing_np_length_empty_support_policy_is_sentinel_zero`

---

## 7. Q7 — Cartridge ownership

### 7.1 Recommended typed plane

Add `p_nucleotide_lengths` to
[`ReferenceEmpiricalModels`](../src/GenAIRR/reference_models.py),
keyed by end label (`"V_3"`, `"D_5"`, `"D_3"`, `"J_5"`),
mapping to the existing `EmpiricalDistributionSpec` shape used
by NP-length and trim distributions:

```python
@dataclass(frozen=True)
class ReferenceEmpiricalModels:
    np_lengths: Dict[str, EmpiricalDistributionSpec] = field(default_factory=dict)
    trims: Dict[str, EmpiricalDistributionSpec] = field(default_factory=dict)
    np_bases: Dict[str, NpBaseModelSpec] = field(default_factory=dict)
    p_nucleotide_lengths: Dict[str, EmpiricalDistributionSpec] = field(default_factory=dict)  # NEW
```

- Empty dict (default) → no P-addition at any end → the
  pipeline omits the P-passes entirely → byte-identical to the
  pre-slice behaviour (zero baseline drift).
- Per-end keys: `"V_3"`, `"D_5"`, `"D_3"`, `"J_5"`. Validator
  rejects unknown labels.
- Per-spec: `EmpiricalDistributionSpec` over int values (the
  P-length). Same validation discipline as `np_lengths` /
  `trims`.

### 7.2 VJ / VDJ asymmetry

In VJ, the `"D_5"` and `"D_3"` keys are rejected at validation
time (the chain has no D segment). Symmetric with how
`trims["D_5"]` etc. are rejected for VJ today.

### 7.3 Legacy `p_nucleotide_length_probs` orphan stays orphan

Confirmed orphan today (§Pre-flight). The slice does **NOT**
auto-lift this field. Reasons:

- Legacy `p_nucleotide_length_probs: Dict[int, float]` is a
  **single dict** with no per-end discrimination. The typed
  plane needs per-end resolution (V_3 vs D_5 vs J_5 etc.) for
  biological correctness. Auto-lift would require either
  copying the same dict to all four ends (biologically wrong —
  P-length distributions are end-specific) or arbitrarily
  picking one end (worse).
- Auto-lift would silently change output bytes vs the pre-
  slice baseline, breaking the same backwards-compatibility
  boundary the Markov slice respected for `NP_transitions`.

Authors who want P-nucleotides today populate
`ReferenceEmpiricalModels.p_nucleotide_lengths` explicitly. An
opt-in auto-lift slice is a separate cartridge-migration
decision.

### 7.4 Manifest extension

The `models` block grows a parallel `p_nucleotide_models`
section:

```python
manifest["models"]["p_nucleotide_models"] = {
    "models": {"V_3": {...}, "D_5": {...}, ...},
    "legacy_fallback": False,
    "legacy_p_nucleotide_length_probs_present": True,
    "supported_ends": ["V_3", "D_5", "D_3", "J_5"],
    "deferred_ends": [],
    "in_plan_signature": True,
    "in_content_hash": False,
}
```

`legacy_p_nucleotide_length_probs_present` surfaces the
bundled cartridges' orphan data so a downstream tool can
decide whether to author a typed spec from it. Same shape as
`np_base_models.legacy_np_transitions_present`.

### 7.5 Pinned

- `pin_present_p_nucleotide_length_probs_is_orphan_today`
- `pin_present_p_nucleotide_length_probs_in_documented_orphan_list`
- `pin_absence_no_p_nucleotide_lengths_field_on_reference_models`
- `pin_absence_no_p_nucleotide_models_in_manifest`

---

## 8. Q8 — Replay

### 8.1 Trace records lengths only

Per the §3 trace model: each `PAdditionPass` records exactly
one `ChoiceValue::Int(length)` at `p.{end}.length`. No
per-base records.

### 8.2 Replay-side base derivation

At replay, `PAdditionPass::execute_with_sampling_mode` reads
the recorded length from the trace cursor, walks back through
the assembled allele segment + recorded trim to find the
post-trim coding flank, computes the palindrome via
`complement_base`, and pushes the deterministic P-bases.

The validator gains a "recorded length must equal the actual
P-region length emitted" invariant — symmetric with the NP
base replay validator's "recorded byte must match the
generator's per-position support" pattern.

### 8.3 Determinism guarantee

Same `(allele, trim, P-length)` triple always produces the
same P-base sequence by construction (`complement_base` is a
pure function, reversal is deterministic). No RNG is consumed
on the P-base side; only the length sample consumes one RNG
word per P-pass.

### 8.4 Plan-signature gate

Two cartridges with the same NP biology but different
P-length distributions produce different plan signatures (§3.5).
Replay across them fails the signature gate before any
choice is consumed — same boundary the NP-length and NP-base
slices established.

### 8.5 Pinned

- `pin_absence_no_p_addition_replay_validator_today`

---

## 9. Pipeline order — VDJ + VJ canonical

For ease of cross-referencing during implementation, the
recommended pass sequence per chain:

### 9.1 VJ

```text
sample_allele.v
sample_allele.j
trim.v_3
trim.j_5
assemble.v
p_addition.v_3                ← NEW
generate_np.np1
p_addition.j_5                ← NEW
assemble.j
[mutate / corrupt / paired_end / fastq export — unchanged]
```

### 9.2 VDJ

```text
sample_allele.v
sample_allele.d
sample_allele.j
trim.v_3
trim.d_5
trim.d_3
trim.j_5
assemble.v
p_addition.v_3                ← NEW
generate_np.np1
p_addition.d_5                ← NEW
[invert_d]
assemble.d
p_addition.d_3                ← NEW
generate_np.np2
p_addition.j_5                ← NEW
assemble.j
[mutate / corrupt / paired_end / fastq export — unchanged]
```

### 9.3 Note on `invert_d`

The D-inversion pass flips D's effective orientation. After
inversion, the "5' end" of D's reference becomes the 3' end of
the assembled segment. The `PAdditionPass(end=D_5)` must read
D's effective_seq under the *current* orientation, so the
inversion has to commit before `p_addition.d_5` reads its
source. Today `invert_d` already commits before `assemble.D`;
the recommended pipeline keeps `p_addition.d_5` between
`invert_d` and `assemble.D` so the orientation is settled.

`p_addition.d_3` runs after `assemble.D`, so it sees D's bases
in the pool in the orientation D was assembled in — same
correctness.

### 9.4 Pinned

- `pin_scaffold_pipeline_order_today_has_no_p_addition`
- `pin_scaffold_invert_d_commits_before_assemble_d`

---

## 10. Edge cases the implementation slice must handle

| Case | Expected behaviour |
|---|---|
| `length == 0` | Pass executes a no-op — no region added, no event emitted. Only the trace record (`Int(0)`) is written. Symmetric with NP-length-0. |
| Distribution returns `length > max_p_nuc_length` (4 typically) | Caller bug — panic at execute time. The Python spec validator rejects authoring such distributions. |
| Trimmed allele coding flank shorter than `length` | The palindrome cannot be derived — fewer source bases exist than the requested mirror length. Treat as a malformed `(allele, trim)` pair: under permissive, clamp `length = min(length, available_coding)` and emit; under strict, raise `PassError::ConstraintSampling`. |
| `p_v_3.length + p_d_5.length + np1.length` overflows pool capacity | Pool capacity is bounded by `u32::MAX`; realistic P-lengths (0–4) make overflow impossible under any plausible cartridge. Defensive `u32` arithmetic catches the case. |
| Empty `ReferenceEmpiricalModels.p_nucleotide_lengths` dict | Lowering skips every P-pass — byte-identical to the pre-slice baseline. Default behaviour. |
| Trim distribution removes ALL of the allele's coding sequence | Existing behaviour — `assemble.{V,D,J}` already handles zero-length assembled segments; P-addition at that end has no source flank to palindrome → emit `length=0` even if the cartridge's distribution says otherwise. The contract reduces to "empty admit set". |
| Cartridge with `p_nucleotide_lengths["V_3"]` but no V_3 trim distribution | Fine — P-addition reads the allele's full untrimmed 3' end. |
| `invert_d` interaction | `p_addition.d_5` runs after `invert_d` committed; the D's effective 5' end reflects the inversion outcome. |
| Replay against a trace recorded under no-P cartridge | Plan signature differs (P-length distribution folds in) → signature gate rejects. |
| Replay against same cartridge | All four `p.*.length` records reconstruct deterministically; the validator's per-end count check passes. |

---

## 11. Plan-signature folding

Each P-pass folds its length distribution + end identifier:

```text
p_addition.v_3(length=[(0:0.5),(1:0.25),(2:0.15),(3:0.07),(4:0.03)])
p_addition.d_5(length=[(0:0.4),...])
p_addition.d_3(length=[...])
p_addition.j_5(length=[...])
```

Existing `fmt_int_dist` ([`engine_rs/src/passes/paramsig.rs:108`](../engine_rs/src/passes/paramsig.rs#L108))
handles this with no new helper.

The pass `name()` carries the end suffix (`p_addition.v_3` /
`.d_5` / `.d_3` / `.j_5`) so it appears once per end in the
plan signature, parallel to how `trim.v_3` etc. already render.

---

## 12. Performance

### 12.1 Per-position cost

- Length sample: one RNG word + one `JunctionStopState` admit
  check.
- Per-base derivation: O(length) `complement_base` lookups +
  one push per base. Identical to the existing NP push loop.
- Length distributions have ≤ 5 entries (0–4 bp); each
  inverse-CDF lookup is O(1) effectively.

Per simulation: 4 P-passes (VDJ) or 2 (VJ), each adding at
most 4 pushes → ≤ 16 extra pool operations + 4 trace records +
4 events. Negligible relative to the existing
`assemble_segment` + `generate_np` workload.

### 12.2 Plan-signature cost

4 extra `p_addition.*` substrings (~60 bytes each) → ~240 extra
bytes per plan signature. One-time at compile, negligible.

---

## 13. Implementation order (recommended)

A single self-contained engine slice. Six sub-steps:

1. **Typed cartridge plane** — Python spec:
   `ReferenceEmpiricalModels.p_nucleotide_lengths: Dict[str,
   EmpiricalDistributionSpec]` + validator. Manifest
   `p_nucleotide_models` block.

2. **Engine types** —
   `engine_rs/src/passes/p_addition/mod.rs` with `PEnd`,
   `PAdditionPass`. `SimulationEvent::PRegionAdded { end:
   PEnd, region }` variant.
   `Nucleotide::flags::P_NUC` flag.
   `ChoiceAddress::PLength { end: PEnd }` variant + canonical
   spellings.

3. **Pass execution** —
   `PAdditionPass::execute_with_sampling_mode` mirrors the
   NP pattern: sample length (admit-mask × cartridge
   distribution), derive deterministic P-bases via
   `complement_base` over reversed allele flank, push into
   pool, emit `PRegionAdded`, record `Int(length)` at
   `p.{end}.length`.

4. **PyO3 bridge** — `push_p_addition(end, length_pairs)`
   factory. `_dataconfig_extract._p_nucleotide_lengths_from_models`
   resolver. `_lower_recombine` adds the four P-pass insertions
   at the §9 positions.

5. **AIRR projection + validator** — four new int fields
   (`p_v_3_length` / `p_d_5_length` / `p_d_3_length` /
   `p_j_5_length`); validator surfaces `PLengthMismatch { end,
   recorded, recomputed }` per end.

6. **Replay validator** — recompute expected length from
   event ledger; assert recorded matches.

Cost estimate:
- ~150 lines Rust (pass + event + address)
- ~70 lines Python (spec + resolver + bridge + manifest)
- ~30 lines AIRR field additions + validator issue kinds
- ~50 lines tests (per-end behaviour + replay + productive
  composition)

Single self-contained slice.

---

## 14. Test surface — what this audit pins

Mirrored in
[`tests/test_p_nucleotide_contract.py`](../tests/test_p_nucleotide_contract.py).

### `pin_scaffold_*` — pre-existing surfaces the slice builds on

1. `Segment` enum has 5 variants today; the slice does NOT
   extend it.
2. `complement_base` is the palindrome primitive — reused.
3. `JunctionStopState::build` precomputes the length admit-set
   that NP-length sampling already consumes — extended for
   P-length.
4. `NP_LENGTH_EMPTY_SUPPORT = EmptySupport::Sentinel(0)` —
   length-0 sentinel policy is the precedent for P-length-0
   no-op.
5. Pipeline order today has no P-addition pass.
6. `invert_d` commits before `assemble.D`.
7. The MCP `p_nucleotides` diagnostic endpoint reads the
   legacy field as a read-only inspection — unchanged.

### `pin_present_*` — current state the slice respects

8. `p_nucleotide_length_probs` exists on `DataConfig` with a
   non-empty default factory.
9. `p_nucleotide_length_probs` is in
   `_DOCUMENTED_ORPHAN_DATACONFIG_FIELDS` (orphan policy in
   force).
10. Productive-only NP path runs clean today — pinning the
    pre-P baseline so the slice's composition demonstrably
    extends it.

### `pin_absence_*` — gaps the slice closes

11. No `PAdditionPass` exists in `engine_rs/src/passes/`.
12. No `PEnd` enum exists in the Rust source.
13. No `PRegionAdded` event variant on `SimulationEvent`.
14. No `PLength` variant on `ChoiceAddress`.
15. No `p.v_3.length` / `p.d_5.length` / `p.d_3.length` /
    `p.j_5.length` addresses parse today.
16. No `P_NUC` flag on `Nucleotide.flags`.
17. No `p_nucleotide_lengths` field on
    `ReferenceEmpiricalModels`.
18. No `p_nucleotide_models` block in `cartridge_manifest()`.
19. No `p_v_3_length` / `p_d_5_length` / `p_d_3_length` /
    `p_j_5_length` AIRR fields on the projected record.
20. No `PLengthMismatch` issue kind in the AIRR validator.

### Doc anchor

21. Audit doc exists and references the contract file; section
    structure intact.

---

## 15. Out of scope

Documented here so a future implementer doesn't accidentally
expand the work.

- **Pre-trim P-addition mode** (true biological order). v1
  ships post-trim P; a separate slice could add a pre-trim
  opt-in with recalibrated trim distributions.
- **Per-base AIRR fields** (`p_v_3: str`, etc.). v1 ships
  lengths only; bases are derivable from `(allele, trim,
  length)`. Additive in v2 if requested.
- **Aggregate `n_p_nucleotides` AIRR field.** Easy to derive
  from the four per-end fields; don't pre-aggregate.
- **Legacy `p_nucleotide_length_probs` auto-lift.** Stays
  deferred — same boundary as `NP_transitions` /
  `NP_first_bases` auto-lift.
- **Per-cartridge or per-locus P-base bias** (departures from
  pure palindrome — e.g. asymmetric hairpin opening). True
  biology is mostly palindromic with rare deviations; v1
  models it as deterministic palindrome only.
- **P-bases inside V or J anchor codons.** P-extensions sit
  3' of the V anchor and 5' of the J anchor (junction-
  internal). They don't displace anchors. The `AnchorPreserved`
  contract doesn't need extension.

---

## 16. Summary table

| Concern | Today | Recommendation |
|---|---|---|
| Biological placement | Not modelled | Post-trim P-addition. Hairpin → trim biology approximated as trim → P, with trim distributions absorbing realised biology |
| Ends receiving P | None | V_3, D_5, D_3, J_5 (VDJ); V_3, J_5 (VJ) |
| Relationship to NP regions | n/a | Separate `Region`s with `Segment` = adjacent V/D/J; provenance via `PRegionAdded { end, region }` event + `Nucleotide::flags::P_NUC` |
| Trace addresses | n/a | `p.v_3.length` / `p.d_5.length` / `p.d_3.length` / `p.j_5.length` — lengths only, bases deterministic |
| Events | `RegionAdded` exists, no P variant | New `PRegionAdded { end: PEnd, region }` variant; new `PEnd` enum |
| AIRR fields | n/a | Four new ints: `p_v_3_length` / `p_d_5_length` / `p_d_3_length` / `p_j_5_length` |
| Productive-only | Composes with NP only | Same `JunctionStopState` admit-mask machinery extends to per-end P-length filtering |
| Cartridge ownership | Orphan dict (`p_nucleotide_length_probs`) | Typed `ReferenceEmpiricalModels.p_nucleotide_lengths: Dict[str, EmpiricalDistributionSpec]` keyed by end |
| Legacy `p_nucleotide_length_probs` auto-lift | n/a | **Stays deferred** — separate cartridge-migration slice |
| Replay | n/a | Lengths in trace, bases derived from `(allele, trim, length)` via `complement_base` |
| Plan-signature folding | n/a | `fmt_int_dist` over per-end length distributions, parallel to existing `trim.v_3(length=...)` shape |
| Performance | n/a | ≤ 16 extra pool ops per VDJ simulation; ~240 extra bytes in plan signature |
| Pre-flight bugs found | **None.** `p_nucleotide_length_probs` is genuinely orphan; the MCP diagnostic endpoint is read-only. Implementation is one focused engine slice with no architectural surprises. | — |

The P-nucleotide modeling story is **architecturally tractable**:
the engine's existing `Region`/`Segment` /`SimulationEvent`
discipline, the `complement_base` primitive, and the
`JunctionStopState` admit-mask machinery all compose cleanly
with a per-end `PAdditionPass`. The slice closes the last
material junction-biology gap left after the N-addition story
shipped (typed length models → typed base models → Markov
generator → release consolidation). Pre-trim P, per-base AIRR
fields, and legacy auto-lift remain explicitly out of scope per
§15.

The recommended next step is the single P-nucleotide
implementation slice per §13. Whether the slice is worth
running depends on the user's downstream requirements
(simulation realism vs incremental engine complexity); the
audit's job is to make that call easy to make.
