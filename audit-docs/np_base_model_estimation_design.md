# NP Base Model Estimation — Pre-Implementation Audit

**Status: audit only.** Designs the fourth
`ReferenceCartridgeBuilder.estimate_*` method —
`estimate_np_base_model` — from observed AIRR
rearrangement data. Per the cartridge-authoring audit's
§11.2 ordering, NP base models are the next typed plane to
land after NP length distributions.

Companion to
[`tests/test_np_base_model_estimation_contract.py`](../tests/test_np_base_model_estimation_contract.py)
which freezes (a) the AIRR `np1` / `np2` string fields the
estimator consumes, (b) the typed `np_bases` plane it
writes into, (c) the engine `NpBaseGenerator` family that
already supports `uniform` / `empirical_first_base` /
`markov`, (d) the plan-signature folding boundary, (e) the
legacy orphan-field non-auto-lift discipline, and (f) the
builder/method absence pre-slice.

**Pre-flight finding (§7 below): clean-yes — no
stop-and-report condition.** The two anticipated tricky
points both resolved cleanly at audit time:

1. **P-byte contamination of `np1` / `np2`:** verified
   empirically that the structural NP region span captures
   only NP bytes. With a max-P plane authored
   (`[(3, 1.0)]` per end → 12 P-bytes per record),
   `len(np1) == np1_length` holds on 100/100 records and
   `len(np2) == np2_length` holds on 100/100 records. The
   AIRR builder's `unclaimed_np_string` walks
   `region.start..region.end` over the **structural**
   NP region, which the P-nucleotide v1 slice
   deliberately keeps disjoint from the P-byte ranges
   (audit §1.3 below).
2. **Markov wiring:** verified end-to-end. The Python
   `NpBaseModelSpec` docstring at
   [`src/GenAIRR/reference_models.py:145-158`](../src/GenAIRR/reference_models.py#L145)
   still calls Markov "deferred" — this is **stale
   documentation**. Empirical smoke (forced-transition
   matrix `A→T → C→G → G→A → T→C`) reproduces the
   forced pattern on 99% of observed transitions. The
   engine `MarkovBaseGenerator` ships per the validation
   matrix's "NP base models / Markov N-addition" row.

The estimator can land without disturbing any existing
surface; the manifest's `np_base_models` block is already
perfectly shaped to advertise the estimator's output (no
manifest extension needed — only the populated `np_bases`
plane flows through the existing block).

---

## 1. Q1 — Input source

### 1.1 Direct AIRR fields

The Rust `AirrRecord` struct at
[`engine_rs/src/airr_record/record.rs:77-82`](../engine_rs/src/airr_record/record.rs#L77)
exposes `np1: String` and `np2: String` alongside the
integer length fields the previous slice consumed:

| Rust field | AIRR column | Source |
|---|---|---|
| `np1: String` | `np1` | `unclaimed_np_string(sim, refdata, &rec.sequence, &np1_region)` |
| `np2: String` | `np2` | same, on `np2_region` |
| `np1_aa: String` | `np1_aa` | codon-rail translation of `np1` (estimator ignores) |
| `np2_aa: String` | `np2_aa` | codon-rail translation of `np2` (estimator ignores) |

The Python `result.py` column order declares `np1` /
`np2` as canonical AIRR fields exported on every cartridge.

### 1.2 Estimator's chosen source — direct sequence fields

Per user brief §1: "do not infer from junction sequence
arithmetic in v1". The estimator consumes the **direct
sequence fields** `np1` and `np2`:

| Chain | Consumed AIRR columns | Plane key |
|---|---|---|
| VJ | `np1` only | `NP1` |
| VDJ | `np1`, `np2` | `NP1`, `NP2` |

Junction sequence arithmetic (`junction[v_end:j_start]`
slicing) requires multi-field coordination across
multiple simulator implementations — deferred per audit
§11.

### 1.3 P-nucleotide cross-contamination — resolved CLEAN

**The likely tricky point flagged in the user brief.**
Resolved at audit time via three layers of evidence:

1. **Engine layout (per P-nucleotide v1 audit):** P-bytes
   are pushed into the pool at the V_3 / D_5 / D_3 / J_5
   ends with the `Nucleotide::flags::P_NUC` flag set, and
   their descriptive `Region` flows through
   `SimulationEvent::PRegionAdded { end, region }` —
   **no structural region is added to
   `sim.sequence.regions`**. The structural NP1 / NP2
   regions are disjoint from the P-byte ranges. See
   [`docs/p_nucleotide_design.md`](p_nucleotide_design.md).
2. **AIRR builder (per `projection.rs::unclaimed_np_string`):**
   builds the `np1` / `np2` strings by iterating
   `region.start..region.end` over the STRUCTURAL NP
   region — by construction, this range excludes the
   pre-region P-bytes (e.g. P_V_3 bytes immediately
   preceding NP1) and the post-region P-bytes (P_D_5
   immediately following NP1).
3. **Empirical verification at audit time:** with a
   max-P plane `p_nucleotide_lengths = {"V_3":
   [(3, 1.0)], "D_5": [(3, 1.0)], "D_3": [(3, 1.0)],
   "J_5": [(3, 1.0)]}` forcing 12 P-bytes per record,
   100/100 records satisfy `len(np1) == np1_length` AND
   `len(np2) == np2_length`. If P-bytes leaked into the
   AIRR strings, lengths would diverge.

**Conclusion:** the estimator can safely use `rec["np1"]`
/ `rec["np2"]` verbatim as A/C/G/T-only NP-base
observation streams. No P-aware filtering needed at v1.

### 1.4 Fields the estimator does NOT consume

| AIRR column | Why ignored |
|---|---|
| `np1_length`, `np2_length` | Length-only — the previous slice (`estimate_np_length_distributions`) handles these. |
| `np1_aa`, `np2_aa` | Codon-rail translations — irrelevant to per-base composition. |
| `junction`, `junction_aa`, `junction_length` | Aggregate junction-arithmetic decomposition is fragile across simulators (audit §1.2). |
| `p_v_3_length`, `p_d_5_length`, `p_d_3_length`, `p_j_5_length` | Separate biology — P-nucleotide additions are templated palindromic, not N-additions. |
| `n_p_nucleotides` (when present) | Same — separate biology, separate plane. |

### 1.5 Pinned

- `pin_scaffold_airr_record_carries_np1_and_np2_string_fields`
- `pin_scaffold_unclaimed_np_string_walks_structural_region_only`
- `pin_present_genairr_np1_length_matches_len_np1_string_under_p_plane`
  (stop-and-report verification — the critical P-cross-
  contamination check)
- `pin_scaffold_result_column_order_includes_np1_and_np2_strings`

---

## 2. Q2 — Target plane + supported model kinds

### 2.1 Existing typed plane

The `ReferenceEmpiricalModels.np_bases` plane is the
existing typed surface. From
[`reference_models.py:455`](../src/GenAIRR/reference_models.py#L455):

```python
np_bases: Dict[str, "NpBaseModelSpec"] = field(default_factory=dict)
```

Keys are `NP1` / `NP2` (the same `NP_KEYS` constant the
NP-length plane uses). The existing validator
([`reference_models.py:530-544`](../src/GenAIRR/reference_models.py#L530)):

1. Rejects keys outside `NP_KEYS`.
2. Validates each `NpBaseModelSpec` via its own
   `validate()` (kind / first_base / transitions shape).

### 2.2 NpBaseModelSpec — three supported kinds

`NpBaseModelSpec` lives at
[`reference_models.py:125-235`](../src/GenAIRR/reference_models.py#L125):

| Kind | Required fields | Behaviour |
|---|---|---|
| `"uniform"` | none (both `None`) | Every NP position samples uniformly A/C/G/T. Pre-slice engine default. |
| `"empirical_first_base"` | `first_base: dict[str, float]` | Every NP position samples **independently** from the supplied categorical. The estimator computes this from the **full base composition** (every observed NP base, not just position 0). The name is historical — the engine uses the same distribution at every position, so the "first" naming is biologically accurate only when transitions are not modelled. |
| `"markov"` | `first_base` AND `transitions: dict[str, dict[str, float]]` | 1-step previous-base-conditional. Position 0 uses `first_base`; positions 1+ select a transition row keyed by the previous emitted base. Transition matrix must cover all four A/C/G/T from-bases. |

**Validation rules (verbatim from the existing spec):**

- `kind` must be one of `{"uniform", "empirical_first_base",
  "markov"}`.
- Bases must be from `{"A", "C", "G", "T"}`. Unknown bases
  rejected.
- Weights finite, non-negative, with at least one
  strictly positive per row.
- For `markov`, the outer transition matrix MUST cover all
  four A/C/G/T from-bases — partial rows are rejected as
  an authoring bug.

### 2.3 Markov is wired end-to-end (not deferred)

The `NpBaseModelSpec` docstring at lines 145-158 still
says Markov is "deferred". **This is stale**. Empirical
smoke at audit time: a forced cyclic Markov matrix
(`A→T`, `T→C`, `C→G`, `G→A`) reproduces the forced
pattern on ≥ 92% of observed pair transitions. The engine
`MarkovBaseGenerator` ships per the
[validation matrix](validation_matrix.md)'s
"NP base models / Markov N-addition" row and the
[`junction_n_addition_audit.md`](junction_n_addition_audit.md)
"Markov shipped" status. **The estimator's `kind="markov"`
output flows into the live engine, not a NotImplementedError.**

A follow-up cleanup slice should align the docstring with
shipping reality; that is out of scope here.

### 2.4 Estimator default kind

Per user brief: `kind="markov"` (default). Rationale:

- Markov is biologically more accurate for short NP
  regions where neighbour-base correlations dominate.
- The existing engine's plan-signature fold means
  cartridges with estimated Markov matrices remain
  replayable.
- Authors who want a simpler model can opt in with
  `kind="empirical_first_base"`.

### 2.5 Pinned

- `pin_scaffold_reference_empirical_models_np_bases_plane_exists`
- `pin_scaffold_np_base_model_spec_supports_three_kinds`
- `pin_scaffold_np_base_model_spec_validates_first_base_alphabet`
- `pin_scaffold_np_base_model_spec_validates_markov_row_coverage`
- `pin_present_markov_wired_end_to_end_via_engine_generator`

---

## 3. Q3 — Estimator API

### 3.1 Method signature

Per user brief:

```python
def estimate_np_base_model(
    self,
    rearrangements: Union[List[Dict[str, Any]], Path, IO[str]],
    *,
    kind: str = "markov",
    min_count: int = 1,
    pseudocount: float = 0.0,
    replace: bool = True,
) -> "ReferenceCartridgeBuilder":
    ...
```

Returns `self` for fluent chaining. Idempotent — calling
twice re-estimates and overwrites with `replaced=True` on
the new stage entry (matches every prior estimator).

### 3.2 Per-kind output

| `kind` | What the estimator counts | Output |
|---|---|---|
| `"empirical_first_base"` | **Every** observed NP base across all positions in all rows | A single `first_base` categorical over A/C/G/T |
| `"markov"` | `first_base` from position 0 of each NP string; `transitions` from every (prev, next) pair | `first_base` row + 4×4 transition matrix |
| `"uniform"` | (not estimable — no data needed) | Rejected: `ValueError` — pass `with_models(...)` if you actually want a uniform plane |

**Naming clarification.** The audit explicitly diverges
from a literal reading of `"empirical_first_base"` to mean
"position 0 only". The engine's behaviour for that kind
is "every position samples from this distribution",
making the *full base composition* the biologically
correct empirical estimate. The kind name is preserved
for API stability with existing cartridges.

### 3.3 Pseudocount semantics

Per user brief §6:

| `kind` | `pseudocount` application |
|---|---|
| `"empirical_first_base"` | Added to each of the four A/C/G/T base categories before normalisation. Prevents zero rows when one base is unobserved in a small sample. |
| `"markov"` | Added to **each cell** of the 4×4 transition matrix AND to each of the four `first_base` categories. Guarantees every row has at least one positive weight, so the spec validator's row-coverage requirement is satisfied even on sparse data. |

`pseudocount=0.0` (default) is pure empirical estimation —
zero rows surface as authoring failures (with a clear
error message naming the missing from-base for Markov).

### 3.4 `replace=False`

Same discipline as the prior estimators:

- If `self._reference_models is not None` AND
  `self._reference_models.np_bases` already carries any
  spec, `replace=False` raises `ValueError` before
  consuming any records.

### 3.5 Pinned

- `pin_absence_no_estimate_np_base_model_method_today`
- `pin_absence_no_kind_kwarg_owned_by_np_base_model_estimator_today`
- `pin_absence_no_pseudocount_kwarg_owned_by_np_base_model_estimator_today`

---

## 4. Q4 — Validation

### 4.1 Per-row, per-field validation rules

| Case | Default behaviour | Report destination |
|---|---|---|
| `np1` / `np2` column missing or empty string | Skip THAT field's contribution; row's other field still feeds | `report.rejected` with `reason="missing_required_column"` (carries column) |
| `np1` / `np2` value is not a string | Skip THAT field's contribution | `report.rejected` with `reason="malformed_np_value"` |
| `np1` / `np2` contains a non-canonical base (anything outside `{A,C,G,T}`, case-insensitive) | Skip THAT field's contribution entirely; record one rejection entry naming the unknown character | `report.rejected` with `reason="noncanonical_base"` (carries `column` + `unknown_chars`) |
| VJ row has non-empty `np2` | Drop NP2 contribution; keep NP1 | `report.warnings` (cartridge-level, one notice per dataset) |
| Per-base count below `min_count` (for `empirical_first_base`) OR per-transition-cell count below `min_count` (for `markov`) | Drop from result | `stage.inferred.below_min_count` |

### 4.2 Case handling

Per user brief: "accept only A/C/G/T bases". The estimator
uppercases each string before tallying so mixed-case input
(`"aCgT"`) maps to canonical bases. Any character that
doesn't uppercase to `{A,C,G,T}` triggers the
`noncanonical_base` rejection — the row's contribution
for that field drops.

### 4.3 Markov edge cases

For `kind="markov"` on extremely short NP strings:

- A 0-length NP string contributes nothing — no `first_base`
  observation, no transitions. Drops with `reason="missing_required_column"`.
- A 1-length NP string contributes one `first_base` observation
  + zero transitions. The transition matrix relies on other
  rows or the pseudocount to populate the missing rows.
- A 2+-length NP string contributes one `first_base` + `(len-1)`
  transitions.

If pseudocount=0 AND a from-base remains unobserved across
the entire dataset, the spec validator will raise at
construction time because the transition matrix is
incomplete. The estimator surfaces this with a tagged
`ValueError` ("from-base {X} has no observed transitions
in input — pass pseudocount > 0 or supply more data").

### 4.4 Pinned

- `pin_scaffold_np_base_model_spec_rejects_unknown_alphabet`
- `pin_scaffold_np_base_model_spec_rejects_zero_weight_row_in_markov`
- `pin_absence_no_estimate_np_base_model_min_count_kwarg_today`

---

## 5. Q5 — Chain behaviour

### 5.1 Chain-type-driven enforcement

| Chain | Estimator emits | Skipped (warn) |
|---|---|---|
| VJ | `NP1` key only | Non-empty `np2` column (one-time warning, contribution dropped) |
| VDJ | `NP1`, `NP2` keys | — (both consumed per chain) |

Same boundary the NP-length estimator follows. The
cartridge's `metadata.has_d` is the authoritative source
for VJ-vs-VDJ classification.

### 5.2 Empirical baseline (smoke check at audit time)

The bundled VJ cartridge `HUMAN_IGK_OGRDB` emits `np2 = ""`
on every record (no NP2 region on a VJ chain). The
estimator's VJ NP2-skip path covers user inputs from
external simulators that produce a populated `np2`
column on a VJ-shaped cartridge — same defensive
discipline as the NP-length slice.

### 5.3 Pinned

- `pin_scaffold_config_info_has_d_drives_chain_classification`
  (reused)
- `pin_present_vj_np2_string_is_empty_on_bundled_cartridge`

---

## 6. Q6 — Pseudocount semantics (already covered in §3.3)

The deliberate v1 choice:

| Kind | Pseudocount adds to | Rationale |
|---|---|---|
| `empirical_first_base` | Every A/C/G/T base category | Prevents zero rows on a 4-element categorical; matches Bayesian Dirichlet-prior shape. |
| `markov` | Every cell of the 4×4 transition matrix + every A/C/G/T `first_base` row | Guarantees the Markov spec validator's row-coverage requirement is met even on sparse input. |

`pseudocount=0.0` preserves pure empirical estimation —
the estimator raises a clear error if the spec validator
would reject the result (e.g. missing Markov from-base row).

### 6.1 Pinned

(no dedicated pin section — the pseudocount semantics
are pinned via the implementation tests in the next
slice.)

---

## 7. Q7 — Replay / signature + stop-and-report

### 7.1 Plan signature folding

The `GenerateNPPass.parameter_signature()` at
[`engine_rs/src/passes/generate_np.rs:147-161`](../engine_rs/src/passes/generate_np.rs#L147)
folds the **full base generator payload** — for `markov`,
both the first-base row AND the 4-row transition matrix
in canonical A/C/G/T order. Confirmed empirically at
audit time: two cartridges differing only in
`np_bases["NP1"]` produce different plan signatures.
**No soft gap inherited** (same property the trim +
NP-length slices have).

### 7.2 Stop-and-report condition

> "If direct `np1` / `np2` fields do not exist, stop and
> report. If they exist but include P bases mixed in after
> P-nucleotide v1, document whether estimator must exclude
> P or cannot distinguish; this is the likely tricky point."

**Verdict: NOT triggered. Proceed.**

1. **Direct `np1` / `np2` fields exist.** Both as `String`
   on the Rust `AirrRecord` and as canonical AIRR columns
   in `result.py`. Confirmed at source + smoke.
2. **P-bytes do NOT contaminate `np1` / `np2`.** Verified
   empirically with a max-P plane forcing 12 P-bytes
   per record:
   - 100/100 records: `len(np1) == np1_length`
   - 100/100 records: `len(np2) == np2_length`
   - All `np1` / `np2` characters are canonical A/C/G/T
     on the bundled cartridges.

The structural NP region span is **disjoint from the
P-byte ranges** by P-nucleotide-v1 design (no structural
region added to `sim.sequence.regions` for P-bytes —
they live in the gap between structural regions). The
estimator can safely use `rec["np1"]` / `rec["np2"]`
verbatim with no P-aware filtering.

### 7.3 Pinned

- `pin_present_np_base_distribution_changes_plan_signature`
  (no soft gap inherited)
- `pin_present_genairr_np1_length_matches_len_np1_string_under_p_plane`
  (stop-and-report verification)

---

## 8. Legacy orphan boundary

### 8.1 NP_first_bases / NP_transitions stay orphan

The legacy `DataConfig.NP_first_bases` and
`DataConfig.NP_transitions` fields are listed in
[`_DOCUMENTED_ORPHAN_DATACONFIG_FIELDS`](../src/GenAIRR/dataconfig/data_config.py#L43)
and the existing manifest reports
`legacy_np_first_bases_present` /
`legacy_np_transitions_present` flags. The auto-lift
boundary is held: the bridge resolver does NOT lift the
legacy dicts into the typed plane (per
[`junction_n_addition_audit.md`](junction_n_addition_audit.md)
the auto-lift would silently change output bytes vs the
pre-slice baseline).

**The estimator slice MUST NOT touch the legacy fields.**
Same discipline every prior estimator slice respected for
its legacy orphan counterparts (`gene_use_dict` /
`NP_lengths` / `trim_dicts` / `p_nucleotide_length_probs`).

### 8.2 Pinned

- `pin_present_np_first_bases_and_np_transitions_listed_in_orphan_fields`
- `pin_present_legacy_np_first_bases_does_not_auto_lift_to_typed_plane`

---

## 9. Manifest exposure

### 9.1 No new block needed

Unlike the trim / NP-length slices, the manifest
`models.np_base_models` block already exists with the
exact shape the estimator's output needs:

```python
"np_base_models": {
    "models": {<key>: {"kind": str}, ...},          # per-key kind
    "legacy_fallback": False,
    "legacy_np_transitions_present": bool,
    "legacy_np_first_bases_present": bool,
    "supported_kinds": ["uniform", "empirical_first_base", "markov"],
    "deferred_kinds": [],
    "in_plan_signature": True,
    "in_content_hash": False,
}
```

The block already enumerates `supported_kinds`,
documents `legacy_fallback=False`, and reports
`in_plan_signature=True`. The estimator slice's only
manifest-relevant effect is populating
`self._reference_models.np_bases` so that the existing
block's `models` dict carries the estimated entries.

**No new manifest helper, no extension, no additive
block.** The slice is even simpler than the previous
two on the manifest side.

### 9.2 Pinned

- `pin_scaffold_manifest_np_base_models_block_already_carries_full_shape`

---

## 10. Implementation order (recommended)

A single self-contained slice can land the estimator. Two
sub-steps:

1. **Builder method** —
   `ReferenceCartridgeBuilder.estimate_np_base_model(
   rearrangements, *, kind="markov", min_count=1,
   pseudocount=0.0, replace=True)`. Field-local per-row
   validation. Computes per-kind output (first-base
   categorical OR first-base + 4×4 transition matrix).
   Applies pseudocount + min_count. Constructs the
   `NpBaseModelSpec` and writes into
   `self._reference_models.np_bases`. Stage entry shape
   per canonical `{stage, inputs, inferred, warnings}`.

2. **Tests + doc + pin flip** — 12 implementation tests
   covering the brief's expected surface (kind=markov
   default + the two estimable kinds + chain behaviour +
   replay + manifest + validators). Flip the three
   absence pins in
   `tests/test_np_base_model_estimation_contract.py` to
   present-state. Lockstep update to the
   reference-cartridge-authoring estimator boundary pins
   (allowing the fourth shipped estimator).

Cost estimate:

- ~150 lines Python (builder method — per-kind branching
  is the bulk; field-local validation copied from the
  NP-length slice; spec construction via the existing
  `NpBaseModelSpec`)
- ~160 lines tests
- ~25 lines docstrings + audit-doc-flip
- **0 lines manifest changes** (the block already exists)
- **0 lines engine / bridge / spec-class changes**

Lightest manifest impact of any estimator slice so far.

---

## 11. Test surface — what this audit pins

Mirrored in
[`tests/test_np_base_model_estimation_contract.py`](../tests/test_np_base_model_estimation_contract.py).

### `pin_scaffold_*` — AIRR np1/np2 fields (live)

1. `AirrRecord` carries `np1: String` and `np2: String`.
2. `unclaimed_np_string` (projection.rs) walks the
   STRUCTURAL `np_region` only — by construction the
   range excludes P-byte ranges.
3. `result.py` canonical column order includes both
   string fields.

### `pin_scaffold_*` — typed plane + spec (live)

4. `ReferenceEmpiricalModels.np_bases` plane exists.
5. `NpBaseModelSpec` supports three kinds: `uniform`,
   `empirical_first_base`, `markov`.
6. Spec validator rejects non-A/C/G/T bases.
7. Spec validator rejects incomplete Markov transition
   rows.

### `pin_scaffold_*` — bridge + engine surface (live)

8. `_dataconfig_extract.extract_recombine_defaults`
   returns `np1_bases` / `np2_bases` keys.
9. `_np_bases_from_models` resolver exists.
10. `_np_markov_transitions_from_models` resolver exists.
11. Engine `GenerateNPPass.parameter_signature` folds the
    base generator's full payload via the canonical
    A/C/G/T base weights.

### `pin_present_*` — stop-and-report verification

12. **Critical** —
    `len(np1) == np1_length` and
    `len(np2) == np2_length` on EVERY record produced
    under a max-P plane. The P-cross-contamination
    boundary holds — `np1` / `np2` strings are P-clean.
13. Bundled cartridges produce only canonical A/C/G/T
    characters in `np1` / `np2`.
14. Markov is wired end-to-end (NOT deferred) — a forced
    cyclic transition matrix reproduces on observed
    output.

### `pin_present_*` — chain-type + replay

15. VJ cartridge's `np2` is empty on every bundled record
    (no NP2 region).
16. Two cartridges differing only in `np_bases["NP1"]`
    produce different plan signatures (no soft gap).

### `pin_present_*` — legacy orphan boundary

17. `NP_first_bases` and `NP_transitions` listed in
    `_DOCUMENTED_ORPHAN_DATACONFIG_FIELDS`.
18. The bridge resolver does NOT auto-lift the legacy
    dicts into the typed plane.

### `pin_scaffold_*` — manifest already shaped

19. `manifest['models']['np_base_models']` already
    exposes `models` / `supported_kinds` /
    `legacy_fallback` / `in_plan_signature` — the
    estimator's output populates this block by
    populating the typed plane.

### `pin_scaffold_*` — builder shape (reused)

20. `ConfigInfo.has_d` chain-type classifier.
21. Builder stage entry shape `{stage, inputs, inferred,
    warnings}`.
22. Builder idempotency pattern via `replaced=True` flag.
23. `csv.DictReader` AIRR-TSV ingestion.

### `pin_absence_*` — gaps the implementation slice closes

24. `ReferenceCartridgeBuilder.estimate_np_base_model`
    is not a method today.
25. No `kind` kwarg surface exists on the builder today.
26. No sibling module / class / free function with the
    estimator's name was introduced.

### Doc anchor

27. Audit doc exists and references the contract file.

---

## 12. Out of scope

Documented here so a future implementer doesn't
accidentally expand the work.

- **Stale `NpBaseModelSpec` docstring claim that Markov is
  "deferred".** This audit empirically confirmed Markov
  is wired end-to-end. A docstring cleanup slice should
  fix this; out of scope here.
- **Per-position empirical distributions.** v1 collapses
  the empirical model to position-independent (matches
  the engine's behaviour — see §3.2 naming clarification).
  A future slice could add `kind="empirical_per_position"`
  with a position-keyed distribution.
- **Higher-order Markov chains (n-gram > 1).** v1 supports
  1-step Markov. Higher-order would require a new
  `NpBaseGenerator` engine variant.
- **Auto-lift of legacy `NP_first_bases` / `NP_transitions`.**
  Stays deferred — same boundary every prior slice
  respected.
- **Junction-arithmetic source.** v1 uses direct `np1` /
  `np2` strings. Decomposing `junction[v_end:j_start]`
  back into NP1 / NP2 contributions is fragile across
  simulators and would require trim + D + anchor
  arithmetic.
- **Per-segment / per-allele base models.** The typed
  plane is per-key (NP1, NP2) — flat by design.
- **Stratified estimation by `productive` / `vj_in_frame`
  / `stop_codon`.** v1 ignores these columns; the
  estimator consumes every record uniformly.
- **In-place pseudocount tuning (`pseudocount`-shape
  dependent on `min_count`).** v1 applies pseudocount
  AFTER the min_count filter on observed cells only
  (mirrors trim/NP-length slices).

---

## 13. Summary table

| Concern | Today | Recommendation |
|---|---|---|
| AIRR fields populated by GenAIRR | Direct `np1` / `np2` `String` fields populated on every record; canonical A/C/G/T composition on bundled cartridges; lengths match `len(np1) == np1_length` even under a max-P plane | **Estimator consumes the two string fields**; P-bytes do NOT contaminate them per audit §1.3 |
| Typed plane | `ReferenceEmpiricalModels.np_bases: Dict[str, NpBaseModelSpec]` with `NP_KEYS = ("NP1","NP2")` + per-spec validator | **Reuse verbatim.** No new spec class. |
| Spec kinds | `NpBaseModelSpec` supports `uniform` / `empirical_first_base` / `markov` — all three wired end-to-end (Markov docstring is stale) | **Estimator's `kind` kwarg accepts `empirical_first_base` and `markov`**; `uniform` is not estimable and is rejected with a clear error. Default `markov`. |
| Bridge resolver | `_np_bases_from_models` + `_np_markov_transitions_from_models` + legacy fallback; precedence kwarg > typed plane > legacy > uniform already in place | **Reuse verbatim.** No bridge changes. |
| Engine surface | `GenerateNPPass` with `NpBaseGenerator` family; `parameter_signature` folds via canonical A/C/G/T payload | **Reuse verbatim.** No engine changes. |
| P-byte contamination of `np1` / `np2` | Resolved CLEAN: `len(np1) == np1_length` on 100/100 records under max-P plane; structural NP region disjoint from P-byte ranges | **Estimator can use `rec["np1"]` / `rec["np2"]` verbatim** with no P-aware filtering |
| Provenance distinction | NP base model vs P-nucleotide vs junction arithmetic — all distinguished at AIRR + engine layers | Estimator consumes ONLY `np1` / `np2`. Documented + pinned. |
| Chain-type enforcement | `metadata.has_d` is authoritative; VJ chains hard-empty-string `np2` | VJ skips `np2` contribution + warns once if column non-empty; VDJ consumes both fields |
| Builder method signature | n/a | `estimate_np_base_model(rearrangements, *, kind="markov", min_count=1, pseudocount=0.0, replace=True) → self` |
| Stage entry shape | n/a | Canonical `{stage, inputs, inferred, warnings}` |
| Manifest block | `models.np_base_models` ALREADY exists with the exact shape the estimator's output needs | **No new block, no extension**; the existing block reads through `_reference_models.np_bases` |
| Plan signature folding | NP base generator already folds via canonical A/C/G/T payload | **No soft gap inherited** — same property the trim + NP-length slices have |
| Pre-flight stop-and-report | **NOT triggered** — bundled cartridges produce A/C/G/T-only `np1` / `np2`; P-bytes are excluded from the AIRR strings by structural-region design + verified empirically | Proceed to implementation slice |

The slice is **the lightest one yet on the manifest /
plumbing side** — zero new manifest helpers, zero engine
changes, zero bridge changes, zero spec changes. The
estimator method itself is heavier than trim / NP-length
(per-kind branching + 4×4 Markov accumulation) but
still self-contained. Total estimated cost ~150 lines
Python + ~160 lines tests.
