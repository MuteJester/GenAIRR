# P-Nucleotide Length Estimation — Pre-Implementation Audit

**Status: audit only.** Designs the fifth
`ReferenceCartridgeBuilder.estimate_*` method —
`estimate_p_nucleotide_lengths` — from observed AIRR
rearrangement data. Per the cartridge-authoring audit's
§11.2 ordering, this is the fifth estimator (sixth and
final remains `estimate_shm_rates`).

Companion to
[`tests/test_p_nucleotide_length_estimation_contract.py`](../tests/test_p_nucleotide_length_estimation_contract.py)
which freezes (a) the AIRR `p_*_length` fields the
estimator consumes, (b) the typed
`p_nucleotide_lengths` plane it writes into, (c) the
chain-type validator + engine PAdditionPass plumbing, (d)
the plan-signature folding, (e) the legacy orphan
boundary, and (f) the builder/method absence pre-slice.

**Pre-flight finding (§7 below): clean-yes — but with a
documented utility caveat.** The plumbing is clean
(direct AIRR fields exist, default to zero without a P
plane, populate reliably with a P plane, plan-signature
folds the distribution, validator chain-type-rejects D
keys on VJ). The caveat is **scope of applicability**:
external AIRR tools do not model P-nucleotide additions,
so external AIRR data does not carry `p_*_length`
columns. The estimator's realistic inputs are:

- GenAIRR's own AIRR records produced from a cartridge
  that **already** carries a P-plane (re-estimation /
  cartridge migration use case).
- Custom datasets where the user explicitly populated
  `p_*_length` columns from their own analysis.

The audit recommends shipping the estimator with a
**strong provenance warning** in the docstring and a
**stage-level warning** when all rows lack the P field
for a given key (indicating the user is running against
P-naïve data).

Empirical smoke at audit time:

| Cartridge | Baseline (no P plane) | With P plane (`V_3=[(0,0.5),(1,0.3),(2,0.2)]`, etc.) |
|---|---|---|
| `HUMAN_IGH_OGRDB` baseline | 0/50 records have any non-zero `p_*_length` | n/a |
| `HUMAN_IGH_OGRDB` + P plane | n/a | `p_v_3_length` 28/50 nonzero, `p_d_5_length` 11/50, `p_d_3_length` 13/50, `p_j_5_length` 11/50 |
| `HUMAN_IGK_OGRDB` baseline | 0/50 records have any non-zero `p_*_length` | n/a |

---

## 1. Q1 — Input source

### 1.1 Direct AIRR fields available today

The Rust `AirrRecord` struct at
[`engine_rs/src/airr_record/record.rs:103-106`](../engine_rs/src/airr_record/record.rs#L103)
exposes four direct integer `p_*_length` fields:

| Rust field | AIRR column | Source |
|---|---|---|
| `p_v_3_length: i64` | `p_v_3_length` | event ledger sum of `SimulationEvent::PRegionAdded { end: V_3 }` region lengths |
| `p_d_5_length: i64` | `p_d_5_length` | event ledger sum on `D_5` |
| `p_d_3_length: i64` | `p_d_3_length` | event ledger sum on `D_3` |
| `p_j_5_length: i64` | `p_j_5_length` | event ledger sum on `J_5` |

`result.py`'s canonical AIRR column order
([`src/GenAIRR/result.py`](../src/GenAIRR/result.py))
declares all four as exported columns on every cartridge.

### 1.2 Default-zero behaviour when no P plane runs

Confirmed empirically at audit time on the bundled
cartridges (no typed P-plane authored, no `Experiment.p_addition(...)`
DSL — neither exists at v1):

- `HUMAN_IGH_OGRDB`: 0/50 records carry any non-zero
  `p_*_length` field across all four ends.
- `HUMAN_IGK_OGRDB`: same — 0/50 records carry any
  non-zero `p_*_length` field.

The defaults flow from the AIRR builder's event-ledger
walk: when no `PAdditionPass` ran, no `PRegionAdded`
events are emitted, so the per-end length sum is zero by
construction. **Records produced without a P-plane are
indistinguishable from records produced with a P-plane
where every length sample happened to be zero.** The
estimator slice MUST emit a stage-level warning when
this ambiguity could matter — see §5.2.

### 1.3 Estimator's chosen source — direct integer fields

Per user brief §2: "v1 should consume only explicit P
fields". The estimator consumes the four direct integer
fields:

| Chain | Consumed AIRR columns | Plane keys |
|---|---|---|
| VJ | `p_v_3_length`, `p_j_5_length` | `V_3`, `J_5` |
| VDJ | `p_v_3_length`, `p_d_5_length`, `p_d_3_length`, `p_j_5_length` | `V_3`, `D_5`, `D_3`, `J_5` |

The mapping is verbatim — AIRR column name → plane key —
modulo the trailing `_length` strip.

### 1.4 Fields the estimator does NOT consume

| AIRR column | Why ignored |
|---|---|
| `junction`, `junction_aa`, `junction_length` | v1 does NOT infer P lengths from junction arithmetic — same boundary every prior length estimator respected (NP-length estimator §5.2 audit-doc reference). |
| `np1`, `np2`, `np1_length`, `np2_length` | NP region is biologically distinct from P region; even if `unclaimed_np_string` is P-clean (verified in NP-base-model audit), deriving P lengths from NP geometry would be a fundamentally different inference problem. |
| `v_trim_3`, `d_trim_5`, `d_trim_3`, `j_trim_5` | Recombination-stage trims, separate biology. The PAdditionPass runs AFTER trim at each end; P-length depends on `(allele, trim, orientation, length)` but the estimator infers length from the AIRR-emitted integer, not from re-deriving via trim arithmetic. |

### 1.5 Pinned

- `pin_scaffold_airr_record_carries_four_p_length_fields`
- `pin_scaffold_python_airr_projection_emits_p_length_fields`
- `pin_scaffold_result_column_order_includes_p_length_fields`
- `pin_present_baseline_no_p_plane_produces_zero_p_lengths`
  (stop-and-report verification: ambiguity boundary held)
- `pin_present_authored_p_plane_produces_nonzero_p_lengths`
  (stop-and-report verification: estimator has signal to
  consume when input has provenance)

---

## 2. Q2 — Estimation validity

### 2.1 External AIRR data limitation

**The fundamental utility constraint.** External AIRR
analysis tools (IgBLAST, MiXCR, repertoire analysis
pipelines) do NOT model P-nucleotide additions —
P-nucleotide biology is a recombination-stage modelling
decision, not an observation-stage analysis decision. So
AIRR records produced by external tools:

- Either omit the `p_*_length` columns entirely (most
  common case).
- Or populate them as `0` per the AIRR-C schema's
  optional-field default.

Running the estimator against such records would produce
either:

- A "no records carried the column" rejection storm (all
  rows in `report.rejected` with
  `reason="missing_required_column"`).
- A degenerate `[(0, 1.0)]` distribution (every observed
  length is zero — biologically meaningless).

### 2.2 v1 boundary — consume only explicit P fields

Per user brief: "v1 should consume only explicit P
fields". The estimator:

- Reads `p_v_3_length` / `p_d_5_length` / `p_d_3_length`
  / `p_j_5_length` directly.
- Does NOT infer from junction arithmetic, NP strings,
  trim lengths, or any other column.
- Surfaces a **stage-level warning** when more than 95%
  of consumed rows reported `0` for a given key (likely
  P-naïve data). Threshold rationale: even a low-rate
  P-addition cartridge yields at least one nonzero P
  length per ~20 records empirically (see audit smoke
  table).

### 2.3 No heuristic inference

Per user brief §5: "do not infer P lengths from junction
arithmetic or NP strings". This is the v1 scope boundary.
A future slice could attempt P-vs-N disambiguation from
the assembled junction (which is biologically
ambiguous — palindromic complement of `N...N` is hard to
distinguish from random N-addition), but that involves
modelling decisions out of scope for the builder-method
shape.

### 2.4 Pinned

- `pin_present_external_airr_tools_do_not_populate_p_length_fields`
  (documented in design doc; pinned by the test
  asserting baseline-cartridge records have all-zero P
  lengths — same observation surface)
- `pin_absence_no_junction_arithmetic_p_inference_today`
- `pin_absence_no_np_string_p_inference_today`

---

## 3. Q3 — Target plane

### 3.1 Existing typed plane

The `ReferenceEmpiricalModels.p_nucleotide_lengths`
plane is the existing typed surface. From
[`reference_models.py:48-49`](../src/GenAIRR/reference_models.py#L48):

```python
P_NUCLEOTIDE_END_KEYS: Tuple[str, ...] = ("V_3", "D_5", "D_3", "J_5")
P_NUCLEOTIDE_END_KEYS_VJ: Tuple[str, ...] = ("V_3", "J_5")
```

The plane is `Dict[str, EmpiricalDistributionSpec]` —
same shape as `trims` and `np_lengths`. The existing
validator
([`reference_models.py:552-575`](../src/GenAIRR/reference_models.py#L552)):

1. Rejects keys outside `P_NUCLEOTIDE_END_KEYS`.
2. **Rejects D-end keys (`D_5` / `D_3`) on VJ chains at
   attach time** — symmetric with `trims`, asymmetric
   with `np_lengths`. This is the strong-boundary
   discipline the estimator inherits without enforcement
   at estimation time.
3. Delegates per-spec validation to
   `EmpiricalDistributionSpec.validate()`.

### 3.2 Chain-aware key set

| Chain | Plane keys allowed | Estimator writes |
|---|---|---|
| VJ | `V_3`, `J_5` | `V_3`, `J_5` |
| VDJ | `V_3`, `D_5`, `D_3`, `J_5` | all four |

### 3.3 Pinned

- `pin_scaffold_reference_empirical_models_p_nucleotide_lengths_plane_exists`
- `pin_scaffold_p_nucleotide_end_keys_constants_hold`
- `pin_scaffold_p_nucleotide_lengths_validator_rejects_unknown_keys`
- `pin_present_p_nucleotide_lengths_validator_rejects_d_keys_on_vj`

---

## 4. Q4 — Validation

### 4.1 Per-row, per-field validation rules

Same shape as trim / NP-length estimators (field-local,
structured rejection):

| Case | Default behaviour | Report destination |
|---|---|---|
| Required `p_*_length` column missing (V_3 / J_5 always; D_5 / D_3 on VDJ) | Skip THAT field's contribution; row's other fields still feed | `report.rejected` with `reason="missing_required_column"` (carries column) |
| Non-integer value | Skip THAT field's contribution | `report.rejected` with `reason="malformed_length_value"` |
| Negative integer | Skip THAT field's contribution | `report.rejected` with `reason="negative_length_value"` |
| VJ row has non-zero `p_d_5_length` / `p_d_3_length` | Drop those columns; keep V_3 / J_5 | `report.warnings` (one notice per dataset per column) |
| Below `min_count` for a length value | Drop value from result | `stage.inferred.below_min_count` (per-key dict) |

### 4.2 Estimator API kwargs

Mirror the trim / NP-length estimator surface:

- `min_count: int = 1` — drop length values whose
  observed count is strictly below the threshold before
  normalisation.
- `pseudocount: float = 0.0` — additive smoothing
  applied to **observed** length values only.
- `replace: bool = True` — idempotency. `False` blocks
  re-entry when typed-plane `p_nucleotide_lengths` is
  already attached.

### 4.3 Pinned

- `pin_absence_no_estimate_p_nucleotide_lengths_method_today`
- `pin_absence_no_estimate_p_nucleotide_lengths_min_count_pseudocount_today`

---

## 5. Q5 — Provenance warning

### 5.1 v1 strong-warning discipline

The estimator's docstring SHOULD lead with a prominent
provenance warning paragraph:

> **Provenance warning.** This estimator infers
> P-nucleotide length distributions ONLY from records
> that carry the canonical AIRR `p_v_3_length` /
> `p_d_5_length` / `p_d_3_length` / `p_j_5_length`
> columns. External AIRR tools (IgBLAST, MiXCR, …) do
> NOT model P-nucleotide additions — their output will
> either omit these columns entirely or populate them as
> zero. Running the estimator against such data produces
> either a rejection storm (all rows in
> `report.rejected` with reason `missing_required_column`)
> or a degenerate `[(0, 1.0)]` distribution.
>
> Realistic inputs:
>
> - GenAIRR's own AIRR records from a cartridge that
>   already carries a typed P-plane
>   (re-estimation / cartridge-migration use case).
> - Custom datasets where you populated `p_*_length`
>   columns from your own analysis.

### 5.2 Per-key auto-warning

The estimator additionally emits a stage-level warning
under `stage.warnings` when **more than 95% of rows that
contributed to a given key reported `0`**:

```python
"warnings": [
    "V_3: 49/50 contributing rows reported p_v_3_length=0. "
    "Input may be P-naïve (external AIRR tool, or "
    "GenAIRR cartridge without a P-plane). The estimated "
    "distribution will be degenerate.",
    ...
]
```

Threshold rationale: even a low-rate P-addition cartridge
produces at least one nonzero P length per ~20 records
in the audit smoke (`p_v_3_length` 28/50 nonzero on the
test cartridge with `V_3=[(0,0.5),(1,0.3),(2,0.2)]`).
Crossing 95% zero is empirically diagnostic of P-naïve
data.

### 5.3 Pinned

- `pin_absence_no_estimate_p_nucleotide_lengths_warns_about_provenance_today`
  (the dedicated warning surface is introduced by the
  implementation slice, pinned via the
  per-key auto-warning shape in the implementation
  tests)

---

## 6. Q6 — Engine consumption

### 6.1 Bridge resolver — already wired

The bridge's
[`_dataconfig_extract.extract_recombine_defaults`](../src/GenAIRR/_dataconfig_extract.py#L241-L252)
returns four `p_*_lengths` keys with the typed-plane
resolver:

```python
"p_v_3_lengths": _p_nucleotide_lengths_from_models(explicit, "V_3"),
"p_d_5_lengths": _p_nucleotide_lengths_from_models(explicit, "D_5"),
"p_d_3_lengths": _p_nucleotide_lengths_from_models(explicit, "D_3"),
"p_j_5_lengths": _p_nucleotide_lengths_from_models(explicit, "J_5"),
```

`_p_nucleotide_lengths_from_models` lives at
[`_dataconfig_extract.py:417`](../src/GenAIRR/_dataconfig_extract.py#L417).
**No bridge changes** — the typed plane already lowers
end-to-end.

### 6.2 Plan signature folding — confirmed at audit time

Empirical check: two cartridges differing only in
`p_nucleotide_lengths["V_3"]` produce **different** plan
signatures. The `PAdditionPass.parameter_signature()`
implementation at
[`engine_rs/src/passes/p_addition.rs:320`](../engine_rs/src/passes/p_addition.rs#L320)
folds the length distribution via `fmt_int_dist`. **No
soft gap inherited.**

### 6.3 Precedence chain

| Priority | Source |
|---|---|
| 1 (highest) | Explicit `Experiment.recombine(p_*_lengths=...)` — these kwargs are NOT exposed at v1 (P-length distribution is a cartridge-only surface per `docs/p_nucleotide_design.md`); included for future-proofing |
| 2 | Typed `ReferenceEmpiricalModels.p_nucleotide_lengths[key]` plane — **the estimator's output** |
| 3 | None / no-op (P-plane absent → no `PAdditionPass` ran, byte-identical to pre-slice baseline) |

The legacy `DataConfig.p_nucleotide_length_probs` orphan
field is NOT lifted into the typed plane. Same boundary
every prior estimator slice respected.

### 6.4 Pinned

- `pin_scaffold_extract_recombine_defaults_returns_four_p_length_keys`
- `pin_scaffold_p_nucleotide_lengths_from_models_resolver_exists`
- `pin_scaffold_engine_p_addition_pass_folds_length_dist_into_signature`
- `pin_present_p_length_distribution_changes_plan_signature`
- `pin_present_legacy_p_nucleotide_length_probs_orphan_boundary_holds`

---

## 7. Q7 — Stop-and-report verification

### 7.1 The brief's stop-and-report condition

> "If direct P fields are absent or not reliable, stop
> and report. Otherwise implementation should mirror
> `estimate_np_length_distributions` almost exactly,
> with stronger doc warnings about external AIRR data."

### 7.2 Verdict — NOT triggered. Proceed.

Direct AIRR `p_*_length` integer fields exist on both
the Rust `AirrRecord` struct and the Python projection.
Their behaviour is fully characterised:

1. **Reliable when source has P-plane.** Empirical smoke
   on `HUMAN_IGH_OGRDB` + max-P plane: 28/50 records
   carry non-zero `p_v_3_length` (with a 50/30/20 mix of
   0/1/2). The four ends report nonzero at rates
   matching the authored P-length distribution.
2. **Hard-zero when source has no P-plane.** 0/50
   records carry any non-zero `p_*_length` on either
   bundled cartridge at baseline. The estimator running
   against such data would produce a degenerate `[(0,
   1.0)]` distribution per key — biologically
   meaningless but algorithmically clean.

### 7.3 Documented utility caveat — strong-warning shape

**Not a stop-and-report.** The narrow utility (external
AIRR tools don't model P) is handled by the per-key
auto-warning at §5.2. The estimator ships; it just
warns loudly when the input lacks signal.

### 7.4 Pinned

- `pin_present_baseline_no_p_plane_produces_zero_p_lengths`
- `pin_present_authored_p_plane_produces_nonzero_p_lengths`

---

## 8. Legacy orphan boundary

### 8.1 `p_nucleotide_length_probs` stays orphan

The legacy `DataConfig.p_nucleotide_length_probs` field
is listed in
[`_DOCUMENTED_ORPHAN_DATACONFIG_FIELDS`](../src/GenAIRR/dataconfig/data_config.py#L43)
and the existing manifest's `p_nucleotide_models` block
reports `legacy_p_nucleotide_length_probs_present` for
inspection. The auto-lift boundary holds — same
discipline every prior estimator slice respected.

### 8.2 Pinned

- `pin_present_legacy_p_nucleotide_length_probs_orphan_boundary_holds`
  (already in §6.4)

---

## 9. Manifest exposure

### 9.1 Existing `p_nucleotide_models` block

The manifest at
[`dataconfig/data_config.py:870-879`](../src/GenAIRR/dataconfig/data_config.py#L870)
already exposes a structured `p_nucleotide_models`
block:

```python
"p_nucleotide_models": {
    "length_keys": [...],                                        # populated subset of P_NUCLEOTIDE_END_KEYS
    "legacy_p_nucleotide_length_probs_present": bool,
    "legacy_fallback": False,
    "supported_ends": ["V_3", "D_5", "D_3", "J_5"],
    "in_plan_signature": True,                                   # confirmed at audit time
    "in_content_hash": False,
}
```

This matches the shape the NP-base-model slice
inherited — no new manifest helper needed. The
estimator's output (populated typed plane) flows
through `length_keys` automatically.

### 9.2 Pinned

- `pin_scaffold_manifest_p_nucleotide_models_block_already_carries_full_shape`

---

## 10. Implementation order (recommended)

A single self-contained slice. Mirrors the NP-length
estimator discipline almost exactly:

1. **Builder method** —
   `ReferenceCartridgeBuilder.estimate_p_nucleotide_lengths(
   rearrangements, *, min_count=1, pseudocount=0.0,
   replace=True)`. Field-local per-row validation. Reads
   `p_v_3_length` always; `p_j_5_length` always;
   `p_d_5_length` / `p_d_3_length` on VDJ; VJ rows with
   nonzero D-end P columns surface as warnings.
   Per-key auto-warning when ≥95% of contributing rows
   reported zero (provenance heuristic). Constructs
   per-key `EmpiricalDistributionSpec` instances and
   writes into `self._reference_models.p_nucleotide_lengths`.

2. **Tests + doc + pin flip** — 12 implementation tests
   covering the brief's expected surface (per-end smoke
   + chain-type enforcement + field-local validation +
   min_count + pseudocount + replay + manifest + the
   per-key auto-warning for P-naïve input). Flip the
   absence pins in
   `tests/test_p_nucleotide_length_estimation_contract.py`
   to present-state. Lockstep update to the
   reference-cartridge-authoring estimator boundary pins.

3. **Doc + provenance warning** — the docstring carries
   the §5.1 "Provenance warning" paragraph verbatim.

Cost estimate:

- ~110 lines Python (builder method — heavier than
  NP-length because of the per-key auto-warning logic
  but otherwise identical shape)
- ~180 lines tests (mirror NP-length implementation
  tests + 1-2 new tests for the per-key auto-warning)
- ~30 lines docstrings + audit-doc-flip
- **0 lines manifest changes** (the block already
  exists)
- **0 lines engine / bridge / spec changes**

Lightest of the empirical-distribution estimator
slices.

---

## 11. Test surface — what this audit pins

Mirrored in
[`tests/test_p_nucleotide_length_estimation_contract.py`](../tests/test_p_nucleotide_length_estimation_contract.py).

### `pin_scaffold_*` — AIRR P-length fields (live)

1. `AirrRecord.p_v_3_length` / `p_d_5_length` /
   `p_d_3_length` / `p_j_5_length` exist as `i64` fields.
2. Python AIRR projection at `_airr_record.py` populates
   all four fields.
3. `result.py`'s canonical column order includes all
   four `p_*_length` columns.

### `pin_scaffold_*` — typed plane + spec (live)

4. `ReferenceEmpiricalModels.p_nucleotide_lengths` plane
   exists.
5. `P_NUCLEOTIDE_END_KEYS` and `P_NUCLEOTIDE_END_KEYS_VJ`
   constants hold.
6. Validator rejects unknown keys.
7. Validator rejects D-end keys on VJ at attach time
   (symmetric with `trims`).
8. `EmpiricalDistributionSpec` validator rejects
   negative values + non-positive weights (reused
   pins).

### `pin_scaffold_*` — bridge + engine surface (live)

9. `extract_recombine_defaults` returns all four
   `p_*_lengths` keys.
10. `_p_nucleotide_lengths_from_models` resolver exists.
11. Engine `PAdditionPass.parameter_signature` folds the
    length distribution via `fmt_int_dist`.

### `pin_present_*` — stop-and-report verification

12. Baseline (no P-plane) cartridges produce 0/50
    records with any non-zero P-length field — the
    estimator running against such data has no signal.
13. Authored P-plane produces nonzero P-length fields
    at meaningful rates (≥ 10/50 per end with a
    50/30/20 V_3 distribution).
14. Plan signature changes when the P-length
    distribution changes (no soft gap).

### `pin_present_*` — legacy orphan boundary

15. `p_nucleotide_length_probs` listed in
    `_DOCUMENTED_ORPHAN_DATACONFIG_FIELDS`.
16. Manifest's `legacy_p_nucleotide_length_probs_present`
    flag reads the orphan field's presence without
    lifting it.

### `pin_scaffold_*` — manifest already shaped

17. `manifest['models']['p_nucleotide_models']` already
    carries `length_keys` / `supported_ends` /
    `legacy_fallback` / `in_plan_signature` / etc. —
    the estimator's output rides the existing block.

### `pin_scaffold_*` — builder shape (reused)

18. `ConfigInfo.has_d` chain-type classifier.
19. Builder stage entry shape.
20. Builder idempotency pattern.
21. `csv.DictReader` AIRR-TSV ingestion.

### `pin_absence_*` — gaps the implementation slice closes

22. `ReferenceCartridgeBuilder.estimate_p_nucleotide_lengths`
    is not a method today.
23. No sibling module / class with the estimator's name.
24. No junction-arithmetic / NP-string heuristic
    inference path today.

### Doc anchor

25. Audit doc exists and references the contract file.

---

## 12. Out of scope

Documented here so a future implementer doesn't
accidentally expand the work.

- **Heuristic inference from junction / NP strings.** Per
  user brief §5 — v1 uses direct fields only. A future
  slice could attempt P-vs-N disambiguation but that
  requires modelling decisions out of scope for the
  builder-method shape.
- **Auto-lift of legacy `p_nucleotide_length_probs`.**
  Stays deferred — same boundary every prior estimator
  slice respected.
- **Per-allele / per-trim conditional P-length
  estimation.** The typed plane is per-end — flat by
  design. A future slice could expose per-allele
  estimation but would need a new typed-spec dataclass
  (the existing `EmpiricalDistributionSpec` is flat).
- **Per-base P-sequence (`p_v_3`, `p_d_5`, …)
  estimation.** v1 P-nucleotide design ships
  lengths-only; per-base P strings are deferred at the
  engine layer (see `docs/p_nucleotide_design.md` §15).
  Estimating them would require the engine surface to
  land first.
- **`Experiment.recombine(p_*_lengths=...)` kwargs.** Not
  exposed at v1 (P-length distribution is a
  cartridge-only surface). The estimator's output flows
  through the cartridge plane, not a per-experiment
  override.
- **Pandas DataFrame as first-class input.** v1 accepts
  list-of-dict / path / file handle.

---

## 13. Summary table

| Concern | Today | Recommendation |
|---|---|---|
| AIRR fields populated by GenAIRR | Four direct `p_*_length` integer fields populated when a typed P-plane is authored; hard-zero on every record without a P-plane | **Estimator consumes the four direct integer fields**; no junction arithmetic, no NP-string derivation |
| External AIRR data | **Does not carry these columns** — external tools (IgBLAST, MiXCR, …) do not model P-nucleotides | **Strong provenance warning in docstring** + per-key auto-warning when ≥ 95% of contributing rows reported zero |
| Typed plane | `ReferenceEmpiricalModels.p_nucleotide_lengths: Dict[str, EmpiricalDistributionSpec]` with `P_NUCLEOTIDE_END_KEYS = ("V_3","D_5","D_3","J_5")` + VJ-D-rejection at attach time | **Reuse verbatim.** No new spec class. |
| Bridge resolver | `_p_nucleotide_lengths_from_models` lives in `_dataconfig_extract.py`; precedence typed plane > None (no legacy fallback) already in place | **Reuse verbatim.** No bridge changes. |
| Engine surface | `PAdditionPass` consumes the length distribution; `parameter_signature` folds via `fmt_int_dist` | **Reuse verbatim.** No engine changes. |
| Provenance distinction | P lengths are distinct from NP lengths + trim lengths + junction at every layer (engine pass / trace addresses / AIRR fields / cartridge plane / DSL) | Estimator consumes ONLY the four direct `p_*_length` columns. Documented + pinned. |
| Chain-type enforcement | `metadata.has_d` is authoritative; spec validator chain-type-rejects D-end keys on VJ at attach time (stronger than `np_bases`) | VJ skips `p_d_5_length` / `p_d_3_length` columns + warns once per column; VDJ consumes all four |
| Builder method signature | n/a | `estimate_p_nucleotide_lengths(rearrangements, *, min_count=1, pseudocount=0.0, replace=True) → self` |
| Stage entry shape | n/a | Canonical `{stage, inputs, inferred, warnings}` |
| Manifest block | `p_nucleotide_models` ALREADY exists with the exact shape the estimator's output needs | **No new block, no extension** |
| Plan signature folding | P-length distributions already fold via `fmt_int_dist` | **No soft gap inherited** |
| Pre-flight stop-and-report | **NOT triggered** — direct fields exist and behave correctly; the only caveat is **narrow utility** for external AIRR data, handled via doc-warning + per-key auto-warning | Proceed to implementation slice |

The slice is **the simplest light-utility slice yet on
the plumbing side** — zero new manifest helpers, zero
engine changes, zero bridge changes, zero spec changes;
the only new code is the builder method body itself
(~110 lines including the per-key auto-warning). Total
estimated cost ~110 lines Python + ~180 lines tests.
The utility caveat is handled at the docstring +
runtime warning level.
