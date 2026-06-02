# Trim Distribution Estimation — Pre-Implementation Audit

**Status: audit only.** Designs the second
`ReferenceCartridgeBuilder.estimate_*` method —
`estimate_trim_distributions` — from observed AIRR
rearrangement data. Per the cartridge-authoring audit's
§11.2 recommended ordering, trim distributions are the
next-most-impactful data-derived plane after allele usage
because every downstream biology slice (NP length,
P-nucleotide, productive-only composition) depends on the
post-trim coding-flank state.

Companion to
[`tests/test_trim_distribution_estimation_contract.py`](../tests/test_trim_distribution_estimation_contract.py)
which freezes (a) the engine's existing trim-distribution
lowering chain that the estimator's output flows into,
(b) the AIRR trim fields GenAIRR populates reliably today,
(c) the end-loss vs recombination-trim boundary, and (d)
the builder/method absence pre-slice.

**Pre-flight finding (§7 below): clean-yes — no
stop-and-report condition.** GenAIRR populates exactly the
four AIRR trim fields the estimator needs (`v_trim_3`,
`d_trim_5`, `d_trim_3`, `j_trim_5`); `v_trim_5` and
`j_trim_3` are hard-zero because `Experiment.recombine`
does NOT trim those ends (confirmed at
[`engine_rs/src/airr_record/builder.rs:82,111`](../engine_rs/src/airr_record/builder.rs#L82)
+ [`src/GenAIRR/_airr_record.py:585`](../src/GenAIRR/_airr_record.py#L585)).
The typed `ReferenceEmpiricalModels.trims` plane already
exists with exactly the right key vocabulary
(`V_3` / `D_5` / `D_3` / `J_5`) and folds into the plan
signature via [`engine_rs/src/passes/trim.rs:281`](../engine_rs/src/passes/trim.rs#L281).
The end-loss surface (`corrupt.end_loss.5` /
`corrupt.end_loss.3` → AIRR `end_loss_5_length` /
`end_loss_3_length`) is rigorously separated from the
recombination-stage trims at every layer.

Smoke verification across 100 records on the bundled
cartridges:

| Cartridge | `v_trim_3` | `d_trim_5` | `d_trim_3` | `j_trim_5` | `v_trim_5` | `j_trim_3` |
|---|---|---|---|---|---|---|
| `HUMAN_IGH_OGRDB` (VDJ) | 71% nonzero | 80% nonzero | 82% nonzero | 86% nonzero | 0/100 | 0/100 |
| `HUMAN_IGK_OGRDB` (VJ) | 83% nonzero | 0/100 | 0/100 | 100% nonzero | 0/100 | 0/100 |

The estimator can land without disturbing any existing
surface; trim distributions become the second plane after
allele usage that recombination consumes by default.

---

## 1. Q1 — Input columns

### 1.1 AIRR trim fields available in the spec

The AIRR-C rearrangement schema defines six trim fields:

```text
v_trim_5, v_trim_3
d_trim_5, d_trim_3
j_trim_5, j_trim_3
```

### 1.2 Which fields GenAIRR actually populates

The engine populates the four AIRR fields corresponding to
the recombination-stage exonuclease passes, and hard-zeros
the two that no current `Experiment.recombine` pass touches:

| AIRR field | Sourced from | Engine pass |
|---|---|---|
| `v_trim_5` | hard `0` ([`builder.rs:82`](../engine_rs/src/airr_record/builder.rs#L82)) | — (no V_5 trim pass in v1) |
| `v_trim_3` | `trace["trim.v_3"]` ([`builder.rs:83`](../engine_rs/src/airr_record/builder.rs#L83)) | `TrimPass(V, Three)` |
| `d_trim_5` | `trace["trim.d_5"]` | `TrimPass(D, Five)` |
| `d_trim_3` | `trace["trim.d_3"]` | `TrimPass(D, Three)` |
| `j_trim_5` | `trace["trim.j_5"]` | `TrimPass(J, Five)` |
| `j_trim_3` | hard `0` ([`builder.rs:111`](../engine_rs/src/airr_record/builder.rs#L111)) | — (no J_3 trim pass in v1) |

The Python `result.py:_AIRR_FIELD_ORDER` block
(lines 83 / 84 / 96 / 97 / 109 / 110) declares all six
fields, and the legacy Python projection at
[`_airr_record.py:561-587`](../src/GenAIRR/_airr_record.py#L561)
hard-zeros `v_trim_5` with the explanatory comment "current
Experiment.recombine doesn't trim V_5".

### 1.3 Which fields the estimator consumes

The estimator consumes exactly the four trim fields the
engine populates — the same four `ReferenceEmpiricalModels.trims`
already accepts (`V_3` / `D_5` / `D_3` / `J_5`). The two
hard-zero AIRR fields (`v_trim_5`, `j_trim_3`) are
**not consumed** even when present in user-supplied AIRR
data:

| AIRR column | Estimator consumes | Reason |
|---|---|---|
| `v_trim_5` | NO | No `V_5` key on the typed plane; no engine pass; would silently drop. |
| `v_trim_3` | YES → key `"V_3"` | Direct map. |
| `d_trim_5` | YES → key `"D_5"` | Direct map (VDJ only). |
| `d_trim_3` | YES → key `"D_3"` | Direct map (VDJ only). |
| `j_trim_5` | YES → key `"J_5"` | Direct map. |
| `j_trim_3` | NO | No `J_3` key on the typed plane; no engine pass. |

Surveys passing AIRR data that does carry non-zero
`v_trim_5` / `j_trim_3` (e.g. records produced by another
simulator with a V_5 trim model) record the dropped
columns in the stage entry's `warnings` list with a single
notice per column — see §6.2.

### 1.4 Pinned

- `pin_scaffold_airr_record_carries_six_trim_fields`
- `pin_scaffold_engine_populates_recombination_trim_fields_only`
- `pin_scaffold_v_trim_5_and_j_trim_3_are_hard_zero`
- `pin_present_genairr_populates_trim_fields_reliably`
  (stop-and-report verification — see §7)

---

## 2. Q2 — Target model surface

### 2.1 Existing typed plane

The `ReferenceEmpiricalModels.trims` plane is the existing
typed surface for trim distributions. From
[`src/GenAIRR/reference_models.py:43`](../src/GenAIRR/reference_models.py#L43):

```python
TRIM_KEYS: Tuple[str, ...] = ("V_3", "D_5", "D_3", "J_5")
TRIM_KEYS_VJ: Tuple[str, ...] = ("V_3", "J_5")
```

The plane is a `Dict[str, EmpiricalDistributionSpec]` where
each value is a flat `[(value, weight), ...]` list. The
existing validator
([`reference_models.py:507-523`](../src/GenAIRR/reference_models.py#L507)):

1. Rejects keys outside `TRIM_KEYS`.
2. Rejects D keys (`D_5` / `D_3`) on VJ cartridges.
3. Delegates per-spec validation (non-negative ints,
   positive finite weights) to
   `EmpiricalDistributionSpec.validate()`.

### 2.2 Why not V_5 / J_3 — confirmation

Both `V_5` and `J_3` are deliberately absent from the typed
plane because:

1. **No engine pass.** The pipeline contains four
   `TrimPass` instances — `(V, Three)` / `(D, Five)` /
   `(D, Three)` / `(J, Five)` — and no others. Adding
   `V_5` / `J_3` keys to the typed plane would introduce
   plane state with no consumer.
2. **No AIRR trace address.** The trace address vocabulary
   is `trim.v_3` / `trim.d_5` / `trim.d_3` / `trim.j_5`;
   no `trim.v_5` / `trim.j_3` address exists at any layer.
3. **Hard zero in projection.** The Rust AIRR builder
   hard-zeros `rec.v_trim_5` and `rec.j_trim_3` regardless
   of trace contents (`builder.rs:82,111`).
4. **Recombination biology.** The two omitted ends
   correspond to allele 5'-end (V_5) and 3'-end (J_3) — in
   B/T cell receptor biology these are anchored to the
   conserved Cys (V) / W-or-F (J) and recombinase
   exonuclease activity at those ends is rare enough that
   the v1 pipeline omits modelling it. Future biology
   work could add V_5 / J_3 trim passes; that would be a
   separate audit covering trace-address expansion + AIRR
   builder un-zeroing + plan-signature folding.

### 2.3 Recombination-stage trims vs end-loss — confirmed separated

The end-loss surface is a completely separate corruption-stage
pass that runs after recombination + assembly:

| Concern | Recombination trim | End loss |
|---|---|---|
| Engine pass | `engine_rs/src/passes/trim.rs::TrimPass` | `engine_rs/src/passes/corrupt/end_loss.rs::EndLossPass` |
| Trace address | `trim.v_3` / `trim.d_5` / `trim.d_3` / `trim.j_5` | `corrupt.end_loss.5` / `corrupt.end_loss.3` |
| AIRR field | `v_trim_3` / `d_trim_5` / `d_trim_3` / `j_trim_5` | `end_loss_5_length` / `end_loss_3_length` |
| Cartridge plane | `ReferenceEmpiricalModels.trims` | n/a — per-experiment DSL kwarg only |
| DSL | `Experiment.recombine(v_3_lengths=...)` etc. + cartridge default | `Experiment.end_loss_5prime(length=...)` / `primer_trim_5prime(...)` alias |
| Biology stage | Recombination (germline exonuclease) | Sequencing / library-prep observation (5'/3' adapter or primer loss) |
| Documentation | This audit | [`docs/primer_trim_end_loss_audit.md`](primer_trim_end_loss_audit.md) §6.1 |

The Rust `AirrRecord` struct's docstring at
[`record.rs:204-207`](../engine_rs/src/airr_record/record.rs#L204)
explicitly distinguishes `end_loss_5_length` from
`v_trim_5`, naming the latter as the
"recombination-stage exonuclease trim".

### 2.4 Pinned

- `pin_scaffold_reference_empirical_models_trims_plane_exists`
- `pin_scaffold_trim_keys_constant_holds`
- `pin_scaffold_trim_keys_vj_constant_holds`
- `pin_scaffold_trims_validator_rejects_unknown_keys`
- `pin_scaffold_trims_validator_rejects_d_keys_on_vj`
- `pin_scaffold_end_loss_pass_distinct_from_trim_pass`
- `pin_scaffold_end_loss_airr_fields_distinct_from_trim_fields`

---

## 3. Q3 — Chain behaviour

### 3.1 Chain-type-driven enforcement

| Chain type | Estimator consumes | Skipped (warn) |
|---|---|---|
| VJ (`BCR_LIGHT_KAPPA` / `BCR_LIGHT_LAMBDA` / `TCR_ALPHA` / `TCR_GAMMA`) | `v_trim_3`, `j_trim_5` | `d_trim_5`, `d_trim_3` (absent on VJ — D-trim contributions ignored with a single one-time warning) |
| VDJ (`BCR_HEAVY` / `TCR_BETA` / `TCR_DELTA`) | `v_trim_3`, `d_trim_5`, `d_trim_3`, `j_trim_5` | — |

The cartridge's `metadata.has_d` (from `ConfigInfo`) is the
authoritative classifier — same source the allele-usage
estimator uses. The `EmpiricalDistributionSpec` plane is
chain-type-validated at attach time
([`reference_models.py:512-516`](../src/GenAIRR/reference_models.py#L512))
so an estimator that mistakenly writes a D key on a VJ
cartridge fails fast at `with_models()` / direct attach.

### 3.2 Comparison with allele-usage estimator

The allele-usage estimator's policy for the same boundary
is "ignore D contribution + warn once" on VJ; the trim
estimator mirrors that policy for consistency.

### 3.3 Pinned

- `pin_scaffold_config_info_has_d_drives_chain_classification`
  (existing pin reused from allele-usage contract)
- `pin_scaffold_vj_trim_keys_subset_of_vdj_trim_keys`

---

## 4. Q4 — Validation

### 4.1 Per-row validation rules

| Case | Default behaviour | Report destination |
|---|---|---|
| Required trim column missing (`v_trim_3` always; `d_trim_5` / `d_trim_3` on VDJ; `j_trim_5` always) | Skip row | `report.rejected` with `reason="missing_required_column"` |
| Required trim column non-integer (`"5.5"`, `"abc"`, `""`) | Skip row | `report.rejected` with `reason="malformed_trim_value"` (carries column + raw value) |
| Required trim column negative integer (`"-1"`) | Skip row | `report.rejected` with `reason="negative_trim_value"` |
| Optional trim column non-zero (`v_trim_5` / `j_trim_3` non-zero on a row) | Drop column contribution, keep row | `report.warnings` (cartridge-level, one notice per dropped column) |
| Below `min_count` threshold for a trim value's final count | Drop from result, surface count | `stage.inferred.below_min_count` (per-key dict) + cartridge-level warnings |

### 4.2 Pseudocount handling

The `pseudocount` kwarg (default `0.0`) adds a uniform
prior to every supported integer value in `[0, max_observed]`
per key. Use cases:

- `pseudocount=0.0` (default) — pure empirical estimation;
  unobserved trim values get probability 0.
- `pseudocount=1.0` — additive smoothing; every integer in
  `[0, max]` becomes drawable with at least nominal weight.
  Useful for small sample sizes where the cartridge's
  trim cap (see [`_dataconfig_extract._trim_cap`](../src/GenAIRR/_dataconfig_extract.py#L93))
  might leave the support sparse after clamping.

Pseudocount addition runs **before** the `min_count` filter
so the filter applies to combined `(observed + pseudocount)`
mass.

### 4.3 Recommended v1 — skip + report counts

Per the user brief: "skip row per missing required key
with structured report entry". The estimator emits a
single `inferred` block per stage entry:

```python
{
    "V_3": {<trim_int>: <normalised_weight>, ...},
    "D_5": {...},
    "D_3": {...},
    "J_5": {...},
    "skipped": {
        "missing_required_column": 14,
        "malformed_trim_value": {"v_trim_3": 2, "j_trim_5": 0, ...},
        "negative_trim_value": {"v_trim_3": 0, "d_trim_5": 1, ...},
    },
    "below_min_count": {"V_3": 0, "D_5": 3, "D_3": 2, "J_5": 0},
    "dropped_columns": {"v_trim_5": 23, "j_trim_3": 18},  # AIRR columns the estimator deliberately doesn't model
}
```

### 4.4 Pinned

- `pin_scaffold_empirical_distribution_spec_rejects_negative_values`
- `pin_scaffold_empirical_distribution_spec_rejects_non_positive_weights`
- `pin_absence_no_estimate_trim_distributions_min_count_kwarg_today`
- `pin_absence_no_estimate_trim_distributions_pseudocount_kwarg_today`

---

## 5. Q5 — Estimator API

### 5.1 Method signature

```python
def estimate_trim_distributions(
    self,
    rearrangements: Union[List[Dict[str, Any]], Path, IO[str]],
    *,
    min_count: float = 1.0,
    pseudocount: float = 0.0,
    replace: bool = True,
) -> "ReferenceCartridgeBuilder":
    ...
```

Returns `self` for fluent chaining. Idempotent — calling
twice re-estimates and overwrites the previous output with
`replaced=True` on the stage entry (same discipline
`infer_v_subregions` and `estimate_allele_usage` use).

### 5.2 Input shape (mirrors allele-usage)

Three accepted forms, mirroring the allele-usage
estimator's convention:

| Shape | Source |
|---|---|
| `list[dict]` (each dict is one rearrangement record) | Native Python pipeline output; MCP server input. |
| `pathlib.Path` / `str` (filesystem path) | AIRR TSV file. |
| open text file handle | In-memory or pre-opened stream. |

Path / file inputs use `csv.DictReader` with
`delimiter='\t'`. **DataFrame inputs are NOT first-class**
in v1 — users convert via `df.to_dict(orient="records")`.

### 5.3 Output — writes into the existing typed plane

The estimator constructs per-key
`EmpiricalDistributionSpec` instances and writes them into
`self._reference_models.trims`. Behaviour matches the
allele-usage discipline:

- If `_reference_models` is `None`, a fresh
  `ReferenceEmpiricalModels(trims={...})` is constructed.
- If `_reference_models` exists, only the affected `trims`
  keys are overwritten — `np_lengths` / `np_bases` /
  `p_nucleotide_lengths` / `allele_usage` are preserved.
- The bridge's existing
  `_trim_from_models` resolver
  ([`_dataconfig_extract.py:477-486`](../src/GenAIRR/_dataconfig_extract.py#L477))
  picks up the typed plane verbatim; the
  `Experiment.recombine(v_3_lengths=...)` kwarg still takes
  precedence (same kwarg > plane > legacy fallback chain
  the trim plane already has today).

### 5.4 Precedence — clean, already in place

The trim distribution lowering chain already implements
the precedence the estimator slice needs. From
[`_dataconfig_extract.py:253-260`](../src/GenAIRR/_dataconfig_extract.py#L253):

```python
"trim_v_3": _trim_from_models(explicit, "V_3", cap=cap_v)
    or extract_trim_distribution(cfg, "V", "3", cap=cap_v),
# ...
```

Resolution order: **(1) typed plane** (highest priority) →
**(2) legacy nested-dict `cfg.trim_dicts`** → **(3) None**
(uniform fallback at the engine layer).

The per-experiment override surface is
[`Experiment.trim(v_3=..., d_5=..., d_3=..., j_5=..., enabled=True)`](../src/GenAIRR/experiment.py#L1509)
— chained AFTER `recombine()` and BEFORE any
mutation / corruption step. The full precedence chain
the estimator inherits:

| Priority | Source |
|---|---|
| 1 (highest) | Explicit `Experiment.trim(v_3=..., …)` after `recombine()` |
| 2 | Typed `ReferenceEmpiricalModels.trims[key]` plane — **the estimator's output** |
| 3 | Legacy nested `DataConfig.trim_dicts` (auto-extracted by `extract_trim_distribution`) |
| 4 (lowest) | Engine no-op / uniform fallback when nothing else is supplied |

`Experiment.recombine` itself does NOT take per-distribution
trim kwargs — overrides flow through `.trim()`. This is
asymmetric with the allele-usage path (which had to add a
precedence shim inside `recombine`) and is **simpler**: the
typed plane is already at priority 2 with no shim required.

**No engine changes required.** No new resolver, no new
precedence shim — the typed plane is already the highest
priority short of the explicit `.trim()` override surface.

### 5.5 Pinned

- `pin_absence_no_estimate_trim_distributions_method_today`
- `pin_scaffold_extract_recombine_defaults_consumes_trims_plane`
- `pin_scaffold_trim_from_models_resolver_exists`
- `pin_scaffold_experiment_recombine_accepts_v_3_lengths_kwarg`

---

## 6. Q6 — Reporting

### 6.1 Stage name + canonical shape

The stage entry follows the canonical `{stage, inputs,
inferred, warnings}` shape pinned by the release smoke
test + the allele-usage estimator. From the user brief —
stage name `"estimate_trim_distributions"`.

### 6.2 Stage entry shape

```python
{
    "stage": "estimate_trim_distributions",
    "inputs": {
        "record_count": 18412,
        "min_count": 1.0,
        "pseudocount": 0.0,
        "source": "<path or 'records'>",
        "replaced": False,
    },
    "inferred": {
        "V_3": {0: 0.3214, 1: 0.2103, 2: 0.1574, ...},
        "D_5": {...},  # empty {} on VJ cartridges
        "D_3": {...},  # empty {} on VJ cartridges
        "J_5": {...},
        "skipped": {
            "missing_required_column": 14,
            "malformed_trim_value": {"v_trim_3": 2, ...},
            "negative_trim_value": {"v_trim_3": 0, ...},
        },
        "below_min_count": {"V_3": 0, "D_5": 3, "D_3": 2, "J_5": 0},
        "dropped_columns": {"v_trim_5": 23, "j_trim_3": 18},
    },
    "warnings": [
        "v_trim_5 column non-zero in 23 record(s); contribution dropped (no V_5 trim pass)",
        "j_trim_3 column non-zero in 18 record(s); contribution dropped (no J_3 trim pass)",
        "VJ cartridge: d_trim_5 / d_trim_3 columns present in records; contribution ignored",
        # ...
    ],
}
```

### 6.3 Per-segment normalisation

Per-key weights are normalised to sum to 1.0. Unlike the
allele-usage estimator (which writes a `{name: weight}`
dict directly), trim weights flow through an
`EmpiricalDistributionSpec([(value, weight), ...])`
container — the bridge's existing `_trim_from_models`
resolver flattens to the `[(int, float), ...]` engine
shape via the spec's `.values` field.

### 6.4 Pinned

- `pin_scaffold_builder_stage_entries_have_inputs_inferred_warnings`
  (existing pin reused from allele-usage contract)
- `pin_scaffold_idempotency_pattern_via_replaced_flag_in_v_subregions`
  (existing pin reused)

---

## 7. Q7 — Stop-and-report verification

### 7.1 The brief's stop-and-report condition

> "If the audit finds GenAIRR does not reliably populate
> recombination trim fields, stop and report. Otherwise
> implementation should be straightforward."

### 7.2 Verdict — NOT triggered. Proceed.

Empirical verification on the two bundled cartridges
(`HUMAN_IGH_OGRDB` + `HUMAN_IGK_OGRDB`) over 100 records
each at fixed seed:

| Cartridge | `v_trim_3` | `d_trim_5` | `d_trim_3` | `j_trim_5` |
|---|---|---|---|---|
| VDJ (`HUMAN_IGH_OGRDB`) | 71 / 100 nonzero | 80 / 100 nonzero | 82 / 100 nonzero | 86 / 100 nonzero |
| VJ (`HUMAN_IGK_OGRDB`) | 83 / 100 nonzero | n/a (D absent) | n/a (D absent) | 100 / 100 nonzero |

Every cartridge-relevant trim field is populated by a
significant fraction of records; the hard-zero fields
(`v_trim_5`, `j_trim_3`) match the documented engine
boundary (no V_5 / J_3 trim pass exists in v1).

### 7.3 No silent corruption

Three structural guarantees the estimator inherits:

1. **Trim distributions already fold into the plan signature.**
   Empirical check at audit time: two cartridges differing
   only in `trims["V_3"]` produce **different** plan
   signatures. Unlike the allele-usage slice's documented
   soft gap 1, the trim plane has NO inherited gap to
   manage.
2. **Trim trace addresses match the typed-plane keys.**
   `trim.v_3` / `trim.d_5` / `trim.d_3` / `trim.j_5`
   precisely correspond to `V_3` / `D_5` / `D_3` / `J_5` —
   no rename / no ambiguity.
3. **Manifest already exposes a `trim_keys` block.**
   The bundled cartridges report `trim_keys=[]` (no typed
   plane attached) and `legacy_trim_dicts_present=True`
   (legacy fallback active). The estimator slice can
   extend the existing block per §8 below.

### 7.4 Pinned

- `pin_present_genairr_populates_trim_fields_reliably`
  (stop-and-report verification: bundled cartridges emit
  the four trim fields at significant rates over a
  100-record smoke at fixed seed)
- `pin_present_trim_distributions_fold_into_plan_signature`
  (no soft gap inherited)

---

## 8. Manifest exposure

### 8.1 Today's `manifest['models']` block

Already exposes (per
[`dataconfig/data_config.py`](../src/GenAIRR/dataconfig/data_config.py)):

```python
{
    "has_reference_models": True / False,
    "np_length_keys": [...],
    "trim_keys": [...],
    "legacy_np_lengths_present": True / False,
    "legacy_trim_dicts_present": True / False,
    "shm": {...},
    "np_base_models": {...},
    "p_nucleotide_models": {...},
    "allele_usage": {...},  # post-allele-usage-slice
}
```

`trim_keys` is the sorted list of `ReferenceEmpiricalModels.trims`
keys; `legacy_trim_dicts_present` reports whether the
legacy nested `cfg.trim_dicts` is non-empty.

### 8.2 Estimator slice manifest extension

The estimator slice extends the existing minimal
`trim_keys` / `legacy_trim_dicts_present` pair into a
richer `trims` block paralleling the `allele_usage` block:

```python
"trims": {
    "available": True / False,          # typed plane non-empty?
    "keys": ["V_3", "D_5", "D_3", "J_5"],  # canonical key vocabulary
    "populated_keys": ["V_3", "J_5"],   # subset actually written
    "legacy_trim_dicts_present": True / False,  # legacy nested dict?
    "legacy_fallback": False,           # we never auto-lift the legacy dict
    "in_plan_signature": True,          # confirmed at audit time
    "source": "ReferenceEmpiricalModels.trims",
}
```

**Backwards compatibility:** the existing top-level
`trim_keys` / `legacy_trim_dicts_present` entries stay
in place. The new `trims` block is **additive**, not
a rename — same migration discipline the
allele-usage slice used.

### 8.3 Pinned

- `pin_scaffold_manifest_exposes_trim_keys_today`
- `pin_scaffold_manifest_exposes_legacy_trim_dicts_present_today`
- `pin_absence_no_trims_block_in_manifest_today`

---

## 9. Implementation order (recommended)

A single self-contained slice can land the estimator.
Mirrors the six-step allele-usage discipline, but
**lighter** because no new typed-spec class is needed
(`EmpiricalDistributionSpec` already exists):

1. **Builder method** — `ReferenceCartridgeBuilder.estimate_trim_distributions(
   rearrangements, *, min_count=1.0, pseudocount=0.0,
   replace=True)`. Reads input via `csv.DictReader` (path
   / file handle) or directly (records). Builds per-key
   integer Counter over `(v_trim_3 / d_trim_5 / d_trim_3
   / j_trim_5)` with chain-type-aware D-skip. Validates
   non-negative integers; rejects per-row with structured
   `report.rejected` entries. Applies `pseudocount`
   + `min_count` filter; normalises per key. Stage entry
   shape per §6.2. Constructs per-key
   `EmpiricalDistributionSpec` and writes into
   `self._reference_models.trims`.

2. **Manifest extension** — extend the existing
   `manifest['models']` block with the new `trims` sub-block
   per §8.2. The existing `trim_keys` /
   `legacy_trim_dicts_present` top-level entries stay
   for backwards compatibility.

3. **Tests** — implementation tests for the estimator
   (per-key smoke + chain-type enforcement + min_count
   + pseudocount + replay-through-built-cartridge
   end-to-end + manifest exposure + report shape).

4. **Doc + contract pin flip** — flip the absence pins in
   `tests/test_trim_distribution_estimation_contract.py`
   to present-state, mirror the allele-usage doc-flip
   discipline. Add the new `trims` block to the
   `_allele_usage_manifest_block` neighbour.

5. **Release smoke** — add
   `tests/test_trim_distribution_estimation_release.py`
   following the allele-usage release smoke pattern: tiny
   inline-FASTA VDJ cartridge → estimator on a synthetic
   AIRR record set → recombination shows the estimated
   bias (a record set heavily favouring trim=0 produces
   visibly shorter trims than legacy fallback).

6. **Validation matrix** — add a row to
   [`docs/validation_matrix.md`](validation_matrix.md)
   matching the allele-usage row's discipline.

Cost estimate:

- ~80 lines Python (builder method only — no new spec
  class, no resolver, no lowering: the precedence chain
  exists end-to-end already)
- ~30 lines manifest extension
- ~180 lines tests (implementation + release smoke)
- ~30 lines docstrings + audit-doc-flip

**Lighter than the allele-usage slice** because:

- No new typed-spec dataclass (estimator writes into
  existing `EmpiricalDistributionSpec`).
- No bridge resolver changes (precedence chain already
  prioritises the typed plane).
- No `Experiment.recombine` precedence shim
  (kwarg-over-plane already works through the existing
  `extract_recombine_defaults` resolver).

---

## 10. Test surface — what this audit pins

Mirrored in
[`tests/test_trim_distribution_estimation_contract.py`](../tests/test_trim_distribution_estimation_contract.py).

### `pin_scaffold_*` — live surfaces the estimator builds on

1. `AirrRecord` (Rust) carries all six AIRR trim fields
   (`v_trim_5` / `v_trim_3` / `d_trim_5` / `d_trim_3` /
   `j_trim_5` / `j_trim_3`).
2. The Rust AIRR builder hard-zeros `v_trim_5` and
   `j_trim_3`, sourcing the other four from
   `trace["trim.<seg>_<end>"]`.
3. The Python projection at `_airr_record.py:585` carries
   the same hard-zero comment for `v_trim_5`.
4. `ReferenceEmpiricalModels.trims` plane exists with key
   constants `TRIM_KEYS = ("V_3", "D_5", "D_3", "J_5")`
   and `TRIM_KEYS_VJ = ("V_3", "J_5")`.
5. `ReferenceEmpiricalModels` validator rejects unknown
   trim keys.
6. `ReferenceEmpiricalModels` validator rejects D-trim
   keys on VJ chains.
7. `EmpiricalDistributionSpec` validator rejects negative
   values + non-positive weights (existing surface the
   estimator uses verbatim).
8. `_dataconfig_extract.extract_recombine_defaults`
   consumes the typed `trims` plane and lowers it into
   `trim_v_3` / `trim_d_5` / `trim_d_3` / `trim_j_5`
   keys.
9. `_dataconfig_extract._trim_from_models` resolver exists
   as the typed-plane → engine-shape adapter.
10. `Experiment.trim(v_3=..., d_5=..., d_3=..., j_5=...,
    enabled=...)` is the per-experiment override surface
    chained AFTER `recombine()`.
11. Engine pass `TrimPass` (with trace addresses
    `trim.v_3` etc.) exists in `engine_rs/src/passes/trim.rs`
    and folds the distribution into the plan signature.
12. `ConfigInfo.has_d` is the authoritative chain-type
    classifier the estimator dispatches on.
13. Builder stage entry shape `{stage, inputs, inferred,
    warnings}` (existing pin from allele-usage contract).
14. Builder idempotency pattern via `replaced=True` flag.
15. `csv.DictReader` (stdlib only) handles AIRR TSV
    inputs.

### `pin_scaffold_*` — end-loss separation

16. `EndLossPass` lives at `engine_rs/src/passes/corrupt/end_loss.rs`
    (separate corruption-pass module).
17. `AirrRecord.end_loss_5_length` / `.end_loss_3_length`
    are distinct from `v_trim_5` / `v_trim_3` etc. (separate
    fields on the record struct).
18. AIRR record docstring at `record.rs:204-207` explicitly
    distinguishes `end_loss_5_length` from `v_trim_5`,
    naming the latter as "recombination-stage exonuclease
    trim".
19. `Experiment.end_loss_5prime` / `primer_trim_5prime`
    alias both wire to `EndLossPass`, NOT `TrimPass`.

### `pin_present_*` — stop-and-report verification

20. `Experiment.on('human_igh').recombine().run_records(n=100)`
    populates `v_trim_3` / `d_trim_5` / `d_trim_3` /
    `j_trim_5` at significant rates (≥ 50% nonzero on a
    fixed-seed smoke).
21. `Experiment.on('human_igk').recombine().run_records(n=100)`
    populates `v_trim_3` / `j_trim_5` at significant rates
    (≥ 50% nonzero on a fixed-seed smoke); D fields are
    zero (correct — VJ chain).
22. `v_trim_5` / `j_trim_3` are zero on every record
    across both bundled cartridges (correct — no V_5 /
    J_3 trim pass).

### `pin_present_*` — documented surface state

23. Trim distributions fold into the plan signature
    (no soft gap inherited): two cartridges differing only
    in `trims["V_3"]` produce different plan signatures.

### `pin_absence_*` — gaps the implementation slice closes

24. `ReferenceCartridgeBuilder.estimate_trim_distributions`
    is not a method today.
25. No `min_count` / `pseudocount` kwarg surface exists
    on the builder today.
26. `cartridge_manifest()["models"]` exposes minimal
    `trim_keys` / `legacy_trim_dicts_present` entries
    but NOT a structured `trims` block.

### Doc anchor

27. Audit doc exists and references the contract file;
    section structure intact.

---

## 11. Out of scope

Documented here so a future implementer doesn't accidentally
expand the work.

- **V_5 / J_3 trim modelling.** The two omitted AIRR fields
  would require adding `TrimPass(V, Five)` and
  `TrimPass(J, Three)` passes to the recombination
  pipeline, extending the trace-address vocabulary,
  un-zeroing the AIRR builder, and folding two new
  distributions into the plan signature. Separate audit
  (biology question, not estimator scope).
- **Auto-lift of legacy `DataConfig.trim_dicts`.** The
  legacy nested-dict format
  (`{family: {gene: {trim: prob}}}`) is per-gene; the
  typed plane is per-segment-and-end (marginalized). Auto-
  lifting would silently change output bytes for the 106
  bundled cartridges. Stays deferred — same boundary the
  Markov / P-nucleotide / allele-usage slices respected
  for their legacy fields.
- **Per-gene / per-allele trim distributions.** The
  typed plane is per-(segment, end) — one distribution
  per key. The legacy `trim_dicts` shape is per-gene;
  a future slice could expose per-allele estimation but
  would need a new typed-spec dataclass (the existing
  `EmpiricalDistributionSpec` is flat by design). Out of
  scope for v1.
- **In-frame / out-of-frame conditioning.** The AIRR-C
  schema allows conditioning trim distributions on the
  rearrangement's productivity. v1 ignores `productive`
  / `vj_in_frame` / `stop_codon` columns; the estimator
  consumes every record uniformly. Stratification would be
  a separate audit covering the joint distribution shape.
- **Pandas DataFrame as first-class input.** v1 accepts
  list-of-dict / path / file handle.
- **End-loss estimation.** Separate audit
  ([`docs/primer_trim_end_loss_audit.md`](primer_trim_end_loss_audit.md))
  and separate biology stage (observation / library prep,
  not recombination).

---

## 12. Summary table

| Concern | Today | Recommendation |
|---|---|---|
| AIRR trim fields populated by GenAIRR | `v_trim_3` / `d_trim_5` / `d_trim_3` / `j_trim_5` populated reliably; `v_trim_5` / `j_trim_3` hard-zero by design (no V_5 / J_3 trim pass) | **Estimator consumes the four populated fields**; surfaces non-zero `v_trim_5` / `j_trim_3` in user data as `dropped_columns` warnings |
| Typed plane | `ReferenceEmpiricalModels.trims: Dict[str, EmpiricalDistributionSpec]` with `TRIM_KEYS = ("V_3","D_5","D_3","J_5")` + VJ-rejects-D-keys validator | **Reuse verbatim.** No new spec class. |
| Bridge resolver | `_trim_from_models` lives in `_dataconfig_extract.py`; precedence kwarg > typed plane > legacy `trim_dicts` > uniform already in place | **Reuse verbatim.** No bridge changes. |
| Engine surface | `TrimPass(segment, end)` for each of the four `(V,Three)` / `(D,Five)` / `(D,Three)` / `(J,Five)` pairs; trace addresses `trim.v_3` etc.; folds into plan signature via `fmt_int_dist` | **Reuse verbatim.** No engine changes. |
| End-loss separation | Completely distinct surface: `EndLossPass` corruption-stage pass + `end_loss_5_length` / `end_loss_3_length` AIRR fields + `corrupt.end_loss.5` / `.3` trace addresses | **Estimator MUST NOT consume `end_loss_*` AIRR fields.** Documented + pinned. |
| Chain type enforcement | `metadata.has_d` is the authoritative classifier; VJ rejects D-trim keys at validator | VJ skips `d_trim_*` columns + warns once; VDJ consumes all four |
| Builder method signature | n/a | `estimate_trim_distributions(rearrangements, *, min_count=1.0, pseudocount=0.0, replace=True) → self` |
| Stage entry shape | n/a | Canonical `{stage, inputs, inferred, warnings}` per §6.2 |
| Manifest block | `trim_keys: []` + `legacy_trim_dicts_present: True` exist today; no structured `trims` sub-block | New `models.trims` block per §8.2 (additive; existing entries stay) |
| Plan signature folding | Trim distributions already fold via `fmt_int_dist` | **No soft gap inherited** — unlike the allele-usage slice |
| Pre-flight stop-and-report | **NOT triggered** — bundled cartridges populate the four trim fields at significant rates; the two hard-zero fields match the documented engine boundary | Proceed to implementation slice |

The estimator is **architecturally trivial**: the typed
plane is already present, the bridge precedence chain
already prioritises it, the engine surface already
consumes it, the plan signature already folds it. The
only missing piece is a single builder method that
populates `_reference_models.trims` from observed AIRR
records following the user-brief signature. The slice
should be substantially smaller than the allele-usage
slice by every metric (no new spec, no resolver changes,
no precedence shim, no soft gap inheritance).
