# NP Length Distribution Estimation — Pre-Implementation Audit

**Status: audit only.** Designs the third
`ReferenceCartridgeBuilder.estimate_*` method —
`estimate_np_length_distributions` — from observed AIRR
rearrangement data. Per the cartridge-authoring audit's
§11.2 recommended ordering, NP length distributions are the
next data-derived plane after trim distributions.

Companion to
[`tests/test_np_length_estimation_contract.py`](../tests/test_np_length_estimation_contract.py)
which freezes (a) the AIRR fields GenAIRR populates for NP
lengths, (b) the typed `ReferenceEmpiricalModels.np_lengths`
plane the estimator writes into, (c) the engine surface
that already consumes the plane and folds it into the plan
signature, (d) the VJ chain boundary (no NP2 region), and
(e) the builder/method absence pre-slice.

**Pre-flight finding (§7 below): clean-yes — no
stop-and-report condition.** GenAIRR populates direct
integer `np1_length` and `np2_length` AIRR fields at high
rates on the bundled cartridges, the typed plane already
exists with the right key vocabulary, the bridge precedence
chain already prioritises it, and the engine pass folds
the length distribution into the plan signature. Empirical
smoke at audit time (100 records per cartridge, fixed seed):

| Cartridge | `np1_length` | `np2_length` |
|---|---|---|
| `HUMAN_IGH_OGRDB` (VDJ) | 98 / 100 nonzero | 93 / 100 nonzero |
| `HUMAN_IGK_OGRDB` (VJ)  | 100 / 100 nonzero | 0 / 100 (no NP2 on VJ) |

Pattern shape is identical to trim estimation: builder
method + report + manifest block, no new spec class, no
engine changes, no bridge changes. The slice is the
lightest one yet because there is no asymmetry in field
shape (single integer per record per key) and the legacy
`cfg.NP_lengths` orphan field stays untouched.

---

## 1. Q1 — Input columns

### 1.1 AIRR fields exposed today

The Rust `AirrRecord` struct at
[`engine_rs/src/airr_record/record.rs:77-82`](../engine_rs/src/airr_record/record.rs#L77)
exposes both the per-NP-region sequence and direct length
fields:

| Rust field | AIRR column | Shape |
|---|---|---|
| `np1: String` | `np1` | nucleotide sequence of the unclaimed NP1 span |
| `np1_aa: String` | `np1_aa` | codon-rail translation of `np1` |
| `np1_length: i64` | `np1_length` | length of the NP1 region (post-claim-reabsorption) |
| `np2: String` | `np2` | sequence of the unclaimed NP2 span |
| `np2_aa: String` | `np2_aa` | codon-rail translation of `np2` |
| `np2_length: i64` | `np2_length` | length of the NP2 region |

The Python projection at
[`src/GenAIRR/_airr_record.py:714,719`](../src/GenAIRR/_airr_record.py#L714)
mirrors the Rust shape — `np1_length` /`np2_length` are
derived from region `(end - start)`, with `0` when the
region is absent (VJ cartridges have no NP2 region).

The Python `result.py` column order at lines 121 / 124
declares `np1_length` / `np2_length` as canonical AIRR
fields exported on every cartridge.

### 1.2 Estimator's chosen source — direct length fields

Per user brief §5: "v1 should use direct fields only". The
estimator consumes the **direct integer length fields**:

| Chain | Consumed AIRR columns | Plane key |
|---|---|---|
| VJ | `np1_length` | `NP1` |
| VDJ | `np1_length`, `np2_length` | `NP1`, `NP2` |

The sequence fields (`np1`, `np2`) are NOT consumed —
neither for length derivation nor for any other purpose.
This avoids the failure mode of disagreement between
`len(np1) != np1_length` (which could arise if a record
came from another simulator that interpreted "NP region"
differently, e.g. without post-claim reabsorption — see
[`engine_rs/src/airr_record/builder.rs:194-196`](../engine_rs/src/airr_record/builder.rs#L194)).

### 1.3 Fields the estimator does NOT consume

| AIRR column | Why ignored |
|---|---|
| `np1`, `np2` | Sequence-derived length is sensitive to the source simulator's NP-vs-D-tail attribution; direct integer fields are canonical. |
| `np1_aa`, `np2_aa` | Codon-rail translations — irrelevant to length estimation. |
| `p_v_3_length`, `p_d_5_length`, `p_d_3_length`, `p_j_5_length` | **Distinct biology stage.** P-nucleotide additions are templated palindromic extensions of coding flanks, NOT N-additions. See provenance distinction §5. |
| `junction_length` | Aggregate length over multiple regions; arithmetic decomposition is fragile across simulators. |

### 1.4 Pinned

- `pin_scaffold_airr_record_carries_np1_and_np2_length_fields`
- `pin_scaffold_airr_record_carries_np1_and_np2_sequence_fields`
- `pin_scaffold_python_airr_projection_emits_np1_length_np2_length`
- `pin_scaffold_result_column_order_includes_np_length_fields`

---

## 2. Q2 — Target model surface

### 2.1 Existing typed plane

The `ReferenceEmpiricalModels.np_lengths` plane is the
existing typed surface for NP-length distributions. From
[`src/GenAIRR/reference_models.py:42`](../src/GenAIRR/reference_models.py#L42):

```python
NP_KEYS: Tuple[str, ...] = ("NP1", "NP2")
```

The plane is `Dict[str, EmpiricalDistributionSpec]` —
same shape as `trims`. The existing validator
([`reference_models.py:494-504`](../src/GenAIRR/reference_models.py#L494)):

1. Rejects keys outside `NP_KEYS`.
2. Delegates per-spec validation (non-negative ints,
   positive finite weights) to
   `EmpiricalDistributionSpec.validate()`.

**Documented asymmetry vs `trims`:** the `np_lengths`
validator does NOT chain-type-reject `NP2` on a VJ
cartridge at attach time. The `trims` validator DOES
reject `D_5` / `D_3` keys on VJ. This is a pre-existing
asymmetry — the estimator should mirror the trim
estimator's pattern (skip NP2 contribution silently +
warn) rather than rely on the spec validator.

### 2.2 Why NP1 only on VJ chains

A VJ chain (no D segment) has a single NP-equivalent
region between V and J — the engine names it `NP1` and
omits `NP2` entirely. The Python AIRR projection emits
`np2_length=0` for every record on a VJ cartridge
(confirmed empirically: 100/100 records zero on
`HUMAN_IGK_OGRDB`). The estimator's VJ behaviour:

| Source `np2_length` on VJ row | Estimator behaviour |
|---|---|
| `0` (engine's own AIRR records) | Field-local: count as `(NP2, 0)` ? No. **Skip silently** — VJ cartridges have no NP2 plane key. |
| Non-zero (from a non-GenAIRR simulator with a different chain definition) | Skip the NP2 contribution + emit a single warning (mirrors the trim estimator's VJ D-column policy). |

The estimator writes ONLY `NP1` to a VJ cartridge's
`np_lengths` plane.

### 2.3 Pinned

- `pin_scaffold_reference_empirical_models_np_lengths_plane_exists`
- `pin_scaffold_np_keys_constant_holds`
- `pin_scaffold_np_lengths_validator_rejects_unknown_keys`
- `pin_present_np_lengths_validator_does_not_reject_np2_on_vj`
  (documented asymmetry vs `trims`)

---

## 3. Q3 — Chain behaviour

### 3.1 Chain-type-driven enforcement

| Chain | Estimator consumes | Estimator emits | Skipped (warn) |
|---|---|---|---|
| VJ | `np1_length` only | `NP1` key | `np2_length` ≠ 0 in source row (one-time warning, contribution dropped) |
| VDJ | `np1_length`, `np2_length` | `NP1`, `NP2` keys | — (both columns required per chain) |

The cartridge's `metadata.has_d` (from `ConfigInfo`) is the
authoritative classifier — same source the trim + allele-
usage estimators use.

### 3.2 Existing engine-side enforcement

`Experiment.recombine(np2_lengths=...)` already raises
`ValueError` on VJ chains:
[`src/GenAIRR/experiment.py:1378-1384`](../src/GenAIRR/experiment.py#L1378).
Message:

```text
np2_lengths is only valid for VDJ chains; the bound refdata is
{chain_type!r} (no D segment, no NP2 region). Drop the
np2_lengths kwarg or bind a VJ refdata.
```

The estimator inherits this safety net at compile time
for any cartridge it builds: a VJ cartridge with an
estimator-written `NP2` plane key would crash at
`Experiment.recombine()` not at `with_models()` /
direct attach. **The estimator MUST NOT write `NP2` on a
VJ cartridge** — the chain-type guard is enforced at
estimation time, not at lowering time.

### 3.3 Field-local skipping on VDJ

A VDJ row missing `np2_length` drops only the NP2
contribution — the NP1 column still feeds its
distribution. Same field-local discipline the trim
estimator follows.

### 3.4 Pinned

- `pin_scaffold_config_info_has_d_drives_chain_classification`
  (reused from trim contract)
- `pin_scaffold_recombine_rejects_np2_lengths_on_vj`
- `pin_scaffold_vj_np2_length_is_hard_zero_in_airr_projection`

---

## 4. Q4 — Validation

### 4.1 Per-row, per-field validation rules

Same shape as the trim estimator (field-local, structured
rejection):

| Case | Default behaviour | Report destination |
|---|---|---|
| Required NP-length column missing (`np1_length` always; `np2_length` on VDJ) | Skip THAT field's contribution; row's other field still feeds | `report.rejected` with `reason="missing_required_column"` (carries column) |
| Non-integer value (`"5.5"`, `"abc"`, `""`) | Skip THAT field's contribution | `report.rejected` with `reason="malformed_length_value"` |
| Negative integer (`"-1"`) | Skip THAT field's contribution | `report.rejected` with `reason="negative_length_value"` |
| VJ row has non-zero `np2_length` | Drop NP2 contribution; keep NP1 | `report.warnings` (cartridge-level, one notice per dataset) |
| Below `min_count` for a length value | Drop value from result | `stage.inferred.below_min_count` (per-key dict) + cartridge-level warning |

### 4.2 Estimator API kwargs

Mirror the trim estimator's signature so the user-facing
surface is consistent:

- `min_count: int = 1` — drop length values whose observed
  count is strictly below the threshold before
  normalisation.
- `pseudocount: float = 0.0` — additive smoothing applied
  to **observed** length values only (no support
  expansion). Same v1 simplicity as the trim estimator.
- `replace: bool = True` — idempotency. `False` blocks
  re-entry when typed-plane `np_lengths` is already
  attached.

### 4.3 Pinned

- `pin_scaffold_empirical_distribution_spec_rejects_negative_values`
  (reused from trim contract)
- `pin_scaffold_empirical_distribution_spec_rejects_non_positive_weights`
  (reused)
- `pin_absence_no_estimate_np_length_distributions_today`
- `pin_absence_no_min_count_pseudocount_kwarg_owned_by_np_length_estimator_today`

---

## 5. Q5 — Provenance distinction

### 5.1 NP lengths vs P-nucleotide lengths — different biology

| Concern | NP region | P-nucleotide |
|---|---|---|
| Biology stage | Recombination — TdT-like N-base addition between V/D/J coding ends | Recombination — palindromic complement of coding flank |
| Engine pass | `engine_rs/src/passes/generate_np.rs::GenerateNPPass` | `engine_rs/src/passes/p_addition.rs::PAdditionPass` |
| Cartridge plane | `ReferenceEmpiricalModels.np_lengths: Dict["NP1"\|"NP2", EmpiricalDistributionSpec]` | `ReferenceEmpiricalModels.p_nucleotide_lengths: Dict["V_3"\|"D_5"\|"D_3"\|"J_5", EmpiricalDistributionSpec]` |
| AIRR length field | `np1_length`, `np2_length` | `p_v_3_length`, `p_d_5_length`, `p_d_3_length`, `p_j_5_length` |
| Documentation | This audit | [`docs/p_nucleotide_design.md`](p_nucleotide_design.md) |

The two surfaces share the `EmpiricalDistributionSpec`
container shape but are **biologically distinct**. The NP-
length estimator MUST NOT consume any `p_*_length` AIRR
column.

### 5.2 NP lengths vs junction_length — do not derive

`junction_length` is the aggregate length spanning the V
Cys → J W/F + 3 — multiple regions (V residue + NP1 +
D + NP2 + J residue). Decomposing it back into NP1 / NP2
contributions requires knowing trim lengths + D length,
and the arithmetic differs between simulators that
attribute the V/J anchor residues differently to the
junction span. **The audit defers this approach.** v1
uses the direct `np1_length` / `np2_length` fields
exclusively.

### 5.3 NP lengths vs np1 / np2 sequence length — do not derive

The Rust builder's note at
[`builder.rs:194-200`](../engine_rs/src/airr_record/builder.rs#L194)
explains that NP regions in the IR are subject to
**post-claim reabsorption** — V/D right-extension claims
can reabsorb NP1 positions back into a structural
segment. This means:

- `np1_length` = the unclaimed-byte count after
  reabsorption (canonical for the estimator).
- `len(np1)` = the same after-reabsorption count on
  GenAIRR records (matches `np1_length`).

A non-GenAIRR record might have `len(np1) != np1_length`
if the producing simulator interpreted NP differently. To
avoid this ambiguity v1 reads only the direct integer
field.

### 5.4 Pinned

- `pin_scaffold_p_nucleotide_length_fields_distinct_from_np_length_fields`
- `pin_scaffold_p_addition_pass_module_distinct_from_generate_np`
- `pin_scaffold_p_nucleotide_lengths_plane_distinct_from_np_lengths`

---

## 6. Q6 — Engine consumption + lowering

### 6.1 Bridge resolver — already wired

The bridge's
[`_dataconfig_extract.extract_recombine_defaults`](../src/GenAIRR/_dataconfig_extract.py#L214-L215)
returns `np1` / `np2` keys with the typed-plane resolver
taking precedence over the legacy nested dict:

```python
"np1": _np_lengths_from_models(explicit, "NP1") or extract_np_lengths(cfg, "NP1"),
"np2": _np_lengths_from_models(explicit, "NP2") or extract_np_lengths(cfg, "NP2"),
```

`_np_lengths_from_models` exists at
[`_dataconfig_extract.py:311-328`](../src/GenAIRR/_dataconfig_extract.py#L311)
and lowers the typed `EmpiricalDistributionSpec` into the
engine's `[(int, float), ...]` shape.

### 6.2 Precedence chain — same as trims

| Priority | Source |
|---|---|
| 1 (highest) | Explicit `Experiment.recombine(np1_lengths=..., np2_lengths=...)` kwarg |
| 2 | Typed `ReferenceEmpiricalModels.np_lengths[key]` plane — **the estimator's output** |
| 3 | Legacy nested `DataConfig.NP_lengths` (auto-extracted by `extract_np_lengths`) |
| 4 (lowest) | Engine uniform fallback (raw-RefDataConfig path with a UserWarning) |

Asymmetry vs trim: `Experiment.recombine` itself takes
`np1_lengths` / `np2_lengths` kwargs (the trim DSL uses
`.trim(v_3=..., d_5=..., d_3=..., j_5=...)` chained
after `recombine`). The estimator inherits the precedence
chain verbatim — no shim required.

### 6.3 Plan signature folding — confirmed at audit time

Empirical check: two cartridges differing only in
`np_lengths["NP1"]` produce **different** plan signatures
at compile time. The `GenerateNPPass.parameter_signature()`
implementation at
[`engine_rs/src/passes/generate_np.rs:147-161`](../engine_rs/src/passes/generate_np.rs#L147)
folds the length distribution via `fmt_int_dist` — same
mechanism the trim plane uses. **No soft gap inherited**
(same property the trim estimator enjoys).

### 6.4 Pinned

- `pin_scaffold_extract_recombine_defaults_consumes_np_lengths_plane`
- `pin_scaffold_np_lengths_from_models_resolver_exists`
- `pin_scaffold_typed_plane_takes_precedence_over_legacy_np_lengths`
- `pin_scaffold_experiment_recombine_accepts_np1_lengths_and_np2_lengths_kwargs`
- `pin_scaffold_engine_generate_np_pass_folds_length_dist_into_signature`
- `pin_present_np_length_distribution_changes_plan_signature`

---

## 7. Q7 — Stop-and-report verification

### 7.1 The brief's stop-and-report condition

> "If direct `np1_length` / `np2_length` fields do not
> exist, stop and report with the recommended source field.
> Do not implement a heuristic estimator in the audit slice."

### 7.2 Verdict — NOT triggered. Proceed.

Direct integer `np1_length` / `np2_length` AIRR fields
exist on both the Rust `AirrRecord` struct and the Python
projection. Empirical verification on the two bundled
cartridges (100 records each at fixed seed):

| Cartridge | `np1_length` populated | `np2_length` populated |
|---|---|---|
| `HUMAN_IGH_OGRDB` (VDJ) | 98 / 100 nonzero | 93 / 100 nonzero |
| `HUMAN_IGK_OGRDB` (VJ) | 100 / 100 nonzero | 0 / 100 (correct — no NP2 on VJ) |

Every cartridge-relevant NP-length field is populated at
high rates; the hard-zero behaviour on VJ matches the
documented engine boundary (no NP2 region exists on a VJ
chain).

### 7.3 No silent corruption

Three structural guarantees the estimator inherits:

1. **NP-length distributions already fold into the plan
   signature.** Empirical check at audit time: two
   cartridges differing only in `np_lengths["NP1"]`
   produce different plan signatures. No soft gap to
   manage.
2. **Typed plane already takes precedence over the legacy
   nested dict.** Bridge plumbing identical to trims —
   nothing new to wire.
3. **Manifest already exposes `np_length_keys` + `legacy_np_lengths_present`.**
   The estimator slice extends this into a structured
   `np_length_models` block (additive, like the trim
   slice).

### 7.4 Pinned

- `pin_present_genairr_populates_np_length_fields_reliably`
- `pin_present_vj_np2_length_is_hard_zero`
- `pin_present_np_length_distributions_fold_into_plan_signature`

---

## 8. Manifest exposure

### 8.1 Today's `manifest['models']` block

Already exposes (per
[`dataconfig/data_config.py`](../src/GenAIRR/dataconfig/data_config.py)):

```python
{
    ...
    "np_length_keys": [...],          # sorted list of typed plane keys
    "legacy_np_lengths_present": True / False,
    ...
}
```

### 8.2 Estimator slice manifest extension

The estimator slice adds a structured `np_length_models`
sub-block paralleling the trim slice's `trim_models`:

```python
"np_length_models": {
    "keys": ["NP1", "NP2"],                  # populated subset of NP_KEYS
    "source": "ReferenceEmpiricalModels.np_lengths",
    "in_plan_signature": True,               # confirmed at audit time
    "legacy_np_lengths_present": True / False,
    "legacy_fallback": False,                # we never auto-lift the legacy dict
}
```

**Backwards compatibility:** the existing top-level
`np_length_keys` / `legacy_np_lengths_present` entries
stay in place. The new `np_length_models` block is
**additive** — same discipline the trim slice used.

### 8.3 Pinned

- `pin_scaffold_manifest_exposes_np_length_keys_today`
- `pin_scaffold_manifest_exposes_legacy_np_lengths_present_today`
- `pin_absence_no_np_length_models_block_in_manifest_today`

---

## 9. Implementation order (recommended)

Even lighter than the trim slice — no new field shape, no
chain-type asymmetry in field availability:

1. **Builder method** —
   `ReferenceCartridgeBuilder.estimate_np_length_distributions(
   rearrangements, *, min_count=1, pseudocount=0.0,
   replace=True)`. Reads `np1_length` always; reads
   `np2_length` only on VDJ. Field-local per-row
   validation: missing / malformed / negative values drop
   that field's contribution only. Constructs per-key
   `EmpiricalDistributionSpec` and writes into
   `self._reference_models.np_lengths`. Stage entry shape
   per §6.2 of trim doc (canonical
   `{stage, inputs, inferred, warnings}`).

2. **Manifest extension** — extend
   `manifest['models']` with the new `np_length_models`
   sub-block per §8.2. The existing `np_length_keys` /
   `legacy_np_lengths_present` top-level entries stay.

3. **Tests** — implementation tests for the estimator
   (per-key smoke + chain-type enforcement + min_count
   + pseudocount + replay-through-built-cartridge end-
   to-end + manifest exposure + report shape).

4. **Doc + contract pin flip** — flip the absence pins in
   `tests/test_np_length_estimation_contract.py` to
   present-state. Mirror the trim slice's lockstep
   updates to the reference-cartridge-authoring
   estimator boundary pins.

5. **Validation matrix** — add a row to
   [`docs/validation_matrix.md`](validation_matrix.md)
   matching the trim slice's discipline.

Cost estimate:

- ~70 lines Python (builder method only — even simpler
  than trim because no four-key fanout / dropped_columns
  bookkeeping)
- ~25 lines manifest extension
- ~140 lines implementation tests
- ~25 lines docstrings + audit-doc-flip

---

## 10. Test surface — what this audit pins

Mirrored in
[`tests/test_np_length_estimation_contract.py`](../tests/test_np_length_estimation_contract.py).

### `pin_scaffold_*` — AIRR fields (live)

1. `AirrRecord.np1_length` and `np2_length` exist as
   `i64` fields on the Rust struct.
2. `AirrRecord.np1` / `np1_aa` / `np2` / `np2_aa`
   exist as `String` fields (the sequence + AA columns
   the estimator does NOT consume).
3. Python projection at `_airr_record.py` populates
   `np1_length` / `np2_length` from region span.
4. `result.py`'s `_AIRR_FIELD_ORDER` declares both fields
   in canonical column order.

### `pin_scaffold_*` — typed plane (live)

5. `ReferenceEmpiricalModels.np_lengths` plane exists.
6. `NP_KEYS = ("NP1", "NP2")` constant holds.
7. The plane validator rejects keys outside `NP_KEYS`.
8. The plane validator does NOT reject `NP2` on a VJ
   cartridge — documented asymmetry vs `trims`.

### `pin_scaffold_*` — bridge + engine surface (live)

9. `_dataconfig_extract.extract_recombine_defaults`
   returns `np1` / `np2` keys.
10. `_np_lengths_from_models` resolver exists.
11. Typed plane takes precedence over legacy `cfg.NP_lengths`.
12. `Experiment.recombine` accepts `np1_lengths` and
    `np2_lengths` kwargs.
13. `Experiment.recombine(np2_lengths=...)` raises
    `ValueError` on VJ chains.
14. `engine_rs/src/passes/generate_np.rs::GenerateNPPass`
    folds the length distribution via `fmt_int_dist` in
    `parameter_signature`.

### `pin_scaffold_*` — provenance distinction

15. `AirrRecord.p_v_3_length` / `p_d_5_length` /
    `p_d_3_length` / `p_j_5_length` are SEPARATE fields
    from `np1_length` / `np2_length`.
16. `engine_rs/src/passes/p_addition.rs` exists as a
    separate engine pass module.
17. `ReferenceEmpiricalModels.p_nucleotide_lengths` is a
    separate plane from `np_lengths`.

### `pin_scaffold_*` — builder shape

18. `ConfigInfo.has_d` chain-type classifier (reused pin).
19. Builder stage entry shape `{stage, inputs, inferred,
    warnings}` (reused pin).
20. Builder idempotency pattern via `replaced=True` flag
    (reused pin).
21. `csv.DictReader` (stdlib only) handles AIRR TSV (reused).

### `pin_present_*` — stop-and-report verification

22. `Experiment.on('human_igh').recombine().run_records(n=100)`
    populates `np1_length` / `np2_length` at ≥ 50 / 100
    each. **Stop-and-report gate.**
23. `Experiment.on('human_igk').recombine().run_records(n=100)`
    populates `np1_length` at ≥ 50 / 100; `np2_length`
    is hard-zero on every record. **Stop-and-report gate.**

### `pin_present_*` — documented surface state

24. NP-length distributions fold into the plan signature
    (no soft gap inherited): two cartridges differing
    only in `np_lengths["NP1"]` produce different plan
    signatures.

### `pin_scaffold_*` — manifest current state

25. `manifest['models']` exposes `np_length_keys` and
    `legacy_np_lengths_present` entries today.

### `pin_absence_*` — gaps the implementation slice closes

26. `ReferenceCartridgeBuilder.estimate_np_length_distributions`
    is not a method today.
27. No sibling `min_count` / `pseudocount` kwarg owned by
    an NP-length estimator method (both currently owned
    by the trim estimator + allele-usage).
28. `manifest['models']` does NOT yet expose a structured
    `np_length_models` sub-block.
29. No sibling module / class / free function with the
    estimator's name was introduced.

### Doc anchor

30. Audit doc exists and references the contract file.

---

## 11. Out of scope

Documented here so a future implementer doesn't accidentally
expand the work.

- **Derivation from junction_length.** Decomposing the
  aggregate junction length into NP1/NP2 contributions
  requires multi-region arithmetic that varies across
  simulators. v1 uses direct `np1_length` / `np2_length`
  fields exclusively.
- **Derivation from `np1` / `np2` sequence length.** Per
  audit §5.3, post-claim reabsorption makes the sequence
  length canonical only on GenAIRR records. v1 reads only
  the direct integer field.
- **Auto-lift of legacy `DataConfig.NP_lengths`.** Stays
  deferred — same boundary the Markov / P-nucleotide /
  allele-usage / trim slices respected for their legacy
  fields.
- **In-frame / out-of-frame stratification.** v1 ignores
  `productive` / `vj_in_frame` / `stop_codon` columns;
  the estimator consumes every record uniformly.
- **Pandas DataFrame as first-class input.** v1 accepts
  list-of-dict / path / file handle. Users convert
  pandas via `df.to_dict(orient="records")`.
- **Per-segment / per-allele NP-length conditioning.**
  The typed plane is per-key (NP1, NP2) — flat by design.
- **Closing the typed-plane validator's `NP2`-on-VJ gap.**
  Documented in §2.2 as a pre-existing asymmetry vs
  `trims`. A separate slice could tighten the validator
  in lockstep with similar tightening elsewhere.

---

## 12. Summary table

| Concern | Today | Recommendation |
|---|---|---|
| AIRR fields populated by GenAIRR | Direct integer `np1_length` and `np2_length` populated reliably; VJ chains hard-zero `np2_length` (no NP2 region) | **Estimator consumes the direct integer fields**; sequence fields (`np1`, `np2`) ignored |
| Typed plane | `ReferenceEmpiricalModels.np_lengths: Dict[str, EmpiricalDistributionSpec]` with `NP_KEYS = ("NP1","NP2")` + key-membership validator (no chain-aware rejection for `NP2` on VJ at attach time) | **Reuse verbatim.** No new spec class. |
| Bridge resolver | `_np_lengths_from_models` + `extract_np_lengths` legacy fallback; precedence kwarg > typed plane > legacy > uniform fallback already in place | **Reuse verbatim.** No bridge changes. |
| Engine surface | `GenerateNPPass` consumes the length distribution; `parameter_signature` folds via `fmt_int_dist` | **Reuse verbatim.** No engine changes. |
| Provenance distinction | NP lengths vs P-nucleotide lengths vs junction arithmetic vs `len(np1)` — all distinguished at AIRR layer | **Estimator consumes ONLY `np1_length` / `np2_length`.** Documented + pinned. |
| Chain-type enforcement | `metadata.has_d` is the authoritative classifier; `Experiment.recombine(np2_lengths=...)` raises on VJ at compile time; spec validator does NOT reject `NP2` on VJ at attach time | VJ skips `np2_length` contribution + warns once; VDJ consumes both fields |
| Builder method signature | n/a | `estimate_np_length_distributions(rearrangements, *, min_count=1, pseudocount=0.0, replace=True) → self` |
| Stage entry shape | n/a | Canonical `{stage, inputs, inferred, warnings}` per trim-audit §6.2 |
| Manifest block | `np_length_keys: []` + `legacy_np_lengths_present: True` exist today; no structured `np_length_models` sub-block | New `models.np_length_models` block per §8.2 (additive) |
| Plan signature folding | NP-length distributions already fold via `fmt_int_dist` | **No soft gap inherited** — same property the trim slice has |
| Pre-flight stop-and-report | **NOT triggered** — bundled cartridges populate `np1_length` / `np2_length` at significant rates; VJ chains correctly hard-zero `np2_length` | Proceed to implementation slice |

The slice is **the lightest one yet**: simpler than trim
because the field shape has no four-key fanout; simpler
than allele-usage because there is no tie-set policy and
no new spec class. The estimator should land in ~70 lines
of Python.
