# Allele Usage Estimation — Audit + v1 Implementation Shipped

**Status: v1 implementation shipped.** The audit's recommended
slice landed end-to-end: `AlleleUsageSpec` typed dataclass +
`ReferenceEmpiricalModels.allele_usage` plane + bridge
resolver in [`_dataconfig_extract.py`](../src/GenAIRR/_dataconfig_extract.py)
+ `Experiment.recombine` precedence shim (kwarg > cartridge
plane > uniform) + `ReferenceCartridgeBuilder.estimate_allele_usage`
builder method with the three `ambiguous` policies + manifest
`models.allele_usage` block. The first data-derived estimator
on the `ReferenceCartridgeBuilder` facade per the
cartridge-authoring audit's §11.2 recommended ordering.

Companion artefacts:

- [`tests/test_allele_usage_estimation_contract.py`](../tests/test_allele_usage_estimation_contract.py)
  — 25 pins, flipped to post-slice present-state (scaffold,
  present, soft-gap-pin, and the five absence pins now
  flipped to present).
- [`tests/test_allele_usage_estimation_implementation.py`](../tests/test_allele_usage_estimation_implementation.py)
  — 15 behaviour tests covering: validation, lowering
  precedence, `gene_use_dict` orphan preservation,
  fractional / truth-first / reject policies, unknown-allele
  rejection, VJ-D warning, VDJ-missing-D rejection,
  `min_count` filter, stage entry shape, manifest, pickle
  round-trip, and contract pin flips.
- [`tests/test_allele_usage_estimation_release.py`](../tests/test_allele_usage_estimation_release.py)
  — release-tier smoke: tiny inline-FASTA VDJ cartridge,
  estimator on a synthetic AIRR record set, recombination
  exhibits the estimated bias, manifest + report carry the
  estimator's output, explicit `recombine` kwarg overrides
  the typed plane.

Two boundaries explicitly stayed deferred per audit §11:

- **Plan-signature soft gap 1 tightening.** The estimator's
  output lowers through the same `v_allele_weights` surface
  that has the documented soft gap. Closing the gap (per
  plan-signature completeness audit §9) is a separate slice
  covering BOTH the per-experiment kwarg AND the
  cartridge-default path in lockstep.
- **`gene_use_dict` auto-lift.** The legacy field stays
  orphan; the estimator does NOT populate it. A separate
  cartridge-migration slice could add explicit opt-in
  auto-lift from `gene_use_dict` → `allele_usage`.

The audit body below is preserved for traceability; sections
marked **[Shipped]** describe how the recommendations actually
landed.

---

## Original audit

Designs the first `ReferenceCartridgeBuilder.estimate_*`
method — `estimate_allele_usage` — from observed AIRR
rearrangement data. Per the cartridge-authoring audit's
§11.2 recommended ordering, this is the simplest and most
useful data-derived estimator to land first.

Companion to
[`tests/test_allele_usage_estimation_contract.py`](../tests/test_allele_usage_estimation_contract.py)
which freezes (a) the engine's existing weighted-allele
selection surface that the estimator's output lowers into,
(b) the cartridge-state gaps the estimator does NOT
silently fill (`gene_use_dict` orphan boundary,
manifest absence, plan-signature soft gap), and (c) the
builder/method absence pre-slice.

**Pre-flight finding (§7 below): clean-yes — no
stop-and-report condition.** The engine already implements
end-to-end weighted allele sampling at the
`Experiment.recombine(v_allele_weights=...)` / refdata-bridge
/ `plan.push_sample_allele(weights=...)` /
`AllelePoolDist::from_weights` chain. The legacy
`DataConfig.gene_use_dict` field is **genuinely orphan** —
listed in `_DOCUMENTED_ORPHAN_DATACONFIG_FIELDS`, never read
by the live pipeline, surfaced only by the MCP read-only
`gene_use` diagnostic endpoint and the dead-code
`DataConfig.validate()` method (which is itself not invoked
anywhere in the live simulation path). The estimator landed
without disturbing any existing surface; the storage decision
took Option B per §3 below.

---

## 1. Q1 — Input shape

### 1.1 Accepted input forms

The estimator's `rearrangements` parameter accepts three
shapes, mirroring the convention the existing MCP summary
helper [`_mcp_summary.py:73-75`](../src/GenAIRR/_mcp_summary.py#L73-L75)
already uses:

| Shape | Source |
|---|---|
| `list[dict]` (each dict is one rearrangement record) | Native Python pipeline output; MCP server input. |
| `pathlib.Path` / `str` (filesystem path) | AIRR TSV file (tab-separated, one row per rearrangement). |
| open text file handle | In-memory or pre-opened stream. |

Path / file inputs are parsed via `csv.DictReader` with
`delimiter='\t'` — AIRR-C TSV is the canonical exchange
format. **DataFrame inputs are NOT first-class** in v1; users
who want to pass a pandas DataFrame call
`df.to_dict(orient="records")` themselves. This matches the
existing repo discipline of treating pandas as an optional
extra.

### 1.2 Required columns

Per AIRR-C rearrangement schema, the only required columns
for allele-usage estimation are:

```text
v_call
d_call    (only for VDJ cartridges; absent / empty OK on VJ)
j_call
```

Any other columns (`junction`, `sequence`, …) are ignored.
Rows missing the required column(s) for the cartridge's
chain type land in `report.rejected` with reason
`"missing_required_column"`.

### 1.3 Comma-separated tie-set handling

`v_call` etc. may carry tie-sets like `"IGHV1*01,IGHV2*01"`
when an aligner couldn't disambiguate. The estimator's
`ambiguous` kwarg picks the policy (§3 below). The existing
`first_call` helper in `_mcp_summary.py` (taking
`value.split(",", 1)[0].strip()`) is the precedent for the
"truth-first" policy.

### 1.4 Pinned

- `pin_present_mcp_summary_first_call_helper_treats_v_call_as_comma_separated`
- `pin_scaffold_csv_dictreader_handles_airr_tsv_via_stdlib_only`

---

## 2. Q2 — Weight semantics + storage surface

### 2.1 Existing engine surface (confirmed wired end-to-end)

Empirical investigation confirms the engine already accepts
weighted allele sampling **as a per-experiment kwarg** through
this chain:

| Layer | Site |
|---|---|
| User-facing DSL | [`Experiment.recombine(v_allele_weights={...}, d_allele_weights={...}, j_allele_weights={...})`](../src/GenAIRR/experiment.py) |
| Validation | [`Experiment._resolve_allele_weights`](../src/GenAIRR/experiment.py) — rejects unknown allele names + non-positive weights + D weights on VJ chains; converts dict → dense pool-aligned `Tuple[float, ...]` with `1.0` default for unlisted alleles |
| Pipeline IR | [`_RecombineStep.weights_{v,d,j}: Optional[Tuple[float, ...]]`](../src/GenAIRR/_pipeline_ir.py) |
| Lowering | [`_compile.py::_lower_recombine`](../src/GenAIRR/_compile.py) → `plan.push_sample_allele("V", refdata, weights=v_weights)` |
| Bridge | [`plan.rs::push_sample_allele`](../engine_rs/src/python/plan.rs) — validates `weights.len() == pool.len()` and each `weight > 0.0 && finite` |
| Engine | [`AllelePoolDist::from_weights(pool, w)`](../engine_rs/src/dist/allele_pool.rs) — categorical distribution wrapped by `SampleAllelePass` |

The estimator's per-segment `{allele_name: weight}` output
maps cleanly onto this surface. The estimator does NOT need
to add a new engine surface.

### 2.2 Storage decision — two clean options

The estimator's output (per-segment weights) needs to live
somewhere on the cartridge so a downstream `Experiment.on(cfg)
.recombine()` consumes it without explicit kwargs. Two
candidates:

**Option A — Per-experiment kwarg only (status quo + builder
report).** The estimator writes results into
`report.stages["estimate_allele_usage"]["inferred"]`. Users
who want the weights applied call
`Experiment.on(cfg).recombine(v_allele_weights=cfg.build_report.estimated_v_weights, ...)`
explicitly. **Pros:** zero new surface; preserves the
documented "soft gap 1" from the plan-signature
completeness audit (which the user explicitly accepted).
**Cons:** the cartridge-attached audit trail doesn't drive
sampling automatically; a user who passes the cartridge to
`Experiment.on(cfg).recombine()` without re-specifying
`v_allele_weights` gets uniform sampling silently.

**Option B — New typed cartridge plane
`ReferenceEmpiricalModels.allele_usage` (recommended).**
The estimator writes results into a new typed plane on the
existing `ReferenceEmpiricalModels` shape, parallel to
`np_lengths` / `trims` / `np_bases` / `p_nucleotide_lengths`.
The `_dataconfig_extract.py` resolver picks up the plane and
threads it through `_RecombineStep.weights_{v,d,j}` so
`Experiment.on(cfg).recombine()` uses the estimated weights
by default. **Pros:** cartridge-driven default; consistent
with the typed-plane discipline the recent biology slices
established; cartridges with builder-estimated weights
remain replayable without explicit kwargs. **Cons:** requires
a new typed `AlleleUsageSpec` validator + extract-recombine-
defaults extension + the plan-signature soft gap surfaces in
the manifest summary block.

**The audit recommends Option B.** Reasons:

1. It matches the documented pattern for every other
   cartridge-driven default (NP lengths / NP bases /
   P-nucleotide lengths / trims).
2. The cartridge becomes self-contained: a user who runs
   `Experiment.on(cfg).recombine()` without explicit
   weights gets the estimated biology, just as they get
   the cartridge's typed NP-length distribution today.
3. The plan-signature soft gap can be addressed separately
   (the plan-signature completeness audit's soft gap 1
   tightening, deferred). For v1 the gap is documented and
   the strict-sampler rejection backstops.

### 2.3 `gene_use_dict` — keep orphan; do NOT auto-lift

The legacy `DataConfig.gene_use_dict` field is **genuinely
orphan today** (confirmed §7 stop-and-report check). The
estimator MUST NOT populate it. Reasons:

- Its dict-of-dict-of-int shape (`{"V": {"IGHV1-2*02":
  0.05, ...}, "D": ..., "J": ...}`) predates the typed
  empirical-models layer and would require auto-lift logic to
  reach the engine.
- Auto-lifting would silently change output bytes for the
  106 bundled cartridges (every one carries a populated
  `gene_use_dict` today).
- Same boundary the Markov / P-nucleotide slices respected
  for `NP_transitions` / `NP_first_bases` /
  `p_nucleotide_length_probs`.

A separate cartridge-migration slice may later add explicit
opt-in auto-lift from `gene_use_dict` → `allele_usage`. The
estimator slice does NOT touch this boundary.

### 2.4 Pinned

- `pin_scaffold_experiment_recombine_accepts_allele_weights_kwargs`
- `pin_scaffold_plan_push_sample_allele_accepts_weights_arg`
- `pin_scaffold_pipeline_ir_recombine_step_has_weights_fields`
- `pin_absence_no_allele_usage_field_on_reference_empirical_models`
- `pin_absence_no_allele_usage_spec_dataclass`

---

## 3. Q3 — Ambiguous (tie-set) calls

### 3.1 Three options analysed

| Policy | Behaviour | Risk |
|---|---|---|
| **fractional** (recommended default) | A row with `v_call="IGHV1*01,IGHV2*01"` contributes `0.5` to each. | Smooths over aligner ambiguity in a principled way; biased toward over-representing alleles that share tie-sets often. |
| **truth-first** | Take only the first comma-separated entry — matches existing `_mcp_summary.first_call` convention. | Loses information but is deterministic; matches what an aligner-consumer would report. |
| **reject** | Drop ambiguous rows entirely; surface count in report. | Most conservative; risks dropping a large fraction of rows when ambiguity is common (which is typical for IGH V allele calls). |

### 3.2 Recommendation

Default `ambiguous="fractional"` — the audit's clean choice
for the smoothing-vs-bias trade-off. Add `ambiguous="truth_first"`
and `ambiguous="reject"` as explicit alternatives so users
who want deterministic policy can choose. **`ambiguous="reject"`
is the recommended strict-mode counterpart for users who want
zero estimator noise.**

### 3.3 Pinned

- `pin_absence_no_ambiguous_policy_enum_today`
- `pin_absence_no_estimate_allele_usage_ambiguous_kwarg_today`

---

## 4. Q4 — Missing / unknown alleles

### 4.1 Decision matrix

| Case | Default behaviour | Report destination |
|---|---|---|
| Row missing required column (`v_call` / `j_call` always; `d_call` on VDJ) | Skip row | `report.rejected` with `reason="missing_required_column"` |
| Row's allele name not found in cartridge's V/D/J pool | Skip row | `report.rejected` with `reason="unknown_allele"` (per-allele) |
| Row's `d_call` empty on VDJ cartridge (no D called by aligner) | Skip row | `report.rejected` with `reason="missing_d_call_on_vdj"` |
| Row's `d_call` non-empty on VJ cartridge (D contamination) | Skip the D contribution, keep V/J | `report.warnings` (cartridge-level), not a per-row rejection |
| Below `min_count` threshold for an allele's final count | Drop from result, surface count | `report.warnings` (cartridge-level): `"3 alleles below min_count=10 dropped"` |

### 4.2 Recommended v1 — skip unknowns + report counts

Per the user brief: "skip unknowns and report counts by
segment". The estimator emits a single `inferred` block with:

```python
{
    "V": {"IGHV1-2*02": 0.0521, ...},
    "D": {...},
    "J": {...},
    "skipped": {
        "missing_required_column": 12,
        "unknown_allele": {"V": 3, "D": 1, "J": 0},
        "missing_d_call_on_vdj": 17,
    },
    "below_min_count": {"V": 2, "D": 0, "J": 0},
}
```

### 4.3 Pinned

- `pin_absence_no_estimate_allele_usage_min_count_kwarg_today`
- `pin_scaffold_report_rejected_entries_carry_stage_segment_reason`

---

## 5. Q5 — Chain behavior

### 5.1 Chain-type-driven enforcement

| Chain type | Required columns | D handling |
|---|---|---|
| VJ (`BCR_LIGHT_KAPPA`, `BCR_LIGHT_LAMBDA`, `TCR_ALPHA`, `TCR_GAMMA`) | `v_call`, `j_call` | `d_call` absent: OK. `d_call` present-and-non-empty: warn and skip the D contribution. |
| VDJ (`BCR_HEAVY`, `TCR_BETA`, `TCR_DELTA`) | `v_call`, `d_call`, `j_call` | Missing `d_call`: skip row (`missing_d_call_on_vdj`). |

The cartridge's `metadata.has_d` (from `ConfigInfo`) is the
authoritative source for VJ-vs-VDJ classification. Same
boundary the existing `Experiment.recombine` enforces for
`np2_lengths` (VJ rejects it loudly) and the existing
`_resolve_allele_weights` enforces for `d_allele_weights`
(VJ rejects it).

### 5.2 Pinned

- `pin_scaffold_config_info_has_d_drives_chain_classification`
- `pin_scaffold_recombine_rejects_np2_on_vj_chain`

---

## 6. Q6 — Builder integration

### 6.1 Method signature

```python
def estimate_allele_usage(
    self,
    rearrangements: Union[List[Dict[str, Any]], Path, IO[str]],
    *,
    min_count: int = 1,
    ambiguous: Literal["fractional", "truth_first", "reject"] = "fractional",
) -> "ReferenceCartridgeBuilder":
    ...
```

Returns `self` for fluent chaining (matches the existing v1
stage methods). Idempotent — calling twice re-estimates and
overwrites the previous output with `replaced=True` on the
stage entry (same discipline `infer_v_subregions` uses).

### 6.2 Stage entry shape

Per the user brief, the appended stage carries the canonical
`{stage, inputs, inferred, warnings}` shape pinned by the
release smoke test:

```python
{
    "stage": "estimate_allele_usage",
    "inputs": {
        "record_count": 18412,
        "ambiguous": "fractional",
        "min_count": 1,
        "source": "<path or 'records'>",
    },
    "inferred": {
        "V": {<allele_name>: <normalised_weight>, ...},
        "D": {<allele_name>: <normalised_weight>, ...},
        "J": {<allele_name>: <normalised_weight>, ...},
        "skipped": {
            "missing_required_column": 12,
            "unknown_allele": {"V": 3, "D": 1, "J": 0},
            "missing_d_call_on_vdj": 17,
        },
        "below_min_count": {"V": 2, "D": 0, "J": 0},
    },
    "warnings": [
        "3 alleles below min_count=10 dropped",
        ...
    ],
}
```

Weights are **normalised so each segment's values sum to
1.0** (matching the existing `_mcp_summary._top_counter`
discipline for "rate of use"). The bridge's
`_resolve_allele_weights` accepts any positive finite
weight (renormalises internally), so the normalisation is
purely about presentation in the build report.

### 6.3 Lowering path (when Option B lands)

Per audit §2.2 recommendation:

1. New typed `AlleleUsageSpec` dataclass in
   `reference_models.py` carrying
   `{"V": Dict[str, float], "D": Dict[str, float], "J":
   Dict[str, float]}` with chain-type-aware validation.
2. New `ReferenceEmpiricalModels.allele_usage: Optional[AlleleUsageSpec] = None`
   field.
3. New `_dataconfig_extract._allele_usage_from_models` resolver
   that returns a dense pool-aligned weight vector per segment
   (using the same `_resolve_allele_weights` shape the
   per-experiment kwarg already produces).
4. `_RecombineStep.weights_{v,d,j}` populated from the
   cartridge default when not overridden by the user's kwarg
   — same precedence discipline `np1_base_pairs` etc. follow.
5. Manifest gains a `models.allele_usage` block:
   ```python
   "allele_usage": {
       "V_count": 285,
       "D_count": 27,
       "J_count": 13,
       "legacy_gene_use_dict_present": True,
       "legacy_fallback": False,
       "in_plan_signature": False,  # documented soft gap
       "in_content_hash": False,
   }
   ```

### 6.4 Pinned

- `pin_absence_no_estimate_allele_usage_method_today`
- `pin_scaffold_builder_stage_entry_shape_includes_inputs_inferred_warnings`
- `pin_scaffold_idempotency_pattern_via_replaced_flag`

---

## 7. Q7 — Lowering to engine + stop-and-report determination

### 7.1 Engine surface inventory (confirmed live)

The engine **already implements** weighted allele sampling
end-to-end through the chain pinned in §2.1. The estimator
slice does NOT need to add a new engine surface; it lowers
into the existing one. This is the cleanest possible
implementation path.

### 7.2 `gene_use_dict` stop-and-report check

The user brief: "If the audit finds `gene_use_dict` is
already silently used but not represented in manifest/
signature, stop and report."

**Verdict: NOT MET — proceed.** The grep-based investigation
confirms `gene_use_dict` is referenced in exactly five live
sites:

1. [`dataconfig/data_config.py:44`](../src/GenAIRR/dataconfig/data_config.py#L44) —
   listed in `_DOCUMENTED_ORPHAN_DATACONFIG_FIELDS`.
2. [`dataconfig/data_config.py:124`](../src/GenAIRR/dataconfig/data_config.py#L124) —
   dataclass field declaration (`gene_use_dict: Dict[str, Any] = field(default_factory=dict)`).
3. [`dataconfig/data_config.py:329-334`](../src/GenAIRR/dataconfig/data_config.py#L329-L334) —
   `DataConfig.validate()` requires the field be non-empty
   AND contain `"V"` / `"J"` keys. **But `DataConfig.validate()`
   is itself dead code** — confirmed not called anywhere in
   `src/GenAIRR/` (the only matching `.validate()` site in the
   pipeline is `rules.validate()` on a `ReferenceRulesSpec`).
4. [`utilities/mcp_helpers.py:641`](../src/GenAIRR/utilities/mcp_helpers.py#L641) —
   MCP read-only `gene_use` diagnostic endpoint:
   `getattr(dc, "gene_use_dict", {}) or {}`. Same shape as
   the `p_nucleotides` / `np_params` endpoints — surfaces the
   dict for inspection, does NOT feed the simulator.

The field is genuinely orphan in the simulation pipeline.
**No silent corruption gap exists.** The estimator slice is
free to:

- Land without touching `gene_use_dict`.
- Recommend Option B storage on a new
  `ReferenceEmpiricalModels.allele_usage` plane.
- Preserve the orphan boundary (no auto-lift).

### 7.3 Pinned

- `pin_present_gene_use_dict_field_exists_and_is_in_orphan_list`
- `pin_present_gene_use_dict_has_no_simulator_pipeline_consumer`
- `pin_scaffold_dataconfig_validate_is_dead_code_today`
- `pin_present_mcp_helpers_gene_use_endpoint_is_read_only`

---

## 8. Manifest exposure + plan-signature behaviour (deferred boundaries)

### 8.1 Manifest

Today's `cartridge_manifest()["models"]` block exposes:

- `np_length_keys` / `trim_keys` (typed plane keys).
- `np_base_models` / `p_nucleotide_models` (typed plane
  blocks with `legacy_*_present` flags + `legacy_fallback`
  + `supported_*` + `in_plan_signature` / `in_content_hash`).
- `shm` (mutation model inventory).

There is **NO existing `allele_usage` or `gene_use` block**
in the manifest. The estimator slice adds the
`allele_usage` block per §6.3.

### 8.2 Plan-signature behaviour — documented soft gap

The plan-signature completeness audit's **soft gap 1**
documents that `Experiment.recombine(v_allele_weights=...)`
does NOT change the plan signature today. Lowering the
estimator's output through the same surface inherits the
soft gap. Per the audit boundary:

- The strict-mode sampler backstops cross-cartridge replay
  when the recorded allele ID is outside the narrowed
  support.
- `refdata_content_hash` covers cartridge pool identity at a
  coarser layer.

The estimator slice does NOT close the soft gap.
A separate tightening slice (audit §9 recommends `fold
SampleAllele distribution support when narrowed`) would
close it for both the per-experiment kwarg AND the
cartridge-driven default in lockstep.

### 8.3 Pinned

- `pin_scaffold_manifest_does_not_currently_expose_allele_usage_block`
- `pin_present_plan_signature_soft_gap_for_allele_weights_holds`

---

## 9. Implementation order — shipped

A single self-contained slice landed the estimator, in six
sub-steps. **[Shipped]** annotations confirm each landed as
designed:

1. **Typed spec** **[Shipped]** — `AlleleUsageSpec`
   dataclass in [`src/GenAIRR/reference_models.py`](../src/GenAIRR/reference_models.py)
   carries `v` / `d` / `j` dict fields (default empty) with
   chain-type-aware `validate(chain_type=...)`. Rejects
   non-empty D entries on VJ chains, non-positive weights,
   non-finite weights, and empty allele names. Per-cartridge
   name validation defers to the bridge layer
   (`_resolve_allele_weights`) for performance.

2. **`ReferenceEmpiricalModels.allele_usage`** **[Shipped]** —
   `Optional[AlleleUsageSpec] = None` field added to the
   existing typed-plane container, validated via the
   `ReferenceEmpiricalModels.validate()` chain.

3. **Resolver + lowering** **[Shipped]** —
   [`_dataconfig_extract.py`](../src/GenAIRR/_dataconfig_extract.py)
   gained `_allele_usage_from_models(models, segment)`
   returning the per-segment `{name: weight}` dict.
   [`Experiment.recombine`](../src/GenAIRR/experiment.py)
   applies the precedence shim before `_resolve_allele_weights`:
   kwarg > cartridge plane > uniform. `_explicit_models`
   empty-check widened to include `allele_usage`.

4. **Builder method** **[Shipped]** —
   [`ReferenceCartridgeBuilder.estimate_allele_usage(
   rearrangements, *, min_count=1.0, ambiguous="fractional",
   replace=True)`](../src/GenAIRR/cartridge_builder.py)
   reads input via `csv.DictReader` (path / file handle) or
   directly (list of dicts). Builds per-segment Counter;
   resolves tie-sets per policy; rejects unknown alleles
   into `report.rejected`. Stage entry shape per §6.2.
   Constructs an `AlleleUsageSpec` and either attaches it
   to the existing `_reference_models` or constructs a fresh
   one. Idempotency via `previously_estimated` + `replace`
   kwarg (mirrors `infer_v_subregions` discipline).

5. **Manifest extension** **[Shipped]** —
   `cartridge_manifest()` gained the `models.allele_usage`
   block per §6.3. The shipped shape carries 7 keys:
   `available` (bool) / `segments` (canonical `["V","D","J"]`)
   / `nonempty_segments` (list reflecting estimator output) /
   `legacy_gene_use_dict_present` /
   `legacy_fallback=False` / `in_plan_signature=False`
   (documented soft gap 1) / `source` (provenance label).

6. **Tests** **[Shipped]** —
   [`tests/test_allele_usage_estimation_implementation.py`](../tests/test_allele_usage_estimation_implementation.py)
   ships 15 behaviour tests (fractional / truth-first /
   reject policies; missing D handling on VJ vs VDJ;
   `min_count` threshold; replay-through-built-cartridge
   end-to-end; manifest exposure; pickle round-trip;
   `gene_use_dict` orphan preservation).
   [`tests/test_allele_usage_estimation_release.py`](../tests/test_allele_usage_estimation_release.py)
   adds release-tier smoke covering the full slice
   composition.

Actual cost — close to the audit estimate:

- ~120 lines Python (spec + resolver + builder method +
  manifest block) — matched estimate
- ~340 lines tests (implementation + release smoke) —
  larger than estimate because the release smoke landed as
  a separate file per consolidation discipline
- ~30 lines docstrings + audit-doc-flip — matched estimate

---

## 10. Test surface — what this audit pins

Mirrored in
[`tests/test_allele_usage_estimation_contract.py`](../tests/test_allele_usage_estimation_contract.py).

### `pin_scaffold_*` — live surfaces the estimator builds on

1. `Experiment.recombine` accepts `v_allele_weights` /
   `d_allele_weights` / `j_allele_weights` kwargs.
2. `_resolve_allele_weights` converts `{name: weight}` dict
   to dense pool-aligned `Tuple[float, ...]`.
3. `_RecombineStep.weights_{v,d,j}` field shape.
4. `_compile.py::_lower_recombine` passes weights to
   `plan.push_sample_allele(weights=...)`.
5. `plan.push_sample_allele` accepts `weights: Option<Vec<f64>>`.
6. `AllelePoolDist::from_weights` is the engine consumer.
7. `_mcp_summary.first_call` parses comma-separated v_call
   by taking the first entry — establishes the convention
   the `ambiguous="truth_first"` policy mirrors.
8. `_mcp_summary` reads `v_call` / `d_call` / `j_call` from
   record dicts — confirms the AIRR column convention.
9. `csv.DictReader` (stdlib only) handles AIRR TSV inputs.
10. `ConfigInfo.has_d` is the authoritative chain-type
    classifier the estimator dispatches on.
11. `ReferenceEmpiricalModels` carries `np_lengths` / `trims`
    / `np_bases` / `p_nucleotide_lengths` typed planes —
    `allele_usage` would slot in as the fifth.
12. `cartridge_manifest()["models"]` block carries
    `np_length_keys` / `trim_keys` / `np_base_models` /
    `p_nucleotide_models` but NO `allele_usage` block.
13. Builder stage entry shape `{stage, inputs, inferred,
    warnings}` (pinned by the release smoke test).
14. Builder idempotency pattern (`replaced=True` flag on
    repeated stage calls — pinned by `infer_v_subregions`).

### `pin_present_*` — stop-and-report verification

15. `gene_use_dict` field exists on `DataConfig` and is in
    `_DOCUMENTED_ORPHAN_DATACONFIG_FIELDS`.
16. `gene_use_dict` has NO simulator-pipeline consumer
    (verified by exhaustive grep against
    `_dataconfig_extract.py` / `_compile.py` /
    `experiment.py` / `_refdata_resolver.py` /
    `engine_rs/src/`).
17. `DataConfig.validate()` exists but is dead code (not
    called anywhere in the live source).
18. MCP `gene_use` endpoint is read-only diagnostic
    (mirrors the `p_nucleotides` orphan-surface pattern).

### `pin_present_*` — documented boundary state

19. Plan-signature completeness audit's soft gap 1 holds:
    two experiments differing only by
    `v_allele_weights` produce equal plan signatures
    (the estimator's output inherits this boundary).

### `pin_absence_*` — gaps the implementation slice closes

20. `ReferenceCartridgeBuilder.estimate_allele_usage` is not
    a method today.
21. `ReferenceEmpiricalModels.allele_usage` field is absent.
22. `AlleleUsageSpec` dataclass is absent.
23. Manifest does not yet expose
    `models.allele_usage` block.
24. No `ambiguous` policy enum / kwarg names exist on the
    builder.

### Doc anchor

25. Audit doc exists and references the contract file;
    section structure intact.

---

## 11. Out of scope

Documented here so a future implementer doesn't accidentally
expand the work.

- **Tightening plan-signature soft gap 1.** The estimator's
  output flows into the same `v_allele_weights` surface
  that has the documented soft gap. Closing the gap (per
  plan-signature completeness audit §9) is a separate
  slice covering BOTH the kwarg AND the cartridge-default
  path in lockstep.
- **Auto-lift of legacy `gene_use_dict`.** Stays
  deferred — same boundary the Markov / P-nucleotide
  slices respected for their legacy fields.
- **`estimate_trim_distributions`** + four other deferred
  estimators per the cartridge-authoring audit §11.2
  ordering. The next-after-this slice; same builder
  facade, same stage shape.
- **Pandas DataFrame as first-class input.** v1 accepts
  list-of-dict / path / file handle; DataFrame users call
  `.to_dict(orient="records")` themselves.
- **Per-gene rollup (`gene_use` separate from
  `allele_use`).** v1 estimates per-allele weights only;
  the bundled cartridges' legacy `gene_use_dict` is
  per-gene + per-allele nested — a translation slice
  would be additive.
- **D-J pairing constraint estimation** (the legacy
  `dj_pairing_map` field). Separate orphan; separate
  audit.

---

## 12. Summary table

| Concern | Today | Recommendation |
|---|---|---|
| Engine weighted-allele surface | Already wired end-to-end via `Experiment.recombine(v_allele_weights=...)` → `_RecombineStep.weights_*` → `plan.push_sample_allele(weights=...)` → `AllelePoolDist::from_weights` | **Reuse verbatim.** No new engine surface. |
| Storage for estimator output | Per-experiment kwarg only; no cartridge-driven default | **Add `ReferenceEmpiricalModels.allele_usage`** (Option B in §2.2); resolver threads into `_RecombineStep.weights_*` with user kwarg taking precedence |
| `gene_use_dict` orphan boundary | Genuinely orphan (no simulator-pipeline consumer); MCP diagnostic + dead-code `validate()` are the only readers | **Keep orphan.** Estimator does NOT auto-lift. |
| AIRR input shape | `_mcp_summary` already parses `v_call` / `d_call` / `j_call` as comma-separated tie-sets | `list[dict]` / path / file handle accepted; `csv.DictReader` for path input |
| Ambiguous tie-set policy | `_mcp_summary` uses truth-first | Builder default `ambiguous="fractional"`; `"truth_first"` + `"reject"` alternatives |
| Missing / unknown alleles | n/a | Skip + per-row `report.rejected` entry with structured reason |
| Chain type enforcement | `metadata.has_d` is the authoritative classifier (existing pattern) | VJ skips D contribution + warns; VDJ requires `d_call` |
| Builder method signature | n/a | `estimate_allele_usage(rearrangements, *, min_count=1, ambiguous="fractional") → self` |
| Stage entry shape | n/a | Canonical `{stage, inputs, inferred, warnings}` per §6.2 |
| Manifest block | No `allele_usage` block today | New `models.allele_usage` block per §6.3 |
| Plan signature folding | `v_allele_weights` is **soft gap 1** in the plan-signature completeness audit | Estimator inherits the gap; tightening is a separate slice |
| Pre-flight stop-and-report | **NOT triggered** — `gene_use_dict` is genuinely orphan; no silent corruption | Proceed to implementation slice |

The estimator is **architecturally tractable**: the engine
surface is already wired, the AIRR input convention is
already established by `_mcp_summary`, the cartridge plane
discipline is well-defined by the recent NP / Markov / P
slices, and the legacy `gene_use_dict` boundary is clean.
The recommended next step is the single estimator
implementation slice per §9.
