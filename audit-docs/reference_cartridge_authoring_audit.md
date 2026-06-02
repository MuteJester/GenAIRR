# Reference Cartridge Authoring / Inference API — Audit + v1 Implementation Shipped

**Status: v1 implementation shipped + first estimator shipped.**
The audit's recommended facade landed:
[`src/GenAIRR/cartridge_builder.py`](../src/GenAIRR/cartridge_builder.py)
exports `ReferenceCartridgeBuilder` + `CartridgeBuildReport`,
re-exported at the top level as `from GenAIRR import
ReferenceCartridgeBuilder, CartridgeBuildReport`. v1 covered
structural authoring (FASTA + identity + V-subregion
derivation + rules/models attachment) plus the build-report
audit trail.

**First estimator landed (post-v1):** `estimate_allele_usage`
shipped per [`docs/allele_usage_estimation_design.md`](allele_usage_estimation_design.md).
The estimator writes per-segment weights into the typed
`ReferenceEmpiricalModels.allele_usage` plane; recombination
samples from it by default unless overridden by an explicit
`Experiment.recombine(v_allele_weights=...)` kwarg. The
remaining five estimators in §11.2 stay deferred.

The three dead-reference sites are **cleaned up** in lockstep:
the `DataConfig.build_report` field docstring now references
the new builder; the private build script raises an explicit
`NotImplementedError` at module load with a porting hint
pointing at the audit doc; the build-cache mirror remains as a
compile artefact (out of the live import path).

Companion artefacts:

- [`tests/test_reference_cartridge_authoring_contract.py`](../tests/test_reference_cartridge_authoring_contract.py)
  — 25 pins, flipped to post-implementation present-state.
- [`tests/test_reference_cartridge_authoring_implementation.py`](../tests/test_reference_cartridge_authoring_implementation.py)
  — 15 behaviour tests covering the user-brief 12-item
  finish-condition plus three error-path bonus tests
  (`d_fasta` on VJ chain, missing `d_fasta` on VDJ chain,
  duplicate allele name).
- [`tests/test_reference_cartridge_authoring_release.py`](../tests/test_reference_cartridge_authoring_release.py)
  — release-tier composition: tiny inline-FASTA VDJ cartridge,
  build, manifest JSON-clean, build-report JSON-clean, compile
  path works under `allow_curatable_refdata()`, report
  contains the expected stages.

The audit body below is preserved for traceability; sections
marked **[Shipped]** describe how the recommendations actually
landed.

---

## Original audit

Surveys the cartridge-authoring landscape after the recent
biology slices (typed length / base / Markov NP, P-nucleotide,
V-subregion rates & counters), maps it against the realistic
input formats users arrive with (FASTA, AIRR TSV, observed
repertoire data), and designs a staged builder API + build-
report surface that makes cartridge creation reproducible and
auditable.

Companion to
[`tests/test_reference_cartridge_authoring_contract.py`](../tests/test_reference_cartridge_authoring_contract.py)
which pins (a) the absence of the proposed builder surface
today, (b) the live inference-adjacent helpers the new
builder will reuse, and (c) the dead/historical references
that need to be cleaned up in lockstep with the
implementation slice.

**Pre-flight finding (§3 below): clean-yes — no hidden live
builder.** A `RandomDataConfigBuilder` /
`CustomDataConfigBuilder` pair existed historically and was
removed from the live source tree. Stale references survive
in three places:

1. The `DataConfig.build_report` field's docstring still says
   "Populated by `RandomDataConfigBuilder.make_from_reference`",
   but `RandomDataConfigBuilder` is no longer importable
   anywhere in `src/GenAIRR/`.
2. `.private/scripts/build_imgt_configs.py` still imports
   `from GenAIRR.dataconfig.make.random import
   RandomDataConfigBuilder` and is **currently broken** —
   running it would `ModuleNotFoundError` at the top.
3. The build-cache mirror at
   `docs/build/lib.linux-x86_64-cpython-312/GenAIRR/dataconfig/make/`
   carries the historical implementation as a compile
   artifact, not as a live module.

This is NOT a "hidden partial builder still in the live
source" condition (the user-brief stop-and-report case).
Designing a new builder from scratch is the right move; the
historical shape can inform the staging but the
implementation slice should start fresh AND clean up the
three dead-reference sites in lockstep.

---

## 1. Current creation paths (live, verified)

| Path | Where | Public? | Audit notes |
|---|---|---|---|
| **Bundled `DataConfig` pickle load** | `GenAIRR.HUMAN_IGH_OGRDB` (and 105 others) via [`src/GenAIRR/data/__init__.py::_load_dataconfig`](../src/GenAIRR/data/__init__.py) | ✓ public (`GenAIRR.list_configs()`) | Lazy `__getattr__` resolves a name → pickle path → `pickle.load()` → `verify_integrity()`. Cartridges produced by an out-of-tree process; **no live in-tree path produces them.** |
| **Pickle-load arbitrary `.pkl` path** | n/a — no public helper; users `pickle.load()` themselves | ✗ (no public surface) | Plain Python pickle works because `DataConfig` is a normal dataclass. Round-trip survives `verify_integrity()` as long as the schema matches. |
| **Manual `DataConfig(...)` construction** | `DataConfig(name="…")` + field-by-field assignment | ✓ public (it's a dataclass) | Used by [`tests/test_cartridge_views.py`](../tests/test_cartridge_views.py) and a handful of integration tests. **No safety net** — a user can construct a malformed `DataConfig` that crashes at compile time. |
| **Manual `RefDataConfig` (Rust-backed) construction** | `ga.RefDataConfig.vj()` / `.vdj()` + `add_v_allele(...)` / `add_d_allele(...)` / `add_j_allele(...)` | ✓ public | Bridge layer; same surface used by [`tests/test_performance_budgets.py`](../tests/test_performance_budgets.py). Carries V/J/D pool only — does NOT carry `ReferenceEmpiricalModels` (NP, P, trim distributions). User must attach those via `cfg.reference_models = ReferenceEmpiricalModels(...)`. |
| **Typed empirical-models authoring** | `ReferenceEmpiricalModels(np_lengths=..., trims=..., np_bases=..., p_nucleotide_lengths=...)` | ✓ public | Validates at construction time via `EmpiricalDistributionSpec.validate` / `NpBaseModelSpec.validate`. Top-level + per-key + per-spec error paths exercised by `test_*_implementation.py`. |
| **Typed rules authoring** | `ReferenceRulesSpec(allowed_bases=..., v_anchor=..., j_anchor=...)` | ✓ public | Programmable interpretation layer; folded into `RefDataConfig.rules` by the bridge. |

---

## 2. Current inference-adjacent helpers (live)

These survive in the live source and the new builder should
reuse them rather than reinvent them.

| Helper | Where | Used by | New-builder reuse |
|---|---|---|---|
| `parse_fasta(file)` | [`utilities/misc.py:113`](../src/GenAIRR/utilities/misc.py#L113) | Top-level `from GenAIRR.utilities import parse_fasta`. Generic `>header\nbases…` parser. | **YES** — `from_fasta` constructor stage feeds this verbatim. |
| `compute_v_region_boundaries(v_allele)` | [`utilities/imgt_regions.py:46`](../src/GenAIRR/utilities/imgt_regions.py#L46) | `_refdata_resolver._resolve_v_subregions` (internal). Maps IMGT-gapped positions → ungapped intervals. | **YES** — `infer_v_subregions` stage feeds this for every V allele with `gapped_seq`. |
| `classify_position(...)` | [`utilities/imgt_regions.py:72`](../src/GenAIRR/utilities/imgt_regions.py#L72) | Inference-time IMGT-region classification for an assembled sequence position. | **YES** — `estimate_v_subregion_rates` stage uses this when AIRR rearrangement data is provided. |
| `_resolve_v_subregions(allele)` (private) | [`_refdata_resolver.py:194`](../src/GenAIRR/_refdata_resolver.py#L194) | Bridge-time V-subregion derivation for `dataconfig_to_refdata`. | Indirectly — the new builder writes subregions to `allele.subregions` so this resolver continues to work. |
| `DataConfig.cartridge_manifest()` | [`dataconfig/data_config.py::cartridge_manifest`](../src/GenAIRR/dataconfig/data_config.py) | Returns JSON-clean per-cartridge summary (identity, catalogue, rules, models, P-nuc, V-subregion, SHM, etc.). | **YES** — the build report's "manifest snapshot" entry uses this verbatim. |
| `DataConfig.verify_integrity()` / `compute_checksum()` | [`data_config.py:207, 261`](../src/GenAIRR/dataconfig/data_config.py) | Bundled-load integrity check. | **YES** — final build step calls `verify_integrity()` so a malformed `DataConfig` surfaces at build time, not at simulate time. |
| `EmpiricalDistributionSpec.validate` / `NpBaseModelSpec.validate` | [`reference_models.py`](../src/GenAIRR/reference_models.py) | Spec-layer validation. | **YES** — every `estimate_*` stage's output flows through these. |

---

## 3. Dead / historical references — **[Cleaned up in the v1 slice]**

### 3.1 Inventory (post-slice state)

| Reference site | Pre-slice state | Post-slice state |
|---|---|---|
| [`dataconfig/data_config.py:159-166`](../src/GenAIRR/dataconfig/data_config.py#L159-L166) — `DataConfig.build_report` field docstring | "Populated by `RandomDataConfigBuilder.make_from_reference`" (named a removed class) | **[Cleaned]** Docstring now references `GenAIRR.cartridge_builder.ReferenceCartridgeBuilder.build`. Pinned by `test_pin_present_build_report_docstring_now_references_new_builder`. |
| [`.private/scripts/build_imgt_configs.py`](../.private/scripts/build_imgt_configs.py) | `from GenAIRR.dataconfig.make.random import RandomDataConfigBuilder` — broken at module load with `ModuleNotFoundError` | **[Cleaned]** Top-level `raise NotImplementedError(...)` with explicit porting hint to `ReferenceCartridgeBuilder` + audit-doc reference; dead import moved into the unreachable function body for porting reference. Pinned by `test_pin_present_private_build_script_now_raises_explicit_legacy_error`. |
| `docs/build/lib.linux-x86_64-cpython-312/GenAIRR/dataconfig/make/...` | Historical implementation of `RandomDataConfigBuilder` / `CustomDataConfigBuilder` as a compile artefact | **Unchanged.** Build-cache mirror is not on the import path; regenerated by the next wheel build. Pinned absent by `test_pin_scaffold_historical_random_builder_module_is_gone` (verifies `ModuleNotFoundError` on import). |

### 3.2 Historical shape (for design inspiration only)

The build artefact shows the old builder pair had this
signature shape, useful as a starting point but not the
target design:

```python
# Historical (REMOVED from live source — do not resurrect verbatim):
builder = RandomDataConfigBuilder(
    convert_to_asc=False, *,
    species=Species.HUMAN, chain_type=ChainType.IGH,
    reference_set="IMGT",
)
config = builder.make(
    v_reference_path="...fasta",
    j_reference_path="...fasta",
    c_reference_path=None,
    d_reference_path="...fasta",
    v_anchor_finder=None, j_anchor_finder=None,
    keep_anchorless=True,
)
```

Limitations the new design should fix:

- **One giant `make(...)`** instead of staged steps — user
  can't override or audit individual inference stages.
- **No build report** — the historical class wrote inferred
  state directly to `self.dataconfig` with no per-stage
  provenance entry. The audit-trail purpose of
  `DataConfig.build_report` was advertised in the docstring
  but never delivered.
- **Pre-typed-models era** — predates
  `ReferenceEmpiricalModels`, `NpBaseModelSpec`,
  `p_nucleotide_lengths`, `ReferenceRulesSpec`. New builder
  emits these typed planes directly.
- **No AIRR-rearrangement empirical inference path** —
  `CustomDataConfigBuilder` did some of this but went straight
  to legacy `NP_lengths` / `trim_dicts` dict-of-dict planes,
  not the typed `ReferenceEmpiricalModels`. Replacing the
  legacy planes' output with typed-plane output is part of
  what makes the new builder genuinely audit-trail-friendly.

### Pinned

- `pin_present_build_report_docstring_references_dead_class`
- `pin_present_private_build_script_imports_dead_module`

---

## 4. Q1 — Creation paths inventory (combined)

The new builder is the **fifth** path:

| # | Path | Use case | Output shape |
|---|---|---|---|
| 1 | Bundled pickle load | Production users with bundled species | `DataConfig` from pickle, validated |
| 2 | Pickle file path | Sharing custom cartridges between teammates | `DataConfig` from arbitrary pickle |
| 3 | Manual `DataConfig(...)` | Tests, programmatic construction | Bare `DataConfig` (user-validated) |
| 4 | Manual `RefDataConfig.vj() / .vdj()` + `add_*_allele` | Tests + perf workloads | Rust-backed refdata only |
| 5 | **`ReferenceCartridgeBuilder` (proposed)** | New cartridge from FASTA / AIRR data | `DataConfig` with `build_report`, typed planes, manifest-clean |

Path 5 is the gap the audit closes.

---

## 5. Q2 — Realistic user inputs

The new builder targets these input shapes:

| Input format | Common provenance | Builder stage(s) |
|---|---|---|
| **V/D/J FASTA** (allele name in header, ungapped sequence body) | IMGT V-QUEST reference directory; OGRDB; hand-curated reference sets | `from_fasta(v_path=..., d_path=..., j_path=...)` constructor |
| **IMGT-gapped FASTA** (allele body carries `.` dots at IMGT positions) | IMGT V-QUEST reference directory; OGRDB | Same constructor; `gapped_seq` populated automatically; `infer_v_subregions` then works |
| **AIRR rearrangement TSV** (one row per observed sequence, AIRR-format columns) | Sequencing pipelines (MiXCR, Igblast, AIRR-compliant tools) | `with_rearrangement_data(path=...)` step, then `estimate_*` methods |
| **OGRDB/IMGT-style allele tables** (CSV with allele name / gene / functional status / anchor) | OGRDB downloads | `from_allele_table(path=...)` constructor variant or pre-`from_fasta` metadata override |
| **Per-sequence metadata** (species, locus, chain type) | User input; AIRR rearrangement column headers | `infer_identity(species=..., chain_type=..., reference_set=...)` step |
| **Observed repertoire data with counts** (sequence_id → count) | Repertoire-sequencing experiments | `with_observed_counts(...)` step input to `estimate_allele_usage` |

---

## 6. Q3 — What can be inferred safely

Per-field inference risk catalogue:

| Inference target | Input requirement | Risk | Builder stage |
|---|---|---|---|
| Identity: `species` / `chain_type` / `reference_set` / `name` / `source` | User-provided (cannot infer from FASTA alone) | None (requires user input) | `infer_identity(species=..., chain_type=..., reference_set=...)` |
| Allele fields: `name` / `gene` / `family` | FASTA header parsing | Low — depends on header convention (IMGT vs custom) | `from_fasta(...)` |
| Allele `ungapped_seq` | FASTA body | None | `from_fasta(...)` |
| Allele `gapped_seq` | FASTA body with IMGT dot-gap convention | Low | `from_fasta(...)` |
| Allele `functional_status` (`F`/`P`/`ORF`) | IMGT header tag or sidecar allele table | Low when present; **default to `None` when absent** rather than guessing | `from_fasta(...)` honours header; `from_allele_table(...)` honours column |
| Allele `anchor` (V Cys / J W/F position) | IMGT-gapped sequence (Cys at IMGT position 105×3 = 312–314) OR explicit override | Medium — gapped derivation only works for IMGT-numbered cartridges | `infer_anchors_from_imgt_gapped()` step; **fallback to None + manifest gap report** when derivation fails |
| V subregions (`FWR1` / `CDR1` / `FWR2` / `CDR2` / `FWR3`) | IMGT-gapped sequence | Low — `compute_v_region_boundaries` does the work | `infer_v_subregions()` step (uses existing helper) |
| Allele usage weights | Observed rearrangement counts | Low when data is well-sized; **emit count + warning when sample size is too small** | `estimate_allele_usage(from_rearrangements=...)` step |
| Trim distributions | Observed rearrangement counts (`v_trim_3` / `d_trim_*` / `j_trim_5` columns) | Medium — needs careful handling of in-frame vs out-of-frame reads and minimum sample size | `estimate_trim_distributions(from_rearrangements=..., min_samples=...)` step |
| NP length distributions | Observed rearrangement counts (`np1_length` / `np2_length` columns) | Medium | `estimate_np_length_distributions(...)` step |
| NP base empirical / Markov models | Observed rearrangement counts (`np1` / `np2` columns) | Medium-high — Markov inference needs care for short NP regions | `estimate_np_base_model(model="empirical_first_base" \| "markov")` step |
| P-nucleotide length distributions | Observed rearrangement counts AFTER decomposing NP into N + P | **High** — distinguishing N from P requires palindrome detection against the post-trim coding flank, which is brittle in practice | `estimate_p_nucleotide_lengths(from_rearrangements=..., method="palindrome_detection")` step; **warn loudly when palindrome ambiguity is detected** |
| SHM rates / segment rates / V-subregion rates | Observed rearrangement counts (`v_mutation_count` etc. columns) + V allele subregion annotations | High — per-segment rate inference depends on getting the segment partitioning right; V-subregion rates additionally depend on accurate subregion annotation | `estimate_shm_rates(from_rearrangements=..., per_segment=True, per_v_subregion=True)` step; **warn when subregion annotations are missing** |

---

## 7. Q4 — What should NOT be inferred automatically

Per-decision opt-in catalogue:

| Decision | Why no auto-inference | Builder surface |
|---|---|---|
| **Biological mechanism choices** (use Markov NP base? P-nucleotides? Productive-only?) | User intent, not data-derivable. Same dataset can support multiple biology models. | Explicit step kwargs (`estimate_np_base_model(kind=...)`); never auto-promote `empirical_first_base` → `markov`. |
| **Productive-only rules** | A pipeline decision, not a cartridge property | NOT a builder concern — `productive_only()` lives on `Experiment`. The builder must NOT bake productive-only into the cartridge. |
| **S5F kernel choice** | Per-experiment parameter, NOT cartridge identity (already documented as `models.shm.in_content_hash=False`) | NOT a builder concern — kernels are loaded at simulation time via `_s5f_loader`. |
| **P-nucleotide vs N-nucleotide decomposition** when palindrome ambiguity is high | Brittle inference; false positives produce wrong cartridge | `estimate_p_nucleotide_lengths` must report decomposition confidence; below threshold → emit zero P-plane + warning |
| **Auto-lift of legacy `NP_transitions` / `NP_first_bases` / `p_nucleotide_length_probs`** | Documented boundary; auto-lift would silently change output bytes | NOT a builder concern — same boundary the Markov / P-nucleotide implementation slices respected. |
| **Anchor inference when `gapped_seq` doesn't follow IMGT convention** | Without IMGT numbering, anchor inference is heuristic at best | Builder emits anchor=None + manifest gap report; cartridge author must override explicitly |

---

## 8. Q5 — Desired API shape

### 8.1 Builder skeleton

```python
from GenAIRR.cartridge import ReferenceCartridgeBuilder

builder = ReferenceCartridgeBuilder.from_fasta(
    v_path="v_alleles.fasta",
    j_path="j_alleles.fasta",
    d_path=None,                  # VJ chain; None for VJ
    # OR
    # v_alleles=[Allele(...), ...],  # pre-parsed allele lists
)

builder.infer_identity(
    species="Homo sapiens",
    chain_type="vdj",
    reference_set="IMGT",
    source="user-curated-2026Q1",
)

builder.infer_anchors_from_imgt_gapped()
builder.infer_v_subregions()
builder.infer_rules(
    allowed_bases="ACGT",
    v_anchor_aa="C",
    j_anchor_aa=("W", "F"),
)

# Rearrangement-driven empirical inference (all optional;
# each step writes a build-report entry):
builder.with_rearrangement_data("rearrangements.tsv")
builder.estimate_allele_usage(min_count=10)
builder.estimate_trim_distributions(min_samples=100)
builder.estimate_np_length_distributions(per_region=True)
builder.estimate_np_base_model(kind="empirical_first_base")
builder.estimate_p_nucleotide_lengths(method="palindrome_detection")

# Finalise.
cfg: ga.DataConfig = builder.build()
report: dict = builder.report()
```

### 8.2 Staged semantics

- Every `from_*` / `infer_*` / `estimate_*` step is
  idempotent — calling twice replaces the previous output
  and writes a `replaced` flag to the build report.
- The build report accumulates one entry per stage:
  `{stage: "infer_v_subregions", inputs: {...},
  assumptions: [...], warnings: [...],
  inferred: {...}, rejected: [...]}`.
- `.build()` runs final validation (calls
  `cfg.verify_integrity()` + `cfg.cartridge_manifest()`,
  rejects malformed states with structured errors) and
  returns a `DataConfig` ready for
  `Experiment.on(cfg)`.

### 8.3 Pinned

- `pin_absence_no_reference_cartridge_builder_module_today`
- `pin_absence_no_from_fasta_public_constructor_today`
- `pin_absence_no_estimate_star_inference_steps_today`

---

## 9. Q6 — Build report shape

```python
{
    "schema_version": 1,
    "builder_version": "x.y.z",   # GenAIRR version
    "created_at": "2026-…",       # ISO timestamp
    "stages": [
        {
            "stage": "from_fasta",
            "inputs": {
                "v_path": "...", "v_alleles_parsed": 285,
                "j_path": "...", "j_alleles_parsed": 13,
                "d_path": "...", "d_alleles_parsed": 27,
            },
            "assumptions": [],
            "warnings": ["3 D alleles lack functional_status tag"],
            "inferred": {"v_allele_count": 285, "j_allele_count": 13, "d_allele_count": 27},
            "rejected": [],
        },
        {
            "stage": "infer_v_subregions",
            "inputs": {"derivation": "imgt_gapped"},
            "assumptions": ["IMGT gapped numbering applies"],
            "warnings": ["12 alleles lack gapped_seq — subregion list empty"],
            "inferred": {"alleles_with_subregions": 273, "alleles_without": 12},
            "rejected": [],
        },
        # ...
        {
            "stage": "estimate_trim_distributions",
            "inputs": {"rearrangement_count": 18412, "min_samples": 100},
            "assumptions": ["AIRR TSV column convention"],
            "warnings": [],
            "inferred": {"V_3": {"support_size": 16, "mean": 2.3},
                          "D_5": {...}, "D_3": {...}, "J_5": {...}},
            "rejected": [],
        },
    ],
    "manifest_snapshot": {
        # `cfg.cartridge_manifest()` output at build time
        "models": {...}, "identity": {...}, "catalogue": {...}, ...
    },
    "checksum_at_build_time": "<sha256>",
}
```

This report is attached as `cfg.build_report` (the existing
field, which is currently `None` on all bundled cartridges)
and exposed via `builder.report()`.

### Pinned

- `pin_present_data_config_build_report_field_exists_but_is_none_today`
- `pin_absence_no_build_report_producer_today`

---

## 10. Q7 — Integration with the existing pipeline

The builder's output is a **plain `DataConfig`**. No
parallel cartridge object, no engine-only side channel. This
preserves:

- Drop-in replaceability with bundled cartridges
  (`Experiment.on(my_cfg)`).
- Manifest discoverability — `cfg.cartridge_manifest()`
  works verbatim and gains a (truthy) `build_report` summary
  block.
- Compile-time validation — the existing
  `dataconfig_to_refdata` bridge validates the cartridge
  shape and the compile-gate test
  ([`test_refdata_validation_compile_gate.py`](../tests/test_refdata_validation_compile_gate.py))
  surfaces malformed cartridges.
- Round-trip preservation — `pickle.dump(cfg)` /
  `pickle.load(path)` works because `DataConfig` is a normal
  dataclass; `build_report` rides through the pickle.

The builder does NOT touch:

- The Rust `RefDataConfig` direct-construction API. That
  surface (`add_v_allele(...)`) stays usable for the
  test/perf workloads it exists for.
- The bundled-pickle loader. The 106 bundled cartridges
  predate the new builder; they continue to load with
  `build_report=None` (treat absence as "legacy cartridge,
  no provenance recorded").

### Pinned

- `pin_scaffold_dataconfig_to_refdata_bridge_validates_at_compile_time`
- `pin_scaffold_pickle_round_trip_preserves_build_report`

---

## 11. Implementation order — **[v1 shipped]**

### 11.1 v1 scope — what landed

Per the user brief, v1 implements only the facade + report
infrastructure + identity / V-subregion stages. The implementation
landed in the cartridge-builder slice:

- [`src/GenAIRR/cartridge_builder.py`](../src/GenAIRR/cartridge_builder.py)
  — ~570 lines: `CartridgeBuildReport` dataclass with
  `to_dict()`; `ReferenceCartridgeBuilder` class with
  `from_fasta` constructor + 4 stages (`infer_identity`,
  `infer_v_subregions`, `with_rules`, `with_models`) +
  `build()` + `report()`. Private `_SafeAnchorMixin` +
  `_BuilderVAllele` / `_BuilderJAllele` subclasses degrade
  gracefully when the native `_native._anchor` C resolver
  isn't built.
- Top-level exports in
  [`src/GenAIRR/__init__.py`](../src/GenAIRR/__init__.py).
- Dead-reference cleanup at the two live sites (§3).

### 11.2 Estimator slice ordering — shipped + deferred

The user brief originally deferred every data-derived
statistical-estimator method. The first one has shipped;
the remaining five stay deferred, pinned by
`test_v1_defers_remaining_statistical_estimators` in the
implementation test file. Each call to a deferred name
raises `AttributeError`. The shipped slice updated the
audit doc + present-state pins in lockstep with each
landing.

**Shipped (1 of 6):**

1. `estimate_allele_usage(rearrangements, *, min_count=1.0,
   ambiguous="fractional", replace=True)` — **[Shipped]**.
   Reads AIRR `v_call` / `d_call` / `j_call` from
   `list[dict]` / path / open text handle; writes
   per-segment weights into the typed
   `ReferenceEmpiricalModels.allele_usage` plane;
   recombination uses estimated weights by default unless
   overridden by `Experiment.recombine(v_allele_weights=...)`.
   Three ambiguity policies (`"fractional"` /
   `"truth_first"` / `"reject"`). See
   [`docs/allele_usage_estimation_design.md`](allele_usage_estimation_design.md).

**Deferred (5 of 6):** in audit-recommended order —

2. `estimate_trim_distributions(from_rearrangements=...,
   min_samples=...)` — next in priority because trim
   distributions are heavily relied on by every downstream
   biology slice (NP length, P-nucleotide, productive-only
   composition). Maps directly to the existing
   `ReferenceEmpiricalModels.trims` plane.
3. `estimate_np_length_distributions(from_rearrangements=...)` —
   straightforward NP1/NP2 histogram from rearrangement data.
4. `estimate_np_base_model(kind="empirical_first_base"|"markov")` —
   recommended kind defaults to `"empirical_first_base"`;
   Markov inference needs care for short NP regions.
5. `estimate_p_nucleotide_lengths(method="palindrome_detection")` —
   highest risk per audit §6 (P-vs-N ambiguity).
6. `estimate_shm_rates(per_segment=True, per_v_subregion=True)` —
   needs the V-subregion annotation to be present (currently
   handled by the v1 `infer_v_subregions` step).

Each estimator's implementation slice MUST:

- Append a new stage entry to `CartridgeBuildReport.stages`
  with the documented `{stage, inputs, inferred, warnings}`
  shape.
- Surface stage-level confidence / sample-size warnings under
  `warnings` rather than `rejected` (the latter is reserved
  for per-allele / per-record rejections at parse time).
- Validate output through the existing
  `ReferenceEmpiricalModels` / `NpBaseModelSpec` /
  `EmpiricalDistributionSpec` typed planes so the resulting
  cartridge survives `dataconfig_to_refdata` bridge
  validation.

### 11.3 Implementation order (recommended) — original audit text

A single self-contained slice can land the builder. Six
sub-steps:

1. **New module** `src/GenAIRR/cartridge/__init__.py` (or
   `src/GenAIRR/cartridge_builder.py` — module name is a
   bikeshed; whichever the user prefers). Exports
   `ReferenceCartridgeBuilder` + a tiny
   `CartridgeBuildReport` dataclass.

2. **Constructors**: `from_fasta(v_path, d_path, j_path)` +
   `from_allele_lists(v_alleles, d_alleles, j_alleles)` +
   `from_dataconfig(cfg)` (a re-entry-point so a partially-
   built cartridge can be "edited" by re-running steps).

3. **Identity stages**: `infer_identity` + `infer_anchors_from_imgt_gapped`
   + `infer_v_subregions` + `infer_rules` (all working off
   pure-FASTA + IMGT-gap-derivable info; no rearrangement
   data required).

4. **Rearrangement-driven stages** (require
   `with_rearrangement_data(path_or_records)` first):
   `estimate_allele_usage` /
   `estimate_trim_distributions` /
   `estimate_np_length_distributions` /
   `estimate_np_base_model` /
   `estimate_p_nucleotide_lengths` /
   `estimate_shm_rates`.

5. **Build + report**: `build()` runs final validation +
   computes the build report (the report is a property that
   accumulates entries as stages run; `build()` snapshots it
   onto `cfg.build_report`).

6. **Dead-reference cleanup**: update the
   `DataConfig.build_report` docstring + update
   `.private/scripts/build_imgt_configs.py` to import from
   the new module.

Cost estimate:

- ~400 lines Python (builder + 10 stages + report dataclass)
- ~120 lines tests (per-stage smoke + idempotency + report
  shape + round-trip)
- ~30 lines docstring + private-script update

Single self-contained slice.

---

## 12. Test surface — what this audit pins

Mirrored in
[`tests/test_reference_cartridge_authoring_contract.py`](../tests/test_reference_cartridge_authoring_contract.py).

### `pin_scaffold_*` — pre-existing live surfaces the new
builder reuses

1. `parse_fasta` is public and parses `>header\nbases` correctly.
2. `compute_v_region_boundaries` derives V subregions from
   `gapped_seq`.
3. `DataConfig.cartridge_manifest()` exists and returns a
   JSON-clean dict.
4. `DataConfig.verify_integrity()` + `compute_checksum()`
   gate corrupted pickles.
5. `RefDataConfig.{vj,vdj}()` + `add_v_allele` /
   `add_d_allele` / `add_j_allele` are the existing manual
   construction path.
6. `ReferenceEmpiricalModels` accepts every typed plane
   (`np_lengths`, `trims`, `np_bases`, `p_nucleotide_lengths`).
7. `ReferenceRulesSpec` is the rules-authoring entry point.
8. `dataconfig_to_refdata` bridge validates malformed
   cartridges at compile time.
9. 106 bundled cartridges load via lazy `__getattr__`.

### `pin_present_*` — dead-reference inventory

10. `DataConfig.build_report` field exists and is `None` on
    every bundled cartridge (no producer in live source).
11. `DataConfig.build_report` docstring still names
    `RandomDataConfigBuilder` — a dead class.
12. `.private/scripts/build_imgt_configs.py` still imports
    from the dead module path.

### `pin_absence_*` — the gaps the slice closes

13. No public `ReferenceCartridgeBuilder` class or
    `cartridge` / `cartridge_builder` module.
14. No public `from_fasta` constructor on any cartridge
    surface.
15. No public `from_airr` / `from_rearrangement_tsv`
    constructor.
16. No `infer_v_subregions` / `infer_anchors_from_imgt_gapped`
    public step methods.
17. `estimate_allele_usage` **[Shipped — post-v1]**;
    `estimate_trim_distributions` /
    `estimate_np_length_distributions` /
    `estimate_np_base_model` / `estimate_p_nucleotide_lengths`
    / `estimate_shm_rates` remain deferred (pinned by
    `test_v1_defers_remaining_statistical_estimators`).
18. No `CartridgeBuildReport` dataclass / API today.
19. No public surface that writes to `DataConfig.build_report`.

### `pin_scaffold_*` — historical builder gone from live source

20. `from GenAIRR.dataconfig.make.random import
    RandomDataConfigBuilder` raises `ModuleNotFoundError`.
21. `from GenAIRR.dataconfig.make import
    CustomDataConfigBuilder` raises `ModuleNotFoundError`.
22. The build-cache mirror at
    `docs/build/lib.../GenAIRR/dataconfig/make/random.py`
    is NOT on the import path (verified by raising
    ModuleNotFoundError).

### Doc anchor

23. Audit doc exists and references the contract file;
    section structure intact.

---

## 13. Out of scope

Documented here so a future implementer doesn't accidentally
expand the work.

- **Cartridge round-trip via JSON / YAML.** Pickle is the
  durable format; the builder doesn't introduce a parallel
  serialisation layer.
- **Cartridge diff / merge tooling.** A future
  `CartridgeDiff` surface could compare two cartridges' build
  reports — separate slice.
- **Live-data download (IMGT / OGRDB / etc.).** The
  historical `.private/scripts/build_imgt_configs.py`
  downloads FASTA from the web. The new builder accepts file
  paths or pre-parsed allele lists; download tooling is a
  separate concern (the private script can stay private,
  updated to import the new builder).
- **GUI / web wrapper.** Out of scope.
- **Auto-bundling.** The builder produces a `DataConfig`
  ready for pickling, but the wheel-build process that ships
  the 106 bundled cartridges is a separate Makefile step.
- **Cross-locus / cross-species merging.** Each builder
  invocation produces ONE cartridge.

---

## 14. Summary table

| Concern | Today | Recommendation |
|---|---|---|
| Public cartridge-creation API | 4 paths (bundled load, pickle load, manual `DataConfig`, manual `RefDataConfig`) | Add 5th path: `ReferenceCartridgeBuilder` (FASTA / AIRR data → `DataConfig`) |
| Inference-adjacent helpers in live source | `parse_fasta`, `compute_v_region_boundaries`, `_resolve_v_subregions`, `cartridge_manifest`, `verify_integrity` | All reused by the new builder verbatim |
| Historical builder (`RandomDataConfigBuilder` + `CustomDataConfigBuilder` + auxiliary builders) | **REMOVED** from live source; 3 dead-reference sites remain | Inventory pinned; cleanup lockstepped with implementation slice |
| `DataConfig.build_report` field | Exists; unconditionally `None` on every bundled cartridge | Builder populates with per-stage entries |
| Build-report shape | Not defined | `CartridgeBuildReport` dataclass per §9 |
| Identity inference (`species`, `chain_type`, `reference_set`) | Manual `ConfigInfo` construction | `infer_identity(species=..., chain_type=...)` stage |
| V-allele subregion inference | `_resolve_v_subregions` (internal) called by bridge | `infer_v_subregions()` stage — public step |
| Allele usage inference | **[Shipped — post-v1]** `estimate_allele_usage(rearrangements, *, min_count, ambiguous, replace)` stage; recombination uses estimated weights by default. See [`docs/allele_usage_estimation_design.md`](allele_usage_estimation_design.md). | — (shipped) |
| Trim distribution inference | None | `estimate_trim_distributions(...)` stage |
| NP length / base / P-nuc inference | None — must author typed planes manually | `estimate_np_length_distributions` + `estimate_np_base_model` + `estimate_p_nucleotide_lengths` stages |
| SHM rate inference | None | `estimate_shm_rates(per_segment=True, per_v_subregion=True)` stage with warnings |
| Stop-and-report condition (hidden live builder) | **NOT MET** — historical builder is gone, three dead references remain | Audit proceeds; cleanup three sites in lockstep with the implementation slice |

The cartridge-authoring story is **architecturally tractable**:
the live source already exposes every inference primitive the
new builder needs (`parse_fasta`,
`compute_v_region_boundaries`, the typed plane validators,
`cartridge_manifest`, `verify_integrity`), and the missing
piece is a staged builder facade that composes them with a
build-report audit trail. The recommended next step is the
single implementation slice per §11, with the three
dead-reference cleanup actions ridden in lockstep.
