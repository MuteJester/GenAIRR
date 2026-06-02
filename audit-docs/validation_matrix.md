# GenAIRR validation matrix

This document is the **navigable index** of GenAIRR's correctness,
reproducibility, distribution, and performance guarantees.

Each row of the matrix below points to:
- the **audit doc** that documents *what* the engine guarantees and
  *why* (with derivation rules, edge cases, drift items);
- the **test file** that pins the guarantee with golden tests so a
  regression fails CI before it ships;
- (where relevant) the **Rust-side test module** that pins the
  kernel-level invariant.

Use it when you:
- Want to know which audit covers a behaviour you're changing.
- Want to know which tests must stay green for a given guarantee.
- Want to know what *isn't* yet validated (look at each audit's §6
  drift section).

---

## 1. The validation matrix

### 1.0 Two-layer integrity model

Engine integrity is checked at **two independent layers**:

| Layer | Question | API | Audience |
|-------|----------|-----|----------|
| **AIRR validator** | Does the *projected* AIRR record agree with an independent re-derivation from the Outcome? | `result.validate_records(refdata)` → `ValidationReport` | User-facing. The downstream contract. |
| **Cache parity** | Does the *cached* `SegmentLiveCall` on the final Simulation equal a from-scratch recompute? | `outcome.check_live_call_cache_parity(refdata)` → `List[dict]` | Engine-side. The integrity guard. |

**Why two layers.** The validator catches projection-layer bugs and indirectly catches cache-layer bugs that leak into projection. The parity harness catches cache-layer bugs directly — closer to the source. Three real refresh-path bugs (the `segment_region_overlaps_dirty` strict-`<`, NP-claim over-extension, and `on_indel_inserted` `<=`) were each caught by one layer or the other; keeping them separate localised the failures to the right code path quickly.

Use BOTH gates side by side in release-tier CI. The canonical example lives in [`tests/test_release_validation.py::test_validator_and_parity_are_independent_layers`](../tests/test_release_validation.py).

| Guarantee | Audit doc | Python tests | Rust kernel |
|-----------|-----------|--------------|-------------|
| **Projected record consistency** | [`airr_record_validator.md`](airr_record_validator.md) | [`test_airr_record_validator.py`](../tests/test_airr_record_validator.py), [`test_release_validation.py`](../tests/test_release_validation.py) | `engine_rs/src/airr_record/validate.rs` |
| **Cached live-call parity** | (this section) | [`test_live_call_cache_parity.py`](../tests/test_live_call_cache_parity.py) | `engine_rs/src/live_call/cache_parity.rs` |

### Standard postcondition pattern

Every realistic pipeline can be capped with a one-line validation
guard. This is the recommended pattern for CI, notebooks, and
release QA:

```python
result = exp.run_records(n=1000, seed=0)
report = result.validate_records(refdata)
assert report, report.summary()
```

`report` is a `ValidationReport` (`from GenAIRR import
ValidationReport`) with `.ok`, `.count`, `.failures`, and
`.summary()` (a histogram of issue kinds across failures).
`.failures` carries `{record_index, sequence_id, issues}` per
failing record; issue dicts are JSON-serializable for CI
artifacts.

Live-tier example in
[`tests/test_release_validation.py`](../tests/test_release_validation.py)
— runs the full IGK/IGL productive stack + a non-productive variant
+ an IGH stack composed with every audit mechanism (D inversion,
receptor revision, paired-end). The earlier "D-tie oracle filtered
out under inversion" exception was removed when the
orientation-aware walker landed; the inverted-D release tests now
run clean without exceptions.

### Contributor checklist

When adding a new mechanism / pass, treat record validation as a
required step in the contributor flow:

1. **Pass emits events** through `SimulationBuilder` (no direct
   `sim.with_*` in production code).
2. **Replay validates** — the trace round-trip test in your audit
   reproduces sequence + AIRR coords + per-pass event count.
3. **Distribution invariant tested** if the mechanism is
   stochastic (MC with ±5σ tolerance or exact enumeration).
4. **Biological stress matrix updated** if the mechanism affects
   productivity / anchors / junction.
5. **`result.validate_records(refdata)` passes** on a
   representative stack. This is the postcondition gate — the
   engine's own truth oracle must agree with the projected AIRR
   record under your new mechanism.

### 1.1 Per-mechanism guarantees

| Guarantee                                  | Audit doc                                                                         | Python tests                                                                          | Rust kernel tests                                                                                       |
|--------------------------------------------|------------------------------------------------------------------------------------|---------------------------------------------------------------------------------------|----------------------------------------------------------------------------------------------------------|
| **Productive validity**                    | [`productive_failure_mode_audit.md`](productive_failure_mode_audit.md)             | [`test_productive_failure_modes.py`](../tests/test_productive_failure_modes.py) (13 tests), [`test_productive_stress_matrix.py`](../tests/test_productive_stress_matrix.py) | `engine_rs/src/contract/productive_junction_frame.rs::tests` (12 unit tests), `engine_rs/src/contract/no_stop_codon_in_junction/tests/` (4 modules), `engine_rs/src/contract/anchor_preserved/tests/` |
| **Indel provenance**                       | [`indel_provenance_audit.md`](indel_provenance_audit.md)                           | [`test_indel_provenance.py`](../tests/test_indel_provenance.py) (26 tests)             | `engine_rs/src/passes/corrupt/indel/tests/constraints.rs` (8 tests incl. distribution oracles)            |
| **End-loss / primer-trim provenance**      | [`primer_trim_end_loss_audit.md`](primer_trim_end_loss_audit.md)                   | [`test_primer_trim_end_loss_provenance.py`](../tests/test_primer_trim_end_loss_provenance.py)            | `engine_rs/src/passes/end_loss/`                                                                         |
| **D-segment inversion** (V(D)J inversion event; `d_inverted` AIRR field) | [`d_inversion_design.md`](d_inversion_design.md) | [`test_invert_d_dsl.py`](../tests/test_invert_d_dsl.py), [`test_d_inversion_contract.py`](../tests/test_d_inversion_contract.py), inverted-D stack in [`test_release_validation.py`](../tests/test_release_validation.py) | `engine_rs/src/passes/invert_d.rs` (20 unit tests: prob + replay + event + metadata), `engine_rs/src/airr_record/tests/projection.rs` (Slice E unit tests for the field) |
| **Receptor revision** (post-recombine V replacement; `receptor_revision_applied` / `original_v_call` AIRR fields) | [`receptor_revision_design.md`](receptor_revision_design.md) | [`test_receptor_revision_dsl.py`](../tests/test_receptor_revision_dsl.py), [`test_receptor_revision_contract.py`](../tests/test_receptor_revision_contract.py), revised-V stack + Bernoulli draw in [`test_release_validation.py`](../tests/test_release_validation.py) / [`test_distribution_invariants.py`](../tests/test_distribution_invariants.py) | `engine_rs/src/passes/receptor_revision.rs` (17 unit tests: prob + replay + event + metadata), `engine_rs/src/live_call/refresh_plan.rs::tests` (Slice B step-mapping + dedup), `engine_rs/src/live_call/refresh_hook.rs::tests` (Slice B AllStructural-equivalent cache parity), `engine_rs/src/airr_record/tests/projection.rs` (Slice E unit tests for the fields + validator) |
| **Paired-end / read layout** (R1 / R2 windows + insert size; eight new AIRR fields `read_layout` / `r1_sequence` / `r2_sequence` / `r1_start/end` / `r2_start/end` / `insert_size`) | [`paired_end_design.md`](paired_end_design.md) | [`test_paired_end_dsl.py`](../tests/test_paired_end_dsl.py), [`test_paired_end_contract.py`](../tests/test_paired_end_contract.py), [`test_paired_end_schema.py`](../tests/test_paired_end_schema.py), paired-end stack + insert-size invariant in [`test_release_validation.py`](../tests/test_release_validation.py) / [`test_distribution_invariants.py`](../tests/test_distribution_invariants.py) | `engine_rs/src/passes/paired_end.rs` (12 unit tests: sampling + replay + relationship validation), `engine_rs/src/airr_record/sequence.rs` (Slice B projection kernel + bounds errors), `engine_rs/src/airr_record/tests/projection.rs` (Slice A/B/C projection + validator tests covering all eight fields + the five `Read*` issue variants) |
| **Targeted SHM** (per-segment rate scalars on the V / D / J / NP buckets; `Experiment.mutate(segment_rates=…)`) | [`shm_segment_rate_design.md`](shm_segment_rate_design.md) (architecture), [`shm_model_audit.md`](shm_model_audit.md) (broader SHM contract) | [`test_shm_segment_rate_contract.py`](../tests/test_shm_segment_rate_contract.py), [`test_segment_rates_implementation.py`](../tests/test_segment_rates_implementation.py), targeted-SHM stack + zero-rate exclusion invariant in [`test_release_validation.py`](../tests/test_release_validation.py) / [`test_distribution_invariants.py`](../tests/test_distribution_invariants.py) | `engine_rs/src/passes/mutate/segment_rates.rs` (`SegmentRateWeights` + `segment_at_position` + 3 unit tests), `engine_rs/src/passes/mutate/s5f/sampling.rs` (S5F profile-build weighting), `engine_rs/src/passes/mutation_transaction/substitution.rs` (uniform fast/constrained/replay paths) |
| **Mutation provenance counters** (per-segment biological-SHM counters partitioning `n_mutations` across V / D / J / NP; new AIRR fields `n_v_mutations` / `n_d_mutations` / `n_j_mutations` / `n_np_mutations`) | [`mutation_provenance_audit.md`](mutation_provenance_audit.md) | [`test_mutation_provenance_contract.py`](../tests/test_mutation_provenance_contract.py), [`test_per_segment_mutation_counters.py`](../tests/test_per_segment_mutation_counters.py), per-segment counter stack + sum-invariant cross-check in [`test_release_validation.py`](../tests/test_release_validation.py) | `engine_rs/src/airr_record/builder.rs` (event-ledger aggregation filtered to `mutate.{uniform,s5f}` pass names, NP1+NP2 rollup), `engine_rs/src/airr_record/validate.rs` (five new issue kinds: `N{V,D,J,Np}MutationsMismatch` + `MutationCountSumMismatch`), `engine_rs/src/python/outcome.rs` (record dict + issue serialisation with documented `details.source` tags) |
| **V-subregion cartridge annotations** (per-V-allele IMGT FWR1/CDR1/FWR2/CDR2/FWR3 intervals derived from `gapped_seq`; surfaced as Python `Allele.subregions`, Rust `VSubregion`/`VSubregionLabel`, manifest `models.shm.v_subregion_support`, and `refdata_content_hash`) | [`v_region_substructure_audit.md`](v_region_substructure_audit.md) (Slice 1 shipped; counters slice deferred) | [`test_v_region_substructure_contract.py`](../tests/test_v_region_substructure_contract.py) (scaffold + present + remaining-absence pins), [`test_v_subregion_annotations.py`](../tests/test_v_subregion_annotations.py) (17 spec tests: coverage, monotonic+bounds, content-hash flip, manifest shape, legacy zero-coverage, user override round-trip, malformed-input rejection), V-subregion sanity check in [`test_release_validation.py`](../tests/test_release_validation.py) | `engine_rs/src/refdata.rs` (`VSubregion` + `VSubregionLabel`, allele field), `engine_rs/src/python/refdata.rs` (`add_v_allele(subregions=…)` + `parse_subregions` bridge-time validator + `PyAllele.subregions` getter), `engine_rs/src/trace_file.rs` (`refdata_content_hash` folds `label:start-end` per allele); Python bridge: [`src/GenAIRR/_refdata_resolver.py`](../src/GenAIRR/_refdata_resolver.py) (`_resolve_v_subregions` derives via `compute_v_region_boundaries` or honours user-set `allele.subregions`) |
| **V-subregion SHM targeting** (per-V-subregion rate scalars on the FWR1 / CDR1 / FWR2 / CDR2 / FWR3 labels — `Experiment.mutate(v_subregion_rates={…})`. Composes multiplicatively with `segment_rates` on V sites; non-V sites and unannotated V alleles receive identity factor `1.0`. Default rates / explicit all-ones / `FWR` / `CDR` aliases with explicit-label override. Rates fold into the plan signature so replay against a different vector fails before consuming choices; rates do NOT enter `refdata_content_hash` — same per-experiment boundary as `segment_rate_support`.) | [`v_subregion_shm_rate_design.md`](v_subregion_shm_rate_design.md) (Slices A + B shipped; counters shipped — separate row), [`shm_segment_rate_design.md`](shm_segment_rate_design.md) (broader SHM rate architecture) | [`test_v_subregion_shm_rate_contract.py`](../tests/test_v_subregion_shm_rate_contract.py) (scaffold + present pins for kwarg / validator / Rust struct / helper / manifest / compile-time rejection), [`test_v_subregion_rates_implementation.py`](../tests/test_v_subregion_rates_implementation.py) (27 spec tests: DSL validation, alias expansion + override, FWR-only / CDR-only, zero-rate CDR / FWR exclusion via baseline-diff classifier, non-V sites unaffected, both uniform + S5F models, productive_only triad preserved, matching + mismatched replay, manifest advertisement), V-subregion-targeting stack + zero-rate exclusion invariant in [`test_release_validation.py`](../tests/test_release_validation.py) / [`test_distribution_invariants.py`](../tests/test_distribution_invariants.py) | `engine_rs/src/passes/mutate/v_subregion_rates.rs` (`VSubregionRateWeights` + `v_subregion_at_position` helper), `engine_rs/src/passes/mutation_transaction/substitution.rs` (`combined_site_factor` unifies segment × V-subregion factor for the constrained + fast paths + replay validation), `engine_rs/src/passes/mutate/s5f/sampling.rs` (`build_profile` calls `combined_site_factor`), `engine_rs/src/passes/mutate/s5f/execution.rs` (S5F replay-time admissibility check), `engine_rs/src/passes/mutate/uniform.rs` (uniform `substitute_position_constrained` call site), `engine_rs/src/passes/paramsig.rs` (`fmt_v_subregion_rates` for plan signature), `engine_rs/src/python/plan.rs` (`push_mutate_*(v_subregion_rates=…)`); Python bridge: [`src/GenAIRR/experiment.py`](../src/GenAIRR/experiment.py) (`_validate_v_subregion_rates` + alias-expansion + `_check_v_subregion_rates_satisfiable`), [`src/GenAIRR/_compile.py`](../src/GenAIRR/_compile.py) (`_lower_mutate`) |
| **V-subregion mutation counters** (per-V-IMGT-subregion biological-SHM counters partitioning `n_v_mutations` across the five canonical labels plus `n_v_unannotated_mutations` for V SHM events that fall outside every subregion interval — V-side CDR3 stretch, unannotated V alleles, or indel-inserted V bases. Six new AIRR fields, all `int`, always present, sum exactly to `n_v_mutations` by construction. Aggregated by walking `outcome.events()`, filtering to the same `mutate.{uniform,s5f}` pass-name allowlist as the per-segment counters, and matching each V `BaseChanged.germline_pos: Option<u16>` against the assigned V allele's `subregions` table. Validator enforces six per-field equalities plus a cross-field sum invariant.) | [`v_subregion_mutation_counters_audit.md`](v_subregion_mutation_counters_audit.md) (slice shipped) | [`test_v_subregion_mutation_counters_contract.py`](../tests/test_v_subregion_mutation_counters_contract.py) (scaffold + present pins for fields / validator issue kinds / manifest block; two-bucket aggregates remain absent), [`test_v_subregion_mutation_counters_implementation.py`](../tests/test_v_subregion_mutation_counters_implementation.py) (13 spec tests: baseline zero, per-cartridge partition invariant, zero-rate CDR / FWR composition, full-stack validation, replay round-trip, mismatch detection through engine-projection contract, unannotated routing on stripped cartridge, corruption isolation, DataFrame columns, manifest advertisement), V-subregion-counters full-stack + two-layer partition invariant in [`test_release_validation.py`](../tests/test_release_validation.py) (`test_v_subregion_mutation_counters_full_stack_partition_and_validate` + extended `test_v_subregion_rates_trace_replay_round_trip_passes_validator` covering all six new fields) | `engine_rs/src/airr_record/record.rs` (six new `n_*_mutations` fields on `AirrRecord`), `engine_rs/src/airr_record/builder.rs` (hoisted V `subregions` lookup + per-event subregion dispatch in the existing `mutate.{uniform,s5f}` filter loop), `engine_rs/src/airr_record/validate.rs` (seven new issue variants: `N{Fwr1,Cdr1,Fwr2,Cdr2,Fwr3,VUnannotated}MutationsMismatch` + `VSubregionMutationCountSumMismatch` with independent recompute), `engine_rs/src/python/outcome.rs` (six new `dict.set_item` calls + serialisation for the seven new issue kinds with `source` tags), `engine_rs/src/airr_record/tests/projection.rs` (four Rust unit tests covering event-time attribution, `germline_pos=None` routing, tampered-counter detection, partition-sum mismatch detection); Python: [`src/GenAIRR/result.py`](../src/GenAIRR/result.py) (canonical column order), [`src/GenAIRR/dataconfig/data_config.py`](../src/GenAIRR/dataconfig/data_config.py) (`v_subregion_counter_support` manifest block) |
| **Paired-end FASTQ export** (`SimulationResult.to_paired_fastq(r1_path, r2_path, *, quality="illumina", overwrite=False, **quality_kwargs)` writes two synchronized FASTQ files — one R1 record + one R2 record per AIRR record — sourced verbatim from the `r1_sequence` / `r2_sequence` projection fields. Headers use `@{sequence_id}/1` and `@{sequence_id}/2` (universally-portable suffix; no `\|`-pipe metadata that some aligners split on). R2 is **already** reverse-complemented at projection time (the `PairedEndWindowMismatch { side: R2 }` validator invariant enforces it); the writer outputs it verbatim and applies no second flip. Reuses the existing pluggable quality models (`ConstantQualityModel` / `IlluminaQualityModel` from `_qmodel.py`) consumed by the single-end `to_fastq`; each read is scored independently so R1 and R2 get their own ramp-up + tail-down shape from position 0. Refuses to clobber existing files unless `overwrite=True`; raises on non-paired records (`read_layout != "paired_end"`), empty R1/R2 windows, and length-mismatched quality arrays. Pure projection — no engine, validator, trace, manifest, or plan-signature touchpoints.) | [`fastq_export_design.md`](fastq_export_design.md) (slice shipped), [`paired_end_design.md`](paired_end_design.md) (the projection plumbing the writer reads from) | [`test_fastq_export_contract.py`](../tests/test_fastq_export_contract.py) (scaffold + present pins for the method + signature; absence pins remain on record-generator surface, manifest block, and trace addresses), [`test_to_paired_fastq.py`](../tests/test_to_paired_fastq.py) (13 spec tests: file record count, header format `@{sequence_id}/1` + `/2`, R1/R2 body equality against AIRR fields, quality length parity, constant + illumina quality models, non-paired rejection, empty R1/R2 rejection, overwrite guard with single-side pre-existence, synchronized R1/R2 sequence_id), paired-end FASTQ export smoke + decodability in [`test_release_validation.py`](../tests/test_release_validation.py) | `src/GenAIRR/result.py` (the `to_paired_fastq` method on `SimulationResult`), `src/GenAIRR/_qmodel.py` (reused unchanged — `ConstantQualityModel`, `IlluminaQualityModel`, `resolve_quality_model`, `phred_to_ascii`). **Zero Rust changes**; pure Python export layer on already-projected AIRR fields. |
| **NP base models / Markov N-addition** (cartridge-authored NP-region base sampling via `NpBaseModelSpec(kind="uniform"\|"empirical_first_base"\|"markov")` on `ReferenceEmpiricalModels.np_bases`. Three concrete generators behind the new `NpBaseGenerator` trait: `UniformNpGenerator` (legacy 4-way ACGT, byte-identical signature to pre-typed-model baseline), `CategoricalNpGenerator` (position-independent weighted A/C/G/T, byte-identical to pre-Markov-slice `CategoricalBase` signature), and `MarkovBaseGenerator` (1-step previous-base-conditional with `first_base: [f64;4]` + `transitions: [[f64;4];4]` in canonical A/C/G/T from/to order). Per-position support is materialised via `base_generator.support(position, previous)`; the `GenerateNPPass` loop threads a `previous: Option<u8>` set to `Some(base)` only after the trace record is committed. **No new trace addresses** — replay reconstructs `previous` from `np.npN.bases[i-1]`. Plan signature folds the full Markov payload (5 rows × A/C/G/T canonical order) so two cartridges with identical `first_base` but different transitions fail the signature gate before any choice is consumed. Productive-only composes through the admit-mask × generator-support intersection; mid-stream `b'N'` sentinel falls back to `first_base`. **Legacy `NP_transitions` / `NP_first_bases` auto-lift remains deferred** — orphan fields stay in `_DOCUMENTED_ORPHAN_DATACONFIG_FIELDS`; manifest reports `legacy_fallback=False`.) | [`junction_n_addition_audit.md`](junction_n_addition_audit.md) (typed model + Markov shipped; legacy auto-lift deferred), [`np_markov_base_generator_design.md`](np_markov_base_generator_design.md) (Markov implementation surface table) | [`test_np_base_model_implementation.py`](../tests/test_np_base_model_implementation.py) (typed model spec + lowering + manifest + replay), [`test_np_markov_base_generator_contract.py`](../tests/test_np_markov_base_generator_contract.py) (20 post-slice pins: trait + generators + previous-base param + bridge kwarg + manifest flip + legacy orphan preservation), [`test_np_markov_base_generator_implementation.py`](../tests/test_np_markov_base_generator_implementation.py) (13 behaviour tests: lowering, deterministic A→T→G→C walk, strong-matrix convergence, replay round-trip, signature-gate mismatch, productive-only triad, byte-identical legacy signatures, manifest flip), [`test_np_markov_release.py`](../tests/test_np_markov_release.py) (release-tier: productive IGH full stack + replay round-trip + dependency invariant) | `engine_rs/src/passes/generate_np/np_base_generator.rs` (trait + `UniformNpGenerator` + `CategoricalNpGenerator` + `MarkovBaseGenerator` + 8 unit tests covering signature shapes, support-by-previous, mid-stream sentinel fallback, row-weight validation, signature divergence), `engine_rs/src/passes/generate_np.rs` (`base_generator: Box<dyn NpBaseGenerator>` field + `with_generator(…)` constructor + `parameter_signature` fold), `engine_rs/src/passes/generate_np/execution.rs` (`let mut previous` loop-local + post-commit update), `engine_rs/src/passes/generate_np/sampling.rs` (`previous: Option<u8>` threaded through `sample_base` + `validate_replayed_np_base` + private `SupportPairsDist` adapter feeding the generic `sample_base_with_admit_mask` / `sample_filtered_with_policy` helpers verbatim), `engine_rs/src/python/plan.rs` (`push_generate_np(..., markov_transitions=…)` kwarg + `parse_canonical_base_weights` row-completeness validator); Python bridge: [`src/GenAIRR/_dataconfig_extract.py`](../src/GenAIRR/_dataconfig_extract.py) (`_np_bases_from_models` + `_np_markov_transitions_from_models`), [`src/GenAIRR/_normalize.py`](../src/GenAIRR/_normalize.py) (`_to_immutable_byte_pair_matrix`), [`src/GenAIRR/_pipeline_ir.py`](../src/GenAIRR/_pipeline_ir.py) (`_RecombineStep.np{1,2}_markov_transitions`), [`src/GenAIRR/_compile.py`](../src/GenAIRR/_compile.py) (`_lower_recombine` thread), [`src/GenAIRR/dataconfig/data_config.py`](../src/GenAIRR/dataconfig/data_config.py) (manifest `np_base_models.supported_kinds = ["uniform", "empirical_first_base", "markov"]`, `deferred_kinds = []`) |
| **P-nucleotide / palindromic addition** (cartridge-authored per-end P-nucleotide length distributions via `ReferenceEmpiricalModels.p_nucleotide_lengths`, keyed by junction side label `"V_3"` / `"D_5"` / `"D_3"` / `"J_5"`. v1 ships **lengths-only**: each `PAdditionPass` records one `Int(length)` choice at `p.{end}.length`; bases derive deterministically from `(allele, trim, orientation, length)` via `complement_base` — 3' extensions (V_3 / D_3) reverse-complement the last `length` post-trim coding bytes; 5' extensions (D_5 / J_5) reverse-complement the first `length`. P-bytes are pushed into the pool with `Nucleotide::flags::P_NUC` set; the descriptive `Region` flows through the new `SimulationEvent::PRegionAdded { end, region }` event only — **no structural region added to `sim.sequence.regions`**, so existing `regions.iter().find(|r| r.segment == ...)` projection sites stay correct (one structural region per biological segment). Pipeline ordering is post-trim: `assemble.V → p_addition.V_3 → generate_np.NP1 → [invert_d] → p_addition.D_5 → assemble.D → p_addition.D_3 → generate_np.NP2 → p_addition.J_5 → assemble.J`. **Critical IR ordering: `p_addition.D_5` runs AFTER `invert_d` commits** so D's effective_seq is read under the post-inversion orientation. AIRR projection surfaces four new int fields (`p_v_3_length` / `p_d_5_length` / `p_d_3_length` / `p_j_5_length`); the validator's `PLengthMismatch { end, reported, event_count }` issue kind catches tampered records via an event-ledger recompute. Plan signature folds each pass's per-end length distribution via `fmt_int_dist`. **Legacy `DataConfig.p_nucleotide_length_probs` is NOT auto-lifted** — orphan field stays in `_DOCUMENTED_ORPHAN_DATACONFIG_FIELDS`, manifest reports `legacy_fallback=False`. Per-base P strings (`p_v_3`, ...), aggregate `n_p_nucleotides`, and pre-trim P mode all remain explicitly out of scope per [`docs/p_nucleotide_design.md`](p_nucleotide_design.md) §15.) | [`p_nucleotide_design.md`](p_nucleotide_design.md) (audit + v1 implementation shipped; pre-trim / per-base strings / aggregate / legacy auto-lift deferred), [`junction_n_addition_audit.md`](junction_n_addition_audit.md) (junction-biology umbrella) | [`test_p_nucleotide_contract.py`](../tests/test_p_nucleotide_contract.py) (21 post-slice pins: PEnd enum / PRegionAdded variant / PLength choice address / PAdditionPass struct / cartridge plane / manifest block / four AIRR fields / PLengthMismatch validator surface / legacy-orphan boundary / IR-correct D_5 lowering position), [`test_p_nucleotide_implementation.py`](../tests/test_p_nucleotide_implementation.py) (10 behaviour tests: byte-identical baseline, VJ V_3/J_5 emission, VDJ all-four-ends emission, `P_NUC` flag emission proof, replay round-trip, productive-only triad, validator clean recompute, manifest exposure, legacy-orphan invariant, plan-signature fold + cross-replay rejection), [`test_p_nucleotide_release.py`](../tests/test_p_nucleotide_release.py) (release-tier: productive IGH full stack with all four ends + D inversion + receptor revision + SHM + corruption + `validate_records` + cache parity; replay round-trip with `invert_d(prob=1.0)`; deterministic palindrome derivation on a synthetic short allele) | `engine_rs/src/passes/p_addition.rs` (`PAdditionPass` + `effective_coding_bytes` + `derive_p_bytes` + 6 unit tests covering palindrome math at each end + zero-length no-op + declared-choice-pattern + forward/reverse orientation), `engine_rs/src/address.rs` (`PEnd { V3, D5, D3, J5 }` enum + `ChoiceAddress::PLength { end }` variant + 4 canonical address spellings + parser + `frozen_address_spellings` test extension + `P_ADDITION_{V_3,D_5,D_3,J_5}` pass-name constants), `engine_rs/src/ir/sim_event.rs` (`PRegionAdded { end: PEnd, region: Region }` variant), `engine_rs/src/ir/builder.rs` (`record_p_region(end, region)` method that broadcasts the event without mutating `sim.sequence.regions`), `engine_rs/src/airr_record/record.rs` (four new `p_*_length: i64` fields), `engine_rs/src/airr_record/builder.rs` (event-ledger walk filtered to the four `p_addition.*` pass names, summing `region.len()` per `PRegionAdded.end`), `engine_rs/src/airr_record/validate.rs` (`PLengthMismatch` variant + per-end recompute), `engine_rs/src/python/plan.rs` (`push_p_addition(end, length_pairs)` PyO3 method + `parse_p_end` helper), `engine_rs/src/python/outcome.rs` (four new `dict.set_item` calls + `PLengthMismatch` serialisation with `details.source: "event-ledger:PRegionAdded"`), `engine_rs/src/airr_record/tests/projection.rs` (2 Rust unit tests: tampered-field detection + clean-projection acceptance), `engine_rs/src/pass/support.rs` (`PassCompileFact::PLengthSupport { end, support }` variant for downstream schedule-analyser introspection), `engine_rs/src/live_call/{dirty_signal_observer,walker_observer,refresh_plan}.rs` + `engine_rs/src/ir/event_log_observer.rs` (no-op match arms documenting that `PRegionAdded` doesn't invalidate live-call state); Python: [`src/GenAIRR/reference_models.py`](../src/GenAIRR/reference_models.py) (`p_nucleotide_lengths: Dict[str, EmpiricalDistributionSpec]` field + VJ-rejects-D-end validation), [`src/GenAIRR/_dataconfig_extract.py`](../src/GenAIRR/_dataconfig_extract.py) (`_p_nucleotide_lengths_from_models` resolver + `_explicit_models` widening), [`src/GenAIRR/_pipeline_ir.py`](../src/GenAIRR/_pipeline_ir.py) (`_RecombineStep.p_{v_3,d_5,d_3,j_5}_lengths` immutable fields), [`src/GenAIRR/_compile.py`](../src/GenAIRR/_compile.py) (conditional `plan.push_p_addition(...)` insertions at the four audited boundaries with `invert_d → p_addition.D_5 → assemble.D` ordering enforced), [`src/GenAIRR/experiment.py`](../src/GenAIRR/experiment.py) (threading through `recombine()`), [`src/GenAIRR/dataconfig/data_config.py`](../src/GenAIRR/dataconfig/data_config.py) (`models.p_nucleotide_models` manifest block: `length_keys` / `legacy_p_nucleotide_length_probs_present` / `legacy_fallback=False` / `supported_ends` / `in_plan_signature=True` / `in_content_hash=False`) |
| **Reference cartridge authoring** (`ReferenceCartridgeBuilder` v1 facade + `CartridgeBuildReport` audit trail. Fluent staged builder: `from_fasta(v_fasta=..., j_fasta=..., d_fasta=..., chain_type=...)` constructor → `infer_identity(species=..., locus=..., reference_set=..., name=..., source=...)` → `infer_v_subregions()` → `with_rules(...)` → `with_models(...)` → `build()` returning a plain `DataConfig` ready for `Experiment.on(cfg)`. v1 ships structural authoring only — FASTA parsing via the existing `parse_fasta` helper; V-subregion derivation via the same `compute_v_region_boundaries` helper used by `dataconfig_to_refdata`; rules / empirical-models attachment validated at attach time. `build()` stamps the canonical `schema_sha256`, attaches the typed `CartridgeBuildReport`, runs `verify_integrity()`. Statistical estimators (`estimate_allele_usage`, `estimate_trim_distributions`, `estimate_np_length_distributions`, `estimate_np_base_model`, `estimate_p_nucleotide_lengths`, `estimate_shm_rates`) are explicitly deferred and raise `AttributeError`. The 3 historical dead-reference sites (`DataConfig.build_report` docstring + `.private/scripts/build_imgt_configs.py` + `docs/build/` mirror) cleaned up in lockstep — docstring rewritten to reference the new builder, private script raises `NotImplementedError` at module load with explicit porting hint to `ReferenceCartridgeBuilder`, build-cache mirror remains as compile artefact only and is verified absent from import path. Output is a plain `DataConfig` — no parallel cartridge object, flows through existing `dataconfig_to_refdata` bridge unchanged. Internal `_SafeAnchorMixin` + `_BuilderVAllele` / `_BuilderJAllele` subclasses degrade gracefully when the native `_native._anchor` C resolver isn't built so v1 succeeds in any environment.) | [`reference_cartridge_authoring_audit.md`](reference_cartridge_authoring_audit.md) (v1 shipped; estimators deferred — §11.2 lists explicit method names + recommended ordering), [`reference_cartridge.md`](reference_cartridge.md) (3-paths comparison: bundled vs manual vs builder + when-not-to-use guide) | [`test_reference_cartridge_authoring_contract.py`](../tests/test_reference_cartridge_authoring_contract.py) (25 post-slice pins covering: live inference-adjacent helpers reused — `parse_fasta` + `compute_v_region_boundaries` + `cartridge_manifest` JSON-cleanliness + `verify_integrity` + `RefDataConfig.{vj,vdj}` + `ReferenceEmpiricalModels` + `ReferenceRulesSpec` + `dataconfig_to_refdata` bridge + bundled-loader's 106 configs; dead-reference cleanup state — `RandomDataConfigBuilder` docstring removed + private build script raises explicit `NotImplementedError` + `GenAIRR.dataconfig.make.*` namespaces stay absent; v1 surface — module landed at `cartridge_builder` + class + dataclass at top level + `from_fasta` on builder only + `infer_*` step methods on builder; v1-deferred surfaces — no `from_airr` / no `estimate_*` at any of the three potential owner locations), [`test_reference_cartridge_authoring_implementation.py`](../tests/test_reference_cartridge_authoring_implementation.py) (15 behaviour tests: top-level imports + `from_fasta().build()` round-trip + compile through Experiment + `infer_identity()` manifest population + `infer_v_subregions()` 5-label coverage + missing-gapped-V warning-not-crash + `with_rules()`/`with_models()` attach + `build_report` pickle round-trip + `.to_dict()` JSON-cleanliness + manual `DataConfig(...)` still supported + dead-reference integration probe + v1 estimator absence + 3 error-path bonus: VJ-rejects-d_fasta, VDJ-requires-d_fasta, duplicate-allele-rejection), [`test_reference_cartridge_authoring_release.py`](../tests/test_reference_cartridge_authoring_release.py) (release-tier composition: tiny inline-FASTA VDJ cartridge built end-to-end + manifest JSON-clean + build-report JSON-clean + compile-and-recombine via `allow_curatable_refdata()` + report contains `from_fasta`/`infer_identity`/`build` stages) | `src/GenAIRR/cartridge_builder.py` (~570 lines: `CartridgeBuildReport` dataclass + `to_dict` + `ReferenceCartridgeBuilder` class + 6 stages + `from_fasta` + `_open_fasta` path/string/file-object dispatcher + `_parse_allele_name` IMGT-pipe-aware header parser + `_parse_segment_fasta` per-segment ingestion with structured rejection + `_SafeAnchorMixin` for graceful native-resolver fallback + `_BuilderVAllele` / `_BuilderJAllele` / `_BuilderDAllele` subclasses + `_resolve_chain_type` / `_resolve_species` / `_locus_label_from_chain_type` helpers), `src/GenAIRR/__init__.py` (top-level exports: `ReferenceCartridgeBuilder`, `CartridgeBuildReport`), `src/GenAIRR/dataconfig/data_config.py` (cleaned docstring on `build_report` field references new builder; defaults stay `None` for bundled + manual paths), `.private/scripts/build_imgt_configs.py` (explicit `raise NotImplementedError` at module load with porting hint); engine bridge: existing `dataconfig_to_refdata` validates builder-produced cartridges unchanged. **Zero Rust changes**; pure Python authoring layer on existing live infrastructure (FASTA parser, V-subregion derivation, typed-plane validators, manifest builder, integrity checker). |
| **Allele usage estimation** (`ReferenceCartridgeBuilder.estimate_allele_usage(rearrangements, *, min_count=1.0, ambiguous="fractional", replace=True)` reads AIRR `v_call` / `d_call` / `j_call` from `list[dict]` / path / open text handle and writes per-segment weights into the typed `ReferenceEmpiricalModels.allele_usage: Optional[AlleleUsageSpec]` plane — `AlleleUsageSpec` carries `v` / `d` / `j` `Dict[str, float]` fields, validates non-empty allele names + finite positive weights + VJ-D-rejection, and lowers through `_dataconfig_extract._allele_usage_from_models` into `Experiment.recombine` with precedence kwarg > cartridge plane > uniform (no Rust changes; reuses the existing `_RecombineStep.weights_*` → `plan.push_sample_allele(weights=...)` → `AllelePoolDist::from_weights` engine surface). Three ambiguity policies: `"fractional"` (splits 1.0 / k credit across the k tie-set entries), `"truth_first"` (uses the first call only — mirrors `_mcp_summary.first_call`), `"reject"` (drops ambiguous rows into `report.rejected` with `reason="ambiguous_rejected"`). Unknown alleles, VJ-D-warning (one-time per stage), VDJ-missing-D (`reason="missing_d_call_on_vdj"`), and below-`min_count` drops all surface in the canonical `{stage, inputs, inferred, warnings}` stage entry on `CartridgeBuildReport`. Per-segment weights normalised to sum to 1.0. Idempotent via `previously_estimated` + `replace` kwarg (mirrors `infer_v_subregions` discipline). Manifest `models.allele_usage` block surfaces `available` / `segments` / `nonempty_segments` / `legacy_gene_use_dict_present` / `legacy_fallback=False` / `in_plan_signature=False` (inherited soft gap 1 — same boundary as the per-experiment `v_allele_weights` kwarg) / `source`. **`gene_use_dict` stays orphan**: the estimator does NOT auto-lift the legacy dict — same boundary the Markov / P-nucleotide slices respected for their legacy fields.) | [`allele_usage_estimation_design.md`](allele_usage_estimation_design.md) (slice shipped; plan-signature soft-gap tightening + `gene_use_dict` auto-lift remain deferred), [`reference_cartridge_authoring_audit.md`](reference_cartridge_authoring_audit.md) §11.2 (estimator ordering — 1 of 6 shipped) | [`test_allele_usage_estimation_contract.py`](../tests/test_allele_usage_estimation_contract.py) (25 post-slice pins: scaffold + present-state surfaces + soft-gap pin held + `gene_use_dict` orphan-preservation pin + the 5 absence pins flipped to present), [`test_allele_usage_estimation_implementation.py`](../tests/test_allele_usage_estimation_implementation.py) (15 behaviour tests: `AlleleUsageSpec` validation + explicit typed plane biases recombination + explicit recombine kwarg overrides plane + `gene_use_dict` orphan-no-auto-lift + fractional splits tie-set credit + truth_first first-call-only + reject drops to `report.rejected` + unknown alleles in rejected + VJ-ignores-D + VDJ-missing-D + `min_count` drops low-support + stage entry shape + manifest block + pickle round-trip + contract-pin-flip integration probe), [`test_allele_usage_estimation_release.py`](../tests/test_allele_usage_estimation_release.py) (release-tier: tiny inline-FASTA VDJ cartridge built via the v1 facade + `estimate_allele_usage` on a synthetic AIRR record set + recombination shows the estimated bias + `manifest["models"]["allele_usage"]["available"]=True` + report carries the `estimate_allele_usage` stage + explicit `recombine(v_allele_weights=...)` override still wins) | **Zero Rust changes**; pure Python authoring layer on the existing engine surface. Reuses `engine_rs/src/dist/allele_pool.rs::AllelePoolDist::from_weights` + `engine_rs/src/python/plan.rs::push_sample_allele(weights=...)` (already wired end-to-end pre-slice — confirmed by `pin_scaffold_*` group in the contract file); Python: [`src/GenAIRR/reference_models.py`](../src/GenAIRR/reference_models.py) (`AlleleUsageSpec` dataclass + chain-type-aware validate + `ReferenceEmpiricalModels.allele_usage` field), [`src/GenAIRR/_dataconfig_extract.py`](../src/GenAIRR/_dataconfig_extract.py) (`_allele_usage_from_models` resolver + `_explicit_models` widening), [`src/GenAIRR/experiment.py`](../src/GenAIRR/experiment.py) (precedence shim before `_resolve_allele_weights`), [`src/GenAIRR/cartridge_builder.py`](../src/GenAIRR/cartridge_builder.py) (`estimate_allele_usage` method + `_split_tie_set` + `_allele_names_for_segment` + `_load_rearrangements` helpers), [`src/GenAIRR/dataconfig/data_config.py`](../src/GenAIRR/dataconfig/data_config.py) (`_allele_usage_manifest_block` + manifest insertion), [`src/GenAIRR/__init__.py`](../src/GenAIRR/__init__.py) (top-level `AlleleUsageSpec` export). |
| **Allele-call ambiguity & disambiguation** | [`allele_call_audit.md`](allele_call_audit.md)                                     | [`test_allele_call_provenance.py`](../tests/test_allele_call_provenance.py) (21 tests) | `engine_rs/src/live_call/tests/` (7 modules: walker, caller, state, bitset, reference_index, …)          |
| **Junction-call fields**                   | [`junction_call_audit.md`](junction_call_audit.md)                                 | [`test_junction_call_provenance.py`](../tests/test_junction_call_provenance.py) (25 tests) | `engine_rs/src/junction.rs::tests` (12 unit tests), `engine_rs/src/airr_record/tests/anchors.rs`         |
| **Strict / replay semantics**              | [`productive_failure_mode_audit.md`](productive_failure_mode_audit.md) §5–§6.2 + [`engine_architecture.md`](engine_architecture.md) §1.3 | [`test_productive_failure_modes.py::test_strict_docstring_*`](../tests/test_productive_failure_modes.py), [`test_productive_failure_modes.py::test_replay_of_permissive_*`](../tests/test_productive_failure_modes.py) | `engine_rs/src/replay.rs`                                                                                |
| **Event ledger consistency**               | [`engine_architecture.md`](engine_architecture.md) §1.4–§1.7 + [`indel_provenance_audit.md`](indel_provenance_audit.md) §1.3 | Implicit in every provenance test (event-count assertions in `test_indel_provenance.py`, `test_allele_call_provenance.py`)  | `engine_rs/src/event.rs::tests` (compile-effect / event-emission consistency policy)                     |
| **Trace replay reproducibility**           | Each audit's "Replay round-trip" section                                            | Parametrised replay round-trip tests in every audit's test file (`test_indel_provenance.py::test_indel_replay_*`, `test_allele_call_provenance.py::test_replay_*`, `test_junction_call_provenance.py::test_replay_*`) | `engine_rs/src/trace_file.rs::tests` (schema, signature, content hash)                                   |
| **Distribution invariants (statistical)**  | [`distribution_invariant_audit.md`](distribution_invariant_audit.md)               | [`test_distribution_invariants.py`](../tests/test_distribution_invariants.py) (7 tests, MC ±5σ) | `engine_rs/src/passes/corrupt/indel/tests/constraints.rs::indel_pass_count_1_..._matches_exact_enumeration`, `..._matches_exact_continuation_weighting` (gold-standard oracle tests) |
| **Performance budgets (regression guard)** | [`performance_baseline.md`](performance_baseline.md)                                | [`test_performance_budgets.py`](../tests/test_performance_budgets.py) (8 budget tests, marker `@pytest.mark.performance`) | None — perf is bounded from Python only today                                                            |
| **Golden trace + AIRR-record compat**      | (no dedicated audit; format stability)                                              | [`test_trace_file_compat.py`](../tests/test_trace_file_compat.py), [`test_trace_file.py`](../tests/test_trace_file.py), fixtures under [`tests/golden/`](../tests/golden/) | `engine_rs/src/trace_file.rs::tests` (v1/v2 schema, content hash validation)                            |

### 1.1 What's NOT in the matrix

The matrix lists **enforced** guarantees. Each audit's §6 catalogues
drift items — known gaps that are documented but not yet pinned. The
canonical drift catalogue lives in:

- `indel_provenance_audit.md` §6
- `allele_call_audit.md` §6
- `junction_call_audit.md` §6
- `productive_failure_mode_audit.md` §6
- `distribution_invariant_audit.md` §6
- `performance_baseline.md` §6

Items appearing in §6 are **API ergonomics / denormalisation
conveniences**, not correctness drift — none of them indicate the
engine is misbehaving. They're catalogued so a future product
decision (or a downstream consumer's request) can scope them quickly.

---

## 2. CI guidance — three tiers

The validation suite splits cleanly into three tiers with different
runtime / coverage tradeoffs. Pick the tier appropriate for your
context.

### 2.1 Fast tier — pre-commit / per-push

**Goal:** catch correctness regressions in ≤ 60 s on a developer
machine.

```bash
# Python tests, skip wall-time budgets
pytest tests/ -m "not performance" -q

# Rust kernel tests (always fast)
cd engine_rs && cargo test --lib
```

What's exercised: every correctness, provenance, replay, and
distribution test from §1, except `test_performance_budgets.py`.
~830 Python tests + ~900 Rust tests.

Use this for `pre-commit`, every push to a feature branch, and the
default `make test` workflow.

### 2.2 Full tier — PR ready-for-review / nightly

**Goal:** include the regression-guard performance budgets so a 10×
slowdown blocks merge.

```bash
pytest tests/ -q                                  # all marks including performance
cd engine_rs && cargo test --all-features
```

Adds: [`test_performance_budgets.py`](../tests/test_performance_budgets.py)
(8 budget tests, ~2 s on dev machine, allow ~10 s on slow CI).

Use this in the PR pipeline (`.github/workflows/test.yml`) so a
budget regression fails CI and prompts a profile-first investigation
per [`performance_baseline.md`](performance_baseline.md) §5.

### 2.3 Release tier — pre-tag

**Goal:** prove the wheel + trace + AIRR-record formats are
byte-compatible with prior versions before publishing.

```bash
# Everything in the full tier, plus:
pytest tests/test_trace_file_compat.py -q   # golden trace round-trip
pytest tests/golden/ -q                      # golden AIRR-record fixtures (if present as a separate suite)
cd engine_rs && cargo test --all-features
maturin build --release                      # build the publishable wheel
```

The golden trace + AIRR-record compat tests pin schema and content
hash stability across versions. If they fail, you're shipping a
breaking change — bump the schema version and update the v1/v2
fixtures rather than silently rebasing the golden bytes.

---

## 3. Adding a new guarantee

When a future slice adds a new audit:

1. Land the audit doc under `docs/<topic>_audit.md` following the
   pattern of the existing audits (§ inventory → derivation rules →
   §6 drift).
2. Land the test file under `tests/test_<topic>_provenance.py` (or
   `test_<topic>_invariants.py`).
3. Add a row to §1 of this matrix.
4. If the slice introduces a new performance hot path, add a
   workload to [`tests/test_performance_budgets.py`](../tests/test_performance_budgets.py)
   and update [`performance_baseline.md`](performance_baseline.md) §3.
5. If the slice introduces a new on-disk format (trace, AIRR field),
   add a golden fixture under [`tests/golden/`](../tests/golden/).

---

## 4. How to read an audit doc

Every GenAIRR audit doc follows the same six-section structure:

1. **What's there today** — mechanism inventory + file:line
   citations into the engine.
2. **AIRR fields affected** / behaviour catalogue — which user-visible
   surfaces this slice touches.
3. **Derivation rules** — the *exact* per-case behaviour.
4. **Behaviour under evidence-changing mechanisms** — how the
   guarantee composes with SHM, indels, PCR errors, etc.
5. **Replay determinism** — confirmation that the guarantee
   reproduces under trace replay.
6. **Drift identified** — known gaps catalogued for follow-up; each
   item has a "possible fix" and a "tradeoff".

Read §1 + §3 for the rules; §6 for what's still open.

---

## 5. Quick reference — `make` targets

The repo's Makefile exposes the tiers as shortcuts:

```
make validate-fast      # Tier 1: correctness only
make validate-full      # Tier 2: correctness + performance budgets
make validate-release   # Tier 3: full + wheel build + golden compat
```

See the Makefile for the exact commands each target runs.
