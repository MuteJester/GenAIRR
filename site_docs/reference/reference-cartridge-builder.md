# ReferenceCartridgeBuilder

<p class="lead">A fluent builder for custom reference cartridges
from raw FASTA. Every call returns the same builder extended by
one stage; at <code>.build()</code> the engine validates the
typed planes and stamps a <code>schema_sha256</code> on the
resulting <code>DataConfig</code>. Every stage's inputs,
warnings, and rejection details are captured in
<code>builder.report()</code> - the cartridge ships with an
auditable trail.</p>

For the broader workflow, see
[Build a reference cartridge](../guides/build-reference-cartridge.md);
for the estimator deep dive, see
[Estimate cartridge models from real data](../guides/estimate-cartridge-models.md).

## Common methods

| Method | Purpose |
|---|---|
| `.from_fasta(...)` | Class-level constructor; parses V/D/J FASTA into the catalogue |
| `.infer_identity(...)` | Resolve species / chain / locus / reference-set / name |
| `.infer_v_subregions()` | Annotate V alleles with FWR/CDR intervals (required for `v_subregion_rates`) |
| `.with_rules(...)` | Attach an explicit `ReferenceRulesSpec` |
| `.with_models(...)` | Attach an explicit `ReferenceEmpiricalModels` bundle |
| `.estimate_allele_usage(records)` | Estimate per-segment allele weights |
| `.estimate_trim_distributions(records)` | Estimate per-key trim distributions |
| `.estimate_np_length_distributions(records)` | Estimate NP1 / NP2 length distributions |
| `.estimate_np_base_model(records, kind=...)` | Estimate NP-base model (Markov or empirical first-base) |
| `.estimate_p_nucleotide_lengths(records)` | Estimate palindromic-insertion length distributions |
| `.build()` | Finalise, stamp checksum, return the `DataConfig` |
| `.report()` | Get the `CartridgeBuildReport` with stages / warnings / rejections |

## Class reference

::: GenAIRR.cartridge_builder.ReferenceCartridgeBuilder
    options:
      show_source: false
      show_root_heading: true
      heading_level: 3
      members:
        - from_fasta
        - infer_identity
        - infer_v_subregions
        - with_rules
        - with_models
        - estimate_allele_usage
        - estimate_trim_distributions
        - estimate_np_length_distributions
        - estimate_np_base_model
        - estimate_p_nucleotide_lengths
        - build
        - report
