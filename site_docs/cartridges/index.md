# Reference Cartridge Authoring

<p class="lead">Build your own reference cartridge: custom
alleles, custom rules, custom empirical biology. Three focused
guides cover the end-to-end loop from FASTA to a build report
you can pin in CI.</p>

## The three guides

1. **[Build a reference cartridge](../guides/build-reference-cartridge.md)**
   the practical builder workflow from FASTA to `build()`.
   Start here for any custom cartridge.
2. **[Estimate models from data](../guides/estimate-cartridge-models.md)**
   turn an AIRR-like rearrangement table into empirical models
   on the cartridge (allele usage, trim, NP length, NP base,
   P-nucleotide length).
3. **[Inspect manifest + build report](../guides/cartridge-manifest-report.md)**
   audit the cartridge's current state and how it was produced;
   the canonical CI-gate surface.

## Background

The conceptual model behind the builder lives at
**[Reference cartridge](../concepts/reference-cartridge.md)**
(four typed planes: identity, catalogue, rules, empirical
models). Read it before starting non-trivial builder work.

For the API-level catalogue of the builder class and the spec
dataclasses, see the API Reference:
**[`ReferenceCartridgeBuilder`](../reference/reference-cartridge-builder.md)**
and **[Reference models and rules](../reference/reference-models.md)**.

For higher-level routing by reader intent, see
**[Choose your path](../learn.md)**.
