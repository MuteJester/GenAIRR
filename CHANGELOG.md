# Changelog

All notable changes to GenAIRR are documented here.

## Unreleased

### Added

- Experimental `CoupledDuplexRepair` and `DuplexJoin` scripts for successive
  V-D/D-J joins on one shared two-strand D substrate, including inversion,
  canonical source ancestry, a global nucleotide ledger and cumulative budgets.
- `CodingEndOutcomeModel` with joint D outcomes, compatible-reference identities,
  generating-history estimation, held-out scoring and finite compatible-history
  likelihood fitting. Calibration records its source and sampling assumptions.
- Experimental `GenAIRR.experimental.DuplexRepair` for authored two-strand
  coding-end opening, resection, N extension, compatible annealing, paired-primer
  template synthesis and ligation, with both strand ancestries and exact replay.
- Trace schema 4 records productivity policy and active contracts. Replay and
  seed reruns share reference-content and address-version validation.

- Opt-in `Experiment.process_coding_ends()` with explicit original coding-end
  boundary confirmation and finite signed net outcomes. Positive values remove
  reference bases; negative values retain original-tip P extensions. The mode
  replaces automatic trim/P settings, preserves per-base origin provenance,
  and supports productive completion constraints. It does not claim calibrated
  microscopic hairpin-opening probabilities.

### Fixed

- The DataConfig bridge no longer assigns IMGT FWR/CDR subregions to V alleles
  whose `gapped_seq` declares no IMGT gaps (e.g. ungapped FASTA). Explicit
  subregions and IMGT-gapped sequences are unchanged; all 106 bundled cartridges
  keep identical native content hashes. Older user DataConfigs built from
  ungapped V FASTA now compile without fabricated subregions, so their native
  content hash changes. `compute_v_region_boundaries()` raises `ValueError` for
  ungapped input.
- `ReferenceCartridgeBuilder.infer_identity()` fills per-allele source, species
  and locus only where the loader left them unset, preserving FASTA-header
  provenance, and reports preserved species labels. New `with_build_date()`
  fixes `metadata.last_updated` for reproducible checksums with or without
  `infer_identity()`.
- The native reference bridge rejects V/D/J sequences longer than 65,535 nt
  with `ValueError` before pool mutation, instead of narrowing lengths to `u16`.
- Reference authoring uses the shared native landmark resolver, retains FASTA
  functionality/source declarations, and keeps missing species unknown.
  Optional checksum-bound evidence transports independent coding phase and
  coding-end confidence without changing legacy reference hashes. A reproducible
  audit marks bundled catalogues with findings as quarantined for default strict
  compilation; unresolved source evidence is explicit and archived data is not
  silently corrected.
- Replace approximate constrained-indel continuation weights with versioned
  natural sequential history sampling and exact weighted fallback. Report work
  exhaustion separately from proven empty support; validate replay count, kind,
  base weights and checked coordinates. Preserve sampling diagnostics separately
  from biological events.
- Use complete productive continuation support for ordinary zero-P VDJ
  recombination as well as P-aware processing, excluding future D/junction stops
  before choosing an allele or length.

- Use versioned observed-span productivity with nullable evidence and stable
  diagnostics. Reference admission overrides no longer legitimize invalid
  landmarks or frames; ambiguous possible stops remain unknown. Default
  anchors require V Cys and J Phe/Trp classes with original provenance.
  Productive-only runs assert accepted final metadata in both modes, and the
  record validator derives an independent assessment. Lineage divisions
  inherit the same constraints instead of silently mutating outside them.
- Filter ordinary and genotype-aware receptor-revision candidates by active
  biological contracts before drawing, preserving their original weights.
  Permissive execution and replay now enforce the same replacement checks.

- Keep standard V/D/J calls and alignment references tied to the generating
  allele; retain post-mutation inference separately as `recalled_*_call`.
  Ground-truth NP boundaries no longer move when random junction bases match
  a neighboring trimmed allele.
- Represent end loss as missing reference coverage, preserving actual
  structural deletions inside the observed interval.

- Preserve mutation-stage molecule lengths in `mutation_rate`, so later read
  loss does not inflate the reported cumulative SHM event burden. Clonal
  descendants retain inherited mutation counts and their event provenance.
- Translate the physical nucleotide junction even when it is out of frame
  or nonproductive; complete codons remain in `junction_aa`.
- Report impossible explicit recombination trims as structured errors instead
  of panicking, including in permissive runs.
- Load every bundled S5F model without pandas or NumPy. The two migrated
  artifacts retain exactly the same mutation-model values.
- Describe D inversion and its configured probability in `Experiment.describe()`.
- Calculate productivity and translation from the retained V coding origin,
  including internal stops and structural frameshifts outside the junction.
- Preserve coding-oriented annotations when emitting reverse-strand observations.
- Preserve P-nucleotide order in alignments and include P bases in combined NP
  strings; align inverted D references in their assembled orientation.
- Keep germline alignments within the called allele, including when its length
  differs from the sampled allele.
- Preserve recombination trims and P-addition provenance through clonal forks
  and receptor revision. Revision selects a different eligible V allele.
- Keep realized loss counters when a read is removed completely.
- Support nonzero P additions with productive-only sampling by retaining
  choices that admit a complete stop-free, in-frame sequence. Fixed P
  lengths are preserved; impossible support fails explicitly.
- Refresh live-call coordinates after no-op and boundary-edit passes so later
  end loss cannot reuse stale D/J positions.
- Use the correct canonical trim-end budgets for cached inverted-D re-calls,
  keeping NP-extension bounds and allele ties consistent with fresh scoring.
- Validate returned record fields and record/outcome counts; prevent custom
  metadata from overwriting generated sequence annotations.

### Compatibility

- Standard `v_call`, `d_call`, and `j_call` now name single generating alleles.
  Consumers of the previous evidence-based calls should read `recalled_v_call`,
  `recalled_d_call`, and `recalled_j_call`; optional `truth_*_call` columns remain
  aliases of the standard calls. Alignments now consistently use the generating
  reference, including when the inferred allele differs in length.
- Donor-genotype re-calls consider only alleles carried across that donor's
  two chromosomes, while retaining their original reference-pool IDs.

- `mutation_rate` sums each mutation stage's realized events divided by that
  stage's input molecule length. It is an event burden, not the mismatch
  fraction in an observed read; lineage records retain net-difference semantics.
- Clonal descendant outcomes now include ancestral trace, revisions, and events
  once, so inherited mutation counters have complete provenance.
- The raw-file checksums of `hh_s5f_60` and `hh_s5f_opposite` changed when their
  pandas objects were replaced by portable dictionaries; kernel values are
  bit-for-bit unchanged.
- `to_tsv()` now defaults to AIRR coordinates, CIGAR placement offsets and
  `T`/`F` boolean encoding. Pass `airr_strict=False` for the previous zero-based
  TSV representation. Python record dictionaries and default DataFrames keep
  zero-based half-open coordinates and segment-local CIGAR operation runs.
- AIRR `np1`/`np2` strings and lengths now include P bases between resolved
  segment alignments; they must not be interpreted as pure TdT N additions.
- Reverse-strand `*_sequence_start/end`, junction and amino-acid annotations
  now use the coding orientation, as AIRR requires. Read-window coordinates
  continue to index the emitted observation.
- Productive-only sampling excludes retained coding stops and protects the
  complete observed coding sequence. Seeded productive results can therefore
  differ from previous versions.
- Productive P-addition sampling restricts each draw's original weights to
  candidates with a productive completion across deterministic P bytes,
  N-addition support, and D orientation. Impossible support raises an error
  even in permissive mode instead of substituting zero P or sentinel N bases.
- NP length/base estimators reject records with explicitly nonzero P-length
  provenance. AIRR NP fields combine N and P; fitting them as N-only models
  while separately fitting P lengths would count P additions twice.

## [2.3.1] - 2026-07-05

Maintenance release. No public API or simulation-behaviour changes.

### Added
- `tools/build_imgt_configs.py` — a first-class maintainer tool that builds
  structural cartridges from IMGT V-QUEST germline FASTA via
  `ReferenceCartridgeBuilder`.
- `__all__` on the documented public authoring modules (`reference_models`,
  `reference_rules`, `cartridge_builder`, `genotype_priors`, `cohort`).

### Changed
- Large behaviour-preserving internal reorganization: several god-files split
  into focused packages/modules (`experiment.py`, `_cartridge_estimators`,
  `_compiled`, `result`, `genotype`, `utilities/visualize`, and the Rust
  `airr_record/validate.rs`), plus dead-code and duplicate-module removal. No
  user-visible import paths changed.

### Documentation
- Clarify that only the human IGH/IGK/IGL and TCRB cartridges carry
  data-derived distributions; every other bundled cartridge uses
  uniform/placeholder parameters (with guidance on fitting real ones).

## [2.3.0] - 2026-06-18

### Added - Per-individual diploid genotypes

A new `Genotype` layer makes V(D)J recombination haplotype-phased: attach a
genotype to an `Experiment` and the V, D and J of each rearrangement are drawn
from a single chromosome, honouring allele presence/absence, zygosity, gene
deletion, copy number, and novel/private alleles.

- **`Genotype` builder** - `from_dataconfig` / `permissive`, `homozygous` /
  `heterozygous` / `delete_gene` / `duplicate_gene`, `add_novel_allele` (private
  alleles injected into an effective reference), `complete_from_reference`,
  `chromosome_weights`, `with_subject`, `to_table` / `to_tsv`.
- **`Experiment.with_genotype(g)`** - phased recombination, with per-record
  `subject_id` / `haplotype` provenance and `result.genotypes`.
- **`Genotype.sample(...)`** - draw a genotype from population priors
  (independent per-gene Hardy-Weinberg: allele frequencies + per-haplotype
  deletion).
- **Cartridge genotype-priors plane** - persist a donor-population germline
  prior on a cartridge (`DataConfig.genotype_priors` / `PopulationGenotypeModel`);
  author or estimate it (`PopulationGenotypeModel.from_genotypes`, builder
  `set_genotype_priors` / `estimate_genotype_priors`); `Genotype.sample(cfg)`
  consumes it with per-input provenance.
- **`Experiment.run_cohort(...)`** - simulate many subjects (each with its own
  genotype) in one call, returning a `CohortResult` with per-subject access and
  combined export.
- **Genotype-aware receptor revision** -
  `receptor_revision(..., same_haplotype=True)` restricts the replacement V to
  the carried alleles on the drawn rearrangement chromosome.

### Documentation

- New genotype guide cluster: overview, sampling & population priors, cohorts,
  and benchmarking genotype inference (with a worked TIgGER / IgDiscover example
  recovering a planted genotype).

### Compatibility

- Purely additive. An experiment with no genotype attached produces
  byte-identical output to 2.2.0 (pinned by a checksum test).

## [2.0.0] — 2026-05-25

### Breaking Changes — Experiment DSL refactor for biological readability

The Experiment DSL was reworked around a single principle: a fluent chain
should read like a wet-lab protocol, not a configuration file. The pre-v2
DSL leaked engine mechanics (`corrupt_*` as an umbrella verb, `respect=`
as a runtime kwarg, `count=` as the unit for somatic hypermutation) and
hid biological assumptions inside `recombine(trim=True)`. v2.0 fixes both
ends. **Old scripts will not run unmodified** — see
[MIGRATING.md](MIGRATING.md) for the full map.

#### New methods

- **`.productive_only()`** — declares the productive-sequence contract
  bundle as a chainable method instead of a `respect=ga.productive()`
  kwarg. Order-independent in the chain.
- **`.restrict_alleles(v=..., d=..., j=...)`** — was `.using(...)`.
  Renamed for semantic precision (narrows sampling support, doesn't
  pin one allele).
- **`.trim(enabled=False)` / `.trim(v_3=..., d_5=..., d_3=..., j_5=...)`** —
  exonuclease-trim override method. Default-on inside `.recombine()`;
  call `.trim()` only to disable or supply custom distributions.
  Chain-position-guarded (raises if called after a mutation/corruption
  step). Was previously hidden behind `.recombine(trim=...)`.
- **`mutate(rate=0.03)`** — per-base SHM rate; engine samples count
  per record from `Poisson(rate × pool_len)`. The canonical
  biology default. `mutate(count=...)` remains as the explicit
  benchmark-friendly form; passing both raises.
- **`exp.describe()`** — renders the chain as a biology-style narrative
  per-step. Use it as a quick sanity check before `compile()`.

#### Renamed methods (full map)

| v1 | v2 |
|---|---|
| `respect=ga.productive()` (kwarg) | `.productive_only()` (method) |
| `.using(v=...)` | `.restrict_alleles(v=...)` |
| `.recombine(trim=True/False)` | `.recombine()` + `.trim(enabled=...)` |
| `mutate(count=15)` | `mutate(rate=0.03)` *(canonical)* or `mutate(count=15)` *(retained)* |
| `corrupt_pcr` | `pcr_amplify` |
| `corrupt_quality` | `sequencing_errors` |
| `corrupt_contaminants` | `contaminate` |
| `corrupt_indels` | `polymerase_indels` |
| `corrupt_ns` | `ambiguous_base_calls` |
| `corrupt_reverse_complement` | `random_strand_orientation` |
| `corrupt_5prime_loss` | `primer_trim_5prime` |
| `corrupt_3prime_loss` | `primer_trim_3prime` |
| `with_clonal_structure(n_clones=N, size=K)` | `expand_clones(n_clones=N, per_clone=K)` |

Both `pcr_amplify` and `sequencing_errors` accept `rate=` (per-base
error probability, runtime Poisson) in addition to the legacy
`count=`. Exactly one must be provided. `mutate(rate=...)` shares
the same shape.

#### Stricter validation

- `recombine(np2_lengths=...)` on a VJ chain is now a `ValueError`
  (was silent ignore).
- Binding a raw `RefDataConfig` (no `DataConfig` backing) emits a
  `UserWarning` when synthetic NP-length defaults or no-op trim are
  injected. Pass explicit `np1_lengths=` / `np2_lengths=` and call
  `.trim(enabled=False)` to silence.

#### Engine

- New `MutationCountSource` enum in `engine_rs/src/passes/mutate/`
  with `Distribution` and `Rate(f64)` variants. The rate variant
  samples count via a Knuth Poisson sampler at execute time
  against the current pool length. `S5FMutationPass::new_rate(...)`
  and `UniformMutationPass::new_rate(...)` constructors expose the
  new mode; the existing count-based constructors are unchanged.

#### Tests + benchmarks

- 631 pytest, 637 cargo lib, 128 cargo integration/property — all green.
- Throughput unchanged on the count path (no engine-side cost when
  rate-mode is unused). Rate-mode adds ~1 Poisson sample per record.

---

## [1.0.0] — 2025-03-10

### Breaking Changes

This is a complete rewrite. The Python simulation engine, pipeline system, and
public API have all been replaced. Code written for 0.6.x will not work without
changes.

- **New Experiment DSL** — Single entry point replaces the old
  `AugmentationPipeline` + step-based system. All simulation is now configured
  through `Experiment.on("config").mutate(...).observe(...).run(n=1000)`.
- **C simulation engine** — All rearrangement, mutation, corruption, and AIRR
  serialization now run in compiled C. Typical throughput is 50–100k
  sequences/second on a single core.
- **Removed packages** — `simulation/`, `sequence/`, `mutation/`, `pipeline/`,
  `steps/`, `validation/`, `container/`, `unbiased/`, `TCR/` have all been
  removed. Their functionality is now in the C backend.
- **No runtime dependencies** — The core package has zero mandatory
  dependencies. Optional extras (`numpy`, `scipy`, `graphviz`, `fastmcp`) are
  available via `pip install GenAIRR[all]`.

### Added

- 106+ built-in DataConfigs covering 23 species (human, mouse, rat, rabbit,
  dog, cat, cow, sheep, pig, horse, gorilla, rhesus, cynomolgus, ferret,
  chicken, salmon, trout, zebrafish, platypus, alpaca, dromedary, goat).
- GDC binary format for fast DataConfig serialization to/from the C engine.
- Selection pressure simulation (CDR/FWR replacement acceptance rates).
- Class switch recombination (CSR) with isotype-specific rates.
- D-gene inversion and receptor revision simulation.
- Pipeline introspection hooks — snapshot the internal sequence state at any
  pipeline stage for debugging and analysis.
- Execution tracing — human-readable log of every internal decision.
- MCP server with 21 tools for AI-assisted sequence analysis.
- Cross-platform support: Linux, macOS, and Windows (pre-built wheels).
- Allele locking — constrain simulation to specific V/D/J alleles.
- Streaming interface via `CompiledSimulator.stream()`.

### Fixed

- S5F silent mutations (substitution could draw the same base).
- TCRB anchorless alleles (anchor=0 treated as valid position).
- Indel mutation string truncation (null bytes in annotation strings).
- Junction bounds with 3' corruption (junction_end could exceed sequence length).

## [0.6.3] — 2024-12-01

Last release of the Python-only engine. See the
[0.6.x documentation](https://genairr.readthedocs.io/) for details.
