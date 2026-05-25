# Changelog

All notable changes to GenAIRR are documented here.

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
