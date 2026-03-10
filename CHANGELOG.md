# Changelog

All notable changes to GenAIRR are documented here.

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
