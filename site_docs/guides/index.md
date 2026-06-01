# Simulation Guides

<p class="lead">Task-focused guides for designing simulations.
Start with the Experiment builder; the rest are deep dives into
specific biological mechanisms or pipeline patterns.</p>

## Start here

- **[The Experiment builder](experiment-builder.md)** — the
  control panel. Every pipeline goes through this builder; this
  is the one guide every user should read.

## Biology mechanisms

- **[Recombination + junction biology](recombination-junction.md)**
  — V(D)J join, trim-and-fill, productivity contract.
- **[Junction N/P additions](junction-additions.md)** — N-base
  composition models, P-nucleotide lengths, layout diagrams.
- **[D inversion + receptor revision](recombination-editing.md)**
  — the two recombination-stage editing mechanisms.
- **[SHM and mutation targeting](shm-targeting.md)** — uniform
  vs S5F, per-segment and per-V-subregion rates, counter
  partitions.
- **[Clonal families](clonal-families.md)** — `expand_clones`,
  the ancestor / descendant phase split, family validation.

## Library + sequencing

- **[Corruption + sequencing artefacts](corruption-sequencing.md)**
  — PCR errors, sequencing errors, indels, end-loss, N
  corruption, strand orientation.
- **[Paired-end reads and FASTQ](paired-end-fastq.md)** — R1 /
  R2 projection + quality models.

## Reproducibility

- **[Trace, replay, reproducibility](trace-replay.md)** — what
  a trace contains, how to replay it, every failure mode.

## Cartridge authoring (separate section)

For building or estimating reference cartridges, see the
dedicated **[Reference Cartridge Authoring](../cartridges/index.md)**
section.

For higher-level routing by reader intent, see
**[Choose your path](../learn.md)**.
