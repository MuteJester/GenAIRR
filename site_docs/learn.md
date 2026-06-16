# Choose your path

<p class="lead">Five paths through the docs, organised by what
you came here to do. Each path is a short curriculum — read in
order and you'll have a working mental model by the end. Skip
to any starting point that matches your current task.</p>

## I want to simulate sequences

The fastest path from `pip install` to a productive heavy-chain
batch.

1. **[Quick start](getting-started/quick-start.md)** — 5 minutes,
   100 productive heavy chains, AIRR records returned.
2. **[Your first AIRR record](getting-started/first-airr-record.md)**
   — field-by-field tour of what the engine emits.
3. **[The Experiment builder](guides/experiment-builder.md)** —
   the full DSL surface: recombination, mutation, constraints,
   clonal expansion, corruption, paired-end, compile reuse.
4. **[Export the results](getting-started/export-results.md)** —
   TSV / CSV / FASTA / FASTQ / DataFrame.

If you only ever read one page after the quick start, make it
the Experiment builder.

## I want to build a reference cartridge

Custom alleles, custom empirical distributions, custom biology.

1. **[Reference cartridge concept](concepts/reference-cartridge.md)**
   — the four-plane model (identity, catalogue, rules, empirical
   models). Read this before any builder work.
2. **[Build a reference cartridge](guides/build-reference-cartridge.md)**
   — the practical builder workflow from FASTA to `build()`.
3. **[Estimate models from data](guides/estimate-cartridge-models.md)**
   — turn an AIRR-like rearrangement table into empirical models
   on the cartridge.
4. **[Inspect manifest + build report](guides/cartridge-manifest-report.md)**
   — audit what's in a cartridge before you ship it or pin it
   in CI.

The bundled cartridges (`HUMAN_IGH_OGRDB`, etc.) work for most
projects; reach for the builder when you have a non-standard
catalogue or want estimated biology to match a specific dataset.

## I want reproducible / validated output

Bit-stable runs, replayable traces, output you can defend.

1. **[Validation hub](validation/index.md)** — the overall
   reproducibility + validation story.
2. **[`validate_records`](validation/validate-records.md)** — the
   per-record AIRR-output correctness gate.
3. **[Trace, replay, reproducibility](guides/trace-replay.md)** —
   what a trace contains, how to replay it, every failure mode,
   strict-mode behaviour.

For deep-publication reproducibility, the canonical pattern is:
seed for fast runs, trace for the records you'll defend.

## I want to understand biology mechanisms

The engine's biological surface — what's modelled, what's not,
where the v1 boundary sits.

1. **[Recombination + junction biology](guides/recombination-junction.md)**
   — V(D)J join, trim-and-fill, productivity contract.
2. **[Junction N/P additions](guides/junction-additions.md)** —
   N-base composition models, P-nucleotide lengths, layout
   diagrams.
3. **[D inversion + receptor revision](guides/recombination-editing.md)**
   — the two recombination-stage editing mechanisms.
4. **[SHM and mutation targeting](guides/shm-targeting.md)** —
   uniform vs S5F, per-segment and per-V-subregion rates,
   counter partitions.
5. **[Clonal simulation overview](guides/clonal-families.md)** —
   choose `clonal_lineage` for BCR trees, `clonal_repertoire` for
   TCR / abundance repertoires, or legacy `expand_clones`.
6. **[Clonal lineage trees](guides/clonal-lineage.md)** — BCR SHM
   trees, selection, final-cell sampling, lineage metadata, and
   tree exports.
7. **[Clonal repertoires](guides/clonal-repertoire.md)** — TCR and
   flat-BCR clone sizes, `duplicate_count`, and AIRR clone-caller
   export.
8. **[Corruption + sequencing artefacts](guides/corruption-sequencing.md)**
   — the observation-stage mechanisms (PCR, sequencing errors,
   indels, end-loss, N corruption, strand).

Each guide opens with the biology, names the v1 boundary
decisions, and links back to the audit doc that defined them.

## I'm contributing to GenAIRR

The contributor doorway — engine architecture, audit-first
workflow, mechanism additions.

1. **[Architecture (Contributor)](architecture/index.md)** — the
   engine mental model, the audit-first workflow, the
   "before you add a new mechanism" checklist, deep links into
   the audit corpus.
2. **[`docs/engine_architecture.md`](https://github.com/MuteJester/GenAIRR/blob/master/docs/engine_architecture.md)**
   — the seven engine invariants. Required reading before any
   kernel work.
3. **[`docs/adding_a_pass.md`](https://github.com/MuteJester/GenAIRR/blob/master/docs/adding_a_pass.md)**
   — the pass-author playbook with the minimal pass template
   and the three required test patterns.
4. **[`docs/validation_matrix.md`](https://github.com/MuteJester/GenAIRR/blob/master/docs/validation_matrix.md)**
   — the navigable map: every guarantee → audit doc → test file
   → kernel invariant.

GenAIRR's release process is **audit-first**: every mechanism
gets specified, pinned, and validated before implementation.
The architecture landing page explains how that discipline
shapes contribution.

---

The above paths cover ~95 % of real reader intents. If yours
isn't here, the [API Reference](reference/index.md) is the
authoritative public-surface catalogue.
