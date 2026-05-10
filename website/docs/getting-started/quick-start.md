---
title: Quick Start
sidebar_label: Quick Start
---

# Quick Start

This page gets you from install to your first simulated sequences in under two minutes.

## Install

```bash
pip install GenAIRR
```

Pre-built wheels are available for **Linux**, **macOS**, and **Windows** (Python 3.9+). No compiler required. 

Installing `GenAIRR` pulls both the pure-Python wrapper and the high-performance Rust simulation kernel `genairr_engine`.

Verify the installation:

```python
import GenAIRR
print(GenAIRR.__version__)
```

## Your first simulation

Generate 1,000 human heavy-chain sequences using standard V(D)J recombination:

```python
import GenAIRR as ga

# 1. Start an experiment on a specific species/chain config
# 2. Add a recombination step
# 3. Run for 1,000 iterations with a fixed seed
result = ga.Experiment.on("human_igh").recombine().run(n=1000, seed=42)
```

`result` is a `SimulationResult` containing 1,000 outcomes. Each outcome produces an AIRR-format record with **~70 ground-truth fields**:

```python
rec = result[0]
rec["sequence"]        # full nucleotide sequence
rec["v_call"]          # e.g. "IGHV3-23*01"
rec["d_call"]          # e.g. "IGHD3-10*01"
rec["j_call"]          # e.g. "IGHJ4*02"
rec["productive"]      # True / False
rec["sequence_length"] # e.g. 365
```

## Add somatic hypermutation

Chain a `.mutate()` step to introduce context-aware SHM using the built-in S5F model:

```python
import GenAIRR as ga

result = (
    ga.Experiment.on("human_igh")
    .recombine()
    .mutate(model="s5f", count=(5, 25))
    .run(n=1000, seed=42)
)

rec = result[0]
rec["n_mutations"]    # e.g. 14
rec["mutation_rate"]  # e.g. 0.038
```

The `count=(5, 25)` argument tells the engine to draw a uniform integer number of mutations between 5 and 25 for every sequence.

## Model sequencing realism

Simulate real-world artifacts like primer trimming, sequencing errors, and ambiguous bases:

```python
import GenAIRR as ga

result = (
    ga.Experiment.on("human_igh")
    .recombine()

    # Somatic hypermutation (5 to 25 mutations)
    .mutate(count=(5, 25))

    # Sequencing artifacts: 5' and 3' end loss
    .corrupt_5prime_loss(length=(5, 30))
    .corrupt_3prime_loss(length=(5, 20))

    # Post-sequencing noise: Indels and N-bases
    .corrupt_indels(count=(0, 2), insertion_prob=0.5)
    .corrupt_ns(count=5)

    .run(n=1000, seed=42)
)
```

Each method appends a "pass" to the simulation pipeline. You only include the stages you need.

## Export results

```python
# Save as AIRR-compliant TSV (standard)
result.to_csv("repertoire.tsv")

# Save as FASTA
result.to_fasta("repertoire.fasta")

# Convert to pandas DataFrame
df = result.to_dataframe()
```

## Switch species or chain

GenAIRR ships with **106 built-in configurations** covering 23 species:

```python
import GenAIRR as ga

# See all available configs
print(ga.list_configs())

# Use any config by name (species_chain)
ga.Experiment.on("mouse_igh").recombine().run(n=500)
ga.Experiment.on("rabbit_tcrb").recombine().run(n=500)
ga.Experiment.on("rhesus_igk").recombine().run(n=500)
```

See [Choosing a Config](/docs/getting-started/choosing-config) for the full list and naming conventions.

## Streaming (memory-efficient)

For large datasets (millions of sequences), use `.stream_records()` to process one record at a time without loading the entire batch into RAM:

```python
import GenAIRR as ga

exp = ga.Experiment.on("human_igh").recombine().mutate(count=10)

# Stream 1,000,000 records
for record in exp.stream_records(n=1_000_000, seed=42):
    process(record) # record is a dict
```

## Reproducibility

GenAIRR is bit-for-bit deterministic. The same seed produces identical results across Linux, macOS, and Windows:

```python
import GenAIRR as ga

exp = ga.Experiment.on("human_igh").recombine()

r1 = exp.run(n=100, seed=42)
r2 = exp.run(n=100, seed=42)

assert r1[0].final_simulation().bases() == r2[0].final_simulation().bases()
```

## Next steps

- [Choosing a Config](/docs/getting-started/choosing-config) — pick the right species and chain type
- [Understanding Output](/docs/getting-started/interpreting-results) — what each of the ~70 fields means
- [Guides](/docs/guides/) — recipes for clonal families, antigen selection, and advanced workflows
