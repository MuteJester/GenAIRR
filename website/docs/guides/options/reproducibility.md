---
title: Reproducibility
sidebar_label: Reproducibility
---

# Reproducibility

Reproducible simulation is essential for scientific work — you need to generate the exact same dataset again when rerunning an experiment, sharing results with a collaborator, or debugging an annotation pipeline. GenAIRR supports fully deterministic simulation through seed-based random number generation.

## How seeds work

The C engine uses a seeded pseudo-random number generator (PRNG). When you pass a `seed` to `.run()` or `.compile()`, the engine initializes its PRNG with that seed before generating any sequences. Since every random decision in the pipeline (allele selection, trimming amounts, NP nucleotide choices, mutation positions, corruption events) flows through this single PRNG, the same seed always produces the same sequence of random numbers and therefore the same output.

```python
from GenAIRR import Experiment

r1 = Experiment.on("human_igh").run(n=1000, seed=42)
r2 = Experiment.on("human_igh").run(n=1000, seed=42)

assert r1[0]["sequence"] == r2[0]["sequence"]
assert r1[0]["v_call"] == r2[0]["v_call"]
assert r1[0]["mutation_rate"] == r2[0]["mutation_rate"]
```

Everything is identical — the sequences, the allele calls, the mutation positions, the junction regions, the trimming amounts. The same seed with the same configuration and the same GenAIRR version always produces bit-for-bit identical results.

## Without a seed

If you omit the seed, GenAIRR uses a time-based seed derived from the system clock, producing different results each run:

```python
r1 = Experiment.on("human_igh").run(n=100)
r2 = Experiment.on("human_igh").run(n=100)
# r1 and r2 will (almost certainly) differ
```

## Seeds with streaming

Seeds work with `.compile()` + `.stream()` too. The PRNG is seeded once when streaming begins, and each subsequent `.stream()` call advances the same PRNG state:

```python
sim = Experiment.on("human_igh").compile(seed=42)

records = []
for rec in sim.stream():
    records.append(rec)
    if len(records) >= 100:
        break
```

## Different seeds, different data

Different seeds produce independent datasets — useful for generating train/test splits or multiple replicates:

```python
train = Experiment.on("human_igh").run(n=5000, seed=1)
test  = Experiment.on("human_igh").run(n=1000, seed=2)
```

## What affects reproducibility

A result is reproducible when the seed **and** the full experiment configuration are identical. Changing any of these will produce different output even with the same seed:

- **Config** — different species/chain, or a different version of the DataConfig
- **Pipeline ops** — adding, removing, or reordering ops (e.g., adding `with_indels()`)
- **Op parameters** — changing `rate(0.02, 0.08)` to `rate(0.03, 0.06)`
- **GenAIRR version** — internal algorithm changes between versions may alter the PRNG call sequence
- **Number of sequences** — `run(n=100, seed=42)` and `run(n=200, seed=42)` share the first 100 sequences, but only if no retry logic is involved (productive mode retries consume additional PRNG calls)

The seed does **not** depend on the platform (Linux, macOS, Windows) or Python version — the C engine uses its own PRNG, independent of Python's `random` module.
