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

Pre-built wheels are available for Linux, macOS, and Windows (Python 3.9+). No compiler required.

Verify the installation:

```python
import GenAIRR
print(GenAIRR.__version__)  # 1.0.0
```

## Your first simulation

Generate 1,000 unmutated human heavy-chain sequences:

```python
from GenAIRR import Experiment

result = Experiment.on("human_igh").run(n=1000, seed=42)
```

That's it. `result` is a `SimulationResult` containing 1,000 AIRR-format records. Each record is a dictionary with 47 ground-truth fields:

```python
rec = result[0]
rec["sequence"]        # full nucleotide sequence
rec["v_call"]          # e.g. "IGHVF10-G50*04"
rec["d_call"]          # e.g. "IGHD2-21*02"
rec["j_call"]          # e.g. "IGHJ4*02"
rec["productive"]      # True / False
rec["sequence_length"] # e.g. 365
```

## Add somatic hypermutation

Chain a `.mutate()` phase to introduce SHM:

```python
from GenAIRR import Experiment
from GenAIRR.ops import rate

result = Experiment.on("human_igh").mutate(rate(0.02, 0.08)).run(n=1000, seed=42)

rec = result[0]
rec["mutation_rate"]  # e.g. 0.054
rec["n_mutations"]    # e.g. 19
rec["mutations"]      # e.g. "13:a>T,30:g>A,..."
```

The `rate(0.02, 0.08)` op tells GenAIRR to draw a mutation rate uniformly between 2% and 8% for each sequence.

## Add sequencing artifacts

Model what happens during a real sequencing experiment by adding more phases:

```python
from GenAIRR import Experiment
from GenAIRR.ops import (
    rate, model,
    with_5prime_loss, with_3prime_loss,
    with_indels, with_ns,
)

result = (
    Experiment.on("human_igh")

    # Somatic hypermutation
    .mutate(
        model("s5f"),
        rate(0.02, 0.08),
    )

    # Sequencing artifacts
    .sequence(
        with_5prime_loss(min_remove=5, max_remove=30),
        with_3prime_loss(min_remove=5, max_remove=20),
    )

    # Post-sequencing noise
    .observe(
        with_indels(prob=0.01),
        with_ns(prob=0.005),
    )

    .run(n=1000, seed=42)
)
```

Each phase is optional — include only what you need. Calling `.run()` on a bare `Experiment.on(...)` gives you clean, unmutated rearrangements.

## Export results

```python
# AIRR-compliant TSV
result.to_csv("repertoire.tsv")

# FASTA
result.to_fasta("repertoire.fasta")

# pandas DataFrame (requires pandas)
df = result.to_dataframe()
```

## Switch species or chain

GenAIRR ships with 106 built-in configurations covering 23 species:

```python
from GenAIRR import Experiment, list_configs

# See all available configs
print(list_configs())

# Use any config by name
Experiment.on("mouse_igh").run(n=500, seed=1)
Experiment.on("rabbit_tcrb").run(n=500, seed=1)
Experiment.on("rhesus_igk").run(n=500, seed=1)
```

Config names follow the pattern `species_chain` (lowercase). See [Choosing a Config](/docs/getting-started/choosing-config) for the full list and naming conventions.

## Streaming (memory-efficient)

For large datasets, use `.compile()` + `.stream()` to process one record at a time without accumulating them in memory:

```python
from GenAIRR import Experiment
from GenAIRR.ops import rate

sim = Experiment.on("human_igh").mutate(rate(0.05, 0.15)).compile(seed=42)

for record in sim.stream():
    print(record["v_call"])  # process one at a time
    break                    # infinite iterator — break when done
```

## Reproducibility

Pass a `seed` to get identical results across runs and platforms:

```python
r1 = Experiment.on("human_igh").run(n=100, seed=42)
r2 = Experiment.on("human_igh").run(n=100, seed=42)
assert r1[0]["sequence"] == r2[0]["sequence"]  # always True
```

## Next steps

- [Choosing a Config](/docs/getting-started/choosing-config) — pick the right species and chain type
- [Understanding Output](/docs/getting-started/interpreting-results) — what each of the 47 fields means
- [Guides](/docs/guides/) — recipes for SHM, artifacts, biological events, export, and more
