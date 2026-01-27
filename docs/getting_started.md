# Installation & Quick Start

This guide covers installing GenAIRR and running your first simulation. By the end, you will have generated a realistic immunoglobulin heavy chain sequence with somatic hypermutation.

## Installation

Install the latest stable release from PyPI:

```bash
pip install GenAIRR
```

GenAIRR requires **Python 3.9+** and depends on NumPy, SciPy, and pandas.

!!! tip "Virtual environment recommended"
    Install GenAIRR inside a virtual environment (`venv` or `conda`) to avoid dependency conflicts with other packages.

To verify the installation:

```python
import GenAIRR
print(GenAIRR.__version__)
```

---

## Core Concepts in 60 Seconds

GenAIRR has four core building blocks:

| Component | Purpose |
|-----------|---------|
| **DataConfig** | Holds germline allele sets, trimming distributions, and empirical data for a species/chain type |
| **Pipeline** | Executes an ordered list of steps against a config to produce a simulated sequence |
| **Steps** | Individual transformations — sequence generation, correction, artifact injection |
| **SimulationContainer** | The output object carrying the sequence, annotations, and metadata |

Their relationship:

```
DataConfig ──► Pipeline(config, steps) ──► SimulationContainer
                    │
                    ├── SimulateSequence
                    ├── Fix...Ambiguity
                    ├── CorrectFor...
                    └── InsertIndels, InsertNs, ...
```

---

## Your First Simulation

### One-liner with `simulate()`

The fastest way to generate a sequence:

```python
from GenAIRR import simulate, HUMAN_IGH_OGRDB, S5F

result = simulate(HUMAN_IGH_OGRDB, S5F(0.01, 0.05))
print(result.sequence)
```

`simulate()` creates a minimal pipeline internally — it runs `SimulateSequence` and `FixVPositionAfterTrimmingIndexAmbiguity`, then returns a `SimulationContainer`.

!!! note "When to use `simulate()` vs `Pipeline`"
    Use `simulate()` for quick exploratory work. Switch to an explicit `Pipeline` when you need artifact simulation, custom step ordering, or measurement steps like `DistillMutationRate`.

**Parameters:**

- `config` — a `DataConfig` instance (e.g., `HUMAN_IGH_OGRDB` for human heavy chain)
- `mutation_model` — an `S5F` or `Uniform` instance specifying the mutation rate range
- `productive` — if `True` (default), only generates in-frame sequences without stop codons
- `n` — number of sequences to generate (default: 1; returns a list when `n > 1`)

### Generating multiple sequences

```python
results = simulate(HUMAN_IGH_OGRDB, S5F(0.01, 0.05), n=100)
print(f"Generated {len(results)} sequences")
```

---

## Building an Explicit Pipeline

For full control, create a `Pipeline` with your choice of steps:

```python
from GenAIRR import Pipeline, steps, HUMAN_IGH_OGRDB, S5F

pipeline = Pipeline(
    config=HUMAN_IGH_OGRDB,
    steps=[
        steps.SimulateSequence(
            S5F(min_mutation_rate=0.01, max_mutation_rate=0.05),
            productive=True
        ),
        steps.FixVPositionAfterTrimmingIndexAmbiguity(),
        steps.FixDPositionAfterTrimmingIndexAmbiguity(),
        steps.FixJPositionAfterTrimmingIndexAmbiguity(),
        steps.CorrectForVEndCut(),
        steps.CorrectForDTrims(),
        steps.DistillMutationRate(),
    ]
)

result = pipeline.execute()
data = result.get_dict()
```

Each call to `pipeline.execute()` produces one `SimulationContainer`. The container holds all annotations:

```python
print(data['sequence'][:60])     # nucleotide sequence
print(data['v_call'])            # V allele used
print(data['d_call'])            # D allele used
print(data['j_call'])            # J allele used
print(data['mutation_rate'])     # fraction of mutated positions
print(data['mutations'])         # dict of {position: "X>Y"} changes
print(data['productive'])        # True if in-frame, no stop codons
```

---

## Available Data Configurations

GenAIRR ships with pre-built configs derived from the OGRDB and IMGT germline databases:

| Import name | Chain type | Species | Segments |
|-------------|-----------|---------|----------|
| `HUMAN_IGH_OGRDB` | Heavy (IGH) | Human | V, D, J |
| `HUMAN_IGH_EXTENDED` | Heavy (IGH) | Human | V, D, J (extended set) |
| `HUMAN_IGK_OGRDB` | Kappa light (IGK) | Human | V, J |
| `HUMAN_IGL_OGRDB` | Lambda light (IGL) | Human | V, J |
| `HUMAN_TCRB_IMGT` | TCR Beta (TRB) | Human | V, D, J |

Import them directly:

```python
from GenAIRR import HUMAN_IGH_OGRDB, HUMAN_IGK_OGRDB, HUMAN_IGL_OGRDB, HUMAN_TCRB_IMGT
```

---

## Available Mutation Models

### S5F — Context-Dependent Mutation

The S5F model captures the context-dependent substitution patterns observed in real somatic hypermutation. Mutation probability at each position depends on the surrounding 5-mer motif.

```python
from GenAIRR import S5F

model = S5F(min_mutation_rate=0.01, max_mutation_rate=0.05)
```

### Uniform — Position-Independent Mutation

Each position has equal probability of mutation, regardless of sequence context. Useful for null-model comparisons.

```python
from GenAIRR.mutation import Uniform

model = Uniform(min_mutation_rate=0.01, max_mutation_rate=0.05)
```

**Choosing mutation rates** — the `min_mutation_rate` and `max_mutation_rate` define a range; each simulated sequence samples a rate uniformly from this range. Typical ranges:

| Cell type | Mutation rate range |
|-----------|-------------------|
| Naive B cells | 0.001 – 0.01 |
| Memory B cells | 0.02 – 0.08 |
| Plasma cells | 0.05 – 0.25 |

---

## Exporting Results

### To pandas DataFrame

```python
import pandas as pd

sequences = [pipeline.execute().get_dict() for _ in range(100)]
df = pd.DataFrame(sequences)
df.to_csv('simulated_sequences.csv', index=False)
```

### To FASTA

```python
with open('sequences.fasta', 'w') as f:
    for i, seq in enumerate(sequences):
        f.write(f">seq_{i:04d}\n{seq['sequence']}\n")
```

---

## Reproducibility

GenAIRR provides seed management for deterministic output:

```python
from GenAIRR import set_seed, get_seed, reset_seed

set_seed(42)
result_a = pipeline.execute()

set_seed(42)
result_b = pipeline.execute()

assert result_a.sequence == result_b.sequence  # identical
```

- `set_seed(n)` — fix the global random state
- `get_seed()` — retrieve the current seed value
- `reset_seed()` — clear the seed, restoring non-deterministic behavior

!!! warning "Remember to reset"
    If you set a seed for a reproducibility check, call `reset_seed()` afterwards to restore non-deterministic behavior for production runs.

---

## Next Steps

- **[Step-by-Step Tutorial](step_by_step_tutorial.md)** — Detailed walkthrough of building a full pipeline with explanations for each step
- **[How the Pipeline Works](genairr_flow.md)** — Architecture deep-dive into DataConfig, Steps, and SimulationContainer
- **[Biological Context](biological_context.md)** — The immunobiology behind GenAIRR's simulation model
- **[Parameter Reference](parameter_reference.md)** — Complete parameter documentation for every step
