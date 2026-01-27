# How the Pipeline Works

This page explains GenAIRR's architecture — how `DataConfig`, `Pipeline`, `Steps`, and `SimulationContainer` fit together, and the order of operations when a pipeline executes.

---

## Architecture Overview

```
                   ┌──────────────┐
                   │  DataConfig   │  Germline alleles, trimming distributions,
                   │  (immutable)  │  NP-region models, correction maps
                   └──────┬───────┘
                          │
                          ▼
              ┌───────────────────────┐
              │       Pipeline        │
              │  config + [steps...]  │
              └───────────┬───────────┘
                          │  .execute()
                          ▼
              ┌───────────────────────┐
              │  SimulationContainer  │ ◄── created empty
              └───────────┬───────────┘
                          │
          ┌───────────────┼───────────────┐
          ▼               ▼               ▼
    ┌──────────┐   ┌──────────┐   ┌──────────┐
    │  Step 1  │──►│  Step 2  │──►│  Step N  │
    └──────────┘   └──────────┘   └──────────┘
          │               │               │
          └───────────────┼───────────────┘
                          ▼
              ┌───────────────────────┐
              │  SimulationContainer  │ ◄── populated with results
              └───────────────────────┘
```

1. The **Pipeline** receives a `DataConfig` and a list of steps
2. On `.execute()`, it creates an empty `SimulationContainer`
3. Each step's `.apply(container)` method is called in order, mutating the container in place
4. The final container is returned to the caller

---

## DataConfig

A `DataConfig` encapsulates everything needed to simulate sequences for a specific species and chain type:

| Attribute | Contents |
|-----------|----------|
| `v_alleles`, `d_alleles`, `j_alleles` | Dictionaries mapping gene names to lists of `Allele` objects |
| `trim_dicts` | Trimming length distributions for each segment end (V_3, D_5, D_3, J_5) |
| `NP_transitions` | Markov chain transition matrices for NP-region nucleotide generation |
| `NP_first_bases` | Starting nucleotide distributions for NP regions |
| `NP_lengths` | Length distributions for NP1 and NP2 regions |
| `correction_maps` | Ambiguity resolution structures (trim similarity maps, N-ambiguity graphs) |

### Built-in Configs

```python
from GenAIRR import HUMAN_IGH_OGRDB, HUMAN_IGK_OGRDB, HUMAN_IGL_OGRDB, HUMAN_TCRB_IMGT
```

These are lazily loaded — the data is only read from disk when first accessed.

### Custom Configs

You can create a `DataConfig` from your own FASTA files using `CustomDataConfigGenerator`. See [Custom DataConfig](custom_data_config.md) for details.

---

## Pipeline

The `Pipeline` class (aliased from `AugmentationPipeline`) is the main entry point:

```python
from GenAIRR import Pipeline, steps, HUMAN_IGH_OGRDB, S5F

pipeline = Pipeline(
    config=HUMAN_IGH_OGRDB,
    steps=[
        steps.SimulateSequence(S5F(min_mutation_rate=0.01, max_mutation_rate=0.05), productive=True),
        steps.FixVPositionAfterTrimmingIndexAmbiguity(),
        steps.CorrectForVEndCut(),
    ]
)
```

**Key behavior:**

- `config` is passed to each step automatically — steps do not need to know about the config at construction time
- Steps execute in list order; earlier steps' modifications are visible to later steps
- Each `.execute()` call produces an independent `SimulationContainer`
- The pipeline is reusable — call `.execute()` as many times as needed

!!! tip "Pipeline reuse"
    Construct the pipeline once, then call `.execute()` in a loop to generate thousands of sequences. There is no per-call overhead beyond the simulation itself.

---

## Steps

Every step inherits from `AugmentationStep` and implements an `.apply(container)` method. Steps fall into four categories:

### 1. Sequence Generation

| Step | Purpose |
|------|---------|
| `SimulateSequence` | Performs V(D)J recombination, applies mutations, populates the container with the initial sequence and annotations |

This is always the first step.

### 2. Position Correction

| Step | Purpose |
|------|---------|
| `FixVPositionAfterTrimmingIndexAmbiguity` | Resolves V segment boundary ambiguity caused by trimming |
| `FixDPositionAfterTrimmingIndexAmbiguity` | Resolves D segment boundary ambiguity (heavy/TRB only) |
| `FixJPositionAfterTrimmingIndexAmbiguity` | Resolves J segment boundary ambiguity |
| `CorrectForVEndCut` | Adjusts V-end position when trimming extends into the coding region |
| `CorrectForDTrims` | Adjusts D boundaries after 5'/3' trimming (heavy/TRB only) |
| `ShortDValidation` | Validates short D segments |

These steps ensure that ground-truth annotations remain accurate after trimming introduces positional ambiguity.

### 3. Measurement

| Step | Purpose |
|------|---------|
| `DistillMutationRate` | Computes and stores the mutation rate. Must be placed **before** artifact steps. |

### 4. Artifact Simulation

| Step | Purpose |
|------|---------|
| `CorruptSequenceBeginning` | Simulates 5' end degradation (truncation, random base addition) |
| `EnforceSequenceLength` | Truncates sequences to a maximum read length |
| `InsertNs` | Replaces random positions with 'N' (ambiguous bases) |
| `InsertIndels` | Introduces insertion/deletion errors |

### Recommended Step Order

```python
steps=[
    # 1. Generation
    steps.SimulateSequence(...),
    # 2. Position corrections
    steps.FixVPositionAfterTrimmingIndexAmbiguity(),
    steps.FixDPositionAfterTrimmingIndexAmbiguity(),   # heavy/TRB only
    steps.FixJPositionAfterTrimmingIndexAmbiguity(),
    # 3. Biological corrections
    steps.CorrectForVEndCut(),
    steps.CorrectForDTrims(),                          # heavy/TRB only
    # 4. Measurement (before artifacts)
    steps.DistillMutationRate(),
    # 5. Artifacts
    steps.CorruptSequenceBeginning(),
    steps.EnforceSequenceLength(),
    steps.InsertNs(),
    # 6. Validation + errors
    steps.ShortDValidation(),
    steps.InsertIndels(),
]
```

---

## SimulationContainer

The `SimulationContainer` is the mutable data object passed through the pipeline. After execution, it contains all sequence data and annotations.

### Accessing Data

```python
result = pipeline.execute()

# As a dictionary
data = result.get_dict()
print(data['sequence'])
print(data['v_call'])

# Dictionary-style access
print(result['mutation_rate'])
result['mutation_rate'] = 0.05
```

### Key Fields

| Field | Type | Description |
|-------|------|-------------|
| `sequence` | `str` | Final nucleotide sequence |
| `v_call` | `str` | V allele name(s) used |
| `d_call` | `str` | D allele name(s) used |
| `j_call` | `str` | J allele name(s) used |
| `v_sequence_start` / `end` | `int` | V segment position in final sequence |
| `d_sequence_start` / `end` | `int` | D segment position in final sequence |
| `j_sequence_start` / `end` | `int` | J segment position in final sequence |
| `v_germline_start` / `end` | `int` | Position within original V gene |
| `d_germline_start` / `end` | `int` | Position within original D gene |
| `j_germline_start` / `end` | `int` | Position within original J gene |
| `junction_start` / `end` | `int` | CDR3 region boundaries |
| `mutations` | `dict` | `{position: "X>Y"}` mutation log |
| `mutation_rate` | `float` | Fraction of mutated positions |
| `productive` | `bool` | In-frame and no stop codons |
| `vj_in_frame` | `bool` | V and J maintain reading frame |
| `stop_codon` | `bool` | Premature stop codon present |
| `indels` | `dict` | Inserted indel information |

### Utility Methods

- `add_mutation(position, change)` — log a mutation
- `shift_positions(offset)` — shift all position annotations (used after insertions/deletions)
- `update_from_dict(d)` — bulk-update attributes from a dictionary
- `get_dict()` — export all fields as a plain dictionary

---

## Writing Custom Steps

Create a custom step by subclassing `AugmentationStep`:

```python
from GenAIRR.steps.StepBase import AugmentationStep
from GenAIRR.container.SimulationContainer import SimulationContainer

class MyCustomStep(AugmentationStep):
    def __init__(self, *, my_param=0.5):
        self.my_param = my_param

    def apply(self, container: SimulationContainer) -> None:
        # Access the config via self.data_config (set automatically by Pipeline)
        # Modify the container in place
        if some_condition:
            container.sequence = modified_sequence
```

Use it in a pipeline:

```python
pipeline = Pipeline(
    config=HUMAN_IGH_OGRDB,
    steps=[
        steps.SimulateSequence(S5F(min_mutation_rate=0.01, max_mutation_rate=0.05), productive=True),
        MyCustomStep(my_param=0.8),
    ]
)
```

The pipeline automatically sets `self.data_config` on each step before calling `.apply()`.

---

## Next Steps

- **[Step-by-Step Tutorial](step_by_step_tutorial.md)** — Build a pipeline from scratch
- **[Parameter Reference](parameter_reference.md)** — Detailed parameter docs for every step
- **[Custom DataConfig](custom_data_config.md)** — Create configs from your own germline data
