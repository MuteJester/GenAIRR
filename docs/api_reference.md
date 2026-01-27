# API Reference

Complete reference for GenAIRR's public API — classes, functions, and their signatures.

---

## Top-Level Imports

```python
from GenAIRR import (
    # Core
    Pipeline,               # AugmentationPipeline alias
    steps,                  # Namespace for all augmentation steps
    simulate,               # One-liner convenience function
    SimulationContainer,    # Result object

    # Mutation models
    S5F,                    # Context-dependent mutation
    Uniform,                # Position-independent mutation

    # Data configurations
    HUMAN_IGH_OGRDB,        # Human heavy chain (OGRDB)
    HUMAN_IGH_EXTENDED,     # Human heavy chain (extended set)
    HUMAN_IGK_OGRDB,        # Human kappa light chain (OGRDB)
    HUMAN_IGL_OGRDB,        # Human lambda light chain (OGRDB)
    HUMAN_TCRB_IMGT,        # Human TCR beta (IMGT)

    # Reproducibility
    set_seed, get_seed, reset_seed,
)
```

---

## `simulate()`

```python
simulate(config, mutation_model, productive=True, n=1)
```

One-liner convenience function. Creates a minimal pipeline (`SimulateSequence` + `FixVPositionAfterTrimmingIndexAmbiguity`) and executes it.

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `config` | `DataConfig` | — | Germline data configuration |
| `mutation_model` | `MutationModel` | — | S5F or Uniform instance |
| `productive` | `bool` | `True` | Only generate in-frame sequences without stop codons |
| `n` | `int` | `1` | Number of sequences to generate |

**Returns:** `SimulationContainer` if `n=1`, `List[SimulationContainer]` if `n > 1`.

```python
from GenAIRR import simulate, HUMAN_IGH_OGRDB, S5F

result = simulate(HUMAN_IGH_OGRDB, S5F(0.01, 0.05))
results = simulate(HUMAN_IGH_OGRDB, S5F(0.01, 0.05), n=100)
```

---

## `Pipeline`

```python
Pipeline(steps, config=None)
```

Alias for `AugmentationPipeline`. Executes a list of steps in order, passing a shared `SimulationContainer` through each.

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `steps` | `List[AugmentationStep]` | — | Ordered list of pipeline steps |
| `config` | `DataConfig` | `None` | Germline data configuration. Required for most steps. |

### Methods

#### `execute()`

```python
pipeline.execute() -> SimulationContainer
```

Creates an empty `SimulationContainer`, runs each step's `.apply()` method on it, and returns the populated container.

#### `plot(filename)`

```python
pipeline.plot(filename) -> None
```

Generates a GraphViz diagram of the pipeline structure and saves it as an image.

---

## Mutation Models

### `S5F`

```python
S5F(min_mutation_rate=0.0, max_mutation_rate=0.0, custom_model=None, productive=False)
```

Context-dependent somatic hypermutation model. Mutation probability at each position depends on the surrounding 5-mer sequence motif.

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `min_mutation_rate` | `float` | `0.0` | Lower bound of mutation rate range |
| `max_mutation_rate` | `float` | `0.0` | Upper bound of mutation rate range |
| `custom_model` | `str` | `None` | Path to a custom S5F model file |
| `productive` | `bool` | `False` | Avoid stop codons and preserve CDR3 anchor residues |

### `Uniform`

```python
Uniform(min_mutation_rate=0.0, max_mutation_rate=0.0, productive=False)
```

Position-independent mutation model. Each base has equal probability of mutation regardless of sequence context.

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `min_mutation_rate` | `float` | `0.0` | Lower bound of mutation rate range |
| `max_mutation_rate` | `float` | `0.0` | Upper bound of mutation rate range |
| `productive` | `bool` | `False` | Avoid stop codons and preserve CDR3 anchor residues |

---

## Steps

All steps inherit from `AugmentationStep` and implement `.apply(container: SimulationContainer) -> None`.

### `SimulateSequence`

```python
steps.SimulateSequence(mutation_model, productive=False, specific_v=None, specific_d=None, specific_j=None)
```

Performs V(D)J recombination and applies mutations. Always the first step in a pipeline.

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `mutation_model` | `MutationModel` | — | S5F or Uniform instance |
| `productive` | `bool` | `False` | Only generate in-frame, stop-codon-free sequences |
| `specific_v` | `VAllele` | `None` | Force a specific V allele |
| `specific_d` | `DAllele` | `None` | Force a specific D allele (heavy/TRB only) |
| `specific_j` | `JAllele` | `None` | Force a specific J allele |

### `FixVPositionAfterTrimmingIndexAmbiguity`

```python
steps.FixVPositionAfterTrimmingIndexAmbiguity()
```

Resolves V segment boundary ambiguity caused by trimming. No parameters.

### `FixDPositionAfterTrimmingIndexAmbiguity`

```python
steps.FixDPositionAfterTrimmingIndexAmbiguity()
```

Resolves D segment boundary ambiguity. Use for heavy chain and TRB only. No parameters.

### `FixJPositionAfterTrimmingIndexAmbiguity`

```python
steps.FixJPositionAfterTrimmingIndexAmbiguity()
```

Resolves J segment boundary ambiguity. No parameters.

### `CorrectForVEndCut`

```python
steps.CorrectForVEndCut()
```

Adjusts V-end position when 3' trimming extends into the coding region. No parameters.

### `CorrectForDTrims`

```python
steps.CorrectForDTrims()
```

Adjusts D segment boundaries after 5'/3' trimming. Heavy chain and TRB only. No parameters.

### `DistillMutationRate`

```python
steps.DistillMutationRate()
```

Computes the mutation rate by comparing the mutated sequence to the germline and stores it in the container. **Place before artifact steps** to get a clean biological mutation rate.

### `CorruptSequenceBeginning`

```python
steps.CorruptSequenceBeginning(
    *,
    probability=0.7,
    event_weights=(0.4, 0.4, 0.2),
    nucleotide_add_coefficient=210,
    nucleotide_remove_coefficient=310,
    nucleotide_add_after_remove_coefficient=50,
    random_sequence_add_probability=1.0
)
```

Simulates 5' end degradation. All parameters are keyword-only.

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `probability` | `float` | `0.7` | Probability of applying corruption |
| `event_weights` | `tuple` | `(0.4, 0.4, 0.2)` | Weights for (add, remove, add-after-remove) events |
| `nucleotide_add_coefficient` | `int` | `210` | Controls distribution of nucleotides added |
| `nucleotide_remove_coefficient` | `int` | `310` | Controls distribution of nucleotides removed |
| `nucleotide_add_after_remove_coefficient` | `int` | `50` | Controls nucleotides added after removal |
| `random_sequence_add_probability` | `float` | `1.0` | Probability of adding random (vs germline-matching) bases |

### `EnforceSequenceLength`

```python
steps.EnforceSequenceLength(*, max_length=576)
```

Truncates sequences exceeding `max_length` nucleotides from the 3' end, simulating sequencing platform read-length limits.

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `max_length` | `int` | `576` | Maximum allowed sequence length |

### `InsertNs`

```python
steps.InsertNs(*, n_ratio=0.02, probability=0.5)
```

Replaces random positions with 'N' to simulate ambiguous base calls.

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `n_ratio` | `float` | `0.02` | Fraction of positions to replace |
| `probability` | `float` | `0.5` | Probability of applying N-insertion to a given sequence |

### `InsertIndels`

```python
steps.InsertIndels(*, probability=0.5, max_indels=5, insertion_probability=0.5, deletion_probability=0.5)
```

Introduces random insertions and deletions.

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `probability` | `float` | `0.5` | Probability of applying indels |
| `max_indels` | `int` | `5` | Maximum number of indel events |
| `insertion_probability` | `float` | `0.5` | Relative weight for insertion events |
| `deletion_probability` | `float` | `0.5` | Relative weight for deletion events |

### `ShortDValidation`

```python
steps.ShortDValidation(short_d_length=5)
```

Validates D segment length and adjusts annotations for short D regions.

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `short_d_length` | `int` | `5` | Threshold below which D is considered "short" |

### `FilterTCRDJAmbiguities`

```python
steps.FilterTCRDJAmbiguities()
```

Filters ambiguous D-J assignments specific to TCR beta chains. No parameters.

---

## `SimulationContainer`

The mutable data object returned by `pipeline.execute()`. Holds the simulated sequence, annotations, and metadata.

### Data Access

```python
result = pipeline.execute()

# Dictionary export
data = result.get_dict()

# Direct attribute access
print(result.sequence)
print(result.v_call)

# Dictionary-style access
print(result['mutation_rate'])
result['mutation_rate'] = 0.05
```

### Methods

| Method | Signature | Description |
|--------|-----------|-------------|
| `get_dict()` | `-> dict` | Export all fields as a plain dictionary |
| `add_mutation(position, change)` | `-> None` | Log a mutation (e.g., `add_mutation(150, "A>T")`) |
| `shift_positions(offset)` | `-> None` | Shift all position annotations by `offset` |
| `update_from_dict(d)` | `-> None` | Bulk-update attributes from a dictionary |

### Output Fields

| Field | Type | Description |
|-------|------|-------------|
| `sequence` | `str` | Final nucleotide sequence |
| `v_call` | `str` | V allele name(s) |
| `d_call` | `str` | D allele name(s) |
| `j_call` | `str` | J allele name(s) |
| `v_sequence_start` / `end` | `int` | V segment position in final sequence |
| `d_sequence_start` / `end` | `int` | D segment position |
| `j_sequence_start` / `end` | `int` | J segment position |
| `v_germline_start` / `end` | `int` | Position within original V gene |
| `d_germline_start` / `end` | `int` | Position within original D gene |
| `j_germline_start` / `end` | `int` | Position within original J gene |
| `junction_start` / `end` | `int` | CDR3 region boundaries |
| `np1_length` | `int` | NP1 region length |
| `np2_length` | `int` | NP2 region length |
| `mutations` | `dict` | `{position: "X>Y"}` mutation log |
| `mutation_rate` | `float` | Fraction of mutated positions |
| `productive` | `bool` | In-frame, no stop codons |
| `vj_in_frame` | `bool` | V and J maintain reading frame |
| `stop_codon` | `bool` | Premature stop codon present |
| `indels` | `dict` | Indel information |

---

## Seed Management

```python
from GenAIRR import set_seed, get_seed, reset_seed
```

| Function | Signature | Description |
|----------|-----------|-------------|
| `set_seed(seed)` | `seed: int` | Fix the global random state for reproducible output |
| `get_seed()` | `-> int` | Retrieve the current seed value |
| `reset_seed()` | `-> None` | Clear the seed, restoring non-deterministic behavior |

---

## Data Configurations

| Import | Chain | Species | Has D segment |
|--------|-------|---------|---------------|
| `HUMAN_IGH_OGRDB` | Heavy (IGH) | Human | Yes |
| `HUMAN_IGH_EXTENDED` | Heavy (IGH) | Human | Yes |
| `HUMAN_IGK_OGRDB` | Kappa light (IGK) | Human | No |
| `HUMAN_IGL_OGRDB` | Lambda light (IGL) | Human | No |
| `HUMAN_TCRB_IMGT` | TCR Beta (TRB) | Human | Yes |

### Accessing Alleles

```python
config = HUMAN_IGH_OGRDB

# List all gene families
list(config.v_alleles.keys())  # ['IGHVF1-G1', 'IGHVF1-G2', ...]
list(config.d_alleles.keys())  # ['IGHD1-1', 'IGHD1-7', ...]
list(config.j_alleles.keys())  # ['IGHJ1', 'IGHJ2', ...]

# Get alleles for a gene
alleles = config.v_alleles['IGHVF1-G1']  # List[VAllele]
allele = alleles[0]
print(allele.name)      # 'IGHVF1-G1*01'
print(allele.family)    # 'IGHVF1'
print(allele.gene)      # 'IGHVF1-G1'
print(allele.length)    # 301
print(allele.sequence)  # nucleotide sequence string
```

---

## Complete Pipeline Examples

### Heavy Chain

```python
from GenAIRR import Pipeline, steps, HUMAN_IGH_OGRDB, S5F

pipeline = Pipeline(
    config=HUMAN_IGH_OGRDB,
    steps=[
        steps.SimulateSequence(S5F(min_mutation_rate=0.01, max_mutation_rate=0.05), productive=True),
        steps.FixVPositionAfterTrimmingIndexAmbiguity(),
        steps.FixDPositionAfterTrimmingIndexAmbiguity(),
        steps.FixJPositionAfterTrimmingIndexAmbiguity(),
        steps.CorrectForVEndCut(),
        steps.CorrectForDTrims(),
        steps.DistillMutationRate(),
        steps.CorruptSequenceBeginning(),
        steps.EnforceSequenceLength(),
        steps.InsertNs(),
        steps.ShortDValidation(),
        steps.InsertIndels(),
    ]
)
```

### Kappa Light Chain

```python
from GenAIRR import Pipeline, steps, HUMAN_IGK_OGRDB, S5F

pipeline = Pipeline(
    config=HUMAN_IGK_OGRDB,
    steps=[
        steps.SimulateSequence(S5F(min_mutation_rate=0.01, max_mutation_rate=0.05), productive=True),
        steps.FixVPositionAfterTrimmingIndexAmbiguity(),
        steps.FixJPositionAfterTrimmingIndexAmbiguity(),
        steps.CorrectForVEndCut(),
        steps.DistillMutationRate(),
        steps.CorruptSequenceBeginning(),
        steps.EnforceSequenceLength(),
        steps.InsertNs(),
        steps.InsertIndels(),
    ]
)
```

### TCR Beta

```python
from GenAIRR import Pipeline, steps, HUMAN_TCRB_IMGT
from GenAIRR.mutation import Uniform

pipeline = Pipeline(
    config=HUMAN_TCRB_IMGT,
    steps=[
        steps.SimulateSequence(Uniform(min_mutation_rate=0.0, max_mutation_rate=0.01), productive=True),
        steps.FixVPositionAfterTrimmingIndexAmbiguity(),
        steps.FixDPositionAfterTrimmingIndexAmbiguity(),
        steps.FixJPositionAfterTrimmingIndexAmbiguity(),
        steps.CorrectForVEndCut(),
        steps.CorrectForDTrims(),
        steps.FilterTCRDJAmbiguities(),
        steps.DistillMutationRate(),
        steps.CorruptSequenceBeginning(),
        steps.EnforceSequenceLength(),
        steps.InsertNs(),
        steps.ShortDValidation(),
        steps.InsertIndels(),
    ]
)
```

### Custom Step

```python
from GenAIRR.steps.StepBase import AugmentationStep
from GenAIRR.container.SimulationContainer import SimulationContainer

class MyStep(AugmentationStep):
    def __init__(self, *, threshold=0.5):
        super().__init__()
        self.threshold = threshold

    def apply(self, container: SimulationContainer) -> None:
        # self.data_config is set automatically by the Pipeline
        # Modify container in place
        pass
```
