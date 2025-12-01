# GenAIRR API Quick Reference

## Essential Imports

```python
from GenAIRR.pipeline import AugmentationPipeline
from GenAIRR.steps import (
    SimulateSequence, FixVPositionAfterTrimmingIndexAmbiguity,
    FixDPositionAfterTrimmingIndexAmbiguity, FixJPositionAfterTrimmingIndexAmbiguity,
    CorrectForVEndCut, CorrectForDTrims, CorruptSequenceBeginning,
    InsertNs, InsertIndels, ShortDValidation, DistillMutationRate
)
from GenAIRR.mutation import S5F, Uniform
from GenAIRR.data import (
    HUMAN_IGH_OGRDB, HUMAN_IGH_EXTENDED,
    HUMAN_IGK_OGRDB, HUMAN_IGL_OGRDB, 
    HUMAN_TCRB_IMGT
)
from GenAIRR.steps.StepBase import AugmentationStep
```

## Data Configurations

| Config | Description | Use Case |
|--------|-------------|----------|
| `HUMAN_IGH_OGRDB` | Human heavy chain immunoglobulin (OGRDB) | BCR heavy chain simulation |
| `HUMAN_IGH_EXTENDED` | Extended human heavy chain immunoglobulin | Extended BCR heavy chain simulation |
| `HUMAN_IGK_OGRDB` | Human kappa light chain (OGRDB) | BCR kappa light chain simulation |
| `HUMAN_IGL_OGRDB` | Human lambda light chain (OGRDB) | BCR lambda light chain simulation |
| `HUMAN_TCRB_IMGT` | Human T-cell receptor beta (IMGT) | TCR-β simulation |

## Basic Pipeline Setup

```python
# Set data configuration
AugmentationStep.set_dataconfig(HUMAN_IGH_OGRDB)

# Create pipeline
pipeline = AugmentationPipeline([
    # Add steps here
])

# Execute simulation
result = pipeline.execute()
sequence_data = result.get_dict()
```

## Core Pipeline Steps

### 1. Sequence Simulation
```python
# S5F mutation model with mutation rate range
SimulateSequence(S5F(min_mutation_rate=0.003, max_mutation_rate=0.25), productive=True)

# Uniform mutation model
SimulateSequence(Uniform(min_mutation_rate=0.0, max_mutation_rate=0.1), productive=True)

# Naive sequence (no mutations)
SimulateSequence(Uniform(min_mutation_rate=0, max_mutation_rate=0), productive=True)

# With specific alleles
SimulateSequence(
    S5F(min_mutation_rate=0.02, max_mutation_rate=0.08),
    productive=True,
    specific_v=v_allele_object,
    specific_d=d_allele_object,
    specific_j=j_allele_object
)
```

**Parameters:**
- `mutation_model`: Instance of mutation model (S5F or Uniform)
- `productive`: Ensure sequence is productive (default: False)
- `specific_v`: Specific V allele object (optional)
- `specific_d`: Specific D allele object (optional, only for chains with D segment)
- `specific_j`: Specific J allele object (optional)

### 2. Position Fixing Steps
```python
FixVPositionAfterTrimmingIndexAmbiguity()
FixDPositionAfterTrimmingIndexAmbiguity()
FixJPositionAfterTrimmingIndexAmbiguity()
```

### 3. Correction Steps
```python
CorrectForVEndCut()
CorrectForDTrims()
```

### 4. Sequence Corruption
```python
# Parameters: corruption_prob, events_proba, max_length, add_coeff, remove_coeff, add_after_remove_coeff
CorruptSequenceBeginning(0.7, [0.4, 0.4, 0.2], 576, 210, 310, 50)
```

### 5. N Insertion
```python
# Parameters: n_ratio, probability
InsertNs(0.02, 0.5)
```

### 6. Indel Insertion
```python
# Parameters: indel_probability, max_indels, insertion_proba, deletion_proba
InsertIndels(0.5, 5, 0.5, 0.5)
```

### 7. Validation Steps
```python
ShortDValidation()  # Validates D region length
DistillMutationRate()  # Calculates final mutation rate
```

## Complete Pipeline Examples

### Heavy Chain (Full Pipeline)
```python
AugmentationStep.set_dataconfig(HUMAN_IGH_OGRDB)

pipeline = AugmentationPipeline([
    SimulateSequence(S5F(min_mutation_rate=0.003, max_mutation_rate=0.25), True),
    FixVPositionAfterTrimmingIndexAmbiguity(),
    FixDPositionAfterTrimmingIndexAmbiguity(),
    FixJPositionAfterTrimmingIndexAmbiguity(),
    CorrectForVEndCut(),
    CorrectForDTrims(),
    CorruptSequenceBeginning(0.7, [0.4, 0.4, 0.2], 576, 210, 310, 50),
    InsertNs(0.02, 0.5),
    ShortDValidation(),
    InsertIndels(0.5, 5, 0.5, 0.5),
    DistillMutationRate()
])
```

### Light Chain (Kappa)
```python
AugmentationStep.set_dataconfig(HUMAN_IGK_OGRDB)

pipeline = AugmentationPipeline([
    SimulateSequence(S5F(min_mutation_rate=0.003, max_mutation_rate=0.25), True),
    FixVPositionAfterTrimmingIndexAmbiguity(),
    FixJPositionAfterTrimmingIndexAmbiguity(),
    CorrectForVEndCut(),
    CorruptSequenceBeginning(0.7, [0.4, 0.4, 0.2], 576, 210, 310, 50),
    InsertNs(0.02, 0.5),
    InsertIndels(0.5, 5, 0.5, 0.5),
    DistillMutationRate()
])
```

### TCR Beta Chain
```python
from GenAIRR.steps import FilterTCRDJAmbiguities

AugmentationStep.set_dataconfig(HUMAN_TCRB_IMGT)

pipeline = AugmentationPipeline([
    SimulateSequence(Uniform(), True),
    FixVPositionAfterTrimmingIndexAmbiguity(),
    FixDPositionAfterTrimmingIndexAmbiguity(),
    FixJPositionAfterTrimmingIndexAmbiguity(),
    CorrectForVEndCut(),
    CorrectForDTrims(),
    FilterTCRDJAmbiguities(),  # TCR-specific filtering
    CorruptSequenceBeginning(0.7, [0.4, 0.4, 0.2], 576, 210, 310, 50),
    InsertNs(0.02, 0.5),
    ShortDValidation(),
    InsertIndels(0.5, 5, 0.5, 0.5),
    DistillMutationRate()
])
```

## Mutation Models

### S5F Model
Context-aware somatic hypermutation model based on empirical data.
```python
S5F(min_mutation_rate=0.003, max_mutation_rate=0.25, custom_model=None, productive=False)
```

**Parameters:**
- `min_mutation_rate`: Minimum mutation rate (default: 0)
- `max_mutation_rate`: Maximum mutation rate (default: 0)
- `custom_model`: Path to custom mutation model file (default: None)
- `productive`: Ensure no stop codons and preserve CDR3 anchors (default: False)

### Uniform Model
Simple uniform random mutation model.
```python
Uniform(min_mutation_rate=0.0, max_mutation_rate=0.0, productive=False)
```

**Parameters:**
- `min_mutation_rate`: Minimum mutation rate (default: 0)
- `max_mutation_rate`: Maximum mutation rate (default: 0)
- `productive`: Ensure no stop codons and preserve CDR3 anchors (default: False)

## Accessing Results

```python
result = pipeline.execute()
data = result.get_dict()

# Key result fields
sequence = data['sequence']          # The simulated sequence
productive = data['productive']      # Whether sequence is productive
v_call = data['v_call']             # V gene assignment
d_call = data['d_call']             # D gene assignment
j_call = data['j_call']             # J gene assignment
mutation_rate = data['mutation_rate'] # Final mutation rate
mutations = data['mutations']        # Mutation positions and changes
```

## Custom Allele Selection

```python
# Access specific alleles by family name (returns list of alleles)
# Note: Allele families are organized by gene family groups (e.g., IGHVF1-G1)
v_allele = HUMAN_IGH_OGRDB.v_alleles['IGHVF1-G1'][0]  # First allele in family
d_allele = HUMAN_IGH_OGRDB.d_alleles['IGHD1-1'][0]    # First allele in family
j_allele = HUMAN_IGH_OGRDB.j_alleles['IGHJ1'][0]      # First allele in family

# View available allele families
print("V allele families:", list(HUMAN_IGH_OGRDB.v_alleles.keys())[:5])
print("D allele families:", list(HUMAN_IGH_OGRDB.d_alleles.keys())[:5])
print("J allele families:", list(HUMAN_IGH_OGRDB.j_alleles.keys())[:5])

# View alleles in a specific family
family_alleles = HUMAN_IGH_OGRDB.v_alleles['IGHVF1-G1']
for allele in family_alleles:
    print(allele.name)  # e.g., 'IGHVF1-G1*01', 'IGHVF1-G1*02', etc.

# Use in simulation
custom_step = SimulateSequence(
    S5F(min_mutation_rate=0.003, max_mutation_rate=0.25),
    productive=True,
    specific_v=v_allele,
    specific_d=d_allele,
    specific_j=j_allele
)
```

## Batch Simulation

```python
# Generate multiple sequences
results = []
for _ in range(100):
    result = pipeline.execute()
    results.append(result.get_dict())

# Convert to DataFrame for analysis
import pandas as pd
df = pd.DataFrame(results)
```

## Common Parameter Patterns

- **Probabilities**: Values between 0.0 and 1.0
- **Mutation rates**: Typically 0.001 to 0.3 (0.1% to 30%)
- **Sequence lengths**: Usually 200-600 nucleotides for full sequences
- **Trimming**: 0-20 nucleotides typically
- **Indel counts**: Usually 1-10 per sequence
