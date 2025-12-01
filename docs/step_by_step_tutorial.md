# Step-by-Step Tutorial: Building Your First Pipeline

This tutorial walks you through creating a GenAIRR simulation pipeline from scratch, explaining each decision and parameter choice.

## Step 1: Environment Setup

```python
# Essential imports
from GenAIRR.pipeline import AugmentationPipeline
from GenAIRR.steps import SimulateSequence
from GenAIRR.steps.StepBase import AugmentationStep
from GenAIRR.mutation import S5F
from GenAIRR.data import HUMAN_IGH_OGRDB
```

**Why these imports?**
- `AugmentationPipeline`: The core pipeline runner
- `SimulateSequence`: The fundamental sequence generation step
- `AugmentationStep`: Base class that manages data configuration
- `S5F`: Biologically realistic mutation model
- `HUMAN_IGH_OGRDB`: Pre-built heavy chain germline database

## Step 2: Configure Data Source

```python
# CRITICAL: Always set this first!
AugmentationStep.set_dataconfig(HUMAN_IGH_OGRDB)
```

**What this does:**
- Loads human heavy chain V, D, J, and C gene segments
- Configures trimming probabilities based on empirical data
- Sets up nucleotide addition patterns for junctions
- Enables all subsequent steps to access this data

**Common mistake:** Forgetting this step leads to AttributeError crashes.

## Step 3: Choose Your Mutation Model

```python
# Option 1: Realistic somatic hypermutation
mutation_model = S5F(min_mutation_rate=0.02, max_mutation_rate=0.08)

# Option 2: Simple uniform mutations  
# from GenAIRR.mutation import Uniform
# mutation_model = Uniform(min_mutation_rate=0.02, max_mutation_rate=0.08)
```

**Parameter guidance:**
- 0.02-0.08 (2-8%): Memory B cells
- 0.001-0.01 (0.1-1%): Naive B cells  
- 0.05-0.25 (5-25%): Plasma cells

**Why S5F?** It models context-dependent mutation patterns seen in real somatic hypermutation.

## Step 4: Create Core Simulation Step

```python
# Generate the basic sequence
simulate_step = SimulateSequence(
    mutation_model,
    productive=True  # Ensures functional sequences
)
```

**The `productive` parameter:**
- `True`: Only generates in-frame sequences without stop codons
- `False`: Allows non-functional sequences (includes ~2/3 of all possible sequences)

## Step 5: Build Your Pipeline

### Minimal Pipeline (Good for Testing)
```python
minimal_pipeline = AugmentationPipeline([
    simulate_step
])

# Test it
result = minimal_pipeline.execute()
print("Basic sequence:", result.sequence[:50] + "...")
```

### Realistic Pipeline (Good for Research)
```python
from GenAIRR.steps import (FixVPositionAfterTrimmingIndexAmbiguity,
                          FixDPositionAfterTrimmingIndexAmbiguity, 
                          FixJPositionAfterTrimmingIndexAmbiguity,
                          CorrectForVEndCut, CorrectForDTrims)

realistic_pipeline = AugmentationPipeline([
    simulate_step,
    FixVPositionAfterTrimmingIndexAmbiguity(),
    FixDPositionAfterTrimmingIndexAmbiguity(), 
    FixJPositionAfterTrimmingIndexAmbiguity(),
    CorrectForVEndCut(),
    CorrectForDTrims()
])
```

### Full Pipeline (Production Ready)
```python
from GenAIRR.steps import (CorruptSequenceBeginning, InsertNs, 
                          InsertIndels, ShortDValidation, DistillMutationRate)

full_pipeline = AugmentationPipeline([
    simulate_step,
    
    # Position corrections
    FixVPositionAfterTrimmingIndexAmbiguity(),
    FixDPositionAfterTrimmingIndexAmbiguity(),
    FixJPositionAfterTrimmingIndexAmbiguity(),
    
    # Biological corrections
    CorrectForVEndCut(),
    CorrectForDTrims(),
    
    # Sequencing artifacts
    CorruptSequenceBeginning(0.7, [0.4, 0.4, 0.2], 576, 210, 310, 50),
    InsertNs(0.02, 0.5),
    
    # Quality control and variants
    ShortDValidation(),
    InsertIndels(0.5, 5, 0.5, 0.5),
    
    # Finalization
    DistillMutationRate()
])
```

## Step 6: Generate and Analyze Sequences

```python
# Generate a single sequence
sequence_container = full_pipeline.execute()
sequence_data = sequence_container.get_dict()

# Examine the output
print("Sequence length:", len(sequence_data['sequence']))
print("Mutation count:", len(sequence_data['mutations']))
print("Is productive:", sequence_data['productive'])
print("V allele used:", sequence_data['v_call'])
print("Mutation rate:", sequence_data['mutation_rate'])
```

### Generate Multiple Sequences
```python
# Generate a small dataset
sequences = []
for i in range(10):
    result = full_pipeline.execute()
    seq_dict = result.get_dict()
    seq_dict['id'] = f"seq_{i:03d}"  # Add unique identifier
    sequences.append(seq_dict)

# Basic statistics
productive_count = sum(1 for s in sequences if s['productive'])
avg_mutations = sum(len(s['mutations']) for s in sequences) / len(sequences)

print(f"Generated {len(sequences)} sequences")
print(f"Productive: {productive_count}/{len(sequences)} ({productive_count/len(sequences)*100:.1f}%)")
print(f"Average mutations per sequence: {avg_mutations:.1f}")
```

## Step 7: Export Your Data

```python
import pandas as pd
import json

# Option 1: Pandas DataFrame
df = pd.DataFrame(sequences)
df.to_csv('simulated_sequences.csv', index=False)

# Option 2: JSON format
with open('simulated_sequences.json', 'w') as f:
    json.dump(sequences, f, indent=2)

# Option 3: FASTA format (sequences only)
with open('simulated_sequences.fasta', 'w') as f:
    for seq in sequences:
        f.write(f">{seq['id']}\n{seq['sequence']}\n")
```

## Common Customizations

### 1. Specific Allele Usage
```python
# Force specific V, D, J combination
# Note: Allele families use IGHVF (family) naming in OGRDB
specific_v = HUMAN_IGH_OGRDB.v_alleles['IGHVF1-G1'][0]  # First allele in family
specific_d = HUMAN_IGH_OGRDB.d_alleles['IGHD1-1'][0]
specific_j = HUMAN_IGH_OGRDB.j_alleles['IGHJ1'][0]

custom_step = SimulateSequence(
    S5F(min_mutation_rate=0.02, max_mutation_rate=0.08), 
    productive=True,
    specific_v=specific_v,
    specific_d=specific_d,
    specific_j=specific_j
)
```

### 2. Light Chain Simulation
```python
from GenAIRR.data import HUMAN_IGK_OGRDB

# Switch to kappa light chain
AugmentationStep.set_dataconfig(HUMAN_IGK_OGRDB)

# Light chain pipeline (no D segment steps)
light_pipeline = AugmentationPipeline([
    SimulateSequence(S5F(min_mutation_rate=0.02, max_mutation_rate=0.08), productive=True),
    FixVPositionAfterTrimmingIndexAmbiguity(),
    FixJPositionAfterTrimmingIndexAmbiguity(),  # No D segment
    CorrectForVEndCut(),
    # Skip CorrectForDTrims and ShortDValidation
    CorruptSequenceBeginning(0.7, [0.4, 0.4, 0.2], 400, 150, 200, 30),
    InsertNs(0.02, 0.5),
    InsertIndels(0.3, 3, 0.5, 0.5),
    DistillMutationRate()
])
```

### 3. Naive vs Memory Comparison
```python
# Naive B cells (low mutation)
naive_step = SimulateSequence(S5F(min_mutation_rate=0.001, max_mutation_rate=0.01), productive=True)

# Memory B cells (moderate mutation)  
memory_step = SimulateSequence(S5F(min_mutation_rate=0.02, max_mutation_rate=0.08), productive=True)

# Generate both types
naive_pipeline = AugmentationPipeline([naive_step, ...])
memory_pipeline = AugmentationPipeline([memory_step, ...])

naive_seqs = [naive_pipeline.execute().get_dict() for _ in range(50)]
memory_seqs = [memory_pipeline.execute().get_dict() for _ in range(50)]
```

## Troubleshooting Your Pipeline

### Pipeline Not Working?
1. **Check DataConfig**: Did you call `AugmentationStep.set_dataconfig()`?
2. **Check imports**: Are all step classes imported?
3. **Check parameters**: Are mutation rates reasonable (0.001-0.3)?

### Sequences Look Wrong?
1. **Low diversity**: Increase mutation rates or check allele usage
2. **Too many non-productive**: Set `productive=True` in SimulateSequence
3. **Sequences too short/long**: Adjust CorruptSequenceBeginning parameters

### Performance Issues?
1. **Slow generation**: Lower mutation rates or remove complex steps
2. **Memory problems**: Generate sequences in batches
3. **Low productivity**: Use `productive=True` to avoid retries

## Next Steps

1. **Read the Parameter Reference** for detailed parameter explanations
2. **Check Best Practices Guide** for optimization tips
3. **Explore Advanced Tutorials** for custom step development
4. **Review Biological Context** to understand what you're simulating

Congratulations! You now have a working GenAIRR pipeline. Start with simple parameters and gradually add complexity as you understand each component.
