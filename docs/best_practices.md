# Best Practices Guide

Guidelines for effective use of GenAIRR in research and production environments.

## Pipeline Design

### 1. Always Set DataConfig First
```python
from GenAIRR.steps.StepBase import AugmentationStep
from GenAIRR.data import HUMAN_IGH_OGRDB

# This should be your first line after imports
AugmentationStep.set_dataconfig(HUMAN_IGH_OGRDB)
```

### 2. Use Realistic Parameter Ranges
Based on biological literature and empirical data:

**Mutation Rates (S5F)**:
- Naive B cells: 0.001 - 0.01 (0.1% - 1%)
- Memory B cells: 0.02 - 0.08 (2% - 8%)
- Plasma cells: 0.05 - 0.25 (5% - 25%)

**Indel Frequencies**:
- Typical: 0.1 - 0.5 probability
- Conservative: 0.05 - 0.2 probability

```python
# Realistic memory B cell simulation
SimulateSequence(S5F(min_mutation_rate=0.02, max_mutation_rate=0.08), productive=True)
```

### 3. Pipeline Order Matters
Follow this recommended step order:

```python
pipeline = AugmentationPipeline([
    # 1. Generate base sequence with mutations
    SimulateSequence(S5F(min_mutation_rate=0.003, max_mutation_rate=0.25), productive=True),
    
    # 2. Fix position ambiguities (order matters!)
    FixVPositionAfterTrimmingIndexAmbiguity(),
    FixDPositionAfterTrimmingIndexAmbiguity(),  # Skip for light chains
    FixJPositionAfterTrimmingIndexAmbiguity(),
    
    # 3. Apply biological corrections
    CorrectForVEndCut(),
    CorrectForDTrims(),  # Skip for light chains
    
    # 4. Add sequencing artifacts
    CorruptSequenceBeginning(0.7, [0.4, 0.4, 0.2], 576, 210, 310, 50),
    InsertNs(0.02, 0.5),
    
    # 5. Validate and filter
    ShortDValidation(),  # Skip for light chains
    
    # 6. Add structural variations
    InsertIndels(0.5, 5, 0.5, 0.5),
    
    # 7. Finalize
    DistillMutationRate()
])
```

## Data Generation

### 1. Batch Processing for Large Datasets
```python
def generate_sequences(n_sequences, batch_size=1000):
    """Generate sequences in batches to manage memory."""
    sequences = []
    
    for i in range(0, n_sequences, batch_size):
        batch_size_actual = min(batch_size, n_sequences - i)
        batch = []
        
        for _ in range(batch_size_actual):
            seq = pipeline.execute()
            batch.append(seq.get_dict())
        
        sequences.extend(batch)
        print(f"Generated {len(sequences)}/{n_sequences} sequences")
    
    return sequences

# Usage
large_dataset = generate_sequences(10000, batch_size=500)
```

### 2. Reproducible Results
```python
import random
import numpy as np

# Set seeds for reproducibility
random.seed(42)
np.random.seed(42)

# Your simulation code here
sequences = [pipeline.execute().get_dict() for _ in range(100)]
```

### 3. Quality Control
Always validate your generated sequences:

```python
def validate_sequences(sequences):
    """Basic quality checks for generated sequences."""
    stats = {
        'total': len(sequences),
        'productive': sum(1 for s in sequences if s.get('productive', False)),
        'avg_length': np.mean([len(s['sequence']) for s in sequences]),
        'avg_mutations': np.mean([len(s.get('mutations', {})) for s in sequences])
    }
    
    print(f"Generated {stats['total']} sequences")
    print(f"Productive: {stats['productive']} ({stats['productive']/stats['total']*100:.1f}%)")
    print(f"Average length: {stats['avg_length']:.1f} bp")
    print(f"Average mutations: {stats['avg_mutations']:.1f}")
    
    return stats

# Use it
stats = validate_sequences(sequences)
```

## Performance Optimization

### 1. Choose Appropriate Chain Types
```python
# For heavy chains (has D segment)
from GenAIRR.data import HUMAN_IGH_OGRDB
AugmentationStep.set_dataconfig(HUMAN_IGH_OGRDB)

# For light chains (no D segment) - faster
from GenAIRR.data import HUMAN_IGK_OGRDB
AugmentationStep.set_dataconfig(HUMAN_IGK_OGRDB)
```

### 2. Optimize Step Parameters
```python
# Faster execution with reasonable quality
pipeline = AugmentationPipeline([
    SimulateSequence(S5F(min_mutation_rate=0.01, max_mutation_rate=0.05), productive=True),  # Lower mutation range
    FixVPositionAfterTrimmingIndexAmbiguity(),
    FixJPositionAfterTrimmingIndexAmbiguity(),
    CorruptSequenceBeginning(0.5, [0.5, 0.5, 0], 400, 150, 200, 30),  # Reduced complexity
    InsertIndels(0.2, 3, 0.5, 0.5),  # Fewer indels
    DistillMutationRate()
])
```

## Research Applications

### 1. Benchmarking Aligners
```python
# Generate ground truth dataset
def create_benchmark_dataset():
    sequences = []
    for mutation_rate in [0.01, 0.05, 0.1, 0.2]:
        step = SimulateSequence(S5F(mutation_rate, mutation_rate), productive=True)
        local_pipeline = AugmentationPipeline([step, ...])
        
        for _ in range(100):
            seq = local_pipeline.execute()
            seq_dict = seq.get_dict()
            seq_dict['true_mutation_rate'] = mutation_rate
            sequences.append(seq_dict)
    
    return sequences
```

### 2. Training ML Models
```python
# Create diverse training data
def create_training_data(n_samples):
    data = []
    
    # Vary key parameters
    mutation_rates = np.random.uniform(0.001, 0.3, n_samples)
    productive_flags = np.random.choice([True, False], n_samples, p=[0.8, 0.2])
    
    for i in range(n_samples):
        step = SimulateSequence(
            S5F(mutation_rates[i] * 0.5, mutation_rates[i]), 
            productive_flags[i]
        )
        
        local_pipeline = AugmentationPipeline([step, ...])
        seq = local_pipeline.execute()
        
        data.append({
            'sequence': seq.sequence,
            'mutation_rate': seq.mutation_rate,
            'productive': seq.productive,
            'label': some_classification_label(seq)
        })
    
    return data
```

## Common Pitfalls to Avoid

1. **Don't forget to set DataConfig** - This is the #1 source of errors
2. **Don't use unrealistic parameters** - Keep mutation rates biological
3. **Don't skip position fix steps** - They ensure metadata accuracy
4. **Don't generate huge datasets in memory** - Use batch processing
5. **Don't ignore productive flag** - Set appropriately for your use case

## Getting Started Checklist

- [ ] Install GenAIRR (`pip install GenAIRR`)
- [ ] Import required modules
- [ ] Set DataConfig with `AugmentationStep.set_dataconfig()`
- [ ] Create pipeline with appropriate steps
- [ ] Test with small dataset first
- [ ] Validate output quality
- [ ] Scale up with batch processing
- [ ] Document your parameters for reproducibility
