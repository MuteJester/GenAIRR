# Best Practices Guide

Guidelines for effective use of GenAIRR in research and production environments.

## Pipeline Design

### 1. Pass Config to Pipeline
```python
from GenAIRR import Pipeline, steps, HUMAN_IGH_OGRDB, S5F

# Config is passed directly to the Pipeline constructor
pipeline = Pipeline(
    config=HUMAN_IGH_OGRDB,
    steps=[...]
)
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
steps.SimulateSequence(S5F(min_mutation_rate=0.02, max_mutation_rate=0.08), productive=True)
```

### 3. Pipeline Order Matters
Follow this recommended step order:

```python
pipeline = Pipeline(
    config=HUMAN_IGH_OGRDB,
    steps=[
        # 1. Generate base sequence with mutations
        steps.SimulateSequence(S5F(min_mutation_rate=0.003, max_mutation_rate=0.25), productive=True),

        # 2. Fix position ambiguities (order matters!)
        steps.FixVPositionAfterTrimmingIndexAmbiguity(),
        steps.FixDPositionAfterTrimmingIndexAmbiguity(),  # Skip for light chains
        steps.FixJPositionAfterTrimmingIndexAmbiguity(),

        # 3. Apply biological corrections
        steps.CorrectForVEndCut(),
        steps.CorrectForDTrims(),  # Skip for light chains

        # 4. Finalize mutation rate
        steps.DistillMutationRate(),

        # 5. Add sequencing artifacts
        steps.CorruptSequenceBeginning(probability=0.7, event_weights=(0.4, 0.4, 0.2)),
        steps.EnforceSequenceLength(max_length=576),
        steps.InsertNs(n_ratio=0.02, probability=0.5),

        # 6. Validate and filter
        steps.ShortDValidation(),  # Skip for light chains

        # 7. Add structural variations
        steps.InsertIndels(probability=0.5, max_indels=5),
    ]
)
```

## Data Generation

### 1. Batch Processing for Large Datasets
```python
def generate_sequences(pipeline, n_sequences, batch_size=1000):
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
large_dataset = generate_sequences(pipeline, 10000, batch_size=500)
```

### 2. Reproducible Results
```python
from GenAIRR import set_seed, get_seed, reset_seed

# Set seed for reproducibility
set_seed(42)

# Your simulation code here
sequences = [pipeline.execute().get_dict() for _ in range(100)]

# Check current seed
print(get_seed())

# Reset to random seed when done
reset_seed()
```

### 3. Quality Control
Always validate your generated sequences:

```python
import numpy as np

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
from GenAIRR import Pipeline, steps, HUMAN_IGH_OGRDB, HUMAN_IGK_OGRDB, S5F

# For heavy chains (has D segment)
heavy_pipeline = Pipeline(
    config=HUMAN_IGH_OGRDB,
    steps=[...]
)

# For light chains (no D segment) - faster
light_pipeline = Pipeline(
    config=HUMAN_IGK_OGRDB,
    steps=[...]
)
```

### 2. Optimize Step Parameters
```python
# Faster execution with reasonable quality
pipeline = Pipeline(
    config=HUMAN_IGK_OGRDB,
    steps=[
        steps.SimulateSequence(S5F(min_mutation_rate=0.01, max_mutation_rate=0.05), productive=True),
        steps.FixVPositionAfterTrimmingIndexAmbiguity(),
        steps.FixJPositionAfterTrimmingIndexAmbiguity(),
        steps.CorrectForVEndCut(),
        steps.DistillMutationRate(),
        steps.CorruptSequenceBeginning(probability=0.5, event_weights=(0.5, 0.5, 0)),
        steps.EnforceSequenceLength(max_length=400),
        steps.InsertIndels(probability=0.2, max_indels=3),
    ]
)
```

## Research Applications

### 1. Benchmarking Aligners
```python
from GenAIRR import Pipeline, steps, HUMAN_IGH_OGRDB, S5F

# Generate ground truth dataset
def create_benchmark_dataset():
    sequences = []
    for mutation_rate in [0.01, 0.05, 0.1, 0.2]:
        pipeline = Pipeline(
            config=HUMAN_IGH_OGRDB,
            steps=[
                steps.SimulateSequence(
                    S5F(min_mutation_rate=mutation_rate, max_mutation_rate=mutation_rate),
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

        for _ in range(100):
            seq = pipeline.execute()
            seq_dict = seq.get_dict()
            seq_dict['true_mutation_rate'] = mutation_rate
            sequences.append(seq_dict)

    return sequences
```

### 2. Training ML Models
```python
import numpy as np
from GenAIRR import Pipeline, steps, HUMAN_IGH_OGRDB, S5F

# Create diverse training data
def create_training_data(n_samples):
    data = []
    mutation_rates = np.random.uniform(0.001, 0.3, n_samples)
    productive_flags = np.random.choice([True, False], n_samples, p=[0.8, 0.2])

    for i in range(n_samples):
        pipeline = Pipeline(
            config=HUMAN_IGH_OGRDB,
            steps=[
                steps.SimulateSequence(
                    S5F(min_mutation_rate=mutation_rates[i] * 0.5, max_mutation_rate=mutation_rates[i]),
                    productive=productive_flags[i]
                ),
                steps.FixVPositionAfterTrimmingIndexAmbiguity(),
                steps.FixDPositionAfterTrimmingIndexAmbiguity(),
                steps.FixJPositionAfterTrimmingIndexAmbiguity(),
                steps.CorrectForVEndCut(),
                steps.CorrectForDTrims(),
                steps.DistillMutationRate(),
            ]
        )
        seq = pipeline.execute()

        data.append({
            'sequence': seq.sequence,
            'mutation_rate': seq.mutation_rate,
            'productive': seq.productive,
        })

    return data
```

## Common Pitfalls to Avoid

!!! danger "Pitfalls that silently corrupt results"
    1. **Unrealistic mutation rates** — keep rates within biological ranges (see table above)
    2. **Skipping position fix steps** — without them, segment boundary annotations are unreliable
    3. **Including D-segment steps for light chains** — causes errors or meaningless annotations
    4. **Measuring mutation rate after artifacts** — place `DistillMutationRate` before `CorruptSequenceBeginning`

!!! warning "Performance pitfalls"
    - **Generating huge datasets in memory** — use batch processing (see above)
    - **`productive=True` with very high mutation rates** — the retry loop can be slow; consider lowering `max_mutation_rate`

## Getting Started Checklist

- [ ] Install GenAIRR (`pip install GenAIRR`)
- [ ] Import required modules (`from GenAIRR import Pipeline, steps, ...`)
- [ ] Create pipeline with `Pipeline(config=..., steps=[...])`
- [ ] Test with small dataset first
- [ ] Validate output quality
- [ ] Scale up with batch processing
- [ ] Use `set_seed()` for reproducibility
