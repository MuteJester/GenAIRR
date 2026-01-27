# Frequently Asked Questions (FAQ)

Common questions about GenAIRR and their answers.

## Getting Started

### Q: What is GenAIRR used for?
**A:** GenAIRR simulates realistic immune receptor sequences (antibodies and T-cell receptors) with full ground truth annotation. It's primarily used for:
- Benchmarking sequence alignment algorithms
- Training machine learning models on immune data
- Studying somatic hypermutation patterns
- Generating synthetic datasets for research

### Q: Do I need biology knowledge to use GenAIRR?

!!! tip "New to immunology?"
    Start with the [Biological Context](biological_context.md) page for a concise primer, then follow the [Step-by-Step Tutorial](step_by_step_tutorial.md).

**A:** Basic understanding helps, but it's not required. Start with:
1. The [Step-by-Step Tutorial](step_by_step_tutorial.md) for hands-on learning
2. The [Biological Context Guide](biological_context.md) for background
3. Use default parameters initially - they're biologically reasonable

### Q: Which Python version does GenAIRR support?
**A:** Python 3.9 or higher. Install with: `pip install GenAIRR`

## Basic Usage

### Q: What's the minimum code to generate a sequence?
**A:** Just 3 lines using the convenience function:
```python
from GenAIRR import simulate, HUMAN_IGH_OGRDB, S5F

result = simulate(HUMAN_IGH_OGRDB, S5F(0.003, 0.25))
print(result.sequence)
```

Or with a pipeline for more control:
```python
from GenAIRR import Pipeline, steps, HUMAN_IGH_OGRDB, S5F

pipeline = Pipeline(
    config=HUMAN_IGH_OGRDB,
    steps=[steps.SimulateSequence(S5F(min_mutation_rate=0.003, max_mutation_rate=0.25), productive=True)]
)
sequence = pipeline.execute()
print(sequence.sequence)
```

### Q: What does "productive=True" mean?
**A:** It ensures the generated sequence is:
- In the correct reading frame
- Free of premature stop codons
- Potentially functional as an antibody

About 1/3 of natural V(D)J recombination events are productive.

### Q: Why do I get an error when creating pipelines?
**A:** Make sure you pass the config to the Pipeline constructor:
```python
pipeline = Pipeline(
    config=HUMAN_IGH_OGRDB,
    steps=[...]
)
```

## Parameters and Configuration

### Q: What mutation rates should I use?
**A:** Depends on the cell type you're modeling:
- **Naive B cells**: 0.001-0.01 (0.1-1%)
- **Memory B cells**: 0.02-0.08 (2-8%)
- **Plasma cells**: 0.05-0.25 (5-25%)

### Q: What's the difference between S5F and Uniform mutation models?
**A:**
- **S5F**: Context-dependent, biologically realistic mutations
- **Uniform**: Simple random mutations at specified rate
- **Recommendation**: Use S5F for research, Uniform for testing

### Q: Can I simulate light chains?
**A:** Yes! Use the appropriate data configuration:
```python
from GenAIRR import Pipeline, steps, HUMAN_IGK_OGRDB, S5F
# or HUMAN_IGL_OGRDB for lambda light chain

pipeline = Pipeline(
    config=HUMAN_IGK_OGRDB,
    steps=[
        steps.SimulateSequence(S5F(min_mutation_rate=0.02, max_mutation_rate=0.08), productive=True),
        steps.FixVPositionAfterTrimmingIndexAmbiguity(),
        steps.FixJPositionAfterTrimmingIndexAmbiguity(),  # No D segment steps
        steps.CorrectForVEndCut(),
        steps.DistillMutationRate(),
    ]
)
```

## Pipeline Design

### Q: Do I need all the pipeline steps?
**A:** No. Start minimal and add complexity:
- **Minimal**: Just `SimulateSequence`
- **Basic**: Add position fix steps
- **Realistic**: Add biological corrections
- **Full**: Add sequencing artifacts

### Q: What order should pipeline steps be in?
**A:** Follow this general order:
1. `SimulateSequence` (always first)
2. Position fixes (`FixVPositionAfterTrimmingIndexAmbiguity`, etc.)
3. Biological corrections (`CorrectForVEndCut`, etc.)
4. Finalization (`DistillMutationRate`)
5. Sequencing artifacts (`CorruptSequenceBeginning`, `EnforceSequenceLength`, `InsertNs`)
6. Quality control (`ShortDValidation`)
7. Structural variants (`InsertIndels`)

### Q: Can I create custom pipeline steps?
**A:** Yes! Inherit from `AugmentationStep` and implement the `apply` method:
```python
from GenAIRR.steps.StepBase import AugmentationStep

class MyCustomStep(AugmentationStep):
    def apply(self, container):
        # Your custom logic here
        container.sequence = container.sequence.upper()
```

## Data and Output

### Q: What data does GenAIRR output?
**A:** Each simulated sequence includes:
- DNA sequence string
- V, D, J allele names used
- Mutation positions and types
- Sequence region boundaries
- Quality metrics (productive, mutation rate, etc.)

### Q: How do I export results to different formats?
**A:**
```python
# Pandas DataFrame
import pandas as pd
df = pd.DataFrame([seq.get_dict() for seq in sequences])

# FASTA format
with open('output.fasta', 'w') as f:
    for i, seq in enumerate(sequences):
        f.write(f">seq_{i}\n{seq.sequence}\n")

# JSON format
import json
with open('output.json', 'w') as f:
    json.dump([seq.get_dict() for seq in sequences], f)
```

### Q: Can I use my own germline database?
**A:** Yes, but it requires creating a custom DataConfig. See the [Custom Data Config](custom_data_config.md) guide.

## Performance and Scaling

### Q: How fast is GenAIRR?
**A:** Speed depends on complexity:
- Simple pipeline: ~100-1000 sequences/second
- Full pipeline: ~10-100 sequences/second
- With high mutation rates or `productive=True`: slower due to retries

### Q: How do I generate large datasets efficiently?
**A:** Use batch processing:
```python
def generate_batch(pipeline, batch_size=1000):
    return [pipeline.execute().get_dict() for _ in range(batch_size)]

# Generate 10,000 sequences in batches
all_sequences = []
for i in range(10):
    batch = generate_batch(pipeline, 1000)
    all_sequences.extend(batch)
    print(f"Generated {len(all_sequences)} sequences")
```

### Q: Why is generation slow when using productive=True?
**A:** The library regenerates sequences until finding productive ones. Solutions:
- Use lower mutation rates
- Accept some non-productive sequences (`productive=False`)
- Use pre-filtered germline databases

## Troubleshooting

### Q: My sequences all look the same!
**A:** Check these:
- Mutation rates aren't zero: `S5F(min_mutation_rate=0.01, max_mutation_rate=0.05)` not `S5F(0, 0)`
- Using different alleles: Check if you're forcing specific alleles
- Random seed: Don't set a fixed seed for production use

### Q: I'm getting very short sequences!
**A:** Adjust corruption and length parameters:
```python
# Less aggressive corruption
steps.CorruptSequenceBeginning(probability=0.3, event_weights=(0.7, 0.3, 0))
steps.EnforceSequenceLength(max_length=400)
```

### Q: How do I reproduce results?
**A:** Use GenAIRR's built-in seed management:
```python
from GenAIRR import set_seed, get_seed, reset_seed

set_seed(42)
# Now generate sequences...
```

## Advanced Usage

### Q: Can I simulate paired heavy/light chains?
**A:** Not directly, but you can generate them separately using different pipelines:
```python
from GenAIRR import Pipeline, steps, HUMAN_IGH_OGRDB, HUMAN_IGK_OGRDB, S5F

heavy_pipeline = Pipeline(
    config=HUMAN_IGH_OGRDB,
    steps=[steps.SimulateSequence(S5F(min_mutation_rate=0.02, max_mutation_rate=0.08), productive=True)]
)

light_pipeline = Pipeline(
    config=HUMAN_IGK_OGRDB,
    steps=[steps.SimulateSequence(S5F(min_mutation_rate=0.02, max_mutation_rate=0.08), productive=True)]
)

heavy = heavy_pipeline.execute()
light = light_pipeline.execute()
```

### Q: How do I model specific diseases or conditions?
**A:** Adjust parameters to reflect biology:
```python
# Autoimmune (higher mutation)
steps.SimulateSequence(S5F(min_mutation_rate=0.05, max_mutation_rate=0.15), productive=True)

# Immunodeficiency (lower diversity - use specific alleles)
steps.SimulateSequence(S5F(min_mutation_rate=0.001, max_mutation_rate=0.02), productive=True, specific_v=common_allele)
```

### Q: Can I add custom mutation patterns?
**A:** Yes, by creating custom mutation models. See the source code of `S5F` and `Uniform` classes as examples.

## Getting Help

### Q: Where can I find more examples?
**A:** Check these resources:
- [Jupyter notebook tutorials](tutorials/Quick Start Guide.ipynb)
- [GitHub repository](https://github.com/MuteJester/GenAIRR)
- [Step-by-step tutorial](step_by_step_tutorial.md)

### Q: How do I report bugs or request features?
**A:**
1. Check existing GitHub issues first
2. Create a minimal reproducible example
3. Include your Python version and GenAIRR version
4. Submit to the GitHub repository

### Q: Is there a community forum or chat?
**A:** Check the GitHub repository for current community resources and discussion channels.
