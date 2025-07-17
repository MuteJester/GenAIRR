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
**A:** Basic understanding helps, but it's not required. Start with:
1. The [Step-by-Step Tutorial](step_by_step_tutorial.md) for hands-on learning
2. The [Biological Context Guide](biological_context.md) for background
3. Use default parameters initially - they're biologically reasonable

### Q: Which Python version does GenAIRR support?
**A:** Python 3.9 or higher. Install with: `pip install GenAIRR`

## Basic Usage

### Q: What's the minimum code to generate a sequence?
**A:** Just 6 lines:
```python
from GenAIRR.pipeline import AugmentationPipeline
from GenAIRR.steps import SimulateSequence
from GenAIRR.mutation import S5F
from GenAIRR.data import HUMAN_IGH_OGRDB
from GenAIRR.steps.StepBase import AugmentationStep

AugmentationStep.set_dataconfig(HUMAN_IGH_OGRDB)
pipeline = AugmentationPipeline([SimulateSequence(S5F(), True)])
sequence = pipeline.execute()
print(sequence.sequence)
```

### Q: What does "productive=True" mean?
**A:** It ensures the generated sequence is:
- In the correct reading frame
- Free of premature stop codons
- Potentially functional as an antibody

About 1/3 of natural V(D)J recombination events are productive.

### Q: Why do I get an AttributeError when creating pipelines?
**A:** You forgot to set the data configuration. Always call this first:
```python
AugmentationStep.set_dataconfig(HUMAN_IGH_OGRDB)
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
from GenAIRR.data import HUMAN_IGK_OGRDB  # Kappa light chain
# or
from GenAIRR.data import HUMAN_IGL_OGRDB  # Lambda light chain

AugmentationStep.set_dataconfig(HUMAN_IGK_OGRDB)
# Note: Skip D-segment steps in your pipeline
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
4. Sequencing artifacts (`CorruptSequenceBeginning`, `InsertNs`, etc.)
5. Quality control (`ShortDValidation`)
6. Structural variants (`InsertIndels`)
7. Finalization (`DistillMutationRate`)

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
- Mutation rates aren't zero: `S5F(0.01, 0.05)` not `S5F(0, 0)`
- Using different alleles: Check if you're forcing specific alleles
- Random seed: Don't set `random.seed()` for production use

### Q: I'm getting very short sequences!
**A:** Adjust `CorruptSequenceBeginning` parameters:
```python
# Less aggressive corruption
CorruptSequenceBeginning(0.3, [0.7, 0.3, 0], 400, 100, 150, 20)
```

### Q: How do I reproduce results?
**A:** Set random seeds before generation:
```python
import random
import numpy as np

random.seed(42)
np.random.seed(42)
# Now generate sequences...
```

## Advanced Usage

### Q: Can I simulate paired heavy/light chains?
**A:** Not directly, but you can generate them separately:
```python
# Heavy chain
AugmentationStep.set_dataconfig(HUMAN_IGH_OGRDB)
heavy = heavy_pipeline.execute()

# Light chain  
AugmentationStep.set_dataconfig(HUMAN_IGK_OGRDB)
light = light_pipeline.execute()

# Pair them logically in your analysis
```

### Q: How do I model specific diseases or conditions?
**A:** Adjust parameters to reflect biology:
```python
# Autoimmune (higher mutation)
SimulateSequence(S5F(0.05, 0.15), True)

# Immunodeficiency (lower diversity - use specific alleles)
SimulateSequence(S5F(0.001, 0.02), True, specific_v=common_allele)
```

### Q: Can I add custom mutation patterns?
**A:** Yes, by creating custom mutation models. See the source code of `S5F` and `Uniform` classes as examples.

## Getting Help

### Q: Where can I find more examples?
**A:** Check these resources:
- [Jupyter notebooks in docs/tutorials/](tutorials/)
- [Examples in the README](../README.md)
- [Step-by-step tutorial](step_by_step_tutorial.md)

### Q: How do I report bugs or request features?
**A:** 
1. Check existing GitHub issues first
2. Create a minimal reproducible example
3. Include your Python version and GenAIRR version
4. Submit to the GitHub repository

### Q: Is there a community forum or chat?
**A:** Check the GitHub repository for current community resources and discussion channels.
