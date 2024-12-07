# GenAIRR: Adaptive Immune Receptor Repertoire Sequence Simulator

GenAIRR is an advanced Python framework for generating synthetic Adaptive Immune Receptor Repertoire (AIRR) sequences. It facilitates benchmarking of alignment algorithms, comprehensive sequence analysis, and exploration of immune receptor diversity.

## Key Features
- **Realistic AIRR Sequence Simulation**: Generate synthetic heavy and light chain sequences with robust customization options.
- **Customizable Augmentation Pipelines**: Build, visualize, and customize sequence simulation pipelines with the `AugmentedPipeline`.
- **Mutation and Indel Simulation**: Incorporates realistic mutation models and indel simulation to capture biological variability.
- **Precision Allele-Specific Handling**: Employ correction maps to manage allele-specific trimming and ambiguities accurately.
- **Support for TCR Sequences**: Simulate TCR-Beta sequences along with BCR sequences for diverse research needs.
- **Custom Noise and Mutation Models**: Incorporates realistic mutation and noise models and allows for usage of custom models.
- **Easily Build and Add Custom Steps**: Implement your own steps and add them easily to your pipeline to further customize your sequnces.


<p align="center">
  <h1 align="center" style="font-size: 3em;">GenAIRR Documentation</h1>
  <p align="center" style="font-size: 1.5em;">
    Access detailed documentation, guides, and examples.
    <br />
    <a href="https://genairr.readthedocs.io/en/latest/"><strong>Explore the Docs »</strong></a>
    <br />
  </p>
</p>

## Acknowledgements
Portions of this project were inspired by [AIRRship](https://github.com/Cowanlab/airrship).

---

## Quick Start Guide

This guide provides a concise overview of setting up GenAIRR, simulating sequences, and customizing your augmentation pipeline.

### Installation

Ensure Python 3.9 is installed and install GenAIRR using pip:

```bash
pip install GenAIRR
```

### Setting Up Your Environment

Begin by importing essential modules and initializing the `AugmentedPipeline`:

```python
from GenAIRR.pipeline import AugmentationPipeline
from GenAIRR.steps import AugmentationStep
from GenAIRR.utilities import DataConfig
from GenAIRR.data import builtin_heavy_chain_data_config, builtin_kappa_chain_data_config
from GenAIRR.pipeline import CHAIN_TYPE_BCR_HEAVY

# Initialize a built-in DataConfig
data_config_builtin = builtin_heavy_chain_data_config()

# Set up the data configuration and chain type for the simulations
AugmentationStep.set_dataconfig(data_config_builtin,chain_type=CHAIN_TYPE_BCR_HEAVY)
```

### Design Your Pipeline
here is an example of the custom BCR heavy chain GenAIRR pipeline:
```python
from GenAIRR.pipeline import AugmentationPipeline
from GenAIRR.steps import SimulateSequence, FixVPositionAfterTrimmingIndexAmbiguity
from GenAIRR.mutation import S5F
from GenAIRR.data import builtin_heavy_chain_data_config
from GenAIRR.steps.StepBase import AugmentationStep
from GenAIRR.pipeline import CHAIN_TYPE_BCR_HEAVY
from GenAIRR.steps import SimulateSequence,FixVPositionAfterTrimmingIndexAmbiguity,FixDPositionAfterTrimmingIndexAmbiguity,FixJPositionAfterTrimmingIndexAmbiguity
from GenAIRR.steps import CorrectForVEndCut,CorrectForDTrims,CorruptSequenceBeginning,InsertNs,InsertIndels,ShortDValidation,DistillMutationRate
from GenAIRR.mutation import S5F

pipeline = AugmentationPipeline([
    SimulateSequence(mutation_model = S5F(min_mutation_rate=0.003,max_mutation_rate=0.25),productive = True),
    FixVPositionAfterTrimmingIndexAmbiguity(),
    FixDPositionAfterTrimmingIndexAmbiguity(),
    FixJPositionAfterTrimmingIndexAmbiguity(),
    CorrectForVEndCut(),
    CorrectForDTrims(),
    CorruptSequenceBeginning(corruption_probability = 0.7,corrupt_events_proba = [0.4,0.4,0.2],max_sequence_length = 576,nucleotide_add_coefficient = 210,
                             nucleotide_remove_coefficient = 310,nucleotide_add_after_remove_coefficient = 50,random_sequence_add_proba = 1,
                             single_base_stream_proba = 0,duplicate_leading_proba = 0,random_allele_proba = 0),
    InsertNs(n_ratio = 0.02,proba = 0.5),
    ShortDValidation(short_d_length= 5),
    InsertIndels(indel_probability = 0.5,max_indels = 5,insertion_proba=0.5,deletion_proba=0.5),
    DistillMutationRate()
    ])
```

### Simulating a single Heavy Chain Sequence

Use the `AugmentedPipeline` to simulate sequences with ease:

```python
# Simulate a heavy chain sequence using the pipeline
simulation_result = pipeline.execute()

print("Simulated Heavy Chain Sequence:", simulation_result.get_dict())
```

### Customizing the Pipeline

GenAIRR allows extensive customization of the augmentation pipeline. Add steps such as mutations, indel simulation, or sequence corruption:

```python
# Customize the pipeline by adding steps
pipeline = AugmentationPipeline([
    SimulateSequence(mutation_model = S5F(min_mutation_rate=0.003,max_mutation_rate=0.25),productive = True),
    FixVPositionAfterTrimmingIndexAmbiguity(),
    FixDPositionAfterTrimmingIndexAmbiguity(),
    FixJPositionAfterTrimmingIndexAmbiguity(),
    CorrectForVEndCut(),
    CorrectForDTrims(),
    CUSTOM_STEP_HERE(), # example
    CorruptSequenceBeginning(corruption_probability = 0.7,corrupt_events_proba = [0.4,0.4,0.2],max_sequence_length = 576,nucleotide_add_coefficient = 210,
                             nucleotide_remove_coefficient = 310,nucleotide_add_after_remove_coefficient = 50,random_sequence_add_proba = 1,
                             single_base_stream_proba = 0,duplicate_leading_proba = 0,random_allele_proba = 0),
    DistillMutationRate(),
    CUSTOM_STEP_HERE()# example
    ])

# Run the customized pipeline
custom_augmented_sequence = pipeline.execute()
print("Custom Augmented Sequence:", custom_augmented_sequence.get_dict())
```
## Generating Naïve Sequences

In immunogenetics, a naïve sequence refers to an antibody sequence that has not undergone the process of somatic hypermutation. GenAIRR allows you to simulate such naïve sequences using the `HeavyChainSequence` class. Let's start by generating a naïve heavy chain sequence.



```python
from GenAIRR.mutation import Uniform
naive_simulation_step = SimulateSequence(mutation_model = Uniform(min_mutation_rate=0,max_mutation_rate=0),productive = True)
# notice mutation rate min and max are set to 0 so mutation will not be generated

pipeline = AugmentationPipeline([
    naive_simulation_step,
    FixVPositionAfterTrimmingIndexAmbiguity(),
    FixDPositionAfterTrimmingIndexAmbiguity(),
    FixJPositionAfterTrimmingIndexAmbiguity(),
    CorrectForVEndCut(),
    CorrectForDTrims(),
    ])

# Create a naive heavy chain sequence
naive_heavy_sequence = HeavyChainSequence.create_random(data_config_builtin)

# Access the generated naive sequence
naive_sequence = naive_heavy_sequence

print("Naïve Heavy Chain Sequence:", naive_sequence)
print('Ungapped Sequence: ')
print(naive_sequence.ungapped_seq)

```



## Mutation Models

GenAIRR offers multiple mutation models to replicate real-world diversity:

- **S5F Mutation Model**: Context-specific model for realistic somatic hypermutation.
- **Uniform Mutation Model**: Applies a uniform mutation rate across sequences.

```python
from GenAIRR.mutation import S5F, Uniform

# Initialize and apply S5F mutation model
s5f_model = S5F(min_mutation_rate=0.01, max_mutation_rate=0.05)
s5f_mutated_sequence, mutations, mutation_rate = s5f_model.apply_mutation(naive_sequence)

print("S5F Mutated Sequence:", s5f_mutated_sequence)

# Initialize and apply Uniform mutation model
uniform_model = Uniform(min_mutation_rate=0.01, max_mutation_rate=0.05)
uniform_mutated_sequence, uniform_mutations, uniform_rate = uniform_model.apply_mutation(naive_sequence)

print("Uniform Mutated Sequence:", uniform_mutated_sequence)
```

## Advanced Use Cases

### Generating Multiple Sequences

Create batches of sequences for data analysis or benchmarking:

```python
sequences = [pipeline.execute() for _ in range(5)]
print("Generated Sequences:", sequences)
```

### Custom Allele Combinations

Simulate sequences with specific V, D, and J alleles:

```python
# Define specific alleles
v_allele = 'IGHVF3-G10*01'
d_allele = 'IGHD5-12*01'
j_allele = 'IGHJ6*03'

modified_sequence_simulation_step = SimulateSequence(mutation_model = S5F(min_mutation_rate=0.003,max_mutation_rate=0.25),
                                                     productive = True,
                                                     specific_v = v_allele,
                                                     specific_d = d_allele,
                                                     specific_j = j_allele)

# Simulate a sequence with the chosen allele combination
pipeline = AugmentationPipeline([
    modified_sequence_simulation_step,
    FixVPositionAfterTrimmingIndexAmbiguity(),
    FixDPositionAfterTrimmingIndexAmbiguity(),
    FixJPositionAfterTrimmingIndexAmbiguity(),
    CorrectForVEndCut(),
    CorrectForDTrims(),
    CorruptSequenceBeginning(corruption_probability = 0.7,corrupt_events_proba = [0.4,0.4,0.2],max_sequence_length = 576,nucleotide_add_coefficient = 210,
                             nucleotide_remove_coefficient = 310,nucleotide_add_after_remove_coefficient = 50,random_sequence_add_proba = 1,
                             single_base_stream_proba = 0,duplicate_leading_proba = 0,random_allele_proba = 0),
    InsertNs(n_ratio = 0.02,proba = 0.5),
    ShortDValidation(short_d_length= 5),
    InsertIndels(indel_probability = 0.5,max_indels = 5,insertion_proba=0.5,deletion_proba=0.5),
    DistillMutationRate()
    ])
custom_sequence = pipeline.execute()
print("Specific Allele Combination Sequence:", custom_sequence.get_dict())
```


## Conclusion

GenAIRR is designed for a broad range of immunogenetics research, from generating large datasets to specific allele simulations and detailed mutation modeling. Explore the full capabilities of GenAIRR in the [official documentation](https://genairr.readthedocs.io/en/latest/).
