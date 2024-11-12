# Getting Started with GenAIRR

**GenAIRR** is a comprehensive and modular framework designed for simulating, analyzing, and enhancing immunoglobulin (Ig) sequence data. It facilitates the creation of sophisticated Ig sequence simulations, supporting both basic and complex scenarios. Whether you are studying immune responses, developing algorithms for sequence analysis, or benchmarking data, GenAIRR offers the tools you need to achieve precise, customized simulations.

## Overview of GenAIRR

GenAIRR's core architecture is modular, allowing flexibility in building, modifying, and visualizing Ig simulation pipelines. Its components can be used independently or combined to create full simulation workflows. These modules are organized into key categories:

1. **Simulation**: Core classes that help construct and simulate sequences, including alleles and mutation models.
2. **Correction and Validation**: Tools and steps for refining and validating sequences to ensure reliable ground truth data.
3. **Pipeline Integration**: The framework for designing complex pipelines using modular steps and easily integrating custom augmenters.

## Tutorials and Examples

To start building simulations and exploring GenAIRR's potential, refer to the following resources:

- **[Introductory Tutorial](tutorials/Quick%20Start%20Guide.ipynb)**: Get started by creating your first GenAIRR simulation pipeline.
- **[Advanced Tutorial](tutorials/Advanced%20Custom%20Generation.ipynb)**: Browse a collection of example models that showcase different capabilities of GenAIRR, providing insights into how to build and extend the framework for specific needs.

## Building a GenAIRR Pipeline

GenAIRR supports customizable pipelines through a set of modular steps, each of which can be modified or extended. Steps such as `SimulateSequence`, `CorrectForVEndCut`, and `InsertIndels` are examples of augmenters that can be combined to create pipelines tailored to various research goals.

### Example of a Basic Pipeline

Below is a simple pipeline that demonstrates how to create a sequence simulation in GenAIRR:

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

# Set up the data configuration and pipeline
AugmentationStep.set_dataconfig(builtin_heavy_chain_data_config(),chain_type=CHAIN_TYPE_BCR_HEAVY)
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


# Run the pipeline and retrieve the simulated sequence
simulation_result = pipeline.execute()
print(simulation_result.get_dict())
```

### Customizing and Adding Steps

GenAIRR's modular nature allows users to insert their own custom steps at any stage of the pipeline. This makes it possible to introduce new augmenters, mutation models, or validation processes. See the [tutorial on extending GenAIRR](tutorials/extend_tutorial) for guidance on implementing custom steps.

## Understanding Key GenAIRR Components

### Sequence and Allele Classes

GenAIRR includes built-in classes that represent different Ig sequence types, such as `HeavyChainSequence` and `LightChainSequence`. These classes encapsulate methods for generating and mutating sequences, ensuring flexible customization for various Ig chain types.

### Correction and Validation Modules

Validation and correction steps are essential to maintain accurate sequence data. GenAIRR includes modules like `CorrectForDTrims` and `FixJPositionAfterTrimmingIndexAmbiguity`, which automatically adjust sequence metadata to reflect trimming or mutation changes.

### Mutation Models

Built-in mutation models such as `S5F` simulate realistic mutation patterns based on empirical data. Custom models can also be defined to match specific research needs.

## Visualization and Analysis

Visualizing simulation pipelines is easy with GenAIRR. The framework includes tools for generating GraphViz representations of pipeline structures, making it straightforward to understand the flow of steps and the transformations applied at each stage.

### Example of Pipeline Visualization

```python
pipeline.plot('pipeline_diagram')
```

This will save a diagram of the pipeline structure as an image, allowing you to visualize the sequential process.

## Best Practices and Recommendations

For optimal use of GenAIRR:

- **Utilize Data Configurations**: Take advantage of built-in data configurations (`builtin_heavy_chain_data_config`, etc.) for empirical distributions and allele references.
- **Verify Step Outputs**: Ensure that each pipeline step produces the expected results by testing with sample data.
- **Explore Customization**: Don't hesitate to implement your own steps if the built-in ones don’t fully meet your needs.

## Further Resources

### Community and Support

- **[GitHub Issues](https://github.com/MuteJester/GenAIRR/issues)**: Report bugs or request features.

