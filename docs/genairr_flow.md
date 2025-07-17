# Structure and Flow of Simulation in GenAIRR

GenAIRR provides a robust framework for simulating immunoglobulin (Ig) sequences, offering flexibility and precision through its modular structure. At the heart of this framework are `DataConfig` instances, `Sequence` classes like `HeavyChainSequence`, and customizable pipelines composed of augmentation steps. This page explores the general structure and flow of simulations in GenAIRR, emphasizing the role of `DataConfig`, how to create and run simulations, and how to develop custom augmentation steps to tailor simulations to specific needs.

## The Importance of `DataConfig`

The `DataConfig` object is the cornerstone of any GenAIRR simulation. It encapsulates all essential distributions, empirical data, and validation mappings required for realistic Ig sequence simulations. A `DataConfig` instance holds:
- **Allelic distributions**: Information about V, D, J alleles and their usage probabilities.
- **Trim and mutation distributions**: Probabilities and distributions for nucleotide trimming and mutation rates.
- **Correction maps**: Objects that help correct and validate generated sequences by resolving ambiguities.
- **Empirical data**: Data derived from real Ig sequence observations for accurate modeling.

### Example: Loading a `DataConfig`
To load a built-in `DataConfig` for heavy chains:
```python
from GenAIRR.data import HUMAN_IGH_OGRDB

heavy_chain_config = HUMAN_IGH_OGRDB
```

This instance can now be used as the foundation for generating sequences and building simulation pipelines.

## Basic Sequence Simulation

With a `DataConfig` instance, users can directly create naive sequences using predefined classes such as `HeavyChainSequence`. This approach is ideal for generating simple, unguided sequences with realistic characteristics.

### Example: Creating a Naive Heavy Chain Sequence
```python
from GenAIRR.sequence import HeavyChainSequence
from GenAIRR.data import HUMAN_IGH_OGRDB

config = HUMAN_IGH_OGRDB
naive_sequence = HeavyChainSequence.create_random(config)

print("Generated Sequence:", naive_sequence.mutated_seq)
```

This code snippet creates a naive heavy chain sequence using the provided configuration. While generating naive sequences can be insightful, more complex simulations often require a sequence to undergo various augmentation steps for better analysis and study.

## Custom Simulation Pipelines

The real power of GenAIRR lies in its ability to create customizable pipelines composed of different augmentation steps. Each pipeline is essentially a sequence of steps that modifies a `SimulationContainer` instance to reflect specific simulation logic.

### Structure of a Pipeline

A pipeline in GenAIRR is simply a list of step instances:
```python
from GenAIRR.steps import SimulateSequence, FixVPositionAfterTrimmingIndexAmbiguity, InsertIndels
from GenAIRR.pipeline import AugmentationPipeline
from GenAIRR.mutation import S5F

pipeline = AugmentationPipeline([
    SimulateSequence(S5F(), True),
    FixVPositionAfterTrimmingIndexAmbiguity(),
    InsertIndels(0.5, 3, 0.5, 0.5)
])
```

### Executing a Pipeline

Each step takes in a `SimulationContainer` instance and modifies it based on the step's logic. Running a pipeline is straightforward:
```python
container = pipeline.execute()
print(container.get_dict())  # View the updated simulation data
```

This flow ensures that the sequence undergoes a series of predefined modifications, reflecting realistic biological processes such as mutations, indels, or trimming corrections.

### Example Steps Explained

- **SimulateSequence**: This step initializes the simulation by creating a sequence using the provided mutation model and ensuring that it meets criteria such as productivity.
- **FixVPositionAfterTrimmingIndexAmbiguity**: Corrects the position of V sequences after trimming to ensure consistency.
- **InsertIndels**: Inserts or deletes bases at random positions to simulate real sequence alterations.

## Developing Custom Steps

GenAIRR's flexibility allows users to create their own augmentation steps to match unique simulation goals. Custom steps inherit from the `AugmentationStep` base class and override the `apply` method to implement specific logic.

### Example: Creating a Custom Step

Below is an example of a custom step that reverses the sequence in the `SimulationContainer`:
```python
from GenAIRR.steps.StepBase import AugmentationStep
from GenAIRR.container.SimulationContainer import SimulationContainer

class ReverseSequenceStep(AugmentationStep):
    def apply(self, container: SimulationContainer) -> None:
        """Reverses the sequence in the SimulationContainer."""
        container.sequence = container.sequence[::-1]  # Reverse the sequence
        container.note += " Sequence reversed."
```

To use this custom step in a pipeline:
```python
pipeline = AugmentationPipeline([
    SimulateSequence(S5F(), True),
    ReverseSequenceStep()
])
```

## The Role of `SimulationContainer`

The `SimulationContainer` class is the central data structure that holds simulation results. Each step in the pipeline reads from and updates this container, ensuring continuity and coherence throughout the process.

### Key Features of `SimulationContainer`

- **Attributes**: Stores sequence data, metadata, mutation logs, and more.
- **Flexible Interaction**: Offers dictionary-like access to attributes and can be updated using dictionaries.
- **Methods for Manipulation**:
  - `add_mutation()`: Log a new mutation.
  - `shift_positions()`: Shift all sequence positions, useful after insertions or deletions.
  - `update_from_dict()`: Update the container's attributes with data from a dictionary.

### Example: Interacting with `SimulationContainer`
```python
container = SimulationContainer()
container.add_mutation(150, "A>T")
container.shift_positions(5)
print(container["mutation_rate"])  # Accessing data
container["mutation_rate"] = 0.05  # Modifying data
```

## Summary

GenAIRR's simulation flow is modular and adaptable, starting with the foundational `DataConfig` instance and building up to complex pipelines with multiple augmentation steps. Whether using predefined sequences like `HeavyChainSequence` or crafting custom simulation pipelines, users can tailor their workflows to fit specific research needs. Each augmentation step modifies the `SimulationContainer`, providing an intuitive and powerful way to simulate Ig sequences accurately and efficiently.
