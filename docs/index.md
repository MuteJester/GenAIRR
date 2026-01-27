# GenAIRR: Modular Ig Sequence Simulation

**GenAIRR** is a Python framework for simulating immunoglobulin (Ig) and adaptive immune receptor sequences. It provides a modular, pipeline-based architecture for generating realistic synthetic AIRR data — from naive B-cell sequences to fully mutated, sequencing-artifact-laden reads.

GenAIRR is designed for researchers and developers who need synthetic immune receptor data for benchmarking alignment tools, training machine learning models, or studying V(D)J recombination and somatic hypermutation.

```python
from GenAIRR import simulate, HUMAN_IGH_OGRDB, S5F

# Generate a single simulated heavy chain sequence
result = simulate(HUMAN_IGH_OGRDB, S5F(0.01, 0.05))
print(result.sequence)
```

---

## Key Features

**Modular Pipeline Architecture**
:   Build simulation workflows from composable steps. Each step modifies a `SimulationContainer` and can be added, removed, or reordered freely.

**Biologically Realistic Mutation Models**
:   Built-in S5F (context-dependent) and Uniform mutation models simulate somatic hypermutation at configurable rates.

**Empirical Germline Data**
:   Pre-built `DataConfig` objects contain V, D, and J allele sets, trimming distributions, and nucleotide addition patterns derived from real repertoire data (OGRDB).

**Sequencing Artifact Simulation**
:   Simulate real-world data imperfections — 5' truncation, N-base insertions, insertions/deletions, and read-length limits.

**Ambiguity Resolution**
:   Automatic correction steps resolve positional ambiguities introduced by trimming, ensuring accurate ground-truth annotations.

**Reproducibility**
:   Seed management (`set_seed`, `get_seed`, `reset_seed`) enables deterministic sequence generation.

---

## Supported Chains

| Config | Chain | Species |
|--------|-------|---------|
| `HUMAN_IGH_OGRDB` | Heavy (IGH) | Human |
| `HUMAN_IGK_OGRDB` | Kappa light (IGK) | Human |
| `HUMAN_IGL_OGRDB` | Lambda light (IGL) | Human |
| `HUMAN_TRB_OGRDB` | T-cell receptor beta (TRB) | Human |

---

## Installation

```bash
pip install GenAIRR
```

Requires Python 3.9+.

---

## Quick Example: Full Pipeline

For complete control over the simulation process, build a pipeline with explicit steps:

```python
from GenAIRR import Pipeline, steps, HUMAN_IGH_OGRDB, S5F

pipeline = Pipeline(
    config=HUMAN_IGH_OGRDB,
    steps=[
        # 1. Generate a mutated sequence
        steps.SimulateSequence(
            S5F(min_mutation_rate=0.01, max_mutation_rate=0.05),
            productive=True
        ),
        # 2. Resolve positional ambiguities
        steps.FixVPositionAfterTrimmingIndexAmbiguity(),
        steps.FixDPositionAfterTrimmingIndexAmbiguity(),
        steps.FixJPositionAfterTrimmingIndexAmbiguity(),
        # 3. Biological corrections
        steps.CorrectForVEndCut(),
        steps.CorrectForDTrims(),
        # 4. Record mutation rate
        steps.DistillMutationRate(),
        # 5. Simulate sequencing artifacts
        steps.CorruptSequenceBeginning(),
        steps.EnforceSequenceLength(),
        steps.InsertNs(),
        # 6. Quality variants
        steps.ShortDValidation(),
        steps.InsertIndels(),
    ]
)

result = pipeline.execute()
data = result.get_dict()
print(data['v_call'], data['mutation_rate'])
```

---

## Documentation Overview

### Getting Started
- **[Installation & Quick Start](getting_started.md)** — Install GenAIRR and run your first simulation
- **[Step-by-Step Tutorial](step_by_step_tutorial.md)** — Build a pipeline from scratch with explanations at each stage

### User Guide
- **[Biological Context](biological_context.md)** — The immunobiology behind GenAIRR's simulation model
- **[How the Pipeline Works](genairr_flow.md)** — Architecture: DataConfig, Pipeline, Steps, and SimulationContainer
- **[Best Practices](best_practices.md)** — Guidelines for realistic and reproducible simulations

### Tutorials (Jupyter Notebooks)
- **[Quick Start Guide](tutorials/Quick Start Guide.ipynb)** — Interactive introduction
- **[Advanced Custom Generation](tutorials/Advanced Custom Generation.ipynb)** — Specific allele selection and chain-type control
- **[Introduction to DataConfig](tutorials/Introduction to the DataConfig Object.ipynb)** — Explore allele sets and distributions
- **[Creating Custom DataConfig](tutorials/Creating Custom DataConfig from FASTA Files.ipynb)** — Build configs from your own FASTA files
- **[Uniform Allele Distribution](tutorials/Generating_Uniform_Allele_Distribution_Dataset.ipynb)** — Generate balanced benchmarking datasets

### Reference
- **[API Reference](api_reference.md)** — Classes, methods, and function signatures
- **[Parameter Reference](parameter_reference.md)** — Detailed parameter descriptions for every step

### Advanced
- **[Custom DataConfig](custom_data_config.md)** — Create DataConfig instances from your own germline databases

### Support
- **[FAQ](faq.md)** — Frequently asked questions
- **[Troubleshooting](troubleshooting.md)** — Common issues and solutions
- **[Migration Guide](migration_guide.md)** — Upgrading from previous versions

---

## Citation

If you use GenAIRR in your research, please cite:

> GenAIRR — *Briefings in Bioinformatics*, 2024.
> [DOI: 10.1093/bib/bbae556](https://academic.oup.com/bib/article/25/6/bbae556/7863770)

## Links

- [GitHub Repository](https://github.com/MuteJester/GenAIRR)
- [Issue Tracker](https://github.com/MuteJester/GenAIRR/issues)
- [PyPI Package](https://pypi.org/project/GenAIRR/)
