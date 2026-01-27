<h1 align="center">GenAIRR</h1>

<p align="center">
  <b>Adaptive Immune Receptor Repertoire Sequence Simulator</b><br/>
  Generate realistic BCR &amp; TCR repertoires with full ground-truth annotations in Python.
</p>

<p align="center">
  <a href="https://pypi.org/project/GenAIRR/"><img src="https://img.shields.io/pypi/v/GenAIRR.svg?logo=pypi&logoColor=white" alt="PyPI version"></a>
  <a href="https://genairr.readthedocs.io/en/latest/"><img src="https://img.shields.io/readthedocs/genairr?logo=readthedocs" alt="Docs"></a>
  <a href="https://github.com/MuteJester/GenAIRR/blob/master/LICENSE"><img src="https://img.shields.io/github/license/MuteJester/GenAIRR" alt="License"></a>
</p>

---

## Why GenAIRR?

Benchmarking sequence aligners, studying somatic hypermutation, or training ML models on immune repertoires requires large, perfectly-annotated datasets — not noisy snippets of real sequencing data.

GenAIRR is a **plug-and-play, fully-extensible simulation engine** that produces realistic immunoglobulin and TCR sequences while giving you complete ground-truth labels for every position, mutation, and gene segment.

---

## Key Features

| Category | Highlights |
| -------- | ---------- |
| **Realistic Simulation** | Context-aware S5F mutations, indels, allele-specific trimming, NP-region modelling |
| **Composable Pipelines** | Chain together built-in & custom steps into simulation pipelines |
| **Multi-Chain Support** | Heavy chain, kappa/lambda light chains, and TCR-beta out of the box |
| **Research-ready Output** | Full ground-truth annotations, JSON/pandas export, deterministic seeds |
| **Docs & Tutorials** | Step-by-step guides, Jupyter notebooks, API reference |

---

## Installation

```bash
# Python >= 3.9
pip install GenAIRR
```

---

## Quick Start

### One-liner

```python
from GenAIRR import simulate, HUMAN_IGH_OGRDB, S5F

result = simulate(HUMAN_IGH_OGRDB, S5F(0.003, 0.25))
print(result.sequence)
```
```
CAGGTGCAGCTGCAGGAGTCGGGCCCAGGACTGGTGAAGCCTTCGGGGACCCTGTCCCTCACCTGCGCTG...
```

Generate multiple sequences at once:
```python
results = simulate(HUMAN_IGH_OGRDB, S5F(0.003, 0.25), n=100)
```

### Pipeline (Full Control)

For complete control over the simulation, use the Pipeline API:

```python
from GenAIRR import Pipeline, steps, HUMAN_IGH_OGRDB, S5F

pipeline = Pipeline(
    config=HUMAN_IGH_OGRDB,
    steps=[
        steps.SimulateSequence(S5F(min_mutation_rate=0.003, max_mutation_rate=0.25), productive=True),
        steps.FixVPositionAfterTrimmingIndexAmbiguity(),
        steps.FixDPositionAfterTrimmingIndexAmbiguity(),
        steps.FixJPositionAfterTrimmingIndexAmbiguity(),
        steps.CorrectForVEndCut(),
        steps.CorrectForDTrims(),
        steps.DistillMutationRate(),
    ]
)

sim = pipeline.execute()
print(sim.get_dict())
```
```python
{
    'sequence': 'CAGGTGCAGCTGCAGGAGTCGGGCCCAGGACTGGTGAAGCCTTCG...',
    'v_call': ['IGHVF3-G8*04'],
    'd_call': ['IGHD6-6*01'],
    'j_call': ['IGHJ4*02'],
    'productive': True,
    'mutation_rate': 0.0027,
    'mutations': {142: 'T>C'},
    'v_sequence_start': 0,
    'v_sequence_end': 293,
    'd_sequence_start': 298,
    'd_sequence_end': 316,
    'j_sequence_start': 323,
    'j_sequence_end': 367,
    # ... and more fields
}
```

Every output includes the full sequence, V/D/J gene calls, mutation positions, region boundaries, and quality metrics — ready for downstream analysis.

---

## Examples

### Full Heavy-Chain Pipeline

A production-ready pipeline that simulates sequences with biological corrections and sequencing artifacts:

```python
from GenAIRR import Pipeline, steps, HUMAN_IGH_OGRDB, S5F

pipeline = Pipeline(
    config=HUMAN_IGH_OGRDB,
    steps=[
        # Core: generate sequence with somatic hypermutation
        steps.SimulateSequence(S5F(min_mutation_rate=0.003, max_mutation_rate=0.25), productive=True),

        # Correct ground-truth positions after trimming ambiguities
        steps.FixVPositionAfterTrimmingIndexAmbiguity(),
        steps.FixDPositionAfterTrimmingIndexAmbiguity(),
        steps.FixJPositionAfterTrimmingIndexAmbiguity(),
        steps.CorrectForVEndCut(),
        steps.CorrectForDTrims(),

        # Calculate final mutation rate
        steps.DistillMutationRate(),

        # Simulate sequencing artifacts
        steps.CorruptSequenceBeginning(),   # 5' end degradation
        steps.EnforceSequenceLength(),      # read-length limit
        steps.InsertNs(),                   # ambiguous base calls
        steps.ShortDValidation(),           # D-region QC
        steps.InsertIndels(),               # sequencing indels
    ]
)

result = pipeline.execute()
```

### Naive Sequence (No Mutations)

```python
from GenAIRR import Pipeline, steps, HUMAN_IGH_OGRDB, Uniform

pipeline = Pipeline(
    config=HUMAN_IGH_OGRDB,
    steps=[steps.SimulateSequence(Uniform(0, 0), productive=True)]
)
naive_seq = pipeline.execute()
```

### Light Chain

```python
from GenAIRR import Pipeline, steps, HUMAN_IGK_OGRDB, S5F

pipeline = Pipeline(
    config=HUMAN_IGK_OGRDB,  # kappa light chain (no D segment)
    steps=[
        steps.SimulateSequence(S5F(min_mutation_rate=0.02, max_mutation_rate=0.08), productive=True),
        steps.FixVPositionAfterTrimmingIndexAmbiguity(),
        steps.FixJPositionAfterTrimmingIndexAmbiguity(),
        steps.CorrectForVEndCut(),
        steps.DistillMutationRate(),
    ]
)
```

### Custom Allele Combination

```python
from GenAIRR import Pipeline, steps, HUMAN_IGH_OGRDB, S5F

pipeline = Pipeline(
    config=HUMAN_IGH_OGRDB,
    steps=[
        steps.SimulateSequence(
            S5F(0.003, 0.25),
            productive=True,
            specific_v=HUMAN_IGH_OGRDB.v_alleles['IGHVF1-G1'][0],
            specific_d=HUMAN_IGH_OGRDB.d_alleles['IGHD1-1'][0],
            specific_j=HUMAN_IGH_OGRDB.j_alleles['IGHJ1'][0]
        )
    ]
)
```

### Batch Generation

```python
import pandas as pd
from GenAIRR import Pipeline, steps, HUMAN_IGH_OGRDB, S5F

pipeline = Pipeline(
    config=HUMAN_IGH_OGRDB,
    steps=[
        steps.SimulateSequence(S5F(min_mutation_rate=0.003, max_mutation_rate=0.25), productive=True),
        steps.FixVPositionAfterTrimmingIndexAmbiguity(),
        steps.FixDPositionAfterTrimmingIndexAmbiguity(),
        steps.FixJPositionAfterTrimmingIndexAmbiguity(),
        steps.CorrectForVEndCut(),
        steps.CorrectForDTrims(),
        steps.DistillMutationRate(),
    ]
)

# Generate 1000 sequences as a DataFrame
df = pd.DataFrame([pipeline.execute().get_dict() for _ in range(1000)])
df.to_csv('simulated_repertoire.csv', index=False)
```

---

## Mutation Models

| Model | Description | When to use |
| ----- | ----------- | ----------- |
| `S5F` | Context-dependent somatic hypermutation based on empirical 5-mer frequencies | Realistic antibody maturation studies |
| `Uniform` | Uniform random mutations | Baselines, ablation experiments |
| **Custom** | Implement `BaseMutationModel` | Your own evolutionary scenarios |

```python
from GenAIRR import S5F, Uniform

# Realistic context-aware SHM
s5f = S5F(min_mutation_rate=0.01, max_mutation_rate=0.05)

# Simple uniform mutations
uniform = Uniform(min_mutation_rate=0.01, max_mutation_rate=0.05)
```

---

## Available Data Configurations

| Config | Chain | Source |
| ------ | ----- | ------ |
| `HUMAN_IGH_OGRDB` | Heavy chain (BCR) | OGRDB |
| `HUMAN_IGH_EXTENDED` | Heavy chain extended | OGRDB |
| `HUMAN_IGK_OGRDB` | Kappa light chain | OGRDB |
| `HUMAN_IGL_OGRDB` | Lambda light chain | OGRDB |
| `HUMAN_TCRB_IMGT` | TCR-beta | IMGT |

---

## Reproducibility

```python
from GenAIRR import set_seed, get_seed, reset_seed

set_seed(42)         # deterministic results
print(get_seed())    # check current seed
reset_seed()         # back to random
```

---

## Documentation

- **[Getting Started](https://genairr.readthedocs.io/en/latest/getting_started.html)** — Overview and first pipeline
- **[Step-by-Step Tutorial](https://genairr.readthedocs.io/en/latest/step_by_step_tutorial.html)** — Build a pipeline from scratch
- **[API Reference](https://genairr.readthedocs.io/en/latest/api_reference.html)** — All classes, parameters, and defaults
- **[Migration Guide](https://genairr.readthedocs.io/en/latest/migration_guide.html)** — Upgrading from older versions
- **[Biological Context](https://genairr.readthedocs.io/en/latest/biological_context.html)** — What biological processes are simulated

---

## Roadmap

- [ ] Selection-aware mutation model
- [ ] Additional germline databases
- [ ] Sphinx auto-generated API docs from docstrings

*See [open issues](https://github.com/MuteJester/GenAIRR/issues).* Feel something's missing? [Open a feature request](https://github.com/MuteJester/GenAIRR/issues/new).

---

## Contributing

Contributions are welcome! Please read our [contributing guide](CONTRIBUTING.md) and check the **good first issue** label.

---

## Citing GenAIRR

If GenAIRR helps your research, please cite:

```
Konstantinovsky T, Peres A, Polak P, Yaari G.
An unbiased comparison of immunoglobulin sequence aligners.
Briefings in Bioinformatics. 2024 Sep 23; 25(6): bbae556.
https://doi.org/10.1093/bib/bbae556
PMID: 39489605 | PMCID: PMC11531861
```

---

## License

Distributed under the GPL-3.0 License. See **[LICENSE](LICENSE)** for details.

---

## Acknowledgements

GenAIRR is inspired by and builds upon work from the immunoinformatics community — especially [AIRRship](https://github.com/Cowanlab/airrship).
