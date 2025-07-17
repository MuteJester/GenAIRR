<!-- -------------------------------------------------------- -->
<!-- Banner / Logo (optional) -->
<!-- Replace the link below with your own logo or SVG banner -->
<p align="center">
  <img src="https://raw.githubusercontent.com/your-org/GenAIRR/main/docs/_static/banner.svg" alt="GenAIRR" width="60%" />
</p>

<h1 align="center">GenAIRR</h1>

<p align="center">
  <b>Adaptive Immune Receptor Repertoire sequence simulator</b><br/>
  Generate realistic BCR & TCR repertoires in a single line of Python.
</p>

<p align="center">
  <a href="https://pypi.org/project/GenAIRR/"><img src="https://img.shields.io/pypi/v/GenAIRR.svg?logo=pypi&logoColor=white" alt="PyPI version"></a>
  <a href="https://genairr.readthedocs.io/en/latest/"><img src="https://img.shields.io/readthedocs/genairr?logo=readthedocs"></a>
</a>
</p>

---

## 📑 Table of Contents
1. [Why GenAIRR?](#-why-genairr)
2. [Key Features](#-key-features)
3. [Installation](#-installation)
4. [Quick Start](#-quick-start)
5. [Examples](#-examples)
6. [Mutation Models](#-mutation-models)
7. [Roadmap](#-roadmap)
8. [Contributing](#-contributing)
9. [Citing GenAIRR](#-citing-genairr)
10. [License](#-license)
11. [Acknowledgements](#-acknowledgements)

---

## 🧐 Why GenAIRR?
<details>
<summary>Click to expand</summary>

*Benchmarking modern aligners, exploring somatic-hypermutation, or stress-testing novel ML pipelines requires large, perfectly-annotated repertoires—not snippets of real data peppered with sequencing error.*  
GenAIRR fills that gap with a **plug-and-play, fully-extensible simulation engine** that produces sequences while giving you full ground-truth labels.

</details>

---

## ✨ Key Features
| Category | Highlights                                                                         |
| -------- |------------------------------------------------------------------------------------|
| **Realistic Simulation** | Context-aware S5F mutations, indels, allele-specific trimming, NP-region modelling |
| **Composable Pipelines** | Chain together built-in & custom `AugmentationStep`s into simulation pipelines     |
| **Multi-Chain Support** | Heavy & light BCRs plus TCR-β out of the box                                       |
| **Research-ready Output** | JSON / pandas export, built-in plotting stubs, deterministic seeds                 |
| **Docs & Tutorials** | Rich API docs, Jupyter notebooks, step-by-step guides                              |

---

## ⚡ Installation
```bash
# Python ≥ 3.9
pip install GenAIRR
# or the bleeding edge
pip install git+https://github.com/MuteJester/GenAIRR.git
````

---

## 🚀 Quick Start

Below is a 60-second tour. See [`/examples`](examples/) for notebooks and CLI usages.

```python
from GenAIRR.pipeline import AugmentationPipeline
from GenAIRR.steps import SimulateSequence, FixVPositionAfterTrimmingIndexAmbiguity
from GenAIRR.mutation import S5F
from GenAIRR.data import HUMAN_IGH_OGRDB
from GenAIRR.steps.StepBase import AugmentationStep

# 1️⃣  Configure built-in germline data
AugmentationStep.set_dataconfig(HUMAN_IGH_OGRDB)

# 2️⃣  Build a minimal pipeline
pipeline = AugmentationPipeline([
    SimulateSequence(S5F(min_mutation_rate=0.003, max_mutation_rate=0.25), True),
    FixVPositionAfterTrimmingIndexAmbiguity()
])

# 3️⃣  Simulate!
sim = pipeline.execute()
print(sim.get_dict())
```

---

## 🧑‍💻 Examples

### 1. Full Heavy-Chain Pipeline

```python
from GenAIRR.steps import (
    FixDPositionAfterTrimmingIndexAmbiguity, FixJPositionAfterTrimmingIndexAmbiguity,
    CorrectForVEndCut, CorrectForDTrims, CorruptSequenceBeginning,
    InsertNs, InsertIndels, ShortDValidation, DistillMutationRate
)

pipeline = AugmentationPipeline([
    SimulateSequence(S5F(min_mutation_rate=0.003, max_mutation_rate=0.25), True),
    FixVPositionAfterTrimmingIndexAmbiguity(),
    FixDPositionAfterTrimmingIndexAmbiguity(),
    FixJPositionAfterTrimmingIndexAmbiguity(),
    CorrectForVEndCut(),
    CorrectForDTrims(),
    CorruptSequenceBeginning(0.7, [0.4, 0.4, 0.2], 576, 210, 310, 50),
    InsertNs(0.02, 0.5),
    ShortDValidation(),
    InsertIndels(0.5, 5, 0.5, 0.5),
    DistillMutationRate()
])
result = pipeline.execute()
```

### 2. Naïve Sequence (no SHM)

```python
from GenAIRR.mutation import Uniform
naive_step = SimulateSequence(Uniform(0, 0), True)
pipeline = AugmentationPipeline([naive_step])
naive_seq = pipeline.execute()
print(naive_seq.sequence)
```

### 3. Custom Allele Combination

```python
custom_step = SimulateSequence(
    S5F(0.003, 0.25),
    True,
    specific_v=HUMAN_IGH_OGRDB.v_alleles['IGHV1-2*02'][0],  # specific V allele
    specific_d=HUMAN_IGH_OGRDB.d_alleles['IGHD3-10*01'][0], # specific D allele  
    specific_j=HUMAN_IGH_OGRDB.j_alleles['IGHJ4*02'][0]     # specific J allele
)
pipeline = AugmentationPipeline([custom_step])
print(pipeline.execute().get_dict())
```
---

## 🔬 Mutation Models

| Model            | Description                             | When to use                   |
| ---------------- | --------------------------------------- | ----------------------------- |
| `S5F`            | Context-specific somatic hyper-mutation | Antibody maturation studies   |
| `Uniform`        | Evenly random mutations                 | Baselines / ablation          |
| **Your Model ➕** | Implement `BaseMutationModel`           | Custom evolutionary scenarios |

```python
from GenAIRR.mutation import S5F
s5f = S5F(min_mutation_rate=0.01, max_mutation_rate=0.05)
mut_seq, muts, rate = s5f.apply_mutation(naive_seq)
```

---

## 🗺️ Roadmap

* [ ] 🚧 **More Complex Mutation Model (With Selection)**
* [ ] 🚧 **More Built-in Data Configs** (e.g., TCR, custom germlines)
* [ ] 🚧 **More Built-in Steps** (e.g., more mutation models, more data augmentation)
* [ ] 🚧 **Deeper Docs** (e.g., more examples, more tutorials)

*See the [open issues](https://github.com/your-org/GenAIRR/issues).*
  Feel something’s missing? [Open a feature request](https://github.com/your-org/GenAIRR/issues/new).

---

## 🤝 Contributing

Contributions are welcome! 💙
Please read our [contributing guide](CONTRIBUTING.md) and check the **good first issue** label.

---

## ✏️ Citing GenAIRR

If GenAIRR helps your research, please cite:

```
Konstantinovsky T, Peres A, Polak P, Yaari G.  
An unbiased comparison of immunoglobulin sequence aligners.
Briefings in Bioinformatics. 2024 Sep 23; 25(6): bbae556.  
https://doi.org/10.1093/bib/bbae556  
PMID: 39489605 | PMCID: PMC11531861
```

---

## 📜 License

Distributed under the GPL3 License. See **[LICENSE](LICENSE)** for details.

---

## 🙏 Acknowledgements

GenAIRR is inspired by and builds upon amazing work from the immunoinformatics community—especially [AIRRship](https://github.com/Cowanlab/airrship).

<!-- End of README -->


