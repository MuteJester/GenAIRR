<h1 align="center">GenAIRR</h1>

<p align="center">
  <b>Synthetic Adaptive Immune Receptor Repertoire Generator</b>
</p>

<p align="center">
  <a href="https://pypi.org/project/GenAIRR/"><img src="https://img.shields.io/pypi/v/GenAIRR.svg?logo=pypi&logoColor=white" alt="PyPI"></a>
  <a href="https://github.com/MuteJester/GenAIRR/actions/workflows/test.yml"><img src="https://github.com/MuteJester/GenAIRR/actions/workflows/test.yml/badge.svg" alt="Tests"></a>
  <a href="https://pypi.org/project/GenAIRR/"><img src="https://img.shields.io/pypi/pyversions/GenAIRR.svg?logo=python&logoColor=white" alt="Python"></a>
  <a href="https://github.com/MuteJester/GenAIRR/blob/master/LICENSE"><img src="https://img.shields.io/github/license/MuteJester/GenAIRR" alt="License"></a>
</p>

<p align="center">
  High-performance BCR and TCR sequence simulation with full ground-truth annotations.<br/>
  C engine &middot; 23 species &middot; zero mandatory dependencies &middot; cross-platform wheels
</p>

<p align="center">
  <a href="https://mutejester.github.io/GenAIRR/"><b>📖 Documentation</b></a>
</p>

---

## Installation

```bash
pip install GenAIRR
```

Pre-built wheels are available for **Linux**, **macOS**, and **Windows** (Python 3.9+). No compiler required.

---

## Quick Start

```python
from GenAIRR import Experiment
from GenAIRR.ops import rate

# Generate 1,000 mutated human heavy-chain sequences
result = Experiment.on("human_igh").mutate(rate(0.02, 0.08)).run(n=1000, seed=42)

# Each record is an AIRR-format dict with full ground truth
rec = result[0]
rec["sequence"]      # full nucleotide sequence
rec["v_call"]        # e.g. "IGHVF10-G50*04"
rec["d_call"]        # e.g. "IGHD2-21*02"
rec["j_call"]        # e.g. "IGHJ4*02"
rec["mutation_rate"]  # e.g. 0.054
rec["productive"]     # True / False
```

### Realistic Sequencing Experiment

```python
from GenAIRR import Experiment
from GenAIRR.ops import (
    with_d_inversion, with_receptor_revision,
    rate, model, with_isotype_rates, with_antigen_selection,
    with_primer_mask, with_umi, with_pcr,
    with_5prime_loss, with_3prime_loss, with_quality_profile,
    with_indels, with_ns,
)

result = (
    Experiment.on("human_igh")

    # V(D)J recombination with biological events
    .recombine(
        with_d_inversion(0.15),
        with_receptor_revision(0.05),
    )

    # Somatic hypermutation with CSR and selection pressure
    .mutate(
        model("s5f"),
        rate(0.01, 0.05),
        with_isotype_rates(),
        with_antigen_selection(0.5),
    )

    # Library preparation
    .prepare(
        with_primer_mask(),
        with_umi(12),
        with_pcr(error_rate=1e-4, cycles=30),
    )

    # Sequencing artifacts
    .sequence(
        with_5prime_loss(min_remove=5, max_remove=30),
        with_3prime_loss(min_remove=5, max_remove=20),
        with_quality_profile(base=0.001, peak=0.02),
    )

    # Post-sequencing noise
    .observe(
        with_indels(prob=0.005),
        with_ns(prob=0.005),
    )

    .run(n=1000, seed=42)
)
```

### Streaming (Memory-Efficient)

```python
from GenAIRR import Experiment
from GenAIRR.ops import rate

sim = Experiment.on("human_igh").mutate(rate(0.05, 0.15)).compile(seed=42)

for record in sim.stream():
    print(record["v_call"])  # one dict at a time, no accumulation
    break                    # infinite iterator — break when done
```

### Export

```python
result.to_csv("repertoire.tsv")          # AIRR TSV
result.to_fasta("repertoire.fasta")      # FASTA
df = result.to_dataframe()               # pandas DataFrame
```

---

## Output Fields

Every record is a dictionary containing the **absolute ground truth** of the
simulated sequence. Unlike real-world aligners that must infer gene assignments
from statistical models and heuristics, GenAIRR knows the exact origin of every
nucleotide. The metadata reflects what a hypothetical perfect aligner would
report if it could examine the sequence with full knowledge of the rearrangement
process &mdash; no ambiguity, no probabilistic gene usage priors, no alignment
scoring trade-offs.

| Field | Type | Description |
|-------|------|-------------|
| `sequence` | str | Full nucleotide sequence (post-corruption if artifacts are enabled) |
| `sequence_length` | int | Length of `sequence` in bases |
| `germline_alignment` | str | Ungapped germline reference aligned to `sequence` (lowercase = germline, uppercase = mutated, `N` = non-templated) |
| **V gene** | | |
| `v_call` | str | V allele name |
| `v_call_true` | str | Simulator ground-truth sampled V allele |
| `v_sequence_start` | int | Start of the V segment in `sequence` (0-based) |
| `v_sequence_end` | int | End of the V segment in `sequence` (exclusive) |
| `v_germline_start` | int | Start position within the V germline allele |
| `v_germline_end` | int | End position within the V germline allele |
| `v_trim_5` | int | Bases trimmed from the 5' end of V |
| `v_trim_3` | int | Bases trimmed from the 3' end of V |
| **D gene** | | |
| `d_call` | str | D allele name (empty for light chains / chains without D) |
| `d_call_true` | str | Simulator ground-truth sampled D allele |
| `d_sequence_start` | int | Start of the D segment in `sequence` |
| `d_sequence_end` | int | End of the D segment in `sequence` |
| `d_germline_start` | int | Start position within the D germline allele |
| `d_germline_end` | int | End position within the D germline allele |
| `d_trim_5` | int | Bases trimmed from the 5' end of D |
| `d_trim_3` | int | Bases trimmed from the 3' end of D |
| **J gene** | | |
| `j_call` | str | J allele name |
| `j_call_true` | str | Simulator ground-truth sampled J allele |
| `j_sequence_start` | int | Start of the J segment in `sequence` |
| `j_sequence_end` | int | End of the J segment in `sequence` |
| `j_germline_start` | int | Start position within the J germline allele |
| `j_germline_end` | int | End position within the J germline allele |
| `j_trim_5` | int | Bases trimmed from the 5' end of J |
| `j_trim_3` | int | Bases trimmed from the 3' end of J |
| **Junction** | | |
| `junction_nt` | str | Junction nucleotide sequence (V-anchor to J-anchor inclusive) |
| `junction_aa` | str | Junction amino acid translation |
| `junction_start` | int | Junction start position in `sequence` |
| `junction_end` | int | Junction end position in `sequence` |
| `junction_length` | int | Junction length in nucleotides |
| **N/P regions** | | |
| `np1_region` | str | NP1 nucleotide sequence (between V and D, or V and J) |
| `np1_length` | int | NP1 length in bases |
| `np2_region` | str | NP2 nucleotide sequence (between D and J) |
| `np2_length` | int | NP2 length in bases |
| **Somatic hypermutation** | | |
| `mutation_rate` | float | Fraction of positions mutated |
| `n_mutations` | int | Total number of point mutations |
| `mutations` | str | Comma-separated list of mutations (e.g. `13:a>T,30:g>A`) |
| **Functionality** | | |
| `productive` | bool | Whether the sequence is productive (in-frame, no stop codons) |
| `vj_in_frame` | bool | Whether V and J segments are in the same reading frame |
| `stop_codon` | bool | Whether the junction contains a stop codon |
| `note` | str | Reason for non-productivity (e.g. `VJ out of frame.`) |
| **Isotype** | | |
| `c_call` | str | Constant region allele (when CSR is enabled) |
| **Artifact annotations** | | |
| `n_pcr_errors` | int | Number of PCR-introduced errors |
| `pcr_errors` | str | Comma-separated list of PCR error positions and substitutions |
| `n_sequencing_errors` | int | Number of sequencing-introduced errors |
| `sequencing_errors` | str | Comma-separated list of sequencing error positions and substitutions |
| `is_reverse_complement` | bool | Whether the sequence was reverse-complemented |
| `is_contaminant` | bool | Whether the sequence is a contaminant spike-in |
| `d_inverted` | bool | Whether D inversion (reverse-complement D event) was applied |
| `receptor_revised` | bool | Whether receptor revision was applied |
| `revision_footprint_length` | int | Preserved old-V footprint length during receptor revision |
| `original_v_allele_name` | str | V allele before receptor revision |

All coordinates are 0-based. Gene segment boundaries account for trimming,
N/P additions, and any 5'/3' corruption, so they point to the exact positions
in the final `sequence` string.

---

## Supported Species & Chains

GenAIRR ships with **106 built-in configurations** covering 23 species (sourced from IMGT and OGRDB).

```python
from GenAIRR import list_configs
print(list_configs())  # all available configs
```

| Species | BCR | TCR |
|---------|-----|-----|
| Human | IGH, IGK, IGL | TCRA, TCRB, TCRD, TCRG |
| Mouse | IGH, IGK, IGL | TCRA, TCRB, TCRD, TCRG |
| Rat | IGH, IGK, IGL | &mdash; |
| Rabbit | IGH, IGK, IGL | TCRA, TCRB, TCRD, TCRG |
| Dog | IGH, IGK, IGL | TCRA, TCRB, TCRD, TCRG |
| Cat | IGK, IGL | TCRA, TCRB, TCRD, TCRG |
| Rhesus | IGH, IGK, IGL | TCRA, TCRB, TCRD, TCRG |

<details>
<summary>All 23 species</summary>

Alpaca, Cat, Chicken, Cow, Cynomolgus, Dog, Dromedary, Ferret, Goat, Gorilla,
Horse, Human, Mouse (generic + C57BL/6J), Pig, Platypus, Rabbit, Rat, Rhesus,
Salmon, Sheep, Trout, Zebrafish.

</details>

```python
Experiment.on("mouse_igh").run(n=500)
Experiment.on("rabbit_tcrb").run(n=500)
Experiment.on("rhesus_igk").run(n=500)
```

---

## Key Features

- **C simulation engine** &mdash; 15,000&ndash;30,000 sequences/second end-to-end on a single core
- **Context-aware S5F somatic hypermutation** &mdash; 5-mer motif-based targeting with empirical substitution profiles
- **Full AIRR-format output** &mdash; V/D/J calls, germline alignments, junction boundaries, mutation annotations
- **Sequencing artifact simulation** &mdash; 5'/3' degradation, indels, N-insertions, quality errors, PCR artifacts
- **Biological processes** &mdash; D-gene inversion, receptor revision, class switch recombination, selection pressure
- **Allele locking** &mdash; constrain simulation to specific V/D/J alleles
- **Deterministic seeds** &mdash; fully reproducible results across runs and platforms
- **Zero mandatory dependencies** &mdash; optional extras for DataConfig building, visualization, and MCP

---

## Reproducibility

```python
# Pass a seed for deterministic, reproducible results
result = Experiment.on("human_igh").run(n=1000, seed=42)
```

---

## Optional Extras

```bash
pip install GenAIRR[all]          # numpy, scipy, graphviz, tqdm, fastmcp
pip install GenAIRR[dataconfig]   # numpy + scipy (custom DataConfig building)
pip install GenAIRR[mcp]          # FastMCP server for AI-assisted analysis
```

---

## Citing GenAIRR

If GenAIRR is useful in your research, please cite:

> Konstantinovsky T, Peres A, Polak P, Yaari G. An unbiased comparison of immunoglobulin sequence aligners. *Briefings in Bioinformatics*. 2024;25(6):bbae556. [doi:10.1093/bib/bbae556](https://doi.org/10.1093/bib/bbae556)

---

## Contributing

Contributions are welcome. See [CONTRIBUTING.md](CONTRIBUTING.md) for development setup and guidelines.

## License

GPL-3.0. See [LICENSE](LICENSE).
