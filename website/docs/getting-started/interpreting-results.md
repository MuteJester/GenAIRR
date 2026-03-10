---
title: Understanding Output
sidebar_label: Understanding Output
---

# Understanding Output

Every GenAIRR record is a dictionary containing the **absolute ground truth** of the simulated sequence. Unlike real-world aligners that must infer gene assignments from statistical models and heuristics, GenAIRR knows the exact origin of every nucleotide — no ambiguity, no probabilistic gene usage priors, no alignment scoring trade-offs.

## Record structure

```python
from GenAIRR import Experiment
from GenAIRR.ops import rate

result = Experiment.on("human_igh").mutate(rate(0.03, 0.06)).run(n=100, seed=42)
rec = result[0]

# rec is a plain dict with 47 keys
print(type(rec))       # <class 'dict'>
print(len(rec))        # 47
print(sorted(rec.keys()))
```

## Field reference

### Sequence

| Field | Type | Description |
|-------|------|-------------|
| `sequence` | str | Full nucleotide sequence (post-corruption if artifacts are enabled) |
| `sequence_length` | int | Length of `sequence` in bases |
| `germline_alignment` | str | Ungapped germline reference aligned to `sequence` |

The `germline_alignment` field encodes what happened at each position:
- **Lowercase** letters — match germline (no mutation)
- **Uppercase** letters — mutated positions (the uppercase letter is the *germline* base that was replaced)
- **`N`** — non-templated nucleotides (NP regions)

```
sequence:           gaggtgcagTtggtggagt...
germline_alignment: gaggtgcagctggtggagt...
                             ^
                    Position 9: germline 'c' was mutated to 'T'
```

### V gene

| Field | Type | Description |
|-------|------|-------------|
| `v_call` | str | V allele name |
| `v_sequence_start` | int | Start of V in `sequence` (0-based) |
| `v_sequence_end` | int | End of V in `sequence` (exclusive) |
| `v_germline_start` | int | Start position within the V germline allele |
| `v_germline_end` | int | End position within the V germline allele |
| `v_trim_5` | int | Bases trimmed from the 5' end of V |
| `v_trim_3` | int | Bases trimmed from the 3' end of V |

### D gene

| Field | Type | Description |
|-------|------|-------------|
| `d_call` | str | D allele name (empty string for light chains / chains without D) |
| `d_sequence_start` | int | Start of D in `sequence` |
| `d_sequence_end` | int | End of D in `sequence` |
| `d_germline_start` | int | Start position within the D germline allele |
| `d_germline_end` | int | End position within the D germline allele |
| `d_trim_5` | int | Bases trimmed from the 5' end of D |
| `d_trim_3` | int | Bases trimmed from the 3' end of D |

For VJ chains (light chains, TCR alpha/gamma), `d_call` is an empty string and all D coordinate fields are 0.

### J gene

| Field | Type | Description |
|-------|------|-------------|
| `j_call` | str | J allele name |
| `j_sequence_start` | int | Start of J in `sequence` |
| `j_sequence_end` | int | End of J in `sequence` |
| `j_germline_start` | int | Start position within the J germline allele |
| `j_germline_end` | int | End position within the J germline allele |
| `j_trim_5` | int | Bases trimmed from the 5' end of J |
| `j_trim_3` | int | Bases trimmed from the 3' end of J |

### Junction

| Field | Type | Description |
|-------|------|-------------|
| `junction_nt` | str | Junction nucleotide sequence (V-anchor to J-anchor inclusive) |
| `junction_aa` | str | Junction amino acid translation |
| `junction_start` | int | Junction start position in `sequence` |
| `junction_end` | int | Junction end position in `sequence` |
| `junction_length` | int | Junction length in nucleotides |

The junction spans from the conserved V-region cysteine (Cys) codon through the conserved J-region tryptophan (Trp, for heavy chains) or phenylalanine (Phe, for light chains) codon, following the IMGT definition.

### N/P regions

| Field | Type | Description |
|-------|------|-------------|
| `np1_region` | str | NP1 nucleotide sequence (between V and D, or V and J for VJ chains) |
| `np1_length` | int | NP1 length in bases |
| `np2_region` | str | NP2 nucleotide sequence (between D and J; empty for VJ chains) |
| `np2_length` | int | NP2 length in bases |

### Somatic hypermutation

| Field | Type | Description |
|-------|------|-------------|
| `mutation_rate` | float | Fraction of positions mutated (e.g. `0.054`) |
| `n_mutations` | int | Total number of point mutations |
| `mutations` | str | Comma-separated list of mutations (e.g. `13:a>T,30:g>A`) |

The `mutations` string uses the format `position:germline>mutant`, where position is 0-based into `sequence`. Lowercase germline base, uppercase mutant base.

### Functionality

| Field | Type | Description |
|-------|------|-------------|
| `productive` | bool | Whether the sequence is productive (in-frame, no stop codons in junction) |
| `vj_in_frame` | bool | Whether V and J segments are in the same reading frame |
| `stop_codon` | bool | Whether the junction contains a stop codon |
| `note` | str | Reason for non-productivity (e.g. `VJ out of frame.`) |

### Isotype

| Field | Type | Description |
|-------|------|-------------|
| `c_call` | str | Constant region allele (populated when CSR is enabled via `with_isotype_rates()`) |

### Artifact annotations

| Field | Type | Description |
|-------|------|-------------|
| `n_pcr_errors` | int | Number of PCR-introduced errors |
| `pcr_errors` | str | Comma-separated list of PCR error positions and substitutions |
| `n_sequencing_errors` | int | Number of sequencing-introduced errors |
| `sequencing_errors` | str | Comma-separated list of sequencing error positions and substitutions |
| `is_reverse_complement` | bool | Whether the sequence was reverse-complemented |
| `is_contaminant` | bool | Whether the sequence is a contaminant spike-in |

## Coordinate system

All coordinates are **0-based**. End positions are **exclusive** (Python slice convention):

```python
rec = result[0]

# Extract the V segment from the sequence
v_segment = rec["sequence"][rec["v_sequence_start"]:rec["v_sequence_end"]]

# Extract the junction
junction = rec["sequence"][rec["junction_start"]:rec["junction_end"]]
assert junction == rec["junction_nt"].replace("N", "n")  # modulo case
```

Gene segment boundaries account for trimming, N/P additions, and any 5'/3' corruption, so they always point to the exact positions in the final `sequence` string.

## VDJ vs VJ chains

Heavy chains (IGH, TCRB, TCRD) use VDJ recombination:

```
V segment ─── NP1 ─── D segment ─── NP2 ─── J segment
```

Light chains (IGK, IGL, TCRA, TCRG) use VJ recombination:

```
V segment ─── NP1 ─── J segment
```

For VJ chains, `d_call` is an empty string, `np2_region` is empty, and all D-related coordinate fields are 0.

```python
from GenAIRR import Experiment

# VDJ chain
heavy = Experiment.on("human_igh").run(n=1, seed=42)[0]
print(heavy["d_call"])      # e.g. "IGHD2-21*02"
print(heavy["np1_length"])  # e.g. 6
print(heavy["np2_length"])  # e.g. 4

# VJ chain
kappa = Experiment.on("human_igk").run(n=1, seed=42)[0]
print(kappa["d_call"])      # ""
print(kappa["np1_length"])  # e.g. 1
print(kappa["np2_length"])  # 0
```

## Next steps

- [Quick Start](/docs/getting-started/quick-start) — install and run your first simulation
- [Choosing a Config](/docs/getting-started/choosing-config) — pick the right species and chain type
- [Guides](/docs/guides/) — recipes for SHM, artifacts, export, and more
