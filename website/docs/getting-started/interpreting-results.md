---
title: Understanding Output
sidebar_label: Understanding Output
---

# Understanding Output

Every GenAIRR record is a dictionary containing the **absolute ground truth** of the simulated sequence. Unlike real-world aligners that must infer gene assignments from statistical models and heuristics, GenAIRR knows the exact origin of every nucleotide — no ambiguity, no probabilistic gene usage priors, no alignment scoring trade-offs.

The new Rust simulation engine produces approximately **70 ground-truth fields** for every sequence, following the AIRR Rearrangement schema while adding simulator-specific annotations.

## Record structure

```python
import GenAIRR as ga

result = ga.Experiment.on("human_igh").recombine().mutate(count=(5, 15)).run(n=100, seed=42)
rec = result[0]

# rec is a plain dict with ~70 keys
print(type(rec))       # <class 'dict'>
print(len(rec))        # 69 (or 72 if expose_provenance=True)
```

## Field reference

### AIRR Metadata & Alignments

| Field | Type | Description |
|-------|------|-------------|
| `sequence_id` | str | Unique identifier for the sequence (e.g. "seq0") |
| `sequence` | str | Full nucleotide sequence (post-mutation and post-corruption) |
| `sequence_aa` | str | Full amino acid translation of the `sequence` |
| `sequence_alignment` | str | The simulated sequence aligned to the germline (includes gaps `-`) |
| `germline_alignment` | str | The germline reference sequence aligned to the simulated sequence |
| `germline_alignment_d_mask` | str | Same as `germline_alignment` but with the D-segment masked by `N`s |
| `sequence_length` | int | Length of `sequence` in bases |
| `rev_comp` | bool | Whether the sequence was reverse-complemented (technical artifact) |
| `locus` | str | The genetic locus (e.g., "IGH", "TRB") |

All alignment strings are returned in **uppercase** for maximum compatibility with AIRR analysis tools.

### V, D, and J Segments

Each segment (V, D, and J) has a consistent set of annotation fields. Replace `[X]` with `v`, `d`, or `j`:

| Field | Type | Description |
|-------|------|-------------|
| `[X]_call` | str | Allele name (e.g. "IGHV3-23*01") |
| `[X]_cigar` | str | CIGAR string representing the alignment (e.g. "301M") |
| `[X]_identity` | float | Fraction of matching bases in the alignment |
| `[X]_sequence_start` | int | Start of the segment in `sequence` (0-based) |
| `[X]_sequence_end` | int | End of the segment in `sequence` (exclusive) |
| `[X]_alignment_start` | int | Start of the segment in `sequence_alignment` |
| `[X]_alignment_end` | int | End of the segment in `sequence_alignment` |
| `[X]_germline_start` | int | Start position within the germline allele |
| `[X]_germline_end` | int | End position within the germline allele |
| `[X]_trim_5` | int | Bases trimmed from the 5' end during recombination |
| `[X]_trim_3` | int | Bases trimmed from the 3' end during recombination |
| `[X]_score` | float | Alignment score (fixed to 1.0 for ground truth) |
| `[X]_support` | float | Statistical support (fixed to 1.0 for ground truth) |

**Note:** For VJ chains (light chains and TCR alpha/gamma), all `d_*` fields are empty or `None`.

### Junction & NP Regions

| Field | Type | Description |
|-------|------|-------------|
| `junction` | str | Junction nucleotide sequence (V-anchor to J-anchor inclusive) |
| `junction_aa` | str | Junction amino acid translation |
| `junction_start` | int | Junction start position in `sequence` |
| `junction_end` | int | Junction end position in `sequence` |
| `junction_length` | int | Junction length in nucleotides |
| `np1` | str | Nucleotide sequence of the NP1 region (V-D or V-J) |
| `np1_aa` | str | Amino acid translation of NP1 |
| `np1_length` | int | Length of NP1 in bases |
| `np2` | str | Nucleotide sequence of the NP2 region (D-J; empty for VJ chains) |
| `np2_aa` | str | Amino acid translation of NP2 |
| `np2_length` | int | Length of NP2 in bases |

### Functionality & Isotype

| Field | Type | Description |
|-------|------|-------------|
| `productive` | bool | Whether the sequence is productive (in-frame, no stop codons) |
| `vj_in_frame` | bool | Whether V and J segments are in the same reading frame |
| `stop_codon` | bool | Whether the sequence contains a stop codon |
| `c_call` | str | Constant region allele (if enabled) |

### GenAIRR Ground-Truth Annotations

These fields are GenAIRR-specific and provide deep insight into the simulation process:

| Field | Type | Description |
|-------|------|-------------|
| `n_mutations` | int | Total number of somatic hypermutations applied |
| `mutation_rate` | float | Fraction of the sequence that was mutated |
| `n_pcr_errors` | int | Number of PCR-introduced base substitutions |
| `n_quality_errors` | int | Number of sequencing-quality base errors |
| `n_indels` | int | Number of sequencing-introduced insertions or deletions |
| `is_contaminant` | bool | Whether the sequence is a non-receptor contaminant |

## Truth vs Evidence-Based Calls

By default, `v_call`, `d_call`, and `j_call` reflect the "evidence-based" assignment—which allele the sequence *looks like* after mutations. In extreme SHM cases, a sequence might look more like a different allele than the one it was originally sampled from.

To export the original sampled alleles (the "absolute" truth), pass `expose_provenance=True` to `run_records()`:

```python
result = exp.run_records(n=100, expose_provenance=True)
rec = result[0]
rec["truth_v_call"]  # The allele originally sampled by the engine
rec["v_call"]        # The allele the final sequence matches most closely
```

## Coordinate System

All coordinates are **0-based**. End positions are **exclusive** (matching Python slice convention):

```python
rec = result[0]

# Extract the V segment from the sequence
v_segment = rec["sequence"][rec["v_sequence_start"]:rec["v_sequence_end"]]

# Extract the junction
junction = rec["sequence"][rec["junction_start"]:rec["junction_end"]]
```

When exporting to TSV/CSV using `result.to_csv(airr_strict=True)`, GenAIRR automatically converts these to **1-based inclusive** coordinates to comply with the AIRR Rearrangement specification.

## Next steps

- [Quick Start](/docs/getting-started/quick-start) — install and run your first simulation
- [Choosing a Config](/docs/getting-started/choosing-config) — pick the right species and chain type
- [Guides](/docs/guides/) — recipes for clonal families, antigen selection, and more
