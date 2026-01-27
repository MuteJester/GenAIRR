# Parameter Reference

Detailed explanation of every configurable parameter in GenAIRR, organized by component.

---

## Mutation Model Parameters

Both `S5F` and `Uniform` accept the same base parameters:

| Parameter | Type | Default | Range | Description |
|-----------|------|---------|-------|-------------|
| `min_mutation_rate` | `float` | `0.0` | 0.0 – 1.0 | Lower bound of the mutation rate range. Each simulated sequence samples a rate uniformly from [`min`, `max`]. |
| `max_mutation_rate` | `float` | `0.0` | 0.0 – 1.0 | Upper bound of the mutation rate range. Setting both to 0 produces unmutated sequences. |
| `productive` | `bool` | `False` | — | When `True`, the mutation process avoids introducing stop codons and preserves CDR3 anchor residues (Cys/Trp). |

**S5F only:**

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `custom_model` | `str` | `None` | Path to a custom 5-mer substitution model file. If `None`, uses the built-in S5F targeting model. |

### Recommended Mutation Rate Ranges

| Biological scenario | `min_mutation_rate` | `max_mutation_rate` |
|---------------------|--------------------|--------------------|
| Unmutated (germline) | 0.0 | 0.0 |
| Naive B cells | 0.001 | 0.01 |
| Memory B cells | 0.02 | 0.08 |
| Plasma cells | 0.05 | 0.25 |
| Germinal center (active) | 0.01 | 0.15 |
| Mixed/broad repertoire | 0.003 | 0.25 |
| TCR (minimal SHM) | 0.0 | 0.01 |

---

## SimulateSequence Parameters

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `mutation_model` | `MutationModel` | *required* | An `S5F` or `Uniform` instance that controls how mutations are applied. |
| `productive` | `bool` | `False` | If `True`, the step retries until it generates an in-frame sequence without stop codons. Roughly 1/3 of random rearrangements are productive. |
| `specific_v` | `VAllele` | `None` | Force a specific V allele instead of random weighted selection. Obtain from `config.v_alleles['GENE_NAME'][index]`. |
| `specific_d` | `DAllele` | `None` | Force a specific D allele. Only applicable for chains with D segments (heavy, TRB). |
| `specific_j` | `JAllele` | `None` | Force a specific J allele. |

**Notes:**
- When `specific_v/d/j` are `None`, alleles are selected randomly using the empirical usage weights in the DataConfig.
- You can specify any combination — e.g., fix V and J but randomize D.

---

## CorruptSequenceBeginning Parameters

All parameters are **keyword-only** (must use `param=value` syntax).

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `probability` | `float` | `0.7` | Probability that corruption is applied to a given sequence. Set to 0.0 to disable. |
| `event_weights` | `tuple[float, float, float]` | `(0.4, 0.4, 0.2)` | Relative weights for three corruption event types: **(1)** add random nucleotides to the 5' end, **(2)** remove nucleotides from the 5' end, **(3)** remove then add. Weights are normalized internally; they do not need to sum to 1. |
| `nucleotide_add_coefficient` | `int` | `210` | Shape parameter controlling the distribution of how many nucleotides are added. Higher values produce longer additions. |
| `nucleotide_remove_coefficient` | `int` | `310` | Shape parameter controlling how many nucleotides are removed. Higher values produce more aggressive trimming. |
| `nucleotide_add_after_remove_coefficient` | `int` | `50` | Shape parameter for the add-after-remove event. Typically smaller than the add coefficient. |
| `random_sequence_add_probability` | `float` | `1.0` | Probability that added nucleotides are random (vs. matching the germline). 1.0 = always random; 0.0 = always germline-matching. |

---

## EnforceSequenceLength Parameters

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `max_length` | `int` | `576` | Maximum allowed sequence length in nucleotides. Sequences longer than this are truncated from the 3' end. Choose based on your sequencing platform's read length. |

**Common values:**

| Platform | Typical max_length |
|----------|--------------------|
| Illumina MiSeq (2x300 paired) | 576 |
| Illumina MiSeq (2x250 paired) | 480 |
| Illumina NextSeq (2x150) | 280 |
| PacBio / Nanopore | 1000+ (or omit this step) |

---

## InsertNs Parameters

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `n_ratio` | `float` | `0.02` | Fraction of sequence positions to replace with 'N'. A value of 0.02 means ~2% of bases become ambiguous. |
| `probability` | `float` | `0.5` | Probability of applying any N-insertion to a given sequence. When the event fires, `n_ratio` fraction of positions are replaced. |

---

## InsertIndels Parameters

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `probability` | `float` | `0.5` | Probability of applying indels to a given sequence. |
| `max_indels` | `int` | `5` | Maximum number of indel events per sequence. The actual count is sampled uniformly from [1, max_indels]. |
| `insertion_probability` | `float` | `0.5` | Relative weight for insertion events. |
| `deletion_probability` | `float` | `0.5` | Relative weight for deletion events. The ratio of `insertion_probability` to `deletion_probability` determines the insertion/deletion balance. |

---

## ShortDValidation Parameters

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `short_d_length` | `int` | `5` | D segments shorter than this threshold (in nucleotides) are flagged as "short" and their annotations may be adjusted. |

---

## Position Correction Steps

The following steps take **no parameters**:

| Step | Applicable chains |
|------|-------------------|
| `FixVPositionAfterTrimmingIndexAmbiguity()` | All |
| `FixDPositionAfterTrimmingIndexAmbiguity()` | Heavy, TRB |
| `FixJPositionAfterTrimmingIndexAmbiguity()` | All |
| `CorrectForVEndCut()` | All |
| `CorrectForDTrims()` | Heavy, TRB |
| `DistillMutationRate()` | All |
| `FilterTCRDJAmbiguities()` | TRB only |

---

## Seed Management

| Function | Parameter | Description |
|----------|-----------|-------------|
| `set_seed(seed)` | `seed: int` | Sets the global random seed across all RNGs (Python `random`, NumPy). Enables deterministic output. |
| `get_seed()` | — | Returns the current seed value, or `None` if not set. |
| `reset_seed()` | — | Clears the seed, returning to non-deterministic behavior. |

---

## Pipeline Constructor

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `steps` | `List[AugmentationStep]` | *required* | Ordered list of steps to execute. |
| `config` | `DataConfig` | `None` | The germline data configuration. Passed to each step automatically. Required for all pipelines that include `SimulateSequence`. |

---

## simulate() Convenience Function

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `config` | `DataConfig` | *required* | Germline data configuration. |
| `mutation_model` | `MutationModel` | *required* | S5F or Uniform instance. |
| `productive` | `bool` | `True` | Only generate productive sequences. |
| `n` | `int` | `1` | Number of sequences. Returns a list when `n > 1`. |
