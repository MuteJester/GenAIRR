# DataConfig In Depth

The `DataConfig` is GenAIRR's central data container. It encapsulates everything needed to simulate sequences for a specific species and chain type: germline allele references, trimming distributions, junction nucleotide models, and ambiguity correction structures.

This page explains what a DataConfig contains, how each part works, and why it is structured the way it is.

---

## What Is a DataConfig?

A DataConfig is a **portable, self-contained reference container**. When you pass a DataConfig to a Pipeline, every step in that pipeline can access alleles, distributions, and correction maps without needing external files or global state.

```python
from GenAIRR import Pipeline, steps, HUMAN_IGH_OGRDB, S5F

pipeline = Pipeline(
    config=HUMAN_IGH_OGRDB,  # ← DataConfig passed here
    steps=[
        steps.SimulateSequence(S5F(min_mutation_rate=0.01, max_mutation_rate=0.05), productive=True),
        # ... steps automatically receive the config
    ]
)
```

This design means:

- **No global state** — different pipelines can use different configs simultaneously
- **Reproducibility** — the config fully specifies the simulation parameters
- **Portability** — a pickled DataConfig can be shared between researchers, ensuring identical simulation conditions

---

## ConfigInfo: Metadata

Every DataConfig carries a `ConfigInfo` metadata object:

| Field | Type | Description |
|-------|------|-------------|
| `species` | `Species` enum | The organism (Human, Mouse, Rat, Rabbit, etc. — 30+ species supported) |
| `chain_type` | `ChainType` enum | The receptor chain (BCR_HEAVY, BCR_LIGHT_KAPPA, BCR_LIGHT_LAMBDA, TCR_ALPHA, TCR_BETA, etc.) |
| `reference_set` | `str` | Source of the germline data (e.g., "OGRDB", "IMGT") |
| `last_updated` | `date` | When the config was last modified |
| `has_d` | `bool` | Whether this chain type includes D segments |

The `has_d` flag is particularly important — it controls whether the simulation generates NP2 junctions, D segment trimming, and D-related correction steps.

!!! tip "Why metadata matters"
    When publishing results generated with GenAIRR, including the DataConfig's metadata (species, chain type, reference set, date) makes your simulation fully reproducible. Another researcher with the same DataConfig will produce statistically identical output given the same seed.

---

## Allele Dictionaries

The core of any DataConfig is its allele collections:

```python
config.v_alleles  # Dict[str, List[VAllele]]
config.d_alleles  # Dict[str, List[DAllele]]
config.j_alleles  # Dict[str, List[JAllele]]
config.c_alleles  # Dict[str, List[CAllele]]
```

Each dictionary maps **gene names** to lists of **allele objects**:

```python
config.v_alleles['IGHVF1-G1']
# → [VAllele(name='IGHVF1-G1*01', ...), VAllele(name='IGHVF1-G1*02', ...), ...]
```

### Allele Object Attributes

Every allele object (VAllele, DAllele, JAllele, CAllele) carries:

| Attribute | Type | Description |
|-----------|------|-------------|
| `name` | `str` | Full allele identifier (e.g., `IGHVF1-G1*01`) |
| `gapped_seq` | `str` | IMGT-gapped reference sequence |
| `ungapped_seq` | `str` | Sequence with gaps removed (used for simulation) |
| `length` | `int` | Length of gapped sequence |
| `ungapped_len` | `int` | Length of ungapped sequence |
| `family` | `str` | Gene family (e.g., `IGHVF1`) |
| `gene` | `str` | Gene name (e.g., `IGHVF1-G1`) |
| `anchor` | `int` | Structural anchor position (see below) |

### Anchor Positions

Anchors are biologically critical positions that define structural landmarks:

- **V alleles**: The conserved **Cysteine** at IMGT position 104 (codon TGT or TGC at positions 309–311 in the gapped sequence). This marks the start of CDR3.
- **J alleles**: The conserved **Trp/Phe** motif in framework region 4. The anchor stores both the position and reading frame (0, 1, or 2).
- **D alleles**: No anchor (D segments have no defined structural landmark).
- **C alleles**: No anchor.

!!! note "OGRDB naming convention"
    Built-in configs use OGRDB-style allele names (`IGHVF1-G1*01`), not IMGT-style (`IGHV1-2*01`). The OGRDB convention groups alleles by sequence similarity using Allele Sequence Clustering (see ASC Tables below).

---

## Gene Usage Probabilities

The `gene_use_dict` attribute stores the probability of selecting each gene during simulation:

```python
config.gene_use_dict  # Dict[str, float]
# e.g., {'IGHVF1-G1': 0.023, 'IGHVF1-G2': 0.018, ...}
```

When `SimulateSequence` selects an allele, it first chooses a gene weighted by these probabilities, then randomly selects one of the alleles within that gene. This captures the biological reality that some genes are used far more frequently than others.

---

## Trimming Distributions

The `trim_dicts` attribute contains probability distributions for how many nucleotides are trimmed from each segment end during V(D)J recombination:

```python
config.trim_dicts  # Dict[str, Dict[str, Dict[str, Dict[int, float]]]]
```

### Structure

The dictionary is hierarchical with four levels:

```
trim_dicts
├── "V_3"              ← Segment + end (V trimmed at 3')
│   ├── "IGHVF1"       ← Allele family
│   │   ├── "IGHVF1-G1"  ← Gene within family
│   │   │   ├── 0: 0.47   ← Trim 0 bases: 47% probability
│   │   │   ├── 1: 0.26   ← Trim 1 base: 26% probability
│   │   │   ├── 2: 0.095
│   │   │   └── ...
│   │   └── "IGHVF1-G2"
│   │       └── ...
│   └── "IGHVF2"
│       └── ...
├── "D_5"              ← D segment trimmed at 5'
├── "D_3"              ← D segment trimmed at 3'
├── "J_5"              ← J segment trimmed at 5'
└── "C_3"              ← C segment trimmed at 3'
```

### Trim Keys

| Key | Meaning |
|-----|---------|
| `V_3` | 3' end of V (V is never trimmed at 5') |
| `D_5` | 5' end of D |
| `D_3` | 3' end of D |
| `J_5` | 5' end of J (J is never trimmed at 3') |
| `C_3` | 3' end of C (C is never trimmed at 5') |

### How Trimming Works

During simulation, GenAIRR:

1. Looks up the selected allele's family and gene
2. Retrieves the trim distribution for that specific gene
3. Samples a trim amount from the distribution
4. Validates that the trim doesn't destroy the allele's anchor or exceed its length
5. Returns the trimmed sequence

The per-gene granularity captures the biological reality that different genes have different trimming profiles — some are trimmed more aggressively than others.

---

## NP Region Model

GenAIRR simulates non-templated (NP) nucleotide additions using a **first-order Markov chain**. Three attributes control this:

### NP_lengths

```python
config.NP_lengths  # Dict[str, Dict[int, float]]
# e.g., {'NP1': {0: 0.052, 1: 0.040, 2: 0.056, ...}, 'NP2': {...}}
```

Probability distribution over possible NP region lengths. Length 0 means no nucleotides are added (segments join directly).

### NP_first_bases

```python
config.NP_first_bases  # Dict[str, Dict[str, float]]
# e.g., {'NP1': {'A': 0.11, 'C': 0.25, 'G': 0.28, 'T': 0.36}, 'NP2': {...}}
```

The probability of each nucleotide being the first base in the NP region.

### NP_transitions

```python
config.NP_transitions  # Dict[str, Dict[int, Dict[str, Dict[str, float]]]]
```

Position-specific transition matrices. For each position in the NP region, given the current nucleotide, what is the probability of each next nucleotide?

```
NP_transitions
├── "NP1"
│   ├── 0                    ← Position in NP region
│   │   ├── "A"              ← Current nucleotide
│   │   │   ├── "A": 0.15   ← P(next=A | current=A, position=0)
│   │   │   ├── "C": 0.35
│   │   │   ├── "G": 0.28
│   │   │   └── "T": 0.22
│   │   ├── "C"
│   │   │   └── ...
│   │   └── ...
│   ├── 1
│   │   └── ...
│   └── ...
└── "NP2"
    └── ...
```

### How NP Generation Works

1. Sample a length from `NP_lengths`
2. Sample the first nucleotide from `NP_first_bases`
3. For each subsequent position:
    - Look up the transition probabilities for the current position and current nucleotide
    - Sample the next nucleotide
4. Return the generated sequence

This Markov model captures the biological observation that NP nucleotide composition is not uniform — TdT (terminal deoxynucleotidyl transferase) has nucleotide preferences that vary by position.

---

## Correction Maps

Correction maps are precomputed data structures that enable GenAIRR's ambiguity resolution framework. They answer the question: "Given modifications to the sequence, which alleles can no longer be distinguished?"

### Trim Similarity Maps

Six maps handle trimming-induced ambiguity:

| Map | Key structure | Purpose |
|-----|---------------|---------|
| `V_3_TRIM_SIMILARITY_MAP` | `allele → trim_amount → [equivalent_alleles]` | V alleles indistinguishable after 3' trimming |
| `V_5_TRIM_SIMILARITY_MAP` | `allele → trim_amount → [equivalent_alleles]` | V alleles indistinguishable after 5' corruption |
| `J_5_TRIM_SIMILARITY_MAP` | `allele → trim_amount → [equivalent_alleles]` | J alleles indistinguishable after 5' trimming |
| `J_3_TRIM_SIMILARITY_MAP` | `allele → trim_amount → [equivalent_alleles]` | J alleles indistinguishable after 3' corruption |
| `C_3_TRIM_SIMILARITY_MAP` | `allele → trim_amount → [equivalent_alleles]` | C alleles indistinguishable after 3' trimming |
| `D_5_3_TRIM_SIMILARITY_MAP` | `allele → (trim_5, trim_3) → [equivalent_alleles]` | D alleles indistinguishable after bilateral trimming |

**How they're built:** For every allele and every possible trim amount, the builder trims the allele's sequence and checks if the trimmed result is a substring of any other allele's sequence. All matches are stored.

**How they're used:** Steps like `CorrectForVEndCut` and `CorrectForDTrims` look up equivalent alleles using these maps, then add them to the container's call lists.

### V_N_AMBIGUITY_CORRECTION_GRAPH

This structure handles a different type of ambiguity: when 'N' bases are inserted at positions that distinguish between V alleles, those alleles become indistinguishable.

**How it works:** The graph stores, for each pair of V alleles, the set of positions where they differ. When Ns are inserted at some of these discriminating positions, the graph identifies which alleles can no longer be told apart. The `InsertNs` step uses this graph to extend `v_call` with N-ambiguous alleles.

---

## ASC Tables

**Allele Sequence Clustering (ASC)** is the method used to organize alleles in OGRDB-based configs. Instead of IMGT's gene-level grouping, ASC uses hierarchical clustering based on sequence similarity:

1. **Compute pairwise distances** between all alleles (Hamming distance, ignoring gaps and Ns)
2. **Cluster at two thresholds**:
    - **Family level** (75% similarity): Groups alleles into broad families (e.g., IGHVF1, IGHVF2)
    - **Gene level** (95% similarity): Groups alleles within families into genes (e.g., IGHVF1-G1, IGHVF1-G2)
3. **Assign names**: `{segment}F{family}-G{gene}*{allele_index}`

The ASC table (`config.asc_tables["V"]`) is a pandas DataFrame recording the clustering results, including which alleles are duplicates and what SNPs distinguish near-identical alleles.

!!! info "Why ASC instead of IMGT naming?"
    ASC groups alleles by actual sequence similarity rather than historical naming conventions. Two alleles with the same IMGT gene name can have quite different sequences, while alleles with different names may be nearly identical. ASC resolves this by letting the sequences speak for themselves.

---

## Built-in Configs

GenAIRR ships with five pre-built DataConfig instances:

| Config | Chain | Species | Source | Alleles (approx.) |
|--------|-------|---------|--------|-------------------|
| `HUMAN_IGH_OGRDB` | Heavy (IGH) | Human | OGRDB | ~350 V, ~40 D, ~12 J |
| `HUMAN_IGH_EXTENDED` | Heavy (IGH) | Human | OGRDB | ~700 V (extended set) |
| `HUMAN_IGK_OGRDB` | Kappa (IGK) | Human | OGRDB | ~150 V, ~10 J |
| `HUMAN_IGL_OGRDB` | Lambda (IGL) | Human | OGRDB | ~180 V, ~14 J |
| `HUMAN_TCRB_IMGT` | TCR Beta (TRB) | Human | IMGT | ~140 V, ~3 D, ~14 J |

These are stored as serialized (pickled) files and **lazily loaded** — the data is only deserialized from disk when you first access the config. Subsequent accesses use the cached object.

```python
from GenAIRR import HUMAN_IGH_OGRDB

# First access: loads from disk (~6 MB pickle)
config = HUMAN_IGH_OGRDB

# Subsequent accesses: returns cached object
config2 = HUMAN_IGH_OGRDB  # Same object, no disk I/O
```

---

## Creating Your Own DataConfig

If you work with a species or allele set not covered by the built-in configs, you can create a custom DataConfig using `CustomDataConfigGenerator`. This requires:

1. **FASTA reference files** for V, D, J (and optionally C) segments
2. **An annotated dataset** (CSV) with columns for sequence, gene calls, trim amounts, and corruption events

The generator derives all distributions (trimming, NP, gene usage) from the provided data.

For full details, see:

- **[Custom DataConfig guide](../custom_data_config.md)** — Attribute-by-attribute documentation
- **[Creating Custom DataConfig tutorial](../tutorials/Creating Custom DataConfig from FASTA Files.ipynb)** — Interactive notebook walkthrough

---

## Inspecting a DataConfig

You can explore any DataConfig interactively:

```python
from GenAIRR import HUMAN_IGH_OGRDB

config = HUMAN_IGH_OGRDB

# Metadata
print(config.metadata)
# ConfigInfo(species=HUMAN, chain_type=BCR_HEAVY, reference_set='OGRDB', ...)

# Allele counts
print(f"V genes: {len(config.v_alleles)}")
print(f"D genes: {len(config.d_alleles)}")
print(f"J genes: {len(config.j_alleles)}")

# List gene families
print(list(config.v_alleles.keys())[:5])

# Inspect an allele
allele = config.v_alleles['IGHVF1-G1'][0]
print(allele.name, allele.ungapped_len, allele.family)

# Check trimming distribution
print(config.trim_dicts['V_3']['IGHVF1']['IGHVF1-G1'])

# Check NP length distribution
print(config.NP_lengths['NP1'])
```
