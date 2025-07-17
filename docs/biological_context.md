# Biological Context Guide

Understanding the biological processes simulated by GenAIRR.

## Immunoglobulin Structure Basics

### Heavy Chain Structure
```
5' V-region — D-region — J-region — C-region 3'
   └─────── Variable Domain ───────┘
```

- **V (Variable)**: Highly diverse, determines antigen specificity
- **D (Diversity)**: Only in heavy chains, adds additional diversity
- **J (Joining)**: Connects variable to constant region
- **C (Constant)**: Determines antibody class (IgG, IgA, etc.)

### Light Chain Structure (Kappa/Lambda)
```
5' V-region — J-region — C-region 3'
   └─ Variable Domain ─┘
```

Light chains lack D segments, making them simpler than heavy chains.

## Biological Processes Simulated

### 1. V(D)J Recombination
**What happens in biology**: During B cell development, one V, D (heavy only), and J gene segment randomly combine to form the initial antibody.

**How GenAIRR simulates it**:
```python
# Random selection of segments
HeavyChainSequence.create_random(config)

# Or controlled selection
SimulateSequence(
    S5F(), True,
    specific_v=v_allele,
    specific_d=d_allele, 
    specific_j=j_allele
)
```

### 2. Somatic Hypermutation (SHM)
**What happens in biology**: After antigen exposure, B cells in germinal centers undergo rapid mutation in their variable regions to improve antigen binding.

**Key characteristics**:
- Occurs primarily in germinal centers
- Targets CDR regions more than framework regions
- Context-dependent (certain sequence motifs mutate more)
- Rate: ~10⁻³ to 10⁻⁴ per base pair per cell division

**How GenAIRR simulates it**:
```python
# S5F model captures context-dependent mutation patterns
S5F(min_mutation_rate=0.02, max_mutation_rate=0.08)  # Memory B cell range

# Uniform model for simpler, random mutations
Uniform(min_mutation_rate=0.01, max_mutation_rate=0.05)
```

### 3. Junctional Diversity
**What happens in biology**: During V(D)J recombination, imprecise joining creates additional diversity through:
- Exonuclease trimming of gene segments
- Random nucleotide addition by TdT enzyme
- P-nucleotide addition from hairpin resolution

**How GenAIRR simulates it**:
```python
# These steps model the biological trimming and joining process
CorrectForVEndCut(),      # V segment 3' trimming
CorrectForDTrims(),       # D segment 5' and 3' trimming
```

### 4. Sequencing Artifacts
**What happens in reality**: NGS introduces various errors and biases:
- 5' end degradation
- Random base calls ('N's)
- PCR/sequencing errors
- Insertion/deletion errors

**How GenAIRR simulates it**:
```python
CorruptSequenceBeginning(0.7, [0.4, 0.4, 0.2], 576, 210, 310, 50),  # 5' degradation
InsertNs(0.02, 0.5),                                                  # Ambiguous bases
InsertIndels(0.5, 5, 0.5, 0.5)                                       # Sequencing errors
```

## Mutation Patterns by Cell Type

### Naive B Cells
- **Location**: Bone marrow, peripheral blood
- **Mutations**: Very few (0.1-1%)
- **Function**: Recently formed, not antigen-experienced

```python
SimulateSequence(S5F(0.001, 0.01), productive=True)
```

### Memory B Cells  
- **Location**: Secondary lymphoid organs
- **Mutations**: Moderate (2-8%)
- **Function**: Long-lived, antigen-experienced

```python
SimulateSequence(S5F(0.02, 0.08), productive=True)
```

### Plasma Cells
- **Location**: Bone marrow, inflamed tissues
- **Mutations**: High (5-25%)
- **Function**: Antibody-secreting effector cells

```python
SimulateSequence(S5F(0.05, 0.25), productive=True)
```

### Germinal Center B Cells
- **Location**: Germinal centers of lymphoid organs
- **Mutations**: Variable, rapidly increasing
- **Function**: Undergoing selection and affinity maturation

```python
SimulateSequence(S5F(0.01, 0.15), productive=True)
```

## Chain Pairing

### BCR Structure
Each B cell receptor consists of:
- 2 identical heavy chains
- 2 identical light chains (either kappa OR lambda)
- Ratio: ~60% kappa, 40% lambda in humans

### TCR Structure
Each T cell receptor consists of:
- 1 alpha chain (similar to light chain)
- 1 beta chain (similar to heavy chain, has D segment)

## Understanding Output Fields

### Sequence Positions
- **v_sequence_start/end**: Where V segment appears in final sequence
- **d_sequence_start/end**: Where D segment appears (heavy/TCR-β only)
- **j_sequence_start/end**: Where J segment appears
- **junction_start/end**: CDR3 region boundaries

### Germline Positions
- **v_germline_start/end**: Position within original V gene
- **d_germline_start/end**: Position within original D gene
- **j_germline_start/end**: Position within original J gene

### Quality Metrics
- **productive**: In-frame, no stop codons
- **vj_in_frame**: V and J segments maintain reading frame
- **stop_codon**: Presence of premature stop codons
- **mutation_rate**: Frequency of mutations in sequence

## Clinical Relevance

### Immune Repertoire Diversity
- Healthy individuals: >10¹¹ unique BCR sequences
- Repertoire changes with age, disease, vaccination
- Clonal expansion indicates immune response

### Disease Applications
- **Cancer**: Tumor-infiltrating lymphocytes
- **Autoimmune**: Self-reactive antibodies
- **Immunodeficiency**: Reduced diversity
- **Vaccines**: Tracking immune response

## Research Applications

### 1. Vaccine Studies
Track how B cell repertoires change post-vaccination:
```python
# Pre-vaccination (naive-like)
pre_vax = SimulateSequence(S5F(0.001, 0.01), True)

# Post-vaccination (activated)
post_vax = SimulateSequence(S5F(0.02, 0.08), True)
```

### 2. Aging Studies
Model how repertoires change with age:
```python
# Young repertoire (high diversity)
young = SimulateSequence(S5F(0.005, 0.02), True)

# Aged repertoire (more mutations, less diversity)
aged = SimulateSequence(S5F(0.02, 0.1), True)
```

### 3. Disease Modeling
Compare healthy vs. disease states:
```python
# Healthy memory response
healthy = SimulateSequence(S5F(0.02, 0.06), True)

# Autoimmune (potentially higher mutation)
autoimmune = SimulateSequence(S5F(0.05, 0.15), True)
```

This biological context helps explain why GenAIRR's simulation steps exist and how to choose appropriate parameters for your research questions.
