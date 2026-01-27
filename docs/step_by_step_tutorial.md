# Step-by-Step Tutorial: Building Your First Pipeline

This tutorial walks you through creating a GenAIRR simulation pipeline from scratch. Each section introduces one concept, explains why it matters, and shows the code.

**Prerequisites:** GenAIRR installed (`pip install GenAIRR`). See [Installation & Quick Start](getting_started.md) if you haven't set up yet.

---

## Step 1: Imports

```python
from GenAIRR import Pipeline, steps, HUMAN_IGH_OGRDB, S5F
```

| Import | What it is |
|--------|-----------|
| `Pipeline` | The pipeline runner — takes a config and a list of steps, executes them in order |
| `steps` | Namespace containing all built-in augmentation steps |
| `HUMAN_IGH_OGRDB` | Pre-built `DataConfig` for human heavy chain (OGRDB germline set) |
| `S5F` | Context-dependent somatic hypermutation model |

---

## Step 2: Choose a Mutation Model

A mutation model determines **how** nucleotide substitutions are introduced into the simulated sequence.

```python
mutation_model = S5F(min_mutation_rate=0.02, max_mutation_rate=0.08)
```

**What this does:** For each simulated sequence, a mutation rate is sampled uniformly from [0.02, 0.08]. The S5F model then applies mutations at that rate, using 5-mer context-dependent substitution probabilities derived from empirical data.

**Why S5F?** Real somatic hypermutation is not uniform — certain motifs (like WRC/GYW hotspots) mutate at much higher rates. S5F captures this pattern. Use `Uniform` if you want a simpler null model:

```python
from GenAIRR.mutation import Uniform
null_model = Uniform(min_mutation_rate=0.02, max_mutation_rate=0.08)
```

**Choosing mutation rate ranges:**

| Scenario | `min_mutation_rate` | `max_mutation_rate` |
|----------|--------------------|--------------------|
| Naive B cells | 0.001 | 0.01 |
| Memory B cells | 0.02 | 0.08 |
| Plasma cells | 0.05 | 0.25 |
| Mixed repertoire | 0.003 | 0.25 |

---

## Step 3: Create the Sequence Generation Step

The `SimulateSequence` step is always the first step in a pipeline. It performs V(D)J recombination and applies mutations.

```python
simulate_step = steps.SimulateSequence(
    mutation_model,
    productive=True
)
```

**Parameters:**

- `mutation_model` — the model created in Step 2
- `productive=True` — restricts output to sequences that are in-frame and lack stop codons. Set to `False` to include non-productive rearrangements (roughly 2/3 of all rearrangements are non-productive in biology).

!!! info "Productivity and performance"
    When `productive=True`, GenAIRR retries up to 25 times to find a valid rearrangement. At very high mutation rates (>0.2), this can slow generation. See [Productive Sequences](concepts/productivity.md) for details.

**What happens internally:**

1. Selects a random V, D, and J allele from the config (weighted by empirical usage)
2. Applies exonuclease trimming at segment junctions (from empirical distributions)
3. Adds N-nucleotides at junctions (P- and N-additions)
4. Concatenates segments into a complete sequence
5. Applies somatic hypermutation at the sampled rate
6. If `productive=True`, repeats until an in-frame sequence is generated

---

## Step 4: Add Position Correction Steps

After trimming, the exact boundaries between V, D, and J segments can become ambiguous — the same nucleotide sequence could correspond to multiple valid trim positions. These correction steps resolve that ambiguity to produce reliable ground-truth annotations.

```python
fix_v = steps.FixVPositionAfterTrimmingIndexAmbiguity()
fix_d = steps.FixDPositionAfterTrimmingIndexAmbiguity()
fix_j = steps.FixJPositionAfterTrimmingIndexAmbiguity()
```

**Why this matters:** If you use simulated data to benchmark an alignment tool, the ground-truth segment boundaries must be unambiguous. Without these corrections, multiple trim positions could explain the same observed sequence, leading to inconsistent ground truth.

!!! warning "Light chains"
    For kappa and lambda light chains, omit `FixDPositionAfterTrimmingIndexAmbiguity` and `CorrectForDTrims` — light chains have no D segment.

---

## Step 5: Add Biological Correction Steps

These steps apply further corrections based on biological constraints:

```python
correct_v = steps.CorrectForVEndCut()
correct_d = steps.CorrectForDTrims()
```

- **`CorrectForVEndCut()`** — adjusts the V-end position when the 3' end of the V segment was trimmed into the coding region
- **`CorrectForDTrims()`** — adjusts D segment boundaries after 5' and 3' trimming. Skip this for light chains.

---

## Step 6: Record the Mutation Rate

```python
distill = steps.DistillMutationRate()
```

`DistillMutationRate` computes and stores the final mutation rate in the `SimulationContainer`.

!!! danger "Step ordering"
    Place `DistillMutationRate` **before** any artifact steps (corruption, N-insertion, indels). Artifact steps modify the sequence in ways that are not biological mutations — measuring afterwards inflates the reported rate.

---

## Step 7: Add Sequencing Artifact Steps (Optional)

Real-world sequencing data contains imperfections. These steps simulate common NGS artifacts:

### 5' End Corruption

```python
corrupt = steps.CorruptSequenceBeginning(
    probability=0.7,
    event_weights=(0.4, 0.4, 0.2)
)
```

Simulates 5' end degradation — the leading portion of the sequence may be truncated, replaced with random nucleotides, or both. The `event_weights` tuple controls the probability of (add random bases, remove bases, add-after-remove).

### Read Length Enforcement

```python
enforce_length = steps.EnforceSequenceLength(max_length=576)
```

Truncates sequences exceeding `max_length` nucleotides, simulating fixed read-length limits of sequencing platforms.

### Ambiguous Base Insertion

```python
insert_ns = steps.InsertNs(n_ratio=0.02, probability=0.5)
```

With probability 0.5, replaces ~2% of bases with 'N' to simulate ambiguous base calls.

### Insertion/Deletion Errors

```python
indels = steps.InsertIndels(
    probability=0.5,
    max_indels=5,
    insertion_probability=0.5,
    deletion_probability=0.5
)
```

Introduces random insertions and deletions to simulate PCR/sequencing errors.

---

## Step 8: Assemble and Run the Pipeline

Combine all steps into a pipeline:

```python
pipeline = Pipeline(
    config=HUMAN_IGH_OGRDB,
    steps=[
        # Sequence generation
        steps.SimulateSequence(
            S5F(min_mutation_rate=0.02, max_mutation_rate=0.08),
            productive=True
        ),
        # Position corrections
        steps.FixVPositionAfterTrimmingIndexAmbiguity(),
        steps.FixDPositionAfterTrimmingIndexAmbiguity(),
        steps.FixJPositionAfterTrimmingIndexAmbiguity(),
        # Biological corrections
        steps.CorrectForVEndCut(),
        steps.CorrectForDTrims(),
        # Record mutation rate (before artifacts)
        steps.DistillMutationRate(),
        # Sequencing artifacts
        steps.CorruptSequenceBeginning(),
        steps.EnforceSequenceLength(),
        steps.InsertNs(),
        # Validation and errors
        steps.ShortDValidation(),
        steps.InsertIndels(),
    ]
)
```

Execute the pipeline:

```python
result = pipeline.execute()
data = result.get_dict()
```

### Inspecting the Output

The `SimulationContainer` holds all simulation metadata:

```python
# Sequence and alleles
print("Sequence:", data['sequence'][:60], "...")
print("V allele:", data['v_call'])
print("D allele:", data['d_call'])
print("J allele:", data['j_call'])

# Positions in the final sequence
print("V region:", data['v_sequence_start'], "-", data['v_sequence_end'])
print("D region:", data['d_sequence_start'], "-", data['d_sequence_end'])
print("J region:", data['j_sequence_start'], "-", data['j_sequence_end'])

# Mutation info
print("Mutation rate:", data['mutation_rate'])
print("Mutations:", data['mutations'])
print("Productive:", data['productive'])
```

---

## Step 9: Generate a Dataset

To produce a batch of sequences:

```python
import pandas as pd

sequences = []
for i in range(100):
    result = pipeline.execute()
    seq_dict = result.get_dict()
    seq_dict['id'] = f"seq_{i:04d}"
    sequences.append(seq_dict)

df = pd.DataFrame(sequences)
print(df[['id', 'v_call', 'd_call', 'j_call', 'mutation_rate']].head())
```

Export:

```python
# CSV
df.to_csv('simulated_sequences.csv', index=False)

# FASTA
with open('simulated_sequences.fasta', 'w') as f:
    for _, row in df.iterrows():
        f.write(f">{row['id']}\n{row['sequence']}\n")
```

---

## Variations

### Light Chain Pipeline

Light chains lack D segments. Omit D-related steps:

```python
from GenAIRR import HUMAN_IGK_OGRDB

kappa_pipeline = Pipeline(
    config=HUMAN_IGK_OGRDB,
    steps=[
        steps.SimulateSequence(
            S5F(min_mutation_rate=0.02, max_mutation_rate=0.08),
            productive=True
        ),
        steps.FixVPositionAfterTrimmingIndexAmbiguity(),
        steps.FixJPositionAfterTrimmingIndexAmbiguity(),
        steps.CorrectForVEndCut(),
        steps.DistillMutationRate(),
        steps.CorruptSequenceBeginning(),
        steps.EnforceSequenceLength(max_length=400),
        steps.InsertNs(),
        steps.InsertIndels(),
    ]
)
```

### Specific Allele Selection

Force particular V, D, J alleles:

```python
v_allele = HUMAN_IGH_OGRDB.v_alleles['IGHVF1-G1'][0]
d_allele = HUMAN_IGH_OGRDB.d_alleles['IGHD1-1'][0]
j_allele = HUMAN_IGH_OGRDB.j_alleles['IGHJ1'][0]

pipeline = Pipeline(
    config=HUMAN_IGH_OGRDB,
    steps=[
        steps.SimulateSequence(
            S5F(min_mutation_rate=0.02, max_mutation_rate=0.08),
            productive=True,
            specific_v=v_allele,
            specific_d=d_allele,
            specific_j=j_allele
        ),
        steps.FixVPositionAfterTrimmingIndexAmbiguity(),
        steps.FixDPositionAfterTrimmingIndexAmbiguity(),
        steps.FixJPositionAfterTrimmingIndexAmbiguity(),
        steps.CorrectForVEndCut(),
        steps.CorrectForDTrims(),
    ]
)
```

### Reproducible Output

```python
from GenAIRR import set_seed

set_seed(42)
result = pipeline.execute()
# Same seed → same sequence every time
```

---

## Troubleshooting

| Problem | Cause | Solution |
|---------|-------|----------|
| `TypeError` on step constructors | Using positional args | Use keyword arguments: `InsertNs(n_ratio=0.02, probability=0.5)` |
| All sequences identical | Seed set without reset | Call `reset_seed()` or use different seeds |
| `KeyError` on allele name | Wrong allele key format | Use OGRDB family keys (e.g., `IGHVF1-G1`), not IMGT names |
| Non-productive sequences | `productive=False` | Set `productive=True` in `SimulateSequence` |

---

## Next Steps

- **[How the Pipeline Works](genairr_flow.md)** — Understand the architecture in detail
- **[Parameter Reference](parameter_reference.md)** — Complete parameter documentation
- **[Best Practices](best_practices.md)** — Guidelines for realistic simulations
- **[Biological Context](biological_context.md)** — The immunobiology behind the simulation
