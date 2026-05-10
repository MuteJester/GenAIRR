---
title: Sequencing Artifacts
sidebar_label: Artifacts
---

# Sequencing Artifacts

Real-world AIRR-seq data is rarely "perfect." From library preparation to the sequencing instrument itself, technical noise is introduced at every step. GenAIRR provides a comprehensive suite of corruption methods to model these effects, ensuring your synthetic data is as realistic as possible.

## Corruption Methods

Unlike biological events, corruption methods represent technical noise. They are prefixed with `.corrupt_*()` and can be chained in any order, though placing them after biological steps (like `.recombine()` and `.mutate()`) is recommended to mirror a real experimental workflow.

```python
import GenAIRR as ga

result = (
    ga.Experiment.on("human_igh")
    .recombine()
    .mutate(count=15)

    # 1. PCR Errors
    .corrupt_pcr(count=(0, 5))

    # 2. Sequencing Quality Errors
    .corrupt_quality(count=(5, 10))

    # 3. End Loss (Primer trimming / read-end decay)
    .corrupt_5prime_loss(length=(5, 25))
    .corrupt_3prime_loss(length=(0, 15))

    # 4. Post-processing noise
    .corrupt_indels(count=2, insertion_prob=0.5)
    .corrupt_ns(count=3)

    .run(n=1000, seed=42)
)
```

### Key Artifacts

| Method | What it models | Effect on IR |
|-------|---------------|--------------|
| `.corrupt_pcr()` | Polymerase substitutions during amplification | Uniform base substitutions |
| `.corrupt_quality()` | Instrument-level base call errors | Position-dependent substitutions |
| `.corrupt_5prime_loss()` | Primer trimming or 5' signal loss | Permanent removal of 5' bases |
| `.corrupt_3prime_loss()` | Read-end signal decay | Permanent removal of 3' bases |
| `.corrupt_indels()` | Sequencing-induced insertions/deletions | Shifts the reading frame |
| `.corrupt_ns()` | Ambiguous base calls (base-caller failure) | Replaces bases with `N` |
| `.corrupt_rev_comp()` | Antisense read capture (~50% of reads) | Full reverse-complementation |
| `.corrupt_contaminants()` | Non-receptor DNA spike-ins | Overwrites sequence with random bases |

## The `count` and `length` Parameters

Just like `.mutate()`, all corruption methods accept a flexible distribution for the number of events:

*   **Fixed:** `count=5` (Exactly 5 events)
*   **Uniform Range:** `count=(0, 10)` (Random integer between 0 and 10)
*   **Empirical:** `count=[(0, 0.8), (5, 0.2)]` (80% chance of no events, 20% chance of 5 events)

## PCR vs. Sequencing Errors

GenAIRR distinguishes between different types of substitutions to provide more granular ground truth:

*   **PCR Errors (`n_pcr_errors`):** Models early-stage amplification bias. These are treated as standard base substitutions and are counted separately in the output.
*   **Quality Errors (`n_quality_errors`):** Models instrument-level noise. In the output `sequence`, these are often represented as **lowercase** bases (if the original was uppercase) to distinguish them from biological mutations.

## End Loss (Trimming)

5' and 3' loss methods permanently remove bases from the simulation's nucleotide pool. This is a critical artifact to model because it hides the beginning of the V gene or the end of the J gene, making alignment and allele calling significantly more challenging for analysis tools.

Because GenAIRR uses a **Persistent IR**, any bases removed by a `.corrupt_*_loss()` step are gone for all subsequent passes in the pipeline.

## Indels and Frame Shifts

The `.corrupt_indels()` method introduces insertions and deletions. Unlike point mutations, indels shift the reading frame, which can introduce premature stop codons or break the conserved anchors (Cys/Trp) required for junction identification.

GenAIRR automatically recomputes the **Codon Rail** and **Junction Boundaries** after every indel event, ensuring that fields like `productive` and `junction_aa` are always accurate relative to the final, corrupted sequence.

## Output Annotations

All artifacts are tracked in the output record:

```python
rec = result[0]
print(rec["n_pcr_errors"])      # Count of PCR-induced substitutions
print(rec["n_quality_errors"])  # Count of sequencing quality errors
print(rec["n_indels"])          # Count of insertions/deletions
print(rec["is_contaminant"])    # True if the sequence was replaced
```

## Next steps

- [Somatic Hypermutation](/docs/guides/options/shm) — Modeling biological variation
- [Understanding Output](/docs/getting-started/interpreting-results) — Detailed field reference
- [Clonal Structure](/docs/guides/basics/experiment-dsl#clonal-structure) — Sharing artifacts within a family
