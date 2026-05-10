---
title: Clonal Families
sidebar_label: Clonal Families
---

# Clonal Families

In a real immune repertoire, multiple B-cells often belong to the same **clonal family** — they share the same initial V(D)J recombination event but have diverged through somatic hypermutation.

GenAIRR provides a first-class way to simulate these hierarchical structures using the `.with_clonal_structure()` method.

## The Clonal Pipeline

To simulate clonal families, you split your simulation pipeline into two distinct phases:

1.  **Parent Phase (Pre-Fork):** Steps that run once per family (e.g., `.recombine()`).
2.  **Descendant Phase (Post-Fork):** Steps that run once per read within the family (e.g., `.mutate()`, `.corrupt_*()`).

```python
import GenAIRR as ga

# Configure a clonal experiment
exp = (
    ga.Experiment.on("human_igh")
    
    # 1. Parent Phase: Shared recombination
    .recombine()
    
    # 2. Fork into 10 families, each with 20 descendants
    .with_clonal_structure(n_clones=10, size=20)
    
    # 3. Descendant Phase: Individual divergence
    .mutate(count=(5, 15))
    .corrupt_pcr(count=2)
)

# Run the experiment
# n is optional here; if provided, it must equal n_clones * size (200)
result = exp.run_records(seed=42)

print(len(result))  # 200
```

## How It Works: The IR Fork

The engine uses the **Persistent IR** to make clonal simulation incredibly efficient.

*   When the pipeline reaches the `.with_clonal_structure()` step, it captures the current state of the sequence (the "parent" IR).
*   It then "forks" this IR into the specified number of descendants.
*   Each descendant starts from the exact same biological sequence (same V/D/J alleles, same trim lengths, same NP regions) and then proceeds through the post-fork steps independently with its own random seed.

## Output Annotations

Every record produced by a clonal experiment includes a `clone_id` field (0-indexed) that identifies its family.

```python
for rec in result:
    print(f"Sequence {rec['sequence_id']} belongs to Clone {rec['clone_id']}")
```

## Advanced Usage: Shared Mutations

If you want descendants to share some mutations (mimicking a lineage tree), you can add a `.mutate()` step *before* the fork:

```python
exp = (
    ga.Experiment.on("human_igh")
    .recombine()
    
    # These mutations are shared by the entire family
    .mutate(count=5)
    
    .with_clonal_structure(n_clones=10, size=20)
    
    # These mutations introduce divergence within the family
    .mutate(count=(5, 15))
)
```

In this example, every sequence in a clone will have the same 5 "ancestral" mutations plus its own unique set of descendant mutations.

## Why Use Clonal Simulation?

*   **Clonal Grouping:** Generate datasets to benchmark clonal grouping and lineage reconstruction algorithms (like Change-O).
*   **Realistic Divergence:** Model how sequencing artifacts and SHM affect sequences that are biologically related.
*   **Performance:** Simulating clonal families is faster than simulating 200 independent sequences because the recombination phase is only executed 10 times instead of 200.

## Next steps

- [Somatic Hypermutation](/docs/guides/options/shm) — Diving deeper into mutation models
- [Experiment DSL](/docs/guides/basics/experiment-dsl) — Full reference for pipeline building
- [Simulation Pipeline](/docs/concepts/simulation-pipeline) — How the engine handles passes and forks
