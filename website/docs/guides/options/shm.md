---
title: Somatic Hypermutation
sidebar_label: SHM
---

# Somatic Hypermutation

Somatic hypermutation (SHM) is the process by which B cells introduce point mutations into the variable regions of their immunoglobulin genes during an immune response. This is the primary mechanism behind antibody affinity maturation.

GenAIRR models SHM through the `.mutate()` method on the `Experiment` object. It provides both context-dependent (S5F) and position-independent (Uniform) mutation models.

## Applying SHM

Use `.mutate()` after your `.recombine()` step to introduce mutations:

```python
import GenAIRR as ga

# Generate human heavy-chain sequences with SHM
result = (
    ga.Experiment.on("human_igh")
    .recombine()
    .mutate(model="s5f", count=(5, 25))
    .run(n=1000, seed=42)
)

rec = result[0]
print(rec["n_mutations"])    # e.g. 14
print(rec["mutation_rate"])  # e.g. 0.038
```

:::note
Somatic hypermutation is a B-cell phenomenon. Calling `.mutate()` on a TCR-configured experiment (like `mouse_tcrb`) will raise a `ValueError` to prevent biological misuse.
:::

## Mutation Models

GenAIRR supports two primary mutation models:

### 1. The S5F Model (`model="s5f"`)
This is the default model. It implements the **Somatic 5-mer Frequency (S5F)** model, where a position's probability of being mutated depends on its 5-nucleotide context (2 bases upstream and 2 bases downstream).

*   **Hotspots:** Captures the biological reality that AID (activation-induced cytidine deaminase) targets specific motifs like **WRC**.
*   **Iterative Application:** The new Rust engine recomputes the context weights after *every* single mutation. This is more accurate than older engines that compute weights once, as a mutation at one position immediately changes the 5-mer context of its four neighbors.

### 2. The Uniform Model (`model="model='uniform'"`)
A position-independent model where every base in the V(D)J segments has an equal probability of being targeted. Each selected position is replaced with a uniformly drawn A/C/G/T base (excluding the original base).

## Controlling the Mutation Count

The `count` argument determines how many mutations are applied per sequence:

*   **Fixed Count:** `count=10` — every sequence receives exactly 10 mutations.
*   **Uniform Range:** `count=(5, 15)` — a random integer in `[5, 15]` is drawn for each sequence.
*   **Empirical Distribution:** `count=[(5, 1.0), (10, 2.0)]` — uses a weighted distribution (here, 10 mutations are twice as likely as 5).

## The Mutation History

Because the new engine uses a **Persistent IR**, every mutation is recorded in the simulation's history. While the final AIRR record reports the aggregate `n_mutations` and `mutation_rate`, the underlying **Addressed Trace** contains the exact position and base for every event:

*   `mutate.s5f.site[0]`: The pool index of the first mutation.
*   `mutate.s5f.base[0]`: The new nucleotide base.

## Germline Alignment

The `germline_alignment` field in the output allows you to visualize the mutations by comparing it to the `sequence` field.

<div className="seq-vis">
  <div><span className="sv-label">sequence</span><span className="sv-v">gaggtgcag</span><span className="sv-mut">T</span><span className="sv-v">tggtggagt</span></div>
  <div><span className="sv-label">germline_alignment</span><span className="sv-v">gaggtgcag</span><span className="sv-v">c</span><span className="sv-v">tggtggagt</span></div>
  <div><span className="sv-label"></span><span style={{opacity:0.5}}>         ↑ position 9: germline 'c' mutated to 'T'</span></div>
</div>

In the `sequence` field, mutated bases are shown in **uppercase** (if they were originally uppercase in the germline) or **lowercase** (if they are technical artifacts). The `germline_alignment` always shows the original germline base.

## Next steps

- [Sequencing Artifacts](/docs/guides/options/artifacts) — Modeling realistic sequencing noise
- [Understanding Output](/docs/getting-started/interpreting-results) — Detailed reference for mutation fields
- [Clonal Structure](/docs/guides/basics/experiment-dsl#clonal-structure) — Simulating families with shared mutations
