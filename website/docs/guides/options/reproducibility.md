---
title: Reproducibility
sidebar_label: Reproducibility
---

# Reproducibility

Reproducible simulation is essential for scientific research — you must be able to generate the exact same dataset when rerunning an experiment, sharing results with collaborators, or debugging an analysis pipeline.

GenAIRR provides **bit-for-bit identical results** across all supported platforms (**Linux**, **macOS**, and **Windows**) through its deterministic Rust simulation kernel.

## How Seeds Work

GenAIRR uses a high-performance, seeded pseudo-random number generator (PRNG) based on the **Xoshiro256++** algorithm. When you pass a `seed` to `.run()`, `.run_records()`, or `.stream_records()`, the engine ensures that every random decision—from allele sampling to mutation site selection—is fully deterministic.

```python
import GenAIRR as ga

# Run two identical experiments with the same seed
exp = ga.Experiment.on("human_igh").recombine().mutate(count=10)

r1 = exp.run(n=1000, seed=42)
r2 = exp.run(n=1000, seed=42)

# Every sequence and its metadata will be identical
assert r1[0].final_simulation().bases() == r2[0].final_simulation().bases()
```

### Deterministic Batches
When you call `.run(n=1000, seed=42)`, the engine uses `seed + i` for each of the 1,000 iterations. This ensures that:
1.  Individual sequences within a batch are independent.
2.  The batch as a whole is reproducible.
3.  Consecutive batches can be "stitched" together by offsetting the seed (e.g., `seed=42` followed by `seed=1042`).

## The Addressed Trace

A major feature of the new engine is the **Addressed Trace**. Every time the simulator makes a random choice, it records that choice in a hierarchical log keyed by a unique "address."

```python
# Inspect the trace of a simulated outcome
outcome = result[0]
trace = outcome.trace

# Look up specific decisions
print(trace.find("sample_allele.v").value)    # The index of the sampled V allele
print(trace.find("np.np1.length").value)      # The length of the NP1 region
```

This trace makes GenAIRR simulations **fully auditable**. You don't just get a sequence; you get a complete record of *why* that sequence was generated, which is invaluable for debugging and verifying biological models.

## Cross-Platform Consistency

Because the simulation logic and PRNG are implemented in Rust and compiled as a native extension, GenAIRR is immune to the platform-specific differences in floating-point math or system-level `rand()` implementations that affect many other simulators.

Whether you run your simulation on a Linux server or a Windows workstation, `seed=42` will produce the exact same sequence of nucleotides and metadata.

## What Affects Reproducibility?

A simulation is reproducible only if the **Seed**, **Reference Data**, and **Pipeline** are identical. Changing any of the following will result in different data:

*   **Reference Configuration:** Using a different species or a newer version of a DataConfig.
*   **Pipeline Passes:** Adding, removing, or reordering steps (e.g., adding a `.corrupt_ns()` step).
*   **Pass Parameters:** Changing a mutation count or an allele weight.
*   **GenAIRR Version:** Internal algorithm changes between major versions may alter the sequence of PRNG calls.

## Next steps

- [Productive Sequences](/docs/guides/options/productive) — Deterministic filtering
- [Simulation Pipeline](/docs/concepts/simulation-pipeline) — How passes use the PRNG
- [Metadata Accuracy](/docs/concepts/metadata-accuracy) — How the trace powers ground truth
