---
title: Biological Events
sidebar_label: Biological Events
---

# Biological Events

In a natural immune system, V(D)J recombination is a stochastic process driven by weighted gene usage frequencies and empirical trimming distributions. However, for benchmarking and controlled experiments, you often need to override these biological defaults.

GenAIRR provides two primary ways to control biological sampling: **Allele Locking** and **Allele Weighting**.

## Allele Locking with `.using()`

The `.using()` method allows you to restrict the simulation to specific V, D, or J alleles. This is essential for benchmarking alignment tools, as it allows you to generate sequences from a known biological origin.

```python
import GenAIRR as ga

# Lock to a specific V allele and a subset of J alleles
result = (
    ga.Experiment.on("human_igh")
    .using(
        v="IGHV3-23*01", 
        j=["IGHJ4*02", "IGHJ6*01"]
    )
    .recombine()
    .run(n=1000, seed=42)
)

# Every sequence in the result will use V3-23*01
assert all(rec["v_call"] == "IGHV3-23*01" for rec in result)
```

### Parameters

| Argument | Type | Description |
|----------|------|-------------|
| `v` | `str` or `list` | V allele name(s) to use. |
| `d` | `str` or `list` | D allele name(s) to use (VDJ chains only). |
| `j` | `str` or `list` | J allele name(s) to use. |

### How it Works
When you provide a single string, the engine is "locked" to that exact allele. If you provide a list, the engine will sample **uniformly** from that subset, ignoring the empirical frequencies defined in the reference configuration.

## Allele Weighting in `.recombine()`

If you want to bias the simulation toward certain alleles without completely excluding others, you can use the weighting parameters in the `.recombine()` method.

```python
# Significantly boost the frequency of one specific V allele
weights = {"IGHV3-23*01": 100.0}

result = (
    ga.Experiment.on("human_igh")
    .recombine(v_allele_weights=weights)
    .run(n=1000)
)
```

By default, every allele has a weight of `1.0`. By setting an allele's weight to `100.0`, you make it 100 times more likely to be sampled than any unlisted allele.

## Biological Constraints (Contracts)

The new engine introduces a concept called **Contracts**, which can enforce biological invariants like productivity during the recombination process. Instead of "filtering" after the fact, contracts prune the engine's random choices in real-time.

```python
from GenAIRR.contract import productive

# Use the productive() contract bundle to ensure in-frame rearrangements
result = (
    ga.Experiment.on("human_igh")
    .recombine()
    .run(n=1000, respect=productive())
)
```

For more details on ensuring realistic, functional sequences, see the [Productive Sequences](/docs/guides/options/productive) guide.

## Future Biological Passes

The new Rust-based engine architecture is modular. While the current version focuses on core V(D)J recombination and SHM, the **PassPlan** system is designed to support additional biological events such as:

*   **D-Gene Inversion:** Reversing the orientation of the D segment.
*   **Receptor Revision:** Secondary V-gene rearrangement.
*   **Class Switch Recombination (CSR):** Constant region switching (e.g., IgM to IgG).

These features are planned for future updates to the GenAIRR engine.

## Next steps

- [Productive Sequences](/docs/guides/options/productive) — Ensuring functional rearrangements
- [Somatic Hypermutation](/docs/guides/options/shm) — Modeling B-cell maturation
- [Metadata Accuracy](/docs/concepts/metadata-accuracy) — How GenAIRR tracks these events
