---
title: Productive Sequences
sidebar_label: Productive
---

# Productive Sequences

In adaptive immunity, a **productive** rearrangement is one that can produce a functional protein. This requires the V and J gene segments to be joined in the same reading frame and the resulting sequence to be free of premature stop codons.

GenAIRR uses an efficient **Contract-Aware Sampling** system to ensure sequences meet these biological requirements.

## What Makes a Sequence Productive?

GenAIRR evaluates productivity based on five fundamental biological rules:

1.  **In-Frame Junction:** The length of the junction (V-anchor to J-anchor) must be a multiple of 3.
2.  **No Stop Codons:** The sequence must not contain any in-frame stop codons (`TAA`, `TAG`, `TGA`).
3.  **VJ Frame Alignment:** The V and J segments must be in the same reading frame.
4.  **Conserved Cysteine:** The junction must start with the conserved 2nd Cysteine of the V region.
5.  **Conserved Anchor:** The junction must end with the conserved Tryptophan (W) or Phenylalanine (F) of the J region.

A sequence is marked as `productive: True` only if all five rules are satisfied.

## Contract-Aware Sampling

Instead of generating a sequence and "checking" if it's productive after the fact (which is inefficient), GenAIRR now uses **Contracts**. When a contract is active, the engine prunes its random choices in real-time.

For example, when sampling the length of an NP region, the engine will *only* consider lengths that will result in an in-frame junction.

### Using the Productive Contract

To ensure your simulation produces productive sequences, pass the `productive()` contract bundle to the `.run()` or `.run_records()` method using the `respect` parameter:

```python
import GenAIRR as ga

# Generate 1,000 guaranteed productive sequences
result = (
    ga.Experiment.on("human_igh")
    .recombine()
    .run(n=1000, seed=42, respect=ga.productive())
)

# Nearly all (if not all) sequences will be productive
productive_count = sum(1 for rec in result if rec["productive"])
print(f"{productive_count}/1000 productive")
```

### Why Contracts are Superior to Retries

1.  **Efficiency:** Retrying a failed rearrangement can be slow, especially for complex pipelines. Contracts ensure the sequence is "correct by construction," eliminating wasted computation.
2.  **Statistical Integrity:** Contracts allow the engine to maintain the correct biological distributions (like trimming lengths) while staying within the allowed search space.
3.  **Deterministic:** Because choices are pruned before they are made, the process remains fully deterministic and recorded in the trace.

## Strict vs. Permissive Sampling

When using contracts, you can control how the engine behaves if a specific step has *no* valid choices (e.g., if a V-allele is so short that no trimming/NP combination can make it in-frame).

*   **Permissive (`strict=False`, default):** If no admissible candidate exists for a specific choice, the engine will fall back to unconstrained sampling for that one draw and continue. This ensures the simulation finishes but may result in a non-productive sequence.
*   **Strict (`strict=True`):** If no admissible candidate exists, the engine raises a `StrictSamplingError`. This is useful for debugging custom reference data or ensuring 100% compliance.

```python
# Force strict compliance
try:
    result = exp.run(n=1000, respect=ga.productive(), strict=True)
except ga.StrictSamplingError as e:
    print(f"Sampling failed at {e.address}: {e}")
```

## Output fields

| Field | Type | Description |
|-------|------|-------------|
| `productive` | bool | `True` if all productivity rules pass. |
| `vj_in_frame` | bool | `True` if V and J are joined in the same reading frame. |
| `stop_codon` | bool | `True` if an in-frame stop codon was detected. |

## When to Use Productive Filtering

*   **Model Training:** Use `respect=ga.productive()` when training ML models on "clean" data, as most experimental repertoires are pre-filtered for productivity.
*   **Tool Benchmarking:** Run simulations *without* the contract to test an annotation tool's ability to correctly identify and flag non-productive sequences.
*   **Repertoire Research:** Use the default (no contract) to study the natural ratio of productive vs. non-productive rearrangements in a species.

## Next steps

- [Biological Events](/docs/guides/options/biology) — Allele locking and weighting
- [Somatic Hypermutation](/docs/guides/options/shm) — Maturation effects
- [Understanding Output](/docs/getting-started/interpreting-results) — Field reference
