---
sidebar_position: 2
title: Simulation Introspection
---

# Simulation Introspection

In the new Rust engine, GenAIRR provides a powerful way to understand the "story" of how a sequence was built through the **Addressed Trace**. Every random decision made during simulation is recorded in a hierarchical log, making the process fully auditable.

## The Addressed Trace

When you run an experiment, GenAIRR returns an `Outcome` object (accessible via `SimulationResult.outcomes`). Each outcome contains a `trace` that you can query using "addresses."

```python
import GenAIRR as ga

# Run a simulation
exp = ga.Experiment.on("human_igh").recombine().mutate(count=15)
result = exp.run(n=1, seed=42)

# Access the trace of the first sequence
outcome = result.outcomes[0]
trace = outcome.trace

# Find specific decisions
v_choice = trace.find("sample_allele.v")
print(f"Sampled V Allele Index: {v_choice.value}")

np1_len = trace.find("np.np1.length")
print(f"NP1 Length: {np1_len.value}")
```

### Common Trace Addresses

| Address | Description |
|---------|-------------|
| `sample_allele.v/d/j` | The index of the sampled V, D, or J allele. |
| `trim.v_3` | Number of bases trimmed from the 3' end of the V segment. |
| `np.np1.length` | The length of the first NP region. |
| `mutate.s5f.count` | The total number of SHM mutations applied. |
| `mutate.s5f.site[i]` | The pool position of the i-th mutation. |
| `corrupt.pcr.count` | The number of PCR substitutions. |

## Auditing and Debugging

The trace is invaluable for verifying biological models or debugging complex simulation pipelines. For example, you can verify exactly why a sequence became non-productive by checking the trim lengths and NP additions recorded in the trace.

```python
# Check all recorded choices
for choice in trace.choices():
    print(f"{choice.address}: {choice.value}")
```

## Replaying Simulations

Because the trace contains every random choice, you can theoretically "replay" a simulation by providing the same trace to the engine. This ensures bit-for-bit identity and allows you to inspect the IR (Intermediate Representation) at any point in the simulation's history.

## Human-Readable Narration

While the raw trace is machine-readable, a formatted, human-readable "narrator" tool is currently available via the **GenAIRR MCP Server**. This tool translates the trace into a biological story:

> "Selected IGHV3-23*01, trimmed 3bp from the 3' end, inserted 6bp of NP1 (ATGCGT), and applied 12 S5F mutations..."

A native Python `narrate()` function for the new engine is planned for a future update.

## Next steps

- [Sequence Visualization](/docs/showcase/visualize-sequence) — Visually dissect simulated sequences
- [Metadata Accuracy](/docs/concepts/metadata-accuracy) — How the trace powers ground truth
- [Simulation Pipeline](/docs/concepts/simulation-pipeline) — How passes record their decisions
