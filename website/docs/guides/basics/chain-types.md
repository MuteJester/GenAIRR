---
title: Chain Types
sidebar_label: Chain Types
---

# Chain Types

GenAIRR supports a wide variety of B-cell receptor (BCR) and T-cell receptor (TCR) chains. The fundamental difference between these chains is their recombination architecture: **VDJ** (three segments) or **VJ** (two segments).

## VDJ Chains (Three Segments)

VDJ chains rearrange a **V** (Variable), **D** (Diversity), and **J** (Joining) segment. These segments are separated by two non-templated junctional regions (**NP1** and **NP2**).

```
5' ── [V Segment] ── NP1 ── [D Segment] ── NP2 ── [J Segment] ── 3'
```

| Chain | Receptor Type | Config Suffix |
|-------|---------------|---------------|
| **IGH** | BCR Heavy Chain | `_igh` |
| **TCRB** | TCR Beta | `_tcrb` |
| **TCRD** | TCR Delta | `_tcrd` |

### Example
```python
import GenAIRR as ga

# Simulate a human heavy chain
result = ga.Experiment.on("human_igh").recombine().run(n=1, seed=42)
rec = result[0]

print(rec["v_call"])  # e.g. "IGHV3-23*01"
print(rec["d_call"])  # e.g. "IGHD3-10*01"
print(rec["j_call"])  # e.g. "IGHJ4*02"
print(rec["np1"])     # Non-templated bases between V and D
print(rec["np2"])     # Non-templated bases between D and J
```

## VJ Chains (Two Segments)

VJ chains rearrange only a **V** and **J** segment, separated by a single non-templated region (**NP1**).

```
5' ── [V Segment] ── NP1 ── [J Segment] ── 3'
```

| Chain | Receptor Type | Config Suffix |
|-------|---------------|---------------|
| **IGK** | BCR Kappa Light Chain | `_igk` |
| **IGL** | BCR Lambda Light Chain | `_igl` |
| **TCRA** | TCR Alpha | `_tcra` |
| **TCRG** | TCR Gamma | `_tcrg` |

### Example
```python
# Simulate a human kappa light chain
result = ga.Experiment.on("human_igk").recombine().run(n=1, seed=42)
rec = result[0]

print(rec["v_call"])  # e.g. "IGKV3-20*01"
print(rec["d_call"])  # "" (Empty for VJ chains)
print(rec["j_call"])  # e.g. "IGKJ1*01"
print(rec["np1"])     # Non-templated bases between V and J
print(rec["np2"])     # "" (Empty for VJ chains)
```

## Output Differences Summary

When working with VJ chains, GenAIRR leaves D-related fields empty or `None` to maintain a consistent AIRR-compliant schema:

| Field | VDJ Chains | VJ Chains |
|-------|------------|-----------|
| `d_call` | Allele name | `""` (Empty string) |
| `d_sequence_start` | Integer | `None` or `0` |
| `np2` | Nucleotide sequence | `""` (Empty string) |
| `np2_length` | Integer | `0` |

## BCR vs TCR Simulation

While the underlying simulation mechanics (sampling, trimming, NP-insertion) are identical for both BCR and TCR chains, GenAIRR uses different **germline reference data** for each.

Additionally, GenAIRR includes biological guards to prevent unrealistic simulations. For example, calling `.mutate()` on a TCR chain will raise an error, as T-cell receptors do not undergo somatic hypermutation in the periphery.

```python
# This will raise a ValueError:
ga.Experiment.on("human_trb").recombine().mutate(count=10)
```

## Next steps

- [Choosing a Config](/docs/getting-started/choosing-config) — Full list of supported species and chains
- [Somatic Hypermutation](/docs/guides/options/shm) — Modeling B-cell maturation
- [Understanding Output](/docs/getting-started/interpreting-results) — Detailed field reference
