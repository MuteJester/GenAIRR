---
title: Productive Sequences
sidebar_label: Productive
---

# Productive Sequences

In immunology, a **productive** rearrangement is one that can produce a functional protein. This requires two conditions: the V and J gene segments must be joined in the same reading frame, and the resulting sequence must not contain any premature stop codons. Only about 1 in 3 random V(D)J rearrangements are productive — the rest are discarded by the cell.

GenAIRR evaluates productivity using a 5-rule validator that mirrors the biological checks a cell performs, and can optionally filter for productive sequences using a retry mechanism.

## What makes a sequence productive

The C engine runs five sequential checks after assembling each rearrangement:

| # | Rule | What it checks | Effect |
|---|------|---------------|--------|
| 1 | **Stop Codon** | Scans the entire sequence for in-frame stop codons (every 3rd base from position 0) | Fatal — stops validation immediately |
| 2 | **Frame Alignment** | Checks that `junction_start % 3 == 0`, `junction_end % 3 == 0`, and `junction_length % 3 == 0` | Sets `vj_in_frame = False` |
| 3 | **Junction Translatable** | Verifies the junction region can be translated (length divisible by 3, valid codons) | Fatal — stops validation immediately |
| 4 | **Conserved Cysteine** | First amino acid of the junction must be Cysteine (C) — the V region's 2nd conserved cysteine at CDR3 start | Non-fatal |
| 5 | **Conserved Anchor** | Last amino acid of the junction must be Tryptophan (W) or Phenylalanine (F) — the J region anchor residue | Non-fatal |

A sequence is `productive = True` only if **all five rules pass**. If any rule fails, the `note` field records the first failure reason.

### Why only ~20% of unfiltered sequences are productive

Random V(D)J recombination can join segments in any of three reading frames. Only one frame (1 in 3 chance) will produce in-frame V-J alignment. Even when in-frame, NP nucleotide additions can introduce stop codons, and the junction must preserve the conserved Cysteine and W/F anchor residues. Together, these constraints mean roughly 20% of unfiltered rearrangements pass all checks — consistent with the biological observation that most V(D)J recombination events produce non-functional receptors.

## Output fields

| Field | Type | Description |
|-------|------|-------------|
| `productive` | bool | `True` if all 5 rules pass |
| `vj_in_frame` | bool | `True` if V and J are in the same reading frame (rules 1+2) |
| `stop_codon` | bool | `True` if a stop codon was found anywhere in the sequence |
| `note` | str | First failure reason (e.g. `"VJ out of frame."`, `"Stop codon present."`) |

## Filtering for productive sequences

Pass `productive=True` to `.run()` to bias the simulation toward productive sequences:

```python
from GenAIRR import Experiment

result = Experiment.on("human_igh").run(n=1000, seed=42, productive=True)

productive_count = sum(1 for rec in result if rec["productive"])
print(f"{productive_count}/{len(result)} productive")  # ~998/1000
```

### How the retry loop works

When `productive=True`, the C engine uses a **retry loop** for the rearrangement phase. It repeats the allele sampling, trimming, and assembly steps up to 25 times until a productive rearrangement is found. Once a productive rearrangement is achieved, the remaining pipeline steps (mutation, corruption, etc.) run once on the successful rearrangement.

This means `productive=True` is **best-effort, not a hard guarantee**. If the engine exhausts all 25 attempts without finding a productive rearrangement, it proceeds with the last non-productive attempt. In practice, this is rare — you'll typically see 998-1000 out of 1000 sequences being productive.

The retry boundary is placed right after the `assess_functionality` step, so only the rearrangement portion is retried. Mutation and corruption steps run exactly once, regardless of how many rearrangement attempts were needed.

## Without filtering

By default, GenAIRR produces both productive and non-productive sequences, reflecting the natural distribution:

```python
result = Experiment.on("human_igh").run(n=1000, seed=42)

productive = sum(1 for rec in result if rec["productive"])
print(f"{productive}/{len(result)} productive")  # ~201/1000
```

The ~20% productive rate is biologically realistic. The two main reasons for non-productivity are:
- **VJ out of frame** (~43% of non-productive) — V and J joined in the wrong reading frame
- **Stop codon present** (~57% of non-productive) — an in-frame stop codon appears in the sequence

## Checking individual records

```python
rec = result[0]

if rec["productive"]:
    print("Productive:", rec["junction_aa"])
else:
    print("Non-productive:", rec["note"])
    print("  In frame:", rec["vj_in_frame"])
    print("  Stop codon:", rec["stop_codon"])
```

## When to use productive filtering

- **Training ML models for annotation**: use `productive=True` — most real datasets are pre-filtered for productive sequences
- **Benchmarking annotation tools**: use both — test your tool's ability to correctly flag non-productive sequences
- **Studying repertoire diversity**: use default (no filter) — the natural productive/non-productive ratio is part of the biology
