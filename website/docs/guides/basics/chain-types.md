---
title: Chain Types
sidebar_label: Chain Types
---

# Chain Types

GenAIRR supports both B-cell receptor (BCR) and T-cell receptor (TCR) chains. The key difference is whether the chain uses VDJ or VJ recombination.

## VDJ chains (three gene segments)

These chains rearrange a V, D, and J gene segment with two N/P junctional regions:

```
5' ── V segment ── NP1 ── D segment ── NP2 ── J segment ── 3'
```

| Chain | Receptor | Config suffix |
|-------|----------|---------------|
| IGH | BCR heavy chain | `_igh` |
| TCRB | TCR beta | `_tcrb` |
| TCRD | TCR delta | `_tcrd` |

```python
from GenAIRR import Experiment

rec = Experiment.on("human_igh").run(n=1, seed=42)[0]
print(rec["d_call"])      # e.g. "IGHD2-21*02"
print(rec["np1_region"])  # e.g. "ATACGTACGC"
print(rec["np2_region"])  # e.g. "GATCATC"
```

## VJ chains (two gene segments)

These chains rearrange only V and J with a single N/P region:

```
5' ── V segment ── NP1 ── J segment ── 3'
```

| Chain | Receptor | Config suffix |
|-------|----------|---------------|
| IGK | BCR kappa light chain | `_igk` |
| IGL | BCR lambda light chain | `_igl` |
| TCRA | TCR alpha | `_tcra` |
| TCRG | TCR gamma | `_tcrg` |

```python
rec = Experiment.on("human_igk").run(n=1, seed=42)[0]
print(rec["d_call"])      # "" (empty — no D gene)
print(rec["np1_region"])  # e.g. "C"
print(rec["np2_region"])  # "" (empty — no NP2)
print(rec["np2_length"])  # 0
```

## Differences in output

For VJ chains, all D-gene fields are empty or zero:

| Field | VDJ value | VJ value |
|-------|-----------|----------|
| `d_call` | allele name | `""` |
| `d_sequence_start` | position | `0` |
| `d_sequence_end` | position | `0` |
| `d_trim_5` / `d_trim_3` | trim amounts | `0` |
| `np2_region` | nucleotides | `""` |
| `np2_length` | length | `0` |

## BCR vs TCR

From GenAIRR's perspective, the simulation mechanics are identical for BCR and TCR chains — the only difference is the germline reference data. However, some ops are specific to certain chain types:

- `with_isotype_rates()` — only meaningful for IGH (class switch recombination)
- `with_d_inversion()` — only meaningful for VDJ chains (IGH, TCRB, TCRD)

```python
# BCR heavy chain with CSR
from GenAIRR.ops import rate, with_isotype_rates

result = (
    Experiment.on("human_igh")
    .mutate(
        rate(0.02, 0.05),
        with_isotype_rates(),
    )
    .run(n=100, seed=42)
)
print(result[0]["c_call"])  # e.g. "IGHG3*17"
```
