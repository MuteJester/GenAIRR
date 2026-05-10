---
title: Choosing a Config
sidebar_label: Choosing a Config
---

# Choosing a Config

Every GenAIRR simulation starts with a **config** — a pre-built germline reference that defines which V, D, and J alleles are available, their usage frequencies, trimming distributions, and NP-region models for a specific species and chain type.

## Listing available configs

```python
import GenAIRR as ga

# Returns a list of all 106 built-in configuration names
print(ga.list_configs())
# ['ALPACA_IGH_IMGT', 'CAT_IGK_IMGT', ..., 'ZEBRAFISH_TCRD_IMGT']
```

GenAIRR ships with **106 built-in configs** covering **23 species and strains**, sourced from IMGT and OGRDB.

## Config naming and aliases

Configs follow the pattern `SPECIES_CHAIN_SOURCE` (e.g., `HUMAN_IGH_OGRDB`). When passed to `Experiment.on()`, you can use a lowercase shorthand. 

For many species, the shorthand points to the `IMGT` version by default. However, for **Human BCR** (IGH, IGK, IGL), the shorthand points to the **OGRDB** versions, which are generally more accurate for population-level simulations.

```python
import GenAIRR as ga

# These use the OGRDB version for Human IGH:
ga.Experiment.on("human_igh")        # Preferred alias
ga.Experiment.on("HUMAN_IGH_OGRDB")  # Explicit name

# This uses the IMGT version:
ga.Experiment.on("HUMAN_IGH_IMGT")
```

For the Mouse C57BL/6J strain:

```python
ga.Experiment.on("mouse_c57bl6j_igh")
```

## Chain types

| Chain | Description | Has D gene? |
|-------|-------------|:-----------:|
| `IGH` | B-cell heavy chain (immunoglobulin) | Yes |
| `IGK` | B-cell kappa light chain | No |
| `IGL` | B-cell lambda light chain | No |
| `TCRA` | T-cell receptor alpha | No |
| `TCRB` | T-cell receptor beta | Yes |
| `TCRD` | T-cell receptor delta | Yes |
| `TCRG` | T-cell receptor gamma | No |

Chains without a D gene produce two-segment (VJ) rearrangements. Chains with a D gene produce three-segment (VDJ) rearrangements with two NP regions (NP1 and NP2).

## Full species coverage

| Species | BCR | TCR |
|---------|-----|-----|
| **Human** | IGH, IGK, IGL | TCRA, TCRB, TCRD, TCRG |
| **Mouse** | IGH, IGK, IGL | TCRA, TCRB, TCRD, TCRG |
| **Mouse (C57BL/6J)** | IGH, IGK, IGL | TCRA, TCRB, TCRD, TCRG |
| **Rabbit** | IGH, IGK, IGL | TCRA, TCRB, TCRD, TCRG |
| **Dog** | IGH, IGK, IGL | TCRA, TCRB, TCRD, TCRG |
| **Gorilla** | IGH, IGK, IGL | TCRA, TCRB, TCRD, TCRG |
| **Rhesus** | IGH, IGK, IGL | TCRA, TCRB, TCRD, TCRG |
| **Ferret** | IGH, IGK, IGL | TCRA, TCRB, TCRD, TCRG |
| **Cow** | IGH, IGK, IGL | TCRA, TCRB, TCRD, TCRG |
| **Cat** | IGK, IGL | TCRA, TCRB, TCRD, TCRG |
| **Sheep** | IGH, IGK, IGL | TCRA, TCRB, TCRD |
| **Pig** | IGH, IGK, IGL | TCRB, TCRG |
| **Rat** | IGH, IGK, IGL | — |
| **Chicken** | IGH, IGL | — |
| **Horse** | IGH, IGK | — |
| **Goat** | IGK, IGL | — |
| **Cynomolgus** | IGH | TCRB |
| **Dromedary** | IGK | TCRB, TCRG |
| **Zebrafish** | IGH | TCRA, TCRD |
| **Trout** | IGH | TCRB |
| **Alpaca** | IGH | — |
| **Platypus** | IGH | — |
| **Salmon** | IGH | — |

:::note
Available chains per species depend on what germline data is published in IMGT and OGRDB. If a species/chain combination isn't listed, GenAIRR doesn't yet have germline reference data for it.
:::

## Quick examples

```python
import GenAIRR as ga

# Human BCR heavy chain (using OGRDB default)
ga.Experiment.on("human_igh").recombine().run(n=1000)

# Mouse kappa light chain
ga.Experiment.on("mouse_igk").recombine().run(n=1000)

# Rabbit TCR beta
ga.Experiment.on("rabbit_tcrb").recombine().run(n=1000)

# Mouse C57BL/6J strain heavy chain
ga.Experiment.on("mouse_c57bl6j_igh").recombine().run(n=1000)
```

## Next steps

- [Quick Start](/docs/getting-started/quick-start) — install and run your first simulation
- [Understanding Output](/docs/getting-started/interpreting-results) — what each of the ~70 fields means
