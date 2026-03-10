---
title: The Experiment DSL
sidebar_label: Experiment DSL
---

# The Experiment DSL

GenAIRR uses a fluent builder called `Experiment` that models the stages of a real sequencing experiment. You configure what you want, then call `.run()` to generate sequences.

## Minimal example

```python
from GenAIRR import Experiment

result = Experiment.on("human_igh").run(n=1000, seed=42)
```

This produces 1,000 unmutated, uncorrupted human heavy-chain rearrangements. Every phase is optional — you only add what you need.

## The five phases

The Experiment DSL has five biological phases. Each corresponds to a stage of a real sequencing experiment:

```python
from GenAIRR import Experiment
from GenAIRR.ops import (
    with_d_inversion, with_receptor_revision, using,
    rate, model, with_isotype_rates, with_antigen_selection,
    with_primer_mask, with_umi, with_pcr,
    with_5prime_loss, with_3prime_loss, with_quality_profile,
    with_reverse_complement,
    with_indels, with_ns, with_contaminants,
)

result = (
    Experiment.on("human_igh")

    # 1. V(D)J recombination
    .recombine(
        with_d_inversion(0.15),
        with_receptor_revision(0.05),
    )

    # 2. Somatic hypermutation
    .mutate(
        model("s5f"),
        rate(0.01, 0.05),
        with_isotype_rates(),
        with_antigen_selection(0.5),
    )

    # 3. Library preparation
    .prepare(
        with_primer_mask(),
        with_umi(12),
        with_pcr(error_rate=1e-4, cycles=30),
    )

    # 4. Sequencing
    .sequence(
        with_5prime_loss(min_remove=5, max_remove=30),
        with_3prime_loss(min_remove=5, max_remove=20),
        with_quality_profile(base=0.001, peak=0.02),
    )

    # 5. Post-sequencing observation
    .observe(
        with_indels(prob=0.005),
        with_ns(prob=0.005),
        with_contaminants(rate=0.01),
    )

    .run(n=1000, seed=42)
)
```

### Phase summary

| Phase | Method | What it models | Available ops |
|-------|--------|----------------|---------------|
| Recombination | `.recombine()` | Biological events during V(D)J joining | `with_d_inversion`, `with_receptor_revision`, `using` |
| Mutation | `.mutate()` | Somatic hypermutation and maturation | `model`, `rate`, `with_isotype_rates`, `with_antigen_selection` |
| Preparation | `.prepare()` | Wet-lab sample processing | `with_primer_mask`, `with_umi`, `with_pcr` |
| Sequencing | `.sequence()` | Instrument-level effects | `with_5prime_loss`, `with_3prime_loss`, `with_quality_profile`, `with_reverse_complement`, `long_read`, `paired_end` |
| Observation | `.observe()` | Post-sequencing noise | `with_indels`, `with_ns`, `with_contaminants` |

## Phase order matters

Phases must be called in biological order: `.recombine()` before `.mutate()` before `.prepare()` before `.sequence()` before `.observe()`. You can skip any phase, but you cannot reorder them:

```python
# Valid — skip recombine and prepare
Experiment.on("human_igh").mutate(rate(0.05, 0.1)).sequence(with_5prime_loss()).run(n=100)

# Valid — just mutation
Experiment.on("human_igh").mutate(rate(0.02, 0.08)).run(n=100)

# Valid — just artifacts, no mutation
Experiment.on("human_igh").observe(with_indels(prob=0.01)).run(n=100)
```

## Ops are validated at construction

Each op function validates its parameters when called, not at run time. This means errors are caught early:

```python
from GenAIRR.ops import rate

rate(0.02, 0.08)   # valid
rate(-0.1, 0.5)    # raises ValueError immediately
rate(0.5, 0.02)    # raises ValueError (min > max)
```

## `.run()` vs `.compile()`

`.run(n=..., seed=...)` is the simple path — it compiles and simulates in one call, returning a `SimulationResult`:

```python
result = Experiment.on("human_igh").mutate(rate(0.02, 0.08)).run(n=1000, seed=42)
print(len(result))  # 1000
```

`.compile(seed=...)` returns a `CompiledSimulator` for streaming or repeated use:

```python
sim = Experiment.on("human_igh").mutate(rate(0.02, 0.08)).compile(seed=42)

# Stream one record at a time (infinite iterator)
for record in sim.stream():
    process(record)
    if done:
        break
```

## Productive filtering

Pass `productive=True` to `.run()` to bias the simulation toward productive sequences (in-frame, no stop codons in the junction):

```python
result = Experiment.on("human_igh").run(n=1000, seed=42, productive=True)

productive_count = sum(1 for rec in result if rec["productive"])
print(f"{productive_count}/{len(result)} productive")  # nearly all
```

:::note
The productive filter uses a retry loop with a maximum attempt limit. The vast majority of sequences will be productive, but a small number may slip through if the retry budget is exhausted. Always check the `productive` field if you need a strict guarantee.
:::

## Next steps

- [Somatic Hypermutation](/docs/guides/options/shm) — mutation models, rates, and selection
- [Sequencing Artifacts](/docs/guides/options/artifacts) — 5'/3' loss, PCR, quality profiles
- [Biological Events](/docs/guides/options/biology) — D-inversion, receptor revision, allele locking
