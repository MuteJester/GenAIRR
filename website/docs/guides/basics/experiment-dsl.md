---
title: The Experiment DSL
sidebar_label: Experiment DSL
---

# The Experiment DSL

GenAIRR uses a fluent builder called `Experiment` that models the stages of a real adaptive immune receptor experiment. You configure the biological and technical stages of your simulation by chaining methods, then call `.run()` to generate your data.

## Minimal example

```python
import GenAIRR as ga

# Generate 1,000 human heavy-chain rearrangements
result = ga.Experiment.on("human_igh").recombine().run(n=1000, seed=42)
```

This produces 1,000 unmutated, uncorrupted human heavy-chain rearrangements. The `.recombine()` step is mandatory if you want V(D)J rearrangement; without it, you'd be simulating an empty sequence.

## Building a Pipeline

The `Experiment` object allows you to build a sophisticated simulation pipeline by chaining methods in their natural biological and technical order.

```python
import GenAIRR as ga

result = (
    ga.Experiment.on("human_igh")
    
    # 1. Biological V(D)J recombination
    # You can restrict alleles via .using()
    .using(v="IGHV3-23*01", j=["IGHJ4*02", "IGHJ6*01"])
    .recombine()

    # 2. Somatic hypermutation (B-cells)
    .mutate(model="s5f", count=(5, 25))

    # 3. Sequencing artifacts (5'/3' loss)
    .corrupt_5prime_loss(length=(5, 30))
    .corrupt_3prime_loss(length=(5, 20))

    # 4. Post-sequencing noise (Indels and N-bases)
    .corrupt_indels(count=(0, 2), insertion_prob=0.5)
    .corrupt_ns(count=5)

    # 5. Add AIRR-compliant metadata to every record
    .with_metadata(sample_id="P1", donor="D001")

    .run(n=1000, seed=42)
)
```

### Key Methods

| Method | Purpose | Key Arguments |
|-------|---------|---------------|
| `.recombine()` | Standard V(D)J joining | `trim` (bool), `v_allele_weights` (dict) |
| `.mutate()` | Somatic hypermutation | `model` ("s5f" or "uniform"), `count` |
| `.corrupt_5prime_loss()` | Primer/signal loss at 5' end | `length` |
| `.corrupt_3prime_loss()` | Read-end degradation | `length` |
| `.corrupt_indels()` | Sequencing insertions/deletions | `count`, `insertion_prob` |
| `.corrupt_ns()` | Ambiguous base calls | `count` |
| `.corrupt_pcr()` | PCR amplification errors | `count` |
| `.corrupt_quality()` | Sequencing quality errors | `count` |
| `.with_metadata()` | Inject AIRR sample columns | `sample_id`, `donor`, etc. |
| `.using()` | Lock specific alleles | `v`, `d`, `j` (names or lists) |

## Customizing Recombination

The `.recombine()` method can be biased by providing specific allele weights. This is useful for simulating specific repertoire biases without completely locking out other alleles.

```python
# Bias toward a specific V gene family
v_weights = {"IGHV3-23*01": 50.0, "IGHV3-7*01": 10.0}

exp = ga.Experiment.on("human_igh").recombine(v_allele_weights=v_weights)
```

## Adding Custom Metadata

You can attach sample-level metadata (like donor ID or sample name) that will be included in every generated AIRR record. This is essential for integrating synthetic data into multi-sample analysis pipelines.

```python
exp = (
    ga.Experiment.on("human_igh")
    .with_metadata(donor="D123", sample_id="S001", custom_tag="Synthetic_v1")
    .recombine()
)
```

## The `count` and `length` Arguments

Most mutation and corruption methods take a `count` or `length` argument that defines how many events to apply per sequence. This argument is highly flexible:

*   **Fixed Integer:** `count=10` — every sequence gets exactly 10 events.
*   **Uniform Range:** `count=(5, 15)` — every sequence gets a uniform random integer between 5 and 15.
*   **Empirical Distribution:** `count=[(5, 1.0), (10, 2.0)]` — a list of `(value, weight)` tuples. In this example, 10 events are twice as likely as 5.

## Clonal Structure

To simulate clonal families (multiple descendants from the same parent rearrangement), use `.with_clonal_structure()`:

```python
exp = (
    ga.Experiment.on("human_igh")
    .recombine() # Runs once per parent
    .with_clonal_structure(n_clones=10, size=20)
    .mutate(count=(5, 15)) # Runs once per descendant
)

# 10 clones * 20 descendants = 200 records
result = exp.run()
```

Steps before `.with_clonal_structure()` are shared by all members of the clone; steps after it introduce divergence.

## Running the Experiment

### `.run(n=..., seed=...)`
The simplest way to get results. Returns a `SimulationResult` containing `n` records.

### `.run_records(n=..., seed=..., expose_provenance=True)`
Returns a list of raw dictionaries. Use `expose_provenance=True` to include `truth_v_call` and other absolute truth fields.

### `.compile(seed=...)`
Returns a `CompiledSimulator` that can be used for streaming:

```python
sim = exp.compile(seed=42)
for outcome in sim.stream(n=10000):
    # Process GenAIRR._engine.Outcome objects
    pass
```

## Reproducibility

GenAIRR is bit-for-bit deterministic. The same pipeline and same seed will produce identical sequences and metadata every time, regardless of whether you are running on Linux, macOS, or Windows.

## Next steps

- [Chain Types](/docs/guides/basics/chain-types) — VDJ vs VJ chains and species coverage
- [Somatic Hypermutation](/docs/guides/options/shm) — Deep dive into SHM models
- [Sequencing Artifacts](/docs/guides/options/artifacts) — Modeling realistic noise
