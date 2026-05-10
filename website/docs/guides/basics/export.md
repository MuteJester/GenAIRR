---
title: Export Formats
sidebar_label: Export Formats
---

# Export Formats

GenAIRR's `SimulationResult` provides several helpers for exporting your simulated sequences and their metadata to industry-standard formats.

## AIRR-Compliant TSV

The most common way to save your results is as an AIRR-compliant tab-separated file.

```python
import GenAIRR as ga

# Run a simulation
result = ga.Experiment.on("human_igh").recombine().run(n=1000, seed=42)

# Save as TSV (default)
result.to_csv("repertoire.tsv")

# Save as CSV
result.to_csv("repertoire.csv", sep=",")
```

This writes a file with **~70 ground-truth fields** per sequence.

### 1-Based Coordinates (AIRR Strict)

By default, GenAIRR uses **0-based half-open** coordinates (Python convention). To export using the **1-based inclusive** coordinates required by the AIRR Rearrangement specification, pass `airr_strict=True`:

```python
# Convert all coordinate fields to 1-based inclusive for AIRR tooling
result.to_csv("repertoire.tsv", airr_strict=True)
```

## FASTA

To export just the nucleotide sequences for use with tools like BLAST or IgBLAST:

```python
result.to_fasta("repertoire.fasta")
```

The sequence ID in the FASTA header will match the `sequence_id` field in the AIRR metadata (e.g., `>seq0`).

## pandas DataFrame

If you have `pandas` installed, you can convert the entire result set into a DataFrame for easy analysis:

```python
df = result.to_dataframe()

print(df.shape)     # (1000, 69)
print(df["v_call"].value_counts())
```

:::note
`to_dataframe()` requires pandas. Install it via `pip install pandas` or `pip install GenAIRR[all]`.
:::

## Accessing Records Directly

`SimulationResult` is a list-like container. You can index it, iterate over it, and check its length. Each record is returned as a plain Python `dict`:

```python
# Get the first record
rec = result[0]
print(rec["sequence"])

# Iterate through results
for rec in result:
    if rec["productive"]:
        print(f"Productive junction: {rec['junction_aa']}")
```

## Memory-Efficient Streaming

For very large datasets (e.g., millions of sequences), you should avoid loading all results into memory. Use `.stream_records()` to lazily generate sequences one at a time:

```python
exp = ga.Experiment.on("human_igh").recombine().mutate(count=10)

# Stream 1,000,000 records without using massive amounts of RAM
for record in exp.stream_records(n=1_000_000, seed=42):
    # 'record' is a dict containing the AIRR fields
    process(record)
```

The stream will yield exactly `n` records (if provided) and then stop.

## Next steps

- [Understanding Output](/docs/getting-started/interpreting-results) — Detailed guide to the ~70 output fields
- [Experiment DSL](/docs/guides/basics/experiment-dsl) — Learn how to customize your simulation pipeline
