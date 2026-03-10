---
title: Export Formats
sidebar_label: Export Formats
---

# Export Formats

GenAIRR's `SimulationResult` can be exported to several formats.

## AIRR TSV

```python
from GenAIRR import Experiment

result = Experiment.on("human_igh").run(n=1000, seed=42)
result.to_csv("repertoire.tsv")
```

This writes a tab-separated file with all 47 fields, following the AIRR Community standard.

## FASTA

```python
result.to_fasta("repertoire.fasta")
```

Writes sequences in FASTA format with the record index as the header.

## pandas DataFrame

```python
df = result.to_dataframe()
print(df.shape)    # (1000, 47)
print(df.columns)  # all 47 field names
```

:::note
`to_dataframe()` requires pandas to be installed. Install it with `pip install pandas` or `pip install GenAIRR[all]`.
:::

## Accessing records directly

`SimulationResult` is list-like — you can index, iterate, and check length:

```python
result = Experiment.on("human_igh").run(n=100, seed=42)

# Index
rec = result[0]
print(rec["v_call"])

# Length
print(len(result))  # 100

# Iterate
for rec in result:
    if rec["productive"]:
        print(rec["junction_aa"])
```

## Streaming (no accumulation)

For very large datasets, use `.compile()` + `.stream()` to avoid holding all records in memory:

```python
from GenAIRR import Experiment
from GenAIRR.ops import rate

sim = Experiment.on("human_igh").mutate(rate(0.05, 0.15)).compile(seed=42)

count = 0
for record in sim.stream():
    # Process one record at a time
    count += 1
    if count >= 10000:
        break
```

The stream is an infinite iterator — you control when to stop.
