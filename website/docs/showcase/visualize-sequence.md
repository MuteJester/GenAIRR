---
sidebar_position: 1
title: Sequence Visualization
---

# Sequence Visualization

`visualize_sequence` generates a **standalone HTML file** with an interactive "exploding view" of a simulated AIRR sequence. It allows you to visually dissect the sequence, showing color-coded segments, mutation sites, and detailed alignments.

## Quick Start

```python
import GenAIRR as ga
from GenAIRR.utilities.visualize import visualize_sequence

# Run a simulation
result = ga.Experiment.on("human_igh").recombine().mutate(count=15).run(n=1, seed=42)

# Generate the visualization
visualize_sequence(result[0], "sequence.html")
```

Open `sequence.html` in any browser to see the full dissection.

## What You Get

The generated HTML contains a professional, publication-ready visualization with several key sections:

### 1. Summary Metrics
At the top, you'll see a high-level overview of the sequence, including V/D/J gene calls, total mutations, mutation rate, CDR3 length, and productivity status.

### 2. The Exploding View
The core of the visualization is a color-coded bar representing the assembled sequence:
*   **V Segment** (Blue)
*   **NP1 Region** (Green)
*   **D Segment** (Red)
*   **NP2 Region** (Green)
*   **J Segment** (Yellow)

Mutation dots appear above the bar, and the **Junction Bracket** below shows the CDR3 amino acid translation.

### 3. Segment Panels
Each biological segment gets its own detailed panel. For the V, D, and J segments, the panel includes:
*   **Length & Identity:** The number of bases and the percentage match to the germline.
*   **Trimming:** Bases nibbled from the 5' or 3' ends.
*   **Mutations:** A density sparkline and the full nucleotide sequence with mutations highlighted.
*   **CIGAR String:** The compressed alignment operations (e.g., `301M`).

### 4. Germline Alignment
If the record includes a `germline_alignment` (which all GenAIRR records do by default), the visualization provides a position-by-position comparison between the simulated sequence and its biological origin.

## API Reference

```python
visualize_sequence(record, path, title=None)
```

| Parameter | Type | Description |
|-----------|------|-------------|
| `record` | `dict` | A single AIRR record (e.g., `result[0]`). |
| `path` | `str` or `Path` | Output file path for the HTML. |
| `title` | `str`, optional | Custom page title (default: "Sequence Dissection"). |

## Batch Visualization

You can easily generate visualizations for an entire batch of sequences:

```python
import GenAIRR as ga
from GenAIRR.utilities.visualize import visualize_sequence

result = ga.Experiment.on("human_igh").recombine().mutate(count=(5, 15)).run(n=10, seed=42)

for i, rec in enumerate(result):
    visualize_sequence(rec, f"visualizations/seq_{i}.html", title=f"Record #{i}")
```

:::tip Zero Dependencies
The generated HTML files are completely self-contained. All styles and scripts are inlined, so you can share the files easily without worrying about external assets.
:::

## Supported Fields

The visualization is optimized for the **~70 ground-truth fields** produced by the new Rust engine, including:
*   `v_identity`, `d_identity`, `j_identity`
*   `v_cigar`, `d_cigar`, `j_cigar`
*   `v_sequence_start/end`, `d_sequence_start/end`, `j_sequence_start/end`
*   `junction_aa`, `junction_start/end`
*   `n_mutations`, `mutation_rate`
*   `germline_alignment`
