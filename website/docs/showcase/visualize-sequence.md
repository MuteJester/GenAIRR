---
sidebar_position: 1
title: Sequence Visualization
---

# Sequence Visualization

`visualize_sequence` generates a **standalone HTML file** with an interactive "exploding view" of a simulated AIRR sequence — no browser dependencies, no JavaScript frameworks, just a single self-contained file you can open anywhere.

## Quick Start

```python
from GenAIRR import Experiment, visualize_sequence
from GenAIRR.ops import rate, model

result = Experiment.on("human_igh").mutate(rate(0.05, 0.10), model("s5f")).run(n=1, seed=42)

visualize_sequence(result[0], "sequence.html")
```

Open `sequence.html` in any browser to see the full dissection.

## What You Get

The generated HTML contains a dark-themed, publication-ready visualization with these sections:

### Summary Metrics

<div style={{background: '#181825', border: '1px solid #2a2a40', borderRadius: '10px', padding: '1rem', marginBottom: '1rem'}}>
<div style={{display: 'grid', gridTemplateColumns: 'repeat(4, 1fr)', gap: '0.5rem'}}>
{[
  {label: 'V-GENE', value: 'IGHVF10-G50*04', color: '#4A90D9'},
  {label: 'D-GENE', value: 'IGHD2-21*02', color: '#E8685A'},
  {label: 'J-GENE', value: 'IGHJ4*02', color: '#E8A838'},
  {label: 'MUTATIONS', value: '28 (7.98%)', color: '#DC2626'},
].map(item => (
  <div key={item.label} style={{background: '#1e1e32', borderRadius: '6px', padding: '0.5rem 0.7rem'}}>
    <span style={{fontSize: '0.6rem', textTransform: 'uppercase', letterSpacing: '0.06em', color: '#8888aa', display: 'block'}}>{item.label}</span>
    <span style={{fontSize: '0.85rem', fontWeight: 600, color: item.color}}>{item.value}</span>
  </div>
))}
</div>
</div>

V/D/J gene calls, junction amino acids, total length, mutation count and rate, CDR3 length, and productive status — all at a glance.

### Color-Coded Assembled Sequence Bar

<div style={{background: '#181825', border: '1px solid #2a2a40', borderRadius: '8px', padding: '0.75rem', marginBottom: '1rem'}}>
<div style={{display: 'flex', height: '36px', borderRadius: '5px', overflow: 'hidden', position: 'relative'}}>
  <div style={{width: '72%', background: '#4A90D9', display: 'flex', alignItems: 'center', justifyContent: 'center'}}>
    <span style={{fontSize: '0.7rem', fontWeight: 700, color: '#fff'}}>V</span>
  </div>
  <div style={{width: '5%', background: '#5DC48C', display: 'flex', alignItems: 'center', justifyContent: 'center'}}>
    <span style={{fontSize: '0.6rem', fontWeight: 700, color: '#fff'}}>NP1</span>
  </div>
  <div style={{width: '8%', background: '#E8685A', display: 'flex', alignItems: 'center', justifyContent: 'center'}}>
    <span style={{fontSize: '0.7rem', fontWeight: 700, color: '#fff'}}>D</span>
  </div>
  <div style={{width: '3%', background: '#5DC48C'}}></div>
  <div style={{width: '12%', background: '#E8A838', display: 'flex', alignItems: 'center', justifyContent: 'center'}}>
    <span style={{fontSize: '0.7rem', fontWeight: 700, color: '#fff'}}>J</span>
  </div>
</div>
<div style={{textAlign: 'center', marginTop: '0.4rem'}}>
  <span style={{fontSize: '0.65rem', color: '#9B6FC4', fontFamily: 'monospace'}}>CDR3: CARDVKK*CGGDLPLL</span>
</div>
</div>

Each segment is proportionally sized. The **junction bracket** below shows the CDR3 amino acid translation. **Mutation dots** appear above the bar at their exact positions.

### Exploded Segment Panels

Each of the five segments (V, NP1, D, NP2, J) gets its own panel showing:

| Segment | Panel Contents |
|---------|---------------|
| **V** | Length, mutation count, 3' trim amount, mutation density sparkline, full nucleotide sequence with mutations highlighted in red |
| **NP1 / NP2** | P-prefix, N-addition, P-suffix lengths, P\|N\|P composition bar, colored nucleotide breakdown |
| **D** | Length, 5'/3' trim amounts, inversion status, trim indicators |
| **J** | Length, mutation count, 5' trim amount, mutation density sparkline, full nucleotide sequence |

### Germline Alignment

When the record includes `germline_alignment`, the visualization shows a position-by-position comparison:

<div style={{background: '#181825', border: '1px solid #2a2a40', borderRadius: '6px', padding: '0.75rem', fontFamily: 'monospace', fontSize: '0.7rem', lineHeight: 1.8, marginBottom: '1rem', overflowX: 'auto'}}>
<div style={{color: '#8888aa'}}>
<span style={{color: '#555', marginRight: '0.5rem'}}>1</span>
<span style={{color: '#555'}}>CAGGTG</span><span style={{color: '#DC2626', fontWeight: 700}}>C</span><span style={{color: '#555'}}>AGCTG</span>
</div>
<div>
<span style={{color: '#555', marginRight: '0.5rem'}}>&nbsp;</span>
<span style={{color: '#555'}}>||||||</span><span style={{color: '#DC2626'}}>*</span><span style={{color: '#555'}}>|||||</span>
</div>
<div style={{color: '#8888aa'}}>
<span style={{color: '#555', marginRight: '0.5rem'}}>&nbsp;</span>
<span style={{color: '#555'}}>CAGGTG</span><span style={{color: '#DC2626', fontWeight: 700}}>T</span><span style={{color: '#555'}}>AGCTG</span>
</div>
</div>

Germline on top, `|` for matches, `*` for mismatches (shown in red).

## API Reference

```python
visualize_sequence(record, path, title=None)
```

| Parameter | Type | Description |
|-----------|------|-------------|
| `record` | `dict` | A single AIRR record from `SimulationResult` (e.g., `result[0]`) |
| `path` | `str` or `Path` | Output file path for the HTML |
| `title` | `str`, optional | Custom page title (default: "Sequence Dissection") |

**Returns:** `Path` — the path to the generated file.

## Batch Visualization

Generate one HTML per sequence:

```python
from GenAIRR import Experiment, visualize_sequence
from GenAIRR.ops import rate, model, with_5prime_loss, with_indels

result = (
    Experiment.on("human_igh")
    .mutate(rate(0.03, 0.08), model("s5f"))
    .sequence(with_5prime_loss())
    .observe(with_indels(prob=0.01))
    .run(n=10, seed=42)
)

for i, rec in enumerate(result):
    visualize_sequence(rec, f"seq_{i:04d}.html", title=f"Sequence #{i + 1}")
```

## Supported AIRR Fields

The visualization reads these fields from the record dict:

| Field | Used For |
|-------|----------|
| `sequence` | Full nucleotide display, segment extraction |
| `germline_alignment` | Germline comparison view |
| `v_call`, `d_call`, `j_call` | Gene labels in summary and panels |
| `v_sequence_start/end`, `d_sequence_start/end`, `j_sequence_start/end` | Segment boundaries for color bar |
| `junction_start`, `junction_end`, `junction_aa` | Junction bracket and CDR3 display |
| `mutation_rate`, `mutations` | Mutation statistics, sparklines, nucleotide highlighting |
| `v_trim_3`, `d_trim_5`, `d_trim_3`, `j_trim_5` | Trim indicators in segment panels |
| `np1_p_prefix`, `np1_n_region`, `np1_p_suffix` | NP1 composition bar |
| `np2_p_prefix`, `np2_n_region`, `np2_p_suffix` | NP2 composition bar |
| `d_inverted` | D-gene inversion badge |
| `corruption_5prime`, `corruption_3prime` | Corruption warning badges |
| `productive` | Productive/non-productive badge |

:::tip Self-Contained Output
The HTML files have zero external dependencies — all CSS is inlined. You can email them, embed them in reports, or host them on any static server.
:::
