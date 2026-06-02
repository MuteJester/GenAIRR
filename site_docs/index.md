---
title: GenAIRR
hide:
  - navigation
  - toc
---

<p class="eyebrow accent">// Synthetic immune repertoires with absolute ground truth</p>

# Simulate immune receptor sequences. Know exactly what the engine did.

<p class="lead">GenAIRR is a high-performance simulator for adaptive immune
receptor repertoires (BCR + TCR). A Rust simulation kernel drives a fluent
Python DSL; every record comes back annotated with by-construction truth
(V/D/J calls, junction, productive flag, identity, mutation counts)
derived from the persistent intermediate representation, not inferred by
an aligner.</p>

<div class="cta-row">
<a class="cta" href="getting-started/quick-start.html"><span class="cta-eyebrow">// 5 min</span><span class="cta-title">First simulation <span class="cta-arrow">&rarr;</span></span><span class="cta-sub">Install, simulate 1,000 productive heavy chains, inspect AIRR records.</span></a><a class="cta" href="concepts/reference-cartridge.html"><span class="cta-eyebrow">// Build</span><span class="cta-title">Reference cartridges <span class="cta-arrow">&rarr;</span></span><span class="cta-sub">Four typed planes. Build from FASTA. Estimate biology from your own AIRR data.</span></a><a class="cta" href="validation/validate-records.html"><span class="cta-eyebrow">// Verify</span><span class="cta-title">Validate AIRR records <span class="cta-arrow">&rarr;</span></span><span class="cta-sub">One call answers: is every reported field internally consistent with the outcome?</span></a><a class="cta" href="reference/index.html"><span class="cta-eyebrow">// Lookup</span><span class="cta-title">API reference <span class="cta-arrow">&rarr;</span></span><span class="cta-sub">Every public symbol: Experiment, SimulationResult, builders, validators.</span></a>
</div>

## What it's for. A DSL for AIRR. An engine that knows the answer.

<div class="feat-grid" markdown>

<div class="feat" markdown>
  <p class="num">// 01</p>
  <h3>Benchmark</h3>
  <p>Compare alignment, clustering, or annotation tools against a ground truth
  the engine emitted by construction, not by an oracle aligner that can be
  wrong.</p>
</div>

<div class="feat" markdown>
  <p class="num">// 02</p>
  <h3>Null models</h3>
  <p>Generate biologically grounded null repertoires with controllable
  productive fraction, SHM rate, V-gene usage, and clonal structure for
  statistical hypothesis testing.</p>
</div>

<div class="feat" markdown>
  <p class="num">// 03</p>
  <h3>Phenomena lab</h3>
  <p>Switch any biological mechanism on or off (P/N additions, D inversion,
  receptor revision, targeted SHM) and observe the downstream effect on
  what aligners report.</p>
</div>

</div>

<div class="engine-strip" markdown>
  <span class="engine-tag">// Engine</span>
  <span>Rust kernel</span>
  <span class="sep">·</span>
  <span>23 species, 106 bundled configs</span>
  <span class="sep">·</span>
  <span>Constraint-aware sampling</span>
  <span class="sep">·</span>
  <span>Byte-stable replay</span>
</div>

## What GenAIRR does. Six capabilities. One composable engine.

<div class="cap-grid" markdown>

<div class="cap-card cap-recombine" markdown>
  <p class="cap-num">// V(D)J · 01</p>
  <h3>Recombine</h3>
  <p>Sample alleles, trim, fill NP1/NP2, assemble. Authored typed empirical
  models per cartridge plane.</p>
</div>

<div class="cap-card cap-mutate" markdown>
  <p class="cap-num">// SHM · 02</p>
  <h3>Mutate</h3>
  <p>Uniform or S5F context-aware SHM. Per-segment + per-V-subregion rate
  targeting (FWR/CDR).</p>
</div>

<div class="cap-card cap-clone" markdown>
  <p class="cap-num">// Lineage · 03</p>
  <h3>Expand clones</h3>
  <p>Parent rearrangement forks N descendants. SHM accumulates
  per-descendant after the fork.</p>
</div>

<div class="cap-card cap-corrupt" markdown>
  <p class="cap-num">// Library · 04</p>
  <h3>Corrupt</h3>
  <p>Primer trimming, structural indels, PCR errors, N-base injection. Each
  knob tunable per workload.</p>
</div>

<div class="cap-card cap-constrain" markdown>
  <p class="cap-num">// Contracts · 05</p>
  <h3>Constrain</h3>
  <p>Productive-only sampling at compile time. The engine never proposes
  out-of-frame or stop-bearing candidates.</p>
</div>

<div class="cap-card cap-replay" markdown>
  <p class="cap-num">// Replay · 06</p>
  <h3>Replay</h3>
  <p>Trace every random draw. Byte-identical replay across runs, platforms,
  and one-knob-changed counterfactuals.</p>
</div>

</div>

```python
import GenAIRR as ga

result = (
    ga.Experiment.on("human_igh")
      .recombine()
      .expand_clones(n_clones=50, per_clone=20)
      .mutate(rate=0.05)
      .pcr_amplify(count=(0, 3))
      .ambiguous_base_calls(count=(0, 2))
      .productive_only().run_records(seed=42)
)

result.to_tsv("repertoire.tsv")
```

## Install. One command. No compiler.

<div class="install-band" markdown>
`pip install GenAIRR`
<span class="mute">Pre-built wheels for Linux (x86_64, aarch64), macOS (Intel + Apple Silicon), and Windows (x64). Python 3.9+. No Rust toolchain needed.</span>
</div>

[See the full quick start →](getting-started/quick-start.md){ .btn .btn-primary }
[Learn by simulating one molecule →](getting-started/index.md){ .btn .btn-ghost }

---

## Choose your path

Five paths through the docs, organised by what you came here to
do. Each one is a short curriculum. Pick the one that matches
your task.

| If you want to ... | Start here |
|---|---|
| **Simulate sequences** | [Quick start](getting-started/quick-start.md) → [The Experiment builder](guides/experiment-builder.md) |
| **Build a reference cartridge** | [Reference cartridge concept](concepts/reference-cartridge.md) → [Build a reference cartridge](guides/build-reference-cartridge.md) |
| **Get reproducible / validated output** | [Validation hub](validation/index.md) → [Trace, replay, reproducibility](guides/trace-replay.md) |
| **Understand the biological mechanisms** | [Recombination + junction biology](guides/recombination-junction.md), [SHM and mutation targeting](guides/shm-targeting.md), [Corruption + sequencing artefacts](guides/corruption-sequencing.md), [Clonal families](guides/clonal-families.md) |
| **Contribute to GenAIRR** | [Architecture (Contributor)](architecture/index.md) |

The **[Choose your path](learn.md)** page expands each of these
into a short ordered reading list.

---

<p class="eyebrow">// Citation</p>

If GenAIRR helps your research, please cite: Konstantinovsky T, Peres A,
Polak P, Yaari G. *An unbiased comparison of immunoglobulin sequence
aligners.* Briefings in Bioinformatics. 2024;25(6):bbae556.
[doi:10.1093/bib/bbae556](https://doi.org/10.1093/bib/bbae556)
