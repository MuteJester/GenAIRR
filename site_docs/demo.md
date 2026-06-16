---
title: GenAIRR Demo
hide:
  - toc
---

<header class="demo-hero" markdown="0">
  <div class="demo-hero-copy">
    <h1>GenAIRR Demo</h1>
    <p class="lead">Generate annotated AIRR records with ground truth,
    replayable traces, and clonal families from one fluent Python DSL.</p>
    <div class="demo-actions">
      <button class="demo-present-button" type="button" data-demo-presentation-toggle aria-pressed="false">
        <span class="demo-present-button-icon" aria-hidden="true"></span>
        <span data-demo-presentation-label>Presentation mode</span>
      </button>
    </div>
  </div>
  <a class="demo-qr" href="https://mutejester.github.io/GenAIRR/demo.html" aria-label="Open the GenAIRR Demo page">
    <img src="assets/demo-qr.svg" alt="QR code for the GenAIRR Demo page">
    <span>Scan demo</span>
  </a>
</header>

<figure class="demo-diagram demo-architecture-flow" markdown="0">
  <div class="demo-diagram-title">
    <span>From DSL to AIRR record</span>
    <small>the part other simulators usually hide</small>
  </div>
  <div class="demo-arch-rail">
    <section class="demo-node demo-node--truth">
      <span class="demo-node-kicker">1. DSL builder</span>
      <div class="demo-arch-codechips">
        <code>Experiment.on(...)</code>
        <code>.recombine()</code>
        <code>.productive_only()</code>
        <code>.mutate()</code>
      </div>
      <small>Declare the pipeline; no record sampled yet.</small>
    </section>
    <div class="demo-flow-arrow" aria-hidden="true"></div>
    <section class="demo-node demo-node--neutral demo-arch-compile">
      <span class="demo-node-kicker">2. <code>compile()</code></span>
      <strong>typed, ordered, signed pass plan</strong>
      <div class="demo-arch-badges">
        <span>contract set</span>
        <span>plan signature</span>
        <span>live-call hooks</span>
      </div>
    </section>
    <div class="demo-flow-arrow" aria-hidden="true"></div>
    <section class="demo-node demo-node--neutral demo-arch-engine">
      <span class="demo-node-kicker">3. executor</span>
      <strong>passes walk the persistent IR</strong>
      <div class="demo-arch-streams">
        <span><b>trace</b> address + value for each draw</span>
        <span><b>events</b> committed IR changes</span>
        <span><b>hooks</b> committed V/D/J live calls</span>
      </div>
    </section>
    <div class="demo-flow-arrow" aria-hidden="true"></div>
    <section class="demo-node demo-node--record demo-arch-record">
      <span class="demo-node-kicker">4. projection</span>
      <strong>AIRR dict / DataFrame row</strong>
      <div class="demo-record-lanes">
        <div><b>observed</b><code>v_call</code> <code>junction_aa</code></div>
        <div><b>truth</b><code>truth_v_call</code> <code>productive</code></div>
        <div><b>counters</b><code>n_mutations</code> <code>n_v_mutations</code></div>
      </div>
    </section>
  </div>
</figure>

---

## 1. Annotated AIRR records

One fluent chain creates records, truth fields, mutation counters,
and a pandas-ready table.

<div class="demo-run-card" markdown="1">
<div class="demo-cell-label">Code</div>

```python
import GenAIRR as ga

result = (
    ga.Experiment.on("HUMAN_IGH_OGRDB")
      .recombine()
      .productive_only()
      .mutate(model="s5f", rate=0.03)
      .run_records(n=1000, seed=42, expose_provenance=True)
)
df = result.to_dataframe()
```

<div class="demo-cell-label demo-cell-label--output">Output</div>

```text
len(result) = 1000
df.shape    = (1000, 102)
elapsed     = 2.0s
```

</div>

<figure class="demo-diagram demo-airr-record" markdown="0">
  <div class="demo-diagram-title">
    <span>What one generated AIRR row contains</span>
    <small>observed call, committed truth, and counters travel together</small>
  </div>
  <div class="demo-airr-card">
    <div class="demo-airr-sequence">
      <span class="demo-node-kicker">row 0 sequence</span>
      <code>CVKDDGNRGYCSGDSCYGHCCALDYWYFDLW</code>
    </div>
    <div class="demo-airr-columns">
      <section>
        <span class="demo-airr-column-title">Observed fields</span>
        <div><b>v_call</b><code>IGHVF10-G38*04</code></div>
        <div><b>d_call</b><code>IGHD2-15*01</code></div>
        <div><b>j_call</b><code>IGHJ2*01</code></div>
      </section>
      <section>
        <span class="demo-airr-column-title">Ground truth</span>
        <div><b>truth_v_call</b><code>IGHVF10-G38*04</code></div>
        <div><b>truth_d_call</b><code>IGHD2-15*01</code></div>
        <div><b>truth_j_call</b><code>IGHJ2*01</code></div>
      </section>
      <section>
        <span class="demo-airr-column-title">Record state</span>
        <div><b>productive</b><code>True</code></div>
        <div><b>n_mutations</b><code>10</code></div>
        <div><b>columns</b><code>102</code></div>
      </section>
    </div>
  </div>
</figure>

---

## 2. Productive by construction

`productive_only()` changes the sampling support before each draw.

<div class="demo-run-card" markdown="1">
<div class="demo-cell-label">Code</div>

```python
unconstrained = (
    ga.Experiment.on("HUMAN_IGH_OGRDB")
      .recombine()
      .run_records(n=1000, seed=42)
)
constrained = (
    ga.Experiment.on("HUMAN_IGH_OGRDB")
      .recombine()
      .productive_only()
      .run_records(n=1000, seed=42)
)

unconstrained_productive = sum(
    1 for r in unconstrained if r["productive"]
) / len(unconstrained)
constrained_productive = sum(
    1 for r in constrained if r["productive"]
) / len(constrained)

print(f"Without productive_only(): {unconstrained_productive:.1%}")
print(f"With    productive_only(): {constrained_productive:.1%}")
```

<div class="demo-cell-label demo-cell-label--output">Output</div>

```text
Without productive_only(): 18.5%
With    productive_only(): 100.0%
```

</div>

<div class="demo-bars" markdown="0">
  <div class="demo-bar-row">
    <div class="demo-bar-label">Without productive_only()</div>
    <div class="demo-bar-track"><div class="demo-bar-fill demo-bar-fill--warn" style="width: 18.5%;"></div></div>
    <div class="demo-bar-value">18.5%</div>
  </div>
  <div class="demo-bar-row">
    <div class="demo-bar-label">With productive_only()</div>
    <div class="demo-bar-track"><div class="demo-bar-fill demo-bar-fill--accent" style="width: 100%;"></div></div>
    <div class="demo-bar-value">100.0%</div>
  </div>
</div>

---

## 3. Auditable replay

One trace file replays the same record, but refuses a changed
pipeline.

<figure class="demo-diagram demo-replay-flow" markdown="0">
  <div class="demo-diagram-title">
    <span>Replay is gated, not just seeded</span>
    <small>same trace, two outcomes depending on signatures</small>
  </div>
  <div class="demo-replay-rail">
    <section class="demo-node demo-node--truth demo-node--compact">
      <span class="demo-node-kicker">Run</span>
      <strong><code>seed=42</code></strong>
    </section>
    <div class="demo-flow-arrow" aria-hidden="true"></div>
    <section class="demo-node demo-node--neutral demo-node--compact">
      <span class="demo-node-kicker">Trace file</span>
      <strong>12,996 bytes</strong>
    </section>
    <div class="demo-flow-arrow" aria-hidden="true"></div>
    <section class="demo-node demo-node--neutral demo-node--compact">
      <span class="demo-node-kicker">Replay</span>
      <strong><code>replay_from_trace_file()</code></strong>
    </section>
    <div class="demo-splitter" aria-hidden="true">
      <span class="demo-splitter-line"></span>
    </div>
    <div class="demo-replay-results">
      <section class="demo-node demo-node--success">
        <span class="demo-node-kicker">same pipeline + cartridge</span>
        <strong>byte-identical Outcome</strong>
      </section>
      <section class="demo-node demo-node--failure">
        <span class="demo-node-kicker">SHM rate changed</span>
        <strong>ValueError before output</strong>
      </section>
    </div>
  </div>
</figure>

<div class="demo-run-card" markdown="1">
<div class="demo-cell-label">Code</div>

```python
exp = (
    ga.Experiment.on("HUMAN_IGH_OGRDB")
      .recombine().productive_only()
      .mutate(model="s5f", rate=0.03)
)
compiled = exp.compile()

outcome = compiled.simulator.run(seed=42)
trace_file = compiled.simulator.trace_file_from(outcome, seed=42)
trace_file.write_to("demo.trace.json")

from GenAIRR._engine import TraceFile
tf = TraceFile.read_from("demo.trace.json")
replayed = compiled.simulator.replay_from_trace_file(tf, strict=False)
```

<div class="demo-cell-label demo-cell-label--output">Output</div>

```text
trace size = 12,996 bytes
byte-identical fields = sequence, V/D/J calls, junction_aa, n_mutations
```

</div>

<div class="demo-run-card" markdown="1">
<div class="demo-cell-label">Changed code</div>

```python
exp_modified = (
    ga.Experiment.on("HUMAN_IGH_OGRDB")
      .recombine().productive_only()
      .mutate(model="s5f", rate=0.05)
)
compiled_modified = exp_modified.compile()
compiled_modified.simulator.replay_from_trace_file(tf, strict=False)
```

<div class="demo-cell-label demo-cell-label--output">Output</div>

```text
ValueError: replay_from_trace_file: pass plan signature mismatch.
```

</div>

---

## 4. Legacy fixed-size clonal families

One parent recombination can fork into many independently mutated
descendants. This demo uses legacy `expand_clones` for the simple
fixed-size star shape; for new clone benchmarks, see
[`clonal_lineage`](guides/clonal-lineage.md) for BCR trees and
[`clonal_repertoire`](guides/clonal-repertoire.md) for TCR / abundance
repertoires.

<div class="demo-run-card" markdown="1">
<div class="demo-cell-label">Code</div>

```python
result = (
    ga.Experiment.on("HUMAN_IGH_OGRDB")
      .recombine().productive_only()
      .expand_clones(n_clones=10, per_clone=50)
      .mutate(model="s5f", rate=0.02)
      .run_records(seed=42, expose_provenance=True)
)
```

<div class="demo-cell-label demo-cell-label--output">Output</div>

```text
records = 500
parents = 10
per clone = 50
clone 0 SHM counts = 4, 14, 8, 7, 5, ...
validate_families() = ok
```

</div>

<figure class="demo-diagram demo-clone-tree" markdown="0">
  <div class="demo-diagram-title">
    <span>What the clonal result looks like</span>
    <small>10 parent outcomes, 500 descendant AIRR records</small>
  </div>
  <div class="demo-clone-grid">
    <section class="demo-clone-parent">
      <span class="demo-node-kicker">clone 0 parent</span>
      <strong>IGHVF10-G38*04 / IGHD2-15*01 / IGHJ2*01</strong>
      <small><code>result.parents[0]</code> holds the pre-fork IR + trace.</small>
    </section>
    <div class="demo-clone-branches" aria-hidden="true">
      <span></span>
    </div>
    <section class="demo-clone-descendants">
      <div class="demo-clone-leaf demo-clone-leaf--hot">
        <b>clone0_desc0</b><span>4 SHM</span>
      </div>
      <div class="demo-clone-leaf demo-clone-leaf--hotter">
        <b>clone0_desc1</b><span>14 SHM</span>
      </div>
      <div class="demo-clone-leaf demo-clone-leaf--warm">
        <b>clone0_desc2</b><span>8 SHM</span>
      </div>
      <div class="demo-clone-leaf demo-clone-leaf--warm">
        <b>clone0_desc3</b><span>7 SHM</span>
      </div>
      <div class="demo-clone-leaf">
        <b>clone0_desc4</b><span>5 SHM</span>
      </div>
      <div class="demo-clone-more">+45 more descendants</div>
    </section>
    <aside class="demo-clone-stack">
      <div><b>clone 1</b><span>50 descendants</span></div>
      <div><b>clone 2</b><span>50 descendants</span></div>
      <div><b>...</b><span>same structure</span></div>
      <div><b>clone 9</b><span>50 descendants</span></div>
    </aside>
  </div>
  <div class="demo-clone-facts">
    <span><b>Every row:</b> <code>clone_id</code> + <code>parent_id</code></span>
    <span><b>Within clone:</b> shared V/D/J truth</span>
    <span><b>After fork:</b> independent SHM</span>
  </div>
</figure>

---

## Next steps

- [**Quick start**](getting-started/quick-start.md)
- [**Experiment builder**](guides/experiment-builder.md)
- [**Architecture**](architecture/index.md)
- [**GitHub repo**](https://github.com/MuteJester/GenAIRR)
