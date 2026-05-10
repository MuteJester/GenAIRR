<h1 align="center">GenAIRR</h1>

<p align="center">
  <b>Synthetic Adaptive Immune Receptor Repertoire Generator</b>
</p>

<p align="center">
  <a href="https://pypi.org/project/GenAIRR/"><img src="https://img.shields.io/pypi/v/GenAIRR.svg?logo=pypi&logoColor=white" alt="PyPI"></a>
  <a href="https://github.com/MuteJester/GenAIRR/actions/workflows/test.yml"><img src="https://github.com/MuteJester/GenAIRR/actions/workflows/test.yml/badge.svg" alt="Tests"></a>
  <a href="https://pypi.org/project/GenAIRR/"><img src="https://img.shields.io/pypi/pyversions/GenAIRR.svg?logo=python&logoColor=white" alt="Python"></a>
  <a href="https://github.com/MuteJester/GenAIRR/blob/master/LICENSE"><img src="https://img.shields.io/github/license/MuteJester/GenAIRR" alt="License"></a>
</p>

<p align="center">
  High-performance BCR and TCR sequence simulation with full ground-truth annotations.<br/>
  Rust kernel &middot; 23 species &middot; constraint-aware sampling &middot; cross-platform wheels
</p>

<p align="center">
  <a href="https://mutejester.github.io/GenAIRR/"><b>📖 Documentation</b></a>
</p>

---

## Installation

```bash
pip install GenAIRR
```

GenAIRR ships as a single wheel that bundles both the Python API and the Rust simulation kernel — no extra packages, no compiler needed. Pre-built wheels are published for **Linux** (x86_64, aarch64), **macOS** (Intel + Apple Silicon), and **Windows** (x64), supporting Python **3.9+**.

Building from source needs a stable Rust toolchain (`rustup install stable`) — see [CONTRIBUTING.md](CONTRIBUTING.md).

---

## Quick Start

```python
import GenAIRR as ga

# Generate 1,000 human heavy-chain sequences via standard V(D)J recombination.
outcomes = ga.Experiment.on("human_igh").recombine().run(n=1000, seed=42)

# Each item is an Outcome wrapping the full pipeline state.
o = outcomes[0]
sim = o.final_simulation()

sim.bases()          # b'gaggtgcagctggtg...' — the assembled sequence
sim.regions()        # [<Region V [0..296) frame_phase=0>,
                     #  <Region NP1 [296..309) frame_phase=2>,
                     #  <Region D [309..336) frame_phase=0>,
                     #  <Region NP2 [336..342) frame_phase=0>,
                     #  <Region J [342..388) frame_phase=0>]
sim.v_allele_id()    # 146  (index into the V pool of the active refdata)
sim.d_allele_id()    # 5
sim.j_allele_id()    # 1

# The trace records every random draw the engine made.
o.trace().find("np.np1.length").value     # 13
o.pass_names()                            # ['sample_allele.v', 'sample_allele.d',
                                          #  'sample_allele.j', 'trim.v_3', ...]
```

`Experiment.on(...)` accepts **a config-name string** (e.g. `"human_igh"`, `"mouse_tcrb"`), **a `DataConfig`** loaded from the bundled species pickles, or **a `RefDataConfig`** for [custom reference data](#custom-reference-data).

> See the full walkthrough in the docs: [Quick Start](https://mutejester.github.io/GenAIRR/docs/getting-started/quick-start) · [Interpreting Results](https://mutejester.github.io/GenAIRR/docs/getting-started/interpreting-results)

---

## Constraint-aware sampling

GenAIRR's signature feature is **constraint-aware sampling**: contracts that prune the candidate distribution at sample time, not retries after the fact. The canonical bundle is `productive()` (in-frame junction + no stop codons + V/J anchors preserved):

```python
import GenAIRR as ga

# Every sequence is productive by construction. No retry loops, no
# post-hoc filtering — the engine only ever picks NP lengths, NP bases,
# and mutation substitutions that satisfy the bundle.
outcomes = (
    ga.Experiment.on("human_igh")
      .recombine()
      .run(n=1000, seed=42, respect=ga.productive())
)
```

> Docs: [Productive sequences](https://mutejester.github.io/GenAIRR/docs/guides/options/productive)

### Strict vs permissive mode

By default, if a contract can't admit any candidate at a sampling step the runtime falls back to unconstrained sampling and the run continues. Pass `strict=True` to surface the failure as an exception instead — useful for catching unsatisfiable plans early during development:

```python
import GenAIRR as ga

try:
    ga.Experiment.on("human_igh").recombine().run(
        n=10, seed=42, respect=ga.productive(), strict=True
    )
except ga.StrictSamplingError as e:
    pass_name, address, reason = e.args
    # pass_name e.g. "generate_np.np1", address e.g. "np.np1.length",
    # reason in {"empty_admissible_support", "support_unavailable", ...}
    print(f"{pass_name} could not satisfy the contract at {address}: {reason}")
```

---

## Reproducibility

```python
import GenAIRR as ga

# Same seed → byte-identical outcomes across runs and platforms.
a = ga.Experiment.on("human_igh").recombine().run(n=100, seed=42)
b = ga.Experiment.on("human_igh").recombine().run(n=100, seed=42)
assert a[0].final_simulation().bases() == b[0].final_simulation().bases()

# `n` runs use seeds [seed, seed+1, ..., seed+n-1] so consecutive
# batches stitch together by offsetting the starting seed.
batch_a = ga.Experiment.on("human_igh").recombine().run(n=100, seed=0)
batch_b = ga.Experiment.on("human_igh").recombine().run(n=100, seed=100)
# batch_a[50] is byte-equal to a one-off run at seed=50.
```

> Docs: [Reproducibility](https://mutejester.github.io/GenAIRR/docs/guides/options/reproducibility)

---

## Compile once, run many times

For a hot loop, `compile()` once and reuse the plan. Contracts (`respect=`) are baked into the compiled plan, so they only need to be passed once:

```python
import GenAIRR as ga

compiled = (
    ga.Experiment.on("human_igk")
      .recombine()
      .compile(respect=ga.productive())
)

# Run 10 batches of 100, seeded so they don't overlap.
for batch in range(10):
    outcomes = compiled.run(n=100, seed=batch * 100)
```

---

## What you get back

Every `outcome` in the returned list is an `Outcome` with:

| Accessor | Returns | Description |
|----------|---------|-------------|
| `outcome.final_simulation()` | `Simulation` | The end-of-pipeline IR snapshot. |
| `outcome.revision(i)` | `Simulation` | The IR after the i-th pass — full step-by-step history. |
| `outcome.revision_after(name)` | `Simulation \| None` | First revision produced by the named pass. |
| `outcome.pass_names()` | `list[str]` | Names of every pass that ran, in order. |
| `outcome.trace()` | `Trace` | Addressed log of every random draw. |

Each `Simulation` exposes `len(sim)` (pool length), `sim.bases() → bytes`, `sim.regions() → list[Region]`, `sim.germline_position(i)`, `sim.v_allele_id() / .d_allele_id() / .j_allele_id()`. Each `Region` carries `segment` (`"V"`/`"D"`/`"J"`/`"NP1"`/`"NP2"`), `start`/`end`/`len()`, `frame_phase`, and `amino_acids() → bytes` (codon-rail translation, including stop markers and ambiguous codons).

`outcome.trace()` supports `find(address)`, `prefix_query(prefix)`, and `prefix_count(prefix)` — every random draw the engine made is keyed by a hierarchical address (`"sample_allele.v"`, `"np.np1.length"`, `"np.np1.bases[3]"`, ...). This is the same trace the engine uses internally for replay determinism, so your downstream tooling sees exactly what the kernel saw.

> Docs: [Simulation Pipeline](https://mutejester.github.io/GenAIRR/docs/concepts/simulation-pipeline) · [Metadata Accuracy](https://mutejester.github.io/GenAIRR/docs/concepts/metadata-accuracy) · [Interpreting Results](https://mutejester.github.io/GenAIRR/docs/getting-started/interpreting-results)

---

## Supported Species & Chains

GenAIRR ships with **106 built-in configurations** covering 23 species (sourced from IMGT and OGRDB).

```python
import GenAIRR as ga
print(ga.list_configs())  # all available configs
```

| Species | BCR | TCR |
|---------|-----|-----|
| Human | IGH, IGK, IGL | TCRA, TCRB, TCRD, TCRG |
| Mouse | IGH, IGK, IGL | TCRA, TCRB, TCRD, TCRG |
| Rat | IGH, IGK, IGL | &mdash; |
| Rabbit | IGH, IGK, IGL | TCRA, TCRB, TCRD, TCRG |
| Dog | IGH, IGK, IGL | TCRA, TCRB, TCRD, TCRG |
| Cat | IGK, IGL | TCRA, TCRB, TCRD, TCRG |
| Rhesus | IGH, IGK, IGL | TCRA, TCRB, TCRD, TCRG |

<details>
<summary>All 23 species</summary>

Alpaca, Cat, Chicken, Cow, Cynomolgus, Dog, Dromedary, Ferret, Goat, Gorilla,
Horse, Human, Mouse (generic + C57BL/6J), Pig, Platypus, Rabbit, Rat, Rhesus,
Salmon, Sheep, Trout, Zebrafish.

</details>

```python
import GenAIRR as ga

ga.Experiment.on("mouse_igh").recombine().run(n=500)
ga.Experiment.on("rabbit_tcrb").recombine().run(n=500)
ga.Experiment.on("rhesus_igk").recombine().run(n=500)
```

> Docs: [Choosing a config](https://mutejester.github.io/GenAIRR/docs/getting-started/choosing-config) · [Chain types](https://mutejester.github.io/GenAIRR/docs/guides/basics/chain-types)

---

## Custom reference data

For non-builtin alleles (custom IMGT pulls, in-house references, etc.) you can build a `RefDataConfig` directly and pass it to `Experiment.on(...)`:

```python
import GenAIRR as ga

cfg = ga.RefDataConfig.vj()
cfg.add_v_allele("v_custom*01", "v_custom", b"GAAGTACAGCTGGTGCAG...", anchor=288)
cfg.add_v_allele("v_custom*02", "v_custom", b"GAAGTACAGCTAGTGCAG...", anchor=288)
cfg.add_j_allele("j_custom*01", "j_custom", b"TGGGGCCAAGGG...",       anchor=10)

outcomes = ga.Experiment.on(cfg).recombine().run(n=100, seed=42)
```

`RefDataConfig.vdj()` builds a heavy-chain-shaped refdata (with a D pool); `add_d_allele(...)` populates it. Anchors are 0-based offsets of the V Cys / J W or F codon's first base, used to keep the junction frame-aligned during recombination.

---

## Key Features

- **Rust simulation kernel** &mdash; persistent IR with full revision history, addressed-trace introspection, `cargo test`-grade unit coverage.
- **Constraint-aware sampling** &mdash; contracts prune candidate distributions at sample time so productive sequences come out of the engine by construction; no retry loops.
- **Strict-mode opt-in** &mdash; surface unsatisfiable plans as `StrictSamplingError` instead of silently relaxing the bundle.
- **Deterministic seeds** &mdash; same seed reproduces every byte of the pool and every entry of the trace, across runs and platforms.
- **Full revision history** &mdash; `outcome.revision(i)` exposes the IR after each pass for fine-grained debugging.
- **Addressed trace** &mdash; every random draw is keyed by a hierarchical string (`"np.np1.bases[3]"`) and survives end-to-end into the returned `Outcome`.
- **23 species, 106 configs** &mdash; built-in IMGT + OGRDB reference pickles ship with the wheel.
- **Zero mandatory Python dependencies** &mdash; one wheel, everything in the box.

---

## Optional Extras

```bash
pip install GenAIRR[all]          # numpy, scipy, graphviz, tqdm
pip install GenAIRR[dataconfig]   # numpy + scipy (custom DataConfig analysis)
pip install GenAIRR[viz]          # graphviz
```

---

## Documentation

The full documentation site is at **[mutejester.github.io/GenAIRR](https://mutejester.github.io/GenAIRR/)**. Useful starting points:

- **Getting started** — [Quick Start](https://mutejester.github.io/GenAIRR/docs/getting-started/quick-start) · [Choosing a Config](https://mutejester.github.io/GenAIRR/docs/getting-started/choosing-config) · [Interpreting Results](https://mutejester.github.io/GenAIRR/docs/getting-started/interpreting-results)
- **Concepts** — [Simulation Pipeline](https://mutejester.github.io/GenAIRR/docs/concepts/simulation-pipeline) · [Metadata Accuracy](https://mutejester.github.io/GenAIRR/docs/concepts/metadata-accuracy)
- **Guides** — [Experiment DSL](https://mutejester.github.io/GenAIRR/docs/guides/basics/experiment-dsl) · [Chain Types](https://mutejester.github.io/GenAIRR/docs/guides/basics/chain-types) · [Export](https://mutejester.github.io/GenAIRR/docs/guides/basics/export)
- **Options** — [Productive](https://mutejester.github.io/GenAIRR/docs/guides/options/productive) · [Reproducibility](https://mutejester.github.io/GenAIRR/docs/guides/options/reproducibility) · [SHM](https://mutejester.github.io/GenAIRR/docs/guides/options/shm) · [Biology](https://mutejester.github.io/GenAIRR/docs/guides/options/biology) · [Artifacts](https://mutejester.github.io/GenAIRR/docs/guides/options/artifacts)

---

## Citing GenAIRR

If GenAIRR is useful in your research, please cite:

> Konstantinovsky T, Peres A, Polak P, Yaari G. An unbiased comparison of immunoglobulin sequence aligners. *Briefings in Bioinformatics*. 2024;25(6):bbae556. [doi:10.1093/bib/bbae556](https://doi.org/10.1093/bib/bbae556)

---

## Contributing

Contributions are welcome. See [CONTRIBUTING.md](CONTRIBUTING.md) for development setup and guidelines.

## License

GPL-3.0. See [LICENSE](LICENSE).
