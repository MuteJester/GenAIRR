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

# Generate 1,000 productive human heavy-chain sequences. Every sequence
# comes back with the full AIRR-format annotation block — gene calls,
# junction, productive flag, identity, mutation counts.
result = (
    ga.Experiment.on("human_igh")
      .recombine()
      .run_records(n=1000, seed=42, respect=ga.productive())
)

# `result` is a SimulationResult — list-like over AIRR record dicts.
# Each dict has the 50+ standard AIRR fields per row.
len(result)                 # 1000
rec = result[0]

rec["sequence"]             # 'gaggtgcagctggtggagtctgggggaggc...' (nucleotide)
rec["sequence_aa"]          # 'EVQLVESGGGLVQPGGSLRLSCSAS...'      (translated)
rec["locus"]                # 'IGH'
rec["v_call"]               # 'IGHVF10-G38*04'   (comma-separated if the call ties)
rec["d_call"]               # 'IGHD2-15*01'
rec["j_call"]               # 'IGHJ2*01'
rec["junction_aa"]          # 'CVKDDGNRGYCSGGSCYGRCCALDYWYFDLW'
rec["productive"]           # True
rec["v_identity"]           # 1.0  (matches/total over the V segment)
rec["n_mutations"]          # 0

# Export in any of the standard formats. TSV/FASTA/FASTQ are dependency-free;
# to_dataframe() needs pandas (pip install GenAIRR[all]).
result.to_tsv("repertoire.tsv")        # AIRR-spec TSV (50+ columns)
result.to_fasta("sequences.fasta")     # FASTA with v_call/j_call in the headers
result.to_fastq("sequences.fastq")     # FASTQ with illumina-shaped quality scores
df = result.to_dataframe()             # one row per record, AIRR columns
```

`Experiment.on(...)` accepts **a config-name string** (e.g. `"human_igh"`, `"mouse_tcrb"`), **a `DataConfig`** loaded from the bundled species pickles, or **a `RefDataConfig`** for [custom reference data](#custom-reference-data). `respect=ga.productive()` is the constraint-aware bundle — covered in the next section. Drop it to allow non-productive sequences (~30% of records will then have stop codons in the junction).

> See the full walkthrough in the docs: [Quick Start](https://mutejester.github.io/GenAIRR/docs/getting-started/quick-start) · [Interpreting Results](https://mutejester.github.io/GenAIRR/docs/getting-started/interpreting-results)

---

## Constraint-aware sampling

GenAIRR's signature feature is **constraint-aware sampling**: contracts that prune the candidate distribution at sample time, not retries after the fact. The canonical bundle is `productive()` (in-frame junction + no stop codons + V/J anchors preserved):

```python
import GenAIRR as ga

# Every sequence is productive by construction. No retry loops, no
# post-hoc filtering — the engine only ever picks NP lengths, NP bases,
# and mutation substitutions that satisfy the bundle.
result = (
    ga.Experiment.on("human_igh")
      .recombine()
      .run_records(n=1000, seed=42, respect=ga.productive())
)
assert all(rec["productive"] for rec in result)
```

> Docs: [Productive sequences](https://mutejester.github.io/GenAIRR/docs/guides/options/productive)

### Strict vs permissive mode

By default, if a contract can't admit any candidate at a sampling step the runtime falls back to unconstrained sampling and the run continues. Pass `strict=True` to surface the failure as an exception instead — useful for catching unsatisfiable plans early during development:

```python
import GenAIRR as ga

try:
    ga.Experiment.on("human_igh").recombine().run_records(
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

# Same seed → byte-identical records across runs and platforms.
a = ga.Experiment.on("human_igh").recombine().run_records(n=100, seed=42)
b = ga.Experiment.on("human_igh").recombine().run_records(n=100, seed=42)
assert a[0]["sequence"] == b[0]["sequence"]

# `n` runs use seeds [seed, seed+1, ..., seed+n-1] so consecutive
# batches stitch together by offsetting the starting seed.
batch_a = ga.Experiment.on("human_igh").recombine().run_records(n=100, seed=0)
batch_b = ga.Experiment.on("human_igh").recombine().run_records(n=100, seed=100)
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
    result = compiled.run_records(n=100, seed=batch * 100)
    result.to_tsv(f"batch_{batch:02d}.tsv")
```

---

## What you get back

`.run_records(...)` returns a `SimulationResult` — a list-like wrapper around a batch of AIRR record dicts:

| Method / attribute | Returns | Description |
|-----|-----|-----|
| `len(result)` | `int` | Number of records in the batch. |
| `result[i]` | `dict` | The i-th AIRR record. Standard 0-based indexing + slicing. |
| `for rec in result:` | iterates `dict`s | Records in `[seed, seed+1, …, seed+n-1]` order. |
| `result.records` | `list[dict]` | The underlying list. Mutate-through is fine. |
| `result.to_tsv(path, *, airr_strict=False)` | — | AIRR-format TSV. `airr_strict=True` converts coordinates to 1-based-inclusive per spec. |
| `result.to_csv(path, *, airr_strict=False)` | — | Comma-separated. Same options as `to_tsv`. |
| `result.to_fasta(path, *, prefix="seq")` | — | FASTA. Headers include `v_call` and `j_call`. |
| `result.to_fastq(path, *, quality="illumina", **kw)` | — | FASTQ. Quality models: `"illumina"` (smoothed trapezoid) or `"constant"`. |
| `result.to_dataframe(*, airr_strict=False)` | `pandas.DataFrame` | One row per record. Requires pandas (`pip install GenAIRR[all]`). |
| `result.outcomes` | `list[Outcome] \| None` | The underlying `Outcome` objects, for advanced introspection (see below). |

Each record dict has 50+ AIRR fields. The most commonly used:

| Field | Example value | Description |
|-----|-----|-----|
| `sequence` | `'gaggtgcagctggtg…'` | Assembled nucleotide sequence (uppercase + lowercase corruption markers). |
| `sequence_aa` | `'EVQLVESGGG…'` | Codon-rail translation. Stops emit `*`, ambiguous codons emit `X`. |
| `locus` | `'IGH'` | Locus code derived from `v_call` / `j_call`. |
| `v_call` / `d_call` / `j_call` | `'IGHV3-23*01'` | Gene calls. Comma-separated tie set when the evidence walker can't disambiguate. |
| `junction` / `junction_aa` | `'TGC…GAC'` / `'CAR…D'` | Junction nucleotide + AA. AA includes the V Cys (anchor) through J W/F+3. |
| `productive` | `True` / `False` / `None` | In-frame junction AND no stop codons AND anchors preserved. `None` when undefined (e.g. junction not present). |
| `v_identity` / `d_identity` / `j_identity` | `0.987` | Match rate over each segment's CIGAR M/D ops. |
| `v_cigar` / `d_cigar` / `j_cigar` | `'17D279M'` | CIGAR strings. Only M/I/D ops are emitted — no soft-clips. |
| `n_mutations` / `n_pcr_errors` / `n_quality_errors` / `n_indels` | `4` / `0` / `2` / `1` | Per-record error counts from the trace. |

The full schema (plus the `*_sequence_start/end`, `*_alignment_start/end`, `*_germline_start/end` coordinate fields, `vj_in_frame`, `stop_codon`, `rev_comp`, and others) is documented at [Interpreting Results](https://mutejester.github.io/GenAIRR/docs/getting-started/interpreting-results).

### Advanced: full pipeline state via `Outcome`

When you need step-by-step IR history or the raw trace of every random draw — debugging an engine bug, building a custom alignment tool, replaying a specific seed — use `.run()` instead of `.run_records()`. It returns a list of `Outcome` objects that carry the full pipeline state:

| Accessor | Returns | Description |
|----------|---------|-------------|
| `outcome.final_simulation()` | `Simulation` | End-of-pipeline IR snapshot. |
| `outcome.revision(i)` | `Simulation` | IR after the i-th pass — full step-by-step history. |
| `outcome.revision_after(name)` | `Simulation \| None` | First revision produced by the named pass. |
| `outcome.pass_names()` | `list[str]` | Names of every pass that ran, in order. |
| `outcome.trace()` | `Trace` | Addressed log of every random draw. |

Each `Simulation` exposes `len(sim)` (pool length), `sim.bases() → bytes`, `sim.regions() → list[Region]`, `sim.germline_position(i)`, `sim.v_allele_id() / .d_allele_id() / .j_allele_id()`. Each `Region` carries `segment` (`"V"`/`"D"`/`"J"`/`"NP1"`/`"NP2"`), `start`/`end`/`len()`, `frame_phase`, and `amino_acids() → bytes` (codon-rail translation, including stop markers and ambiguous codons).

`outcome.trace()` supports `find(address)`, `prefix_query(prefix)`, and `prefix_count(prefix)` — every random draw is keyed by a hierarchical address (`"sample_allele.v"`, `"np.np1.length"`, `"np.np1.bases[3]"`, …). This is the same trace the engine uses internally for replay determinism.

`.run_records(...)` also exposes these via `result.outcomes[i]` — so you can have both the AIRR records *and* the deep introspection from a single call.

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

ga.Experiment.on("mouse_igh").recombine().run_records(n=500)
ga.Experiment.on("rabbit_tcrb").recombine().run_records(n=500)
ga.Experiment.on("rhesus_igk").recombine().run_records(n=500)
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

result = ga.Experiment.on(cfg).recombine().run_records(n=100, seed=42)
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
