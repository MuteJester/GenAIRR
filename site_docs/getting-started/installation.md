# Installation

<p class="lead">GenAIRR ships as a single pre-built wheel that bundles
both the Python API and the Rust simulation kernel. No compiler is
required on any supported platform.</p>

## Requirements

- **Python 3.9 or newer** (3.9, 3.10, 3.11, 3.12, and 3.13 are
  officially supported and CI-tested).
- One of the following platforms with a pre-built wheel:
    - Linux x86_64 or aarch64
    - macOS Intel or Apple Silicon
    - Windows x64

If your platform isn't covered, the source build needs a stable
Rust toolchain (`rustup install stable`); see
[CONTRIBUTING.md](https://github.com/MuteJester/GenAIRR/blob/master/CONTRIBUTING.md).

## Install

```bash
pip install GenAIRR
```

That's the whole install. Nothing else needs to be on your system —
no Rust toolchain, no compiler, no external services. A virtual
environment (`venv` or `conda`) is recommended to avoid clashing
with other repertoire-analysis packages.

## Smoke test

After install, a one-liner confirms everything is in place:

```python
import GenAIRR as ga

print(ga.__version__)
```

If this prints a version string, the Python API loaded and the
Rust kernel was importable. To run a 5-record simulation as a
deeper smoke:

```python
result = ga.Experiment.on("human_igh").recombine().run_records(n=5, seed=0)
assert len(result) == 5
print(result[0]["v_call"], result[0]["junction_aa"])
```

## Optional extras

A few features live behind opt-in extras to keep the base wheel
small:

```bash
pip install GenAIRR[all]          # numpy, scipy, graphviz, tqdm, fastmcp
pip install GenAIRR[dataconfig]   # numpy + scipy (custom DataConfig analysis)
pip install GenAIRR[viz]          # graphviz
pip install GenAIRR[mcp]          # fastmcp (for the MCP server)
```

`pandas` is needed for `result.to_dataframe()`; it ships under
`[all]`. The core simulator works without any extras — every export
format except `to_dataframe` is implemented in pure Python on the
standard library.

## Documentation build (optional, for contributors)

To build this documentation site locally:

```bash
pip install -r docs/requirements-docs.txt
make docs-serve
```

The site then renders at `http://localhost:8000` with live reload.
`make docs-build` runs a `--strict` build that fails on broken
links — the same command CI uses.

---

## Next step

→ [First simulation](quick-start.md) — generate 1,000 productive
heavy-chain sequences and inspect what comes back.
