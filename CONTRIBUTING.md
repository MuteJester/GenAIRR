# Contributing to GenAIRR

**GenAIRR** welcomes contributions from everyone — bug reports, documentation improvements, new features, and code patches. No contribution is too small.

Large changes should be discussed first via [GitHub Issues](https://github.com/MuteJester/GenAIRR/issues).

## Development Setup

```bash
git clone https://github.com/MuteJester/GenAIRR.git
cd GenAIRR
python -m venv .venv && source .venv/bin/activate
pip install maturin pytest pandas
maturin develop --release -m engine_rs/Cargo.toml   # builds the Rust kernel
pip install -e .                                    # editable Python install
pytest tests/ -x -q                                 # run the test suite
```

Or, if you prefer the Makefile:

```bash
make build       # rebuilds Rust + reinstalls Python
make test        # pytest
make test-rust   # cargo test
```

**Requirements**: Python 3.9+ and a stable Rust toolchain (`rustup install stable`). No C compiler or CMake needed — the simulation kernel is pure Rust.

## Project Structure

```
.
├── engine_rs/                 # Rust simulation kernel (cargo workspace)
│   ├── Cargo.toml
│   └── src/
│       ├── ir.rs              # Persistent IR
│       ├── pass.rs            # Pass trait + runtime
│       ├── passes/            # Concrete passes (recombination, mutation, corruption)
│       ├── contract/          # Contract trait + filter/verify implementations
│       ├── dist.rs            # Distribution trait + concrete distributions
│       └── python/            # PyO3 wrappers exposed as `genairr_engine`
└── src/GenAIRR/
    ├── experiment.py          # Experiment DSL (public entry point)
    ├── dataconfig/            # DataConfig + species metadata
    ├── data/                  # 106+ built-in DataConfig pickles
    └── utilities/             # Helpers, IMGT regions
```

## Running Tests

```bash
pytest tests/ -x -q                  # Python tests (DSL + integration)
cd engine_rs && cargo test           # Rust tests (unit + integration)
```

## Release Process

GenAIRR ships two PyPI packages:

- **`genairr_engine`** — Rust kernel built per OS/arch via maturin (abi3-py39, so one wheel per platform covers Python 3.9+).
- **`GenAIRR`** — pure Python (single universal wheel + sdist).

GitHub Actions handles CI and publishing:

1. **Tests** run on every push/PR (Linux, macOS, Windows × Python 3.10–3.12) plus a separate `cargo test` job.
2. **Wheels** are built when a tag is pushed: maturin builds engine wheels for Linux x86_64+aarch64, macOS x86_64+arm64, Windows x64; setuptools builds a universal Python wheel.
3. **Publishing** to PyPI happens automatically on GitHub Release or version tag.

### Creating a Release

1. Update the version in `setup.py` and `engine_rs/Cargo.toml` (keep them in sync).
2. Update `CHANGELOG.md`.
3. Commit, tag, and push:
   ```bash
   git tag v1.0.0
   git push origin v1.0.0
   ```
4. Create a GitHub Release from the tag.

### Prerequisites

- `PYPI_API_TOKEN` secret must be configured in GitHub repository settings.
