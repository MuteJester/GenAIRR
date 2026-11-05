# Contributing to GenAIRR

**GenAIRR** welcomes contributions from everyone — bug reports, documentation improvements, new features, and code patches. No contribution is too small.

Large changes should be discussed first via [GitHub Issues](https://github.com/MuteJester/GenAIRR/issues).

## Development Setup

```bash
git clone https://github.com/MuteJester/GenAIRR.git
cd GenAIRR
python -m venv .venv && source .venv/bin/activate
pip install -e ".[test,mcp]"
maturin develop --release -m engine_rs/Cargo.toml   # builds the Rust kernel
pip install -e .                                    # editable Python install
make test-fast                                 # run the test suite
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
│       └── python/            # PyO3 wrappers exposed as `GenAIRR._engine`
└── src/GenAIRR/
    ├── experiment.py          # Experiment DSL (public entry point)
    ├── dataconfig/            # DataConfig + species metadata
    ├── data/                  # the 106 built-in DataConfig pickles
    └── utilities/             # Helpers, IMGT regions
```

## Running Tests

```bash
make test-fast                  # Python tests (DSL + integration)
cd engine_rs && cargo test           # Rust tests (unit + integration)
```

## Release Process

GenAIRR ships as **one** PyPI package. `maturin` builds the Rust kernel
(`engine_rs/`) and the Python wrappers (`src/GenAIRR/`) into a single wheel per
OS/arch; the kernel is exposed as the private submodule `GenAIRR._engine`. The
extension is `abi3-py39`, so one wheel per platform covers Python 3.9+.

GitHub Actions handles CI and publishing:

1. **Tests** run on every push/PR (Linux, macOS, Windows × Python 3.10–3.13)
   plus a separate `cargo test` job.
2. **Wheels** are built when a tag is pushed, for Linux x86_64+aarch64, macOS
   x86_64+arm64 and Windows x64, alongside a maturin sdist. Each built wheel is
   installed and smoke-run on Python 3.9 and 3.13 on all three platforms before
   anything is published.
3. **Publishing** to PyPI happens automatically on GitHub Release or version tag.

### Creating a Release

1. Update `CHANGELOG.md` and `MIGRATING.md`. If the release changes what a seed
   produces, bump `DETERMINISM_VERSION` in `engine_rs/src/determinism.rs`, add
   its row to the table in that file's header, and regenerate the goldens.
2. Tag and push — the release workflow derives the version from the tag name and
   rewrites `pyproject.toml` and `engine_rs/Cargo.toml` itself, so there is no
   version to edit by hand:
   ```bash
   git tag v2.4.0
   git push origin v2.4.0
   ```
3. Create a GitHub Release from the tag.

### Prerequisites

- `PYPI_API_TOKEN` secret must be configured in GitHub repository settings.


See [tests/README.md](tests/README.md) for where a new test goes, how to
pick a tier, and how to register a biological mechanism.
