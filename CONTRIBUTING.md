# Contributing to GenAIRR

**GenAIRR** welcomes contributions from everyone — bug reports, documentation improvements, new features, and code patches. No contribution is too small.

Large changes should be discussed first via [GitHub Issues](https://github.com/MuteJester/GenAIRR/issues).

## Development Setup

```bash
git clone https://github.com/MuteJester/GenAIRR.git
cd GenAIRR
python -m venv .venv && source .venv/bin/activate
pip install -e ".[all]"    # builds the C backend automatically
pytest tests/ -x -q        # run the test suite
```

**Requirements**: Python 3.9+, a C compiler (gcc/clang/MSVC), and CMake 3.14+ (auto-installed by pip as a build dependency).

## Project Structure

```
src/GenAIRR/
├── experiment.py        # Experiment DSL (public entry point)
├── protocol.py          # Compilation to C engine
├── steps.py             # Step descriptors
├── _native/             # C backend (ctypes bindings)
│   ├── __init__.py      # CSimulator wrapper
│   └── csrc/            # C source code
├── dataconfig/          # DataConfig, builders, GDC I/O
├── data/                # 106+ built-in DataConfig pickles
└── utilities/           # Helpers, IMGT regions
```

## Running Tests

```bash
pytest tests/ -x -q          # all tests
pytest tests/test_c_backend.py -x -q   # C backend only
```

## Release Process

GenAIRR uses GitHub Actions for CI and automated publishing:

1. **Tests** run on every push/PR (Linux, macOS, Windows × Python 3.10–3.12).
2. **Wheels** are built via cibuildwheel for all platforms when a tag is pushed.
3. **Publishing** to PyPI happens automatically on GitHub Release or version tag.

### Creating a Release

1. Update the version in `setup.py` and `src/GenAIRR/_native/csrc/include/genairr/genairr.h`
2. Update `CHANGELOG.md`
3. Commit, tag, and push:
   ```bash
   git tag v1.0.0
   git push origin v1.0.0
   ```
4. Create a GitHub Release from the tag.

### Prerequisites

- `PYPI_API_TOKEN` secret must be configured in GitHub repository settings.
