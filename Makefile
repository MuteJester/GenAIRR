# GenAIRR — developer Makefile.
#
# GenAIRR ships as a single wheel that bundles the Rust simulation
# kernel (built into `GenAIRR/_engine`) with the Python wrappers under
# `src/GenAIRR/`. `make build` runs maturin in develop mode so the Rust
# extension lands in your active venv and the Python sources are
# editable.

PYTHON ?= python
PIP    ?= $(PYTHON) -m pip

.PHONY: help build clean rebuild test test-rust test-python

help:
	@echo "GenAIRR developer targets:"
	@echo "  make build       Build Rust extension + install GenAIRR (editable)"
	@echo "  make rebuild     clean + build"
	@echo "  make clean       Remove build artifacts"
	@echo "  make test        Run the Python test suite (pytest)"
	@echo "  make test-rust   Run the Rust unit + integration tests (cargo test)"
	@echo "  make test-python Same as 'make test'"

# --- Build ---------------------------------------------------------

build:
	$(PIP) install --upgrade "maturin>=1.5,<2.0"
	maturin develop --release

rebuild: clean build

# --- Clean ---------------------------------------------------------

clean:
	rm -rf build/ dist/ src/GenAIRR.egg-info/
	rm -rf engine_rs/target/

# --- Tests ---------------------------------------------------------

test: test-python

test-python:
	$(PYTHON) -m pytest tests/

test-rust:
	cd engine_rs && cargo test --all-features
