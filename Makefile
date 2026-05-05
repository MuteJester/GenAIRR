# GenAIRR — developer Makefile.
#
# Two packages live in this repo:
#   - engine_rs/   Rust simulation kernel, built via maturin into the
#                  active venv as the `genairr_engine` extension wheel.
#   - src/GenAIRR/ Pure-Python wrappers around `genairr_engine`,
#                  installed editable.
#
# `make build` rebuilds both. Run it whenever Rust sources change;
# pure-Python edits hot-reload through the editable install.

PYTHON ?= python
PIP    ?= $(PYTHON) -m pip

.PHONY: help build clean rebuild test test-rust test-python

help:
	@echo "GenAIRR developer targets:"
	@echo "  make build       Rebuild Rust engine + reinstall Python package"
	@echo "  make rebuild     clean + build"
	@echo "  make clean       Remove build artifacts"
	@echo "  make test        Run the Python test suite (pytest)"
	@echo "  make test-rust   Run the Rust unit + integration tests (cargo test)"
	@echo "  make test-python Same as 'make test'"

# --- Build ---------------------------------------------------------

build:
	$(PIP) install --upgrade maturin
	maturin develop --release -m engine_rs/Cargo.toml
	$(PIP) install -e . --no-build-isolation

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
