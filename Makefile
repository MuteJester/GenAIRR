# GenAIRR — developer Makefile.
#
# GenAIRR ships as a single wheel that bundles the Rust simulation
# kernel (built into `GenAIRR/_engine`) with the Python wrappers under
# `src/GenAIRR/`. `make build` runs maturin in develop mode so the Rust
# extension lands in your active venv and the Python sources are
# editable.

PYTHON ?= python
PIP    ?= $(PYTHON) -m pip

.PHONY: help build clean rebuild test test-rust test-python \
        validate-fast validate-full validate-release \
        docs-install docs-build docs-serve docs-clean

help:
	@echo "GenAIRR developer targets:"
	@echo "  make build              Build Rust extension + install GenAIRR (editable)"
	@echo "  make rebuild            clean + build"
	@echo "  make clean              Remove build artifacts"
	@echo "  make test               Run the Python test suite (pytest, default fast)"
	@echo "  make test-rust          Run the Rust unit + integration tests (cargo test)"
	@echo "  make test-python        Same as 'make test'"
	@echo ""
	@echo "Validation tiers (see docs/validation_matrix.md):"
	@echo "  make validate-fast      Tier 1: correctness only (pytest -m \"not performance\" + cargo test --lib)"
	@echo "  make validate-full      Tier 2: correctness + performance budgets"
	@echo "  make validate-release   Tier 3: full + wheel build + golden trace compat"
	@echo ""
	@echo "Documentation (MkDocs Material — see docs/docs_mkdocs_migration_plan.md):"
	@echo "  make docs-install       Install docs-build deps from docs/requirements-docs.txt"
	@echo "  make docs-build         mkdocs build --strict (fails on broken links / nav drift)"
	@echo "  make docs-serve         mkdocs serve (live-reload local preview on :8000)"
	@echo "  make docs-clean         Remove the built site/ directory"

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

# --- Validation tiers ---------------------------------------------
#
# See docs/validation_matrix.md for the full breakdown. The three
# tiers trade runtime for coverage:
#
#   fast    — pre-commit / per-push    correctness + replay + dist invariants
#   full    — PR ready-for-review      fast + performance budgets (regression guard)
#   release — pre-tag                  full + wheel build + golden trace compat

validate-fast:
	$(PYTHON) -m pytest tests/ -m "not performance" -q
	cd engine_rs && cargo test --lib

validate-full:
	$(PYTHON) -m pytest tests/ -q
	cd engine_rs && cargo test --all-features

validate-release: validate-full
	$(PYTHON) -m pytest tests/test_trace_file_compat.py -q
	maturin build --release

# --- Docs (MkDocs Material — Phase 1 onwards) ---------------------
#
# The docs source lives at site_docs/. The contributor-only audit/
# design markdowns at docs/*.md stay rendered on GitHub and are
# linked from the future Architecture section of the site.
#
# `docs-build` uses `--strict` so broken nav links / missing files
# / dead anchors fail the build. Run this in CI once Phase 6
# cutover lands.

DOCS_REQ ?= docs/requirements-docs.txt

docs-install:
	$(PIP) install -r $(DOCS_REQ)

docs-build:
	$(PYTHON) -m mkdocs build --strict

docs-serve:
	$(PYTHON) -m mkdocs serve

docs-clean:
	rm -rf site/
