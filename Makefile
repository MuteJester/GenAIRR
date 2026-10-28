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
	@echo "  make docs-install       Install docs-build deps from requirements-docs.txt"
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
	$(PYTHON) -m pytest tests/components/io/test_trace_file_compat.py -q
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

DOCS_REQ ?= requirements-docs.txt

docs-install:
	$(PIP) install -r $(DOCS_REQ)

docs-build:
	$(PYTHON) -m mkdocs build --strict

docs-serve:
	$(PYTHON) -m mkdocs serve

docs-clean:
	rm -rf site/

# ── test selections ─────────────────────────────────────────────────
# Names say what they do. `test-fast` is the inner loop: it drops the
# eight slowest modules and the wall-time budgets and distributes across
# cores. The marker alone only buys about 2x -- most of the speed is
# xdist -- so this is "parallel minus the heavy tail", not a serial
# fast subset.
.PHONY: test-fast test-parallel test-budgets test-lint

test-fast:
	pytest tests/ -q -m "not slow and not performance" -n auto

test-parallel:
	pytest tests/ -q -m "not performance" -n auto

# Serially and on their own: under xdist the other workers saturate the
# cores and the measured time reflects contention, not the engine.
test-budgets:
	pytest tests/ -q -m performance

# .private/scripts is named explicitly because ruff respects .gitignore, so a
# repo-wide run skips it -- which left the gate that enforces the migration
# invariant unlinted. The path is skipped silently when the directory is absent,
# as it is in any clone.
test-lint:
	ruff check tests/
	@test -d .private/scripts && ruff check .private/scripts/ || true

# The cross-platform reproducibility references. Regenerate only when the
# output format changes on purpose, then commit them: they are what every CI
# platform is checked against.
.PHONY: goldens
goldens:
	python tools/regenerate_golden_files.py
