# GenAIRR — developer Makefile (T2-7).
#
# The canonical dev workflow is now `make build`, which runs setup.py's
# CMakeBuild and copies the freshly-built libgenairr.so into
# src/GenAIRR/_native/ — the single location _find_library() searches.
# Running `cmake --build src/GenAIRR/_native/csrc/build` directly no
# longer suffices, because that path is no longer on the search list
# (it used to silently shadow installed libraries).

PYTHON      ?= python
PIP         ?= $(PYTHON) -m pip
CSRC_DIR    := src/GenAIRR/_native/csrc
NATIVE_DIR  := src/GenAIRR/_native

.PHONY: help build clean test test-c rebuild

help:
	@echo "GenAIRR developer targets:"
	@echo "  make build      Rebuild C backend + reinstall (canonical dev rebuild)"
	@echo "  make clean      Remove build trees and bundled .so files"
	@echo "  make rebuild    clean + build"
	@echo "  make test       Run the Python test suite (pytest)"
	@echo "  make test-c     Configure + build the C tests, then run ctest"

# --- Build ---------------------------------------------------------

build:
	$(PIP) install -e . --no-build-isolation

rebuild: clean build

# --- Clean ---------------------------------------------------------
#
# Removes *all* in-tree C build artifacts. The bundled .so glob covers
# every SOVERSION suffix (libgenairr.so, libgenairr.so.0,
# libgenairr.so.1.0.0, libgenairr.dylib, genairr.dll) so a stale build
# from a previous version cannot lurk after a clean.

clean:
	rm -rf $(CSRC_DIR)/build $(CSRC_DIR)/build_*
	rm -f $(NATIVE_DIR)/libgenairr.so* $(NATIVE_DIR)/libgenairr.dylib* \
	      $(NATIVE_DIR)/genairr.dll
	rm -rf build/ dist/ src/GenAIRR.egg-info/

# --- Tests ---------------------------------------------------------

test:
	$(PYTHON) -m pytest tests/

test-c:
	cmake -S $(CSRC_DIR) -B $(CSRC_DIR)/build -DCMAKE_BUILD_TYPE=Release
	cmake --build $(CSRC_DIR)/build
	cd $(CSRC_DIR)/build && ctest --output-on-failure
