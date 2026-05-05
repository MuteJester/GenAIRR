"""Smoke tests for the Rust simulation kernel.

Verifies that the Rust kernel at ``engine_rs/`` builds, installs, and
is importable as the ``genairr_engine`` Python module — the smallest
end-to-end test of the Cargo → maturin → Python import chain.

If this file fails, every later test will fail — fix the build
infrastructure (typically ``maturin develop -m engine_rs/Cargo.toml``)
before chasing other test failures.
"""
from __future__ import annotations

import importlib


def test_genairr_engine_module_is_importable():
    """The maturin-built Rust extension must be importable as
    ``genairr_engine``. If this fails, run:
        .venv/bin/maturin develop --manifest-path engine_rs/Cargo.toml
    """
    mod = importlib.import_module("genairr_engine")
    assert mod is not None


def test_genairr_engine_version_function_exists():
    """The smoke-test ``version()`` function from
    ``engine_rs/src/lib.rs`` must be callable from Python. Confirms
    the PyO3 binding wiring."""
    import genairr_engine

    assert hasattr(genairr_engine, "version")
    assert callable(genairr_engine.version)


def test_genairr_engine_version_string_is_well_formed():
    """The version string must match Cargo.toml's package.version. The
    test does not pin the exact value — it only asserts a non-empty
    SemVer-shaped string so the assertion survives version bumps
    without manual edits."""
    import genairr_engine

    v = genairr_engine.version()
    assert isinstance(v, str)
    assert len(v) > 0
    parts = v.split(".")
    assert len(parts) >= 2, f"version {v!r} is not SemVer-shaped"
    assert all(p.isdigit() for p in parts[:2]), (
        f"version {v!r} major.minor must be numeric"
    )
