"""V6 engine — Phase A.1 skeleton smoke tests.

Verifies that the Rust kernel at `engine_rs/` builds, installs, and is
importable as the `genairr_engine` Python module. This is the smallest
end-to-end test of the Cargo → maturin → Python import chain.

If this file fails, every later V6 test will fail — fix the build
infrastructure before chasing other test failures.
"""
from __future__ import annotations

import importlib


def test_genairr_engine_module_is_importable():
    """The maturin-built Rust extension must be importable as
    `genairr_engine`. If this fails, run:
        .venv/bin/maturin develop --manifest-path engine_rs/Cargo.toml
    """
    mod = importlib.import_module("genairr_engine")
    assert mod is not None


def test_genairr_engine_version_function_exists():
    """The smoke-test `version()` function from `engine_rs/src/lib.rs`
    must be callable from Python. Confirms the PyO3 binding wiring."""
    import genairr_engine

    assert hasattr(genairr_engine, "version")
    assert callable(genairr_engine.version)


def test_genairr_engine_version_string_is_well_formed():
    """The version string must match Cargo.toml's package.version
    (currently 0.0.1 during the V6 skeleton phase). The test does not
    pin the exact value — it only asserts a non-empty SemVer-shaped
    string so the assertion survives version bumps without manual edits."""
    import genairr_engine

    v = genairr_engine.version()
    assert isinstance(v, str)
    assert len(v) > 0
    parts = v.split(".")
    assert len(parts) >= 2, f"version {v!r} is not SemVer-shaped"
    assert all(p.isdigit() for p in parts[:2]), (
        f"version {v!r} major.minor must be numeric"
    )


def test_genairr_engine_does_not_collide_with_v5_native():
    """The V6 Rust kernel must coexist with the V5 C extension at
    `GenAIRR._native` without collision. Both must be importable in
    the same Python session."""
    # V5 C extension still works.
    from GenAIRR._native import CSimulator  # noqa: F401

    # V6 Rust extension also works.
    import genairr_engine

    # They are distinct modules — no shared state.
    assert genairr_engine.__name__ != "GenAIRR._native"
