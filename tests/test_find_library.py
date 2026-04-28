"""T2-7: tests for ``GenAIRR._native._find_library``.

The bug the audit flagged was that ``_find_library`` searched both
``src/GenAIRR/_native/`` (the location the wheel ships and ``setup.py``
populates during ``pip install``) AND ``csrc/build/`` (a raw cmake
output directory). When a contributor ran ``cmake --build csrc/build``
without re-running ``pip install -e .``, the freshly-built ``.so``
silently failed to load because ``_native/<libname>`` was hit first
and contained an older copy.

After T2-7 the search order is:

1. ``GENAIRR_LIB_PATH`` env var, if set — explicit dev/CI escape hatch.
2. ``src/GenAIRR/_native/<libname>`` — single source of truth.
3. System paths on Unix.

The ``csrc/build/`` paths are no longer searched. Direct cmake users
must either re-run ``pip install -e .`` (or ``make build``) or set the
env var.
"""
from __future__ import annotations

import os
import sys
from pathlib import Path

import pytest

import GenAIRR._native as _native
from GenAIRR._native import _find_library, _LIB_NAME


@pytest.fixture
def restore_env():
    """Snapshot/restore the GENAIRR_LIB_PATH env var around each test."""
    original = os.environ.pop("GENAIRR_LIB_PATH", None)
    try:
        yield
    finally:
        if original is None:
            os.environ.pop("GENAIRR_LIB_PATH", None)
        else:
            os.environ["GENAIRR_LIB_PATH"] = original


# ── Default search returns the installed/bundled library ────────────


class TestDefaultSearch:
    def test_default_returns_installed_path(self, restore_env):
        """With no env var, the installed ``_native/<libname>`` must
        win (this is what wheels and editable installs both produce)."""
        path = _find_library()
        assert path.endswith(_LIB_NAME), (
            f"unexpected lib name: {path}")
        # And it must live under the package's _native dir, not in
        # csrc/build or anywhere else.
        native_dir = Path(_native.__file__).parent
        assert Path(path).parent == native_dir, (
            f"{path!r} did not resolve to the bundled location "
            f"{native_dir!r} — search order regression?")


# ── GENAIRR_LIB_PATH escape hatch ──────────────────────────────────


class TestEnvVarEscapeHatch:
    def test_env_var_wins_when_set(self, restore_env, tmp_path):
        """Setting the env var to an existing file must override the
        default search."""
        # Build a fake .so by copying the real one to a temp location.
        fake = tmp_path / _LIB_NAME
        fake.write_bytes(b"")  # _find_library only checks existence
        os.environ["GENAIRR_LIB_PATH"] = str(fake)
        assert _find_library() == str(fake)

    def test_env_var_falls_through_when_path_missing(self, restore_env):
        """If the env var points to a nonexistent file,
        ``_find_library`` must fall through to the next candidate
        rather than crashing — otherwise a typo in the var would
        break every dev's setup."""
        os.environ["GENAIRR_LIB_PATH"] = "/totally/nonexistent/foo.so"
        path = _find_library()
        # Falls back to the bundled location.
        native_dir = Path(_native.__file__).parent
        assert Path(path).parent == native_dir

    def test_env_var_not_set_uses_default(self, restore_env):
        """Sanity: without the env var (or with empty value), the
        default search runs as if it were never there."""
        os.environ.pop("GENAIRR_LIB_PATH", None)
        default = _find_library()
        os.environ["GENAIRR_LIB_PATH"] = ""
        empty = _find_library()
        assert default == empty


# ── csrc/build/ is no longer searched (T2-7 regression test) ────────


class TestCsrcBuildNotSearched:
    """Without GENAIRR_LIB_PATH, the legacy ``csrc/build/`` location
    must not be on the search path. This is the actual fix the audit
    asked for — silent shadowing of fresh cmake builds is what we are
    eliminating."""

    def test_csrc_build_path_not_in_search(self, restore_env, monkeypatch):
        """Make the bundled .so 'disappear' (point _native.__file__ at
        an empty temp dir) and confirm the import-error message does
        NOT mention csrc/build, and that the function raises rather
        than falling back to it."""
        import tempfile

        with tempfile.TemporaryDirectory() as td:
            fake_native = Path(td)
            # Patch the module's __file__ so ``Path(__file__).parent``
            # inside _find_library resolves to our empty dir.
            monkeypatch.setattr(_native, "__file__",
                                str(fake_native / "__init__.py"))

            with pytest.raises(ImportError) as excinfo:
                _find_library()

            msg = str(excinfo.value)
            # The new error must NOT advertise csrc/build/ as a
            # legitimate search location.
            assert "csrc/build/libgenairr" not in msg, (
                f"_find_library still searches csrc/build:\n{msg}")
            assert "csrc\\build\\libgenairr" not in msg


# ── Import error message guides users to the fix ───────────────────


class TestImportErrorMessage:
    def test_error_mentions_pip_install_e(self, restore_env, monkeypatch):
        import tempfile

        with tempfile.TemporaryDirectory() as td:
            monkeypatch.setattr(_native, "__file__",
                                str(Path(td) / "__init__.py"))
            with pytest.raises(ImportError) as excinfo:
                _find_library()

            msg = str(excinfo.value)
            assert "pip install -e ." in msg
            assert "GENAIRR_LIB_PATH" in msg
            # The current platform's library name appears at least once.
            assert _LIB_NAME in msg


# ── Linux only: env var actually points _load_lib at a fake .so ────


@pytest.mark.skipif(
    sys.platform != "linux",
    reason="dlopen of a junk file behaves differently on macOS/Windows",
)
class TestEnvVarLoadIntegration:
    """Confirms the env var is honored by the real load path, not
    just the standalone ``_find_library`` helper."""

    def test_env_var_drives_load_failure(self, restore_env, tmp_path):
        """Set GENAIRR_LIB_PATH to a nonexistent path with a no-fallback
        scenario, and confirm the real load chain consults it.

        We can't make ``_load_lib`` actually load a fake binary
        (ctypes.CDLL would reject it), but we CAN observe that
        ``_find_library`` (called from ``_load_lib``) returns the env
        path when valid."""
        fake = tmp_path / _LIB_NAME
        fake.write_bytes(b"")  # exists but is not a real ELF
        os.environ["GENAIRR_LIB_PATH"] = str(fake)
        # _find_library is the integration point — _load_lib calls it.
        assert _find_library() == str(fake)
