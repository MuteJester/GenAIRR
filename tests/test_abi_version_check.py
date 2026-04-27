"""T1-12: ABI version-check sanity tests.

The ``_load_lib()`` path in ``GenAIRR._native`` asserts that the
loaded C library reports the same version string as the Python
package's ``__version__``. A mismatch typically means a stale
``libgenairr.so`` was found on the search path (e.g. an old editable
install left behind, or a system-wide build of a different version) —
loading a mismatched binary yields silent ABI corruption.

This test directly exercises the check via ``_check_abi_version``
since faking a real library mismatch would require hot-swapping the
shared object on disk.
"""
from __future__ import annotations

import pytest

import GenAIRR
from GenAIRR._native import _check_abi_version, get_version


class _FakeLib:
    """Minimal stand-in for ctypes.CDLL — we only need genairr_version."""

    def __init__(self, version: str):
        self._version = version.encode("ascii")

        class _Fn:
            def __init__(self, v):
                self._v = v

            def __call__(self):
                return self._v

        self.genairr_version = _Fn(self._version)


def test_happy_path_no_raise_on_match():
    """The real loaded C library matches Python __version__ — must
    not raise. This is the production happy-path."""
    real_c = get_version()
    assert real_c == GenAIRR.__version__, (
        f"Setup mismatch: C lib reports {real_c!r}, Python is "
        f"{GenAIRR.__version__!r}. Bump CMakeLists.txt and "
        f"GENAIRR_C_VERSION_STRING in lockstep with setup.py."
    )


def test_check_raises_on_mismatch():
    """Direct call to _check_abi_version with a fake C library that
    reports a different version must raise RuntimeError with a clear
    message redirecting to a reinstall."""
    fake = _FakeLib("0.1.0-stale")
    with pytest.raises(RuntimeError) as excinfo:
        _check_abi_version(fake)   # type: ignore[arg-type]
    msg = str(excinfo.value)
    assert "0.1.0-stale" in msg
    assert GenAIRR.__version__ in msg
    assert "reinstall" in msg.lower()


def test_check_no_raise_when_versions_match():
    """Direct call with a fake C library reporting the SAME version
    as Python — must not raise."""
    fake = _FakeLib(GenAIRR.__version__)
    _check_abi_version(fake)   # type: ignore[arg-type]
