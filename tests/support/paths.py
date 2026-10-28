"""Where the suite's fixture data lives.

Every module resolves data through here rather than counting ``..`` from its
own location, so moving a test module between directories cannot break its
data access. The anchor is this file, whose position is fixed.
"""
from __future__ import annotations

from pathlib import Path

_TESTS_ROOT = Path(__file__).resolve().parent.parent


def tests_root() -> Path:
    """The ``tests`` package directory."""
    return _TESTS_ROOT


def repo_root() -> Path:
    """The repository checkout root, the directory holding ``pyproject.toml``."""
    return _TESTS_ROOT.parent


def data(*parts: str) -> Path:
    """A path under the shared fixture-data directory."""
    return _TESTS_ROOT.joinpath("data", *parts)


def golden(*parts: str) -> Path:
    """A path under the golden-snapshot directory."""
    return _TESTS_ROOT.joinpath("golden", *parts)
