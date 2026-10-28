"""Structural guarantees of the test suite itself.

These tests are about the shape of ``tests/``, not about GenAIRR. They exist
because the properties they check are the ones that quietly decay: a fixture
borrowed from a neighbouring test module, a data path built by counting
``..`` from a module's own location, a helper that drags pytest into code
meant to be usable from a script.
"""
from __future__ import annotations

import ast
from pathlib import Path

from tests.support import paths

TESTS_ROOT = Path(__file__).resolve().parents[2]

# Modules that still borrow from another test module. Empty, and it stays
# empty: shared setup belongs in tests/support. A new entry fails the test
# below, so this cannot quietly refill.
KNOWN_CROSS_IMPORTS: dict[str, set[str]] = {}


def _python_files(root: Path):
    return sorted(p for p in root.rglob("*.py") if "__pycache__" not in p.parts)


def _imported_modules(path: Path):
    """Every dotted module name imported by ``path``, including inside functions."""
    for node in ast.walk(ast.parse(path.read_text(encoding="utf-8"))):
        if isinstance(node, ast.ImportFrom) and node.module:
            yield node.module
        elif isinstance(node, ast.Import):
            for alias in node.names:
                yield alias.name


def _relative_test_imports(path: Path) -> set[str]:
    """Test modules reached by a relative import, e.g. ``from .test_x import y``.

    These are easy to miss: they carry no ``tests.`` prefix, so a check that
    looks only for the absolute form reports none while several remain.
    """
    found = set()
    for node in ast.walk(ast.parse(path.read_text(encoding="utf-8"))):
        if isinstance(node, ast.ImportFrom) and node.level and node.module \
                and node.module.split(".")[-1].startswith("test_"):
            found.add(node.module.split(".")[-1])
    return found


def _test_module_imports(path: Path) -> set[str]:
    found = _relative_test_imports(path)
    for name in _imported_modules(path):
        parts = name.split(".")
        if parts[0] == "tests" and len(parts) > 1 and parts[1].startswith("test_"):
            found.add(parts[1])
    return found


def _modules_importing_test_modules() -> dict[str, set[str]]:
    out = {}
    for path in _python_files(TESTS_ROOT):
        borrowed = _test_module_imports(path)
        if borrowed:
            out[str(path.relative_to(TESTS_ROOT))] = borrowed
    return out


def _string_literals(source: str) -> list[str]:
    """String constants in ``source``, excluding docstrings."""
    tree = ast.parse(source)
    docstrings = set()
    for node in ast.walk(tree):
        if isinstance(node, (ast.Module, ast.FunctionDef, ast.AsyncFunctionDef, ast.ClassDef)):
            body = getattr(node, "body", [])
            if body and isinstance(body[0], ast.Expr) and isinstance(body[0].value, ast.Constant) \
                    and isinstance(body[0].value.value, str):
                docstrings.add(id(body[0].value))
    return [n.value for n in ast.walk(tree)
            if isinstance(n, ast.Constant) and isinstance(n.value, str) and id(n) not in docstrings]


# ── path resolution ──────────────────────────────────────────────────────────

def test_data_resolves_regardless_of_caller_depth():
    assert paths.data("IGHV.fasta").is_file()


def test_golden_resolves_regardless_of_caller_depth():
    assert paths.golden("seed42_n10_human_igh.tsv").is_file()


def test_repo_root_is_the_directory_holding_pyproject():
    assert (paths.repo_root() / "pyproject.toml").is_file()


def test_paths_module_uses_no_hardcoded_separator():
    """Path pieces are joined with the ``/`` operator, never spelled into a literal.

    A literal "a/b" works on Windows by accident of Python's tolerance; a
    literal "a\\b" does not work anywhere else. Neither belongs here.
    """
    offenders = [s for s in _string_literals((TESTS_ROOT / "support" / "paths.py").read_text())
                 if "/" in s or "\\" in s]
    assert offenders == [], f"path separators spelled into literals: {offenders}"


# ── structural rules ─────────────────────────────────────────────────────────

def test_no_test_module_imports_another_test_module():
    """A test module is not a library. Shared setup belongs in tests/support.

    Borrowing a fixture from a neighbouring test module couples the two
    modules' collection order, makes the lender impossible to move, and hides
    where a fixture came from.
    """
    found = _modules_importing_test_modules()
    appeared = {k: v for k, v in found.items() if k not in KNOWN_CROSS_IMPORTS}
    assert appeared == {}, (
        f"new test-to-test imports: {appeared}. Put the shared setup in tests/support.")
    resolved = sorted(set(KNOWN_CROSS_IMPORTS) - set(found))
    assert resolved == [], (
        f"these no longer borrow from a test module: {resolved}. "
        "Remove them from KNOWN_CROSS_IMPORTS so the record stays honest.")


def test_builders_do_not_import_pytest():
    """tests/support/builders must stay usable from a script or the REPL."""
    builders = TESTS_ROOT / "support" / "builders"
    offenders = [str(p.relative_to(TESTS_ROOT)) for p in _python_files(builders)
                 if any(m == "pytest" or m.startswith("pytest.") for m in _imported_modules(p))]
    assert offenders == [], f"builders importing pytest: {offenders}"


def test_no_two_test_paths_differ_only_by_case():
    """Windows is in the CI matrix and its filesystem is case-insensitive."""
    lowered = [str(p.relative_to(TESTS_ROOT)).lower() for p in _python_files(TESTS_ROOT)]
    duplicates = {p for p in lowered if lowered.count(p) > 1}
    assert duplicates == set(), f"paths colliding under a case-insensitive filesystem: {duplicates}"
