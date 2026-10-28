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
    """Test modules that ``path`` borrows from, by any import form.

    Any segment after "tests" may be the module, not just the second: once the
    tree nested, every path became tests.components.<area>.test_x, and a check
    that looked only at parts[1] saw "components" and reported nothing.
    """
    found = _relative_test_imports(path)
    for name in _imported_modules(path):
        parts = name.split(".")
        if parts[0] != "tests":
            continue
        found.update(p for p in parts[1:] if p.startswith("test_"))
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


def test_no_test_modules_sit_at_the_suite_root():
    """Every test module belongs to an area: components, biology, acceptance or
    performance. The root holds conftest.py and the package marker, nothing else.
    """
    stray = sorted(p.name for p in TESTS_ROOT.glob("test_*.py"))
    assert stray == [], f"test modules at the suite root: {stray}"


def test_every_test_directory_is_a_package():
    """pytest can collect without __init__.py, but the suite is imported as a
    package (from tests.support import ...), and a missing marker makes two
    modules with the same basename collide."""
    missing = []
    for directory in sorted({p.parent for p in TESTS_ROOT.rglob("test_*.py")}):
        if "__pycache__" in directory.parts:
            continue
        if not (directory / "__init__.py").is_file():
            missing.append(str(directory.relative_to(TESTS_ROOT)))
    assert missing == [], f"test directories without __init__.py: {missing}"


#: Tiers every product test module declares, at module level.
TIERS = ("pin", "behavior", "composition")


def _mark_names(value: ast.expr) -> set[str]:
    """The marker names in a pytestmark value.

    Matching the name matters: a substring search over the unparsed expression
    read "pin" out of ``pytest.mark.skip(reason="will pin later")``, so the
    reason string could satisfy the tier requirement.
    """
    items = value.elts if isinstance(value, (ast.List, ast.Tuple)) else [value]
    names = set()
    for item in items:
        target = item.func if isinstance(item, ast.Call) else item
        if isinstance(target, ast.Attribute):
            names.add(target.attr)
    return names


def _tier_of(module: Path) -> str | None:
    tree = ast.parse(module.read_text(encoding="utf-8"))
    for node in tree.body:
        targets = node.targets if isinstance(node, ast.Assign) else []
        if any(isinstance(x, ast.Name) and x.id == "pytestmark" for x in targets):
            found = _mark_names(node.value) & set(TIERS)
            return next(iter(sorted(found)), None)
    return None


def _product_modules():
    """Test modules that test GenAIRR, as opposed to the suite's own structure."""
    return sorted(p for p in TESTS_ROOT.rglob("test_*.py")
                  if "__pycache__" not in p.parts and "support" not in p.parts)


def test_every_test_module_declares_a_tier():
    """pin, behavior or composition, declared once per module.

    The distinction was carried by the _contract / _implementation / _release
    filenames and is load-bearing during a refactor: pins say what must not
    shift, composition says the feature still works against the whole stack.
    Declared at module level so no test body changes and `-m pin` selects
    across the entire suite.

    Modules under support/ are exempt: they test the suite's own structure,
    not GenAIRR, and the tiers do not describe them.
    """
    missing = [str(m.relative_to(TESTS_ROOT)) for m in _product_modules() if _tier_of(m) is None]
    assert missing == [], f"{len(missing)} module(s) with no tier marker: {missing[:12]}"


# ── the detectors' own self-tests ────────────────────────────────────────────
# These exist because both detectors below silently stopped working: the
# cross-import check was written against the flat layout and never updated for
# the nested one it now protects, and the tier check matched a substring of the
# whole marker expression.

def test_the_cross_import_detector_catches_every_import_form(tmp_path):
    """Absolute, relative, and `import ... as`.

    The absolute form is what every IDE auto-import and every contributor
    produces, and it is the one that stopped being detected: the check asked
    whether parts[1] began with "test_", which held while modules sat directly
    under the suite root and stopped holding once they moved into an area.
    """
    module = tmp_path / "t.py"
    module.write_text(
        "from tests.components.experiment.test_run import something\n"
        "from .test_sibling import helper\n"
        "import tests.components.duplex.test_duplex_observation as dup\n",
        encoding="utf-8")
    assert _test_module_imports(module) == {
        "test_run", "test_sibling", "test_duplex_observation"}


def test_the_cross_import_detector_ignores_support_and_product_imports(tmp_path):
    module = tmp_path / "t.py"
    module.write_text(
        "import GenAIRR as ga\n"
        "from tests.support.fixtures.source_api import from_files\n"
        "from tests.support import paths\n",
        encoding="utf-8")
    assert _test_module_imports(module) == set()


def test_a_skip_marker_is_not_a_tier(tmp_path):
    """A substring search over the unparsed marker read "pin" out of a skip
    reason, so a module could skip its entire contents and satisfy the tier
    requirement with the reason string doing the work."""
    module = tmp_path / "t.py"
    module.write_text('import pytest\npytestmark = pytest.mark.skip(reason="will pin later")\n',
                      encoding="utf-8")
    assert _tier_of(module) is None


def test_an_xfail_marker_is_not_a_tier(tmp_path):
    module = tmp_path / "t.py"
    module.write_text(
        'import pytest\npytestmark = pytest.mark.xfail(reason="composition is broken")\n',
        encoding="utf-8")
    assert _tier_of(module) is None


def test_a_real_tier_marker_is_read(tmp_path):
    for declaration in ("pytestmark = pytest.mark.pin",
                        "pytestmark = [pytest.mark.pin, pytest.mark.slow]",
                        "pytestmark = pytest.mark.behavior"):
        module = tmp_path / "t.py"
        module.write_text(f"import pytest\n{declaration}\n", encoding="utf-8")
        assert _tier_of(module) in TIERS, declaration


#: Modules whose whole contents are conditionally skipped, with why. A
#: conditional skip is legitimate, but one nobody has written down means an
#: area silently never runs while every count and tier selection still looks
#: right. Recorded so that is visible. The list only shrinks.
CONDITIONALLY_SKIPPED = {
    "acceptance/test_docs_mkdocs_migration_contract.py":
        "skipif on docs/docs_mkdocs_migration_plan.md, which .gitignore excludes "
        "-- so its 26 tests never run in CI, only in a contributor checkout",
}


def _module_pytestmark(module: Path):
    for node in ast.parse(module.read_text(encoding="utf-8")).body:
        targets = node.targets if isinstance(node, ast.Assign) else []
        if any(isinstance(x, ast.Name) and x.id == "pytestmark" for x in targets):
            return node.value
    return None


def test_no_module_unconditionally_skips_or_xfails_its_whole_contents():
    """An unconditional module-level skip or xfail removes a whole area from
    the suite while every count and every tier selection still looks right."""
    offenders = []
    for module in _product_modules():
        value = _module_pytestmark(module)
        if value is not None and _mark_names(value) & {"skip", "xfail"}:
            offenders.append(str(module.relative_to(TESTS_ROOT)))
    assert offenders == [], f"modules skipping or xfailing everything: {offenders}"


def test_every_conditionally_skipped_module_is_recorded():
    """So that "this area never runs in CI" is a line someone can read."""
    found = set()
    for module in _product_modules():
        value = _module_pytestmark(module)
        if value is not None and _mark_names(value) & {"skipif", "xfail"}:
            found.add(str(module.relative_to(TESTS_ROOT)).replace("\\", "/"))
    unrecorded = sorted(found - set(CONDITIONALLY_SKIPPED))
    assert unrecorded == [], (
        f"modules that may skip their whole contents, unrecorded: {unrecorded}. "
        "Add them to CONDITIONALLY_SKIPPED with the reason, so an area that "
        "never runs is visible.")
    resolved = sorted(set(CONDITIONALLY_SKIPPED) - found)
    assert resolved == [], f"no longer conditionally skipped, remove: {resolved}"
