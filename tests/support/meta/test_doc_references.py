"""Documentation that points at a test file must point at a file that exists.

The audit docs reference contract modules by path, both as inline code and as
markdown links, and several test modules assert that their own audit doc still
names them. That lockstep is a real guard against documentation drift, and
relocating a test module silently breaks it: the link rots while the
basename-only assertion keeps passing.

Dated verification records (JSON) are excluded. They record what was true on a
date and are not kept current.
"""
from __future__ import annotations

import re

from tests.support import paths

REPO = paths.repo_root()
# Only tracked documentation. /docs/ is gitignored, so scanning it would make
# this test's result depend on what a particular checkout happens to hold.
DOC_ROOTS = ("audit-docs",)
TEST_PATH = re.compile(r"tests/[A-Za-z0-9_/]*test_[A-Za-z0-9_]*\.py")

# References to modules that do not exist and never did: absent on master, and
# with no add-commit anywhere in the history. They name tests that were planned
# and not written, so the right target is not something to guess at. Recorded
# here so this guard catches drift the restructure introduces without being
# blocked by drift that predates it. The list only shrinks.
KNOWN_STALE = {
    "tests/test_docs_executable_examples.py",
    "tests/test_docs_links.py",
    "tests/test_paired_end.py",
    "tests/test_performance_baseline.py",
    "tests/test_trim_distribution_estimation_release.py",
}


def _prose_docs():
    for root in DOC_ROOTS:
        directory = REPO / root
        if directory.is_dir():
            yield from sorted(directory.rglob("*.md"))


def test_every_referenced_test_path_exists():
    missing = {}
    for doc in _prose_docs():
        gone = sorted({m for m in TEST_PATH.findall(doc.read_text(encoding="utf-8"))
                       if not (REPO / m).is_file() and m not in KNOWN_STALE})
        if gone:
            missing[str(doc.relative_to(REPO))] = gone
    assert missing == {}, (
        "documentation references test files that do not exist at that path:\n"
        + "\n".join(f"  {doc}: {refs}" for doc, refs in sorted(missing.items()))
        + "\n\nA relocated test module needs its doc references repointed in the same "
          "commit as the move.")


def test_the_known_stale_record_does_not_rot():
    """An entry that starts resolving must leave the record.

    Otherwise the list grows stale in the other direction and stops describing
    anything real.
    """
    resolved = sorted(ref for ref in KNOWN_STALE if (REPO / ref).is_file())
    assert resolved == [], (
        f"these now exist and should be removed from KNOWN_STALE: {resolved}")


def test_no_test_module_names_a_test_path_that_does_not_exist():
    """A test module that names another module by path must name a real one.

    Several modules carry such paths: acceptance criteria mapped to node ids,
    pins asserting a companion file exists, audit-doc lockstep assertions.
    Relocating a module leaves these pointing at nothing, and most of them are
    not themselves asserted -- so the break stays silent until something else
    reads one.
    """
    tests_root = paths.tests_root()
    stale = {}
    for module in sorted(tests_root.rglob("*.py")):
        if "__pycache__" in module.parts:
            continue
        gone = sorted({m for m in TEST_PATH.findall(module.read_text(encoding="utf-8"))
                       if not (REPO / m).is_file() and m not in KNOWN_STALE})
        if gone:
            stale[str(module.relative_to(tests_root))] = gone
    assert stale == {}, (
        "test modules naming test paths that do not exist:\n"
        + "\n".join(f"  {mod}: {refs}" for mod, refs in sorted(stale.items())))


#: Tracked, non-documentation places that name test modules by path: build
#: configuration, CI, and comments in the package that point a reader at the
#: test covering a rule. Dated verification records (JSON) are excluded, as
#: above: they describe what was true on a date.
REFERRING_FILES = (
    ".github/workflows/test.yml",
    ".github/workflows/build-wheels.yml",
    "Makefile",
    "mkdocs.yml",
    "pyproject.toml",
    "tools/build_published_reference.py",
    "tools/README.md",
)
REFERRING_TREES = ("src/GenAIRR",)


def _referring_sources():
    for name in REFERRING_FILES:
        path = REPO / name
        if path.is_file():
            yield path
    for tree in REFERRING_TREES:
        root = REPO / tree
        if root.is_dir():
            yield from sorted(p for p in root.rglob("*.py") if "__pycache__" not in p.parts)


def test_build_config_and_package_comments_name_real_test_paths():
    """A moved test module leaves these pointing at nothing.

    None of them is exercised by a test, so the staleness is silent: a Makefile
    target that names a vanished path simply collects nothing, and a comment
    that cites the wrong file misleads the next reader.
    """
    stale = {}
    for source in _referring_sources():
        gone = sorted({m for m in TEST_PATH.findall(source.read_text(encoding="utf-8"))
                       if not (REPO / m).is_file() and m not in KNOWN_STALE})
        if gone:
            stale[str(source.relative_to(REPO))] = gone
    assert stale == {}, (
        "build configuration or package comments name test paths that do not exist:\n"
        + "\n".join(f"  {src}: {refs}" for src, refs in sorted(stale.items())))
