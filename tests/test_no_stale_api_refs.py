"""T1-10: regression test that prevents stale API references from
re-entering the codebase.

This is a fast, project-wide grep that fails if any source file,
docstring, or example code contains references to long-removed
methods or classes. New stale references should be caught at PR
review by this test rather than at runtime by an `AttributeError`.

If the user-facing API legitimately changes again, update the
``_DEAD_API_PATTERNS`` mapping with the new dead names and a
short rationale.
"""
from __future__ import annotations

import re
from pathlib import Path

import pytest


# Map: dead-API regex → human-readable reason. The regex must be
# precise enough that the test file itself (this one) does not match
# any pattern (we filter the test file out below regardless).
_DEAD_API_PATTERNS = {
    # Removed when the v2.0 Experiment DSL replaced the old fluent
    # AugmentedSimulator interface. Anything that mentions these
    # method names in source or docstrings is stale and will fail
    # at runtime.
    r"\.somatic_hypermutation\s*\(":
        "removed v2.0 API; use .mutate(rate(...), model('s5f'))",
    r"\.using_s5f\s*\(":
        "removed v2.0 API; use .mutate(model('s5f'))",
    # Old top-level types that don't exist any more.
    r"\bAugmentedSimulator\b":
        "removed v2.0 API; use Experiment.compile()",
    r"\bGenAIRR_Simulator\b":
        "renamed; use CompiledSimulator",
    # Old Protocol class — replaced by _compile_to_c.
    r"\bProtocol\s*\(":
        "removed v2.0 API; Experiment + _compile_to_c is the only path",
}

# Files / directories the scan walks. Source + tests + docs + README.
_REPO_ROOT = Path(__file__).resolve().parent.parent
_SCAN_ROOTS = [
    _REPO_ROOT / "src" / "GenAIRR",
    _REPO_ROOT / "tests",
    _REPO_ROOT / "docs",
    _REPO_ROOT / "README.md",
]
# Skip generated files, build outputs, and this test file (which
# legitimately mentions the dead names in string literals).
_SKIP_DIR_NAMES = {"__pycache__", "build", "build_san", "build_profile",
                   ".pytest_cache", "_old_docs", ".git"}
_THIS_FILE = Path(__file__).resolve()
# File extensions we care about — Python sources and Markdown / docs.
_TEXT_EXTS = {".py", ".md", ".rst", ".txt"}


def _iter_text_files():
    for root in _SCAN_ROOTS:
        if not root.exists():
            continue
        if root.is_file():
            yield root
            continue
        for p in root.rglob("*"):
            if any(part in _SKIP_DIR_NAMES for part in p.parts):
                continue
            if p.is_file() and p.suffix in _TEXT_EXTS:
                if p.resolve() == _THIS_FILE:
                    continue
                yield p


@pytest.mark.parametrize("pattern,reason", list(_DEAD_API_PATTERNS.items()))
def test_no_stale_api_reference(pattern: str, reason: str) -> None:
    rx = re.compile(pattern)
    hits: list[str] = []
    for path in _iter_text_files():
        try:
            text = path.read_text(encoding="utf-8")
        except (UnicodeDecodeError, OSError):
            continue
        for lineno, line in enumerate(text.splitlines(), start=1):
            if rx.search(line):
                rel = path.relative_to(_REPO_ROOT)
                hits.append(f"  {rel}:{lineno}: {line.strip()}")
    assert not hits, (
        f"Stale API reference matching {pattern!r} found "
        f"({reason}):\n" + "\n".join(hits)
    )
