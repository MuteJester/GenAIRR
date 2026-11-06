"""What CI promises to run, pinned.

The other CI assertions in this suite live beside the thing they protect --
the inventory gate's step in ``test_inventory_gate.py``, the campaign's seed in
``test_pipeline_registry.py``. The trigger has no single owner because it
decides whether *any* of them run, so it is pinned here.

Text assertions rather than a YAML parse: ``pyyaml`` is not in the ``test``
extra, and a new dependency for one test costs more than reading the lines.
This is also what the two existing CI tests do.
"""
from __future__ import annotations

import pytest

from tests.support import paths

pytestmark = pytest.mark.behavior


def _trigger_block() -> list[str]:
    """The lines of the workflow's top-level ``on:`` mapping.

    From ``on:`` to the next line that begins at column zero, which is the
    next top-level key.
    """
    text = (paths.repo_root() / ".github" / "workflows" / "test.yml").read_text(
        encoding="utf-8")
    lines = text.splitlines()
    start = next(i for i, line in enumerate(lines) if line.rstrip() == "on:")
    block = []
    for line in lines[start + 1:]:
        if line and not line[0].isspace():
            break
        block.append(line)
    return block


def test_ci_runs_on_every_pull_request_whatever_its_base():
    """A branch filter on `pull_request` makes a stacked PR invisible, not
    pending.

    `.github/workflows/test.yml` filtered `pull_request` to `[master]`, so a
    PR based on another branch reported "no checks reported" -- no pending
    run to wait for and no failure to read, just silence. On 2026-11-01 that
    silence let #45 merge into `engine/one-productivity-assessment` two hours
    after that branch had itself merged to master, so its merge commit was
    never an ancestor of master and the work did not ship. Nothing failed;
    nothing had run.

    `push` keeps its filter: that one exists to avoid a duplicate run on
    every branch push, and a push to a side branch is not a request for
    review.
    """
    block = _trigger_block()
    pull_request = [i for i, line in enumerate(block)
                    if line.strip() == "pull_request:"]
    assert len(pull_request) == 1, (
        f"expected exactly one `pull_request:` trigger, found {len(pull_request)}")
    start = pull_request[0]
    indent = len(block[start]) - len(block[start].lstrip())
    body = []
    for line in block[start + 1:]:
        if line.strip() and (len(line) - len(line.lstrip())) <= indent:
            break
        body.append(line.strip())
    offending = [line for line in body if line.startswith("branches")]
    assert offending == [], (
        f"`pull_request` is filtered by branch ({offending}), so a pull "
        f"request based on anything else gets no checks at all -- not a "
        f"failure, not a pending run, nothing to notice")


def _release_workflow() -> list[str]:
    """The lines of the release workflow."""
    return (paths.repo_root() / ".github" / "workflows" / "build-wheels.yml").read_text(
        encoding="utf-8").splitlines()


def test_the_release_workflow_installs_the_declared_test_extra():
    """A hand-written dependency list in the release job drifts silently, and
    the drift only shows up when you tag.

    `build-wheels.yml` installed `pytest pandas` by name and then the built
    wheel with the `[mcp]` extra alone. That matched the suite's needs until a
    test module imported something else. `hypothesis` arrived in
    `tests/support/pipelines/registry.py` on 2026-10-31; from then on every
    release test job died at collection with `ModuleNotFoundError`, which
    `pytest -x` reports as a single failure 0.85s in. Nothing else noticed,
    because `test.yml` installs the wheel with `[test,mcp]` and was green
    throughout -- so the break was invisible until `v2.4.0` was pushed, and
    the first evidence was a failed release.

    The extras in `pyproject.toml` are the one place that knows what the
    suite needs. The release job must install them rather than restate them.
    """
    lines = _release_workflow()

    wheel_installs = [line for line in lines
                      if "pip install" in line and "${WHEEL}" in line]
    assert wheel_installs, (
        "the release workflow no longer installs the built wheel by path; "
        "this pin cannot see what extras the test job gets")
    missing = [line.strip() for line in wheel_installs if "[test" not in line]
    assert missing == [], (
        "the release job installs the built wheel without the [test] extra:\n  "
        + "\n  ".join(missing)
        + "\n\nThe extra is what carries pytest, pandas, hypothesis, "
        "pytest-xdist and airr. Without it the suite fails at collection "
        "when a test module imports any of them -- and only when you tag.")

    restated = [line.strip() for line in lines
                if "pip install" in line
                and any(f" {pkg}" in line
                        for pkg in ("pytest", "pandas", "hypothesis"))
                and "${WHEEL}" not in line]
    assert restated == [], (
        "the release workflow names test dependencies by hand:\n  "
        + "\n  ".join(restated)
        + "\n\nInstall them through the wheel's [test] extra instead, so "
        "pyproject.toml stays the only place that knows the list.")
