"""Tests for `tools/inventory_gate.py` -- the gate over this suite's own diffs.

The gate hashes every test body and refuses a change the author has not
declared. It ran for three weeks out of a private scripts directory, where no
test covered it and its guards were trusted rather than exercised; moving it
into ``tools/`` and into CI makes it a thing the suite depends on, so each
guard is shown here to fire on the change it exists to catch.

Every case builds a two-file suite in a temp directory and captures it, which
is the gate's real entry point -- ``capture`` shells out to pytest to collect
node ids, and a stub would not exercise that.
"""
from __future__ import annotations

import importlib.util
import json
import textwrap

import pytest

from tests.support import paths

pytestmark = pytest.mark.behavior

_TOOL_PATH = paths.repo_root() / "tools" / "inventory_gate.py"


def _load_tool():
    """Load the standalone tool module (it lives outside the package)."""
    spec = importlib.util.spec_from_file_location("inventory_gate", _TOOL_PATH)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


gate = _load_tool()


_SUITE = {
    "test_alpha.py": '''
        def test_one():
            assert 1 + 1 == 2

        def test_two():
            value = [3, 1, 2]
            assert sorted(value) == [1, 2, 3]
    ''',
    "test_beta.py": '''
        LIMIT = 7

        def _helper(value):
            assert value < LIMIT
            return value

        def test_three():
            assert _helper(4) == 4
    ''',
}


def _write_suite(root, files=None):
    """Write a miniature suite under ``root/tests`` and return that directory."""
    tests = root / "tests"
    tests.mkdir(parents=True, exist_ok=True)
    for name, body in (files or _SUITE).items():
        (tests / name).write_text(textwrap.dedent(body).lstrip(), encoding="utf-8")
    return tests


def test_a_capture_is_the_same_wherever_the_suite_is_read_from(tmp_path):
    """The defect this found: keys were the path given on the command line.

    A baseline captured with ``--tests-dir`` pointing at a git worktree wrote
    absolute keys, and the member and constant checks fall back to comparing
    bare qualnames when a key's path differs -- so such a baseline silently
    degraded both to name-only for the whole suite, and a helper moving
    between modules would have passed. CI captures its baseline from a
    worktree, so this is the arrangement the gate actually runs in.
    """
    here = gate.capture(str(_write_suite(tmp_path / "a")))
    elsewhere = gate.capture(str(_write_suite(tmp_path / "b")))
    assert set(here["members"]) == set(elsewhere["members"])
    assert here["members"] == elsewhere["members"]
    assert here["constants"] == elsewhere["constants"]
    assert here["collected"] == elsewhere["collected"]
    # Keyed by the suite root and the module inside it, not by where the tree
    # was read from. Spelled in parts: the fixture modules are temp files, and
    # a literal path to one reads to the doc-reference guard as a reference to
    # a test module that does not exist.
    assert all(key.startswith("tests/") for key in here["members"])
    assert any(key.endswith("test_beta.py::_helper") for key in here["members"])


def test_an_edited_assertion_is_a_violation_until_the_author_names_it(tmp_path):
    """The whole point of the gate: a body that changed is reported."""
    tests = _write_suite(tmp_path)
    baseline = gate.capture(str(tests))
    (tests / "test_alpha.py").write_text(textwrap.dedent('''
        def test_one():
            assert 1 + 1 == 2

        def test_two():
            value = [3, 1, 2]
            assert sorted(value)
    ''').lstrip(), encoding="utf-8")
    current = gate.capture(str(tests))

    violations = gate.check(baseline, current)
    assert any("test_two" in v for v in violations), violations

    used: dict = {}
    assert gate.check(baseline, current, revised={"test_two"}, used=used) == []
    assert used["revised"] == {"test_two"}


def test_a_deleted_test_is_not_excused_by_naming_it_revised(tmp_path):
    """``--revised`` accepts an edit, never a removal: a deleted test is the
    gate's loudest finding and has to cost more than the same flag an edit
    costs."""
    tests = _write_suite(tmp_path)
    baseline = gate.capture(str(tests))
    (tests / "test_beta.py").write_text("LIMIT = 7\n", encoding="utf-8")
    current = gate.capture(str(tests))

    violations = gate.check(baseline, current, revised={"test_three", "_helper"})
    assert any("test functions no longer present" in v for v in violations), violations


def test_a_replacement_is_accepted_only_once_it_exists_and_is_collected(tmp_path):
    """``--replaced`` is the one way a removal passes, and it is checked rather
    than trusted: the successor has to be a test function the current tree
    collects, which is what makes the declaration evidence and not a note."""
    tests = _write_suite(tmp_path)
    baseline = gate.capture(str(tests))
    (tests / "test_beta.py").write_text(textwrap.dedent('''
        LIMIT = 7

        def _helper(value):
            assert value < LIMIT
            return value

        def test_three_renamed_and_rewritten():
            assert _helper(4) == 4
            assert _helper(0) == 0
    ''').lstrip(), encoding="utf-8")
    current = gate.capture(str(tests))

    missing = gate.check(baseline, current,
                         replaced={"test_three": "test_three_that_does_not_exist"})
    assert any("not a test function in the current tree" in v for v in missing), missing

    used: dict = {}
    assert gate.check(baseline, current,
                      replaced={"test_three": "test_three_renamed_and_rewritten"},
                      used=used) == []
    assert used["replaced"] == {"test_three"}


def test_a_declaration_with_no_reason_is_refused(tmp_path):
    """The reason is not machine-checkable, which is exactly why it is
    mandatory: it is what a reviewer reads, and a file of bare names is how an
    allowlist accumulates."""
    path = tmp_path / "declarations.json"
    path.write_text(json.dumps({"revised": {"test_one": "  "}}), encoding="utf-8")
    with pytest.raises(ValueError, match="no reason"):
        gate.load_declarations(path)

    path.write_text(json.dumps({"replaced": {"test_one": {"why": "because"}}}),
                    encoding="utf-8")
    with pytest.raises(ValueError, match="successor"):
        gate.load_declarations(path)

    path.write_text(json.dumps({"revisedd": {}}), encoding="utf-8")
    with pytest.raises(ValueError, match="unknown section"):
        gate.load_declarations(path)


def test_a_declaration_that_explains_nothing_is_reported_stale(tmp_path):
    """A declaration outlives the branch that needed it. Once the change is in
    the merge base it excuses nothing, and would wave through the next edit to
    the same name."""
    path = tmp_path / "declarations.json"
    path.write_text(json.dumps({
        "revised": {"test_two": "edited here", "test_long_gone": "merged weeks ago"},
    }), encoding="utf-8")
    stale = gate.stale_declarations(path, {"revised": {"test_two"}})
    assert stale == ["revised['test_long_gone']"]
    assert gate.stale_declarations(path, {"revised": {"test_two", "test_long_gone"}}) == []


def test_an_allow_class_with_no_rule_is_refused(tmp_path):
    """An empty rule used to mean "explains everything", which is the opposite
    of what a transform class is for."""
    tests = _write_suite(tmp_path)
    snapshot = gate.capture(str(tests))
    with pytest.raises(ValueError, match="unknown transform class"):
        gate.check(snapshot, snapshot, allow={"whatever-i-felt-like"})


def test_a_short_baseline_raises_rather_than_reporting_clean(tmp_path):
    """A baseline missing a section turns that section's check off. Every such
    check used ``baseline.get(section)``, so an older capture ran a fraction of
    the gate and printed ``clean``."""
    tests = _write_suite(tmp_path)
    snapshot = gate.capture(str(tests))
    short = {key: value for key, value in snapshot.items() if key != "members"}
    with pytest.raises(ValueError, match="missing"):
        gate.check(short, snapshot)


def test_ci_runs_the_gate_against_the_base_branch():
    """A gate nobody runs is a gate nobody has.

    It lived in a private scripts directory for three weeks and was run by
    hand, which means every branch's compliance depended on the author
    remembering. Pinned on the base branch specifically: `--since HEAD` or a
    committed baseline would both pass trivially.
    """
    workflow = paths.repo_root() / ".github" / "workflows" / "test.yml"
    text = workflow.read_text(encoding="utf-8")
    assert "tools/inventory_gate.py check" in text, (
        "no CI step runs the inventory gate")
    gate_lines = [line for line in text.splitlines()
                  if "inventory_gate.py" in line or "--since" in line]
    assert any("--since" in line and "base_ref" in line for line in gate_lines), (
        f"the gate does not run against the base branch: {gate_lines}")
    assert any("--declarations" in line for line in gate_lines), (
        "the gate runs without the declarations file, so every declared "
        "change would read as a violation")


def test_the_declarations_file_parses_and_every_entry_carries_a_reason():
    """The file CI depends on, checked here rather than only in CI: a syntax
    error or a reasonless entry fails the gate step with no suite output to
    read."""
    declarations = paths.repo_root() / "tools" / "inventory-declarations.json"
    allow, revised, replaced = gate.load_declarations(declarations)
    assert allow <= set(gate.ALLOWED_TRANSFORMS), (
        f"declared transform classes the gate does not know: "
        f"{sorted(allow - set(gate.ALLOWED_TRANSFORMS))}")
    # Not asserted non-empty: an empty declarations file is the steady state,
    # and is what a branch that changed no existing test should carry.
    assert isinstance(revised, set) and isinstance(replaced, dict)
