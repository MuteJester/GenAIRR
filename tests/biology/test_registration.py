"""Registering a verification, and reporting what is covered.

A test declares which mechanism it verifies and what claim it establishes. The
report then aggregates those declarations, so coverage reflects what was proved
rather than where the file sits.
"""
from __future__ import annotations

import pytest

from tests.biology import catalog
from tests.biology.catalog import coverage_rows, verifies

pytestmark = pytest.mark.pin


def test_verifies_records_the_mechanism_and_the_claim():
    @verifies("junction.n_addition", claim="np1 is exactly the untemplated bases",
              fields=("np1",))
    def sample():
        return 42

    assert sample.genairr_verifies == ("junction.n_addition",
                                       "np1 is exactly the untemplated bases", ("np1",))


def test_verifies_returns_the_same_function_object():
    """No wrapper: decorating must not change what the test does, and must not
    interpose anything between pytest and the body."""
    def sample():
        return 42
    decorated = verifies("junction.n_addition", claim="x", fields=("np1",))(sample)
    assert decorated is sample
    assert decorated() == 42


def test_verifies_marks_the_test_as_biology():
    @verifies("junction.n_addition", claim="x", fields=("np1",))
    def sample():
        pass
    assert any(m.name == "biology" for m in sample.pytestmark)


def test_an_unknown_mechanism_id_raises_at_decoration_time():
    """A typo must surface. Silently recording an unknown id would let a
    mechanism look verified by a test that does not exist, or count nothing at
    all -- and either way no one would be told."""
    with pytest.raises(KeyError, match="junction.n_additon"):
        @verifies("junction.n_additon", claim="typo", fields=("np1",))
        def sample():
            pass


def test_an_empty_claim_is_refused():
    """The claim is the point of the registration: it says what was established,
    which is what a reader of the report needs."""
    with pytest.raises(ValueError):
        @verifies("junction.n_addition", claim="   ", fields=("np1",))
        def sample():
            pass


def test_coverage_reports_a_row_per_mechanism():
    rows = coverage_rows({})
    assert [r.mechanism_id for r in rows] == [m.id for m in catalog.all_mechanisms()]


def test_a_mechanism_with_no_test_is_unverified():
    rows = {r.mechanism_id: r for r in coverage_rows({})}
    assert rows["junction.n_addition"].status == "unverified"
    assert rows["junction.n_addition"].python_tests == ()


def test_a_registered_test_counts_wherever_it_lives():
    """Registration, not location: a test under components/ counts the same as
    one under biology/."""
    rows = {r.mechanism_id: r
            for r in coverage_rows({"junction.n_addition": ("components/engine/test_x.py::test_y",)})}
    assert rows["junction.n_addition"].python_tests == ("components/engine/test_x.py::test_y",)


def test_one_test_against_many_fields_is_partial_not_verified():
    """The honest distinction. junction.n_addition owns six fields; a test that
    checks one of them has not verified the mechanism."""
    mechanism = catalog.by_id("junction.n_addition")
    assert len(mechanism.truth_fields) > 1
    rows = {r.mechanism_id: r for r in coverage_rows(
        {"junction.n_addition": ("test_x.py::test_y",)},
        covered_by_mechanism={"junction.n_addition": ("np1",)})}
    row = rows["junction.n_addition"]
    assert row.status == "partial"
    assert "np2" in row.uncovered_fields


def test_covering_every_field_is_verified():
    mechanism = catalog.by_id("junction.n_addition")
    rows = {r.mechanism_id: r for r in coverage_rows(
        {"junction.n_addition": ("test_x.py::test_y",)},
        covered_by_mechanism={"junction.n_addition": mechanism.truth_fields})}
    assert rows["junction.n_addition"].status == "verified"
    assert rows["junction.n_addition"].uncovered_fields == ()


def test_a_field_the_mechanism_does_not_own_is_refused():
    with pytest.raises(ValueError, match="does not own"):
        @verifies("junction.n_addition", claim="x", fields=("v_call",))
        def sample():
            pass


def test_fields_must_be_given():
    with pytest.raises(ValueError, match="must name the truth fields"):
        @verifies("junction.n_addition", claim="x", fields=())
        def sample():
            pass


def test_status_is_one_of_the_declared_values():
    """no-truth-exported is the fourth: a mechanism whose truth is not in the
    record cannot be established by any test, and saying "unverified" would
    invite someone to go and write one."""
    assert {r.status for r in coverage_rows({})} <= {
        "verified", "partial", "unverified", "no-truth-exported"}


def test_a_typo_is_a_collection_error_in_a_real_module(tmp_path):
    """Not merely a raise in isolation: pytest must refuse to collect.

    This is the whole safety of the scheme. If a mistyped id were tolerated,
    the mechanism it was meant to cover would read as unverified while the
    author believed it covered -- and the report would be confidently wrong.
    """
    import subprocess
    import sys
    from tests.support import paths

    module = tmp_path / "test_typo.py"
    module.write_text(
        "from tests.biology.catalog import verifies\n\n\n"
        '@verifies("junction.n_additon", claim="typo", fields=("np1",))\n'
        "def test_x():\n"
        "    pass\n",
        encoding="utf-8")
    result = subprocess.run(
        [sys.executable, "-m", "pytest", str(module), "--collect-only", "-q", "-p", "no:cacheprovider"],
        capture_output=True, text=True, cwd=paths.repo_root())
    assert result.returncode != 0, result.stdout
    assert "junction.n_additon" in result.stdout + result.stderr


def _coverage_line(tmp_path, *extra_args):
    """Run a module with two registrations and return its report line."""
    import subprocess
    import sys
    from tests.support import paths

    module = tmp_path / "test_registered.py"
    module.write_text(
        "from tests.biology.catalog import verifies\n\n\n"
        '@verifies("junction.p_nucleotide", claim="the 5\' V P length is recorded",\n'
        '          fields=("p_v_3_length",))\n'
        "def test_one():\n"
        "    assert True\n\n\n"
        '@verifies("junction.p_nucleotide", claim="the 3\' J P length is recorded",\n'
        '          fields=("p_j_5_length",))\n'
        "def test_two():\n"
        "    assert True\n",
        encoding="utf-8")
    result = subprocess.run(
        [sys.executable, "-m", "pytest", str(module), "-q", "--biology-coverage",
         "-p", "tests.biology.catalog._plugin", "-p", "no:cacheprovider", *extra_args],
        capture_output=True, text=True, cwd=paths.repo_root())
    lines = [ln for ln in result.stdout.splitlines() if "junction.p_nucleotide" in ln]
    assert lines, f"no report line produced:\n{result.stdout}\n{result.stderr}"
    return lines[0].split()


def test_registrations_survive_distributed_collection(tmp_path):
    """The report must say the same thing under -n as serial.

    Registrations are gathered during collection, which happens in the xdist
    workers, while the report is written by the controller. CI runs -n auto, so
    without shipping them across that boundary the report would say
    "unverified" the moment anyone annotated a test -- the inverted lie, by
    default, under the project's own command.
    """
    serial = _coverage_line(tmp_path)
    distributed = _coverage_line(tmp_path, "-n", "2")
    assert serial == distributed, f"serial {serial} != distributed {distributed}"


def test_two_registrations_on_one_mechanism_are_both_counted(tmp_path):
    fields = _coverage_line(tmp_path)
    assert "py:2" in fields, fields
    assert "2/4" in fields, fields


def test_a_failing_or_skipped_registration_does_not_count_as_coverage(tmp_path):
    """An xfail or skip marks a mechanism as NOT verified.

    A registration on a test that cannot pass would otherwise credit the
    mechanism with a field the test does not actually establish -- which is the
    same lie as a prose-sourced green row, arriving by a different route. Known
    defects are pinned as strict xfails, so this matters immediately.
    """
    import subprocess
    import sys
    from tests.support import paths

    module = tmp_path / "test_xfailed.py"
    module.write_text(
        "import pytest\n"
        "from tests.biology.catalog import verifies\n\n\n"
        '@pytest.mark.xfail(strict=True, reason="known defect")\n'
        '@verifies("junction.p_nucleotide", claim="x", fields=("p_v_3_length",))\n'
        "def test_known_defect():\n"
        "    assert False\n",
        encoding="utf-8")
    result = subprocess.run(
        [sys.executable, "-m", "pytest", str(module), "-q", "--biology-coverage",
         "-p", "tests.biology.catalog._plugin", "-p", "no:cacheprovider"],
        capture_output=True, text=True, cwd=paths.repo_root())
    line = next((ln for ln in result.stdout.splitlines()
                 if "junction.p_nucleotide" in ln), "")
    assert line, result.stdout
    assert "unverified" in line and "py:0" in line, (
        f"an xfailed registration must not count: {line!r}")


def test_a_registration_on_a_failing_test_does_not_count(tmp_path):
    """The plain-failure route to a false green.

    The xfail route was closed by inspecting markers, but the report was built
    at collection time, so a registered test whose body simply fails still
    credited its mechanism. A mechanism cannot be verified by a test that the
    same run watched fail.
    """
    import subprocess
    import sys
    from tests.support import paths

    module = tmp_path / "test_failing.py"
    module.write_text(
        "from tests.biology.catalog import verifies\n\n\n"
        '@verifies("assay.contamination", claim="x", fields=("is_contaminant",))\n'
        "def test_it_fails():\n"
        "    assert False\n",
        encoding="utf-8")
    result = subprocess.run(
        [sys.executable, "-m", "pytest", str(module), "-q", "--biology-coverage",
         "-p", "tests.biology.catalog._plugin", "-p", "no:cacheprovider"],
        capture_output=True, text=True, cwd=paths.repo_root())
    # the report line, identified by its status column -- the traceback also
    # names the mechanism, in the decorator it echoes back
    line = next((ln for ln in result.stdout.splitlines()
                 if "assay.contamination" in ln
                 and any(s in ln for s in ("unverified", "partial", "verified",
                                           "no-truth-exported"))), "")
    assert line, result.stdout
    # "unverified" contains "verified", so check the status token, not a substring
    status = line.split()[1]
    assert status == "unverified" and "py:0" in line, (
        f"a failing registration credited its mechanism: {line!r}")


def test_the_coverage_plumbing_is_not_silently_dead(tmp_path):
    """The report must reflect registrations that passed.

    When credit started requiring an observed pass, the hook that records
    outcomes was not wired into conftest, so every row read `unverified` while
    the whole suite was green -- a report that had quietly stopped measuring
    anything. Nothing noticed, because nothing asserted the report was
    non-degenerate.
    """
    import subprocess
    import sys
    from tests.support import paths

    module = tmp_path / "test_counts.py"
    module.write_text(
        "from tests.biology.catalog import verifies\n\n\n"
        '@verifies("assay.contamination", claim="it is flagged",\n'
        '          fields=("is_contaminant",))\n'
        "def test_passes():\n"
        "    assert True\n",
        encoding="utf-8")
    result = subprocess.run(
        [sys.executable, "-m", "pytest", str(module), "-q", "--biology-coverage",
         "-p", "tests.biology.catalog._plugin", "-p", "no:cacheprovider"],
        capture_output=True, text=True, cwd=paths.repo_root())
    line = next((ln for ln in result.stdout.splitlines()
                 if "assay.contamination" in ln and "fields" in ln
                 and "@verifies" not in ln), "")
    assert line, result.stdout
    assert "py:1" in line, (
        f"a passing registration was not counted, so the report is measuring "
        f"nothing: {line!r}")
    assert line.split()[1] == "verified", line
