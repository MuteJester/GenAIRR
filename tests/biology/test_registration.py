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


def test_status_is_one_of_three_values():
    assert {r.status for r in coverage_rows({})} <= {"verified", "partial", "unverified"}


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
