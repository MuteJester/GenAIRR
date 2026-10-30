"""Verification that already exists on the Rust side.

engine_rs/tests/metadata_truth/ carries the same charter as this layer -- its
README opens "GenAIRR's flagship guarantee is that every simulated sequence
carries metadata that exactly reflects what a perfect aligner / biologist /
oracle would extract". Eighty-six tests across eight category files.

Not recording it would make the report declare mechanisms unverified that are
in fact well covered, which is the same lie as declaring them verified, only
inverted -- and it would create a standing incentive to reimplement in Python
what the engine already proves about itself.
"""
from __future__ import annotations

import pytest

from tests.biology import catalog
from tests.biology.catalog import all_rust_sources, coverage_rows, rust_coverage
from tests.support import paths

pytestmark = pytest.mark.pin


def test_every_recorded_rust_source_exists():
    """The anti-drift guard. If the Rust suite is reorganised, a stale record
    fails here instead of silently crediting coverage that has moved."""
    missing = [s.path for s in all_rust_sources() if not (paths.repo_root() / s.path).is_file()]
    assert missing == [], f"recorded Rust sources that no longer exist: {missing}"


def test_every_recorded_source_names_a_known_mechanism():
    unknown = sorted({m for s in all_rust_sources() for m in s.mechanism_ids
                      if m not in {x.id for x in catalog.all_mechanisms()}})
    assert unknown == [], f"Rust records naming unknown mechanisms: {unknown}"


def test_recorded_fields_are_owned_by_the_mechanism_claimed():
    wrong = []
    for source in all_rust_sources():
        for mechanism_id, fields in source.fields.items():
            owned = set(catalog.by_id(mechanism_id).truth_fields)
            wrong.extend(f"{source.path}:{mechanism_id}:{f}" for f in fields if f not in owned)
    assert wrong == [], f"Rust records claiming fields a mechanism does not own: {wrong}"


def test_rust_records_attach_to_mechanisms_and_are_visible():
    """They are recorded and shown; what they must not do is set a status --
    that is pinned by test_rust_assertions_do_not_set_a_mechanism_status."""
    rows = {r.mechanism_id: r for r in coverage_rows({}, rust_by_mechanism=rust_coverage())}
    attached = {mid for mid, row in rows.items() if row.rust_asserted}
    assert attached == set(rust_coverage()), (
        "every recorded Rust claim should appear against its mechanism")


def test_the_projection_concern_is_recorded_separately():
    """Alignment-truth fields are the projection of what every mechanism did --
    a trim, an N addition, a mutation and an indel all land in the same CIGAR.
    Letting one mechanism claim them would assert ownership it does not have.
    """
    assert catalog.PROJECTION_SOURCES
    for source in catalog.PROJECTION_SOURCES:
        assert (paths.repo_root() / source.path).is_file(), source.path
        assert source.mechanism_ids == (), "the projection is not one mechanism's"


def test_rust_assertions_do_not_set_a_mechanism_status():
    """Rust credit is hand-transcribed from module docstrings and is not
    audited, so it must not move the status.

    The transcription is demonstrably loose: indel_effects.rs references
    n_v_indels and n_d_indels zero times, corruption_stack.rs references
    end_loss zero times, and cross_segment.rs references v_germline_start and
    j_germline_end zero times -- yet all three were credited as covering those
    fields, and corruption_stack.rs's own heading reads "Invariants this module
    SHOULD test". A green row sourced from prose is worse than a red one,
    because it gets quoted.
    """
    rows = {r.mechanism_id: r for r in coverage_rows({}, rust_by_mechanism=rust_coverage())}
    credited = [r for r in rows.values() if r.rust_asserted]
    assert credited, "the Rust records credit nothing, which cannot be right"
    unearned = [r.mechanism_id for r in credited if r.status in ("verified", "partial")]
    assert unearned == [], f"Rust assertions must not set a status: {unearned}"


def test_the_python_column_reports_only_what_registered():
    """With no product test registered, the honest Python number is zero."""
    rows = coverage_rows({}, rust_by_mechanism=rust_coverage())
    assert all(r.python_tests == () for r in rows)
    assert all(r.status in ("unverified", "no-truth-exported") for r in rows)


def test_a_python_registration_is_what_moves_a_status():
    mechanism = catalog.by_id("junction.n_addition")
    rows = {r.mechanism_id: r for r in coverage_rows(
        {"junction.n_addition": ("x::y",)},
        rust_by_mechanism=rust_coverage(),
        covered_by_mechanism={"junction.n_addition": mechanism.truth_fields},
        basis_by_mechanism={"junction.n_addition":
                            {f: "reference" for f in mechanism.truth_fields}})}
    assert rows["junction.n_addition"].status == "verified"
