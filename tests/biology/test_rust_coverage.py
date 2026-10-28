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


def test_rust_coverage_lifts_mechanisms_off_unverified():
    covered = {}
    for source in all_rust_sources():
        for mechanism_id, fields in source.fields.items():
            covered.setdefault(mechanism_id, set()).update(fields)
    rows = coverage_rows({}, rust_by_mechanism=rust_coverage(),
                         covered_by_mechanism={k: tuple(v) for k, v in covered.items()})
    lifted = [r.mechanism_id for r in rows if r.status != "unverified"]
    assert lifted, "the Rust suite covers nothing, which cannot be right"


def test_the_projection_concern_is_recorded_separately():
    """Alignment-truth fields are the projection of what every mechanism did --
    a trim, an N addition, a mutation and an indel all land in the same CIGAR.
    Letting one mechanism claim them would assert ownership it does not have.
    """
    assert catalog.PROJECTION_SOURCES
    for source in catalog.PROJECTION_SOURCES:
        assert (paths.repo_root() / source.path).is_file(), source.path
        assert source.mechanism_ids == (), "the projection is not one mechanism's"
