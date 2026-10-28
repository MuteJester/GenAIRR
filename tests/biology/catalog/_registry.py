"""Declaring that a test verifies a mechanism, and reporting the result."""
from __future__ import annotations

from dataclasses import dataclass

import pytest

from ._stages import MECHANISMS


def _mechanism_or_raise(mechanism_id: str):
    for mechanism in MECHANISMS:
        if mechanism.id == mechanism_id:
            return mechanism
    raise KeyError(
        f"unknown mechanism {mechanism_id!r}. Known ids: "
        f"{sorted(m.id for m in MECHANISMS)}")


def verifies(mechanism_id: str, *, claim: str, fields: tuple[str, ...]):
    """Declare that this test verifies ``mechanism_id``, and what it establishes.

    ``fields`` names the truth fields this test actually checks. It is required
    and validated against the mechanism, because the report stands or falls on
    it: a test that checks np1 while the mechanism is answerable for six fields
    is a partial result, and defaulting to "all of them" would let one test
    report the mechanism verified.

    The id and the fields are checked now, at decoration time, so a typo fails
    during collection. Recording an unknown id silently would let a mechanism
    look verified by a test that does not exist, or count nothing at all, and
    either way nobody would be told.

    The function is returned unchanged -- marked, not wrapped. Nothing is
    interposed between pytest and the body, and annotating an existing test does
    not alter what it does.
    """
    mechanism = _mechanism_or_raise(mechanism_id)
    if not claim.strip():
        raise ValueError(
            f"a registration for {mechanism_id!r} needs a claim: what does this "
            "test establish about the recorded truth?")
    if not fields:
        raise ValueError(
            f"a registration for {mechanism_id!r} must name the truth fields it "
            f"checks; the mechanism owns {list(mechanism.truth_fields)}")
    unknown = sorted(set(fields) - set(mechanism.truth_fields))
    if unknown:
        raise ValueError(
            f"{mechanism_id!r} does not own {unknown}; it owns "
            f"{list(mechanism.truth_fields)}")

    def decorate(function):
        function.genairr_verifies = (mechanism_id, claim.strip(), tuple(fields))
        return pytest.mark.biology(function)

    return decorate


@dataclass(frozen=True)
class CoverageRow:
    """One mechanism's verification status, and where that verification lives."""

    mechanism_id: str
    stage: str
    status: str                       # verified | partial | unverified
    truth_fields: tuple[str, ...]
    python_tests: tuple[str, ...] = ()
    rust_asserted: tuple[str, ...] = ()
    covered_fields: tuple[str, ...] = ()

    @property
    def uncovered_fields(self) -> tuple[str, ...]:
        """Truth fields this mechanism owns that nothing has checked."""
        return tuple(f for f in self.truth_fields if f not in self.covered_fields)


def coverage_rows(python_by_mechanism: dict[str, tuple[str, ...]],
                  rust_by_mechanism: dict[str, tuple[str, ...]] | None = None,
                  covered_by_mechanism: dict[str, tuple[str, ...]] | None = None,
                  ) -> tuple[CoverageRow, ...]:
    """A row per declared mechanism, in developmental order.

    Status is set by registered Python tests only. Rust credit rides alongside
    in ``rust_asserted`` and deliberately cannot move it: those records are
    transcribed by hand from the category files' module docstrings, several of
    which are headed "Invariants this module SHOULD test", and spot-checking
    found fields credited to files that never mention them. Until a Rust test
    can declare its own mechanism, that is an assertion by a reader, not a
    measurement, and must not read as green.

    ``partial`` means registered tests exist but do not reach every truth field
    the mechanism owns -- one test against a mechanism answerable for twelve
    fields is a much weaker claim than "verified".
    """
    from . import all_mechanisms
    rust_by_mechanism = rust_by_mechanism or {}
    covered_by_mechanism = covered_by_mechanism or {}
    rows = []
    for mechanism in all_mechanisms():
        python = tuple(python_by_mechanism.get(mechanism.id, ()))
        rust = tuple(rust_by_mechanism.get(mechanism.id, ()))
        covered = tuple(f for f in mechanism.truth_fields
                        if f in set(covered_by_mechanism.get(mechanism.id, ())))
        if not python:
            status = "unverified"
        elif set(covered) >= set(mechanism.truth_fields):
            status = "verified"
        else:
            status = "partial"
        rows.append(CoverageRow(mechanism_id=mechanism.id, stage=mechanism.stage,
                                status=status, truth_fields=mechanism.truth_fields,
                                python_tests=python, rust_asserted=rust,
                                covered_fields=covered))
    return tuple(rows)
