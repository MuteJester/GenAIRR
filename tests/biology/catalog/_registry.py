"""Declaring that a test verifies a mechanism, and reporting the result."""
from __future__ import annotations

from dataclasses import dataclass, field

import pytest

from ._stages import MECHANISMS


def _mechanism_or_raise(mechanism_id: str):
    for mechanism in MECHANISMS:
        if mechanism.id == mechanism_id:
            return mechanism
    raise KeyError(
        f"unknown mechanism {mechanism_id!r}. Known ids: "
        f"{sorted(m.id for m in MECHANISMS)}")


#: How a registration establishes its claim, strongest first. Required on
#: every registration, because the report previously counted a mechanism
#: "verified" whenever a passing test named its fields -- and nothing asked
#: what the body asserted. Two mechanisms sat green for months while their
#: fields were measurably wrong: `n_pcr_errors` and `n_quality_errors`
#: reported the number of errors requested rather than applied, and the tests
#: holding them green asserted `n_mutations == 0` and `count > 0`, which any
#: value satisfies.
#:
#: - ``reference``  the expected value comes from outside the engine's output
#:                  entirely: a cartridge sequence, an allele length, a codon
#:                  table, a declared configuration parameter. The strongest
#:                  available, because a wrong value cannot also move the
#:                  thing it is checked against.
#: - ``provenance`` an independent oracle derives the expected value from the
#:                  engine's recorded per-base provenance. A different code
#:                  path from the field under test, but the same producer, so
#:                  it cannot see an error both of them share.
#: - ``cross-run``  the expected value comes from a second run of the library
#:                  -- same seed, one pass added or removed -- and a diff of
#:                  the two outputs. Independent of the ledger and trace the
#:                  field is read from.
#: - ``consistent`` two exported views of the same recorded number are
#:                  compared: a partition that must sum, two spans that must
#:                  agree. Catches a transcription slip; cannot catch a value
#:                  both views agree on.
#: - ``exercised``  the mechanism runs and something is asserted that a wrong
#:                  value would also satisfy: a field is present, a count is
#:                  non-zero, a flag varies. Establishes that the code path is
#:                  reachable and nothing about the value.
VERIFICATION_BASES = (
    "reference",
    "provenance",
    "cross-run",
    "consistent",
    "exercised",
)

#: The bases that establish a field's value. A mechanism all of whose fields
#: are covered only by the weaker two is reported ``consistent``, not
#: ``verified``.
STRONG_BASES = frozenset({"reference", "provenance", "cross-run"})


def verifies(mechanism_id: str, *, claim: str, fields: tuple[str, ...],
             basis: str):
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

    ``basis`` says how the claim is established -- see
    :data:`VERIFICATION_BASES`. It is required and has no default on purpose:
    a default would be applied to the registrations nobody thought about,
    which are exactly the ones whose basis matters.

    The function is returned unchanged -- marked, not wrapped. Nothing is
    interposed between pytest and the body, and annotating an existing test does
    not alter what it does.
    """
    if basis not in VERIFICATION_BASES:
        raise ValueError(
            f"a registration for {mechanism_id!r} declares basis {basis!r}; "
            f"it must be one of {list(VERIFICATION_BASES)}. The basis is how "
            "the test establishes its claim, not how important the claim is.")
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
    # Validated against the denominator the report actually uses. Against
    # truth_fields a registration could claim a field excluded from coverage
    # and print a numerator above its denominator -- "6/3 fields of 6".
    excluded = sorted(set(fields) - set(mechanism.verifiable_fields))
    if excluded:
        raise ValueError(
            f"{mechanism_id!r} cannot be verified on {excluded}: the catalog "
            f"records {'them' if len(excluded) > 1 else 'it'} as not being any "
            "mechanism's truth, so a test cannot establish it. If that "
            "judgement is wrong, remove the entry from UNCLAIMED_FIELDS and "
            "update the pin in test_catalog_integrity -- which is what "
            "productivity_reasons and insert_size needed.")

    def decorate(function):
        function.genairr_verifies = (mechanism_id, claim.strip(), tuple(fields), basis)
        return pytest.mark.biology(function)

    return decorate


@dataclass(frozen=True)
class CoverageRow:
    """One mechanism's verification status, and where that verification lives."""

    mechanism_id: str
    stage: str
    #: verified | consistent | gap-limited | partial | unverified |
    #: no-truth-exported
    status: str
    truth_fields: tuple[str, ...]
    python_tests: tuple[str, ...] = ()
    rust_asserted: tuple[str, ...] = ()
    covered_fields: tuple[str, ...] = ()
    truth_gap: str = ""
    #: The subset of ``truth_fields`` a test could honestly claim, which is
    #: what coverage is measured against. See ``Mechanism.verifiable_fields``.
    verifiable_fields: tuple[str, ...] = ()
    #: field -> the strongest basis any passing registration brings to it.
    #: See :data:`VERIFICATION_BASES`.
    basis_by_field: dict[str, str] = field(default_factory=dict)

    @property
    def weakly_covered_fields(self) -> tuple[str, ...]:
        """Covered fields no registration establishes the value of."""
        return tuple(f for f in self.covered_fields
                     if self.basis_by_field.get(f) not in STRONG_BASES)

    @property
    def uncovered_fields(self) -> tuple[str, ...]:
        """Verifiable fields this mechanism owns that nothing has checked."""
        return tuple(f for f in self.verifiable_fields if f not in self.covered_fields)


def coverage_rows(python_by_mechanism: dict[str, tuple[str, ...]],
                  rust_by_mechanism: dict[str, tuple[str, ...]] | None = None,
                  covered_by_mechanism: dict[str, tuple[str, ...]] | None = None,
                  basis_by_mechanism: dict[str, dict[str, str]] | None = None,
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

    ``consistent`` means every field is reached and at least one of them is
    reached only by a registration that cannot establish its value: a
    comparison between two exported views of the same number, or an assertion
    a wrong value would also satisfy. It sits between ``partial`` and
    ``verified`` because the fields ARE all covered -- the gap is in how.
    Before this distinction existed, `assay.pcr_error` and
    `assay.quality_error` were reported ``verified`` while their counters
    reported the number of errors requested rather than the number applied;
    the tests holding them green asserted `n_mutations == 0` and `count > 0`.
    """
    from . import all_mechanisms
    rust_by_mechanism = rust_by_mechanism or {}
    covered_by_mechanism = covered_by_mechanism or {}
    rows = []
    for mechanism in all_mechanisms():
        python = tuple(python_by_mechanism.get(mechanism.id, ()))
        rust = tuple(rust_by_mechanism.get(mechanism.id, ()))
        covered = tuple(f for f in mechanism.verifiable_fields
                        if f in set(covered_by_mechanism.get(mechanism.id, ())))
        bases = dict((basis_by_mechanism or {}).get(mechanism.id, {}))
        verifiable = mechanism.verifiable_fields
        if mechanism.truth_gap and not python:
            # Distinct from unverified on purpose: "no one has written the test"
            # and "there is nothing to test against" are different problems, and
            # only one of them is fixed by writing a test.
            status = "no-truth-exported"
        elif not python:
            status = "unverified"
        elif not set(covered) >= set(verifiable):
            status = "partial"
        elif any(bases.get(f) not in STRONG_BASES for f in covered):
            status = "consistent"
        elif mechanism.truth_gap:
            # Every exported field is covered, and some of what the mechanism
            # means is still not exported. `verified` in that column reads as
            # "this mechanism is covered", which is the one thing a declared
            # gap says is not true, so the ceiling for a gapped mechanism is
            # its own status rather than the green one.
            status = "gap-limited"
        else:
            status = "verified"
        rows.append(CoverageRow(mechanism_id=mechanism.id, stage=mechanism.stage,
                                status=status, truth_fields=mechanism.truth_fields,
                                python_tests=python, rust_asserted=rust,
                                covered_fields=covered, truth_gap=mechanism.truth_gap,
                                verifiable_fields=verifiable, basis_by_field=bases))
    return tuple(rows)
