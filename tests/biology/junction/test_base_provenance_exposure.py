"""Per-base provenance, read from the engine rather than inferred from case.

Python could already read half of each base's origin: `simulation.bases()`
and `simulation.germline_position(i)` are exposed, so "this base has no
germline coordinate" -- untemplated -- is genuine provenance. Which *kind*
of untemplated base it is was not exposed, and the oracle layer inferred it
from letter case instead:

    tests/support/oracles/reconstruction.py
    "Case is the channel: lowercase is a P nucleotide, uppercase an N
     addition."

Case is overloaded three ways, so that inference is a proxy that breaks:

- a P nucleotide is written lowercase;
- `sequencing_errors` lowercases every base it corrupts, so a corrupted
  TdT addition reads as a P nucleotide;
- somatic hypermutation writes uppercase, so a mutated P nucleotide reads
  as a TdT addition.

Measured before this module existed: over 200 human IGH records,
`recombine()` alone produced 382 untemplated runs, every one correctly
classified as pure N. Adding `sequencing_errors(count=8)` turned 54 of
them into `mixed` -- 54 runs apparently containing palindromic bases that
the recombination never wrote.

That mattered beyond tidiness. An oracle has to derive a fact
independently of the field it checks, and for the TdT/P/insertion split of
the junction the only independent derivation available in Python was the
one that breaks. So the tests on it could claim `consistent` at best. With
the flags exposed, they can claim `provenance`.
"""
from __future__ import annotations

import pytest

from tests.biology.catalog import verifies
from tests.support.builders import experiment as experiment_builder

pytestmark = pytest.mark.behavior

RECORDS = 60
SEED = 4

def _popcount(mask: bytes) -> int:
    return sum(1 for byte in mask if byte)


@verifies("junction.p_nucleotide",
          claim="the engine reports which positions are palindromic "
                "additions, and the count of those positions equals the sum "
                "of the four p_*_length fields -- independently of letter "
                "case, so it still holds on a read whose bases have been "
                "lowercased by sequencing error or uppercased by "
                "hypermutation",
          fields=("p_v_3_length", "p_d_5_length", "p_d_3_length", "p_j_5_length"),
          basis="provenance")
@pytest.mark.parametrize("noise", [
    pytest.param((), id="clean"),
    pytest.param((("sequencing_errors", {"count": 8}),), id="lowercasing"),
    pytest.param((("mutate", {"model": "s5f", "count": 8}),), id="uppercasing"),
])
def test_the_p_mask_counts_the_recorded_p_nucleotides(noise):
    """The mask against the record's own P lengths, under both collisions.

    The `clean` case would pass against a case-based oracle too. The other
    two are the point: each applies a pass that rewrites the channel case
    was carrying, and the mask has to be unmoved by it.
    """
    exp = experiment_builder.with_p_nucleotides_at_every_end()
    for name, kwargs in noise:
        exp = getattr(exp, name)(**kwargs)

    compared = 0
    for outcome in exp.run(n=RECORDS, seed=SEED):
        record = outcome.airr_record(exp.refdata)
        simulation = outcome.final_simulation()
        # Named one by one rather than looped over a constant: the
        # catalogue checks that a registration reads the fields it claims,
        # and a field reached through a variable is invisible to it.
        recorded = (record["p_v_3_length"] + record["p_d_5_length"]
                    + record["p_d_3_length"] + record["p_j_5_length"])
        assert _popcount(simulation.p_nucleotide_mask()) == recorded, (
            f"the P mask marks {_popcount(simulation.p_nucleotide_mask())} "
            f"positions but the record reports {recorded} P nucleotides")
        compared += recorded
    assert compared, "no record carried a P nucleotide, so nothing was checked"


def test_the_masks_are_one_byte_per_position():
    """Shape, so a caller may zip them against `bases()` without checking.

    A mask shorter than the pool would silently align the wrong positions
    with the wrong flags, which is worse than an error.
    """
    exp = experiment_builder.with_p_nucleotides_at_every_end()
    for outcome in exp.run(n=5, seed=SEED):
        simulation = outcome.final_simulation()
        expected = len(simulation.bases())
        for name in ("p_nucleotide_mask", "tdt_addition_mask",
                     "indel_inserted_mask", "junction_mask",
                     "inverted_mask"):
            mask = getattr(simulation, name)()
            assert len(mask) == expected, (
                f"{name} returned {len(mask)} bytes for a pool of {expected}")
            assert set(mask) <= {0, 1}, (
                f"{name} returned values other than 0/1: {sorted(set(mask))}")


def test_p_and_tdt_are_disjoint_and_cover_the_untemplated_bases():
    """The partition the junction counters will be built on.

    A base with no germline coordinate was written by the recombination
    rather than copied from a germline allele, so it is a palindromic
    addition, a TdT addition, or a later polymerase insertion -- and
    exactly one of those. If the three overlapped or left a gap, a
    partition exported from them would not sum to the region length, so
    this is the identity Phase 4's counters depend on.
    """
    exp = experiment_builder.with_p_nucleotides_at_every_end().polymerase_indels(count=2)
    checked = 0
    for outcome in exp.run(n=RECORDS, seed=SEED):
        simulation = outcome.final_simulation()
        p = simulation.p_nucleotide_mask()
        n = simulation.tdt_addition_mask()
        inserted = simulation.indel_inserted_mask()
        for index in range(len(simulation.bases())):
            flags = (bool(p[index]), bool(n[index]), bool(inserted[index]))
            assert sum(flags) <= 1, (
                f"position {index} carries more than one origin: "
                f"p={flags[0]} tdt={flags[1]} inserted={flags[2]}")
            untemplated = simulation.germline_position(index) is None
            if untemplated:
                assert any(flags), (
                    f"position {index} has no germline coordinate but no "
                    "origin flag either, so it belongs to no bucket")
                checked += 1
            else:
                assert not flags[2], (
                    f"position {index} is marked indel-inserted but carries "
                    "a germline coordinate")
    assert checked, "no untemplated base was seen, so the partition is untested"
