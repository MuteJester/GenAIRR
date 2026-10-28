"""What the junction's untemplated content is, and what it is made of."""
from __future__ import annotations

import pytest

from tests.biology.catalog import verifies
from tests.support.builders import experiment as experiment_builder
from tests.support.oracles import reconstruction

pytestmark = pytest.mark.behavior


@verifies("junction.n_addition",
          claim="np1 and np2 are exactly the bases carrying no germline coordinate, "
                "in order, and their recorded lengths match",
          fields=("np1", "np2", "np1_length", "np2_length"))
def test_np_fields_are_the_bases_with_no_germline_coordinate():
    """A zero-length NP region contributes no run, so the comparison is against
    the non-empty ones. Records with np1_length == 0 do occur and are the reason
    this is not a positional comparison."""
    exp = experiment_builder.recombined()
    saw_empty = False
    for outcome in exp.run(n=50, seed=0):
        record = outcome.airr_record(exp.refdata)
        runs = reconstruction.untemplated_runs(outcome.final_simulation())
        saw_empty |= record["np1_length"] == 0 or record["np2_length"] == 0
        disagreement = reconstruction.np_fields_disagreement(record, runs)
        assert disagreement is None, disagreement
    assert saw_empty, "no record had a zero-length NP region; widen the sample"


@verifies("junction.p_nucleotide",
          claim="the V 3' P nucleotides are the reverse complement of the V's "
                "own post-trim tip, and their recorded length matches the bases "
                "present",
          fields=("p_v_3_length",))
def test_p_nucleotides_are_palindromes_of_the_segment_tip():
    """P is templated -- by the segment's own end, reversed and complemented.

    The provenance oracle records it with no germline coordinate, which is why
    it cannot verify this on its own: the distinction between P and N survives
    only in letter case.
    """
    exp = experiment_builder.with_p_nucleotides()
    checked = 0
    for outcome in exp.run(n=25, seed=1):
        simulation = outcome.final_simulation()
        record = outcome.airr_record(exp.refdata)
        bases = simulation.bases().decode()
        v_end = record["v_sequence_end"]
        p_length = record["p_v_3_length"]
        if not p_length or v_end is None:
            continue
        disagreement = reconstruction.p_tip_disagreement(bases, v_end, p_length)
        assert disagreement is None, (
            f"the V 3' P run is {disagreement[0]!r} but its tip's reverse "
            f"complement is {disagreement[1]!r}")
        checked += 1
    assert checked, "no record carried a V 3' P nucleotide; the fixture is wrong"


@verifies("junction.p_nucleotide",
          claim="P and N are distinguishable in the output, so a consumer can tell "
                "templated palindromic bases from untemplated ones",
          fields=("p_v_3_length", "p_d_5_length", "p_d_3_length", "p_j_5_length"))
def test_p_and_n_are_distinguishable_in_the_untemplated_runs():
    exp = experiment_builder.with_p_nucleotides()
    seen = set()
    for outcome in exp.run(n=25, seed=1):
        record = outcome.airr_record(exp.refdata)
        total_p = sum(record[f] for f in ("p_v_3_length", "p_d_5_length",
                                          "p_d_3_length", "p_j_5_length"))
        runs = reconstruction.untemplated_runs(outcome.final_simulation())
        lowercase = sum(sum(1 for c in r.bases if c.islower()) for r in runs)
        assert lowercase == total_p, (
            f"{lowercase} lowercase untemplated bases but {total_p} P nucleotides "
            "recorded; the P/N channel does not agree with the P lengths")
        seen.update(r.kind for r in runs)
    assert {"P", "N"} & seen, f"no run was classified P or N: {seen}"
