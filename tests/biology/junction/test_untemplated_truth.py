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
          fields=("np1", "np2", "np1_length", "np2_length"),
          basis="provenance")
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
          fields=("p_v_3_length",),
          basis="provenance")
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
          fields=("p_v_3_length", "p_d_5_length", "p_d_3_length", "p_j_5_length"),
          basis="consistent")
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


@verifies("junction.n_addition",
          claim="np1 and np2 are the stretch of the molecule between the "
                "segment bounds -- np1 == sequence[v_sequence_end:"
                "d_sequence_start], and its length -- which is what the "
                "fields promise and what holds when a polymerase insertion "
                "lands in the junction. What it does NOT establish is the "
                "TdT-added amount, which the interval also contains P "
                "nucleotides and any such insertion on top of; see the "
                "mechanism's truth_gap",
          fields=("np1", "np2", "np1_length", "np2_length"),
          basis="consistent")
@pytest.mark.parametrize("with_indels", [False, True], ids=["plain", "indels"])
def test_np_fields_are_the_interval_between_the_segments(with_indels):
    """The np fields' documented definition, checked on both sides of an indel
    pass because that is where it diverges from the mechanism's own summary.

    The run-splitting oracle cannot be used here: an inserted base carries no
    germline coordinate, so it forms an untemplated run of its own and the
    oracle reports a disagreement on 194 of 200 records whose np fields are
    correct. This compares the fields against the sequence at the coordinates
    the record itself publishes, which is why the basis is `consistent` -- two
    exported views of one molecule. It is enough to catch an np field drifting
    from its own interval, and deliberately not claimed as more.
    """
    exp = experiment_builder.recombined()
    if with_indels:
        exp = exp.polymerase_indels(count=(1, 3), insertion_prob=1.0)
    checked = nonempty = 0
    for index, outcome in enumerate(exp.run(n=100, seed=5)):
        record = outcome.airr_record(exp.refdata)
        sequence = record["sequence"]
        for np_field, length_field, start, end in (
                ("np1", "np1_length", "v_sequence_end", "d_sequence_start"),
                ("np2", "np2_length", "d_sequence_end", "j_sequence_start")):
            interval = sequence[record[start]:record[end]]
            assert record[np_field].upper() == interval.upper(), (
                f"record {index}: {np_field} is {record[np_field]!r} but the "
                f"molecule between {start} and {end} is {interval!r}")
            assert record[length_field] == len(record[np_field]), (
                f"record {index}: {length_field} is {record[length_field]} for "
                f"a {len(record[np_field])}-base {np_field}")
            checked += 1
            nonempty += bool(interval)
    assert checked == 200, f"only {checked} intervals were compared"
    assert nonempty > 100, (
        f"only {nonempty} intervals carried any base; the identity is "
        f"satisfied trivially by empty strings")
