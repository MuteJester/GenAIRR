"""The assembled sequence must come from the alleles it says it came from."""
from __future__ import annotations

import pytest

from tests.biology.catalog import verifies
from tests.support.builders import experiment as experiment_builder
from tests.support.oracles import reconstruction

pytestmark = pytest.mark.behavior

LOCI = ("human_igh", "human_igk", "human_tcrb")


@verifies("recombination.assembly",
          claim="every base carrying a germline coordinate is the reference "
                "allele's base at that coordinate -- its COMPLEMENT for a "
                "segment assembled in reverse orientation, which half these "
                "parametrizations exercise -- so the bases copied are the "
                "allele's own. An allele-identity claim about the assembled "
                "sequence and nothing about the record's coordinate fields: "
                "the engine lays each byte down beside its coordinate from the "
                "same two numbers, so shifting both agrees, measured 0 "
                "mismatches of 366",
          fields=("sequence",),
          basis="reference")
@pytest.mark.parametrize("invert_d", [False, True])
@pytest.mark.parametrize("locus", LOCI)
def test_every_templated_base_matches_the_reference_at_its_claimed_coordinate(locus, invert_d):
    """Checked against the reference tables, not against the engine's own record.

    This is what catches a wrong allele assembled with internally consistent
    provenance. It does NOT catch a germline coordinate shifted by one: the
    engine lays each byte down beside its coordinate from the same two numbers,
    so a shift moves both together and every comparison here still agrees --
    measured, zero mismatches over 366 such positions. An earlier version of
    this docstring claimed otherwise. The coordinates are checked by the
    ``coordinates`` oracle, which compares them to the reference allele's own
    length instead of to its bases.

    The inverting parametrization is the point of the ``inverted`` argument
    below. Without it the argument was always an empty set: a previous commit
    added the plumbing and claimed "the assembly truth test now runs the
    inverting stack", and `d_inverted` was False on 25 of 25 records for every
    locus. The condition that let the original defect survive -- no registered
    test running an inverting stack -- was still true after the fix.
    """
    if invert_d and locus == "human_igk":
        pytest.skip("a VJ locus has no D to invert")
    exp = experiment_builder.recombined(locus=locus)
    if invert_d:
        exp = exp.invert_d(prob=1.0)
    checked = 0
    for outcome in exp.run(n=25, seed=0):
        simulation = outcome.final_simulation()
        record = outcome.airr_record(exp.refdata)
        inverted = {"D"} if record.get("d_inverted") else set()
        check = reconstruction.templated_check(simulation, exp.refdata,
                                               inverted=inverted)
        assert check.mismatches == (), (
            f"{locus}: recorded germline disagrees with the reference at "
            f"{len(check.mismatches)} position(s); first: {check.mismatches[0]}")
        # An empty mismatch tuple also describes a record where nothing was
        # compared, so the count is asserted rather than assumed.
        assert check.clean, f"{locus}: no templated position was compared"
        if invert_d:
            assert record["d_inverted"] is True, "the inverting stack did not invert"
            assert check.inverted_compared > 0, (
                "no inverted position was compared, so the inverted case "
                "passed by being skipped -- the original defect")
        assert check.checked > 50, (
            f"{locus}: only {check.checked} positions compared; a rearrangement "
            "carries hundreds of templated bases")
        # Tie the IR-level check to the exported field, so the claim on
        # `sequence` is something the body establishes rather than implies.
        assert record["sequence"].upper() == simulation.bases().decode().upper(), (
            f"{locus}: the exported sequence is not the assembled bases")
        checked += 1
    assert checked == 25, f"{locus}: only {checked} records were checked"


@verifies("recombination.trimming",
          claim="the V span present equals the allele length minus v_trim_3, "
                "checked against the reference allele rather than the engine's "
                "arithmetic, and the two non-biological trim fields are zero. "
                "D and J trims are not established here",
          fields=("v_trim_3", "v_trim_5", "j_trim_3"),
          basis="reference")
def test_recorded_trims_account_for_the_germline_that_is_missing():
    """Checked against the reference allele's length, not against the engine's
    own arithmetic.

    v_trim_5 and j_trim_3 must be zero: RAG cleaves only the RSS-adjacent
    coding end, so there is no 5' V or 3' J exonuclease loss to record. A
    non-zero value would mean the model had acquired biology it does not claim.
    """
    exp = experiment_builder.recombined()
    for outcome in exp.run(n=30, seed=4):
        record = outcome.airr_record(exp.refdata)
        simulation = outcome.final_simulation()

        assert record["v_trim_5"] == 0, "RAG does not cleave the V 5' end"
        assert record["j_trim_3"] == 0, "RAG does not cleave the J 3' end"

        v_allele = exp.refdata.v_allele(simulation.v_allele_id())
        v_sequence = v_allele.seq() if callable(v_allele.seq) else v_allele.seq
        present = record["v_germline_end"] - record["v_germline_start"]
        assert present == len(v_sequence) - record["v_trim_3"], (
            f"V spans {present} germline bases but the allele is "
            f"{len(v_sequence)} long with v_trim_3={record['v_trim_3']}")


@verifies("recombination.d_inversion",
          claim="d_inverted records the orientation the D segment was joined in, "
                "and when set the assembled D is the reverse complement of its "
                "germline span",
          fields=("d_inverted",),
          basis="reference")
def test_an_inverted_d_is_recorded_and_is_actually_reverse_complemented():
    """The flag is load-bearing: the germline coordinates stay forward-reference
    while the sequence carries the reverse complement, so anything
    reconstructing germline without reading this gets the wrong strand.
    """
    exp = experiment_builder.recombined().invert_d(prob=1.0)
    checked = 0
    for outcome in exp.run(n=30, seed=6):
        record = outcome.airr_record(exp.refdata)
        assert record["d_inverted"] is True
        start, end = record["d_sequence_start"], record["d_sequence_end"]
        if start is None or end is None or end <= start:
            continue
        simulation = outcome.final_simulation()
        d_allele = exp.refdata.d_allele(simulation.d_allele_id())
        germline = (d_allele.seq() if callable(d_allele.seq) else d_allele.seq).decode()
        span = germline[record["d_germline_start"]:record["d_germline_end"]]
        assembled = record["sequence"][start:end].upper()
        assert assembled == reconstruction.reverse_complement(span), (
            f"inverted D is {assembled!r} but the reverse complement of its "
            f"germline span is {reconstruction.reverse_complement(span)!r}")
        checked += 1
    assert checked, "no record carried an inverted D; the fixture establishes nothing"
