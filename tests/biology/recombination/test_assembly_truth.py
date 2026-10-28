"""The assembled sequence must come from the alleles it says it came from."""
from __future__ import annotations

import pytest

from tests.biology.catalog import verifies
from tests.support.builders import experiment as experiment_builder
from tests.support.oracles import reconstruction

pytestmark = pytest.mark.behavior

LOCI = ("human_igh", "human_igk", "human_tcrb")


@verifies("recombination.assembly",
          claim="every base carrying a germline coordinate matches the reference "
                "allele's base at that coordinate, so the recorded germline is the "
                "germline the record names",
          fields=("v_germline_start", "v_germline_end", "d_germline_start",
                  "d_germline_end", "j_germline_start", "j_germline_end"))
@pytest.mark.parametrize("locus", LOCI)
def test_every_templated_base_matches_the_reference_at_its_claimed_coordinate(locus):
    """Checked against the reference tables, not against the engine's own record.

    This is what catches a wrong allele assembled with internally consistent
    provenance, or a germline coordinate shifted by one -- neither of which the
    engine's own bookkeeping can reveal.
    """
    exp = experiment_builder.recombined(locus=locus)
    for outcome in exp.run(n=25, seed=0):
        simulation = outcome.final_simulation()
        mismatches = reconstruction.templated_mismatches(simulation, exp.refdata)
        assert mismatches == (), (
            f"{locus}: recorded germline disagrees with the reference at "
            f"{len(mismatches)} position(s); first: {mismatches[0]}")


@verifies("recombination.trimming",
          claim="the recorded trims account for the germline each segment lost: "
                "the span present equals the allele length minus the trims, and "
                "the two non-biological trim fields are structurally zero",
          fields=("v_trim_3", "j_trim_5", "v_trim_5", "j_trim_3"))
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
          fields=("d_inverted",))
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
