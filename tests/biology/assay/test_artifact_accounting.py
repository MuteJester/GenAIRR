"""Each assay artifact is counted as itself and not as biology.

This is the distinction that matters most to a benchmark: a PCR
misincorporation and a somatic hypermutation both look like a difference from
germline, and a tool scored against GenAIRR truth is only fair if the record
says which was which.
"""
from __future__ import annotations

import pytest

from tests.biology.catalog import verifies
from tests.support.builders import experiment as experiment_builder

pytestmark = pytest.mark.behavior


@verifies("assay.pcr_error",
          claim="amplification error is counted in n_pcr_errors and not folded "
                "into the somatic mutation count",
          fields=("n_pcr_errors",))
def test_pcr_error_is_counted_separately_from_hypermutation():
    exp = experiment_builder.recombined().pcr_amplify(count=3)
    seen = 0
    for outcome in exp.run(n=40, seed=5):
        record = outcome.airr_record(exp.refdata)
        assert record["n_mutations"] == 0, (
            "no SHM was configured, so a non-zero mutation count means an "
            "amplification error was counted as somatic")
        seen += record["n_pcr_errors"]
    assert seen > 0, "no PCR error was produced; the fixture establishes nothing"


@verifies("assay.quality_error",
          claim="base-calling error is counted in n_quality_errors, separately "
                "from both amplification error and somatic mutation",
          fields=("n_quality_errors",))
def test_quality_error_is_counted_separately():
    exp = experiment_builder.recombined().sequencing_errors(count=3)
    seen = 0
    for outcome in exp.run(n=40, seed=5):
        record = outcome.airr_record(exp.refdata)
        assert record["n_mutations"] == 0 and record["n_pcr_errors"] == 0
        seen += record["n_quality_errors"]
    assert seen > 0, "no quality error was produced; the fixture establishes nothing"


@verifies("assay.ambiguous_calls",
          claim="an uncalled position appears as N and is not counted as a "
                "mutation -- the claim the mechanism's truth_gap says has no "
                "field of its own to carry it",
          fields=("sequence",))
def test_ambiguous_positions_are_not_counted_as_mutations():
    """The mechanism declares a truth gap because no N-count field exists. The
    claim can still be established against the sequence itself."""
    exp = experiment_builder.recombined().ambiguous_base_calls(count=8)
    saw_n = False
    for outcome in exp.run(n=40, seed=5):
        record = outcome.airr_record(exp.refdata)
        ns = record["sequence"].upper().count("N")
        saw_n |= ns > 0
        assert record["n_mutations"] == 0, (
            f"{ns} ambiguous positions were counted as {record['n_mutations']} "
            "mutations; an N is an absence of information, not a substitution")
    assert saw_n, "no ambiguous base was produced; the fixture establishes nothing"


@verifies("assay.orientation",
          claim="rev_comp records which strand was read, and the recorded "
                "sequence is the strand the coordinates are relative to",
          fields=("rev_comp",))
def test_orientation_is_recorded_and_consistent():
    """AIRR puts positions with respect to the reverse complement of `sequence`
    when rev_comp is true, so the flag is what makes every coordinate in the
    record interpretable."""
    exp = experiment_builder.recombined().random_strand_orientation(prob=1.0)
    for outcome in exp.run(n=20, seed=5):
        record = outcome.airr_record(exp.refdata)
        assert record["rev_comp"] is True
        assert record["v_sequence_start"] is not None
        assert 0 <= record["v_sequence_start"] < len(record["sequence"])


@verifies("assay.contamination",
          claim="a read that did not come from the sampled repertoire is flagged, "
                "so a benchmark can exclude it rather than scoring against it",
          fields=("is_contaminant",))
def test_contamination_is_flagged():
    exp = experiment_builder.recombined().contaminate(prob=1.0)
    flags = {outcome.airr_record(exp.refdata)["is_contaminant"]
             for outcome in exp.run(n=20, seed=5)}
    assert flags == {True}, f"contaminated reads were not all flagged: {flags}"


@verifies("assay.paired_end",
          claim="each read's recorded span is the length of the read, and lies "
                "within the molecule, so a tool reconstructing the fragment from "
                "R1 and R2 is scored against coordinates that hold",
          fields=("r1_sequence", "r1_start", "r1_end", "r2_sequence",
                  "r2_start", "r2_end"))
def test_each_read_is_the_stretch_its_offsets_claim():
    """Fields are read by name rather than through an f-string, so a static
    check can confirm the claim matches what the body actually inspects."""
    exp = experiment_builder.recombined().paired_end(r1_length=120, insert_size=260)
    checked = 0
    for outcome in exp.run(n=25, seed=7):
        record = outcome.airr_record(exp.refdata)
        molecule = len(record["sequence"])
        for sequence, start, end, which in (
                (record["r1_sequence"], record["r1_start"], record["r1_end"], "r1"),
                (record["r2_sequence"], record["r2_start"], record["r2_end"], "r2")):
            if not sequence or start is None or end is None:
                continue
            assert end - start == len(sequence), (
                f"{which}: offsets span {end - start} but the read is "
                f"{len(sequence)} bases")
            assert 0 <= start <= end <= molecule, (
                f"{which}: offsets {start}..{end} fall outside a "
                f"{molecule}-base molecule")
            checked += 1
    assert checked, "no paired read was produced; the fixture establishes nothing"


@verifies("assay.end_loss",
          claim="5' end loss is recorded in its own field, so sequence missing "
                "because the read did not reach it is distinguishable from "
                "germline trimming",
          fields=("end_loss_5_length",))
def test_end_loss_is_recorded_separately_from_trimming():
    exp = experiment_builder.recombined().end_loss_5prime(length=12)
    lost = 0
    for outcome in exp.run(n=30, seed=7):
        record = outcome.airr_record(exp.refdata)
        assert record["v_trim_5"] == 0, (
            "end loss was folded into the trim fields; a read that did not "
            "reach the 5' end is not exonuclease loss")
        lost += record["end_loss_5_length"]
    assert lost > 0, "no end loss was produced; the fixture establishes nothing"
