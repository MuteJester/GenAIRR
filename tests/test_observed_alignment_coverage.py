"""Explicit read clipping changes coverage, not structural deletion burden."""

import csv

import pytest

import GenAIRR as ga
from ._p_productive_support import reverse_complement, synthetic_cartridge


@pytest.mark.parametrize("reverse", [False, True])
@pytest.mark.parametrize("loss5,loss3", [(2, 0), (0, 2), (2, 2)])
def test_read_clipping_projects_observed_reference_bounds(loss5, loss3, reverse, tmp_path):
    cfg = synthetic_cartridge({})
    exp = ga.Experiment.on(cfg).recombine(np1_lengths=[(0, 1)]).trim(enabled=False)
    if loss5:
        exp.end_loss_5prime(length=loss5)
    if loss3:
        exp.end_loss_3prime(length=loss3)
    if reverse:
        exp.random_strand_orientation(prob=1.0)
    result = exp.run_records(n=1, seed=0)
    record = result[0]
    expected = "AAACCCTGTTTTGGG"[loss5:15 - loss3]
    assert record["sequence"].upper() == (reverse_complement(expected) if reverse else expected)
    assert record["sequence_alignment"] == expected
    assert record["germline_alignment"] == expected
    assert record["v_cigar"] == f"{9 - loss5}M"
    assert record["j_cigar"] == f"{6 - loss3}M"
    assert (record["v_germline_start"], record["v_germline_end"]) == (loss5, 9)
    assert (record["j_germline_start"], record["j_germline_end"]) == (0, 6 - loss3)
    assert record["v_identity"] == record["j_identity"] == 1.0
    assert record["v_trim_5"] == record["j_trim_3"] == 0
    assert (record["end_loss_5_length"], record["end_loss_3_length"]) == (loss5, loss3)
    assert result.validate_records(exp.refdata).ok

    output = tmp_path / "observed-coverage.tsv"
    result.to_tsv(str(output))
    with output.open(newline="") as stream:
        exported = next(csv.DictReader(stream, delimiter="\t"))
    assert int(exported["v_germline_start"]) == loss5 + 1
    assert int(exported["j_germline_end"]) == 6 - loss3
    assert exported["v_cigar"] == (f"{loss5}N" if loss5 else "") + f"{9 - loss5}M"
    assert exported["j_cigar"] == f"{9 - loss5}S{6 - loss3}M"
