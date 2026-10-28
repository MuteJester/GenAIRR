"""Alignment truth follows source provenance even when SHM changes inference."""
import re

import pytest

import GenAIRR as ga

pytestmark = pytest.mark.behavior


@pytest.mark.parametrize("inverted", [False, True])
def test_high_shm_alignment_reconstructs_every_generating_reference_position(inverted):
    exp = ga.Experiment.on("human_igh").recombine()
    if inverted:
        exp.invert_d(prob=1)
    compiled = exp.mutate(rate=0.15).polymerase_indels(count=3).compile()
    result = compiled.run_records(n=400, seed=11, expose_provenance=True)
    recalls_changed = 0
    deletions = 0
    for outcome, record in zip(result.outcomes, result):
        sim = outcome.final_simulation()
        regions = {r.segment.lower(): r for r in sim.regions()}
        for segment in ("v", "d", "j"):
            call = record[f"{segment}_call"]
            assert call == record[f"truth_{segment}_call"]
            assert "," not in call
            recalls_changed += call not in record[f"recalled_{segment}_call"].split(",")
            allele_id = getattr(sim, f"{segment}_allele_id")()
            allele = getattr(compiled.refdata, f"{segment}_allele")(allele_id)
            reference = bytes(allele.seq()).decode().upper()
            start, end = record[f"{segment}_germline_start"], record[f"{segment}_germline_end"]
            assert start == record[f"{segment}_trim_5"]
            assert end == len(reference) - record[f"{segment}_trim_3"]
            expected_positions = list(range(start, end))
            if segment == "d" and inverted:
                expected_positions.reverse()
                reference = reference.translate(str.maketrans("ACGT", "TGCA"))[::-1]
                expected_ref = reference[len(reference) - end:len(reference) - start]
            else:
                expected_ref = reference[start:end]
            a, b = record[f"{segment}_alignment_start"], record[f"{segment}_alignment_end"]
            assert record["germline_alignment"][a:b].replace("-", "").upper() == expected_ref.upper()
            region = regions[segment]
            assert record[f"{segment}_sequence_start"] == region.start
            assert record[f"{segment}_sequence_end"] == region.end
            query = region.start
            ref_index = 0
            for amount, op in re.findall(r"(\d+)([MID])", record[f"{segment}_cigar"]):
                for _ in range(int(amount)):
                    if op == "M":
                        assert sim.germline_position(query) == expected_positions[ref_index]
                    elif op == "I":
                        assert sim.germline_position(query) is None
                    else:
                        deletions += 1
                    query += op in "MI"
                    ref_index += op in "MD"
            assert query == region.end
            assert ref_index == len(expected_positions)
    assert recalls_changed > 0, "exercise actual post-SHM reference disagreements"
    assert deletions > 0, "exercise source positions removed by indel corruption"
    assert result.validate_records(compiled.refdata).ok


def test_unmutated_np_lengths_are_sampled_lengths_despite_boundary_matches():
    exp = ga.Experiment.on("human_igh").recombine(
        np1_lengths=[(3, 1)], np2_lengths=[(3, 1)],
    )
    result = exp.run_records(n=1000, seed=11)
    assert all(r["np1_length"] == r["np2_length"] == 3 for r in result)
    for outcome, record in zip(result.outcomes, result):
        sim = outcome.final_simulation()
        for region in sim.regions():
            if region.segment in {"NP1", "NP2"}:
                assert record[region.segment.lower()].upper() == bytes(sim.bases())[region.start:region.end].decode().upper()
