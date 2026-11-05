"""The protein affinity selection actually scores.

`sim_to_aa` concatenated the per-region codon rail, and a codon is kept only if
it fits wholly inside one region (`ir/region.rs:112`). Every region boundary in
a rearrangement is a junction boundary, so the residues selection could not see
were exactly the ones it exists to select on: measured over 200 `human_igh`
records, 2.68 codons dropped per record and 536 of 536 of them straddling
V/NP1, NP1/D, D/NP2 or NP2/J.

The oracle here is `sequence_aa`, which the record projection derives from
`coding_frame_start`/`coding_frame_end` -- a different code path from the one
under test, which is what makes it evidence rather than a tautology.
"""
import pytest

import GenAIRR as ga
import GenAIRR.data as gdata

pytestmark = pytest.mark.behavior


def _lineage(target=None, *, seed=1, rate=0.0):
    return (ga.Experiment.on(gdata.HUMAN_IGH_OGRDB).recombine()
            .clonal_lineage(n_clones=1, max_generations=4, n_max=50, n_sample=3,
                            rate=rate, lambda_base=2.0, selection_strength=0.0,
                            target_aa=target)
            .compile().run_records(seed=seed).records)


@pytest.fixture(scope="module")
def unmutated():
    """A lineage grown with no SHM, so every cell's protein is the founder's."""
    records = _lineage()
    protein = records[0]["sequence_aa"]
    # Guard the premise twice. A target may not contain a stop (the DSL refuses
    # it), and the whole point of rate=0 is that there is one protein in the
    # family -- if cells differed, scoring a cell against the founder's protein
    # would not have to give 1.0.
    assert "*" not in protein and "X" not in protein, (
        f"this founder is not translatable end to end: {protein}")
    assert len({r["sequence_aa"] for r in records}) == 1, "cells differ despite rate=0"
    return protein


def test_a_cell_scored_against_its_own_protein_is_a_perfect_match(unmutated):
    """The sharpest statement of the defect: it reported 3.57e-67.

    With `rate=0` a sampled cell *is* the founder, so scored against the
    founder's own protein its affinity can only be 1.0. It was not, because the
    cell was translated by dropping its junction-boundary codons while the
    target was the real protein, leaving the two misaligned from CDR3 onward.
    """
    records = _lineage(target=unmutated)
    assert records[0]["sequence_aa"] == unmutated, "the founder moved; nothing was compared"
    affinity = records[0]["lineage_affinity"]
    assert affinity == pytest.approx(1.0), (
        f"a cell identical to the target scored {affinity!r}")


def test_every_cell_of_an_unmutated_family_is_a_perfect_match(unmutated):
    """Not just the first record: no cell in the family may be scored short."""
    for record in _lineage(target=unmutated):
        assert record["lineage_affinity"] == pytest.approx(1.0), (
            f"cell {record.get('lineage_node_id')} scored "
            f"{record['lineage_affinity']!r}")


def test_a_target_the_wrong_length_is_refused(unmutated):
    """A misaligned target produced plausible numbers that meant nothing.

    Everything past the shorter of the two sequences was ignored, so a truncated
    or over-long target scored silently rather than being reported.
    """
    with pytest.raises(ValueError, match=r"target_aa"):
        _lineage(target=unmutated[:-5])
    with pytest.raises(ValueError, match=r"target_aa"):
        _lineage(target=unmutated + "AAAAA")
