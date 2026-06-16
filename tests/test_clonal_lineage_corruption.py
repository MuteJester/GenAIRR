"""Library-prep / sequencing artefact passes applied to clonal_lineage reads.

``clonal_lineage`` produces per-observed-cell AIRR records. Prior to this
slice the compiler rejected ALL steps following the lineage fork, so the
reads were pristine — no sequencing / PCR / indel noise, a regression vs
``expand_clones``. These tests pin the corruption-aware path: each observed
cell now gets independent corruption, and its record reports BOTH the
founder recombination provenance (trims, v/d/j) AND the SHM counts AND the
corruption counts (``n_quality_errors`` etc.).
"""
import pytest

import GenAIRR as ga


def test_clonal_lineage_with_corruption_applies_artefacts():
    exp = (
        ga.Experiment.on("human_igh")
        .recombine()
        .clonal_lineage(
            n_clones=3, max_generations=4, n_max=200, n_sample=20, rate=0.02
        )
        .sequencing_errors(rate=0.05)
    )
    result = exp.run_records(seed=0)
    assert len(result.records) > 0
    # Corruption was actually applied to at least one observed cell.
    assert any(r.get("n_quality_errors", 0) > 0 for r in result.records)
    for r in result.records:
        # Founder recombination provenance survives the merge.
        assert r["v_call"]
        assert r["clone_id"] in (0, 1, 2)
        assert "lineage_node_id" in r
        # SHM per-segment counts stay self-consistent after corruption merge.
        assert r["n_mutations"] == (
            r["n_v_mutations"]
            + r["n_d_mutations"]
            + r["n_j_mutations"]
            + r["n_np_mutations"]
        )


def test_clonal_lineage_rejects_mutate_after():
    with pytest.raises(Exception):
        (
            ga.Experiment.on("human_igh")
            .recombine()
            .clonal_lineage(n_clones=2, n_sample=5)
            .mutate(rate=0.05)
            .compile()
        )


def test_clonal_lineage_rejects_paired_end_after():
    with pytest.raises(Exception):
        (
            ga.Experiment.on("human_igh")
            .recombine()
            .clonal_lineage(n_clones=2, n_sample=5)
            .paired_end(r1_length=150, insert_size=300)
            .compile()
        )


def test_clonal_lineage_without_corruption_still_works():
    result = (
        ga.Experiment.on("human_igh")
        .recombine()
        .clonal_lineage(n_clones=2, n_sample=10, rate=0.02)
        # seed chosen so the families survive: sampling draws from the living
        # final generation, so an all-extinct seed (e.g. 0) yields zero records.
        .run_records(seed=1)
    )
    assert len(result.records) > 0
    # No corruption pass → no quality errors stamped.
    assert all(r.get("n_quality_errors", 0) == 0 for r in result.records)
