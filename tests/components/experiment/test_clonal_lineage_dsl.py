import pytest
import GenAIRR as ga

pytestmark = pytest.mark.behavior


def _exp(**kw):
    base = dict(n_clones=3, max_generations=6, n_max=200, n_sample=20,
                rate=0.05, lambda_base=1.6)
    base.update(kw)
    return ga.Experiment.on("human_igh").recombine().clonal_lineage(**base)


def test_clonal_lineage_runs_and_tags_records():
    result = _exp().run_records(seed=0)
    # Sampling draws from the LIVING final-generation population, so an extinct
    # clone (founder drew 0 offspring) contributes no records. By default the
    # founder-survival guard retries extinct clones with fresh deterministic
    # sub-seeds, so every requested clone survives and all clone_ids are present.
    assert len(result.records) > 0
    cids = {r["clone_id"] for r in result.records}
    assert cids == {0, 1, 2}
    for r in result.records:
        assert r["v_call"]            # real recombination provenance
        assert r["sequence"]
        assert "lineage_node_id" in r
        assert "lineage_generation" in r
        assert "lineage_abundance" in r
        assert "lineage_affinity" in r
        # mutation counts are self-consistent (pool-derived)
        per_seg = (r["n_v_mutations"] + r["n_d_mutations"]
                   + r["n_j_mutations"] + r["n_np_mutations"])
        assert r["n_mutations"] == per_seg


def test_clonal_lineage_exposes_per_clone_trees():
    result = _exp(n_clones=2).run_records(seed=1)
    trees = result.lineage_trees
    assert trees is not None and len(trees) == 2
    for t in trees:
        t.validate()                  # raises if malformed
        assert t.to_newick().endswith(";")
        assert t.to_fasta().startswith(">node")


def test_clonal_lineage_selection_populates_affinities_in_range():
    """The DSL plumbs `selection_strength` through and affinities come back.

    This test used to assert that selection raises the mean affinity, with
    `target_aa="A" * 100`. That target was shorter than the receptors it scored,
    so the distance compared the first 100 positions and ignored the rest, and
    the resulting affinities underflowed to near zero for every cell -- which is
    what made the comparison stable.

    With the translation corrected the distances are real, and the claim is no
    longer stably testable here: at beta=0.001 the mean over 7 to 17 sampled
    cells is dominated by sampling noise (0.637, 0.642, 0.634, 0.631, 0.628,
    0.651 across strengths 0.05 to 0.7), and raising the strength enough to
    separate the arms drives fitness to zero for diverged cells and extinguishes
    the family. The claim is kept where the founder can be pinned and the target
    given the right length: `test_affinity_selection_raises_mean_affinity` in
    tests/components/clonal/test_lineage_engine.py.

    What is left here is what only the DSL can break: the keyword reaching the
    engine, and affinities arriving populated and in range.
    """
    selected = _exp(n_clones=2, selection_strength=0.5, beta=0.001,
                    mature_substitutions=5, max_generations=10,
                    n_sample=40).run_records(seed=5)
    affinities = [r["lineage_affinity"] for r in selected.records]
    # Guard the premise: no records means the assertion below is vacuous.
    assert affinities, "no cells were sampled"
    assert all(0.0 < a <= 1.0 for a in affinities), (
        f"affinity outside (0, 1]: {sorted(affinities)[:3]}")
    # A neutral run with no target builds no model at all and reports 0.0, so a
    # populated affinity is itself evidence the selection keyword arrived.
    neutral = _exp(n_clones=2, selection_strength=0.0, max_generations=10,
                   n_sample=40).run_records(seed=5)
    assert all(r["lineage_affinity"] == 0.0 for r in neutral.records)


def test_clonal_lineage_validation_rejects_bad_args():
    with pytest.raises((ValueError, TypeError)):
        ga.Experiment.on("human_igh").recombine().clonal_lineage(n_clones=0)
