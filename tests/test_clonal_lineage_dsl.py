import pytest
import GenAIRR as ga


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


def test_clonal_lineage_selection_raises_mean_affinity():
    neutral = _exp(n_clones=2, selection_strength=0.0, target_aa="A" * 100,
                   beta=0.001, max_generations=10, n_sample=40).run_records(seed=5)
    selected = _exp(n_clones=2, selection_strength=50.0, target_aa="A" * 100,
                    beta=0.001, max_generations=10, n_sample=40).run_records(seed=5)
    def mean_aff(res):
        a = [r["lineage_affinity"] for r in res.records]
        return sum(a) / max(1, len(a))
    assert mean_aff(selected) > mean_aff(neutral)


def test_clonal_lineage_validation_rejects_bad_args():
    with pytest.raises((ValueError, TypeError)):
        ga.Experiment.on("human_igh").recombine().clonal_lineage(n_clones=0)
