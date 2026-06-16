import pytest
import GenAIRR as ga
from GenAIRR import _engine
from GenAIRR._s5f_loader import load_builtin_s5f_kernel

S5F_MODEL = "hh_s5f"


def _founder():
    compiled = ga.Experiment.on("human_igh").recombine().compile()
    outcome = compiled.run(n=1, seed=0)[0]
    return outcome.final_simulation()


def _kernel():
    return load_builtin_s5f_kernel(S5F_MODEL)


def test_simulate_lineage_produces_valid_tree():
    founder = _founder()
    mut, sub = _kernel()
    tree = _engine.simulate_lineage(founder, mut, sub, 0.05, 1.5, 0.0, 8, 500, 30, 2024)
    tree.validate()
    assert len(tree) >= 1
    nodes = tree.nodes()
    assert nodes[0].parent_id is None
    assert nodes[0].generation == 0
    assert any(n.mutations_from_parent > 0 for n in nodes) or len(tree) == 1


def test_simulate_lineage_exports_are_wellformed():
    founder = _founder()
    mut, sub = _kernel()
    tree = _engine.simulate_lineage(founder, mut, sub, 0.05, 1.5, 0.0, 8, 500, 30, 7)
    nwk = tree.to_newick()
    assert nwk.endswith(";")
    assert nwk.count("(") == nwk.count(")")
    assert tree.to_fasta().startswith(">node")
    tsv = tree.to_node_table_tsv()
    assert tsv.splitlines()[0].startswith("node_id\t")
    assert len(tsv.splitlines()) == len(tree) + 1


def test_simulate_lineage_is_deterministic():
    founder = _founder()
    mut, sub = _kernel()
    a = _engine.simulate_lineage(founder, mut, sub, 0.05, 1.5, 0.0, 8, 500, 30, 99)
    b = _engine.simulate_lineage(founder, mut, sub, 0.05, 1.5, 0.0, 8, 500, 30, 99)
    assert a.to_newick() == b.to_newick()
    assert a.to_fasta() == b.to_fasta()


def test_simulate_lineage_rejects_bad_kernel():
    founder = _founder()
    with pytest.raises(ValueError):
        _engine.simulate_lineage(founder, [0.1] * 10, [0.1] * 10, 0.05, 1.5, 0.0, 8, 500, 30, 0)


def test_simulate_lineage_rejects_nonfinite_or_negative_kernel():
    founder = _founder()
    _mut, sub = _kernel()
    # NaN in the mutability table must raise (not panic across the FFI boundary).
    with pytest.raises(ValueError):
        _engine.simulate_lineage(
            founder, [float("nan")] * 1024, sub, 0.05, 1.5, 0.0, 8, 500, 30, 0
        )
    # Negative value in the mutability table must raise too.
    with pytest.raises(ValueError):
        _engine.simulate_lineage(
            founder, [-1.0] + [0.0] * 1023, sub, 0.05, 1.5, 0.0, 8, 500, 30, 0
        )


def test_simulate_lineage_rejects_zero_n_sample():
    founder = _founder()
    mut, sub = _kernel()
    with pytest.raises(ValueError):
        _engine.simulate_lineage(founder, mut, sub, 0.05, 1.5, 0.0, 8, 500, 0, 0)


def test_simulate_lineage_rejects_nonfinite_lambda():
    founder = _founder()
    mut, sub = _kernel()
    # NaN/inf/negative lambda_base must raise, not silently return a founder-only tree.
    for bad in (float("nan"), float("inf"), -1.0):
        with pytest.raises(ValueError):
            _engine.simulate_lineage(founder, mut, sub, 0.05, bad, 0.0, 8, 500, 30, 0)
    # lambda_mut is validated too.
    with pytest.raises(ValueError):
        _engine.simulate_lineage(founder, mut, sub, 0.05, 1.5, float("nan"), 8, 500, 30, 0)


def test_simulate_lineage_rejects_excessive_max_generations():
    founder = _founder()
    mut, sub = _kernel()
    with pytest.raises(ValueError):
        _engine.simulate_lineage(founder, mut, sub, 0.05, 1.5, 0.0, 5000, 500, 30, 0)


def test_affinity_node_getter_and_neutral_default():
    founder = _founder()
    mut, sub = _kernel()
    # Default call (no affinity args) is the neutral path; affinity getter exists.
    tree = _engine.simulate_lineage(founder, mut, sub, 0.05, 1.5, 0.0, 8, 500, 30, 2024)
    nodes = tree.nodes()
    assert all(hasattr(n, "affinity") for n in nodes)
    # neutral path leaves affinity at 0.0
    assert all(n.affinity == 0.0 for n in nodes)


def test_affinity_selection_raises_mean_affinity():
    founder = _founder()
    mut, sub = _kernel()
    founder_aa_len = 100  # target length need not match; distance uses min length
    # Use "A" (alanine) target with small beta so exp(-beta*d) stays in (0,1)
    # and doesn't underflow. "W"*100 + beta=1.0 collapses to denorm for IgH
    # sequences, making neutral and selected indistinguishable.
    target = "A" * founder_aa_len  # a fixed explicit target both runs share
    common = dict(target_aa=target, beta=0.001, mature_substitutions=5)
    # neutral: selection_strength = 0 (affinities populated but not selected on)
    neutral = _engine.simulate_lineage(
        founder, mut, sub, 0.1, 1.6, 0.0, 12, 800, 60, 4242,
        selection_strength=0.0, **common,
    )
    # strong selection toward the same target, same seed
    selected = _engine.simulate_lineage(
        founder, mut, sub, 0.1, 1.6, 0.0, 12, 800, 60, 4242,
        selection_strength=50.0, **common,
    )

    def mean_aff(tree):
        ns = tree.nodes()
        return sum(n.affinity for n in ns) / max(1, len(ns))

    # selection should enrich high-affinity cells → higher mean affinity
    assert mean_aff(selected) > mean_aff(neutral), (
        f"selected {mean_aff(selected)} !> neutral {mean_aff(neutral)}"
    )
    # affinities are populated (in (0,1]) when a model is active
    assert all(0.0 < n.affinity <= 1.0 + 1e-9 for n in selected.nodes())


def test_affinity_auto_target_is_deterministic():
    founder = _founder()
    mut, sub = _kernel()
    a = _engine.simulate_lineage(founder, mut, sub, 0.1, 1.5, 0.0, 10, 500, 40, 7,
                                 selection_strength=5.0)
    b = _engine.simulate_lineage(founder, mut, sub, 0.1, 1.5, 0.0, 10, 500, 40, 7,
                                 selection_strength=5.0)
    assert a.to_newick() == b.to_newick()
    assert [n.affinity for n in a.nodes()] == [n.affinity for n in b.nodes()]


def test_simulate_family_outcomes_extinct_founder_yields_no_observed():
    """lambda_base=0 => founder never divides => extinct => zero observed cells.

    Sampling must draw from the LIVING final-generation population, not all tree
    leaves. An extinct clone (no living cells) is not observed, so it contributes
    no node/observed Outcomes.
    """
    c = ga.Experiment.on("human_igh").recombine().compile()
    founder = c.run(n=1, seed=0)[0]  # full Outcome (not just the Simulation)
    refdata = c.refdata
    mut, sub = _kernel()
    fam = _engine.simulate_family_outcomes(
        founder, refdata, mut, sub, 0.05, 0.0, 0.0, 6, 300, 30, 2024
    )
    assert len(fam.observed_outcomes()) == 0
    assert all(o is None for o in fam.node_outcomes())
    assert all(not n.observed and n.abundance == 0 for n in fam.tree().nodes())


def test_simulate_family_outcomes_healthy_lambda_produces_observed():
    """A healthy lambda_base produces observed cells whose abundances sum to n_sample.

    Scans seeds because a single founder can go extinct even at a healthy rate.
    """
    c = ga.Experiment.on("human_igh").recombine().compile()
    founder = c.run(n=1, seed=0)[0]  # full Outcome (not just the Simulation)
    refdata = c.refdata
    mut, sub = _kernel()
    n_sample = 30
    for seed in range(20):
        fam = _engine.simulate_family_outcomes(
            founder, refdata, mut, sub, 0.05, 1.6, 0.0, 6, 300, n_sample, seed
        )
        observed = [n for n in fam.tree().nodes() if n.observed]
        if not observed:
            continue  # extinct for this seed
        assert len(fam.observed_outcomes()) == len(observed)
        assert sum(n.abundance for n in observed) == n_sample
        max_gen = max(n.generation for n in fam.tree().nodes())
        assert all(n.generation == max_gen for n in observed)
        break
    else:
        raise AssertionError("no seed in 0..20 produced a surviving family")


def test_affinity_rejects_bad_params():
    founder = _founder()
    mut, sub = _kernel()
    with pytest.raises(ValueError):
        _engine.simulate_lineage(founder, mut, sub, 0.05, 1.5, 0.0, 8, 500, 30, 0,
                                 selection_strength=float("nan"))
    with pytest.raises(ValueError):
        _engine.simulate_lineage(founder, mut, sub, 0.05, 1.5, 0.0, 8, 500, 30, 0,
                                 beta=-1.0)
