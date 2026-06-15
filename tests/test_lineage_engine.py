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
