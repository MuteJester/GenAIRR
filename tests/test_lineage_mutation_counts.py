import GenAIRR as ga
from GenAIRR import _engine
from GenAIRR._s5f_loader import load_builtin_s5f_kernel


def _founder_and_refdata():
    compiled = ga.Experiment.on("human_igh").recombine().compile()
    return compiled.run(n=1, seed=0)[0], compiled.refdata


def test_lineage_record_mutation_counts_are_consistent():
    founder, refdata = _founder_and_refdata()
    mut, sub = load_builtin_s5f_kernel("hh_s5f")
    fam = _engine.simulate_family_outcomes(founder, refdata, mut, sub, 0.1, 1.6, 0.0, 8, 300, 30, 7)
    recs = fam.airr_records(refdata)   # NEW method (Task 1b)
    assert len(recs) >= 1
    saw_mutated = False
    for r in recs:
        per_seg = r["n_v_mutations"] + r["n_d_mutations"] + r["n_j_mutations"] + r["n_np_mutations"]
        assert r["n_mutations"] == per_seg, (
            f"n_mutations {r['n_mutations']} != sum per-segment {per_seg}"
        )
        v_sub = (r["n_fwr1_mutations"] + r["n_cdr1_mutations"] + r["n_fwr2_mutations"]
                 + r["n_cdr2_mutations"] + r["n_fwr3_mutations"] + r["n_v_unannotated_mutations"])
        assert r["n_v_mutations"] == v_sub, (
            f"n_v_mutations {r['n_v_mutations']} != sum subregion {v_sub}"
        )
        if r["n_mutations"] > 0:
            saw_mutated = True
    assert saw_mutated, "expected at least one mutated node record"
