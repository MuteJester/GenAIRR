"""
Correctness test: per-node Outcomes from simulate_family_outcomes must be
self-consistent so that native build_airr_record + validate_record produce
correct mutation-count fields and pass mutation-count invariants.

Bug fixed: node Outcomes were built with an empty event ledger but a
nonzero Simulation.mutation_count, so validate_record would false-fail the
mutation-count sum invariant on mutated nodes.

Note: AlleleCallTieSetMismatch issues on D can still appear for SHM-mutated
nodes (SHM changes D-segment match scores, making formerly-unique calls
ambiguous). These are a biological reality of the lineage path, not caused
by the event-ledger bug, and are not checked here.
"""
import GenAIRR as ga
from GenAIRR import _engine
from GenAIRR._s5f_loader import load_builtin_s5f_kernel

# Mutation-count issue kinds that this fix must eliminate.
_MUTATION_COUNT_KINDS = {
    "NMutationsMismatch",
    "NVMutationsMismatch",
    "NDMutationsMismatch",
    "NJMutationsMismatch",
    "NNpMutationsMismatch",
    "MutationCountSumMismatch",
    "NFwr1MutationsMismatch",
    "NCdr1MutationsMismatch",
    "NFwr2MutationsMismatch",
    "NCdr2MutationsMismatch",
    "NFwr3MutationsMismatch",
    "NVUnannotatedMutationsMismatch",
    "VSubregionMutationCountSumMismatch",
}


def _founder_refdata():
    c = ga.Experiment.on("human_igh").recombine().compile()
    return c.run(n=1, seed=0)[0], c.refdata


def test_node_outcomes_no_mutation_count_issues():
    """After the fix, no mutation-count validation issues on any lineage node."""
    founder, refdata = _founder_refdata()
    mut, sub = load_builtin_s5f_kernel("hh_s5f")
    fam = _engine.simulate_family_outcomes(founder, mut, sub, 0.05, 1.6, 0.0, 6, 300, 30, 7)
    saw_mut = False
    for o in fam.observed_outcomes():
        issues = o.validate_record(refdata)
        mut_issues = [i for i in issues if i["kind"] in _MUTATION_COUNT_KINDS]
        assert len(mut_issues) == 0, f"mutation-count validate_record issues: {mut_issues}"
        rec = o.airr_record(refdata)
        if rec["n_mutations"] > 0:
            saw_mut = True
            per_seg = (rec["n_v_mutations"] + rec["n_d_mutations"]
                       + rec["n_j_mutations"] + rec["n_np_mutations"])
            assert rec["n_mutations"] == per_seg, (
                f"n_mutations {rec['n_mutations']} != per-segment sum {per_seg}"
            )
            v_sub = (rec["n_fwr1_mutations"] + rec["n_cdr1_mutations"]
                     + rec["n_fwr2_mutations"] + rec["n_cdr2_mutations"]
                     + rec["n_fwr3_mutations"] + rec["n_v_unannotated_mutations"])
            assert rec["n_v_mutations"] == v_sub, (
                f"n_v_mutations {rec['n_v_mutations']} != subregion sum {v_sub}"
            )
    assert saw_mut, "expected at least one mutated observed node"


def test_airr_records_mutation_counts_consistent():
    """PyFamilyOutcome.airr_records() produces consistent mutation counts."""
    founder, refdata = _founder_refdata()
    mut, sub = load_builtin_s5f_kernel("hh_s5f")
    fam = _engine.simulate_family_outcomes(founder, mut, sub, 0.05, 1.6, 0.0, 6, 300, 30, 7)
    recs = fam.airr_records(refdata)
    assert len(recs) >= 1
    saw_mut = False
    for rec in recs:
        per_seg = (rec["n_v_mutations"] + rec["n_d_mutations"]
                   + rec["n_j_mutations"] + rec["n_np_mutations"])
        assert rec["n_mutations"] == per_seg, (
            f"n_mutations {rec['n_mutations']} != per-segment sum {per_seg}"
        )
        if rec["n_mutations"] > 0:
            saw_mut = True
    assert saw_mut, "expected at least one mutated node record"
