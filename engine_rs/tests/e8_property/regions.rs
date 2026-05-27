use super::common::{
    assembled_v_sim, assert_codon_rails_consistent, assert_frame_phases_consistent,
    assert_region_ranges_valid, uniform_s5f, vj_plan, vj_refdata, SEED_RANGE,
};
use genairr_engine::dist::{EmpiricalLengthDist, UniformBase};
use genairr_engine::ir::Simulation;
use genairr_engine::pass::testing::PassRuntime;
use genairr_engine::pass::PassPlan;
use genairr_engine::passes::{IndelPass, S5FMutationPass};

#[test]
fn property_region_ranges_valid_after_indel() {
    let mut plan = PassPlan::new();
    plan.push(Box::new(IndelPass::new(
        Box::new(EmpiricalLengthDist::from_pairs(vec![(5, 1.0)])),
        0.5,
        Box::new(UniformBase),
    )));

    for seed in 0..SEED_RANGE {
        let outcome = PassRuntime::execute(&plan, assembled_v_sim(), seed);
        assert_region_ranges_valid(outcome.final_simulation(), &format!("indel seed={}", seed));
    }
}

#[test]
fn property_region_ranges_valid_after_full_vj_pipeline() {
    let refdata = vj_refdata();
    let recomb = vj_plan(&refdata);

    for seed in 0..SEED_RANGE {
        let recomb_outcome =
            PassRuntime::execute_with_refdata(&recomb, Simulation::default(), seed, &refdata);
        let assembled = recomb_outcome.final_simulation().clone();

        let mut corrupt = PassPlan::new();
        corrupt.push(Box::new(S5FMutationPass::new(
            uniform_s5f(),
            Box::new(EmpiricalLengthDist::from_pairs(vec![(2, 1.0)])),
        )));
        corrupt.push(Box::new(IndelPass::new(
            Box::new(EmpiricalLengthDist::from_pairs(vec![(2, 1.0)])),
            0.5,
            Box::new(UniformBase),
        )));
        let corr_outcome = PassRuntime::execute(&corrupt, assembled, seed);

        assert_region_ranges_valid(
            corr_outcome.final_simulation(),
            &format!("vj+corruption seed={}", seed),
        );
        assert_frame_phases_consistent(
            corr_outcome.final_simulation(),
            &format!("vj+corruption seed={}", seed),
        );
        assert_codon_rails_consistent(
            corr_outcome.final_simulation(),
            &format!("vj+corruption seed={}", seed),
        );
    }
}

#[test]
fn property_frame_phases_valid_after_indel_in_vj_pipeline() {
    let refdata = vj_refdata();
    let recomb = vj_plan(&refdata);

    let mut indel = PassPlan::new();
    indel.push(Box::new(IndelPass::new(
        Box::new(EmpiricalLengthDist::from_pairs(vec![(5, 1.0)])),
        0.5,
        Box::new(UniformBase),
    )));

    for seed in 0..SEED_RANGE {
        let recomb_outcome =
            PassRuntime::execute_with_refdata(&recomb, Simulation::default(), seed, &refdata);
        let outcome = PassRuntime::execute(
            &indel,
            recomb_outcome.final_simulation().clone(),
            seed ^ 0x9e37_79b9,
        );

        assert_frame_phases_consistent(
            outcome.final_simulation(),
            &format!("vj+indel seed={}", seed),
        );
        assert_codon_rails_consistent(
            outcome.final_simulation(),
            &format!("vj+indel seed={}", seed),
        );
    }
}
