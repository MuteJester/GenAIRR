use super::common::{
    assembled_v_sim, assert_codon_rails_consistent, assert_region_ranges_valid, uniform_s5f,
    SEED_RANGE,
};
use genairr_engine::dist::{EmpiricalLengthDist, UniformBase};
use genairr_engine::pass::testing::PassRuntime;
use genairr_engine::pass::PassPlan;
use genairr_engine::passes::{
    ContaminantPass, IndelPass, PCRErrorPass, QualityErrorPass, S5FMutationPass,
    UniformMutationPass,
};

#[test]
fn property_codon_rail_consistency_uniform_mutation() {
    let mut plan = PassPlan::new();
    plan.push(Box::new(UniformMutationPass::new(
        Box::new(EmpiricalLengthDist::from_pairs(vec![(3, 1.0)])),
        Box::new(UniformBase),
    )));

    for seed in 0..SEED_RANGE {
        let outcome = PassRuntime::execute(&plan, assembled_v_sim(), seed);
        assert_codon_rails_consistent(
            outcome.final_simulation(),
            &format!("UniformMutation seed={}", seed),
        );
    }
}

#[test]
fn property_codon_rail_consistency_s5f() {
    let mut plan = PassPlan::new();
    plan.push(Box::new(S5FMutationPass::new(
        uniform_s5f(),
        Box::new(EmpiricalLengthDist::from_pairs(vec![(3, 1.0)])),
    )));

    for seed in 0..SEED_RANGE {
        let outcome = PassRuntime::execute(&plan, assembled_v_sim(), seed);
        assert_codon_rails_consistent(
            outcome.final_simulation(),
            &format!("S5FMutation seed={}", seed),
        );
    }
}

#[test]
fn property_codon_rail_consistency_pcr() {
    let mut plan = PassPlan::new();
    plan.push(Box::new(PCRErrorPass::new(
        Box::new(EmpiricalLengthDist::from_pairs(vec![(3, 1.0)])),
        Box::new(UniformBase),
    )));

    for seed in 0..SEED_RANGE {
        let outcome = PassRuntime::execute(&plan, assembled_v_sim(), seed);
        assert_codon_rails_consistent(
            outcome.final_simulation(),
            &format!("PCRError seed={}", seed),
        );
    }
}

#[test]
fn property_codon_rail_consistency_quality() {
    let mut plan = PassPlan::new();
    plan.push(Box::new(QualityErrorPass::new(
        Box::new(EmpiricalLengthDist::from_pairs(vec![(3, 1.0)])),
        Box::new(UniformBase),
    )));

    for seed in 0..SEED_RANGE {
        let outcome = PassRuntime::execute(&plan, assembled_v_sim(), seed);
        assert_codon_rails_consistent(
            outcome.final_simulation(),
            &format!("QualityError seed={}", seed),
        );
    }
}

#[test]
fn property_codon_rail_consistency_contaminant() {
    let mut plan = PassPlan::new();
    plan.push(Box::new(ContaminantPass::new(0.5, Box::new(UniformBase))));

    for seed in 0..SEED_RANGE {
        let outcome = PassRuntime::execute(&plan, assembled_v_sim(), seed);
        assert_codon_rails_consistent(
            outcome.final_simulation(),
            &format!("Contaminant seed={}", seed),
        );
    }
}

#[test]
fn property_codon_rail_consistency_indel() {
    let mut plan = PassPlan::new();
    plan.push(Box::new(IndelPass::new(
        Box::new(EmpiricalLengthDist::from_pairs(vec![(3, 1.0)])),
        0.5,
        Box::new(UniformBase),
    )));

    for seed in 0..SEED_RANGE {
        let outcome = PassRuntime::execute(&plan, assembled_v_sim(), seed);
        assert_codon_rails_consistent(outcome.final_simulation(), &format!("Indel seed={}", seed));
    }
}

#[test]
fn property_codon_rail_consistency_full_corruption_stack() {
    let mut plan = PassPlan::new();
    plan.push(Box::new(S5FMutationPass::new(
        uniform_s5f(),
        Box::new(EmpiricalLengthDist::from_pairs(vec![(2, 1.0)])),
    )));
    plan.push(Box::new(PCRErrorPass::new(
        Box::new(EmpiricalLengthDist::from_pairs(vec![(1, 1.0)])),
        Box::new(UniformBase),
    )));
    plan.push(Box::new(IndelPass::new(
        Box::new(EmpiricalLengthDist::from_pairs(vec![(2, 1.0)])),
        0.5,
        Box::new(UniformBase),
    )));
    plan.push(Box::new(QualityErrorPass::new(
        Box::new(EmpiricalLengthDist::from_pairs(vec![(1, 1.0)])),
        Box::new(UniformBase),
    )));
    plan.push(Box::new(ContaminantPass::new(0.1, Box::new(UniformBase))));

    for seed in 0..SEED_RANGE {
        let outcome = PassRuntime::execute(&plan, assembled_v_sim(), seed);
        assert_codon_rails_consistent(
            outcome.final_simulation(),
            &format!("full-corruption seed={}", seed),
        );
        assert_region_ranges_valid(
            outcome.final_simulation(),
            &format!("full-corruption seed={}", seed),
        );
    }
}
