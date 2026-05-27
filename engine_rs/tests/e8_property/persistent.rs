use super::common::{assembled_v_sim, SEED_RANGE};
use genairr_engine::dist::{EmpiricalLengthDist, UniformBase};
use genairr_engine::ir::compute_codon_rail;
use genairr_engine::pass::testing::PassRuntime;
use genairr_engine::pass::PassPlan;
use genairr_engine::passes::{IndelPass, UniformMutationPass};

#[test]
fn property_persistent_ir_uniform_mutation() {
    for seed in 0..SEED_RANGE {
        let sim = assembled_v_sim();
        let pre_len = sim.pool.len();
        let pre_region_end = sim.sequence.regions[0].end.index();
        let pre_amino = compute_codon_rail(&sim.sequence.regions[0], &sim.pool)
            .amino_acids
            .clone();

        let mut plan = PassPlan::new();
        plan.push(Box::new(UniformMutationPass::new(
            Box::new(EmpiricalLengthDist::from_pairs(vec![(3, 1.0)])),
            Box::new(UniformBase),
        )));
        let _ = PassRuntime::execute(&plan, sim.clone(), seed);

        assert_eq!(sim.pool.len(), pre_len);
        assert_eq!(sim.sequence.regions[0].end.index(), pre_region_end);
        assert_eq!(
            compute_codon_rail(&sim.sequence.regions[0], &sim.pool).amino_acids,
            pre_amino
        );
    }
}

#[test]
fn property_persistent_ir_indel_pass() {
    for seed in 0..SEED_RANGE {
        let sim = assembled_v_sim();
        let pre_len = sim.pool.len();
        let pre_region_end = sim.sequence.regions[0].end.index();
        let pre_amino = compute_codon_rail(&sim.sequence.regions[0], &sim.pool)
            .amino_acids
            .clone();

        let mut plan = PassPlan::new();
        plan.push(Box::new(IndelPass::new(
            Box::new(EmpiricalLengthDist::from_pairs(vec![(4, 1.0)])),
            0.5,
            Box::new(UniformBase),
        )));
        let _ = PassRuntime::execute(&plan, sim.clone(), seed);

        assert_eq!(sim.pool.len(), pre_len);
        assert_eq!(sim.sequence.regions[0].end.index(), pre_region_end);
        assert_eq!(
            compute_codon_rail(&sim.sequence.regions[0], &sim.pool).amino_acids,
            pre_amino
        );
    }
}
