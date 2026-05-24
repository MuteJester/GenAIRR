use super::common::{assembled_v_sim, SEED_RANGE};
use genairr_engine::dist::{EmpiricalLengthDist, UniformBase};
use genairr_engine::ir::NucHandle;
use genairr_engine::pass::PassPlan;
use genairr_engine::pass::testing::PassRuntime;
use genairr_engine::passes::UniformMutationPass;
use genairr_engine::trace::ChoiceValue;

#[test]
fn property_uniform_mutation_trace_faithfulness() {
    let count = 5usize;
    let mut plan = PassPlan::new();
    plan.push(Box::new(UniformMutationPass::new(
        Box::new(EmpiricalLengthDist::from_pairs(vec![(count as i64, 1.0)])),
        Box::new(UniformBase),
    )));

    for seed in 0..SEED_RANGE {
        let outcome = PassRuntime::execute(&plan, assembled_v_sim(), seed);
        let final_sim = outcome.final_simulation();

        let mut last: std::collections::HashMap<u32, u8> = std::collections::HashMap::new();
        for i in 0..count {
            let site = match outcome
                .trace
                .find(&format!("mutate.uniform.site[{}]", i))
                .unwrap()
                .value
            {
                ChoiceValue::Int(s) => s as u32,
                _ => unreachable!(),
            };
            let base = match outcome
                .trace
                .find(&format!("mutate.uniform.base[{}]", i))
                .unwrap()
                .value
            {
                ChoiceValue::Base(b) => b,
                _ => unreachable!(),
            };
            last.insert(site, base);
        }

        for (&site, &expected) in last.iter() {
            assert_eq!(
                final_sim.pool.get(NucHandle::new(site)).unwrap().base,
                expected,
                "seed {}: pool base at site {} doesn't match traced {} ",
                seed,
                site,
                expected as char
            );
        }
    }
}
