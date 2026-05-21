use super::*;

#[test]
fn s5f_mutation_pass_pool_reflects_recorded_mutations() {
    let mut plan = PassPlan::new();
    plan.push(Box::new(S5FMutationPass::new(
        s5f_uniform_kernel(),
        Box::new(EmpiricalLengthDist::from_pairs(vec![(5, 1.0)])),
    )));
    let outcome = PassRuntime::execute(&plan, s5f_test_sim(), 99);
    let final_sim = outcome.final_simulation();

    let mut last_at_site: std::collections::HashMap<u32, u8> = std::collections::HashMap::new();
    for i in 0..5 {
        let s = match outcome
            .trace
            .find(&format!("mutate.s5f.site[{}]", i))
            .unwrap()
            .value
        {
            ChoiceValue::Int(n) => n as u32,
            _ => unreachable!(),
        };
        let b = match outcome
            .trace
            .find(&format!("mutate.s5f.base[{}]", i))
            .unwrap()
            .value
        {
            ChoiceValue::Base(b) => b,
            _ => unreachable!(),
        };
        last_at_site.insert(s, b);
    }
    for (&site, &expected_base) in last_at_site.iter() {
        let actual = final_sim.pool.get(NucHandle::new(site)).unwrap().base;
        assert_eq!(
            actual, expected_base,
            "trace says site {} got base {}, but pool has {}",
            site, expected_base as char, actual as char
        );
    }
}

#[test]
fn s5f_mutation_pass_refreshes_codon_rail() {
    let mut plan = PassPlan::new();
    plan.push(Box::new(S5FMutationPass::new(
        s5f_uniform_kernel(),
        Box::new(EmpiricalLengthDist::from_pairs(vec![(8, 1.0)])),
    )));
    let outcome = PassRuntime::execute(&plan, s5f_test_sim(), 42);
    let final_sim = outcome.final_simulation();

    let stored_aa = &final_sim.sequence.regions[0].amino_acids;
    let fresh = final_sim.sequence.regions[0].with_codon_rail_recomputed(&final_sim.pool);
    assert_eq!(stored_aa, &fresh.amino_acids);
    assert_eq!(
        final_sim.sequence.regions[0].stop_codon_positions,
        fresh.stop_codon_positions
    );
}
