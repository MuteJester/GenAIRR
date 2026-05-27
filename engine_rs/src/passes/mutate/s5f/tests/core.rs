use super::*;

#[test]
fn s5f_mutation_pass_zero_count_is_noop() {
    let mut plan = PassPlan::new();
    plan.push(Box::new(S5FMutationPass::new(
        s5f_uniform_kernel(),
        Box::new(EmpiricalLengthDist::from_pairs(vec![(0, 1.0)])),
    )));

    let sim = s5f_test_sim();
    let outcome = PassRuntime::execute(&plan, sim.clone(), 0);

    assert_eq!(outcome.final_simulation().pool.len(), sim.pool.len());
    for i in 0..sim.pool.len() {
        assert_eq!(
            outcome
                .final_simulation()
                .pool
                .get(NucHandle::new(i as u32))
                .unwrap()
                .base,
            sim.pool.get(NucHandle::new(i as u32)).unwrap().base
        );
    }
    assert_eq!(outcome.trace.len(), 1);
    assert_eq!(
        outcome.trace.find("mutate.s5f.count").unwrap().value,
        ChoiceValue::Int(0)
    );
}

#[test]
fn s5f_mutation_pass_short_pool_emits_no_mutations() {
    let mut sim = Simulation::new();
    for (i, b) in b"AAAA".iter().enumerate() {
        let (next, _) = sim.with_nucleotide_pushed(Nucleotide::germline(*b, i as u16, Segment::V));
        sim = next;
    }

    let mut plan = PassPlan::new();
    plan.push(Box::new(S5FMutationPass::new(
        s5f_uniform_kernel(),
        Box::new(EmpiricalLengthDist::from_pairs(vec![(5, 1.0)])),
    )));
    let outcome = PassRuntime::execute(&plan, sim, 0);

    assert_eq!(
        outcome.trace.find("mutate.s5f.count").unwrap().value,
        ChoiceValue::Int(5)
    );
    for i in 0..5 {
        assert!(outcome
            .trace
            .find(&format!("mutate.s5f.site[{}]", i))
            .is_none());
    }
}

#[test]
fn s5f_mutation_pass_zero_mutability_kernel_emits_no_mutations() {
    let mut plan = PassPlan::new();
    plan.push(Box::new(S5FMutationPass::new(
        s5f_zero_kernel(),
        Box::new(EmpiricalLengthDist::from_pairs(vec![(10, 1.0)])),
    )));
    let outcome = PassRuntime::execute(&plan, s5f_test_sim(), 0);

    assert_eq!(
        outcome.trace.find("mutate.s5f.count").unwrap().value,
        ChoiceValue::Int(10)
    );
    for i in 0..10 {
        assert!(
            outcome
                .trace
                .find(&format!("mutate.s5f.site[{}]", i))
                .is_none(),
            "expected no mutation at index {}",
            i
        );
    }
}

#[test]
fn s5f_mutation_pass_applies_n_mutations_with_uniform_kernel() {
    let mut plan = PassPlan::new();
    plan.push(Box::new(S5FMutationPass::new(
        s5f_uniform_kernel(),
        Box::new(EmpiricalLengthDist::from_pairs(vec![(7, 1.0)])),
    )));
    let outcome = PassRuntime::execute(&plan, s5f_test_sim(), 1234);

    assert_eq!(outcome.trace.len(), 15);
    for i in 0..7 {
        let site = outcome
            .trace
            .find(&format!("mutate.s5f.site[{}]", i))
            .unwrap();
        let base = outcome
            .trace
            .find(&format!("mutate.s5f.base[{}]", i))
            .unwrap();
        match site.value {
            ChoiceValue::Int(s) => {
                assert!(s >= 2 && s < 18, "site {} out of range [2, 18)", s);
            }
            _ => panic!("wrong variant"),
        }
        match base.value {
            ChoiceValue::Base(b) => assert!(matches!(b, b'A' | b'C' | b'G' | b'T')),
            _ => panic!("wrong variant"),
        }
    }
}

#[test]
fn s5f_mutation_pass_is_deterministic_under_same_seed() {
    let plan = || {
        let mut p = PassPlan::new();
        p.push(Box::new(S5FMutationPass::new(
            s5f_uniform_kernel(),
            Box::new(EmpiricalLengthDist::from_pairs(vec![(6, 1.0)])),
        )));
        p
    };

    let oa = PassRuntime::execute(&plan(), s5f_test_sim(), 0xc0ff_ee);
    let ob = PassRuntime::execute(&plan(), s5f_test_sim(), 0xc0ff_ee);
    assert_eq!(oa.trace.choices(), ob.trace.choices());
    for i in 0..oa.final_simulation().pool.len() {
        let h = NucHandle::new(i as u32);
        assert_eq!(
            oa.final_simulation().pool.get(h).unwrap().base,
            ob.final_simulation().pool.get(h).unwrap().base
        );
    }
}

#[test]
fn s5f_mutation_pass_targets_mutability_hotspot() {
    let mut mu = vec![0.0; S5F_NUM_CONTEXTS];
    mu[0] = 1.0;
    let mut sub = vec![0.0; S5F_SUBSTITUTION_LEN];
    sub[3] = 1.0;
    let kernel = S5FKernel::new(mu, sub);

    let mut sim = Simulation::new();
    for i in 0..16 {
        let (next, _) =
            sim.with_nucleotide_pushed(Nucleotide::germline(b'A', i as u16, Segment::V));
        sim = next;
    }

    let mut plan = PassPlan::new();
    plan.push(Box::new(S5FMutationPass::new(
        kernel,
        Box::new(EmpiricalLengthDist::from_pairs(vec![(3, 1.0)])),
    )));
    let outcome = PassRuntime::execute(&plan, sim, 0);

    for i in 0..3 {
        let base_addr = format!("mutate.s5f.base[{}]", i);
        if let Some(rec) = outcome.trace.find(&base_addr) {
            match rec.value {
                ChoiceValue::Base(b) => assert_eq!(
                    b, b'T',
                    "mutation {} produced base {} (expected T)",
                    i, b as char
                ),
                _ => panic!("wrong variant"),
            }
        }
    }
}

#[test]
fn s5f_mutation_pass_declared_choices() {
    let pass = S5FMutationPass::new(
        s5f_uniform_kernel(),
        Box::new(EmpiricalLengthDist::from_pairs(vec![(3, 1.0)])),
    );
    let declared = pass.declared_choices();
    assert_eq!(declared.len(), 3);
    assert!(declared.contains(&"mutate.s5f.count".to_string()));
    assert!(declared.contains(&"mutate.s5f.site[0..n]".to_string()));
    assert!(declared.contains(&"mutate.s5f.base[0..n]".to_string()));

    assert_eq!(
        pass.declared_choice_patterns(),
        vec![
            address::ChoiceAddressPattern::MutateS5fCount,
            address::ChoiceAddressPattern::MutateS5fSite,
            address::ChoiceAddressPattern::MutateS5fBase,
        ]
    );
}
