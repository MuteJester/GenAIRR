use super::*;

fn d_fixture(lengths: &[usize], five: i64, three: i64) -> (RefDataConfig, PassPlan) {
    let mut cfg = RefDataConfig::empty(ChainType::Vdj);
    for (index, length) in lengths.iter().enumerate() {
        let _ = cfg.d_pool.push(Allele {
            name: format!("D{index}*01"),
            gene: format!("D{index}"),
            seq: vec![b'C'; *length],
            segment: Segment::D,
            anchor: None,
            functional_status: None,
            subregions: vec![],
        });
    }
    let mut plan = PassPlan::new();
    plan.push(Box::new(SampleAllelePass::new(
        Segment::D,
        Box::new(AllelePoolDist::uniform(&cfg.d_pool)),
    )));
    for (end, value) in [(TrimEnd::Five, five), (TrimEnd::Three, three)] {
        plan.push(Box::new(TrimPass::new(
            Segment::D,
            end,
            Box::new(EmpiricalLengthDist::from_pairs(vec![(value, 1.0)])),
        )));
    }
    plan.push(Box::new(AssembleSegmentPass::new(Segment::D)));
    (cfg, plan)
}

#[test]
fn impossible_combined_trim_is_a_compile_error_under_both_policies() {
    let (cfg, plan) = d_fixture(&[3], 2, 2);
    for policy in [ExecutionPolicy::Strict, ExecutionPolicy::Permissive] {
        let error = match compile_test_fixture(&plan, Some(&cfg), None, policy) {
            Err(error) => error,
            Ok(_) => panic!("combined trim cannot fit"),
        };
        assert!(error.errors.iter().any(|error| matches!(&error.kind,
            CompileErrorKind::InvalidParameterSupport { reason, .. }
                if reason.contains("trim_exceeds_allele_length")
        )));
    }
}

#[test]
fn realized_overtrim_is_a_structured_error_without_clamping_in_permissive_mode() {
    let (cfg, plan) = d_fixture(&[3, 12], 2, 2);
    for policy in [ExecutionPolicy::Strict, ExecutionPolicy::Permissive] {
        let compiled = compile_test_fixture(&plan, Some(&cfg), None, policy).unwrap();
        let mut saw_valid = false;
        let mut saw_invalid = false;
        for seed in 0..24 {
            match compiled.execute_transactional(Simulation::new(), seed) {
                Ok(outcome) => {
                    saw_valid = true;
                    assert_eq!(outcome.final_simulation().pool.len(), 8);
                }
                Err(abort) => {
                    saw_invalid = true;
                    assert!(
                        matches!(&abort.error, PassError::InvalidPlanState { reason, .. }
                        if reason == "trim_exceeds_allele_length")
                    );
                    let d = abort
                        .committed
                        .final_simulation()
                        .assignments
                        .get(Segment::D)
                        .unwrap();
                    assert_eq!((d.trim_5, d.trim_3), (2, 2));
                }
            }
        }
        assert!(saw_valid && saw_invalid);
    }
}

#[test]
fn trimming_exactly_the_complete_d_allele_is_allowed() {
    let (cfg, plan) = d_fixture(&[3], 1, 2);
    let compiled =
        compile_test_fixture(&plan, Some(&cfg), None, ExecutionPolicy::Permissive).unwrap();
    let outcome = compiled.run_one(0).unwrap();
    assert_eq!(outcome.final_simulation().pool.len(), 0);
    let d = outcome
        .final_simulation()
        .assignments
        .get(Segment::D)
        .unwrap();
    assert_eq!((d.trim_5, d.trim_3), (1, 2));
}
