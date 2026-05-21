use super::*;

#[test]
#[should_panic(expected = "insertion_prob must be in [0.0, 1.0]")]
fn indel_pass_rejects_negative_insertion_prob() {
    let _ = IndelPass::new(fixed_count(0), -0.5, Box::new(UniformBase));
}

#[test]
#[should_panic(expected = "insertion_prob must be in [0.0, 1.0]")]
fn indel_pass_rejects_insertion_prob_above_one() {
    let _ = IndelPass::new(fixed_count(0), 1.5, Box::new(UniformBase));
}

#[test]
#[should_panic(expected = "insertion_prob must be in [0.0, 1.0]")]
fn indel_pass_rejects_nan_insertion_prob() {
    let _ = IndelPass::new(fixed_count(0), f64::NAN, Box::new(UniformBase));
}

#[test]
fn indel_pass_zero_count_is_noop() {
    let mut plan = PassPlan::new();
    plan.push(Box::new(IndelPass::new(
        fixed_count(0),
        0.5,
        Box::new(UniformBase),
    )));
    let outcome = PassRuntime::execute(&plan, indel_test_sim(), 0);

    assert_eq!(outcome.trace.len(), 1);
    assert_eq!(
        outcome.trace.find("corrupt.indel.count").unwrap().value,
        ChoiceValue::Int(0)
    );
    assert_eq!(outcome.final_simulation().pool.len(), 12);
}

#[test]
#[should_panic(expected = "count distribution returned negative")]
fn indel_pass_negative_count_panics() {
    use crate::dist::UniformInt;
    let mut plan = PassPlan::new();
    plan.push(Box::new(IndelPass::new(
        Box::new(UniformInt::new(-3, -2)),
        0.5,
        Box::new(UniformBase),
    )));
    let _ = PassRuntime::execute(&plan, indel_test_sim(), 0);
}

#[test]
fn indel_pass_strict_errors_on_negative_count() {
    use crate::dist::UniformInt;
    let mut plan = PassPlan::new();
    plan.push(Box::new(IndelPass::new(
        Box::new(UniformInt::new(-3, -2)),
        0.5,
        Box::new(UniformBase),
    )));

    let err = PassRuntime::execute_strict_with_context(&plan, indel_test_sim(), 0, None, None)
        .unwrap_err();

    assert_eq!(err.pass_name(), "corrupt.indel");
    assert_eq!(err.address(), "corrupt.indel.count");
    assert!(matches!(
        err,
        PassError::InvalidDistributionOutput { value: -3, .. }
    ));
}

#[test]
fn indel_pass_insertion_prob_one_grows_pool() {
    let mut plan = PassPlan::new();
    plan.push(Box::new(IndelPass::new(
        fixed_count(4),
        1.0,
        Box::new(UniformBase),
    )));
    let outcome = PassRuntime::execute(&plan, indel_test_sim(), 7);
    assert_eq!(outcome.final_simulation().pool.len(), 12 + 4);

    for i in 0..4 {
        assert_eq!(
            outcome
                .trace
                .find(&format!("corrupt.indel.kind[{}]", i))
                .unwrap()
                .value,
            ChoiceValue::Bool(true)
        );
    }
}

#[test]
fn indel_pass_insertion_prob_zero_shrinks_pool() {
    let mut plan = PassPlan::new();
    plan.push(Box::new(IndelPass::new(
        fixed_count(4),
        0.0,
        Box::new(UniformBase),
    )));
    let outcome = PassRuntime::execute(&plan, indel_test_sim(), 11);
    assert_eq!(outcome.final_simulation().pool.len(), 12 - 4);

    for i in 0..4 {
        assert_eq!(
            outcome
                .trace
                .find(&format!("corrupt.indel.kind[{}]", i))
                .unwrap()
                .value,
            ChoiceValue::Bool(false)
        );
    }
}

#[test]
fn indel_pass_records_canonical_addresses() {
    let mut plan = PassPlan::new();
    plan.push(Box::new(IndelPass::new(
        fixed_count(3),
        1.0,
        Box::new(UniformBase),
    )));
    let outcome = PassRuntime::execute(&plan, indel_test_sim(), 5);
    assert_eq!(outcome.trace.len(), 10);
    for i in 0..3 {
        assert!(outcome
            .trace
            .find(&format!("corrupt.indel.kind[{}]", i))
            .is_some());
        assert!(outcome
            .trace
            .find(&format!("corrupt.indel.site[{}]", i))
            .is_some());
        assert!(outcome
            .trace
            .find(&format!("corrupt.indel.base[{}]", i))
            .is_some());
    }
}

#[test]
fn indel_pass_deletion_path_does_not_record_base() {
    let mut plan = PassPlan::new();
    plan.push(Box::new(IndelPass::new(
        fixed_count(3),
        0.0,
        Box::new(UniformBase),
    )));
    let outcome = PassRuntime::execute(&plan, indel_test_sim(), 5);
    assert_eq!(outcome.trace.len(), 7);
    for i in 0..3 {
        assert!(outcome
            .trace
            .find(&format!("corrupt.indel.base[{}]", i))
            .is_none());
    }
}

#[test]
fn indel_pass_empty_pool_deletion_records_minus_one_site() {
    let mut plan = PassPlan::new();
    plan.push(Box::new(IndelPass::new(
        fixed_count(2),
        0.0,
        Box::new(UniformBase),
    )));
    let outcome = PassRuntime::execute(&plan, Simulation::new(), 0);

    assert_eq!(outcome.final_simulation().pool.len(), 0);
    for i in 0..2 {
        assert_eq!(
            outcome
                .trace
                .find(&format!("corrupt.indel.site[{}]", i))
                .unwrap()
                .value,
            ChoiceValue::Int(-1)
        );
    }
}

#[test]
fn indel_pass_half_probability_produces_mixed_kinds() {
    let plan = || {
        let mut p = PassPlan::new();
        p.push(Box::new(IndelPass::new(
            fixed_count(4),
            0.5,
            Box::new(UniformBase),
        )));
        p
    };

    let mut saw_insertion = false;
    let mut saw_deletion = false;
    for seed in 0..50u64 {
        let outcome = PassRuntime::execute(&plan(), indel_test_sim(), seed);
        for i in 0..4 {
            match outcome
                .trace
                .find(&format!("corrupt.indel.kind[{}]", i))
                .unwrap()
                .value
            {
                ChoiceValue::Bool(true) => saw_insertion = true,
                ChoiceValue::Bool(false) => saw_deletion = true,
                _ => unreachable!(),
            }
        }
        if saw_insertion && saw_deletion {
            break;
        }
    }
    assert!(saw_insertion, "expected at least one insertion");
    assert!(saw_deletion, "expected at least one deletion");
}

#[test]
fn indel_pass_declared_choices() {
    let pass = IndelPass::new(fixed_count(0), 0.5, Box::new(UniformBase));
    let declared = pass.declared_choices();
    assert!(declared.contains(&"corrupt.indel.count".to_string()));
    assert!(declared.contains(&"corrupt.indel.kind[0..n]".to_string()));
    assert!(declared.contains(&"corrupt.indel.site[0..n]".to_string()));
    assert!(declared.contains(&"corrupt.indel.base[0..n]".to_string()));
}
