use super::*;

#[test]
fn indel_pass_productive_filters_insertion_to_structurally_safe_site() {
    let (cfg, sim) = make_substitution_productive_vj_fixture();
    let contracts = productive();

    let mut plan = PassPlan::new();
    plan.push(Box::new(IndelPass::new(
        fixed_count(1),
        1.0,
        Box::new(UniformBase),
    )));

    let outcome = PassRuntime::execute_with_context(&plan, sim, 0, Some(&cfg), Some(&contracts));

    assert_eq!(
        outcome.trace.find("corrupt.indel.kind[0]").unwrap().value,
        ChoiceValue::Bool(true)
    );
    assert_eq!(
        outcome.trace.find("corrupt.indel.site[0]").unwrap().value,
        ChoiceValue::Int(6)
    );
    assert!(contracts
        .verify(outcome.final_simulation(), Some(&cfg))
        .is_ok());
}

#[test]
fn indel_pass_permissive_falls_back_when_deletions_cannot_satisfy_contracts() {
    let (cfg, sim) = make_substitution_productive_vj_fixture();
    let contracts = productive();

    let mut plan = PassPlan::new();
    plan.push(Box::new(IndelPass::new(
        fixed_count(1),
        0.0,
        Box::new(UniformBase),
    )));

    let outcome = PassRuntime::execute_with_context(&plan, sim, 0, Some(&cfg), Some(&contracts));

    assert_eq!(
        outcome.trace.find("corrupt.indel.kind[0]").unwrap().value,
        ChoiceValue::Bool(false)
    );
    let violations = contracts
        .verify(outcome.final_simulation(), Some(&cfg))
        .unwrap_err();
    assert!(
        violations.iter().any(|v| {
            v.contract_name == "productive_junction_frame"
                || v.contract_name == "anchor_preserved.v"
                || v.contract_name == "anchor_preserved.j"
        }),
        "expected productive-bundle structural violation, got {:?}",
        violations
    );
}

#[test]
fn indel_pass_strict_errors_when_deletion_filter_empty() {
    let (cfg, sim) = make_substitution_productive_vj_fixture();
    let contracts = productive();

    let mut plan = PassPlan::new();
    plan.push(Box::new(IndelPass::new(
        fixed_count(1),
        0.0,
        Box::new(UniformBase),
    )));

    let err = PassRuntime::execute_strict_with_context(&plan, sim, 0, Some(&cfg), Some(&contracts))
        .unwrap_err();

    assert_eq!(err.pass_name(), "corrupt.indel");
    assert_eq!(err.address(), "corrupt.indel.site[0]");
    assert_eq!(
        err.constraint_reason(),
        Some(FilteredSampleError::EmptyAdmissibleSupport)
    );
}
