use super::*;
use crate::contract::productive;
use crate::dist::FilteredSampleError;

#[test]
fn s5f_mutation_productive_filters_base_that_would_create_stop() {
    let kernel = s5f_stop_filter_kernel(true);
    let (cfg, sim) = s5f_productive_vj_fixture();
    let seed = find_seed_for_s5f_unconstrained_base(&kernel, &sim, b'A');
    let contracts = productive();

    let constrained = PassRuntime::execute_with_context(
        &s5f_single_mutation_plan(kernel.clone()),
        sim.clone(),
        seed,
        Some(&cfg),
        Some(&contracts),
    );

    assert_eq!(
        constrained.trace.find("mutate.s5f.site[0]").unwrap().value,
        ChoiceValue::Int(2)
    );
    assert_eq!(
        constrained.trace.find("mutate.s5f.base[0]").unwrap().value,
        ChoiceValue::Base(b'C')
    );
    assert!(contracts
        .verify(constrained.final_simulation(), Some(&cfg))
        .is_ok());

    let unconstrained = PassRuntime::execute_with_context(
        &s5f_single_mutation_plan(kernel),
        sim,
        seed,
        Some(&cfg),
        None,
    );
    assert_eq!(
        unconstrained
            .trace
            .find("mutate.s5f.base[0]")
            .unwrap()
            .value,
        ChoiceValue::Base(b'A')
    );
    let violations = productive()
        .verify(unconstrained.final_simulation(), Some(&cfg))
        .unwrap_err();
    assert!(violations
        .iter()
        .any(|v| v.contract_name == "no_stop_codon_in_junction"));
}

#[test]
fn s5f_mutation_strict_errors_when_base_filter_empty() {
    let (cfg, sim) = s5f_productive_vj_fixture();
    let contracts = productive();

    let err = PassRuntime::execute_strict_with_context(
        &s5f_single_mutation_plan(s5f_stop_filter_kernel(false)),
        sim,
        0,
        Some(&cfg),
        Some(&contracts),
    )
    .unwrap_err();

    assert_eq!(err.pass_name(), "mutate.s5f");
    assert_eq!(err.address(), "mutate.s5f.base[0]");
    assert_eq!(
        err.constraint_reason(),
        Some(FilteredSampleError::EmptyAdmissibleSupport)
    );
}
