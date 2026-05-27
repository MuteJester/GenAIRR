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
fn s5f_mutation_permissive_empty_support_skips_slot_without_extra_rng_or_fallback() {
    // v3.0 architectural rule: under active contracts, the pass must
    // NEVER propose an action outside the contract-admissible support.
    // This regression test pins the rule in place: with an empty
    // filtered support in permissive mode, S5F must skip the slot
    // (no `site[i]`/`base[i]` trace records, no pool mutation, no
    // RNG consumed past the candidate filter), NOT fall through to
    // the unconstrained kernel-row sampler.
    let (cfg, sim) = s5f_productive_vj_fixture();
    // `s5f_stop_filter_kernel(false)` targets the TACAA context at
    // site 2; the only candidate base is `A` which creates a TAA
    // stop. Under `productive()` the mask zeroes site 2's `A` bit,
    // leaving the filtered support empty.
    let contracts = crate::contract::productive();
    let outcome = PassRuntime::execute_with_context(
        &s5f_single_mutation_plan(s5f_stop_filter_kernel(false)),
        sim.clone(),
        0,
        Some(&cfg),
        Some(&contracts),
    );

    // No per-event trace entries: the slot was skipped, not
    // re-proposed via fallback.
    assert!(outcome.trace.find("mutate.s5f.site[0]").is_none());
    assert!(outcome.trace.find("mutate.s5f.base[0]").is_none());
    // Mutation count sidecar reflects realized (zero) mutations,
    // not the requested count.
    assert_eq!(outcome.final_simulation().mutation_count, 0);
    // Pool byte-identical to input.
    for i in 0..sim.pool.len() {
        let h = NucHandle::new(i as u32);
        assert_eq!(
            outcome.final_simulation().pool.get(h).unwrap().base,
            sim.pool.get(h).unwrap().base,
            "pool byte at {} changed despite empty admissible support",
            i
        );
    }
    // Bundle still verifies — no contract was violated.
    assert!(contracts
        .verify(outcome.final_simulation(), Some(&cfg))
        .is_ok());
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
