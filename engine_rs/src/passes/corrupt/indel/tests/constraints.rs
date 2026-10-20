use super::*;
use crate::contract::ContractSet;
use crate::ir::flag;
use crate::refdata::RefDataConfig;

struct EnumeratedEvent {
    insertion: bool,
    delta: i32,
    natural_weight: f64,
    post: Simulation,
}

/// Independently enumerate the small fixture's uniform natural event law.
/// Do not call the sampler's candidate enumeration: these tests must
/// detect errors in its probability weights as well as in its support.
fn enumerate_events(
    sim: &Simulation,
    _cfg: &RefDataConfig,
    _contracts: &ContractSet,
) -> Vec<EnumeratedEvent> {
    let len = sim.pool.len() as u32;
    assert!(len > 0);
    let mut events = Vec::new();
    for site in 0..=len {
        for insertion in [true, false] {
            if !insertion && site == len {
                continue;
            }
            // Deleting original anchor provenance is irreversible. Inserting
            // inside an anchor may be undone by the second event.
            if !insertion {
                let n = &sim.pool.as_slice()[site as usize];
                if matches!(n.segment, Segment::V | Segment::J)
                    && n.germline_pos.get().is_some_and(|g| g < 3)
                {
                    continue;
                }
            }
            let anchor = |segment| {
                let region = sim
                    .sequence
                    .regions
                    .iter()
                    .find(|r| r.segment == segment)
                    .unwrap();
                crate::coding::anchor_position(sim, region, 0).unwrap()
            };
            let (v, j) = (anchor(Segment::V), anchor(Segment::J));
            let delta: i8 = if insertion {
                i8::from(v < site && site <= j)
            } else {
                -i8::from(v <= site && site < j)
            };
            if insertion {
                for base in b"ACGT" {
                    let nuc = Nucleotide::synthetic(
                        *base,
                        IndelPass::insertion_segment(sim, site),
                        flag::INDEL_INSERTED,
                    );
                    events.push(EnumeratedEvent {
                        insertion,
                        delta: i32::from(delta),
                        natural_weight: 0.5 / f64::from(len + 1) / 4.0,
                        post: sim.with_indel_inserted(site, nuc),
                    });
                }
            } else {
                events.push(EnumeratedEvent {
                    insertion,
                    delta: i32::from(delta),
                    natural_weight: 0.5 / f64::from(len),
                    post: sim.with_indel_deleted(site),
                });
            }
        }
    }
    events
}

/// Exhaust the two-event proposal tree and condition on full-bundle acceptance.
/// Every history receives the product of its natural event probabilities,
/// evaluated on the actual intermediate molecule.
fn exact_two_event_first_delta(
    sim: &Simulation,
    cfg: &RefDataConfig,
    contracts: &ContractSet,
) -> [f64; 3] {
    let first_events = enumerate_events(sim, cfg, contracts);
    let mut accepted = [0.0; 3];
    for first in &first_events {
        let second_events: Vec<_> = enumerate_events(&first.post, cfg, contracts)
            .into_iter()
            .collect();
        for second in second_events {
            if contracts.verify(&second.post, Some(cfg)).is_ok() {
                accepted[first.delta.rem_euclid(3) as usize] +=
                    first.natural_weight * second.natural_weight;
            }
        }
    }
    let total: f64 = accepted.iter().sum();
    assert!(total > 0.0);
    accepted.map(|weight| weight / total)
}

#[test]
fn interleaved_independent_probes_match_complete_three_event_oracle() {
    let (cfg, sim) = balanced_dp_distribution_fixture();
    let contracts = productive();
    let mut masses = [0.0; 3];
    for first in enumerate_events(&sim, &cfg, &contracts) {
        for second in enumerate_events(&first.post, &cfg, &contracts) {
            for third in enumerate_events(&second.post, &cfg, &contracts) {
                if contracts.verify(&third.post, Some(&cfg)).is_ok() {
                    masses[first.delta.rem_euclid(3) as usize] +=
                        first.natural_weight * second.natural_weight * third.natural_weight;
                }
            }
        }
    }
    let total: f64 = masses.iter().sum();
    assert!(total > 0.0);
    let expected = masses.map(|m| m / total);
    let mut counts = [0u32; 3];
    let mut probes = 0;
    let mut plan = PassPlan::new();
    plan.push(Box::new(
        IndelPass::new(fixed_count(3), 0.5, Box::new(UniformBase)).with_sampling_options(
            IndelSamplingOptions {
                rejection_attempts: 1,
                exact_node_budget: 100_000,
                ..Default::default()
            },
        ),
    ));
    const N: u32 = 1024;
    for seed in 0..N {
        let out = PassRuntime::execute_strict_with_context(
            &plan,
            sim.clone(),
            u64::from(seed),
            Some(&cfg),
            Some(&contracts),
        )
        .unwrap();
        let ChoiceValue::Bool(insertion) = out.trace.find("corrupt.indel.kind[0]").unwrap().value
        else {
            panic!()
        };
        let ChoiceValue::Int(site) = out.trace.find("corrupt.indel.site[0]").unwrap().value else {
            panic!()
        };
        let delta = if insertion {
            i32::from(0 < site && site <= 6)
        } else {
            -i32::from(site < 6)
        };
        counts[delta.rem_euclid(3) as usize] += 1;
        let report = &out.trace.diagnostics()[0];
        if report.strategy == "exact_search_independent_probe" {
            probes += 1;
            assert!(report.explored_states >= 19);
            assert!(report.attempts > 1);
        }
    }
    assert!(
        probes > N / 2,
        "the test must exercise independent restart acceptance"
    );
    for (count, probability) in counts.into_iter().zip(expected) {
        let tolerance = 5.0 * (f64::from(N) * probability * (1.0 - probability)).sqrt();
        assert!((f64::from(count) - f64::from(N) * probability).abs() < tolerance);
    }
}

#[test]
fn indel_pass_productive_filters_insertion_to_structurally_safe_site() {
    // With V=TGT at [0..3) and J=TGG at [3..6), insertion
    // inside an anchor codon is forbidden. Insertion immediately
    // before either intact anchor can pass anchor preservation,
    // but fails the complete coding frame or junction frame.
    // Only post-pool insertion survives the full productive bundle:
    // {site=6} × 4 bases.
    let (cfg, sim) = make_substitution_productive_vj_fixture();
    let contracts = productive();

    let mut plan = PassPlan::new();
    plan.push(Box::new(IndelPass::new(
        fixed_count(1),
        1.0,
        Box::new(UniformBase),
    )));

    let outcome = PassRuntime::execute_with_context(&plan, sim, 0, Some(&cfg), Some(&contracts));

    // An insertion was committed (kind=true) at the only
    // bundle-admissible site for this fixture.
    assert_eq!(
        outcome.trace.find("corrupt.indel.kind[0]").unwrap().value,
        ChoiceValue::Bool(true)
    );
    assert_eq!(
        outcome.trace.find("corrupt.indel.site[0]").unwrap().value,
        ChoiceValue::Int(6)
    );
    // The complete tuple preserves the coding frame and introduces no stop.
    assert!(contracts
        .verify(outcome.final_simulation(), Some(&cfg))
        .is_ok());
}

#[test]
fn indel_pass_permissive_reduces_to_noop_when_deletions_cannot_satisfy_contracts() {
    // When every (kind, site) candidate is rejected by the active
    // contract bundle, permissive mode honors the contract by
    // emitting `IndelEvent::NoOp` instead of falling through to
    // unconstrained sampling. The trace still records the slot
    // (kind=false, site=-1) so consumers can audit the reduction,
    // and the final simulation satisfies the contract bundle.
    let (cfg, sim) = make_substitution_productive_vj_fixture();
    let contracts = productive();

    let mut plan = PassPlan::new();
    plan.push(Box::new(IndelPass::new(
        fixed_count(1),
        0.0,
        Box::new(UniformBase),
    )));

    let outcome = PassRuntime::execute_with_context(&plan, sim, 0, Some(&cfg), Some(&contracts));

    // The slot was consumed as a NoOp — kind=false (matches the
    // deletion sentinel) and site=-1 (the "no-site" marker).
    assert_eq!(
        outcome.trace.find("corrupt.indel.kind[0]").unwrap().value,
        ChoiceValue::Bool(false)
    );
    assert_eq!(
        outcome.trace.find("corrupt.indel.site[0]").unwrap().value,
        ChoiceValue::Int(-1)
    );
    // The contract bundle holds after permissive reduce-and-skip.
    assert!(
        contracts
            .verify(outcome.final_simulation(), Some(&cfg))
            .is_ok(),
        "permissive NoOp must leave the final simulation contract-admissible"
    );
}

#[test]
fn indel_pass_tuple_sampler_balances_frame_under_productive_count_2() {
    // Two edits can compensate each other; all accepted complete tuples
    // must satisfy the same final productive bundle.
    let (cfg, sim) = make_productive_indel_balance_fixture();
    let contracts = productive();

    // Run the pass under productive() at a handful of seeds to
    // pin the architectural property: every committed tuple leaves the
    // complete bundle satisfied.
    for seed in 0..16u64 {
        let mut plan = PassPlan::new();
        plan.push(Box::new(IndelPass::new(
            fixed_count(2),
            0.5,
            Box::new(UniformBase),
        )));
        let outcome = PassRuntime::execute_with_context(
            &plan,
            sim.clone(),
            seed,
            Some(&cfg),
            Some(&contracts),
        );
        assert!(
            contracts
                .verify(outcome.final_simulation(), Some(&cfg))
                .is_ok(),
            "seed {} produced a contract violation",
            seed
        );
    }
}

#[test]
fn indel_pass_tuple_sampler_balances_frame_under_productive_count_3() {
    // Three edits can compensate their coding-frame effects. Intermediate
    // frames are not the acceptance criterion; the complete state is.
    let (cfg, sim) = make_productive_indel_balance_fixture();
    let contracts = productive();
    for seed in 0..16u64 {
        let mut plan = PassPlan::new();
        plan.push(Box::new(IndelPass::new(
            fixed_count(3),
            0.5,
            Box::new(UniformBase),
        )));
        let outcome = PassRuntime::execute_with_context(
            &plan,
            sim.clone(),
            seed,
            Some(&cfg),
            Some(&contracts),
        );
        assert!(
            contracts
                .verify(outcome.final_simulation(), Some(&cfg))
                .is_ok(),
            "count=3 seed {} produced a contract violation",
            seed
        );
    }
}

/// Compare the single-event law against all frame-neutral proposals that
/// satisfy the complete productive bundle. J-tail events can be junction-
/// neutral while violating the whole coding frame, so class counts alone
/// are insufficient. Condition empirical counts on an event being committed:
/// complete proof of empty support permits a NoOp, not an invalid deletion.
#[test]
fn indel_pass_count_1_under_productive_matches_exact_enumeration() {
    let (cfg, sim) = balanced_dp_distribution_fixture();
    let contracts = productive();
    let accepted_events: Vec<_> = enumerate_events(&sim, &cfg, &contracts)
        .into_iter()
        .filter(|e| e.delta == 0 && contracts.verify(&e.post, Some(&cfg)).is_ok())
        .collect();
    let accepted_mass: f64 = accepted_events.iter().map(|e| e.natural_weight).sum();
    let insertion_mass: f64 = accepted_events
        .iter()
        .filter(|e| e.insertion)
        .map(|e| e.natural_weight)
        .sum();
    let expected_p = insertion_mass / accepted_mass;
    // A deletion of the final J base also changes the coding endpoint
    // phase. Only insertion beyond the coding region survives count=1.
    assert_eq!(expected_p, 1.0);

    const N: u32 = 4000;
    let mut insertions = 0u32;
    let mut deletions = 0u32;
    let mut noops = 0u32;
    for seed in 0..N as u64 {
        let mut plan = PassPlan::new();
        plan.push(Box::new(IndelPass::new(
            fixed_count(1),
            0.5,
            Box::new(UniformBase),
        )));
        let outcome = PassRuntime::execute_with_context(
            &plan,
            sim.clone(),
            seed,
            Some(&cfg),
            Some(&contracts),
        );
        assert!(
            contracts
                .verify(outcome.final_simulation(), Some(&cfg))
                .is_ok(),
            "seed {} violated contracts",
            seed
        );
        if outcome.trace.find("corrupt.indel.site[0]").unwrap().value == ChoiceValue::Int(-1) {
            noops += 1;
            continue;
        }
        match outcome.trace.find("corrupt.indel.kind[0]").unwrap().value {
            ChoiceValue::Bool(true) => insertions += 1,
            ChoiceValue::Bool(false) => deletions += 1,
            _ => panic!("unexpected kind variant"),
        }
    }

    let total = insertions + deletions;
    assert_eq!(total + noops, N);
    assert!(
        total > 3000,
        "too few committed events: {total}, NoOps: {noops}"
    );
    // Degenerate support has zero binomial variance: any sampled deletion
    // is a violation, rather than a statistical fluctuation.
    assert_eq!(deletions, 0, "accepted={total}, NoOps={noops}");
}

/// Distribution-invariant test for count=2: every accepted tuple
/// must have net frame delta ≡ 0 mod 3 (the coding-frame
/// guarantee). Verified by sampling many tuples and asserting
/// the post-tuple junction length is divisible by 3 against the
/// pre-tuple length.
#[test]
fn indel_pass_count_2_under_productive_preserves_frame_invariant() {
    use crate::junction::compute_junction;

    let (cfg, sim) = balanced_dp_distribution_fixture();
    let contracts = productive();
    let pre_junction = compute_junction(&sim, &cfg).expect("junction defined");
    assert!(pre_junction.is_in_frame());

    for seed in 0..256u64 {
        let mut plan = PassPlan::new();
        plan.push(Box::new(IndelPass::new(
            fixed_count(2),
            0.5,
            Box::new(UniformBase),
        )));
        let outcome = PassRuntime::execute_with_context(
            &plan,
            sim.clone(),
            seed,
            Some(&cfg),
            Some(&contracts),
        );
        let post = outcome.final_simulation();
        // This fixture has nonempty support, so each requested tuple must
        // commit and preserve the junction reading frame.
        match compute_junction(post, &cfg) {
            Some(j) => assert!(
                j.is_in_frame(),
                "seed {} produced out-of-frame junction (length {})",
                seed,
                j.length
            ),
            None => panic!("seed {} lost the junction window", seed),
        }
        assert!(
            contracts.verify(post, Some(&cfg)).is_ok(),
            "seed {} violated bundle",
            seed
        );
    }
}

/// Check natural-history weights against exact two-event enumeration, including
/// the full post-state contract check. The accepted marginal differs from the
/// first-proposal marginal because coding stops and J-tail frameshifts reject
/// some otherwise junction-frame-balanced tuples.
#[test]
fn indel_pass_count_2_matches_exact_natural_sequential_weighting() {
    let (cfg, sim) = balanced_dp_distribution_fixture();
    let contracts = productive();
    let expected = exact_two_event_first_delta(&sim, &cfg, &contracts);

    const N: u32 = 4000;
    let mut delta_plus = 0u32;
    let mut delta_zero = 0u32;
    let mut delta_minus = 0u32;
    let mut noops = 0u32;
    for seed in 0..N as u64 {
        let mut plan = PassPlan::new();
        plan.push(Box::new(IndelPass::new(
            fixed_count(2),
            0.5,
            Box::new(UniformBase),
        )));
        let outcome = PassRuntime::execute_with_context(
            &plan,
            sim.clone(),
            seed,
            Some(&cfg),
            Some(&contracts),
        );
        assert!(
            contracts
                .verify(outcome.final_simulation(), Some(&cfg))
                .is_ok(),
            "seed {seed} violated contracts"
        );

        // Extract step-0 (kind, site) from the trace. Classify
        // delta against the *initial* sim — that's the working
        // sim seen by step 0.
        let kind_0 = match outcome.trace.find("corrupt.indel.kind[0]").unwrap().value {
            ChoiceValue::Bool(b) => b,
            _ => panic!("unexpected kind variant"),
        };
        let site_0 = match outcome.trace.find("corrupt.indel.site[0]").unwrap().value {
            ChoiceValue::Int(s) => s,
            _ => panic!("unexpected site variant"),
        };
        // Per the IndelEvent trace convention: kind=true is
        // insertion, kind=false is deletion. NoOp records site=−1,
        // which means permissive empty-support — skip it.
        if site_0 == -1 {
            noops += 1;
            continue;
        }
        let delta = if kind_0 {
            i32::from(0 < site_0 && site_0 <= 6)
        } else {
            -i32::from(site_0 < 6)
        };
        match delta {
            0 => delta_zero += 1,
            1 => delta_plus += 1,
            -1 => delta_minus += 1,
            _ => unreachable!(),
        }
    }
    let total = delta_plus + delta_zero + delta_minus;
    assert_eq!(total + noops, N);
    assert!(
        total > 3000,
        "too few committed tuples: {total}, NoOps: {noops}"
    );
    for (residue, count) in [delta_zero, delta_plus, delta_minus]
        .into_iter()
        .enumerate()
    {
        let expected_p = expected[residue];
        let observed_p = f64::from(count) / f64::from(total);
        let tolerance = 5.0 * (expected_p * (1.0 - expected_p) / f64::from(total)).sqrt();
        assert!(
            (observed_p - expected_p).abs() < tolerance,
            "p(δ mod 3={residue}) = {observed_p:.4}, expected {expected_p:.4} ± \
             {tolerance:.4} (5σ); counts: + {delta_plus} 0 {delta_zero} − \
             {delta_minus}; NoOps={noops}"
        );
    }
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
    assert_eq!(err.address(), "corrupt.indel.count");
    assert_eq!(
        err.constraint_reason(),
        Some(FilteredSampleError::EmptyAdmissibleSupport)
    );
}
