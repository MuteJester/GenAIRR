use super::*;
use crate::duplex_repair::model::tests::fixture;
use crate::family_replay::{with_budget_scope, FamilyReplayLimits};
use crate::lineage::branching::{execute_mutation_choices, LineageMutationContext};
use crate::lineage::{
    affinity::{duplex_sim_to_aa, make_mature_target, AffinityModel},
    family,
};
use crate::passes::S5FMutationPass;
use crate::s5f::S5FKernel;

fn setup() -> (crate::refdata::RefDataConfig, Outcome, AttemptConfig) {
    let (reference, model) = fixture();
    let founder = Outcome {
        revisions: vec![model.prepare_row(0, &reference).unwrap().seed],
        pass_names: vec![],
        trace: Trace::new(),
        events: vec![],
    };
    let bit = |x| F64Bits::new(x).unwrap();
    let config = AttemptConfig {
        rate: bit(0.2),
        lambda_base: bit(3.0),
        lambda_mut: bit(0.0),
        max_generations: 3,
        n_max: 30,
        n_sample: 12,
        seed_label: "12".into(),
        selection_strength: bit(1.0),
        beta: bit(0.1),
        target_aa: None,
        mature_substitutions: 3,
        mutability: vec![bit(1.0); 1024],
        substitution: vec![bit(1.0); 4096],
        contracts: None,
        strict: false,
    };
    (reference, founder, config)
}
fn run(
    reference: &crate::refdata::RefDataConfig,
    founder: &Outcome,
    config: AttemptConfig,
    expected: Option<&FamilyAttemptReceipt>,
) -> Result<(LineageTree, Vec<Option<Outcome>>, FamilyAttemptReceipt), String> {
    with_budget_scope(&mut FamilyReplayBudget::default(), |budget| {
        family::simulate_family_receipt(founder, reference, config, expected, budget)
            .map_err(|e| e.to_string())
    })
}
#[test]
fn complete_attempt_roundtrips_targets_growth_edges_samples_and_actual_states() {
    let (reference, founder, config) = setup();
    let (tree, edges, receipt) = run(&reference, &founder, config.clone(), None).unwrap();
    let wire = serde_json::to_string(&receipt).unwrap();
    let decoded: FamilyAttemptReceipt = serde_json::from_str(&wire).unwrap();
    let (replayed, replayed_edges, reconstructed) =
        run(&reference, &founder, config, Some(&decoded)).unwrap();
    assert_eq!(tree.nodes, replayed.nodes);
    assert_eq!(
        serde_json::to_value(receipt).unwrap(),
        serde_json::to_value(reconstructed).unwrap()
    );
    for (a, b) in edges.iter().zip(&replayed_edges) {
        if let (Some(a), Some(b)) = (a, b) {
            assert_eq!(a.events, b.events);
            assert_eq!(a.trace.choices(), b.trace.choices());
            assert_eq!(
                a.final_simulation().duplex_provenance,
                b.final_simulation().duplex_provenance
            );
        }
    }
}
#[test]
fn capture_preserves_existing_generation_target_and_sampling_streams() {
    let (reference, founder, config) = setup();
    let aa = duplex_sim_to_aa(founder.final_simulation(), &reference).unwrap();
    let target = make_mature_target(
        &aa,
        config.mature_substitutions,
        &mut Rng::new(12 ^ 0x7461_7267_6574_0001),
    );
    let affinity = AffinityModel::new(
        target.clone(),
        vec![1.0; aa.len()],
        config.beta.get(),
        config.selection_strength.get(),
        &aa,
    );
    let context = LineageMutationContext {
        refdata: &reference,
        contracts: None,
        strict: false,
    };
    let mutator = S5FMutationPass::new_rate(
        S5FKernel::new(vec![1.0; 1024], vec![1.0; 4096]),
        config.rate.get(),
    );
    let (ordinary, _, _) = family::simulate_family_edges_with_context(
        founder.final_simulation(),
        &config.params().unwrap(),
        &mutator,
        Some(&affinity),
        &context,
    )
    .unwrap();
    let (captured, _, receipt) = run(&reference, &founder, config, None).unwrap();
    assert_eq!(ordinary.nodes, captured.nodes);
    assert_eq!(
        receipt.target.resolved_target,
        Some(String::from_utf8(target).unwrap())
    );
}
#[test]
fn malformed_semantic_choices_fail_after_recomputable_wire_changes() {
    let (reference, founder, config) = setup();
    let (_, _, receipt) = run(&reference, &founder, config.clone(), None).unwrap();
    let mut variants = vec![];
    let mut x = receipt.clone();
    x.schema_version = 2;
    variants.push(x);
    let mut x = receipt.clone();
    x.founder_state_digest = "wrong".into();
    variants.push(x);
    let mut x = receipt.clone();
    x.target.choices[0].position = u32::MAX;
    variants.push(x);
    let mut x = receipt.clone();
    x.growth[0].parent_id = 7;
    variants.push(x);
    let mut x = receipt.clone();
    x.growth[0].rate = F64Bits::new(0.0).unwrap();
    variants.push(x);
    let mut x = receipt.clone();
    x.growth.push(x.growth[0].clone());
    variants.push(x);
    if !receipt.edges.is_empty() {
        let mut x = receipt.clone();
        x.edges[0].parent_id = 9;
        variants.push(x);
        let mut x = receipt.clone();
        x.edges[0].seed_label = "01".into();
        variants.push(x);
        let mut x = receipt.clone();
        x.edges[0]
            .trace
            .record("unexpected", crate::trace::ChoiceValue::Int(0));
        variants.push(x);
        let mut x = receipt.clone();
        x.edges[0].events_digest = "wrong".into();
        variants.push(x);
    }
    let mut x = receipt.clone();
    x.final_live.push(0);
    variants.push(x);
    let mut x = receipt.clone();
    x.completed_generations += 1;
    variants.push(x);
    let mut x = receipt.clone();
    x.nodes[0].abundance += 1;
    variants.push(x);
    if !receipt.samples.is_empty() {
        let mut x = receipt.clone();
        x.samples[0].live_index = u32::MAX;
        variants.push(x);
    }
    for x in variants {
        assert!(run(&reference, &founder, config.clone(), Some(&x)).is_err());
    }
}
#[test]
fn attempt_budgets_include_nodes_choices_history_and_raw_poisson_work() {
    let (reference, founder, config) = setup();
    for which in 0..4 {
        let mut limits = FamilyReplayLimits::default();
        match which {
            0 => limits.max_nodes = 0,
            1 => limits.max_choices = 0,
            2 => limits.max_history_entries = 0,
            _ => limits.max_work = 1,
        };
        let mut budget = FamilyReplayBudget::new(limits);
        let result = with_budget_scope(&mut budget, |budget| {
            family::simulate_family_receipt(&founder, &reference, config.clone(), None, budget)
                .map_err(|e| e.to_string())
        });
        assert!(result.is_err());
    }
    let mut rng = Rng::new(44);
    let mut words = 0;
    let result = crate::lineage::poisson::poisson_sample_with_work(&mut rng, 64.0, || {
        words += 1;
        if words > 1 {
            Err("work fence".into())
        } else {
            Ok(())
        }
    });
    assert!(result.is_err());
    assert_eq!(rng.words_consumed(), 1);
}

fn literal_capacity_receipt(
    reference: &crate::refdata::RefDataConfig,
    founder: &Outcome,
    mut config: AttemptConfig,
) -> FamilyAttemptReceipt {
    config.rate = F64Bits::new(0.0).unwrap();
    config.selection_strength = F64Bits::new(0.0).unwrap();
    config.n_max = 4;
    config.n_sample = 3;
    config.max_generations = 2;
    let mut budget = FamilyReplayBudget::default();
    let driver = FamilyReplayDriver::new(config.clone(), founder, None, &mut budget).unwrap();
    // Build static commitments through the real zero-mutation pass. The
    // topology and sample arithmetic below are literal independent controls.
    let context = LineageMutationContext {
        refdata: reference,
        contracts: None,
        strict: false,
    };
    let mutator = S5FMutationPass::new_rate(S5FKernel::new(vec![1.0; 1024], vec![1.0; 4096]), 0.0);
    let mut trace = Trace::new();
    trace.record("mutate.s5f.count", crate::trace::ChoiceValue::Int(0));
    let edge =
        execute_mutation_choices(founder.final_simulation(), &trace, &mutator, &context).unwrap();
    let state = semantic_state_digest(founder.final_simulation()).unwrap();
    let events = outcome_events_digest(&edge).unwrap();
    let mut result = driver.receipt;
    result.growth = vec![
        ParentProposal {
            generation: 1,
            parent_id: 0,
            rate: F64Bits::new(2.25).unwrap(),
            proposed_offspring: 3,
        },
        ParentProposal {
            generation: 2,
            parent_id: 1,
            rate: F64Bits::new(0.75).unwrap(),
            proposed_offspring: 4,
        },
        ParentProposal {
            generation: 2,
            parent_id: 2,
            rate: F64Bits::new(0.75).unwrap(),
            proposed_offspring: 0,
        },
        ParentProposal {
            generation: 2,
            parent_id: 3,
            rate: F64Bits::new(0.75).unwrap(),
            proposed_offspring: 1,
        },
    ];
    result.edges = (1..=7)
        .map(|id| DivisionReceipt {
            node_id: id,
            parent_id: if id <= 3 { 0 } else { 1 },
            generation: if id <= 3 { 1 } else { 2 },
            seed_label: id.to_string(),
            trace: trace.clone(),
            parent_state_digest: state.clone(),
            child_state_digest: state.clone(),
            events_digest: events.clone(),
        })
        .collect();
    result.final_live = vec![4, 5, 6, 7];
    result.samples = vec![
        SampleChoice {
            live_index: 0,
            node_id: 4,
        },
        SampleChoice {
            live_index: 3,
            node_id: 7,
        },
        SampleChoice {
            live_index: 0,
            node_id: 4,
        },
    ];
    result.completed_generations = 2;
    result.peak_live = 4;
    result.termination = "max_generations".into();
    result.nodes = (0..=7)
        .map(|id| NodeSummary {
            node_id: id,
            parent_id: if id == 0 {
                None
            } else {
                Some(if id <= 3 { 0 } else { 1 })
            },
            generation: if id == 0 {
                0
            } else if id <= 3 {
                1
            } else {
                2
            },
            mutations_from_parent: 0,
            affinity: F64Bits::new(0.0).unwrap(),
            abundance: if id == 4 {
                2
            } else if id == 7 {
                1
            } else {
                0
            },
            observed: id == 4 || id == 7,
        })
        .collect();
    result
}
#[test]
fn literal_capacity_replay_consumes_suppressed_proposal_without_inventing_child() {
    let (reference, founder, config) = setup();
    let receipt = literal_capacity_receipt(&reference, &founder, config);
    let (tree, _, _) = run(&reference, &founder, receipt.config.clone(), Some(&receipt)).unwrap();
    assert_eq!(tree.nodes.len(), 8);
    assert_eq!(tree.nodes[4].parent_id, Some(1));
    assert_eq!(tree.nodes[4].abundance, 2);
    assert_eq!(tree.nodes[7].abundance, 1);
    let mut missing = receipt.clone();
    missing.growth.pop();
    assert!(run(&reference, &founder, missing.config.clone(), Some(&missing)).is_err());
    let mut extra = receipt.clone();
    extra.edges.push(extra.edges[0].clone());
    assert!(run(&reference, &founder, extra.config.clone(), Some(&extra)).is_err());
}
#[test]
fn extinct_attempt_replays_no_sampling_and_preserves_root() {
    let (reference, founder, mut config) = setup();
    config.lambda_base = F64Bits::new(0.0).unwrap();
    config.selection_strength = F64Bits::new(0.0).unwrap();
    let (tree, _, receipt) = run(&reference, &founder, config.clone(), None).unwrap();
    assert_eq!(tree.len(), 1);
    assert_eq!(receipt.termination, "extinct");
    assert!(receipt.samples.is_empty());
    assert!(receipt.final_live.is_empty());
    run(&reference, &founder, config, Some(&receipt)).unwrap();
}

#[test]
fn literal_target_repeated_positions_and_reversion_have_positive_support() {
    let (reference, founder, mut config) = setup();
    config.lambda_base = F64Bits::new(0.0).unwrap();
    let (_, _, mut receipt) = run(&reference, &founder, config.clone(), None).unwrap();
    // Founder translation is KCTFK; K -> A -> K -> R at the same site.
    receipt.target.choices = ["A", "K", "R"]
        .into_iter()
        .map(|aa| TargetChoice {
            position: 0,
            residue: aa.into(),
        })
        .collect();
    receipt.target.resolved_target = Some("RCTFK".into());
    // BLOSUM K,K=5 and K,R=2, so distance is 3 at beta=.1.
    receipt.nodes[0].affinity = F64Bits::new((-0.30000000000000004f64).exp()).unwrap();
    run(&reference, &founder, config.clone(), Some(&receipt)).unwrap();
    receipt.target.choices[1].residue = "A".into(); // A -> A has zero replacement support.
    assert!(run(&reference, &founder, config, Some(&receipt)).is_err());
}

#[test]
fn family_s5f_replay_checks_normalizable_support_without_product_underflow() {
    let (reference, founder, _) = setup();
    let mut trace = Trace::new();
    trace.record("mutate.s5f.count", crate::trace::ChoiceValue::Int(1));
    trace.record("mutate.s5f.site[0]", crate::trace::ChoiceValue::Int(5));
    trace.record("mutate.s5f.base[0]", crate::trace::ChoiceValue::Base(b'A'));
    for (mu, row, allowed) in [
        (f64::MAX, 1.0, false),
        (1.0, f64::MAX, false),
        (1e-200, 1e-200, true),
    ] {
        let mutator =
            S5FMutationPass::new_rate(S5FKernel::new(vec![mu; 1024], vec![row; 4096]), 0.2);
        let context = LineageMutationContext {
            refdata: &reference,
            contracts: None,
            strict: false,
        };
        let result = crate::family_replay::with_replay_scope(|| {
            execute_mutation_choices(founder.final_simulation(), &trace, &mutator, &context)
                .map_err(|e| e.to_string())
        });
        assert_eq!(result.is_ok(), allowed, "mu={mu} row={row}");
    }
}

#[test]
fn permissive_family_pass_dispatch_returns_structured_invalid_choice_errors() {
    use crate::dist::AllelePoolDist;
    use crate::pass::{Pass, PassContext};
    use crate::passes::{InvertDPass, SampleAllelePass};
    use crate::replay::TraceCursor;
    use crate::trace::ChoiceValue;

    let (reference, founder, _) = setup();
    let mut invalid_s5f = Trace::new();
    invalid_s5f.record("mutate.s5f.count", ChoiceValue::Int(1));
    invalid_s5f.record("mutate.s5f.site[0]", ChoiceValue::Int(5));
    invalid_s5f.record("mutate.s5f.base[0]", ChoiceValue::Base(b'A'));
    let mut invalid_allele = Trace::new();
    invalid_allele.record("sample_allele.v", ChoiceValue::AlleleId(99));
    let mut invalid_orientation = Trace::new();
    invalid_orientation.record("sample_allele.d.inverted", ChoiceValue::Bool(true));
    let cases: Vec<(Box<dyn Pass>, Trace)> = vec![
        (
            Box::new(S5FMutationPass::new_rate(
                S5FKernel::new(vec![f64::MAX; 1024], vec![1.0; 4096]),
                0.2,
            )),
            invalid_s5f,
        ),
        (
            Box::new(SampleAllelePass::new(
                crate::ir::Segment::V,
                Box::new(AllelePoolDist::uniform(&reference.v_pool)),
            )),
            invalid_allele,
        ),
        (Box::new(InvertDPass::new(0.0)), invalid_orientation),
    ];
    for (pass, choices) in cases {
        let mut cursor = TraceCursor::new(choices.choices());
        let mut trace = Trace::new();
        let mut rng = Rng::new(1);
        let result = crate::family_replay::with_replay_scope(|| {
            let mut context = PassContext {
                trace: &mut trace,
                rng: &mut rng,
                pass_index: 0,
                refdata: Some(&reference),
                contracts: None,
                feasibility: None,
                reference_index: None,
                replay_cursor: Some(&mut cursor),
                event_log_sink: None,
            };
            pass.execute_with_policy(founder.final_simulation(), &mut context, false)
                .map_err(|e| e.to_string())
        });
        assert!(result.is_err(), "{} accepted invalid choice", pass.name());
        assert_eq!(rng.words_consumed(), 0, "{} used fallback RNG", pass.name());
    }
}

#[test]
fn compiled_family_replay_preserves_permissive_empty_contract_noop_policy() {
    use crate::compiled::{ExecutionPolicy, OwnedCompiledSimulator};
    use crate::contract::ContractSet;
    use crate::dist::EmpiricalLengthDist;
    use crate::pass::PassPlan;
    let (reference, founder, _) = setup();
    let mut plan = PassPlan::new();
    plan.push(Box::new(S5FMutationPass::new(
        S5FKernel::new(vec![1.0; 1024], vec![0.0; 4096]),
        Box::new(EmpiricalLengthDist::from_pairs([(1, 1.0)])),
    )));
    let compiled = OwnedCompiledSimulator::compile(
        plan,
        Some(reference),
        Some(ContractSet::new()),
        ExecutionPolicy::Permissive,
    )
    .unwrap();
    let initial = founder.final_simulation().clone();
    let fresh = compiled.run_one_from(initial.clone(), 7).unwrap();
    assert_eq!(fresh.trace.len(), 1); // requested count; no admissible event.
    let replayed = crate::family_replay::with_replay_scope(|| {
        compiled
            .replay_family_from(
                initial.clone(),
                fresh.trace.choices(),
                ExecutionPolicy::Permissive,
            )
            .map_err(|e| e.to_string())
    })
    .unwrap();
    assert_eq!(fresh.events, replayed.events);
    assert_eq!(
        semantic_state_digest(fresh.final_simulation()).unwrap(),
        semantic_state_digest(replayed.final_simulation()).unwrap()
    );
    // Outside the family scope, historical source-trace dispatch stays strict.
    assert!(compiled
        .replay_family_from(initial, fresh.trace.choices(), ExecutionPolicy::Permissive)
        .is_err());
}

#[test]
fn bounded_capture_rejects_s5f_numerical_noops_and_fences_empty_support_loops() {
    use crate::contract::ContractSet;
    use crate::dist::EmpiricalLengthDist;
    use crate::pass::{Pass, PassContext};
    let (reference, founder, _) = setup();
    let empty = ContractSet::new();
    for (mu, row, count, max_work) in [
        (f64::MAX, 1.0, 1, 1000),
        (1.0, f64::MAX, 1, 1000),
        (1.0, 0.0, 1000, 60),
    ] {
        let mutator = S5FMutationPass::new(
            S5FKernel::new(vec![mu; 1024], vec![row; 4096]),
            Box::new(EmpiricalLengthDist::from_pairs([(count, 1.0)])),
        );
        let mut trace = Trace::new();
        let mut rng = Rng::new(7);
        let mut budget = FamilyReplayBudget::new(FamilyReplayLimits {
            max_work,
            ..FamilyReplayLimits::default()
        });
        let result = with_budget_scope(&mut budget, |_| {
            let mut context = PassContext {
                trace: &mut trace,
                rng: &mut rng,
                pass_index: 0,
                refdata: Some(&reference),
                contracts: Some(&empty),
                feasibility: None,
                reference_index: None,
                replay_cursor: None,
                event_log_sink: None,
            };
            mutator
                .execute_with_policy(founder.final_simulation(), &mut context, false)
                .map_err(|e| e.to_string())
        });
        assert!(result.is_err(), "mu={mu} row={row} count={count}");
        if count > 1 {
            assert!(result.unwrap_err().contains("resource limit: work"));
        }
    }
}

#[test]
fn family_replay_requires_selected_site_for_unconstrained_empty_row() {
    use crate::contract::ContractSet;
    use crate::dist::EmpiricalLengthDist;
    use crate::trace::ChoiceValue;
    let (reference, founder, _) = setup();
    let mutator = S5FMutationPass::new(
        S5FKernel::new(vec![1.0; 1024], vec![0.0; 4096]),
        Box::new(EmpiricalLengthDist::from_pairs([(1, 1.0)])),
    );
    let mut count_only = Trace::new();
    count_only.record("mutate.s5f.count", ChoiceValue::Int(1));
    let mut site_only = count_only.clone();
    site_only.record("mutate.s5f.site[0]", ChoiceValue::Int(5));
    let context = LineageMutationContext {
        refdata: &reference,
        contracts: None,
        strict: false,
    };
    for (trace, allowed) in [(&site_only, true), (&count_only, false)] {
        let result = crate::family_replay::with_replay_scope(|| {
            execute_mutation_choices(founder.final_simulation(), trace, &mutator, &context)
                .map_err(|e| e.to_string())
        });
        assert_eq!(result.is_ok(), allowed);
    }
    let empty = ContractSet::new();
    let with_contracts = LineageMutationContext {
        contracts: Some(&empty),
        ..context
    };
    crate::family_replay::with_replay_scope(|| {
        execute_mutation_choices(
            founder.final_simulation(),
            &count_only,
            &mutator,
            &with_contracts,
        )
        .map_err(|e| e.to_string())
    })
    .unwrap();
}
