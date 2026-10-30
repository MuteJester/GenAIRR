use super::branching::{replay_mutation_edge, LineageMutationContext};
use super::*;
use crate::contract::{AnchorPreserved, ContractSet};
use crate::duplex_repair::model::tests::fixture;
use crate::duplex_repair::observation::project_duplex_simulation;
use crate::ir::{DuplexEditAction, Segment, Simulation};
use crate::pass::Outcome;
use crate::passes::S5FMutationPass;
use crate::rng::Rng;
use crate::s5f::S5FKernel;

fn kernel() -> S5FKernel {
    // Context-neutral equal destination probabilities. No measured mutational
    // preferences are asserted -- and note this kernel admits the centre base
    // as a destination, which none of the four bundled kernels does: measured
    // over 960 records at 15 attempts each, every bundled kernel realised
    // every attempt, so no S5F write in a publicly reachable configuration is
    // a no-op. That is a property of the shipped tables rather than of the
    // sampler, which takes every positive row entry -- which is why this
    // test's own kernel can produce one.
    S5FKernel::new(vec![1.0; 1024], vec![1.0; 4096])
}

fn params() -> BranchingParams {
    BranchingParams {
        lambda_base: 3.0,
        lambda_mut: 0.0,
        max_generations: 3,
        n_max: 40,
        n_sample: 30,
        seed: 12,
    }
}

fn grown(
    rate: f64,
    strict: bool,
) -> (
    crate::refdata::RefDataConfig,
    Simulation,
    LineageTree,
    Vec<Simulation>,
    Vec<Option<Outcome>>,
) {
    let (refdata, model) = fixture();
    let founder = model.prepare_row(0, &refdata).unwrap().seed;
    let productive = crate::contract::productive();
    let context = LineageMutationContext {
        refdata: &refdata,
        contracts: strict.then_some(&productive),
        strict,
    };
    let (tree, sims, edges) = super::family::simulate_family_edges_with_context(
        &founder,
        &params(),
        &S5FMutationPass::new_rate(kernel(), rate),
        None,
        &context,
    )
    .unwrap();
    (refdata, founder, tree, sims, edges)
}

#[test]
fn duplex_divisions_preserve_shared_founder_and_independent_actual_event_history() {
    let (refdata, founder, tree, sims, edges) = grown(0.2, false);
    assert!(tree.len() > 3);
    let repair = &founder.duplex_provenance.as_ref().unwrap().founder;
    for node in tree.nodes.iter().skip(1) {
        let sim = &sims[node.id as usize];
        let p = sim.duplex_provenance.as_ref().unwrap();
        assert!(std::sync::Arc::ptr_eq(repair, &p.founder));
        let parent = node.parent_id.unwrap() as usize;
        let before = sims[parent].duplex_provenance.as_ref().unwrap();
        assert!(p.edits.starts_with(&before.edits));
        assert_eq!(
            sim.mutation_count - sims[parent].mutation_count,
            node.mutations_from_parent
        );
        // The ledger is a physical-ancestry record and holds every write,
        // including one whose drawn base equalled the base already there.
        // `mutation_count` holds the writes that CHANGED something, so the
        // two agree only once the no-op entries are filtered out -- 7 ledger
        // substitutions against a count of 6 on this family otherwise. The
        // ledger is deliberately not redefined to match the counter: it
        // records what happened to the molecule, which is its job.
        assert_eq!(
            p.edits
                .iter()
                .filter(|e| match e.action {
                    DuplexEditAction::Substitution {
                        old_base, new_base, ..
                    } => !old_base.eq_ignore_ascii_case(&new_base),
                    _ => false,
                })
                .count(),
            sim.mutation_count as usize
        );
        let projected = project_duplex_simulation(sim, &refdata).unwrap();
        assert_eq!(projected["record"]["n_mutations"], sim.mutation_count);
        let per_segment: u64 = ["v", "d", "j", "np"]
            .iter()
            .map(|s| {
                projected["record"][format!("n_{s}_mutations")]
                    .as_u64()
                    .unwrap()
            })
            .sum();
        assert_eq!(per_segment, sim.mutation_count as u64);
        let edge = edges[node.id as usize].as_ref().unwrap();
        assert_eq!(edge.revisions.len(), 2);
        assert_eq!(edge.pass_names, ["mutate.s5f"]);
        // Every write fires a BaseChanged event, including one whose drawn
        // base equalled the base already there. `mutations_from_parent` is a
        // difference of mutation counts, which exclude those, so the event
        // stream is filtered to the writes that changed something -- 7 events
        // against 6 mutations on this family otherwise. The event stream is
        // deliberately not narrowed: it records what the pass did.
        let actual = edge.events[0]
            .simulation_events
            .iter()
            .filter(|e| match e {
                crate::ir::SimulationEvent::BaseChanged {
                    old_base, new_base, ..
                } => !old_base.eq_ignore_ascii_case(new_base),
                _ => false,
            })
            .count();
        assert_eq!(actual, node.mutations_from_parent as usize);
    }
    assert!(founder.duplex_provenance.as_ref().unwrap().edits.is_empty());
    assert_eq!(founder.mutation_count, 0);
}

#[test]
fn duplex_recorded_edges_replay_exactly_in_both_sampling_policies() {
    for strict in [false, true] {
        let (refdata, _, tree, sims, edges) = grown(0.2, strict);
        let contracts = crate::contract::productive();
        let context = LineageMutationContext {
            refdata: &refdata,
            contracts: strict.then_some(&contracts),
            strict,
        };
        let mutator = S5FMutationPass::new_rate(kernel(), 0.2);
        for node in tree.nodes.iter().skip(1) {
            let parent = &sims[node.parent_id.unwrap() as usize];
            let replayed = replay_mutation_edge(
                parent,
                edges[node.id as usize].as_ref().unwrap(),
                &mutator,
                &context,
            )
            .unwrap();
            assert_eq!(
                replayed.final_simulation().duplex_provenance,
                sims[node.id as usize].duplex_provenance
            );
            assert_eq!(
                replayed.trace.choices(),
                edges[node.id as usize].as_ref().unwrap().trace.choices()
            );
            if strict {
                assert_eq!(
                    crate::productivity::assess(replayed.final_simulation(), &refdata).productive,
                    Some(true)
                );
            }
        }
    }
}

#[test]
fn duplex_edge_replay_rejects_changed_parent_trace_and_events() {
    let (refdata, _, tree, sims, edges) = grown(0.2, false);
    let node = &tree.nodes[1];
    let parent = &sims[node.parent_id.unwrap() as usize];
    let recorded = edges[1].as_ref().unwrap();
    let context = LineageMutationContext {
        refdata: &refdata,
        contracts: None,
        strict: false,
    };
    let mutator = S5FMutationPass::new_rate(kernel(), 0.2);
    let changed_parent = parent.with_mutation_count(parent.mutation_count + 1);
    assert!(replay_mutation_edge(&changed_parent, recorded, &mutator, &context).is_err());
    let mut changed = recorded.clone();
    changed.events[0].simulation_events.clear();
    assert!(replay_mutation_edge(parent, &changed, &mutator, &context).is_err());
    let mut changed = recorded.clone();
    changed
        .trace
        .record("unexpected", crate::trace::ChoiceValue::Int(4));
    assert!(replay_mutation_edge(parent, &changed, &mutator, &context).is_err());
    let different =
        S5FMutationPass::new_rate(S5FKernel::new(vec![0.0; 1024], vec![0.0; 4096]), 0.2);
    assert!(replay_mutation_edge(parent, recorded, &different, &context).is_err());
    let mut changed = recorded.clone();
    changed.trace = crate::trace::Trace::new();
    changed
        .trace
        .record("mutate.s5f.count", crate::trace::ChoiceValue::Int(3));
    assert!(replay_mutation_edge(parent, &changed, &mutator, &context).is_err());
    changed
        .trace
        .record("mutate.s5f.site[0]", crate::trace::ChoiceValue::Int(5));
    assert!(replay_mutation_edge(parent, &changed, &mutator, &context).is_err());
    let mut changed = recorded.clone();
    changed.trace = crate::trace::Trace::new();
    changed.trace.record("mutate.s5f.count", crate::trace::ChoiceValue::Int(3));
    changed.trace.record("mutate.s5f.site[0]", crate::trace::ChoiceValue::Int(i64::from(u32::MAX)));
    assert!(replay_mutation_edge(parent, &changed, &mutator, &context).is_err());
}

#[test]
fn duplex_zero_mutability_count_only_division_replays_in_both_policies() {
    let (refdata, model) = fixture();
    let founder = model.prepare_row(0, &refdata).unwrap().seed;
    let mutator = S5FMutationPass::new_rate(S5FKernel::new(vec![0.0; 1024], vec![1.0; 4096]), 0.2);
    for strict in [false, true] {
        let context = LineageMutationContext {
            refdata: &refdata,
            contracts: None,
            strict,
        };
        let (tree, sims, edges) = family::simulate_family_edges_with_context(
            &founder,
            &params(),
            &mutator,
            None,
            &context,
        )
        .unwrap();
        for node in tree.nodes.iter().skip(1) {
            let edge = edges[node.id as usize].as_ref().unwrap();
            assert_eq!(edge.trace.choices().len(), 1);
            assert_eq!(edge.final_simulation().mutation_count, 0);
            replay_mutation_edge(
                &sims[node.parent_id.unwrap() as usize],
                edge,
                &mutator,
                &context,
            )
            .unwrap();
        }
    }
}

#[test]
fn duplex_zero_destination_rows_record_sites_and_replay_without_rng() {
    let (refdata, model) = fixture();
    let founder = model.prepare_row(0, &refdata).unwrap().seed;
    let mutator = S5FMutationPass::new_rate(S5FKernel::new(vec![1.0; 1024], vec![0.0; 4096]), 0.2);
    let context = LineageMutationContext {
        refdata: &refdata,
        contracts: None,
        strict: false,
    };
    let (tree, sims, edges) =
        family::simulate_family_edges_with_context(&founder, &params(), &mutator, None, &context)
            .unwrap();
    for node in tree.nodes.iter().skip(1) {
        let edge = edges[node.id as usize].as_ref().unwrap();
        let crate::trace::ChoiceValue::Int(count) = edge.trace.choices()[0].value else {
            panic!("count is not typed integer")
        };
        assert_eq!(edge.trace.choices().len(), 1 + count as usize);
        assert!(edge.trace.choices()[1..]
            .iter()
            .all(|c| c.address.starts_with("mutate.s5f.site[")));
        assert_eq!(edge.final_simulation().mutation_count, 0);
        replay_mutation_edge(
            &sims[node.parent_id.unwrap() as usize],
            edge,
            &mutator,
            &context,
        )
        .unwrap();
    }
}

#[test]
fn duplex_permissive_empty_contract_support_replays_but_strict_rejects() {
    let (refdata, model) = fixture();
    let founder = model.prepare_row(0, &refdata).unwrap().seed;
    let mutability = {
        let mut values = vec![0.0; 1024];
        values[S5FKernel::encode_context(b'T', b'G', b'T', b'A', b'C').unwrap() as usize] = 1.0;
        values
    };
    let substitutions: Vec<_> = (0..1024).flat_map(|_| [1.0, 0.0, 0.0, 0.0]).collect();
    let mutator = S5FMutationPass::new_rate(S5FKernel::new(mutability, substitutions), 0.2);
    let contracts = crate::contract::productive();
    let context = LineageMutationContext {
        refdata: &refdata,
        contracts: Some(&contracts),
        strict: false,
    };
    let (tree, sims, edges) =
        family::simulate_family_edges_with_context(&founder, &params(), &mutator, None, &context)
            .unwrap();
    assert!(tree.len() > 1);
    for node in tree.nodes.iter().skip(1) {
        let edge = edges[node.id as usize].as_ref().unwrap();
        assert_eq!(edge.trace.choices().len(), 1);
        assert_eq!(edge.final_simulation().mutation_count, 0);
        replay_mutation_edge(
            &sims[node.parent_id.unwrap() as usize],
            edge,
            &mutator,
            &context,
        )
        .unwrap();
    }
    let strict = LineageMutationContext {
        strict: true,
        ..context
    };
    assert!(family::simulate_family_edges_with_context(
        &founder,
        &params(),
        &mutator,
        None,
        &strict
    )
    .is_err());
    assert!(replay_mutation_edge(&founder, edges[1].as_ref().unwrap(), &mutator, &strict).is_err());
}

#[test]
fn equal_dna_in_different_cells_keeps_separate_sampled_ancestry() {
    let (_, _, tree, _, _) = grown(0.0, false);
    let observed: Vec<_> = tree.nodes.iter().filter(|n| n.observed).collect();
    assert!(observed.len() > 1);
    assert_eq!(
        observed.iter().map(|n| n.abundance).sum::<u32>(),
        params().n_sample
    );
    assert!(observed.windows(2).all(|w| w[0].genotype == w[1].genotype));
    assert!(observed
        .iter()
        .all(|n| n.generation == params().max_generations));
}

#[test]
fn duplex_native_founder_validation_rejects_missing_reference_and_legacy_contract() {
    let (refdata, model) = fixture();
    let founder = model.prepare_row(0, &refdata).unwrap().seed;
    let mutator = S5FMutationPass::new_rate(kernel(), 0.0);
    assert!(branching::validate_duplex_founder(&founder, None, Some(&mutator)).is_err());
    let contracts = ContractSet::new().with(Box::new(AnchorPreserved::new(Segment::V)));
    let context = LineageMutationContext {
        refdata: &refdata,
        contracts: Some(&contracts),
        strict: false,
    };
    assert!(branching::validate_duplex_founder(&founder, Some(&context), Some(&mutator)).is_err());
    let mut wrong = refdata.clone();
    wrong.v_pool.push(
        wrong
            .v_pool
            .get(crate::refdata::AlleleId::new(0))
            .unwrap()
            .clone(),
    );
    let context = LineageMutationContext {
        refdata: &wrong,
        contracts: None,
        strict: false,
    };
    assert!(branching::validate_duplex_founder(&founder, Some(&context), Some(&mutator)).is_err());
}

#[test]
fn duplex_sampling_uses_same_literal_uniform_cell_draws() {
    let (_, _, mut tree, _, _) = grown(0.0, false);
    for n in &mut tree.nodes {
        n.observed = false;
        n.abundance = 0;
    }
    let living: Vec<_> = tree
        .nodes
        .iter()
        .filter(|n| n.generation == params().max_generations)
        .map(|n| n.id)
        .collect();
    let mut expected = vec![0; tree.len()];
    let mut oracle = Rng::new(128);
    for _ in 0..100 {
        expected[living[oracle.range_u32(living.len() as u32) as usize] as usize] += 1;
    }
    let mut rng = Rng::new(128);
    sampling::sample_cells(&mut tree, &living, 100, &mut rng);
    assert_eq!(
        tree.nodes.iter().map(|n| n.abundance).collect::<Vec<_>>(),
        expected
    );
    assert_eq!(rng.words_consumed(), oracle.words_consumed());
}

#[test]
fn ordinary_duplex_source_replays_zero_mutability_and_zero_destination_sites() {
    use crate::compiled::{CompileOptions, ExecutionPolicy, OwnedCompiledSimulator};
    use crate::dist::AllelePoolDist;
    use crate::pass::PassPlan;
    use crate::passes::{DuplexRecombinePass, SampleAllelePass};
    for (mutability, substitution) in [
        (vec![0.0; 1024], vec![1.0; 4096]),
        (vec![1.0; 1024], vec![0.0; 4096]),
        (
            vec![1.0; 1024],
            (0..4096)
                .map(|i| if i / 4 % 2 == 0 { 0.0 } else { 1.0 })
                .collect(),
        ),
    ] {
        let (refdata, model) = fixture();
        let mut plan = PassPlan::new();
        for segment in [Segment::V, Segment::D, Segment::J] {
            plan.push(Box::new(SampleAllelePass::new(
                segment,
                Box::new(AllelePoolDist::uniform(refdata.pool_for(segment).unwrap())),
            )));
        }
        plan.push(Box::new(
            DuplexRecombinePass::new(std::sync::Arc::new(model)).unwrap(),
        ));
        plan.push(Box::new(S5FMutationPass::new_rate(
            S5FKernel::new(mutability, substitution),
            0.2,
        )));
        let compiled = OwnedCompiledSimulator::compile_with_options(
            plan,
            Some(refdata),
            None,
            ExecutionPolicy::Permissive,
            CompileOptions::skip_refdata_validation(),
        )
        .unwrap();
        for seed in 1..16 {
            let original = compiled.run_one(seed).unwrap();
            let replayed = compiled
                .replay_from_trace_records(original.trace.choices(), 999, ExecutionPolicy::Strict)
                .unwrap();
            assert_eq!(
                replayed.final_simulation().duplex_provenance,
                original.final_simulation().duplex_provenance
            );
            assert_eq!(
                replayed.final_simulation().mutation_count,
                original.final_simulation().mutation_count
            );
            assert_eq!(replayed.events, original.events);
        }
    }
}

#[test]
fn finite_duplex_trace_binds_numeric_s5f_kernel_without_changing_legacy_signature() {
    use crate::dist::AllelePoolDist;
    use crate::pass::{Pass, PassPlan};
    use crate::passes::{DuplexRecombinePass, SampleAllelePass};
    let (refdata, model) = fixture();
    let make_plan = |weight: f64, duplex: bool| {
        let mut plan = PassPlan::new();
        if duplex {
            for segment in [Segment::V, Segment::D, Segment::J] {
                plan.push(Box::new(SampleAllelePass::new(
                    segment,
                    Box::new(AllelePoolDist::uniform(refdata.pool_for(segment).unwrap())),
                )));
            }
            plan.push(Box::new(
                DuplexRecombinePass::new(std::sync::Arc::new(model.clone())).unwrap(),
            ));
        }
        plan.push(Box::new(
            S5FMutationPass::new_rate(S5FKernel::new(vec![weight; 1024], vec![1.0; 4096]), 0.2)
                .with_kernel_name("same_name"),
        ));
        plan
    };
    let a = make_plan(1.0, true);
    let b = make_plan(2.0, true);
    let trace = crate::trace_file::TraceFile::build(&a, &refdata, 1, crate::trace::Trace::new());
    assert!(trace.validate_for_replay(&b, &refdata).is_err());
    assert_ne!(
        crate::trace_file::pass_plan_signature(&a),
        crate::trace_file::pass_plan_signature(&b)
    );
    assert_eq!(
        crate::trace_file::pass_plan_signature(&make_plan(1.0, false)),
        crate::trace_file::pass_plan_signature(&make_plan(2.0, false))
    );
    let p = S5FMutationPass::new_rate(kernel(), 0.2).with_kernel_name("same_name");
    assert!(p
        .sampling_law_signature()
        .contains("numeric_kernel=sha256:"));
    assert!(!p.parameter_signature().contains("numeric_kernel"));
}

#[test]
fn exact_back_mutations_keep_events_counts_and_cell_identity() {
    let (refdata, model) = fixture();
    let founder = model.prepare_row(0, &refdata).unwrap().seed;
    let a = S5FKernel::encode_context(b'T', b'G', b'T', b'A', b'C').unwrap() as usize;
    let b = S5FKernel::encode_context(b'T', b'G', b'C', b'A', b'C').unwrap() as usize;
    let mut mutability = vec![0.0; 1024];
    mutability[a] = 1.0;
    mutability[b] = 1.0;
    let mut substitutions = vec![0.0; 4096];
    substitutions[a * 4 + 1] = 1.0; // T -> C
    substitutions[b * 4 + 3] = 1.0; // C -> T
    let mutator = S5FMutationPass::new(
        S5FKernel::new(mutability, substitutions),
        Box::new(crate::dist::EmpiricalLengthDist::from_pairs(vec![(2, 1.0)])),
    );
    let context = LineageMutationContext {
        refdata: &refdata,
        contracts: None,
        strict: true,
    };
    let (tree, sims, edges) =
        family::simulate_family_edges_with_context(&founder, &params(), &mutator, None, &context)
            .unwrap();
    for node in tree.nodes.iter().skip(1) {
        let sim = &sims[node.id as usize];
        assert_eq!(node.genotype, tree.nodes[0].genotype);
        assert_eq!(node.mutations_from_parent, 2);
        assert_eq!(sim.mutation_count, node.generation * 2);
        assert_eq!(sim.mutation_burden, f64::from(node.generation) * 2.0 / 15.0);
        replay_mutation_edge(
            &sims[node.parent_id.unwrap() as usize],
            edges[node.id as usize].as_ref().unwrap(),
            &mutator,
            &context,
        )
        .unwrap();
    }
    assert!(tree.nodes.iter().filter(|n| n.observed).count() > 1);
}
