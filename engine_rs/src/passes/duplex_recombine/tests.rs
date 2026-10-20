use super::*;
use crate::assignment::AlleleInstance;
use crate::compiled::{CompileOptions, CompiledSimulator, ExecutionPolicy, OwnedCompiledSimulator};
use crate::dist::AllelePoolDist;
use crate::duplex_repair::model::tests::fixture;
use crate::pass::PassPlan;
use crate::passes::SampleAllelePass;
use crate::refdata::RefDataConfig;
use crate::rng::Rng;
use crate::trace::Trace;

fn plan(refdata: &RefDataConfig, model: FiniteDuplexModel) -> PassPlan {
    let mut plan = PassPlan::new();
    for segment in [Segment::V, Segment::D, Segment::J] {
        plan.push(Box::new(SampleAllelePass::new(
            segment,
            Box::new(AllelePoolDist::uniform(refdata.pool_for(segment).unwrap())),
        )));
    }
    plan.push(Box::new(DuplexRecombinePass::new(Arc::new(model)).unwrap()));
    plan
}

fn initial(seed: &Simulation) -> Simulation {
    let mut sim = Simulation::new();
    for segment in [Segment::V, Segment::D, Segment::J] {
        let source = seed.assignments.get(segment).unwrap();
        let mut instance = AlleleInstance::new(source.allele_id);
        instance.orientation = source.orientation;
        instance.haplotype = Some(1);
        sim = sim.with_allele_assigned(segment, instance);
    }
    let mut carried = crate::genotype::Haplotype::new();
    for segment in [Segment::V, Segment::D, Segment::J] {
        carried.set(
            segment,
            crate::refdata::GeneId::new(0),
            vec![crate::genotype::GeneCopy {
                allele: sim.assignments.get(segment).unwrap().allele_id,
                copies: 1,
                weight: 1.0,
            }],
        );
    }
    let genotype =
        crate::genotype::Genotype::new([carried.clone(), carried], [1.0, 1.0], None, String::new());
    sim.with_live_call_support(Arc::new(crate::live_call::LiveCallSupport::from_genotype(
        &genotype,
    )))
}

#[test]
fn weighted_authored_rows_remain_distinct_even_when_products_are_identical() {
    let (refdata, model) = fixture();
    let compiled = OwnedCompiledSimulator::compile_with_options(
        plan(&refdata, model.clone()),
        Some(refdata),
        None,
        ExecutionPolicy::Strict,
        CompileOptions::skip_refdata_validation(),
    )
    .unwrap();
    let n = 2048;
    let mut counts = [0; 3];
    for seed in 1..=n {
        let out = compiled.run_one(seed).unwrap();
        let sim = out.final_simulation();
        let p = sim.duplex_provenance.as_ref().unwrap();
        let selection = p.founder.selection.as_ref().unwrap();
        counts[selection.row_index] += 1;
        assert_eq!(selection.history_id, model.rows[selection.row_index].id);
        assert_eq!(p.founder.annotated.product.sequence, "AAATGTACTTTTAAA");
        assert_eq!(
            out.trace.choices().last().unwrap().address,
            "duplex_repair.v1.history"
        );
    }
    let expected = n as f64 * 0.75;
    let sigma = (n as f64 * 0.75 * 0.25).sqrt();
    assert!(
        (counts[1] as f64 - expected).abs() < 5.0 * sigma,
        "{counts:?}"
    );
    assert_eq!(counts[2], 0);
}

#[test]
fn replay_reexecutes_recorded_history_and_rejects_invalid_support() {
    let (refdata, model) = fixture();
    let compiled = OwnedCompiledSimulator::compile_with_options(
        plan(&refdata, model),
        Some(refdata),
        None,
        ExecutionPolicy::Permissive,
        CompileOptions::skip_refdata_validation(),
    )
    .unwrap();
    let original = compiled.run_one(20).unwrap();
    let replayed = compiled
        .replay_from_trace_records(original.trace.choices(), 999, ExecutionPolicy::Strict)
        .unwrap();
    assert_eq!(
        original.final_simulation().duplex_provenance,
        replayed.final_simulation().duplex_provenance
    );
    for invalid in [-1, 2, 3, i64::MAX] {
        let mut choices = original.trace.choices().to_vec();
        choices.last_mut().unwrap().value = ChoiceValue::Int(invalid);
        let error = compiled
            .replay_from_trace_records(&choices, 20, ExecutionPolicy::Strict)
            .unwrap_err();
        assert_eq!(error.pass_name(), "duplex.recombine");
    }
}

#[test]
fn atomic_installation_replays_all_body_caches_and_keeps_haplotype() {
    let (refdata, model) = fixture();
    let source = DuplexRecombinePass::new(Arc::new(model.clone())).unwrap();
    let plan = plan(&refdata, model.clone());
    let compiled = CompiledSimulator::compile_with_options(
        &plan,
        Some(&refdata),
        None,
        ExecutionPolicy::Strict,
        CompileOptions::skip_refdata_validation(),
    )
    .unwrap();
    let prepared = model.prepare_row(0, &refdata).unwrap();
    let initial = initial(&prepared.seed);
    let mut rng = Rng::new(1);
    let mut trace = Trace::new();
    let mut events = Vec::new();
    let mut ctx = PassContext {
        rng: &mut rng,
        trace: &mut trace,
        pass_index: 3,
        refdata: Some(&refdata),
        contracts: None,
        feasibility: compiled.feasibility(),
        reference_index: None,
        replay_cursor: None,
        event_log_sink: Some(&mut events),
    };
    let result = source.execute_checked(&initial, &mut ctx).unwrap();
    let replayed = crate::ir::event_log_observer::replay_events(initial.clone(), &events);
    assert_eq!(result.pool.as_slice(), replayed.pool.as_slice());
    assert_eq!(result.assignments, replayed.assignments);
    assert_eq!(result.sequence.regions, replayed.sequence.regions);
    assert_eq!(result.duplex_provenance, replayed.duplex_provenance);
    assert_eq!(result.live_call_support, initial.live_call_support);
    for segment in [Segment::V, Segment::D, Segment::J] {
        assert_eq!(result.assignments.get(segment).unwrap().haplotype, Some(1));
    }
    result
        .duplex_provenance
        .as_ref()
        .unwrap()
        .validate_pool(result.pool.as_slice())
        .unwrap();
    let mut cursor = crate::replay::TraceCursor::new(trace.choices());
    let mut replay_trace = Trace::new();
    let mut replay_rng = Rng::new(999);
    let mut replay_ctx = PassContext {
        rng: &mut replay_rng,
        trace: &mut replay_trace,
        pass_index: 3,
        refdata: Some(&refdata),
        contracts: None,
        feasibility: compiled.feasibility(),
        reference_index: None,
        replay_cursor: Some(&mut cursor),
        event_log_sink: None,
    };
    let regenerated = source.execute_checked(&initial, &mut replay_ctx).unwrap();
    assert_eq!(result.duplex_provenance, regenerated.duplex_provenance);
    assert_eq!(replay_rng.words_consumed(), 0);
    cursor.assert_drained().unwrap();
}

#[test]
fn failed_atomic_installation_does_not_emit_events_or_partial_query() {
    let (refdata, model) = fixture();
    let prepared = model.prepare_row(0, &refdata).unwrap();
    let initial =
        initial(&prepared.seed).with_trim(Segment::V, crate::assignment::TrimEnd::Three, 1);
    let mut builder = SimulationBuilder::from_simulation(initial.clone());
    builder.attach_event_log_observer();
    assert!(builder.install_duplex_seed(&prepared.seed).is_err());
    assert!(builder.seal_event_log_observer().is_empty());
    assert_eq!(builder.peek().pool.as_slice(), initial.pool.as_slice());
    assert_eq!(builder.peek().assignments, initial.assignments);
    assert!(builder.peek().duplex_provenance.is_none());
}

#[test]
fn history_address_is_a_singleton_round_trip() {
    let address = ChoiceAddress::DuplexRepairHistory;
    assert_eq!(ChoiceAddress::parse(&address.to_string()), Some(address));
    let pattern = ChoiceAddressPattern::DuplexRepairHistory;
    assert_eq!(
        ChoiceAddressPattern::parse(&pattern.to_string()),
        Some(pattern)
    );
}
