use super::*;
use crate::compiled::CompiledSimulator;
use crate::dist::{EmpiricalLengthDist, UniformBase};
use crate::duplex_repair::bindings::DuplexReferenceBindings;
use crate::duplex_repair::{repair_coupled_annotated, JoinName, JoinScript, RepairOperation, Side};
use crate::ir::Segment;
use crate::pass::Pass;
use crate::passes::*;
use crate::refdata::{Allele, ChainType};
use std::sync::atomic::{AtomicBool, Ordering};
use std::sync::Arc;

fn fixed(count: i64) -> Box<EmpiricalLengthDist> {
    Box::new(EmpiricalLengthDist::from_pairs(vec![(count, 1.0)]))
}

fn seed() -> (RefDataConfig, Simulation) {
    seed_with_orientation(false)
}

fn seed_with_orientation(inverted: bool) -> (RefDataConfig, Simulation) {
    let mut cfg = RefDataConfig::empty(ChainType::Vdj);
    for (segment, name, sequence, anchor) in [
        (Segment::V, "V*01", "AAATGT", Some(3)),
        (Segment::D, "D*01", "CCC", None),
        (Segment::J, "J*01", "TTTAAA", Some(0)),
    ] {
        let allele = Allele {
            name: name.into(),
            gene: name.split('*').next().unwrap().into(),
            seq: sequence.as_bytes().to_vec(),
            segment,
            anchor,
            functional_status: None,
            subregions: vec![],
        };
        let _ = match segment {
            Segment::V => cfg.v_pool.push(allele),
            Segment::D => cfg.d_pool.push(allele),
            Segment::J => cfg.j_pool.push(allele),
            _ => unreachable!(),
        };
    }
    let joins: Vec<_> = [(JoinName::VD, 6), (JoinName::DJ, 3)]
        .into_iter()
        .map(|(join, right_start)| JoinScript {
            join,
            operations: vec![
                RepairOperation::Open {
                    side: Side::Left,
                    offset: 0,
                },
                RepairOperation::Open {
                    side: Side::Right,
                    offset: 0,
                },
                RepairOperation::Synapse { right_start },
                RepairOperation::Ligate {},
            ],
        })
        .collect();
    let annotated =
        repair_coupled_annotated("AAATGT", "CCC", "TTTAAA", inverted, &joins, 256, 100_000)
            .unwrap();
    let bindings = DuplexReferenceBindings::new(
        &cfg,
        ["V*01", "D*01", "J*01"],
        ["AAATGT", "CCC", "TTTAAA"],
        inverted,
        true,
    )
    .unwrap();
    let sim = crate::duplex_repair::seed::seed_from_coupled(annotated, bindings, &cfg).unwrap();
    (cfg, sim)
}

#[test]
fn duplex_insertions_keep_an_exact_partition_at_every_query_boundary() {
    use crate::duplex_repair::observation::{layout::validate_layout, project_duplex_simulation};
    use crate::ir::{NucFlags, Nucleotide, SimulationBuilder};
    let (cfg, initial) = seed();
    for at in 0..=initial.pool.len() as u32 {
        let segment = if at < 6 {
            Segment::V
        } else if at < 9 {
            Segment::D
        } else {
            Segment::J
        };
        let mut builder = SimulationBuilder::from_simulation(initial.clone());
        builder.attach_event_log_observer();
        builder.insert_indel(at, Nucleotide::synthetic(b'G', segment, NucFlags::empty()));
        let events = builder.seal_event_log_observer();
        let changed = builder.seal();
        validate_layout(&changed).unwrap();
        let replayed = crate::ir::event_log_observer::replay_events(initial.clone(), &events);
        validate_layout(&replayed).unwrap();
        assert_eq!(
            changed
                .sequence
                .regions
                .iter()
                .map(|r| (r.segment, r.start, r.end))
                .collect::<Vec<_>>(),
            replayed
                .sequence
                .regions
                .iter()
                .map(|r| (r.segment, r.start, r.end))
                .collect::<Vec<_>>()
        );
        let projected = project_duplex_simulation(&changed, &cfg).unwrap();
        let record = &projected["record"];
        assert_eq!(
            record["sequence_alignment"]
                .as_str()
                .unwrap()
                .replace('-', ""),
            record["sequence"].as_str().unwrap()
        );
        if at == initial.pool.len() as u32 {
            assert_eq!(record["j_cigar"], "6M1I");
            assert_eq!(record["j_sequence_end"], 16);
        }
    }
}

#[test]
fn duplex_insertion_into_empty_query_keeps_zero_width_sources_in_order() {
    use crate::duplex_repair::observation::layout::validate_layout;
    use crate::ir::{NucFlags, Nucleotide, SimulationBuilder};
    let (cfg, initial) = seed();
    for observation_loss in [false, true] {
        let empty = if observation_loss {
            let mut plan = PassPlan::new();
            plan.push(Box::new(EndLossPass::new(
                LossEnd::Five,
                fixed(initial.pool.len() as i64),
            )));
            let compiled =
                CompiledSimulator::compile(&plan, Some(&cfg), None, ExecutionPolicy::Strict)
                    .unwrap();
            execute_transactional(compiled.inputs(), initial.clone(), 1)
                .unwrap()
                .final_simulation()
                .clone()
        } else {
            let mut builder = SimulationBuilder::from_simulation(initial.clone());
            for _ in 0..initial.pool.len() {
                builder.delete_indel(0);
            }
            builder.seal()
        };
        let mut builder = SimulationBuilder::from_simulation(empty);
        builder.insert_indel(
            0,
            Nucleotide::synthetic(b'A', Segment::V, NucFlags::empty()),
        );
        let changed = builder.seal();
        validate_layout(&changed).unwrap();
        assert_eq!(
            changed
                .sequence
                .regions
                .iter()
                .map(|r| (r.segment, r.start.index(), r.end.index()))
                .collect::<Vec<_>>(),
            [(Segment::V, 0, 1), (Segment::D, 1, 1), (Segment::J, 1, 1)]
        );
    }
}

#[test]
fn duplex_coverage_cache_must_match_audited_loss_and_preserves_repeated_inverted_clips() {
    use crate::duplex_repair::observation::project_duplex_simulation;
    use crate::ir::{NucHandle, SimulationBuilder};
    for inverted in [false, true] {
        let (cfg, initial) = seed_with_orientation(inverted);
        for loss_end in [LossEnd::Five, LossEnd::Three] {
            let mut first = PassPlan::new();
            first.push(Box::new(EndLossPass::new(loss_end, fixed(3))));
            let first =
                CompiledSimulator::compile(&first, Some(&cfg), None, ExecutionPolicy::Strict)
                    .unwrap();
            let parent = execute_transactional(first.inputs(), initial.clone(), 1).unwrap();
            let mut second = PassPlan::new();
            second.push(Box::new(EndLossPass::new(loss_end, fixed(5))));
            let second =
                CompiledSimulator::compile(&second, Some(&cfg), None, ExecutionPolicy::Strict)
                    .unwrap();
            let child =
                execute_transactional(second.inputs(), parent.final_simulation().clone(), 2)
                    .unwrap();
            let projected = project_duplex_simulation(child.final_simulation(), &cfg).unwrap();
            assert_eq!(projected["record"]["sequence_length"], 7);
            assert_eq!(projected["record"]["d_cigar"], "1M");
            let d = child
                .final_simulation()
                .assignments
                .get(Segment::D)
                .unwrap();
            let five = matches!(loss_end, LossEnd::Five);
            assert_eq!(
                (d.observed_ref_start, d.observed_ref_end),
                if five != inverted {
                    (Some(2), None)
                } else {
                    (None, Some(1))
                }
            );
            let mut builder = SimulationBuilder::from_simulation(child.final_simulation().clone());
            let mut forged = *d;
            forged.observed_ref_start = Some(0);
            forged.observed_ref_end = Some(0);
            builder.assign_allele(Segment::D, forged);
            assert!(project_duplex_simulation(&builder.seal(), &cfg).is_err());
        }
        let mut builder = SimulationBuilder::from_simulation(initial.clone());
        let mut forged = *initial.assignments.get(Segment::V).unwrap();
        forged.observation_loss_5 = 1;
        builder.assign_allele(Segment::V, forged);
        assert!(project_duplex_simulation(&builder.seal(), &cfg).is_err());
        let raw = initial.with_observation_deleted(0);
        assert!(project_duplex_simulation(&raw, &cfg)
            .unwrap_err()
            .contains("audited terminal"));
        // Ordinary literal substitutions retain the complete coverage history.
        assert!(project_duplex_simulation(
            &initial.with_base_changed(NucHandle::new(0), b'C'),
            &cfg
        )
        .is_ok());
    }
}

#[test]
fn duplex_preflight_and_projection_reject_forged_region_partition() {
    let (cfg, initial) = seed();
    let plan = PassPlan::new();
    let compiled =
        CompiledSimulator::compile(&plan, Some(&cfg), None, ExecutionPolicy::Strict).unwrap();
    for defect in 0..4 {
        let mut forged = initial.clone();
        let mut sequence = (*forged.sequence).clone();
        match defect {
            0 => {
                sequence.regions.remove(1);
            }
            1 => {
                sequence.regions[1].start = crate::ir::NucHandle::new(7);
            }
            2 => {
                sequence.regions.swap(0, 1);
            }
            3 => {
                sequence.regions[1].segment = Segment::Np1;
            }
            _ => unreachable!(),
        }
        forged.sequence = Arc::new(sequence);
        assert!(execute_transactional(compiled.inputs(), forged.clone(), 1).is_err());
        assert!(
            crate::duplex_repair::observation::project_duplex_simulation(&forged, &cfg).is_err()
        );
    }
}

struct SpoofedPass(Arc<AtomicBool>);

impl Pass for SpoofedPass {
    fn name(&self) -> &str {
        crate::address::CORRUPT_REV_COMP
    }

    fn execute(&self, sim: &Simulation, _: &mut PassContext) -> Simulation {
        self.0.store(true, Ordering::SeqCst);
        sim.clone()
    }
}

#[test]
fn duplex_whitelist_uses_concrete_types_not_names_or_effects() {
    let allowed: Vec<Box<dyn Pass>> = vec![
        Box::new(UniformMutationPass::new(fixed(0), Box::new(UniformBase))),
        Box::new(S5FMutationPass::new(
            crate::s5f::S5FKernel::new(vec![1.0; 1024], vec![0.25; 4096]),
            fixed(0),
        )),
        Box::new(PCRErrorPass::new(fixed(0), Box::new(UniformBase))),
        Box::new(QualityErrorPass::new(fixed(0), Box::new(UniformBase))),
        Box::new(NCorruptionPass::new(fixed(0))),
        Box::new(IndelPass::new(fixed(0), 0.5, Box::new(UniformBase))),
        Box::new(EndLossPass::new(LossEnd::Five, fixed(0))),
        Box::new(RevCompPass::new(0.0)),
    ];
    for pass in allowed {
        assert!(duplex_pass_supported(pass.as_ref()), "{}", pass.name());
    }
    let denied: Vec<Box<dyn Pass>> = vec![
        Box::new(SpoofedPass(Arc::new(AtomicBool::new(false)))),
        Box::new(ContaminantPass::new(0.0, Box::new(UniformBase))),
        Box::new(EchoPass::new(b'A', 0, Segment::V)),
        Box::new(AssembleSegmentPass::new(Segment::V)),
        Box::new(PairedEndSamplingPass::new(PairedEndLayoutSpec::new(
            fixed(3),
            fixed(3),
            fixed(6),
        ))),
    ];
    for pass in denied {
        assert!(!duplex_pass_supported(pass.as_ref()), "{}", pass.name());
    }
}

#[test]
fn duplex_preflights_entire_continuation_in_both_modes() {
    let (cfg, initial) = seed();
    let marker = Arc::new(AtomicBool::new(false));
    let mut plan = PassPlan::new();
    plan.push(Box::new(RevCompPass::new(1.0)));
    plan.push(Box::new(SpoofedPass(marker.clone())));
    for policy in [ExecutionPolicy::Permissive, ExecutionPolicy::Strict] {
        let compiled = CompiledSimulator::compile(&plan, Some(&cfg), None, policy).unwrap();
        let error = execute_transactional(compiled.inputs(), initial.clone(), 3).unwrap_err();
        assert!(matches!(error.error, PassError::InvalidPlanState { .. }));
        assert!(error
            .error
            .to_string()
            .contains("unsupported for a duplex seed"));
        assert_eq!(error.committed.revisions.len(), 1);
        assert!(error.committed.pass_names.is_empty());
        assert!(error.committed.events.is_empty());
        assert!(error.committed.trace.is_empty());
        assert_eq!(
            error.committed.final_simulation().duplex_provenance,
            initial.duplex_provenance
        );
        assert!(!marker.load(Ordering::SeqCst));
    }
}

#[test]
fn duplex_empty_plan_preserves_seed_and_ordinary_custom_passes_still_work() {
    let (cfg, initial) = seed();
    let plan = PassPlan::new();
    let compiled =
        CompiledSimulator::compile(&plan, Some(&cfg), None, ExecutionPolicy::Permissive).unwrap();
    let output = compiled.run_one_from(initial.clone(), 2).unwrap();
    assert_eq!(
        output.final_simulation().duplex_provenance,
        initial.duplex_provenance
    );
    let marker = Arc::new(AtomicBool::new(false));
    let mut ordinary = PassPlan::new();
    ordinary.push(Box::new(SpoofedPass(marker.clone())));
    let compiled =
        CompiledSimulator::compile(&ordinary, None, None, ExecutionPolicy::Permissive).unwrap();
    compiled.run_one(2).unwrap();
    assert!(marker.load(Ordering::SeqCst));
}

#[test]
fn duplex_missing_or_mismatched_reference_rejected_even_for_noop() {
    let (mut cfg, initial) = seed();
    let plan = PassPlan::new();
    let compiled =
        CompiledSimulator::compile(&plan, None, None, ExecutionPolicy::Permissive).unwrap();
    let error = compiled.run_one_from(initial.clone(), 0).unwrap_err();
    assert!(error.to_string().contains("bound reference data"));
    cfg.identity.source = Some("different reference provenance".into());
    let compiled =
        CompiledSimulator::compile(&plan, Some(&cfg), None, ExecutionPolicy::Permissive).unwrap();
    let error = compiled.run_one_from(initial, 0).unwrap_err();
    assert!(error.to_string().contains("binding identity"));
}

#[test]
fn duplex_legacy_preservation_contract_rejected_before_any_work() {
    let (cfg, initial) = seed();
    let plan = PassPlan::new();
    let contracts =
        ContractSet::new().with(Box::new(crate::contract::AnchorPreserved::new(Segment::V)));
    for policy in [ExecutionPolicy::Permissive, ExecutionPolicy::Strict] {
        let compiled =
            CompiledSimulator::compile(&plan, Some(&cfg), Some(&contracts), policy).unwrap();
        let error = compiled.run_one_from(initial.clone(), 0).unwrap_err();
        assert!(error.to_string().contains("legacy AnchorPreserved"));
    }
}

#[test]
fn duplex_rejects_rewritten_generating_assignment_on_reentry() {
    let (cfg, initial) = seed();
    let mut builder = crate::ir::SimulationBuilder::from_simulation(initial.clone());
    let mut assignment = *initial.assignments.get(Segment::V).unwrap();
    assignment.trim_3 += 1;
    builder.assign_allele(Segment::V, assignment);
    let altered = builder.seal();
    let plan = PassPlan::new();
    let compiled =
        CompiledSimulator::compile(&plan, Some(&cfg), None, ExecutionPolicy::Permissive).unwrap();
    let error = compiled.run_one_from(altered, 0).unwrap_err();
    assert!(error.to_string().contains("founder trims changed"));
}

#[test]
fn duplex_supported_reverse_pass_keeps_ancestry_and_records_presentation() {
    let (cfg, initial) = seed();
    let mut plan = PassPlan::new();
    plan.push(Box::new(RevCompPass::new(1.0)));
    let compiled =
        CompiledSimulator::compile(&plan, Some(&cfg), None, ExecutionPolicy::Strict).unwrap();
    let output = compiled.run_one_from(initial.clone(), 4).unwrap();
    let provenance = output
        .final_simulation()
        .duplex_provenance
        .as_ref()
        .unwrap();
    assert!(provenance.rev_comp);
    assert_eq!(
        provenance.founder,
        initial.duplex_provenance.as_ref().unwrap().founder
    );
    assert_eq!(
        output.final_simulation().pool.as_slice(),
        initial.pool.as_slice()
    );
}

#[test]
fn duplex_legacy_native_record_builder_fails_closed() {
    let (cfg, initial) = seed();
    let plan = PassPlan::new();
    let compiled =
        CompiledSimulator::compile(&plan, Some(&cfg), None, ExecutionPolicy::Permissive).unwrap();
    let outcome = compiled.run_one_from(initial, 2).unwrap();
    let error = crate::airr_record::try_build_airr_record(&outcome, &cfg, "duplex").unwrap_err();
    assert!(error.contains("unknown duplex P/NP"));
    let legacy = std::panic::catch_unwind(|| {
        crate::airr_record::build_airr_record(&outcome, &cfg, "duplex")
    });
    assert!(legacy.is_err());
}
