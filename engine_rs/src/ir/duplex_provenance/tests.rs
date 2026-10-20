use super::*;
use crate::duplex_repair::{repair_coupled_annotated, JoinName, JoinScript, RepairOperation, Side};
use crate::ir::{NucFlags, NucHandle, Simulation, SimulationBuilder};
use crate::refdata::{Allele, ChainType, RefDataConfig};

fn founder_with_reference(inverted: bool) -> (Simulation, RefDataConfig) {
    let mut cfg = RefDataConfig::empty(ChainType::Vdj);
    for (segment, name, seq, anchor) in [
        (Segment::V, "V*01", "AAATGT", Some(3)),
        (Segment::D, "D*01", "ACT", None),
        (Segment::J, "J*01", "TTTAAA", Some(0)),
    ] {
        let a = Allele {
            name: name.into(),
            gene: name.split('*').next().unwrap().into(),
            seq: seq.as_bytes().to_vec(),
            segment,
            anchor,
            functional_status: None,
            subregions: vec![],
        };
        match segment {
            Segment::V => cfg.v_pool.push(a),
            Segment::D => cfg.d_pool.push(a),
            Segment::J => cfg.j_pool.push(a),
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
        repair_coupled_annotated("AAATGT", "ACT", "TTTAAA", inverted, &joins, 32, 100).unwrap();
    let bindings = DuplexReferenceBindings::new(
        &cfg,
        ["V*01", "D*01", "J*01"],
        ["AAATGT", "ACT", "TTTAAA"],
        inverted,
        true,
    )
    .unwrap();
    let sim = crate::duplex_repair::seed::seed_from_coupled(annotated, bindings, &cfg).unwrap();
    (sim, cfg)
}

fn founder(inverted: bool) -> Simulation {
    founder_with_reference(inverted).0
}

#[test]
fn branch_edits_preserve_founder_and_distinguish_loss_from_deletion() {
    let original = founder(false);
    let changed = original
        .with_base_changed(NucHandle::new(1), b'C')
        .with_indel_inserted(
            2,
            Nucleotide::synthetic(b'G', Segment::V, NucFlags::empty()),
        )
        .with_indel_deleted(4)
        .with_observation_deleted(0);
    let p = changed.duplex_provenance.as_ref().unwrap();
    p.validate_pool(changed.pool.as_slice()).unwrap();
    assert_eq!(
        p.removed.iter().map(|r| r.cause).collect::<Vec<_>>(),
        [
            DuplexRemovalCause::StructuralDeletion,
            DuplexRemovalCause::ObservationLoss
        ]
    );
    assert_eq!(p.body_reference_position(Segment::V, 0), Some(1));
    assert_eq!(p.body_reference_position(Segment::V, 1), None);
    assert!(p.edits.iter().all(|e| e.source_pass.is_none()));
    assert!(Arc::ptr_eq(
        &p.founder,
        &original.duplex_provenance.as_ref().unwrap().founder
    ));
    assert!(original
        .duplex_provenance
        .as_ref()
        .unwrap()
        .edits
        .is_empty());
    assert_eq!(original.pool.get(NucHandle::new(1)).unwrap().base, b'A');
    assert!(matches!(
        p.founder.p_metrics.v3,
        DuplexMetric::Unavailable { .. }
    ));
}

#[test]
fn builder_event_replay_preserves_source_cause_and_latest_orientation() {
    let original = founder(false);
    let mut b = SimulationBuilder::from_simulation(original.clone());
    b.attach_event_log_observer();
    b.set_duplex_edit_context(Some("mutate.uniform"));
    b.change_base(NucHandle::new(1), b'A'); // actual same-byte event remains visible
    b.set_duplex_edit_context(Some("corrupt.indel"));
    b.insert_indel(
        3,
        Nucleotide::synthetic(b'C', Segment::V, NucFlags::empty()),
    );
    b.delete_indel(2);
    b.set_duplex_edit_context(Some("corrupt.end_loss"));
    b.delete_observation_base(0);
    b.record_reverse_complement_flag(true);
    b.record_reverse_complement_flag(false);
    let events = b.seal_event_log_observer();
    let changed = b.seal();
    let replayed = crate::ir::event_log_observer::replay_events(original, &events);
    assert_eq!(changed.duplex_provenance, replayed.duplex_provenance);
    let p = changed.duplex_provenance.as_ref().unwrap();
    p.validate_pool(changed.pool.as_slice()).unwrap();
    assert_eq!(p.edits[0].source_pass.as_deref(), Some("mutate.uniform"));
    assert_eq!(p.edits[3].source_pass.as_deref(), Some("corrupt.end_loss"));
    assert_eq!(p.edits.len(), 6);
    assert!(!p.rev_comp);
    assert!(events
        .iter()
        .any(|e| matches!(e, crate::ir::SimulationEvent::BaseDeleted { .. })));
}

#[test]
fn inverted_reference_cache_tracks_canonical_positions_after_edits() {
    let original = founder(true);
    let p = original.duplex_provenance.as_ref().unwrap();
    assert_eq!(p.body_reference_position(Segment::D, 6), Some(2));
    assert_eq!(p.body_reference_position(Segment::D, 8), Some(0));
    p.validate_pool(original.pool.as_slice()).unwrap();
    let changed = original
        .with_indel_deleted(7)
        .with_base_changed(NucHandle::new(6), b'N');
    changed
        .duplex_provenance
        .as_ref()
        .unwrap()
        .validate_pool(changed.pool.as_slice())
        .unwrap();
    assert_eq!(
        changed
            .duplex_provenance
            .as_ref()
            .unwrap()
            .body_reference_position(Segment::D, 7),
        Some(0)
    );
}

#[test]
fn invalid_cache_and_missing_or_duplicated_ledger_entries_are_rejected() {
    let original = founder(false);
    let p = original.duplex_provenance.as_ref().unwrap();
    let mut forged = (**p).clone();
    forged.links[0] = DuplexQueryLink::FounderColumn(1);
    assert!(forged.validate_mapping(original.pool.len()).is_err());
    let changed = original.with_indel_deleted(1);
    let mut forged = (**changed.duplex_provenance.as_ref().unwrap()).clone();
    forged.removed.clear();
    assert!(forged.validate_mapping(changed.pool.len()).is_err());
    let corrupt = original.with_nucleotide_changed(
        NucHandle::new(0),
        Nucleotide::synthetic(b'A', Segment::V, NucFlags::empty()),
    );
    assert!(corrupt
        .duplex_provenance
        .as_ref()
        .unwrap()
        .validate_pool(corrupt.pool.as_slice())
        .is_err());
}

#[test]
fn replacement_records_all_removed_origins_and_empty_extension_is_noop() {
    let original = founder(false);
    let (changed, _, _) = original.with_segment_replaced(
        Segment::D,
        vec![Nucleotide::synthetic(b'C', Segment::D, NucFlags::empty())],
    );
    let p = changed.duplex_provenance.as_ref().unwrap();
    p.validate_pool(changed.pool.as_slice()).unwrap();
    assert_eq!(p.removed.len(), 3);
    assert!(p
        .removed
        .iter()
        .all(|r| r.cause == DuplexRemovalCause::Replacement));
    let (same, _) = changed.with_nucleotides_extended([]);
    assert_eq!(same.duplex_provenance, changed.duplex_provenance);
}

#[test]
fn fully_clipped_query_preserves_sampled_error_attempts_with_no_applied_edits() {
    use crate::dist::{EmpiricalLengthDist, UniformBase};
    use crate::pass::{testing::PassRuntime, PassPlan};
    use crate::passes::{EndLossPass, LossEnd, PCRErrorPass, QualityErrorPass};
    let (initial, cfg) = founder_with_reference(false);
    let fixed = |n| Box::new(EmpiricalLengthDist::from_pairs(vec![(n, 1.0)]));
    let mut plan = PassPlan::new();
    plan.push(Box::new(EndLossPass::new(
        LossEnd::Five,
        fixed(initial.pool.len() as i64),
    )));
    plan.push(Box::new(PCRErrorPass::new(fixed(2), Box::new(UniformBase))));
    plan.push(Box::new(QualityErrorPass::new(
        fixed(3),
        Box::new(UniformBase),
    )));
    plan.push(Box::new(PCRErrorPass::new(fixed(0), Box::new(UniformBase))));
    let outcome =
        PassRuntime::execute_strict_with_context(&plan, initial.clone(), 81, Some(&cfg), None)
            .unwrap();
    let sim = outcome.final_simulation();
    assert!(sim.pool.is_empty());
    let p = sim.duplex_provenance.as_ref().unwrap();
    p.validate_pool(sim.pool.as_slice()).unwrap();
    assert_eq!(p.sampled_attempts(DuplexAttemptStage::Pcr), 2);
    assert_eq!(p.sampled_attempts(DuplexAttemptStage::Quality), 3);
    assert_eq!(p.attempt_stages.len(), 3);
    assert!(p
        .attempt_stages
        .iter()
        .all(|s| s.input_pool_len == 0 && s.edit_ordinal == 15));
    assert!(p.edits.iter().all(|e| matches!(
        e.action,
        DuplexEditAction::Deletion {
            cause: DuplexRemovalCause::ObservationLoss,
            ..
        }
    )));
    let events: Vec<_> = outcome
        .events()
        .iter()
        .flat_map(|e| e.simulation_events.clone())
        .collect();
    let replayed = crate::ir::event_log_observer::replay_events(initial, &events);
    assert_eq!(replayed.duplex_provenance, sim.duplex_provenance);
}

#[test]
fn attempt_stages_accumulate_once_across_repeated_passes_and_cloned_continuations() {
    use crate::dist::{EmpiricalLengthDist, UniformBase};
    use crate::pass::{testing::PassRuntime, PassPlan};
    use crate::passes::{PCRErrorPass, QualityErrorPass};
    let (initial, cfg) = founder_with_reference(false);
    let fixed = |n| Box::new(EmpiricalLengthDist::from_pairs(vec![(n, 1.0)]));
    let mut first = PassPlan::new();
    first.push(Box::new(PCRErrorPass::new(fixed(2), Box::new(UniformBase))));
    first.push(Box::new(QualityErrorPass::new(
        fixed(0),
        Box::new(UniformBase),
    )));
    let first_outcome =
        PassRuntime::execute_strict_with_context(&first, initial, 2, Some(&cfg), None).unwrap();
    let parent = first_outcome.final_simulation().clone();
    let mut next = PassPlan::new();
    next.push(Box::new(PCRErrorPass::new(fixed(3), Box::new(UniformBase))));
    next.push(Box::new(QualityErrorPass::new(
        fixed(1),
        Box::new(UniformBase),
    )));
    let outcome =
        PassRuntime::execute_strict_with_context(&next, parent.clone(), 3, Some(&cfg), None)
            .unwrap();
    let child = outcome.final_simulation();
    let p = child.duplex_provenance.as_ref().unwrap();
    p.validate_pool(child.pool.as_slice()).unwrap();
    assert_eq!(p.sampled_attempts(DuplexAttemptStage::Pcr), 5);
    assert_eq!(p.sampled_attempts(DuplexAttemptStage::Quality), 1);
    assert_eq!(
        p.attempt_stages
            .iter()
            .map(|s| (s.stage_ordinal, s.edit_ordinal))
            .collect::<Vec<_>>(),
        [(0, 0), (1, 2), (2, 2), (3, 5)]
    );
    let pp = parent.duplex_provenance.as_ref().unwrap();
    assert_eq!(pp.sampled_attempts(DuplexAttemptStage::Pcr), 2);
    assert_eq!(pp.attempt_stages.len(), 2);
    let events: Vec<_> = outcome
        .events()
        .iter()
        .flat_map(|e| e.simulation_events.clone())
        .collect();
    let replayed = crate::ir::event_log_observer::replay_events(parent.clone(), &events);
    assert_eq!(replayed.duplex_provenance, child.duplex_provenance);
    let mut forged = (**p).clone();
    forged.attempt_stages[0].input_pool_len = 0;
    assert!(forged.validate_pool(child.pool.as_slice()).is_err());
}

#[test]
fn unrecorded_pool_bytes_and_forged_source_bytes_are_rejected() {
    let original = founder(false);
    let mut forged = original.clone();
    forged.pool = forged.pool.with_base_changed(NucHandle::new(1), b'C');
    assert!(forged
        .duplex_provenance
        .as_ref()
        .unwrap()
        .validate_pool(forged.pool.as_slice())
        .unwrap_err()
        .contains("recorded observation history"));
    let recorded = original.with_base_changed(NucHandle::new(1), b'n');
    recorded
        .duplex_provenance
        .as_ref()
        .unwrap()
        .validate_pool(recorded.pool.as_slice())
        .unwrap();
    let mut forged_history = (**recorded.duplex_provenance.as_ref().unwrap()).clone();
    if let DuplexEditAction::Substitution { old_base, .. } = &mut forged_history.edits[0].action {
        *old_base = b'G';
    } else {
        panic!("expected recorded substitution");
    }
    assert!(forged_history
        .validate_pool(recorded.pool.as_slice())
        .unwrap_err()
        .contains("source byte"));
    let deleted = recorded.with_indel_deleted(1);
    deleted
        .duplex_provenance
        .as_ref()
        .unwrap()
        .validate_pool(deleted.pool.as_slice())
        .unwrap();
    let mut forged_history = (**deleted.duplex_provenance.as_ref().unwrap()).clone();
    if let DuplexEditAction::Deletion { base, .. } = &mut forged_history.edits[1].action {
        *base = b'A';
    } else {
        panic!("expected recorded deletion");
    }
    assert!(forged_history
        .validate_pool(deleted.pool.as_slice())
        .unwrap_err()
        .contains("deleted byte"));
}

#[test]
fn replay_accounts_for_mutated_insertions_and_replacement_bases() {
    let original = founder(false);
    let changed = original
        .with_indel_inserted(
            3,
            Nucleotide::synthetic(b'A', Segment::V, NucFlags::empty()),
        )
        .with_base_changed(NucHandle::new(3), b't')
        .with_observation_deleted(0);
    changed
        .duplex_provenance
        .as_ref()
        .unwrap()
        .validate_pool(changed.pool.as_slice())
        .unwrap();
    let (replaced, _, _) = changed.with_segment_replaced(
        Segment::D,
        vec![Nucleotide::synthetic(b'G', Segment::D, NucFlags::empty())],
    );
    let position = replaced
        .sequence
        .regions
        .iter()
        .find(|r| r.segment == Segment::D)
        .unwrap()
        .start;
    let replaced = replaced.with_base_changed(position, b'N');
    replaced
        .duplex_provenance
        .as_ref()
        .unwrap()
        .validate_pool(replaced.pool.as_slice())
        .unwrap();
    let removed = replaced.with_indel_deleted(position.index());
    removed
        .duplex_provenance
        .as_ref()
        .unwrap()
        .validate_pool(removed.pool.as_slice())
        .unwrap();
}
