//! Persistent ancestry tests independent of the stochastic processing kernel.

use super::*;
use crate::address::PEnd;
use crate::airr_record::{build_airr_record, validate_airr_record, RecordValidationIssue};
use crate::assignment::SegmentOrientation;
use crate::ir::{flag, SimulationBuilder};

fn reference(segment: Segment, name: &str) -> Allele {
    Allele {
        name: name.into(),
        gene: name.into(),
        seq: b"AACGTC".to_vec(),
        segment,
        anchor: None,
        functional_status: None,
        subregions: Vec::new(),
    }
}

fn fixture() -> (RefDataConfig, Simulation) {
    let mut cfg = RefDataConfig::empty(ChainType::Vj);
    let original = cfg.v_pool.push(reference(Segment::V, "original*01"));
    let _ = cfg.v_pool.push(reference(Segment::V, "replacement*01"));
    let initial = Simulation::new()
        .with_allele_assigned(Segment::V, AlleleInstance::new(original).with_p_length_3(2));
    let mut builder = SimulationBuilder::from_simulation(initial);
    builder.record_coding_end_choice(Segment::V, original, 0, -2, 7);
    for (position, &base) in b"AACGTC".iter().enumerate() {
        builder.push_nucleotide(Nucleotide::germline(base, position as u16, Segment::V));
    }
    builder.add_region(Region::new(
        Segment::V,
        NucHandle::new(0),
        NucHandle::new(6),
    ));
    builder.push_coding_end_p(b'G', Segment::V, PEnd::V3, 5);
    builder.push_coding_end_p(b'A', Segment::V, PEnd::V3, 4);
    (cfg, builder.seal())
}

fn metadata_issues(sim: Simulation, cfg: &RefDataConfig) -> Vec<RecordValidationIssue> {
    let outcome = outcome_from_sim(sim);
    let record = build_airr_record(&outcome, cfg, "ancestry");
    validate_airr_record(&record, &outcome, cfg)
        .into_iter()
        .filter(|issue| {
            matches!(
                issue,
                RecordValidationIssue::CodingEndProcessingMismatch { .. }
            )
        })
        .collect()
}

#[test]
fn coding_end_provenance_survives_mask_indel_deletion_and_persistent_fork() {
    let (cfg, parent) = fixture();
    let child = parent
        .with_base_changed(NucHandle::new(6), b'N')
        .with_indel_inserted(
            0,
            Nucleotide::synthetic(b'T', Segment::V, flag::INDEL_INSERTED),
        );
    let record = build_airr_record(&outcome_from_sim(child.clone()), &cfg, "masked");
    let metadata = record.coding_end_processing.unwrap();
    assert_eq!(metadata["p_origins"][0]["query_position"], 7);
    assert_eq!(metadata["p_origins"][0]["current_base"], "N");
    assert_eq!(metadata["p_origins"][0]["original_base"], "G");
    assert!(child
        .pool
        .get(NucHandle::new(7))
        .unwrap()
        .flags
        .contains(flag::P_NUC));
    assert_eq!(
        parent.coding_end_provenance.as_ref().unwrap().p_origins[0].query_position,
        Some(6)
    );

    let clipped = child.with_indel_deleted(7);
    let outcome = outcome_from_sim(clipped.clone());
    let metadata = build_airr_record(&outcome, &cfg, "clipped")
        .coding_end_processing
        .unwrap();
    assert!(metadata["p_origins"][0]["query_position"].is_null());
    assert!(metadata["p_origins"][0]["current_base"].is_null());
    assert_eq!(metadata["p_origins"][1]["query_position"], 7);
    assert_eq!(metadata["retained_p_at_join"]["V_3"], 2);
    assert_eq!(metadata["observed_p"]["V_3"], 1);
    assert!(metadata_issues(clipped, &cfg).is_empty());
}

#[test]
fn coding_end_provenance_keeps_old_allele_after_receptor_replacement() {
    let (cfg, parent) = fixture();
    let revised = parent.with_allele_assigned(
        Segment::V,
        AlleleInstance::new(AlleleId::new(1)).with_p_length_3(2),
    );
    let (revised, _, _) = revised.with_segment_replaced(
        Segment::V,
        vec![
            Nucleotide::germline(b'A', 0, Segment::V),
            Nucleotide::germline(b'A', 1, Segment::V),
        ],
    );
    let metadata = build_airr_record(&outcome_from_sim(revised.clone()), &cfg, "revised")
        .coding_end_processing
        .unwrap();
    assert_eq!(metadata["choices"][0]["original_allele"], "original*01");
    assert_eq!(metadata["p_origins"][0]["original_allele"], "original*01");
    assert_eq!(metadata["p_origins"][0]["query_position"], 2);
    assert_eq!(metadata["p_origins"][1]["query_position"], 3);
    assert!(metadata_issues(revised, &cfg).is_empty());
}

#[test]
fn coding_end_provenance_distinguishes_canonical_and_assembled_inverted_d_ends() {
    let mut cfg = RefDataConfig::empty(ChainType::Vdj);
    let id = cfg.d_pool.push(reference(Segment::D, "d*01"));
    let initial = Simulation::new().with_allele_assigned(
        Segment::D,
        AlleleInstance::new(id).with_orientation(SegmentOrientation::ReverseComplement),
    );
    let mut builder = SimulationBuilder::from_simulation(initial);
    builder.record_coding_end_choice(Segment::D, id, -2, -1, 0);
    builder.push_coding_end_p(b'C', Segment::D, PEnd::D5, 5);
    builder.push_coding_end_p(b'A', Segment::D, PEnd::D3, 1);
    builder.push_coding_end_p(b'A', Segment::D, PEnd::D3, 0);
    let sim = builder.seal();
    let metadata = build_airr_record(&outcome_from_sim(sim.clone()), &cfg, "inverted")
        .coding_end_processing
        .unwrap();
    assert_eq!(metadata["choices"][0]["retained_p_5"], 2);
    assert_eq!(metadata["retained_p_at_join"]["D_5"], 1);
    assert_eq!(metadata["retained_p_at_join"]["D_3"], 2);
    assert_eq!(metadata["p_origins"][0]["canonical_end"], "3");
    assert_eq!(metadata["p_origins"][1]["canonical_end"], "5");
    assert!(metadata_issues(sim, &cfg).is_empty());
}

#[test]
fn coding_end_provenance_events_replay_choices_and_origin_positions() {
    let (_, fixture) = fixture();
    let initial = Simulation::new()
        .with_allele_assigned(Segment::V, *fixture.assignments.get(Segment::V).unwrap());
    let mut builder = SimulationBuilder::from_simulation(initial.clone());
    builder.attach_event_log_observer();
    builder.record_coding_end_choice(Segment::V, AlleleId::new(0), 0, -2, 7);
    builder.push_nucleotide(Nucleotide::germline(b'C', 5, Segment::V));
    builder.push_coding_end_p(b'G', Segment::V, PEnd::V3, 5);
    builder.push_coding_end_p(b'A', Segment::V, PEnd::V3, 4);
    builder.change_base(NucHandle::new(1), b'N');
    builder.insert_indel(
        0,
        Nucleotide::synthetic(b'T', Segment::V, flag::INDEL_INSERTED),
    );
    builder.delete_indel(3);
    let events = builder.seal_event_log_observer();
    let expected = builder.seal();
    let replayed = crate::ir::event_log_observer::replay_events(initial, &events);
    assert_eq!(replayed.pool.as_slice(), expected.pool.as_slice());
    assert_eq!(
        replayed.coding_end_provenance,
        expected.coding_end_provenance
    );
}

#[test]
fn coding_end_metadata_validator_detects_corrupt_count_and_original_tip() {
    let (cfg, mut sim) = fixture();
    let outcome = outcome_from_sim(sim.clone());
    let mut record = build_airr_record(&outcome, &cfg, "tampered");
    record.coding_end_processing.as_mut().unwrap()["observed_p"]["V_3"] = 99.into();
    assert!(validate_airr_record(&record, &outcome, &cfg).iter().any(|issue| {
        matches!(issue, RecordValidationIssue::CodingEndProcessingMismatch { field, .. } if field == "observed_p.V_3")
    }));
    std::sync::Arc::make_mut(sim.coding_end_provenance.as_mut().unwrap()).p_origins[0]
        .template_position = 2;
    assert!(!metadata_issues(sim, &cfg).is_empty());
}

#[test]
fn coding_end_provenance_complete_observation_loss_retains_history() {
    let (cfg, mut sim) = fixture();
    while !sim.pool.is_empty() {
        sim = sim.with_indel_deleted(0);
    }
    let metadata = build_airr_record(&outcome_from_sim(sim.clone()), &cfg, "empty")
        .coding_end_processing
        .unwrap();
    assert_eq!(metadata["retained_p_at_join"]["V_3"], 2);
    assert_eq!(metadata["observed_p"]["V_3"], 0);
    assert_eq!(metadata["p_origins"].as_array().unwrap().len(), 2);
    assert!(metadata_issues(sim, &cfg).is_empty());
}

#[test]
fn coding_end_validator_detects_missing_emission_even_with_matching_pool_and_counters() {
    let (cfg, parent) = fixture();
    let mut sim = parent
        .with_region_added(Region::new(
            Segment::J,
            NucHandle::new(8),
            NucHandle::new(8),
        ))
        .with_indel_deleted(6)
        .with_allele_assigned(
            Segment::V,
            AlleleInstance::new(AlleleId::new(0)).with_p_length_3(1),
        );
    // A fabricated ledger, pool and legacy count all agree on one P base,
    // while the immutable sampled outcome requires two original-tip bases.
    std::sync::Arc::make_mut(sim.coding_end_provenance.as_mut().unwrap())
        .p_origins
        .remove(0);
    assert!(metadata_issues(sim, &cfg).iter().any(|issue| {
        matches!(issue, RecordValidationIssue::CodingEndProcessingMismatch { reason, .. }
            if reason.contains("every selected"))
    }));
}

#[test]
fn coding_end_metadata_marks_partial_assembly_and_preserves_completion_after_clipping() {
    let (cfg, partial) = fixture();
    let record = build_airr_record(&outcome_from_sim(partial.clone()), &cfg, "partial");
    assert_eq!(
        record.coding_end_processing.unwrap()["recombination_complete"],
        false
    );
    let complete = partial.with_region_added(Region::new(
        Segment::J,
        NucHandle::new(8),
        NucHandle::new(8),
    ));
    let record = build_airr_record(&outcome_from_sim(complete.clone()), &cfg, "complete");
    assert_eq!(
        record.coding_end_processing.unwrap()["recombination_complete"],
        true
    );
    let missing_d =
        complete.with_allele_assigned(Segment::D, AlleleInstance::new(AlleleId::new(0)));
    let record = build_airr_record(&outcome_from_sim(missing_d), &cfg, "missing-d");
    assert_eq!(
        record.coding_end_processing.unwrap()["recombination_complete"],
        false
    );
    let mut clipped = complete;
    while !clipped.pool.is_empty() {
        clipped = clipped.with_indel_deleted(0);
    }
    let outcome = outcome_from_sim(clipped);
    let mut record = build_airr_record(&outcome, &cfg, "clipped");
    assert_eq!(
        record.coding_end_processing.as_ref().unwrap()["recombination_complete"],
        true
    );
    record.coding_end_processing.as_mut().unwrap()["recombination_complete"] = false.into();
    assert!(validate_airr_record(&record, &outcome, &cfg)
        .iter()
        .any(|issue| {
            matches!(issue, RecordValidationIssue::CodingEndProcessingMismatch { field, .. }
            if field == "recombination_complete")
        }));
}
