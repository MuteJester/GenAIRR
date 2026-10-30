//! Literal descendants of a chemically restored V anchor. These tests do not
//! infer placement from the final sequence or from the nucleotide cache.
use super::*;
use crate::contract::{productive, AnchorPreserved, CodingSequenceIntegrity, Contract};
use crate::duplex_repair::bindings::DuplexReferenceBindings;
use crate::duplex_repair::{
    repair_coupled_annotated, JoinName, JoinScript, RepairOperation, Side, Strand,
};
use crate::ir::{flag, NucFlags};
use crate::refdata::{Allele, ChainType};

fn restored_fixture() -> (RefDataConfig, Simulation) {
    restored_fixture_with_reference("AAATGT", false, None)
}

fn restored_fixture_with_reference(
    v: &str,
    verified: bool,
    phase: Option<u8>,
) -> (RefDataConfig, Simulation) {
    let mut cfg = RefDataConfig::empty(ChainType::Vdj);
    for (segment, name, sequence, anchor) in [
        (Segment::V, "V*01", v, Some(3)),
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
        match segment {
            Segment::V => cfg.v_pool.push(allele),
            Segment::D => cfg.d_pool.push(allele),
            Segment::J => cfg.j_pool.push(allele),
            _ => unreachable!(),
        };
    }
    if verified {
        use crate::refdata::{AlleleEvidence, CodingBoundaryStatus};
        cfg.set_allele_evidence(
            Segment::V,
            AlleleId::new(0),
            AlleleEvidence {
                schema_version: 1,
                sequence_sha256: crate::refdata::evidence::sequence_sha256(v.as_bytes()),
                coding_phase: phase,
                anchor_status: AnchorEvidenceStatus::Verified,
                boundary_5: CodingBoundaryStatus::Complete,
                boundary_3: CodingBoundaryStatus::Complete,
                source: Some("independently positioned fixture landmark".into()),
                method: "authored_coordinate".into(),
                notes: vec![],
            },
        )
        .unwrap();
    }
    let joins = [
        JoinScript {
            join: JoinName::VD,
            operations: vec![
                RepairOperation::Open {
                    side: Side::Left,
                    offset: 0,
                },
                RepairOperation::Open {
                    side: Side::Right,
                    offset: 0,
                },
                RepairOperation::Resect {
                    side: Side::Left,
                    strand: Strand::Top,
                    count: 1,
                },
                RepairOperation::Synapse { right_start: 6 },
                RepairOperation::Fill {
                    strand: Strand::Top,
                },
                RepairOperation::Ligate {},
            ],
        },
        JoinScript {
            join: JoinName::DJ,
            operations: vec![
                RepairOperation::Open {
                    side: Side::Left,
                    offset: 0,
                },
                RepairOperation::Open {
                    side: Side::Right,
                    offset: 0,
                },
                RepairOperation::Synapse { right_start: 3 },
                RepairOperation::Ligate {},
            ],
        },
    ];
    let annotated = repair_coupled_annotated(v, "CCC", "TTTAAA", false, &joins, 256, 1000).unwrap();
    assert_eq!(annotated.product.sequence, format!("{v}CCCTTTAAA"));
    assert_eq!(annotated.product.columns[5].top.kind, "templated_synthesis");
    let bindings = DuplexReferenceBindings::new(
        &cfg,
        ["V*01", "D*01", "J*01"],
        [v, "CCC", "TTTAAA"],
        false,
        true,
    )
    .unwrap();
    let sim = crate::duplex_repair::seed::seed_from_coupled(annotated, bindings, &cfg).unwrap();
    (cfg, sim)
}

fn without_cached_restoration(sim: &Simulation) -> Simulation {
    // Deliberately remove the logical cache entry. The certified paired-column
    // placement, not a synthetic nucleotide's apparent origin, is authoritative.
    sim.with_nucleotide_changed(
        NucHandle::new(5),
        Nucleotide::synthetic(b'T', Segment::V, NucFlags::empty()),
    )
}

fn clip(
    sim: Simulation,
    cfg: &RefDataConfig,
    end: crate::passes::LossEnd,
    count: i64,
) -> Simulation {
    use crate::compiled::{CompileOptions, CompiledSimulator, ExecutionPolicy};
    use crate::dist::EmpiricalLengthDist;
    use crate::pass::PassPlan;
    use crate::passes::EndLossPass;
    let mut plan = PassPlan::new();
    plan.push(Box::new(EndLossPass::new(
        end,
        Box::new(EmpiricalLengthDist::from_pairs(vec![(count, 1.0)])),
    )));
    let prior_removals = sim.duplex_provenance.as_ref().unwrap().removed.len();
    let compiled = CompiledSimulator::compile_with_options(
        &plan,
        Some(cfg),
        None,
        ExecutionPolicy::Strict,
        CompileOptions::allow_curatable_refdata(),
    )
    .unwrap();
    let outcome = compiled.run_one_from(sim, 0).unwrap();
    let clipped = outcome.final_simulation().clone();
    assert!(clipped
        .duplex_provenance
        .as_ref()
        .unwrap()
        .removed
        .iter()
        .skip(prior_removals)
        .all(|removed| removed.cause == crate::ir::DuplexRemovalCause::ObservationLoss,));
    clipped
}

#[test]
fn restored_anchor_uses_certificates_and_keeps_legacy_policy_unchanged() {
    let (cfg, initial) = restored_fixture();
    let sim = without_cached_restoration(&initial);
    let a = assess(&sim, &cfg);
    assert_eq!(policy_name(&sim), DUPLEX_POLICY);
    assert_eq!(a.productive, Some(true));
    assert_eq!(a.junction, Some((3, 12)));
    assert_eq!(a.stop_codon, Some(false));
    assert!(productive().verify(&sim, Some(&cfg)).is_ok());
    let mut legacy = sim.clone();
    legacy.duplex_provenance = None;
    assert_eq!(policy_name(&legacy), POLICY);
    assert_eq!(assess(&legacy, &cfg).productive, Some(false));
    assert!(assess(&legacy, &cfg)
        .reasons
        .contains(&"anchor_provenance_disrupted"));
}

#[test]
fn restored_anchor_substitutions_are_constrained_at_the_certified_position() {
    let (cfg, initial) = restored_fixture();
    let sim = without_cached_restoration(&initial);
    let anchor = AnchorPreserved::canonical(Segment::V);
    // Restored TGT's final T can become C (Cys), but A would form TGA.
    assert!(anchor.admits_fixed_base_at(&sim, Some(&cfg), NucHandle::new(5), b'C'));
    assert!(!anchor.admits_fixed_base_at(&sim, Some(&cfg), NucHandle::new(5), b'A'));
    assert_eq!(
        assess(&sim.with_base_changed(NucHandle::new(5), b'C'), &cfg).productive,
        Some(true)
    );
    let stopped = assess(&sim.with_base_changed(NucHandle::new(5), b'A'), &cfg);
    assert_eq!(stopped.productive, Some(false));
    assert_eq!(stopped.stop_codon, Some(true));
    // AAA before the anchor remains in the same measured coding span.
    assert!(!CodingSequenceIntegrity.admits_fixed_base_at(
        &sim,
        Some(&cfg),
        NucHandle::new(0),
        b'T',
    ));
}

#[test]
fn verified_noncanonical_source_does_not_unconstrain_a_now_canonical_live_anchor() {
    let (cfg, original) = restored_fixture_with_reference("AAAGGG", true, Some(0));
    assert_eq!(assess(&original, &cfg).productive, Some(false));
    let corrected = original
        .with_base_changed(NucHandle::new(3), b'T')
        .with_base_changed(NucHandle::new(5), b'T');
    assert_eq!(assess(&corrected, &cfg).productive, Some(true));
    assert!(productive().verify(&corrected, Some(&cfg)).is_ok());
    let contract = AnchorPreserved::canonical(Segment::V);
    assert!(!contract.admits_fixed_base_at(&corrected, Some(&cfg), NucHandle::new(3), b'A'));
    assert!(!contract.admits_fixed_base_at(&corrected, Some(&cfg), NucHandle::new(5), b'A'));
    assert!(contract.admits_fixed_base_at(&corrected, Some(&cfg), NucHandle::new(5), b'C'));
    assert_ne!(
        contract.admissible_bases_at(&corrected, Some(&cfg), NucHandle::new(3)),
        crate::contract::BaseMask::UNCONSTRAINED
    );
}

#[test]
fn known_rejected_restored_residue_remains_false_without_an_invented_phase() {
    let (cfg, sim) = restored_fixture_with_reference("AAAGGG", true, None);
    let a = assess(&sim, &cfg);
    assert_eq!(a.junction, Some((3, 12)));
    assert_eq!(a.productive, Some(false));
    assert_eq!(a.vj_in_frame, None);
    assert_eq!(a.stop_codon, None);
    assert!(a.reasons.contains(&"anchor_residue_rejected"));
    assert!(a.reasons.contains(&"coding_phase_unresolved"));
    let corrected = sim
        .with_base_changed(NucHandle::new(3), b'T')
        .with_base_changed(NucHandle::new(5), b'T');
    assert_eq!(assess(&corrected, &cfg).productive, None);
    assert_eq!(assess(&corrected, &cfg).stop_codon, None);
}

#[test]
fn a_new_matching_base_cannot_recreate_the_removed_founder_landmark() {
    let (cfg, sim) = restored_fixture();
    let replacement = sim
        .with_indel_deleted(5)
        .with_indel_inserted(5, Nucleotide::germline(b'T', 5, Segment::V));
    assert_eq!(
        replacement
            .pool
            .as_slice()
            .iter()
            .map(|n| n.base)
            .collect::<Vec<_>>(),
        sim.pool
            .as_slice()
            .iter()
            .map(|n| n.base)
            .collect::<Vec<_>>(),
    );
    // Even a misleading reference-cache entry cannot grant body placement.
    assert_eq!(replacement.body_reference_position(Segment::V, 5), None);
    assert_eq!(assess(&replacement, &cfg).productive, Some(false));
    assert!(productive().verify(&replacement, Some(&cfg)).is_err());
}

#[test]
fn first_v_deletion_and_read_loss_preserve_different_coding_origins() {
    let (cfg, sim) = restored_fixture();
    let deleted = sim.with_indel_deleted(0);
    assert_eq!(assess(&deleted, &cfg).vj_in_frame, Some(false));
    assert_eq!(assess(&deleted, &cfg).productive, Some(false));
    let clipped = clip(sim, &cfg, crate::passes::LossEnd::Five, 1);
    assert_eq!(clipped.pool.as_slice(), deleted.pool.as_slice());
    assert_eq!(
        clipped.assignment(Segment::V).unwrap().observation_loss_5,
        1
    );
    assert_eq!(assess(&clipped, &cfg).vj_in_frame, Some(true));
    assert_eq!(assess(&clipped, &cfg).productive, Some(true));
    assert_eq!(policy_name(&clipped), DUPLEX_POLICY);
}

#[test]
fn observing_only_part_of_the_restored_anchor_yields_unknown() {
    let (cfg, sim) = restored_fixture();
    let clipped = clip(sim, &cfg, crate::passes::LossEnd::Five, 4);
    let a = assess(&clipped, &cfg);
    assert_eq!(a.productive, None);
    assert_eq!(a.junction, None);
    assert!(a.reasons.contains(&"anchor_not_observed"));
}

#[test]
fn deleting_all_d_body_does_not_invent_a_missing_vj_landmark() {
    let (cfg, mut sim) = restored_fixture();
    for _ in 0..3 {
        sim = sim.with_indel_deleted(6);
    }
    let sequence: Vec<_> = sim.pool.as_slice().iter().map(|n| n.base).collect();
    assert_eq!(sequence, b"AAATGTTTTAAA");
    assert_eq!(sim.region(Segment::D), Some((6, 6)));
    let a = assess(&sim, &cfg);
    assert_eq!(a.junction, Some((3, 9)));
    assert_eq!(a.vj_in_frame, Some(true));
    assert_eq!(a.productive, Some(true));
}

#[test]
fn j_terminal_deletion_is_not_reinterpreted_as_observation_loss() {
    let (cfg, sim) = restored_fixture();
    let deleted = sim.with_indel_deleted(14);
    assert_eq!(assess(&deleted, &cfg).productive, Some(false));
    let clipped = clip(sim, &cfg, crate::passes::LossEnd::Three, 1);
    assert_eq!(clipped.pool.as_slice(), deleted.pool.as_slice());
    assert_eq!(
        clipped.assignment(Segment::J).unwrap().observation_loss_3,
        1
    );
    assert_eq!(assess(&clipped, &cfg).productive, Some(true));
}

#[test]
fn clipping_an_inserted_j_suffix_restores_the_observed_reference_endpoint() {
    let (cfg, sim) = restored_fixture();
    let inserted = sim.with_indel_inserted(
        15,
        Nucleotide::synthetic(b'A', Segment::J, flag::INDEL_INSERTED),
    );
    assert_eq!(assess(&inserted, &cfg).productive, Some(false));
    let clipped = clip(inserted, &cfg, crate::passes::LossEnd::Three, 1);
    assert_eq!(clipped.pool.as_slice(), sim.pool.as_slice());
    assert_eq!(clipped.j_reference_end_for_phase(6), Some(6));
    assert_eq!(assess(&clipped, &cfg).productive, Some(true));
}

#[test]
fn clipping_past_a_prior_j_deletion_limits_assessment_to_observed_reference_coverage() {
    let (cfg, sim) = restored_fixture();
    let deleted = sim.with_indel_deleted(14);
    assert_eq!(assess(&deleted, &cfg).productive, Some(false));
    let clipped = clip(deleted, &cfg, crate::passes::LossEnd::Three, 1);
    let ordinary_clip = clip(sim, &cfg, crate::passes::LossEnd::Three, 2);
    assert_eq!(clipped.pool.as_slice(), ordinary_clip.pool.as_slice());
    assert_eq!(
        clipped.assignment(Segment::J).unwrap().observed_ref_end,
        Some(4)
    );
    assert_eq!(assess(&clipped, &cfg).productive, Some(true));
}

#[test]
fn temporary_internal_insertions_only_accept_a_fully_restored_codon() {
    let (cfg, sim) = restored_fixture();
    let inserted = sim.with_indel_inserted(
        4,
        Nucleotide::synthetic(b'C', Segment::V, flag::INDEL_INSERTED),
    );
    assert_eq!(assess(&inserted, &cfg).productive, Some(false));
    let restored = inserted.with_indel_deleted(4);
    assert_eq!(assess(&restored, &cfg).productive, Some(true));
    assert_eq!(restored.body_reference_position(Segment::V, 5), Some(5));
}

#[test]
fn legacy_anchor_preservation_is_explicitly_unsupported_for_duplex() {
    let (cfg, sim) = restored_fixture();
    let contract = AnchorPreserved::new(Segment::V);
    let violation = contract.verify(&sim, Some(&cfg)).unwrap_err();
    assert!(violation.reason.contains("unsupported for duplex ancestry"));
    assert!(!contract.admits_fixed_base_at(&sim, Some(&cfg), NucHandle::new(0), b'A'));
    assert_eq!(
        contract.admissible_bases_at(&sim, Some(&cfg), NucHandle::new(0)),
        crate::contract::BaseMask::EMPTY
    );
    assert!(AnchorPreserved::canonical(Segment::V)
        .verify(&sim, Some(&cfg))
        .is_ok());
}

#[test]
fn exact_indel_support_preserves_the_restored_duplex_landmark() {
    use crate::compiled::{CompileOptions, CompiledSimulator, ExecutionPolicy};
    use crate::dist::{EmpiricalLengthDist, UniformBase};
    use crate::pass::PassPlan;
    use crate::passes::corrupt::indel::IndelSamplingOptions;
    use crate::passes::IndelPass;
    let (cfg, sim) = restored_fixture();
    let mut plan = PassPlan::new();
    plan.push(Box::new(
        IndelPass::new(
            Box::new(EmpiricalLengthDist::from_pairs(vec![(2, 1.0)])),
            0.5,
            Box::new(UniformBase),
        )
        .with_sampling_options(IndelSamplingOptions {
            rejection_attempts: 0,
            exact_node_budget: 100_000,
            legacy_replay_sentinels: false,
        }),
    ));
    let contracts = productive();
    let compiled = CompiledSimulator::compile_with_options(
        &plan,
        Some(&cfg),
        Some(&contracts),
        ExecutionPolicy::Strict,
        CompileOptions::allow_curatable_refdata(),
    )
    .unwrap();
    let outcome = compiled
        .run_one_from(sim, 17)
        .expect("complete duplex-preserving two-event histories exist");
    let final_state = outcome.final_simulation();
    assert_eq!(policy_name(final_state), DUPLEX_POLICY);
    assert_eq!(assess(final_state, &cfg).productive, Some(true));
    assert!(productive().verify(final_state, Some(&cfg)).is_ok());
}

#[test]
fn nullable_record_oracle_uses_certificates_and_observed_reference_end_independently() {
    let (cfg, sim) = restored_fixture();
    let inserted = sim.with_indel_inserted(
        15,
        Nucleotide::synthetic(b'A', Segment::J, flag::INDEL_INSERTED),
    );
    let clipped = clip(inserted, &cfg, crate::passes::LossEnd::Three, 1);
    for state in [&sim, &clipped] {
        let mut record = crate::duplex_repair::observation::project_duplex_simulation(state, &cfg)
            .unwrap()["record"]
            .clone();
        assert!(crate::airr_record::validate_duplex_productivity(&record, state, &cfg).is_empty());
        // The independent oracle must not depend on an apparent original
        // identity placed in the nucleotide cache.
        let no_cache = without_cached_restoration(state);
        assert!(
            crate::airr_record::validate_duplex_productivity(&record, &no_cache, &cfg).is_empty()
        );
        record["productive"] = serde_json::json!(false);
        let issues = crate::airr_record::validate_duplex_productivity(&record, state, &cfg);
        assert_eq!(issues.len(), 1);
    }
}

#[test]
fn coding_integrity_mask_agrees_with_asking_about_each_base_separately() {
    // `admissible_bases_at` answers for all four bases at once and
    // `admits_fixed_base_at` for one; they have to give the same answer.
    //
    // Worth stating what this test was and became: before the per-base loop in
    // `admissible_bases_at` was hoisted, both hooks called one shared
    // `allows_base` and this assertion was a tautology. The hoisted version
    // computes the coding frame and the two unchanged codon bytes once, so the
    // two hooks are now genuinely different code and this is the guard that
    // they did not diverge -- including on the three early-exit paths, where
    // the mask has to be "all four admitted" rather than "none".
    let (cfg, initial) = restored_fixture();
    let sim = without_cached_restoration(&initial);
    let mut sites_inside_a_codon = 0;
    for index in 0..sim.pool.len() {
        let site = NucHandle::new(index as u32);
        let combined = CodingSequenceIntegrity.admissible_bases_at(&sim, Some(&cfg), site);
        let mut separate = 0u8;
        for bit in 0..4 {
            if CodingSequenceIntegrity.admits_fixed_base_at(
                &sim,
                Some(&cfg),
                site,
                crate::contract::BaseMask::base_for_bit(bit),
            ) {
                separate |= 1 << bit;
            }
        }
        assert_eq!(
            combined.0, separate,
            "site {index}: mask {:#06b} from admissible_bases_at, {separate:#06b} from \
             four admits_fixed_base_at calls",
            combined.0
        );
        if combined.0 != 0b1111 {
            sites_inside_a_codon += 1;
        }
    }
    assert!(
        sites_inside_a_codon > 0,
        "every site admitted all four bases, so this fixture exercises only the \
         early-exit paths and the codon comparison is untested"
    );
}
