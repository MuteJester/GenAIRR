//! Independently specified regressions from the biological metadata audit.

use super::*;
use crate::airr_record::build_airr_record;
use crate::assignment::SegmentOrientation;
use crate::ir::{complement_base, flag};
use crate::live_call::{
    AlleleBitSet, EvidenceScore, HypothesisFlags, PlacementHypothesis, SegmentCalls,
    SegmentLiveCall,
};

fn allele(segment: Segment, name: &str, sequence: &[u8]) -> Allele {
    Allele {
        name: name.into(),
        gene: name.into(),
        seq: sequence.to_vec(),
        segment,
        anchor: None,
        functional_status: None,
        subregions: Vec::new(),
    }
}

fn push(sim: Simulation, nucleotide: Nucleotide) -> Simulation {
    sim.with_nucleotide_pushed(nucleotide).0
}

fn region(sim: Simulation, segment: Segment, start: u32, end: u32) -> Simulation {
    sim.with_region_added(Region::new(
        segment,
        NucHandle::new(start),
        NucHandle::new(end),
    ))
}

#[test]
fn p_additions_remain_between_their_regions_in_alignment_pool_order() {
    let mut cfg = RefDataConfig::empty(ChainType::Vj);
    let v = cfg.v_pool.push(allele(Segment::V, "v*01", b"AAA"));
    let j = cfg.j_pool.push(allele(Segment::J, "j*01", b"TTT"));
    let mut sim = Simulation::new()
        .with_allele_assigned(Segment::V, AlleleInstance::new(v))
        .with_allele_assigned(Segment::J, AlleleInstance::new(j));
    for position in 0..3 {
        sim = push(sim, Nucleotide::germline(b'A', position, Segment::V));
    }
    sim = region(sim, Segment::V, 0, 3);
    sim = push(sim, Nucleotide::synthetic(b'T', Segment::Np1, flag::P_NUC));
    for _ in 0..2 {
        sim = push(sim, Nucleotide::synthetic(b'C', Segment::Np1, flag::N_NUC));
    }
    sim = region(sim, Segment::Np1, 4, 6);
    sim = push(sim, Nucleotide::synthetic(b'A', Segment::Np1, flag::P_NUC));
    for position in 0..3 {
        sim = push(sim, Nucleotide::germline(b'T', position, Segment::J));
    }
    sim = region(sim, Segment::J, 7, 10);

    let rec = build_airr_record(&outcome_from_sim(sim), &cfg, "p-order");
    assert_eq!(rec.sequence, "AAATCCATTT");
    assert_eq!(rec.sequence_alignment, "AAATCCATTT");
    assert_eq!(rec.germline_alignment, "AAANNNNTTT");
    assert_eq!(rec.j_alignment_start, Some(7));
    assert_eq!(rec.j_alignment_end, Some(10));
    assert_eq!(rec.j_sequence_start, Some(7));
    assert_eq!(rec.j_cigar, "3M");
}

fn inverted_d_sim(cfg: &mut RefDataConfig, trim_5: u16, trim_3: u16) -> (Simulation, AlleleId) {
    let reference = b"AACGTA";
    let d = cfg.d_pool.push(allele(Segment::D, "d*01", reference));
    let mut sim = Simulation::new().with_allele_assigned(
        Segment::D,
        AlleleInstance::new(d)
            .with_trim_5(trim_5)
            .with_trim_3(trim_3)
            .with_orientation(SegmentOrientation::ReverseComplement),
    );
    for position in (trim_5 as usize..reference.len() - trim_3 as usize).rev() {
        sim = push(
            sim,
            Nucleotide::germline(
                complement_base(reference[position]),
                position as u16,
                Segment::D,
            )
            .with_flags(flag::INVERTED),
        );
    }
    let end = sim.pool.len() as u32;
    (region(sim, Segment::D, 0, end), d)
}

#[test]
fn inverted_d_projects_the_oriented_reference_without_fictitious_deletions() {
    let mut cfg = RefDataConfig::empty(ChainType::Vdj);
    let (sim, _) = inverted_d_sim(&mut cfg, 0, 0);
    let rec = build_airr_record(&outcome_from_sim(sim), &cfg, "inverted");
    assert_eq!(rec.sequence, "TACGTT");
    assert_eq!(rec.sequence_alignment, "TACGTT");
    assert_eq!(rec.germline_alignment, "TACGTT");
    assert_eq!(rec.germline_alignment_d_mask, "NNNNNN");
    assert_eq!(rec.d_cigar, "6M");
    assert_eq!(rec.d_identity, Some(1.0));
    assert_eq!(rec.d_germline_start, Some(0));
    assert_eq!(rec.d_germline_end, Some(6));
}

#[test]
fn inverted_d_with_asymmetric_trims_keeps_original_allele_coordinates() {
    let mut cfg = RefDataConfig::empty(ChainType::Vdj);
    let (sim, _) = inverted_d_sim(&mut cfg, 1, 2);
    // Retained canonical reference [1,4) is ACG, oriented query CGT.
    let rec = build_airr_record(&outcome_from_sim(sim), &cfg, "inverted-trim");
    assert_eq!(rec.sequence_alignment, "CGT");
    assert_eq!(rec.germline_alignment, "CGT");
    assert_eq!(rec.d_cigar, "3M");
    assert_eq!(rec.d_identity, Some(1.0));
    assert_eq!(rec.d_germline_start, Some(1));
    assert_eq!(rec.d_germline_end, Some(4));
}

#[test]
fn inverted_d_internal_indels_are_rendered_in_query_order() {
    let mut cfg = RefDataConfig::empty(ChainType::Vdj);
    let d = cfg.d_pool.push(allele(Segment::D, "d*01", b"AACGTA"));
    let mut sim = Simulation::new().with_allele_assigned(
        Segment::D,
        AlleleInstance::new(d).with_orientation(SegmentOrientation::ReverseComplement),
    );
    // Oriented reference TACGTT: insert A after its first T and delete G.
    for (base, position) in [
        (b'T', Some(5)),
        (b'A', None),
        (b'A', Some(4)),
        (b'C', Some(3)),
        (b'T', Some(1)),
        (b'T', Some(0)),
    ] {
        sim = push(
            sim,
            match position {
                Some(position) => {
                    Nucleotide::germline(base, position, Segment::D).with_flags(flag::INVERTED)
                }
                None => Nucleotide::synthetic(base, Segment::D, flag::INDEL_INSERTED),
            },
        );
    }
    sim = region(sim, Segment::D, 0, 6);
    let rec = build_airr_record(&outcome_from_sim(sim), &cfg, "inverted-indels");
    assert_eq!(rec.sequence_alignment, "TAAC-TT");
    assert_eq!(rec.germline_alignment, "T-ACGTT");
    assert_eq!(rec.d_cigar, "1M1I2M1D2M");
    assert_eq!(rec.d_identity, Some(5.0 / 7.0));
    assert_eq!(rec.d_germline_start, Some(0));
    assert_eq!(rec.d_germline_end, Some(6));
}

#[test]
fn inverted_d_recall_extensions_do_not_reassign_generated_np_bases() {
    let mut cfg = RefDataConfig::empty(ChainType::Vdj);
    let d = cfg.d_pool.push(allele(Segment::D, "d*01", b"AACGTA"));
    let mut sim = Simulation::new().with_allele_assigned(
        Segment::D,
        AlleleInstance::new(d)
            .with_trim_5(1)
            .with_trim_3(1)
            .with_orientation(SegmentOrientation::ReverseComplement),
    );
    sim = push(sim, Nucleotide::synthetic(b'T', Segment::Np1, flag::N_NUC));
    sim = region(sim, Segment::Np1, 0, 1);
    for (base, position) in [(b'A', 4), (b'C', 3), (b'G', 2), (b'T', 1)] {
        sim = push(sim, Nucleotide::germline(base, position, Segment::D));
    }
    sim = region(sim, Segment::D, 1, 5);
    sim = push(sim, Nucleotide::synthetic(b'T', Segment::Np2, flag::N_NUC));
    sim = region(sim, Segment::Np2, 5, 6);
    let hypothesis = PlacementHypothesis::new(
        Segment::D,
        0,
        6,
        0,
        6,
        AlleleBitSet::from_ids(cfg.d_pool.len(), [d]),
        EvidenceScore::exact(6, 0),
        HypothesisFlags::EMPTY,
    );
    let call = SegmentLiveCall::from_hypotheses(Segment::D, cfg.d_pool.len(), vec![hypothesis], 1);
    sim = sim.with_segment_calls(SegmentCalls::empty().with_segment_call(call));
    let rec = build_airr_record(&outcome_from_sim(sim), &cfg, "inverted-extensions");
    assert_eq!(rec.sequence_alignment, "TACGTT");
    assert_eq!(rec.germline_alignment, "NACGTN");
    assert_eq!(rec.d_cigar, "4M");
    assert_eq!(rec.d_identity, Some(1.0));
    assert_eq!(rec.d_sequence_start, Some(1));
    assert_eq!(rec.d_sequence_end, Some(5));
    assert_eq!(rec.d_germline_start, Some(1));
    assert_eq!(rec.d_germline_end, Some(5));
}

#[test]
fn shorter_recalled_allele_does_not_replace_generating_alignment_reference() {
    let mut cfg = RefDataConfig::empty(ChainType::Vj);
    let original = cfg
        .v_pool
        .push(allele(Segment::V, "original*01", b"AAACCCGGG"));
    let called = cfg.v_pool.push(allele(Segment::V, "short*01", b"TTTAAA"));
    let mut sim = Simulation::new().with_allele_assigned(Segment::V, AlleleInstance::new(original));
    for (position, &base) in b"TTTAAAGGG".iter().enumerate() {
        let mut nucleotide =
            Nucleotide::germline(b"AAACCCGGG"[position], position as u16, Segment::V);
        nucleotide.base = base;
        sim = push(sim, nucleotide);
    }
    sim = region(sim, Segment::V, 0, 9);
    let hypothesis = PlacementHypothesis::new(
        Segment::V,
        0,
        9,
        0,
        9,
        AlleleBitSet::from_ids(cfg.v_pool.len(), [called]),
        EvidenceScore::exact(6, 0),
        HypothesisFlags::EMPTY,
    );
    let call = SegmentLiveCall::from_hypotheses(Segment::V, cfg.v_pool.len(), vec![hypothesis], 1);
    sim = sim.with_segment_calls(SegmentCalls::empty().with_segment_call(call));
    let rec = build_airr_record(&outcome_from_sim(sim.clone()), &cfg, "short-reference");
    assert_eq!(rec.v_call, "original*01");
    assert_eq!(rec.recalled_v_call, "short*01");
    assert_eq!(rec.sequence_alignment, "TTTAAAGGG");
    assert_eq!(rec.germline_alignment, "AAACCCGGG");
    assert_eq!(rec.v_cigar, "9M");
    assert_eq!(rec.v_germline_start, Some(0));
    assert_eq!(rec.v_germline_end, Some(9));
    assert_eq!(rec.v_identity, Some(3.0 / 9.0));

    // A terminal source deletion beyond the six-base re-called allele
    // must remain visible, together with the independent internal deletion.
    let deleted = sim.with_indel_deleted(8).with_indel_deleted(3);
    let rec = build_airr_record(
        &outcome_from_sim(deleted),
        &cfg,
        "short-recall-with-deletions",
    );
    assert_eq!(rec.v_call, "original*01");
    assert_eq!(rec.recalled_v_call, "short*01");
    assert_eq!(rec.sequence_alignment, "TTT-AAGG-");
    assert_eq!(rec.germline_alignment, "AAACCCGGG");
    assert_eq!(rec.v_cigar, "3M1D4M1D");
    assert_eq!(rec.v_germline_end, Some(9));
    assert_eq!(rec.v_identity, Some(2.0 / 9.0));
}

#[test]
fn longer_recalled_allele_does_not_invent_deletion_of_an_unsampled_tail() {
    let mut cfg = RefDataConfig::empty(ChainType::Vj);
    let original = cfg
        .v_pool
        .push(allele(Segment::V, "original*01", b"AAACCC"));
    let called = cfg.v_pool.push(allele(Segment::V, "long*01", b"TTTAAAGGG"));
    let mut sim = Simulation::new().with_allele_assigned(Segment::V, AlleleInstance::new(original));
    for (position, &base) in b"TTTAAA".iter().enumerate() {
        let mut nucleotide = Nucleotide::germline(b"AAACCC"[position], position as u16, Segment::V);
        nucleotide.base = base;
        sim = push(sim, nucleotide);
    }
    sim = region(sim, Segment::V, 0, 6);
    let hypothesis = PlacementHypothesis::new(
        Segment::V,
        0,
        6,
        0,
        6,
        AlleleBitSet::from_ids(cfg.v_pool.len(), [called]),
        EvidenceScore::exact(6, 0),
        HypothesisFlags::EMPTY,
    );
    let call = SegmentLiveCall::from_hypotheses(Segment::V, cfg.v_pool.len(), vec![hypothesis], 1);
    sim = sim.with_segment_calls(SegmentCalls::empty().with_segment_call(call));
    let rec = build_airr_record(&outcome_from_sim(sim), &cfg, "long-reference");
    assert_eq!(rec.v_call, "original*01");
    assert_eq!(rec.recalled_v_call, "long*01");
    assert_eq!(rec.sequence_alignment, "TTTAAA");
    assert_eq!(rec.germline_alignment, "AAACCC");
    assert_eq!(rec.v_cigar, "6M");
    assert_eq!(rec.v_germline_start, Some(0));
    assert_eq!(rec.v_germline_end, Some(6));
    assert_eq!(rec.v_identity, Some(0.0));
}

#[test]
fn np_base_matching_both_trimmed_germline_neighbors_retains_its_generated_origin() {
    let mut cfg = RefDataConfig::empty(ChainType::Vj);
    let v = cfg.v_pool.push(allele(Segment::V, "v*01", b"AAAC"));
    let j = cfg.j_pool.push(allele(Segment::J, "j*01", b"CGGG"));
    let mut sim = Simulation::new()
        .with_allele_assigned(Segment::V, AlleleInstance::new(v).with_trim_3(1))
        .with_allele_assigned(Segment::J, AlleleInstance::new(j).with_trim_5(1));
    for position in 0..3 {
        sim = push(sim, Nucleotide::germline(b'A', position, Segment::V));
    }
    sim = region(sim, Segment::V, 0, 3);
    sim = push(sim, Nucleotide::synthetic(b'C', Segment::Np1, flag::N_NUC));
    sim = region(sim, Segment::Np1, 3, 4);
    for position in 1..4 {
        sim = push(sim, Nucleotide::germline(b'G', position, Segment::J));
    }
    sim = region(sim, Segment::J, 4, 7);
    let mut calls = SegmentCalls::empty();
    for (segment, id, start, end) in [(Segment::V, v, 0, 4), (Segment::J, j, 3, 7)] {
        let hypothesis = PlacementHypothesis::new(
            segment,
            start,
            end,
            0,
            4,
            AlleleBitSet::from_ids(1, [id]),
            EvidenceScore::exact(4, 0),
            HypothesisFlags::EMPTY,
        );
        calls = calls.with_segment_call(SegmentLiveCall::from_hypotheses(
            segment,
            1,
            vec![hypothesis],
            1,
        ));
    }
    sim = sim.with_segment_calls(calls);
    let outcome = outcome_from_sim(sim);
    let rec = build_airr_record(&outcome, &cfg, "np-truth");
    assert_eq!(rec.sequence_alignment, "AAACGGG");
    assert_eq!(rec.germline_alignment, "AAANGGG");
    assert_eq!(rec.v_cigar, "3M");
    assert_eq!(rec.j_cigar, "3M");
    assert_eq!(
        (rec.v_sequence_end, rec.j_sequence_start),
        (Some(3), Some(4))
    );
    assert_eq!(rec.np1, "C");
    assert_eq!(rec.np1_length, 1);
    assert_eq!((rec.v_identity, rec.j_identity), (Some(1.0), Some(1.0)));
    // Record projection does not remove the evidence-based inference.
    let live = &outcome.final_simulation().segment_calls;
    assert_eq!(live.get(Segment::V).unwrap().hypotheses[0].seq_end, 4);
    assert_eq!(live.get(Segment::J).unwrap().hypotheses[0].seq_start, 3);
}

fn lose_end(sim: Simulation, end: crate::passes::LossEnd, count: i64) -> Outcome {
    use crate::dist::EmpiricalLengthDist;
    use crate::pass::{testing::PassRuntime, PassPlan};
    use crate::passes::EndLossPass;
    let mut plan = PassPlan::new();
    plan.push(Box::new(EndLossPass::new(
        end,
        Box::new(EmpiricalLengthDist::from_pairs(vec![(count, 1.0)])),
    )));
    PassRuntime::execute(&plan, sim, 0)
}

fn full_v_for_coverage() -> (RefDataConfig, Simulation) {
    let mut cfg = RefDataConfig::empty(ChainType::Vj);
    let id = cfg.v_pool.push(allele(Segment::V, "v*01", b"AACCGGTT"));
    let mut sim = Simulation::new().with_allele_assigned(Segment::V, AlleleInstance::new(id));
    for (position, &base) in b"AACCGGTT".iter().enumerate() {
        sim = push(sim, Nucleotide::germline(base, position as u16, Segment::V));
    }
    (cfg, region(sim, Segment::V, 0, 8))
}

#[test]
fn explicit_read_loss_changes_coverage_without_deletion_columns_or_identity_penalty() {
    use crate::passes::LossEnd;
    for segment in [Segment::V, Segment::D, Segment::J] {
        for reverse in [false, true] {
            if reverse && segment != Segment::D {
                continue;
            }
            for loss_end in [LossEnd::Five, LossEnd::Three] {
                let mut cfg = RefDataConfig::empty(ChainType::Vdj);
                let source = b"ACGTAACC";
                let entry = allele(segment, "coverage*01", source);
                let id = match segment {
                    Segment::V => cfg.v_pool.push(entry),
                    Segment::D => cfg.d_pool.push(entry),
                    Segment::J => cfg.j_pool.push(entry),
                    _ => unreachable!(),
                };
                let mut instance = AlleleInstance::new(id).with_trim_5(1).with_trim_3(1);
                if reverse {
                    instance = instance.with_orientation(SegmentOrientation::ReverseComplement);
                }
                let mut sim = Simulation::new().with_allele_assigned(segment, instance);
                let positions: Vec<_> = if reverse {
                    (1..7).rev().collect()
                } else {
                    (1..7).collect()
                };
                for position in positions {
                    let base = if reverse {
                        complement_base(source[position])
                    } else {
                        source[position]
                    };
                    sim = push(sim, Nucleotide::germline(base, position as u16, segment));
                }
                sim = region(sim, segment, 0, 6);
                let outcome = lose_end(sim, loss_end, 2);
                let rec = build_airr_record(&outcome, &cfg, "observed-coverage");
                let (cigar, identity, start, end, trim5, trim3) = match segment {
                    Segment::V => (
                        &rec.v_cigar,
                        rec.v_identity,
                        rec.v_germline_start,
                        rec.v_germline_end,
                        rec.v_trim_5,
                        rec.v_trim_3,
                    ),
                    Segment::D => (
                        &rec.d_cigar,
                        rec.d_identity,
                        rec.d_germline_start,
                        rec.d_germline_end,
                        rec.d_trim_5,
                        rec.d_trim_3,
                    ),
                    Segment::J => (
                        &rec.j_cigar,
                        rec.j_identity,
                        rec.j_germline_start,
                        rec.j_germline_end,
                        rec.j_trim_5,
                        rec.j_trim_3,
                    ),
                    _ => unreachable!(),
                };
                let clips_start = (loss_end == LossEnd::Five) != reverse;
                assert_eq!(
                    (start, end),
                    if clips_start {
                        (Some(3), Some(7))
                    } else {
                        (Some(1), Some(5))
                    }
                );
                assert_eq!(cigar, "4M");
                assert_eq!(identity, Some(1.0));
                assert_eq!(
                    (trim5, trim3),
                    (1, 1),
                    "read loss must not rewrite recombination provenance"
                );
                assert_eq!(rec.sequence_alignment, rec.sequence.to_ascii_uppercase());
                assert_eq!(rec.germline_alignment, rec.sequence_alignment);
                // The observation bounds survive a fork with no founding trace.
                let fork = build_airr_record(
                    &outcome_from_sim(outcome.final_simulation().clone()),
                    &cfg,
                    "fork",
                );
                assert_eq!(fork.sequence_alignment, rec.sequence_alignment);
                assert_eq!(fork.germline_alignment, rec.germline_alignment);
            }
        }
    }
}

#[test]
fn structural_terminal_deletions_remain_inside_observed_coverage() {
    use crate::passes::LossEnd;
    let (cfg, original) = full_v_for_coverage();
    let structural = original.with_indel_deleted(0).with_indel_deleted(0);
    let rec = build_airr_record(&outcome_from_sim(structural), &cfg, "structural-only");
    assert_eq!(rec.v_cigar, "2D6M");
    assert_eq!(
        (rec.v_germline_start, rec.v_germline_end),
        (Some(0), Some(8))
    );
    assert_eq!(rec.v_identity, Some(0.75));

    // Clip one base, then delete the new terminal base biologically.
    let clipped = lose_end(original.clone(), LossEnd::Five, 1);
    let deleted = clipped.final_simulation().with_indel_deleted(0);
    let rec = build_airr_record(&outcome_from_sim(deleted), &cfg, "loss-then-deletion");
    assert_eq!(rec.v_cigar, "1D6M");
    assert_eq!(rec.sequence_alignment, "-CCGGTT");
    assert_eq!(rec.germline_alignment, "ACCGGTT");
    assert_eq!(
        (rec.v_germline_start, rec.v_germline_end),
        (Some(1), Some(8))
    );
    assert_eq!(rec.v_identity, Some(6.0 / 7.0));

    // Reverse operation order: the prior deletion is outside the new read
    // coverage after clipping the next surviving reference base.
    let deleted = original.clone().with_indel_deleted(0);
    let clipped = lose_end(deleted, LossEnd::Five, 1);
    let rec = build_airr_record(&clipped, &cfg, "deletion-then-loss");
    assert_eq!(rec.v_cigar, "6M");
    assert_eq!(
        (rec.v_germline_start, rec.v_germline_end),
        (Some(2), Some(8))
    );
    assert_eq!(rec.v_identity, Some(1.0));

    let clipped = lose_end(original, LossEnd::Three, 1);
    let deleted = clipped.final_simulation().with_indel_deleted(6);
    let rec = build_airr_record(
        &outcome_from_sim(deleted),
        &cfg,
        "three-prime-loss-then-deletion",
    );
    assert_eq!(rec.v_cigar, "6M1D");
    assert_eq!(
        (rec.v_germline_start, rec.v_germline_end),
        (Some(0), Some(7))
    );
    assert_eq!(rec.v_identity, Some(6.0 / 7.0));
}

#[test]
fn clipping_inserted_bases_does_not_consume_germline_coordinates() {
    use crate::passes::LossEnd;
    let (cfg, original) = full_v_for_coverage();
    let inserted = original.with_indel_inserted(
        0,
        Nucleotide::synthetic(b'T', Segment::V, flag::INDEL_INSERTED),
    );
    let clipped = lose_end(inserted, LossEnd::Five, 1);
    let rec = build_airr_record(&clipped, &cfg, "insert-then-loss");
    assert_eq!(rec.sequence_alignment, "AACCGGTT");
    assert_eq!(rec.germline_alignment, "AACCGGTT");
    assert_eq!(rec.v_cigar, "8M");
    assert_eq!(
        (rec.v_germline_start, rec.v_germline_end),
        (Some(0), Some(8))
    );
    assert_eq!(rec.v_identity, Some(1.0));
}

#[test]
fn clipping_across_a_whole_segment_does_not_resurrect_its_alignment() {
    use crate::passes::LossEnd;
    let (mut cfg, mut sim) = full_v_for_coverage();
    let j = cfg.j_pool.push(allele(Segment::J, "j*01", b"TTTGGG"));
    sim = sim.with_allele_assigned(Segment::J, AlleleInstance::new(j));
    for (position, &base) in b"TTTGGG".iter().enumerate() {
        sim = push(sim, Nucleotide::germline(base, position as u16, Segment::J));
    }
    sim = region(sim, Segment::J, 8, 14);
    let clipped = lose_end(sim, LossEnd::Five, 10);
    let rec = build_airr_record(&clipped, &cfg, "clipped-v-and-j-prefix");
    assert_eq!(rec.sequence_alignment, "TGGG");
    assert_eq!(rec.germline_alignment, "TGGG");
    assert_eq!(rec.v_cigar, "");
    assert_eq!(rec.v_identity, None);
    assert_eq!((rec.v_germline_start, rec.v_germline_end), (None, None));
    assert_eq!(rec.j_cigar, "4M");
    assert_eq!(
        (rec.j_germline_start, rec.j_germline_end),
        (Some(2), Some(6))
    );
    assert_eq!(rec.j_identity, Some(1.0));
}

#[test]
fn clipping_crosses_already_deleted_d_without_leaving_unobserved_deletion_columns() {
    use crate::passes::LossEnd;
    for reverse in [false, true] {
        let mut cfg = RefDataConfig::empty(ChainType::Vdj);
        let v = cfg.v_pool.push(allele(Segment::V, "v*01", b"AAA"));
        let d = cfg.d_pool.push(allele(Segment::D, "d*01", b"CCG"));
        let j = cfg.j_pool.push(allele(Segment::J, "j*01", b"TTT"));
        let mut d_assignment = AlleleInstance::new(d);
        if reverse {
            d_assignment = d_assignment.with_orientation(SegmentOrientation::ReverseComplement);
        }
        let mut sim = Simulation::new()
            .with_allele_assigned(Segment::V, AlleleInstance::new(v))
            .with_allele_assigned(Segment::D, d_assignment)
            .with_allele_assigned(Segment::J, AlleleInstance::new(j));
        for (segment, source) in [
            (Segment::V, b"AAA"),
            (Segment::D, b"CCG"),
            (Segment::J, b"TTT"),
        ] {
            let start = sim.pool.len() as u32;
            for index in 0..3 {
                let position = if segment == Segment::D && reverse {
                    2 - index
                } else {
                    index
                };
                let base = if segment == Segment::D && reverse {
                    complement_base(source[position])
                } else {
                    source[position]
                };
                sim = push(sim, Nucleotide::germline(base, position as u16, segment));
            }
            sim = region(sim, segment, start, start + 3);
        }
        let deleted = sim
            .with_indel_deleted(3)
            .with_indel_deleted(3)
            .with_indel_deleted(3);
        let before = build_airr_record(&outcome_from_sim(deleted.clone()), &cfg, "deleted-d");
        assert_eq!(before.d_cigar, "3D");
        for end in [LossEnd::Five, LossEnd::Three] {
            let clipped = lose_end(deleted.clone(), end, 4);
            let rec = build_airr_record(&clipped, &cfg, "unobserved-deleted-d");
            assert_eq!(
                rec.sequence_alignment,
                if end == LossEnd::Five { "TT" } else { "AA" }
            );
            assert_eq!(rec.germline_alignment, rec.sequence_alignment);
            assert_eq!(rec.d_cigar, "");
            assert_eq!((rec.d_germline_start, rec.d_germline_end), (None, None));
            assert_eq!(rec.d_identity, None);
        }
    }
}
