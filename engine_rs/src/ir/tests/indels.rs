use super::*;

// ── Indel IR API tests ─────────────────────────────────────────

#[test]
fn pool_with_inserted_grows_and_shifts() {
    let mut p = NucleotidePool::new();
    let _ = p.push(Nucleotide::germline(b'A', 0, Segment::V));
    let _ = p.push(Nucleotide::germline(b'C', 1, Segment::V));
    let _ = p.push(Nucleotide::germline(b'G', 2, Segment::V));

    let p2 = p.with_inserted(1, Nucleotide::germline(b'X', 99, Segment::Np1));
    assert_eq!(p2.len(), 4);
    assert_eq!(p2.get(NucHandle::new(0)).unwrap().base, b'A');
    assert_eq!(p2.get(NucHandle::new(1)).unwrap().base, b'X');
    assert_eq!(p2.get(NucHandle::new(2)).unwrap().base, b'C'); // shifted up
    assert_eq!(p2.get(NucHandle::new(3)).unwrap().base, b'G');

    // Original pool unchanged (persistent IR).
    assert_eq!(p.len(), 3);
}

#[test]
fn pool_with_inserted_at_end_appends() {
    let mut p = NucleotidePool::new();
    let _ = p.push(Nucleotide::germline(b'A', 0, Segment::V));
    let p2 = p.with_inserted(1, Nucleotide::germline(b'C', 1, Segment::V));
    assert_eq!(p2.len(), 2);
    assert_eq!(p2.get(NucHandle::new(1)).unwrap().base, b'C');
}

#[test]
#[should_panic(expected = "with_inserted: at")]
fn pool_with_inserted_out_of_bounds_panics() {
    let p = NucleotidePool::new();
    let _ = p.with_inserted(5, Nucleotide::germline(b'A', 0, Segment::V));
}

#[test]
fn pool_with_deleted_shrinks_and_shifts() {
    let mut p = NucleotidePool::new();
    let _ = p.push(Nucleotide::germline(b'A', 0, Segment::V));
    let _ = p.push(Nucleotide::germline(b'C', 1, Segment::V));
    let _ = p.push(Nucleotide::germline(b'G', 2, Segment::V));

    let p2 = p.with_deleted(1);
    assert_eq!(p2.len(), 2);
    assert_eq!(p2.get(NucHandle::new(0)).unwrap().base, b'A');
    assert_eq!(p2.get(NucHandle::new(1)).unwrap().base, b'G'); // shifted down
    assert_eq!(p.len(), 3); // original unchanged
}

#[test]
#[should_panic(expected = "with_deleted: at")]
fn pool_with_deleted_out_of_bounds_panics() {
    let p = NucleotidePool::new();
    let _ = p.with_deleted(0);
}

#[test]
fn sequence_indel_adjustment_region_entirely_before() {
    let s = Sequence::new().with_region_added(Region::new(
        Segment::V,
        NucHandle::new(0),
        NucHandle::new(5),
    ));
    // Insertion at position 7 (after region.end = 5).
    let s2 = s.with_indel_adjusted(7, 1);
    // Region unchanged.
    assert_eq!(s2.regions[0].start, NucHandle::new(0));
    assert_eq!(s2.regions[0].end, NucHandle::new(5));
}

#[test]
fn sequence_indel_adjustment_region_entirely_after_insertion() {
    let s = Sequence::new().with_region_added(Region::new(
        Segment::V,
        NucHandle::new(5),
        NucHandle::new(10),
    ));
    // Insertion at position 3 — both start and end shift up by 1.
    let s2 = s.with_indel_adjusted(3, 1);
    assert_eq!(s2.regions[0].start, NucHandle::new(6));
    assert_eq!(s2.regions[0].end, NucHandle::new(11));
}

#[test]
fn sequence_indel_adjustment_region_spanning_insertion_grows() {
    let s = Sequence::new().with_region_added(Region::new(
        Segment::V,
        NucHandle::new(0),
        NucHandle::new(10),
    ));
    // Insertion at position 5 (inside the region).
    let s2 = s.with_indel_adjusted(5, 1);
    // Start unchanged, end grew.
    assert_eq!(s2.regions[0].start, NucHandle::new(0));
    assert_eq!(s2.regions[0].end, NucHandle::new(11));
}

#[test]
fn sequence_indel_adjustment_region_spanning_deletion_shrinks() {
    let s = Sequence::new().with_region_added(Region::new(
        Segment::V,
        NucHandle::new(0),
        NucHandle::new(10),
    ));
    // Deletion at position 5 (inside the region).
    let s2 = s.with_indel_adjusted(5, -1);
    assert_eq!(s2.regions[0].start, NucHandle::new(0));
    assert_eq!(s2.regions[0].end, NucHandle::new(9));
}

#[test]
fn sequence_indel_adjustment_region_after_deletion_shifts() {
    let s = Sequence::new().with_region_added(Region::new(
        Segment::V,
        NucHandle::new(5),
        NucHandle::new(10),
    ));
    // Deletion at position 3 (before the region).
    let s2 = s.with_indel_adjusted(3, -1);
    assert_eq!(s2.regions[0].start, NucHandle::new(4));
    assert_eq!(s2.regions[0].end, NucHandle::new(9));
}

#[test]
fn sequence_frame_phase_recompute_chains_from_region_lengths() {
    let s = Sequence::new()
        .with_region_added(
            Region::new(Segment::V, NucHandle::new(0), NucHandle::new(1)).with_frame_phase(2),
        )
        .with_region_added(
            Region::new(Segment::Np1, NucHandle::new(1), NucHandle::new(3)).with_frame_phase(2),
        )
        .with_region_added(
            Region::new(Segment::J, NucHandle::new(3), NucHandle::new(6)).with_frame_phase(2),
        );

    let recomputed = s.with_frame_phases_recomputed();

    assert_eq!(recomputed.regions[0].frame_phase, 0);
    assert_eq!(recomputed.regions[1].frame_phase, 1);
    assert_eq!(recomputed.regions[2].frame_phase, 0);
}

#[test]
fn sequence_indel_adjustment_recomputes_downstream_frame_phase_after_insertion() {
    let s = Sequence::new()
        .with_region_added(
            Region::new(Segment::V, NucHandle::new(0), NucHandle::new(2)).with_frame_phase(0),
        )
        .with_region_added(
            Region::new(Segment::J, NucHandle::new(2), NucHandle::new(8)).with_frame_phase(2),
        );

    let adjusted = s.with_indel_adjusted(1, 1);

    assert_eq!(adjusted.regions[0].start, NucHandle::new(0));
    assert_eq!(adjusted.regions[0].end, NucHandle::new(3));
    assert_eq!(adjusted.regions[0].frame_phase, 0);
    assert_eq!(adjusted.regions[1].start, NucHandle::new(3));
    assert_eq!(adjusted.regions[1].end, NucHandle::new(9));
    assert_eq!(adjusted.regions[1].frame_phase, 0);
}

#[test]
fn sequence_indel_adjustment_recomputes_downstream_frame_phase_after_deletion() {
    let s = Sequence::new()
        .with_region_added(
            Region::new(Segment::V, NucHandle::new(0), NucHandle::new(4)).with_frame_phase(0),
        )
        .with_region_added(
            Region::new(Segment::J, NucHandle::new(4), NucHandle::new(10)).with_frame_phase(1),
        );

    let adjusted = s.with_indel_adjusted(1, -1);

    assert_eq!(adjusted.regions[0].start, NucHandle::new(0));
    assert_eq!(adjusted.regions[0].end, NucHandle::new(3));
    assert_eq!(adjusted.regions[0].frame_phase, 0);
    assert_eq!(adjusted.regions[1].start, NucHandle::new(3));
    assert_eq!(adjusted.regions[1].end, NucHandle::new(9));
    assert_eq!(adjusted.regions[1].frame_phase, 0);
}

#[test]
fn simulation_with_indel_inserted_refreshes_downstream_frame_phase_and_rail() {
    let mut sim = Simulation::new();
    for (i, b) in b"AAATGCCC".iter().enumerate() {
        let segment = if i < 2 { Segment::V } else { Segment::J };
        let (next, _) = sim.with_nucleotide_pushed(Nucleotide::germline(*b, i as u16, segment));
        sim = next;
    }
    let v_region =
        Region::new(Segment::V, NucHandle::new(0), NucHandle::new(2)).with_frame_phase(0);
    let j_region =
        Region::new(Segment::J, NucHandle::new(2), NucHandle::new(8)).with_frame_phase(2);
    sim = sim.with_region_added(v_region).with_region_added(j_region);

    let mutated = sim.with_indel_inserted(
        1,
        Nucleotide::synthetic(b'G', Segment::V, crate::ir::flag::INDEL_INSERTED),
    );

    let j = &mutated.sequence.regions[1];
    assert_eq!(j.start, NucHandle::new(3));
    assert_eq!(j.end, NucHandle::new(9));
    assert_eq!(j.frame_phase, 0);
    let rail = crate::ir::compute_codon_rail(j, &mutated.pool);
    assert_eq!(rail.amino_acids, b"MP");
}

#[test]
fn simulation_with_indel_deleted_refreshes_downstream_frame_phase_and_rail() {
    let mut sim = Simulation::new();
    for (i, b) in b"AAAAATGCCC".iter().enumerate() {
        let segment = if i < 4 { Segment::V } else { Segment::J };
        let (next, _) = sim.with_nucleotide_pushed(Nucleotide::germline(*b, i as u16, segment));
        sim = next;
    }
    let v_region =
        Region::new(Segment::V, NucHandle::new(0), NucHandle::new(4)).with_frame_phase(0);
    let j_region =
        Region::new(Segment::J, NucHandle::new(4), NucHandle::new(10)).with_frame_phase(1);
    sim = sim.with_region_added(v_region).with_region_added(j_region);

    let mutated = sim.with_indel_deleted(1);

    let j = &mutated.sequence.regions[1];
    assert_eq!(j.start, NucHandle::new(3));
    assert_eq!(j.end, NucHandle::new(9));
    assert_eq!(j.frame_phase, 0);
    let rail = crate::ir::compute_codon_rail(j, &mutated.pool);
    assert_eq!(rail.amino_acids, b"MP");
}

#[test]
fn simulation_with_indel_inserted_full_round_trip() {
    let mut sim = Simulation::new();
    for (i, b) in b"ATGCCCGGG".iter().enumerate() {
        let (next, _) = sim.with_nucleotide_pushed(Nucleotide::germline(*b, i as u16, Segment::V));
        sim = next;
    }
    let region = Region::new(Segment::V, NucHandle::new(0), NucHandle::new(9));
    sim = sim.with_region_added(region);
    // Initial codon rail: ATG CCC GGG → M P G.
    assert_eq!(
        compute_codon_rail(&sim.sequence.regions[0], &sim.pool).amino_acids,
        b"MPG"
    );

    // Insert a 'G' at position 3 — the region grows by 1, and
    // codons shift: ATG | GCC CGG G → M, A, R, + leftover G.
    let mutated = sim.with_indel_inserted(
        3,
        Nucleotide::synthetic(b'G', Segment::V, crate::ir::flag::INDEL_INSERTED),
    );
    assert_eq!(mutated.pool.len(), 10);
    assert_eq!(mutated.sequence.regions[0].len(), 10);
    // ATG GCC CGG → M A R (3 codons, 10 = 3*3 + 1 leftover)
    assert_eq!(
        compute_codon_rail(&mutated.sequence.regions[0], &mutated.pool).amino_acids,
        b"MAR"
    );
    // Original sim untouched.
    assert_eq!(sim.pool.len(), 9);
    assert_eq!(
        compute_codon_rail(&sim.sequence.regions[0], &sim.pool).amino_acids,
        b"MPG"
    );
}

#[test]
fn simulation_with_indel_deleted_full_round_trip() {
    let mut sim = Simulation::new();
    for (i, b) in b"ATGAAACCCGGG".iter().enumerate() {
        let (next, _) = sim.with_nucleotide_pushed(Nucleotide::germline(*b, i as u16, Segment::V));
        sim = next;
    }
    let region = Region::new(Segment::V, NucHandle::new(0), NucHandle::new(12));
    sim = sim.with_region_added(region);
    // Codon rail: ATG AAA CCC GGG → M K P G.
    assert_eq!(
        compute_codon_rail(&sim.sequence.regions[0], &sim.pool).amino_acids,
        b"MKPG"
    );

    // Delete position 3 (the first 'A' of the second codon).
    // Region shrinks to 11 bases.
    let mutated = sim.with_indel_deleted(3);
    assert_eq!(mutated.pool.len(), 11);
    assert_eq!(mutated.sequence.regions[0].len(), 11);
    // Bases now: ATG AAC CCG GG → 3 full codons + 2 leftover.
    // ATG=M, AAC=N, CCG=P. Three amino acids.
    assert_eq!(
        compute_codon_rail(&mutated.sequence.regions[0], &mutated.pool).amino_acids,
        b"MNP"
    );
}

#[test]
fn simulation_indel_branching_revisions_are_independent() {
    let mut sim = Simulation::new();
    for (i, b) in b"ATGAAA".iter().enumerate() {
        let (next, _) = sim.with_nucleotide_pushed(Nucleotide::germline(*b, i as u16, Segment::V));
        sim = next;
    }

    // Branch A: insert at 3.
    let branch_a = sim.with_indel_inserted(
        3,
        Nucleotide::synthetic(b'C', Segment::V, crate::ir::flag::INDEL_INSERTED),
    );
    // Branch B: delete at 3.
    let branch_b = sim.with_indel_deleted(3);

    assert_eq!(branch_a.pool.len(), 7);
    assert_eq!(branch_b.pool.len(), 5);
    // Original unchanged.
    assert_eq!(sim.pool.len(), 6);
}

#[test]
fn simulation_with_nucleotide_changed_also_refreshes_codon_rail() {
    // Same fix applies to `with_nucleotide_changed`, which can
    // change germline_pos / segment / flags in addition to base.
    let mut sim = Simulation::new();
    for (i, b) in b"AAACCC".iter().enumerate() {
        let (next, _) = sim.with_nucleotide_pushed(Nucleotide::germline(*b, i as u16, Segment::V));
        sim = next;
    }
    let region = Region::new(Segment::V, NucHandle::new(0), NucHandle::new(6));
    sim = sim.with_region_added(region);
    assert_eq!(
        compute_codon_rail(&sim.sequence.regions[0], &sim.pool).amino_acids,
        b"KP"
    );

    // Replace position 3 with a new nucleotide whose base is 'T'.
    // Codon at [3..6) was CCC = P, becomes TCC = S.
    let mutated =
        sim.with_nucleotide_changed(NucHandle::new(3), Nucleotide::germline(b'T', 3, Segment::V));
    assert_eq!(
        compute_codon_rail(&mutated.sequence.regions[0], &mutated.pool).amino_acids,
        b"KS"
    );
}

#[test]
fn simulation_assignments_chain_with_pool_changes() {
    // Assignments survive pool / sequence changes — they live in
    // their own slot, not in the pool.
    use crate::assignment::AlleleInstance;
    use crate::refdata::AlleleId;

    let s0 =
        Simulation::new().with_allele_assigned(Segment::V, AlleleInstance::new(AlleleId::new(3)));

    let (s1, _h) = s0.with_nucleotide_pushed(Nucleotide::germline(b'A', 0, Segment::V));

    // Assignment survives across the pool change.
    assert_eq!(
        s1.assignments.get(Segment::V).copied().unwrap().allele_id,
        AlleleId::new(3)
    );
}

#[test]
fn simulation_branching_revisions_are_independent() {
    // Two divergent branches of mutation from the same parent must
    // not affect each other.
    let s0 = Simulation::new();
    let (s1, h) = s0.with_nucleotide_pushed(Nucleotide::germline(b'A', 0, Segment::V));

    let branch_left = s1.with_base_changed(h, b'C');
    let branch_right = s1.with_base_changed(h, b'G');

    // Common ancestor unchanged.
    assert_eq!(s1.pool.get(h).unwrap().base, b'A');
    // Left branch: A → C.
    assert_eq!(branch_left.pool.get(h).unwrap().base, b'C');
    // Right branch: A → G.
    assert_eq!(branch_right.pool.get(h).unwrap().base, b'G');
}
