use super::*;

#[test]
fn simulation_assignments_default_is_empty() {
    let sim = Simulation::new();
    assert!(sim.assignments.v.is_none());
    assert!(sim.assignments.d.is_none());
    assert!(sim.assignments.j.is_none());
    assert!(sim.assignments.c.is_none());
}

#[test]
fn simulation_with_allele_assigned_populates_slot_persistently() {
    use crate::assignment::AlleleInstance;
    use crate::refdata::AlleleId;

    let s0 = Simulation::new();
    let s1 = s0.with_allele_assigned(Segment::V, AlleleInstance::new(AlleleId::new(7)));

    // s0 unchanged.
    assert!(s0.assignments.v.is_none());
    // s1 has the V allele.
    assert_eq!(s1.assignments.v.unwrap().allele_id, AlleleId::new(7));
    assert_eq!(s1.assignments.v.unwrap().trim_5, 0);
    assert_eq!(s1.assignments.v.unwrap().trim_3, 0);
}

#[test]
fn simulation_with_trim_updates_assignment_persistently() {
    use crate::assignment::{AlleleInstance, TrimEnd};
    use crate::refdata::AlleleId;

    let s0 =
        Simulation::new().with_allele_assigned(Segment::V, AlleleInstance::new(AlleleId::new(0)));
    let s1 = s0.with_trim(Segment::V, TrimEnd::Three, 5);

    // s0 retains zero trim.
    assert_eq!(s0.assignments.v.unwrap().trim_3, 0);
    // s1 has the trim applied.
    assert_eq!(s1.assignments.v.unwrap().trim_3, 5);
    // 5' end untouched.
    assert_eq!(s1.assignments.v.unwrap().trim_5, 0);
}

#[test]
fn simulation_with_base_changed_refreshes_overlapping_region_codon_rail() {
    // Audit finding: `with_base_changed` previously updated the
    // pool but cloned the sequence verbatim, leaving
    // `Region.amino_acids` stale against the new pool. This
    // would silently desync any future SHM/uniform pass that
    // mutates a base inside an assembled region.
    //
    // Pin the fix: change a base inside a Region and verify the
    // Region's codon rail reflects the new base.

    // Build a sim with a 9-base V region whose codon rail says KPG.
    let mut sim = Simulation::new();
    for (i, b) in b"AAACCCGGG".iter().enumerate() {
        let (next, _) = sim.with_nucleotide_pushed(Nucleotide::germline(*b, i as u16, Segment::V));
        sim = next;
    }
    let region = Region::new(Segment::V, NucHandle::new(0), NucHandle::new(9))
        .with_codon_rail_recomputed(&sim.pool);
    sim = sim.with_region_added(region);

    assert_eq!(sim.sequence.regions[0].amino_acids, b"KPG");

    // Mutate position 5 (C → A): codon at [3..6) was CCC = P, now becomes CCA = P (no change).
    // Mutate position 4 (C → A): codon CCC → CAC = H. Should change amino_acids[1] to H.
    let mutated = sim.with_base_changed(NucHandle::new(4), b'A');
    assert_eq!(mutated.sequence.regions[0].amino_acids, b"KHG");
    // Pool reflects the change.
    assert_eq!(mutated.pool.get(NucHandle::new(4)).unwrap().base, b'A');
    // Original sim's region still says KPG (persistent IR).
    assert_eq!(sim.sequence.regions[0].amino_acids, b"KPG");
}

#[test]
fn simulation_with_base_changed_outside_any_region_leaves_rails_alone() {
    // Edge case: change a base that's NOT inside any region.
    // Codon rails should remain whatever they were; the change
    // affects only the pool.
    let mut sim = Simulation::new();
    for (i, b) in b"AAACCCGGG".iter().enumerate() {
        let (next, _) = sim.with_nucleotide_pushed(Nucleotide::germline(*b, i as u16, Segment::V));
        sim = next;
    }
    // Add ONE region covering [0, 6) with codon rail.
    let region = Region::new(Segment::V, NucHandle::new(0), NucHandle::new(6))
        .with_codon_rail_recomputed(&sim.pool);
    sim = sim.with_region_added(region);

    assert_eq!(sim.sequence.regions[0].amino_acids, b"KP");

    // Mutate position 7 (outside the region's [0, 6) range).
    let mutated = sim.with_base_changed(NucHandle::new(7), b'A');
    // Region's codon rail unchanged.
    assert_eq!(mutated.sequence.regions[0].amino_acids, b"KP");
}

#[test]
fn simulation_with_base_changed_refreshes_correct_region_among_many() {
    // Multi-region case: only the region covering the changed
    // handle should have its rail recomputed; others stay
    // intact (and should be the same Region values, by Eq).
    let mut sim = Simulation::new();
    for (i, b) in b"AAACCCGGG".iter().enumerate() {
        let (next, _) = sim.with_nucleotide_pushed(Nucleotide::germline(*b, i as u16, Segment::V));
        sim = next;
    }
    // Two regions: [0, 3) with codon AAA → K, and [6, 9) with codon GGG → G.
    let r0 = Region::new(Segment::V, NucHandle::new(0), NucHandle::new(3))
        .with_codon_rail_recomputed(&sim.pool);
    let r1 = Region::new(Segment::V, NucHandle::new(6), NucHandle::new(9))
        .with_codon_rail_recomputed(&sim.pool);
    sim = sim.with_region_added(r0).with_region_added(r1);

    assert_eq!(sim.sequence.regions[0].amino_acids, b"K");
    assert_eq!(sim.sequence.regions[1].amino_acids, b"G");

    // Mutate position 7 (inside r1: GGG → GAG = E).
    let mutated = sim.with_base_changed(NucHandle::new(7), b'A');
    // r0 unchanged.
    assert_eq!(mutated.sequence.regions[0].amino_acids, b"K");
    // r1 recomputed: GAG = E.
    assert_eq!(mutated.sequence.regions[1].amino_acids, b"E");
}
