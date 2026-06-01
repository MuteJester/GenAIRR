use super::*;

#[test]
fn simulation_assignments_default_is_empty() {
    let sim = Simulation::new();
    assert!(sim.assignments.get(Segment::V).is_none());
    assert!(sim.assignments.get(Segment::D).is_none());
    assert!(sim.assignments.get(Segment::J).is_none());
    assert_eq!(sim.assignments.iter().count(), 0);
}

#[test]
fn simulation_with_allele_assigned_populates_slot_persistently() {
    use crate::assignment::AlleleInstance;
    use crate::refdata::AlleleId;

    let s0 = Simulation::new();
    let s1 = s0.with_allele_assigned(Segment::V, AlleleInstance::new(AlleleId::new(7)));

    // s0 unchanged.
    assert!(s0.assignments.get(Segment::V).is_none());
    // s1 has the V allele.
    assert_eq!(
        s1.assignments.get(Segment::V).copied().unwrap().allele_id,
        AlleleId::new(7)
    );
    assert_eq!(s1.assignments.get(Segment::V).copied().unwrap().trim_5, 0);
    assert_eq!(s1.assignments.get(Segment::V).copied().unwrap().trim_3, 0);
}

#[test]
fn simulation_with_trim_updates_assignment_persistently() {
    use crate::assignment::{AlleleInstance, TrimEnd};
    use crate::refdata::AlleleId;

    let s0 =
        Simulation::new().with_allele_assigned(Segment::V, AlleleInstance::new(AlleleId::new(0)));
    let s1 = s0.with_trim(Segment::V, TrimEnd::Three, 5);

    // s0 retains zero trim.
    assert_eq!(s0.assignments.get(Segment::V).copied().unwrap().trim_3, 0);
    // s1 has the trim applied.
    assert_eq!(s1.assignments.get(Segment::V).copied().unwrap().trim_3, 5);
    // 5' end untouched.
    assert_eq!(s1.assignments.get(Segment::V).copied().unwrap().trim_5, 0);
}

// The three tests that previously verified `with_base_changed`
// refreshing the per-region codon rail are now obsolete: the rail
// isn't stored on Region anymore, so there's no field to keep in
// sync. The equivalent invariant — "computing the rail on demand
// against the new pool reflects the mutation" — is covered by the
// canonical codon-rail tests in `codon_rail.rs`.

// ── D-inversion Slice A: Simulation::with_allele_orientation ──────

#[test]
fn simulation_with_allele_orientation_updates_target_persistently() {
    use crate::assignment::{AlleleInstance, SegmentOrientation};
    use crate::refdata::AlleleId;

    let s0 = Simulation::new()
        .with_allele_assigned(Segment::D, AlleleInstance::new(AlleleId::new(2)));
    let s1 = s0.with_allele_orientation(Segment::D, SegmentOrientation::ReverseComplement);

    // s0 retains forward orientation.
    assert_eq!(
        s0.assignments.get(Segment::D).copied().unwrap().orientation,
        SegmentOrientation::Forward,
    );
    // s1 has the orientation flipped.
    assert_eq!(
        s1.assignments.get(Segment::D).copied().unwrap().orientation,
        SegmentOrientation::ReverseComplement,
    );
    // allele_id and trims untouched.
    let updated = s1.assignments.get(Segment::D).copied().unwrap();
    assert_eq!(updated.allele_id, AlleleId::new(2));
    assert_eq!(updated.trim_5, 0);
    assert_eq!(updated.trim_3, 0);
}

#[test]
fn simulation_with_allele_orientation_only_touches_target_segment() {
    use crate::assignment::{AlleleInstance, SegmentOrientation};
    use crate::refdata::AlleleId;

    let s = Simulation::new()
        .with_allele_assigned(Segment::V, AlleleInstance::new(AlleleId::new(0)))
        .with_allele_assigned(Segment::D, AlleleInstance::new(AlleleId::new(1)))
        .with_allele_assigned(Segment::J, AlleleInstance::new(AlleleId::new(2)));
    let s1 = s.with_allele_orientation(Segment::D, SegmentOrientation::ReverseComplement);

    assert_eq!(
        s1.assignments.get(Segment::V).copied().unwrap().orientation,
        SegmentOrientation::Forward,
    );
    assert_eq!(
        s1.assignments.get(Segment::D).copied().unwrap().orientation,
        SegmentOrientation::ReverseComplement,
    );
    assert_eq!(
        s1.assignments.get(Segment::J).copied().unwrap().orientation,
        SegmentOrientation::Forward,
    );
}

#[test]
fn simulation_with_allele_orientation_preserves_other_simulation_state() {
    use crate::assignment::{AlleleInstance, SegmentOrientation};
    use crate::refdata::AlleleId;

    // Pin that flipping orientation is a *pure* assignments update —
    // scalar state on the simulation (mutation_count, pool length)
    // is preserved across the persistent revision. Slice A's
    // byte-identical guarantee depends on this.
    let s = Simulation::new()
        .with_allele_assigned(Segment::D, AlleleInstance::new(AlleleId::new(0)));
    let pool_len_before = s.pool.len();
    let regions_before = s.sequence.regions.len();
    let mutation_before = s.mutation_count;

    let s1 = s.with_allele_orientation(Segment::D, SegmentOrientation::ReverseComplement);
    assert_eq!(s1.pool.len(), pool_len_before);
    assert_eq!(s1.sequence.regions.len(), regions_before);
    assert_eq!(s1.mutation_count, mutation_before);
}
