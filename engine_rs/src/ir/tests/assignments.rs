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
    assert_eq!(s1.assignments.get(Segment::V).copied().unwrap().allele_id, AlleleId::new(7));
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
