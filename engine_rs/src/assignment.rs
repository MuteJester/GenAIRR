//! Per-simulation allele state.
//!
//! ## What lives here
//!
//! `AlleleInstance` and `AlleleAssignments` — the per-simulation
//! record of "which V/D/J alleles got sampled and what trim was
//! applied." Reference data (the immutable allele pool) lives in
//! `refdata.rs`; this module is the mutable IR side.
//!
//! ## Relationship to the entity hierarchy
//!
//! Per §3 of the design doc, `AlleleInstance` is conceptually owned
//! by a `Region` (V/D/J — never NP). In this implementation we keep
//! the instances in a flat `AlleleAssignments` slot on the
//! `Simulation` rather than embedded in each `Region`, because:
//!
//! - Allele sampling happens *before* assembly — pre-assembly there
//!   are no `Region`s yet, but the sampled allele has to live
//!   somewhere.
//! - At any moment the simulation has at most one V instance, one D
//!   instance (VDJ only), one J instance, and optionally one C. A
//!   four-slot record is the natural shape.
//! - `Region` looks up its instance by segment via
//!   `simulation.assignments.get(region.segment)` rather than
//!   carrying a redundant copy.
//!
//! ## Scope
//!
//! Just the instance / assignments / trim-end types and their
//! persistent `with_*` API. No sampling pass, no assembly — those
//! arrive in C.5–C.8.

use crate::ir::Segment;
use crate::refdata::AlleleId;

// ──────────────────────────────────────────────────────────────────
// AlleleInstance — per-simulation state of one sampled allele
// ──────────────────────────────────────────────────────────────────

/// One sampled allele plus the trim choices applied to it.
///
/// `allele_id` indexes into the appropriate `AllelePool` of the
/// `RefDataConfig` the simulation was built against. The trim fields
/// give the number of bases removed from each end of the reference
/// sequence; the post-trim retained slice is what assembly (C.8)
/// will copy into the nucleotide pool.
///
/// `trim_5` and `trim_3` are sized as `u16` since allele lengths in
/// the reference data sit comfortably within that range (~300 bases
/// for V, ~100 for J, etc.). A trim larger than the reference
/// length is a caller bug; assembly will catch it.
#[derive(Copy, Clone, Eq, PartialEq, Debug)]
pub struct AlleleInstance {
    pub allele_id: AlleleId,
    pub trim_5: u16,
    pub trim_3: u16,
}

impl AlleleInstance {
    /// Construct a fresh instance with both trims at zero. Trim
    /// passes (C.6) update the trim fields via the persistent API.
    pub fn new(allele_id: AlleleId) -> Self {
        Self {
            allele_id,
            trim_5: 0,
            trim_3: 0,
        }
    }

    /// Return a new instance with `trim_5` updated; receiver
    /// unchanged. Persistent IR contract (D1).
    #[must_use]
    pub fn with_trim_5(self, trim_5: u16) -> Self {
        Self { trim_5, ..self }
    }

    /// Return a new instance with `trim_3` updated; receiver
    /// unchanged.
    #[must_use]
    pub fn with_trim_3(self, trim_3: u16) -> Self {
        Self { trim_3, ..self }
    }
}

// ──────────────────────────────────────────────────────────────────
// TrimEnd — which side of an allele a trim applies to
// ──────────────────────────────────────────────────────────────────

/// Which end of a segment a trim operation applies to.
///
/// `Five` = 5' end of the segment (start side; e.g., J 5' trim
/// removes bases from the start of the J reference). `Three` = 3'
/// end (terminal side; e.g., V 3' trim removes bases from the end
/// of the V reference).
///
/// Biology: V is conventionally trimmed on the 3' end (the side
/// that joins NP1); J on the 5' end (joining NP1 for VJ chains or
/// NP2 for VDJ chains); D on both ends.
#[derive(Copy, Clone, Eq, PartialEq, Hash, Debug)]
pub enum TrimEnd {
    Five,
    Three,
}

// ──────────────────────────────────────────────────────────────────
// AlleleAssignments — the four optional slots on a simulation
// ──────────────────────────────────────────────────────────────────

/// Per-simulation allele assignments, one optional slot per
/// segment role.
///
/// Pre-recombination this is all-`None`. Each `SampleAllele*` pass
/// (C.5) populates one slot. NP segments do not appear here — they
/// have no allele.
#[derive(Copy, Clone, Eq, PartialEq, Debug, Default)]
pub struct AlleleAssignments {
    pub v: Option<AlleleInstance>,
    pub d: Option<AlleleInstance>,
    pub j: Option<AlleleInstance>,
    pub c: Option<AlleleInstance>,
}

impl AlleleAssignments {
    pub fn new() -> Self {
        Self::default()
    }

    /// Look up the assigned instance for a segment role. Returns
    /// `None` for unassigned slots and for NP segments (which never
    /// have alleles).
    pub fn get(&self, segment: Segment) -> Option<&AlleleInstance> {
        match segment {
            Segment::V => self.v.as_ref(),
            Segment::D => self.d.as_ref(),
            Segment::J => self.j.as_ref(),
            Segment::Np1 | Segment::Np2 => None,
        }
    }

    /// Whether a particular segment has an assigned instance.
    pub fn has(&self, segment: Segment) -> bool {
        self.get(segment).is_some()
    }

    // ── Persistent update API ──────────────────────────────────────

    /// Return new assignments with `instance` set on the given
    /// segment slot; receiver unchanged. Panics if `segment` is an
    /// NP segment — those don't have alleles.
    pub fn with_assigned(&self, segment: Segment, instance: AlleleInstance) -> Self {
        let mut next = *self;
        match segment {
            Segment::V => next.v = Some(instance),
            Segment::D => next.d = Some(instance),
            Segment::J => next.j = Some(instance),
            Segment::Np1 | Segment::Np2 => panic!(
                "AlleleAssignments::with_assigned: cannot assign an allele \
                 to NP segment {:?}",
                segment
            ),
        }
        next
    }

    /// Return new assignments with the given segment's instance
    /// updated by `update`. Panics if no instance is currently
    /// assigned to that segment, or if `segment` is NP.
    pub fn with_updated<F>(&self, segment: Segment, update: F) -> Self
    where
        F: FnOnce(AlleleInstance) -> AlleleInstance,
    {
        let current = self.get(segment).copied().unwrap_or_else(|| {
            panic!(
                "AlleleAssignments::with_updated: no instance assigned to \
                 segment {:?} — assign one first via with_assigned",
                segment
            )
        });
        self.with_assigned(segment, update(current))
    }

    /// Return new assignments with the given segment's trim updated.
    /// Convenience wrapper over `with_updated`.
    pub fn with_trim(&self, segment: Segment, end: TrimEnd, value: u16) -> Self {
        self.with_updated(segment, |inst| match end {
            TrimEnd::Five => inst.with_trim_5(value),
            TrimEnd::Three => inst.with_trim_3(value),
        })
    }
}

// ──────────────────────────────────────────────────────────────────
// Tests
// ──────────────────────────────────────────────────────────────────

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn allele_instance_starts_with_zero_trims() {
        let inst = AlleleInstance::new(AlleleId::new(7));
        assert_eq!(inst.allele_id, AlleleId::new(7));
        assert_eq!(inst.trim_5, 0);
        assert_eq!(inst.trim_3, 0);
    }

    #[test]
    fn allele_instance_with_trim_5_isolates() {
        let a = AlleleInstance::new(AlleleId::new(0));
        let b = a.with_trim_5(3);
        assert_eq!(a.trim_5, 0);
        assert_eq!(b.trim_5, 3);
        assert_eq!(b.trim_3, 0); // unchanged
        assert_eq!(b.allele_id, AlleleId::new(0)); // unchanged
    }

    #[test]
    fn allele_instance_with_trim_3_isolates() {
        let a = AlleleInstance::new(AlleleId::new(0));
        let b = a.with_trim_3(5);
        assert_eq!(a.trim_3, 0);
        assert_eq!(b.trim_3, 5);
        assert_eq!(b.trim_5, 0); // unchanged
    }

    #[test]
    fn allele_instance_chained_trims_compose() {
        let a = AlleleInstance::new(AlleleId::new(0))
            .with_trim_5(2)
            .with_trim_3(7);
        assert_eq!(a.trim_5, 2);
        assert_eq!(a.trim_3, 7);
    }

    #[test]
    fn assignments_starts_empty() {
        let a = AlleleAssignments::new();
        assert!(a.v.is_none());
        assert!(a.d.is_none());
        assert!(a.j.is_none());
        assert!(a.c.is_none());
        assert!(!a.has(Segment::V));
        assert!(!a.has(Segment::J));
        assert!(a.get(Segment::V).is_none());
    }

    #[test]
    fn assignments_get_returns_none_for_np_segments() {
        let a = AlleleAssignments::new();
        assert!(a.get(Segment::Np1).is_none());
        assert!(a.get(Segment::Np2).is_none());
        assert!(!a.has(Segment::Np1));
    }

    #[test]
    fn assignments_with_assigned_isolates_revisions() {
        let a = AlleleAssignments::new();
        let inst = AlleleInstance::new(AlleleId::new(3));
        let b = a.with_assigned(Segment::V, inst);

        // Old assignments untouched.
        assert!(a.v.is_none());
        // New assignments has the V instance.
        assert_eq!(b.v, Some(inst));
        assert!(b.d.is_none());
        assert!(b.j.is_none());
        assert!(b.has(Segment::V));
    }

    #[test]
    fn assignments_with_assigned_works_for_all_real_segments() {
        let a = AlleleAssignments::new();
        let inst_v = AlleleInstance::new(AlleleId::new(1));
        let inst_d = AlleleInstance::new(AlleleId::new(2));
        let inst_j = AlleleInstance::new(AlleleId::new(3));

        let b = a
            .with_assigned(Segment::V, inst_v)
            .with_assigned(Segment::D, inst_d)
            .with_assigned(Segment::J, inst_j);

        assert_eq!(b.v, Some(inst_v));
        assert_eq!(b.d, Some(inst_d));
        assert_eq!(b.j, Some(inst_j));
    }

    #[test]
    #[should_panic(expected = "cannot assign an allele to NP segment")]
    fn assignments_with_assigned_rejects_np1() {
        let a = AlleleAssignments::new();
        let _ = a.with_assigned(Segment::Np1, AlleleInstance::new(AlleleId::new(0)));
    }

    #[test]
    #[should_panic(expected = "cannot assign an allele to NP segment")]
    fn assignments_with_assigned_rejects_np2() {
        let a = AlleleAssignments::new();
        let _ = a.with_assigned(Segment::Np2, AlleleInstance::new(AlleleId::new(0)));
    }

    #[test]
    #[should_panic(expected = "no instance assigned to segment")]
    fn assignments_with_updated_panics_on_unassigned() {
        let a = AlleleAssignments::new();
        let _ = a.with_updated(Segment::V, |i| i.with_trim_3(5));
    }

    #[test]
    fn assignments_with_updated_modifies_existing_instance() {
        let a = AlleleAssignments::new()
            .with_assigned(Segment::V, AlleleInstance::new(AlleleId::new(7)));
        let b = a.with_updated(Segment::V, |i| i.with_trim_3(4));

        // Old assignments retains zero trim.
        assert_eq!(a.v.unwrap().trim_3, 0);
        // New assignments has updated trim, same allele_id.
        assert_eq!(b.v.unwrap().trim_3, 4);
        assert_eq!(b.v.unwrap().allele_id, AlleleId::new(7));
    }

    #[test]
    fn assignments_with_trim_updates_correct_end() {
        let a = AlleleAssignments::new()
            .with_assigned(Segment::V, AlleleInstance::new(AlleleId::new(0)))
            .with_assigned(Segment::J, AlleleInstance::new(AlleleId::new(1)));

        let b = a
            .with_trim(Segment::V, TrimEnd::Three, 6)
            .with_trim(Segment::J, TrimEnd::Five, 2);

        // V got 3' trim, J got 5' trim. Other ends and segments unchanged.
        assert_eq!(b.v.unwrap().trim_5, 0);
        assert_eq!(b.v.unwrap().trim_3, 6);
        assert_eq!(b.j.unwrap().trim_5, 2);
        assert_eq!(b.j.unwrap().trim_3, 0);

        // Original assignments untouched.
        assert_eq!(a.v.unwrap().trim_3, 0);
        assert_eq!(a.j.unwrap().trim_5, 0);
    }

    #[test]
    fn assignments_branching_revisions_are_independent() {
        let base = AlleleAssignments::new()
            .with_assigned(Segment::V, AlleleInstance::new(AlleleId::new(5)));
        let branch_a = base.with_trim(Segment::V, TrimEnd::Three, 1);
        let branch_b = base.with_trim(Segment::V, TrimEnd::Three, 9);

        assert_eq!(base.v.unwrap().trim_3, 0);
        assert_eq!(branch_a.v.unwrap().trim_3, 1);
        assert_eq!(branch_b.v.unwrap().trim_3, 9);
    }

    #[test]
    fn assignments_assign_replaces_previous_instance() {
        let a = AlleleAssignments::new()
            .with_assigned(Segment::V, AlleleInstance::new(AlleleId::new(1)));
        let b = a.with_assigned(Segment::V, AlleleInstance::new(AlleleId::new(2)));

        assert_eq!(a.v.unwrap().allele_id, AlleleId::new(1));
        assert_eq!(b.v.unwrap().allele_id, AlleleId::new(2));
    }
}
