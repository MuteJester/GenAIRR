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

use crate::ir::{PerSegment, Segment};
use crate::refdata::AlleleId;

// ──────────────────────────────────────────────────────────────────
// SegmentOrientation — Forward vs ReverseComplement for one assignment
// ──────────────────────────────────────────────────────────────────

/// Orientation of a sampled segment relative to its reference
/// allele. Lives on [`AlleleInstance`] — same allele can be sampled
/// forward in one simulation and reverse-complemented in the next.
///
/// `Forward` is the established biology default. `ReverseComplement`
/// models V(D)J inversion (RSS-symmetric inversion event for D
/// segments, biologically real with low prevalence; rare but
/// documented for V/J).
///
/// **Slice A scope (IR primitives only).** Adding this enum and the
/// `AlleleInstance.orientation` field is data-model preparation; no
/// pass reads or writes it yet. Assembly emission, trace recording,
/// AIRR projection, and the DSL surface arrive in later slices.
/// See [`docs/d_inversion_design.md`](../../../docs/d_inversion_design.md)
/// for the full plan.
#[derive(Copy, Clone, Eq, PartialEq, Hash, Debug, Default)]
pub enum SegmentOrientation {
    /// Reference orientation. Default for every `AlleleInstance`.
    #[default]
    Forward,
    /// Reverse-complemented relative to the reference allele.
    ReverseComplement,
}

impl SegmentOrientation {
    /// Whether this orientation reverses the reference allele.
    /// `true` only for [`Self::ReverseComplement`].
    pub fn is_reverse(self) -> bool {
        matches!(self, SegmentOrientation::ReverseComplement)
    }
}

// ──────────────────────────────────────────────────────────────────
// AlleleInstance — per-simulation state of one sampled allele
// ──────────────────────────────────────────────────────────────────

/// One sampled allele plus the trim and orientation choices applied
/// to it.
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
///
/// `orientation` carries Slice A's prepared D-inversion data model.
/// Default [`SegmentOrientation::Forward`]; no pass reads or writes
/// it yet (the D-inversion implementation phases land in later
/// slices — see [`docs/d_inversion_design.md`](../../../docs/d_inversion_design.md)).
///
/// `receptor_revision_original_id` is the pre-revision V allele
/// when receptor revision rewrote this slot, or `None` otherwise.
/// It is the IR-side source of truth for the AIRR
/// `receptor_revision_applied` / `original_v_call` projection
/// (see [`docs/clonal_plan_split_design.md`](../../../docs/clonal_plan_split_design.md)
/// / Bug D — trace-sourced provenance silently dropped on
/// post-fork clonal descendants because the descendant trace
/// doesn't carry pre-fork choices, but the assignments do).
/// Only the V slot ever has a non-`None` value today — D / J
/// receptor revision is out of scope.
#[derive(Copy, Clone, Eq, PartialEq, Debug)]
pub struct AlleleInstance {
    pub allele_id: AlleleId,
    pub trim_5: u16,
    pub trim_3: u16,
    pub orientation: SegmentOrientation,
    /// Pre-revision allele id when receptor revision applied; `None`
    /// when revision did not run on this slot. Preserves the
    /// *first* original id across hypothetical chained revisions
    /// (today the DSL prevents chaining, but the field stays
    /// stable under any future relaxation).
    pub receptor_revision_original_id: Option<AlleleId>,
    /// The diploid rearrangement chromosome (0 or 1) this assignment was
    /// sampled from, when phased genotype recombination ran; `None` on the
    /// flat (no-genotype) path. Immutable provenance of the *original*
    /// rearrangement — receptor revision preserves it even when a
    /// cross-haplotype replacement (`same_haplotype=false`) comes from the
    /// other chromosome.
    pub haplotype: Option<u8>,
}

impl AlleleInstance {
    /// Construct a fresh instance with both trims at zero, forward
    /// orientation, and no receptor-revision provenance. Trim
    /// passes (C.6) update the trim fields; the future D-inversion
    /// pass updates `orientation` via [`Self::with_orientation`];
    /// receptor revision installs provenance via
    /// [`Self::with_receptor_revision_original_id`].
    pub fn new(allele_id: AlleleId) -> Self {
        Self {
            allele_id,
            trim_5: 0,
            trim_3: 0,
            orientation: SegmentOrientation::Forward,
            receptor_revision_original_id: None,
            haplotype: None,
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

    /// Return a new instance with `orientation` updated; receiver
    /// unchanged. Mirror of [`Self::with_trim_5`] / [`Self::with_trim_3`].
    #[must_use]
    pub fn with_orientation(self, orientation: SegmentOrientation) -> Self {
        Self { orientation, ..self }
    }

    /// Return a new instance with the receptor-revision pre-revision
    /// allele id set to `original_id`; receiver unchanged. Mirrors
    /// the existing persistent setters. Pass-level helper — the
    /// receptor-revision pass installs this after capturing the
    /// pre-revision V identity; the AIRR builder reads it back via
    /// the field directly.
    #[must_use]
    pub fn with_receptor_revision_original_id(self, original_id: AlleleId) -> Self {
        Self {
            receptor_revision_original_id: Some(original_id),
            ..self
        }
    }

    /// Return a new instance with the rearrangement `haplotype` (0 or 1)
    /// set; receiver unchanged. Panics if `hap` is not 0 or 1 — a diploid
    /// genotype has exactly two haplotypes.
    #[must_use]
    pub fn with_haplotype(self, hap: u8) -> Self {
        assert!(hap <= 1, "haplotype must be 0 or 1, got {}", hap);
        Self {
            haplotype: Some(hap),
            ..self
        }
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

/// Per-simulation allele assignments, indexed by `Segment`.
///
/// Pre-recombination every slot is empty. Each `SampleAllele*` pass
/// populates one slot via [`Self::with_assigned`]. NP segments do
/// not appear here — they have no allele.
///
/// Backed by [`PerSegment`] rather than parallel `Option<...>`
/// fields, which makes adding a new assignable segment (e.g. C-region
/// for constant regions, or a paired-chain second segment) a
/// one-line change to [`Segment::assignable`] — no schema migration
/// here.
#[derive(Copy, Clone, Eq, PartialEq, Debug, Default)]
pub struct AlleleAssignments {
    slots: PerSegment<AlleleInstance>,
}

impl AlleleAssignments {
    pub fn new() -> Self {
        Self::default()
    }

    /// Look up the assigned instance for a segment role. Returns
    /// `None` for unassigned slots and for NP segments (which never
    /// have alleles).
    pub fn get(&self, segment: Segment) -> Option<&AlleleInstance> {
        if !Segment::assignable().contains(&segment) {
            return None;
        }
        self.slots.get(segment)
    }

    /// Whether a particular segment has an assigned instance.
    pub fn has(&self, segment: Segment) -> bool {
        self.get(segment).is_some()
    }

    /// Iterate every populated `(segment, instance)` entry in
    /// canonical segment order (V → D → J). Skips empty slots.
    pub fn iter(&self) -> impl Iterator<Item = (Segment, &AlleleInstance)> + '_ {
        self.slots.iter()
    }

    // ── Persistent update API ──────────────────────────────────────

    /// Return new assignments with `instance` set on the given
    /// segment slot; receiver unchanged. Panics if `segment` is not
    /// in [`Segment::assignable`] (NP segments, or any future
    /// non-allele-bearing segment).
    pub fn with_assigned(&self, segment: Segment, instance: AlleleInstance) -> Self {
        assert!(
            Segment::assignable().contains(&segment),
            "AlleleAssignments::with_assigned: cannot assign an allele \
             to non-assignable segment {:?}",
            segment
        );
        let mut next = *self;
        next.slots.set(segment, instance);
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

    /// Return new assignments with the given segment's orientation
    /// updated. Convenience wrapper over `with_updated`. Panics if
    /// no instance is currently assigned to that segment, with the
    /// same wording as the trim setters so a caller misusing the
    /// API sees a uniform diagnostic.
    pub fn with_orientation(
        &self,
        segment: Segment,
        orientation: SegmentOrientation,
    ) -> Self {
        self.with_updated(segment, |inst| inst.with_orientation(orientation))
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
    fn allele_instance_haplotype_defaults_none_and_with_haplotype_isolates() {
        let a = AlleleInstance::new(AlleleId::new(0));
        assert_eq!(a.haplotype, None);
        let b = a.with_haplotype(1);
        assert_eq!(a.haplotype, None); // receiver unchanged (persistent)
        assert_eq!(b.haplotype, Some(1));
    }

    #[test]
    #[should_panic(expected = "haplotype must be 0 or 1")]
    fn with_haplotype_rejects_out_of_range() {
        let _ = AlleleInstance::new(AlleleId::new(0)).with_haplotype(2);
    }

    #[test]
    fn assignments_starts_empty() {
        let a = AlleleAssignments::new();
        assert!(a.get(Segment::V).is_none());
        assert!(a.get(Segment::D).is_none());
        assert!(a.get(Segment::J).is_none());
        assert!(!a.has(Segment::V));
        assert!(!a.has(Segment::J));
        assert_eq!(a.iter().count(), 0);
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
        assert!(a.get(Segment::V).is_none());
        // New assignments has the V instance.
        assert_eq!(b.get(Segment::V).copied(), Some(inst));
        assert!(b.get(Segment::D).is_none());
        assert!(b.get(Segment::J).is_none());
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

        assert_eq!(b.get(Segment::V).copied(), Some(inst_v));
        assert_eq!(b.get(Segment::D).copied(), Some(inst_d));
        assert_eq!(b.get(Segment::J).copied(), Some(inst_j));
    }

    #[test]
    #[should_panic(expected = "cannot assign an allele to non-assignable segment")]
    fn assignments_with_assigned_rejects_np1() {
        let a = AlleleAssignments::new();
        let _ = a.with_assigned(Segment::Np1, AlleleInstance::new(AlleleId::new(0)));
    }

    #[test]
    #[should_panic(expected = "cannot assign an allele to non-assignable segment")]
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
        assert_eq!(a.get(Segment::V).unwrap().trim_3, 0);
        // New assignments has updated trim, same allele_id.
        assert_eq!(b.get(Segment::V).unwrap().trim_3, 4);
        assert_eq!(b.get(Segment::V).unwrap().allele_id, AlleleId::new(7));
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
        assert_eq!(b.get(Segment::V).unwrap().trim_5, 0);
        assert_eq!(b.get(Segment::V).unwrap().trim_3, 6);
        assert_eq!(b.get(Segment::J).unwrap().trim_5, 2);
        assert_eq!(b.get(Segment::J).unwrap().trim_3, 0);

        // Original assignments untouched.
        assert_eq!(a.get(Segment::V).unwrap().trim_3, 0);
        assert_eq!(a.get(Segment::J).unwrap().trim_5, 0);
    }

    #[test]
    fn assignments_branching_revisions_are_independent() {
        let base = AlleleAssignments::new()
            .with_assigned(Segment::V, AlleleInstance::new(AlleleId::new(5)));
        let branch_a = base.with_trim(Segment::V, TrimEnd::Three, 1);
        let branch_b = base.with_trim(Segment::V, TrimEnd::Three, 9);

        assert_eq!(base.get(Segment::V).unwrap().trim_3, 0);
        assert_eq!(branch_a.get(Segment::V).unwrap().trim_3, 1);
        assert_eq!(branch_b.get(Segment::V).unwrap().trim_3, 9);
    }

    #[test]
    fn assignments_assign_replaces_previous_instance() {
        let a = AlleleAssignments::new()
            .with_assigned(Segment::V, AlleleInstance::new(AlleleId::new(1)));
        let b = a.with_assigned(Segment::V, AlleleInstance::new(AlleleId::new(2)));

        assert_eq!(a.get(Segment::V).unwrap().allele_id, AlleleId::new(1));
        assert_eq!(b.get(Segment::V).unwrap().allele_id, AlleleId::new(2));
    }

    // ── SegmentOrientation + AlleleInstance.orientation ──────────────

    #[test]
    fn segment_orientation_default_is_forward() {
        // Pin the `Default` impl: the enum's #[default] is the
        // load-bearing knob that keeps every existing
        // `AlleleInstance::new(...)` byte-identical.
        assert_eq!(SegmentOrientation::default(), SegmentOrientation::Forward);
    }

    #[test]
    fn segment_orientation_is_reverse_only_for_reverse_complement() {
        assert!(!SegmentOrientation::Forward.is_reverse());
        assert!(SegmentOrientation::ReverseComplement.is_reverse());
    }

    #[test]
    fn new_allele_instance_defaults_to_forward_orientation() {
        let inst = AlleleInstance::new(AlleleId::new(7));
        assert_eq!(inst.orientation, SegmentOrientation::Forward);
        assert!(!inst.orientation.is_reverse());
    }

    #[test]
    fn with_orientation_persists_without_touching_other_fields() {
        let a = AlleleInstance::new(AlleleId::new(3))
            .with_trim_5(2)
            .with_trim_3(5);
        let b = a.with_orientation(SegmentOrientation::ReverseComplement);
        // Receiver untouched (persistent IR contract).
        assert_eq!(a.orientation, SegmentOrientation::Forward);
        assert_eq!(b.orientation, SegmentOrientation::ReverseComplement);
        // Trims survive the orientation update.
        assert_eq!(b.trim_5, 2);
        assert_eq!(b.trim_3, 5);
        assert_eq!(b.allele_id, AlleleId::new(3));
    }

    #[test]
    fn with_orientation_round_trip_back_to_forward() {
        let a = AlleleInstance::new(AlleleId::new(0))
            .with_orientation(SegmentOrientation::ReverseComplement)
            .with_orientation(SegmentOrientation::Forward);
        assert_eq!(a.orientation, SegmentOrientation::Forward);
    }

    #[test]
    fn with_trim_setters_preserve_orientation() {
        let a = AlleleInstance::new(AlleleId::new(0))
            .with_orientation(SegmentOrientation::ReverseComplement);
        let b = a.with_trim_5(4);
        let c = b.with_trim_3(6);
        // Orientation survives both trim updates.
        assert_eq!(b.orientation, SegmentOrientation::ReverseComplement);
        assert_eq!(c.orientation, SegmentOrientation::ReverseComplement);
        assert_eq!(c.trim_5, 4);
        assert_eq!(c.trim_3, 6);
    }

    // ── AlleleAssignments::with_orientation ──────────────────────────

    #[test]
    fn assignments_with_orientation_updates_only_target_segment() {
        let a = AlleleAssignments::new()
            .with_assigned(Segment::V, AlleleInstance::new(AlleleId::new(0)))
            .with_assigned(Segment::D, AlleleInstance::new(AlleleId::new(1)))
            .with_assigned(Segment::J, AlleleInstance::new(AlleleId::new(2)));
        let b = a.with_orientation(Segment::D, SegmentOrientation::ReverseComplement);

        // D flipped; V and J unchanged.
        assert_eq!(
            b.get(Segment::D).unwrap().orientation,
            SegmentOrientation::ReverseComplement,
        );
        assert_eq!(
            b.get(Segment::V).unwrap().orientation,
            SegmentOrientation::Forward,
        );
        assert_eq!(
            b.get(Segment::J).unwrap().orientation,
            SegmentOrientation::Forward,
        );
        // Receiver still has every segment in Forward.
        assert_eq!(
            a.get(Segment::D).unwrap().orientation,
            SegmentOrientation::Forward,
        );
    }

    #[test]
    fn assignments_with_orientation_preserves_trims_and_allele_id() {
        let v = AlleleInstance::new(AlleleId::new(11))
            .with_trim_5(3)
            .with_trim_3(7);
        let a = AlleleAssignments::new().with_assigned(Segment::V, v);
        let b = a.with_orientation(Segment::V, SegmentOrientation::ReverseComplement);
        let updated = b.get(Segment::V).unwrap();
        assert_eq!(updated.orientation, SegmentOrientation::ReverseComplement);
        assert_eq!(updated.allele_id, AlleleId::new(11));
        assert_eq!(updated.trim_5, 3);
        assert_eq!(updated.trim_3, 7);
    }

    #[test]
    #[should_panic(expected = "no instance assigned to segment")]
    fn assignments_with_orientation_panics_on_unassigned_segment() {
        // Mirror of `with_updated_panics_on_unassigned` — orientation
        // setter must surface the same wording so a caller misusing
        // either setter sees a uniform diagnostic.
        let a = AlleleAssignments::new();
        let _ = a.with_orientation(Segment::D, SegmentOrientation::ReverseComplement);
    }

    #[test]
    fn assignments_orientation_branching_revisions_are_independent() {
        let base = AlleleAssignments::new()
            .with_assigned(Segment::D, AlleleInstance::new(AlleleId::new(5)));
        let branch_a = base.with_orientation(Segment::D, SegmentOrientation::Forward);
        let branch_b =
            base.with_orientation(Segment::D, SegmentOrientation::ReverseComplement);

        // Base unchanged.
        assert_eq!(
            base.get(Segment::D).unwrap().orientation,
            SegmentOrientation::Forward,
        );
        // Branches independent.
        assert_eq!(
            branch_a.get(Segment::D).unwrap().orientation,
            SegmentOrientation::Forward,
        );
        assert_eq!(
            branch_b.get(Segment::D).unwrap().orientation,
            SegmentOrientation::ReverseComplement,
        );
    }
}
