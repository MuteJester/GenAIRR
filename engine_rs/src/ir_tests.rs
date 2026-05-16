    use super::*;
    use std::mem::size_of;

    #[test]
    fn nucleotide_size_is_small() {
        // We want Nucleotide to stay tight enough for cache-friendly
        // arena scans. 8 bytes today; lifting the cap to 16 is fine
        // but should be a deliberate decision documented in the design
        // doc, not silent struct growth.
        assert!(
            size_of::<Nucleotide>() <= 16,
            "Nucleotide grew past 16 bytes ({} bytes) — review IR layout",
            size_of::<Nucleotide>()
        );
    }

    #[test]
    fn handle_types_are_zero_cost() {
        // Newtype wrappers around u32 must not add overhead.
        assert_eq!(size_of::<NucHandle>(), size_of::<u32>());
        assert_eq!(size_of::<RegionHandle>(), size_of::<u32>());
    }

    #[test]
    fn handle_index_round_trip() {
        let h = NucHandle::new(42);
        assert_eq!(h.index(), 42);
        assert_eq!(h.as_usize(), 42);
    }

    #[test]
    fn segment_equality_is_structural() {
        assert_eq!(Segment::V, Segment::V);
        assert_ne!(Segment::V, Segment::J);
        assert_ne!(Segment::Np1, Segment::Np2);
    }

    #[test]
    fn nuc_flags_set_clear_test() {
        let empty = NucFlags::empty();
        assert!(!empty.contains(flag::P_NUC));

        let with_p = empty.with(flag::P_NUC);
        assert!(with_p.contains(flag::P_NUC));
        assert!(!with_p.contains(flag::N_NUC));

        let with_p_and_junction = with_p.with(flag::JUNCTION);
        assert!(with_p_and_junction.contains(flag::P_NUC));
        assert!(with_p_and_junction.contains(flag::JUNCTION));

        let cleared = with_p_and_junction.without(flag::P_NUC);
        assert!(!cleared.contains(flag::P_NUC));
        assert!(cleared.contains(flag::JUNCTION));
    }

    #[test]
    fn nucleotide_germline_constructor() {
        let n = Nucleotide::germline(b'A', 12, Segment::V);
        assert_eq!(n.base, b'A');
        assert_eq!(n.germline, b'A');
        assert_eq!(n.germline_pos, GermlinePos::pos(12));
        assert_eq!(n.segment, Segment::V);
        assert_eq!(n.flags, NucFlags::empty());
    }

    #[test]
    fn nucleotide_synthetic_constructor_has_no_germline_pos() {
        let n = Nucleotide::synthetic(b'a', Segment::Np1, flag::N_NUC);
        assert_eq!(n.base, b'a');
        assert_eq!(n.germline, b'a');
        assert!(n.germline_pos.is_none());
        assert_eq!(n.segment, Segment::Np1);
        assert!(n.flags.contains(flag::N_NUC));
    }

    #[test]
    fn pool_push_returns_sequential_handles() {
        let mut pool = NucleotidePool::new();
        assert!(pool.is_empty());

        let h0 = pool.push(Nucleotide::germline(b'A', 0, Segment::V));
        let h1 = pool.push(Nucleotide::germline(b'C', 1, Segment::V));
        let h2 = pool.push(Nucleotide::germline(b'G', 2, Segment::V));

        assert_eq!(h0.index(), 0);
        assert_eq!(h1.index(), 1);
        assert_eq!(h2.index(), 2);
        assert_eq!(pool.len(), 3);
        assert!(!pool.is_empty());
    }

    #[test]
    fn pool_get_by_handle_returns_correct_nucleotide() {
        let mut pool = NucleotidePool::new();
        let h = pool.push(Nucleotide::germline(b'T', 5, Segment::J));

        let got = pool.get(h).expect("handle should be valid");
        assert_eq!(got.base, b'T');
        assert_eq!(got.germline_pos, GermlinePos::pos(5));
        assert_eq!(got.segment, Segment::J);
    }

    #[test]
    fn pool_get_out_of_bounds_returns_none() {
        let pool = NucleotidePool::new();
        assert!(pool.get(NucHandle::new(0)).is_none());
        assert!(pool.get(NucHandle::new(99)).is_none());
    }

    #[test]
    fn region_range_and_length() {
        let r = Region::new(Segment::V, NucHandle::new(0), NucHandle::new(10));
        assert_eq!(r.segment, Segment::V);
        assert_eq!(r.len(), 10);
        assert!(!r.is_empty());

        let empty = Region::new(Segment::Np1, NucHandle::new(5), NucHandle::new(5));
        assert!(empty.is_empty());
        assert_eq!(empty.len(), 0);
    }

    #[test]
    fn sequence_starts_empty() {
        let s = Sequence::new();
        assert_eq!(s.region_count(), 0);
        assert!(s.regions.is_empty());
    }

    #[test]
    fn simulation_default_construction_is_clean() {
        let s = Simulation::new();
        assert!(s.pool.is_empty());
        assert_eq!(s.sequence.region_count(), 0);
        assert!(s.live_calls.is_none());
    }

    #[test]
    fn simulation_with_capacity_pre_allocates() {
        // Reserved capacity is a perf hint, not a guarantee of `len`.
        let s = Simulation::with_capacity(500);
        assert!(s.pool.is_empty());
        assert_eq!(s.sequence.region_count(), 0);
        assert!(s.live_calls.is_none());
    }

    #[test]
    fn simulation_live_call_sidecar_is_persistent_and_dormant() {
        let live = crate::live_call::LiveCallState::empty().with_dirty_window(
            crate::live_call::DirtyWindow::new(
                1,
                4,
                crate::live_call::DirtyReason::BaseEdited { site: 2 },
            ),
        );
        let s0 = Simulation::new().with_live_calls(live.clone());
        let (s1, handle) = s0.with_nucleotide_pushed(Nucleotide::germline(b'A', 0, Segment::V));
        let s2 = s1.with_base_changed(handle, b'C');
        let s3 = s2.with_region_added(Region::new(
            Segment::V,
            NucHandle::new(0),
            NucHandle::new(1),
        ));

        assert_eq!(s0.live_calls.as_ref(), Some(&live));
        assert_eq!(s1.live_calls.as_ref(), Some(&live));
        assert_eq!(s2.live_calls.as_ref(), Some(&live));
        assert_eq!(s3.live_calls.as_ref(), Some(&live));

        // Dormant means core structural edits do not yet interpret or
        // mutate live calls; they only preserve the sidecar.
        assert_eq!(s3.pool.len(), 1);
        assert_eq!(s3.sequence.region_count(), 1);
    }

    #[test]
    fn simulation_clone_is_deep_at_value_level() {
        // simple Clone is a real deep copy.
        let mut a = Simulation::new();
        a.pool.push(Nucleotide::germline(b'A', 0, Segment::V));
        a.sequence.regions.push(Region::new(
            Segment::V,
            NucHandle::new(0),
            NucHandle::new(1),
        ));

        let b = a.clone();
        assert_eq!(b.pool.len(), 1);
        assert_eq!(b.sequence.region_count(), 1);

        // Mutating `a` must not affect `b` (Clone is a deep copy in A.2).
        a.pool.push(Nucleotide::germline(b'C', 1, Segment::V));
        assert_eq!(a.pool.len(), 2);
        assert_eq!(b.pool.len(), 1);
    }

    // ── Persistent update API tests ─────────────────────────────────

    #[test]
    fn pool_with_pushed_returns_new_pool_old_unchanged() {
        let p0 = NucleotidePool::new();
        let (p1, h) = p0.with_pushed(Nucleotide::germline(b'A', 0, Segment::V));

        // Old pool: untouched.
        assert_eq!(p0.len(), 0);
        // New pool: has the nucleotide at handle 0.
        assert_eq!(p1.len(), 1);
        assert_eq!(h.index(), 0);
        assert_eq!(p1.get(h).unwrap().base, b'A');
        // Old pool still does not see it.
        assert!(p0.get(h).is_none());
    }

    #[test]
    fn pool_with_base_changed_isolates_revisions() {
        let p0 = NucleotidePool::new();
        let (p1, h) = p0.with_pushed(Nucleotide::germline(b'A', 0, Segment::V));
        let p2 = p1.with_base_changed(h, b'G');

        // Three independent revisions of the pool exist simultaneously.
        assert!(p0.get(h).is_none());
        assert_eq!(p1.get(h).unwrap().base, b'A');
        assert_eq!(p2.get(h).unwrap().base, b'G');

        // Other fields preserved on the changed revision.
        let n2 = p2.get(h).unwrap();
        assert_eq!(n2.germline, b'A'); // germline unchanged
        assert_eq!(n2.germline_pos, GermlinePos::pos(0));
        assert_eq!(n2.segment, Segment::V);
    }

    #[test]
    fn pool_with_nucleotide_changed_replaces_only_target() {
        let mut p = NucleotidePool::new();
        let h0 = p.push(Nucleotide::germline(b'A', 0, Segment::V));
        let h1 = p.push(Nucleotide::germline(b'C', 1, Segment::V));
        let h2 = p.push(Nucleotide::germline(b'G', 2, Segment::V));

        let p2 = p.with_nucleotide_changed(h1, Nucleotide::germline(b'T', 99, Segment::J));

        // Old pool retains everything.
        assert_eq!(p.get(h0).unwrap().base, b'A');
        assert_eq!(p.get(h1).unwrap().base, b'C');
        assert_eq!(p.get(h1).unwrap().segment, Segment::V);
        assert_eq!(p.get(h2).unwrap().base, b'G');

        // New pool: target replaced, neighbours unchanged.
        assert_eq!(p2.get(h0).unwrap().base, b'A');
        assert_eq!(p2.get(h1).unwrap().base, b'T');
        assert_eq!(p2.get(h1).unwrap().germline_pos, GermlinePos::pos(99));
        assert_eq!(p2.get(h1).unwrap().segment, Segment::J);
        assert_eq!(p2.get(h2).unwrap().base, b'G');
    }

    #[test]
    fn region_with_end_extended_keeps_other_fields() {
        let r0 = Region::new(Segment::V, NucHandle::new(0), NucHandle::new(5));
        let r1 = r0.with_end_extended(NucHandle::new(10));

        assert_eq!(r0.end.index(), 5);
        assert_eq!(r1.end.index(), 10);
        assert_eq!(r1.start.index(), 0);
        assert_eq!(r1.segment, Segment::V);
    }

    #[test]
    fn region_with_frame_phase_isolated() {
        let r0 = Region::new(Segment::V, NucHandle::new(0), NucHandle::new(9));
        assert_eq!(r0.frame_phase, 0);

        let r1 = r0.with_frame_phase(2);
        assert_eq!(r0.frame_phase, 0);
        assert_eq!(r1.frame_phase, 2);
        assert_eq!(r1.start, r0.start);
        assert_eq!(r1.end, r0.end);
        assert_eq!(r1.segment, r0.segment);
    }

    #[test]
    fn sequence_with_region_added_isolated() {
        let s0 = Sequence::new();
        let s1 = s0.with_region_added(Region::new(
            Segment::V,
            NucHandle::new(0),
            NucHandle::new(5),
        ));
        let s2 = s1.with_region_added(Region::new(
            Segment::J,
            NucHandle::new(5),
            NucHandle::new(8),
        ));

        assert_eq!(s0.region_count(), 0);
        assert_eq!(s1.region_count(), 1);
        assert_eq!(s2.region_count(), 2);
        assert_eq!(s2.regions[0].segment, Segment::V);
        assert_eq!(s2.regions[1].segment, Segment::J);
    }

    #[test]
    fn sequence_with_region_replaced_isolated() {
        let s = Sequence::new()
            .with_region_added(Region::new(
                Segment::V,
                NucHandle::new(0),
                NucHandle::new(3),
            ))
            .with_region_added(Region::new(
                Segment::J,
                NucHandle::new(3),
                NucHandle::new(6),
            ));

        let s2 = s.with_region_replaced(
            0,
            Region::new(Segment::V, NucHandle::new(0), NucHandle::new(10)),
        );

        // Old sequence: V region had length 3.
        assert_eq!(s.regions[0].len(), 3);
        // New sequence: V region has length 10, J unchanged.
        assert_eq!(s2.regions[0].len(), 10);
        assert_eq!(s2.regions[1].segment, Segment::J);
        assert_eq!(s2.regions[1].len(), 3);
    }

    #[test]
    fn simulation_persistent_chain_preserves_history() {
        // The point of D1: every revision in a chain remains accessible.
        let s0 = Simulation::new();
        let (s1, h_a) = s0.with_nucleotide_pushed(Nucleotide::germline(b'A', 0, Segment::V));
        let (s2, h_c) = s1.with_nucleotide_pushed(Nucleotide::germline(b'C', 1, Segment::V));
        let s3 = s2.with_base_changed(h_a, b'T');
        let s4 = s3.with_region_added(Region::new(
            Segment::V,
            NucHandle::new(0),
            NucHandle::new(2),
        ));

        // s0: empty.
        assert_eq!(s0.pool.len(), 0);
        assert_eq!(s0.sequence.region_count(), 0);

        // s1: one nucleotide (A).
        assert_eq!(s1.pool.len(), 1);
        assert_eq!(s1.pool.get(h_a).unwrap().base, b'A');
        assert_eq!(s1.sequence.region_count(), 0);

        // s2: two nucleotides (A, C), no regions.
        assert_eq!(s2.pool.len(), 2);
        assert_eq!(s2.pool.get(h_a).unwrap().base, b'A');
        assert_eq!(s2.pool.get(h_c).unwrap().base, b'C');
        assert_eq!(s2.sequence.region_count(), 0);

        // s3: A→T mutation, C unchanged.
        assert_eq!(s3.pool.len(), 2);
        assert_eq!(s3.pool.get(h_a).unwrap().base, b'T');
        assert_eq!(s3.pool.get(h_c).unwrap().base, b'C');
        // Original germline preserved on the mutated nucleotide.
        assert_eq!(s3.pool.get(h_a).unwrap().germline, b'A');
        assert_eq!(s3.sequence.region_count(), 0);

        // s4: region added.
        assert_eq!(s4.pool.len(), 2);
        assert_eq!(s4.sequence.region_count(), 1);
        assert_eq!(s4.sequence.regions[0].segment, Segment::V);

        // Earlier revisions still see their own pool state.
        assert_eq!(s2.pool.get(h_a).unwrap().base, b'A'); // pre-mutation
        assert_eq!(s3.sequence.region_count(), 0); // pre-region-add
    }

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

        let s0 = Simulation::new()
            .with_allele_assigned(Segment::V, AlleleInstance::new(AlleleId::new(0)));
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
            let (next, _) =
                sim.with_nucleotide_pushed(Nucleotide::germline(*b, i as u16, Segment::V));
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
            let (next, _) =
                sim.with_nucleotide_pushed(Nucleotide::germline(*b, i as u16, Segment::V));
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
            let (next, _) =
                sim.with_nucleotide_pushed(Nucleotide::germline(*b, i as u16, Segment::V));
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
        let v_region = Region::new(Segment::V, NucHandle::new(0), NucHandle::new(2))
            .with_frame_phase(0)
            .with_codon_rail_recomputed(&sim.pool);
        let j_region = Region::new(Segment::J, NucHandle::new(2), NucHandle::new(8))
            .with_frame_phase(2)
            .with_codon_rail_recomputed(&sim.pool);
        sim = sim.with_region_added(v_region).with_region_added(j_region);

        let mutated = sim.with_indel_inserted(
            1,
            Nucleotide::synthetic(b'G', Segment::V, crate::ir::flag::INDEL_INSERTED),
        );

        let j = &mutated.sequence.regions[1];
        assert_eq!(j.start, NucHandle::new(3));
        assert_eq!(j.end, NucHandle::new(9));
        assert_eq!(j.frame_phase, 0);
        assert_eq!(j.amino_acids, b"MP");
    }

    #[test]
    fn simulation_with_indel_deleted_refreshes_downstream_frame_phase_and_rail() {
        let mut sim = Simulation::new();
        for (i, b) in b"AAAAATGCCC".iter().enumerate() {
            let segment = if i < 4 { Segment::V } else { Segment::J };
            let (next, _) = sim.with_nucleotide_pushed(Nucleotide::germline(*b, i as u16, segment));
            sim = next;
        }
        let v_region = Region::new(Segment::V, NucHandle::new(0), NucHandle::new(4))
            .with_frame_phase(0)
            .with_codon_rail_recomputed(&sim.pool);
        let j_region = Region::new(Segment::J, NucHandle::new(4), NucHandle::new(10))
            .with_frame_phase(1)
            .with_codon_rail_recomputed(&sim.pool);
        sim = sim.with_region_added(v_region).with_region_added(j_region);

        let mutated = sim.with_indel_deleted(1);

        let j = &mutated.sequence.regions[1];
        assert_eq!(j.start, NucHandle::new(3));
        assert_eq!(j.end, NucHandle::new(9));
        assert_eq!(j.frame_phase, 0);
        assert_eq!(j.amino_acids, b"MP");
    }

    #[test]
    fn simulation_with_indel_inserted_full_round_trip() {
        let mut sim = Simulation::new();
        for (i, b) in b"ATGCCCGGG".iter().enumerate() {
            let (next, _) =
                sim.with_nucleotide_pushed(Nucleotide::germline(*b, i as u16, Segment::V));
            sim = next;
        }
        let region = Region::new(Segment::V, NucHandle::new(0), NucHandle::new(9))
            .with_codon_rail_recomputed(&sim.pool);
        sim = sim.with_region_added(region);
        // Initial codon rail: ATG CCC GGG → M P G.
        assert_eq!(sim.sequence.regions[0].amino_acids, b"MPG");

        // Insert a 'G' at position 3 — the region grows by 1, and
        // codons shift: ATG | GCC CGG G → M, A, R, + leftover G.
        let mutated = sim.with_indel_inserted(
            3,
            Nucleotide::synthetic(b'G', Segment::V, crate::ir::flag::INDEL_INSERTED),
        );
        assert_eq!(mutated.pool.len(), 10);
        assert_eq!(mutated.sequence.regions[0].len(), 10);
        // ATG GCC CGG → M A R (3 codons, 10 = 3*3 + 1 leftover)
        assert_eq!(mutated.sequence.regions[0].amino_acids, b"MAR");
        // Original sim untouched.
        assert_eq!(sim.pool.len(), 9);
        assert_eq!(sim.sequence.regions[0].amino_acids, b"MPG");
    }

    #[test]
    fn simulation_with_indel_deleted_full_round_trip() {
        let mut sim = Simulation::new();
        for (i, b) in b"ATGAAACCCGGG".iter().enumerate() {
            let (next, _) =
                sim.with_nucleotide_pushed(Nucleotide::germline(*b, i as u16, Segment::V));
            sim = next;
        }
        let region = Region::new(Segment::V, NucHandle::new(0), NucHandle::new(12))
            .with_codon_rail_recomputed(&sim.pool);
        sim = sim.with_region_added(region);
        // Codon rail: ATG AAA CCC GGG → M K P G.
        assert_eq!(sim.sequence.regions[0].amino_acids, b"MKPG");

        // Delete position 3 (the first 'A' of the second codon).
        // Region shrinks to 11 bases.
        let mutated = sim.with_indel_deleted(3);
        assert_eq!(mutated.pool.len(), 11);
        assert_eq!(mutated.sequence.regions[0].len(), 11);
        // Bases now: ATG AAC CCG GG → 3 full codons + 2 leftover.
        // ATG=M, AAC=N, CCG=P. Three amino acids.
        assert_eq!(mutated.sequence.regions[0].amino_acids, b"MNP");
    }

    #[test]
    fn simulation_indel_branching_revisions_are_independent() {
        let mut sim = Simulation::new();
        for (i, b) in b"ATGAAA".iter().enumerate() {
            let (next, _) =
                sim.with_nucleotide_pushed(Nucleotide::germline(*b, i as u16, Segment::V));
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
            let (next, _) =
                sim.with_nucleotide_pushed(Nucleotide::germline(*b, i as u16, Segment::V));
            sim = next;
        }
        let region = Region::new(Segment::V, NucHandle::new(0), NucHandle::new(6))
            .with_codon_rail_recomputed(&sim.pool);
        sim = sim.with_region_added(region);
        assert_eq!(sim.sequence.regions[0].amino_acids, b"KP");

        // Replace position 3 with a new nucleotide whose base is 'T'.
        // Codon at [3..6) was CCC = P, becomes TCC = S.
        let mutated = sim
            .with_nucleotide_changed(NucHandle::new(3), Nucleotide::germline(b'T', 3, Segment::V));
        assert_eq!(mutated.sequence.regions[0].amino_acids, b"KS");
    }

    #[test]
    fn simulation_assignments_chain_with_pool_changes() {
        // Assignments survive pool / sequence changes — they live in
        // their own slot, not in the pool.
        use crate::assignment::AlleleInstance;
        use crate::refdata::AlleleId;

        let s0 = Simulation::new()
            .with_allele_assigned(Segment::V, AlleleInstance::new(AlleleId::new(3)));

        let (s1, _h) = s0.with_nucleotide_pushed(Nucleotide::germline(b'A', 0, Segment::V));

        // Assignment survives across the pool change.
        assert_eq!(s1.assignments.v.unwrap().allele_id, AlleleId::new(3));
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

    // ── Codon-rail metadata tests ───────────────────────────────────

    /// Helper: build a pool from a base string with all nucleotides in
    /// segment V, germline_pos = position within the string. Returns
    /// the pool and a Region covering the whole pool with frame_phase=0.
    fn pool_from_string(s: &str) -> (NucleotidePool, Region) {
        let mut pool = NucleotidePool::new();
        for (i, b) in s.bytes().enumerate() {
            pool.push(Nucleotide::germline(b, i as u16, Segment::V));
        }
        let region = Region::new(
            Segment::V,
            NucHandle::new(0),
            NucHandle::new(s.len() as u32),
        );
        (pool, region)
    }

    #[test]
    fn translate_codon_canonical_examples() {
        assert_eq!(translate_codon(b'A', b'T', b'G'), b'M'); // Met (start)
        assert_eq!(translate_codon(b'T', b'T', b'T'), b'F'); // Phe
        assert_eq!(translate_codon(b'T', b'G', b'G'), b'W'); // Trp
        assert_eq!(translate_codon(b'C', b'A', b'C'), b'H'); // His
        assert_eq!(translate_codon(b'A', b'A', b'A'), b'K'); // Lys
        assert_eq!(translate_codon(b'G', b'G', b'C'), b'G'); // Gly
    }

    #[test]
    fn translate_codon_stops() {
        assert_eq!(translate_codon(b'T', b'A', b'A'), AMINO_STOP);
        assert_eq!(translate_codon(b'T', b'A', b'G'), AMINO_STOP);
        assert_eq!(translate_codon(b'T', b'G', b'A'), AMINO_STOP);
        // TGG is *not* a stop — it is Trp (W).
        assert_eq!(translate_codon(b'T', b'G', b'G'), b'W');
    }

    #[test]
    fn translate_codon_lowercase_is_case_insensitive() {
        assert_eq!(translate_codon(b'a', b't', b'g'), b'M');
        assert_eq!(translate_codon(b'T', b'a', b'a'), AMINO_STOP);
    }

    #[test]
    fn translate_codon_ambiguous_is_x() {
        assert_eq!(translate_codon(b'A', b'T', b'N'), AMINO_AMBIGUOUS);
        assert_eq!(translate_codon(b'N', b'A', b'T'), AMINO_AMBIGUOUS);
        assert_eq!(translate_codon(b'A', b'-', b'C'), AMINO_AMBIGUOUS);
    }

    #[test]
    fn translate_codon_uracil_treated_as_thymine() {
        assert_eq!(translate_codon(b'A', b'U', b'G'), b'M');
        assert_eq!(translate_codon(b'U', b'A', b'A'), AMINO_STOP);
    }

    #[test]
    fn region_codon_rail_starts_empty_until_recomputed() {
        let r = Region::new(Segment::V, NucHandle::new(0), NucHandle::new(9));
        assert!(r.amino_acids.is_empty());
        assert!(r.stop_codon_positions.is_empty());
        assert_eq!(r.stop_codon_count(), 0);
    }

    #[test]
    fn region_with_codon_rail_recomputed_translates_in_frame() {
        // ATG GGG CAC → M G H, no stops.
        let (pool, region) = pool_from_string("ATGGGGCAC");
        let r = region.with_codon_rail_recomputed(&pool);

        assert_eq!(r.amino_acids, b"MGH");
        assert_eq!(r.stop_codon_positions, vec![]);
        assert_eq!(r.stop_codon_count(), 0);
        // Receiver unchanged (persistent contract).
        assert!(region.amino_acids.is_empty());
    }

    #[test]
    fn region_with_codon_rail_recomputed_picks_up_stops() {
        // ATG TAA TGG → M * W, one stop at position 3.
        let (pool, region) = pool_from_string("ATGTAATGG");
        let r = region.with_codon_rail_recomputed(&pool);

        assert_eq!(r.amino_acids, b"M*W");
        assert_eq!(r.stop_codon_positions, vec![NucHandle::new(3)]);
        assert_eq!(r.stop_codon_count(), 1);
    }

    #[test]
    fn region_with_codon_rail_recomputed_drops_partial_codon_at_end() {
        // ATG GG → M, plus 2 incomplete bases. No second amino acid.
        let (pool, region) = pool_from_string("ATGGG");
        let r = region.with_codon_rail_recomputed(&pool);

        assert_eq!(r.amino_acids, b"M");
        assert_eq!(r.stop_codon_count(), 0);
    }

    #[test]
    fn region_with_codon_rail_respects_frame_phase_one() {
        // frame_phase=1 means position 0 is the 2nd base of a codon
        // started in a (notional) previous region. Skip 2 bases, then
        // translate.
        //
        // Bases:   X X A T G C C C
        // Phase:   1 2 0 1 2 0 1 2     (0 = first base of codon)
        // Codons fully in this region: ATG, CCC → M P
        let (pool, region) = pool_from_string("XXATGCCC");
        let r = region.with_frame_phase(1).with_codon_rail_recomputed(&pool);
        assert_eq!(r.amino_acids, b"MP");
    }

    #[test]
    fn region_with_codon_rail_respects_frame_phase_two() {
        // frame_phase=2 means position 0 is the 3rd base of a codon.
        // Skip 1 base, then translate.
        //
        // Bases:   X T A C G G G
        // Phase:   2 0 1 2 0 1 2
        // Codons fully in this region: TAC, GGG → Y G
        let (pool, region) = pool_from_string("XTACGGG");
        let r = region.with_frame_phase(2).with_codon_rail_recomputed(&pool);
        assert_eq!(r.amino_acids, b"YG");
    }

    #[test]
    fn region_with_codon_rail_handles_ambiguous_bases() {
        // Codon containing N → X.
        let (pool, region) = pool_from_string("ATGNAATGG");
        let r = region.with_codon_rail_recomputed(&pool);
        // ATG = M, NAA = X (ambiguous), TGG = W
        assert_eq!(r.amino_acids, b"MXW");
    }

    #[test]
    fn region_with_codon_rail_recomputes_after_base_mutation() {
        // Mutation in this region must produce a region revision with
        // updated amino acids.
        let (pool0, region0) = pool_from_string("ATGTACTGG");
        let r0 = region0.with_codon_rail_recomputed(&pool0);
        assert_eq!(r0.amino_acids, b"MYW");

        // Mutate the second base of the second codon: TAC → TGC = C.
        let pool1 = pool0.with_base_changed(NucHandle::new(4), b'G');
        let r1 = region0.with_codon_rail_recomputed(&pool1);
        assert_eq!(r1.amino_acids, b"MCW");

        // Old region revision unaffected.
        assert_eq!(r0.amino_acids, b"MYW");
    }

    #[test]
    fn region_with_codon_rail_detects_stop_introduced_by_mutation() {
        // TAC → TAA (Y → stop) by mutating one base.
        let (pool0, region0) = pool_from_string("ATGTACGGG");
        let r0 = region0.with_codon_rail_recomputed(&pool0);
        assert_eq!(r0.stop_codon_count(), 0);
        assert_eq!(r0.amino_acids, b"MYG");

        // Mutate position 5: TAC → TAA. This creates a stop codon.
        let pool1 = pool0.with_base_changed(NucHandle::new(5), b'A');
        let r1 = region0.with_codon_rail_recomputed(&pool1);
        assert_eq!(r1.amino_acids, b"M*G");
        assert_eq!(r1.stop_codon_count(), 1);
        assert_eq!(r1.stop_codon_positions, vec![NucHandle::new(3)]);

        // Old revision unchanged.
        assert_eq!(r0.stop_codon_count(), 0);
    }

    #[test]
    fn region_codon_rail_empty_region_produces_empty_metadata() {
        let pool = NucleotidePool::new();
        let region = Region::new(Segment::Np1, NucHandle::new(0), NucHandle::new(0));
        let r = region.with_codon_rail_recomputed(&pool);
        assert!(r.amino_acids.is_empty());
        assert!(r.stop_codon_positions.is_empty());
    }

    #[test]
    fn region_codon_rail_frame_phase_2_on_one_base_region() {
        // frame_phase = 2 means the region's first base is the third
        // base of a codon started in a previous region. `skip = 1`,
        // so we'd need at least 4 bases (1 skipped + 3 for a fresh
        // codon) to emit anything. With one base, output is empty.
        let (pool, region) = pool_from_string("X");
        let r = region.with_frame_phase(2).with_codon_rail_recomputed(&pool);
        assert!(r.amino_acids.is_empty());
        assert!(r.stop_codon_positions.is_empty());
    }

    #[test]
    fn region_codon_rail_frame_phase_1_on_two_base_region() {
        // frame_phase = 1 → skip 2 → no bases left after skip.
        // Output is empty.
        let (pool, region) = pool_from_string("XY");
        let r = region.with_frame_phase(1).with_codon_rail_recomputed(&pool);
        assert!(r.amino_acids.is_empty());
    }

    #[test]
    fn region_codon_rail_malformed_end_less_than_start_is_safe() {
        // Defensive contract: a region constructed with end < start
        // (which `len()` already saturates to 0) produces empty
        // codon rail metadata, not a panic. Builders should not
        // produce such regions, but the recompute should be robust.
        let mut pool = NucleotidePool::new();
        for i in 0..6 {
            pool.push(Nucleotide::germline(b'A', i, Segment::V));
        }
        let region = Region::new(Segment::V, NucHandle::new(5), NucHandle::new(2));
        assert_eq!(region.len(), 0); // saturating_sub clamps to 0

        let r = region.with_codon_rail_recomputed(&pool);
        assert!(r.amino_acids.is_empty());
        assert!(r.stop_codon_positions.is_empty());
    }

    #[test]
    fn region_codon_rail_skip_overruns_end_is_safe() {
        // 1-base region with frame_phase=1 → skip=2 → start_idx > end_idx
        // after the skip. Same defensive case as above, different
        // path through the code.
        let (pool, region) = pool_from_string("X");
        let r = region.with_frame_phase(1).with_codon_rail_recomputed(&pool);
        assert!(r.amino_acids.is_empty());
        assert!(r.stop_codon_positions.is_empty());
    }

    // ── Stress test for the persistent IR + codon rail ──────────────

    /// Tiny deterministic PRNG (xorshift32) for the stress test. We
    /// avoid bringing in `rand` here so the test has zero external
    /// dependencies and reproduces identically across machines.
    struct Xorshift32(u32);

    impl Xorshift32 {
        fn new(seed: u32) -> Self {
            // xorshift32 cannot have a zero seed.
            Self(if seed == 0 { 0xdead_beef } else { seed })
        }
        fn next(&mut self) -> u32 {
            let mut x = self.0;
            x ^= x << 13;
            x ^= x >> 17;
            x ^= x << 5;
            self.0 = x;
            x
        }
    }

    /// Apply 1000 random base mutations through the persistent API,
    /// keeping every IR revision in memory. Verifies four properties
    /// that together pin the IR's architectural commitments:
    ///
    ///   1. Persistent contract (D1): the initial revision retains its
    ///      original amino-acid sequence after all 1000 mutations have
    ///      been applied to descendant revisions. No mutation leaks
    ///      back upstream.
    ///   2. Persistent + entity-attached metadata (D5): every revision's
    ///      stored `amino_acids` field equals the result of recomputing
    ///      the codon rail from that revision's own pool. Stale
    ///      metadata is structurally impossible.
    ///   3. Branching independence: at any mid-revision, the stored
    ///      amino acids match the revision's pool — they were captured
    ///      at construction time and were not perturbed by later writes.
    ///   4. The persistent API is not silently coalescing changes —
    ///      most random mutations produce a visible delta in the new
    ///      revision's pool versus the previous revision's pool.
    #[test]
    fn stress_persistent_ir_with_codon_rail() {
        const N_NUCS: u32 = 90; // 30 codons
        const N_MUTATIONS: usize = 1000;

        // Build the initial pool: 90 'A' nucleotides → 30 codons of AAA → 30 K's.
        let mut pool0 = NucleotidePool::new();
        for i in 0..N_NUCS {
            pool0.push(Nucleotide::germline(b'A', i as u16, Segment::V));
        }
        let region0 = Region::new(Segment::V, NucHandle::new(0), NucHandle::new(N_NUCS))
            .with_codon_rail_recomputed(&pool0);
        let sim0 = Simulation {
            pool: pool0,
            sequence: Sequence::new().with_region_added(region0),
            assignments: crate::assignment::AlleleAssignments::new(),
            live_calls: None,
        };

        // Property check at revision 0.
        assert_eq!(sim0.sequence.regions[0].amino_acids.len(), 30);
        assert!(sim0.sequence.regions[0]
            .amino_acids
            .iter()
            .all(|&aa| aa == b'K'));

        // Snapshot the initial amino acid sequence to verify property
        // (1) — that revision 0 is never perturbed.
        let initial_amino_acids = sim0.sequence.regions[0].amino_acids.clone();

        // Apply 1000 random mutations. Every mutation produces a new
        // Simulation revision; every previous revision stays alive.
        let mut history: Vec<Simulation> = Vec::with_capacity(N_MUTATIONS + 1);
        history.push(sim0);

        let mut rng = Xorshift32::new(0x517c_c1ed); // arbitrary fixed seed
        let bases = [b'A', b'C', b'G', b'T'];
        let mut visible_mutations = 0usize;

        for _ in 0..N_MUTATIONS {
            let prev = history.last().unwrap();
            let pos = rng.next() % N_NUCS;
            let new_base = bases[(rng.next() & 0b11) as usize];

            let prev_base = prev.pool.get(NucHandle::new(pos)).unwrap().base;
            if prev_base != new_base {
                visible_mutations += 1;
            }

            // Persistent update: pool, then region (with codon rail
            // recomputed against the new pool), then sequence.
            let new_pool = prev.pool.with_base_changed(NucHandle::new(pos), new_base);
            let new_region = prev.sequence.regions[0].with_codon_rail_recomputed(&new_pool);
            let new_sequence = prev.sequence.with_region_replaced(0, new_region);

            history.push(Simulation {
                pool: new_pool,
                sequence: new_sequence,
                assignments: crate::assignment::AlleleAssignments::new(),
                live_calls: prev.live_calls.clone(),
            });
        }

        // Sanity: history size is correct.
        assert_eq!(history.len(), N_MUTATIONS + 1);

        // Property (1): revision 0's amino acids are exactly what they
        // were before any mutation happened.
        assert_eq!(
            history[0].sequence.regions[0].amino_acids,
            initial_amino_acids
        );
        assert!(history[0].sequence.regions[0]
            .amino_acids
            .iter()
            .all(|&aa| aa == b'K'));

        // Property (2): every revision's stored amino acids equal a
        // fresh recomputation from its own pool. This is the entity-
        // attached metadata invariant — staleness is impossible.
        for (i, rev) in history.iter().enumerate() {
            let recomputed = rev.sequence.regions[0].with_codon_rail_recomputed(&rev.pool);
            assert_eq!(
                rev.sequence.regions[0].amino_acids, recomputed.amino_acids,
                "revision {} amino acids drift from pool",
                i
            );
            assert_eq!(
                rev.sequence.regions[0].stop_codon_positions, recomputed.stop_codon_positions,
                "revision {} stop codon positions drift from pool",
                i
            );
        }

        // Property (3): mid-revisions hold their own state, not the
        // tail's state. We pick a few specific revisions to confirm.
        for &i in &[1usize, 250, 500, 750, 999] {
            assert!(i < history.len());
            let rev = &history[i];
            let recomputed = rev.sequence.regions[0].with_codon_rail_recomputed(&rev.pool);
            assert_eq!(rev.sequence.regions[0].amino_acids, recomputed.amino_acids);
        }

        // Property (4): most mutations produced a visible base delta.
        // With a uniform 4-base alphabet, the expected fraction of
        // self-substitutions is 1/4, so visible should be ~75% × 1000.
        // Use a generous lower bound to absorb statistical jitter.
        assert!(
            visible_mutations >= 600,
            "expected ≥600 visible mutations, got {}",
            visible_mutations
        );
    }

    // ── GermlinePos newtype invariants ───────────────────────────────

    #[test]
    #[should_panic(expected = "u16::MAX is reserved for NONE")]
    fn germline_pos_pos_rejects_max() {
        // u16::MAX is reserved for the NONE sentinel — constructing a
        // positioned GermlinePos with that value must panic, otherwise
        // a caller could silently create a value that compares equal
        // to NONE under derive(PartialEq).
        let _ = GermlinePos::pos(u16::MAX);
    }

    #[test]
    fn germline_pos_none_projects_to_option_none() {
        assert_eq!(GermlinePos::NONE.get(), None);
        assert!(GermlinePos::NONE.is_none());
        assert!(!GermlinePos::NONE.is_some());
    }

    #[test]
    fn germline_pos_pos_round_trips() {
        let p = GermlinePos::pos(42);
        assert_eq!(p.get(), Some(42));
        assert!(p.is_some());
        assert!(!p.is_none());
    }

    #[test]
    fn nucleotide_size_unchanged() {
        // Pin the layout guarantee from `#[repr(transparent)]` on
        // GermlinePos: Nucleotide stays 6 bytes after the migration.
        // If this fails after a future change, the newtype lost its
        // niche / transparency and the migration regressed memory
        // cost on the hottest data structure in the engine.
        // (Plan's spec said 8 bytes; actual baseline was 6 — both
        // pre- and post-migration measured 6 with the current field
        // set: base/germline/germline_pos/segment/flags. The point
        // of the assertion is unchangedness, not the magic number.)
        assert_eq!(std::mem::size_of::<Nucleotide>(), 6);
    }
