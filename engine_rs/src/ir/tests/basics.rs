use super::*;

#[test]
fn nucleotide_size_is_small() {
    // We want Nucleotide to stay tight enough for cache-friendly
    // arena scans. 6 bytes today; lifting the cap to 16 is fine
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
    assert_eq!(s.segment_calls.version, 0);
    assert!(s.dirty_log.is_empty());
    assert_eq!(s.mutation_count, 0);
}

#[test]
fn simulation_with_capacity_pre_allocates() {
    // Reserved capacity is a perf hint, not a guarantee of `len`.
    let s = Simulation::with_capacity(500);
    assert!(s.pool.is_empty());
    assert_eq!(s.sequence.region_count(), 0);
    assert_eq!(s.segment_calls.version, 0);
    assert!(s.dirty_log.is_empty());
}

#[test]
fn simulation_live_call_sidecar_is_persistent_and_dormant() {
    // Sidecars (segment_calls / dirty_log / mutation_count) propagate
    // unchanged through structural edits — none of `with_*` interpret
    // or mutate them.
    let mut log = crate::live_call::DirtyLog::empty();
    log.push(crate::live_call::DirtyWindow::new(
        1,
        4,
        crate::live_call::DirtyReason::BaseEdited { site: 2 },
    ));
    let s0 = Simulation::new()
        .with_dirty_log(log.clone())
        .with_mutation_count(7);
    let (s1, handle) = s0.with_nucleotide_pushed(Nucleotide::germline(b'A', 0, Segment::V));
    let s2 = s1.with_base_changed(handle, b'C');
    let s3 = s2.with_region_added(Region::new(
        Segment::V,
        NucHandle::new(0),
        NucHandle::new(1),
    ));

    for sim in [&s0, &s1, &s2, &s3] {
        assert_eq!(&*sim.dirty_log, &log);
        assert_eq!(sim.mutation_count, 7);
    }

    // Dormant means core structural edits do not yet interpret or
    // mutate live calls; they only preserve the sidecar.
    assert_eq!(s3.pool.len(), 1);
    assert_eq!(s3.sequence.region_count(), 1);
}

#[test]
fn simulation_clone_is_deep_at_value_level() {
    // Clone is structurally shared via Arc but copy-on-write —
    // mutating `a` after a clone still leaves `b` semantically
    // independent. This test pins that observable invariant.
    let mut a = Simulation::new();
    a.pool.push(Nucleotide::germline(b'A', 0, Segment::V));
    std::sync::Arc::make_mut(&mut a.sequence)
        .regions
        .push(Region::new(
            Segment::V,
            NucHandle::new(0),
            NucHandle::new(1),
        ));

    let b = a.clone();
    assert_eq!(b.pool.len(), 1);
    assert_eq!(b.sequence.region_count(), 1);

    // Mutating `a` must not affect `b` (copy-on-write semantics).
    a.pool.push(Nucleotide::germline(b'C', 1, Segment::V));
    assert_eq!(a.pool.len(), 2);
    assert_eq!(b.pool.len(), 1);
}

// ── Persistent update API tests ─────────────────────────────────
