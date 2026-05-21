use super::*;

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

// ── Batched-push equivalence (Phase 1 IR-allocation refactor) ───

/// `with_nucleotides_extended` must produce a simulation pool that
/// is bit-for-bit identical to one produced by N repeated
/// `with_nucleotide_pushed` calls — same bases, same handles, same
/// metadata. This pins the equivalence the assembly pass relies on.
#[test]
fn simulation_with_nucleotides_extended_matches_repeated_push() {
    let bases: Vec<Nucleotide> = (0..10)
        .map(|i| Nucleotide::germline(b"ACGTACGTAC"[i] as u8, i as u16, Segment::V))
        .collect();

    // Reference: per-base push.
    let mut by_push = Simulation::new();
    let mut push_handles: Vec<NucHandle> = Vec::new();
    for n in &bases {
        let (next, h) = by_push.with_nucleotide_pushed(*n);
        by_push = next;
        push_handles.push(h);
    }

    // Under test: one batched extend.
    let (by_extend, range) = Simulation::new().with_nucleotides_extended(bases.iter().copied());

    // Same pool length.
    assert_eq!(by_push.pool.len(), by_extend.pool.len());
    assert_eq!(by_extend.pool.len(), 10);

    // Range covers exactly the pushed bases [0, 10).
    assert_eq!(range, 0..10);

    // Every nucleotide is identical, including provenance metadata.
    for i in 0..10u32 {
        let h = NucHandle::new(i);
        assert_eq!(by_push.pool.get(h), by_extend.pool.get(h));
    }

    // Handles returned by push form the same contiguous range.
    assert_eq!(push_handles.first().map(|h| h.index()), Some(range.start));
    assert_eq!(push_handles.last().map(|h| h.index() + 1), Some(range.end));
}

/// `with_nucleotides_extended` against a non-empty pool must continue
/// numbering handles from `pool.len()` — same offset semantics as a
/// per-base push loop.
#[test]
fn simulation_with_nucleotides_extended_chains_after_existing_pool() {
    let (s1, _h) =
        Simulation::new().with_nucleotide_pushed(Nucleotide::germline(b'G', 0, Segment::V));
    assert_eq!(s1.pool.len(), 1);

    let extra = vec![
        Nucleotide::germline(b'A', 1, Segment::V),
        Nucleotide::germline(b'T', 2, Segment::V),
    ];
    let (s2, range) = s1.with_nucleotides_extended(extra);

    assert_eq!(s2.pool.len(), 3);
    assert_eq!(range, 1..3);
    assert_eq!(s2.pool.get(NucHandle::new(0)).unwrap().base, b'G');
    assert_eq!(s2.pool.get(NucHandle::new(1)).unwrap().base, b'A');
    assert_eq!(s2.pool.get(NucHandle::new(2)).unwrap().base, b'T');
}

/// Empty extension must be a no-op clone — pool length unchanged,
/// returned range empty at the current pool boundary.
#[test]
fn simulation_with_nucleotides_extended_empty_iter_is_noop() {
    let (s1, _h) =
        Simulation::new().with_nucleotide_pushed(Nucleotide::germline(b'C', 0, Segment::V));
    let len_before = s1.pool.len() as u32;
    let (s2, range) = s1.with_nucleotides_extended(std::iter::empty());

    assert_eq!(s2.pool.len() as u32, len_before);
    assert_eq!(range, len_before..len_before);
}
