use super::*;
use crate::ir::compute_codon_rail;

// ── Codon-rail computation tests ─────────────────────────────────

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
fn compute_codon_rail_translates_in_frame() {
    // ATG GGG CAC → M G H, no stops.
    let (pool, region) = pool_from_string("ATGGGGCAC");
    let rail = compute_codon_rail(&region, &pool);

    assert_eq!(rail.amino_acids, b"MGH");
    assert_eq!(rail.stop_codon_positions, vec![]);
    assert_eq!(rail.stop_codon_count(), 0);
}

#[test]
fn compute_codon_rail_picks_up_stops() {
    // ATG TAA TGG → M * W, one stop at position 3.
    let (pool, region) = pool_from_string("ATGTAATGG");
    let rail = compute_codon_rail(&region, &pool);

    assert_eq!(rail.amino_acids, b"M*W");
    assert_eq!(rail.stop_codon_positions, vec![NucHandle::new(3)]);
    assert_eq!(rail.stop_codon_count(), 1);
}

#[test]
fn compute_codon_rail_drops_partial_codon_at_end() {
    // ATG GG → M, plus 2 incomplete bases. No second amino acid.
    let (pool, region) = pool_from_string("ATGGG");
    let rail = compute_codon_rail(&region, &pool);

    assert_eq!(rail.amino_acids, b"M");
    assert_eq!(rail.stop_codon_count(), 0);
}

#[test]
fn compute_codon_rail_respects_frame_phase_one() {
    // frame_phase=1 means position 0 is the 2nd base of a codon
    // started in a (notional) previous region. Skip 2 bases, then
    // translate.
    //
    // Bases:   X X A T G C C C
    // Phase:   1 2 0 1 2 0 1 2     (0 = first base of codon)
    // Codons fully in this region: ATG, CCC → M P
    let (pool, region) = pool_from_string("XXATGCCC");
    let r = region.with_frame_phase(1);
    let rail = compute_codon_rail(&r, &pool);
    assert_eq!(rail.amino_acids, b"MP");
}

#[test]
fn compute_codon_rail_respects_frame_phase_two() {
    // frame_phase=2 means position 0 is the 3rd base of a codon.
    // Skip 1 base, then translate.
    //
    // Bases:   X T A C G G G
    // Phase:   2 0 1 2 0 1 2
    // Codons fully in this region: TAC, GGG → Y G
    let (pool, region) = pool_from_string("XTACGGG");
    let r = region.with_frame_phase(2);
    let rail = compute_codon_rail(&r, &pool);
    assert_eq!(rail.amino_acids, b"YG");
}

#[test]
fn compute_codon_rail_handles_ambiguous_bases() {
    // Codon containing N → X.
    let (pool, region) = pool_from_string("ATGNAATGG");
    let rail = compute_codon_rail(&region, &pool);
    // ATG = M, NAA = X (ambiguous), TGG = W
    assert_eq!(rail.amino_acids, b"MXW");
}

#[test]
fn compute_codon_rail_re_evaluates_after_base_mutation() {
    // Mutation in this region must produce a different rail when
    // re-computed against the mutated pool.
    let (pool0, region) = pool_from_string("ATGTACTGG");
    let rail0 = compute_codon_rail(&region, &pool0);
    assert_eq!(rail0.amino_acids, b"MYW");

    // Mutate the second base of the second codon: TAC → TGC = C.
    let pool1 = pool0.with_base_changed(NucHandle::new(4), b'G');
    let rail1 = compute_codon_rail(&region, &pool1);
    assert_eq!(rail1.amino_acids, b"MCW");

    // Old rail unaffected (the rail is a free-function return value;
    // no stale state to worry about).
    assert_eq!(rail0.amino_acids, b"MYW");
}

#[test]
fn compute_codon_rail_detects_stop_introduced_by_mutation() {
    // TAC → TAA (Y → stop) by mutating one base.
    let (pool0, region) = pool_from_string("ATGTACGGG");
    let rail0 = compute_codon_rail(&region, &pool0);
    assert_eq!(rail0.stop_codon_count(), 0);
    assert_eq!(rail0.amino_acids, b"MYG");

    // Mutate position 5: TAC → TAA. This creates a stop codon.
    let pool1 = pool0.with_base_changed(NucHandle::new(5), b'A');
    let rail1 = compute_codon_rail(&region, &pool1);
    assert_eq!(rail1.amino_acids, b"M*G");
    assert_eq!(rail1.stop_codon_count(), 1);
    assert_eq!(rail1.stop_codon_positions, vec![NucHandle::new(3)]);
}

#[test]
fn compute_codon_rail_empty_region_produces_empty_metadata() {
    let pool = NucleotidePool::new();
    let region = Region::new(Segment::Np1, NucHandle::new(0), NucHandle::new(0));
    let rail = compute_codon_rail(&region, &pool);
    assert!(rail.amino_acids.is_empty());
    assert!(rail.stop_codon_positions.is_empty());
}

#[test]
fn compute_codon_rail_frame_phase_2_on_one_base_region() {
    // frame_phase = 2 means the region's first base is the third
    // base of a codon started in a previous region. `skip = 1`,
    // so we'd need at least 4 bases (1 skipped + 3 for a fresh
    // codon) to emit anything. With one base, output is empty.
    let (pool, region) = pool_from_string("X");
    let r = region.with_frame_phase(2);
    let rail = compute_codon_rail(&r, &pool);
    assert!(rail.amino_acids.is_empty());
    assert!(rail.stop_codon_positions.is_empty());
}

#[test]
fn compute_codon_rail_frame_phase_1_on_two_base_region() {
    // frame_phase = 1 → skip 2 → no bases left after skip.
    // Output is empty.
    let (pool, region) = pool_from_string("XY");
    let r = region.with_frame_phase(1);
    let rail = compute_codon_rail(&r, &pool);
    assert!(rail.amino_acids.is_empty());
}

#[test]
fn compute_codon_rail_malformed_end_less_than_start_is_safe() {
    // Defensive contract: a region constructed with end < start
    // (which `len()` already saturates to 0) produces empty
    // codon rail metadata, not a panic.
    let mut pool = NucleotidePool::new();
    for i in 0..6 {
        pool.push(Nucleotide::germline(b'A', i, Segment::V));
    }
    let region = Region::new(Segment::V, NucHandle::new(5), NucHandle::new(2));
    assert_eq!(region.len(), 0);

    let rail = compute_codon_rail(&region, &pool);
    assert!(rail.amino_acids.is_empty());
    assert!(rail.stop_codon_positions.is_empty());
}

#[test]
fn compute_codon_rail_skip_overruns_end_is_safe() {
    // 1-base region with frame_phase=1 → skip=2 → start_idx > end_idx
    // after the skip. Same defensive case as above, different
    // path through the code.
    let (pool, region) = pool_from_string("X");
    let r = region.with_frame_phase(1);
    let rail = compute_codon_rail(&r, &pool);
    assert!(rail.amino_acids.is_empty());
    assert!(rail.stop_codon_positions.is_empty());
}
