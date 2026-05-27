//! Tests for the shared allele-call scoring kernel.
//!
//! Covers the documented match semantics + tie-set construction
//! plus a walker-vs-validator agreement check on a synthetic fixture.
//! The walker uses the inverted-index path (`ReferenceMatchIndex`);
//! the validator uses `score_alleles_in_region`. Both route their
//! match decisions through `scoring::classify_base`, so they must
//! return the same tie-set on any input.

use crate::ir::{flag, NucHandle, Nucleotide, Region, Segment, Simulation};
use crate::live_call::scoring::{
    classify_base, matches_observed, score_alleles_in_region, tie_set_at_max_score,
    tie_set_ids_at_max_score, BaseKind,
};
use crate::live_call::AlleleBitSet;
use crate::refdata::{Allele, AlleleId, AllelePool};

// ──────────────────────────────────────────────────────────────────
// Match semantics — classify_base / matches_observed
// ──────────────────────────────────────────────────────────────────

#[test]
fn canonical_bases_classify_to_uppercase_canonical() {
    for b in *b"ACGT" {
        assert_eq!(classify_base(b), BaseKind::Canonical(b));
    }
}

#[test]
fn lowercase_bases_classify_to_uppercase_canonical_not_wildcard() {
    // The sequencing-error convention writes mutated bases as
    // lowercase. The scoring kernel must treat them as their
    // uppercase counterparts, NOT as wildcards. This is the load-
    // bearing rule from the allele-call audit §1.3.
    for (lower, upper) in [(b'a', b'A'), (b'c', b'C'), (b'g', b'G'), (b't', b'T')] {
        assert_eq!(classify_base(lower), BaseKind::Canonical(upper));
    }
}

#[test]
fn n_and_lowercase_n_classify_to_wildcard() {
    assert_eq!(classify_base(b'N'), BaseKind::Wildcard);
    assert_eq!(classify_base(b'n'), BaseKind::Wildcard);
}

#[test]
fn non_iupac_bases_classify_to_invalid() {
    for b in [b'.', b'-', b'X', b'?', 0u8, b' '] {
        assert_eq!(classify_base(b), BaseKind::Invalid);
    }
}

#[test]
fn matches_observed_canonical_equal_matches() {
    assert!(matches_observed(b'A', b'A'));
    assert!(matches_observed(b'c', b'C')); // lowercase as canonical
    assert!(matches_observed(b'G', b'g'));
}

#[test]
fn matches_observed_canonical_unequal_does_not_match() {
    assert!(!matches_observed(b'A', b'C'));
    assert!(!matches_observed(b't', b'A')); // lowercase t != A
}

#[test]
fn matches_observed_n_is_wildcard_against_canonical() {
    for g in *b"ACGT" {
        assert!(matches_observed(b'N', g), "N should wildcard-match {}", g as char);
        assert!(matches_observed(b'n', g));
    }
}

#[test]
fn matches_observed_invalid_never_matches() {
    assert!(!matches_observed(b'.', b'A'));
    assert!(!matches_observed(b'A', b'.'));
    assert!(!matches_observed(0, b'A'));
}

// ──────────────────────────────────────────────────────────────────
// Region scorer — score_alleles_in_region
// ──────────────────────────────────────────────────────────────────

fn allele(name: &str, seq: &[u8]) -> Allele {
    Allele {
        name: name.to_string(),
        gene: name.split('*').next().unwrap_or(name).to_string(),
        seq: seq.to_vec(),
        segment: Segment::V,
        anchor: None,
    }
}

/// Build a Simulation whose pool contains exactly `seq`, with each
/// byte tagged `germline_pos == its pool index` and `segment = V`.
/// Returns the sim plus the V region [0, seq.len()).
fn sim_with_v_region(seq: &[u8]) -> (Simulation, Region) {
    let mut sim = Simulation::new();
    for (i, &b) in seq.iter().enumerate() {
        let nuc = Nucleotide::germline(b, i as u16, Segment::V);
        let (next, _h) = sim.with_nucleotide_pushed(nuc);
        sim = next;
    }
    let region = Region::new(
        Segment::V,
        NucHandle::new(0),
        NucHandle::new(seq.len() as u32),
    );
    sim = sim.with_region_added(region.clone());
    (sim, region)
}

#[test]
fn identical_alleles_produce_a_tie() {
    // Two byte-identical alleles must always have equal scores
    // against any observed region — neither can ever lose.
    let alleles = [allele("v1*01", b"AAACCCGGG"), allele("v1*02", b"AAACCCGGG")];
    let (sim, region) = sim_with_v_region(b"AAACCCGGG");
    let scores = score_alleles_in_region(
        sim.pool.as_slice(),
        region.start.index(),
        region.end.index(),
        &alleles,
    );
    assert_eq!(scores[0], scores[1]);
    assert_eq!(scores[0], 9);
    let tied = tie_set_ids_at_max_score(&scores);
    assert_eq!(tied, vec![AlleleId::new(0), AlleleId::new(1)]);
}

#[test]
fn distinguishing_snp_picks_one_allele() {
    // v1 and v2 differ at pos 4 (C vs T). Observed sequence has C
    // at pos 4 → v1 scores 9, v2 scores 8.
    let alleles = [allele("v1*01", b"AAACCCGGG"), allele("v2*01", b"AAACTCGGG")];
    let (sim, region) = sim_with_v_region(b"AAACCCGGG");
    let scores = score_alleles_in_region(
        sim.pool.as_slice(),
        region.start.index(),
        region.end.index(),
        &alleles,
    );
    assert_eq!(scores, vec![9, 8]);
    let tied = tie_set_ids_at_max_score(&scores);
    assert_eq!(tied, vec![AlleleId::new(0)]);
}

#[test]
fn n_at_distinguishing_position_widens_tie_set() {
    // Same v1/v2 fixture. Observed has N at pos 4 (the
    // distinguishing position). Wildcard matches both alleles
    // there → both score 9 → tie.
    let alleles = [allele("v1*01", b"AAACCCGGG"), allele("v2*01", b"AAACTCGGG")];
    let (sim, region) = sim_with_v_region(b"AAACNCGGG");
    let scores = score_alleles_in_region(
        sim.pool.as_slice(),
        region.start.index(),
        region.end.index(),
        &alleles,
    );
    assert_eq!(scores, vec![9, 9]);
}

#[test]
fn lowercase_mutation_at_distinguishing_position_does_not_widen() {
    // Observed has lowercase 't' at pos 4 (the distinguishing
    // position). 't' is upper-cased to T, so it matches v2 (T at
    // pos 4) and mismatches v1 (C at pos 4). The call switches,
    // not widens.
    let alleles = [allele("v1*01", b"AAACCCGGG"), allele("v2*01", b"AAACTCGGG")];
    let (sim, region) = sim_with_v_region(b"AAACtCGGG");
    let scores = score_alleles_in_region(
        sim.pool.as_slice(),
        region.start.index(),
        region.end.index(),
        &alleles,
    );
    assert_eq!(scores, vec![8, 9]);
    let tied = tie_set_ids_at_max_score(&scores);
    assert_eq!(tied, vec![AlleleId::new(1)]); // v2 wins
}

#[test]
fn mismatching_base_does_not_create_false_positive_evidence() {
    // Observed has 'G' at pos 4 (matches neither v1's C nor v2's
    // T). Both alleles score 8. Tie at lower max — but the lower
    // max IS the new max, so the tie-set widens to include both.
    // This pins the rule: a mismatched base reduces every allele's
    // score equally; it does not add evidence to anyone.
    let alleles = [allele("v1*01", b"AAACCCGGG"), allele("v2*01", b"AAACTCGGG")];
    let (sim, region) = sim_with_v_region(b"AAACGCGGG");
    let scores = score_alleles_in_region(
        sim.pool.as_slice(),
        region.start.index(),
        region.end.index(),
        &alleles,
    );
    assert_eq!(scores, vec![8, 8]);
}

#[test]
fn indel_insert_bytes_contribute_no_evidence() {
    // Build a sim with an insertion at position 2: pool[2] has
    // germline_pos == None (synthetic indel insert). The scoring
    // kernel must skip it.
    let alleles = [allele("v1*01", b"AAACCCGGG")];
    let mut sim = Simulation::new();
    let germline_seq = b"AAACCCGGG";
    // First 2 germline bytes
    for (i, &b) in germline_seq[..2].iter().enumerate() {
        let (next, _) = sim.with_nucleotide_pushed(Nucleotide::germline(b, i as u16, Segment::V));
        sim = next;
    }
    // One synthetic insertion (germline_pos = None)
    let (next, _) = sim.with_nucleotide_pushed(Nucleotide::synthetic(
        b'T',
        Segment::V,
        flag::INDEL_INSERTED,
    ));
    sim = next;
    // Remaining germline bytes
    for (i, &b) in germline_seq[2..].iter().enumerate() {
        let gp = (i + 2) as u16;
        let (next, _) = sim.with_nucleotide_pushed(Nucleotide::germline(b, gp, Segment::V));
        sim = next;
    }
    let region = Region::new(Segment::V, NucHandle::new(0), NucHandle::new(10));
    sim = sim.with_region_added(region.clone());

    let scores = score_alleles_in_region(
        sim.pool.as_slice(),
        region.start.index(),
        region.end.index(),
        &alleles,
    );
    // 9 germline bytes that all match → score 9. The indel-insert
    // byte at pool[2] is skipped despite being a T that wouldn't
    // match the allele's reference at pos 2.
    assert_eq!(scores, vec![9]);
}

// ──────────────────────────────────────────────────────────────────
// Tie-set construction
// ──────────────────────────────────────────────────────────────────

#[test]
fn tie_set_at_max_score_returns_none_for_zero_evidence() {
    let scores = vec![0u32, 0, 0];
    assert!(tie_set_at_max_score(&scores, 3).is_none());
    assert!(tie_set_ids_at_max_score(&scores).is_empty());
}

#[test]
fn tie_set_at_max_score_picks_all_alleles_at_max() {
    let scores = vec![5u32, 7, 7, 6];
    let set = tie_set_at_max_score(&scores, 4).expect("non-zero scores");
    let mut expected = AlleleBitSet::empty(4);
    expected.insert(AlleleId::new(1));
    expected.insert(AlleleId::new(2));
    assert_eq!(set, expected);
}

#[test]
fn tie_set_ids_at_max_score_returns_sorted_by_index() {
    let scores = vec![3u32, 5, 5, 4, 5];
    let ids = tie_set_ids_at_max_score(&scores);
    assert_eq!(
        ids,
        vec![AlleleId::new(1), AlleleId::new(2), AlleleId::new(4)]
    );
}

// ──────────────────────────────────────────────────────────────────
// Walker-validator agreement
// ──────────────────────────────────────────────────────────────────

#[test]
fn walker_and_kernel_agree_on_synthetic_fixture() {
    // Build a region whose observed pool exactly matches v1's
    // germline. The walker (via ReferenceMatchIndex /
    // walker_observer) and the kernel's score_alleles_in_region
    // must agree on the per-allele match counts: same alleles,
    // same scores, same max-score tie-set.
    let alleles = vec![
        allele("v1*01", b"AAACCCGGG"),
        allele("v2*01", b"AAACTCGGG"),
        allele("v3*01", b"AAATCCGGG"),
    ];
    let (sim, region) = sim_with_v_region(b"AAACCCGGG");

    // Kernel path.
    let kernel_scores = score_alleles_in_region(
        sim.pool.as_slice(),
        region.start.index(),
        region.end.index(),
        &alleles,
    );

    // Walker path: build the same per-allele count by routing
    // each pool byte through the walker's match check (via
    // ReferenceMatchIndex), incrementing scores the same way the
    // walker_observer does.
    let allele_pool = AllelePool::from_vec(alleles.clone());
    let segment_index = crate::live_call::reference_index::SegmentRefIndex::build(
        Segment::V,
        &allele_pool,
        crate::live_call::DEFAULT_REFERENCE_KMER_LEN,
    );
    let mut walker_scores = vec![0u32; alleles.len()];
    for i in region.start.index()..region.end.index() {
        let nuc = &sim.pool.as_slice()[i as usize];
        let Some(gp) = nuc.germline_pos.get() else {
            continue;
        };
        if let Some(evidence) = segment_index.compatible_alleles_at(gp as usize, nuc.base) {
            evidence.allele_ids.for_each_id(|id| {
                walker_scores[id.as_usize()] += 1;
            });
        }
    }

    assert_eq!(
        kernel_scores, walker_scores,
        "kernel and walker disagree on synthetic fixture"
    );

    // Tie-set agreement too.
    let kernel_tie = tie_set_ids_at_max_score(&kernel_scores);
    let walker_tie = tie_set_ids_at_max_score(&walker_scores);
    assert_eq!(kernel_tie, walker_tie);
}

#[test]
fn extension_narrows_tie_set_basics() {
    use crate::live_call::scoring::extension_narrows_tie_set;
    // No evidence yet — extensions can't narrow from zero.
    assert!(!extension_narrows_tie_set(&[0, 0, 0], &[true, true, false]));

    // Two alleles tied at 5; only one matches at the next byte.
    // Strictly narrows: the matched one moves to 6, the missed one
    // stays at 5 → new tie set is the matched-only allele.
    assert!(extension_narrows_tie_set(&[5, 5, 3], &[true, false, false]));

    // All current max-tied alleles match — extension would just
    // shift the tie set up by 1, not narrow it.
    assert!(!extension_narrows_tie_set(&[5, 5, 3], &[true, true, false]));

    // No current max-tied allele matches — extension would widen
    // (pre_max+1 unreachable from current ties).
    assert!(!extension_narrows_tie_set(&[5, 5, 3], &[false, false, true]));
}

#[test]
fn alleles_compatible_at_matches_per_allele_byte() {
    use crate::live_call::scoring::alleles_compatible_at;
    let alleles = vec![
        allele("v1*01", b"AAACCCGGG"),
        allele("v2*01", b"AAACTCGGG"),
    ];
    // ref_pos 4: v1='C', v2='T'. Observed 'C' matches v1 only.
    assert_eq!(alleles_compatible_at(&alleles, 4, b'C'), vec![true, false]);
    // Observed 'T' matches v2 only.
    assert_eq!(alleles_compatible_at(&alleles, 4, b'T'), vec![false, true]);
    // Observed 'N' wildcards both.
    assert_eq!(alleles_compatible_at(&alleles, 4, b'N'), vec![true, true]);
    // ref_pos past allele length: no match.
    assert_eq!(alleles_compatible_at(&alleles, 100, b'A'), vec![false, false]);
}

#[test]
fn extension_aware_oracle_narrows_tie_under_right_extension() {
    //! V alleles tie on the structural region but differ at a
    //! position that's trimmed off the 3' end. NP1 follows V with
    //! a byte that matches one allele's trimmed-off position.
    //! Under trim_3=1 the right extension recovers that byte's
    //! evidence and narrows the tie.
    use crate::ir::Region;
    use crate::live_call::scoring::score_alleles_with_extensions;

    // v1 has 'A' at pos 8, v2 has 'C' at pos 8. Both have 'T' at
    // pos 9. Identical at positions 0-7.
    let alleles = vec![
        allele("v1*01", b"AAACCCGGAT"),
        allele("v2*01", b"AAACCCGGCT"),
    ];

    // V structural region: 8 bytes at germline_pos 0..8 (trim_3=2
    // → positions 8 and 9 trimmed). NP1 follows immediately, with
    // its first byte = 'A' (which matches v1[8]='A' but not
    // v2[8]='C').
    let mut sim = Simulation::new();
    for (i, &b) in b"AAACCCGG".iter().enumerate() {
        let (next, _) = sim.with_nucleotide_pushed(Nucleotide::germline(
            b,
            i as u16,
            Segment::V,
        ));
        sim = next;
    }
    let (next, _) = sim.with_nucleotide_pushed(Nucleotide::germline(
        b'A', 0, Segment::Np1,
    ));
    sim = next;
    let v_region = Region::new(Segment::V, NucHandle::new(0), NucHandle::new(8));
    let np_region = Region::new(Segment::Np1, NucHandle::new(8), NucHandle::new(9));
    sim = sim.with_region_added(v_region.clone());
    sim = sim.with_region_added(np_region);

    // Without extension cap: both alleles tie at 8 (identical on
    // structural positions 0..7).
    let structural_only =
        score_alleles_with_extensions(&sim, Segment::V, &alleles, &v_region, 0, 0);
    assert_eq!(structural_only, vec![8, 8]);

    // With trim_3=1: right extension steps 1 byte into NP1 at
    // ref_pos=8. v1[8]='A' matches the observed 'A'; v2[8]='C'
    // does not. Extension narrows the tie set → v1 to 9, v2 stays
    // at 8.
    let with_extension =
        score_alleles_with_extensions(&sim, Segment::V, &alleles, &v_region, 0, 1);
    assert_eq!(with_extension, vec![9, 8]);
    let tie = tie_set_ids_at_max_score(&with_extension);
    assert_eq!(tie, vec![AlleleId::new(0)]);
}

#[test]
fn extension_aware_oracle_skips_extension_when_no_narrowing() {
    //! When the NP byte matches BOTH alleles, the extension
    //! doesn't narrow the tie set and the walker (per its
    //! conservative policy) skips it. The oracle must do the same.
    use crate::ir::Region;
    use crate::live_call::scoring::score_alleles_with_extensions;

    // v1[8] = v2[8] = 'A'. NP1 byte = 'A' matches both → no
    // narrowing.
    let alleles = vec![
        allele("v1*01", b"AAACCCGGAT"),
        allele("v2*01", b"GAACCCGGAT"),  // differs at pos 0 only
    ];

    let mut sim = Simulation::new();
    for (i, &b) in b"AAACCCGG".iter().enumerate() {
        let (next, _) = sim.with_nucleotide_pushed(Nucleotide::germline(
            b,
            i as u16,
            Segment::V,
        ));
        sim = next;
    }
    let (next, _) = sim.with_nucleotide_pushed(Nucleotide::germline(
        b'A', 0, Segment::Np1,
    ));
    sim = next;
    let v_region = Region::new(Segment::V, NucHandle::new(0), NucHandle::new(8));
    let np_region = Region::new(Segment::Np1, NucHandle::new(8), NucHandle::new(9));
    sim = sim.with_region_added(v_region.clone());
    sim = sim.with_region_added(np_region);

    // v2[0]='G' differs from observed 'A' at structural pos 0 →
    // v1 scores 8, v2 scores 7. Already not tied. Extension at
    // ref_pos 8 would match v1 (A=A); v2 also matches (A=A) but
    // v2 is not at max → not a narrowing concern.
    //
    // Pre-max = 8 (v1 only). matched at pos 8 = [v1, v2] (both A).
    // Per the narrowing rule, the byte must:
    //   - have at least one max-tied allele matching (v1: yes)
    //   - have at least one max-tied allele NOT matching (no; v1
    //     is the only max-tied and it matches)
    // → no narrowing → extension skipped.
    let scores =
        score_alleles_with_extensions(&sim, Segment::V, &alleles, &v_region, 0, 1);
    assert_eq!(
        scores,
        vec![8, 7],
        "extension should be skipped because it doesn't narrow"
    );
}
