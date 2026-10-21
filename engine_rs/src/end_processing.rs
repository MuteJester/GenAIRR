//! Version 1 finite outcomes at the original coding ends.
//!
//! Signed values are canonical-reference net deletions: positive removes
//! reference bases, negative adds a palindrome from the original tip.
//! This is an effective sequence model, not a reconstruction of hairpin
//! opening, strand resection, or polymerase history.

use crate::assignment::SegmentOrientation;
use crate::ir::complement_base;

/// One authored row of a joint two-ended outcome distribution.
///
/// Row order is replay-significant. Physical legality depends on the
/// selected allele and is checked by [`materialize`].
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct EndOutcome {
    pub five: i64,
    pub three: i64,
    pub weight: f64,
}

/// Deterministic output and template ancestry for one processed segment.
///
/// Reference bounds remain canonical even when the body is reversed.
/// P strings and their parallel template-position vectors are in assembled
/// orientation. P template positions are ancestry, not alignment positions.
#[derive(Clone, Debug, Eq, PartialEq)]
pub struct ProcessedSegment {
    pub retained_start: usize,
    pub retained_end: usize,
    pub body: Vec<u8>,
    pub p5: Vec<u8>,
    pub p3: Vec<u8>,
    pub p5_template_positions: Vec<u32>,
    pub p3_template_positions: Vec<u32>,
}

/// Materialize a joint original-end outcome without random draws or state.
///
/// Complete reference deletion is valid. Opposite positive deletions may
/// meet but cannot cross. Each extension is bounded by the original allele
/// length, independently of the opposite end's deletion.
pub fn materialize(
    seq: &[u8],
    five: i64,
    three: i64,
    orientation: SegmentOrientation,
) -> Result<ProcessedSegment, String> {
    let length = seq.len();
    if u32::try_from(length).is_err() {
        return Err("coding_end_reference_exceeds_coordinate_capacity".into());
    }
    let magnitude = |value: i64| {
        usize::try_from(value.unsigned_abs())
            .ok()
            .filter(|&amount| amount <= length)
            .ok_or_else(|| "coding_end_outcome_exceeds_original_allele".to_string())
    };
    let five_amount = magnitude(five)?;
    let three_amount = magnitude(three)?;
    let retained_start = if five > 0 { five_amount } else { 0 };
    let retained_end = length - if three > 0 { three_amount } else { 0 };
    if retained_start > retained_end {
        return Err("coding_end_deletions_overlap".into());
    }
    let p5_len = if five < 0 { five_amount } else { 0 };
    let p3_len = if three < 0 { three_amount } else { 0 };
    let p5_template_positions: Vec<u32> = (0..p5_len).rev().map(|i| i as u32).collect();
    let p3_template_positions: Vec<u32> =
        (length - p3_len..length).rev().map(|i| i as u32).collect();
    if p5_template_positions
        .iter()
        .chain(&p3_template_positions)
        .any(|&position| {
            !matches!(
                seq[position as usize].to_ascii_uppercase(),
                b'A' | b'C' | b'G' | b'T'
            )
        })
    {
        return Err("coding_end_ambiguous_original_template".into());
    }
    let copy = |positions: &[u32]| {
        positions
            .iter()
            .map(|&i| complement_base(seq[i as usize]))
            .collect::<Vec<_>>()
    };
    let mut processed = ProcessedSegment {
        retained_start,
        retained_end,
        body: seq[retained_start..retained_end].to_vec(),
        p5: copy(&p5_template_positions),
        p3: copy(&p3_template_positions),
        p5_template_positions,
        p3_template_positions,
    };
    if orientation.is_reverse() {
        let reverse_complement =
            |bytes: &[u8]| bytes.iter().rev().map(|&b| complement_base(b)).collect();
        processed.body = reverse_complement(&processed.body);
        let canonical_p5 = std::mem::take(&mut processed.p5);
        processed.p5 = reverse_complement(&processed.p3);
        processed.p3 = reverse_complement(&canonical_p5);
        std::mem::swap(
            &mut processed.p5_template_positions,
            &mut processed.p3_template_positions,
        );
        processed.p5_template_positions.reverse();
        processed.p3_template_positions.reverse();
    }
    Ok(processed)
}

#[cfg(test)]
mod tests {
    use super::*;
    use SegmentOrientation::{Forward, ReverseComplement};

    #[test]
    fn original_tip_templates_keep_canonical_ancestry() {
        let p = materialize(b"AACGTC", -2, -3, Forward).unwrap();
        assert_eq!((p.retained_start, p.retained_end), (0, 6));
        assert_eq!(p.body, b"AACGTC");
        assert_eq!(p.p5, b"TT");
        assert_eq!(p.p3, b"GAC");
        assert_eq!(p.p5_template_positions, [1, 0]);
        assert_eq!(p.p3_template_positions, [5, 4, 3]);
    }

    #[test]
    fn opposite_end_deletion_does_not_move_the_original_template() {
        let p = materialize(b"AACGTC", 4, -3, Forward).unwrap();
        assert_eq!((p.retained_start, p.retained_end), (4, 6));
        assert_eq!(p.body, b"TC");
        assert_eq!(p.p3, b"GAC");
        assert_eq!(p.p3_template_positions, [5, 4, 3]);
    }

    #[test]
    fn inversion_swaps_ends_and_reverses_template_positions() {
        let p = materialize(b"AACGTC", -2, -3, ReverseComplement).unwrap();
        assert_eq!(p.body, b"GACGTT");
        assert_eq!(p.p5, b"GTC");
        assert_eq!(p.p3, b"AA");
        assert_eq!(p.p5_template_positions, [3, 4, 5]);
        assert_eq!(p.p3_template_positions, [0, 1]);
        let p = materialize(b"AACGTC", 2, -2, ReverseComplement).unwrap();
        assert_eq!((p.retained_start, p.retained_end), (2, 6));
        assert_eq!(p.body, b"GACG");
        assert_eq!(p.p5, b"TC");
        assert_eq!(p.p5_template_positions, [4, 5]);
        assert!(p.p3.is_empty());
    }

    #[test]
    fn complete_deletion_is_distinct_from_crossing_deletions() {
        for orientation in [Forward, ReverseComplement] {
            let empty = materialize(b"ACG", 1, 2, orientation).unwrap();
            assert_eq!((empty.retained_start, empty.retained_end), (1, 1));
            assert!(empty.body.is_empty());
            assert!(materialize(b"ACG", 2, 2, orientation).is_err());
        }
    }

    #[test]
    fn signed_extremes_and_oversized_extensions_fail_without_panicking() {
        for (five, three) in [(i64::MIN, 0), (0, i64::MAX), (-4, 0), (0, -4)] {
            assert!(materialize(b"ACG", five, three, Forward).is_err());
        }
        assert!(materialize(b"", 0, 0, Forward).unwrap().body.is_empty());
        assert!(materialize(b"", -1, 0, Forward).is_err());
    }

    #[test]
    fn ambiguous_copied_tips_are_rejected_without_banning_unrelated_ambiguity() {
        for orientation in [Forward, ReverseComplement] {
            for reference in [b"RACGT".as_slice(), b"NACGT", b"UACGT", b"-ACGT"] {
                assert!(materialize(reference, -1, 0, orientation).is_err());
            }
            assert!(materialize(b"ACGTN", 0, -1, orientation).is_err());
            assert!(materialize(b"ACGTN", -2, 0, orientation).is_ok());
            assert!(materialize(b"ACNGT", -2, -2, orientation).is_ok());
            assert!(materialize(b"NNNNN", 1, 1, orientation).is_ok());
        }
        let lower = materialize(b"acgt", -2, -2, Forward).unwrap();
        assert_eq!(lower.p5, b"gt");
        assert_eq!(lower.p3, b"ac");
    }

    #[test]
    fn exhaustive_small_outcomes_reconstruct_the_oriented_whole_molecule() {
        let reference = b"AACGTC";
        for five in -6..=6i64 {
            for three in -6..=6i64 {
                let canonical = materialize(reference, five, three, Forward);
                assert_eq!(canonical.is_ok(), five.max(0) + three.max(0) <= 6);
                let Ok(canonical) = canonical else { continue };
                let reverse = materialize(reference, five, three, ReverseComplement).unwrap();
                let forward_bytes = [canonical.p5, canonical.body, canonical.p3].concat();
                let reverse_bytes = [reverse.p5.clone(), reverse.body, reverse.p3.clone()].concat();
                let expected: Vec<_> = forward_bytes
                    .iter()
                    .rev()
                    .map(|base| match base {
                        b'A' => b'T',
                        b'C' => b'G',
                        b'G' => b'C',
                        b'T' => b'A',
                        _ => unreachable!(),
                    })
                    .collect();
                assert_eq!(reverse_bytes, expected, "{five}, {three}");
                // In reverse orientation the two complement operations cancel:
                // each extension byte equals its canonical template base.
                for (bytes, positions) in [
                    (&reverse.p5, &reverse.p5_template_positions),
                    (&reverse.p3, &reverse.p3_template_positions),
                ] {
                    for (&base, &position) in bytes.iter().zip(positions) {
                        assert_eq!(base, reference[position as usize]);
                    }
                }
            }
        }
    }
}
