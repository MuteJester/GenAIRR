//! C1: Structural record invariants.

use super::*;

use super::super::sequence::pool_bases;
use crate::ir::Segment;

pub(super) fn check_structural(
    record: &AirrRecord,
    sim: &crate::ir::Simulation,
    issues: &mut Vec<RecordValidationIssue>,
) {
    // sequence_length matches actual sequence string bytes.
    let actual_bytes = record.sequence.len();
    if record.sequence_length as usize != actual_bytes {
        issues.push(RecordValidationIssue::SequenceLengthMismatch {
            reported: record.sequence_length,
            actual_bytes,
        });
    }

    // sequence content matches the pool's bases (case-folded).
    // Compare case-insensitively because sequencing errors lowercase
    // their substitutions, and the AIRR record may carry rev-comp.
    let pool_seq = String::from_utf8(pool_bases(sim)).unwrap_or_default();
    if !record.rev_comp && !record.sequence.eq_ignore_ascii_case(&pool_seq) {
        let prefix_len = record.sequence.len().min(pool_seq.len()).min(40);
        issues.push(RecordValidationIssue::SequenceContentMismatch {
            reported_prefix: record.sequence.chars().take(prefix_len).collect(),
            actual_prefix: pool_seq.chars().take(prefix_len).collect(),
        });
    }

    // Coordinate ordering for V/D/J segments and germline ranges.
    check_coord_pair(
        Segment::V,
        record.v_sequence_start,
        record.v_sequence_end,
        false,
        issues,
    );
    check_coord_pair(
        Segment::V,
        record.v_germline_start,
        record.v_germline_end,
        true,
        issues,
    );
    check_coord_pair(
        Segment::D,
        record.d_sequence_start,
        record.d_sequence_end,
        false,
        issues,
    );
    check_coord_pair(
        Segment::D,
        record.d_germline_start,
        record.d_germline_end,
        true,
        issues,
    );
    check_coord_pair(
        Segment::J,
        record.j_sequence_start,
        record.j_sequence_end,
        false,
        issues,
    );
    check_coord_pair(
        Segment::J,
        record.j_germline_start,
        record.j_germline_end,
        true,
        issues,
    );

    // CIGAR span sanity: M+I ops must equal the sequence-side span.
    for (seg, start, end, cigar) in [
        (
            Segment::V,
            record.v_sequence_start,
            record.v_sequence_end,
            &record.v_cigar,
        ),
        (
            Segment::D,
            record.d_sequence_start,
            record.d_sequence_end,
            &record.d_cigar,
        ),
        (
            Segment::J,
            record.j_sequence_start,
            record.j_sequence_end,
            &record.j_cigar,
        ),
    ] {
        if cigar.is_empty() {
            continue;
        }
        match parse_cigar(cigar) {
            Ok(ops) => {
                let query_span: usize = ops
                    .iter()
                    .filter(|(_, op)| *op == b'M' || *op == b'I')
                    .map(|(n, _)| *n)
                    .sum();
                if let (Some(s), Some(e)) = (start, end) {
                    let sequence_span = (e - s).max(0) as usize;
                    if query_span != sequence_span {
                        issues.push(RecordValidationIssue::CigarSpanMismatch {
                            segment: seg,
                            cigar_query_span: query_span,
                            sequence_span,
                        });
                    }
                }
            }
            Err(reason) => issues.push(RecordValidationIssue::CigarReadsInvalid {
                segment: seg,
                reason,
            }),
        }
    }
}

fn check_coord_pair(
    segment: Segment,
    start: Option<i64>,
    end: Option<i64>,
    is_germline: bool,
    issues: &mut Vec<RecordValidationIssue>,
) {
    if let (Some(s), Some(e)) = (start, end) {
        // Half-open: end > start is required (or both 0 for empty).
        if s < 0 || e < s {
            let issue = if is_germline {
                RecordValidationIssue::GermlineCoordinatesOutOfOrder {
                    segment,
                    start: s,
                    end: e,
                }
            } else {
                RecordValidationIssue::SegmentCoordinatesOutOfOrder {
                    segment,
                    start: s,
                    end: e,
                }
            };
            issues.push(issue);
        }
    }
}

fn parse_cigar(cigar: &str) -> Result<Vec<(usize, u8)>, String> {
    let mut ops = Vec::new();
    let mut digits = String::new();
    for ch in cigar.chars() {
        if ch.is_ascii_digit() {
            digits.push(ch);
        } else {
            let n: usize = digits
                .parse()
                .map_err(|_| format!("non-numeric op length in CIGAR {cigar:?}"))?;
            if !matches!(ch as u8, b'M' | b'I' | b'D' | b'S' | b'N' | b'P' | b'X' | b'=') {
                return Err(format!("unrecognized CIGAR op {ch:?} in {cigar:?}"));
            }
            ops.push((n, ch as u8));
            digits.clear();
        }
    }
    if !digits.is_empty() {
        return Err(format!("trailing digits in CIGAR {cigar:?}"));
    }
    Ok(ops)
}
