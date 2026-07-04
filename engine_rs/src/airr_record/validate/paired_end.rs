//! C6: Paired-end / read-layout invariants.
//!
//! Slice A: enforces the no-layout default invariant (`read_layout
//! == ""` ⇒ every paired-end field is at its default).
//!
//! Slice B: extends the dispatch with per-layout geometry checks.
//! When `read_layout == "paired_end"`, the validator re-derives
//! R1/R2 windows from the rules in `docs/paired_end_design.md` §8
//! and surfaces any divergence as one of the four reserved variants
//! (`ReadWindowOutOfBounds`, `ReadSequenceMismatch`,
//! `ReadInsertSizeMismatch`). Unknown non-empty layouts surface as
//! `ReadLayoutMismatch`. `"single_end"` is documented as reserved
//! (§2.2) and treated as a no-op for now.

use super::*;

pub(super) fn check_paired_end_defaults(
    record: &AirrRecord,
    issues: &mut Vec<RecordValidationIssue>,
) {
    match record.read_layout.as_str() {
        // Slice A: no-layout default invariant.
        "" => check_paired_end_default_values(record, issues),
        // Slice B: geometry against the projection rules.
        "paired_end" => check_paired_end_geometry(record, issues),
        // Reserved (per §2.2). A future slice may attach a
        // single-read geometry check; today we don't validate
        // anything beyond accepting the layout string.
        "single_end" => {}
        // Anything else is an unsupported layout value. Surface
        // the structured mismatch so a typo (`"pair_end"`) or a
        // refactor that introduced a new layout without wiring
        // it through the validator fails closed.
        _ => issues.push(RecordValidationIssue::ReadLayoutMismatch {
            reported: record.read_layout.clone(),
            expected: r#"one of: "", "paired_end", "single_end""#.to_string(),
        }),
    }
}

fn check_paired_end_default_values(
    record: &AirrRecord,
    issues: &mut Vec<RecordValidationIssue>,
) {
    if !record.r1_sequence.is_empty() {
        issues.push(RecordValidationIssue::PairedEndFieldWithoutLayout {
            field: PairedEndField::R1Sequence,
        });
    }
    if !record.r2_sequence.is_empty() {
        issues.push(RecordValidationIssue::PairedEndFieldWithoutLayout {
            field: PairedEndField::R2Sequence,
        });
    }
    if record.r1_start.is_some() {
        issues.push(RecordValidationIssue::PairedEndFieldWithoutLayout {
            field: PairedEndField::R1Start,
        });
    }
    if record.r1_end.is_some() {
        issues.push(RecordValidationIssue::PairedEndFieldWithoutLayout {
            field: PairedEndField::R1End,
        });
    }
    if record.r2_start.is_some() {
        issues.push(RecordValidationIssue::PairedEndFieldWithoutLayout {
            field: PairedEndField::R2Start,
        });
    }
    if record.r2_end.is_some() {
        issues.push(RecordValidationIssue::PairedEndFieldWithoutLayout {
            field: PairedEndField::R2End,
        });
    }
    if record.insert_size != 0 {
        issues.push(RecordValidationIssue::PairedEndFieldWithoutLayout {
            field: PairedEndField::InsertSize,
        });
    }
}

fn check_paired_end_geometry(
    record: &AirrRecord,
    issues: &mut Vec<RecordValidationIssue>,
) {
    let seq_len = record.sequence_length;
    let seq = &record.sequence;

    // R1 window: bounds + byte equality.
    match resolve_window(record.r1_start, record.r1_end, seq_len) {
        WindowResolution::Valid { start, end } => {
            let expected = seq[start as usize..end as usize].to_string();
            if record.r1_sequence != expected {
                issues.push(RecordValidationIssue::ReadSequenceMismatch {
                    side: PairedEndRead::R1,
                    reported: record.r1_sequence.clone(),
                    expected,
                });
            }
        }
        WindowResolution::OutOfBounds { start, end } => {
            issues.push(RecordValidationIssue::ReadWindowOutOfBounds {
                side: PairedEndRead::R1,
                start,
                end,
                sequence_length: seq_len,
            });
        }
    }

    // R2 window: bounds + reverse-complement byte equality + insert
    // size consistency. The audit pins
    // `insert_size == r2_end` (§8) — the only insert-size-mismatch
    // surface today.
    match resolve_window(record.r2_start, record.r2_end, seq_len) {
        WindowResolution::Valid { start, end } => {
            let r2_inner = &seq[start as usize..end as usize];
            let expected = super::super::sequence::reverse_complement(r2_inner);
            if record.r2_sequence != expected {
                issues.push(RecordValidationIssue::ReadSequenceMismatch {
                    side: PairedEndRead::R2,
                    reported: record.r2_sequence.clone(),
                    expected,
                });
            }
            if record.insert_size != end {
                issues.push(RecordValidationIssue::ReadInsertSizeMismatch {
                    reported: record.insert_size,
                    expected: end,
                });
            }
        }
        WindowResolution::OutOfBounds { start, end } => {
            issues.push(RecordValidationIssue::ReadWindowOutOfBounds {
                side: PairedEndRead::R2,
                start,
                end,
                sequence_length: seq_len,
            });
            // With R2 bounds unresolved we can't recompute the
            // expected insert size; skip the insert-size check
            // rather than fabricate a sentinel. The window
            // out-of-bounds issue is the actionable signal.
        }
    }
}

/// Result of resolving a paired-end window `(start, end)` against
/// the projected `sequence_length`. `OutOfBounds` collapses missing
/// coords (`None`) and out-of-range coords (negative, swapped,
/// past the end) into one variant whose `start`/`end` fields use
/// the sentinel `-1` for missing values — surfaces cleanly in the
/// structured-issue payload without inventing a new variant.
enum WindowResolution {
    Valid { start: i64, end: i64 },
    OutOfBounds { start: i64, end: i64 },
}

fn resolve_window(
    start: Option<i64>,
    end: Option<i64>,
    sequence_length: i64,
) -> WindowResolution {
    let start_val = start.unwrap_or(-1);
    let end_val = end.unwrap_or(-1);
    match (start, end) {
        (Some(s), Some(e)) if s >= 0 && e >= s && e <= sequence_length => {
            WindowResolution::Valid { start: s, end: e }
        }
        _ => WindowResolution::OutOfBounds {
            start: start_val,
            end: end_val,
        },
    }
}
