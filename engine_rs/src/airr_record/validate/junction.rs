//! Junction and nullable productivity: independent evidence reconstruction.
use super::*;

pub(super) fn check_junction(
    record: &AirrRecord,
    sim: &crate::ir::Simulation,
    refdata: &RefDataConfig,
    issues: &mut Vec<RecordValidationIssue>,
) {
    let expected = super::productivity_oracle::derive(sim, refdata);
    for (field, reported, recomputed) in [
        ("productive", record.productive, expected.productive),
        ("vj_in_frame", record.vj_in_frame, expected.frame),
        ("stop_codon", record.stop_codon, expected.stop),
    ] {
        if reported == recomputed {
            continue;
        }
        match (field, recomputed) {
            ("productive", Some(recomputed)) => {
                issues.push(RecordValidationIssue::ProductiveMismatch {
                    reported,
                    recomputed,
                    reason: if expected.frame == Some(false) {
                        ProductiveDecidedBy::OutOfFrame
                    } else if expected.stop == Some(true) {
                        ProductiveDecidedBy::CodingStopCodon
                    } else if recomputed {
                        ProductiveDecidedBy::InFrameAndAnchorsPreserved
                    } else {
                        ProductiveDecidedBy::VAnchorAaChanged
                    },
                })
            }
            ("vj_in_frame", Some(recomputed)) => {
                issues.push(RecordValidationIssue::VjInFrameMismatch {
                    reported,
                    recomputed,
                })
            }
            ("stop_codon", Some(recomputed)) => {
                issues.push(RecordValidationIssue::StopCodonMismatch {
                    reported,
                    recomputed,
                })
            }
            _ => issues.push(RecordValidationIssue::ProductivityAssessmentMismatch {
                field: field.into(),
                reported: format!("{reported:?}"),
                expected: format!("{recomputed:?}"),
            }),
        }
    }
    for (field, reported, recomputed) in [
        (
            "productivity_policy",
            record.productivity_policy.clone(),
            crate::productivity::POLICY.to_string(),
        ),
        (
            "productivity_scope",
            record.productivity_scope.clone(),
            crate::productivity::SCOPE.to_string(),
        ),
        (
            "productivity_reasons",
            record.productivity_reasons.join(","),
            expected.reasons.join(","),
        ),
    ] {
        if reported != recomputed {
            issues.push(RecordValidationIssue::ProductivityAssessmentMismatch {
                field: field.into(),
                reported,
                expected: recomputed,
            });
        }
    }
    let sequence: String = sim.pool.as_slice().iter().map(|n| n.base as char).collect();
    let content = expected
        .junction
        .map(|(s, e)| sequence[s..e].to_string())
        .unwrap_or_default();
    let aa = crate::codon::translate_seq(&content);
    if !record.junction.eq_ignore_ascii_case(&content) {
        issues.push(RecordValidationIssue::JunctionContentMismatch {
            reported: record.junction.clone(),
            recomputed: content,
        });
    }
    if record.junction_aa != aa {
        issues.push(RecordValidationIssue::JunctionAaMismatch {
            reported: record.junction_aa.clone(),
            recomputed: aa,
        });
    }
    for (field, reported, recomputed) in [
        (
            "junction_start",
            record.junction_start,
            expected.junction.map(|(s, _)| s as i64),
        ),
        (
            "junction_end",
            record.junction_end,
            expected.junction.map(|(_, e)| e as i64),
        ),
        (
            "junction_length",
            record.junction_length,
            expected.junction.map(|(s, e)| (e - s) as i64),
        ),
    ] {
        if reported != recomputed {
            if field == "junction_length" && recomputed.is_some() {
                issues.push(RecordValidationIssue::JunctionLengthMismatch {
                    reported,
                    recomputed: recomputed.unwrap() as u32,
                });
            } else {
                issues.push(RecordValidationIssue::ProductivityAssessmentMismatch {
                    field: field.into(),
                    reported: format!("{reported:?}"),
                    expected: format!("{recomputed:?}"),
                });
            }
        }
    }
}
