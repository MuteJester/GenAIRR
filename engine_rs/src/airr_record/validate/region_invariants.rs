//! C5: Region / live-call structural invariants.
//!
//! Per the §5 architectural audit:
//!   - Each Segment (V/D/J) appears at most once in
//!     sim.sequence.regions. Live-call and AIRR projection pick the
//!     "latest" and "first" respectively; identical-by-construction
//!     today, drift-able if invariant breaks.
//!   - SegmentLiveCall.hypotheses has length <= 1 in production runs.
//!     Projection silently uses hypotheses[0]; multi-hypothesis would
//!     be lossy.

use super::*;

use crate::ir::Segment;

pub(super) fn check_region_and_hypothesis_invariants(
    sim: &crate::ir::Simulation,
    issues: &mut Vec<RecordValidationIssue>,
) {
    for seg in [Segment::V, Segment::D, Segment::J] {
        let count = sim
            .sequence
            .regions
            .iter()
            .filter(|r| r.segment == seg)
            .count();
        if count > 1 {
            issues.push(RecordValidationIssue::MultipleRegionsForSegment {
                segment: seg,
                count,
            });
        }
    }

    for seg in [Segment::V, Segment::D, Segment::J] {
        if let Some(call) = sim.segment_calls.get(seg) {
            if call.hypotheses.len() > 1 {
                issues.push(RecordValidationIssue::MultipleHypothesesInLiveCall {
                    segment: seg,
                    count: call.hypotheses.len(),
                });
            }
        }
    }
}
