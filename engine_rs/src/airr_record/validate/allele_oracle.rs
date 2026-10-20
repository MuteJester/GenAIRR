//! C4: Allele-call oracle.
//!
//! Standard calls must name the final generating assignment. Evidence-based
//! recalled calls are independently re-scored from the observed segment bytes,
//! with truth-first ordering only when the source remains in the tie set.

use super::*;

use crate::ir::Segment;
use crate::live_call::scoring::{
    allele_pool_for_segment, score_alleles_with_extensions, tie_set_ids_at_max_score,
};

pub(super) fn check_allele_oracle(
    record: &AirrRecord,
    outcome: &Outcome,
    refdata: &RefDataConfig,
    issues: &mut Vec<RecordValidationIssue>,
) {
    let sim = outcome.final_simulation();
    for (segment, reported_call, recalled_call) in [
        (Segment::V, &record.v_call, &record.recalled_v_call),
        (Segment::D, &record.d_call, &record.recalled_d_call),
        (Segment::J, &record.j_call, &record.recalled_j_call),
    ] {
        let expected = sim
            .assignments
            .get(segment)
            .and_then(|a| allele_pool_for_segment(refdata, segment)?.get(a.allele_id.as_usize()))
            .map(|a| a.name.clone())
            .unwrap_or_default();
        if *reported_call != expected {
            issues.push(RecordValidationIssue::AlleleCallTieSetMismatch {
                segment,
                reported: reported_call
                    .split(',')
                    .filter(|s| !s.is_empty())
                    .map(str::to_owned)
                    .collect(),
                recomputed: if expected.is_empty() {
                    vec![]
                } else {
                    vec![expected]
                },
            });
        }
        oracle_check_segment(sim, refdata, segment, recalled_call, issues);
    }
}

fn oracle_check_segment(
    sim: &crate::ir::Simulation,
    refdata: &RefDataConfig,
    segment: Segment,
    reported_call: &str,
    issues: &mut Vec<RecordValidationIssue>,
) {
    // Direct/native partial outcomes may omit inference entirely. Missing
    // inference must export an empty field, rather than a generating fallback.
    // Cache parity is checked separately against a fresh walk when requested.
    if sim
        .segment_calls
        .get(segment)
        .map_or(true, |call| call.allele_call.is_empty())
    {
        if !reported_call.is_empty() {
            issues.push(RecordValidationIssue::RecalledAlleleCallMismatch {
                segment,
                reported: reported_call.to_owned(),
                recomputed: String::new(),
            });
        }
        return;
    }
    let Some(region) = sim.sequence.regions.iter().find(|r| r.segment == segment) else {
        return;
    };
    let Some(allele_pool) = allele_pool_for_segment(refdata, segment) else {
        return;
    };
    if allele_pool.is_empty() {
        return;
    }
    let assignment = sim.assignments.get(segment);
    let truth_id = assignment.map(|a| a.allele_id);
    let trim_5_cap = assignment.map(|a| a.trim_5 as u32).unwrap_or(0);
    let trim_3_cap = assignment.map(|a| a.trim_3 as u32).unwrap_or(0);
    // Orientation drives the per-byte comparison rule via the
    // shared `matches_observed_with_orientation` primitive: under
    // `ReverseComplement` the observed byte is pre-complemented
    // before matching the allele's germline byte at the same
    // `germline_pos`. Defaults to `Forward` when the segment is
    // unassigned. See `scoring::observed_in_germline_orientation`
    // for the rationale.
    let orientation = assignment
        .map(|a| a.orientation)
        .unwrap_or(crate::assignment::SegmentOrientation::Forward);

    // Independent rescore via the shared scoring kernel: structural
    // region + NP-region extensions under the assigned allele's trim
    // caps. Mirrors `live_call::walker::call_from_region` so the
    // oracle and the walker agree on the tie-set under arbitrary
    // trim.
    let mut scores = score_alleles_with_extensions(
        sim,
        segment,
        allele_pool,
        region,
        trim_5_cap,
        trim_3_cap,
        orientation,
    );
    if let Some(allowed) = sim.live_call_alleles(segment) {
        for (id, score) in scores.iter_mut().enumerate() {
            if !allowed.iter().any(|a| a.as_usize() == id) {
                *score = 0;
            }
        }
    }
    let tied_ids = tie_set_ids_at_max_score(&scores);
    if tied_ids.is_empty() {
        return; // No germline evidence; oracle abstains.
    }
    let tied_indices: Vec<usize> = tied_ids.iter().map(|id| id.as_usize()).collect();

    // Expected CSV order: truth allele first when in tie-set,
    // otherwise ascending by allele id (already sorted by the
    // kernel's iteration order).
    let truth_idx = truth_id
        .map(|id| id.as_usize())
        .filter(|&i| i < allele_pool.len());

    let mut expected_order = tied_indices.clone();
    if let Some(t) = truth_idx {
        if let Some(pos) = expected_order.iter().position(|&i| i == t) {
            let truth_first = expected_order.remove(pos);
            expected_order.insert(0, truth_first);
        }
    }

    let expected_names: Vec<String> = expected_order
        .iter()
        .map(|&i| allele_pool[i].name.clone())
        .collect();
    let recomputed = expected_names.join(",");
    if reported_call != recomputed {
        issues.push(RecordValidationIssue::RecalledAlleleCallMismatch {
            segment,
            reported: reported_call.to_string(),
            recomputed,
        });
    }
}

// Match semantics live in crate::live_call::scoring; this module
// uses them via score_alleles_in_region / tie_set_ids_at_max_score.
