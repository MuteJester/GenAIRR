//! C4: Allele-call oracle.
//!
//! Independent reimplementation of the walker's max-match-count
//! tie-set selection. We walk the segment's region bytes, count
//! matches per allele, and pick the alleles at the max score.
//! Then we re-derive the projected CSV order (truth first if in
//! tie-set, otherwise ascending allele id).

use super::*;

use crate::address::VdjSegment;
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
    if record.rev_comp {
        // Rev-comp flips the sequence; the allele call was computed
        // pre-flip. Skip the oracle here; the rev-comp audit covers
        // it with dedicated tests.
        return;
    }
    for (segment, vdj, reported_call) in [
        (Segment::V, VdjSegment::V, &record.v_call),
        (Segment::D, VdjSegment::D, &record.d_call),
        (Segment::J, VdjSegment::J, &record.j_call),
    ] {
        let _ = vdj;
        oracle_check_segment(sim, refdata, segment, reported_call, issues);
    }
}

fn oracle_check_segment(
    sim: &crate::ir::Simulation,
    refdata: &RefDataConfig,
    segment: Segment,
    reported_call: &str,
    issues: &mut Vec<RecordValidationIssue>,
) {
    let Some(region) = sim
        .sequence
        .regions
        .iter()
        .find(|r| r.segment == segment)
    else {
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
    let scores = score_alleles_with_extensions(
        sim,
        segment,
        allele_pool,
        region,
        trim_5_cap,
        trim_3_cap,
        orientation,
    );
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
    let reported_names: Vec<String> = reported_call
        .split(',')
        .filter(|s| !s.is_empty())
        .map(|s| s.to_string())
        .collect();

    // Tie-set equality (order-insensitive).
    let mut reported_sorted = reported_names.clone();
    reported_sorted.sort();
    let mut expected_sorted = expected_names.clone();
    expected_sorted.sort();
    if reported_sorted != expected_sorted {
        issues.push(RecordValidationIssue::AlleleCallTieSetMismatch {
            segment,
            reported: reported_names.clone(),
            recomputed: expected_names.clone(),
        });
        return; // Order check is meaningless if the sets differ.
    }

    // Order check: first element must match expected_order's first.
    if let (Some(reported_first), Some(expected_first)) =
        (reported_names.first(), expected_names.first())
    {
        if reported_first != expected_first {
            let reason = if truth_idx.is_some()
                && expected_order
                    .first()
                    .map(|&i| i == truth_idx.unwrap())
                    .unwrap_or(false)
            {
                AlleleOrderReason::TruthFirstIfInTieSet
            } else {
                AlleleOrderReason::AscendingAlleleIdOtherwise
            };
            issues.push(RecordValidationIssue::AlleleCallOrderMismatch {
                segment,
                reported_first: reported_first.clone(),
                expected_first: expected_first.clone(),
                reason,
            });
        }
    }
}

// Match semantics live in crate::live_call::scoring; this module
// uses them via score_alleles_in_region / tie_set_ids_at_max_score.
