use super::super::{HypothesisFlags, SegmentRefIndex};
use crate::ir::{NucHandle, Region, Segment, Simulation};
use crate::refdata::AlleleId;

/// Conservative narrowing check for the boundary extension walkers.
/// The walkers extend a V/D/J live-call hypothesis one byte at a
/// time into an adjacent NP region (or across a structural boundary
/// into a neighbor segment). The current policy is *conservative*:
/// extend only when the byte strictly narrows the current max-score
/// tie set (as opposed to the greedy "extend whenever the byte
/// matches some allele" alternative).
///
/// "Narrowing" is precisely:
/// - among alleles currently tied at `pre_max`,
/// - at least one is in `matched_ids` (so post_max = pre_max + 1
///   exists),
/// - and at least one is NOT in `matched_ids` (so the post tie set
///   excludes them, strictly shrinking).
///
/// Cases ruled out:
/// - All current max-score alleles match → tie set size unchanged
///   after extension (everyone moves up by one). Don't extend.
/// - No current max-score allele matches → would *widen* the tie set
///   to include former-non-max alleles that bumped up to pre_max+
///   (some byte they happened to match before). Don't extend.
/// - Score floor (pre_max == 0) → nothing to narrow against (no
///   informative bytes have landed yet). Don't extend; the structural
///   walk is what populates initial scores.
///
/// Returns `true` iff the extension is allowed to proceed.
fn extension_narrows_tie_set(scores: &[u32], matched_ids: &super::super::AlleleBitSet) -> bool {
    let pre_max = scores.iter().copied().max().unwrap_or(0);
    if pre_max == 0 {
        return false;
    }
    let mut any_matched = false;
    let mut any_missed = false;
    for (i, &score) in scores.iter().enumerate() {
        if score != pre_max {
            continue;
        }
        let id = AlleleId::new(i as u32);
        if matched_ids.contains(id) {
            any_matched = true;
        } else {
            any_missed = true;
        }
        if any_matched && any_missed {
            return true;
        }
    }
    false
}

pub(crate) struct ExtensionWalkState<'a> {
    pub scores: &'a mut [u32],
    pub informative_matches: &'a mut u32,
    pub wildcard_matches: &'a mut u32,
    pub ref_start: &'a mut u32,
    pub ref_end: &'a mut u32,
    pub seq_start: &'a mut u32,
    pub seq_end: &'a mut u32,
    pub flags: &'a mut HypothesisFlags,
}

pub(crate) fn walk_left_extension(
    sim: &Simulation,
    segment_index: &SegmentRefIndex,
    np_region: &Region,
    trim_cap_5: Option<u32>,
    state: &mut ExtensionWalkState<'_>,
) {
    if np_region.end.index() != *state.seq_start {
        return;
    }

    let np_start = np_region.start.index();
    let mut extended_seq_start = *state.seq_start;
    let mut extended_ref_start = *state.ref_start;
    let mut extended_inform = *state.informative_matches;
    let mut extended_wildcard = *state.wildcard_matches;
    let mut extension_added = false;
    let mut crossed_into_overlap = false;
    let mut steps_taken: u32 = 0;

    let lower_bound: u32 = sim
        .sequence
        .regions
        .iter()
        .find(|r| {
            matches!(r.segment, Segment::V | Segment::D | Segment::J) && r.end == np_region.start
        })
        .map(|r| r.start.index())
        .unwrap_or(0);

    for seq_pos in (lower_bound..np_region.end.index()).rev() {
        if let Some(cap) = trim_cap_5 {
            if steps_taken >= cap {
                break;
            }
        }
        if extended_ref_start == 0 {
            break;
        }
        let candidate_ref_pos = extended_ref_start - 1;
        let nucleotide = sim
            .pool
            .get(NucHandle::new(seq_pos))
            .expect("extension range must point into the nucleotide pool");
        let Some(evidence) =
            segment_index.compatible_alleles_at(candidate_ref_pos as usize, nucleotide.base)
        else {
            break;
        };
        if evidence.allele_ids.is_empty() {
            break;
        }
        // only extend when this byte strictly narrows the
        // current max-score tie set. If every current candidate
        // matches the byte (or none of them do), the extension would
        // not resolve any ambiguity — leave the structural boundary
        // alone and stop walking.
        if !extension_narrows_tie_set(state.scores, &evidence.allele_ids) {
            break;
        }
        evidence.allele_ids.for_each_id(|id| {
            let slot = &mut state.scores[id.as_usize()];
            *slot = slot.saturating_add(1);
        });
        extended_seq_start = seq_pos;
        extended_ref_start = candidate_ref_pos;
        if evidence.informative {
            extended_inform = extended_inform.saturating_add(1);
        } else {
            extended_wildcard = extended_wildcard.saturating_add(1);
        }
        if seq_pos < np_start {
            crossed_into_overlap = true;
        }
        extension_added = true;
        steps_taken = steps_taken.saturating_add(1);
    }

    if extension_added {
        *state.seq_start = extended_seq_start;
        *state.ref_start = extended_ref_start;
        *state.informative_matches = extended_inform;
        *state.wildcard_matches = extended_wildcard;
        state.flags.insert(HypothesisFlags::BOUNDARY_ELASTIC);
        if crossed_into_overlap {
            state.flags.insert(HypothesisFlags::OVERLAPS_OTHER_SEGMENT);
        }
    }
}

pub(crate) fn walk_right_extension(
    sim: &Simulation,
    segment_index: &SegmentRefIndex,
    np_region: &Region,
    trim_cap_3: Option<u32>,
    state: &mut ExtensionWalkState<'_>,
) {
    if np_region.start.index() != *state.seq_end {
        return;
    }

    let np_end = np_region.end.index();
    let pool_len = sim.pool.len() as u32;
    let mut extended_seq_end = *state.seq_end;
    let mut extended_ref_pos = *state.ref_end;
    let mut extended_inform = *state.informative_matches;
    let mut extended_wildcard = *state.wildcard_matches;
    let mut extension_added = false;
    let mut crossed_into_overlap = false;
    let mut steps_taken: u32 = 0;

    let upper_bound: u32 = sim
        .sequence
        .regions
        .iter()
        .find(|r| {
            matches!(r.segment, Segment::V | Segment::D | Segment::J) && r.start == np_region.end
        })
        .map(|r| r.end.index())
        .unwrap_or(pool_len)
        .min(pool_len);

    for seq_pos in np_region.start.index()..upper_bound {
        if let Some(cap) = trim_cap_3 {
            if steps_taken >= cap {
                break;
            }
        }
        let nucleotide = sim
            .pool
            .get(NucHandle::new(seq_pos))
            .expect("extension range must point into the nucleotide pool");
        let Some(evidence) =
            segment_index.compatible_alleles_at(extended_ref_pos as usize, nucleotide.base)
        else {
            break;
        };
        if evidence.allele_ids.is_empty() {
            break;
        }
        // see commentary in `walk_left_extension`.
        if !extension_narrows_tie_set(state.scores, &evidence.allele_ids) {
            break;
        }
        evidence.allele_ids.for_each_id(|id| {
            let slot = &mut state.scores[id.as_usize()];
            *slot = slot.saturating_add(1);
        });
        extended_seq_end = seq_pos.saturating_add(1);
        extended_ref_pos = extended_ref_pos.saturating_add(1);
        if evidence.informative {
            extended_inform = extended_inform.saturating_add(1);
        } else {
            extended_wildcard = extended_wildcard.saturating_add(1);
        }
        if extended_seq_end > np_end {
            crossed_into_overlap = true;
        }
        extension_added = true;
        steps_taken = steps_taken.saturating_add(1);
    }

    if extension_added {
        *state.seq_end = extended_seq_end;
        *state.ref_end = extended_ref_pos;
        *state.informative_matches = extended_inform;
        *state.wildcard_matches = extended_wildcard;
        state.flags.insert(HypothesisFlags::BOUNDARY_ELASTIC);
        if crossed_into_overlap {
            state.flags.insert(HypothesisFlags::OVERLAPS_OTHER_SEGMENT);
        }
    }
}
