use super::super::{HypothesisFlags, SegmentRefIndex};
use crate::ir::{NucHandle, Region, Segment, Simulation};

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
