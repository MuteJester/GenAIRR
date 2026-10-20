//! Persistent ancestry of the original-coding-end retained-outcome model.
//!
//! Pool positions are observations, not identities: deletions mark an origin
//! absent while retaining its original allele, template position and choice.

use crate::address::PEnd;
use crate::assignment::SegmentOrientation;
use crate::refdata::AlleleId;

use super::Segment;

#[derive(Clone, Debug, Eq, PartialEq)]
pub struct CodingEndChoice {
    pub choice_id: usize,
    pub segment: Segment,
    pub allele_id: AlleleId,
    pub five: i64,
    pub three: i64,
    pub outcome_index: usize,
}

#[derive(Clone, Debug, Eq, PartialEq)]
pub struct CodingEndPOrigin {
    pub choice_id: usize,
    pub source_segment: Segment,
    pub allele_id: AlleleId,
    /// Assembled junction side; D ends swap under inversion.
    pub end: PEnd,
    pub orientation: SegmentOrientation,
    /// Zero-based position in the original, forward reference allele.
    pub template_position: u32,
    pub original_base: u8,
    /// Position in the coding-oriented pool, or absent after removal.
    pub query_position: Option<u32>,
}

#[derive(Clone, Debug, Default, Eq, PartialEq)]
pub struct CodingEndProvenance {
    pub choices: Vec<CodingEndChoice>,
    pub p_origins: Vec<CodingEndPOrigin>,
}

impl CodingEndProvenance {
    pub(crate) fn shifted_for_replacement(
        &self,
        start: u32,
        end: u32,
        replacement_len: u32,
    ) -> Self {
        let mut next = self.clone();
        let delta = i64::from(replacement_len) - i64::from(end - start);
        for origin in &mut next.p_origins {
            let Some(position) = origin.query_position else {
                continue;
            };
            if position >= start && position < end {
                origin.query_position = None;
            } else if position >= end {
                origin.query_position = Some((i64::from(position) + delta) as u32);
            }
        }
        next
    }
}
