//! Live V/D/J call evidence state.
//!
//! Data structures plus small, testable operations that back the
//! dynamic allele-call layer.

use crate::ir::Segment;

mod bitset;
mod call;
mod model;
mod reference_index;
mod walker;

pub use bitset::AlleleBitSet;
pub use call::{assembled_segment_live_call, with_assembled_segment_live_call};
pub use model::{
    BoundarySummary, BoundaryValue, DirtyReason, DirtyWindow, EvidenceScore, HypothesisFlags,
    LiveCallConfidence, LiveCallState, PlacementHypothesis, SegmentLiveCall,
};
pub use reference_index::{
    BaseBitsets, BaseEvidence, BoundaryIndex, IndexedAllele, KmerHit, KmerIndex,
    ReferenceMatchIndex, SegmentRefIndex, DEFAULT_REFERENCE_KMER_LEN,
};

fn assert_live_segment(segment: Segment) {
    assert!(
        matches!(segment, Segment::V | Segment::D | Segment::J),
        "live allele calls are only defined for V/D/J segments, got {:?}",
        segment
    );
}

#[cfg(test)]
mod tests;
