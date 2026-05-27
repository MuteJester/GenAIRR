//! Live V/D/J call evidence state.
//!
//! Data structures plus small, testable operations that back the
//! dynamic allele-call layer.

use crate::ir::Segment;

mod bitset;
mod call;
pub(crate) mod dirty_signal_observer;
mod model;
mod reference_index;
mod refresh_hook;
mod refresh_plan;
pub(crate) mod scoring;
mod walker;
pub(crate) mod walker_observer;

pub use bitset::AlleleBitSet;
pub use call::{assembled_segment_live_call, with_assembled_segment_live_call};
pub use model::{
    BoundarySummary, BoundaryValue, DirtyLog, DirtyReason, DirtyWindow, EvidenceScore,
    HypothesisFlags, LiveCallConfidence, PlacementHypothesis, SegmentCalls, SegmentLiveCall,
};
pub use reference_index::{
    BaseBitsets, BaseEvidence, BoundaryIndex, IndexedAllele, KmerHit, KmerIndex,
    ReferenceMatchIndex, SegmentRefIndex, DEFAULT_REFERENCE_KMER_LEN,
};
pub use refresh_hook::LiveCallRefreshHook;
pub use refresh_plan::{LiveCallRefreshPlan, LiveCallRefreshStep};

fn assert_live_segment(segment: Segment) {
    assert!(
        matches!(segment, Segment::V | Segment::D | Segment::J),
        "live allele calls are only defined for V/D/J segments, got {:?}",
        segment
    );
}

#[cfg(test)]
mod tests;
