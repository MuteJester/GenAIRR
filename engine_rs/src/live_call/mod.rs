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
mod walker;
pub(crate) mod walker_observer;

pub use bitset::AlleleBitSet;
pub use call::{assembled_segment_live_call, with_assembled_segment_live_call};
pub use model::{
    BoundarySummary, BoundaryValue, DirtyReason, DirtyWindow, EvidenceScore, HypothesisFlags,
    LiveCallConfidence, LiveCallState, PlacementHypothesis, SegmentLiveCall,
};
pub use refresh_hook::LiveCallRefreshHook;
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
