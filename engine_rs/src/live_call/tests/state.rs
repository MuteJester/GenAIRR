use super::super::{
    AlleleBitSet, BoundaryValue, DirtyLog, DirtyReason, DirtyWindow, EvidenceScore,
    HypothesisFlags, LiveCallConfidence, PlacementHypothesis, SegmentCalls, SegmentLiveCall,
};
use super::id;
use crate::ir::Segment;

#[test]
fn segment_live_call_summarizes_hypothesis_boundaries_and_alleles() {
    let h1 = PlacementHypothesis::new(
        Segment::V,
        10,
        30,
        0,
        20,
        AlleleBitSet::from_ids(6, [id(1), id(2)]),
        EvidenceScore::exact(20, 0),
        HypothesisFlags::EMPTY,
    );
    let h2 = PlacementHypothesis::new(
        Segment::V,
        10,
        33,
        0,
        23,
        AlleleBitSet::from_ids(6, [id(2), id(4)]),
        EvidenceScore::exact(22, 1),
        HypothesisFlags::BOUNDARY_ELASTIC,
    );

    let call = SegmentLiveCall::from_hypotheses(Segment::V, 6, vec![h1, h2], 7);

    assert_eq!(call.allele_call.to_ids(), vec![id(1), id(2), id(4)]);
    assert_eq!(call.confidence, LiveCallConfidence::ExactAmbiguous);
    assert_eq!(call.boundary_summary.seq_start, BoundaryValue::Single(10));
    assert_eq!(
        call.boundary_summary.seq_end,
        BoundaryValue::Ambiguous(vec![30, 33])
    );
    assert_eq!(call.boundary_summary.ref_start, BoundaryValue::Single(0));
    assert_eq!(
        call.boundary_summary.ref_end,
        BoundaryValue::Ambiguous(vec![20, 23])
    );
    assert_eq!(call.evidence_version, 7);
}

#[test]
fn segment_calls_updates_are_persistent() {
    let state = SegmentCalls::empty();
    let call = SegmentLiveCall::unresolved(Segment::J, 4);
    let with_call = state.with_segment_call(call);

    assert!(state.get(Segment::J).is_none());
    assert_eq!(state.version, 0);
    assert!(with_call.get(Segment::J).is_some());
    assert_eq!(with_call.version, 1);
}

#[test]
fn dirty_log_is_a_separate_sidecar() {
    let mut log = DirtyLog::empty();
    assert!(log.is_empty());
    log.push(DirtyWindow::new(
        5,
        8,
        DirtyReason::StructuralIndel { site: 6, delta: 1 },
    ));
    assert_eq!(log.len(), 1);
    log.clear();
    assert!(log.is_empty());
}

#[test]
#[should_panic(expected = "only defined for V/D/J")]
fn placement_hypothesis_rejects_np_segments() {
    let _ = PlacementHypothesis::new(
        Segment::Np1,
        0,
        1,
        0,
        1,
        AlleleBitSet::full(1),
        EvidenceScore::exact(1, 0),
        HypothesisFlags::EMPTY,
    );
}
