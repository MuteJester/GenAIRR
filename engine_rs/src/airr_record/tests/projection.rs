use super::super::build_airr_record;
use super::*;
use crate::live_call::{
    AlleleBitSet, EvidenceScore, HypothesisFlags, PlacementHypothesis, SegmentCalls,
    SegmentLiveCall,
};

#[test]
fn airr_call_projection_falls_back_to_origin_assignment_without_live_call() {
    let (cfg, sim, _v0, _v1) = call_projection_fixture();
    let outcome = outcome_from_sim(sim);
    let rec = build_airr_record(&outcome, &cfg, "fallback");

    assert_eq!(rec.v_call, "IGHV1-1*01");
    assert_eq!(rec.locus, "IGH");
}

#[test]
fn airr_call_projection_prefers_live_call_allele_set() {
    let (cfg, sim, v0, v1) = call_projection_fixture();
    let hypothesis = PlacementHypothesis::new(
        Segment::V,
        0,
        6,
        0,
        6,
        AlleleBitSet::from_ids(cfg.v_pool.len(), [v0, v1]),
        EvidenceScore::exact(6, 0),
        HypothesisFlags::EMPTY,
    );
    let call = SegmentLiveCall::from_hypotheses(Segment::V, cfg.v_pool.len(), vec![hypothesis], 1);
    let live = SegmentCalls::empty().with_segment_call(call);
    let outcome = outcome_from_sim(sim.with_segment_calls(live));

    let rec = build_airr_record(&outcome, &cfg, "live");

    assert_eq!(rec.v_call, "IGHV1-1*01,IGHV1-1*02");
    assert_eq!(rec.locus, "IGH");
}

#[test]
fn airr_call_projection_falls_back_to_truth_when_live_call_is_unsupported() {
    let (cfg, sim, _v0, _v1) = call_projection_fixture();
    let call = SegmentLiveCall::from_hypotheses(Segment::V, cfg.v_pool.len(), Vec::new(), 1);
    let live = SegmentCalls::empty().with_segment_call(call);
    let outcome = outcome_from_sim(sim.with_segment_calls(live));

    let rec = build_airr_record(&outcome, &cfg, "unsupported");

    assert_eq!(rec.v_call, "IGHV1-1*01");
}
