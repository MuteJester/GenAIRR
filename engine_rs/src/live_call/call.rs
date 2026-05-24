use super::{assert_live_segment, ReferenceMatchIndex, SegmentLiveCall};
use crate::ir::{Region, Segment, Simulation};

use super::walker::call_from_region;

/// Return a simulation with the live call for an assembled V/D/J
/// structural region refreshed from exact current-base evidence.
///
/// Discovery is restricted to the existing structural region. It
/// does not extend into NP evidence, repair indel-shifted
/// hypotheses, or alter AIRR projection.
///
/// **Fast path (streaming-walker observer).** When the
/// preceding `AssembleSegmentPass` ran with a `WalkerObserverState`
/// attached, it produced the segment's `SegmentLiveCall` inline with
/// the per-base push loop and stashed it on `sim.live_calls` via
/// `LiveCallState::with_segment_call_observed` *without* bumping the
/// `version`. We detect that pre-staged call here by checking
/// `evidence_version == base_state.version + 1` and absorb it by
/// performing the version bump that the observer path deliberately
/// skipped. The slow path (full `call_from_region` re-walk) remains
/// for callers that did not run the observer (e.g. test fixtures
/// that build `Simulation` directly without going through
/// `AssembleSegmentPass`).
pub fn with_assembled_segment_live_call(
    sim: &Simulation,
    reference_index: &ReferenceMatchIndex,
    segment: Segment,
) -> Simulation {
    assert_live_segment(segment);
    let base_state = sim.live_calls.clone().unwrap_or_default();
    let evidence_version = base_state.version.saturating_add(1);

    // Fast path: did `AssembleSegmentPass` already populate the
    // segment's call via the streaming walker observer?
    if let Some(existing) = base_state.get(segment) {
        if existing.evidence_version == evidence_version {
            // Observer-staged call detected. Absorb it by bumping
            // the version so subsequent passes see a consistent
            // monotonic version trajectory.
            return sim.with_live_calls(base_state.with_segment_call(existing.clone()));
        }
    }

    let Some(call) = assembled_segment_live_call(sim, reference_index, segment, evidence_version)
    else {
        return sim.clone();
    };
    sim.with_live_calls(base_state.with_segment_call(call))
}

/// Build the exact live call implied by the latest structural region
/// for `segment`.
pub fn assembled_segment_live_call(
    sim: &Simulation,
    reference_index: &ReferenceMatchIndex,
    segment: Segment,
    evidence_version: u64,
) -> Option<SegmentLiveCall> {
    assert_live_segment(segment);
    let segment_index = reference_index.get(segment)?;
    let region = latest_region_for_segment(sim, segment)?;
    let left_extension = left_extension_region_for(sim, segment, region);
    let right_extension = right_extension_region_for(sim, segment, region);
    Some(call_from_region(
        sim,
        segment_index,
        region,
        left_extension,
        right_extension,
        evidence_version,
    ))
}

fn latest_region_for_segment(sim: &Simulation, segment: Segment) -> Option<&Region> {
    sim.sequence
        .regions
        .iter()
        .rev()
        .find(|region| region.segment == segment)
}

/// Pick the NP region (if any) whose bases sit immediately to the right
/// of `segment`'s assembled region.
///
/// - V -> NP1.
/// - D -> NP2.
/// - J never has a right neighbour in either VJ or VDJ chain layouts.
fn right_extension_region_for<'a>(
    sim: &'a Simulation,
    segment: Segment,
    region: &Region,
) -> Option<&'a Region> {
    match segment {
        Segment::V => sim
            .sequence
            .regions
            .iter()
            .find(|r| r.segment == Segment::Np1 && r.start == region.end),
        Segment::D => sim
            .sequence
            .regions
            .iter()
            .find(|r| r.segment == Segment::Np2 && r.start == region.end),
        Segment::J => None,
        Segment::Np1 | Segment::Np2 => None,
    }
}

/// Pick the NP region (if any) whose bases sit immediately to the left
/// of `segment`'s assembled region.
///
/// - J -> adjacent NP region: NP1 in VJ chains (V -> NP1 -> J), NP2 in
///   VDJ chains (... -> NP2 -> J). Both cases collapse to "the NP region
///   whose `end` equals `region.start`", so we don't need to know the
///   chain type.
/// - D -> NP1.
/// - V never has a left neighbour in our DSL (V is always the first
///   region).
fn left_extension_region_for<'a>(
    sim: &'a Simulation,
    segment: Segment,
    region: &Region,
) -> Option<&'a Region> {
    match segment {
        Segment::J => sim
            .sequence
            .regions
            .iter()
            .find(|r| matches!(r.segment, Segment::Np1 | Segment::Np2) && r.end == region.start),
        Segment::D => sim
            .sequence
            .regions
            .iter()
            .find(|r| r.segment == Segment::Np1 && r.end == region.start),
        Segment::V => None,
        Segment::Np1 | Segment::Np2 => None,
    }
}
