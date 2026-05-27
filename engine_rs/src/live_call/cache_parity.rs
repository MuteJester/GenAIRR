//! Live-call cache equivalence harness.
//!
//! The runtime path stages and commits V/D/J `SegmentLiveCall`s via
//! the walker observer's incremental state, with the
//! `LiveCallRefreshHook` re-running `assembled_segment_live_call`
//! when events demand it. The end result is *supposed* to equal a
//! from-scratch recomputation against the final `Simulation`, but
//! the incremental machinery has historically harboured boundary
//! bugs (e.g. the `segment_region_overlaps_dirty` strict-`<` that
//! the AIRR record validator exposed).
//!
//! This module compares the cached call on `sim.segment_calls`
//! against a fresh `assembled_segment_live_call` over the same sim
//! + refdata, per V/D/J segment, and reports any divergence. Use
//! it as an explicit cache-equivalence guard so stale-cache bugs
//! fail closer to the source instead of through the AIRR projection.
//!
//! The fresh recomputation builds its own `ReferenceMatchIndex` so
//! the harness is self-contained — callers don't need to thread the
//! simulator's index through.

use crate::ir::{Segment, Simulation};
use crate::refdata::{AlleleId, RefDataConfig};

use super::call::assembled_segment_live_call;
use super::ReferenceMatchIndex;

/// One per-segment parity result.
#[derive(Clone, Debug, PartialEq)]
pub struct SegmentParity {
    pub segment: Segment,
    pub tie_set_matches: bool,
    pub cached_tie_set: Vec<AlleleId>,
    pub fresh_tie_set: Vec<AlleleId>,
    pub cached_present: bool,
    pub fresh_present: bool,
    /// Hypothesis bounds equality. `None` when either side lacks a
    /// resolved hypothesis (e.g. unresolved live call). When both
    /// sides have a hypothesis, `Some(true)` iff every bound
    /// (`seq_start`, `seq_end`, `ref_start`, `ref_end`) matches.
    pub hypothesis_bounds_match: Option<bool>,
    pub cached_hypothesis: Option<HypothesisBounds>,
    pub fresh_hypothesis: Option<HypothesisBounds>,
}

#[derive(Copy, Clone, Debug, Eq, PartialEq)]
pub struct HypothesisBounds {
    pub seq_start: u32,
    pub seq_end: u32,
    pub ref_start: u32,
    pub ref_end: u32,
}

/// Check live-call cache parity for every assignable segment (V/D/J).
///
/// Returns one `SegmentParity` per segment that has either a cached
/// call OR a non-trivial fresh recomputation. Segments with no
/// assignment and no region are skipped (nothing to compare).
pub fn check_segment_calls_parity(
    sim: &Simulation,
    refdata: &RefDataConfig,
) -> Vec<SegmentParity> {
    let reference_index = ReferenceMatchIndex::build(refdata);
    let mut out = Vec::new();
    for segment in [Segment::V, Segment::D, Segment::J] {
        let cached = sim.segment_calls.get(segment);
        let evidence_version = sim.segment_calls.version.saturating_add(1);
        let fresh = assembled_segment_live_call(sim, &reference_index, segment, evidence_version);

        let cached_present = cached.is_some();
        let fresh_present = fresh.is_some();
        if !cached_present && !fresh_present {
            // Neither side has a call — segment isn't assembled.
            // Skip rather than emit a vacuous match record.
            continue;
        }

        let cached_tie_set: Vec<AlleleId> = cached
            .map(|c| c.allele_call.iter_ids().collect())
            .unwrap_or_default();
        let fresh_tie_set: Vec<AlleleId> = fresh
            .as_ref()
            .map(|c| c.allele_call.iter_ids().collect())
            .unwrap_or_default();

        let tie_set_matches = cached_present == fresh_present && cached_tie_set == fresh_tie_set;

        let cached_hypothesis = cached.and_then(hypothesis_bounds);
        let fresh_hypothesis = fresh.as_ref().and_then(|c| hypothesis_bounds(c));
        let hypothesis_bounds_match = match (cached_hypothesis, fresh_hypothesis) {
            (Some(a), Some(b)) => Some(a == b),
            _ => None,
        };

        out.push(SegmentParity {
            segment,
            tie_set_matches,
            cached_tie_set,
            fresh_tie_set,
            cached_present,
            fresh_present,
            hypothesis_bounds_match,
            cached_hypothesis,
            fresh_hypothesis,
        });
    }
    out
}

fn hypothesis_bounds(call: &super::SegmentLiveCall) -> Option<HypothesisBounds> {
    let h = call.hypotheses.first()?;
    Some(HypothesisBounds {
        seq_start: h.seq_start,
        seq_end: h.seq_end,
        ref_start: h.ref_start,
        ref_end: h.ref_end,
    })
}
