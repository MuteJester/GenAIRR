use super::{assert_live_segment, AlleleBitSet};
use crate::ir::Segment;

/// Evidence quality for one placement hypothesis.
#[derive(Copy, Clone, Debug, Default, Eq, PartialEq)]
pub struct EvidenceScore {
    pub conflicts: u32,
    pub informative_matches: u32,
    pub wildcard_matches: u32,
    pub span_len: u32,
}

impl EvidenceScore {
    pub const fn exact(informative_matches: u32, wildcard_matches: u32) -> Self {
        Self {
            conflicts: 0,
            informative_matches,
            wildcard_matches,
            span_len: informative_matches + wildcard_matches,
        }
    }

    pub const fn is_exact(self) -> bool {
        self.conflicts == 0
    }
}

/// Flags describing how a hypothesis was discovered.
#[derive(Copy, Clone, Debug, Default, Eq, PartialEq)]
pub struct HypothesisFlags(u32);

impl HypothesisFlags {
    pub const EMPTY: Self = Self(0);
    pub const BOUNDARY_ELASTIC: Self = Self(1 << 0);
    pub const OVERLAPS_OTHER_SEGMENT: Self = Self(1 << 1);
    pub const HAS_WILDCARD_EVIDENCE: Self = Self(1 << 2);

    pub const fn contains(self, flag: Self) -> bool {
        (self.0 & flag.0) == flag.0
    }

    pub fn insert(&mut self, flag: Self) {
        self.0 |= flag.0;
    }
}

/// One exact-equivalence placement hypothesis for a V/D/J segment.
#[derive(Clone, Debug, Eq, PartialEq)]
pub struct PlacementHypothesis {
    pub segment: Segment,
    pub seq_start: u32,
    pub seq_end: u32,
    pub ref_start: u32,
    pub ref_end: u32,
    pub allele_ids: AlleleBitSet,
    pub score: EvidenceScore,
    pub flags: HypothesisFlags,
}

impl PlacementHypothesis {
    #[allow(clippy::too_many_arguments)]
    pub fn new(
        segment: Segment,
        seq_start: u32,
        seq_end: u32,
        ref_start: u32,
        ref_end: u32,
        allele_ids: AlleleBitSet,
        score: EvidenceScore,
        flags: HypothesisFlags,
    ) -> Self {
        assert_live_segment(segment);
        assert!(seq_start <= seq_end, "hypothesis seq range is inverted");
        assert!(ref_start <= ref_end, "hypothesis ref range is inverted");
        Self {
            segment,
            seq_start,
            seq_end,
            ref_start,
            ref_end,
            allele_ids,
            score,
            flags,
        }
    }
}

#[derive(Clone, Debug, Eq, PartialEq)]
pub enum BoundaryValue {
    Unresolved,
    Single(u32),
    Ambiguous(Vec<u32>),
}

impl BoundaryValue {
    pub fn from_values(values: impl IntoIterator<Item = u32>) -> Self {
        let mut values: Vec<u32> = values.into_iter().collect();
        values.sort_unstable();
        values.dedup();

        match values.len() {
            0 => Self::Unresolved,
            1 => Self::Single(values[0]),
            _ => Self::Ambiguous(values),
        }
    }
}

#[derive(Clone, Debug, Eq, PartialEq)]
pub struct BoundarySummary {
    pub seq_start: BoundaryValue,
    pub seq_end: BoundaryValue,
    pub ref_start: BoundaryValue,
    pub ref_end: BoundaryValue,
}

impl BoundarySummary {
    pub fn unresolved() -> Self {
        Self {
            seq_start: BoundaryValue::Unresolved,
            seq_end: BoundaryValue::Unresolved,
            ref_start: BoundaryValue::Unresolved,
            ref_end: BoundaryValue::Unresolved,
        }
    }

    pub fn from_hypotheses(hypotheses: &[PlacementHypothesis]) -> Self {
        Self {
            seq_start: BoundaryValue::from_values(hypotheses.iter().map(|h| h.seq_start)),
            seq_end: BoundaryValue::from_values(hypotheses.iter().map(|h| h.seq_end)),
            ref_start: BoundaryValue::from_values(hypotheses.iter().map(|h| h.ref_start)),
            ref_end: BoundaryValue::from_values(hypotheses.iter().map(|h| h.ref_end)),
        }
    }
}

#[derive(Copy, Clone, Debug, Eq, PartialEq)]
pub enum LiveCallConfidence {
    Unresolved,
    Unsupported,
    ExactSingle,
    ExactAmbiguous,
}

/// Live evidence state for one segment.
#[derive(Clone, Debug, Eq, PartialEq)]
pub struct SegmentLiveCall {
    pub segment: Segment,
    pub hypotheses: Vec<PlacementHypothesis>,
    pub allele_call: AlleleBitSet,
    pub boundary_summary: BoundarySummary,
    pub confidence: LiveCallConfidence,
    pub evidence_version: u64,
}

impl SegmentLiveCall {
    pub fn unresolved(segment: Segment, allele_universe_len: usize) -> Self {
        assert_live_segment(segment);
        Self {
            segment,
            hypotheses: Vec::new(),
            allele_call: AlleleBitSet::empty(allele_universe_len),
            boundary_summary: BoundarySummary::unresolved(),
            confidence: LiveCallConfidence::Unresolved,
            evidence_version: 0,
        }
    }

    pub fn from_hypotheses(
        segment: Segment,
        allele_universe_len: usize,
        hypotheses: Vec<PlacementHypothesis>,
        evidence_version: u64,
    ) -> Self {
        assert_live_segment(segment);
        let mut allele_call = AlleleBitSet::empty(allele_universe_len);
        for hypothesis in &hypotheses {
            assert_eq!(
                hypothesis.segment, segment,
                "SegmentLiveCall cannot contain {:?} hypothesis for {:?}",
                hypothesis.segment, segment
            );
            allele_call.union_with(&hypothesis.allele_ids);
        }

        let confidence = if hypotheses.is_empty() {
            LiveCallConfidence::Unsupported
        } else if allele_call.len() == 1 {
            LiveCallConfidence::ExactSingle
        } else {
            LiveCallConfidence::ExactAmbiguous
        };

        let boundary_summary = BoundarySummary::from_hypotheses(&hypotheses);
        Self {
            segment,
            hypotheses,
            allele_call,
            boundary_summary,
            confidence,
            evidence_version,
        }
    }
}

/// Dirty live-call recomputation interval.
#[derive(Clone, Debug, Eq, PartialEq)]
pub struct DirtyWindow {
    pub start: u32,
    pub end: u32,
    pub reason: DirtyReason,
}

impl DirtyWindow {
    pub fn new(start: u32, end: u32, reason: DirtyReason) -> Self {
        assert!(start <= end, "dirty window range is inverted");
        Self { start, end, reason }
    }
}

#[derive(Clone, Debug, Eq, PartialEq)]
pub enum DirtyReason {
    AlleleSampled(Segment),
    TrimChanged(Segment),
    SegmentAssembled(Segment),
    NpGenerated(Segment),
    BaseEdited { site: u32 },
    StructuralIndel { site: u32, delta: i32 },
    ContaminantReplacement,
}

/// Dormant top-level live-call sidecar.
#[derive(Clone, Debug, Default, Eq, PartialEq)]
pub struct LiveCallState {
    pub v: Option<SegmentLiveCall>,
    pub d: Option<SegmentLiveCall>,
    pub j: Option<SegmentLiveCall>,
    pub dirty_windows: Vec<DirtyWindow>,
    pub version: u64,
}

impl LiveCallState {
    pub fn empty() -> Self {
        Self::default()
    }

    pub fn get(&self, segment: Segment) -> Option<&SegmentLiveCall> {
        match segment {
            Segment::V => self.v.as_ref(),
            Segment::D => self.d.as_ref(),
            Segment::J => self.j.as_ref(),
            Segment::Np1 | Segment::Np2 => None,
        }
    }

    pub fn with_segment_call(&self, call: SegmentLiveCall) -> Self {
        let mut next = self.clone();
        match call.segment {
            Segment::V => next.v = Some(call),
            Segment::D => next.d = Some(call),
            Segment::J => next.j = Some(call),
            Segment::Np1 | Segment::Np2 => unreachable!("SegmentLiveCall rejects NP segments"),
        }
        next.version = next.version.saturating_add(1);
        next
    }

    pub fn with_dirty_window(&self, window: DirtyWindow) -> Self {
        let mut next = self.clone();
        next.dirty_windows.push(window);
        next.version = next.version.saturating_add(1);
        next
    }
}
