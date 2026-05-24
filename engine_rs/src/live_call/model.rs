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

/// V/D/J live-call evidence + monotonic version chain.
///
/// **Storage shape — typestate for the fast-path absorb.** Each
/// V/D/J segment has TWO slots: a `committed` call (the consistent
/// state) and an optional `staged` call (a pending observer-produced
/// call waiting for absorption). The streaming walker observer writes
/// to `staged`; the post-pass refresh hook commits via
/// [`Self::commit_staged`].
///
/// This replaces the previous "`evidence_version == base.version + 1`"
/// arithmetic test that detected staged calls. The version-arithmetic
/// contract was a fragile convention: it required
/// [`crate::ir::SimulationBuilder::seal_with_committed_live_calls`]
/// to stamp staged calls with exactly the version offset the absorb
/// would test for, and required the post-pass refresh to iterate
/// V → D → J in lockstep. Both invariants lived in prose. The
/// staged/committed split makes "this call was produced by an
/// observer and hasn't been committed yet" a structural property of
/// the type — the compiler enforces what the convention used to.
///
/// `version` still bumps on every commit and is surfaced to
/// downstream consumers (AIRR provenance, change-detection tests),
/// but it's no longer load-bearing for the absorb path.
#[derive(Clone, Debug, Default, Eq, PartialEq)]
pub struct SegmentCalls {
    /// Per-segment slots holding the consistent (committed) call for
    /// each V/D/J segment. Use [`Self::get`] / [`Self::iter`] to
    /// read; modifications flow through [`Self::with_segment_call`].
    committed: crate::ir::PerSegment<SegmentLiveCall>,
    /// Per-segment slots holding observer-staged calls awaiting
    /// commitment. A staged call shadows its committed counterpart on
    /// read — it's the freshest evidence the engine has seen for that
    /// segment.
    staged: crate::ir::PerSegment<SegmentLiveCall>,
    pub version: u64,
}

impl SegmentCalls {
    pub fn empty() -> Self {
        Self::default()
    }

    /// Return the freshest call for `segment`: the staged one if any
    /// (pending commit), otherwise the committed one. `None` if
    /// neither slot is populated, or `segment` is non-assignable.
    pub fn get(&self, segment: Segment) -> Option<&SegmentLiveCall> {
        if !Segment::assignable().contains(&segment) {
            return None;
        }
        self.staged
            .get(segment)
            .or_else(|| self.committed.get(segment))
    }

    /// Iterate populated `(segment, &call)` entries in canonical
    /// V → D → J order. Yields the freshest (staged-shadowed) view
    /// of each segment. Skips empty slots.
    pub fn iter(&self) -> impl Iterator<Item = (Segment, &SegmentLiveCall)> + '_ {
        Segment::assignable()
            .iter()
            .copied()
            .filter_map(move |seg| self.get(seg).map(|c| (seg, c)))
    }

    /// `true` when `segment` has a staged-but-uncommitted call.
    pub fn has_staged(&self, segment: Segment) -> bool {
        self.staged.contains(segment)
    }

    /// Commit `call` directly: write to the committed slot, clear any
    /// pending stage on the same segment, bump `version`.
    pub fn with_segment_call(&self, call: SegmentLiveCall) -> Self {
        assert!(
            Segment::assignable().contains(&call.segment),
            "SegmentCalls::with_segment_call: non-assignable segment {:?}",
            call.segment
        );
        let mut next = self.clone();
        next.staged.take(call.segment);
        next.committed.set(call.segment, call);
        next.version = next.version.saturating_add(1);
        next
    }

    /// Stash a pre-computed segment call in the staged slot. Does
    /// **not** bump `version` and does **not** touch the committed
    /// slot — the post-pass refresh hook commits via
    /// [`Self::commit_staged`].
    ///
    /// Used only by the streaming walker observer during
    /// `AssembleSegmentPass` / mutation passes. Direct in-place
    /// mutator (no clone) to amortise across the V → D → J staging
    /// loop in `SimulationBuilder::seal_with_committed_live_calls`.
    pub fn stage_segment_call(&mut self, call: SegmentLiveCall) {
        assert!(
            Segment::assignable().contains(&call.segment),
            "SegmentCalls::stage_segment_call: non-assignable segment {:?}",
            call.segment
        );
        self.staged.set(call.segment, call);
    }

    /// If `segment` has a staged call, commit it (move staged →
    /// committed, bump `version`) and return the new `SegmentCalls`.
    /// Returns `None` when nothing is staged for that segment so the
    /// caller can fall through to a from-scratch recompute path.
    pub fn commit_staged(&self, segment: Segment) -> Option<Self> {
        let staged_call = self.staged.get(segment)?.clone();
        Some(self.with_segment_call(staged_call))
    }
}

/// Inter-pass dirty-window message log. Producer:
/// [`crate::live_call::dirty_signal_observer::DirtySignalObserver`].
/// Consumer: [`crate::live_call::LiveCallRefreshHook`]. The
/// post-pass refresh in `compiled/execute.rs` reads the windows to
/// narrow which V/D/J segments need recomputation, then drains via
/// [`crate::ir::Simulation::with_dirty_log`].
///
/// Was previously a `Vec<DirtyWindow>` field inside `LiveCallState`;
/// pulled into its own sidecar so it stops version-coupling to the
/// segment calls and stops cluttering call sites that only want to
/// touch the mutation counter or the V/D/J evidence.
#[derive(Clone, Debug, Default, Eq, PartialEq)]
pub struct DirtyLog {
    pub windows: Vec<DirtyWindow>,
}

impl DirtyLog {
    pub fn empty() -> Self {
        Self::default()
    }

    pub fn is_empty(&self) -> bool {
        self.windows.is_empty()
    }

    pub fn len(&self) -> usize {
        self.windows.len()
    }

    pub fn push(&mut self, window: DirtyWindow) {
        self.windows.push(window);
    }

    pub fn extend<I: IntoIterator<Item = DirtyWindow>>(&mut self, windows: I) {
        self.windows.extend(windows);
    }

    pub fn clear(&mut self) {
        self.windows.clear();
    }
}
