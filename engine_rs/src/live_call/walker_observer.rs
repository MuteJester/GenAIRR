//! Streaming live-call walker attached to `SimulationBuilder`.
//!
//! Today's `call_from_region` (see [`super::walker`]) runs as a
//! post-pass after `AssembleSegmentPass`: it re-walks the just-pushed
//! structural region position-by-position to accumulate per-allele
//! score totals, then runs the optional NP-region extension walks
//! and packages the result as a `SegmentLiveCall`.
//!
//! `WalkerObserverState` performs the **same scoring state machine**
//! incrementally while `AssembleSegmentPass` is still pushing bases
//! into the pool. Each `on_base_pushed` call mirrors one iteration of
//! `call_from_region`'s structural loop ([walker.rs:43-102]). At
//! `seal` time the observer hands back a `SealedWalkerState` that the
//! assembly pass feeds into the existing post-seal NP-extension walks
//! (`walk_left_extension`, `walk_right_extension`), producing exactly
//! the same `SegmentLiveCall` the from-scratch path would have built.
//!
//! ## Why
//!
//! Profiling on `human_tcrb` recombine+productive shows `call_from_region`
//! at 31.96% self-time, 37.27% inclusive — the single biggest cost after
//! the IR persistent-clone work. The structural loop *re-reads* the bases
//! the assembly pass just finished pushing; doing the scoring inline with
//! the push saves one full per-base walk per assembled segment.
//!
//! ## Bit-identical invariant
//!
//! The observer state machine matches `call_from_region`'s loop byte-for-byte:
//!
//! - Cross-segment nucleotide → `SealedWalkerState::Unsupported` (mirrors
//!   `unsupported_call` at walker.rs:58-60).
//! - `germline_pos.is_none()` (indel-insert) → skipped silently (mirrors
//!   walker.rs:66-68).
//! - Backwards `ref_pos` motion → `SealedWalkerState::Unsupported` (mirrors
//!   walker.rs:75-77).
//! - Forward `ref_pos` jump or first ref_pos → updates bookkeeping, runs
//!   `compatible_alleles_at`, OR-accumulates `scores` via `for_each_id`,
//!   bumps the informative/wildcard tally.
//! - Every base was indel-insert (no ref_pos ever set) → `SealedWalkerState::Unresolved`
//!   (mirrors walker.rs:104-109).
//!
//! The from-scratch walker remains the property-test oracle: the
//! `walker_observer_matches_call_from_region` property test in
//! [`super::tests`] asserts byte-for-byte equality across the full
//! existing test corpus.

use super::reference_index::SegmentRefIndex;
use super::walker::extensions::{walk_left_extension, walk_right_extension, ExtensionWalkState};
use super::{AlleleBitSet, EvidenceScore, HypothesisFlags, PlacementHypothesis, SegmentLiveCall};
use crate::ir::{
    GermlinePos, NucFlags, NucHandle, Nucleotide, Region, Segment, Simulation, SimulationEvent,
    SimulationEventSink,
};
use crate::refdata::AlleleId;

/// State a streaming observer accumulates as `SimulationBuilder`
/// pushes one base after another into the pool for an assembled
/// structural region.
///
/// One observer instance covers one segment's assembly. The
/// `SimulationBuilder` holds it as `Option<WalkerObserverState<'idx>>`
/// where `'idx` is the lifetime of the borrowed `SegmentRefIndex`.
pub(crate) struct WalkerObserverState<'idx> {
    /// Per-allele match scores accumulated as the segment's bases
    /// are pushed (matches the `scores` local in walker.rs:43).
    scores: Vec<u32>,
    /// Count of canonical (A/C/G/T) match positions seen so far
    /// (matches `informative_matches` in walker.rs:44).
    informative_matches: u32,
    /// Count of wildcard (`N`) match positions seen so far
    /// (matches `wildcard_matches` in walker.rs:45).
    wildcard_matches: u32,
    /// First germline position observed, set on the first base whose
    /// `germline_pos` is non-NONE (matches walker.rs:46, 79).
    ref_start: Option<u32>,
    /// Next expected germline position; advanced on every base with
    /// real provenance (matches walker.rs:47, 81).
    next_ref_pos: Option<u32>,
    /// Segment whose region this observer is attached to. Used to
    /// reject cross-segment nucleotides that would set `malformed`.
    segment: Segment,
    /// Number of alleles in `segment`'s reference pool. Carried so
    /// the seal step can build the right-sized `AlleleBitSet` and
    /// `Vec<u32>` for the empty-evidence / unsupported branches.
    allele_universe_len: usize,
    /// First pool index of the assembled region, captured at attach
    /// time before any push. The observer can't observe the
    /// post-push pool length so the assembly pass anchors this
    /// here. Matches `region.start.index()` in walker.rs:112.
    seq_start: u32,
    /// Highest pool index the observer has seen (= seq_start + N
    /// after N pushes). For assembly this grows by 1 per
    /// `on_base_pushed`; for `from_existing_region` it's initialised
    /// to `region.end.index()` so the bulk-seal API can read it
    /// without a separate caller-provided value.
    seq_end_seen: u32,
    /// `true` once the observer has seen a cross-segment nucleotide
    /// or a backwards `ref_pos` jump. Seal returns
    /// `SealedWalkerState::Unsupported` regardless of how many
    /// good bases preceded the bad one.
    malformed: bool,
    /// Compiled reference index for `segment`, borrowed from
    /// `PassContext::reference_index`. The `compatible_alleles_at`
    /// lookup on each push reads exactly the same data structure
    /// the from-scratch walker uses, so scoring is bit-identical.
    segment_index: &'idx SegmentRefIndex,
    /// set when an `on_indel_*` event arrives that the
    /// observer cannot patch in place — currently any deletion
    /// inside the region (which may move `ref_start` /
    /// `next_ref_pos` boundaries) and a future-safe default for
    /// edge cases. The builder rebuilds the observer from the
    /// post-indel simulation at seal time via
    /// [`Self::rebuild_if_stale`] when this flag is set.
    needs_rebuild: bool,
}

/// Result of sealing the streaming observer.
///
/// The three variants mirror the three exit paths from
/// `call_from_region` at walker.rs:104-228:
///
/// - `Unsupported` → walker.rs:59 / walker.rs:76, the malformed-IR
///   `unsupported_call` path.
/// - `Unresolved` → walker.rs:108, every base in the region was an
///   indel-insert with no germline provenance.
/// - `Resolved` → walker.rs:110-227, the normal path that runs the
///   extension walks and builds the `PlacementHypothesis`.
pub(crate) enum SealedWalkerState {
    Unsupported {
        segment: Segment,
        allele_universe_len: usize,
    },
    Unresolved {
        segment: Segment,
        allele_universe_len: usize,
    },
    Resolved(ResolvedWalkerState),
}

/// Owned per-allele scoring state ready to feed into the post-seal
/// NP-extension walks and the final `SegmentLiveCall` construction.
///
/// All fields shadow the locals in `call_from_region` at walker.rs:110-114
/// so the existing `ExtensionWalkState` (which takes `&'a mut` borrows
/// into each of these) can borrow into this struct unchanged.
pub(crate) struct ResolvedWalkerState {
    pub scores: Vec<u32>,
    pub informative_matches: u32,
    pub wildcard_matches: u32,
    pub ref_start: u32,
    pub ref_end: u32,
    pub seq_start: u32,
    pub seq_end: u32,
    pub flags: HypothesisFlags,
    pub segment: Segment,
    pub allele_universe_len: usize,
}

impl<'idx> WalkerObserverState<'idx> {
    /// Initialize a fresh observer for one segment's assembly.
    ///
    /// `seq_start` must be `sim.pool.len()` at the moment of attach —
    /// that is, the soon-to-be-assembled region's `region.start.index()`.
    /// The assembly pass anchors this once at attach time and the
    /// observer never inspects the pool again; everything it needs
    /// to know about positions arrives via `on_base_pushed`.
    pub(crate) fn new(segment_index: &'idx SegmentRefIndex, seq_start: u32) -> Self {
        let allele_universe_len = segment_index.allele_count();
        Self {
            scores: vec![0; allele_universe_len],
            informative_matches: 0,
            wildcard_matches: 0,
            ref_start: None,
            next_ref_pos: None,
            segment: segment_index.segment,
            allele_universe_len,
            seq_start,
            seq_end_seen: seq_start,
            malformed: false,
            segment_index,
            needs_rebuild: false,
        }
    }

    /// rebuild a walker observer's internal state from an
    /// already-assembled `Region` plus the current pool. Walks the
    /// region's bytes through `on_base_pushed` so the resulting
    /// score vector matches what `call_from_region`'s structural
    /// loop would produce on the same pool (extensions excluded —
    /// those re-run in `finalize_with_extensions` at seal time).
    ///
    /// Used by post-assembly mutation passes (S5F, PCR, …) that
    /// want to attach walker observers to a `SimulationBuilder`
    /// after assembly has finished, so subsequent `change_base`
    /// events can update scores incrementally instead of triggering
    /// the from-scratch `PassEffect::EditBases` post-pass refresh.
    ///
    /// Cost: O(region_len) re-walk at rebuild time. Same as a
    /// from-scratch `call_from_region` for the structural region.
    /// Net break-even per record when ≥1 mutation lands in the
    /// segment; net win when ≥2 mutations land (one rebuild + N
    /// O(matched_alleles) deltas vs. N full rebuilds).
    pub(crate) fn from_existing_region(
        segment_index: &'idx SegmentRefIndex,
        sim: &Simulation,
        region: &Region,
    ) -> Self {
        let mut obs = Self::new(segment_index, region.start.index());
        for seq_pos in region.start.index()..region.end.index() {
            let handle = NucHandle::new(seq_pos);
            // Safe: region range is by construction inside the pool.
            // If a handle were out of bounds the from-scratch walker
            // would panic too (walker.rs:49-53), so behaviour matches.
            let nuc = sim
                .pool
                .get(handle)
                .expect("from_existing_region: region range must point into the pool");
            obs.on_base_pushed(handle, nuc);
        }
        // After re-walking the region, the seq_end is region.end.
        // `on_base_pushed` doesn't update `seq_end_seen` (the field is
        // dedicated to bulk-seal-time read), so we set it explicitly.
        obs.seq_end_seen = region.end.index();
        obs
    }

    /// Highest pool index the observer has tracked. For
    /// `from_existing_region` this equals the region's end at attach
    /// time; for assembly use it's updated lazily by call sites
    /// (the single-observer seal path passes `seq_end` explicitly).
    /// Used by [`SimulationBuilder::seal_all_walker_observers`] to
    /// derive each observer's seal arg without external bookkeeping.
    pub(crate) fn seq_end_hint(&self) -> u32 {
        self.seq_end_seen
    }

    /// Replace the observer with a fresh
    /// `from_existing_region` rebuild against the post-mutation
    /// simulation if `needs_rebuild` is set. The rebuild walks the
    /// region currently assigned to this observer's segment in
    /// `sim`, which after an indel pass already reflects every
    /// `with_indel_inserted` / `with_indel_deleted` shift.
    pub(crate) fn rebuild_if_stale(self, sim: &Simulation) -> Self {
        if !self.needs_rebuild {
            return self;
        }
        let Some(region) = sim
            .sequence
            .regions
            .iter()
            .find(|r| r.segment == self.segment)
            .cloned()
        else {
            return self;
        };
        Self::from_existing_region(self.segment_index, sim, &region)
    }

    /// Notify the observer that one base just got committed to the
    /// pool. Mirrors one iteration of the for-loop in
    /// `call_from_region` at walker.rs:49-102.
    ///
    /// `_handle` is the index the nucleotide will occupy in the
    /// pool — the walker's score state is per-allele, not per-pool-
    /// position, so the parameter is ignored here. It exists on the
    /// `IrEventObserver` trait shape for 's codon-rail
    /// observer (which keys per-region amino-acid state by pool
    /// position).
    ///
    /// Returns immediately on the malformed-IR paths so subsequent
    /// pushes don't keep mutating `scores`. The malformed flag stays
    /// set and seal time produces `SealedWalkerState::Unsupported`.
    pub(crate) fn on_base_pushed(&mut self, _handle: NucHandle, nucleotide: &Nucleotide) {
        if self.malformed {
            return;
        }

        // Single segment check covers both branches below: cross-segment
        // nucleotides — synthetic (indel-insert) or otherwise — are a
        // malformed IR for this region and produce `unsupported_call`.
        // (walker.rs:58-60)
        if nucleotide.segment != self.segment {
            self.malformed = true;
            return;
        }

        // an indel-inserted nucleotide ends up inside V/D/J's
        // region with `germline_pos == GermlinePos::NONE`. It carries
        // no allele evidence (no germline byte to compare), so we
        // skip it without failing the call. (walker.rs:66-68)
        let Some(ref_pos) = nucleotide.germline_pos.get().map(|p| p as u32) else {
            return;
        };

        // ref_pos may *jump forward* if a base was deleted between
        // this position and the previous one. We allow gap-up but
        // still reject backwards motion, which would indicate a
        // genuinely-broken IR. (walker.rs:70-80)
        match self.next_ref_pos {
            Some(expected) if ref_pos < expected => {
                self.malformed = true;
                return;
            }
            Some(_) => {}
            None => self.ref_start = Some(ref_pos),
        }
        self.next_ref_pos = Some(ref_pos.saturating_add(1));

        // Score increment: every allele whose germline at this ref_pos
        // matches the observed (current) base picks up +1. (walker.rs:83-101)
        let Some(evidence) = self
            .segment_index
            .compatible_alleles_at(ref_pos as usize, nucleotide.base)
        else {
            return;
        };
        let scores = &mut self.scores;
        evidence.allele_ids.for_each_id(|id| {
            let slot = &mut scores[id.as_usize()];
            *slot = slot.saturating_add(1);
        });
        if evidence.informative {
            self.informative_matches = self.informative_matches.saturating_add(1);
        } else {
            self.wildcard_matches = self.wildcard_matches.saturating_add(1);
        }
    }

    /// Consume the observer at the end of region assembly.
    ///
    /// `seq_end` is the post-push pool length (the soon-to-be-built
    /// region's `region.end.index()`); the assembly pass passes it
    /// in because the observer doesn't see the builder's seal.
    pub(crate) fn seal(self, seq_end: u32) -> SealedWalkerState {
        if self.malformed {
            return SealedWalkerState::Unsupported {
                segment: self.segment,
                allele_universe_len: self.allele_universe_len,
            };
        }

        let Some(ref_start) = self.ref_start else {
            // Every position in the region was an indel-insert (no
            // germline_pos). Defer to the unresolved state — there is
            // no structural evidence to anchor a hypothesis.
            // (walker.rs:104-109)
            return SealedWalkerState::Unresolved {
                segment: self.segment,
                allele_universe_len: self.allele_universe_len,
            };
        };
        let ref_end = self
            .next_ref_pos
            .expect("non-empty region should set ref_end");

        SealedWalkerState::Resolved(ResolvedWalkerState {
            scores: self.scores,
            informative_matches: self.informative_matches,
            wildcard_matches: self.wildcard_matches,
            ref_start,
            ref_end,
            seq_start: self.seq_start,
            seq_end,
            flags: HypothesisFlags::EMPTY,
            segment: self.segment,
            allele_universe_len: self.allele_universe_len,
        })
    }
}

impl<'idx> WalkerObserverState<'idx> {
    /// Incremental score delta on a base change.
    ///
    /// When the byte at `handle` flips from `old_n.base` to
    /// `new_base`, the walker's per-allele scores need to be
    /// updated by *exactly* the set of alleles whose germline at
    /// `old_n.germline_pos` matched the old base (decrement) and
    /// the set whose germline matches the new base (increment).
    /// Alleles whose germline matched neither are unchanged.
    ///
    /// This is the architectural payoff of a base change no
    /// longer triggering a full from-scratch walker rebuild via
    /// `PassEffect::EditBases` + `call_from_region`. The observer
    /// updates its score vector in O(matched_alleles) time.
    ///
    /// **Scope filter.** Ignores events for nucleotides outside
    /// the observer's segment (the walker scores only its own
    /// segment). Also ignores events whose `germline_pos` is
    /// `None` — those are indel-inserted bases with no allele
    /// evidence (matches the `on_base_pushed` skip at
    /// walker.rs:66-68).
    pub(crate) fn on_base_changed(
        &mut self,
        _handle: NucHandle,
        old_n: &Nucleotide,
        new_base: u8,
    ) {
        if self.malformed {
            return;
        }
        if old_n.segment != self.segment {
            return;
        }
        let Some(ref_pos) = old_n.germline_pos.get().map(|p| p as u32) else {
            return;
        };
        if old_n.base == new_base {
            return;
        }

        // Decrement scores for alleles that matched the OLD base
        // at this ref_pos.
        if let Some(old_evidence) = self
            .segment_index
            .compatible_alleles_at(ref_pos as usize, old_n.base)
        {
            let scores = &mut self.scores;
            old_evidence.allele_ids.for_each_id(|id| {
                let slot = &mut scores[id.as_usize()];
                *slot = slot.saturating_sub(1);
            });
            if old_evidence.informative {
                self.informative_matches = self.informative_matches.saturating_sub(1);
            } else {
                self.wildcard_matches = self.wildcard_matches.saturating_sub(1);
            }
        }

        // Increment scores for alleles that match the NEW base at
        // this ref_pos.
        if let Some(new_evidence) = self
            .segment_index
            .compatible_alleles_at(ref_pos as usize, new_base)
        {
            let scores = &mut self.scores;
            new_evidence.allele_ids.for_each_id(|id| {
                let slot = &mut scores[id.as_usize()];
                *slot = slot.saturating_add(1);
            });
            if new_evidence.informative {
                self.informative_matches = self.informative_matches.saturating_add(1);
            } else {
                self.wildcard_matches = self.wildcard_matches.saturating_add(1);
            }
        }
    }

    /// Handle an indel insertion event.
    ///
    /// Inserted nucleotides are synthetic (`germline_pos == NONE`)
    /// — they contribute zero to every allele's score. So the
    /// walker can absorb internal insertions losslessly: bump
    /// `seq_end_seen` by 1. External insertions (before our
    /// region) shift both bounds.
    pub(crate) fn on_indel_inserted(&mut self, at: u32) {
        if self.malformed || self.needs_rebuild {
            return;
        }
        if at < self.seq_start {
            self.seq_start = self.seq_start.saturating_add(1);
            self.seq_end_seen = self.seq_end_seen.saturating_add(1);
        } else if at < self.seq_end_seen {
            // Strict `<`: an insertion at `at == seq_end_seen` lands
            // JUST PAST the walker's tracked range. `Sequence::with_indel_inserted`
            // doesn't grow `region.end` in this case either (the
            // `region.end > at` shift rule is strict), so the new
            // byte stays outside the segment's region. Bumping
            // `seq_end_seen` here would diverge the cached
            // hypothesis from a from-scratch recompute by one
            // position. (Caught by the live-call cache-parity
            // harness on IGH J under indel events.)
            self.seq_end_seen = self.seq_end_seen.saturating_add(1);
        }
    }

    /// Handle an indel deletion event.
    ///
    /// External deletions (before our region) shift both bounds.
    /// Internal deletions invalidate observer state: removing a
    /// boundary germline byte changes `ref_start` / `next_ref_pos`,
    /// and the in-place score decrement is correct but doesn't fix
    /// the boundary tracking. Mark `needs_rebuild` and defer to
    /// the seal-time `from_existing_region` rebuild.
    pub(crate) fn on_indel_deleted(&mut self, at: u32) {
        if self.malformed || self.needs_rebuild {
            return;
        }
        if at < self.seq_start {
            self.seq_start = self.seq_start.saturating_sub(1);
            self.seq_end_seen = self.seq_end_seen.saturating_sub(1);
        } else if at < self.seq_end_seen {
            self.needs_rebuild = true;
        }
    }
}

impl SimulationEventSink for WalkerObserverState<'_> {
    /// Dispatch incoming [`SimulationEvent`]s to the inherent
    /// handlers above. Reserved variants are ignored.
    ///
    /// The walker only inspects three fields from a nucleotide:
    /// `base`, `segment`, and `germline_pos`. For `BasePushed` and
    /// `BaseChanged` we reconstruct a transient `Nucleotide` from
    /// the event payload so the inherent methods can keep their
    /// `&Nucleotide` signatures (which are still used by
    /// `from_existing_region` reading directly from the pool).
    /// The `germline` byte and `flags` are filler — neither is
    /// read by the walker's scoring loop.
    fn on_event(&mut self, event: &SimulationEvent) {
        match *event {
            SimulationEvent::BasePushed {
                handle,
                base,
                segment,
                germline_pos,
                flags,
            } => {
                let n = Nucleotide {
                    base,
                    germline: base,
                    germline_pos: match germline_pos {
                        Some(g) => GermlinePos::pos(g),
                        None => GermlinePos::NONE,
                    },
                    segment,
                    flags,
                };
                self.on_base_pushed(handle, &n);
            }
            SimulationEvent::BaseChanged {
                handle,
                old_base,
                new_base,
                segment,
                germline_pos,
            } => {
                let old_n = Nucleotide {
                    base: old_base,
                    germline: old_base,
                    germline_pos: match germline_pos {
                        Some(g) => GermlinePos::pos(g),
                        None => GermlinePos::NONE,
                    },
                    segment,
                    flags: NucFlags::empty(),
                };
                self.on_base_changed(handle, &old_n, new_base);
            }
            SimulationEvent::IndelInserted { at, .. } => {
                self.on_indel_inserted(at);
            }
            SimulationEvent::IndelDeleted { at, .. } => {
                self.on_indel_deleted(at);
            }
            // Reserved (non-pool) variants describe sidecar state
            // — they don't move the walker's score vector or
            // boundary tracking.
            SimulationEvent::BaseDeleted { .. }
            | SimulationEvent::AssignmentChanged { .. }
            | SimulationEvent::TrimChanged { .. }
            | SimulationEvent::RegionAdded { .. }
            | SimulationEvent::RegionReplaced { .. }
            | SimulationEvent::ReverseComplementFlagRecorded { .. }
            | SimulationEvent::MutationCountChanged { .. } => {}
        }
    }
}

impl ResolvedWalkerState {
    /// Borrow the resolved state's fields as an `ExtensionWalkState`
    /// the existing `walk_left_extension` / `walk_right_extension`
    /// helpers can mutate in place. The layout matches walker.rs:157-166
    /// and walker.rs:172-181 exactly.
    pub(crate) fn as_extension_state(&mut self) -> ExtensionWalkState<'_> {
        ExtensionWalkState {
            scores: &mut self.scores,
            informative_matches: &mut self.informative_matches,
            wildcard_matches: &mut self.wildcard_matches,
            ref_start: &mut self.ref_start,
            ref_end: &mut self.ref_end,
            seq_start: &mut self.seq_start,
            seq_end: &mut self.seq_end,
            flags: &mut self.flags,
        }
    }

    /// Materialize the final `SegmentLiveCall` from the (extension-walked)
    /// resolved state. Mirrors walker.rs:185-227 byte-for-byte.
    pub(crate) fn into_segment_live_call(
        mut self,
        segment_index: &SegmentRefIndex,
        evidence_version: u64,
    ) -> SegmentLiveCall {
        if self.wildcard_matches > 0 {
            self.flags.insert(HypothesisFlags::HAS_WILDCARD_EVIDENCE);
        }

        // ── Build the tie-set ────────────────────────────────────────────
        // The candidate allele_call is the strict tie-set at max score —
        // alleles with strictly more matches win, alleles tied at the
        // maximum share the call. If no allele matched any position
        // (max_score == 0), every allele is equally consistent with the
        // absence of evidence — return the full pool. (walker.rs:189-209)
        let max_score = self.scores.iter().copied().max().unwrap_or(0);
        let allele_call = if max_score == 0 {
            segment_index.all_alleles.clone()
        } else {
            let mut bitset = AlleleBitSet::empty(self.allele_universe_len);
            for (idx, &score) in self.scores.iter().enumerate() {
                if score == max_score {
                    bitset.insert(AlleleId::new(idx as u32));
                }
            }
            bitset
        };

        let hypothesis = PlacementHypothesis::new(
            self.segment,
            self.seq_start,
            self.seq_end,
            self.ref_start,
            self.ref_end,
            allele_call,
            EvidenceScore::exact(self.informative_matches, self.wildcard_matches),
            self.flags,
        );

        SegmentLiveCall::from_hypotheses(
            self.segment,
            self.allele_universe_len,
            vec![hypothesis],
            evidence_version,
        )
    }
}

impl SealedWalkerState {
    /// Segment this sealed state belongs to. Used by the bulk-seal
    /// flow on `SimulationBuilder` to dispatch each sealed state to
    /// the matching V/D/J slot in `LiveCallState`.
    pub(crate) fn segment(&self) -> Segment {
        match self {
            SealedWalkerState::Unsupported { segment, .. } => *segment,
            SealedWalkerState::Unresolved { segment, .. } => *segment,
            SealedWalkerState::Resolved(r) => r.segment,
        }
    }

    /// Finalize this sealed state into a `SegmentLiveCall`, pulling
    /// the structural region's `seq_start` / `seq_end` from the
    /// state itself (for `Resolved`) instead of requiring the caller
    /// to thread them. Equivalent to
    /// [`Self::finalize_with_extensions`] with those args; provided
    /// as a convenience for callers (like S5F) that
    /// drain many sealed states and don't want to re-derive the
    /// coordinates per state.
    pub(crate) fn into_live_call(
        self,
        sim_in_progress: &Simulation,
        segment_index: &SegmentRefIndex,
        evidence_version: u64,
    ) -> SegmentLiveCall {
        let (seq_start, seq_end) = match &self {
            SealedWalkerState::Resolved(r) => (r.seq_start, r.seq_end),
            // Unused for the Unresolved / Unsupported branches, which
            // don't run the extension walks.
            _ => (0, 0),
        };
        self.finalize_with_extensions(
            sim_in_progress,
            segment_index,
            evidence_version,
            seq_start,
            seq_end,
        )
    }

    /// Run the post-seal NP-extension walks (against `sim_in_progress`)
    /// and materialize the final `SegmentLiveCall`.
    ///
    /// `sim_in_progress` is the simulation as the assembly pass has
    /// pushed the structural region but **not yet** added the
    /// `Region` to `sequence.regions`. The extension walks read
    /// adjacent NP-region nucleotides from `sim_in_progress.pool`
    /// (which they index by `NucHandle`, never by region membership)
    /// and look up adjacent NP regions via `sim_in_progress.sequence.regions`.
    /// At assembly time only the upstream NP region exists (e.g. NP1 is
    /// adjacent to D's region.start when D is assembled), so the right
    /// extension is typically a no-op — the eventual post-NP refresh
    /// path through `call_from_region` picks it up later.
    ///
    /// Mirrors walker.rs:155-228 exactly, including the trim-cap and
    /// boundary-region lookups.
    pub(crate) fn finalize_with_extensions(
        self,
        sim_in_progress: &Simulation,
        segment_index: &SegmentRefIndex,
        evidence_version: u64,
        seq_start: u32,
        seq_end: u32,
    ) -> SegmentLiveCall {
        match self {
            SealedWalkerState::Unsupported {
                segment,
                allele_universe_len,
            } => SegmentLiveCall::from_hypotheses(
                segment,
                allele_universe_len,
                Vec::new(),
                evidence_version,
            ),
            SealedWalkerState::Unresolved {
                segment,
                allele_universe_len,
            } => SegmentLiveCall::unresolved(segment, allele_universe_len),
            SealedWalkerState::Resolved(mut resolved) => {
                let seg = resolved.segment;
                // Trim-bounded extension caps. Mirrors walker.rs:129-130.
                let trim_cap_5 = sim_in_progress
                    .assignments
                    .get(seg)
                    .map(|a| a.trim_5 as u32);
                let trim_cap_3 = sim_in_progress
                    .assignments
                    .get(seg)
                    .map(|a| a.trim_3 as u32);

                // Construct a *synthetic* in-progress region matching
                // the bases the observer just saw. It's not added to
                // `sim_in_progress.sequence.regions` yet — the
                // extension-region lookups search by adjacent (NP)
                // region match, not by self-presence in the regions
                // list, so this synthetic region only needs the right
                // start/end coordinates to anchor the neighbour search.
                let region = Region::new(seg, NucHandle::new(seq_start), NucHandle::new(seq_end));
                let left_extension = left_extension_region_for(sim_in_progress, seg, &region);
                let right_extension = right_extension_region_for(sim_in_progress, seg, &region);

                if let Some(np_region) = left_extension {
                    let mut state = resolved.as_extension_state();
                    walk_left_extension(
                        sim_in_progress,
                        segment_index,
                        np_region,
                        trim_cap_5,
                        &mut state,
                    );
                }
                if let Some(np_region) = right_extension {
                    let mut state = resolved.as_extension_state();
                    walk_right_extension(
                        sim_in_progress,
                        segment_index,
                        np_region,
                        trim_cap_3,
                        &mut state,
                    );
                }

                resolved.into_segment_live_call(segment_index, evidence_version)
            }
        }
    }
}

/// Pick the NP region (if any) whose bases sit immediately to the right
/// of `segment`'s assembled region. Mirrors `right_extension_region_for`
/// in `live_call/call.rs` byte-for-byte.
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
/// of `segment`'s assembled region. Mirrors `left_extension_region_for`
/// in `live_call/call.rs` byte-for-byte.
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
