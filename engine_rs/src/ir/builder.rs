//! `SimulationBuilder` — owns an in-progress `Simulation` by value and
//! supports per-base appends without per-step refcount drama.
//!
//! ## Why this exists
//!
//! `Simulation::with_nucleotide_pushed` is a persistent operation: it
//! takes `&self`, clones the inner pool's Arc, and CoWs the inner `Vec`
//! through `Arc::make_mut`. When that clone is immediately reassigned
//! to a new local (`current = current.with_nucleotide_pushed(...)`),
//! the prior `current` is *briefly co-owned* during the assignment,
//! so `Arc::make_mut` falls into the deep-clone branch every iteration
//! of the NP draw loop.
//!
//! `SimulationBuilder` solves this by:
//!
//! 1. Taking ownership of a `Simulation` by value at construction.
//! 2. Calling `Arc::make_mut` once up front on the pool's inner `Vec`
//!    so the Arc refcount becomes 1.
//! 3. Holding the now-unique `Simulation` across many `push_nucleotide`
//!    calls — each of which is a direct `Vec::push` (the Arc refcount
//!    stays at 1, so subsequent `make_mut` calls hit the cheap path).
//!
//! ## What this *does not* do
//!
//! - It does not support `change_base_in_place`, `with_indel_*`, region
//!   updates, or any other operation beyond per-base append. The NP pass
//!   is the sole intended consumer for now.
//! - It does not break the persistent IR public contract — the only
//!   `&mut self` operation it exposes returns a `NucHandle` rather than
//!   a new `Simulation`. Callers seal the builder with `seal()` to obtain
//!   the immutable snapshot when done.
//!
//! ## Contract visibility
//!
//! `peek()` returns `&Simulation` whose `pool` reflects every base
//! pushed so far. Contracts reading `sim.pool` (e.g. `read_slot` in
//! `JunctionStopState`) see exactly the same committed prefix they
//! would under the old per-step `with_nucleotide_pushed` loop. This is
//! what preserves bit-identical regression output.

use super::event_log_observer::EventLogObserver;
use super::sim_event::{SimulationEvent, SimulationEventSink};
use super::{NucHandle, Nucleotide, Segment, Simulation};
use crate::contract::admit_mask_observer::ProductiveAdmitMaskObserver;
use crate::contract::JunctionStopState;
use crate::live_call::dirty_signal_observer::DirtySignalObserver;
use crate::live_call::walker_observer::{SealedWalkerState, WalkerObserverState};
use crate::live_call::SegmentRefIndex;

/// Mutable, by-value owner of an in-progress `Simulation` whose pool
/// supports cheap per-base append. See the module docs for the design
/// rationale.
///
/// The lifetime `'idx` is the lifetime of borrows held by attached
/// observers (a walker's `SegmentRefIndex`, an admit-mask's
/// `JunctionStopState`). Builders without any borrow-holding
/// observer can safely use `'static`.
///
/// ## Observer architecture
///
/// Every derived-state consumer on the builder rides one channel —
/// [`SimulationEventSink`]. Each observer kind has its own typed
/// slot on this struct (rather than living behind `dyn` in a
/// heterogeneous `Vec`) so the seal path can return that observer's
/// own typed output without downcasting and the broadcast loop
/// stays a series of static dispatches. The previous dual-channel
/// world (`IrEventObserver` + `SimulationEventSink`) collapsed in
/// the admit-mask migration; only the typed sinks remain.
pub struct SimulationBuilder<'idx> {
    simulation: Simulation,
    /// Recorder sink for the typed event stream — `Some` when a
    /// caller has opted in via [`Self::attach_event_log_observer`].
    /// Drained at seal time via
    /// [`Self::seal_event_log_observer`] into a flat
    /// `Vec<SimulationEvent>`.
    event_log_sink: Option<EventLogObserver>,
    /// Refresher sink — accumulates per-edit `DirtyWindow`s for
    /// segment-overlap gating in `apply_live_call_updates`.
    dirty_signal_sink: Option<DirtySignalObserver>,
    /// Derived-state sinks for live-call walking. Multiple
    /// concurrent walkers are allowed: mutation passes attach one
    /// per V/D/J region with a matching segment index.
    walker_sinks: Vec<WalkerObserverState<'idx>>,
    /// Derived-state sink for the productive-admit-mask
    /// fast path. NP-sampling attaches one per draw loop and
    /// queries the cached mask via [`Self::current_admit_mask`].
    admit_mask_sink: Option<ProductiveAdmitMaskObserver<'idx>>,
}

impl<'idx> SimulationBuilder<'idx> {
    /// Take ownership of `sim` and unique-ify its pool's inner `Vec`
    /// so subsequent `push_nucleotide` calls don't trigger a deep
    /// clone via `Arc::make_mut`.
    pub fn from_simulation(sim: Simulation) -> Self {
        let mut builder = Self {
            simulation: sim,
            event_log_sink: None,
            dirty_signal_sink: None,
            walker_sinks: Vec::new(),
            admit_mask_sink: None,
        };
        // Force unique ownership of the pool's nucleotide vector. After
        // this, every subsequent `Arc::make_mut` inside `pool.push`
        // is a refcount==1 fast path: it just hands back a
        // `&mut Vec<Nucleotide>` without copying.
        builder.simulation.pool.make_unique();
        builder
    }

    /// Broadcast a [`SimulationEvent`] to every attached
    /// [`SimulationEventSink`]. Each typed slot is checked in turn;
    /// missing slots are silent no-ops. This is the only fanout
    /// path — the legacy `IrEventObserver` trait was retired in
    /// the admit-mask slice.
    #[inline]
    fn broadcast_event(&mut self, event: SimulationEvent) {
        if let Some(sink) = self.event_log_sink.as_mut() {
            sink.on_event(&event);
        }
        if let Some(sink) = self.dirty_signal_sink.as_mut() {
            sink.on_event(&event);
        }
        for sink in self.walker_sinks.iter_mut() {
            sink.on_event(&event);
        }
        if let Some(sink) = self.admit_mask_sink.as_mut() {
            sink.on_event(&event);
        }
    }

    /// Read-only view of the in-progress simulation. Contracts and
    /// other observers see every base committed via `push_nucleotide`
    /// up to this point.
    pub fn peek(&self) -> &Simulation {
        &self.simulation
    }

    /// Attach a streaming live-call walker observer for one segment's
    /// assembly. `seq_start` is the pool index of the first base of
    /// the soon-to-be-assembled region — typically `peek().pool.len()`
    /// at the moment of attach.
    ///
    /// Only one walker observer may be attached at a time; calling
    /// this while another is active replaces it (the old state is
    /// dropped without being sealed, mirroring the test-fixture
    /// "no observer" path).
    pub(crate) fn attach_walker_observer(
        &mut self,
        segment_index: &'idx SegmentRefIndex,
        seq_start: u32,
    ) {
        self.walker_sinks
            .push(WalkerObserverState::new(segment_index, seq_start));
    }

    /// Attach a walker observer **rebuilt** from an already-assembled
    /// region.
    ///
    /// Walks the region's bytes in the current simulation through
    /// scoring (cost: O(region_len), same as a from-scratch
    /// `call_from_region` for that region) and installs the resulting
    /// state. Subsequent `change_base` events update the scores
    /// incrementally via `on_base_changed`.
    ///
    /// Used by post-assembly mutation passes (S5F, PCR, ...) so the
    /// post-pass `PassEffect::EditBases` refresh hits the live-call
    /// fast path instead of re-walking from scratch.
    pub(crate) fn attach_walker_observer_for_region(
        &mut self,
        segment_index: &'idx SegmentRefIndex,
        region: &super::Region,
    ) {
        self.walker_sinks
            .push(WalkerObserverState::from_existing_region(
                segment_index,
                &self.simulation,
                region,
            ));
    }

    /// Append one nucleotide and return its handle. Amortized O(1):
    /// no Arc clone, no `Vec` reallocation if capacity holds.
    ///
    /// Each attached [`SimulationEventSink`] is notified with a
    /// [`SimulationEvent::BasePushed`] *before* the pool's `push`
    /// consumes `n`. The handle passed to sinks equals the index
    /// `n` will occupy: `NucHandle::new(pool.len() as u32)` at the
    /// moment of call.
    pub fn push_nucleotide(&mut self, n: Nucleotide) -> NucHandle {
        let handle = NucHandle::new(self.simulation.pool.len() as u32);
        self.broadcast_event(SimulationEvent::BasePushed {
            handle,
            base: n.base,
            segment: n.segment,
            germline_pos: n.germline_pos.get(),
            flags: n.flags,
        });
        // The pool's `push` already does `Arc::make_mut` internally;
        // after `from_simulation`'s prepayment, that call is the cheap
        // refcount==1 path on every iteration. The returned handle
        // equals the one we computed above (Pool::push uses pool.len()
        // before the push to mint the handle, same as we did).
        self.simulation.pool.push(n)
    }

    /// Insert `n` into the pool at position `at`, shifting every
    /// existing handle `>= at` up by 1. Notifies every attached
    /// observer via `on_indel_inserted` before the underlying
    /// `Simulation::with_indel_inserted` applies the shift.
    ///
    /// Walker observers track the shift in place for external indels
    /// and mark themselves `needs_rebuild` for internal ones —
    /// `seal_all_walker_observers` rebuilds stale observers from the
    /// post-mutation pool.
    ///
    /// Indels are rare (typical workloads have 0-5 per record), so
    /// this delegates to the persistent `Simulation::with_indel_*`
    /// API rather than mutating the inner pool/sequence in place.
    /// The Arc-clones inside `with_indel_inserted` are negligible at
    /// this frequency.
    ///
    /// **Pass authors:** do not call this directly — route insertions
    /// through [`crate::passes::mutation_transaction::MutationTransaction::insert_base`],
    /// which validates the handle, fires observer events, and applies
    /// the reference-index-aware seal at `commit`. The only legitimate
    /// non-TX callers are the event-log replay path (`ir::event_log_observer`)
    /// and observer unit tests (`live_call::tests`).
    pub(crate) fn insert_indel(&mut self, at: u32, n: Nucleotide) {
        self.broadcast_event(SimulationEvent::IndelInserted {
            at,
            base: n.base,
            segment: n.segment,
            flags: n.flags,
        });
        // The per-region codon rail is not maintained in the hot
        // path; consumers compute on demand. The persistent layer's
        // recompute would be wasted work.
        self.simulation = self.simulation.with_indel_inserted(at, n);
    }

    /// Delete the nucleotide at position `at`, shifting every existing
    /// handle `> at` down by 1. Notifies observers of `on_indel_deleted`
    /// before the underlying pool mutation applies the shift.
    /// Returns the removed nucleotide (caller may ignore).
    ///
    /// **Pass authors:** do not call this directly — route deletions
    /// through [`crate::passes::mutation_transaction::MutationTransaction::delete_base`].
    /// Same restriction and rationale as [`Self::insert_indel`].
    pub(crate) fn delete_indel(&mut self, at: u32) -> Option<Nucleotide> {
        let removed = *self.simulation.pool.get(NucHandle::new(at))?;
        self.broadcast_event(SimulationEvent::IndelDeleted {
            at,
            removed_base: removed.base,
            segment: removed.segment,
        });
        self.simulation = self.simulation.with_indel_deleted(at);
        Some(removed)
    }

    // ──────────────────────────────────────────────────────────────
    // Non-pool consequences
    //
    // Each method below mirrors one persistent `Simulation::with_*`
    // mutation but additionally emits the matching
    // [`SimulationEvent`] before applying it. The method is the
    // only sanctioned way for a pass to make the change once
    // event emission is wired in — calling `sim.with_*` directly
    // bypasses every attached sink.
    //
    // No method here adds policy on top of the existing `with_*`:
    // they compute the pre-mutation `old` for sinks that want a
    // diff, fan out one event, then delegate. Stricter validation
    // (replay support / contract / feasibility) stays at the
    // calling pass.
    // ──────────────────────────────────────────────────────────────

    /// Install (or replace) the allele assignment for `segment` and
    /// emit [`SimulationEvent::AssignmentChanged`].
    pub(crate) fn assign_allele(
        &mut self,
        segment: Segment,
        instance: crate::assignment::AlleleInstance,
    ) {
        let old = self.simulation.assignments.get(segment).copied();
        self.broadcast_event(SimulationEvent::AssignmentChanged {
            segment,
            old,
            new: instance,
        });
        self.simulation = self.simulation.with_allele_assigned(segment, instance);
    }

    /// Update one trim length on `segment` and emit
    /// [`SimulationEvent::TrimChanged`]. `old` is the prior value
    /// for the *same* (segment, end) if the segment already had an
    /// allele assignment, otherwise `None`.
    pub(crate) fn update_trim(
        &mut self,
        segment: Segment,
        end: crate::assignment::TrimEnd,
        value: u16,
    ) {
        let old = self.simulation.assignments.get(segment).map(|inst| match end {
            crate::assignment::TrimEnd::Five => inst.trim_5,
            crate::assignment::TrimEnd::Three => inst.trim_3,
        });
        self.broadcast_event(SimulationEvent::TrimChanged {
            segment,
            end,
            old,
            new: value,
        });
        self.simulation = self.simulation.with_trim(segment, end, value);
    }

    /// Append a new structural region to `Simulation::sequence` and
    /// emit [`SimulationEvent::RegionAdded`].
    #[allow(dead_code)]
    pub(crate) fn add_region(&mut self, region: super::Region) {
        self.broadcast_event(SimulationEvent::RegionAdded {
            region: region.clone(),
        });
        self.simulation = self.simulation.with_region_added(region);
    }

    /// Replace the existing region of `replacement.segment` with
    /// `replacement` and emit [`SimulationEvent::RegionReplaced`]
    /// carrying the prior region.
    ///
    /// Panics if no region exists for `replacement.segment` —
    /// matches `Simulation::with_region_replaced_for_segment`'s
    /// expectation that the segment has a region to replace.
    #[allow(dead_code)]
    pub(crate) fn replace_region(&mut self, replacement: super::Region) {
        let old = self
            .simulation
            .sequence
            .regions
            .iter()
            .find(|r| r.segment == replacement.segment)
            .cloned()
            .expect("replace_region: no existing region for segment");
        self.broadcast_event(SimulationEvent::RegionReplaced {
            old,
            new: replacement.clone(),
        });
        self.simulation = self
            .simulation
            .with_region_replaced_for_segment(replacement);
    }

    /// Set the simulation's `n_mutations` sidecar to `new_count`
    /// and emit [`SimulationEvent::MutationCountChanged`] with the
    /// signed delta pre-computed.
    #[allow(dead_code)]
    pub(crate) fn set_mutation_count(&mut self, new_count: u32) {
        let old = self.simulation.mutation_count;
        let delta = (new_count as i32).saturating_sub(old as i32);
        self.broadcast_event(SimulationEvent::MutationCountChanged {
            old,
            new: new_count,
            delta,
        });
        self.simulation = self.simulation.with_mutation_count(new_count);
    }

    /// Record a reverse-complement decision. `RevCompPass` is a
    /// flag-only pass — the IR is **not** mutated; the AIRR
    /// projection picks up the trace flag and post-flips the
    /// sequence at projection time. The event is emitted so
    /// downstream sinks can observe the decision without scanning
    /// the trace.
    #[allow(dead_code)]
    pub(crate) fn record_reverse_complement_flag(&mut self, applied: bool) {
        self.broadcast_event(SimulationEvent::ReverseComplementFlagRecorded { applied });
    }

    /// Attach a structured event-log observer that captures every
    /// emitted [`SimulationEvent`] into a flat `Vec`. Drain via
    /// [`Self::seal_event_log_observer`] to obtain the captured
    /// stream. Used for trace replay, MCP audit, metrics.
    ///
    /// Rides the [`SimulationEventSink`] channel — see
    /// [`crate::ir::sim_event`] for the trace-vs-events split.
    #[allow(dead_code)]
    pub(crate) fn attach_event_log_observer(&mut self) {
        self.event_log_sink = Some(EventLogObserver::new());
    }

    /// Drain the attached event-log observer and return the
    /// captured events. Panics if no observer is attached.
    #[allow(dead_code)]
    pub(crate) fn seal_event_log_observer(&mut self) -> Vec<SimulationEvent> {
        let sink = self
            .event_log_sink
            .take()
            .expect("seal_event_log_observer called without a prior attach_event_log_observer");
        sink.seal()
    }

    /// Attach a dirty-signal observer that records a `DirtyWindow`
    /// for every base change / indel event. The post-pass
    /// `apply_live_call_updates` dispatch reads the captured windows
    /// to refresh only segments whose region overlaps a dirty
    /// position, rather than re-sweeping V/D/J unconditionally.
    ///
    /// Mutation passes (S5F, PCR, quality, contaminant, ncorrupt,
    /// uniform-mutation, indel, end-loss) opt in to this observer
    /// alongside their walker observers. Pure-assembly paths don't
    /// attach it: assembly is dispatched via
    /// `PassEffect::AssembleSegment` / `AppendRegion`, which have
    /// their own refresh logic in `apply_live_call_updates`.
    ///
    /// Rides the [`SimulationEventSink`] channel.
    pub(crate) fn attach_dirty_signal_observer(&mut self) {
        self.dirty_signal_sink = Some(DirtySignalObserver::new());
    }

    /// Change the base byte at `handle` in place and notify every
    /// attached observer of an `on_base_changed` event.
    ///
    /// Returns the `Nucleotide` as it was *before* the change.
    ///
    /// Companion to [`Self::push_nucleotide`]: observers see the old
    /// nucleotide via the event parameter (so they can compute score
    /// deltas, retranslate codons, etc.) *before* the pool reflects
    /// the new value. The new base is applied immediately after the
    /// notifications return.
    ///
    /// **Pass authors:** do not call this directly — route base
    /// substitutions through [`crate::passes::mutation_transaction::MutationTransaction::substitute_base`]
    /// (contract-filtered, trace-recording) or [`crate::passes::mutation_transaction::MutationTransaction::substitute_base_fixed`]
    /// (unconditional). The TX validates the handle, picks the
    /// reference-index-aware seal path, and emits structured
    /// `PassError`s for out-of-range handles instead of panicking
    /// through this method's `.expect()`. Legitimate non-TX callers
    /// are limited to the event-log replay path and observer unit tests.
    pub(crate) fn change_base(&mut self, handle: NucHandle, new_base: u8) -> Nucleotide {
        // Read the old nucleotide *before* mutating the pool — the
        // event observers need the old value with its segment and
        // germline_pos still populated.
        let old = *self
            .simulation
            .pool
            .get(handle)
            .expect("change_base: handle must point into the pool");

        self.broadcast_event(SimulationEvent::BaseChanged {
            handle,
            old_base: old.base,
            new_base,
            segment: old.segment,
            germline_pos: old.germline_pos.get(),
        });

        let prev = self.simulation.pool.change_base_in_place(handle, new_base);
        debug_assert_eq!(
            prev.base, old.base,
            "change_base: pool returned different old base than peek saw"
        );
        old
    }

    /// Drain the attached walker observer and return its sealed
    /// state. `seq_end` is the pool index just past the last pushed
    /// base — typically `peek().pool.len()` at the moment of seal.
    ///
    /// Panics if there isn't exactly one walker observer attached —
    /// call [`Self::seal_all_walker_observers`] for the multi-region
    /// case (mutation passes attaching one observer per existing
    /// V/D/J region).
    pub(crate) fn seal_walker_observer(&mut self, seq_end: u32) -> SealedWalkerState {
        let walkers = self.take_walker_observers();
        assert!(
            walkers.len() == 1,
            "seal_walker_observer expects exactly one attached observer; \
             got {}. Use seal_all_walker_observers for multi-region passes.",
            walkers.len()
        );
        walkers.into_iter().next().unwrap().seal(seq_end)
    }

    /// Drain every attached walker observer, returning one sealed
    /// state per attachment. Each observer is sealed with its own
    /// `seq_end_hint` — for observers built via
    /// `attach_walker_observer_for_region` this is the original
    /// region's end. Each sealed state carries its own segment +
    /// (for `Resolved`) `seq_start` / `seq_end`, so the caller can
    /// finalize it into a `SegmentLiveCall` and write back without
    /// further bookkeeping.
    ///
    /// Observers whose `needs_rebuild` flag is set (by an internal
    /// indel event) are rebuilt via `from_existing_region` against
    /// the current (post-mutation) simulation *before* sealing, so
    /// the sealed state reflects the post-indel pool rather than
    /// stale pre-indel scores.
    pub(crate) fn seal_all_walker_observers(&mut self) -> Vec<SealedWalkerState> {
        let walkers = self.take_walker_observers();
        let sim = &self.simulation;
        walkers
            .into_iter()
            .map(|obs| obs.rebuild_if_stale(sim))
            .map(|obs| {
                let seq_end = obs.seq_end_hint();
                obs.seal(seq_end)
            })
            .collect()
    }

    /// Drain every attached walker observer out of its typed sink
    /// slot, returning them in attach order. Other sinks stay
    /// attached. Used internally by
    /// [`Self::seal_walker_observer`] /
    /// [`Self::seal_all_walker_observers`].
    fn take_walker_observers(&mut self) -> Vec<WalkerObserverState<'idx>> {
        std::mem::take(&mut self.walker_sinks)
    }

    /// Take the attached dirty-signal sink (if any) out of its
    /// typed slot. The legacy observer registry is no longer
    /// involved; the sink rides the `SimulationEventSink` channel.
    /// Used internally by
    /// [`Self::seal_with_committed_live_calls`] to drain dirty
    /// windows; also exposed to in-crate tests that need to
    /// inspect the captured stream directly.
    pub(crate) fn take_dirty_signal_observer(&mut self) -> Option<DirtySignalObserver> {
        self.dirty_signal_sink.take()
    }

    /// Standard observer-attach pattern for the mutation passes
    /// (S5F, uniform, PCR, quality, ncorrupt, contaminant, indel,
    /// end-loss).
    ///
    /// When a `ReferenceMatchIndex` is available, attaches the
    /// dirty-signal observer (consumed by `apply_live_call_updates`
    /// in the compiled path) and a walker observer for every V/D/J
    /// region with a matching segment index. Without a
    /// `ReferenceMatchIndex` (test harnesses like
    /// `pass::testing::PassRuntime::execute_with_refdata`), this
    /// method is a no-op.
    ///
    /// Codon-rail data is *not* stored on `Region` — see the struct
    /// comment in [`crate::ir::region`]. No production consumer reads
    /// a cached rail (`NoStopCodonInJunction` reads pool bytes
    /// directly; AIRR projection re-translates from the raw
    /// sequence). Callers that need the rail call
    /// [`crate::ir::compute_codon_rail`] against the current pool at
    /// point of use.
    pub(crate) fn attach_standard_mutation_observers(
        &mut self,
        reference_index: Option<&'idx crate::live_call::ReferenceMatchIndex>,
    ) {
        let Some(ref_index) = reference_index else {
            return;
        };
        self.attach_dirty_signal_observer();
        let regions: Vec<super::Region> =
            self.simulation.sequence.regions.iter().cloned().collect();
        for region in regions
            .iter()
            .filter(|r| ref_index.get(r.segment).is_some())
        {
            let segment_index = ref_index.get(region.segment).unwrap();
            self.attach_walker_observer_for_region(segment_index, region);
        }
    }

    /// Seal walker observers, finalize each into a `SegmentLiveCall`,
    /// and stage the resulting live calls on the simulation with
    /// **monotonically increasing evidence_versions in V→D→J order**.
    ///
    /// The version sequence matches the order in which
    /// `compiled/execute.rs::apply_live_call_updates` processes
    /// `PassEffect::EditBases` (V first, then D, then J). With each
    /// staged call's `evidence_version` set to `base.version + 1`,
    /// `base.version + 2`, `base.version + 3` respectively, the live-
    /// call fast path fires for every segment in the post-pass
    /// refresh: each segment absorbs its staged call via
    /// [`with_segment_call`](crate::live_call::LiveCallState::with_segment_call)
    /// (which bumps the version) and the from-scratch
    /// `call_from_region` rebuild is skipped entirely.
    ///
    /// Mutation passes attach observers via
    /// [`Self::attach_standard_mutation_observers`] and finalize via
    /// this method.
    pub(crate) fn seal_with_committed_live_calls(
        mut self,
        ref_index: &crate::live_call::ReferenceMatchIndex,
    ) -> Simulation {
        let walker_sealed = self.seal_all_walker_observers();
        // Drain dirty signals captured during the pass so
        // `apply_live_call_updates` can use them to skip refreshes
        // for segments that received no edits.
        let dirty_windows = self
            .take_dirty_signal_observer()
            .map(|o| o.seal())
            .unwrap_or_default();
        let current = self.seal();

        // Take ownership of the SegmentCalls once, stage in place
        // across V → D → J, and reattach via a single
        // `with_segment_calls` at the end.
        //
        // Stage walker calls ONLY for segments whose region overlaps
        // a dirty window. Clean segments inherit their prior live
        // call unchanged.
        //
        // Each staged call's `evidence_version` is stamped with the
        // global version it'll have AFTER commit (base+1 for V,
        // base+2 for D, base+3 for J in the V → D → J commit order).
        // The post-pass `LiveCallRefreshHook` calls
        // `with_assembled_segment_live_call` per segment; the absorb
        // path detects staged calls structurally via
        // `SegmentCalls::commit_staged` — the `evidence_version`
        // stamp is now informational provenance, not the absorb's
        // control signal.
        let mut calls: crate::live_call::SegmentCalls = (*current.segment_calls).clone();
        let base_version = calls.version;
        let has_dirty_observer_signal = !dirty_windows.is_empty();
        let mut walker_sealed = walker_sealed;
        let mut staged_count: u64 = 0;
        for target_segment in Segment::assignable() {
            let position = match walker_sealed
                .iter()
                .position(|s| s.segment() == *target_segment)
            {
                Some(p) => p,
                None => continue,
            };
            if has_dirty_observer_signal
                && !segment_region_overlaps_dirty(&current, *target_segment, &dirty_windows)
            {
                walker_sealed.swap_remove(position);
                continue;
            }
            let sealed = walker_sealed.swap_remove(position);
            let segment_index = match ref_index.get(*target_segment) {
                Some(s) => s,
                None => continue,
            };
            let evidence_version = base_version.saturating_add(staged_count.saturating_add(1));
            staged_count = staged_count.saturating_add(1);
            let live_call = sealed.into_live_call(&current, segment_index, evidence_version);

            calls.stage_segment_call(live_call);
        }

        // Stash dirty windows on their own sidecar. Segment-call
        // updates and dirty-log updates are now independent — the
        // mutation pass touches both via two focused writes instead
        // of one clone-modify-reattach on the shared grab bag.
        let mut dirty_log: crate::live_call::DirtyLog = (*current.dirty_log).clone();
        dirty_log.extend(dirty_windows);

        current.with_segment_calls(calls).with_dirty_log(dirty_log)
    }

    /// Attach a productive-admit-mask observer that caches "which of
    /// `{A,C,G,T}` admit at the current NP slot" via the precomputed
    /// `JunctionStopState`. The sampler queries the mask via
    /// [`Self::current_admit_mask`] before each NP-base sample.
    ///
    /// Rides the [`SimulationEventSink`] channel.
    pub(crate) fn attach_admit_mask_observer(
        &mut self,
        state: &'idx JunctionStopState,
        np_segment: Segment,
    ) {
        self.admit_mask_sink = Some(ProductiveAdmitMaskObserver::new(state, np_segment));
    }

    /// 4-bit admit mask for the current NP slot, or `None` when no
    /// admit-mask observer is attached. Bit layout: bit 0 → `A`,
    /// bit 1 → `C`, bit 2 → `G`, bit 3 → `T`. See
    /// [`ProductiveAdmitMaskObserver::current_admit_mask`].
    pub(crate) fn current_admit_mask(&self) -> Option<u8> {
        self.admit_mask_sink
            .as_ref()
            .map(|am| am.current_admit_mask(&self.simulation))
    }

    /// Consume the builder and return the sealed `Simulation`.
    ///
    /// At seal time the inner pool still holds an `Arc<Vec<Nucleotide>>`
    /// with refcount==1; subsequent persistent operations on the sealed
    /// `Simulation` (e.g. `with_region_added`) clone the outer pool
    /// struct (cheap) without forcing a deep clone of the vector.
    ///
    /// If a walker observer is still attached, it is dropped without
    /// being sealed — equivalent to the "no observer" path.
    pub fn seal(self) -> Simulation {
        self.simulation
    }
}

/// Does the segment's assembled region (if any) overlap any of the
/// given dirty windows? Used by `seal_with_committed_live_calls` to
/// gate which walker observers stage their calls.
fn segment_region_overlaps_dirty(
    sim: &Simulation,
    segment: Segment,
    windows: &[crate::live_call::DirtyWindow],
) -> bool {
    let Some(region) = sim.sequence.regions.iter().find(|r| r.segment == segment) else {
        return false;
    };
    let region_start = region.start.index();
    let region_end = region.end.index();
    // Upper bound is INCLUSIVE: a dirty window whose start == region.end
    // represents an IndelDeleted at exactly the boundary byte the
    // deletion just shrunk out of the region. The pre-deletion position
    // was inside the region; the dirty signal must trigger a refresh
    // even though post-deletion the region.end equals the dirty site.
    // (Strict `<` here was the bug behind the IGK J residual: end-loss
    // 3' deletions at the J boundary silently slipped past the overlap
    // check, leaving the stale pre-deletion live call committed.)
    //
    // Lower bound stays strict (`>`): a dirty window at position
    // region.start - 1 represents a byte JUST BEFORE the region whose
    // deletion shifts the whole region down without changing its
    // content, so no refresh is needed.
    windows
        .iter()
        .any(|w| w.start <= region_end && w.end > region_start)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::ir::{flag, Segment};

    fn make_n(base: u8) -> Nucleotide {
        Nucleotide::synthetic(base, Segment::Np1, flag::N_NUC)
    }

    #[test]
    fn builder_pushes_are_visible_through_peek() {
        let mut b = SimulationBuilder::from_simulation(Simulation::new());
        assert_eq!(b.peek().pool.len(), 0);
        let h0 = b.push_nucleotide(make_n(b'A'));
        let h1 = b.push_nucleotide(make_n(b'C'));
        assert_eq!(h0, NucHandle::new(0));
        assert_eq!(h1, NucHandle::new(1));
        assert_eq!(b.peek().pool.len(), 2);
        assert_eq!(b.peek().pool.get(h0).unwrap().base, b'A');
        assert_eq!(b.peek().pool.get(h1).unwrap().base, b'C');
    }

    #[test]
    fn seal_yields_simulation_with_all_committed_bases() {
        let mut b = SimulationBuilder::from_simulation(Simulation::new());
        for &c in b"GATTACA" {
            b.push_nucleotide(make_n(c));
        }
        let sim = b.seal();
        assert_eq!(sim.pool.len(), 7);
        let bases: Vec<u8> = sim.pool.as_slice().iter().map(|n| n.base).collect();
        assert_eq!(&bases, b"GATTACA");
    }

    #[test]
    fn builder_uniqueifies_shared_pool_on_construction() {
        // Pre-share the pool's inner Vec, then build. The shared
        // reference must not see our pushes — make_mut at construction
        // must have split the buffer.
        let mut sim = Simulation::new();
        sim.pool.push(make_n(b'A'));
        let shared_pool = sim.pool.clone();
        assert_eq!(shared_pool.len(), 1);

        let mut b = SimulationBuilder::from_simulation(sim);
        b.push_nucleotide(make_n(b'T'));
        b.push_nucleotide(make_n(b'C'));

        // Shared snapshot is unchanged.
        assert_eq!(shared_pool.len(), 1);
        assert_eq!(shared_pool.get(NucHandle::new(0)).unwrap().base, b'A');

        // Builder's view reflects the new pushes.
        let sealed = b.seal();
        assert_eq!(sealed.pool.len(), 3);
        let bases: Vec<u8> = sealed.pool.as_slice().iter().map(|n| n.base).collect();
        assert_eq!(&bases, b"ATC");
    }

    // ──────────────────────────────────────────────────────────────
    // Non-pool consequence builder methods
    //
    // Each test attaches an `EventLogObserver`, calls one builder
    // method, seals, and asserts the recorded event matches the
    // expected payload byte-for-byte. With no observer attached
    // the method is a pure mutation — that side is exercised by
    // every existing pass test.
    // ──────────────────────────────────────────────────────────────

    use crate::assignment::{AlleleInstance, TrimEnd};
    use crate::ir::Region;
    use crate::refdata::AlleleId;

    fn drain_log(mut builder: SimulationBuilder<'static>) -> (Vec<SimulationEvent>, Simulation) {
        let events = builder.seal_event_log_observer();
        (events, builder.seal())
    }

    #[test]
    fn assign_allele_emits_assignment_changed_with_old_none_on_first_install() {
        let mut b = SimulationBuilder::from_simulation(Simulation::new());
        b.attach_event_log_observer();
        let inst = AlleleInstance::new(AlleleId::new(3));
        b.assign_allele(Segment::V, inst);
        let (events, sealed) = drain_log(b);

        assert_eq!(events.len(), 1);
        assert_eq!(
            events[0],
            SimulationEvent::AssignmentChanged {
                segment: Segment::V,
                old: None,
                new: inst,
            }
        );
        // Persistent mutation actually applied.
        assert_eq!(sealed.assignments.get(Segment::V).copied(), Some(inst));
    }

    #[test]
    fn assign_allele_emits_old_some_on_replacement() {
        let first = AlleleInstance::new(AlleleId::new(0));
        let replacement = AlleleInstance::new(AlleleId::new(5));
        let mut b = SimulationBuilder::from_simulation(
            Simulation::new().with_allele_assigned(Segment::V, first),
        );
        b.attach_event_log_observer();
        b.assign_allele(Segment::V, replacement);
        let (events, sealed) = drain_log(b);

        assert_eq!(events.len(), 1);
        assert_eq!(
            events[0],
            SimulationEvent::AssignmentChanged {
                segment: Segment::V,
                old: Some(first),
                new: replacement,
            }
        );
        assert_eq!(
            sealed.assignments.get(Segment::V).copied(),
            Some(replacement)
        );
    }

    #[test]
    fn update_trim_emits_trim_changed_carrying_old_value() {
        let inst = AlleleInstance::new(AlleleId::new(1));
        let mut b = SimulationBuilder::from_simulation(
            Simulation::new().with_allele_assigned(Segment::V, inst),
        );
        b.attach_event_log_observer();
        b.update_trim(Segment::V, TrimEnd::Three, 4);
        let (events, sealed) = drain_log(b);

        assert_eq!(events.len(), 1);
        assert_eq!(
            events[0],
            SimulationEvent::TrimChanged {
                segment: Segment::V,
                end: TrimEnd::Three,
                old: Some(0),
                new: 4,
            }
        );
        assert_eq!(sealed.assignments.get(Segment::V).unwrap().trim_3, 4);
    }

    #[test]
    fn add_region_emits_region_added_carrying_full_region() {
        let region = Region::new(Segment::V, NucHandle::new(0), NucHandle::new(10));
        let mut b = SimulationBuilder::from_simulation(Simulation::new());
        b.attach_event_log_observer();
        b.add_region(region.clone());
        let (events, sealed) = drain_log(b);

        assert_eq!(events.len(), 1);
        assert_eq!(
            events[0],
            SimulationEvent::RegionAdded {
                region: region.clone()
            }
        );
        assert_eq!(sealed.sequence.regions.len(), 1);
        assert_eq!(sealed.sequence.regions[0], region);
    }

    #[test]
    fn replace_region_emits_region_replaced_carrying_old_and_new() {
        let v_old = Region::new(Segment::V, NucHandle::new(0), NucHandle::new(10));
        let v_new = Region::new(Segment::V, NucHandle::new(0), NucHandle::new(13));
        let base = Simulation::new().with_region_added(v_old.clone());
        let mut b = SimulationBuilder::from_simulation(base);
        b.attach_event_log_observer();
        b.replace_region(v_new.clone());
        let (events, sealed) = drain_log(b);

        assert_eq!(events.len(), 1);
        assert_eq!(
            events[0],
            SimulationEvent::RegionReplaced {
                old: v_old,
                new: v_new.clone(),
            }
        );
        assert_eq!(sealed.sequence.regions[0], v_new);
    }

    #[test]
    fn set_mutation_count_emits_with_signed_delta() {
        let base = Simulation::new().with_mutation_count(7);
        let mut b = SimulationBuilder::from_simulation(base);
        b.attach_event_log_observer();
        b.set_mutation_count(3);
        let (events, sealed) = drain_log(b);

        assert_eq!(events.len(), 1);
        assert_eq!(
            events[0],
            SimulationEvent::MutationCountChanged {
                old: 7,
                new: 3,
                delta: -4,
            }
        );
        assert_eq!(sealed.mutation_count, 3);
    }

    #[test]
    fn record_reverse_complement_flag_emits_event_without_mutating_sim() {
        // Rev-comp is flag-only: the IR is unchanged but the event
        // must fire so downstream sinks can observe the decision.
        let base = Simulation::new();
        let mut b = SimulationBuilder::from_simulation(base.clone());
        b.attach_event_log_observer();
        b.record_reverse_complement_flag(true);
        let (events, sealed) = drain_log(b);

        assert_eq!(events.len(), 1);
        assert_eq!(
            events[0],
            SimulationEvent::ReverseComplementFlagRecorded { applied: true }
        );
        // Pool / sequence / assignments are byte-identical: no IR mutation.
        assert_eq!(sealed.pool.len(), base.pool.len());
        assert_eq!(sealed.sequence.regions.len(), base.sequence.regions.len());
        assert_eq!(sealed.mutation_count, base.mutation_count);
    }

    #[test]
    fn no_event_recorded_when_no_log_observer_is_attached() {
        // The builder methods are pure mutations when no sink is
        // attached — the per-method broadcast is a no-op.
        let mut b = SimulationBuilder::from_simulation(Simulation::new());
        b.assign_allele(Segment::V, AlleleInstance::new(AlleleId::new(0)));
        b.add_region(Region::new(Segment::V, NucHandle::new(0), NucHandle::new(1)));
        b.set_mutation_count(2);
        b.record_reverse_complement_flag(false);
        let sealed = b.seal();
        // Mutations still applied; nothing to assert about events
        // because there's no sink to inspect.
        assert!(sealed.assignments.get(Segment::V).is_some());
        assert_eq!(sealed.sequence.regions.len(), 1);
        assert_eq!(sealed.mutation_count, 2);
    }
}
