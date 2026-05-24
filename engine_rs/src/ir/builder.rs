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

use super::event_log_observer::{EventLogObserver, IrEvent};
use super::{NucHandle, Nucleotide, Segment, Simulation};
use crate::contract::admit_mask_observer::ProductiveAdmitMaskObserver;
use crate::contract::JunctionStopState;
use crate::live_call::dirty_signal_observer::DirtySignalObserver;
use crate::live_call::walker_observer::{SealedWalkerState, WalkerObserverState};
use crate::live_call::SegmentRefIndex;

/// Reactive subscriber for IR mutation events emitted by
/// `SimulationBuilder`.
///
/// Today's events:
///
/// - [`on_base_pushed`](Self::on_base_pushed) — one nucleotide is
///   about to be appended to the pool.
/// - [`on_base_changed`](Self::on_base_changed) — the base byte at an
///   existing handle is about to be replaced with a new value (the
///   segment and germline_pos of the nucleotide do not change).
///
/// Indel events (`on_indel_inserted` / `on_indel_deleted`) fire on
/// the same broadcast surface — observers can either patch their
/// state in place or mark themselves stale for rebuild at seal time.
///
/// **Visibility.** `pub(crate)` deliberately: every observer lives in
/// the same crate (walker, productive admit-mask, event log, dirty
/// signal). External crates are not expected to implement this; if
/// they ever do, the trait can be widened to `pub` then.
///
/// **Ordering invariant.** Events fire *before* the underlying pool
/// mutation is applied. Observers receive the nucleotide values
/// directly via the event parameters; they must not rely on
/// `sim.pool` reflecting the change yet. Implementers that need to
/// inspect adjacent bases should hold their own state (built up via
/// `on_base_pushed` over the course of region assembly), not look
/// the new one up after the fact.
///
/// **Why `on_base_changed` carries the full `old_n: &Nucleotide`
/// rather than just `old_base: u8`.** The walker observer's
/// incremental score-delta needs the old nucleotide's `germline_pos`
/// to look up `compatible_alleles_at(ref_pos, old_base)`. It could
/// recover that from its own state — but passing it through the
/// event keeps the observer stateless w.r.t. position mapping and
/// keeps the trait honest about what's needed.
pub(crate) trait IrEventObserver {
    /// Notify the observer that one nucleotide is about to be
    /// committed to the pool at index `handle`.
    fn on_base_pushed(&mut self, handle: NucHandle, n: &Nucleotide);

    /// Notify the observer that the base byte at `handle` is about
    /// to change from `old_n.base` to `new_base`. The nucleotide's
    /// segment, germline_pos, and flags are unchanged — only the
    /// base byte is replaced.
    ///
    /// Default implementation: no-op. Observers that don't care
    /// about edits (e.g. the productive admit-mask, which only
    /// tracks the next-NP-slot counter) can ignore this event by
    /// not overriding the default.
    #[allow(unused_variables)]
    fn on_base_changed(&mut self, handle: NucHandle, old_n: &Nucleotide, new_base: u8) {}

    /// Notify the observer that a nucleotide is about to be
    /// **inserted** into the pool at position `at`. Every pool
    /// handle `>= at` shifts up by 1 as a result. `n` is the
    /// nucleotide being inserted (typically a synthetic
    /// indel-insert with `germline_pos == GermlinePos::NONE`).
    ///
    /// Default implementation: no-op. Observers that maintain
    /// handle-keyed state (e.g. the walker) inherit the no-op and
    /// become **stale** after an indel — the post-pass
    /// `PassEffect::StructuralIndel` dispatch in
    /// `compiled/execute.rs::apply_live_call_updates` triggers a
    /// full from-scratch refresh that overwrites the stale state.
    /// Indels are rare in production workloads (typical: 0-5 per
    /// record), so the from-scratch refresh is acceptable.
    #[allow(unused_variables)]
    fn on_indel_inserted(&mut self, at: u32, n: &Nucleotide) {}

    /// Notify the observer that the nucleotide at position `at` is
    /// about to be **deleted** from the pool. Every pool handle
    /// `> at` shifts down by 1; the handle at `at` is invalidated.
    /// `removed` is the nucleotide being removed.
    ///
    /// Default: no-op. Same staleness story as `on_indel_inserted`.
    #[allow(unused_variables)]
    fn on_indel_deleted(&mut self, at: u32, removed: &Nucleotide) {}
}

/// One slot in the [`SimulationBuilder`]'s observer registry.
///
/// Wraps each concrete observer kind in a typed enum variant so the
/// builder can hold a single `Vec<ObserverSlot<'idx>>` instead of
/// four parallel `Option<T>` / `Vec<T>` fields. Each mutation method
/// broadcasts through one loop; this enum's `IrEventObserver` impl
/// fans the event out to the wrapped observer's matching trait
/// method.
///
/// Adding a fifth observer kind = one new variant + one match arm
/// per event method (today: `on_base_pushed`, `on_base_changed`,
/// `on_indel_inserted`, `on_indel_deleted`). Mutation methods stay
/// untouched.
///
/// Conceptually the variants split into three **tiers**, mirroring
/// game-engine signal/slot ordering conventions (Bevy
/// `SystemSet::Recorder`/`Derived`/`Refresher`):
///
/// - **Derived** — observers that maintain incremental derived
///   state in reaction to events ([`Walker`](Self::Walker) tracks
///   per-allele scores; [`AdmitMask`](Self::AdmitMask) tracks the
///   next-NP-slot counter). Most common case.
/// - **Recorder** — passive observers that record the event stream
///   ([`EventLog`](Self::EventLog)). Today's only consumer is debug
///   replay and the MCP audit log.
/// - **Refresher** — observers that record what changed so a
///   downstream consumer can react ([`DirtySignal`](Self::DirtySignal)
///   accumulates [`DirtyWindow`]s for `apply_live_call_updates`).
///
/// Tier is not enforced in code today — observers fire in
/// `Vec`-insertion order. The trait contract guarantees that no
/// observer reads another's state during a broadcast, so the
/// ordering is purely a recording-order convention.
pub(crate) enum ObserverSlot<'idx> {
    /// Derived live-call walker observer. The builder may hold one
    /// (assembly path) or up to three (mutation-pass path attaching
    /// one observer per existing V/D/J region).
    Walker(WalkerObserverState<'idx>),
    /// Derived productive-admit-mask observer queried per NP slot
    /// via [`SimulationBuilder::current_admit_mask`]. No seal output
    /// — pure ephemeral query state.
    AdmitMask(ProductiveAdmitMaskObserver<'idx>),
    /// Recorder that captures every IR mutation event into a
    /// `Vec<IrEvent>`. Drained at seal time via
    /// [`SimulationBuilder::seal_event_log_observer`].
    EventLog(EventLogObserver),
    /// Refresher that accumulates `DirtyWindow`s for every base
    /// change / indel event. Drained inside
    /// [`SimulationBuilder::seal_with_committed_live_calls`].
    DirtySignal(DirtySignalObserver),
}

impl<'idx> IrEventObserver for ObserverSlot<'idx> {
    fn on_base_pushed(&mut self, handle: NucHandle, n: &Nucleotide) {
        match self {
            Self::Walker(o) => o.on_base_pushed(handle, n),
            Self::AdmitMask(o) => o.on_base_pushed(handle, n),
            Self::EventLog(o) => o.on_base_pushed(handle, n),
            Self::DirtySignal(o) => o.on_base_pushed(handle, n),
        }
    }

    fn on_base_changed(&mut self, handle: NucHandle, old_n: &Nucleotide, new_base: u8) {
        match self {
            Self::Walker(o) => o.on_base_changed(handle, old_n, new_base),
            Self::AdmitMask(o) => o.on_base_changed(handle, old_n, new_base),
            Self::EventLog(o) => o.on_base_changed(handle, old_n, new_base),
            Self::DirtySignal(o) => o.on_base_changed(handle, old_n, new_base),
        }
    }

    fn on_indel_inserted(&mut self, at: u32, n: &Nucleotide) {
        match self {
            Self::Walker(o) => o.on_indel_inserted(at, n),
            Self::AdmitMask(o) => o.on_indel_inserted(at, n),
            Self::EventLog(o) => o.on_indel_inserted(at, n),
            Self::DirtySignal(o) => o.on_indel_inserted(at, n),
        }
    }

    fn on_indel_deleted(&mut self, at: u32, removed: &Nucleotide) {
        match self {
            Self::Walker(o) => o.on_indel_deleted(at, removed),
            Self::AdmitMask(o) => o.on_indel_deleted(at, removed),
            Self::EventLog(o) => o.on_indel_deleted(at, removed),
            Self::DirtySignal(o) => o.on_indel_deleted(at, removed),
        }
    }
}

/// Mutable, by-value owner of an in-progress `Simulation` whose pool
/// supports cheap per-base append. See the module docs for the design
/// rationale.
///
/// The lifetime `'idx` is the lifetime of the borrowed
/// `SegmentRefIndex` used by an attached walker observer (if any).
/// Builders without an observer can safely use `'static`.
pub struct SimulationBuilder<'idx> {
    simulation: Simulation,
    /// Unified observer registry. Replaces the previous four parallel
    /// `Option<T>` / `Vec<T>` fields. Observers fire in insertion
    /// order during broadcasts. See [`ObserverSlot`] for the variant
    /// list and the tier convention.
    observers: Vec<ObserverSlot<'idx>>,
}

impl<'idx> SimulationBuilder<'idx> {
    /// Take ownership of `sim` and unique-ify its pool's inner `Vec`
    /// so subsequent `push_nucleotide` calls don't trigger a deep
    /// clone via `Arc::make_mut`.
    pub fn from_simulation(sim: Simulation) -> Self {
        let mut builder = Self {
            simulation: sim,
            observers: Vec::new(),
        };
        // Force unique ownership of the pool's nucleotide vector. After
        // this, every subsequent `Arc::make_mut` inside `pool.push`
        // is a refcount==1 fast path: it just hands back a
        // `&mut Vec<Nucleotide>` without copying.
        builder.simulation.pool.make_unique();
        builder
    }

    /// Broadcast `event_fn` to every attached observer in
    /// insertion order. Replaces the previous four-arm `if let
    /// Some(obs) = ...` fanout that each mutation method
    /// duplicated.
    #[inline]
    fn broadcast<F>(&mut self, mut event_fn: F)
    where
        F: FnMut(&mut ObserverSlot<'idx>),
    {
        for obs in &mut self.observers {
            event_fn(obs);
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
        self.observers.push(ObserverSlot::Walker(
            WalkerObserverState::new(segment_index, seq_start),
        ));
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
        self.observers.push(ObserverSlot::Walker(
            WalkerObserverState::from_existing_region(
                segment_index,
                &self.simulation,
                region,
            ),
        ));
    }

    /// Append one nucleotide and return its handle. Amortized O(1):
    /// no Arc clone, no `Vec` reallocation if capacity holds.
    ///
    /// Each attached observer is notified via
    /// [`IrEventObserver::on_base_pushed`] *before* the pool's `push`
    /// consumes `n`. The handle passed to observers equals the index
    /// `n` will occupy: `NucHandle::new(pool.len() as u32)` at the
    /// moment of call.
    pub fn push_nucleotide(&mut self, n: Nucleotide) -> NucHandle {
        let handle = NucHandle::new(self.simulation.pool.len() as u32);
        self.broadcast(|obs| obs.on_base_pushed(handle, &n));
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
    pub(crate) fn insert_indel(&mut self, at: u32, n: Nucleotide) {
        self.broadcast(|obs| obs.on_indel_inserted(at, &n));
        // The per-region codon rail is not maintained in the hot
        // path; consumers compute on demand. The persistent layer's
        // recompute would be wasted work.
        self.simulation = self.simulation.with_indel_inserted(at, n);
    }

    /// Delete the nucleotide at position `at`, shifting every existing
    /// handle `> at` down by 1. Notifies observers of `on_indel_deleted`
    /// before the underlying pool mutation applies the shift.
    /// Returns the removed nucleotide (caller may ignore).
    pub(crate) fn delete_indel(&mut self, at: u32) -> Option<Nucleotide> {
        let removed = *self.simulation.pool.get(NucHandle::new(at))?;
        self.broadcast(|obs| obs.on_indel_deleted(at, &removed));
        self.simulation = self.simulation.with_indel_deleted(at);
        Some(removed)
    }

    /// Attach a structured event-log observer that captures every
    /// IR mutation event into a `Vec<IrEvent>`.
    /// Drain via [`Self::seal_event_log_observer`] to obtain the
    /// captured stream. Used for trace replay, MCP audit, metrics.
    #[allow(dead_code)]
    pub(crate) fn attach_event_log_observer(&mut self) {
        self.observers
            .push(ObserverSlot::EventLog(EventLogObserver::new()));
    }

    /// Drain the attached event-log observer and return the
    /// captured events. Panics if no observer is attached.
    #[allow(dead_code)]
    pub(crate) fn seal_event_log_observer(&mut self) -> Vec<IrEvent> {
        let idx = self
            .observers
            .iter()
            .position(|o| matches!(o, ObserverSlot::EventLog(_)))
            .expect("seal_event_log_observer called without a prior attach_event_log_observer");
        match self.observers.remove(idx) {
            ObserverSlot::EventLog(o) => o.seal(),
            _ => unreachable!("position() located the EventLog variant"),
        }
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
    pub(crate) fn attach_dirty_signal_observer(&mut self) {
        self.observers
            .push(ObserverSlot::DirtySignal(DirtySignalObserver::new()));
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
    pub(crate) fn change_base(&mut self, handle: NucHandle, new_base: u8) -> Nucleotide {
        // Read the old nucleotide *before* mutating the pool — the
        // event observers need the old value with its segment and
        // germline_pos still populated.
        let old = *self
            .simulation
            .pool
            .get(handle)
            .expect("change_base: handle must point into the pool");

        self.broadcast(|obs| obs.on_base_changed(handle, &old, new_base));

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

    /// Find-and-remove every `ObserverSlot::Walker` from the
    /// registry, returning the wrapped `WalkerObserverState`s in
    /// attach order. Non-walker observers (admit-mask, event-log,
    /// dirty-signal) stay attached and continue receiving events
    /// until their own seal / drop. Used internally by
    /// [`Self::seal_walker_observer`] / [`Self::seal_all_walker_observers`].
    fn take_walker_observers(&mut self) -> Vec<WalkerObserverState<'idx>> {
        let mut walkers = Vec::new();
        let mut keep = Vec::with_capacity(self.observers.len());
        for slot in std::mem::take(&mut self.observers) {
            match slot {
                ObserverSlot::Walker(w) => walkers.push(w),
                other => keep.push(other),
            }
        }
        self.observers = keep;
        walkers
    }

    /// Find-and-remove the attached `ObserverSlot::DirtySignal` (if
    /// any) and return the wrapped observer. Other observers stay
    /// in the registry. Used internally by
    /// [`Self::seal_with_committed_live_calls`] to drain dirty
    /// windows.
    fn take_dirty_signal_observer(&mut self) -> Option<DirtySignalObserver> {
        let idx = self
            .observers
            .iter()
            .position(|o| matches!(o, ObserverSlot::DirtySignal(_)))?;
        match self.observers.remove(idx) {
            ObserverSlot::DirtySignal(o) => Some(o),
            _ => unreachable!("position() located the DirtySignal variant"),
        }
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
    /// `PassRuntime::execute_with_refdata`), this method is a no-op.
    ///
    /// The per-region codon rail (`amino_acids` /
    /// `stop_codon_positions`) is intentionally NOT maintained: no
    /// production consumer reads it (`NoStopCodonInJunction` reads
    /// pool directly, AIRR re-translates from raw sequence). Callers
    /// that need the rail compute it on demand via
    /// `Region::with_codon_rail_recomputed`.
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
        let dirty_windows = self.take_dirty_signal_observer().map(|o| o.seal()).unwrap_or_default();
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
            let position = match walker_sealed.iter().position(|s| s.segment() == *target_segment) {
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
    pub(crate) fn attach_admit_mask_observer(
        &mut self,
        state: &'idx JunctionStopState,
        np_segment: Segment,
    ) {
        self.observers.push(ObserverSlot::AdmitMask(
            ProductiveAdmitMaskObserver::new(state, np_segment),
        ));
    }

    /// 4-bit admit mask for the current NP slot, or `None` when no
    /// admit-mask observer is attached. Bit layout: bit 0 → `A`,
    /// bit 1 → `C`, bit 2 → `G`, bit 3 → `T`. See
    /// [`ProductiveAdmitMaskObserver::current_admit_mask`].
    pub(crate) fn current_admit_mask(&self) -> Option<u8> {
        self.observers.iter().find_map(|o| match o {
            ObserverSlot::AdmitMask(am) => Some(am.current_admit_mask(&self.simulation)),
            _ => None,
        })
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
    let Some(region) = sim
        .sequence
        .regions
        .iter()
        .find(|r| r.segment == segment)
    else {
        return false;
    };
    let region_start = region.start.index();
    let region_end = region.end.index();
    windows
        .iter()
        .any(|w| w.start < region_end && w.end > region_start)
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
}
