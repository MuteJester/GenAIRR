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

use super::codon_rail_observer::{CodonRailObserverState, SealedCodonRail};
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
/// Still to come (Phase 5+): `on_indel_inserted`, `on_indel_deleted`,
/// `on_region_added`, `on_assignment_made`. Each will be additive —
/// existing implementers inherit no-op default methods and keep
/// working.
///
/// **Visibility.** `pub(crate)` deliberately: every observer lives in
/// the same crate (walker, codon rail, productive admit-mask, future
/// telemetry). External crates are not expected to implement this;
/// if they ever do, the trait can be widened to `pub` then.
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
    /// handle-keyed state (walker, codon rail) inherit the
    /// no-op and become **stale** after an indel — the post-pass
    /// `PassEffect::StructuralIndel` dispatch in
    /// `compiled/execute.rs::apply_live_call_updates` triggers a
    /// full from-scratch refresh that overwrites the stale state.
    ///
    /// A future phase can replace the no-op default with proper
    /// handle-shift propagation in observers, suppressing the
    /// post-pass refresh — but indels are rare in production
    /// workloads (typical: 0-5 per record), so the from-scratch
    /// refresh is acceptable for now.
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

// ── FUTURE TRAIT SHAPE ───────────────────────────────────────────────────────
//
// Phase 1 wires *one* concrete observer (the live-call walker) into
// `SimulationBuilder` via a single optional slot. The shape that *will*
// generalise in Phase 1.5+ is sketched here so the trait extraction
// has the design data in front of it:
//
//   trait BaseAppendObserver {
//       fn on_base_pushed(&mut self, nuc: &Nucleotide);
//   }
//
// Prospective consumers and the union of fields each will need from
// the event payload:
//
//   * Walker (this phase, `WalkerObserverState`): needs
//     `nuc.segment`, `nuc.base`, `nuc.germline_pos`. Accumulates a
//     per-allele score vector and informative/wildcard tallies. See
//     `live_call/walker_observer.rs`.
//
//   * Codon rail (Phase 2): needs `nuc.base` plus the post-push
//     `seq_pos` (= pool.len()-1 at observation time). Maintains a
//     per-region `amino_acids: Vec<u8>` + `stop_codon_positions:
//     Vec<NucHandle>` incrementally — see `ir/region.rs:64-100` for
//     the current from-scratch rebuild and `ir/simulation.rs:203-210`
//     for the per-handle refresh pattern that drives it today.
//
//   * Productive-contract admit-mask (Phase 3): needs `seq_pos`,
//     `nuc.base`, the frame phase of the upcoming-codon slot, and a
//     junction-relative position. Maintains a 4-bit mask of
//     "admissible next bases" per upcoming NP slot. Replaces the
//     per-candidate `NoStopCodonInJunction::admits_with_context`
//     query with an O(1) mask read.
//
//   * Telemetry / trace (Phase 5+): needs all of the above plus the
//     event ordering. Likely fans out via a `DynObserverSet`
//     adapter implementing the same trait.
//
// Two events the current API doesn't yet emit but that downstream
// observers will need:
//
//   * `on_base_changed(handle, old_base, new_base)` — fired from
//     `Simulation::with_base_changed`; load-bearing for S5F/PCR/quality
//     observers in Phase 4.
//   * `on_indel_inserted` / `on_indel_deleted` — Phase 4 again, when
//     `with_indel_*` start emitting and observers stop relying on the
//     coarse `PassEffect::StructuralIndel` post-pass refresh.
//
// We do NOT publish the trait in Phase 1 because we only have one
// concrete implementer; designing in the abstract with only one data
// point is speculative. The walker ships as a private slot; the trait
// extraction happens in Phase 1.5 once the codon rail port begins.

/// Mutable, by-value owner of an in-progress `Simulation` whose pool
/// supports cheap per-base append. See the module docs for the design
/// rationale.
///
/// The lifetime `'idx` is the lifetime of the borrowed
/// `SegmentRefIndex` used by an attached walker observer (if any).
/// Builders without an observer can safely use `'static` (the borrow
/// slot stays `None`).
pub struct SimulationBuilder<'idx> {
    simulation: Simulation,
    /// Streaming live-call walker observers, one per segment under
    /// observation. During assembly there is exactly one (the V/D/J
    /// segment currently being built). During post-assembly mutation
    /// passes (S5F, PCR, …) there may be up to three — one per
    /// existing V/D/J region so each can receive `on_base_changed`
    /// events for bases in its own range. Events broadcast to every
    /// attached observer; observers self-filter by segment.
    walker_observers: Vec<WalkerObserverState<'idx>>,
    /// Streaming codon-rail observers, one per region under
    /// observation. During assembly there is exactly one (the
    /// region currently being built). During post-assembly
    /// mutation passes (S5F, PCR, …) there may be several — one
    /// per existing V/D/J region so each can receive
    /// `on_base_changed` events for bases in its own range.
    /// Each event is broadcast to every observer; observers
    /// self-filter by segment in `on_base_pushed` and by both
    /// segment and pool-range in `on_base_changed`.
    codon_rail_observers: Vec<CodonRailObserverState>,
    /// Optional streaming productive-admit-mask observer. Populated
    /// by `attach_admit_mask_observer` and queried per NP slot via
    /// `current_admit_mask`; advanced lazily on each
    /// `push_nucleotide` whose segment matches the observer's NP
    /// segment. No `seal_*` method — the observer carries no
    /// committable output, only ephemeral query state.
    admit_mask_observer: Option<ProductiveAdmitMaskObserver<'idx>>,
    /// Phase 13: optional structured event-log observer. Captures
    /// every IR mutation event into a `Vec<IrEvent>` for trace
    /// replay / MCP audit / metrics. Pure recorder; doesn't react.
    event_log_observer: Option<EventLogObserver>,
    /// Phase 15: optional dirty-signal observer. Records a
    /// `DirtyWindow` for every base change / indel event so the
    /// post-pass `apply_live_call_updates` dispatch can refresh only
    /// segments whose region overlaps a dirty position instead of
    /// re-sweeping V/D/J unconditionally.
    dirty_signal_observer: Option<DirtySignalObserver>,
}

impl<'idx> SimulationBuilder<'idx> {
    /// Take ownership of `sim` and unique-ify its pool's inner `Vec`
    /// so subsequent `push_nucleotide` calls don't trigger a deep
    /// clone via `Arc::make_mut`.
    pub fn from_simulation(sim: Simulation) -> Self {
        let mut builder = Self {
            simulation: sim,
            walker_observers: Vec::new(),
            codon_rail_observers: Vec::new(),
            admit_mask_observer: None,
            event_log_observer: None,
            dirty_signal_observer: None,
        };
        // Force unique ownership of the pool's nucleotide vector. After
        // this, every subsequent `Arc::make_mut` inside `pool.push`
        // is a refcount==1 fast path: it just hands back a
        // `&mut Vec<Nucleotide>` without copying.
        builder.simulation.pool.make_unique();
        builder
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
        self.walker_observers
            .push(WalkerObserverState::new(segment_index, seq_start));
    }

    /// Phase 9: attach a walker observer **rebuilt** from an
    /// already-assembled region.
    ///
    /// Walks the region's bytes in the current simulation through
    /// scoring (cost: O(region_len), same as a from-scratch
    /// `call_from_region` for that region) and installs the
    /// resulting state. Subsequent `change_base` events update the
    /// scores incrementally via `on_base_changed`.
    ///
    /// Companion to
    /// [`Self::attach_codon_rail_observer_for_region`]: together
    /// they let post-assembly mutation passes (S5F, PCR, …) attach
    /// the full observer set to a `SimulationBuilder` so the
    /// post-pass `PassEffect::EditBases` refresh becomes redundant.
    /// The actual suppression / version-bumping coordination is
    /// done by the mutation pass at seal time; this helper just
    /// installs one observer.
    pub(crate) fn attach_walker_observer_for_region(
        &mut self,
        segment_index: &'idx SegmentRefIndex,
        region: &super::Region,
    ) {
        self.walker_observers
            .push(WalkerObserverState::from_existing_region(
                segment_index,
                &self.simulation,
                region,
            ));
    }

    /// Attach a streaming codon-rail observer for one region's
    /// assembly. `seq_start` is the pool index of the first base of
    /// the soon-to-be-assembled region (= `peek().pool.len()` at
    /// attach time). `frame_phase` is the region's frame phase
    /// (chained from prior regions' cumulative-length-mod-3).
    ///
    /// Only one codon-rail observer may be attached at a time;
    /// calling this while another is active replaces it (the old
    /// state is dropped without being sealed, matching the
    /// fixture-friendly "no observer" path).
    pub(crate) fn attach_codon_rail_observer(
        &mut self,
        segment: Segment,
        seq_start: u32,
        frame_phase: u8,
    ) {
        self.codon_rail_observers
            .push(CodonRailObserverState::new(segment, seq_start, frame_phase));
    }

    /// Phase 6: attach a codon-rail observer reconstructed from an
    /// already-assembled region. Used by post-assembly mutation
    /// passes (S5F, PCR, …) that want to receive incremental
    /// `on_base_changed` updates without re-running assembly.
    ///
    /// The observer is initialised with the region's existing
    /// `amino_acids` and `stop_codon_positions` plus the current
    /// byte values in `[region.start, region.end)`. Subsequent
    /// `change_base` calls update the rail incrementally; sealing
    /// returns the updated `(amino_acids, stop_codon_positions)`.
    pub(crate) fn attach_codon_rail_observer_for_region(&mut self, region: &super::Region) {
        let start = region.start.index() as usize;
        let end = region.end.index() as usize;
        let slice: Vec<u8> = self
            .simulation
            .pool
            .as_slice()
            .get(start..end)
            .expect("attach_codon_rail_observer_for_region: region range outside pool")
            .iter()
            .map(|n| n.base)
            .collect();
        self.codon_rail_observers
            .push(CodonRailObserverState::from_existing_region(region, &slice));
    }

    /// Append one nucleotide and return its handle. Amortized O(1):
    /// no Arc clone, no `Vec` reallocation if capacity holds.
    ///
    /// Each attached observer is notified via
    /// [`BaseAppendObserver::on_base_pushed`] *before* the pool's
    /// `push` consumes `n`. The handle passed to observers equals the
    /// index `n` will occupy: `NucHandle::new(pool.len() as u32)` at
    /// the moment of call.
    pub fn push_nucleotide(&mut self, n: Nucleotide) -> NucHandle {
        let handle = NucHandle::new(self.simulation.pool.len() as u32);
        for obs in &mut self.walker_observers {
            obs.on_base_pushed(handle, &n);
        }
        for obs in &mut self.codon_rail_observers {
            obs.on_base_pushed(handle, &n);
        }
        if let Some(obs) = &mut self.admit_mask_observer {
            obs.on_base_pushed(handle, &n);
        }
        if let Some(obs) = &mut self.event_log_observer {
            obs.on_base_pushed(handle, &n);
        }
        if let Some(obs) = &mut self.dirty_signal_observer {
            obs.on_base_pushed(handle, &n);
        }
        // The pool's `push` already does `Arc::make_mut` internally;
        // after `from_simulation`'s prepayment, that call is the cheap
        // refcount==1 path on every iteration. The returned handle
        // equals the one we computed above (Pool::push uses pool.len()
        // before the push to mint the handle, same as we did).
        self.simulation.pool.push(n)
    }

    /// Phase 12: insert `n` into the pool at position `at`, shifting
    /// every existing handle `>= at` up by 1. Notifies every
    /// attached observer via `on_indel_inserted` before the
    /// underlying `Simulation::with_indel_inserted` applies the
    /// shift.
    ///
    /// Indels are structurally different from base-changes: they
    /// invalidate handle-keyed observer state. The walker and
    /// codon-rail observers inherit the no-op default impl of
    /// `on_indel_inserted`, so their state becomes stale —
    /// downstream consumers must fall back to from-scratch
    /// refresh, which the post-pass `PassEffect::StructuralIndel`
    /// dispatch in `compiled/execute.rs` already does.
    ///
    /// Indels are rare (typical workloads have 0-5 per record), so
    /// this delegates to the persistent `Simulation::with_indel_*`
    /// API rather than mutating the inner pool/sequence in place.
    /// The Arc-clones inside `with_indel_inserted` are negligible
    /// at this frequency. A future phase can add in-place handling
    /// (and proper observer handle-shift propagation) if profiling
    /// shows indels need it.
    #[allow(dead_code)]
    pub(crate) fn insert_indel(&mut self, at: u32, n: Nucleotide) {
        for obs in &mut self.walker_observers {
            obs.on_indel_inserted(at, &n);
        }
        for obs in &mut self.codon_rail_observers {
            obs.on_indel_inserted(at, &n);
        }
        if let Some(obs) = &mut self.admit_mask_observer {
            obs.on_indel_inserted(at, &n);
        }
        if let Some(obs) = &mut self.event_log_observer {
            obs.on_indel_inserted(at, &n);
        }
        if let Some(obs) = &mut self.dirty_signal_observer {
            obs.on_indel_inserted(at, &n);
        }
        self.simulation = self.simulation.with_indel_inserted(at, n);
    }

    /// Phase 12: delete the nucleotide at position `at`, shifting
    /// every existing handle `> at` down by 1. Notifies observers
    /// of `on_indel_deleted` before the underlying
    /// `Simulation::with_indel_deleted` applies the shift.
    /// Returns the removed nucleotide (caller may ignore).
    #[allow(dead_code)]
    pub(crate) fn delete_indel(&mut self, at: u32) -> Option<Nucleotide> {
        let removed = *self.simulation.pool.get(NucHandle::new(at))?;
        for obs in &mut self.walker_observers {
            obs.on_indel_deleted(at, &removed);
        }
        for obs in &mut self.codon_rail_observers {
            obs.on_indel_deleted(at, &removed);
        }
        if let Some(obs) = &mut self.admit_mask_observer {
            obs.on_indel_deleted(at, &removed);
        }
        if let Some(obs) = &mut self.event_log_observer {
            obs.on_indel_deleted(at, &removed);
        }
        if let Some(obs) = &mut self.dirty_signal_observer {
            obs.on_indel_deleted(at, &removed);
        }
        self.simulation = self.simulation.with_indel_deleted(at);
        Some(removed)
    }

    /// Phase 13: attach a structured event-log observer that
    /// captures every IR mutation event into a `Vec<IrEvent>`.
    /// Drain via [`Self::seal_event_log_observer`] to obtain the
    /// captured stream. Used for trace replay, MCP audit, metrics.
    #[allow(dead_code)]
    pub(crate) fn attach_event_log_observer(&mut self) {
        self.event_log_observer = Some(EventLogObserver::new());
    }

    /// Drain the attached event-log observer and return the
    /// captured events. Panics if no observer is attached.
    #[allow(dead_code)]
    pub(crate) fn seal_event_log_observer(&mut self) -> Vec<IrEvent> {
        let observer = self
            .event_log_observer
            .take()
            .expect("seal_event_log_observer called without a prior attach_event_log_observer");
        observer.seal()
    }

    /// Phase 15: attach a dirty-signal observer that records a
    /// `DirtyWindow` for every base change / indel event. The post-
    /// pass `apply_live_call_updates` dispatch reads the captured
    /// windows to refresh only segments whose region overlaps a
    /// dirty position, replacing the prior unconditional V/D/J sweep.
    ///
    /// Mutation passes (S5F, PCR, quality, contaminant, ncorrupt,
    /// uniform-mutation, indel) opt in to this observer alongside
    /// their walker and codon-rail observers. Pure-assembly paths
    /// don't attach it: assembly is dispatched via
    /// `PassEffect::AssembleSegment` / `AppendRegion`, which have
    /// their own refresh logic in `apply_live_call_updates`.
    pub(crate) fn attach_dirty_signal_observer(&mut self) {
        self.dirty_signal_observer = Some(DirtySignalObserver::new());
    }

    /// Phase 4: change the base byte at `handle` in place and notify
    /// every attached observer of an `on_base_changed` event.
    ///
    /// Returns the `Nucleotide` as it was *before* the change.
    ///
    /// Companion to [`Self::push_nucleotide`]: observers see the
    /// old nucleotide via the event parameter (so they can compute
    /// score deltas, retranslate codons, etc.) *before* the pool
    /// reflects the new value. The new base is applied immediately
    /// after the notifications return.
    ///
    /// Currently unused by production passes — mutation passes
    /// (S5F, PCR, quality, contaminant) still operate on the
    /// persistent `Simulation::with_base_changed` API and route
    /// through the post-pass `PassEffect::EditBases` coarse
    /// refresh. Phase 5+ will port them to this builder method so
    /// the walker / codon-rail observers can update incrementally
    /// instead.
    #[allow(dead_code)]
    pub(crate) fn change_base(&mut self, handle: NucHandle, new_base: u8) -> Nucleotide {
        // Read the old nucleotide *before* mutating the pool — the
        // event observers need the old value with its segment and
        // germline_pos still populated.
        let old = *self
            .simulation
            .pool
            .get(handle)
            .expect("change_base: handle must point into the pool");

        for obs in &mut self.walker_observers {
            obs.on_base_changed(handle, &old, new_base);
        }
        for obs in &mut self.codon_rail_observers {
            obs.on_base_changed(handle, &old, new_base);
        }
        if let Some(obs) = &mut self.admit_mask_observer {
            obs.on_base_changed(handle, &old, new_base);
        }
        if let Some(obs) = &mut self.event_log_observer {
            obs.on_base_changed(handle, &old, new_base);
        }
        if let Some(obs) = &mut self.dirty_signal_observer {
            obs.on_base_changed(handle, &old, new_base);
        }

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
        assert!(
            self.walker_observers.len() == 1,
            "seal_walker_observer expects exactly one attached observer; \
             got {}. Use seal_all_walker_observers for multi-region passes.",
            self.walker_observers.len()
        );
        self.walker_observers.pop().unwrap().seal(seq_end)
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
    /// Phase 16: observers whose `needs_rebuild` flag is set (by an
    /// internal indel event) are rebuilt via `from_existing_region`
    /// against the current (post-mutation) simulation *before*
    /// sealing, so the sealed state reflects the post-indel pool
    /// rather than stale pre-indel scores.
    pub(crate) fn seal_all_walker_observers(&mut self) -> Vec<SealedWalkerState> {
        let sim = &self.simulation;
        std::mem::take(&mut self.walker_observers)
            .into_iter()
            .map(|obs| obs.rebuild_if_stale(sim))
            .map(|obs| {
                let seq_end = obs.seq_end_hint();
                obs.seal(seq_end)
            })
            .collect()
    }

    /// Drain the attached codon-rail observer and return its sealed
    /// state. Returns `(amino_acids, stop_codon_positions)` that the
    /// caller can write directly into a new `Region`, skipping the
    /// post-pass `with_codon_rail_recomputed` rebuild.
    ///
    /// Panics if there isn't exactly one observer attached — call
    /// [`Self::seal_all_codon_rail_observers`] for the multi-region
    /// case (mutation passes attaching one observer per existing
    /// V/D/J region).
    pub(crate) fn seal_codon_rail_observer(&mut self) -> SealedCodonRail {
        assert!(
            self.codon_rail_observers.len() == 1,
            "seal_codon_rail_observer expects exactly one attached observer; \
             got {}. Use seal_all_codon_rail_observers for multi-region passes.",
            self.codon_rail_observers.len()
        );
        self.codon_rail_observers.pop().unwrap().seal()
    }

    /// Drain every attached codon-rail observer, returning one sealed
    /// rail per attachment in attach order. Used by post-assembly
    /// mutation passes (S5F, PCR, …) that attached one observer per
    /// existing V/D/J region.
    ///
    /// Returns each sealed rail paired with the segment + region
    /// start handle the observer was attached to, so the caller can
    /// write the rail back into the matching `Region` in the
    /// `Simulation.sequence.regions` vec.
    pub(crate) fn seal_all_codon_rail_observers(&mut self) -> Vec<(Segment, NucHandle, SealedCodonRail)> {
        let sim = &self.simulation;
        std::mem::take(&mut self.codon_rail_observers)
            .into_iter()
            .map(|obs| obs.rebuild_if_stale(sim))
            .map(|obs| {
                let segment = obs.observed_segment();
                let start = obs.observed_start_handle();
                (segment, start, obs.seal())
            })
            .collect()
    }

    /// Convenience: attach a codon-rail observer for every existing
    /// region in the current simulation. Used by post-assembly
    /// mutation passes (S5F, PCR, quality, contaminant, ncorrupt,
    /// uniform-mutation) that need to receive `on_base_changed`
    /// events for bytes in *any* region — observers self-filter
    /// by their own segment.
    pub(crate) fn attach_codon_rail_observers_for_all_regions(&mut self) {
        // Clone the regions snapshot because the subsequent attach
        // method borrows `&self.simulation.pool` and we'll then take
        // `&mut self` to install observers.
        let regions: Vec<super::Region> =
            self.simulation.sequence.regions.iter().cloned().collect();
        for region in &regions {
            self.attach_codon_rail_observer_for_region(region);
        }
    }

    /// Phase 22: standard observer-attach pattern for the mutation
    /// passes (S5F, uniform, PCR, quality, ncorrupt, contaminant,
    /// indel, end-loss). Always attaches a codon-rail observer for
    /// every existing region. When a `ReferenceMatchIndex` is
    /// available, additionally attaches the dirty-signal observer
    /// (consumed by `apply_live_call_updates` in the compiled path)
    /// and a walker observer for every V/D/J region with a matching
    /// segment index.
    ///
    /// Without a `ReferenceMatchIndex`, only the codon-rail observer
    /// attaches. Test harnesses like `PassRuntime::execute_with_refdata`
    /// take this branch — they don't run `apply_live_call_updates`,
    /// so dirty-signal + walker observers would be dead state.
    /// Phase 22 collapsed the previous "attach dirty signal observer
    /// for symmetry" branch in `indel.rs` / `end_loss.rs` because
    /// nothing consumed it.
    pub(crate) fn attach_standard_mutation_observers(
        &mut self,
        reference_index: Option<&'idx crate::live_call::ReferenceMatchIndex>,
    ) {
        self.attach_codon_rail_observers_for_all_regions();
        let Some(ref_index) = reference_index else {
            return;
        };
        self.attach_dirty_signal_observer();
        // Clone the regions snapshot because the subsequent attach
        // method borrows `&self.simulation` and we then take
        // `&mut self` to install observers.
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

    /// Phase 10: seal **both** codon-rail observers and walker
    /// observers, write back codon rails to their regions, finalize
    /// each walker observer into a `SegmentLiveCall`, and stage the
    /// resulting live calls on the simulation with **monotonically
    /// increasing evidence_versions in V→D→J order**.
    ///
    /// The version sequence matches the order in which
    /// `compiled/execute.rs::apply_live_call_updates` processes
    /// `PassEffect::EditBases` (V first, then D, then J). With each
    /// staged call's `evidence_version` set to `base.version + 1`,
    /// `base.version + 2`, `base.version + 3` respectively, all three
    /// Phase-1 fast paths fire in the post-pass refresh: each
    /// segment absorbs its staged call via
    /// [`with_segment_call`](crate::live_call::LiveCallState::with_segment_call)
    /// (which bumps the version) and the from-scratch `call_from_region`
    /// rebuild is skipped entirely.
    ///
    /// Mutation passes that attach both observer types call this
    /// instead of `seal_with_committed_codon_rails`. Passes that only
    /// attach codon-rail observers (the Phase 8 set) continue to use
    /// the simpler helper.
    pub(crate) fn seal_with_committed_live_calls(
        mut self,
        ref_index: &crate::live_call::ReferenceMatchIndex,
    ) -> Simulation {
        // Drain observer kinds while we still have the builder.
        let codon_rail_sealed = self.seal_all_codon_rail_observers();
        let walker_sealed = self.seal_all_walker_observers();
        // Phase 15: drain dirty signals captured during the pass so
        // `apply_live_call_updates` can use them to skip refreshes
        // for segments that received no edits.
        let dirty_windows = self
            .dirty_signal_observer
            .take()
            .map(|obs| obs.seal())
            .unwrap_or_default();
        let mut current = self.seal();

        // Write back codon rails (same as the codon-rail-only helper).
        for (segment, start, rail) in codon_rail_sealed {
            if let Some(existing) = current
                .sequence
                .regions
                .iter()
                .find(|r| r.segment == segment && r.start == start)
                .cloned()
            {
                current = current.with_region_replaced_for_segment(super::Region::from_sealed_codon_rail(
                    existing.segment,
                    existing.start,
                    existing.end,
                    existing.frame_phase,
                    rail,
                ));
            }
        }

        // Phase 23: take ownership of the LiveCallState once,
        // mutate in place across V→D→J staging + dirty-window stash,
        // and reattach via a single `with_live_calls` at the end.
        // This eliminates the 7-clone churn the previous
        // clone-per-segment pattern produced.
        //
        // Phase 15: stage walker calls ONLY for segments whose region
        // overlaps a dirty window. Clean segments inherit their prior
        // live call unchanged, so the consumer
        // (`apply_live_call_updates`) can safely skip them — this is
        // what makes the dirty-window optimisation correctness-
        // preserving without breaking the V→D→J monotonic version
        // chain that the Phase-1 fast path relies on.
        //
        // When no observer was attached (test paths feeding an empty
        // `dirty_windows`), every segment falls through to staging so
        // the existing version-chain behaviour is preserved.
        let mut state: crate::live_call::LiveCallState = current
            .live_calls
            .as_ref()
            .map(|s| (**s).clone())
            .unwrap_or_default();
        let base_version = state.version;
        let has_dirty_observer_signal = !dirty_windows.is_empty();
        let order = [
            crate::ir::Segment::V,
            crate::ir::Segment::D,
            crate::ir::Segment::J,
        ];
        let mut walker_sealed = walker_sealed;
        let mut staged_count: u64 = 0;
        for target_segment in order.iter() {
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

            state.stage_segment_call(live_call);
        }

        // Stash dirty windows on the same state instance.
        state.dirty_windows.extend(dirty_windows);

        current.with_live_calls(state)
    }

    /// Convenience: seal every attached codon-rail observer and
    /// commit each sealed rail back into the matching region in the
    /// final `Simulation`. Consumes the builder.
    ///
    /// For each observer, looks up the region in
    /// `sim.sequence.regions` by `(segment, start)` (which is stable
    /// across mutation passes that don't shift handles) and
    /// replaces it with one whose `amino_acids` /
    /// `stop_codon_positions` come from the observer. If no
    /// matching region is found (e.g. an indel pass shifted things
    /// — out of scope for Phase 7-9), the rail is silently dropped
    /// and the region keeps its pre-mutation rail.
    pub(crate) fn seal_with_committed_codon_rails(mut self) -> Simulation {
        let sealed_rails = self.seal_all_codon_rail_observers();
        let mut current = self.seal();
        for (segment, start, rail) in sealed_rails {
            if let Some(existing) = current
                .sequence
                .regions
                .iter()
                .find(|r| r.segment == segment && r.start == start)
                .cloned()
            {
                current = current.with_region_replaced_for_segment(super::Region::from_sealed_codon_rail(
                    existing.segment,
                    existing.start,
                    existing.end,
                    existing.frame_phase,
                    rail,
                ));
            }
        }
        current
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
        self.admit_mask_observer = Some(ProductiveAdmitMaskObserver::new(state, np_segment));
    }

    /// 4-bit admit mask for the current NP slot, or `None` when no
    /// admit-mask observer is attached. Bit layout: bit 0 → `A`,
    /// bit 1 → `C`, bit 2 → `G`, bit 3 → `T`. See
    /// [`ProductiveAdmitMaskObserver::current_admit_mask`].
    pub(crate) fn current_admit_mask(&self) -> Option<u8> {
        self.admit_mask_observer
            .as_ref()
            .map(|obs| obs.current_admit_mask(&self.simulation))
    }

    /// Consume the builder and return the sealed `Simulation`.
    ///
    /// At seal time the inner pool still holds an `Arc<Vec<Nucleotide>>`
    /// with refcount==1; subsequent persistent operations on the sealed
    /// `Simulation` (e.g. `with_region_added`) clone the outer pool
    /// struct (cheap) without forcing a deep clone of the vector.
    ///
    /// If a walker observer is still attached, it is dropped without
    /// being sealed — this matches the pre-Phase-1 "no observer" path.
    pub fn seal(self) -> Simulation {
        self.simulation
    }
}

/// Phase 15 helper: does the segment's assembled region (if any)
/// overlap any of the given dirty windows? Used by
/// `seal_with_committed_live_calls` to gate which walker observers
/// stage their calls.
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
