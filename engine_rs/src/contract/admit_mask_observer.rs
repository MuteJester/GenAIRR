//! Streaming productive-admit-mask observer attached to
//! `SimulationBuilder`.
//!
//! `NoStopCodonInJunction` enforces "no stop codon in the junction"
//! by being **queried per candidate**: for each NP-base draw, the
//! sampler asks the contract four times (once per `{A,C,G,T}`
//! candidate). Each query is O(1) against the precomputed
//! [`JunctionStopState`] — but the per-candidate dispatch traffic,
//! plus `sample_filtered_result`'s temporary `Vec<(T, f64)>`
//! materialisation around it, was previously a dominant hot path.
//!
//! `ProductiveAdmitMaskObserver` rides the
//! [`SimulationEventSink`](crate::ir::SimulationEventSink) channel
//! — the same one
//! [`crate::ir::event_log_observer::EventLogObserver`],
//! [`crate::live_call::dirty_signal_observer::DirtySignalObserver`],
//! and [`crate::live_call::walker_observer::WalkerObserverState`]
//! migrated onto. It tracks one piece of derived state — the
//! next NP-slot index — and uses it to compute a per-slot
//! admissibility cache (rather than per-allele scoring or
//! per-position translation).
//!
//! ## Reactive state model
//!
//! The observer borrows the precomputed [`JunctionStopState`] and
//! tracks one integer — the next NP-slot index to be drawn. Every
//! time the NP pass commits a base belonging to the observer's
//! `np_segment`, the observer increments its counter. The 4-bit
//! admit mask for the *current* slot is computed lazily, on
//! demand, via [`Self::current_admit_mask`] — at which point the
//! simulation pool already reflects every previously-committed NP
//! byte, so the mask reads consistent state.
//!
//! Bit layout of the mask:
//! - bit 0 → `A`
//! - bit 1 → `C`
//! - bit 2 → `G`
//! - bit 3 → `T`
//!
//! ## Which events drive state
//!
//! - [`SimulationEvent::BasePushed`] with `segment == np_segment` →
//!   bump `current_np_index` by 1. This is the only state-changing
//!   event.
//! - Every other variant (`BasePushed` of a different segment,
//!   `BaseChanged`, `IndelInserted`, `IndelDeleted`, and the
//!   reserved future-emission variants): no-op. NP-region sampling
//!   doesn't run alongside mutation passes, so we never expect to
//!   see those during a real attach window, but the sink contract
//!   is to silently match-and-ignore rather than panic.
//!
//! ## Where this fires
//!
//! Attached in `GenerateNPPass::execute_with_sampling_mode` right
//! after `JunctionStopState::build`. Per NP base sample the
//! sampler calls `current_admit_mask(builder.peek())` once,
//! routes around `sample_filtered_result`'s generic predicate
//! loop, and inverse-CDFs over the admitted subset of the base
//! distribution's support directly (one Vec allocation per draw
//! instead of two; one mask-lookup per candidate instead of one
//! full contract dispatch).

use crate::ir::{Segment, Simulation, SimulationEvent, SimulationEventSink};

use super::junction_stop_state::JunctionStopState;

/// Streaming observer that caches "which of {A,C,G,T} admit at the
/// current NP slot" via the precomputed [`JunctionStopState`].
///
/// One observer instance covers one NP region's sampling. The
/// `SimulationBuilder` holds it as
/// `Option<ProductiveAdmitMaskObserver<'state>>`, where `'state` is
/// the lifetime of the borrowed `JunctionStopState`.
pub(crate) struct ProductiveAdmitMaskObserver<'state> {
    /// Borrowed precomputed contract state. Built once per record by
    /// `JunctionStopState::build` and shared across the entire NP
    /// region's sampling.
    state: &'state JunctionStopState,
    /// Which NP segment this observer is tracking — either
    /// `Segment::Np1` or `Segment::Np2`. Used to decide whether an
    /// inbound `BasePushed` event belongs to *our* NP region or to
    /// (say) an adjacent assembly that doesn't advance our slot
    /// counter.
    np_segment: Segment,
    /// Index of the next NP slot to be sampled. Starts at 0 and
    /// increments by 1 on every `BasePushed` whose nucleotide is in
    /// our `np_segment`.
    current_np_index: u32,
}

impl<'state> ProductiveAdmitMaskObserver<'state> {
    /// Initialize a fresh observer at NP-slot index 0.
    pub(crate) fn new(state: &'state JunctionStopState, np_segment: Segment) -> Self {
        Self {
            state,
            np_segment,
            current_np_index: 0,
        }
    }

    /// Compute the 4-bit admit mask for the current NP slot.
    ///
    /// `sim` must be the current in-progress simulation as observed
    /// by `builder.peek()` — its pool reflects every previously-
    /// committed NP byte. The mask is computed by asking the
    /// underlying [`JunctionStopState`] whether each canonical base
    /// admits at the current slot.
    ///
    /// Bit layout: `(mask >> i) & 1 == 1` iff base `BASES[i]` admits,
    /// where `BASES = [b'A', b'C', b'G', b'T']`.
    pub(crate) fn current_admit_mask(&self, sim: &Simulation) -> u8 {
        let mut mask = 0u8;
        for (i, &base) in BASES.iter().enumerate() {
            if self
                .state
                .admits_np_candidate(sim, self.np_segment, self.current_np_index, base)
                .is_ok()
            {
                mask |= 1 << i;
            }
        }
        mask
    }

    /// Current NP-slot index — the position the next sampled base
    /// will land at.
    #[allow(dead_code)]
    pub(crate) fn current_np_index(&self) -> u32 {
        self.current_np_index
    }
}

impl SimulationEventSink for ProductiveAdmitMaskObserver<'_> {
    fn on_event(&mut self, event: &SimulationEvent) {
        match *event {
            // Only bases of our own NP segment advance the slot
            // counter. Bases from other segments (assembly pushes
            // inside the same builder lifetime, in principle) are
            // ignored.
            SimulationEvent::BasePushed { segment, .. } if segment == self.np_segment => {
                self.current_np_index = self.current_np_index.saturating_add(1);
            }
            // Everything else: silently ignore. NP-region sampling
            // doesn't run alongside base-change or indel events, but
            // the sink contract is to match-and-ignore rather than
            // panic if a foreign event arrives.
            _ => {}
        }
    }
}

/// Canonical base layout for the admit mask. Index `i` of this
/// array corresponds to bit `i` of the mask.
const BASES: [u8; 4] = [b'A', b'C', b'G', b'T'];

#[cfg(test)]
mod tests {
    use super::*;
    use crate::ir::{NucFlags, NucHandle};

    /// Build a minimal `JunctionStopState`. The sink-channel parity
    /// tests below never call `current_admit_mask`, so this state's
    /// contents are irrelevant — only the slot counter behaviour is
    /// under test here. The full admit-mask vs. contract-query
    /// parity is covered by `tests::contract_query_equivalence` in
    /// the integration suite.
    fn empty_state() -> JunctionStopState {
        JunctionStopState::empty_for_observer_tests()
    }

    fn push_event(segment: Segment) -> SimulationEvent {
        SimulationEvent::BasePushed {
            handle: NucHandle::new(0),
            base: b'A',
            segment,
            germline_pos: None,
            flags: NucFlags::empty(),
        }
    }

    #[test]
    fn fresh_observer_starts_at_slot_zero() {
        let state = empty_state();
        let obs = ProductiveAdmitMaskObserver::new(&state, Segment::Np1);
        assert_eq!(obs.current_np_index(), 0);
    }

    #[test]
    fn np_segment_push_bumps_slot_counter() {
        let state = empty_state();
        let mut obs = ProductiveAdmitMaskObserver::new(&state, Segment::Np1);
        obs.on_event(&push_event(Segment::Np1));
        obs.on_event(&push_event(Segment::Np1));
        obs.on_event(&push_event(Segment::Np1));
        assert_eq!(obs.current_np_index(), 3);
    }

    #[test]
    fn non_np_segment_push_does_not_advance() {
        // Pushes for V/D/J/Np2 must not move the Np1 observer's slot.
        let state = empty_state();
        let mut obs = ProductiveAdmitMaskObserver::new(&state, Segment::Np1);
        for &seg in &[Segment::V, Segment::D, Segment::J, Segment::Np2] {
            obs.on_event(&push_event(seg));
        }
        assert_eq!(obs.current_np_index(), 0);
    }

    #[test]
    fn base_changed_indel_inserted_indel_deleted_are_noops() {
        // Parity with the legacy `IrEventObserver` impl, which
        // inherited the default no-op trait methods for these
        // variants. The admit-mask never reacts to non-push events.
        let state = empty_state();
        let mut obs = ProductiveAdmitMaskObserver::new(&state, Segment::Np1);
        obs.on_event(&SimulationEvent::BaseChanged {
            handle: NucHandle::new(0),
            old_base: b'A',
            new_base: b'G',
            segment: Segment::Np1,
            germline_pos: None,
        });
        obs.on_event(&SimulationEvent::IndelInserted {
            at: 0,
            base: b'N',
            segment: Segment::Np1,
            flags: NucFlags::empty(),
        });
        obs.on_event(&SimulationEvent::IndelDeleted {
            at: 0,
            removed_base: b'A',
            segment: Segment::V,
        });
        assert_eq!(obs.current_np_index(), 0);
    }

    #[test]
    fn reserved_variants_are_silently_ignored() {
        // Future-emission variants must not corrupt the counter
        // until they're consciously wired in. Covers every
        // reserved variant the enum can carry today — the catch-
        // all `_ => {}` arm in the sink impl means it's safe by
        // construction, but pinning it here keeps the no-op
        // guarantee visible alongside the live `BasePushed` arm.
        use crate::assignment::{AlleleInstance, TrimEnd};
        use crate::ir::Region;
        use crate::refdata::AlleleId;
        let state = empty_state();
        let mut obs = ProductiveAdmitMaskObserver::new(&state, Segment::Np1);
        obs.on_event(&SimulationEvent::BaseDeleted {
            at: 0,
            removed_base: b'A',
        });
        obs.on_event(&SimulationEvent::AssignmentChanged {
            segment: Segment::V,
            old: None,
            new: AlleleInstance::new(AlleleId::new(0)),
        });
        obs.on_event(&SimulationEvent::TrimChanged {
            segment: Segment::V,
            end: TrimEnd::Five,
            old: Some(0),
            new: 2,
        });
        obs.on_event(&SimulationEvent::RegionAdded {
            region: Region::new(Segment::V, NucHandle::new(0), NucHandle::new(3)),
        });
        obs.on_event(&SimulationEvent::RegionReplaced {
            old: Region::new(Segment::V, NucHandle::new(0), NucHandle::new(3)),
            new: Region::new(Segment::V, NucHandle::new(0), NucHandle::new(5)),
        });
        obs.on_event(&SimulationEvent::ReverseComplementFlagRecorded { applied: true });
        obs.on_event(&SimulationEvent::MutationCountChanged {
            old: 0,
            new: 1,
            delta: 1,
        });
        assert_eq!(obs.current_np_index(), 0);
    }
}
