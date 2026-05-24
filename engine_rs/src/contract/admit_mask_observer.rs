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
//! `ProductiveAdmitMaskObserver` is an [`IrEventObserver`]
//! implementer that exercises the trait against a per-slot
//! **admissibility caching** consumer (rather than per-allele
//! scoring or per-position translation).
//!
//! ## Reactive state model
//!
//! The observer borrows the precomputed [`JunctionStopState`] and
//! tracks one integer — the next NP-slot index to be drawn. Every
//! time the assembly pass (or here, the NP pass) commits a base
//! belonging to the observer's `np_segment`, the observer
//! increments its counter. The 4-bit admit mask for the
//! *current* slot is computed lazily, on demand, via
//! [`Self::current_admit_mask`] — at which point the simulation
//! pool already reflects every previously-committed NP byte, so
//! the mask reads consistent state.
//!
//! Bit layout of the mask:
//! - bit 0 → `A`
//! - bit 1 → `C`
//! - bit 2 → `G`
//! - bit 3 → `T`
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

use crate::ir::builder::IrEventObserver;
use crate::ir::{NucHandle, Nucleotide, Segment, Simulation};

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
    /// inbound `on_base_pushed` event belongs to *our* NP region or
    /// to (say) an adjacent assembly that doesn't advance our slot
    /// counter.
    np_segment: Segment,
    /// Index of the next NP slot to be sampled. Starts at 0 and
    /// increments by 1 on every `on_base_pushed` whose nucleotide is
    /// in our `np_segment`.
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

impl IrEventObserver for ProductiveAdmitMaskObserver<'_> {
    fn on_base_pushed(&mut self, _handle: NucHandle, n: &Nucleotide) {
        // Only bases of our own NP segment advance the slot counter.
        // Bases from other segments (assembly pushes inside the same
        // builder lifetime, in principle) are ignored — this is
        // why the trait's `n: &Nucleotide` parameter is essential
        // here (the walker ignores it; the codon rail uses only
        // `n.base`; we use `n.segment`).
        if n.segment == self.np_segment {
            self.current_np_index = self.current_np_index.saturating_add(1);
        }
    }
}

/// Canonical base layout for the admit mask. Index `i` of this
/// array corresponds to bit `i` of the mask.
const BASES: [u8; 4] = [b'A', b'C', b'G', b'T'];
