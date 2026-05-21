//! `JunctionStopState` — precomputed per-record state that makes the
//! per-candidate `NoStopCodonInJunction` check O(1).
//!
//! `JunctionStopState` precomputes the junction layout that stays
//! constant during an NP draw loop: V-tail bytes, optional D body,
//! J-head bytes, and NP slots. Query time then re-translates only the
//! codon touched by the current candidate.

use super::ContractViolation;
use model::{JunctionSlot, StaticViolation};

mod build;
mod layout;
mod model;
mod query;

/// Per-record precomputed state for the `NoStopCodonInJunction`
/// fast path. Build it once before an `NP` draw loop; query
/// `admits_np_candidate` per candidate base.
#[derive(Clone, Debug)]
pub struct JunctionStopState {
    /// Junction offsets in [0, slots.len()). `slots[i]` describes
    /// what byte source lives at junction offset `i`.
    slots: Vec<JunctionSlot>,
    /// Concatenation of all known fixed bytes (V-tail, D body,
    /// J-head). Indexed by `SlotSource::Fixed { fixed_index }`.
    fixed_bytes: Vec<u8>,
    /// First fixed-only codon containing a stop, if any. Causes
    /// every candidate query to return `Err` — matching the slow
    /// path's behavior where the rebuild-and-walk path would
    /// detect this stop regardless of the candidate.
    static_violation: Option<StaticViolation>,
}

#[cfg(test)]
mod tests;
