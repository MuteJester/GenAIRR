use crate::ir::Segment;

pub(super) const CONTRACT_NAME: &str = "no_stop_codon_in_junction";

/// What kind of byte lives at a given junction offset.
#[derive(Copy, Clone, Debug, Eq, PartialEq)]
pub(super) enum SlotSource {
    /// A known fixed byte (V-tail, D body, or J-head). Value is in
    /// the parent state's `fixed_bytes[fixed_index]`.
    Fixed { fixed_index: u32 },
    /// An NP byte. Resolves to:
    /// - the candidate if `(segment, np_index)` matches the query,
    /// - else `sim.pool[pool_pos]` if that handle is < `np_committed_pool_len`,
    /// - else `None` (future not-yet-drawn).
    Np {
        segment: Segment,
        np_index: u32,
        /// Pool position where this NP byte will land if/when drawn.
        /// Computed at build time from the NP region's pool start
        /// (already committed) or the projected pool length (not yet
        /// committed).
        pool_pos: u32,
    },
}

/// One row of the precomputed junction layout.
#[derive(Copy, Clone, Debug)]
pub(super) struct JunctionSlot {
    pub(super) source: SlotSource,
}

/// A codon found at build time to already contain a stop using
/// only fixed bytes. The first such codon's offset and bases is
/// stored so query can return a matching error.
#[derive(Copy, Clone, Debug)]
pub(super) struct StaticViolation {
    pub(super) offset: u32,
    pub(super) codon: [u8; 3],
}
