/// The biological segment a nucleotide or region belongs to.
///
/// Order matters for assembly: V -> NP1 -> D -> NP2 -> J for VDJ chains;
/// V -> NP1 -> J for VJ chains (NP1 spans the V-J interval and there
/// are no D / NP2 entries).
///
/// Discriminant values are explicit so they can be used to index into
/// [`super::PerSegment`] in O(1). Variants are kept tightly packed
/// (0..[`Self::COUNT`]) so the indexing is safe without bounds
/// checks.
#[derive(Copy, Clone, Eq, PartialEq, Hash, Debug)]
#[repr(u8)]
pub enum Segment {
    V = 0,
    Np1 = 1,
    D = 2,
    Np2 = 3,
    J = 4,
}

impl Segment {
    /// Total number of `Segment` variants. Sized to match the highest
    /// discriminant + 1. Used by [`super::PerSegment`] as the array
    /// length for slot storage.
    pub const COUNT: usize = 5;

    /// Canonical list of "assignable" segments — the V/D/J biology
    /// slots that can carry an allele instance or a live call. NP
    /// segments are excluded (they hold non-templated bases, not
    /// allele-backed segments).
    ///
    /// Centralised here so future extensions (C-region for constant
    /// regions, paired-chain segments, etc.) are one-line additions
    /// to this constant rather than scattered `[V, D, J]` array
    /// literals across the codebase.
    pub const fn assignable() -> &'static [Segment] {
        &[Segment::V, Segment::D, Segment::J]
    }

    /// Stable index into a `PerSegment`-style array. Equals the
    /// enum discriminant.
    #[inline]
    pub const fn index(self) -> usize {
        self as usize
    }
}
