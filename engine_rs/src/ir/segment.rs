/// The biological segment a nucleotide or region belongs to.
///
/// Order matters for assembly: V -> NP1 -> D -> NP2 -> J for VDJ chains;
/// V -> NP1 -> J for VJ chains (NP1 spans the V-J interval and there
/// are no D / NP2 entries).
#[derive(Copy, Clone, Eq, PartialEq, Hash, Debug)]
pub enum Segment {
    V,
    Np1,
    D,
    Np2,
    J,
}
