/// The biological chain configuration.
///
/// `Vj` chains (light: kappa, lambda, TCR alpha/gamma) recombine V
/// and J only — there is no D segment, NP1 spans the V→J interval,
/// and NP2 / D pools are unused.
///
/// `Vdj` chains (heavy: IGH, TCR beta/delta) recombine V, D, and J
/// with NP1 between V and D and NP2 between D and J.
///
/// Other axes (allele frequency models, isotype constants, etc.)
/// are independent of `ChainType` and live elsewhere in the config.
#[derive(Copy, Clone, Eq, PartialEq, Hash, Debug)]
pub enum ChainType {
    Vj,
    Vdj,
}

impl ChainType {
    /// Whether this chain has a D segment (and therefore an NP2 region).
    pub const fn has_d(self) -> bool {
        matches!(self, ChainType::Vdj)
    }
}
