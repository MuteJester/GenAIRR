use super::{Allele, AlleleId, AllelePool, ChainType};
use crate::ir::Segment;

/// Top-level immutable reference configuration for one simulation
/// target (e.g., human IGH).
///
/// In production, this is loaded from a `.v6dat` file (D10) at
/// `Experiment.compile()` time. In tests, build directly via
/// `RefDataConfig::builder()` or by constructing the fields
/// manually.
///
/// **Invariant:** for `chain_type == ChainType::Vj`, the `d_pool` is
/// expected to be empty. Construction does not enforce this — callers
/// of the builder API in tests are expected to honor it; the
/// assembly pass in C.8 will explicitly skip D for VJ chains.
#[derive(Clone, Debug)]
pub struct RefDataConfig {
    pub chain_type: ChainType,
    pub v_pool: AllelePool,
    pub d_pool: AllelePool,
    pub j_pool: AllelePool,
    pub c_pool: AllelePool,
}

impl RefDataConfig {
    /// Empty config for the given chain type. Use the builder /
    /// direct field assignment to populate the pools.
    pub fn empty(chain_type: ChainType) -> Self {
        Self {
            chain_type,
            v_pool: AllelePool::new(),
            d_pool: AllelePool::new(),
            j_pool: AllelePool::new(),
            c_pool: AllelePool::new(),
        }
    }

    /// Resolve an allele by its segment role and id. Returns `None`
    /// if the segment isn't in this config's pools (e.g., asking for
    /// a D allele on a VJ chain) or if the id is out of bounds.
    pub fn get(&self, segment: Segment, id: AlleleId) -> Option<&Allele> {
        match segment {
            Segment::V => self.v_pool.get(id),
            Segment::D => self.d_pool.get(id),
            Segment::J => self.j_pool.get(id),
            // NP regions don't have alleles; return None defensively
            // rather than panicking — a caller asking for an NP
            // allele is misusing the API.
            Segment::Np1 | Segment::Np2 => None,
        }
    }

    /// Pool for the given segment, or `None` for NP segments.
    pub fn pool_for(&self, segment: Segment) -> Option<&AllelePool> {
        match segment {
            Segment::V => Some(&self.v_pool),
            Segment::D => Some(&self.d_pool),
            Segment::J => Some(&self.j_pool),
            Segment::Np1 | Segment::Np2 => None,
        }
    }
}
