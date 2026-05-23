use super::{Allele, AlleleId};

/// All alleles available for one segment role (V, D, J, or C).
///
/// Indexed by `AlleleId`. The pool is the reference data; *which*
/// allele a particular simulation sampled is recorded by an
/// `AlleleInstance` referring to a stable `AlleleId`.
#[derive(Clone, Debug, Default)]
pub struct AllelePool {
    alleles: Vec<Allele>,
}

impl AllelePool {
    pub fn new() -> Self {
        Self::default()
    }

    /// Construct from an existing `Vec<Allele>`. Caller is responsible
    /// for ensuring all alleles have the right segment role for the
    /// pool being built.
    pub fn from_vec(alleles: Vec<Allele>) -> Self {
        Self { alleles }
    }

    /// Number of alleles in the pool.
    pub fn len(&self) -> usize {
        self.alleles.len()
    }

    /// Whether the pool contains zero alleles.
    pub fn is_empty(&self) -> bool {
        self.alleles.is_empty()
    }

    /// Append an allele. Returns the issued `AlleleId` so callers can
    /// reference the just-added allele in subsequent setup. This is
    /// the construction-time API; reference data is sealed before any
    /// simulation runs against it.
    #[must_use = "AllelePool::push returns the issued AlleleId; use \
                  `_ = pool.push(...)` if it isn't needed"]
    pub fn push(&mut self, allele: Allele) -> AlleleId {
        let id = AlleleId::new(self.alleles.len() as u32);
        self.alleles.push(allele);
        id
    }

    /// Look up an allele by id. Returns `None` if `id` is out of
    /// bounds — defensive to allow `RefDataConfig` to expose
    /// fallible lookup at the boundary even though correct callers
    /// always have valid handles.
    pub fn get(&self, id: AlleleId) -> Option<&Allele> {
        self.alleles.get(id.as_usize())
    }

    /// Read-only slice of all alleles in pool order. Iterated by
    /// `AllelePoolDist` to build cumulative frequency tables.
    pub fn as_slice(&self) -> &[Allele] {
        &self.alleles
    }

    /// Iterator over `(AlleleId, &Allele)` pairs in pool order.
    pub fn iter(&self) -> impl Iterator<Item = (AlleleId, &Allele)> {
        self.alleles
            .iter()
            .enumerate()
            .map(|(i, a)| (AlleleId::new(i as u32), a))
    }

    /// Find the first allele whose name matches exactly. O(N).
    /// Used in tests and tooling; production sampling goes through
    /// `AllelePoolDist` (C.3) and never name-resolves here.
    pub fn find_by_name(&self, name: &str) -> Option<(AlleleId, &Allele)> {
        self.iter().find(|(_, a)| a.name == name)
    }
}
