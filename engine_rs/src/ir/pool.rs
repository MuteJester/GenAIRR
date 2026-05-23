use std::sync::Arc;

use super::{NucHandle, Nucleotide};

/// Arena of all nucleotides in a `Simulation`. `NucHandle` indexes
/// into this pool.
#[derive(Clone, Debug, Default)]
pub struct NucleotidePool {
    nucleotides: Arc<Vec<Nucleotide>>,
}

impl NucleotidePool {
    /// New empty pool.
    pub fn new() -> Self {
        Self::default()
    }

    /// New empty pool pre-allocated for `cap` nucleotides.
    pub fn with_capacity(cap: usize) -> Self {
        Self {
            nucleotides: Arc::new(Vec::with_capacity(cap)),
        }
    }

    /// Number of nucleotides currently in the pool.
    pub fn len(&self) -> usize {
        self.nucleotides.len()
    }

    /// Whether the pool contains zero nucleotides.
    pub fn is_empty(&self) -> bool {
        self.nucleotides.is_empty()
    }

    /// Read-only access to a nucleotide by handle.
    pub fn get(&self, h: NucHandle) -> Option<&Nucleotide> {
        self.nucleotides.get(h.as_usize())
    }

    /// Read-only access to the underlying slice.
    pub fn as_slice(&self) -> &[Nucleotide] {
        &self.nucleotides
    }

    /// Ensure the pool's inner `Arc<Vec<Nucleotide>>` is unique
    /// (refcount == 1) so subsequent in-place mutations don't deep-clone.
    ///
    /// Used by `SimulationBuilder::from_simulation` to prepay the
    /// cost-of-write *once* before a hot per-base loop, after which
    /// every internal `Arc::make_mut` in `push` / `extend` is the cheap
    /// refcount==1 path.
    pub(crate) fn make_unique(&mut self) {
        let _ = Arc::make_mut(&mut self.nucleotides);
    }

    /// Append a nucleotide to the pool, returning its handle.
    pub fn push(&mut self, n: Nucleotide) -> NucHandle {
        let inner = Arc::make_mut(&mut self.nucleotides);
        let h = NucHandle::new(inner.len() as u32);
        inner.push(n);
        h
    }

    /// Append every nucleotide produced by `iter` into the pool in place.
    pub fn extend<I: IntoIterator<Item = Nucleotide>>(&mut self, iter: I) {
        let inner = Arc::make_mut(&mut self.nucleotides);
        inner.extend(iter);
    }

    #[must_use = "with_pushed returns (NucleotidePool, NucHandle); both must \
                  be used or destructured. Drop the handle explicitly with \
                  `let (next, _) = ...` if it isn't needed."]
    pub fn with_pushed(&self, n: Nucleotide) -> (Self, NucHandle) {
        let mut next = self.clone();
        let h = next.push(n);
        (next, h)
    }

    #[must_use]
    pub fn with_pushed_many<I: IntoIterator<Item = Nucleotide>>(
        &self,
        iter: I,
    ) -> (Self, std::ops::Range<u32>) {
        let mut next = self.clone();
        let start = next.nucleotides.len() as u32;
        next.extend(iter);
        let end = next.nucleotides.len() as u32;
        (next, start..end)
    }

    /// Return a new pool with the nucleotide at `handle` replaced.
    pub fn with_nucleotide_changed(&self, handle: NucHandle, new_n: Nucleotide) -> Self {
        let mut next = self.clone();
        let inner = Arc::make_mut(&mut next.nucleotides);
        inner[handle.as_usize()] = new_n;
        next
    }

    /// Return a new pool with only the `base` field changed.
    pub fn with_base_changed(&self, handle: NucHandle, new_base: u8) -> Self {
        let mut next = self.clone();
        let inner = Arc::make_mut(&mut next.nucleotides);
        inner[handle.as_usize()].base = new_base;
        next
    }

    /// Mutate the `base` byte at `handle` in place, returning the
    /// nucleotide as it was *before* the change (caller borrows the
    /// previous value to feed event observers etc.).
    ///
    /// Companion to [`Self::push`]: uses `Arc::make_mut` internally,
    /// which is the refcount==1 cheap path once `make_unique` has
    /// been called on the owning `SimulationBuilder`.
    pub(crate) fn change_base_in_place(&mut self, handle: NucHandle, new_base: u8) -> Nucleotide {
        let inner = Arc::make_mut(&mut self.nucleotides);
        let idx = handle.as_usize();
        let old = inner[idx];
        inner[idx].base = new_base;
        old
    }

    /// Return a new pool with a fresh nucleotide inserted at position `at`.
    pub fn with_inserted(&self, at: u32, n: Nucleotide) -> Self {
        assert!(
            (at as usize) <= self.nucleotides.len(),
            "NucleotidePool::with_inserted: at ({}) > pool length ({})",
            at,
            self.nucleotides.len()
        );
        let mut next = self.clone();
        let inner = Arc::make_mut(&mut next.nucleotides);
        inner.insert(at as usize, n);
        next
    }

    /// Return a new pool with the nucleotide at `at` removed.
    pub fn with_deleted(&self, at: u32) -> Self {
        assert!(
            (at as usize) < self.nucleotides.len(),
            "NucleotidePool::with_deleted: at ({}) >= pool length ({})",
            at,
            self.nucleotides.len()
        );
        let mut next = self.clone();
        let inner = Arc::make_mut(&mut next.nucleotides);
        inner.remove(at as usize);
        next
    }
}
