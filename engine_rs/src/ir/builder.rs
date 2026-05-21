//! `SimulationBuilder` — owns an in-progress `Simulation` by value and
//! supports per-base appends without per-step refcount drama.
//!
//! ## Why this exists
//!
//! `Simulation::with_nucleotide_pushed` is a persistent operation: it
//! takes `&self`, clones the inner pool's Arc, and CoWs the inner `Vec`
//! through `Arc::make_mut`. When that clone is immediately reassigned
//! to a new local (`current = current.with_nucleotide_pushed(...)`),
//! the prior `current` is *briefly co-owned* during the assignment,
//! so `Arc::make_mut` falls into the deep-clone branch every iteration
//! of the NP draw loop.
//!
//! `SimulationBuilder` solves this by:
//!
//! 1. Taking ownership of a `Simulation` by value at construction.
//! 2. Calling `Arc::make_mut` once up front on the pool's inner `Vec`
//!    so the Arc refcount becomes 1.
//! 3. Holding the now-unique `Simulation` across many `push_nucleotide`
//!    calls — each of which is a direct `Vec::push` (the Arc refcount
//!    stays at 1, so subsequent `make_mut` calls hit the cheap path).
//!
//! ## What this *does not* do
//!
//! - It does not support `change_base_in_place`, `with_indel_*`, region
//!   updates, or any other operation beyond per-base append. The NP pass
//!   is the sole intended consumer for now.
//! - It does not break the persistent IR public contract — the only
//!   `&mut self` operation it exposes returns a `NucHandle` rather than
//!   a new `Simulation`. Callers seal the builder with `seal()` to obtain
//!   the immutable snapshot when done.
//!
//! ## Contract visibility
//!
//! `peek()` returns `&Simulation` whose `pool` reflects every base
//! pushed so far. Contracts reading `sim.pool` (e.g. `read_slot` in
//! `JunctionStopState`) see exactly the same committed prefix they
//! would under the old per-step `with_nucleotide_pushed` loop. This is
//! what preserves bit-identical regression output.

use super::{NucHandle, Nucleotide, Simulation};

/// Mutable, by-value owner of an in-progress `Simulation` whose pool
/// supports cheap per-base append. See the module docs for the design
/// rationale.
pub struct SimulationBuilder {
    simulation: Simulation,
}

impl SimulationBuilder {
    /// Take ownership of `sim` and unique-ify its pool's inner `Vec`
    /// so subsequent `push_nucleotide` calls don't trigger a deep
    /// clone via `Arc::make_mut`.
    pub fn from_simulation(sim: Simulation) -> Self {
        let mut builder = Self { simulation: sim };
        // Force unique ownership of the pool's nucleotide vector. After
        // this, every subsequent `Arc::make_mut` inside `pool.push`
        // is a refcount==1 fast path: it just hands back a
        // `&mut Vec<Nucleotide>` without copying.
        builder.simulation.pool.make_unique();
        builder
    }

    /// Read-only view of the in-progress simulation. Contracts and
    /// other observers see every base committed via `push_nucleotide`
    /// up to this point.
    pub fn peek(&self) -> &Simulation {
        &self.simulation
    }

    /// Append one nucleotide and return its handle. Amortized O(1):
    /// no Arc clone, no `Vec` reallocation if capacity holds.
    pub fn push_nucleotide(&mut self, n: Nucleotide) -> NucHandle {
        // The pool's `push` already does `Arc::make_mut` internally;
        // after `from_simulation`'s prepayment, that call is the cheap
        // refcount==1 path on every iteration.
        self.simulation.pool.push(n)
    }

    /// Consume the builder and return the sealed `Simulation`.
    ///
    /// At seal time the inner pool still holds an `Arc<Vec<Nucleotide>>`
    /// with refcount==1; subsequent persistent operations on the sealed
    /// `Simulation` (e.g. `with_region_added`) clone the outer pool
    /// struct (cheap) without forcing a deep clone of the vector.
    pub fn seal(self) -> Simulation {
        self.simulation
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::ir::{flag, Segment};

    fn make_n(base: u8) -> Nucleotide {
        Nucleotide::synthetic(base, Segment::Np1, flag::N_NUC)
    }

    #[test]
    fn builder_pushes_are_visible_through_peek() {
        let mut b = SimulationBuilder::from_simulation(Simulation::new());
        assert_eq!(b.peek().pool.len(), 0);
        let h0 = b.push_nucleotide(make_n(b'A'));
        let h1 = b.push_nucleotide(make_n(b'C'));
        assert_eq!(h0, NucHandle::new(0));
        assert_eq!(h1, NucHandle::new(1));
        assert_eq!(b.peek().pool.len(), 2);
        assert_eq!(b.peek().pool.get(h0).unwrap().base, b'A');
        assert_eq!(b.peek().pool.get(h1).unwrap().base, b'C');
    }

    #[test]
    fn seal_yields_simulation_with_all_committed_bases() {
        let mut b = SimulationBuilder::from_simulation(Simulation::new());
        for &c in b"GATTACA" {
            b.push_nucleotide(make_n(c));
        }
        let sim = b.seal();
        assert_eq!(sim.pool.len(), 7);
        let bases: Vec<u8> = sim.pool.as_slice().iter().map(|n| n.base).collect();
        assert_eq!(&bases, b"GATTACA");
    }

    #[test]
    fn builder_uniqueifies_shared_pool_on_construction() {
        // Pre-share the pool's inner Vec, then build. The shared
        // reference must not see our pushes — make_mut at construction
        // must have split the buffer.
        let mut sim = Simulation::new();
        sim.pool.push(make_n(b'A'));
        let shared_pool = sim.pool.clone();
        assert_eq!(shared_pool.len(), 1);

        let mut b = SimulationBuilder::from_simulation(sim);
        b.push_nucleotide(make_n(b'T'));
        b.push_nucleotide(make_n(b'C'));

        // Shared snapshot is unchanged.
        assert_eq!(shared_pool.len(), 1);
        assert_eq!(shared_pool.get(NucHandle::new(0)).unwrap().base, b'A');

        // Builder's view reflects the new pushes.
        let sealed = b.seal();
        assert_eq!(sealed.pool.len(), 3);
        let bases: Vec<u8> = sealed.pool.as_slice().iter().map(|n| n.base).collect();
        assert_eq!(&bases, b"ATC");
    }
}
