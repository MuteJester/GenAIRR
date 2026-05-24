//! [`PerSegment`] — array-backed typed map from [`Segment`] to `T`.
//!
//! Replaces hand-rolled `v: Option<T>, d: Option<T>, j: Option<T>`
//! field clusters (and their match-dispatch getters) with a single
//! container that:
//!
//! - Stores values in a fixed-size `[Option<T>; Segment::COUNT]`
//!   array indexed by [`Segment::index`]. Same memory footprint as
//!   the parallel `Option<T>` fields; same O(1) access cost; no
//!   hashing or comparisons.
//! - Iterates populated entries in canonical segment order
//!   (`Segment::assignable()` first, then NP segments).
//! - Lets future segment kinds (e.g. C-region for constant regions,
//!   paired-chain heavy/light) become available without touching
//!   every consumer site.
//!
//! ## Why not `BTreeMap<Segment, T>`?
//!
//! With N ≤ 5 entries the map's hashing/comparison overhead dwarfs
//! the actual work. The array is identical in performance to the
//! previous per-segment `Option` fields while keeping the
//! schema-free extensibility benefit at the consumer API surface.

use super::Segment;

/// Array-backed map from [`Segment`] to `T`. See module docs.
#[derive(Clone, Debug, Eq, PartialEq)]
pub struct PerSegment<T> {
    slots: [Option<T>; Segment::COUNT],
}

impl<T> PerSegment<T> {
    /// An empty map with no populated slots.
    pub const fn new() -> Self {
        Self {
            slots: [const { None }; Segment::COUNT],
        }
    }

    /// Read the value associated with `seg`, if any.
    #[inline]
    pub fn get(&self, seg: Segment) -> Option<&T> {
        self.slots[seg.index()].as_ref()
    }

    /// Mutable access to the value associated with `seg`, if any.
    #[inline]
    #[allow(dead_code)]
    pub fn get_mut(&mut self, seg: Segment) -> Option<&mut T> {
        self.slots[seg.index()].as_mut()
    }

    /// Insert (or replace) the value at `seg`. Returns the previous
    /// value, if any.
    #[inline]
    pub fn set(&mut self, seg: Segment, value: T) -> Option<T> {
        self.slots[seg.index()].replace(value)
    }

    /// Remove and return the value at `seg`, if any.
    #[inline]
    #[allow(dead_code)]
    pub fn take(&mut self, seg: Segment) -> Option<T> {
        self.slots[seg.index()].take()
    }

    /// `true` when `seg` has a populated value.
    #[inline]
    pub fn contains(&self, seg: Segment) -> bool {
        self.slots[seg.index()].is_some()
    }

    /// Iterate populated `(segment, &value)` entries in
    /// [`Segment::assignable`] order (V → D → J), then NP segments
    /// if any are ever populated. Skips empty slots.
    pub fn iter(&self) -> impl Iterator<Item = (Segment, &T)> + '_ {
        const ALL_SEGMENTS: [Segment; Segment::COUNT] = [
            Segment::V,
            Segment::D,
            Segment::J,
            Segment::Np1,
            Segment::Np2,
        ];
        ALL_SEGMENTS
            .iter()
            .filter_map(move |&seg| self.get(seg).map(|v| (seg, v)))
    }
}

impl<T> Default for PerSegment<T> {
    fn default() -> Self {
        Self::new()
    }
}

// `Copy` when `T: Copy`. We can't blanket-impl Copy on the generic
// because `Option<T>: Copy` requires `T: Copy`; rustc derives this
// for us when we use derive(Copy), but PerSegment uses an array
// init that the derive macro can't yet handle generically. Manual
// `Copy` for the common case where `T` is `Copy` (e.g.
// `AlleleInstance`).
impl<T: Copy> Copy for PerSegment<T> {}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn new_is_empty_across_all_segments() {
        let p: PerSegment<u32> = PerSegment::new();
        for &seg in &[
            Segment::V,
            Segment::Np1,
            Segment::D,
            Segment::Np2,
            Segment::J,
        ] {
            assert!(!p.contains(seg));
            assert!(p.get(seg).is_none());
        }
    }

    #[test]
    fn set_and_get_round_trip() {
        let mut p: PerSegment<u32> = PerSegment::new();
        assert_eq!(p.set(Segment::V, 42), None);
        assert_eq!(p.get(Segment::V), Some(&42));
        assert!(p.contains(Segment::V));
        assert!(!p.contains(Segment::D));
    }

    #[test]
    fn set_replaces_and_returns_previous() {
        let mut p: PerSegment<u32> = PerSegment::new();
        p.set(Segment::V, 1);
        assert_eq!(p.set(Segment::V, 2), Some(1));
        assert_eq!(p.get(Segment::V), Some(&2));
    }

    #[test]
    fn take_returns_and_clears() {
        let mut p: PerSegment<u32> = PerSegment::new();
        p.set(Segment::D, 7);
        assert_eq!(p.take(Segment::D), Some(7));
        assert!(!p.contains(Segment::D));
        assert_eq!(p.take(Segment::D), None);
    }

    #[test]
    fn iter_skips_empties_and_yields_assignable_first() {
        let mut p: PerSegment<u32> = PerSegment::new();
        p.set(Segment::J, 3);
        p.set(Segment::V, 1);
        p.set(Segment::Np1, 99);
        let collected: Vec<(Segment, u32)> = p.iter().map(|(s, v)| (s, *v)).collect();
        // V comes before J (assignable order), and Np1 comes after
        // the assignables.
        assert_eq!(
            collected,
            vec![(Segment::V, 1), (Segment::J, 3), (Segment::Np1, 99)]
        );
    }

    #[test]
    fn copy_works_for_copyable_payloads() {
        let mut p: PerSegment<u32> = PerSegment::new();
        p.set(Segment::V, 5);
        let q = p; // Copy
        assert_eq!(p.get(Segment::V), Some(&5));
        assert_eq!(q.get(Segment::V), Some(&5));
    }
}
