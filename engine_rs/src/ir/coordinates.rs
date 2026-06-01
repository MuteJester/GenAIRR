//! Coordinate algebra kernel.
//!
//! Centralises all boundary arithmetic for half-open coordinate
//! ranges. Every `<` / `<=` / `>` / `>=` decision on a segment or
//! window boundary is expressed once, here, and consumed by the
//! refresh path, indel projector, dirty-window overlap, and AIRR
//! walker through typed primitives.
//!
//! **Why this lives in one file.** Boundary inequalities were the
//! root cause of every live-call refresh divergence the cache parity
//! harness has surfaced. Two examples:
//!
//! - `segment_region_overlaps_dirty` originally used strict `<` on
//!   the upper bound; under 3' end-loss the dirty window's `start`
//!   landed at exactly `region.end` (the byte deletion just shrunk
//!   out of the region) and the refresh was silently skipped. The
//!   fix widened to `<=` — captured here as [`PoolRange::overlaps_dirty`].
//!
//! - `WalkerObserverState::on_indel_inserted` originally used `<=`
//!   on the inside-range check; this diverged from
//!   `Sequence::with_indel_adjusted`'s strict-`<` shift rule by one
//!   position at the right edge. The fix flipped to `<` — captured
//!   here as [`PoolRange::after_insertion`].
//!
//! Both bugs had the same shape: two call sites both encoded the
//! same boundary rule, and they drifted. Collapsing them through this
//! kernel makes future drift impossible by construction.
//!
//! ## Coordinate spaces
//!
//! - [`PoolRange`] — `[start, end)` over the [`NucleotidePool`](crate::ir::NucleotidePool).
//!   `u32`, matches the underlying [`NucHandle`](crate::ir::NucHandle)
//!   index space. Used for `Region` bounds, dirty windows, and
//!   live-call walker tracking.
//! - [`RefRange`] — `[start, end)` over an allele's reference
//!   sequence. `i64`, matches existing AIRR walker arithmetic
//!   (`-1` is the synthetic-nucleotide "no ref position" sentinel).
//!
//! Both ranges are half-open: `start == end` is empty, `start > end`
//! is invalid and triggers a debug assertion at construction.

use core::cmp::{max, min};

// ──────────────────────────────────────────────────────────────────
// PoolRange
// ──────────────────────────────────────────────────────────────────

/// Half-open `[start, end)` window into the nucleotide pool.
#[derive(Copy, Clone, Debug, Eq, PartialEq)]
pub struct PoolRange {
    pub start: u32,
    pub end: u32,
}

impl PoolRange {
    pub fn new(start: u32, end: u32) -> Self {
        debug_assert!(start <= end, "pool range inverted: {start}..{end}");
        Self { start, end }
    }

    pub fn empty_at(pos: u32) -> Self {
        Self { start: pos, end: pos }
    }

    pub fn len(&self) -> u32 {
        self.end.saturating_sub(self.start)
    }

    pub fn is_empty(&self) -> bool {
        self.start >= self.end
    }

    /// Half-open containment: `pos in [start, end)`.
    pub fn contains(&self, pos: u32) -> bool {
        pos >= self.start && pos < self.end
    }

    /// Standard `[start, end)` overlap. Empty ranges never overlap.
    pub fn overlaps(&self, other: PoolRange) -> bool {
        self.start < other.end && other.start < self.end
    }

    /// Overlap with a dirty window.
    ///
    /// **Asymmetric on purpose.** Unlike [`Self::overlaps`], the
    /// upper bound on `self` is inclusive of `dirty.start`: a dirty
    /// window starting AT `self.end` represents an event at exactly
    /// the boundary position the region just shrunk out of (an
    /// `IndelDeleted` at `region.end - 1` post-shift). The pre-event
    /// position was inside the region; the dirty signal must trigger
    /// a refresh even though post-event `region.end == dirty.start`.
    ///
    /// The lower bound stays strict: a dirty event at `self.start - 1`
    /// shifts the region down without changing content, so no refresh
    /// is needed.
    ///
    /// History: strict-`<` on the upper bound stranded stale live
    /// calls under 3' end-loss on IGK J; the AIRR validator caught it.
    pub fn overlaps_dirty(&self, dirty: PoolRange) -> bool {
        dirty.start <= self.end && dirty.end > self.start
    }

    /// `self` and `other` share an exact boundary, with `other`
    /// sitting immediately to the right of `self` (i.e. `self.end ==
    /// other.start`). Used to find an NP region adjacent to a
    /// segment's right edge.
    pub fn touches_right(&self, other: PoolRange) -> bool {
        self.end == other.start
    }

    /// Symmetric to [`Self::touches_right`]: `other` ends exactly at
    /// `self`'s start (`other.end == self.start`).
    pub fn touches_left(&self, other: PoolRange) -> bool {
        other.end == self.start
    }

    /// Adjusted range after a single-nucleotide insertion at pool
    /// position `at`.
    ///
    /// Strict shift semantics:
    ///
    /// - `at < self.start` — insertion is before the range: shift
    ///   both bounds up by 1.
    /// - `at >= self.start && at < self.end` — insertion lands
    ///   inside the range: stretch `end` up by 1; `start` stays put.
    /// - `at >= self.end` — insertion is at or past the right edge:
    ///   range unchanged.
    ///
    /// `at == self.end` is **not** inside: the new byte lands one
    /// position past the range. Bumping `end` here would diverge
    /// cached walker state from a from-scratch recompute by one
    /// position. (See module docstring history.)
    pub fn after_insertion(&self, at: u32) -> Self {
        if at < self.start {
            Self::new(self.start + 1, self.end + 1)
        } else if at < self.end {
            Self::new(self.start, self.end + 1)
        } else {
            *self
        }
    }

    /// Adjusted range after a single-nucleotide deletion at pool
    /// position `at`.
    ///
    /// - `at < self.start` — deletion is before the range: shift
    ///   both bounds down by 1.
    /// - `at >= self.start && at < self.end` — deletion lands inside
    ///   the range: shrink `end` by 1; `start` stays.
    /// - `at >= self.end` — range unchanged.
    pub fn after_deletion(&self, at: u32) -> Self {
        if at < self.start {
            Self::new(self.start.saturating_sub(1), self.end.saturating_sub(1))
        } else if at < self.end {
            Self::new(self.start, self.end.saturating_sub(1))
        } else {
            *self
        }
    }

    /// Adjusted range after a bulk segment replacement.
    ///
    /// `replaced` is the old pool span being excised; `new_len` is
    /// the length of the byte sequence being installed in its place
    /// (so the new occupied span is
    /// `[replaced.start, replaced.start + new_len)`). `self` is some
    /// OTHER region in the same simulation whose coordinates need to
    /// be re-anchored after the swap.
    ///
    /// Receptor-revision Slice A is the use case: a V replacement
    /// shifts the downstream Np1 / D / Np2 / J regions in lockstep
    /// by `new_len - replaced.len()`.
    ///
    /// Cases:
    ///
    /// - Fully upstream (`self.end <= replaced.start`): unchanged.
    /// - Fully downstream (`self.start >= replaced.end`): both
    ///   bounds shifted by `new_len - replaced.len()`.
    /// - Otherwise (intersection or contained-by): undefined.
    ///   The receptor-revision builder enforces the one-region-per-
    ///   segment invariant before calling, so the only intersecting
    ///   region is the replacement itself — and the caller hands the
    ///   *new* region to its consumer instead of transforming the
    ///   old one through this kernel.
    pub fn after_segment_replacement(&self, replaced: PoolRange, new_len: u32) -> Self {
        if self.end <= replaced.start {
            return *self;
        }
        debug_assert!(
            self.start >= replaced.end,
            "after_segment_replacement: range {self:?} intersects replaced span {replaced:?}",
        );
        let old_len = replaced.len();
        if new_len >= old_len {
            let shift = new_len - old_len;
            Self::new(self.start + shift, self.end + shift)
        } else {
            let shift = old_len - new_len;
            Self::new(self.start - shift, self.end - shift)
        }
    }

    pub fn intersection(&self, other: PoolRange) -> Option<PoolRange> {
        let s = max(self.start, other.start);
        let e = min(self.end, other.end);
        if s < e {
            Some(Self::new(s, e))
        } else {
            None
        }
    }
}

// ──────────────────────────────────────────────────────────────────
// RefRange
// ──────────────────────────────────────────────────────────────────

/// Half-open `[start, end)` window into an allele's reference
/// sequence. `i64` matches AIRR walker arithmetic.
#[derive(Copy, Clone, Debug, Eq, PartialEq)]
pub struct RefRange {
    pub start: i64,
    pub end: i64,
}

impl RefRange {
    pub fn new(start: i64, end: i64) -> Self {
        debug_assert!(start <= end, "ref range inverted: {start}..{end}");
        Self { start, end }
    }

    /// Half-open containment: `pos in [start, end)`.
    pub fn contains(&self, pos: i64) -> bool {
        pos >= self.start && pos < self.end
    }

    /// Return a range that covers `self` plus the single position `pos`.
    pub fn extend_to_cover(&self, pos: i64) -> Self {
        Self::new(self.start.min(pos), self.end.max(pos + 1))
    }
}

// ──────────────────────────────────────────────────────────────────
// tests
// ──────────────────────────────────────────────────────────────────

#[cfg(test)]
mod tests {
    use super::*;

    // ── PoolRange ────────────────────────────────────────────────

    #[test]
    fn pool_range_contains_is_half_open() {
        let r = PoolRange::new(3, 7);
        assert!(!r.contains(2));
        assert!(r.contains(3));
        assert!(r.contains(6));
        assert!(!r.contains(7));
        assert!(!r.contains(8));
    }

    #[test]
    fn pool_range_overlaps_strict_on_both_ends() {
        let r = PoolRange::new(3, 7);
        assert!(!r.overlaps(PoolRange::new(0, 3)));
        assert!(r.overlaps(PoolRange::new(0, 4)));
        assert!(r.overlaps(PoolRange::new(6, 10)));
        assert!(!r.overlaps(PoolRange::new(7, 10)));
        assert!(!r.overlaps(PoolRange::new(3, 3)));
    }

    #[test]
    fn pool_range_overlaps_dirty_inclusive_upper_strict_lower() {
        let r = PoolRange::new(3, 7);
        assert!(!r.overlaps_dirty(PoolRange::new(0, 3)));
        assert!(r.overlaps_dirty(PoolRange::new(0, 4)));
        assert!(r.overlaps_dirty(PoolRange::new(7, 7)));
        assert!(r.overlaps_dirty(PoolRange::new(7, 10)));
        assert!(!r.overlaps_dirty(PoolRange::new(8, 10)));
    }

    #[test]
    fn pool_range_after_insertion_shift_rules() {
        let r = PoolRange::new(3, 7);
        assert_eq!(r.after_insertion(0), PoolRange::new(4, 8));
        assert_eq!(r.after_insertion(2), PoolRange::new(4, 8));
        assert_eq!(r.after_insertion(3), PoolRange::new(3, 8));
        assert_eq!(r.after_insertion(6), PoolRange::new(3, 8));
        assert_eq!(r.after_insertion(7), PoolRange::new(3, 7));
        assert_eq!(r.after_insertion(10), PoolRange::new(3, 7));
    }

    #[test]
    fn pool_range_after_deletion_shift_rules() {
        let r = PoolRange::new(3, 7);
        assert_eq!(r.after_deletion(0), PoolRange::new(2, 6));
        assert_eq!(r.after_deletion(2), PoolRange::new(2, 6));
        assert_eq!(r.after_deletion(3), PoolRange::new(3, 6));
        assert_eq!(r.after_deletion(6), PoolRange::new(3, 6));
        assert_eq!(r.after_deletion(7), PoolRange::new(3, 7));
        assert_eq!(r.after_deletion(10), PoolRange::new(3, 7));
    }

    #[test]
    fn pool_range_after_segment_replacement_same_length_is_identity() {
        let replaced = PoolRange::new(5, 10); // len 5
        let downstream = PoolRange::new(10, 14);
        let upstream = PoolRange::new(0, 5);
        assert_eq!(downstream.after_segment_replacement(replaced, 5), downstream);
        assert_eq!(upstream.after_segment_replacement(replaced, 5), upstream);
    }

    #[test]
    fn pool_range_after_segment_replacement_grow_shifts_downstream_right() {
        let replaced = PoolRange::new(5, 10); // len 5 → len 8 (+3)
        let np1 = PoolRange::new(10, 12);
        let d = PoolRange::new(12, 15);
        assert_eq!(np1.after_segment_replacement(replaced, 8), PoolRange::new(13, 15));
        assert_eq!(d.after_segment_replacement(replaced, 8), PoolRange::new(15, 18));
    }

    #[test]
    fn pool_range_after_segment_replacement_shrink_shifts_downstream_left() {
        let replaced = PoolRange::new(5, 10); // len 5 → len 2 (-3)
        let np1 = PoolRange::new(10, 12);
        let j = PoolRange::new(15, 20);
        assert_eq!(np1.after_segment_replacement(replaced, 2), PoolRange::new(7, 9));
        assert_eq!(j.after_segment_replacement(replaced, 2), PoolRange::new(12, 17));
    }

    #[test]
    fn pool_range_after_segment_replacement_upstream_untouched_even_when_grown() {
        let replaced = PoolRange::new(10, 20); // len 10 → len 30 (+20)
        let upstream = PoolRange::new(0, 10);
        assert_eq!(upstream.after_segment_replacement(replaced, 30), upstream);
    }

    #[test]
    #[should_panic(expected = "intersects replaced span")]
    fn pool_range_after_segment_replacement_overlap_panics_in_debug() {
        let replaced = PoolRange::new(5, 10);
        let overlapping = PoolRange::new(7, 12);
        let _ = overlapping.after_segment_replacement(replaced, 5);
    }

    #[test]
    fn pool_range_touches() {
        let a = PoolRange::new(3, 7);
        let b = PoolRange::new(7, 9);
        let c = PoolRange::new(1, 3);
        assert!(a.touches_right(b));
        assert!(!a.touches_right(c));
        assert!(a.touches_left(c));
        assert!(!a.touches_left(b));
    }

    #[test]
    fn pool_range_intersection() {
        let a = PoolRange::new(3, 7);
        assert_eq!(a.intersection(PoolRange::new(5, 10)), Some(PoolRange::new(5, 7)));
        assert_eq!(a.intersection(PoolRange::new(0, 3)), None);
        assert_eq!(a.intersection(PoolRange::new(7, 10)), None);
    }

    // ── RefRange ─────────────────────────────────────────────────

    #[test]
    fn ref_range_contains_is_half_open() {
        let r = RefRange::new(0, 10);
        assert!(r.contains(0));
        assert!(r.contains(9));
        assert!(!r.contains(10));
        assert!(!r.contains(-1));
    }

    #[test]
    fn ref_range_extend_to_cover_grows_both_sides() {
        let r = RefRange::new(5, 8);
        assert_eq!(r.extend_to_cover(3), RefRange::new(3, 8));
        assert_eq!(r.extend_to_cover(10), RefRange::new(5, 11));
        assert_eq!(r.extend_to_cover(6), RefRange::new(5, 8));
    }
}
