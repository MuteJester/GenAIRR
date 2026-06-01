//! `SegmentRateWeights` — per-biological-segment SHM rate scalars.
//!
//! Threads the four-bucket rate vector (V / D / J / NP) from the
//! Python DSL through both `UniformMutationPass` and
//! `S5FMutationPass`. Sites within a segment are multiplied by the
//! segment's rate at site-selection time; zero-rate sites drop
//! out of the support before contract admissibility runs (per the
//! audit's §4 constrain-before-propose ordering).
//!
//! Default: all four buckets at `1.0` — flat substrate, byte-
//! identical to the pre-slice engine.
//!
//! See `docs/shm_segment_rate_design.md` for the architecture
//! contract.

use crate::ir::Segment;

/// Per-biological-segment SHM rate scalars.
///
/// Fields are in canonical DSL order (V / D / J / NP) so positional
/// plumbing through PyO3 stays straightforward. NP folds both
/// `Segment::Np1` and `Segment::Np2`.
///
/// Construction invariants the caller must honour:
/// - Every field is finite and `>= 0.0`.
/// - At least one field is strictly positive.
///
/// The Python DSL boundary in `experiment.py::_validate_segment_rates`
/// is the only call site that constructs non-default values; it
/// enforces these invariants before the tuple crosses the bridge.
#[derive(Copy, Clone, Debug, PartialEq)]
pub struct SegmentRateWeights {
    pub v: f64,
    pub d: f64,
    pub j: f64,
    pub np: f64,
}

impl SegmentRateWeights {
    /// Flat default — all four buckets at `1.0`. Produces byte-
    /// identical sampling distributions to the pre-slice engine
    /// (both uniform and S5F passes detect this via
    /// [`Self::is_default`] and take the existing fast path).
    pub const fn default() -> Self {
        Self {
            v: 1.0,
            d: 1.0,
            j: 1.0,
            np: 1.0,
        }
    }

    /// Construct from the canonical DSL tuple ordering. Caller
    /// must have validated the inputs (no NaN/inf/negative); this
    /// constructor is a thin packaging step.
    pub const fn from_tuple(v: f64, d: f64, j: f64, np: f64) -> Self {
        Self { v, d, j, np }
    }

    /// Whether this rate vector matches the flat default. Used by
    /// both mutation passes to short-circuit the per-position
    /// segment lookup and preserve the pre-slice fast path.
    ///
    /// Implementation note: we compare via exact equality with
    /// `1.0` because the DSL boundary normalises to `1.0` exactly
    /// when the user omits a bucket. A user who passes
    /// `{"V": 1.0}` explicitly hits the same path as omitting
    /// the kwarg entirely.
    #[inline]
    pub fn is_default(&self) -> bool {
        self.v == 1.0 && self.d == 1.0 && self.j == 1.0 && self.np == 1.0
    }

    /// Per-segment rate lookup. Maps the engine's five-variant
    /// `Segment` enum onto the four-bucket DSL model: both `Np1`
    /// and `Np2` return the NP rate.
    #[inline]
    pub fn rate_for(&self, segment: Segment) -> f64 {
        match segment {
            Segment::V => self.v,
            Segment::D => self.d,
            Segment::J => self.j,
            Segment::Np1 | Segment::Np2 => self.np,
        }
    }
}

impl Default for SegmentRateWeights {
    fn default() -> Self {
        Self::default()
    }
}

/// Locate the [`Segment`] a position belongs to by scanning the
/// assembled sequence's regions in biological order. Returns
/// `None` for positions outside every region (e.g. pre-assembly,
/// post-end-loss).
///
/// Linear over the region count (typically ≤ 5 for a VDJ pool)
/// — the audit's pre-flight check rules out a more elaborate
/// lookup data structure for this slice. The S5F pass's profile-
/// building loop calls this once per non-zero-mutability position;
/// the uniform pass's site-selection loop calls it once per pool
/// position (only when `segment_rates` is non-default).
#[inline]
pub fn segment_at_position(
    sequence: &crate::ir::Sequence,
    pos: u32,
) -> Option<Segment> {
    for region in &sequence.regions {
        if region.start.index() <= pos && pos < region.end.index() {
            return Some(region.segment);
        }
    }
    None
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn default_is_all_ones() {
        let w = SegmentRateWeights::default();
        assert_eq!(w.v, 1.0);
        assert_eq!(w.d, 1.0);
        assert_eq!(w.j, 1.0);
        assert_eq!(w.np, 1.0);
        assert!(w.is_default());
    }

    #[test]
    fn rate_for_maps_np_to_both_np_variants() {
        let w = SegmentRateWeights::from_tuple(0.5, 0.25, 0.75, 0.0);
        assert_eq!(w.rate_for(Segment::V), 0.5);
        assert_eq!(w.rate_for(Segment::D), 0.25);
        assert_eq!(w.rate_for(Segment::J), 0.75);
        assert_eq!(w.rate_for(Segment::Np1), 0.0);
        assert_eq!(w.rate_for(Segment::Np2), 0.0);
    }

    #[test]
    fn is_default_false_for_any_non_one() {
        assert!(!SegmentRateWeights::from_tuple(0.5, 1.0, 1.0, 1.0).is_default());
        assert!(!SegmentRateWeights::from_tuple(1.0, 0.5, 1.0, 1.0).is_default());
        assert!(!SegmentRateWeights::from_tuple(1.0, 1.0, 0.5, 1.0).is_default());
        assert!(!SegmentRateWeights::from_tuple(1.0, 1.0, 1.0, 0.5).is_default());
    }
}
