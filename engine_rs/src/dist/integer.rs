//! `UniformInt` — uniform distribution over a half-open integer range `[min, max)`.

use super::Distribution;
use crate::rng::Rng;

/// Uniform distribution over the half-open integer range `[min, max)`.
///
/// `min` and `max` are inclusive / exclusive bounds respectively.
/// Returns `i64` to match `ChoiceValue::Int`, the trace
/// representation for integer choices.
///
/// **Bounds:** the underlying RNG primitive operates on `u32`, so the
/// total span `max - min` is bounded to `u32::MAX`. The constructor
/// rejects spans exceeding this — better an explicit panic at
/// `UniformInt::new` than a silently biased sample later. The domain
/// (NP lengths ≤ 50, trim amounts ≤ 30, indel counts ≤ 10) is far
/// beneath the cap; if a future use needs larger ranges, the right
/// fix is to widen the underlying primitive (`range_u64` based on
/// `next_u64`), not to silently truncate here.
#[derive(Clone, Debug)]
pub struct UniformInt {
    min: i64,
    max: i64,      // exclusive
    span_u32: u32, // cached for sampling — validated at construction
}

impl UniformInt {
    /// Construct a uniform integer distribution over `[min, max)`.
    ///
    /// Panics if either:
    /// - `max <= min` (degenerate range — caller bug), or
    /// - `max - min > u32::MAX` (range exceeds the underlying RNG
    ///   primitive's capacity — see struct doc).
    pub fn new(min: i64, max: i64) -> Self {
        assert!(
            max > min,
            "UniformInt: max ({}) must be strictly greater than min ({})",
            max,
            min
        );
        let span = (max - min) as u64;
        assert!(
            span <= u32::MAX as u64,
            "UniformInt: span ({}) exceeds u32::MAX. Use a smaller range \
             or widen the RNG primitive — silently truncating would \
             produce a biased sample.",
            span
        );
        Self {
            min,
            max,
            span_u32: span as u32,
        }
    }

    /// Inclusive lower bound.
    pub fn min(&self) -> i64 {
        self.min
    }

    /// Exclusive upper bound.
    pub fn max(&self) -> i64 {
        self.max
    }

    /// Number of distinct integer values in this distribution's
    /// support (`max - min`).
    pub fn span(&self) -> u64 {
        self.span_u32 as u64
    }
}

impl Distribution for UniformInt {
    type Output = i64;

    fn sample(&self, rng: &mut Rng) -> i64 {
        // span_u32 is validated to be > 0 (max > min) and <= u32::MAX
        // at construction time, so range_u32 receives a valid argument
        // and no truncation can occur.
        self.min + rng.range_u32(self.span_u32) as i64
    }

    fn support(&self) -> Option<Vec<(i64, f64)>> {
        // Enumerate [min, max) with uniform weights. For very large
        // ranges (>1024) we return None to avoid blowing up memory;
        // the simulation domain only ever uses small ranges so this
        // cap is forgiving for actual use cases.
        const MAX_ENUMERATED: u32 = 1024;
        if self.span_u32 > MAX_ENUMERATED {
            return None;
        }
        let mut pairs = Vec::with_capacity(self.span_u32 as usize);
        for i in 0..self.span_u32 {
            pairs.push((self.min + i as i64, 1.0));
        }
        Some(pairs)
    }
}
