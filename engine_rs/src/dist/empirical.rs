//! `EmpiricalLengthDist` — categorical histogram over integer values.

use super::Distribution;
use crate::rng::Rng;

/// Empirical distribution over integer values, weighted by an
/// arbitrary positive-real histogram.
///
/// Used for NP lengths, trim amounts, indel counts — anything where
/// a real-world distribution from the literature or per-allele
/// empirical data needs to be reproduced. Empirical per-allele trim
/// distributions are layered on top of this primitive.
///
/// **Sampling:** inverse-CDF via binary search over cumulative
/// weights. `O(log N)` per draw where `N` is the number of distinct
/// values. Determinism is preserved exactly — `next_f64` consumes
/// one RNG word per sample regardless of histogram shape.
///
/// **Construction discipline:**
/// - Values may be negative (rare but allowed).
/// - Weights must be strictly positive (a zero weight is rejected;
///   if a value should never be drawn, omit it).
/// - The histogram must contain at least one entry.
/// - Duplicate values are not coalesced — the constructor accepts
///   them as-is so callers can choose semantics. Drawing from
///   `[(0, 1.0), (0, 2.0)]` produces 0 with probability 1.0
///   regardless of which row was hit; behaviorally indistinguish-
///   able from `[(0, 3.0)]`.
#[derive(Clone, Debug)]
pub struct EmpiricalLengthDist {
    values: Vec<i64>,
    /// Cumulative sums of the input weights, in the same order as
    /// `values`. `cumulative[i]` is the upper bound of the interval
    /// for `values[i]` on the unit-line `[0, total)`.
    cumulative: Vec<f64>,
    total: f64,
}

impl EmpiricalLengthDist {
    /// Construct from an iterator of `(value, weight)` pairs.
    ///
    /// Panics if:
    /// - the iterator is empty (no valid sample possible),
    /// - any weight is non-positive (NaN, infinite, ≤ 0), or
    /// - the total weight is non-positive after summation.
    pub fn from_pairs<I>(pairs: I) -> Self
    where
        I: IntoIterator<Item = (i64, f64)>,
    {
        let pairs: Vec<(i64, f64)> = pairs.into_iter().collect();
        assert!(
            !pairs.is_empty(),
            "EmpiricalLengthDist: histogram must have at least one entry"
        );

        let mut values = Vec::with_capacity(pairs.len());
        let mut cumulative = Vec::with_capacity(pairs.len());
        let mut running = 0.0f64;

        for (i, (v, w)) in pairs.into_iter().enumerate() {
            assert!(
                w > 0.0 && w.is_finite(),
                "EmpiricalLengthDist: weight at index {} must be a finite \
                 positive number, got {}",
                i,
                w
            );
            running += w;
            values.push(v);
            cumulative.push(running);
        }

        assert!(
            running > 0.0 && running.is_finite(),
            "EmpiricalLengthDist: total weight must be a finite positive \
             number, got {}",
            running
        );

        Self {
            values,
            cumulative,
            total: running,
        }
    }

    /// Construct from `Vec<i64>` (values) and `Vec<f64>` (weights).
    /// The two vectors must have the same length. Same panic
    /// conditions as `from_pairs`.
    pub fn from_values_and_weights(values: Vec<i64>, weights: Vec<f64>) -> Self {
        assert_eq!(
            values.len(),
            weights.len(),
            "EmpiricalLengthDist: values and weights must have equal lengths \
             (got {} and {})",
            values.len(),
            weights.len()
        );
        Self::from_pairs(values.into_iter().zip(weights))
    }

    /// Number of distinct entries in the histogram.
    pub fn len(&self) -> usize {
        self.values.len()
    }

    /// Whether the histogram is empty. Always `false` for any
    /// successfully-constructed `EmpiricalLengthDist` (the
    /// constructor panics on empty input); kept for API symmetry.
    pub fn is_empty(&self) -> bool {
        self.values.is_empty()
    }

    /// Total (unnormalized) weight of the histogram.
    pub fn total_weight(&self) -> f64 {
        self.total
    }

    /// Read-only slice of the value column. Iteration order matches
    /// the order passed to the constructor.
    pub fn values(&self) -> &[i64] {
        &self.values
    }
}

impl Distribution for EmpiricalLengthDist {
    type Output = i64;

    fn sample(&self, rng: &mut Rng) -> i64 {
        // Draw a uniform on [0, total) and binary-search for the
        // smallest cumulative bound that strictly exceeds it. This
        // is the inverse-CDF method for a categorical distribution.
        //
        // `partition_point(|&c| c <= r)` returns the first index
        // where `c > r` — exactly the value we want.
        let r = rng.next_f64() * self.total;
        let idx = self.cumulative.partition_point(|&c| c <= r);
        // Defensive clamp: float rounding could nudge `r` above the
        // last cumulative bound by 1 ULP. Saturate to the last index.
        let idx = idx.min(self.values.len() - 1);
        self.values[idx]
    }

    fn support(&self) -> Option<Vec<(i64, f64)>> {
        // Reconstruct (value, weight) pairs from the parallel arrays.
        // Weights are the per-bucket pre-cumulative weights.
        let mut pairs = Vec::with_capacity(self.values.len());
        let mut prev_cum = 0.0;
        for (i, &v) in self.values.iter().enumerate() {
            let w = self.cumulative[i] - prev_cum;
            pairs.push((v, w));
            prev_cum = self.cumulative[i];
        }
        Some(pairs)
    }
}
