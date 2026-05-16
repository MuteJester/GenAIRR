//! Distributions — typed sources of random samples.
//!
//! ## Design (D4)
//!
//! Distributions are heap-allocated trait objects:
//! `Box<dyn Distribution<Output = T>>`. The associated `Output` type
//! lets each call site keep its return type — an NP length distribution
//! returns `i64`, a base distribution returns `u8`, an allele-pool
//! distribution returns an `AlleleId`. Trait-object dispatch costs
//! one indirection per `sample` call (~5–10 ns), dwarfed by the
//! per-pass overhead from D2 — acceptable.
//!
//! ## Scope
//!
//! The trait + two concrete implementations:
//!
//! - [`UniformBase`] — uniform over `{b'A', b'C', b'G', b'T'}`. Used
//!   by sampling passes when no empirical model is configured.
//! - [`UniformInt`] — uniform over a half-open integer range
//!   `[min, max)`. Foundation for trim / NP-length sampling.

use crate::refdata::{AlleleId, AllelePool};
use crate::rng::Rng;

// ──────────────────────────────────────────────────────────────────
// Trait
// ──────────────────────────────────────────────────────────────────

/// A typed source of random samples.
///
/// Implementors expose `Output = T` so call sites preserve the natural
/// return type. The `Box<dyn Distribution<Output = T>>` form is
/// dyn-compatible because `Output` is fixed at the box site.
pub trait Distribution {
    /// The natural output type of this distribution.
    type Output;

    /// Draw one sample from `rng`. Implementors must consume a
    /// deterministic number of RNG words per call so reproducibility
    /// (D11) holds — same seed + same call sequence ⇒ same output.
    fn sample(&self, rng: &mut Rng) -> Self::Output;

    /// Optional: enumerate the discrete support of this distribution
    /// as `(value, weight)` pairs.
    ///
    /// Used by **constraint-aware sampling** (D.6): a sampling pass
    /// asks the contract set which candidate values are admissible,
    /// filters the support to those, and renormalizes for sampling.
    /// This is what makes "no retries, just sample admissible
    /// values" work — filter-before-sample rather than
    /// retry-and-reject.
    ///
    /// `None` (the default) means the support is too large or
    /// continuous to enumerate practically. Sampling passes that
    /// receive `None` fall back to unconstrained sampling — which
    /// means contracts can't filter that draw. Concrete categorical
    /// distributions (`EmpiricalLengthDist`, `UniformBase`,
    /// `UniformInt` for small ranges, `AllelePoolDist`) override
    /// to return their full support.
    ///
    /// Weights need not be normalized; the caller renormalizes.
    fn support(&self) -> Option<Vec<(Self::Output, f64)>> {
        None
    }
}

// ──────────────────────────────────────────────────────────────────
// sample_filtered — draw from a distribution restricted to a predicate
// ──────────────────────────────────────────────────────────────────

/// Why a filtered sample could not be produced.
#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub enum FilteredSampleError {
    /// The distribution does not expose enumerable support.
    SupportUnavailable,
    /// Enumerable support existed, but every candidate was rejected.
    EmptyAdmissibleSupport,
    /// Filtered support had invalid weights.
    InvalidFilteredSupport,
}

/// Draw one sample from `dist` restricted to values satisfying
/// `predicate`, preserving the reason if no filtered sample exists.
///
/// Returns `Err` when:
/// - the distribution doesn't expose a `support()` (continuous /
///   too-large categorical),
/// - the filtered support is empty (no admissible candidates), or
/// - the filtered support weights are invalid.
///
/// **The architectural primitive for constraint-aware sampling.**
/// Sampling passes (D.6) call this with `predicate` wrapping a
/// `ContractSet::admits` check. The pass receives `Ok(value)` if
/// any admissible value exists, or a structured error explaining why
/// a strict runtime cannot sample under the requested constraints.
///
/// Determinism: consumes exactly one RNG word only when a filtered
/// sample can be drawn. Failure paths consume no RNG words.
pub fn sample_filtered_result<T, D, F>(
    rng: &mut Rng,
    dist: &D,
    predicate: F,
) -> Result<T, FilteredSampleError>
where
    T: Clone,
    D: Distribution<Output = T> + ?Sized,
    F: Fn(&T) -> bool,
{
    let support = dist
        .support()
        .ok_or(FilteredSampleError::SupportUnavailable)?;
    let mut filtered: Vec<(T, f64)> = support.into_iter().filter(|(v, _)| predicate(v)).collect();

    if filtered.is_empty() {
        return Err(FilteredSampleError::EmptyAdmissibleSupport);
    }

    // Renormalize and sample by inverse CDF.
    let total: f64 = filtered.iter().map(|(_, w)| w).sum();
    if !total.is_finite() || total <= 0.0 {
        return Err(FilteredSampleError::InvalidFilteredSupport);
    }
    let r = rng.next_f64() * total;
    let mut cum = 0.0;
    for (v, w) in filtered.iter() {
        cum += w;
        if r < cum {
            return Ok(v.clone());
        }
    }
    // Defensive last-index fallback for ULP rounding.
    Ok(filtered.pop().unwrap().0)
}

/// Permissive filtered sampling helper.
///
/// Callers get a filtered value when possible and `None` when
/// filtering cannot be performed. Strict runtime paths should call
/// [`sample_filtered_result`] so they can surface structured errors
/// instead of silently falling back.
pub fn sample_filtered<T, D, F>(rng: &mut Rng, dist: &D, predicate: F) -> Option<T>
where
    T: Clone,
    D: Distribution<Output = T> + ?Sized,
    F: Fn(&T) -> bool,
{
    sample_filtered_result(rng, dist, predicate).ok()
}

// ──────────────────────────────────────────────────────────────────
// UniformBase — uniform over A, C, G, T
// ──────────────────────────────────────────────────────────────────

/// Uniform distribution over the four canonical DNA bases
/// `{b'A', b'C', b'G', b'T'}` (uppercase). Each base has probability
/// 0.25.
///
/// This is the default fallback when no empirical TdT / NP-base
/// model is configured. Empirical transition matrices override at
/// runtime through the same trait surface.
#[derive(Clone, Debug, Default)]
pub struct UniformBase;

impl Distribution for UniformBase {
    type Output = u8;

    fn sample(&self, rng: &mut Rng) -> u8 {
        const BASES: [u8; 4] = [b'A', b'C', b'G', b'T'];
        BASES[rng.range_u32(4) as usize]
    }

    fn support(&self) -> Option<Vec<(u8, f64)>> {
        Some(vec![(b'A', 1.0), (b'C', 1.0), (b'G', 1.0), (b'T', 1.0)])
    }
}

// ──────────────────────────────────────────────────────────────────
// UniformInt — uniform over [min, max)
// ──────────────────────────────────────────────────────────────────

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

// ──────────────────────────────────────────────────────────────────
// EmpiricalLengthDist — categorical histogram over integer values
// ──────────────────────────────────────────────────────────────────

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

// ──────────────────────────────────────────────────────────────────
// AllelePoolDist — categorical distribution over AlleleIds in a pool
// ──────────────────────────────────────────────────────────────────

/// Distribution that samples an `AlleleId` from a specific
/// `AllelePool`, weighted by per-allele frequencies.
///
/// Construction binds the distribution to the pool *size* at
/// build-time, so every sampled `AlleleId` is guaranteed to be in
/// bounds for that pool. This is the discipline that lets future
/// passes call `pool.get(dist.sample(rng)).unwrap()` without
/// runtime fallibility — the unwrap is structurally safe.
///
/// Internally identical to `EmpiricalLengthDist` (cumulative
/// weights + binary search) but typed to return `AlleleId` and
/// constructed from a pool reference rather than free-form pairs.
#[derive(Clone, Debug)]
pub struct AllelePoolDist {
    cumulative: Vec<f64>,
    /// `Some(ids)` when the distribution covers a strict subset of
    /// the pool: cumulative bucket `i` corresponds to `ids[i]`.
    /// `None` when the distribution covers the entire pool with
    /// identity index→ID mapping (the default for `uniform` and
    /// `from_weights`).
    ids: Option<Vec<AlleleId>>,
    total: f64,
}

impl AllelePoolDist {
    /// Uniform distribution over all alleles in `pool`. Panics if
    /// the pool is empty (no valid sample possible).
    pub fn uniform(pool: &AllelePool) -> Self {
        assert!(
            !pool.is_empty(),
            "AllelePoolDist::uniform: pool must contain at least one allele"
        );
        Self::from_weights(pool, vec![1.0; pool.len()])
    }

    /// Uniform distribution restricted to the alleles named in
    /// `allowed_ids`. The resulting distribution samples each listed
    /// allele with equal probability and ignores the rest of the
    /// pool. Used by the `Experiment.using(...)` allele-lock API.
    ///
    /// Panics if:
    /// - `allowed_ids` is empty,
    /// - any ID is out of range for `pool` (i.e. `id.index() >= pool.len()`),
    /// - `allowed_ids` contains duplicates.
    pub fn restricted_uniform(pool: &AllelePool, allowed_ids: Vec<AlleleId>) -> Self {
        assert!(
            !allowed_ids.is_empty(),
            "AllelePoolDist::restricted_uniform: allowed_ids must be non-empty"
        );
        let pool_len = pool.len() as u32;
        let mut seen = std::collections::HashSet::with_capacity(allowed_ids.len());
        for id in &allowed_ids {
            assert!(
                id.index() < pool_len,
                "AllelePoolDist::restricted_uniform: AlleleId({}) is out of range \
                 for pool of size {}",
                id.index(),
                pool_len
            );
            assert!(
                seen.insert(*id),
                "AllelePoolDist::restricted_uniform: duplicate AlleleId({}) in \
                 allowed_ids",
                id.index()
            );
        }

        let n = allowed_ids.len();
        let mut cumulative = Vec::with_capacity(n);
        for i in 0..n {
            cumulative.push((i + 1) as f64);
        }
        Self {
            cumulative,
            ids: Some(allowed_ids),
            total: n as f64,
        }
    }

    /// Empirical distribution over `pool` with one weight per
    /// allele. `weights[i]` is the unnormalized frequency of the
    /// allele at `AlleleId(i)`.
    ///
    /// Panics if:
    /// - `weights.len() != pool.len()` (size mismatch),
    /// - any weight is non-positive / NaN / infinite, or
    /// - the total weight is non-positive after summation.
    pub fn from_weights(pool: &AllelePool, weights: Vec<f64>) -> Self {
        assert_eq!(
            weights.len(),
            pool.len(),
            "AllelePoolDist::from_weights: weights ({}) must match pool size ({})",
            weights.len(),
            pool.len()
        );
        assert!(
            !weights.is_empty(),
            "AllelePoolDist::from_weights: pool must contain at least one allele"
        );

        let mut cumulative = Vec::with_capacity(weights.len());
        let mut running = 0.0f64;
        for (i, w) in weights.into_iter().enumerate() {
            assert!(
                w > 0.0 && w.is_finite(),
                "AllelePoolDist::from_weights: weight at index {} must be a \
                 finite positive number, got {}",
                i,
                w
            );
            running += w;
            cumulative.push(running);
        }

        assert!(
            running > 0.0 && running.is_finite(),
            "AllelePoolDist::from_weights: total weight must be a finite \
             positive number, got {}",
            running
        );

        Self {
            cumulative,
            ids: None,
            total: running,
        }
    }

    /// Number of alleles in the distribution's support (equal to
    /// the pool size at construction time).
    pub fn len(&self) -> usize {
        self.cumulative.len()
    }

    /// Whether the distribution is empty. Always `false` for any
    /// successfully-constructed `AllelePoolDist`.
    pub fn is_empty(&self) -> bool {
        self.cumulative.is_empty()
    }

    /// Total (unnormalized) weight.
    pub fn total_weight(&self) -> f64 {
        self.total
    }
}

impl Distribution for AllelePoolDist {
    type Output = AlleleId;

    fn sample(&self, rng: &mut Rng) -> AlleleId {
        let r = rng.next_f64() * self.total;
        let idx = self.cumulative.partition_point(|&c| c <= r);
        // Same defensive last-index clamp as EmpiricalLengthDist.
        let idx = idx.min(self.cumulative.len() - 1);
        match &self.ids {
            Some(ids) => ids[idx],
            None => AlleleId::new(idx as u32),
        }
    }

    fn support(&self) -> Option<Vec<(AlleleId, f64)>> {
        let mut pairs = Vec::with_capacity(self.cumulative.len());
        let mut prev_cum = 0.0;
        for (i, &cum) in self.cumulative.iter().enumerate() {
            let w = cum - prev_cum;
            let id = match &self.ids {
                Some(ids) => ids[i],
                None => AlleleId::new(i as u32),
            };
            pairs.push((id, w));
            prev_cum = cum;
        }
        Some(pairs)
    }
}

// ──────────────────────────────────────────────────────────────────
// Tests
// ──────────────────────────────────────────────────────────────────


#[cfg(test)]
#[path = "dist_tests.rs"]
mod tests;
