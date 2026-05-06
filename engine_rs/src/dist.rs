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
//! per-pass overhead from D2 — acceptable at this phase.
//!
//! ## Phase B.3 scope
//!
//! Just the trait + two concrete implementations:
//!
//! - [`UniformBase`] — uniform over `{b'A', b'C', b'G', b'T'}`. Used
//!   by Phase B.5's first sampling pass and by future NP-base passes
//!   when no empirical model is configured.
//! - [`UniformInt`] — uniform over a half-open integer range
//!   `[min, max)`. Foundation for trim / NP-length sampling once
//!   real empirical distributions land in Phase E.
//!
//! `filter()` and `support_summary()` (D4 trait surface) are
//! deliberately omitted here — they get added when contracts arrive
//! in Phase D and have something to filter against.

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
    /// values" work — the architectural pivot from V5's retry-and-
    /// reject to V6's filter-before-sample.
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
/// This preserves the original Phase D.6 behaviour: callers get a
/// filtered value when possible and `None` when filtering cannot be
/// performed. Strict runtime paths should call
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
/// model is configured. Real biology will override with empirical
/// transition matrices in Phase E; the trait surface is unchanged.
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
/// empirical data needs to be reproduced. Phase E will add empirical
/// per-allele trim distributions on top of this primitive.
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
mod tests {
    use super::*;

    // ── UniformBase ────────────────────────────────────────────────

    #[test]
    fn uniform_base_only_emits_canonical_bases() {
        let mut rng = Rng::new(1);
        let dist = UniformBase;
        for _ in 0..1000 {
            let b = dist.sample(&mut rng);
            assert!(
                matches!(b, b'A' | b'C' | b'G' | b'T'),
                "UniformBase emitted unexpected byte 0x{:02x}",
                b
            );
        }
    }

    #[test]
    fn uniform_base_covers_all_four_bases() {
        let mut rng = Rng::new(123);
        let dist = UniformBase;
        let mut seen = [false; 4];
        for _ in 0..1000 {
            let b = dist.sample(&mut rng);
            let idx = match b {
                b'A' => 0,
                b'C' => 1,
                b'G' => 2,
                b'T' => 3,
                _ => panic!("unexpected base"),
            };
            seen[idx] = true;
        }
        for (i, &was_seen) in seen.iter().enumerate() {
            assert!(was_seen, "base index {} never appeared in 1000 draws", i);
        }
    }

    #[test]
    fn uniform_base_distribution_is_roughly_uniform() {
        let mut rng = Rng::new(0xfeed);
        let dist = UniformBase;
        let mut counts = [0u32; 4];
        let n = 10_000u32;
        for _ in 0..n {
            let idx = match dist.sample(&mut rng) {
                b'A' => 0,
                b'C' => 1,
                b'G' => 2,
                b'T' => 3,
                _ => unreachable!(),
            };
            counts[idx] += 1;
        }
        // Expected count per bucket = 2500. Generous tolerance.
        for &c in &counts {
            assert!(
                (2200..=2800).contains(&c),
                "UniformBase bucket count {} outside [2200, 2800]",
                c
            );
        }
    }

    #[test]
    fn uniform_base_same_seed_same_stream() {
        let mut a = Rng::new(7);
        let mut b = Rng::new(7);
        let dist = UniformBase;
        for _ in 0..100 {
            assert_eq!(dist.sample(&mut a), dist.sample(&mut b));
        }
    }

    // ── UniformInt ─────────────────────────────────────────────────

    #[test]
    #[should_panic(expected = "UniformInt: max")]
    fn uniform_int_new_rejects_max_le_min() {
        let _ = UniformInt::new(5, 5);
    }

    #[test]
    #[should_panic(expected = "UniformInt: max")]
    fn uniform_int_new_rejects_inverted_range() {
        let _ = UniformInt::new(10, 5);
    }

    #[test]
    #[should_panic(expected = "UniformInt: span")]
    fn uniform_int_new_rejects_oversized_span() {
        // Span exceeds u32::MAX → must panic explicitly, not silently
        // truncate. Pre-cleanup the constructor accepted this and
        // sample() biased the result; now it's caught at construction.
        let _ = UniformInt::new(0, (u32::MAX as i64) + 2);
    }

    #[test]
    fn uniform_int_new_accepts_max_span() {
        // Right at the boundary — span == u32::MAX is allowed.
        let dist = UniformInt::new(0, u32::MAX as i64);
        assert_eq!(dist.span(), u32::MAX as u64);
    }

    #[test]
    fn uniform_int_stays_in_bounds() {
        let mut rng = Rng::new(42);
        let dist = UniformInt::new(3, 10); // [3, 10) → 3..=9
        for _ in 0..1000 {
            let v = dist.sample(&mut rng);
            assert!(v >= 3, "UniformInt(3,10) produced {} (< 3)", v);
            assert!(v < 10, "UniformInt(3,10) produced {} (>= 10)", v);
        }
    }

    #[test]
    fn uniform_int_covers_full_range() {
        let mut rng = Rng::new(99);
        let dist = UniformInt::new(0, 5);
        let mut seen = [false; 5];
        for _ in 0..1000 {
            let v = dist.sample(&mut rng);
            seen[v as usize] = true;
        }
        for (i, &was_seen) in seen.iter().enumerate() {
            assert!(was_seen, "UniformInt(0,5) never produced {}", i);
        }
    }

    #[test]
    fn uniform_int_negative_min_works() {
        let mut rng = Rng::new(11);
        let dist = UniformInt::new(-3, 4); // [-3, 4) → -3, -2, -1, 0, 1, 2, 3
        for _ in 0..1000 {
            let v = dist.sample(&mut rng);
            assert!(v >= -3, "{}", v);
            assert!(v < 4, "{}", v);
        }
    }

    #[test]
    fn uniform_int_span_one_returns_min() {
        let mut rng = Rng::new(1);
        let dist = UniformInt::new(7, 8);
        for _ in 0..100 {
            assert_eq!(dist.sample(&mut rng), 7);
        }
    }

    #[test]
    fn uniform_int_same_seed_same_stream() {
        let mut a = Rng::new(0xbabe);
        let mut b = Rng::new(0xbabe);
        let dist = UniformInt::new(0, 100);
        for _ in 0..100 {
            assert_eq!(dist.sample(&mut a), dist.sample(&mut b));
        }
    }

    #[test]
    fn uniform_int_accessors_round_trip() {
        let dist = UniformInt::new(-5, 12);
        assert_eq!(dist.min(), -5);
        assert_eq!(dist.max(), 12);
        assert_eq!(dist.span(), 17);
    }

    // ── Trait object usage ─────────────────────────────────────────

    #[test]
    fn box_dyn_distribution_dispatches_correctly() {
        // The whole point of D4: heterogeneous distributions stored
        // as trait objects. This compiles only if the trait is
        // actually dyn-compatible with the chosen Output type.
        let mut rng = Rng::new(1);

        let base_dist: Box<dyn Distribution<Output = u8>> = Box::new(UniformBase);
        let int_dist: Box<dyn Distribution<Output = i64>> = Box::new(UniformInt::new(0, 10));

        let b = base_dist.sample(&mut rng);
        let i = int_dist.sample(&mut rng);

        assert!(matches!(b, b'A' | b'C' | b'G' | b'T'));
        assert!((0..10).contains(&i));
    }

    #[test]
    fn vec_of_homogeneous_dyn_distributions() {
        // Storing many distributions in a heterogeneous container.
        let dists: Vec<Box<dyn Distribution<Output = u8>>> = vec![
            Box::new(UniformBase),
            Box::new(UniformBase),
            Box::new(UniformBase),
        ];
        let mut rng = Rng::new(0);
        for d in &dists {
            let b = d.sample(&mut rng);
            assert!(matches!(b, b'A' | b'C' | b'G' | b'T'));
        }
    }

    #[derive(Clone, Debug)]
    struct NoSupportDist;

    impl Distribution for NoSupportDist {
        type Output = i64;

        fn sample(&self, _rng: &mut Rng) -> i64 {
            0
        }
    }

    #[test]
    fn sample_filtered_result_reports_support_unavailable() {
        let mut rng = Rng::new(1);
        let err = sample_filtered_result(&mut rng, &NoSupportDist, |_| true).unwrap_err();
        assert_eq!(err, FilteredSampleError::SupportUnavailable);
    }

    #[test]
    fn sample_filtered_result_reports_empty_admissible_support() {
        let mut rng = Rng::new(1);
        let dist = UniformInt::new(0, 4);
        let err = sample_filtered_result(&mut rng, &dist, |_| false).unwrap_err();
        assert_eq!(err, FilteredSampleError::EmptyAdmissibleSupport);
    }

    #[test]
    fn sample_filtered_result_samples_from_admissible_subset() {
        let mut rng = Rng::new(1);
        let dist = UniformInt::new(0, 10);
        for _ in 0..100 {
            let value = sample_filtered_result(&mut rng, &dist, |v| *v >= 7).unwrap();
            assert!((7..10).contains(&value));
        }
    }

    #[test]
    fn sample_filtered_permissive_collapses_filter_errors_to_none() {
        let mut rng = Rng::new(1);
        let dist = UniformInt::new(0, 4);
        assert_eq!(sample_filtered(&mut rng, &dist, |_| false), None);
        assert_eq!(sample_filtered(&mut rng, &NoSupportDist, |_| true), None);
    }

    // ── EmpiricalLengthDist ────────────────────────────────────────

    #[test]
    #[should_panic(expected = "must have at least one entry")]
    fn empirical_rejects_empty_histogram() {
        let _ = EmpiricalLengthDist::from_pairs(std::iter::empty());
    }

    #[test]
    #[should_panic(expected = "weight at index")]
    fn empirical_rejects_zero_weight() {
        let _ = EmpiricalLengthDist::from_pairs(vec![(5, 0.0), (6, 1.0)]);
    }

    #[test]
    #[should_panic(expected = "weight at index")]
    fn empirical_rejects_negative_weight() {
        let _ = EmpiricalLengthDist::from_pairs(vec![(5, -1.0)]);
    }

    #[test]
    #[should_panic(expected = "weight at index")]
    fn empirical_rejects_nan_weight() {
        let _ = EmpiricalLengthDist::from_pairs(vec![(5, f64::NAN)]);
    }

    #[test]
    #[should_panic(expected = "weight at index")]
    fn empirical_rejects_infinite_weight() {
        let _ = EmpiricalLengthDist::from_pairs(vec![(5, f64::INFINITY)]);
    }

    #[test]
    fn empirical_single_value_always_returned() {
        let dist = EmpiricalLengthDist::from_pairs(vec![(42, 1.0)]);
        let mut rng = Rng::new(0xc0ff_ee);
        for _ in 0..1000 {
            assert_eq!(dist.sample(&mut rng), 42);
        }
    }

    #[test]
    fn empirical_construction_accessors_round_trip() {
        let dist = EmpiricalLengthDist::from_pairs(vec![(0, 1.0), (1, 2.0), (2, 1.0)]);
        assert_eq!(dist.len(), 3);
        assert!(!dist.is_empty());
        assert_eq!(dist.values(), &[0, 1, 2]);
        assert!((dist.total_weight() - 4.0).abs() < 1e-12);
    }

    #[test]
    fn empirical_from_values_and_weights_matches_from_pairs() {
        let pairs = EmpiricalLengthDist::from_pairs(vec![(0, 1.0), (1, 2.0), (2, 3.0)]);
        let split =
            EmpiricalLengthDist::from_values_and_weights(vec![0, 1, 2], vec![1.0, 2.0, 3.0]);

        let mut a = Rng::new(7);
        let mut b = Rng::new(7);
        for _ in 0..100 {
            assert_eq!(pairs.sample(&mut a), split.sample(&mut b));
        }
    }

    #[test]
    #[should_panic(expected = "values and weights must have equal lengths")]
    fn empirical_from_values_and_weights_rejects_length_mismatch() {
        let _ = EmpiricalLengthDist::from_values_and_weights(vec![0, 1], vec![1.0]);
    }

    #[test]
    fn empirical_stays_in_value_set() {
        let dist = EmpiricalLengthDist::from_pairs(vec![(3, 1.0), (7, 2.0), (11, 1.0), (-2, 0.5)]);
        let mut rng = Rng::new(17);
        for _ in 0..1000 {
            let v = dist.sample(&mut rng);
            assert!(matches!(v, 3 | 7 | 11 | -2));
        }
    }

    #[test]
    fn empirical_covers_full_value_set() {
        let dist = EmpiricalLengthDist::from_pairs(vec![(0, 1.0), (1, 1.0), (2, 1.0), (3, 1.0)]);
        let mut rng = Rng::new(99);
        let mut seen = [false; 4];
        for _ in 0..1000 {
            let v = dist.sample(&mut rng);
            seen[v as usize] = true;
        }
        for (i, &was_seen) in seen.iter().enumerate() {
            assert!(was_seen, "value {} never appeared", i);
        }
    }

    #[test]
    fn empirical_weights_are_respected() {
        // 90/10 weight split — over 10,000 draws the heavy value
        // should appear roughly 9,000 times. Generous tolerance
        // catches statistical jitter.
        let dist = EmpiricalLengthDist::from_pairs(vec![(0, 0.9), (1, 0.1)]);
        let mut rng = Rng::new(0xabba);
        let mut zero_count = 0;
        let n = 10_000;
        for _ in 0..n {
            if dist.sample(&mut rng) == 0 {
                zero_count += 1;
            }
        }
        assert!(
            (8500..=9500).contains(&zero_count),
            "expected zero_count ~9000 of 10000, got {}",
            zero_count
        );
    }

    #[test]
    fn empirical_uniform_weights_produce_roughly_uniform_distribution() {
        let dist =
            EmpiricalLengthDist::from_pairs(vec![(0, 1.0), (1, 1.0), (2, 1.0), (3, 1.0), (4, 1.0)]);
        let mut rng = Rng::new(0xdead);
        let mut counts = [0u32; 5];
        let n = 10_000u32;
        for _ in 0..n {
            counts[dist.sample(&mut rng) as usize] += 1;
        }
        // Expected ~2000 per bucket. Generous tolerance.
        for (i, &c) in counts.iter().enumerate() {
            assert!(
                (1700..=2300).contains(&c),
                "bucket {} count {} outside [1700, 2300]",
                i,
                c
            );
        }
    }

    #[test]
    fn empirical_same_seed_same_stream() {
        let dist = EmpiricalLengthDist::from_pairs(vec![(10, 1.0), (20, 2.0), (30, 3.0)]);
        let mut a = Rng::new(0xfeed);
        let mut b = Rng::new(0xfeed);
        for _ in 0..100 {
            assert_eq!(dist.sample(&mut a), dist.sample(&mut b));
        }
    }

    #[test]
    fn empirical_negative_values_supported() {
        let dist = EmpiricalLengthDist::from_pairs(vec![(-5, 1.0), (0, 1.0), (5, 1.0)]);
        let mut rng = Rng::new(13);
        for _ in 0..100 {
            let v = dist.sample(&mut rng);
            assert!(matches!(v, -5 | 0 | 5));
        }
    }

    #[test]
    fn empirical_works_through_box_dyn() {
        let dist: Box<dyn Distribution<Output = i64>> =
            Box::new(EmpiricalLengthDist::from_pairs(vec![(7, 1.0)]));
        let mut rng = Rng::new(0);
        assert_eq!(dist.sample(&mut rng), 7);
    }

    // ── AllelePoolDist ─────────────────────────────────────────────

    use crate::ir::Segment;
    use crate::refdata::Allele;

    /// Build a pool of `n` named alleles for testing. The base
    /// sequences are all single-byte `b'A'` since we don't care
    /// about content here, only sampling behavior.
    fn make_pool(n: usize) -> AllelePool {
        let mut p = AllelePool::new();
        for i in 0..n {
            let _ = p.push(Allele {
                name: format!("a{}*01", i),
                gene: format!("a{}", i),
                seq: b"A".to_vec(),
                segment: Segment::V,
                anchor: None,
            });
        }
        p
    }

    #[test]
    #[should_panic(expected = "pool must contain at least one allele")]
    fn allele_pool_dist_uniform_rejects_empty_pool() {
        let p = AllelePool::new();
        let _ = AllelePoolDist::uniform(&p);
    }

    #[test]
    #[should_panic(expected = "weights")]
    fn allele_pool_dist_from_weights_rejects_size_mismatch() {
        let p = make_pool(3);
        let _ = AllelePoolDist::from_weights(&p, vec![1.0, 2.0]);
    }

    #[test]
    #[should_panic(expected = "weight at index")]
    fn allele_pool_dist_from_weights_rejects_zero_weight() {
        let p = make_pool(2);
        let _ = AllelePoolDist::from_weights(&p, vec![0.0, 1.0]);
    }

    #[test]
    #[should_panic(expected = "weight at index")]
    fn allele_pool_dist_from_weights_rejects_nan() {
        let p = make_pool(1);
        let _ = AllelePoolDist::from_weights(&p, vec![f64::NAN]);
    }

    #[test]
    fn allele_pool_dist_uniform_construction_round_trip() {
        let p = make_pool(5);
        let dist = AllelePoolDist::uniform(&p);
        assert_eq!(dist.len(), 5);
        assert!(!dist.is_empty());
        assert!((dist.total_weight() - 5.0).abs() < 1e-12);
    }

    #[test]
    fn allele_pool_dist_uniform_covers_all_alleles() {
        let p = make_pool(10);
        let dist = AllelePoolDist::uniform(&p);
        let mut rng = Rng::new(99);
        let mut seen = vec![false; 10];
        for _ in 0..2000 {
            let id = dist.sample(&mut rng);
            assert!(id.as_usize() < 10, "out-of-bounds AlleleId");
            seen[id.as_usize()] = true;
        }
        for (i, &was_seen) in seen.iter().enumerate() {
            assert!(was_seen, "AlleleId({}) never sampled", i);
        }
    }

    #[test]
    fn allele_pool_dist_uniform_is_roughly_uniform() {
        let p = make_pool(4);
        let dist = AllelePoolDist::uniform(&p);
        let mut rng = Rng::new(0xfade);
        let mut counts = [0u32; 4];
        let n = 10_000;
        for _ in 0..n {
            counts[dist.sample(&mut rng).as_usize()] += 1;
        }
        // Expected ~2500 per bucket.
        for (i, &c) in counts.iter().enumerate() {
            assert!(
                (2200..=2800).contains(&c),
                "bucket {} count {} outside [2200, 2800]",
                i,
                c
            );
        }
    }

    #[test]
    fn allele_pool_dist_weights_are_respected() {
        // 80/10/10 weighted pool. Heavy allele should dominate.
        let p = make_pool(3);
        let dist = AllelePoolDist::from_weights(&p, vec![0.8, 0.1, 0.1]);
        let mut rng = Rng::new(0xfeed);
        let mut counts = [0u32; 3];
        let n = 10_000;
        for _ in 0..n {
            counts[dist.sample(&mut rng).as_usize()] += 1;
        }
        // Allele 0 ~80% (8000); alleles 1 and 2 ~10% each (1000).
        assert!(
            (7500..=8500).contains(&counts[0]),
            "heavy bucket count {} outside [7500, 8500]",
            counts[0]
        );
        assert!(
            (700..=1300).contains(&counts[1]) && (700..=1300).contains(&counts[2]),
            "light bucket counts {:?} outside [700, 1300]",
            &counts[1..]
        );
    }

    #[test]
    fn allele_pool_dist_same_seed_same_stream() {
        let p = make_pool(7);
        let dist = AllelePoolDist::uniform(&p);
        let mut a = Rng::new(0xcafe);
        let mut b = Rng::new(0xcafe);
        for _ in 0..100 {
            assert_eq!(dist.sample(&mut a), dist.sample(&mut b));
        }
    }

    #[test]
    fn allele_pool_dist_single_allele_always_returned() {
        let p = make_pool(1);
        let dist = AllelePoolDist::uniform(&p);
        let mut rng = Rng::new(1);
        for _ in 0..100 {
            assert_eq!(dist.sample(&mut rng), AlleleId::new(0));
        }
    }

    #[test]
    fn allele_pool_dist_works_through_box_dyn() {
        let p = make_pool(3);
        let dist: Box<dyn Distribution<Output = AlleleId>> = Box::new(AllelePoolDist::uniform(&p));
        let mut rng = Rng::new(0);
        let id = dist.sample(&mut rng);
        assert!(id.as_usize() < 3);
    }

    #[test]
    fn allele_pool_dist_sampled_ids_resolve_in_pool() {
        // Integration sanity: sampled AlleleIds round-trip through
        // the pool. This is the structural guarantee the dist gives:
        // sample → get → Some(allele).
        let p = make_pool(5);
        let dist = AllelePoolDist::uniform(&p);
        let mut rng = Rng::new(42);
        for _ in 0..50 {
            let id = dist.sample(&mut rng);
            let resolved = p.get(id);
            assert!(resolved.is_some(), "sampled id {:?} did not resolve", id);
        }
    }

    // ── restricted_uniform tests ────────────────────────────────────

    #[test]
    fn allele_pool_dist_restricted_uniform_only_samples_allowed_ids() {
        let p = make_pool(10);
        let allowed = vec![AlleleId::new(2), AlleleId::new(5), AlleleId::new(7)];
        let dist = AllelePoolDist::restricted_uniform(&p, allowed.clone());
        assert_eq!(dist.len(), 3);
        let mut rng = Rng::new(0xbeef);
        for _ in 0..500 {
            let id = dist.sample(&mut rng);
            assert!(
                allowed.contains(&id),
                "sampled id {:?} not in allowed set",
                id
            );
        }
    }

    #[test]
    fn allele_pool_dist_restricted_uniform_covers_all_allowed() {
        let p = make_pool(8);
        let allowed = vec![AlleleId::new(1), AlleleId::new(4), AlleleId::new(6)];
        let dist = AllelePoolDist::restricted_uniform(&p, allowed.clone());
        let mut rng = Rng::new(0xc0de);
        let mut seen = vec![false; allowed.len()];
        for _ in 0..2000 {
            let id = dist.sample(&mut rng);
            let pos = allowed.iter().position(|&a| a == id).unwrap();
            seen[pos] = true;
        }
        assert!(seen.iter().all(|&s| s), "not every allowed id was sampled");
    }

    #[test]
    fn allele_pool_dist_restricted_uniform_single_id_always_returns_it() {
        let p = make_pool(20);
        let dist = AllelePoolDist::restricted_uniform(&p, vec![AlleleId::new(13)]);
        let mut rng = Rng::new(7);
        for _ in 0..50 {
            assert_eq!(dist.sample(&mut rng), AlleleId::new(13));
        }
    }

    #[test]
    fn allele_pool_dist_restricted_uniform_support_round_trip() {
        let p = make_pool(6);
        let allowed = vec![AlleleId::new(0), AlleleId::new(3), AlleleId::new(5)];
        let dist = AllelePoolDist::restricted_uniform(&p, allowed.clone());
        let support = dist.support().expect("support is Some");
        let support_ids: Vec<AlleleId> = support.iter().map(|(id, _)| *id).collect();
        assert_eq!(support_ids, allowed);
        for (_, weight) in &support {
            assert!((*weight - 1.0).abs() < 1e-12);
        }
    }

    #[test]
    #[should_panic(expected = "allowed_ids must be non-empty")]
    fn allele_pool_dist_restricted_uniform_rejects_empty() {
        let p = make_pool(3);
        let _ = AllelePoolDist::restricted_uniform(&p, vec![]);
    }

    #[test]
    #[should_panic(expected = "out of range")]
    fn allele_pool_dist_restricted_uniform_rejects_out_of_range_id() {
        let p = make_pool(3);
        let _ = AllelePoolDist::restricted_uniform(&p, vec![AlleleId::new(5)]);
    }

    #[test]
    #[should_panic(expected = "duplicate")]
    fn allele_pool_dist_restricted_uniform_rejects_duplicate_ids() {
        let p = make_pool(3);
        let _ = AllelePoolDist::restricted_uniform(
            &p,
            vec![AlleleId::new(1), AlleleId::new(1)],
        );
    }
}
