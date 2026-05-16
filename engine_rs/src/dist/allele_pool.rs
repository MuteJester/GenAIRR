//! `AllelePoolDist` — categorical distribution over `AlleleId`s in a pool.

use super::Distribution;
use crate::refdata::{AlleleId, AllelePool};
use crate::rng::Rng;

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
