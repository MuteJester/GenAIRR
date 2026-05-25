//! `MutationCountSource` — the count abstraction used by
//! [`super::UniformMutationPass`] and [`super::S5FMutationPass`].
//!
//! Two variants:
//!
//! - **`Distribution`** — the classic shape: an explicit
//!   `Distribution<Output = i64>` (typically empirical or a fixed
//!   value) sampled once per pass execution. The pool length is
//!   ignored. This is the path the v1 `mutate(count=...)` DSL took
//!   and the only path before v2.0.
//! - **`Rate`** — a per-base mutation rate. Per execution, the
//!   pass samples `count ~ Poisson(rate * pool_len)`. This is the
//!   biologically natural way to specify SHM intensity: a rate of
//!   `0.03` reads as "3 % of bases mutated, on average," which
//!   matches how immunologists report SHM in the literature.
//!
//! The rate-mode samples per-record against the *current pool
//! length*, not against a refdata-mean length, so trimmed records,
//! VJ vs VDJ chains, and per-record stochastic length differences
//! all see a rate that's a true fraction of their own length. This
//! is intentional — the architecture review explicitly rejected
//! compile-time conversion via mean refdata length on the grounds
//! that it "silently breaks for non-IGH and trimmed records."

use crate::dist::Distribution;
use crate::rng::Rng;

pub enum MutationCountSource {
    /// Empirical / explicit count distribution. Sampled once per
    /// pass execution; the result is the literal mutation count.
    Distribution(Box<dyn Distribution<Output = i64>>),
    /// Per-base mutation rate (e.g. `0.03` for 3 % SHM). Per pass
    /// execution, the count is drawn from `Poisson(rate * pool_len)`
    /// against the current pool length.
    Rate(f64),
}

impl MutationCountSource {
    /// Sample the mutation count for one execution of the pass.
    ///
    /// `pool_len` is consulted only in the `Rate` variant. Returns
    /// an `i64` to match the existing `Distribution<Output = i64>`
    /// shape — the caller does the same negative / overflow checks
    /// as before.
    pub fn sample(&self, rng: &mut Rng, pool_len: u32) -> i64 {
        match self {
            Self::Distribution(d) => d.sample(rng),
            Self::Rate(rate) => {
                let lambda = (*rate) * (pool_len as f64);
                poisson_sample(rng, lambda) as i64
            }
        }
    }

    /// True iff this source samples against the pool length. Used
    /// by the trace / describe layers to label the count event.
    #[allow(dead_code)]
    pub fn is_rate(&self) -> bool {
        matches!(self, Self::Rate(_))
    }
}

/// Sample a `Poisson(lambda)`-distributed value using Knuth's
/// algorithm. Numerically stable for the regime expected by this
/// engine (`lambda` between 0 and ~50, since SHM rates are <10%
/// and sequence lengths are <500bp). For `lambda <= 0`, returns 0.
///
/// Reference: Knuth, *The Art of Computer Programming*, vol 2,
/// 3.4.1. Mean / variance are both `lambda`. Runtime is `O(lambda)`
/// expected, so the algorithm is unsuitable for large `lambda` —
/// but for the genuinely-immunological regime where rate × len is
/// at most a few dozen, it's the right pick.
fn poisson_sample(rng: &mut Rng, lambda: f64) -> u32 {
    if !lambda.is_finite() || lambda <= 0.0 {
        return 0;
    }
    let exp_neg_lambda = (-lambda).exp();
    let mut k: u32 = 0;
    let mut p: f64 = 1.0;
    loop {
        k += 1;
        p *= rng.next_f64();
        if p <= exp_neg_lambda {
            return k - 1;
        }
        // Safety net for pathological `lambda` (large): cap at a
        // very generous bound so we never spin forever. The real
        // engine guard is the upstream `count_raw > u32::MAX` check
        // in the pass; here we just stop the loop.
        if k > 1_000_000 {
            return k;
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::dist::EmpiricalLengthDist;

    fn fixed(n: i64) -> MutationCountSource {
        MutationCountSource::Distribution(Box::new(EmpiricalLengthDist::from_pairs(vec![
            (n, 1.0),
        ])))
    }

    #[test]
    fn distribution_variant_ignores_pool_len() {
        let mut rng = Rng::new(42);
        let src = fixed(7);
        assert_eq!(src.sample(&mut rng, 1), 7);
        assert_eq!(src.sample(&mut rng, 1_000_000), 7);
    }

    #[test]
    fn rate_variant_returns_zero_on_zero_length() {
        let mut rng = Rng::new(42);
        let src = MutationCountSource::Rate(0.05);
        assert_eq!(src.sample(&mut rng, 0), 0);
    }

    #[test]
    fn rate_variant_returns_zero_for_zero_rate() {
        let mut rng = Rng::new(42);
        let src = MutationCountSource::Rate(0.0);
        assert_eq!(src.sample(&mut rng, 300), 0);
    }

    #[test]
    fn rate_variant_mean_matches_lambda_in_expectation() {
        // Run many samples and check the empirical mean is close
        // to lambda. lambda = 0.03 * 300 = 9 → expect 1000-sample
        // mean within ~10 % of 9.
        let mut rng = Rng::new(42);
        let src = MutationCountSource::Rate(0.03);
        let pool_len: u32 = 300;
        let n = 1000usize;
        let total: i64 = (0..n).map(|_| src.sample(&mut rng, pool_len)).sum();
        let mean = total as f64 / n as f64;
        let expected = 0.03 * (pool_len as f64);
        let tolerance = expected * 0.15; // 15 % envelope on the Monte Carlo estimate
        assert!(
            (mean - expected).abs() < tolerance,
            "empirical mean {mean} differs from expected {expected} by more than {tolerance}"
        );
    }

    #[test]
    fn rate_variant_is_deterministic_under_same_seed() {
        let src = MutationCountSource::Rate(0.05);
        let mut a = Rng::new(123);
        let mut b = Rng::new(123);
        for _ in 0..50 {
            assert_eq!(src.sample(&mut a, 250), src.sample(&mut b, 250));
        }
    }

    #[test]
    fn is_rate_discriminates_variants() {
        assert!(MutationCountSource::Rate(0.05).is_rate());
        assert!(!fixed(8).is_rate());
    }
}
