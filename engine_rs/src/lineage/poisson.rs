//! Self-contained Poisson sampler over the engine `Rng` (Knuth's algorithm).

use crate::rng::Rng;

/// Draw a Poisson(`lambda`) variate using Knuth's multiplicative method.
///
/// Deterministic for a given `Rng` state. Non-finite (NaN / ±inf) or
/// non-positive `lambda` always returns 0 — the `is_finite` guard matches the
/// engine's other Poisson sampler (`passes::count_source`) so an infinite
/// lambda cannot drive the loop to a pathological underflow-termination.
/// Suitable for the small lambdas used in clonal branching (offspring ~1.5,
/// mutations ~<1). Not optimized for very large lambda.
pub fn poisson_sample(rng: &mut Rng, lambda: f64) -> u32 {
    if !(lambda.is_finite() && lambda > 0.0) {
        return 0;
    }
    let l = (-lambda).exp();
    let mut k: u32 = 0;
    let mut p: f64 = 1.0;
    loop {
        k += 1;
        p *= rng.next_f64();
        if p <= l {
            return k - 1;
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::rng::Rng;

    #[test]
    fn poisson_zero_lambda_is_always_zero() {
        let mut rng = Rng::new(7);
        for _ in 0..100 {
            assert_eq!(poisson_sample(&mut rng, 0.0), 0);
        }
    }

    #[test]
    fn poisson_is_deterministic_for_a_seed() {
        let mut a = Rng::new(123);
        let mut b = Rng::new(123);
        let xs: Vec<u32> = (0..50).map(|_| poisson_sample(&mut a, 1.5)).collect();
        let ys: Vec<u32> = (0..50).map(|_| poisson_sample(&mut b, 1.5)).collect();
        assert_eq!(xs, ys);
    }

    #[test]
    fn poisson_nonfinite_lambda_is_zero() {
        let mut rng = Rng::new(3);
        assert_eq!(poisson_sample(&mut rng, f64::NAN), 0);
        assert_eq!(poisson_sample(&mut rng, f64::INFINITY), 0);
        assert_eq!(poisson_sample(&mut rng, -1.0), 0);
    }

    #[test]
    fn poisson_mean_is_near_lambda() {
        let mut rng = Rng::new(99);
        let n = 20_000u32;
        let lambda = 1.5;
        let total: u64 = (0..n).map(|_| poisson_sample(&mut rng, lambda) as u64).sum();
        let mean = total as f64 / n as f64;
        assert!((mean - lambda).abs() < 0.1, "mean {mean} not near {lambda}");
    }
}
