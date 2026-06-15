//! Clone-size distributions for repertoire composition. T-cell (and the
//! singleton tail of B-cell) repertoires have heavy-tailed clone-size
//! distributions; this module samples integer clone sizes >= 1.

use crate::rng::Rng;

/// A heavy-tailed clone-size distribution. Sizes are integers in `[1, x_max]`.
#[derive(Clone, Debug)]
pub enum CloneSizeDist {
    /// Truncated discrete power law P(size=k) ∝ k^(-exponent), k in [1, x_max].
    /// `exponent` ~2–3 is typical for TCR repertoires.
    PowerLaw { exponent: f64, x_max: u32 },
    /// Log-normal sizes: round(exp(mu + sigma*Z)), clamped to [1, x_max].
    LogNormal { mu: f64, sigma: f64, x_max: u32 },
}

impl Default for CloneSizeDist {
    fn default() -> Self {
        CloneSizeDist::PowerLaw {
            exponent: 2.0,
            x_max: 100_000,
        }
    }
}

/// Standard-normal variate via Box–Muller using two uniforms from `rng`.
fn next_standard_normal(rng: &mut Rng) -> f64 {
    let mut u1 = rng.next_f64();
    if u1 < 1e-300 {
        u1 = 1e-300;
    }
    let u2 = rng.next_f64();
    (-2.0 * u1.ln()).sqrt() * (std::f64::consts::TAU * u2).cos()
}

/// Sample one clone size (>= 1) from `dist`, deterministic for the `rng` state.
pub fn sample_clone_size(rng: &mut Rng, dist: &CloneSizeDist) -> u32 {
    match *dist {
        CloneSizeDist::PowerLaw { exponent, x_max } => {
            let x_max = x_max.max(1) as f64;
            let u = rng.next_f64();
            let x = if (exponent - 1.0).abs() < 1e-9 {
                x_max.powf(u)
            } else {
                let a = 1.0 - exponent;
                (1.0 + u * (x_max.powf(a) - 1.0)).powf(1.0 / a)
            };
            (x.round() as i64).clamp(1, x_max as i64) as u32
        }
        CloneSizeDist::LogNormal { mu, sigma, x_max } => {
            let z = next_standard_normal(rng);
            let x = (mu + sigma * z).exp();
            (x.round() as i64).clamp(1, x_max.max(1) as i64) as u32
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::rng::Rng;

    #[test]
    fn sizes_are_at_least_one_and_within_max() {
        let mut rng = Rng::new(1);
        let d = CloneSizeDist::PowerLaw { exponent: 2.0, x_max: 1000 };
        for _ in 0..1000 {
            let s = sample_clone_size(&mut rng, &d);
            assert!(s >= 1 && s <= 1000, "size {s} out of [1,1000]");
        }
    }

    #[test]
    fn power_law_is_deterministic() {
        let d = CloneSizeDist::PowerLaw { exponent: 2.0, x_max: 1000 };
        let mut a = Rng::new(42);
        let mut b = Rng::new(42);
        let xs: Vec<u32> = (0..100).map(|_| sample_clone_size(&mut a, &d)).collect();
        let ys: Vec<u32> = (0..100).map(|_| sample_clone_size(&mut b, &d)).collect();
        assert_eq!(xs, ys);
    }

    #[test]
    fn power_law_is_heavy_tailed_mostly_small_some_large() {
        let mut rng = Rng::new(7);
        let d = CloneSizeDist::PowerLaw { exponent: 2.0, x_max: 100_000 };
        let n = 20_000;
        let sizes: Vec<u32> = (0..n).map(|_| sample_clone_size(&mut rng, &d)).collect();
        let ones = sizes.iter().filter(|&&s| s == 1).count();
        let big = sizes.iter().filter(|&&s| s >= 100).count();
        assert!(ones > n / 3, "expected many singletons, got {ones}/{n}");
        assert!(big > 0, "expected some large clones in the tail");
        let max = *sizes.iter().max().unwrap();
        assert!(max > 50, "tail did not reach large sizes (max {max})");
    }

    #[test]
    fn lognormal_is_deterministic_and_in_range() {
        let d = CloneSizeDist::LogNormal { mu: 1.0, sigma: 1.0, x_max: 10_000 };
        let mut a = Rng::new(9);
        let mut b = Rng::new(9);
        let xs: Vec<u32> = (0..100).map(|_| sample_clone_size(&mut a, &d)).collect();
        let ys: Vec<u32> = (0..100).map(|_| sample_clone_size(&mut b, &d)).collect();
        assert_eq!(xs, ys);
        for &s in &xs {
            assert!(s >= 1 && s <= 10_000);
        }
    }
}
