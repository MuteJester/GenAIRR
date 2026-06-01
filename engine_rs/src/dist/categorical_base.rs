//! `CategoricalBase` — weighted categorical over the four canonical DNA bases.
//!
//! Sibling of [`super::UniformBase`] but with per-base weights. The
//! audit's first-base empirical NP model (the v1 cartridge-owned NP
//! base distribution surface) consumes this type — every NP position
//! samples independently from the same weighted categorical.
//!
//! True Markov / previous-base-conditional sampling cannot be
//! expressed through the `Distribution<Output = u8>` trait because
//! `fn sample(&self, rng) -> u8` takes no context. That requires a
//! pass-level architectural change and is deferred to a follow-up
//! slice — see
//! [`docs/junction_n_addition_audit.md`](../../../../docs/junction_n_addition_audit.md)
//! §Markov-deferred.

use super::Distribution;
use crate::rng::Rng;

/// Canonical DNA bases the categorical supports, in canonical order.
/// The order here is the audit's pinned convention: A → C → G → T.
const CANONICAL_BASES: [u8; 4] = [b'A', b'C', b'G', b'T'];

/// Weighted categorical over the four canonical DNA bases
/// `{b'A', b'C', b'G', b'T'}`. Weights need not be normalised — the
/// distribution renormalises at sample time.
///
/// Construction rejects malformed input at the bridge layer (the
/// Python `NpBaseModelSpec` is the authoritative validator; this
/// `from_weights` re-checks defensively for custom Rust callers).
#[derive(Clone, Debug)]
pub struct CategoricalBase {
    /// Per-base weights in canonical A/C/G/T order. Each entry is
    /// non-negative finite; at least one entry is strictly positive.
    weights: [f64; 4],
    /// Cumulative weight after each base in canonical order
    /// (`cumulative[3]` is the total weight). Pre-computed for the
    /// inverse-CDF sampling path.
    cumulative: [f64; 4],
}

impl CategoricalBase {
    /// Construct from per-base weights in canonical order.
    ///
    /// Panics if any weight is NaN, negative, or infinite, or if all
    /// four weights are zero — these are construction-time invariants
    /// the Python DSL enforces before reaching the bridge; the panic
    /// is the defensive last line.
    pub fn from_weights(weights: [f64; 4]) -> Self {
        for (i, w) in weights.iter().enumerate() {
            assert!(
                w.is_finite() && *w >= 0.0,
                "CategoricalBase weight {i} must be finite and non-negative, got {w}",
            );
        }
        let total: f64 = weights.iter().sum();
        assert!(
            total > 0.0,
            "CategoricalBase must have at least one positive weight; total was {total}",
        );
        let mut cumulative = [0.0; 4];
        let mut running = 0.0;
        for i in 0..4 {
            running += weights[i];
            cumulative[i] = running;
        }
        Self { weights, cumulative }
    }

    /// Construct from a `(base, weight)` pair list. Bases not in the
    /// canonical A/C/G/T set are an error.
    pub fn from_pairs(pairs: Vec<(u8, f64)>) -> Self {
        let mut weights = [0.0f64; 4];
        for (base, weight) in pairs {
            let idx = match base {
                b'A' => 0,
                b'C' => 1,
                b'G' => 2,
                b'T' => 3,
                _ => panic!(
                    "CategoricalBase: unsupported base byte {base} (0x{base:02X}); \
                     only canonical A/C/G/T are allowed",
                ),
            };
            // Duplicate bases accumulate (matches the audit's
            // empirical-distribution discipline).
            weights[idx] += weight;
        }
        Self::from_weights(weights)
    }
}

impl Distribution for CategoricalBase {
    type Output = u8;

    fn sample(&self, rng: &mut Rng) -> u8 {
        // Inverse-CDF over the four canonical bases. Same RNG
        // consumption shape as the existing `EmpiricalLengthDist`
        // (`empirical.rs::sample`) so the categorical drops into
        // the same reproducibility discipline.
        let total = self.cumulative[3];
        let r = rng.next_f64() * total;
        // partition_point semantics on a 4-entry array — find the
        // first index where cumulative > r. Linear scan is faster
        // than a binary search at this scale.
        for (i, &cum) in self.cumulative.iter().enumerate() {
            if r < cum {
                return CANONICAL_BASES[i];
            }
        }
        // Float-rounding defensive: clamp to the last base.
        CANONICAL_BASES[3]
    }

    fn support(&self) -> Option<Vec<(u8, f64)>> {
        Some(
            CANONICAL_BASES
                .iter()
                .zip(self.weights.iter())
                .map(|(&b, &w)| (b, w))
                .collect(),
        )
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn uniform_weights_recover_uniform_base_support() {
        // Uniform weights produce the same support shape as
        // UniformBase — a future slice could collapse the two if
        // desired, but they stay separate today.
        let cat = CategoricalBase::from_weights([1.0, 1.0, 1.0, 1.0]);
        let support = cat.support().unwrap();
        assert_eq!(
            support,
            vec![(b'A', 1.0), (b'C', 1.0), (b'G', 1.0), (b'T', 1.0)]
        );
    }

    #[test]
    fn from_pairs_accumulates_duplicates() {
        let cat = CategoricalBase::from_pairs(vec![(b'A', 0.5), (b'A', 0.5), (b'G', 1.0)]);
        let support = cat.support().unwrap();
        assert_eq!(support[0], (b'A', 1.0));
        assert_eq!(support[1], (b'C', 0.0));
        assert_eq!(support[2], (b'G', 1.0));
        assert_eq!(support[3], (b'T', 0.0));
    }

    #[test]
    fn biased_weights_skew_sampling() {
        // With weight 1000 on 'A' and 1 elsewhere, ~997/1000 draws
        // should be 'A'.
        let cat = CategoricalBase::from_weights([1000.0, 1.0, 1.0, 1.0]);
        let mut rng = Rng::new(4242);
        let mut a_count = 0;
        for _ in 0..1000 {
            if cat.sample(&mut rng) == b'A' {
                a_count += 1;
            }
        }
        // Strict-but-wide bound: 990-1000 range gives the test air
        // for RNG variance without losing the bias signal.
        assert!(
            a_count >= 990,
            "biased categorical produced only {a_count}/1000 As — sampler bias broken",
        );
    }

    #[test]
    #[should_panic(expected = "must have at least one positive weight")]
    fn all_zero_weights_panic() {
        let _ = CategoricalBase::from_weights([0.0; 4]);
    }

    #[test]
    #[should_panic(expected = "unsupported base byte")]
    fn non_canonical_base_panics_in_from_pairs() {
        let _ = CategoricalBase::from_pairs(vec![(b'N', 1.0)]);
    }
}
