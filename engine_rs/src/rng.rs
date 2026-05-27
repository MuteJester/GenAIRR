//! Deterministic pseudo-random number generator.
//!
//! ## Why we don't use `rand`
//!
//! For reproducibility we need a PRNG whose output is locked down
//! to *our* implementation, not whatever happens to be the current
//! default in the `rand` crate. We also want zero external
//! dependencies in the kernel until something genuinely needs one.
//!
//! ## SplitMix64
//!
//! SplitMix64 is the canonical "small splittable PRNG" — used as
//! the seeding stream of every Java 8+ `SplittableRandom` and the
//! initialization step of `xoshiro` family generators. It is:
//!
//! - **Deterministic given a seed** — same seed always produces the
//!   same stream.
//! - **Cheap** — one multiply + a few shifts + xors per word.
//! - **Acceptable quality** — passes BigCrush at 64 bits. Not
//!   cryptographic; we don't need cryptographic guarantees here.
//!
//! Reference: Vigna, "Further scramblings of Marsaglia's xorshift
//! generators" (2014); Steele et al. "Fast Splittable Pseudorandom
//! Number Generators" (2014, OOPSLA).
//!
//! ## Future
//!
//! When we need parallel batches across rayon worker threads we'll
//! add `Rng::split() -> Rng` for sub-stream generation. Not needed
//! yet.

/// Deterministic PRNG for the engine.
///
/// State is a single u64. Each `next_*` call advances the stream by
/// one word. Reseeding via `Rng::new(seed)` is idempotent — the same
/// seed always produces the same stream, byte for byte.
#[derive(Clone, Debug)]
pub struct Rng {
    state: u64,
    words_consumed: u64,
}

impl Rng {
    /// Construct an `Rng` from a 64-bit seed. Seed `0` is replaced
    /// with a fixed non-zero constant so the stream cannot stall.
    pub fn new(seed: u64) -> Self {
        Self {
            state: if seed == 0 {
                0x9e37_79b9_7f4a_7c15
            } else {
                seed
            },
            words_consumed: 0,
        }
    }

    /// Number of 64-bit PRNG words consumed from this stream.
    ///
    /// This is an observability hook for deterministic tests and
    /// replay diagnostics. It does not affect the generated stream.
    pub fn words_consumed(&self) -> u64 {
        self.words_consumed
    }

    /// Advance the stream and return the next 64-bit word.
    /// SplitMix64 finalizer — see Vigna 2014.
    pub fn next_u64(&mut self) -> u64 {
        self.words_consumed = self.words_consumed.saturating_add(1);
        self.state = self.state.wrapping_add(0x9e37_79b9_7f4a_7c15);
        let mut z = self.state;
        z = (z ^ (z >> 30)).wrapping_mul(0xbf58_476d_1ce4_e5b9);
        z = (z ^ (z >> 27)).wrapping_mul(0x94d0_49bb_1331_11eb);
        z ^ (z >> 31)
    }

    /// Next 32-bit word. Takes the high 32 bits of `next_u64`,
    /// which is the recommended way to extract a 32-bit value
    /// from SplitMix64 (the high bits have better statistical
    /// properties than the low bits).
    pub fn next_u32(&mut self) -> u32 {
        (self.next_u64() >> 32) as u32
    }

    /// Sample a uniformly-distributed integer in `[0, max)`.
    /// Returns 0 when `max` is 0. Uses Lemire's nearly-divisionless
    /// rejection method for unbiased range mapping.
    ///
    /// Reference: Lemire, "Fast Random Integer Generation in an
    /// Interval" (2019, ACM TOMS).
    pub fn range_u32(&mut self, max: u32) -> u32 {
        if max == 0 {
            return 0;
        }
        // Fast path covers the vast majority of draws; the rejection
        // loop only fires on a tiny biased slice near u32::MAX.
        let max = max as u64;
        let mut x = self.next_u32() as u64;
        let mut m = x.wrapping_mul(max);
        let mut l = m as u32;
        if (l as u64) < max {
            let t = (u32::MAX as u64).wrapping_sub(max - 1) % max;
            while (l as u64) < t {
                x = self.next_u32() as u64;
                m = x.wrapping_mul(max);
                l = m as u32;
            }
        }
        (m >> 32) as u32
    }

    /// Sample a uniformly-distributed `f64` in `[0.0, 1.0)`.
    /// Constructs the float from the high 53 bits of `next_u64`,
    /// which is the standard portable approach for getting a
    /// double-precision uniform.
    pub fn next_f64(&mut self) -> f64 {
        // 53 high bits, then divide by 2^53. Equivalent to
        // (n as f64) * (1.0 / (1u64 << 53) as f64) without
        // intermediate float widening surprises.
        let bits = self.next_u64() >> 11;
        (bits as f64) * (1.0 / ((1u64 << 53) as f64))
    }
}

// ──────────────────────────────────────────────────────────────────
// Tests
// ──────────────────────────────────────────────────────────────────

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn same_seed_produces_same_stream() {
        let mut a = Rng::new(42);
        let mut b = Rng::new(42);
        for _ in 0..1000 {
            assert_eq!(a.next_u64(), b.next_u64());
        }
    }

    #[test]
    fn different_seeds_produce_different_streams() {
        let mut a = Rng::new(42);
        let mut b = Rng::new(43);
        // The first word should already differ for any reasonable PRNG.
        assert_ne!(a.next_u64(), b.next_u64());
    }

    #[test]
    fn zero_seed_is_replaced_with_nonzero_constant() {
        // A zero seed used to stall xorshift; SplitMix64 doesn't
        // suffer from that, but we still normalize to a non-zero
        // constant so the stream is well-defined for `seed=0`.
        let mut a = Rng::new(0);
        let mut b = Rng::new(0);
        // Both produce the same stream (same seed → same stream),
        // and that stream is non-trivial (first word is not zero).
        let first = a.next_u64();
        assert_ne!(first, 0);
        assert_eq!(first, b.next_u64());
    }

    #[test]
    fn next_u32_consumes_one_u64_per_call() {
        // Round-trip: pulling a u32 should advance the state. Two
        // consecutive u32 draws give different values for any seed.
        let mut r = Rng::new(7);
        assert_eq!(r.words_consumed(), 0);
        let a = r.next_u32();
        assert_eq!(r.words_consumed(), 1);
        let b = r.next_u32();
        assert_eq!(r.words_consumed(), 2);
        assert_ne!(a, b);
    }

    #[test]
    fn range_u32_zero_max_returns_zero() {
        let mut r = Rng::new(1);
        assert_eq!(r.range_u32(0), 0);
        // Calling range_u32(0) does not advance state or panic.
        assert_eq!(r.range_u32(0), 0);
        assert_eq!(r.words_consumed(), 0);
    }

    #[test]
    fn rng_tracks_consumed_words_without_changing_stream() {
        let mut counted = Rng::new(123);
        let mut plain = Rng::new(123);

        assert_eq!(counted.words_consumed(), 0);
        assert_eq!(counted.next_u64(), plain.next_u64());
        assert_eq!(counted.words_consumed(), 1);
        assert_eq!(counted.next_f64(), plain.next_f64());
        assert_eq!(counted.words_consumed(), 2);
        assert_eq!(counted.range_u32(16), plain.range_u32(16));
        assert_eq!(counted.words_consumed(), 3);
    }

    #[test]
    fn range_u32_stays_in_bounds() {
        let mut r = Rng::new(1234);
        for _ in 0..10_000 {
            let x = r.range_u32(7);
            assert!(x < 7, "range_u32(7) returned out-of-bounds {}", x);
        }
        for _ in 0..10_000 {
            let x = r.range_u32(1);
            assert_eq!(x, 0);
        }
    }

    #[test]
    fn range_u32_covers_full_range() {
        // Sample a small range many times and confirm every value
        // appears at least once. Catches degenerate biases.
        let mut r = Rng::new(99);
        let mut seen = [false; 5];
        for _ in 0..10_000 {
            let x = r.range_u32(5) as usize;
            seen[x] = true;
        }
        for (i, was_seen) in seen.iter().enumerate() {
            assert!(was_seen, "range_u32(5) never produced {}", i);
        }
    }

    #[test]
    fn next_f64_in_unit_interval() {
        let mut r = Rng::new(0xc0ff_ee);
        for _ in 0..10_000 {
            let x = r.next_f64();
            assert!((0.0..1.0).contains(&x), "next_f64 produced {}", x);
        }
    }

    #[test]
    fn next_f64_distribution_is_roughly_uniform() {
        // Bin into 10 bins of width 0.1, draw 10,000 samples, each
        // bin should hit roughly 1,000. Allow generous slack so this
        // test is rock-solid against statistical jitter.
        let mut r = Rng::new(0xfeed_face);
        let mut bins = [0u32; 10];
        let n = 10_000u32;
        for _ in 0..n {
            let x = r.next_f64();
            let bin = (x * 10.0).floor() as usize;
            assert!(bin < 10);
            bins[bin] += 1;
        }
        // Each bin should be in [700, 1300] — Chernoff bound makes
        // this comfortable for a uniform PRNG.
        for (i, &c) in bins.iter().enumerate() {
            assert!(
                (700..=1300).contains(&c),
                "bin {} count {} outside [700, 1300]",
                i,
                c
            );
        }
    }

    #[test]
    fn rng_clone_is_independent() {
        let mut a = Rng::new(123);
        let _ = a.next_u64();
        let mut b = a.clone();
        // Clone captures current state, then evolves independently.
        let av = a.next_u64();
        let bv = b.next_u64();
        assert_eq!(av, bv); // both start from the same state
        let av2 = a.next_u64();
        let bv2 = b.next_u64();
        assert_eq!(av2, bv2); // continue lock-step
    }
}
