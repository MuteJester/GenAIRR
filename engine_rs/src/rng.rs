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
//! ## Substreams
//!
//! One stream shared by a whole plan makes every pass's draws depend
//! on how many words the passes in front of it happened to consume.
//! [`fold_in`] is the primitive that ends that coupling: a substream
//! seed is *computed* from a run seed and an identifying word rather
//! than taken from the parent stream, so no amount of drawing in one
//! substream can move another, and no ordering of the derivations can
//! either. [`crate::rng::substream`] builds the typed layer on top of
//! it; this module stays ignorant of what the identifying word means.

pub mod substream;

/// Guard: every production generator is born in a place that says
/// why it is not a substream. Test-time only; nothing here ships.
#[cfg(test)]
mod discipline;

/// The SplitMix64 increment — 2^64 divided by the golden ratio. An
/// odd constant, so repeatedly adding it walks the whole 64-bit state
/// space before returning.
const GOLDEN_GAMMA: u64 = 0x9e37_79b9_7f4a_7c15;

/// The SplitMix64 finalizer — see Vigna 2014. A bijection on `u64`
/// that makes every output bit depend on every input bit, which is
/// what both the generator and [`fold_in`] need and the only reason
/// either is sound.
#[inline]
fn mix64(mut z: u64) -> u64 {
    z = (z ^ (z >> 30)).wrapping_mul(0xbf58_476d_1ce4_e5b9);
    z = (z ^ (z >> 27)).wrapping_mul(0x94d0_49bb_1331_11eb);
    z ^ (z >> 31)
}

/// Fold an identifying word into a seed, yielding a seed whose stream
/// is unrelated to the original's.
///
/// This is seed derivation, not splitting: it reads no state and
/// advances nothing, so `fold_in(seed, a)` is the same stream however
/// many other words were folded into the same seed first, and in
/// whatever order. That is the property per-pass substreams need —
/// sequential splitting would reproduce the coupling it removes,
/// moving it from the stream position into the split count.
///
/// Both arguments pass through the finalizer before they meet, so
/// neither a zero seed (the library default, and a fixed point of
/// [`mix64`]) nor a zero data word degenerates. Each step is a
/// bijection, so distinct inputs give distinct seeds — the collision
/// bound is exact, not probabilistic.
#[inline]
pub fn fold_in(seed: u64, data: u64) -> u64 {
    mix64(mix64(seed) ^ mix64(data.wrapping_add(GOLDEN_GAMMA)))
}

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
            state: if seed == 0 { GOLDEN_GAMMA } else { seed },
            words_consumed: 0,
        }
    }

    /// A generator for a path whose draws cannot matter.
    ///
    /// A few places need an `Rng` to satisfy a signature while
    /// drawing nothing that survives: a replay whose every choice
    /// comes from a recorded trace, where a stray draw is rejected by
    /// a consumption fence; a probe whose result is discarded; a
    /// proof that a sampling support is empty. Those sites used to
    /// write `Rng::new(1)` or `Rng::new(0)`, which reads like a seed
    /// somebody chose for a reason and leaves the next reader to work
    /// out that nothing depends on it.
    ///
    /// The state behind it is arbitrary and no caller may depend on
    /// it. A site that would behave differently under a different
    /// placeholder is drawing something that matters, and wants a
    /// substream.
    ///
    /// It takes no seed, so it cannot hand a pass a stream of its own
    /// by the back door -- which is why the guard in
    /// [`discipline`](crate::rng) watches seed-taking construction and
    /// not this.
    pub fn placeholder() -> Self {
        Self::new(1)
    }

    /// Number of 64-bit PRNG words consumed from this stream.
    ///
    /// This is an observability hook for deterministic tests and
    /// replay diagnostics. It does not affect the generated stream.
    pub fn words_consumed(&self) -> u64 {
        self.words_consumed
    }

    /// Advance the stream and return the next 64-bit word.
    pub fn next_u64(&mut self) -> u64 {
        crate::family_replay::primitive_random_word();
        self.words_consumed = self.words_consumed.saturating_add(1);
        self.state = self.state.wrapping_add(GOLDEN_GAMMA);
        mix64(self.state)
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
    fn a_placeholder_is_a_working_stream_with_no_hidden_state() {
        // The three sites that use it draw nothing that survives, but
        // a draw must not panic or stall: a fence rejecting a stray
        // draw has to see the draw happen first.
        let mut a = Rng::placeholder();
        let mut b = Rng::placeholder();
        let first: Vec<u64> = (0..4).map(|_| a.next_u64()).collect();
        let second: Vec<u64> = (0..4).map(|_| b.next_u64()).collect();
        assert_eq!(first, second, "two placeholders are the same stream");
        assert!(
            first.iter().any(|&w| w != 0),
            "a stalled stream would make a fence's word count meaningless"
        );
        assert_eq!(a.words_consumed(), 4);
    }

    #[test]
    fn a_placeholder_cannot_carry_a_runs_seed() {
        // Why `discipline` watches `Rng::new` and ignores this: a
        // placeholder takes no seed, so it cannot be used to give a
        // pass a stream of its own while dodging the guard.
        let from_placeholder = Rng::placeholder().next_u64();
        let distinct = (0..512u64)
            .map(|seed| Rng::new(seed).next_u64())
            .filter(|&w| w == from_placeholder)
            .count();
        assert!(
            distinct <= 1,
            "a placeholder's stream is one fixed stream, not a seed's"
        );
    }

    #[test]
    fn fold_in_is_a_pure_function_of_its_two_words() {
        // The whole reproducibility claim rests on this: a substream's
        // seed is computed, never accumulated, so it cannot depend on
        // how many other substreams were derived first.
        assert_eq!(fold_in(7, 11), fold_in(7, 11));
        assert_ne!(fold_in(7, 11), fold_in(11, 7));
    }

    #[test]
    fn folding_consecutive_data_words_gives_unrelated_streams() {
        // `data` is a label hash in production, but the labels of two
        // passes with the same name differ by 1, so neighbouring words
        // are the case that has to work.
        let first_words: Vec<u64> = (0..4096)
            .map(|data| Rng::new(fold_in(42, data)).next_u64())
            .collect();
        let mut sorted = first_words.clone();
        sorted.sort_unstable();
        sorted.dedup();
        assert_eq!(
            sorted.len(),
            first_words.len(),
            "two of 4096 consecutive data words collided into one stream"
        );
    }

    #[test]
    fn folding_into_a_zero_seed_still_depends_on_the_data() {
        // `seed=0` is the library's default seed, and the SplitMix64
        // finalizer maps 0 to 0 — so a derivation that mixed the seed
        // and then xored the data in would degenerate on exactly the
        // seed most runs use. Asserted on the real default rather than
        // on an arbitrary one.
        let first_words: Vec<u64> = (0..4096)
            .map(|data| Rng::new(fold_in(0, data)).next_u64())
            .collect();
        assert!(
            first_words.iter().all(|word| *word != 0),
            "a derived stream starts at zero"
        );
        let mut sorted = first_words.clone();
        sorted.sort_unstable();
        sorted.dedup();
        assert_eq!(sorted.len(), first_words.len(), "seed 0 collapsed two labels");
    }

    #[test]
    fn folding_the_same_data_into_consecutive_seeds_gives_unrelated_streams() {
        // Batch runs seed record i with `seed + i`, so consecutive
        // seeds holding one pass's label is the other case that has to
        // work: record 7's `mutate` must not echo record 8's.
        let first_words: Vec<u64> = (0..4096)
            .map(|seed| Rng::new(fold_in(seed, 0xdead_beef)).next_u64())
            .collect();
        let mut sorted = first_words.clone();
        sorted.sort_unstable();
        sorted.dedup();
        assert_eq!(
            sorted.len(),
            first_words.len(),
            "two of 4096 consecutive seeds collided for one label"
        );
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
