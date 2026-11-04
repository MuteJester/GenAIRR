//! The named streams a clonal family draws from.
//!
//! A family's run makes three kinds of choice, and each wants its
//! draws unmoved by the others: **growth** (which cells divide, how
//! many times, and each child's seed), **sampling** (which of the
//! surviving cells are observed) and **target** (the affinity target
//! a maturing family is selected toward). They are roles in one run
//! seeded once, exactly the situation [`crate::rng::substream`]
//! exists for, so each takes a named substream of the family's seed.
//!
//! ## What this replaced, and why it was wrong
//!
//! These three used to be spelled by salting the seed with XOR:
//! `Rng::new(params.seed ^ SAMPLE_SEED_SALT)` and
//! `Rng::new(params.seed ^ 0x7461_7267_6574_0001)`. That does not
//! separate anything. XOR is invertible, so a salted stream *is* the
//! plain stream of another seed, and for any two salts there is an
//! offset `d = salt_a ^ salt_b` such that the first role's stream in
//! a run seeded `s` is bit-identical to the second role's stream in
//! the run seeded `s ^ d`. For the two salts this module replaced,
//! `d = 0x2720_3f37_2931_0000`: a family's sampling draws were
//! another family's target draws, word for word. The test below
//! holds the measurement, and holds it against the old scheme too,
//! so the defect cannot come back unnoticed.
//!
//! Folding a label in through [`crate::rng::fold_in`] has no such
//! offset, because the fold is not invertible and mixes both words.

use crate::rng::substream::{StreamRoot, SubstreamLabel};
use crate::rng::Rng;

/// Which cells divide, how many times, and each child's seed.
pub const GROWTH: &str = "lineage.growth";

/// Which of the surviving cells are observed.
pub const SAMPLING: &str = "lineage.sampling";

/// The affinity target a maturing family is selected toward.
pub const TARGET: &str = "lineage.target";

/// The stream `role` draws from in a family seeded `seed`.
///
/// Call it at the point of use rather than holding the result: a
/// substream is a pure function of the seed and the role, so two
/// calls with the same arguments are the same stream, and nothing
/// about the order they are made in can change either.
pub fn stream(seed: u64, role: &str) -> Rng {
    StreamRoot::new(seed).substream(SubstreamLabel::named(role, 0))
}

#[cfg(test)]
mod tests {
    use super::{stream, GROWTH, SAMPLING, TARGET};
    use crate::rng::Rng;

    /// The salts this module replaced, kept only so the test below
    /// can show what they did.
    const OLD_SAMPLE_SALT: u64 = 0x5341_4D50_4C45_0001;
    const OLD_TARGET_SALT: u64 = 0x7461_7267_6574_0001;

    fn first_words(mut rng: Rng, n: usize) -> Vec<u64> {
        (0..n).map(|_| rng.next_u64()).collect()
    }

    #[test]
    fn a_role_is_the_same_stream_every_time() {
        for seed in [0u64, 1, 42, u64::MAX] {
            for role in [GROWTH, SAMPLING, TARGET] {
                assert_eq!(
                    first_words(stream(seed, role), 8),
                    first_words(stream(seed, role), 8),
                    "seed {seed} role {role}"
                );
            }
        }
    }

    #[test]
    fn the_three_roles_are_three_streams() {
        for seed in 0..64u64 {
            let g = first_words(stream(seed, GROWTH), 4);
            let s = first_words(stream(seed, SAMPLING), 4);
            let t = first_words(stream(seed, TARGET), 4);
            assert_ne!(g, s, "seed {seed}");
            assert_ne!(g, t, "seed {seed}");
            assert_ne!(s, t, "seed {seed}");
        }
    }

    /// The defect this module exists to remove, stated as the
    /// measurement that found it: under XOR salting there is one
    /// seed offset at which a family's sampling stream is another
    /// family's target stream, word for word.
    #[test]
    fn xor_salting_made_one_familys_sampling_another_familys_target() {
        let offset = OLD_SAMPLE_SALT ^ OLD_TARGET_SALT;
        assert_eq!(offset, 0x2720_3f37_2931_0000);

        let seed = 1234u64;
        assert_eq!(
            first_words(Rng::new(seed ^ OLD_SAMPLE_SALT), 8),
            first_words(Rng::new((seed ^ offset) ^ OLD_TARGET_SALT), 8),
            "the old scheme's two streams coincide at this offset -- if this \
             ever stops holding, the regression test below is testing nothing"
        );
    }

    #[test]
    fn no_offset_makes_one_familys_sampling_another_familys_target() {
        // The offset that broke the old scheme, and then a sweep: no
        // pair of seeds within it lines the two roles up at all.
        let offset = OLD_SAMPLE_SALT ^ OLD_TARGET_SALT;
        for seed in [0u64, 1, 42, 1234, u64::MAX] {
            assert_ne!(
                first_words(stream(seed, SAMPLING), 4),
                first_words(stream(seed ^ offset, TARGET), 4),
                "seed {seed} at the historical offset"
            );
        }

        let sampling: Vec<Vec<u64>> = (0..128u64)
            .map(|seed| first_words(stream(seed, SAMPLING), 2))
            .collect();
        let targets: std::collections::HashSet<Vec<u64>> = (0..128u64)
            .map(|seed| first_words(stream(seed, TARGET), 2))
            .collect();
        for (seed, words) in sampling.iter().enumerate() {
            assert!(
                !targets.contains(words),
                "seed {seed}'s sampling stream is some seed's target stream"
            );
        }
    }

    #[test]
    fn a_role_is_named_not_positioned() {
        // Roles are independent of each other's existence: naming
        // one cannot shift another, whatever order they are asked
        // for. This is the property the pass plan relies on too.
        let seed = 99u64;
        let expected = first_words(stream(seed, TARGET), 4);
        let _ = first_words(stream(seed, GROWTH), 100);
        let _ = first_words(stream(seed, SAMPLING), 100);
        assert_eq!(first_words(stream(seed, TARGET), 4), expected);
    }
}
