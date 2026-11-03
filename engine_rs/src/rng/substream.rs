//! Per-pass random substreams.
//!
//! One run draws from many streams, one per pass, each derived from
//! the run's seed and the pass's identity. Nothing is taken from a
//! parent stream, so a pass's draws depend on the pass and the seed
//! and on nothing else in the plan — which is what makes a seeded
//! ablation ("the same repertoire, with and without one pass") a
//! controlled comparison about the molecule rather than only about
//! the rearrangement.
//!
//! ## Why identity and not position
//!
//! A pass's stream is keyed by `(run seed, pass name, n-th pass of
//! that name)`, never by its index in the plan. Keying on the index
//! would reproduce the coupling this module exists to remove: insert
//! or delete one pass and every later index shifts, so every later
//! stream moves. Two `pcr_amplify` steps in one plan stay
//! distinguishable because the occurrence counter separates them,
//! and `mutate` draws the same words whether or not an `invert_d`
//! sits in front of it.
//!
//! ## Why a root that cannot draw
//!
//! [`StreamRoot`] has no `next_*` method. The run loop holds one and
//! can therefore only derive from it, which makes "no pass reaches a
//! shared stream" a property of the types rather than a convention
//! that holds until someone adds a draw. Derivation takes `&self`
//! and advances nothing, so the order in which substreams are
//! derived cannot matter either.

use super::{fold_in, Rng};

/// FNV-1a 64-bit offset basis and prime. FNV because the hash has to
/// be written down and never change: the pass names are short ASCII
/// labels, the only requirement is that distinct names give distinct
/// words, and a hash whose definition fits in four lines cannot drift
/// under us the way a library default can.
const FNV_OFFSET: u64 = 0xcbf2_9ce4_8422_2325;
const FNV_PRIME: u64 = 0x0000_0100_0000_01b3;

#[inline]
fn fnv1a64(bytes: &[u8]) -> u64 {
    let mut hash = FNV_OFFSET;
    for byte in bytes {
        hash ^= *byte as u64;
        hash = hash.wrapping_mul(FNV_PRIME);
    }
    hash
}

/// A pass's stream identity, as one word.
///
/// Derived from the pass's name and which occurrence of that name it
/// is, both of which are properties of the plan. Computed once when
/// the pass joins a [`crate::pass::Schedule`], not once per record.
#[derive(Copy, Clone, Debug, Eq, PartialEq, Hash)]
pub struct SubstreamLabel(u64);

impl SubstreamLabel {
    /// The label for the `occurrence`-th pass named `name`, counting
    /// from zero in the order the passes were added to the plan.
    pub fn for_pass(name: &str, occurrence: u32) -> Self {
        Self(fold_in(fnv1a64(name.as_bytes()), occurrence as u64))
    }
}

/// One run's seed, in the one form from which substreams are taken.
///
/// Deliberately not a generator: see the module docs.
#[derive(Copy, Clone, Debug)]
pub struct StreamRoot {
    seed: u64,
}

impl StreamRoot {
    /// The root for a run seeded with `seed`.
    pub fn new(seed: u64) -> Self {
        Self { seed }
    }

    /// The stream belonging to `label` under this run's seed.
    pub fn substream(&self, label: SubstreamLabel) -> Rng {
        Rng::new(self.substream_seed(label))
    }

    /// The seed of that stream. Exposed for oracles that step the
    /// generator by hand to check which word a choice lands on: they
    /// need somewhere to start, and starting from the run seed stopped
    /// being right when passes stopped sharing a stream.
    pub fn substream_seed(&self, label: SubstreamLabel) -> u64 {
        fold_in(self.seed, label.0)
    }
}

// ──────────────────────────────────────────────────────────────────
// Tests
// ──────────────────────────────────────────────────────────────────

#[cfg(test)]
mod tests {
    use super::*;

    /// Real pass names, so the properties below are asserted about the
    /// strings the engine actually hands this module — several of
    /// which share long prefixes and differ in one trailing byte.
    const PLAN: &[&str] = &[
        "sample_allele.V",
        "sample_allele.D",
        "sample_allele.J",
        "trim.V",
        "trim.D",
        "trim.J",
        "generate_np.1",
        "generate_np.2",
        "assemble_segment.V",
        "assemble_segment.D",
        "assemble_segment.J",
        "invert_d",
        "receptor_revision",
        "mutate.uniform",
        "mutate.s5f",
        "corrupt.pcr",
        "corrupt.quality",
        "corrupt.ns",
        "corrupt.indel",
        "corrupt.end_loss",
    ];

    fn first_word(root: StreamRoot, name: &str, occurrence: u32) -> u64 {
        root.substream(SubstreamLabel::for_pass(name, occurrence))
            .next_u64()
    }

    #[test]
    fn the_same_pass_under_the_same_seed_draws_the_same_words() {
        let root = StreamRoot::new(4242);
        let mut a = root.substream(SubstreamLabel::for_pass("mutate.uniform", 0));
        let mut b = root.substream(SubstreamLabel::for_pass("mutate.uniform", 0));
        for _ in 0..1000 {
            assert_eq!(a.next_u64(), b.next_u64());
        }
    }

    #[test]
    fn two_passes_of_the_same_name_draw_different_words() {
        // Two `pcr_amplify` steps in one plan must not apply the same
        // errors at the same sites, which is what sharing a stream
        // would do.
        let root = StreamRoot::new(4242);
        assert_ne!(
            first_word(root, "corrupt.pcr", 0),
            first_word(root, "corrupt.pcr", 1)
        );
    }

    #[test]
    fn every_pass_in_a_plan_draws_a_different_first_word() {
        // Including the names that differ only in their last byte:
        // `sample_allele.V` against `.D` against `.J`, and
        // `trim.*` likewise.
        for seed in [0, 1, 42, u64::MAX] {
            let root = StreamRoot::new(seed);
            let mut words: Vec<u64> = PLAN
                .iter()
                .map(|name| first_word(root, name, 0))
                .collect();
            let before = words.len();
            words.sort_unstable();
            words.dedup();
            assert_eq!(
                words.len(), before,
                "two passes of the plan share a stream under seed {seed}"
            );
        }
    }

    #[test]
    fn the_default_seed_separates_the_passes_like_any_other() {
        // `seed=0` is the default everywhere in the Python API and a
        // fixed point of the SplitMix64 finalizer, so it is the seed a
        // careless derivation degenerates on. Asserted as its own test
        // rather than folded into the loop above, because it is the
        // one most runs use.
        let root = StreamRoot::new(0);
        let words: Vec<u64> = PLAN.iter().map(|name| first_word(root, name, 0)).collect();
        assert!(words.iter().all(|word| *word != 0));
        let mut sorted = words.clone();
        sorted.sort_unstable();
        sorted.dedup();
        assert_eq!(sorted.len(), words.len());
    }

    #[test]
    fn one_pass_draws_different_words_in_consecutive_records() {
        // A batch seeds record `i` with `seed + i`, so a pass that
        // repeated itself across neighbouring seeds would repeat
        // itself across neighbouring records.
        let words: Vec<u64> = (0..4096)
            .map(|seed| first_word(StreamRoot::new(seed), "mutate.uniform", 0))
            .collect();
        let mut sorted = words.clone();
        sorted.sort_unstable();
        sorted.dedup();
        assert_eq!(sorted.len(), words.len());
    }

    #[test]
    fn drawing_from_one_substream_cannot_move_another() {
        // The property the whole phase is for, stated at the level it
        // is true at: exhausting one pass's stream leaves every other
        // pass's stream exactly where it was, and so does deriving
        // them in a different order.
        let root = StreamRoot::new(4242);
        let mut spent = root.substream(SubstreamLabel::for_pass("corrupt.pcr", 0));
        for _ in 0..10_000 {
            spent.next_u64();
        }
        let after = first_word(root, "mutate.uniform", 0);

        let fresh = StreamRoot::new(4242);
        let untouched = first_word(fresh, "mutate.uniform", 0);
        assert_eq!(
            after, untouched,
            "ten thousand draws in one pass's stream moved another pass's"
        );
    }

    #[test]
    fn a_root_is_not_advanced_by_deriving_from_it() {
        // `substream` takes `&self`; this is the behavioural statement
        // of that signature, so a later change to interior state shows
        // up here rather than in a regenerated golden file.
        let root = StreamRoot::new(7);
        let labels: Vec<SubstreamLabel> = PLAN
            .iter()
            .map(|name| SubstreamLabel::for_pass(name, 0))
            .collect();
        let forward: Vec<u64> = labels
            .iter()
            .map(|label| root.substream(*label).next_u64())
            .collect();
        let mut backward: Vec<u64> = labels
            .iter()
            .rev()
            .map(|label| root.substream(*label).next_u64())
            .collect();
        backward.reverse();
        assert_eq!(forward, backward);
    }
}
