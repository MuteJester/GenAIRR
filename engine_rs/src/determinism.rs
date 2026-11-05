//! Which sampling model drew a dataset.
//!
//! The engine's reproducibility promise is narrow and worth stating
//! exactly: **one seed and one determinism version give one answer**,
//! on every platform the library builds on. The seed is the caller's
//! half of that; this constant is the engine's.
//!
//! Two runs at the same version agree byte for byte. Two runs at
//! different versions make no promise at all -- not "mostly the same",
//! not "statistically equivalent". A changed draw order moves every
//! record in the batch, so there is no useful middle claim to make.
//!
//! ## When to increment
//!
//! Whenever a released change alters what a given seed produces: a new
//! or removed draw, a different order of draws, a different mapping
//! from words to values, a different substream key. Not for a change
//! that only adds fields, renames a column or speeds something up with
//! the same draws.
//!
//! The number covers the library's sampling and not only this crate's.
//! Several runs split into sub-simulations in the Python layer -- a
//! clonal run into clones, retries, reads and nodes, a cohort into
//! subjects, a genotype draw into viability attempts -- and that layer
//! derives a seed for each (`GenAIRR._seed_schedule`). How a seed
//! becomes a sub-simulation's seed decides what a seed produces just as
//! directly as a pass does, so a change to it is a change to this
//! number. The clonal schedule carries a name of its own because a
//! family-replay bundle has to replay under the scheme that drew it,
//! but the name is not a second version: a reader holding a dataset
//! reads this one.
//!
//! Incrementing is not a failure to be avoided. Pinning a number while
//! the draws move underneath it is strictly worse than bumping it,
//! because it is the number that makes the difference discoverable.
//!
//! ## Why the numbering starts at 2
//!
//! Version 1 is the unnumbered model that predates the stamp:
//! everything GenAIRR released up to and including 2.3.1, whose
//! datasets carry no such column at all. Giving the first *stamped*
//! model the number 1 would have made 1 mean both "the first versioned
//! model" and "we do not know", which is the one thing a provenance
//! field must never do.
//!
//! | version | model |
//! |---------|-------|
//! | (absent) / 1 | one random stream shared by the whole plan; a substitution could redraw the base already present |
//! | 2 | one random substream per pass, and one per clonal-family role (growth, sampling, target) instead of a seed salted by XOR; a clone, a survival retry and a read take keyed substreams of the run seed instead of offsets added to it; a cohort subject and a genotype viability attempt take keyed substreams instead of words from a `random.Random` seeded with the run seed; a substitution draws a base the site does not hold |
//! | 3 | clonal lineage population regulation is Beverton-Holt, `lambda/(1 + (lambda-1)N/n_max)`, in place of the discrete logistic map `lambda(1 - N/n_max)`; `n_max` is the population a lineage settles at rather than the pivot of a law settling at `n_max(1 - 1/lambda)`, the law is monotone so it no longer oscillates above `lambda = 3`, and a generation may overshoot `n_max` instead of being truncated at it -- the truncation landed the population on `n_max` exactly, where the old law returned a rate of exactly zero, so reaching capacity extinguished the lineage |
//! | 4 | affinity selection scores the protein the cell expresses. The cell's amino-acid sequence is the continuous translation of its coding span (`productivity::coding_span`, resolved from the V allele's native reference phase) instead of the per-region codon rail, which dropped every codon straddling a region boundary -- all of them junction codons, so all of them CDR3; a stop or an unresolved residue is charged the maximal substitution cost instead of cancelling to zero against every target |

/// The engine's sampling-model version. Stamped on every exported
/// record and readable as `GenAIRR.DETERMINISM_VERSION`.
pub const DETERMINISM_VERSION: u32 = 4;
