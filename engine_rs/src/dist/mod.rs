//! Distributions — typed sources of random samples.
//!
//! ## Design (D4)
//!
//! Distributions are heap-allocated trait objects:
//! `Box<dyn Distribution<Output = T>>`. The associated `Output` type
//! lets each call site keep its return type — an NP length distribution
//! returns `i64`, a base distribution returns `u8`, an allele-pool
//! distribution returns an `AlleleId`. Trait-object dispatch costs
//! one indirection per `sample` call (~5–10 ns), dwarfed by the
//! per-pass overhead from D2 — acceptable.
//!
//! ## Scope
//!
//! The trait + concrete implementations:
//!
//! - [`UniformBase`] — uniform over `{b'A', b'C', b'G', b'T'}`. Used
//!   by sampling passes when no empirical model is configured.
//! - [`UniformInt`] — uniform over a half-open integer range
//!   `[min, max)`. Foundation for trim / NP-length sampling.
//! - [`EmpiricalLengthDist`] — categorical histogram over integer values.
//! - [`AllelePoolDist`] — categorical distribution over `AlleleId`s.
//!
//! Plus the filtered-sampling primitives in [`filtered`].

use crate::rng::Rng;

// ──────────────────────────────────────────────────────────────────
// Trait
// ──────────────────────────────────────────────────────────────────

/// A typed source of random samples.
///
/// Implementors expose `Output = T` so call sites preserve the natural
/// return type. The `Box<dyn Distribution<Output = T>>` form is
/// dyn-compatible because `Output` is fixed at the box site.
pub trait Distribution {
    /// The natural output type of this distribution.
    type Output;

    /// Draw one sample from `rng`. Implementors must consume a
    /// deterministic number of RNG words per call so reproducibility
    /// (D11) holds — same seed + same call sequence ⇒ same output.
    fn sample(&self, rng: &mut Rng) -> Self::Output;

    /// Optional: enumerate the discrete support of this distribution
    /// as `(value, weight)` pairs.
    ///
    /// Used by **constraint-aware sampling** (D.6): a sampling pass
    /// asks the contract set which candidate values are admissible,
    /// filters the support to those, and renormalizes for sampling.
    /// This is what makes "no retries, just sample admissible
    /// values" work — filter-before-sample rather than
    /// retry-and-reject.
    ///
    /// `None` (the default) means the support is too large or
    /// continuous to enumerate practically. Constraint-aware passes
    /// cannot safely narrow such a distribution: in strict mode they
    /// surface `SupportUnavailable`; in permissive mode they use the
    /// pass's explicit no-op / sentinel behavior. Concrete
    /// categorical distributions (`EmpiricalLengthDist`, `UniformBase`,
    /// `UniformInt` for small ranges, `AllelePoolDist`) override this
    /// to return their full support.
    ///
    /// Weights need not be normalized; the caller renormalizes.
    fn support(&self) -> Option<Vec<(Self::Output, f64)>> {
        None
    }
}

// ──────────────────────────────────────────────────────────────────
// Submodules
// ──────────────────────────────────────────────────────────────────

mod allele_pool;
mod categorical_base;
mod empirical;
mod filtered;
mod integer;
mod uniform;

pub use allele_pool::AllelePoolDist;
pub use categorical_base::CategoricalBase;
pub use empirical::EmpiricalLengthDist;
pub use filtered::{
    sample_base_with_admit_mask, sample_filtered_result, sample_filtered_with_policy, EmptySupport,
    FilteredSampleError,
};
pub use integer::UniformInt;
pub use uniform::UniformBase;

#[cfg(test)]
mod tests;
