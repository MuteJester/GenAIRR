//! Mutation passes — SHM models that substitute pool bases.
//!
//! Currently houses two mutation models:
//! - [`UniformMutationPass`] — position-independent uniform-base
//!   substitution.
//! - [`S5FMutationPass`] — context-dependent S5F kernel substitution
//!   (Yaari et al. 2013).

pub mod s5f;
pub mod segment_rates;
pub mod uniform;
pub mod v_subregion_rates;

pub use crate::passes::count_source::CountSource;
pub use s5f::S5FMutationPass;
pub use segment_rates::{segment_at_position, SegmentRateWeights};
pub use uniform::UniformMutationPass;
pub use v_subregion_rates::{v_subregion_at_position, VSubregionRateWeights};
