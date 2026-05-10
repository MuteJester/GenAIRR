//! Mutation passes — SHM models that substitute pool bases.
//!
//! Currently houses two mutation models:
//! - [`UniformMutationPass`] — position-independent uniform-base
//!   substitution.
//! - [`S5FMutationPass`] — context-dependent S5F kernel substitution
//!   (Yaari et al. 2013).

pub mod s5f;
pub mod uniform;

pub use s5f::S5FMutationPass;
pub use uniform::UniformMutationPass;
