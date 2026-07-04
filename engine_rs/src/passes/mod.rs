//! Concrete pass implementations.
//!
//! Each biology operation lives in its own submodule; this module
//! re-exports the pass types so external code keeps the flat
//! `crate::passes::FooPass` import surface.
//!
//! Submodule layout:
//! - [`echo`] — `EchoPass`, the deterministic transform reference.
//! - [`sample_base`] — `SampleBasePass`, the sampling reference.
//! - [`sample_allele`] — V/D/J allele sampling.
//! - [`trim`] — recombination trim sampling.
//! - [`assemble_segment`] — copy a germline allele slice into
//!   the pool.
//! - [`generate_np`] — TdT-like N-nucleotide region generation
//!  .
//! - [`mutate`] — SHM passes (uniform + S5F).
//! - [`corrupt`] — observation-stage perturbations (PCR error,
//!   quality error, contamination, indels).

pub mod assemble_segment;
pub mod corrupt;
pub(crate) mod count_source;
pub mod echo;
pub mod generate_np;
pub mod invert_d;
pub mod mutate;
pub(crate) mod mutation_transaction;
pub mod p_addition;
pub mod paired_end;
pub(crate) mod paramsig;
pub mod receptor_revision;
pub mod sample_allele;
pub mod sample_base;
pub mod sample_genotype;
pub mod trim;

#[cfg(test)]
pub(crate) mod test_support;

pub use assemble_segment::AssembleSegmentPass;
pub use corrupt::{
    ContaminantPass, EndLossPass, IndelPass, LossEnd, NCorruptionPass, PCRErrorPass,
    QualityErrorPass, RevCompPass,
};
pub use echo::EchoPass;
pub use generate_np::GenerateNPPass;
pub use invert_d::InvertDPass;
pub use mutate::{S5FMutationPass, UniformMutationPass};
pub use p_addition::PAdditionPass;
pub use paired_end::{PairedEndLayoutSpec, PairedEndSamplingPass};
pub use receptor_revision::ReceptorRevisionPass;
pub use sample_allele::SampleAllelePass;
pub use sample_base::SampleBasePass;
pub use trim::TrimPass;

// ──────────────────────────────────────────────────────────────────
// Cross-cutting integration tests (D.6 productive bundle)
// ──────────────────────────────────────────────────────────────────

#[cfg(test)]
mod tests;
