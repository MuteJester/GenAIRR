//! Metadata-truth test suite.
//!
//! Each submodule targets one class of metadata-correctness invariant.
//! See `README.md` for the catalog.

pub mod common;

mod airr_projection;
mod corruption_stack;
mod cross_segment;
mod indel_effects;
mod junction_anchor;
mod mutation_effects;
mod np_extension;
mod trim_ambiguity;
