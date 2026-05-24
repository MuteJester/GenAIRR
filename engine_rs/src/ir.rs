//! Intermediate-representation (IR) type definitions for the
//! simulation engine.
//!
//! This module defines the typed data shapes from §3 of the design
//! document. It deliberately contains **no biological logic** — only
//! the structural definitions of the entities. Logic lives in
//! sibling modules (passes, contracts, queries).
//!
//! ## Architectural commitments reflected here
//!
//! - **Arena allocation (D-§9):** `NucleotidePool` is a flat arena.
//!   Cross-entity references are typed `u32` newtype handles, not pointers.
//! - **Persistent IR (D1):** all top-level entities derive `Clone` and have
//!   no interior mutability.
//! - **Entity-attached metadata (D5):** derived state lives on the entity
//!   that owns it (live-call state, dirty log, mutation count). Codon-rail
//!   data is *not* stored on `Region` — it's a pure on-demand projection
//!   produced by [`compute_codon_rail`] at the Python boundary; see the
//!   struct comment in [`region`].
//!
//! ## Performance / cost model
//!
//! D1 commits to a persistent IR contract: every `with_*` method takes
//! `&self`, returns a new value, and leaves the receiver intact. The current
//! implementation uses `Arc` plus copy-on-write for the nucleotide pool,
//! sequence, and live-call state so ordinary `Clone` calls are cheap while
//! preserving the public persistent API.

pub(crate) mod builder;
pub(crate) mod event_log_observer;
mod handle;
mod nucleotide;
mod per_segment;
mod pool;
mod region;
mod segment;
mod sequence;
mod simulation;

pub use crate::codon::{translate_codon, AMINO_AMBIGUOUS, AMINO_STOP};
pub use builder::SimulationBuilder;
pub use handle::{NucHandle, RegionHandle};
pub use nucleotide::{encode_base, flag, GermlinePos, NucFlags, Nucleotide};
pub use per_segment::PerSegment;
pub use pool::NucleotidePool;
pub use region::{compute_codon_rail, CodonRail, Region};
pub use segment::Segment;
pub use sequence::Sequence;
pub use simulation::Simulation;

#[cfg(test)]
mod tests;
