//! Clonal lineage simulation: grow a real mutation tree from one founder
//! `Simulation` via a generation-synchronous birth–death process.
//!
//! ```ignore
//! use genairr_engine::lineage::{simulate_family, BranchingParams};
//! use genairr_engine::passes::UniformMutationPass;
//! // build a founder Simulation + a mutator (S5F in production), then:
//! let params = BranchingParams {
//!     lambda_base: 1.5, lambda_mut: 0.0, max_generations: 12,
//!     n_max: 1000, n_sample: 60, seed: 42,
//! };
//! let tree = simulate_family(&founder, &params, &mutator);
//! assert!(tree.validate().is_ok());
//! ```
//!
//! Neutral mode only (no affinity selection — a later plan). The mutator is any
//! `Pass`; production wires the S5F pass, tests use `UniformMutationPass` (no
//! reference cartridge required).

pub mod tree;
pub mod poisson;
pub mod branching;
pub mod sampling;
pub mod family;

pub use tree::{LineageNode, LineageTree};
pub use branching::{grow_lineage, grow_topology, BranchingParams};
pub use sampling::sample_and_collapse;
pub use family::simulate_family;
