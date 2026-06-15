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
//! The mutator is any `Pass`; production wires the S5F pass, tests use
//! `UniformMutationPass` (no reference cartridge required). Affinity-driven
//! selection (BLOSUM-weighted distance to a target antigen → fitness-modulated
//! offspring) lives in [`affinity`] and is used via [`simulate_family_with_affinity`].
//!
//! Ground-truth export (Newick / FASTA / node-table TSV) lives in [`export`].
//! Heavy-tailed clone-size distributions + repertoire composition (the TCR core
//! and the singleton tail of BCR repertoires) live in [`clone_size`].

pub mod tree;
pub mod poisson;
pub mod branching;
pub mod sampling;
pub mod family;
pub mod export;
pub mod affinity;
pub mod clone_size;

pub use clone_size::{sample_clone_size, sample_repertoire_sizes, CloneSizeDist};

pub use affinity::{sim_to_aa, AffinityModel};
pub use tree::{LineageNode, LineageTree};
pub use branching::{grow_lineage, grow_lineage_with_affinity, grow_topology, BranchingParams};
pub use sampling::sample_and_collapse;
pub use family::{simulate_family, simulate_family_with_affinity};
pub use export::{to_fasta, to_newick, to_node_table_tsv};
