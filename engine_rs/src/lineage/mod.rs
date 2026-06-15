//! Clonal lineage simulation: grow a real mutation tree from one founder
//! `Simulation` via a generation-synchronous birth–death process.

pub mod tree;
pub mod poisson;
pub mod branching;
pub mod sampling;
pub mod family;

pub use tree::{LineageNode, LineageTree};
pub use branching::{grow_lineage, grow_topology, BranchingParams};
pub use sampling::sample_and_collapse;
pub use family::simulate_family;
