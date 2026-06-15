//! Clonal lineage simulation: grow a real mutation tree from one founder
//! `Simulation` via a generation-synchronous birth–death process.

pub mod tree;
pub mod poisson;
pub mod branching;

pub use tree::{LineageNode, LineageTree};
pub use branching::{grow_lineage, grow_topology, BranchingParams};
