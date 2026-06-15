//! Generation-synchronous clonal branching loop.

use crate::ir::Simulation;
use crate::rng::Rng;

use super::poisson::poisson_sample;
use super::tree::{LineageNode, LineageTree};

/// Parameters controlling one clonal family's growth (neutral mode).
#[derive(Clone, Debug)]
pub struct BranchingParams {
    /// Base expected offspring per cell per generation (neutral λ).
    pub lambda_base: f64,
    /// Expected mutations introduced per cell division (used in a later task).
    pub lambda_mut: f64,
    /// Maximum number of generations to grow.
    pub max_generations: u32,
    /// Population cap; growth stops adding once the live set reaches this.
    pub n_max: u32,
    /// Number of cells to sample at the end (used in a later task).
    pub n_sample: u32,
    /// Master seed for the whole family (deterministic).
    pub seed: u64,
}

/// Extract a node's genotype (pool bases) from its `Simulation`.
fn genotype_of(sim: &Simulation) -> Vec<u8> {
    sim.pool.as_slice().iter().map(|n| n.base).collect()
}

/// Grow the lineage TOPOLOGY only: children are exact clones of their parent
/// `Simulation` (no mutation). Returns the tree; live cells are whatever
/// remained at the final generation. A later task layers mutation on top.
pub fn grow_topology(founder: &Simulation, params: &BranchingParams) -> LineageTree {
    let mut nodes: Vec<LineageNode> = Vec::new();
    let mut sims: Vec<Simulation> = Vec::new();
    let mut rng = Rng::new(params.seed);

    nodes.push(LineageNode {
        id: 0,
        parent_id: None,
        generation: 0,
        genotype: genotype_of(founder),
        mutations_from_parent: 0,
        abundance: 0,
        observed: false,
    });
    sims.push(founder.clone());

    let mut live: Vec<u32> = vec![0];
    let mut next_id: u32 = 1;

    for gen in 1..=params.max_generations {
        if live.is_empty() {
            break;
        }
        let mut next_live: Vec<u32> = Vec::new();
        for &parent_id in &live {
            let k = poisson_sample(&mut rng, params.lambda_base);
            for _ in 0..k {
                if next_id >= params.n_max {
                    break;
                }
                let parent_sim = sims[parent_id as usize].clone();
                let child_sim = parent_sim; // exact clone for now; mutated in a later task
                nodes.push(LineageNode {
                    id: next_id,
                    parent_id: Some(parent_id),
                    generation: gen,
                    genotype: genotype_of(&child_sim),
                    mutations_from_parent: 0,
                    abundance: 0,
                    observed: false,
                });
                sims.push(child_sim);
                next_live.push(next_id);
                next_id += 1;
            }
        }
        live = next_live;
    }

    LineageTree { nodes }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::ir::{Nucleotide, Region, Segment, NucHandle};
    use crate::ir::Simulation;

    /// Minimal founder: a 4-base V region "AAAA".
    fn founder() -> Simulation {
        let mut sim = Simulation::new();
        for (i, b) in b"AAAA".iter().enumerate() {
            let (next, _) = sim.with_nucleotide_pushed(
                Nucleotide::germline(*b, i as u16, Segment::V),
            );
            sim = next;
        }
        sim.with_region_added(Region::new(Segment::V, NucHandle::new(0), NucHandle::new(4)))
    }

    fn neutral_params() -> BranchingParams {
        BranchingParams {
            lambda_base: 1.5,
            lambda_mut: 0.0,
            max_generations: 5,
            n_max: 1000,
            n_sample: 10,
            seed: 42,
        }
    }

    #[test]
    fn grows_root_and_children_with_monotonic_generations() {
        let tree = grow_topology(&founder(), &neutral_params());
        assert_eq!(tree.nodes.iter().filter(|n| n.parent_id.is_none()).count(), 1);
        let root = tree.root();
        assert_eq!(root.id, 0);
        assert_eq!(root.generation, 0);
        for n in &tree.nodes {
            if let Some(pid) = n.parent_id {
                let parent = tree.get(pid).unwrap();
                assert_eq!(n.generation, parent.generation + 1);
            }
        }
        assert!(tree.len() > 1, "tree did not grow: {} nodes", tree.len());
    }

    #[test]
    fn growth_is_deterministic_for_a_seed() {
        let a = grow_topology(&founder(), &neutral_params());
        let b = grow_topology(&founder(), &neutral_params());
        assert_eq!(a.len(), b.len());
        for (x, y) in a.nodes.iter().zip(b.nodes.iter()) {
            assert_eq!(x.id, y.id);
            assert_eq!(x.parent_id, y.parent_id);
            assert_eq!(x.generation, y.generation);
        }
    }
}
