//! One-call clonal family simulation: grow a lineage and sample it.

use crate::ir::Simulation;
use crate::pass::Pass;
use crate::rng::Rng;

use super::affinity::AffinityModel;
use super::branching::{grow_lineage_retaining_sims, BranchingParams};
use super::sampling::sample_and_collapse;
use super::tree::LineageTree;

/// Salt mixed into `params.seed` to give sampling an independent RNG stream
/// from growth, while staying deterministic for a given seed.
const SAMPLE_SEED_SALT: u64 = 0x5341_4D50_4C45_0001; // "SAMPLE\0\1"

/// Grow a clonal family from `founder` and sample observed cells.
///
/// Deterministic for `params.seed`. Growth and sampling use separate RNG
/// streams (the sampling stream is `params.seed ^ SAMPLE_SEED_SALT`), so the
/// whole family reproduces byte-for-byte across runs.
pub fn simulate_family(
    founder: &Simulation,
    params: &BranchingParams,
    mutator: &dyn Pass,
) -> LineageTree {
    let (mut tree, _sims, live) = grow_lineage_retaining_sims(founder, params, mutator, None);
    let mut sample_rng = Rng::new(params.seed ^ SAMPLE_SEED_SALT);
    sample_and_collapse(&mut tree, &live, params.n_sample, &mut sample_rng);
    tree
}

/// Grow + sample a clonal family with affinity selection.
pub fn simulate_family_with_affinity(
    founder: &Simulation,
    params: &BranchingParams,
    mutator: &dyn Pass,
    model: &AffinityModel,
) -> LineageTree {
    let (mut tree, _sims, live) =
        grow_lineage_retaining_sims(founder, params, mutator, Some(model));
    let mut sample_rng = Rng::new(params.seed ^ SAMPLE_SEED_SALT);
    sample_and_collapse(&mut tree, &live, params.n_sample, &mut sample_rng);
    tree
}

/// Grow + sample a family, returning the tree AND the per-node Simulation arena.
/// Index in the arena equals the node id (arena[node.id] is the Simulation for
/// that node). Only observed (sampled) nodes are useful for AIRR projection; the
/// arena is full-tree and never trimmed.
pub fn simulate_family_sims(
    founder: &Simulation,
    params: &BranchingParams,
    mutator: &dyn Pass,
    affinity: Option<&AffinityModel>,
) -> (LineageTree, Vec<Simulation>) {
    let (mut tree, sims, live) = grow_lineage_retaining_sims(founder, params, mutator, affinity);
    let mut sample_rng = Rng::new(params.seed ^ SAMPLE_SEED_SALT);
    sample_and_collapse(&mut tree, &live, params.n_sample, &mut sample_rng);
    (tree, sims)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::dist::{EmpiricalLengthDist, UniformBase};
    use crate::ir::{Nucleotide, NucHandle, Region, Segment, Simulation};
    use crate::lineage::BranchingParams;
    use crate::passes::UniformMutationPass;

    fn founder() -> Simulation {
        let mut sim = Simulation::new();
        for (i, b) in b"AAAAAAAA".iter().enumerate() {
            let (next, _) = sim.with_nucleotide_pushed(
                Nucleotide::germline(*b, i as u16, Segment::V));
            sim = next;
        }
        sim.with_region_added(Region::new(Segment::V, NucHandle::new(0), NucHandle::new(8)))
    }

    #[test]
    fn simulate_family_returns_valid_sampled_tree() {
        let params = BranchingParams {
            lambda_base: 1.5, lambda_mut: 0.0, max_generations: 6,
            n_max: 300, n_sample: 20, seed: 2024,
        };
        let mutator = UniformMutationPass::new(
            Box::new(EmpiricalLengthDist::from_pairs(vec![(1, 1.0)])),
            Box::new(UniformBase),
        );
        let tree = simulate_family(&founder(), &params, &mutator);

        assert!(tree.validate().is_ok(), "invalid tree: {:?}", tree.validate());
        let total_abundance: u32 = tree.nodes.iter().map(|n| n.abundance).sum();
        // abundance equals n_sample unless the family went extinct early
        assert!(total_abundance == params.n_sample || total_abundance == 0);

        let tree2 = simulate_family(&founder(), &params,
            &UniformMutationPass::new(
                Box::new(EmpiricalLengthDist::from_pairs(vec![(1, 1.0)])),
                Box::new(UniformBase)));
        assert_eq!(tree.len(), tree2.len());
    }

    #[test]
    fn extinct_founder_yields_no_observed_cells() {
        // lambda_base = 0 => founder never divides => the living final population
        // is empty (extinct). An extinct clone must NOT be sampled: no observed node.
        let params = BranchingParams {
            lambda_base: 0.0, lambda_mut: 0.0, max_generations: 6,
            n_max: 300, n_sample: 20, seed: 2024,
        };
        let mutator = UniformMutationPass::new(
            Box::new(EmpiricalLengthDist::from_pairs(vec![(1, 1.0)])),
            Box::new(UniformBase),
        );
        let (tree, _sims) = simulate_family_sims(&founder(), &params, &mutator, None);
        assert!(
            tree.nodes.iter().all(|n| !n.observed && n.abundance == 0),
            "extinct clone was sampled: {} observed nodes",
            tree.nodes.iter().filter(|n| n.observed).count()
        );
    }

    #[test]
    fn observed_cells_are_all_at_the_final_living_generation() {
        // In a grown family, every OBSERVED node must be a member of the final
        // living generation (the maximum generation reached) — never an
        // early-terminal/dead cell that drew 0 offspring at a lower generation.
        // A single founder can go extinct, so scan seeds for one that grows.
        let mutator = || UniformMutationPass::new(
            Box::new(EmpiricalLengthDist::from_pairs(vec![(1, 1.0)])),
            Box::new(UniformBase),
        );
        let mut grew = false;
        for seed in 0..50u64 {
            let params = BranchingParams {
                lambda_base: 1.6, lambda_mut: 0.0, max_generations: 6,
                n_max: 300, n_sample: 20, seed,
            };
            let (tree, _sims) = simulate_family_sims(&founder(), &params, &mutator(), None);
            let observed: Vec<&_> = tree.nodes.iter().filter(|n| n.observed).collect();
            if observed.is_empty() {
                continue; // extinct for this seed
            }
            grew = true;
            let max_gen = tree.nodes.iter().map(|n| n.generation).max().unwrap();
            for n in &observed {
                assert_eq!(
                    n.generation, max_gen,
                    "observed node {} at generation {} is not at the final generation {} (seed {seed})",
                    n.id, n.generation, max_gen
                );
            }
        }
        assert!(grew, "no seed in 0..50 produced a surviving family");
    }

    #[test]
    fn simulate_family_sims_arena_length_matches_tree_and_observed_nodes_have_nonempty_pool() {
        let params = BranchingParams {
            lambda_base: 1.5, lambda_mut: 0.0, max_generations: 6,
            n_max: 300, n_sample: 20, seed: 9999,
        };
        let mutator = UniformMutationPass::new(
            Box::new(EmpiricalLengthDist::from_pairs(vec![(1, 1.0)])),
            Box::new(UniformBase),
        );
        let (tree, sims) = simulate_family_sims(&founder(), &params, &mutator, None);

        // Arena length must equal tree length (one entry per node).
        assert_eq!(sims.len(), tree.len(),
            "arena len {} != tree len {}", sims.len(), tree.len());

        // Every observed node's sim must have a non-empty pool.
        for node in tree.nodes.iter().filter(|n| n.observed) {
            let sim = &sims[node.id as usize];
            assert!(!sim.pool.is_empty(),
                "observed node {} has empty pool", node.id);
        }
    }
}
