//! One-call clonal family simulation: grow a lineage and sample it.

use crate::ir::Simulation;
use crate::pass::Pass;
use crate::rng::Rng;

use super::affinity::AffinityModel;
use super::branching::{grow_lineage, grow_lineage_with_affinity, BranchingParams};
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
    let mut tree = grow_lineage(founder, params, mutator);
    let mut sample_rng = Rng::new(params.seed ^ SAMPLE_SEED_SALT);
    sample_and_collapse(&mut tree, params.n_sample, &mut sample_rng);
    tree
}

/// Grow + sample a clonal family with affinity selection.
pub fn simulate_family_with_affinity(
    founder: &Simulation,
    params: &BranchingParams,
    mutator: &dyn Pass,
    model: &AffinityModel,
) -> LineageTree {
    let mut tree = grow_lineage_with_affinity(founder, params, mutator, model);
    let mut sample_rng = Rng::new(params.seed ^ SAMPLE_SEED_SALT);
    sample_and_collapse(&mut tree, params.n_sample, &mut sample_rng);
    tree
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
}
