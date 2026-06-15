//! Generation-synchronous clonal branching loop.

use crate::ir::Simulation;
use crate::pass::{Pass, PassContext};
use crate::rng::Rng;
use crate::trace::Trace;

use super::affinity::{sim_to_aa, AffinityModel};
use super::poisson::poisson_sample;
use super::tree::{LineageNode, LineageTree};

/// Parameters controlling one clonal family's growth (neutral mode).
#[derive(Clone, Debug)]
pub struct BranchingParams {
    /// Base expected offspring per cell per generation (neutral λ).
    pub lambda_base: f64,
    /// Expected mutations per cell division. NOTE: currently inert — per-division
    /// mutation count is driven by the injected `Pass` mutator; this field is
    /// reserved for a later task and is not read yet.
    pub lambda_mut: f64,
    /// Maximum number of generations to grow.
    pub max_generations: u32,
    /// Carrying capacity: the live population per generation is capped at this,
    /// and it is the pivot for logistic offspring-rate damping (the effective
    /// rate falls to zero as the live set approaches `n_max`).
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

/// Apply one division's mutations: run `mutator` against `parent` with a fresh,
/// deterministic per-division RNG seeded by `child_seed`. Returns the mutated child.
fn mutate_child(parent: &Simulation, mutator: &dyn Pass, child_seed: u64) -> Simulation {
    let mut trace = Trace::new();
    let mut rng = Rng::new(child_seed);
    let mut ctx = PassContext {
        replay_cursor: None,
        trace: &mut trace,
        rng: &mut rng,
        pass_index: 0,
        refdata: None,
        contracts: None,
        feasibility: None,
        reference_index: None,
        event_log_sink: None,
    };
    mutator.execute(parent, &mut ctx)
}

/// Core generation-synchronous growth. With `mutator = None`, children are exact
/// clones (and NO extra RNG is consumed). With `Some(m)`, each division draws a
/// deterministic sub-seed and applies `m`. With `affinity = Some(model)`, each
/// cell's offspring rate is modulated by the model's fitness; `None` leaves the
/// rate unchanged (byte-identical to the pre-affinity path). Returns
/// (tree, peak_live_population, sims_arena) where sims_arena[node.id] is the
/// Simulation for that node.
fn grow_core(
    founder: &Simulation,
    params: &BranchingParams,
    mutator: Option<&dyn Pass>,
    affinity: Option<&AffinityModel>,
) -> (LineageTree, usize, Vec<Simulation>) {
    let mut nodes: Vec<LineageNode> = Vec::new();
    let mut sims: Vec<Simulation> = Vec::new();
    let mut rng = Rng::new(params.seed);

    let root_affinity = affinity
        .map(|m| m.affinity_value(&sim_to_aa(founder)))
        .unwrap_or(0.0);

    nodes.push(LineageNode {
        id: 0,
        parent_id: None,
        generation: 0,
        genotype: genotype_of(founder),
        mutations_from_parent: 0,
        affinity: root_affinity,
        abundance: 0,
        observed: false,
    });
    sims.push(founder.clone());

    let mut live: Vec<u32> = vec![0];
    let mut next_id: u32 = 1;
    let mut peak_live: usize = live.len();

    for gen in 1..=params.max_generations {
        if live.is_empty() {
            break; // lineage went extinct (every cell drew 0 offspring)
        }
        // Logistic carrying-capacity damping: as the live population approaches
        // n_max, the effective offspring rate falls toward zero (plateau).
        let live_frac = (live.len() as f64) / (params.n_max as f64);
        let eff_lambda = params.lambda_base * (1.0 - live_frac).max(0.0);

        let mut next_live: Vec<u32> = Vec::new();
        'generation: for &parent_id in &live {
            let cell_lambda = match affinity {
                Some(m) => eff_lambda * m.fitness_from_affinity(nodes[parent_id as usize].affinity),
                None => eff_lambda,
            };
            let k = poisson_sample(&mut rng, cell_lambda);
            for _ in 0..k {
                // Hard cap: a Poisson draw can still overshoot near saturation.
                if next_live.len() >= params.n_max as usize {
                    break 'generation;
                }
                let parent_mut_count = sims[parent_id as usize].mutation_count;
                let child_sim = match mutator {
                    Some(m) => {
                        let child_seed = rng.next_u64();
                        mutate_child(&sims[parent_id as usize], m, child_seed)
                    }
                    None => sims[parent_id as usize].clone(),
                };
                let muts = child_sim
                    .mutation_count
                    .saturating_sub(parent_mut_count);
                let child_affinity = affinity
                    .map(|m| m.affinity_value(&sim_to_aa(&child_sim)))
                    .unwrap_or(0.0);
                nodes.push(LineageNode {
                    id: next_id,
                    parent_id: Some(parent_id),
                    generation: gen,
                    genotype: genotype_of(&child_sim),
                    mutations_from_parent: muts,
                    affinity: child_affinity,
                    abundance: 0,
                    observed: false,
                });
                sims.push(child_sim);
                next_live.push(next_id);
                next_id += 1;
            }
        }
        live = next_live;
        if live.len() > peak_live {
            peak_live = live.len();
        }
    }

    (LineageTree { nodes }, peak_live, sims)
}

/// Grow a full clonal lineage with per-division mutation via `mutator`.
/// Deterministic for `params.seed`.
pub fn grow_lineage(
    founder: &Simulation,
    params: &BranchingParams,
    mutator: &dyn Pass,
) -> LineageTree {
    grow_core(founder, params, Some(mutator), None).0
}

/// Grow the lineage TOPOLOGY only (children are exact clones of their parent
/// `Simulation`; no mutation — a later task layers mutation on top).
pub fn grow_topology(founder: &Simulation, params: &BranchingParams) -> LineageTree {
    grow_core(founder, params, None, None).0
}

/// Grow a clonal lineage with per-division mutation AND affinity selection.
/// Deterministic for `params.seed`.
pub fn grow_lineage_with_affinity(
    founder: &Simulation,
    params: &BranchingParams,
    mutator: &dyn Pass,
    model: &AffinityModel,
) -> LineageTree {
    grow_core(founder, params, Some(mutator), Some(model)).0
}

/// Grow + mutate a lineage and ALSO return the per-node `Simulation` arena
/// (index == node id), for building per-node AIRR `Outcome`s.
pub fn grow_lineage_retaining_sims(
    founder: &Simulation,
    params: &BranchingParams,
    mutator: &dyn Pass,
    affinity: Option<&AffinityModel>,
) -> (LineageTree, Vec<Simulation>) {
    let (tree, _peak, sims) = grow_core(founder, params, Some(mutator), affinity);
    (tree, sims)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::dist::{EmpiricalLengthDist, UniformBase};
    use crate::ir::{Nucleotide, NucHandle, Region, Segment, Simulation};
    use crate::passes::UniformMutationPass;

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

    #[test]
    fn carrying_capacity_bounds_live_population() {
        let params = BranchingParams {
            lambda_base: 4.0,      // would explode unbounded
            lambda_mut: 0.0,
            max_generations: 40,
            n_max: 200,
            n_sample: 10,
            seed: 7,
        };
        let (_tree, peak_live, _sims) = grow_core(&founder(), &params, None, None);
        // hard cap: live population never exceeds n_max
        assert!(peak_live <= params.n_max as usize,
            "peak live {peak_live} exceeded n_max {}", params.n_max);
        // damping plateau (equilibrium ~ (1 - 1/lambda_base) * n_max = 150) is
        // comfortably above half capacity — confirms it grew, not died early
        assert!(peak_live > (params.n_max as usize) / 2,
            "peak live {peak_live} did not approach capacity");
    }

    /// A mutator that applies exactly 2 substitutions per division.
    fn two_mut_mutator() -> UniformMutationPass {
        UniformMutationPass::new(
            Box::new(EmpiricalLengthDist::from_pairs(vec![(2, 1.0)])),
            Box::new(UniformBase),
        )
    }

    #[test]
    fn children_accumulate_mutations_and_branch_lengths() {
        let params = BranchingParams {
            lambda_base: 1.2, lambda_mut: 0.0, max_generations: 4,
            n_max: 500, n_sample: 10, seed: 11,
        };
        let mutator = two_mut_mutator();
        let tree = grow_lineage(&founder(), &params, &mutator);

        assert_eq!(tree.root().mutations_from_parent, 0);

        let mut saw_mutated_child = false;
        for n in &tree.nodes {
            if n.parent_id.is_some() {
                // <= 2 (not == 2): on the 4-base "AAAA" founder the two uniform
                // substitutions can draw the same position, yielding 1 net change.
                assert!(n.mutations_from_parent <= 2);
                if n.mutations_from_parent > 0 {
                    saw_mutated_child = true;
                    let parent = tree.get(n.parent_id.unwrap()).unwrap();
                    assert_ne!(n.genotype, parent.genotype);
                }
            }
        }
        assert!(saw_mutated_child, "no child accumulated any mutation");
    }

    #[test]
    fn mutated_growth_is_deterministic() {
        let params = BranchingParams {
            lambda_base: 1.2, lambda_mut: 0.0, max_generations: 4,
            n_max: 500, n_sample: 10, seed: 11,
        };
        let a = grow_lineage(&founder(), &params, &two_mut_mutator());
        let b = grow_lineage(&founder(), &params, &two_mut_mutator());
        assert_eq!(a.len(), b.len());
        for (x, y) in a.nodes.iter().zip(b.nodes.iter()) {
            assert_eq!(x.genotype, y.genotype);
            assert_eq!(x.mutations_from_parent, y.mutations_from_parent);
        }
    }

    #[test]
    fn affinity_run_populates_node_affinities_and_root_matches() {
        use crate::lineage::affinity::{sim_to_aa, AffinityModel};
        let f = founder();
        let founder_aa = sim_to_aa(&f);
        let w = vec![1.0; founder_aa.len().max(1)];
        let model = AffinityModel::new(b"W".to_vec(), w, 1.0, 1.0, &founder_aa);
        let params = BranchingParams { lambda_base:1.2, lambda_mut:0.0, max_generations:4, n_max:500, n_sample:10, seed:11 };
        let tree = grow_lineage_with_affinity(&f, &params, &two_mut_mutator(), &model);
        assert!((tree.root().affinity - model.affinity_value(&founder_aa)).abs() < 1e-9);
        for n in &tree.nodes {
            assert!(n.affinity > 0.0 && n.affinity <= 1.0 + 1e-9, "affinity out of range: {}", n.affinity);
        }
    }

    #[test]
    fn affinity_strength_zero_matches_neutral_topology() {
        use crate::lineage::affinity::{sim_to_aa, AffinityModel};
        let f = founder();
        let founder_aa = sim_to_aa(&f);
        let w = vec![1.0; founder_aa.len().max(1)];
        // selection_strength = 0 => fitness identically 1 => identical topology/genotypes to neutral
        let model = AffinityModel::new(b"W".to_vec(), w, 1.0, 0.0, &founder_aa);
        let params = BranchingParams { lambda_base:1.2, lambda_mut:0.0, max_generations:4, n_max:500, n_sample:10, seed:11 };
        let with_aff = grow_lineage_with_affinity(&f, &params, &two_mut_mutator(), &model);
        let neutral = grow_lineage(&f, &params, &two_mut_mutator());
        assert_eq!(with_aff.len(), neutral.len());
        for (a, b) in with_aff.nodes.iter().zip(neutral.nodes.iter()) {
            assert_eq!(a.genotype, b.genotype);
            assert_eq!(a.parent_id, b.parent_id);
            assert_eq!(a.generation, b.generation);
            assert_eq!(a.mutations_from_parent, b.mutations_from_parent);
        }
    }

    #[test]
    fn affinity_growth_is_deterministic() {
        use crate::lineage::affinity::{sim_to_aa, AffinityModel};
        let f = founder();
        let founder_aa = sim_to_aa(&f);
        let mk = || AffinityModel::new(b"W".to_vec(), vec![1.0; founder_aa.len().max(1)], 1.0, 2.0, &founder_aa);
        let params = BranchingParams { lambda_base:1.5, lambda_mut:0.0, max_generations:5, n_max:500, n_sample:10, seed:21 };
        let a = grow_lineage_with_affinity(&f, &params, &two_mut_mutator(), &mk());
        let b = grow_lineage_with_affinity(&f, &params, &two_mut_mutator(), &mk());
        assert_eq!(a.len(), b.len());
        for (x, y) in a.nodes.iter().zip(b.nodes.iter()) {
            assert_eq!(x.genotype, y.genotype);
            assert!((x.affinity - y.affinity).abs() < 1e-12);
        }
    }
}
