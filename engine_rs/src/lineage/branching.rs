//! Generation-synchronous clonal branching loop.

use crate::contract::{ContractSet, ContractViolation};
use crate::event::{EventRecord, StateSummary, TraceSpan};
use crate::ir::Simulation;
use crate::pass::{Outcome, Pass, PassContext, PassError};
use crate::refdata::RefDataConfig;
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
    /// Carrying capacity: the live population a lineage settles at. It is the
    /// pivot for Beverton-Holt offspring-rate regulation -- the effective rate
    /// is `lambda_base` in an empty environment and exactly 1 at `n_max`, so
    /// the population levels off there rather than being truncated. A single
    /// generation may overshoot it; the next one is regulated back.
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

/// Execution context inherited by every division. It conditions mutation
/// support directly; sampled children are never repaired or post-filtered.
pub struct LineageMutationContext<'a> {
    pub refdata: &'a RefDataConfig,
    pub contracts: Option<&'a ContractSet>,
    pub strict: bool,
}

fn require_productive(
    sim: &Simulation,
    context: Option<&LineageMutationContext<'_>>,
) -> Result<(), PassError> {
    let Some(context) = context else {
        return Ok(());
    };
    if context
        .contracts
        .is_some_and(ContractSet::requires_productive_emission)
    {
        let assessment = crate::productivity::assess(sim, context.refdata);
        if assessment.productive != Some(true) {
            return Err(PassError::contract_violation(
                "lineage.productive",
                vec![ContractViolation::new(
                    crate::productivity::policy_name(sim),
                    crate::productivity::joined(&assessment.reasons),
                )],
            ));
        }
    }
    Ok(())
}

/// Apply one division using the same context and execution policy as the
/// founder. The unconstrained path preserves the previous RNG/call sequence.
fn mutate_child(
    parent: &Simulation,
    mutator: &dyn Pass,
    child_seed: u64,
    context: Option<&LineageMutationContext<'_>>,
    retain_edge: bool,
) -> Result<(Simulation, Option<Outcome>), PassError> {
    let mut trace = Trace::new();
    // The child's seed *is* this pass's stream, not a plan-wide one that
    // several passes draw from in turn: one division runs one mutator, and
    // each child's seed is a fresh word. There is no position in a shared
    // stream here for another pass to move, so there is nothing for
    // `rng::substream` to separate — see `compiled::execute`, which does
    // need it because a plan runs many passes against one seed.
    let mut rng = Rng::new(child_seed);
    let mut events = Vec::new();
    let mut ctx = PassContext {
        replay_cursor: None,
        trace: &mut trace,
        rng: &mut rng,
        pass_index: 0,
        refdata: context.map(|c| c.refdata),
        contracts: context.and_then(|c| c.contracts),
        feasibility: None,
        reference_index: None,
        event_log_sink: retain_edge.then_some(&mut events),
    };
    let any: &dyn std::any::Any = mutator;
    let child = if let Some(s5f) = any.downcast_ref::<crate::passes::S5FMutationPass>() {
        s5f.execute_for_lineage(parent, &mut ctx, context.is_some_and(|c| c.strict))?
    } else {
        match context {
            Some(context) => mutator.execute_with_policy(parent, &mut ctx, context.strict)?,
            None => mutator.execute(parent, &mut ctx),
        }
    };
    if let Some(context) = context.filter(|c| c.strict) {
        if let Some(contracts) = context.contracts {
            contracts
                .verify(&child, Some(context.refdata))
                .map_err(|violations| {
                    PassError::contract_violation("lineage.mutation", violations)
                })?;
        }
    }
    require_productive(&child, context)?;
    let edge = retain_edge.then(|| {
        let event = EventRecord::pass_committed(
            0,
            mutator.name(),
            mutator.effects(),
            TraceSpan::new(0, trace.len()),
            StateSummary::from_simulation(parent),
            StateSummary::from_simulation(&child),
            events,
        );
        Outcome {
            revisions: vec![parent.clone(), child.clone()],
            pass_names: vec![mutator.name().to_string()],
            trace,
            events: vec![event],
        }
    });
    Ok((child, edge))
}

/// Validate the complete founder before division, including external native
/// callers that do not enter through the Experiment compiler.
pub(crate) fn validate_duplex_founder(
    founder: &Simulation,
    context: Option<&LineageMutationContext<'_>>,
    mutator: Option<&dyn Pass>,
) -> Result<(), PassError> {
    if founder.duplex_provenance.is_none() {
        return Ok(());
    }
    let fail = |reason: String| PassError::invalid_plan_state("lineage.duplex", reason);
    let context =
        context.ok_or_else(|| fail("duplex lineage requires its bound reference data".into()))?;
    if context.contracts.is_some_and(|contracts| {
        contracts.iter().any(|c| {
            matches!(
                c.kind(),
                crate::contract::ContractKind::Custom
                    | crate::contract::ContractKind::AnchorPreserved { .. }
            )
        })
    }) {
        return Err(fail(
            "duplex lineage does not support custom or legacy AnchorPreserved contracts".into(),
        ));
    }
    if mutator.is_some_and(|m| {
        use std::any::TypeId;
        ![
            TypeId::of::<crate::passes::S5FMutationPass>(),
            TypeId::of::<crate::passes::UniformMutationPass>(),
        ]
        .contains(&m.type_id())
    }) {
        return Err(fail(
            "duplex lineage division requires an audited substitution mutator".into(),
        ));
    }
    crate::duplex_repair::observation::project_duplex_simulation(founder, context.refdata)
        .map_err(fail)?;
    Ok(())
}

/// The protein a cell expresses, for affinity scoring.
///
/// Both paths now translate the coding span continuously. The plain path used
/// to concatenate the per-region codon rail (`sim_to_aa`), which drops every
/// codon straddling a region boundary -- and every boundary in a rearrangement
/// is a junction boundary, so it dropped CDR3. A cell scored against its own
/// protein reported an affinity of 3.57e-67 instead of 1.
///
/// The reference is not optional here: `coding_span` resolves the frame from
/// the V allele's native phase, its 5' trim and its observation loss, so a
/// translation that starts at pool position zero can be in the wrong frame
/// entirely. That is why an affinity model travels with its `RefDataConfig`
/// rather than borrowing one from the mutation context -- the context also
/// selects the mutator's execution policy, and giving the neutral path one
/// just to reach its reference would move its draws.
fn affinity_aa(sim: &Simulation, refdata: &RefDataConfig) -> Result<Vec<u8>, PassError> {
    super::affinity::coding_sim_to_aa(sim, refdata)
        .map_err(|reason| PassError::invalid_plan_state("lineage.affinity", reason))
}

/// Core generation-synchronous growth. With `mutator = None`, children are exact
/// clones (and NO extra RNG is consumed). With `Some(m)`, each division draws a
/// deterministic sub-seed and applies `m`. With `affinity = Some(model)`, each
/// cell's offspring rate is modulated by the model's fitness; `None` leaves the
/// rate unchanged (byte-identical to the pre-affinity path). Returns
/// (tree, peak_live_population, sims_arena, final_live) where sims_arena[node.id]
/// is the Simulation for that node, and `final_live` is the set of node ids that
/// are ALIVE when growth stopped — the cells produced in the final completed
/// generation (every cell is replaced by its offspring each generation, so the
/// living population at sampling time is exactly this set). It is empty when the
/// family went extinct (a cell drew 0 offspring and left no progeny).
fn grow_core(
    founder: &Simulation,
    params: &BranchingParams,
    mutator: Option<&dyn Pass>,
    affinity: Option<(&AffinityModel, &RefDataConfig)>,
    context: Option<&LineageMutationContext<'_>>,
    mut replay: Option<&mut super::replay_receipt::FamilyReplayDriver<'_>>,
) -> Result<
    (
        LineageTree,
        usize,
        Vec<Simulation>,
        Vec<u32>,
        Vec<Option<Outcome>>,
    ),
    PassError,
> {
    validate_duplex_founder(founder, context, mutator)?;
    require_productive(founder, context)?;
    if let Some(driver) = replay.as_deref_mut() {
        driver
            .budget
            .charge_nodes(1)
            .map_err(super::replay_receipt::error)?;
        driver
            .budget
            .charge_bases(founder.pool.len() as u64)
            .map_err(super::replay_receipt::error)?;
    }
    let retain_edges = founder.duplex_provenance.is_some();
    let mut edges = if retain_edges { vec![None] } else { Vec::new() };
    let mut nodes: Vec<LineageNode> = Vec::new();
    let mut sims: Vec<Simulation> = Vec::new();
    let mut rng = super::streams::stream(params.seed, super::streams::GROWTH);

    let root_affinity = match affinity {
        Some((m, refdata)) => m.affinity_value(&affinity_aa(founder, refdata)?),
        None => 0.0,
    };

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
    let mut completed_generations = 0;

    for gen in 1..=params.max_generations {
        if live.is_empty() {
            break; // lineage went extinct (every cell drew 0 offspring)
        }
        // Beverton-Holt regulation. The effective offspring rate falls from
        // `lambda_base` toward 1 as the live population approaches `n_max`, and
        // below 1 above it, so `n_max` is the population a lineage settles at.
        //
        // The discrete logistic law this replaces, `lambda * (1 - N/n_max)`,
        // settled at `n_max * (1 - 1/lambda)` rather than at `n_max`; above
        // `lambda = 3` it had no stable population at all and oscillated; and
        // paired with a hard cap that landed the population on `n_max` exactly
        // it returned an effective rate of exactly zero there, so reaching
        // capacity extinguished the lineage with certainty. Beverton-Holt is
        // monotone in N, so it cannot oscillate at any lambda, and it is
        // exactly 1 at N = n_max: a population at capacity replaces itself.
        //
        // Sub-replacement growth is excluded rather than damped. It is already
        // shrinking, so there is nothing to regulate, and the denominator would
        // turn negative past `n_max / (1 - lambda)`.
        let eff_lambda = if params.lambda_base <= 1.0 {
            params.lambda_base
        } else {
            let live_frac = (live.len() as f64) / (params.n_max as f64);
            params.lambda_base / (1.0 + (params.lambda_base - 1.0) * live_frac)
        };

        let mut next_live: Vec<u32> = Vec::new();
        for &parent_id in &live {
            let cell_lambda = match affinity {
                Some((m, _)) => {
                    eff_lambda * m.fitness_from_affinity(nodes[parent_id as usize].affinity)
                }
                None => eff_lambda,
            };
            let k = if let Some(driver) = replay.as_deref_mut() {
                driver.offspring(gen, parent_id, cell_lambda, &mut rng)?
            } else {
                poisson_sample(&mut rng, cell_lambda)
            };
            for _ in 0..k {
                // No truncation. A generation may overshoot `n_max`; regulation
                // returns an effective rate below 1 there and the next
                // generation shrinks. Discarding the overshoot instead would
                // have to choose which cells die, and the cheap choice -- stop
                // at the first parent that fills the quota -- kept the broods of
                // early-indexed parents and dropped every later parent's brood
                // entire, letting vector position decide survival.
                let parent_mut_count = sims[parent_id as usize].mutation_count;
                let (child_sim, edge) = match mutator {
                    Some(m) => {
                        if let Some(driver) = replay.as_deref_mut() {
                            let parent = &sims[parent_id as usize];
                            let (child_seed, recorded) =
                                driver.begin_edge(next_id, parent_id, gen, parent, &mut rng)?;
                            let (child, edge) = if let Some(trace) = recorded {
                                let any: &dyn std::any::Any = m;
                                let s5f = any
                                    .downcast_ref::<crate::passes::S5FMutationPass>()
                                    .ok_or_else(|| {
                                        super::replay_receipt::error(
                                            "unsupported replay division mutator",
                                        )
                                    })?;
                                let edge = execute_mutation_choices(
                                    parent,
                                    &trace,
                                    s5f,
                                    context.ok_or_else(|| {
                                        super::replay_receipt::error(
                                            "missing replay division context",
                                        )
                                    })?,
                                )?;
                                (edge.final_simulation().clone(), Some(edge))
                            } else {
                                mutate_child(parent, m, child_seed, context, true)?
                            };
                            driver.commit_edge(
                                next_id,
                                parent_id,
                                gen,
                                child_seed,
                                edge.as_ref().ok_or_else(|| {
                                    super::replay_receipt::error(
                                        "division did not retain its event history",
                                    )
                                })?,
                            )?;
                            (child, edge)
                        } else {
                            let child_seed = rng.next_u64();
                            mutate_child(
                                &sims[parent_id as usize],
                                m,
                                child_seed,
                                context,
                                retain_edges,
                            )?
                        }
                    }
                    None => (sims[parent_id as usize].clone(), None),
                };
                let muts = child_sim.mutation_count.saturating_sub(parent_mut_count);
                let child_affinity = match affinity {
                    Some((m, refdata)) => m.affinity_value(&affinity_aa(&child_sim, refdata)?),
                    None => 0.0,
                };
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
                if retain_edges {
                    edges.push(edge);
                }
                next_live.push(next_id);
                next_id = next_id
                    .checked_add(1)
                    .ok_or_else(|| super::replay_receipt::error("lineage node ID overflow"))?;
            }
        }
        live = next_live;
        completed_generations = gen;
        if live.len() > peak_live {
            peak_live = live.len();
        }
    }
    if let Some(driver) = replay.as_deref_mut() {
        driver.finish_growth(&live, completed_generations, peak_live)?;
        if driver.is_replay() && rng.words_consumed() != 0 {
            return Err(super::replay_receipt::error(
                "semantic growth replay consumed RNG",
            ));
        }
    }

    // `live` now holds the cells alive when growth stopped: either the last
    // completed generation's offspring (ran to max_generations) or empty (the
    // population went extinct — every live cell drew 0 offspring). This is the
    // living population sampling must draw from.
    Ok((LineageTree { nodes }, peak_live, sims, live, edges))
}

/// Grow a full clonal lineage with per-division mutation via `mutator`.
/// Deterministic for `params.seed`.
pub fn grow_lineage(
    founder: &Simulation,
    params: &BranchingParams,
    mutator: &dyn Pass,
) -> LineageTree {
    grow_core(founder, params, Some(mutator), None, None, None)
        .expect("unconstrained lineage")
        .0
}

/// Grow the lineage TOPOLOGY only (children are exact clones of their parent
/// `Simulation`; no mutation — a later task layers mutation on top).
pub fn grow_topology(founder: &Simulation, params: &BranchingParams) -> LineageTree {
    grow_core(founder, params, None, None, None, None)
        .expect("unconstrained lineage")
        .0
}

/// Like `grow_topology`, but also returns the final living node ids (cells alive
/// when growth stopped). Empty when the family went extinct.
#[cfg(test)]
pub fn grow_topology_with_live(
    founder: &Simulation,
    params: &BranchingParams,
) -> (LineageTree, Vec<u32>) {
    let (tree, _peak, _sims, live, _edges) =
        grow_core(founder, params, None, None, None, None).expect("unconstrained lineage");
    (tree, live)
}

/// Grow a clonal lineage with per-division mutation AND affinity selection.
/// Deterministic for `params.seed`.
pub fn grow_lineage_with_affinity(
    founder: &Simulation,
    params: &BranchingParams,
    mutator: &dyn Pass,
    model: &AffinityModel,
    refdata: &RefDataConfig,
) -> LineageTree {
    grow_core(founder, params, Some(mutator), Some((model, refdata)), None, None)
        .expect("unconstrained lineage")
        .0
}

/// Grow + mutate a lineage and ALSO return the per-node `Simulation` arena
/// (index == node id) and the final living node ids (cells alive when growth
/// stopped; empty when the family went extinct), for building per-node AIRR
/// `Outcome`s and sampling the living population.
pub fn grow_lineage_retaining_sims(
    founder: &Simulation,
    params: &BranchingParams,
    mutator: &dyn Pass,
    affinity: Option<(&AffinityModel, &RefDataConfig)>,
) -> (LineageTree, Vec<Simulation>, Vec<u32>) {
    let (tree, _peak, sims, live, _edges) =
        grow_core(founder, params, Some(mutator), affinity, None, None)
            .expect("unconstrained lineage");
    (tree, sims, live)
}

pub fn grow_lineage_retaining_sims_with_context(
    founder: &Simulation,
    params: &BranchingParams,
    mutator: &dyn Pass,
    affinity: Option<(&AffinityModel, &RefDataConfig)>,
    context: &LineageMutationContext<'_>,
) -> Result<(LineageTree, Vec<Simulation>, Vec<u32>), PassError> {
    let (tree, _peak, sims, live, _edges) = grow_core(
        founder,
        params,
        Some(mutator),
        affinity,
        Some(context),
        None,
    )?;
    Ok((tree, sims, live))
}

/// Duplex-specific growth retains each incoming mutation transaction for
/// exact branch replay without copying the whole ancestral path per cell.
pub fn grow_lineage_retaining_edges_with_context(
    founder: &Simulation,
    params: &BranchingParams,
    mutator: &dyn Pass,
    affinity: Option<(&AffinityModel, &RefDataConfig)>,
    context: &LineageMutationContext<'_>,
) -> Result<(LineageTree, Vec<Simulation>, Vec<u32>, Vec<Option<Outcome>>), PassError> {
    let (tree, _peak, sims, live, edges) = grow_core(
        founder,
        params,
        Some(mutator),
        affinity,
        Some(context),
        None,
    )?;
    Ok((tree, sims, live, edges))
}

pub(crate) fn grow_lineage_with_receipt(
    founder: &Simulation,
    params: &BranchingParams,
    mutator: &dyn Pass,
    affinity: Option<(&AffinityModel, &RefDataConfig)>,
    context: &LineageMutationContext<'_>,
    driver: &mut super::replay_receipt::FamilyReplayDriver<'_>,
) -> Result<(LineageTree, Vec<Simulation>, Vec<u32>, Vec<Option<Outcome>>), PassError> {
    let (tree, _, sims, live, edges) = grow_core(
        founder,
        params,
        Some(mutator),
        affinity,
        Some(context),
        Some(driver),
    )?;
    Ok((tree, sims, live, edges))
}

/// Replay one stored division. All choices come from the recorded trace;
/// even a hidden fallback draw is rejected by the RNG consumption fence.
pub(crate) fn execute_mutation_choices(
    parent: &Simulation,
    choices: &Trace,
    mutator: &crate::passes::S5FMutationPass,
    context: &LineageMutationContext<'_>,
) -> Result<Outcome, PassError> {
    let fail = |reason: String| PassError::invalid_plan_state("lineage.replay", reason);
    let mut cursor = crate::replay::TraceCursor::new(choices.choices());
    let mut trace = Trace::new();
    let mut rng = Rng::placeholder();
    let mut events = Vec::new();
    let mut ctx = PassContext {
        replay_cursor: Some(&mut cursor),
        trace: &mut trace,
        rng: &mut rng,
        pass_index: 0,
        refdata: Some(context.refdata),
        contracts: context.contracts,
        feasibility: None,
        reference_index: None,
        event_log_sink: Some(&mut events),
    };
    let child = mutator.execute_for_lineage(parent, &mut ctx, context.strict)?;
    drop(ctx);
    cursor.assert_drained().map_err(|e| fail(e.to_string()))?;
    if rng.words_consumed() != 0 {
        return Err(fail("division replay consumed RNG".into()));
    }
    if context.strict {
        if let Some(contracts) = context.contracts {
            contracts
                .verify(&child, Some(context.refdata))
                .map_err(|v| PassError::contract_violation("lineage.replay", v))?;
        }
    }
    require_productive(&child, Some(context))?;
    let event = EventRecord::pass_committed(
        0,
        mutator.name(),
        mutator.effects(),
        TraceSpan::new(0, choices.len()),
        StateSummary::from_simulation(parent),
        StateSummary::from_simulation(&child),
        events,
    );
    Ok(Outcome {
        revisions: vec![parent.clone(), child],
        pass_names: vec![mutator.name().into()],
        trace: choices.clone(),
        events: vec![event],
    })
}

pub fn replay_mutation_edge(
    parent: &Simulation,
    recorded: &Outcome,
    mutator: &crate::passes::S5FMutationPass,
    context: &LineageMutationContext<'_>,
) -> Result<Outcome, PassError> {
    let fail = |reason: String| PassError::invalid_plan_state("lineage.replay", reason);
    let boundary = Outcome {
        revisions: vec![parent.clone()],
        pass_names: vec![],
        trace: Trace::new(),
        events: vec![],
    };
    boundary
        .with_continuation(recorded)
        .map_err(|s| fail(s.into()))?;
    if recorded.pass_names != [mutator.name()]
        || recorded.revisions.len() != 2
        || recorded.events.len() != 1
    {
        return Err(fail(
            "recorded division has an invalid pass/revision/event layout".into(),
        ));
    }
    let replayed = execute_mutation_choices(parent, &recorded.trace, mutator, context)?;
    let child = replayed.final_simulation();
    let event = &replayed.events[0];
    let expected = recorded.final_simulation();
    if child.pool.as_slice() != expected.pool.as_slice()
        || child.sequence.regions != expected.sequence.regions
        || child.assignments != expected.assignments
        || child.duplex_provenance != expected.duplex_provenance
        || child.coding_end_provenance != expected.coding_end_provenance
        || child.live_call_support != expected.live_call_support
        || child.mutation_count != expected.mutation_count
        || child.mutation_burden.to_bits() != expected.mutation_burden.to_bits()
        || event != &recorded.events[0]
    {
        return Err(fail(
            "division replay differs from stored state or event evidence".into(),
        ));
    }
    Ok(replayed)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::dist::{EmpiricalLengthDist, UniformBase};
    use crate::ir::{NucHandle, Nucleotide, Region, Segment, Simulation};
    use crate::passes::UniformMutationPass;

    /// Minimal founder: a 4-base V region "AAAA".
    fn founder() -> Simulation {
        let mut sim = Simulation::new();
        for (i, b) in b"AAAA".iter().enumerate() {
            let (next, _) =
                sim.with_nucleotide_pushed(Nucleotide::germline(*b, i as u16, Segment::V));
            sim = next;
        }
        sim.with_region_added(Region::new(
            Segment::V,
            NucHandle::new(0),
            NucHandle::new(4),
        ))
    }

    /// A founder with a resolved V assignment, and the reference it was built
    /// from. Affinity scoring needs both: the coding frame comes from the V
    /// allele's native reference phase, so a translation cannot be derived from
    /// the pool alone. The synthetic `founder()` has no assignment and no
    /// reference, which is why these tests used to score the per-region rail.
    fn affinity_founder() -> (crate::refdata::RefDataConfig, Simulation) {
        let (reference, model) = crate::duplex_repair::model::tests::fixture();
        let mut sim = model.prepare_row(0, &reference).unwrap().seed;
        // That fixture is the reachable reference-plus-assignment pair, but its
        // founder carries duplex provenance, which routes growth down the
        // duplex path and demands a mutation context. These tests exercise the
        // plain path, so the provenance is dropped.
        sim.duplex_provenance = None;
        (reference, sim)
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

    /// How many seeds a premise search may try before the search is
    /// itself the thing that is broken.
    const SEED_SCAN: u64 = 200;

    /// `base` with the first seed whose family satisfies `premise`.
    ///
    /// A branching process with `lambda_base` near 1 goes extinct at
    /// the founder often -- about a fifth of the time at 1.5, a third
    /// at 1.2 -- so a test about a family that grew, or about a child
    /// that mutated, has to go and find one. A seed that happened to
    /// grow when the test was written stops growing the moment the
    /// growth stream changes, and the failure then reads as a defect
    /// in growth rather than as a test naming a seed instead of its
    /// premise.
    fn seed_where(
        base: &BranchingParams,
        premise: impl Fn(&BranchingParams) -> bool,
    ) -> BranchingParams {
        (0..SEED_SCAN)
            .map(|seed| BranchingParams {
                seed,
                ..base.clone()
            })
            .find(|p| premise(p))
            .unwrap_or_else(|| {
                panic!("no seed below {SEED_SCAN} produced a family satisfying the premise")
            })
    }

    #[test]
    fn lineage_divisions_inherit_productive_support_and_final_guarantee() {
        use crate::contract::test_support::{
            make_assembled_sim_from_refdata, make_vj_with_anchor_codons,
        };
        let cfg = make_vj_with_anchor_codons(b"TGT", b"TGG");
        let founder = make_assembled_sim_from_refdata(&cfg);
        let mutator = UniformMutationPass::new(
            Box::new(EmpiricalLengthDist::from_pairs(vec![(15, 1.0)])),
            Box::new(UniformBase),
        );
        let contracts = crate::contract::productive();
        let params = BranchingParams {
            lambda_base: 2.0,
            max_generations: 3,
            n_max: 40,
            n_sample: 20,
            seed: 7,
            lambda_mut: 0.0,
        };
        for strict in [false, true] {
            let context = LineageMutationContext {
                refdata: &cfg,
                contracts: Some(&contracts),
                strict,
            };
            let (tree, sims, _) = grow_lineage_retaining_sims_with_context(
                &founder, &params, &mutator, None, &context,
            )
            .unwrap();
            assert!(tree.nodes.len() > 1);
            assert!(sims
                .iter()
                .all(|s| crate::productivity::assess(s, &cfg).productive == Some(true)));
            assert!(sims.iter().skip(1).any(|s| s.mutation_count > 0));
        }
        // Supplying reference context alone must not change the unconstrained law.
        let (_, old, _) = grow_lineage_retaining_sims(&founder, &params, &mutator, None);
        let context = LineageMutationContext {
            refdata: &cfg,
            contracts: None,
            strict: false,
        };
        let (_, new, _) =
            grow_lineage_retaining_sims_with_context(&founder, &params, &mutator, None, &context)
                .unwrap();
        assert_eq!(old.len(), new.len());
        for (a, b) in old.iter().zip(&new) {
            assert_eq!(a.pool.as_slice(), b.pool.as_slice());
        }
    }

    #[test]
    fn grows_root_and_children_with_monotonic_generations() {
        let params = seed_where(&neutral_params(), |p| {
            grow_topology(&founder(), p).len() > 1
        });
        let tree = grow_topology(&founder(), &params);
        assert_eq!(
            tree.nodes.iter().filter(|n| n.parent_id.is_none()).count(),
            1
        );
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
        // On a seed whose family went extinct at the founder this
        // compares one node with one node and claims determinism.
        let params = seed_where(&neutral_params(), |p| {
            grow_topology(&founder(), p).len() > 1
        });
        let a = grow_topology(&founder(), &params);
        let b = grow_topology(&founder(), &params);
        assert!(a.len() > 1, "nothing grew, so nothing was compared");
        assert_eq!(a.len(), b.len());
        for (x, y) in a.nodes.iter().zip(b.nodes.iter()) {
            assert_eq!(x.id, y.id);
            assert_eq!(x.parent_id, y.parent_id);
            assert_eq!(x.generation, y.generation);
        }
    }

    #[test]
    fn regulation_settles_the_population_at_capacity() {
        let params = BranchingParams {
            lambda_base: 4.0, // would explode unbounded
            lambda_mut: 0.0,
            max_generations: 40,
            n_max: 200,
            n_sample: 10,
            seed: 7,
        };
        let (_tree, peak_live, _sims, live, _edges) =
            grow_core(&founder(), &params, None, None, None, None).unwrap();
        // This is the assertion that was missing. For 40 generations the test
        // asserted a bound on a lineage that had been extinct since generation
        // 11 and called the result a plateau: the hard cap landed the
        // population on n_max exactly, where the old logistic law returned an
        // effective rate of exactly zero.
        assert!(
            !live.is_empty(),
            "the population was regulated and then died: peak {peak_live}, \
             nothing alive after {} generations",
            params.max_generations
        );
        // Beverton-Holt's fixed point is n_max itself. The logistic law it
        // replaces settled at n_max * (1 - 1/lambda) = 150 here, so this bound
        // also distinguishes the two laws rather than merely bounding growth.
        let n = params.n_max as f64;
        assert!(
            (live.len() as f64) > 0.8 * n && (live.len() as f64) < 1.2 * n,
            "final live {} is not at capacity {}",
            live.len(),
            params.n_max
        );
        // A generation may overshoot capacity, by design -- truncating it would
        // have to pick which cells die. Regulation returns a rate below 1 above
        // n_max, so the overshoot stays small and self-correcting.
        assert!(
            peak_live < 2 * params.n_max as usize,
            "peak live {peak_live} ran away from capacity {}",
            params.n_max
        );
    }

    /// The plateau is `n_max` whatever `lambda_base` is. Under the logistic law
    /// this replaces it was `n_max * (1 - 1/lambda)` -- 100 at lambda 2, 150 at
    /// lambda 4 -- so `n_max` named a capacity the population never reached.
    #[test]
    fn the_plateau_does_not_depend_on_lambda_base() {
        let settled = |lambda_base: f64| {
            let params = BranchingParams {
                lambda_base,
                lambda_mut: 0.0,
                max_generations: 40,
                n_max: 200,
                n_sample: 1,
                seed: 7,
            };
            let (_t, _peak, _s, live, _e) =
                grow_core(&founder(), &params, None, None, None, None).unwrap();
            live.len()
        };
        let slow = settled(2.0);
        let fast = settled(4.0);
        // Guard the premise: two dead lineages would agree at zero.
        assert!(slow > 0 && fast > 0, "a lineage died: slow {slow}, fast {fast}");
        for (name, n) in [("lambda 2", slow), ("lambda 4", fast)] {
            assert!(
                n > 160 && n < 240,
                "{name} settled at {n}, not at the capacity 200"
            );
        }
    }

    /// Reaching the carrying capacity must not be fatal. The parameter's name
    /// and the damping comment both promised a population that levels off at
    /// `n_max`; the logistic law plus a hard cap delivered a population that
    /// landed on `n_max` exactly, drew an effective rate of exactly zero, and
    /// died in one generation with certainty -- 67 cap hits, 67 extinctions.
    #[test]
    fn no_lineage_that_reaches_capacity_dies_from_reaching_it() {
        let mut reached = 0;
        let mut died = 0;
        for seed in 1..40u64 {
            let params = BranchingParams {
                lambda_base: 4.0,
                lambda_mut: 0.0,
                max_generations: 40,
                n_max: 50,
                n_sample: 1,
                seed,
            };
            let (_t, peak, _s, live, _e) =
                grow_core(&founder(), &params, None, None, None, None).unwrap();
            if peak >= params.n_max as usize {
                reached += 1;
                if live.is_empty() {
                    died += 1;
                }
            }
        }
        // Guard the premise: seeds that never reach capacity say nothing about
        // what reaching it costs.
        assert!(reached > 30, "only {reached} of 39 seeds reached capacity");
        assert_eq!(died, 0, "{died} of {reached} lineages died after reaching capacity");
    }

    #[test]
    fn final_live_is_empty_when_founder_never_divides() {
        // lambda_base = 0 => the founder draws 0 offspring => after generation 1
        // the live set is empty and the loop breaks. The final living population
        // is empty (the clone went extinct).
        let params = BranchingParams {
            lambda_base: 0.0,
            ..neutral_params()
        };
        let (_tree, live) = grow_topology_with_live(&founder(), &params);
        assert!(
            live.is_empty(),
            "expected extinct (empty) final-live set, got {:?}",
            live
        );
    }

    #[test]
    fn final_live_is_the_max_generation_set() {
        // In a grown family, every node in the final-live set is at the maximum
        // generation reached, and the set is non-empty. A single founder has a
        // real chance of immediate extinction at lambda_base ~1.5, so scan seeds
        // for one that grows (the invariant under test is about non-extinct runs).
        let mut grew = false;
        for seed in 0..50u64 {
            let params = BranchingParams {
                seed,
                ..neutral_params()
            };
            let (tree, live) = grow_topology_with_live(&founder(), &params);
            if live.is_empty() {
                continue; // extinct for this seed; try another
            }
            grew = true;
            let max_gen = tree.nodes.iter().map(|n| n.generation).max().unwrap();
            for &id in &live {
                assert_eq!(
                    tree.get(id).unwrap().generation,
                    max_gen,
                    "final-live node {id} not at max generation {max_gen} (seed {seed})"
                );
            }
        }
        assert!(grew, "no seed in 0..50 produced a surviving family");
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
        let base = BranchingParams {
            lambda_base: 1.2,
            lambda_mut: 0.0,
            max_generations: 4,
            n_max: 500,
            n_sample: 10,
            seed: 11,
        };
        let mutator = two_mut_mutator();
        let params = seed_where(&base, |p| {
            let tree = grow_lineage(&founder(), p, &mutator);
            tree.nodes.iter().any(|n| {
                n.parent_id
                    .and_then(|pid| tree.get(pid))
                    .is_some_and(|parent| n.genotype != parent.genotype)
            })
        });
        let tree = grow_lineage(&founder(), &params, &mutator);

        assert_eq!(tree.root().mutations_from_parent, 0);

        // The implication runs one way only. A changed genotype has
        // to be accounted for by a recorded mutation. The converse is
        // false: on a four-base "AAAA" founder the child's two
        // substitutions can land on the same position and compose
        // back to the base they started from -- A to X, then X to A,
        // two real draws and no net change. (Neither draw is an
        // identity: a substitution never draws the base the site
        // already holds, measured here as 1464 one-substitution
        // children with not one unchanged sequence.) Asserting the
        // converse is what pinned this test to one seed.
        let mut saw_changed_child = false;
        for n in &tree.nodes {
            let Some(parent_id) = n.parent_id else {
                continue;
            };
            // <= 2 (not == 2): the two uniform substitutions can draw
            // the same position, yielding one net change.
            assert!(n.mutations_from_parent <= 2);
            let parent = tree.get(parent_id).unwrap();
            if n.genotype != parent.genotype {
                saw_changed_child = true;
                assert!(
                    n.mutations_from_parent > 0,
                    "a changed genotype with no recorded mutation"
                );
            }
        }
        assert!(saw_changed_child, "no child's genotype changed");
    }

    #[test]
    fn mutated_growth_is_deterministic() {
        let params = BranchingParams {
            lambda_base: 1.2,
            lambda_mut: 0.0,
            max_generations: 4,
            n_max: 500,
            n_sample: 10,
            seed: 11,
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
        use crate::lineage::affinity::{coding_sim_to_aa, AffinityModel};
        let (reference, f) = affinity_founder();
        let founder_aa = coding_sim_to_aa(&f, &reference).expect("founder has a coding span");
        let w = vec![1.0; founder_aa.len().max(1)];
        let model = AffinityModel::new(b"W".to_vec(), w, 1.0, 1.0, &founder_aa);
        let params = BranchingParams {
            lambda_base: 1.2,
            lambda_mut: 0.0,
            max_generations: 4,
            n_max: 500,
            n_sample: 10,
            seed: 11,
        };
        let tree = grow_lineage_with_affinity(&f, &params, &two_mut_mutator(), &model, &reference);
        assert!((tree.root().affinity - model.affinity_value(&founder_aa)).abs() < 1e-9);
        for n in &tree.nodes {
            assert!(
                n.affinity > 0.0 && n.affinity <= 1.0 + 1e-9,
                "affinity out of range: {}",
                n.affinity
            );
        }
    }

    #[test]
    fn affinity_strength_zero_matches_neutral_topology() {
        use crate::lineage::affinity::{coding_sim_to_aa, AffinityModel};
        let (reference, f) = affinity_founder();
        let founder_aa = coding_sim_to_aa(&f, &reference).expect("founder has a coding span");
        let w = vec![1.0; founder_aa.len().max(1)];
        // selection_strength = 0 => fitness identically 1 => identical topology/genotypes to neutral
        let model = AffinityModel::new(b"W".to_vec(), w, 1.0, 0.0, &founder_aa);
        let params = BranchingParams {
            lambda_base: 1.2,
            lambda_mut: 0.0,
            max_generations: 4,
            n_max: 500,
            n_sample: 10,
            seed: 11,
        };
        let with_aff = grow_lineage_with_affinity(&f, &params, &two_mut_mutator(), &model, &reference);
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
        use crate::lineage::affinity::{coding_sim_to_aa, AffinityModel};
        let (reference, f) = affinity_founder();
        let founder_aa = coding_sim_to_aa(&f, &reference).expect("founder has a coding span");
        let mk = || {
            AffinityModel::new(
                b"W".to_vec(),
                vec![1.0; founder_aa.len().max(1)],
                1.0,
                2.0,
                &founder_aa,
            )
        };
        let params = BranchingParams {
            lambda_base: 1.5,
            lambda_mut: 0.0,
            max_generations: 5,
            n_max: 500,
            n_sample: 10,
            seed: 21,
        };
        let a = grow_lineage_with_affinity(&f, &params, &two_mut_mutator(), &mk(), &reference);
        let b = grow_lineage_with_affinity(&f, &params, &two_mut_mutator(), &mk(), &reference);
        assert_eq!(a.len(), b.len());
        for (x, y) in a.nodes.iter().zip(b.nodes.iter()) {
            assert_eq!(x.genotype, y.genotype);
            assert!((x.affinity - y.affinity).abs() < 1e-12);
        }
    }

}
