//! PyO3 bindings for the clonal lineage engine: `LineageNode`, `LineageTree`,
//! and the `simulate_lineage` entry point.

use pyo3::prelude::*;
use pyo3::types::PyDict;

use crate::airr_record::build_airr_record;
use crate::event::{EventRecord, StateSummary, TraceSpan};
use crate::ir::{NucHandle, SimulationEvent};
use crate::ir::Simulation;
use crate::lineage::export::{to_fasta, to_newick, to_node_table_tsv};
use crate::lineage::tree::{LineageNode, LineageTree};
use crate::lineage::{simulate_family, simulate_family_sims, simulate_family_with_affinity, sim_to_aa, AffinityModel, BranchingParams};
use crate::lineage::affinity::make_mature_target;
use crate::live_call::{refresh_live_calls, ReferenceMatchIndex};
use crate::pass::{Outcome, PassCompileEffect};
use crate::rng::Rng;
use crate::passes::S5FMutationPass;
use crate::s5f::S5FKernel;

use super::outcome::{PyOutcome, airr_record_to_pydict};
use super::refdata::PyRefDataConfig;
use super::simulation::PySimulation;

/// One node of a clonal lineage tree (read-only view).
#[pyclass(name = "LineageNode", module = "GenAIRR._engine", frozen)]
pub struct PyLineageNode {
    pub(crate) inner: LineageNode,
}

#[pymethods]
impl PyLineageNode {
    #[getter]
    fn id(&self) -> u32 {
        self.inner.id
    }
    #[getter]
    fn parent_id(&self) -> Option<u32> {
        self.inner.parent_id
    }
    #[getter]
    fn generation(&self) -> u32 {
        self.inner.generation
    }
    #[getter]
    fn mutations_from_parent(&self) -> u32 {
        self.inner.mutations_from_parent
    }
    #[getter]
    fn abundance(&self) -> u32 {
        self.inner.abundance
    }
    #[getter]
    fn observed(&self) -> bool {
        self.inner.observed
    }
    #[getter]
    fn affinity(&self) -> f64 {
        self.inner.affinity
    }
    /// Nucleotide sequence (pool bases) as a string.
    #[getter]
    fn sequence(&self) -> String {
        String::from_utf8_lossy(&self.inner.genotype).into_owned()
    }

    fn __repr__(&self) -> String {
        format!(
            "LineageNode(id={}, parent_id={:?}, generation={}, abundance={}, observed={})",
            self.inner.id,
            self.inner.parent_id,
            self.inner.generation,
            self.inner.abundance,
            self.inner.observed
        )
    }
}

/// A clonal lineage tree (read-only view) with ground-truth export.
#[pyclass(name = "LineageTree", module = "GenAIRR._engine", frozen)]
pub struct PyLineageTree {
    pub(crate) inner: LineageTree,
}

impl PyLineageTree {
    pub(crate) fn new(inner: LineageTree) -> Self {
        Self { inner }
    }
}

#[pymethods]
impl PyLineageTree {
    fn __len__(&self) -> usize {
        self.inner.len()
    }

    /// All nodes (founder + descendants), in arena order (ascending id).
    fn nodes(&self) -> Vec<PyLineageNode> {
        self.inner
            .nodes
            .iter()
            .cloned()
            .map(|inner| PyLineageNode { inner })
            .collect()
    }

    /// Validate structural invariants; raises ValueError if malformed.
    fn validate(&self) -> PyResult<()> {
        self.inner
            .validate()
            .map_err(pyo3::exceptions::PyValueError::new_err)
    }

    /// Standard rooted Newick string (branch length = per-edge mutation count).
    fn to_newick(&self) -> String {
        to_newick(&self.inner)
    }

    /// FASTA of every node (ancestral + observed) sequence.
    fn to_fasta(&self) -> String {
        to_fasta(&self.inner)
    }

    /// Tab-separated node table.
    fn to_node_table_tsv(&self) -> String {
        to_node_table_tsv(&self.inner)
    }
}

/// Grow + sample a clonal lineage family from `founder` using an S5F mutator
/// built from the supplied kernel tables. Returns the ground-truth tree.
///
/// `mutability` must have 1024 entries and `substitution` 4096 (the S5F 5-mer
/// kernel); all values must be finite and non-negative. `rate` is the per-base
/// SHM rate in [0, 1]. Determinism is keyed on `seed`. `lambda_mut` is reserved
/// for a future plan and currently has no effect (mutations are driven by the
/// S5F `rate`).
#[pyfunction]
#[pyo3(signature = (
    founder, mutability, substitution, rate,
    lambda_base, lambda_mut, max_generations, n_max, n_sample, seed,
    selection_strength=0.0, beta=1.0, target_aa=None, mature_substitutions=5
))]
#[allow(clippy::too_many_arguments)]
pub(crate) fn simulate_lineage(
    founder: &PySimulation,
    mutability: Vec<f64>,
    substitution: Vec<f64>,
    rate: f64,
    lambda_base: f64,
    lambda_mut: f64,
    max_generations: u32,
    n_max: u32,
    n_sample: u32,
    seed: u64,
    selection_strength: f64,
    beta: f64,
    target_aa: Option<String>,
    mature_substitutions: u32,
) -> PyResult<PyLineageTree> {
    use pyo3::exceptions::PyValueError;

    if mutability.len() != 1024 {
        return Err(PyValueError::new_err(format!(
            "simulate_lineage: mutability must have 1024 entries, got {}",
            mutability.len()
        )));
    }
    if substitution.len() != 4096 {
        return Err(PyValueError::new_err(format!(
            "simulate_lineage: substitution must have 4096 entries, got {}",
            substitution.len()
        )));
    }
    if !(rate.is_finite() && (0.0..=1.0).contains(&rate)) {
        return Err(PyValueError::new_err(format!(
            "simulate_lineage: rate must be in [0.0, 1.0], got {rate}"
        )));
    }
    if mutability.iter().any(|&m| !m.is_finite() || m < 0.0) {
        return Err(PyValueError::new_err(
            "simulate_lineage: mutability values must be finite and non-negative",
        ));
    }
    if substitution.iter().any(|&s| !s.is_finite() || s < 0.0) {
        return Err(PyValueError::new_err(
            "simulate_lineage: substitution values must be finite and non-negative",
        ));
    }
    if n_max == 0 {
        return Err(PyValueError::new_err("simulate_lineage: n_max must be > 0"));
    }
    if n_sample == 0 {
        return Err(PyValueError::new_err("simulate_lineage: n_sample must be > 0"));
    }
    if !(lambda_base.is_finite() && lambda_base >= 0.0) {
        return Err(PyValueError::new_err(format!(
            "simulate_lineage: lambda_base must be finite and >= 0, got {lambda_base}"
        )));
    }
    if !(lambda_mut.is_finite() && lambda_mut >= 0.0) {
        return Err(PyValueError::new_err(format!(
            "simulate_lineage: lambda_mut must be finite and >= 0, got {lambda_mut}"
        )));
    }
    // Cap generations: `to_newick` recurses to a depth equal to the generation
    // count, so an unbounded value could overflow the stack across the FFI
    // boundary. 1000 is far beyond any biological germinal-center reaction.
    if max_generations > 1000 {
        return Err(PyValueError::new_err(format!(
            "simulate_lineage: max_generations must be <= 1000, got {max_generations}"
        )));
    }

    // Lengths and value ranges are now guaranteed, so S5FKernel::new cannot panic.
    let kernel = S5FKernel::new(mutability, substitution);
    let mutator = S5FMutationPass::new_rate(kernel, rate);
    let params = BranchingParams {
        lambda_base,
        lambda_mut,
        max_generations,
        n_max,
        n_sample,
        seed,
    };

    // Validate affinity params.
    if !(beta.is_finite() && beta >= 0.0) {
        return Err(PyValueError::new_err(format!(
            "simulate_lineage: beta must be finite and >= 0, got {beta}"
        )));
    }
    if !(selection_strength.is_finite() && selection_strength >= 0.0) {
        return Err(PyValueError::new_err(format!(
            "simulate_lineage: selection_strength must be finite and >= 0, got {selection_strength}"
        )));
    }

    // Neutral fast path: no selection and no explicit target → byte-identical to
    // the pre-affinity behavior.
    if selection_strength == 0.0 && target_aa.is_none() {
        let tree = simulate_family(&founder.inner, &params, &mutator);
        return Ok(PyLineageTree::new(tree));
    }

    // Build the affinity model. Uniform per-position weights for v1 (the model
    // accepts arbitrary weights; CDR3 region-weighting is a follow-up).
    let founder_aa = sim_to_aa(&founder.inner);
    let target = match target_aa {
        Some(s) => s.into_bytes(),
        None => {
            // Deterministic auto "mature" target derived from the seed.
            let mut trng = Rng::new(seed ^ 0x7461_7267_6574_0001); // "target\0\1"
            make_mature_target(&founder_aa, mature_substitutions, &mut trng)
        }
    };
    let weights = vec![1.0; founder_aa.len()];
    let model = AffinityModel::new(target, weights, beta, selection_strength, &founder_aa);
    let tree = simulate_family_with_affinity(&founder.inner, &params, &mutator, &model);
    Ok(PyLineageTree::new(tree))
}

/// A grown clonal family: the lineage tree plus per-node AIRR-projectable
/// `Outcome`s (only observed nodes carry one).
#[pyclass(name = "FamilyOutcome", module = "GenAIRR._engine")]
pub struct PyFamilyOutcome {
    tree: LineageTree,
    node_outcomes: Vec<Option<Outcome>>, // index == node id
}

#[pymethods]
impl PyFamilyOutcome {
    /// The ground-truth lineage tree.
    fn tree(&self) -> PyLineageTree {
        PyLineageTree::new(self.tree.clone())
    }

    /// Per-node Outcomes aligned with `tree().nodes()`; None for unsampled nodes.
    fn node_outcomes(&self) -> Vec<Option<PyOutcome>> {
        self.node_outcomes
            .iter()
            .cloned()
            .map(|o: Option<Outcome>| o.map(PyOutcome::new))
            .collect()
    }

    /// Observed nodes' Outcomes only (abundance > 0), in node-id order.
    fn observed_outcomes(&self) -> Vec<PyOutcome> {
        self.node_outcomes
            .iter()
            .cloned()
            .flatten()
            .map(PyOutcome::new)
            .collect()
    }

    /// Build AIRR Rearrangement record dicts for every observed node,
    /// in node-id order (aligned with `observed_outcomes()`).
    ///
    /// Each node `Outcome` now carries a synthesized `MUTATE_S5F`
    /// `EventRecord` encoding every pool position where
    /// `base != germline` as a `BaseChanged` event, and
    /// `mutation_count` is set to the number of such events.
    /// `build_airr_record` therefore produces correct per-segment
    /// and V-subregion mutation counts natively, without any manual
    /// field overwrite.
    fn airr_records<'py>(
        &self,
        py: Python<'py>,
        refdata: &PyRefDataConfig,
    ) -> PyResult<Vec<Bound<'py, PyDict>>> {
        let mut results = Vec::new();
        for (node_id, outcome) in self.node_outcomes.iter().enumerate() {
            let Some(outcome) = outcome else { continue };
            let sequence_id = format!("node{}", node_id);
            let rec = build_airr_record(outcome, refdata.inner(), &sequence_id);
            let dict = airr_record_to_pydict(py, &rec)?;
            results.push(dict);
        }
        Ok(results)
    }
}

// ──────────────────────────────────────────────────────────────────
// SHM event synthesis helper
// ──────────────────────────────────────────────────────────────────

/// Synthesize a single `EventRecord` that encodes all net SHM
/// mutations visible in `sim.pool` (positions where `base !=
/// germline`) as `SimulationEvent::BaseChanged` events.
///
/// The record is tagged with `address::MUTATE_S5F` as its pass
/// name so the AIRR builder's event-loop filter picks it up and
/// the per-segment / per-subregion mutation counts it produces are
/// identical to the pool-scan. The returned count equals the number
/// of synthesized events; callers should set `sim.mutation_count`
/// to this value before wrapping `sim` in an `Outcome`.
fn synthesize_shm_event_record(sim: &Simulation) -> (EventRecord, u32) {
    let mut events: Vec<SimulationEvent> = Vec::new();
    for (i, nuc) in sim.pool.as_slice().iter().enumerate() {
        if nuc.base == nuc.germline {
            continue;
        }
        events.push(SimulationEvent::BaseChanged {
            handle: NucHandle::new(i as u32),
            old_base: nuc.germline,
            new_base: nuc.base,
            segment: nuc.segment,
            germline_pos: nuc.germline_pos.get(),
        });
    }
    let net = events.len() as u32;
    let summary = StateSummary::from_simulation(sim);
    let record = EventRecord::pass_committed(
        0,
        crate::address::MUTATE_S5F,
        vec![PassCompileEffect::EditBases],
        TraceSpan::new(0, 0),
        summary,
        summary,
        events,
    );
    (record, net)
}


/// Grow + sample a clonal lineage family from `founder` (a full `Outcome`) using
/// an S5F mutator, returning a `FamilyOutcome` with per-node AIRR-projectable
/// `Outcome`s for every observed (sampled) node.
#[pyfunction]
#[pyo3(signature = (
    founder, refdata, mutability, substitution, rate,
    lambda_base, lambda_mut, max_generations, n_max, n_sample, seed,
    selection_strength=0.0, beta=1.0, target_aa=None, mature_substitutions=5
))]
#[allow(clippy::too_many_arguments)]
pub(crate) fn simulate_family_outcomes(
    founder: &PyOutcome,
    refdata: &PyRefDataConfig,
    mutability: Vec<f64>,
    substitution: Vec<f64>,
    rate: f64,
    lambda_base: f64,
    lambda_mut: f64,
    max_generations: u32,
    n_max: u32,
    n_sample: u32,
    seed: u64,
    selection_strength: f64,
    beta: f64,
    target_aa: Option<String>,
    mature_substitutions: u32,
) -> PyResult<PyFamilyOutcome> {
    use pyo3::exceptions::PyValueError;

    if mutability.len() != 1024 {
        return Err(PyValueError::new_err(format!(
            "simulate_family_outcomes: mutability must have 1024 entries, got {}",
            mutability.len()
        )));
    }
    if substitution.len() != 4096 {
        return Err(PyValueError::new_err(format!(
            "simulate_family_outcomes: substitution must have 4096 entries, got {}",
            substitution.len()
        )));
    }
    if !(rate.is_finite() && (0.0..=1.0).contains(&rate)) {
        return Err(PyValueError::new_err(format!(
            "simulate_family_outcomes: rate must be in [0.0, 1.0], got {rate}"
        )));
    }
    if mutability.iter().any(|&m| !m.is_finite() || m < 0.0) {
        return Err(PyValueError::new_err(
            "simulate_family_outcomes: mutability values must be finite and non-negative",
        ));
    }
    if substitution.iter().any(|&s| !s.is_finite() || s < 0.0) {
        return Err(PyValueError::new_err(
            "simulate_family_outcomes: substitution values must be finite and non-negative",
        ));
    }
    if n_max == 0 {
        return Err(PyValueError::new_err("simulate_family_outcomes: n_max must be > 0"));
    }
    if n_sample == 0 {
        return Err(PyValueError::new_err("simulate_family_outcomes: n_sample must be > 0"));
    }
    if !(lambda_base.is_finite() && lambda_base >= 0.0) {
        return Err(PyValueError::new_err(format!(
            "simulate_family_outcomes: lambda_base must be finite and >= 0, got {lambda_base}"
        )));
    }
    if !(lambda_mut.is_finite() && lambda_mut >= 0.0) {
        return Err(PyValueError::new_err(format!(
            "simulate_family_outcomes: lambda_mut must be finite and >= 0, got {lambda_mut}"
        )));
    }
    if max_generations > 1000 {
        return Err(PyValueError::new_err(format!(
            "simulate_family_outcomes: max_generations must be <= 1000, got {max_generations}"
        )));
    }
    if !(beta.is_finite() && beta >= 0.0) {
        return Err(PyValueError::new_err(format!(
            "simulate_family_outcomes: beta must be finite and >= 0, got {beta}"
        )));
    }
    if !(selection_strength.is_finite() && selection_strength >= 0.0) {
        return Err(PyValueError::new_err(format!(
            "simulate_family_outcomes: selection_strength must be finite and >= 0, got {selection_strength}"
        )));
    }

    let founder_sim = founder.inner.final_simulation().clone();
    let founder_trace = founder.inner.trace.clone();

    let kernel = S5FKernel::new(mutability, substitution);
    let mutator = S5FMutationPass::new_rate(kernel, rate);
    let params = BranchingParams {
        lambda_base,
        lambda_mut,
        max_generations,
        n_max,
        n_sample,
        seed,
    };

    // Optional affinity model (same logic as simulate_lineage's non-neutral path).
    let model: Option<AffinityModel> = if selection_strength == 0.0 && target_aa.is_none() {
        None
    } else {
        let founder_aa = sim_to_aa(&founder_sim);
        let target = match target_aa {
            Some(s) => s.into_bytes(),
            None => {
                let mut trng = Rng::new(seed ^ 0x7461_7267_6574_0001);
                make_mature_target(&founder_aa, mature_substitutions, &mut trng)
            }
        };
        let weights = vec![1.0; founder_aa.len()];
        Some(AffinityModel::new(target, weights, beta, selection_strength, &founder_aa))
    };

    let (tree, sims) = simulate_family_sims(&founder_sim, &params, &mutator, model.as_ref());

    // The branching loop mutates child sims with a bare `Pass::execute`
    // (no `LiveCallRefreshHook`), so each node sim's `segment_calls`
    // sidecar still reflects the FOUNDER's alignment, not its own
    // mutated pool. Recompute the V/D/J live calls from each observed
    // node's pool before building its Outcome; otherwise the
    // synthesized record carries stale v/d/j calls and validate_record
    // false-fails with AlleleCallTieSetMismatch under heavy SHM.
    let reference_index = ReferenceMatchIndex::build(refdata.inner());

    let node_outcomes: Vec<Option<Outcome>> = tree
        .nodes
        .iter()
        .map(|n| {
            if n.observed {
                let refreshed = refresh_live_calls(sims[n.id as usize].clone(), &reference_index);
                let s = &refreshed;
                // Synthesize a MUTATE_S5F EventRecord encoding every
                // pool position where base != germline as a
                // BaseChanged event.  This makes the Outcome
                // self-consistent: build_airr_record reads these
                // events to derive per-segment / per-subregion
                // mutation counts, and the net count is stored in
                // mutation_count so the NMutationsMismatch validator
                // check passes.
                let (event_record, net) = synthesize_shm_event_record(s);
                let s2 = s.with_mutation_count(net);
                Some(Outcome {
                    revisions: vec![s2],
                    pass_names: Vec::new(),
                    trace: founder_trace.clone(),
                    events: vec![event_record],
                })
            } else {
                None
            }
        })
        .collect();

    Ok(PyFamilyOutcome { tree, node_outcomes })
}
