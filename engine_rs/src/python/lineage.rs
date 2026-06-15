//! PyO3 bindings for the clonal lineage engine: `LineageNode`, `LineageTree`,
//! and the `simulate_lineage` entry point.

use pyo3::prelude::*;
use pyo3::types::PyDict;

use crate::airr_record::build_airr_record;
use crate::ir::Segment;
use crate::ir::Simulation;
use crate::lineage::export::{to_fasta, to_newick, to_node_table_tsv};
use crate::lineage::tree::{LineageNode, LineageTree};
use crate::lineage::{simulate_family, simulate_family_sims, simulate_family_with_affinity, sim_to_aa, AffinityModel, BranchingParams};
use crate::lineage::affinity::make_mature_target;
use crate::pass::Outcome;
use crate::refdata::RefDataConfig;
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
    /// Per-segment and V-subregion mutation counts are recomputed from
    /// the node's final simulation pool (net mutations from germline:
    /// positions where `base != germline`), overwriting the zero
    /// counts that `build_airr_record` would produce from the empty
    /// event ledger carried by lineage node `Outcome`s.
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

            // Recompute mutation counts from the pool (net mutations
            // from germline) to fix the empty-event-ledger bug.
            let sim = outcome.final_simulation();
            let counts = pool_mutation_counts(sim, refdata.inner());

            dict.set_item("n_mutations", counts.n_mutations)?;
            dict.set_item("n_v_mutations", counts.n_v_mutations)?;
            dict.set_item("n_d_mutations", counts.n_d_mutations)?;
            dict.set_item("n_j_mutations", counts.n_j_mutations)?;
            dict.set_item("n_np_mutations", counts.n_np_mutations)?;
            dict.set_item("n_fwr1_mutations", counts.n_fwr1_mutations)?;
            dict.set_item("n_cdr1_mutations", counts.n_cdr1_mutations)?;
            dict.set_item("n_fwr2_mutations", counts.n_fwr2_mutations)?;
            dict.set_item("n_cdr2_mutations", counts.n_cdr2_mutations)?;
            dict.set_item("n_fwr3_mutations", counts.n_fwr3_mutations)?;
            dict.set_item("n_v_unannotated_mutations", counts.n_v_unannotated_mutations)?;

            // Recompute mutation_rate = n_mutations / sequence_length.
            let seq_len = rec.sequence_length;
            let mutation_rate = if seq_len > 0 {
                counts.n_mutations as f64 / seq_len as f64
            } else {
                0.0
            };
            dict.set_item("mutation_rate", mutation_rate)?;

            results.push(dict);
        }
        Ok(results)
    }
}

// ──────────────────────────────────────────────────────────────────
// Pool-derived mutation counter helper
// ──────────────────────────────────────────────────────────────────

/// Per-segment and V-subregion mutation counts derived by scanning
/// the simulation pool for positions where `base != germline`.
/// This is the "net mutations from germline" quantity that lineage
/// tools use, computed without relying on the event ledger.
struct PoolMutCounts {
    n_mutations: i64,
    n_v_mutations: i64,
    n_d_mutations: i64,
    n_j_mutations: i64,
    n_np_mutations: i64,
    n_fwr1_mutations: i64,
    n_cdr1_mutations: i64,
    n_fwr2_mutations: i64,
    n_cdr2_mutations: i64,
    n_fwr3_mutations: i64,
    n_v_unannotated_mutations: i64,
}

/// Scan `sim.pool` and count positions where `base != germline`,
/// bucketing by segment and (for V) by V-subregion. Mirrors the
/// event-loop logic in `airr_record/builder.rs` exactly: same
/// `v_subregions` table lookup, same `germline_pos.get()` → Option<u16>
/// coordinate, same fallthrough to `n_v_unannotated_mutations`.
fn pool_mutation_counts(sim: &Simulation, refdata: &RefDataConfig) -> PoolMutCounts {
    // Hoist the V-allele subregion table out of the per-nucleotide
    // loop (invariant within a record) — same pattern as builder.rs.
    let v_subregions: Option<&[crate::refdata::VSubregion]> = sim
        .assignments
        .get(Segment::V)
        .and_then(|inst| refdata.v_pool.get(inst.allele_id))
        .map(|allele| allele.subregions.as_slice());

    let mut n_v_mutations = 0i64;
    let mut n_d_mutations = 0i64;
    let mut n_j_mutations = 0i64;
    let mut n_np_mutations = 0i64;
    let mut n_fwr1_mutations = 0i64;
    let mut n_cdr1_mutations = 0i64;
    let mut n_fwr2_mutations = 0i64;
    let mut n_cdr2_mutations = 0i64;
    let mut n_fwr3_mutations = 0i64;
    let mut n_v_unannotated_mutations = 0i64;

    for nuc in sim.pool.as_slice() {
        // A position is mutated iff its current base differs from germline.
        // No case folding — the engine stores mutations as lowercase bytes
        // (SHM traces lowercase), so a raw byte comparison is correct and
        // mirrors `base != germline` in the event-loop path.
        if nuc.base == nuc.germline {
            continue;
        }
        match nuc.segment {
            Segment::V => {
                n_v_mutations += 1;
                // Dispatch into the V-subregion partition using the same
                // logic as builder.rs:321-347. `germline_pos.get()` yields
                // the allele-relative Option<u16> coordinate that the
                // VSubregion [start, end) intervals are keyed on.
                let label = v_subregions.and_then(|subs| {
                    nuc.germline_pos.get().and_then(|pos| {
                        subs.iter()
                            .find(|s| s.start <= pos && pos < s.end)
                            .map(|s| s.label)
                    })
                });
                match label {
                    Some(crate::refdata::VSubregionLabel::Fwr1) => n_fwr1_mutations += 1,
                    Some(crate::refdata::VSubregionLabel::Cdr1) => n_cdr1_mutations += 1,
                    Some(crate::refdata::VSubregionLabel::Fwr2) => n_fwr2_mutations += 1,
                    Some(crate::refdata::VSubregionLabel::Cdr2) => n_cdr2_mutations += 1,
                    Some(crate::refdata::VSubregionLabel::Fwr3) => n_fwr3_mutations += 1,
                    None => n_v_unannotated_mutations += 1,
                }
            }
            Segment::D => n_d_mutations += 1,
            Segment::J => n_j_mutations += 1,
            Segment::Np1 | Segment::Np2 => n_np_mutations += 1,
        }
    }

    let n_mutations = n_v_mutations + n_d_mutations + n_j_mutations + n_np_mutations;
    PoolMutCounts {
        n_mutations,
        n_v_mutations,
        n_d_mutations,
        n_j_mutations,
        n_np_mutations,
        n_fwr1_mutations,
        n_cdr1_mutations,
        n_fwr2_mutations,
        n_cdr2_mutations,
        n_fwr3_mutations,
        n_v_unannotated_mutations,
    }
}

/// Grow + sample a clonal lineage family from `founder` (a full `Outcome`) using
/// an S5F mutator, returning a `FamilyOutcome` with per-node AIRR-projectable
/// `Outcome`s for every observed (sampled) node.
#[pyfunction]
#[pyo3(signature = (
    founder, mutability, substitution, rate,
    lambda_base, lambda_mut, max_generations, n_max, n_sample, seed,
    selection_strength=0.0, beta=1.0, target_aa=None, mature_substitutions=5
))]
#[allow(clippy::too_many_arguments)]
pub(crate) fn simulate_family_outcomes(
    founder: &PyOutcome,
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

    let node_outcomes: Vec<Option<Outcome>> = tree
        .nodes
        .iter()
        .map(|n| {
            if n.observed {
                Some(Outcome {
                    revisions: vec![sims[n.id as usize].clone()],
                    pass_names: Vec::new(),
                    trace: founder_trace.clone(),
                    events: Vec::new(),
                })
            } else {
                None
            }
        })
        .collect();

    Ok(PyFamilyOutcome { tree, node_outcomes })
}
