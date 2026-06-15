//! PyO3 bindings for the clonal lineage engine: `LineageNode`, `LineageTree`,
//! and the `simulate_lineage` entry point.

use pyo3::prelude::*;

use crate::lineage::export::{to_fasta, to_newick, to_node_table_tsv};
use crate::lineage::tree::{LineageNode, LineageTree};
use crate::lineage::{simulate_family, BranchingParams};
use crate::passes::S5FMutationPass;
use crate::s5f::S5FKernel;

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
    lambda_base, lambda_mut, max_generations, n_max, n_sample, seed
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
    let tree = simulate_family(&founder.inner, &params, &mutator);
    Ok(PyLineageTree::new(tree))
}
