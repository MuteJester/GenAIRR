//! PyO3 bindings for the clonal lineage engine: `LineageNode`, `LineageTree`,
//! and the `simulate_lineage` entry point.

use pyo3::prelude::*;

use crate::lineage::export::{to_fasta, to_newick, to_node_table_tsv};
use crate::lineage::tree::{LineageNode, LineageTree};

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
