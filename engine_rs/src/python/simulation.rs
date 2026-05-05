//! `PySimulation` — read-only view of [`crate::ir::Simulation`].

use pyo3::prelude::*;
use pyo3::types::PyBytes;

use crate::ir::{NucHandle, Simulation};

use super::region::PyRegion;

/// A read-only snapshot of a simulation IR revision.
///
/// Exposes the assembled nucleotide pool (as `bytes`) and the
/// region structure. V/D/J/NP/C metadata sits behind accessor
/// methods so the Python surface stays narrow and stable.
#[pyclass(name = "Simulation", module = "genairr_engine", frozen)]
pub struct PySimulation {
    pub(crate) inner: Simulation,
}

impl PySimulation {
    pub(crate) fn new(inner: Simulation) -> Self {
        Self { inner }
    }
}

#[pymethods]
impl PySimulation {
    /// Number of nucleotides in the pool.
    fn __len__(&self) -> usize {
        self.inner.pool.len()
    }

    /// Whether the pool contains zero nucleotides.
    fn is_empty(&self) -> bool {
        self.inner.pool.is_empty()
    }

    /// Concatenated pool bases as a `bytes` object. The result is
    /// ASCII-encoded — uppercase bases mark germline-derived
    /// positions, lowercase mark mutated / corrupted positions.
    fn bases<'py>(&self, py: Python<'py>) -> Bound<'py, PyBytes> {
        let bytes: Vec<u8> = self.inner.pool.as_slice().iter().map(|n| n.base).collect();
        PyBytes::new_bound(py, &bytes)
    }

    /// Original germline base at each pool position. For synthetic
    /// nucleotides (NP, P-nucs, indel insertions) the germline byte
    /// is whatever the constructor stored — typically `b'\0'` or the
    /// inserted base itself; downstream code reading germline
    /// provenance should use [`PySimulation.germline_position`]
    /// instead.
    fn germline_bases<'py>(&self, py: Python<'py>) -> Bound<'py, PyBytes> {
        let bytes: Vec<u8> = self
            .inner
            .pool
            .as_slice()
            .iter()
            .map(|n| n.germline)
            .collect();
        PyBytes::new_bound(py, &bytes)
    }

    /// Position-in-source-allele for the nucleotide at the given
    /// pool index. Returns `None` for synthetic nucleotides (NP /
    /// indel-inserted) that have no germline provenance.
    fn germline_position(&self, index: u32) -> Option<u16> {
        let n = self.inner.pool.get(NucHandle::new(index))?;
        if n.germline_pos == crate::ir::Nucleotide::NO_GERMLINE_POS {
            None
        } else {
            Some(n.germline_pos)
        }
    }

    /// Number of regions assembled into this simulation.
    fn region_count(&self) -> usize {
        self.inner.sequence.region_count()
    }

    /// All assembled regions in 5'→3' order. Each entry is a
    /// [`PyRegion`].
    fn regions(&self) -> Vec<PyRegion> {
        self.inner
            .sequence
            .regions
            .iter()
            .cloned()
            .map(PyRegion::new)
            .collect()
    }

    /// Allele identifier assigned to the V slot, or `None` if no
    /// V allele has been sampled yet.
    fn v_allele_id(&self) -> Option<u32> {
        self.inner.assignments.v.map(|i| i.allele_id.index())
    }

    /// Allele identifier assigned to the D slot, or `None`.
    fn d_allele_id(&self) -> Option<u32> {
        self.inner.assignments.d.map(|i| i.allele_id.index())
    }

    /// Allele identifier assigned to the J slot, or `None`.
    fn j_allele_id(&self) -> Option<u32> {
        self.inner.assignments.j.map(|i| i.allele_id.index())
    }

    fn __repr__(&self) -> String {
        format!(
            "<Simulation pool_len={} regions={}>",
            self.inner.pool.len(),
            self.inner.sequence.region_count(),
        )
    }
}
