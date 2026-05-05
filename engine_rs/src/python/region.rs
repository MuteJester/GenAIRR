//! `PyRegion` — read-only view of [`crate::ir::Region`].

use pyo3::prelude::*;
use pyo3::types::PyBytes;

use crate::ir::{Region, Segment};

/// A read-only view of one region of an assembled receptor sequence.
///
/// Carries the segment role (V/D/J/Np1/Np2/C), the half-open
/// `[start, end)` pool index range, and the codon-rail metadata
/// (`frame_phase` + `amino_acids`).
#[pyclass(name = "Region", module = "genairr_engine", frozen)]
pub struct PyRegion {
    pub(crate) inner: Region,
}

impl PyRegion {
    pub(crate) fn new(inner: Region) -> Self {
        Self { inner }
    }
}

#[pymethods]
impl PyRegion {
    /// Biological role of this region, returned as a short uppercase
    /// label (`"V"`, `"D"`, `"J"`, `"NP1"`, `"NP2"`).
    #[getter]
    fn segment(&self) -> &'static str {
        match self.inner.segment {
            Segment::V => "V",
            Segment::D => "D",
            Segment::J => "J",
            Segment::Np1 => "NP1",
            Segment::Np2 => "NP2",
        }
    }

    /// Inclusive lower bound of the pool range.
    #[getter]
    fn start(&self) -> u32 {
        self.inner.start.index()
    }

    /// Exclusive upper bound of the pool range.
    #[getter]
    fn end(&self) -> u32 {
        self.inner.end.index()
    }

    /// Number of nucleotides covered by this region (`end - start`).
    fn __len__(&self) -> usize {
        self.inner.len() as usize
    }

    /// Whether the region covers zero nucleotides (e.g., zero-length NP1).
    fn is_empty(&self) -> bool {
        self.inner.is_empty()
    }

    /// Position within the codon frame at the first nucleotide of
    /// this region. `0` = codon-aligned start, `1` = second codon
    /// base, `2` = third.
    #[getter]
    fn frame_phase(&self) -> u8 {
        self.inner.frame_phase
    }

    /// Codon-rail amino acid bytes for this region. `b'*'` marks a
    /// stop codon, `b'X'` marks an ambiguous codon containing
    /// non-{A,C,G,T,U} bases.
    fn amino_acids<'py>(&self, py: Python<'py>) -> Bound<'py, PyBytes> {
        PyBytes::new_bound(py, &self.inner.amino_acids)
    }

    fn __repr__(&self) -> String {
        format!(
            "<Region {} [{}..{}) frame_phase={}>",
            self.segment(),
            self.start(),
            self.end(),
            self.frame_phase(),
        )
    }
}
