//! `PySimulation` — read-only view of [`crate::ir::Simulation`].

use pyo3::prelude::*;
use pyo3::types::PyBytes;

use crate::ir::{NucHandle, Segment, Simulation};

use super::region::PyRegion;

/// A read-only snapshot of a simulation IR revision.
///
/// Exposes the assembled nucleotide pool (as `bytes`) and the
/// region structure. V/D/J/NP/C metadata sits behind accessor
/// methods so the Python surface stays narrow and stable.
#[pyclass(name = "Simulation", module = "GenAIRR._engine", frozen)]
pub struct PySimulation {
    pub(crate) inner: Simulation,
}

impl PySimulation {
    pub(crate) fn new(inner: Simulation) -> Self {
        Self { inner }
    }

    /// One byte per pool position, `1` where `flag` is set.
    ///
    /// Outside the `#[pymethods]` block deliberately: the mask accessors
    /// that call it are the public surface, and a caller passing a raw
    /// `NucFlags` would be building on the bit layout rather than on a
    /// named provenance question.
    fn flag_mask<'py>(
        &self,
        py: Python<'py>,
        flag: crate::ir::NucFlags,
    ) -> Bound<'py, PyBytes> {
        let bytes: Vec<u8> = self
            .inner
            .pool
            .as_slice()
            .iter()
            .map(|n| u8::from(n.flags.contains(flag)))
            .collect();
        PyBytes::new_bound(py, &bytes)
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
    /// For duplex seeds this is certified alignment placement, potentially
    /// supported by the opposite strand; it is not physical birth ancestry.
    /// The duplex projector exports chemical origins separately.
    fn germline_position(&self, index: u32) -> Option<u16> {
        let n = self.inner.pool.get(NucHandle::new(index))?;
        n.germline_pos.get()
    }

    /// Positions that are palindromic (P) additions, one byte per pool
    /// position: `1` where the flag is set, `0` elsewhere.
    ///
    /// Exposed because letter case is not a substitute for it. A P
    /// nucleotide is written lowercase, but so is every base
    /// `sequencing_errors` corrupts, and hypermutation writes uppercase --
    /// so inferring P-versus-N from case misreads a corrupted TdT addition
    /// as palindromic and a mutated palindrome as TdT. Measured over 200
    /// human IGH records, `sequencing_errors(count=8)` turned 54 of 382
    /// pure-N untemplated runs into apparently mixed ones.
    ///
    /// A mask rather than raw flag bits: the bit layout includes
    /// `INDEL_OUTSIDE_REGIONS`, which exists to work around a counter
    /// attribution problem and is not something a caller should build on.
    fn p_nucleotide_mask<'py>(&self, py: Python<'py>) -> Bound<'py, PyBytes> {
        self.flag_mask(py, crate::ir::flag::P_NUC)
    }

    /// Positions that are TdT (N) additions. See
    /// [`PySimulation.p_nucleotide_mask`] for why this is a flag rather
    /// than a reading of the base's case.
    fn tdt_addition_mask<'py>(&self, py: Python<'py>) -> Bound<'py, PyBytes> {
        self.flag_mask(py, crate::ir::flag::N_NUC)
    }

    /// Positions written by a polymerase indel pass, which have no germline
    /// provenance and are not part of the recombination. These sit inside
    /// the np interval when the insertion landed between two segments, which
    /// is why the np length alone does not report the TdT amount.
    fn indel_inserted_mask<'py>(&self, py: Python<'py>) -> Bound<'py, PyBytes> {
        self.flag_mask(py, crate::ir::flag::INDEL_INSERTED)
    }

    /// Positions inside the junction, V anchor codon through J anchor codon.
    fn junction_mask<'py>(&self, py: Python<'py>) -> Bound<'py, PyBytes> {
        self.flag_mask(py, crate::ir::flag::JUNCTION)
    }

    /// Positions copied from a D segment that was assembled in inverted
    /// orientation.
    fn inverted_mask<'py>(&self, py: Python<'py>) -> Bound<'py, PyBytes> {
        self.flag_mask(py, crate::ir::flag::INVERTED)
    }

    /// Number of regions assembled into this simulation.
    fn region_count(&self) -> usize {
        self.inner.sequence.region_count()
    }

    /// All assembled regions in 5'→3' order. Each entry is a
    /// [`PyRegion`].
    ///
    /// The codon rail (amino acids + stop-codon positions) is not
    /// stored on [`crate::ir::Region`] — we compute it on demand here
    /// at the Python boundary so `Region.amino_acids()` continues to
    /// work for offline tooling, without paying the cost on every
    /// internal pass.
    fn regions(&self) -> Vec<PyRegion> {
        use crate::ir::compute_codon_rail;
        let pool = &self.inner.pool;
        self.inner
            .sequence
            .regions
            .iter()
            .map(|r| PyRegion::new(r.clone(), compute_codon_rail(r, pool)))
            .collect()
    }

    /// Allele identifier assigned to the V slot, or `None` if no
    /// V allele has been sampled yet.
    fn v_allele_id(&self) -> Option<u32> {
        self.inner
            .assignments
            .get(Segment::V)
            .copied()
            .map(|i| i.allele_id.index())
    }

    /// Allele identifier assigned to the D slot, or `None`.
    fn d_allele_id(&self) -> Option<u32> {
        self.inner
            .assignments
            .get(Segment::D)
            .copied()
            .map(|i| i.allele_id.index())
    }

    /// Allele identifier assigned to the J slot, or `None`.
    fn j_allele_id(&self) -> Option<u32> {
        self.inner
            .assignments
            .get(Segment::J)
            .copied()
            .map(|i| i.allele_id.index())
    }

    fn __repr__(&self) -> String {
        format!(
            "<Simulation pool_len={} regions={}>",
            self.inner.pool.len(),
            self.inner.sequence.region_count(),
        )
    }
}
