//! `PyAllele` + `PyRefDataConfig` — the Python-facing reference data.
//!
//! `PyRefDataConfig` is a *mutable builder* on the Python side: you
//! construct it for a chain type (`vj` or `vdj`) and grow each
//! segment pool by adding alleles. Once built, it's passed to a
//! runner (e.g. `run_vj_recombination`) which clones the inner
//! Rust `RefDataConfig` for the duration of the simulation.
//!
//! `PyAllele` is the corresponding read-only view returned by
//! pool accessors. F.3 keeps the Python surface narrow: enough
//! to drive recombination and inspect what was supplied, no more.

use pyo3::exceptions::{PyIndexError, PyValueError};
use pyo3::prelude::*;
use pyo3::types::PyBytes;

use crate::ir::Segment;
use crate::refdata::{Allele, AlleleId, AllelePool, ChainType, RefDataConfig};

/// A read-only view of one V / D / J allele.
#[pyclass(name = "Allele", module = "genairr_engine", frozen)]
pub struct PyAllele {
    pub(crate) inner: Allele,
}

impl PyAllele {
    pub(crate) fn new(inner: Allele) -> Self {
        Self { inner }
    }
}

#[pymethods]
impl PyAllele {
    /// IMGT allele identifier (e.g. `"IGHV1-2*01"`).
    #[getter]
    fn name(&self) -> &str {
        &self.inner.name
    }

    /// IMGT gene name (e.g. `"IGHV1-2"`).
    #[getter]
    fn gene(&self) -> &str {
        &self.inner.gene
    }

    /// Reference nucleotide sequence as ASCII bytes.
    fn seq<'py>(&self, py: Python<'py>) -> Bound<'py, PyBytes> {
        PyBytes::new_bound(py, &self.inner.seq)
    }

    /// Length of the reference sequence in bases.
    fn __len__(&self) -> usize {
        self.inner.len() as usize
    }

    /// Segment role (`"V"`, `"D"`, or `"J"`).
    #[getter]
    fn segment(&self) -> &'static str {
        match self.inner.segment {
            Segment::V => "V",
            Segment::D => "D",
            Segment::J => "J",
            Segment::Np1 | Segment::Np2 => unreachable!("alleles are only V/D/J"),
        }
    }

    /// Position of the conserved anchor codon (V Cys / J W/F) in
    /// the allele coordinate system, or `None` for anchorless
    /// pseudogenes / partial alleles.
    #[getter]
    fn anchor(&self) -> Option<u16> {
        self.inner.anchor
    }

    fn __repr__(&self) -> String {
        match self.inner.anchor {
            Some(a) => format!(
                "<Allele {} {} {}bp anchor={}>",
                self.segment(),
                self.inner.name,
                self.inner.len(),
                a
            ),
            None => format!(
                "<Allele {} {} {}bp anchor=None>",
                self.segment(),
                self.inner.name,
                self.inner.len()
            ),
        }
    }
}

/// Mutable builder for the V/D/J reference pools used by every
/// recombination plan.
///
/// Construction:
/// ```python
/// cfg = RefDataConfig.vj()
/// cfg.add_v_allele("v1*01", "v1", b"AAACCCGGG", anchor=6)
/// cfg.add_j_allele("j1*01", "j1", b"TTTAAA", anchor=0)
/// ```
#[pyclass(name = "RefDataConfig", module = "genairr_engine")]
pub struct PyRefDataConfig {
    pub(crate) inner: RefDataConfig,
}

impl PyRefDataConfig {
    /// Borrowed view of the inner Rust config — used by runners that
    /// need to thread it into a `PassRuntime`. Internal API only.
    pub(crate) fn inner(&self) -> &RefDataConfig {
        &self.inner
    }
}

#[pymethods]
impl PyRefDataConfig {
    /// Construct an empty refdata for the given chain type.
    /// `chain_type` is `"vj"` (light chain) or `"vdj"` (heavy chain),
    /// case-insensitive.
    #[new]
    fn new(chain_type: &str) -> PyResult<Self> {
        let ct = parse_chain_type(chain_type)?;
        Ok(Self {
            inner: RefDataConfig::empty(ct),
        })
    }

    /// Construct an empty VJ-chain (light) refdata.
    #[staticmethod]
    fn vj() -> Self {
        Self {
            inner: RefDataConfig::empty(ChainType::Vj),
        }
    }

    /// Construct an empty VDJ-chain (heavy) refdata.
    #[staticmethod]
    fn vdj() -> Self {
        Self {
            inner: RefDataConfig::empty(ChainType::Vdj),
        }
    }

    /// Chain type as a lowercase string (`"vj"` or `"vdj"`).
    #[getter]
    fn chain_type(&self) -> &'static str {
        match self.inner.chain_type {
            ChainType::Vj => "vj",
            ChainType::Vdj => "vdj",
        }
    }

    /// Whether the chain type carries a D segment (i.e. is heavy).
    fn has_d(&self) -> bool {
        self.inner.chain_type.has_d()
    }

    /// Append a V allele. Returns the new allele id (`u32`).
    #[pyo3(signature = (name, gene, seq, *, anchor=None))]
    fn add_v_allele(
        &mut self,
        name: &str,
        gene: &str,
        seq: &[u8],
        anchor: Option<u16>,
    ) -> u32 {
        push_allele(&mut self.inner.v_pool, Segment::V, name, gene, seq, anchor)
    }

    /// Append a D allele. Returns the new allele id (`u32`). Panics
    /// if called on a VJ chain — D pools are only valid for heavy chains.
    #[pyo3(signature = (name, gene, seq, *, anchor=None))]
    fn add_d_allele(
        &mut self,
        name: &str,
        gene: &str,
        seq: &[u8],
        anchor: Option<u16>,
    ) -> PyResult<u32> {
        if !self.inner.chain_type.has_d() {
            return Err(PyValueError::new_err(
                "cannot add D allele to a VJ-chain RefDataConfig",
            ));
        }
        Ok(push_allele(
            &mut self.inner.d_pool,
            Segment::D,
            name,
            gene,
            seq,
            anchor,
        ))
    }

    /// Append a J allele. Returns the new allele id (`u32`).
    #[pyo3(signature = (name, gene, seq, *, anchor=None))]
    fn add_j_allele(
        &mut self,
        name: &str,
        gene: &str,
        seq: &[u8],
        anchor: Option<u16>,
    ) -> u32 {
        push_allele(&mut self.inner.j_pool, Segment::J, name, gene, seq, anchor)
    }

    /// Number of V alleles in the pool.
    fn v_pool_size(&self) -> usize {
        self.inner.v_pool.len()
    }

    /// Number of D alleles in the pool. Always `0` for VJ chains.
    fn d_pool_size(&self) -> usize {
        self.inner.d_pool.len()
    }

    /// Number of J alleles in the pool.
    fn j_pool_size(&self) -> usize {
        self.inner.j_pool.len()
    }

    /// V allele at the given id, or raises `IndexError`.
    fn v_allele(&self, id: u32) -> PyResult<PyAllele> {
        self.inner
            .v_pool
            .get(AlleleId::new(id))
            .cloned()
            .map(PyAllele::new)
            .ok_or_else(|| PyIndexError::new_err(format!("V allele id {} out of range", id)))
    }

    /// D allele at the given id, or raises `IndexError`.
    fn d_allele(&self, id: u32) -> PyResult<PyAllele> {
        self.inner
            .d_pool
            .get(AlleleId::new(id))
            .cloned()
            .map(PyAllele::new)
            .ok_or_else(|| PyIndexError::new_err(format!("D allele id {} out of range", id)))
    }

    /// J allele at the given id, or raises `IndexError`.
    fn j_allele(&self, id: u32) -> PyResult<PyAllele> {
        self.inner
            .j_pool
            .get(AlleleId::new(id))
            .cloned()
            .map(PyAllele::new)
            .ok_or_else(|| PyIndexError::new_err(format!("J allele id {} out of range", id)))
    }

    fn __repr__(&self) -> String {
        format!(
            "<RefDataConfig {} V={} D={} J={}>",
            self.chain_type(),
            self.v_pool_size(),
            self.d_pool_size(),
            self.j_pool_size(),
        )
    }
}

fn parse_chain_type(s: &str) -> PyResult<ChainType> {
    match s.to_ascii_lowercase().as_str() {
        "vj" => Ok(ChainType::Vj),
        "vdj" => Ok(ChainType::Vdj),
        other => Err(PyValueError::new_err(format!(
            "unknown chain type {:?}; expected 'vj' or 'vdj'",
            other
        ))),
    }
}

fn push_allele(
    pool: &mut AllelePool,
    segment: Segment,
    name: &str,
    gene: &str,
    seq: &[u8],
    anchor: Option<u16>,
) -> u32 {
    let id = pool.push(Allele {
        name: name.to_string(),
        gene: gene.to_string(),
        seq: seq.to_vec(),
        segment,
        anchor,
    });
    id.index()
}
