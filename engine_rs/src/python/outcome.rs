//! `PyOutcome` — read-only view of [`crate::pass::Outcome`].

use pyo3::exceptions::PyIndexError;
use pyo3::prelude::*;
use pyo3::types::PyDict;

use crate::airr_record::{build_airr_record, AirrRecord};
use crate::pass::Outcome;

use super::refdata::PyRefDataConfig;
use super::simulation::PySimulation;
use super::trace::PyTrace;

/// The result of executing a `PassPlan`: the IR revision history,
/// the per-revision pass names, and the addressed-choice trace.
#[pyclass(name = "Outcome", module = "genairr_engine", frozen)]
pub struct PyOutcome {
    pub(crate) inner: Outcome,
}

impl PyOutcome {
    pub(crate) fn new(inner: Outcome) -> Self {
        Self { inner }
    }
}

#[pymethods]
impl PyOutcome {
    /// Number of revisions (= initial + one per executed pass).
    fn revision_count(&self) -> usize {
        self.inner.revisions.len()
    }

    /// Revision at the given index. `0` is the initial simulation
    /// (input to the first pass); `i+1` is the output of the `i`-th
    /// pass. Raises `IndexError` when `index` is out of range.
    fn revision(&self, index: usize) -> PyResult<PySimulation> {
        self.inner
            .revisions
            .get(index)
            .cloned()
            .map(PySimulation::new)
            .ok_or_else(|| {
                PyIndexError::new_err(format!(
                    "revision index {} out of range (have {})",
                    index,
                    self.inner.revisions.len(),
                ))
            })
    }

    /// The final IR revision after every pass has run.
    fn final_simulation(&self) -> PySimulation {
        PySimulation::new(
            self.inner
                .revisions
                .last()
                .cloned()
                .expect("Outcome guarantees at least one revision"),
        )
    }

    /// Pass names in execution order. `pass_names()[i]` is the name
    /// of the pass that produced `revision(i + 1)`. Length equals the
    /// number of executed passes.
    fn pass_names(&self) -> Vec<String> {
        self.inner.pass_names.clone()
    }

    /// First revision produced by the pass with the given `name`, or
    /// `None` if no pass with that name ran.
    fn revision_after(&self, name: &str) -> Option<PySimulation> {
        self.inner
            .revision_after(name)
            .cloned()
            .map(PySimulation::new)
    }

    /// The addressed-choice trace recorded during the run.
    fn trace(&self) -> PyTrace {
        PyTrace::new(self.inner.trace.clone())
    }

    /// Build a fully-populated AIRR Rearrangement record dict from
    /// this outcome and the reference data it ran against.
    ///
    /// Returns a Python `dict` matching the field shape the Phase H
    /// Python builder produced — same ~50 columns, same coordinate
    /// convention (0-based half-open). The Python `SimulationResult`
    /// layer adds the AIRR-strict 1-based-inclusive transformation
    /// at TSV/CSV/DataFrame export time.
    ///
    /// `sequence_id` defaults to an empty string when omitted.
    #[pyo3(signature = (refdata, *, sequence_id = ""))]
    fn airr_record<'py>(
        &self,
        py: Python<'py>,
        refdata: &PyRefDataConfig,
        sequence_id: &str,
    ) -> PyResult<Bound<'py, PyDict>> {
        let rec = build_airr_record(&self.inner, refdata.inner(), sequence_id);
        airr_record_to_pydict(py, &rec)
    }

    fn __repr__(&self) -> String {
        format!(
            "<Outcome revisions={} passes={} trace_len={}>",
            self.inner.revisions.len(),
            self.inner.pass_names.len(),
            self.inner.trace.len(),
        )
    }
}

/// Convert an `AirrRecord` into a Python `dict` matching the field
/// names + types the Python `_airr_record` builder used to emit.
///
/// `Option<i64>` becomes `int | None`, `Option<f64>` becomes `float
/// | None`, `Option<bool>` becomes `bool | None`. Strings and the
/// few `bool`/`i64`/`f64` non-optional fields go through directly.
fn airr_record_to_pydict<'py>(
    py: Python<'py>,
    rec: &AirrRecord,
) -> PyResult<Bound<'py, PyDict>> {
    let dict = PyDict::new_bound(py);

    // AIRR metadata
    dict.set_item("sequence_id", &rec.sequence_id)?;
    dict.set_item("sequence", &rec.sequence)?;
    dict.set_item("sequence_aa", &rec.sequence_aa)?;
    dict.set_item("sequence_alignment", &rec.sequence_alignment)?;
    dict.set_item("germline_alignment", &rec.germline_alignment)?;
    dict.set_item("germline_alignment_d_mask", &rec.germline_alignment_d_mask)?;
    dict.set_item("sequence_length", rec.sequence_length)?;
    dict.set_item("rev_comp", rec.rev_comp)?;
    dict.set_item("locus", &rec.locus)?;

    // V
    dict.set_item("v_call", &rec.v_call)?;
    dict.set_item("v_cigar", &rec.v_cigar)?;
    set_opt_f64(&dict, "v_score", rec.v_score)?;
    set_opt_f64(&dict, "v_identity", rec.v_identity)?;
    set_opt_f64(&dict, "v_support", rec.v_support)?;
    set_opt_i64(&dict, "v_sequence_start", rec.v_sequence_start)?;
    set_opt_i64(&dict, "v_sequence_end", rec.v_sequence_end)?;
    set_opt_i64(&dict, "v_alignment_start", rec.v_alignment_start)?;
    set_opt_i64(&dict, "v_alignment_end", rec.v_alignment_end)?;
    set_opt_i64(&dict, "v_germline_start", rec.v_germline_start)?;
    set_opt_i64(&dict, "v_germline_end", rec.v_germline_end)?;
    dict.set_item("v_trim_5", rec.v_trim_5)?;
    dict.set_item("v_trim_3", rec.v_trim_3)?;

    // D
    dict.set_item("d_call", &rec.d_call)?;
    dict.set_item("d_cigar", &rec.d_cigar)?;
    set_opt_f64(&dict, "d_score", rec.d_score)?;
    set_opt_f64(&dict, "d_identity", rec.d_identity)?;
    set_opt_f64(&dict, "d_support", rec.d_support)?;
    set_opt_i64(&dict, "d_sequence_start", rec.d_sequence_start)?;
    set_opt_i64(&dict, "d_sequence_end", rec.d_sequence_end)?;
    set_opt_i64(&dict, "d_alignment_start", rec.d_alignment_start)?;
    set_opt_i64(&dict, "d_alignment_end", rec.d_alignment_end)?;
    set_opt_i64(&dict, "d_germline_start", rec.d_germline_start)?;
    set_opt_i64(&dict, "d_germline_end", rec.d_germline_end)?;
    dict.set_item("d_trim_5", rec.d_trim_5)?;
    dict.set_item("d_trim_3", rec.d_trim_3)?;

    // J
    dict.set_item("j_call", &rec.j_call)?;
    dict.set_item("j_cigar", &rec.j_cigar)?;
    set_opt_f64(&dict, "j_score", rec.j_score)?;
    set_opt_f64(&dict, "j_identity", rec.j_identity)?;
    set_opt_f64(&dict, "j_support", rec.j_support)?;
    set_opt_i64(&dict, "j_sequence_start", rec.j_sequence_start)?;
    set_opt_i64(&dict, "j_sequence_end", rec.j_sequence_end)?;
    set_opt_i64(&dict, "j_alignment_start", rec.j_alignment_start)?;
    set_opt_i64(&dict, "j_alignment_end", rec.j_alignment_end)?;
    set_opt_i64(&dict, "j_germline_start", rec.j_germline_start)?;
    set_opt_i64(&dict, "j_germline_end", rec.j_germline_end)?;
    dict.set_item("j_trim_5", rec.j_trim_5)?;
    dict.set_item("j_trim_3", rec.j_trim_3)?;

    // C
    dict.set_item("c_call", &rec.c_call)?;

    // Junction
    dict.set_item("junction", &rec.junction)?;
    dict.set_item("junction_aa", &rec.junction_aa)?;
    set_opt_i64(&dict, "junction_start", rec.junction_start)?;
    set_opt_i64(&dict, "junction_end", rec.junction_end)?;
    set_opt_i64(&dict, "junction_length", rec.junction_length)?;

    // NP regions
    dict.set_item("np1", &rec.np1)?;
    dict.set_item("np1_aa", &rec.np1_aa)?;
    dict.set_item("np1_length", rec.np1_length)?;
    dict.set_item("np2", &rec.np2)?;
    dict.set_item("np2_aa", &rec.np2_aa)?;
    dict.set_item("np2_length", rec.np2_length)?;

    // Functionality
    set_opt_bool(&dict, "productive", rec.productive)?;
    set_opt_bool(&dict, "vj_in_frame", rec.vj_in_frame)?;
    set_opt_bool(&dict, "stop_codon", rec.stop_codon)?;

    // SHM + corruption
    dict.set_item("n_mutations", rec.n_mutations)?;
    dict.set_item("mutation_rate", rec.mutation_rate)?;
    dict.set_item("n_pcr_errors", rec.n_pcr_errors)?;
    dict.set_item("n_quality_errors", rec.n_quality_errors)?;
    dict.set_item("n_indels", rec.n_indels)?;
    dict.set_item("is_contaminant", rec.is_contaminant)?;

    Ok(dict)
}

#[inline]
fn set_opt_i64(dict: &Bound<'_, PyDict>, key: &str, val: Option<i64>) -> PyResult<()> {
    match val {
        Some(v) => dict.set_item(key, v),
        None => dict.set_item(key, dict.py().None()),
    }
}

#[inline]
fn set_opt_f64(dict: &Bound<'_, PyDict>, key: &str, val: Option<f64>) -> PyResult<()> {
    match val {
        Some(v) => dict.set_item(key, v),
        None => dict.set_item(key, dict.py().None()),
    }
}

#[inline]
fn set_opt_bool(dict: &Bound<'_, PyDict>, key: &str, val: Option<bool>) -> PyResult<()> {
    match val {
        Some(v) => dict.set_item(key, v),
        None => dict.set_item(key, dict.py().None()),
    }
}
