//! `PyOutcome` — read-only view of [`crate::pass::Outcome`].

use pyo3::exceptions::PyIndexError;
use pyo3::prelude::*;

use crate::pass::Outcome;

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

    fn __repr__(&self) -> String {
        format!(
            "<Outcome revisions={} passes={} trace_len={}>",
            self.inner.revisions.len(),
            self.inner.pass_names.len(),
            self.inner.trace.len(),
        )
    }
}
