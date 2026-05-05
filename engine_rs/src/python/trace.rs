//! `PyChoiceRecord` + `PyTrace` — read-only views of the addressed trace.

use pyo3::prelude::*;
use pyo3::types::PyBytes;

use crate::trace::{ChoiceRecord, ChoiceValue, Trace};

/// Convert a Rust [`ChoiceValue`] to a Python value.
///
/// - `Int(n)` and `AlleleId(id)` map to native ints.
/// - `Base(b)` and `Bases(v)` map to `bytes`.
/// - `Bool(b)` maps to `bool`.
fn choice_value_to_py(py: Python<'_>, v: &ChoiceValue) -> PyObject {
    match v {
        ChoiceValue::Int(n) => n.into_py(py),
        ChoiceValue::Base(b) => PyBytes::new_bound(py, &[*b]).into_py(py),
        ChoiceValue::Bases(bs) => PyBytes::new_bound(py, bs).into_py(py),
        ChoiceValue::AlleleId(id) => id.into_py(py),
        ChoiceValue::Bool(b) => b.into_py(py),
    }
}

/// One entry in the addressed trace: an `address` string and the
/// `value` that was sampled there.
#[pyclass(name = "ChoiceRecord", module = "genairr_engine", frozen)]
pub struct PyChoiceRecord {
    pub(crate) inner: ChoiceRecord,
}

impl PyChoiceRecord {
    pub(crate) fn new(inner: ChoiceRecord) -> Self {
        Self { inner }
    }
}

#[pymethods]
impl PyChoiceRecord {
    /// Hierarchical-string address (D3) at which this choice was made.
    #[getter]
    fn address(&self) -> &str {
        &self.inner.address
    }

    /// Sampled value, returned as the natural Python type
    /// (`int` / `bytes` / `bool`).
    #[getter]
    fn value(&self, py: Python<'_>) -> PyObject {
        choice_value_to_py(py, &self.inner.value)
    }

    fn __repr__(&self) -> String {
        format!("<ChoiceRecord {}={:?}>", self.inner.address, self.inner.value)
    }
}

/// The append-only record of every addressed choice made during one
/// simulation run.
#[pyclass(name = "Trace", module = "genairr_engine", frozen)]
pub struct PyTrace {
    pub(crate) inner: Trace,
}

impl PyTrace {
    pub(crate) fn new(inner: Trace) -> Self {
        Self { inner }
    }
}

#[pymethods]
impl PyTrace {
    /// Number of recorded choices.
    fn __len__(&self) -> usize {
        self.inner.len()
    }

    /// Whether no choices have been recorded.
    fn is_empty(&self) -> bool {
        self.inner.is_empty()
    }

    /// All recorded choices, in chronological order.
    fn choices(&self) -> Vec<PyChoiceRecord> {
        self.inner
            .choices()
            .iter()
            .cloned()
            .map(PyChoiceRecord::new)
            .collect()
    }

    /// First choice recorded at the given address, or `None`.
    fn find(&self, address: &str) -> Option<PyChoiceRecord> {
        self.inner.find(address).cloned().map(PyChoiceRecord::new)
    }

    /// All choices whose address starts with `prefix`, in chronological
    /// order. Useful for queries like `prefix_query("mutate.s5f.")`.
    fn prefix_query(&self, prefix: &str) -> Vec<PyChoiceRecord> {
        self.inner
            .prefix_query(prefix)
            .cloned()
            .map(PyChoiceRecord::new)
            .collect()
    }

    /// Number of choices whose address starts with `prefix`.
    fn prefix_count(&self, prefix: &str) -> usize {
        self.inner.prefix_count(prefix)
    }

    fn __repr__(&self) -> String {
        format!("<Trace len={}>", self.inner.len())
    }
}
