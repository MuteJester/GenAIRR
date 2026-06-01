//! `PyTraceFile` — the Python-facing durable trace artifact.
//!
//! Mirrors [`crate::trace_file::TraceFile`] one-to-one. The Python
//! surface lets callers persist a run's trace to JSON, reload it on
//! a future invocation, and rerun the engine to reproduce the
//! original `Outcome`.
//!
//! See [`crate::trace_file`] for the wire-format design and the
//! distinction between "rerun-by-seed" (this slice, Option A) and
//! true trace-injected replay (Option B, follow-up).

use pyo3::exceptions::PyValueError;
use pyo3::prelude::*;

use crate::trace_file::{TraceFile, TraceFileError, TRACE_FILE_SCHEMA_VERSION};

use super::trace::PyTrace;

#[pyclass(name = "TraceFile", module = "GenAIRR._engine", frozen)]
pub struct PyTraceFile {
    pub(crate) inner: TraceFile,
}

impl PyTraceFile {
    pub(crate) fn new(inner: TraceFile) -> Self {
        Self { inner }
    }
}

fn trace_file_error_to_pyerr(e: TraceFileError) -> PyErr {
    PyValueError::new_err(e.to_string())
}

#[pymethods]
impl PyTraceFile {
    /// Stable on-disk schema version this engine supports.
    ///
    /// Callers can compare against
    /// `PyTraceFile.SCHEMA_VERSION` to detect forward-incompatible
    /// trace files before loading.
    #[classattr]
    const SCHEMA_VERSION: u32 = TRACE_FILE_SCHEMA_VERSION;

    /// On-disk schema version of this specific trace file.
    #[getter]
    fn schema_version(&self) -> u32 {
        self.inner.schema_version
    }

    /// Engine crate version that produced this trace file.
    #[getter]
    fn engine_version(&self) -> &str {
        &self.inner.engine_version
    }

    /// Seed used to produce the original simulation.
    #[getter]
    fn seed(&self) -> u64 {
        self.inner.seed
    }

    /// Canonical signature of the pass plan that produced the trace.
    #[getter]
    fn pass_plan_signature(&self) -> &str {
        &self.inner.pass_plan_signature
    }

    /// Canonical signature of the refdata that produced the trace.
    #[getter]
    fn refdata_signature(&self) -> &str {
        &self.inner.refdata_signature
    }

    /// Refdata content hash (sha256, v2+ only). `None` on v1 traces.
    #[getter]
    fn refdata_content_hash(&self) -> Option<String> {
        self.inner.refdata_content_hash.clone()
    }

    /// The recorded choice trace.
    #[getter]
    fn trace(&self) -> PyTrace {
        PyTrace::new(self.inner.trace.clone())
    }

    /// Serialise to a pretty-printed JSON string.
    fn to_json(&self) -> PyResult<String> {
        self.inner.to_json_pretty().map_err(trace_file_error_to_pyerr)
    }

    /// Deserialise from a JSON string. Rejects unsupported schema versions.
    #[staticmethod]
    fn from_json(s: &str) -> PyResult<Self> {
        TraceFile::from_json(s)
            .map(Self::new)
            .map_err(trace_file_error_to_pyerr)
    }

    /// Write to disk as pretty-printed JSON.
    fn write_to(&self, path: &str) -> PyResult<()> {
        self.inner.write_to(path).map_err(trace_file_error_to_pyerr)
    }

    /// Read from disk and parse.
    #[staticmethod]
    fn read_from(path: &str) -> PyResult<Self> {
        TraceFile::read_from(path)
            .map(Self::new)
            .map_err(trace_file_error_to_pyerr)
    }

    /// Confirm every recorded address parses to a built-in
    /// [`crate::address::ChoiceAddress`]. Raises `ValueError` on
    /// the first non-parsing address found. `allow_custom_addresses=True`
    /// skips the check (the trace may carry addresses from user-defined
    /// passes).
    #[pyo3(signature = (*, allow_custom_addresses = false))]
    fn validate_addresses(&self, allow_custom_addresses: bool) -> PyResult<()> {
        self.inner
            .validate_addresses(allow_custom_addresses)
            .map_err(trace_file_error_to_pyerr)
    }

    fn __repr__(&self) -> String {
        format!(
            "<TraceFile schema_version={} engine_version={:?} seed={} trace_len={}>",
            self.inner.schema_version,
            self.inner.engine_version,
            self.inner.seed,
            self.inner.trace.len(),
        )
    }
}
