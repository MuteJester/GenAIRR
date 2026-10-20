use pyo3::exceptions::PyValueError;
use pyo3::prelude::*;

#[pyfunction]
#[pyo3(signature = (left, right, operations_json, max_operations=256, max_bases=100_000))]
pub fn duplex_repair(
    left: &str,
    right: &str,
    operations_json: &str,
    max_operations: usize,
    max_bases: usize,
) -> PyResult<String> {
    let operations: Vec<crate::duplex_repair::RepairOperation> =
        serde_json::from_str(operations_json)
            .map_err(|e| PyValueError::new_err(format!("invalid duplex repair history: {e}")))?;
    let product = crate::duplex_repair::repair(left, right, &operations, max_operations, max_bases)
        .map_err(PyValueError::new_err)?;
    serde_json::to_string(&product).map_err(|e| PyValueError::new_err(e.to_string()))
}
