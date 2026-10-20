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

#[pyfunction]
#[pyo3(signature = (v, d, j, joins_json, *, d_inverted=false, max_operations=256, max_bases=100_000))]
pub fn coupled_duplex_repair(
    v: &str,
    d: &str,
    j: &str,
    joins_json: &str,
    d_inverted: bool,
    max_operations: usize,
    max_bases: usize,
) -> PyResult<String> {
    let joins: Vec<crate::duplex_repair::JoinScript> = serde_json::from_str(joins_json)
        .map_err(|e| PyValueError::new_err(format!("invalid coupled duplex history: {e}")))?;
    let product = crate::duplex_repair::repair_coupled(
        v,
        d,
        j,
        d_inverted,
        &joins,
        max_operations,
        max_bases,
    )
    .map_err(PyValueError::new_err)?;
    serde_json::to_string(&product).map_err(|e| PyValueError::new_err(e.to_string()))
}
