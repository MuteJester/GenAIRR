//! GenAIRR engine — Rust kernel.
//!
//! This crate implements the simulation architecture described in
//! `.private/engine_v6_living_design_2026-05-05.md`.

use pyo3::prelude::*;

pub mod airr_record;
pub mod assignment;
pub mod codon;
pub mod compiled;
pub mod contract;
pub mod dist;
pub mod event;
pub mod feasibility;
pub mod ir;
pub mod junction;
pub mod live_call;
pub mod pass;
pub mod passes;
pub mod python;
pub mod refdata;
pub mod rng;
pub mod s5f;
pub mod trace;

/// Crate version string. Used as a smoke-test target to confirm the
/// Cargo → maturin → Python import chain is functional end to end.
#[pyfunction]
fn version() -> &'static str {
    env!("CARGO_PKG_VERSION")
}

/// PyO3 module entry point. Python sees this module as `GenAIRR._engine`.
#[pymodule]
fn _engine(m: &Bound<'_, PyModule>) -> PyResult<()> {
    m.add_function(wrap_pyfunction!(version, m)?)?;
    python::register(m)?;
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;

    /// The version string is non-empty and matches the Cargo metadata.
    #[test]
    fn version_is_cargo_pkg_version() {
        assert!(!version().is_empty());
        assert_eq!(version(), env!("CARGO_PKG_VERSION"));
    }
}
