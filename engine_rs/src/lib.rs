//! GenAIRR V6 engine — Rust kernel.
//!
//! This crate implements the V6 simulation architecture described in
//! `.private/engine_v6_living_design_2026-05-05.md`. It is intentionally
//! decoupled from the V5 C engine at `src/GenAIRR/_native/csrc/` and
//! coexists with it during the migration phases.
//!
//! Architecture entry points (added incrementally as phases land):
//! - Phase A: project skeleton + entity types + persistent IR
//! - Phase B: pass manager + addressed traces
//! - Phase C: recombination vertical slice
//! - Phase D-G: see the design doc.

use pyo3::prelude::*;

pub mod airr_record;
pub mod assignment;
pub mod contract;
pub mod dist;
pub mod ir;
pub mod junction;
pub mod pass;
pub mod passes;
pub mod python;
pub mod refdata;
pub mod rng;
pub mod s5f;
pub mod trace;

/// Crate version string. Smoke-test target for the Phase A.1 scaffolding —
/// proves the Cargo → maturin → Python import chain is functional end to end.
#[pyfunction]
fn version() -> &'static str {
    env!("CARGO_PKG_VERSION")
}

/// PyO3 module entry point. Python sees this module as `genairr_engine`.
#[pymodule]
fn genairr_engine(m: &Bound<'_, PyModule>) -> PyResult<()> {
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
