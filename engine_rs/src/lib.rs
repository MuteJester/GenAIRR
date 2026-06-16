//! GenAIRR engine — Rust kernel.
//!
//! This crate implements the simulation architecture described in
//! `.private/engine_v6_living_design_2026-05-05.md`.
//!
//! # Contributor-facing architecture
//!
//! Two companion docs at the repository root:
//!
//! - **`docs/engine_architecture.md`** — the architectural spine
//!   (read once, refer back to when you need to justify a design
//!   choice).
//! - **`docs/adding_a_pass.md`** — copy-pasteable template,
//!   required test patterns using
//!   [`crate::passes::test_support`] helpers, and a per-mechanism
//!   crib sheet of which existing pass to model the new one on
//!   (read each time you start a new pass).
//!
//! The architecture doc codifies:
//!
//! - The seven engine invariants (contracts narrow support before
//!   proposals; trace = choices; replay = validated proposal
//!   consumption; simulation events = consequences; compile
//!   effects = scheduling facts; live-call refresh follows events;
//!   built-in mutating passes route through [`crate::ir::SimulationBuilder`]).
//! - A six-step checklist for adding a new pass.
//! - The anti-patterns CI lockdowns and policy-conformance tests
//!   catch (direct `sim.with_*` in production passes, replay
//!   force-apply, declaring an effect without emitting matching
//!   events, etc.).
//!
//! The code enforces most rules — `boundary_lockdown`, the
//! event-emission policy test, the divergence tests — but the doc
//! teaches the right way before CI yells.

use pyo3::prelude::*;

pub mod address;
pub mod airr_record;
pub mod assignment;
pub mod codon;
pub mod compiled;
pub mod contract;
pub mod dist;
pub mod event;
pub mod feasibility;
pub mod genotype;
pub mod ir;
pub mod junction;
pub mod lineage;
pub mod live_call;
pub mod pass;
pub mod passes;
pub mod python;
pub mod refdata;
pub mod replay;
pub mod rng;
pub mod s5f;
pub mod trace;
pub mod trace_file;

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
