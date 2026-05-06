//! Python bindings for the V6 engine.
//!
//! Registers PyO3 wrapper types around the Rust IR and runtime
//! output: read-only views of `Simulation`, `Region`, `Trace`,
//! `ChoiceRecord`, and `Outcome`. Submodules:
//!
//! - [`region`] — `PyRegion`, a read-only view of [`crate::ir::Region`].
//! - [`simulation`] — `PySimulation`, a read-only view of
//!   [`crate::ir::Simulation`].
//! - [`trace`] — `PyChoiceRecord` + `PyTrace`, read-only views of
//!   [`crate::trace::ChoiceRecord`] and [`crate::trace::Trace`].
//! - [`event`] — `PyEventRecord` + `PyStateSummary`, read-only views
//!   of committed event ledger records.
//! - [`outcome`] — `PyOutcome`, a read-only view of
//!   [`crate::pass::Outcome`].
//! - [`runner`] — module-level Python entry points that build and
//!   execute plans, returning `PyOutcome`.
//!
//! Wrappers are thin newtype structs that own the inner Rust value.
//! They never expose `&mut` — the Python surface is strictly read-only
//! at this phase. Mutating bindings (e.g., for plan construction) land
//! in later F.x steps.

use pyo3::prelude::*;

pub mod contract;
pub mod event;
pub mod outcome;
pub mod plan;
pub mod refdata;
pub mod region;
pub mod runner;
pub mod simulation;
pub mod trace;

/// Register every PyO3 type and module-level function on `m`.
/// Called once from `lib.rs` during module initialisation.
pub(crate) fn register(m: &Bound<'_, PyModule>) -> PyResult<()> {
    m.add_class::<region::PyRegion>()?;
    m.add_class::<simulation::PySimulation>()?;
    m.add_class::<trace::PyChoiceRecord>()?;
    m.add_class::<trace::PyTrace>()?;
    m.add_class::<event::PyStateSummary>()?;
    m.add_class::<event::PyEventRecord>()?;
    m.add_class::<outcome::PyOutcome>()?;
    m.add_class::<refdata::PyAllele>()?;
    m.add_class::<refdata::PyRefDataConfig>()?;
    m.add_class::<plan::PyPassPlan>()?;
    contract::register(m)?;
    runner::register(m)?;
    Ok(())
}
