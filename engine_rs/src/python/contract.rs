//! Python bindings for the contract subsystem.
//!
//! Exposes:
//! - [`PyContractSet`] — a read-only handle around the Rust
//!   [`ContractSet`] with introspection (`len`, `names`, `is_empty`).
//! - The [`productive()`] module-level factory that returns the
//!   canonical productive-sequence bundle.
//! - The [`StrictSamplingError`] Python exception class raised from
//!   strict-mode runs when a candidate cannot be sampled, a pass
//!   precondition is not satisfied, or a post-pass contract fence
//!   fails.
//!
//! Custom contract construction from Python (e.g. building a
//! `ContractSet` allele-by-allele) is intentionally **not** in this
//! phase: every contract our users currently want is captured by the
//! `productive()` bundle. Bespoke contract authoring lands later if
//! the need arises.

use pyo3::create_exception;
use pyo3::exceptions::PyException;
use pyo3::prelude::*;

use crate::contract::{productive as rust_productive, ContractSet};
use crate::dist::FilteredSampleError;
use crate::pass::PassError;

create_exception!(
    _engine,
    StrictSamplingError,
    PyException,
    "Structured sampling, replay, or contract validation failure.\n\
     Args are (pass_name, address, stable_reason). The reason \
     sampling_budget_exceeded means support remains unresolved, never empty; \
     it is raised in strict and permissive execution, with requested_count, \
     attempts, explored_states, state_budget, and strategy attributes.\n\
     Replay validates positive configured support and the final contracts. \
     Historical all-no-op indel traces require explicit legacy compatibility.\n\
     This exception extends Exception directly, not ValueError. Compile-time \
     static support failures remain ValueError."
);

/// Convert a [`PassError`] into a Python [`StrictSamplingError`]
/// with the pass name, address, and a stable lowercase reason
/// string in `args`.
pub(crate) fn pass_error_to_pyerr(err: PassError) -> PyErr {
    match err {
        PassError::SamplingBudgetExceeded {
            pass_name,
            address,
            requested_count,
            attempts,
            explored_states,
            state_budget,
            strategy,
        } => {
            let error =
                StrictSamplingError::new_err((pass_name, address, "sampling_budget_exceeded"));
            Python::with_gil(|py| {
                let value = error.value_bound(py);
                let _ = value.setattr("requested_count", requested_count);
                let _ = value.setattr("attempts", attempts);
                let _ = value.setattr("explored_states", explored_states);
                let _ = value.setattr("state_budget", state_budget);
                let _ = value.setattr("strategy", strategy);
            });
            error
        }
        PassError::ConstraintSampling {
            pass_name,
            address,
            reason,
        } => {
            let reason_str = match reason {
                FilteredSampleError::SupportUnavailable => "support_unavailable",
                FilteredSampleError::EmptyAdmissibleSupport => "empty_admissible_support",
                FilteredSampleError::InvalidFilteredSupport => "invalid_filtered_support",
            };
            StrictSamplingError::new_err((pass_name, address, reason_str.to_string()))
        }
        PassError::MissingRefData { pass_name } => {
            StrictSamplingError::new_err((pass_name, String::new(), "missing_refdata".to_string()))
        }
        PassError::MissingAssignment { pass_name, segment } => StrictSamplingError::new_err((
            pass_name,
            String::new(),
            format!("missing_assignment.{:?}", segment).to_ascii_lowercase(),
        )),
        PassError::MissingAllele {
            pass_name,
            segment,
            allele_id,
        } => StrictSamplingError::new_err((
            pass_name,
            String::new(),
            format!("missing_allele.{:?}.{}", segment, allele_id).to_ascii_lowercase(),
        )),
        PassError::InvalidDistributionOutput {
            pass_name,
            address,
            reason,
            ..
        } => StrictSamplingError::new_err((
            pass_name,
            address,
            format!("invalid_distribution_output.{}", reason),
        )),
        PassError::InvalidPlanState { pass_name, reason } => StrictSamplingError::new_err((
            pass_name,
            String::new(),
            format!("invalid_plan_state.{}", reason),
        )),
        PassError::ContractViolation {
            pass_name,
            violations,
        } => {
            let contract_names = violations
                .iter()
                .map(|violation| violation.contract_name.as_str())
                .collect::<Vec<_>>()
                .join(",");
            StrictSamplingError::new_err((
                pass_name,
                String::new(),
                format!("contract_violation.{}", contract_names),
            ))
        }
        PassError::Replay { pass_name, reason } => {
            // `ReplayError`'s `Display` already produces a structured
            // diagnostic; surface it as the reason payload with a
            // stable `replay.<kind>` prefix so downstream code can
            // dispatch on it without parsing the human-readable text.
            let reason_tag = match reason {
                crate::replay::ReplayError::Exhausted { .. } => "replay.exhausted",
                crate::replay::ReplayError::AddressMismatch { .. } => "replay.address_mismatch",
                crate::replay::ReplayError::ValueKindMismatch { .. } => {
                    "replay.value_kind_mismatch"
                }
                crate::replay::ReplayError::UnusedTrailingRecords { .. } => {
                    "replay.unused_trailing_records"
                }
            };
            StrictSamplingError::new_err((
                pass_name,
                String::new(),
                format!("{}: {}", reason_tag, reason),
            ))
        }
    }
}

/// A bundle of contracts that all must hold for the simulation.
///
/// Returned by factory functions like [`productive()`]. Pass to
/// runners as the `respect=` argument; the compiled simulator
/// threads it through each pass so sampling can filter candidates
/// and strict execution can verify post-pass state.
#[pyclass(name = "ContractSet", module = "GenAIRR._engine", unsendable)]
pub struct PyContractSet {
    pub(crate) inner: ContractSet,
}

impl PyContractSet {
    pub(crate) fn inner(&self) -> &ContractSet {
        &self.inner
    }
}

#[pymethods]
impl PyContractSet {
    /// Number of contracts in the bundle.
    fn __len__(&self) -> usize {
        self.inner.len()
    }

    /// Whether the bundle contains zero contracts.
    fn is_empty(&self) -> bool {
        self.inner.is_empty()
    }

    /// Stable names of the contracts in the bundle, in order.
    fn names(&self) -> Vec<String> {
        self.inner.iter().map(|c| c.name().to_string()).collect()
    }

    fn __repr__(&self) -> String {
        format!(
            "<ContractSet len={} names={:?}>",
            self.inner.len(),
            self.names(),
        )
    }
}

/// The canonical productive-sequence contract bundle.
///
/// Composes:
/// 1. `ProductiveJunctionFrame` — junction length divisible by 3
/// 2. `NoStopCodonInJunction` — no stops inside the junction
/// 3. Canonical V anchor — original landmark retained, encoding Cys
/// 4. Canonical J anchor — original landmark retained, encoding Phe/Trp
/// 5. `CodingSequenceIntegrity` — stop-free retained coding sequence and
///    preservation of the original V coding phase through the J tail
///
/// Pass this to ``Experiment.run(respect=...)`` (or to
/// ``GenAIRR._engine.run(plan, seed, respect=...)``) to constrain
/// every NP-base draw, every length sample, and every mutation /
/// corruption substitution to admissible values. In strict mode, the
/// compiled simulator also verifies the bundle after each pass. Both modes
/// require a final accepted `canonical_observed_vdj_v1` assessment.
#[pyfunction]
pub fn productive() -> PyContractSet {
    PyContractSet {
        inner: rust_productive(),
    }
}

pub(crate) fn register(m: &Bound<'_, PyModule>) -> PyResult<()> {
    m.add_class::<PyContractSet>()?;
    m.add_function(wrap_pyfunction!(productive, m)?)?;
    m.add(
        "StrictSamplingError",
        m.py().get_type_bound::<StrictSamplingError>(),
    )?;
    Ok(())
}
