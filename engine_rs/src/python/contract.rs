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
    genairr_engine,
    StrictSamplingError,
    PyException,
    "Raised by strict-mode runs when a pass cannot execute safely.\n\
     \n\
     Exception args are a 3-tuple ``(pass_name, address, reason)`` where:\n\
     - ``pass_name`` (str) is the name of the failing pass (e.g. \
       ``\"generate_np.np1\"``);\n\
     - ``address`` (str) is the trace address when applicable (e.g. \
       ``\"np.np1.length\"``);\n\
     - ``reason`` (str) is a stable lowercase error reason."
);

/// Convert a [`PassError`] into a Python [`StrictSamplingError`]
/// with the pass name, address, and a stable lowercase reason
/// string in `args`.
pub(crate) fn pass_error_to_pyerr(err: PassError) -> PyErr {
    match err {
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
    }
}

/// A bundle of contracts that all must hold for the simulation.
///
/// Returned by factory functions like [`productive()`]. Pass to
/// runners as the `respect=` argument; the compiled simulator
/// threads it through each pass so sampling can filter candidates
/// and strict execution can verify post-pass state.
#[pyclass(name = "ContractSet", module = "genairr_engine", unsendable)]
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
/// 3. `AnchorPreserved::V` — V Cys codon retained after V trim
/// 4. `AnchorPreserved::J` — J W/F codon retained after J trim
///
/// Pass this to ``Experiment.run(respect=...)`` (or to
/// ``genairr_engine.run(plan, seed, respect=...)``) to constrain
/// every NP-base draw, every length sample, and every mutation /
/// corruption substitution to admissible values. In strict mode, the
/// compiled simulator also verifies the bundle after each pass.
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
