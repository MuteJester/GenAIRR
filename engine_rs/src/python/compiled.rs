//! Python-facing owning compiled simulator.
//!
//! `PassPlan` is a builder/IR. `CompiledSimulator` is the durable
//! execution artifact: it owns the plan plus cloned reference data and
//! contracts, so repeated `run()` calls do not recompile or reborrow
//! from Python-side builder objects.

use pyo3::exceptions::PyValueError;
use pyo3::prelude::*;

use crate::compiled::{CompiledSimulator, ExecutionPolicy, OwnedCompiledSimulator};
use crate::contract::ContractSet;
use crate::refdata::RefDataConfig;

use super::contract::pass_error_to_pyerr;
use super::simulation::PySimulation;
use super::outcome::PyOutcome;
use super::plan::PyPassPlan;

#[pyclass(name = "CompiledSimulator", module = "GenAIRR._engine", unsendable)]
pub struct PyCompiledSimulator {
    inner: OwnedCompiledSimulator,
}

impl PyCompiledSimulator {
    pub(crate) fn compile_from_plan(
        plan_builder: &mut PyPassPlan,
        refdata: Option<RefDataConfig>,
        contracts: Option<ContractSet>,
        policy: ExecutionPolicy,
    ) -> PyResult<Self> {
        let (report, feasibility) = {
            let plan = plan_builder.inner()?;
            let borrowed =
                CompiledSimulator::compile(plan, refdata.as_ref(), contracts.as_ref(), policy)
                    .map_err(|err| PyValueError::new_err(err.to_string()))?;
            (borrowed.report().clone(), borrowed.feasibility().cloned())
        };

        let plan = plan_builder.take_inner()?;
        Ok(Self {
            inner: OwnedCompiledSimulator::from_validated_parts(
                plan,
                refdata,
                contracts,
                policy,
                report,
                feasibility,
            ),
        })
    }
}

fn policy_from_strict_override(
    default_policy: ExecutionPolicy,
    strict: Option<bool>,
) -> ExecutionPolicy {
    match strict {
        Some(true) => ExecutionPolicy::Strict,
        Some(false) => ExecutionPolicy::Permissive,
        None => default_policy,
    }
}

#[pymethods]
impl PyCompiledSimulator {
    /// Runtime failure policy selected at compile time.
    fn policy(&self) -> &'static str {
        match self.inner.policy() {
            ExecutionPolicy::Permissive => "permissive",
            ExecutionPolicy::Strict => "strict",
        }
    }

    /// Stable names of the compiled pass sequence.
    fn pass_names(&self) -> Vec<String> {
        self.inner.report().pass_names()
    }

    /// Declared stochastic choices as `(pass_index, pass_name, address)`.
    fn declared_choices(&self) -> Vec<(usize, String, String)> {
        self.inner
            .report()
            .declared_choices
            .iter()
            .map(|choice| {
                (
                    choice.pass_index,
                    choice.pass_name.clone(),
                    choice.address.clone(),
                )
            })
            .collect()
    }

    /// Names of active contracts captured by this compiled simulator.
    fn active_contracts(&self) -> Vec<String> {
        self.inner.report().active_contracts.clone()
    }

    /// Non-fatal compile diagnostics.
    fn warnings(&self) -> Vec<String> {
        self.inner
            .report()
            .warnings
            .iter()
            .map(|warning| warning.message.clone())
            .collect()
    }

    /// Run one simulation using this compiled artifact.
    #[pyo3(signature = (seed, *, strict=None))]
    fn run(&self, seed: u64, strict: Option<bool>) -> PyResult<PyOutcome> {
        let policy = policy_from_strict_override(self.inner.policy(), strict);
        let outcome = self
            .inner
            .run_one_with_policy(seed, policy)
            .map_err(pass_error_to_pyerr)?;
        Ok(PyOutcome::new(outcome))
    }

    /// Run `n` deterministic simulations using seeds `seed + i`.
    #[pyo3(signature = (n, seed, *, strict=None))]
    fn run_batch(&self, n: usize, seed: u64, strict: Option<bool>) -> PyResult<Vec<PyOutcome>> {
        let policy = policy_from_strict_override(self.inner.policy(), strict);
        let mut outcomes = Vec::with_capacity(n);
        for result in self.inner.run_batch_with_policy(n, seed, policy) {
            outcomes.push(PyOutcome::new(result.map_err(pass_error_to_pyerr)?));
        }
        Ok(outcomes)
    }

    /// Run one simulation starting from `initial` (clonal expansion).
    /// The plan executes against the supplied parent IR rather than
    /// `Simulation::new()`, letting Python orchestrate "fork from a
    /// parent" semantics for clonal-family generation.
    #[pyo3(signature = (initial, seed, *, strict=None))]
    fn run_from(
        &self,
        initial: &PySimulation,
        seed: u64,
        strict: Option<bool>,
    ) -> PyResult<PyOutcome> {
        let policy = policy_from_strict_override(self.inner.policy(), strict);
        let outcome = self
            .inner
            .run_one_from_with_policy(initial.inner.clone(), seed, policy)
            .map_err(pass_error_to_pyerr)?;
        Ok(PyOutcome::new(outcome))
    }

    fn __repr__(&self) -> String {
        format!(
            "<CompiledSimulator policy={} passes={} contracts={}>",
            self.policy(),
            self.inner.plan().len(),
            self.inner
                .contracts()
                .map(|contracts| contracts.len())
                .unwrap_or(0)
        )
    }
}
