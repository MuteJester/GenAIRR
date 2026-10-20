//! Python-facing owning compiled simulator.
//!
//! `PassPlan` is a builder/IR. `CompiledSimulator` is the durable
//! execution artifact: it owns the plan plus cloned reference data and
//! contracts, so repeated `run()` calls do not recompile or reborrow
//! from Python-side builder objects.

use pyo3::exceptions::PyValueError;
use pyo3::prelude::*;

use crate::compiled::{CompileOptions, CompiledSimulator, ExecutionPolicy, OwnedCompiledSimulator};
use crate::contract::ContractSet;
use crate::refdata::RefDataConfig;
use crate::trace_file::TraceFile;

use super::contract::pass_error_to_pyerr;
use super::outcome::PyOutcome;
use super::plan::PyPassPlan;
use super::simulation::PySimulation;
use super::trace_file::PyTraceFile;

#[pyclass(name = "CompiledSimulator", module = "GenAIRR._engine", unsendable)]
pub struct PyCompiledSimulator {
    inner: OwnedCompiledSimulator,
    duplex_execution: Option<std::sync::Arc<super::outcome::DuplexExecutionBinding>>,
}

impl PyCompiledSimulator {
    pub(crate) fn compile_from_plan(
        plan_builder: &mut PyPassPlan,
        refdata: Option<RefDataConfig>,
        contracts: Option<ContractSet>,
        policy: ExecutionPolicy,
        allow_curatable_refdata: bool,
    ) -> PyResult<Self> {
        let options = if allow_curatable_refdata {
            CompileOptions::allow_curatable_refdata()
        } else {
            CompileOptions::default()
        };
        let (report, feasibility) = {
            let plan = plan_builder.inner()?;
            let borrowed = CompiledSimulator::compile_with_options(
                plan,
                refdata.as_ref(),
                contracts.as_ref(),
                policy,
                options,
            )
            .map_err(|err| PyValueError::new_err(err.to_string()))?;
            (borrowed.report().clone(), borrowed.feasibility().cloned())
        };

        let plan = plan_builder.take_inner()?;
        let duplex_execution = crate::trace_file::duplex_source(&plan).map(|_| {
            let mut active_contracts: Vec<_> = contracts
                .as_ref()
                .map(|set| set.iter().map(|c| c.name().to_owned()).collect())
                .unwrap_or_default();
            active_contracts.sort();
            std::sync::Arc::new(super::outcome::DuplexExecutionBinding {
                plan_signature: crate::trace_file::pass_plan_signature(&plan),
                reference_hash: crate::trace_file::refdata_content_hash(
                    refdata
                        .as_ref()
                        .expect("duplex source compilation requires references"),
                ),
                active_contracts,
            })
        });
        Ok(Self {
            duplex_execution,
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

    /// Run one **fresh** simulation using this compiled artifact.
    ///
    /// The optional `strict` flag overrides the policy compiled into
    /// this artifact for this call only.
    ///
    /// **Strict mode (`strict=True`):** raises
    /// [`StrictSamplingError`](crate::python::contract) when a pass's
    /// contract-narrowed candidate set becomes empty at sample time.
    /// The exception's `args` tuple is `(pass_name, address, reason)`.
    ///
    /// **Permissive mode (`strict=False`, default):** when a sampler's
    /// admissible support is empty, the pass falls back to its
    /// declared `EmptySupport` policy and records a documented sentinel
    /// value to the trace (e.g. indel `site = -1` NoOp, NP length `0`,
    /// NP base `N`, trim `0`) or skips the slot.
    ///
    /// Productive initial assignment and P-aware recombination raise on
    /// empty support even in permissive mode. Replay consumes recorded
    /// proposals without replacement draws; migrated passes validate their
    /// admissibility. See [`Self::replay_from_trace_file`] for details.
    #[pyo3(signature = (seed, *, strict=None))]
    fn run(&self, seed: u64, strict: Option<bool>) -> PyResult<PyOutcome> {
        let policy = policy_from_strict_override(self.inner.policy(), strict);
        let outcome = self
            .inner
            .run_one_with_policy(seed, policy)
            .map_err(pass_error_to_pyerr)?;
        Ok(PyOutcome::new(outcome).with_duplex_execution(self.duplex_execution.clone()))
    }

    /// Run `n` deterministic simulations using seeds `seed + i`.
    #[pyo3(signature = (n, seed, *, strict=None))]
    fn run_batch(&self, n: usize, seed: u64, strict: Option<bool>) -> PyResult<Vec<PyOutcome>> {
        let policy = policy_from_strict_override(self.inner.policy(), strict);
        let mut outcomes = Vec::with_capacity(n);
        for result in self.inner.run_batch_with_policy(n, seed, policy) {
            outcomes.push(
                PyOutcome::new(result.map_err(pass_error_to_pyerr)?)
                    .with_duplex_execution(self.duplex_execution.clone()),
            );
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

    /// Bundle the given `outcome`'s trace into a durable
    /// [`PyTraceFile`] paired with this simulator's plan and refdata
    /// signatures plus the `seed` it was produced from.
    ///
    /// The resulting file can be serialised to JSON via
    /// `trace_file.write_to(path)` and reloaded later with
    /// `TraceFile.read_from(path)` + this simulator's
    /// [`Self::rerun_from_trace_file`].
    ///
    /// Raises `ValueError` when the simulator has no refdata
    /// attached (trace files require both signatures to round-trip
    /// safely).
    fn trace_file_from(&self, outcome: &PyOutcome, seed: u64) -> PyResult<PyTraceFile> {
        if self.duplex_execution.is_some() || super::outcome::has_duplex_provenance(&outcome.inner)
        {
            if self.duplex_execution.is_none() || outcome.duplex_execution != self.duplex_execution
            {
                return Err(PyValueError::new_err(
                    "duplex trace source binding mismatch: an initial repair snapshot or an Outcome generated by this exact compiled sampling law is required",
                ));
            }
            let source = crate::trace_file::duplex_source(self.inner.plan());
            let generated_here = outcome.inner.revisions.first().is_some_and(|s| {
                s.duplex_provenance.is_none()
                    && s.pool.is_empty()
                    && s.sequence.regions.is_empty()
                    && [
                        crate::ir::Segment::V,
                        crate::ir::Segment::D,
                        crate::ir::Segment::J,
                    ]
                    .into_iter()
                    .all(|segment| s.assignments.get(segment).is_none())
                    && s.mutation_count == 0
                    && s.mutation_burden == 0.0
            }) && outcome.inner.pass_names == self.inner.report().pass_names()
                && source.is_some_and(|pass| {
                    outcome
                        .inner
                        .final_simulation()
                        .duplex_provenance
                        .as_ref()
                        .and_then(|p| p.founder.selection.as_ref())
                        .is_some_and(|s| s.model_digest == pass.model().digest())
                });
            if !generated_here {
                return Err(PyValueError::new_err(
                    "duplex continuation trace files are unsupported: the legacy trace schema does not bind the initial repair snapshot; use a finite duplex source plan for durable replay",
                ));
            }
        }
        let refdata = self.inner.refdata().ok_or_else(|| {
            PyValueError::new_err(
                "trace_file_from: simulator has no refdata; trace files require both \
                 plan and refdata signatures to round-trip",
            )
        })?;
        let mut tf = TraceFile::build(
            self.inner.plan(),
            refdata,
            seed,
            outcome.inner.trace.clone(),
        );
        tf.active_contracts = Some(self.active_contracts());
        Ok(PyTraceFile::new(tf))
    }

    /// Trace-injected replay (Option B): run the plan against an
    /// empty initial IR, consuming the recorded values in
    /// `trace_file.trace` at every sampling site that has been
    /// migrated to the consume-trace path. Sites that have not yet
    /// been migrated still draw from the RNG seeded at
    /// `trace_file.seed`.
    ///
    /// This is the architectural sibling of
    /// [`Self::rerun_from_trace_file`]: same input bundle, same
    /// signature checks, but the cursor — not the RNG — is the
    /// source of truth for migrated sites. As more passes migrate
    /// (per the staged plan: allele → trim/NP-length → SHM count
    /// loops → structural), the fraction of values that come from
    /// the cursor grows; once every site is migrated, the engine is
    /// **deterministic in trace**, independent of RNG implementation.
    ///
    /// Recorded choices are proposals: each migrated pass checks address,
    /// value kind and its supported replay rules. P-aware productive
    /// recombination additionally requires positive natural support and a
    /// productive completion for P/N lengths, N bases and D orientation,
    /// regardless of `strict`. Invalid proposals raise without resampling.
    /// Existing no-op/sentinel traces remain valid only where the relevant
    /// pass's replay policy permits that sentinel. To rerun fresh sampling
    /// from the original seed instead, use [`Self::rerun_from_trace_file`].
    ///
    /// Raises `ValueError` on signature mismatch *or* on a
    /// cursor-vs-plan disagreement raised by a migrated pass
    /// (address mismatch, value-kind mismatch, exhausted trace,
    /// trailing unused records).
    #[pyo3(signature = (trace_file, *, strict=None))]
    fn replay_from_trace_file(
        &self,
        trace_file: &PyTraceFile,
        strict: Option<bool>,
    ) -> PyResult<PyOutcome> {
        let refdata = self.inner.refdata().ok_or_else(|| {
            PyValueError::new_err("replay_from_trace_file: simulator has no refdata")
        })?;
        trace_file
            .inner
            .validate_for_replay(self.inner.plan(), refdata)
            .map_err(|e| PyValueError::new_err(format!("replay_from_trace_file: {e}")))?;
        trace_file
            .inner
            .validate_contracts(self.inner.contracts())
            .map_err(|e| PyValueError::new_err(e.to_string()))?;

        let policy = policy_from_strict_override(self.inner.policy(), strict);
        let outcome = self
            .inner
            .replay_from_trace_records(
                trace_file.inner.trace.choices(),
                trace_file.inner.seed,
                policy,
            )
            .map_err(pass_error_to_pyerr)?;
        Ok(PyOutcome::new(outcome).with_duplex_execution(self.duplex_execution.clone()))
    }

    /// Rerun the engine from `trace_file.seed` against this
    /// simulator's plan + refdata, after verifying both signatures
    /// match the trace file's recorded signatures. The replay
    /// remains **seed-based** — engine determinism reproduces the
    /// original `Outcome`. To consume recorded choices instead of
    /// redrawing from the RNG, use [`Self::replay_from_trace_file`].
    ///
    /// **`strict` semantics:** because this method re-runs the
    /// sampler (rather than consuming recorded values), `strict=True`
    /// here behaves the same as a fresh
    /// [`Self::run`] with `strict=True` — it will raise
    /// [`StrictSamplingError`](crate::python::contract) if the seed
    /// produces an empty-support situation under the active
    /// contracts. This is the key difference from
    /// [`Self::replay_from_trace_file`], which consumes recorded
    /// proposals and validates each pass's applicable replay rules.
    ///
    /// Raises `ValueError` on signature mismatch (plan changed,
    /// refdata changed, simulator has no refdata, …).
    #[pyo3(signature = (trace_file, *, strict=None))]
    fn rerun_from_trace_file(
        &self,
        trace_file: &PyTraceFile,
        strict: Option<bool>,
    ) -> PyResult<PyOutcome> {
        let refdata = self.inner.refdata().ok_or_else(|| {
            PyValueError::new_err("rerun_from_trace_file: simulator has no refdata")
        })?;
        trace_file
            .inner
            .validate_against(self.inner.plan(), refdata)
            .map_err(|e| PyValueError::new_err(format!("rerun_from_trace_file: {e}")))?;
        trace_file
            .inner
            .validate_contracts(self.inner.contracts())
            .map_err(|e| PyValueError::new_err(e.to_string()))?;

        let policy = policy_from_strict_override(self.inner.policy(), strict);
        let outcome = self
            .inner
            .run_one_with_policy(trace_file.inner.seed, policy)
            .map_err(pass_error_to_pyerr)?;
        Ok(PyOutcome::new(outcome).with_duplex_execution(self.duplex_execution.clone()))
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
