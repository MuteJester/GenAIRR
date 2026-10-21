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
    family_codec: std::cell::RefCell<Option<std::sync::Arc<FamilyPlanBinding>>>,
    family_allow_curatable: bool,
}

#[derive(Debug)]
pub(crate) struct FamilyPlanBinding {
    pub spec: crate::family_replay::plan::ReplayPlanSpecV1,
    pub json: String,
    pub digest: String,
    pub reference_digest: String,
}

#[derive(serde::Serialize, serde::Deserialize)]
#[serde(deny_unknown_fields)]
struct FamilyTraceReceiptV1 {
    schema_version: u32,
    plan_digest: String,
    reference_digest: String,
    productivity_policy: String,
    strict: bool,
    source: bool,
    seed_label: String,
    initial_state_digest: String,
    final_state_digest: String,
    events_digest: String,
    trace: crate::trace::Trace,
    source_trace: Option<crate::trace_file::TraceFile>,
}

impl PyCompiledSimulator {
    fn prepare_family_codec(&self) -> Result<std::sync::Arc<FamilyPlanBinding>, String> {
        if let Some(binding) = self.family_codec.borrow().as_ref() {
            return Ok(binding.clone());
        }
        let refdata = self
            .inner
            .refdata()
            .ok_or("family replay requires reference data")?;
        let spec = crate::family_replay::plan::ReplayPlanSpecV1::capture(
            self.inner.plan(),
            refdata,
            self.inner.contracts(),
            self.inner.policy() == ExecutionPolicy::Strict,
            self.family_allow_curatable,
        )?;
        let (rebuilt, rebuilt_reference, _) = spec.reconstruct()?;
        if crate::trace_file::refdata_content_hash(&rebuilt_reference)
            != crate::trace_file::refdata_content_hash(refdata)
            || self
                .inner
                .plan()
                .passes()
                .iter()
                .zip(rebuilt.passes())
                .any(|(a, b)| {
                    a.name() != b.name() || a.sampling_law_signature() != b.sampling_law_signature()
                })
        {
            return Err(
                "native descriptor does not preserve the complete execution context".into(),
            );
        }
        let original_order = self
            .inner
            .plan()
            .compile(true)
            .map_err(|e| format!("family source schedule: {e:?}"))?;
        let rebuilt_order = rebuilt
            .compile(true)
            .map_err(|e| format!("family rebuilt schedule: {e:?}"))?;
        if original_order != rebuilt_order {
            return Err("family replay v1 cannot encode an execution schedule changed by explicit ordering edges".into());
        }
        let digest = spec.digest()?;
        let reference_digest = spec.reference.digest()?;
        let json = serde_json::to_string(&spec).map_err(|e| e.to_string())?;
        let binding = std::sync::Arc::new(FamilyPlanBinding {
            spec,
            json,
            digest,
            reference_digest,
        });
        *self.family_codec.borrow_mut() = Some(binding.clone());
        Ok(binding)
    }
    fn family_context(
        &self,
        policy: ExecutionPolicy,
        seed: u64,
    ) -> Option<super::outcome::FamilyExecutionContext> {
        self.family_codec
            .borrow()
            .as_ref()
            .map(|binding| super::outcome::FamilyExecutionContext {
                binding: binding.clone(),
                strict: policy == ExecutionPolicy::Strict,
                seed,
            })
    }
    fn stamp_family(&self, outcome: PyOutcome, policy: ExecutionPolicy, seed: u64) -> PyOutcome {
        outcome.with_family_execution(self.family_context(policy, seed))
    }
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
            family_codec: std::cell::RefCell::new(None),
            family_allow_curatable: allow_curatable_refdata,
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
    /// Prepare a closed executable descriptor before any captured draw.
    #[pyo3(signature=(budget=None))]
    fn family_replay_spec(
        &self,
        budget: Option<&super::family_replay::PyFamilyReplayBudget>,
    ) -> PyResult<String> {
        let prepare = |_: &mut crate::family_replay::FamilyReplayBudget| {
            self.prepare_family_codec().map(|b| b.json.clone())
        };
        if let Some(b) = budget {
            let mut guard = b
                .inner
                .lock()
                .map_err(|_| PyValueError::new_err("budget lock poisoned"))?;
            crate::family_replay::with_budget_scope(&mut guard, prepare)
        } else {
            let mut b = crate::family_replay::FamilyReplayBudget::new(
                crate::family_replay::FamilyReplayLimits::default(),
            );
            crate::family_replay::with_budget_scope(&mut b, prepare)
        }
        .map_err(PyValueError::new_err)
    }

    #[staticmethod]
    #[pyo3(signature=(text,budget=None))]
    fn from_family_replay_spec(
        text: &str,
        budget: Option<&super::family_replay::PyFamilyReplayBudget>,
    ) -> PyResult<Self> {
        use crate::family_replay::{
            plan::ReplayPlanSpecV1, strict_json, FamilyReplayBudget, FamilyReplayLimits,
        };
        let limits = budget
            .map(|b| {
                b.inner
                    .lock()
                    .map(|g| g.limits.clone())
                    .map_err(|_| PyValueError::new_err("budget lock poisoned"))
            })
            .transpose()?
            .unwrap_or_default();
        let spec: ReplayPlanSpecV1 =
            strict_json::from_str(text, &limits).map_err(PyValueError::new_err)?;
        let build = |b: &mut FamilyReplayBudget| -> Result<Self, String> {
            let bases = [
                &spec.reference.v,
                &spec.reference.d,
                &spec.reference.j,
                &spec.reference.c,
            ]
            .into_iter()
            .flatten()
            .try_fold(0u64, |n, a| {
                n.checked_add(a.sequence.len() as u64)
                    .ok_or("reference size overflow")
            })?;
            b.charge_bases(bases)?;
            b.charge_work(spec.passes.len() as u64)?;
            let (plan, reference, contracts) = spec.reconstruct()?;
            let mut builder = PyPassPlan::from_native(plan);
            let compiled = Self::compile_from_plan(
                &mut builder,
                Some(reference),
                contracts,
                if spec.strict {
                    ExecutionPolicy::Strict
                } else {
                    ExecutionPolicy::Permissive
                },
                spec.allow_curatable_refdata,
            )
            .map_err(|e| e.to_string())?;
            let actual = compiled.prepare_family_codec()?;
            if actual.spec != spec {
                return Err("native plan descriptor did not round-trip exactly".into());
            }
            Ok(compiled)
        };
        if let Some(b) = budget {
            let mut guard = b
                .inner
                .lock()
                .map_err(|_| PyValueError::new_err("budget lock poisoned"))?;
            crate::family_replay::with_budget_scope(&mut guard, build)
                .map_err(PyValueError::new_err)
        } else {
            let mut b = FamilyReplayBudget::new(FamilyReplayLimits::default());
            crate::family_replay::with_budget_scope(&mut b, build).map_err(PyValueError::new_err)
        }
    }

    fn family_replay_refdata(&self) -> PyResult<super::refdata::PyRefDataConfig> {
        Ok(super::refdata::PyRefDataConfig {
            inner: self
                .inner
                .refdata()
                .ok_or_else(|| PyValueError::new_err("missing reference data"))?
                .clone(),
        })
    }

    #[pyo3(signature=(seed,strict,budget,initial=None))]
    fn run_for_family_capture(
        &self,
        seed: u64,
        strict: bool,
        budget: &super::family_replay::PyFamilyReplayBudget,
        initial: Option<&PySimulation>,
    ) -> PyResult<PyOutcome> {
        let source = crate::trace_file::duplex_source(self.inner.plan()).is_some();
        if source != initial.is_none()
            || initial.is_some_and(|s| s.inner.duplex_provenance.is_none())
        {
            return Err(PyValueError::new_err(
                "family source/continuation initial state mismatch",
            ));
        }
        let policy = if strict {
            ExecutionPolicy::Strict
        } else {
            ExecutionPolicy::Permissive
        };
        let mut guard = budget
            .inner
            .lock()
            .map_err(|_| PyValueError::new_err("budget lock poisoned"))?;
        let outcome = crate::family_replay::with_budget_scope(&mut guard, |b| {
            self.prepare_family_codec()?;
            let initial = initial.map(|s| s.inner.clone()).unwrap_or_default();
            b.charge_bases(initial.pool.len() as u64)?;
            let out = self
                .inner
                .run_one_from_with_policy(initial, seed, policy)
                .map_err(|e| e.to_string())?;
            b.charge_choices(out.trace.len() as u64)?;
            b.charge_history(out.final_simulation())?;
            Ok(out)
        })
        .map_err(PyValueError::new_err)?;
        let wrapped = PyOutcome::new(outcome).with_duplex_execution(if source {
            self.duplex_execution.clone()
        } else {
            None
        });
        Ok(self.stamp_family(wrapped, policy, seed))
    }

    #[pyo3(signature=(outcome,strict,budget=None))]
    fn family_trace_receipt(
        &self,
        outcome: &PyOutcome,
        strict: bool,
        budget: Option<&super::family_replay::PyFamilyReplayBudget>,
    ) -> PyResult<String> {
        use crate::family_replay::{
            outcome_events_digest, semantic_state_digest, FamilyReplayBudget, FamilyReplayLimits,
        };
        let capture = |b: &mut FamilyReplayBudget| -> Result<String, String> {
            let binding = self.prepare_family_codec()?;
            let actual = outcome
                .family_execution
                .as_ref()
                .ok_or("outcome lacks captured execution context")?;
            if actual.strict != strict
                || actual.binding.digest != binding.digest
                || actual.binding.reference_digest != binding.reference_digest
            {
                return Err("family outcome execution binding/policy mismatch".into());
            }
            let source = crate::trace_file::duplex_source(self.inner.plan()).is_some();
            b.charge_choices(outcome.inner.trace.len() as u64)?;
            b.charge_history(outcome.inner.final_simulation())?;
            let source_trace = if source {
                Some(
                    self.trace_file_from(outcome, actual.seed)
                        .map_err(|e| e.to_string())?
                        .inner,
                )
            } else {
                None
            };
            let initial = outcome
                .inner
                .revisions
                .first()
                .ok_or("outcome has no initial state")?;
            if !source && initial.duplex_provenance.is_none() {
                return Err("family continuation lacks bound duplex initial state".into());
            }
            if outcome.inner.pass_names != self.inner.report().pass_names() {
                return Err("family stage pass layout mismatch".into());
            }
            let receipt = FamilyTraceReceiptV1 {
                schema_version: 1,
                plan_digest: binding.digest.clone(),
                reference_digest: binding.reference_digest.clone(),
                productivity_policy: crate::productivity::DUPLEX_POLICY.into(),
                strict,
                source,
                seed_label: actual.seed.to_string(),
                initial_state_digest: semantic_state_digest(initial)?,
                final_state_digest: semantic_state_digest(outcome.inner.final_simulation())?,
                events_digest: outcome_events_digest(&outcome.inner)?,
                trace: outcome.inner.trace.clone(),
                source_trace,
            };
            serde_json::to_string(&receipt).map_err(|e| e.to_string())
        };
        if let Some(b) = budget {
            let mut guard = b
                .inner
                .lock()
                .map_err(|_| PyValueError::new_err("budget lock poisoned"))?;
            crate::family_replay::with_budget_scope(&mut guard, capture)
        } else {
            let mut b = FamilyReplayBudget::new(FamilyReplayLimits::default());
            crate::family_replay::with_budget_scope(&mut b, capture)
        }
        .map_err(PyValueError::new_err)
    }

    #[pyo3(signature=(initial,receipt_json,budget=None))]
    fn replay_family_trace(
        &self,
        initial: Option<&PySimulation>,
        receipt_json: &str,
        budget: Option<&super::family_replay::PyFamilyReplayBudget>,
    ) -> PyResult<PyOutcome> {
        use crate::family_replay::{
            outcome_events_digest, semantic_state_digest, strict_json, FamilyReplayBudget,
            FamilyReplayLimits,
        };
        let binding = self.prepare_family_codec().map_err(PyValueError::new_err)?;
        let limits = budget
            .map(|b| {
                b.inner
                    .lock()
                    .map(|g| g.limits.clone())
                    .map_err(|_| PyValueError::new_err("budget lock poisoned"))
            })
            .transpose()?
            .unwrap_or_default();
        let value = strict_json::parse(receipt_json, &limits).map_err(PyValueError::new_err)?;
        let receipt: FamilyTraceReceiptV1 = serde_json::from_value(value.clone())
            .map_err(|e| PyValueError::new_err(e.to_string()))?;
        if serde_json::to_value(&receipt).map_err(|e| PyValueError::new_err(e.to_string()))?
            != value
        {
            return Err(PyValueError::new_err("unknown/noncanonical receipt fields"));
        }
        if receipt.schema_version != 1
            || receipt.plan_digest != binding.digest
            || receipt.reference_digest != binding.reference_digest
            || receipt.productivity_policy != crate::productivity::DUPLEX_POLICY
        {
            return Err(PyValueError::new_err("family receipt context mismatch"));
        }
        let source = crate::trace_file::duplex_source(self.inner.plan()).is_some();
        if source != receipt.source
            || source != initial.is_none()
            || (!source && initial.is_some_and(|s| s.inner.duplex_provenance.is_none()))
        {
            return Err(PyValueError::new_err(
                "family receipt initial/source mismatch",
            ));
        }
        if let Some(trace) = &receipt.source_trace {
            if !source {
                return Err(PyValueError::new_err("source trace on continuation"));
            }
            trace
                .validate_for_replay(self.inner.plan(), self.inner.refdata().unwrap())
                .map_err(|e| PyValueError::new_err(e.to_string()))?;
            trace
                .validate_contracts(self.inner.contracts())
                .map_err(|e| PyValueError::new_err(e.to_string()))?;
            if serde_json::to_value(&trace.trace).unwrap()
                != serde_json::to_value(&receipt.trace).unwrap()
            {
                return Err(PyValueError::new_err("source and receipt choices disagree"));
            }
        } else if source {
            return Err(PyValueError::new_err("missing source v4 trace"));
        }
        let seed = receipt
            .seed_label
            .parse::<u64>()
            .map_err(|_| PyValueError::new_err("invalid seed label"))?;
        if seed.to_string() != receipt.seed_label {
            return Err(PyValueError::new_err("noncanonical seed label"));
        }
        if receipt
            .source_trace
            .as_ref()
            .is_some_and(|t| t.seed != seed)
        {
            return Err(PyValueError::new_err("source seed labels disagree"));
        }
        let policy = if receipt.strict {
            ExecutionPolicy::Strict
        } else {
            ExecutionPolicy::Permissive
        };
        let replay = |b: &mut FamilyReplayBudget| {
            let start = initial.map(|s| s.inner.clone()).unwrap_or_default();
            if semantic_state_digest(&start)? != receipt.initial_state_digest {
                return Err("family initial state commitment mismatch".into());
            }
            b.charge_choices(receipt.trace.len() as u64)?;
            b.charge_bases(start.pool.len() as u64)?;
            let out = crate::family_replay::with_replay_scope(|| {
                self.inner
                    .replay_family_from(start.clone(), receipt.trace.choices(), policy)
                    .map_err(|e| e.to_string())
            })?;
            b.charge_history(out.final_simulation())?;
            if semantic_state_digest(out.final_simulation())? != receipt.final_state_digest
                || outcome_events_digest(&out)? != receipt.events_digest
            {
                return Err("family replay state/event commitment mismatch".into());
            }
            Ok(out)
        };
        let outcome = if let Some(b) = budget {
            let mut guard = b
                .inner
                .lock()
                .map_err(|_| PyValueError::new_err("budget lock poisoned"))?;
            crate::family_replay::with_budget_scope(&mut guard, replay)
        } else {
            let mut b = FamilyReplayBudget::new(FamilyReplayLimits::default());
            crate::family_replay::with_budget_scope(&mut b, replay)
        }
        .map_err(PyValueError::new_err)?;
        Ok(self.stamp_family(
            PyOutcome::new(outcome).with_duplex_execution(if source {
                self.duplex_execution.clone()
            } else {
                None
            }),
            policy,
            seed,
        ))
    }
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
        Ok(self.stamp_family(
            PyOutcome::new(outcome).with_duplex_execution(self.duplex_execution.clone()),
            policy,
            seed,
        ))
    }

    /// Run `n` deterministic simulations using seeds `seed + i`.
    #[pyo3(signature = (n, seed, *, strict=None))]
    fn run_batch(&self, n: usize, seed: u64, strict: Option<bool>) -> PyResult<Vec<PyOutcome>> {
        let policy = policy_from_strict_override(self.inner.policy(), strict);
        let mut outcomes = Vec::with_capacity(n);
        for result in self.inner.run_batch_with_policy(n, seed, policy) {
            outcomes.push(
                PyOutcome::new(result.map_err(pass_error_to_pyerr)?)
                    .with_duplex_execution(self.duplex_execution.clone())
                    .with_family_execution(
                        self.family_context(policy, seed.wrapping_add(outcomes.len() as u64)),
                    ),
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
        Ok(self.stamp_family(PyOutcome::new(outcome), policy, seed))
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
        Ok(self.stamp_family(
            PyOutcome::new(outcome).with_duplex_execution(self.duplex_execution.clone()),
            policy,
            trace_file.inner.seed,
        ))
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
        Ok(self.stamp_family(
            PyOutcome::new(outcome).with_duplex_execution(self.duplex_execution.clone()),
            policy,
            trace_file.inner.seed,
        ))
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
