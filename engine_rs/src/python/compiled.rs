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
use crate::trace_file::{
    pass_plan_signature, pass_plan_signature_names_only, refdata_signature, TraceFile,
};

use super::contract::pass_error_to_pyerr;
use super::outcome::PyOutcome;
use super::plan::PyPassPlan;
use super::simulation::PySimulation;
use super::trace_file::PyTraceFile;

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
    /// **Strict semantics apply only to fresh sampling.** Trace replay
    /// via [`Self::replay_from_trace_file`] consumes recorded values
    /// verbatim and does not re-validate them against the contract
    /// bundle — even with `strict=True`. See that method's docstring
    /// for details.
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
        let refdata = self.inner.refdata().ok_or_else(|| {
            PyValueError::new_err(
                "trace_file_from: simulator has no refdata; trace files require both \
                 plan and refdata signatures to round-trip",
            )
        })?;
        let tf = TraceFile::build(
            self.inner.plan(),
            refdata,
            seed,
            outcome.inner.trace.clone(),
        );
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
    /// **The `strict` flag has limited effect on replay.** Replay
    /// consumes each recorded value as a proposal at its sampling
    /// slot and validates *address consistency* and *value-kind
    /// consistency* with the live plan — but does NOT re-run the
    /// contract bundle's admissibility check against the recorded
    /// value. As a consequence:
    ///
    /// - A trace recorded by a permissive run that hit empty support
    ///   and wrote a sentinel value (indel `site = -1`, NP length `0`,
    ///   NP base `N`, trim `0`) **replays cleanly under `strict=True`**
    ///   — the sentinel is consumed verbatim, the original outcome
    ///   reproduces, and no `StrictSamplingError` fires. This is by
    ///   design: replay's contract is "rebuild this exact outcome,"
    ///   not "re-execute the sampler with these inputs."
    /// - To get strict-mode-on-fresh-sampling semantics, call
    ///   [`Self::run`] with `strict=True` and the trace's original
    ///   seed instead.
    ///
    /// See `docs/productive_failure_mode_audit.md` §5 / §6.2 for the
    /// failure matrix and the rationale behind this divergence.
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
        // Signature checks reuse the same flow as rerun_from_trace_file.
        // v1/v2 fixtures recorded the legacy names-only signature; compare
        // like-for-like by rebuilding the live plan in the matching shape.
        // v3+ fixtures embed each pass's compile-time parameter digest, so
        // the live plan is rendered through `pass_plan_signature`.
        let live_plan_sig = if trace_file.inner.schema_version < 3 {
            pass_plan_signature_names_only(self.inner.plan())
        } else {
            pass_plan_signature(self.inner.plan())
        };
        if live_plan_sig != trace_file.inner.pass_plan_signature {
            return Err(PyValueError::new_err(format!(
                "replay_from_trace_file: pass plan signature mismatch.\n  \
                 expected: {}\n  got:      {}",
                trace_file.inner.pass_plan_signature, live_plan_sig,
            )));
        }
        let refdata = self.inner.refdata().ok_or_else(|| {
            PyValueError::new_err(
                "replay_from_trace_file: simulator has no refdata; cannot \
                 verify refdata signature against the trace file",
            )
        })?;
        let live_refdata_sig = refdata_signature(refdata);
        if live_refdata_sig != trace_file.inner.refdata_signature {
            return Err(PyValueError::new_err(format!(
                "replay_from_trace_file: refdata signature mismatch.\n  \
                 expected: {}\n  got:      {}",
                trace_file.inner.refdata_signature, live_refdata_sig,
            )));
        }
        // Content-hash gate (v2+ traces). Catches cartridge differences
        // the structural signature can't see — e.g. identical pool
        // sizes but different rules, identity, or curation provenance.
        // v1 traces (which lack the field) are accepted without this
        // check; the signature alone is the contract there.
        if let Some(recorded_hash) = trace_file.inner.refdata_content_hash.as_ref() {
            let live_hash =
                crate::trace_file::refdata_content_hash(refdata);
            if live_hash != *recorded_hash {
                return Err(PyValueError::new_err(format!(
                    "replay_from_trace_file: refdata content hash mismatch.\n  \
                     expected: {}\n  got:      {}\n  \
                     hint: trace was produced against a different cartridge \
                     (rules / identity / curation may differ).",
                    recorded_hash, live_hash,
                )));
            }
        }

        let policy = policy_from_strict_override(self.inner.policy(), strict);
        let outcome = self
            .inner
            .replay_from_trace_records(
                trace_file.inner.trace.choices(),
                trace_file.inner.seed,
                policy,
            )
            .map_err(pass_error_to_pyerr)?;
        Ok(PyOutcome::new(outcome))
    }

    /// Rerun the engine from `trace_file.seed` against this
    /// simulator's plan + refdata, after verifying both signatures
    /// match the trace file's recorded signatures. The replay
    /// remains **seed-based** — engine determinism reproduces the
    /// original `Outcome`. True trace-injected replay (consuming
    /// recorded values instead of redrawing from the RNG) is a
    /// follow-up; the on-disk schema already supports it.
    ///
    /// **`strict` semantics:** because this method re-runs the
    /// sampler (rather than consuming recorded values), `strict=True`
    /// here behaves the same as a fresh
    /// [`Self::run`] with `strict=True` — it will raise
    /// [`StrictSamplingError`](crate::python::contract) if the seed
    /// produces an empty-support situation under the active
    /// contracts. This is the key difference from
    /// [`Self::replay_from_trace_file`], which consumes recorded
    /// values and does not re-evaluate contract admissibility.
    ///
    /// Raises `ValueError` on signature mismatch (plan changed,
    /// refdata changed, simulator has no refdata, …).
    #[pyo3(signature = (trace_file, *, strict=None))]
    fn rerun_from_trace_file(
        &self,
        trace_file: &PyTraceFile,
        strict: Option<bool>,
    ) -> PyResult<PyOutcome> {
        // Signature checks: plan always, refdata when present.
        // Same v1/v2 vs v3 comparator-shape discipline as
        // `replay_from_trace_file` — see that method for the
        // schema-version-aware rationale.
        let live_plan_sig = if trace_file.inner.schema_version < 3 {
            pass_plan_signature_names_only(self.inner.plan())
        } else {
            pass_plan_signature(self.inner.plan())
        };
        if live_plan_sig != trace_file.inner.pass_plan_signature {
            return Err(PyValueError::new_err(format!(
                "rerun_from_trace_file: pass plan signature mismatch.\n  \
                 expected: {}\n  got:      {}",
                trace_file.inner.pass_plan_signature, live_plan_sig,
            )));
        }
        let refdata = self.inner.refdata().ok_or_else(|| {
            PyValueError::new_err(
                "rerun_from_trace_file: simulator has no refdata; cannot \
                 verify refdata signature against the trace file",
            )
        })?;
        let live_refdata_sig = refdata_signature(refdata);
        if live_refdata_sig != trace_file.inner.refdata_signature {
            return Err(PyValueError::new_err(format!(
                "rerun_from_trace_file: refdata signature mismatch.\n  \
                 expected: {}\n  got:      {}",
                trace_file.inner.refdata_signature, live_refdata_sig,
            )));
        }

        let policy = policy_from_strict_override(self.inner.policy(), strict);
        let outcome = self
            .inner
            .run_one_with_policy(trace_file.inner.seed, policy)
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
