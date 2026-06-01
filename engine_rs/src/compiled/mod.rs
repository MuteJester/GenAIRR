//! Compiled simulator boundary.
//!
//! `PassPlan` is the engine's intermediate representation: an ordered
//! sequence of concrete simulation steps. It should not be the final
//! execution artifact. A `CompiledSimulator` binds that IR to the
//! reference data, active contracts, execution policy, and compile-time
//! analysis report needed to run safely and introspectably.

use crate::address;
use crate::assignment::TrimEnd;
use crate::contract::ContractSet;
use crate::feasibility::FeasibilityContext;
use crate::ir::{Segment, Simulation};
use crate::live_call::{LiveCallRefreshHook, ReferenceMatchIndex};
use crate::pass::{EffectHook, NodeId, Outcome, PassError, PassPlan};
use crate::refdata::{RefDataConfig, RefDataValidationMode};

// Test-only re-imports — the `compiled::tests` submodules pull these
// names through `use super::*`. Gated under `cfg(test)` so non-test
// builds don't carry unused imports.
#[cfg(test)]
use crate::event::TraceSpan;
#[cfg(test)]
use crate::pass::{PassContext, PassEffect};
#[cfg(test)]
use crate::refdata::AlleleId;

/// Compile-time options. `Default` is the production-safe set
/// (strict refdata validation enforced).
///
/// Three knob settings, in order of strictness:
///
/// 1. `Default::default()` — strict refdata validation: every
///    structural problem AND every pseudogene-shaped curatable
///    issue blocks compile.
/// 2. [`CompileOptions::allow_curatable_refdata`] — strict on
///    Fatal issues, lenient on Curatable ones. Use when sampling
///    from a real catalogue (bundled mouse_igh, human_tcrb) that
///    includes pseudogene/ORF alleles you want to keep in the pool.
/// 3. [`CompileOptions::skip_refdata_validation`] — test-only
///    escape hatch; bypasses the validator entirely for synthetic
///    fixtures that intentionally build partial refdata.
#[derive(Copy, Clone, Debug, Eq, PartialEq)]
pub struct CompileOptions {
    /// Whether to run [`RefDataConfig::validate`] at the start of
    /// compile. Always `true` in production paths; `false` only for
    /// unit-test fixtures.
    pub validate_refdata: bool,
    /// Mode for refdata validation when it runs. See
    /// [`RefDataValidationMode`]. Default `Strict`.
    pub refdata_validation_mode: RefDataValidationMode,
}

impl Default for CompileOptions {
    fn default() -> Self {
        Self {
            validate_refdata: true,
            refdata_validation_mode: RefDataValidationMode::Strict,
        }
    }
}

impl CompileOptions {
    /// Production opt-in for catalogues that include pseudogene/ORF
    /// alleles. Curatable-severity issues pass; Fatal issues still
    /// reject. Intended for users explicitly sampling from a raw
    /// IMGT-style catalogue without first filtering it down to
    /// functional alleles.
    pub fn allow_curatable_refdata() -> Self {
        Self {
            validate_refdata: true,
            refdata_validation_mode: RefDataValidationMode::AllowCuratable,
        }
    }

    /// Test-only escape hatch: skip refdata structural validation
    /// entirely. Use ONLY for unit-test fixtures that are
    /// deliberately partial (anchorless synthetic V alleles,
    /// single-segment pools, truncated bytes designed to exercise a
    /// specific engine path). Production callers should never use
    /// this — pseudogene-tolerant simulations should set
    /// `refdata_validation_mode = AllowCuratable` instead, so Fatal
    /// issues still gate.
    #[allow(dead_code)]
    pub(crate) fn skip_refdata_validation() -> Self {
        Self {
            validate_refdata: false,
            refdata_validation_mode: RefDataValidationMode::Strict,
        }
    }
}

/// Runtime failure policy selected at compile time.
#[derive(Copy, Clone, Debug, Eq, PartialEq)]
pub enum ExecutionPolicy {
    /// Preserve permissive runtime behavior: constrained draws that
    /// have no admissible support use the pass's explicit no-op or
    /// sentinel behavior instead of surfacing a structured error.
    /// The recombination-stage allele sampler is the documented
    /// exception: it may fall back to the natural allele draw because
    /// downstream assembly requires an assigned allele.
    Permissive,
    /// Fail loudly with structured errors when a pass cannot execute
    /// safely or cannot sample an admissible candidate.
    Strict,
}

mod report;
pub use report::{CompileReport, CompileWarning, DeclaredChoice, PassSummary};

mod error;
pub use error::{CompileError, CompileErrorKind, CompileErrors};

pub(crate) mod execute;
use execute::{execute_transactional, ExecutionAbort, ExecutionInputs};

/// A borrowed compiled simulation artifact ready for deterministic execution.
///
/// This is the low-level boundary used when the caller owns the
/// `PassPlan`, reference data, and contracts elsewhere. Use
/// [`OwnedCompiledSimulator`] when the compiled artifact itself must
/// outlive the builder inputs, such as Python compile-once/run-many
/// usage.
pub struct CompiledSimulator<'a> {
    plan: &'a PassPlan,
    /// Topologically-sorted execution order produced by
    /// `Schedule::compile`. `execute_transactional` iterates this
    /// instead of `plan.passes()` so the canonical pipeline order is
    /// derived from the dep graph, not the user's push order.
    sorted_order: Vec<NodeId>,
    /// Effect hooks invoked after each pass commits. Today the only
    /// registered hook is the live-call refresh; new derived-state
    /// consumers register here without touching the executor.
    hooks: Vec<Box<dyn EffectHook>>,
    refdata: Option<&'a RefDataConfig>,
    contracts: Option<&'a ContractSet>,
    policy: ExecutionPolicy,
    report: CompileReport,
    feasibility: Option<FeasibilityContext>,
    reference_index: Option<ReferenceMatchIndex>,
}

/// An owning compiled simulation artifact.
///
/// This is the stronger compiled simulator boundary: `PassPlan` is
/// treated as an intermediate representation consumed at compile time,
/// while the resulting simulator owns every input needed for repeated
/// deterministic execution.
pub struct OwnedCompiledSimulator {
    plan: PassPlan,
    /// Topologically-sorted execution order from `Schedule::compile`.
    sorted_order: Vec<NodeId>,
    /// Effect hooks invoked after each pass commits.
    hooks: Vec<Box<dyn EffectHook>>,
    refdata: Option<RefDataConfig>,
    contracts: Option<ContractSet>,
    policy: ExecutionPolicy,
    report: CompileReport,
    feasibility: Option<FeasibilityContext>,
    reference_index: Option<ReferenceMatchIndex>,
}

impl<'a> CompiledSimulator<'a> {
    /// Compile a pass plan into an executable simulator.
    ///
    /// Validation is intentionally semantic and typed: it consumes
    /// `Pass::requirements()` / `Pass::effects()` rather than parsing
    /// pass names. This is the seam where deeper contract precondition
    /// checks and resource validation should be added.
    pub fn compile(
        plan: &'a PassPlan,
        refdata: Option<&'a RefDataConfig>,
        contracts: Option<&'a ContractSet>,
        policy: ExecutionPolicy,
    ) -> Result<Self, CompileErrors> {
        Self::compile_with_options(plan, refdata, contracts, policy, CompileOptions::default())
    }

    /// Compile with explicit options. See [`CompileOptions`]. The
    /// only production-relevant knob is whether refdata validation
    /// runs; the bare [`Self::compile`] enforces it.
    pub fn compile_with_options(
        plan: &'a PassPlan,
        refdata: Option<&'a RefDataConfig>,
        contracts: Option<&'a ContractSet>,
        policy: ExecutionPolicy,
        options: CompileOptions,
    ) -> Result<Self, CompileErrors> {
        // Gate 0: structural refdata validation. Done BEFORE the
        // schedule compile (and before pass precondition checks) so
        // a malformed reference universe surfaces as one aggregated
        // diagnostic, not as a cascade of downstream missing-allele
        // or precondition errors masking the real problem.
        if options.validate_refdata {
            if let Some(rd) = refdata {
                if let Err(errs) = rd.validate_with_mode(options.refdata_validation_mode) {
                    return Err(CompileErrors::from_refdata_validation(errs));
                }
            }
        }
        let sorted_order = match plan.compile(refdata.is_some()) {
            Ok(order) => order,
            Err(schedule_err) => {
                let err = error::schedule_error_into_compile_error(schedule_err, |idx| {
                    plan.passes()[idx].name().to_string()
                });
                return Err(CompileErrors { errors: vec![err] });
            }
        };
        let (mut report, errors, feasibility) = analyze_plan(plan, refdata, contracts, policy);
        if !errors.is_empty() {
            return Err(CompileErrors { errors });
        }
        report.reorder_to_execution(&sorted_order);
        let reference_index = refdata.map(ReferenceMatchIndex::build);

        Ok(Self {
            plan,
            sorted_order,
            hooks: default_hooks(),
            refdata,
            contracts,
            policy,
            report,
            feasibility,
            reference_index,
        })
    }

    pub fn report(&self) -> &CompileReport {
        &self.report
    }

    pub fn policy(&self) -> ExecutionPolicy {
        self.policy
    }

    pub fn plan(&self) -> &PassPlan {
        self.plan
    }

    pub fn refdata(&self) -> Option<&RefDataConfig> {
        self.refdata
    }

    pub fn contracts(&self) -> Option<&ContractSet> {
        self.contracts
    }

    pub(crate) fn feasibility(&self) -> Option<&FeasibilityContext> {
        self.feasibility.as_ref()
    }

    pub fn reference_index(&self) -> Option<&ReferenceMatchIndex> {
        self.reference_index.as_ref()
    }

    fn inputs(&self) -> ExecutionInputs<'_> {
        self.inputs_with_policy(self.policy)
    }

    fn inputs_with_policy(&self, policy: ExecutionPolicy) -> ExecutionInputs<'_> {
        ExecutionInputs {
            plan: self.plan,
            sorted_order: &self.sorted_order,
            hooks: &self.hooks,
            refdata: self.refdata,
            contracts: self.contracts,
            feasibility: self.feasibility.as_ref(),
            reference_index: self.reference_index.as_ref(),
            policy,
            replay_records: None,
        }
    }

    /// Run one simulation from an empty initial IR.
    pub fn run_one(&self, seed: u64) -> Result<Outcome, PassError> {
        self.run_one_from(Simulation::new(), seed)
    }

    /// Run one simulation with an execution policy override.
    ///
    /// The compiled artifact still owns the analyzed plan, reference
    /// data, and active contracts. Only failure semantics differ.
    pub fn run_one_with_policy(
        &self,
        seed: u64,
        policy: ExecutionPolicy,
    ) -> Result<Outcome, PassError> {
        self.run_one_from_with_policy(Simulation::new(), seed, policy)
    }

    /// Run one simulation from a caller-supplied initial IR.
    pub fn run_one_from(&self, initial: Simulation, seed: u64) -> Result<Outcome, PassError> {
        self.execute_transactional(initial, seed)
            .map_err(|abort| abort.error)
    }

    /// Run one simulation from a caller-supplied initial IR with an
    /// execution policy override.
    pub fn run_one_from_with_policy(
        &self,
        initial: Simulation,
        seed: u64,
        policy: ExecutionPolicy,
    ) -> Result<Outcome, PassError> {
        execute_transactional(self.inputs_with_policy(policy), initial, seed)
            .map_err(|abort| abort.error)
    }

    /// Run a deterministic seed-stitched batch. The i-th simulation
    /// uses `seed + i`, matching the current Python `CompiledExperiment`
    /// contract. Parallel seed splitting can replace this behind the
    /// same method later.
    pub fn run_batch(&self, n: usize, seed: u64) -> Vec<Result<Outcome, PassError>> {
        (0..n)
            .map(|i| self.run_one(seed.wrapping_add(i as u64)))
            .collect()
    }

    /// Run a deterministic seed-stitched batch with an execution
    /// policy override.
    pub fn run_batch_with_policy(
        &self,
        n: usize,
        seed: u64,
        policy: ExecutionPolicy,
    ) -> Vec<Result<Outcome, PassError>> {
        (0..n)
            .map(|i| self.run_one_with_policy(seed.wrapping_add(i as u64), policy))
            .collect()
    }

    /// Compiled execution owns the step lifecycle so the compiled
    /// artifact, not `PassPlan`, is the authority that decides when
    /// trace deltas and state revisions are safe to commit.
    fn execute_transactional(
        &self,
        initial: Simulation,
        seed: u64,
    ) -> Result<Outcome, ExecutionAbort> {
        execute_transactional(self.inputs(), initial, seed)
    }
}

impl OwnedCompiledSimulator {
    /// Compile and take ownership of a pass plan plus its execution
    /// dependencies.
    pub fn compile(
        plan: PassPlan,
        refdata: Option<RefDataConfig>,
        contracts: Option<ContractSet>,
        policy: ExecutionPolicy,
    ) -> Result<Self, CompileErrors> {
        Self::compile_with_options(plan, refdata, contracts, policy, CompileOptions::default())
    }

    /// Owning compile with explicit options. See [`CompileOptions`].
    pub fn compile_with_options(
        plan: PassPlan,
        refdata: Option<RefDataConfig>,
        contracts: Option<ContractSet>,
        policy: ExecutionPolicy,
        options: CompileOptions,
    ) -> Result<Self, CompileErrors> {
        // Gate 0: structural refdata validation (see borrowed-form
        // counterpart in `CompiledSimulator::compile_with_options`).
        if options.validate_refdata {
            if let Some(rd) = refdata.as_ref() {
                if let Err(errs) = rd.validate_with_mode(options.refdata_validation_mode) {
                    return Err(CompileErrors::from_refdata_validation(errs));
                }
            }
        }
        let sorted_order = match plan.compile(refdata.is_some()) {
            Ok(order) => order,
            Err(schedule_err) => {
                let err = error::schedule_error_into_compile_error(schedule_err, |idx| {
                    plan.passes()[idx].name().to_string()
                });
                return Err(CompileErrors { errors: vec![err] });
            }
        };
        let (mut report, errors, feasibility) =
            analyze_plan(&plan, refdata.as_ref(), contracts.as_ref(), policy);
        if !errors.is_empty() {
            return Err(CompileErrors { errors });
        }
        report.reorder_to_execution(&sorted_order);
        let reference_index = refdata.as_ref().map(ReferenceMatchIndex::build);

        Ok(Self {
            plan,
            sorted_order,
            hooks: default_hooks(),
            refdata,
            contracts,
            policy,
            report,
            feasibility,
            reference_index,
        })
    }

    pub(crate) fn from_validated_parts(
        plan: PassPlan,
        refdata: Option<RefDataConfig>,
        contracts: Option<ContractSet>,
        policy: ExecutionPolicy,
        report: CompileReport,
        feasibility: Option<FeasibilityContext>,
    ) -> Self {
        // Callers reach here only after the plan already passed
        // `Schedule::compile` upstream (otherwise the validated-parts
        // promise is broken). Recompute the sorted order here rather
        // than threading it through the API.
        let sorted_order = plan
            .compile(refdata.is_some())
            .expect("from_validated_parts requires a schedule that already topo-sorted cleanly");
        let reference_index = refdata.as_ref().map(ReferenceMatchIndex::build);
        Self {
            plan,
            sorted_order,
            hooks: default_hooks(),
            refdata,
            contracts,
            policy,
            report,
            feasibility,
            reference_index,
        }
    }

    pub fn report(&self) -> &CompileReport {
        &self.report
    }

    pub fn policy(&self) -> ExecutionPolicy {
        self.policy
    }

    pub fn plan(&self) -> &PassPlan {
        &self.plan
    }

    pub fn refdata(&self) -> Option<&RefDataConfig> {
        self.refdata.as_ref()
    }

    pub fn contracts(&self) -> Option<&ContractSet> {
        self.contracts.as_ref()
    }

    pub fn reference_index(&self) -> Option<&ReferenceMatchIndex> {
        self.reference_index.as_ref()
    }

    fn inputs(&self) -> ExecutionInputs<'_> {
        self.inputs_with_policy(self.policy)
    }

    fn inputs_with_policy(&self, policy: ExecutionPolicy) -> ExecutionInputs<'_> {
        ExecutionInputs {
            plan: &self.plan,
            sorted_order: &self.sorted_order,
            hooks: &self.hooks,
            refdata: self.refdata.as_ref(),
            contracts: self.contracts.as_ref(),
            feasibility: self.feasibility.as_ref(),
            reference_index: self.reference_index.as_ref(),
            policy,
            replay_records: None,
        }
    }

    /// Run one simulation from an empty initial IR.
    pub fn run_one(&self, seed: u64) -> Result<Outcome, PassError> {
        self.run_one_from(Simulation::new(), seed)
    }

    /// Run one simulation with an execution policy override.
    pub fn run_one_with_policy(
        &self,
        seed: u64,
        policy: ExecutionPolicy,
    ) -> Result<Outcome, PassError> {
        self.run_one_from_with_policy(Simulation::new(), seed, policy)
    }

    /// Run one simulation from a caller-supplied initial IR.
    pub fn run_one_from(&self, initial: Simulation, seed: u64) -> Result<Outcome, PassError> {
        self.execute_transactional(initial, seed)
            .map_err(|abort| abort.error)
    }

    /// Run one simulation from a caller-supplied initial IR with an
    /// execution policy override.
    pub fn run_one_from_with_policy(
        &self,
        initial: Simulation,
        seed: u64,
        policy: ExecutionPolicy,
    ) -> Result<Outcome, PassError> {
        execute_transactional(self.inputs_with_policy(policy), initial, seed)
            .map_err(|abort| abort.error)
    }

    /// Trace-injected replay (Option B): run the plan against an
    /// empty initial IR, consuming the supplied `records` instead of
    /// drawing from the RNG at sampling sites that have been
    /// migrated to the trace-injected path.
    ///
    /// The `seed` is still used for unmigrated sites (those still
    /// branch on `ctx.rng`); migrated sites read from the cursor.
    /// `policy` controls the strict / permissive failure mode for
    /// unmigrated sites as usual.
    ///
    /// Returns `PassError::Replay` on cursor-vs-plan disagreement —
    /// address mismatch, value-kind mismatch, exhausted trace, or
    /// trailing unused records.
    pub fn replay_from_trace_records(
        &self,
        records: &[crate::trace::ChoiceRecord],
        seed: u64,
        policy: ExecutionPolicy,
    ) -> Result<Outcome, PassError> {
        let mut inputs = self.inputs_with_policy(policy);
        inputs.replay_records = Some(records);
        execute_transactional(inputs, Simulation::new(), seed).map_err(|abort| abort.error)
    }

    /// Run a deterministic seed-stitched batch using the same policy
    /// as [`CompiledSimulator::run_batch`].
    pub fn run_batch(&self, n: usize, seed: u64) -> Vec<Result<Outcome, PassError>> {
        (0..n)
            .map(|i| self.run_one(seed.wrapping_add(i as u64)))
            .collect()
    }

    /// Run a deterministic seed-stitched batch with an execution
    /// policy override.
    pub fn run_batch_with_policy(
        &self,
        n: usize,
        seed: u64,
        policy: ExecutionPolicy,
    ) -> Vec<Result<Outcome, PassError>> {
        (0..n)
            .map(|i| self.run_one_with_policy(seed.wrapping_add(i as u64), policy))
            .collect()
    }

    fn execute_transactional(
        &self,
        initial: Simulation,
        seed: u64,
    ) -> Result<Outcome, ExecutionAbort> {
        execute_transactional(self.inputs(), initial, seed)
    }
}

mod analyze;
use analyze::analyze_plan;

mod feasibility_builder;

fn sample_allele_address(segment: Segment) -> &'static str {
    match segment {
        Segment::V | Segment::D | Segment::J => address::sample_allele_vdj(segment),
        Segment::Np1 | Segment::Np2 => address::SAMPLE_ALLELE_INVALID,
    }
}

fn trim_address(segment: Segment, end: TrimEnd) -> &'static str {
    match segment {
        Segment::V | Segment::D | Segment::J => address::trim_vdj(segment, end),
        Segment::Np1 | Segment::Np2 => address::TRIM_INVALID,
    }
}

fn np_length_address(segment: Segment) -> &'static str {
    match segment {
        Segment::Np1 | Segment::Np2 => address::np_length_region(segment),
        Segment::V | Segment::D | Segment::J => address::NP_INVALID_LENGTH,
    }
}

/// Default effect-hook registration for a freshly-compiled
/// simulator. Today this is just the live-call refresh; future
/// derived-state consumers (e.g. a junction cache) get added here.
fn default_hooks() -> Vec<Box<dyn EffectHook>> {
    vec![Box::new(LiveCallRefreshHook::new())]
}

#[cfg(test)]
mod tests;
