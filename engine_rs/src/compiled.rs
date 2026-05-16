//! Compiled simulator boundary.
//!
//! `PassPlan` is the engine's intermediate representation: an ordered
//! sequence of concrete simulation steps. It should not be the final
//! execution artifact. A `CompiledSimulator` binds that IR to the
//! reference data, active contracts, execution policy, and compile-time
//! analysis report needed to run safely and introspectably.

use crate::assignment::TrimEnd;
use crate::contract::{ContractKind, ContractSet};
use crate::event::{EventRecord, StateSummary, TraceSpan};
use crate::feasibility::{ChoiceDomain, FeasibilityContext, VjProductiveFeasibility};
use crate::ir::{Segment, Simulation};
use crate::live_call::{with_assembled_segment_live_call, ReferenceMatchIndex};
use crate::pass::{
    AlleleIdSupport, IntegerSupport, Outcome, PassCompileFact, PassContext, PassEffect, PassError,
    PassPlan, PassRequirement,
};
use crate::refdata::{AlleleId, RefDataConfig};
use crate::rng::Rng;
use crate::trace::Trace;
use std::collections::{HashMap, HashSet};
use std::fmt;

/// Runtime failure policy selected at compile time.
#[derive(Copy, Clone, Debug, Eq, PartialEq)]
pub enum ExecutionPolicy {
    /// Preserve permissive runtime behavior: constrained draws may
    /// fall back to unconstrained sampling when no support is available.
    Permissive,
    /// Fail loudly with structured errors when a pass cannot execute
    /// safely or cannot sample an admissible candidate.
    Strict,
}

/// One pass-level compile summary.
#[derive(Clone, Debug, Eq, PartialEq)]
pub struct PassSummary {
    pub index: usize,
    pub name: String,
    pub declared_choices: Vec<String>,
    pub requirements: Vec<PassRequirement>,
    pub effects: Vec<PassEffect>,
}

/// One declared stochastic choice, annotated with the pass that owns it.
#[derive(Clone, Debug, Eq, PartialEq)]
pub struct DeclaredChoice {
    pub pass_index: usize,
    pub pass_name: String,
    pub address: String,
}

/// Compile-time analysis output for a plan.
#[derive(Clone, Debug, Eq, PartialEq)]
pub struct CompileReport {
    pub policy: ExecutionPolicy,
    pub pass_summaries: Vec<PassSummary>,
    pub declared_choices: Vec<DeclaredChoice>,
    pub active_contracts: Vec<String>,
    pub warnings: Vec<CompileWarning>,
}

impl CompileReport {
    pub fn pass_names(&self) -> Vec<String> {
        self.pass_summaries
            .iter()
            .map(|summary| summary.name.clone())
            .collect()
    }
}

/// Non-fatal compile-time diagnostic.
#[derive(Clone, Debug, Eq, PartialEq)]
pub struct CompileWarning {
    pub pass_index: Option<usize>,
    pub pass_name: Option<String>,
    pub message: String,
}

/// Fatal compile-time error.
#[derive(Clone, Debug, Eq, PartialEq)]
pub struct CompileError {
    pub pass_index: Option<usize>,
    pub pass_name: Option<String>,
    pub kind: CompileErrorKind,
}

/// Stable class of compile-time plan error.
#[derive(Clone, Debug, Eq, PartialEq)]
pub enum CompileErrorKind {
    MissingRefData,
    MissingAssignment {
        segment: Segment,
    },
    InvalidPassOrder {
        reason: String,
    },
    InvalidParameterSupport {
        address: String,
        reason: String,
    },
    ContractPrecondition {
        contract_name: String,
        reason: String,
    },
}

/// All fatal compile-time errors found in one plan scan.
#[derive(Clone, Debug, Eq, PartialEq)]
pub struct CompileErrors {
    pub errors: Vec<CompileError>,
}

impl CompileErrors {
    pub fn is_empty(&self) -> bool {
        self.errors.is_empty()
    }
}

impl fmt::Display for CompileErrors {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        if self.errors.is_empty() {
            return write!(f, "plan compilation failed with no errors");
        }

        write!(f, "plan compilation failed")?;
        for error in &self.errors {
            let location = match (error.pass_index, error.pass_name.as_deref()) {
                (Some(index), Some(name)) => format!("pass {index} ({name})"),
                (Some(index), None) => format!("pass {index}"),
                (None, Some(name)) => name.to_string(),
                (None, None) => "plan".to_string(),
            };
            match &error.kind {
                CompileErrorKind::MissingRefData => {
                    write!(f, "; {location} requires reference data")?
                }
                CompileErrorKind::MissingAssignment { segment } => write!(
                    f,
                    "; {location} requires a prior {:?} allele assignment",
                    segment
                )?,
                CompileErrorKind::InvalidPassOrder { reason } => {
                    write!(f, "; {location} has invalid pass order: {reason}")?
                }
                CompileErrorKind::InvalidParameterSupport { address, reason } => write!(
                    f,
                    "; {location} has invalid parameter support at {address}: {reason}",
                )?,
                CompileErrorKind::ContractPrecondition {
                    contract_name,
                    reason,
                } => write!(
                    f,
                    "; contract {contract_name} precondition failed: {reason}",
                )?,
            }
        }
        Ok(())
    }
}

impl std::error::Error for CompileErrors {}

/// A borrowed compiled simulation artifact ready for deterministic execution.
///
/// This is the low-level boundary used when the caller owns the
/// `PassPlan`, reference data, and contracts elsewhere. Use
/// [`OwnedCompiledSimulator`] when the compiled artifact itself must
/// outlive the builder inputs, such as Python compile-once/run-many
/// usage.
pub struct CompiledSimulator<'a> {
    plan: &'a PassPlan,
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
    refdata: Option<RefDataConfig>,
    contracts: Option<ContractSet>,
    policy: ExecutionPolicy,
    report: CompileReport,
    feasibility: Option<FeasibilityContext>,
    reference_index: Option<ReferenceMatchIndex>,
}

#[derive(Copy, Clone)]
struct ExecutionInputs<'a> {
    plan: &'a PassPlan,
    refdata: Option<&'a RefDataConfig>,
    contracts: Option<&'a ContractSet>,
    feasibility: Option<&'a FeasibilityContext>,
    reference_index: Option<&'a ReferenceMatchIndex>,
    policy: ExecutionPolicy,
}

#[derive(Debug)]
struct ExecutionAbort {
    error: PassError,
    /// Already-committed prefix at the point of failure. Kept private:
    /// public runners expose only the structured pass error, while
    /// internal tests can assert transaction commit boundaries.
    #[cfg_attr(not(test), allow(dead_code))]
    committed: Outcome,
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
        let (report, errors, feasibility) = analyze_plan(plan, refdata, contracts, policy);
        if !errors.is_empty() {
            return Err(CompileErrors { errors });
        }
        let reference_index = refdata.map(ReferenceMatchIndex::build);

        Ok(Self {
            plan,
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
            refdata: self.refdata,
            contracts: self.contracts,
            feasibility: self.feasibility.as_ref(),
            reference_index: self.reference_index.as_ref(),
            policy,
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
        let (report, errors, feasibility) =
            analyze_plan(&plan, refdata.as_ref(), contracts.as_ref(), policy);
        if !errors.is_empty() {
            return Err(CompileErrors { errors });
        }
        let reference_index = refdata.as_ref().map(ReferenceMatchIndex::build);

        Ok(Self {
            plan,
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
        let reference_index = refdata.as_ref().map(ReferenceMatchIndex::build);
        Self {
            plan,
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
            refdata: self.refdata.as_ref(),
            contracts: self.contracts.as_ref(),
            feasibility: self.feasibility.as_ref(),
            reference_index: self.reference_index.as_ref(),
            policy,
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

fn execute_transactional(
    inputs: ExecutionInputs<'_>,
    initial: Simulation,
    seed: u64,
) -> Result<Outcome, ExecutionAbort> {
    let mut trace = Trace::new();
    let mut rng = Rng::new(seed);

    let mut revisions: Vec<Simulation> = Vec::with_capacity(inputs.plan.len() + 1);
    let mut pass_names: Vec<String> = Vec::with_capacity(inputs.plan.len());
    let mut events: Vec<EventRecord> = Vec::with_capacity(inputs.plan.len());

    revisions.push(initial);

    for pass in inputs.plan.passes() {
        let prev = revisions.last().expect("non-empty by construction");
        let pass_index = pass_names.len();
        let pass_name = pass.name().to_string();
        let effects = pass.effects();
        let trace_start = trace.len();
        let pre = StateSummary::from_simulation(prev);
        let mut trace_delta = Trace::new();
        let mut ctx = PassContext {
            trace: &mut trace_delta,
            rng: &mut rng,
            pass_index,
            refdata: inputs.refdata,
            contracts: inputs.contracts,
            feasibility: inputs.feasibility,
        };

        let mut next = match inputs.policy {
            ExecutionPolicy::Permissive => pass.execute(prev, &mut ctx),
            ExecutionPolicy::Strict => match pass.execute_checked(prev, &mut ctx) {
                Ok(next) => next,
                Err(error) => {
                    return Err(abort(error, revisions, pass_names, trace, events));
                }
            },
        };
        drop(ctx);

        if let Some(reference_index) = inputs.reference_index {
            next = apply_live_call_updates(next, &effects, reference_index);
        }

        if inputs.policy == ExecutionPolicy::Strict {
            if let Some(contracts) = inputs.contracts {
                if let Err(violations) = contracts.verify(&next, inputs.refdata) {
                    return Err(abort(
                        PassError::contract_violation(pass_name, violations),
                        revisions,
                        pass_names,
                        trace,
                        events,
                    ));
                }
            }
        }

        let trace_end = trace_start + trace_delta.len();
        let post = StateSummary::from_simulation(&next);
        let event = EventRecord::pass_committed(
            pass_index,
            pass_name.clone(),
            effects,
            TraceSpan::new(trace_start, trace_end),
            pre,
            post,
        );

        trace.append_delta(trace_delta);
        pass_names.push(pass_name);
        revisions.push(next);
        events.push(event);
    }

    Ok(Outcome {
        revisions,
        pass_names,
        trace,
        events,
    })
}

fn apply_live_call_updates(
    mut sim: Simulation,
    effects: &[PassEffect],
    reference_index: &ReferenceMatchIndex,
) -> Simulation {
    for effect in effects {
        match effect {
            PassEffect::AssembleSegment(segment) => {
                sim = with_assembled_segment_live_call(&sim, reference_index, *segment);
                // assembling a downstream segment introduces
                // new bases that an earlier segment's right-extension
                // walker can reach into when its allele suffix happens
                // to match. Retrigger the upstream segment's refresh
                // so any cross-boundary overlap surfaces as
                // OVERLAPS_OTHER_SEGMENT on the upstream hypothesis.
                //
                // - Assembling D → refresh V (V right overlaps into D).
                // - Assembling J → refresh D (D right overlaps into J).
                //
                // No symmetric upstream hook for V right-into-J (we
                // never assemble J without D in VDJ chains, and the
                // walker would have to traverse D's region first;
                // the targeted hops stay conservative).
                match segment {
                    Segment::D => {
                        sim = with_assembled_segment_live_call(
                            &sim,
                            reference_index,
                            Segment::V,
                        );
                    }
                    Segment::J => {
                        sim = with_assembled_segment_live_call(
                            &sim,
                            reference_index,
                            Segment::D,
                        );
                    }
                    _ => {}
                }
            }
            // any base edit (SHM, uniform mutation, PCR, quality
            // / N injection, contaminant overwrite) can change which
            // alleles the assembled bases support. We do a conservative
            // refresh — every assembled V/D/J segment is recomputed from
            // the current pool. Segments without an assembled region are
            // a no-op inside `with_assembled_segment_live_call`.
            //
            // A dirty-window optimization (scanning the trace delta for
            // edited positions and skipping segments outside the dirty
            // range) was evaluated. The trace-scan cost (~5 µs/record)
            // cancelled the refresh savings — net change was 0-2
            // µs/record on typical loads, sometimes a small loss. The
            // simpler full-refresh policy is kept for clarity.
            PassEffect::EditBases => {
                for segment in [Segment::V, Segment::D, Segment::J] {
                    sim = with_assembled_segment_live_call(&sim, reference_index, segment);
                }
            }
            // an NP region appearing right-adjacent to V
            // can extend V's right boundary if the NP bases happen to
            // continue exactly into a V allele's suffix. The walker
            // inside `with_assembled_segment_live_call` does the
            // extension automatically — we just need to retrigger the
            // V refresh now that NP1 exists.
            PassEffect::AppendRegion(Segment::Np1) => {
                sim = with_assembled_segment_live_call(&sim, reference_index, Segment::V);
            }
            // NP2 appears AFTER D is assembled (the typical
            // VDJ pipeline order is `... → assemble(D) → np2 → ...`),
            // so D's right boundary cannot pick up NP2 bases at
            // assembly time. Retrigger D refresh once NP2 exists so
            // the walker can extend D's right side into NP2 if the
            // bases match.
            //
            // J left-extension does NOT need a separate hook here:
            // J is assembled after every NP region exists, so its
            // own `AssembleSegment(J)` refresh already sees the full
            // NP context.
            PassEffect::AppendRegion(Segment::Np2) => {
                sim = with_assembled_segment_live_call(&sim, reference_index, Segment::D);
            }
            // structural indels (insertions / deletions) shift
            // the pool layout under V/D/J coding regions. Indel-inserted
            // nucleotides have NO_GERMLINE_POS so the walker skips them;
            // deletions show up as forward jumps in `germline_pos` which
            // the walker now tolerates. Refresh every assembled V/D/J
            // segment from the post-indel pool so the live calls reflect
            // the new evidence layout.
            PassEffect::StructuralIndel => {
                for segment in [Segment::V, Segment::D, Segment::J] {
                    sim = with_assembled_segment_live_call(&sim, reference_index, segment);
                }
            }
            _ => {}
        }
    }
    sim
}

fn abort(
    error: PassError,
    revisions: Vec<Simulation>,
    pass_names: Vec<String>,
    trace: Trace,
    events: Vec<EventRecord>,
) -> ExecutionAbort {
    ExecutionAbort {
        error,
        committed: Outcome {
            revisions,
            pass_names,
            trace,
            events,
        },
    }
}

#[derive(Clone, Debug)]
struct Located<T> {
    value: T,
    pass_index: usize,
    pass_name: String,
}

impl<T> Located<T> {
    fn new(value: T, pass_index: usize, pass_name: &str) -> Self {
        Self {
            value,
            pass_index,
            pass_name: pass_name.to_string(),
        }
    }

    fn error(&self, kind: CompileErrorKind) -> CompileError {
        CompileError {
            pass_index: Some(self.pass_index),
            pass_name: Some(self.pass_name.clone()),
            kind,
        }
    }
}

#[derive(Default)]
struct CompileFactIndex {
    allele_supports: HashMap<Segment, Located<AlleleIdSupport>>,
    trim_supports: HashMap<(Segment, TrimEnd), Located<IntegerSupport>>,
    np_length_supports: HashMap<Segment, Located<IntegerSupport>>,
    assigned_segments: HashSet<Segment>,
}

impl CompileFactIndex {
    fn record_fact(
        &mut self,
        fact: PassCompileFact,
        pass_index: usize,
        pass_name: &str,
        refdata: Option<&RefDataConfig>,
        errors: &mut Vec<CompileError>,
    ) {
        match fact {
            PassCompileFact::AlleleSampleSupport { segment, support } => {
                let located = Located::new(support, pass_index, pass_name);
                validate_allele_support(segment, &located, refdata, errors);
                self.allele_supports.insert(segment, located);
            }
            PassCompileFact::TrimSupport {
                segment,
                end,
                support,
            } => {
                let located = Located::new(support, pass_index, pass_name);
                validate_trim_support(segment, end, &located, errors);
                self.trim_supports.insert((segment, end), located);
            }
            PassCompileFact::NpLengthSupport { segment, support } => {
                let located = Located::new(support, pass_index, pass_name);
                validate_np_length_support(segment, &located, errors);
                self.np_length_supports.insert(segment, located);
            }
        }
    }
}

fn pass_scoped_error(index: usize, name: &str, kind: CompileErrorKind) -> CompileError {
    CompileError {
        pass_index: Some(index),
        pass_name: Some(name.to_string()),
        kind,
    }
}

fn plan_scoped_error(kind: CompileErrorKind) -> CompileError {
    CompileError {
        pass_index: None,
        pass_name: None,
        kind,
    }
}

fn invalid_parameter(address: impl Into<String>, reason: impl Into<String>) -> CompileErrorKind {
    CompileErrorKind::InvalidParameterSupport {
        address: address.into(),
        reason: reason.into(),
    }
}

fn invalid_pass_order(reason: impl Into<String>) -> CompileErrorKind {
    CompileErrorKind::InvalidPassOrder {
        reason: reason.into(),
    }
}

fn contract_precondition(
    contract_name: impl Into<String>,
    reason: impl Into<String>,
) -> CompileErrorKind {
    CompileErrorKind::ContractPrecondition {
        contract_name: contract_name.into(),
        reason: reason.into(),
    }
}

fn sample_allele_address(segment: Segment) -> &'static str {
    match segment {
        Segment::V => "sample_allele.v",
        Segment::D => "sample_allele.d",
        Segment::J => "sample_allele.j",
        Segment::Np1 | Segment::Np2 => "sample_allele.<invalid>",
    }
}

fn trim_address(segment: Segment, end: TrimEnd) -> &'static str {
    match (segment, end) {
        (Segment::V, TrimEnd::Five) => "trim.v_5",
        (Segment::V, TrimEnd::Three) => "trim.v_3",
        (Segment::D, TrimEnd::Five) => "trim.d_5",
        (Segment::D, TrimEnd::Three) => "trim.d_3",
        (Segment::J, TrimEnd::Five) => "trim.j_5",
        (Segment::J, TrimEnd::Three) => "trim.j_3",
        _ => "trim.<invalid>",
    }
}

fn np_length_address(segment: Segment) -> &'static str {
    match segment {
        Segment::Np1 => "np.np1.length",
        Segment::Np2 => "np.np2.length",
        _ => "np.<invalid>.length",
    }
}

fn validate_allele_support(
    segment: Segment,
    support: &Located<AlleleIdSupport>,
    refdata: Option<&RefDataConfig>,
    errors: &mut Vec<CompileError>,
) {
    match &support.value {
        AlleleIdSupport::Invalid(reason) => {
            errors.push(support.error(invalid_parameter(
                sample_allele_address(segment),
                reason.clone(),
            )));
        }
        AlleleIdSupport::Enumerated(ids) => {
            if ids.is_empty() {
                errors.push(support.error(invalid_parameter(
                    sample_allele_address(segment),
                    "empty_support",
                )));
                return;
            }
            let Some(refdata) = refdata else {
                return;
            };
            let Some(pool) = refdata.pool_for(segment) else {
                return;
            };
            for id in ids {
                if id.as_usize() >= pool.len() {
                    errors.push(support.error(invalid_parameter(
                        sample_allele_address(segment),
                        format!("allele_id_{}_out_of_range", id.index()),
                    )));
                }
            }
        }
        AlleleIdSupport::Unavailable => {}
    }
}

fn validate_trim_support(
    segment: Segment,
    end: TrimEnd,
    support: &Located<IntegerSupport>,
    errors: &mut Vec<CompileError>,
) {
    let address = trim_address(segment, end);
    match &support.value {
        IntegerSupport::Invalid(reason) => {
            errors.push(support.error(invalid_parameter(address, reason.clone())));
        }
        IntegerSupport::Enumerated(values) => {
            if values.is_empty() {
                errors.push(support.error(invalid_parameter(address, "empty_support")));
                return;
            }
            for value in values {
                if *value < 0 {
                    errors.push(support.error(invalid_parameter(address, "negative_trim")));
                } else if *value > u16::MAX as i64 {
                    errors.push(support.error(invalid_parameter(address, "trim_exceeds_u16")));
                }
            }
        }
        IntegerSupport::Unavailable => {}
    }
}

fn validate_np_length_support(
    segment: Segment,
    support: &Located<IntegerSupport>,
    errors: &mut Vec<CompileError>,
) {
    let address = np_length_address(segment);
    match &support.value {
        IntegerSupport::Invalid(reason) => {
            errors.push(support.error(invalid_parameter(address, reason.clone())));
        }
        IntegerSupport::Enumerated(values) => {
            if values.is_empty() {
                errors.push(support.error(invalid_parameter(address, "empty_support")));
                return;
            }
            for value in values {
                if *value < 0 {
                    errors.push(support.error(invalid_parameter(address, "negative_length")));
                } else if *value > u32::MAX as i64 {
                    errors.push(support.error(invalid_parameter(address, "length_exceeds_u32")));
                }
            }
        }
        IntegerSupport::Unavailable => {}
    }
}

fn validate_contract_preconditions(
    contracts: Option<&ContractSet>,
    refdata: Option<&RefDataConfig>,
    facts: &CompileFactIndex,
    errors: &mut Vec<CompileError>,
) {
    let Some(contracts) = contracts else {
        return;
    };

    if contracts.contains_kind(ContractKind::ProductiveJunctionFrame) {
        validate_productive_frame_preconditions(refdata, facts, errors);
    }
}

fn validate_productive_frame_preconditions(
    refdata: Option<&RefDataConfig>,
    facts: &CompileFactIndex,
    errors: &mut Vec<CompileError>,
) {
    let contract_name = "productive_junction_frame";
    let Some(refdata) = refdata else {
        errors.push(plan_scoped_error(contract_precondition(
            contract_name,
            "productive frame validation requires reference data",
        )));
        return;
    };

    let Some(v_residues) = anchor_tail_residues(Segment::V, refdata, facts, errors) else {
        return;
    };
    let Some(j_residues) = anchor_head_residues(Segment::J, refdata, facts, errors) else {
        return;
    };

    let uses_d = facts.assigned_segments.contains(&Segment::D);
    if uses_d {
        let Some(np1_residues) = optional_np_length_residues(Segment::Np1, facts, errors) else {
            return;
        };
        let Some(d_residues) = d_length_residues(refdata, facts, errors) else {
            return;
        };
        let Some(np2_residues) = required_np_length_residues(Segment::Np2, facts, errors) else {
            return;
        };

        if residues_admit_frame(&[
            &v_residues,
            &np1_residues,
            &d_residues,
            &np2_residues,
            &j_residues,
        ]) {
            return;
        }

        errors.push(plan_scoped_error(contract_precondition(
            contract_name,
            "NP2 length support has no in-frame mass for the active V/D/J anchor and trim supports",
        )));
    } else {
        let Some(np1_residues) = required_np_length_residues(Segment::Np1, facts, errors) else {
            return;
        };

        if residues_admit_frame(&[&v_residues, &np1_residues, &j_residues]) {
            return;
        }

        errors.push(plan_scoped_error(contract_precondition(
            contract_name,
            "NP1 length support has no in-frame mass for the active V/J anchor and trim supports",
        )));
    }
}

fn allele_ids_for_segment(
    segment: Segment,
    refdata: &RefDataConfig,
    facts: &CompileFactIndex,
) -> Vec<AlleleId> {
    let Some(pool) = refdata.pool_for(segment) else {
        return Vec::new();
    };

    match facts
        .allele_supports
        .get(&segment)
        .map(|entry| &entry.value)
    {
        Some(AlleleIdSupport::Enumerated(ids)) => ids
            .iter()
            .copied()
            .filter(|id| id.as_usize() < pool.len())
            .collect(),
        Some(AlleleIdSupport::Invalid(_)) => Vec::new(),
        Some(AlleleIdSupport::Unavailable) | None => pool.iter().map(|(id, _)| id).collect(),
    }
}

fn trim_values(
    segment: Segment,
    end: TrimEnd,
    facts: &CompileFactIndex,
    errors: &mut Vec<CompileError>,
) -> Option<Vec<u32>> {
    let Some(entry) = facts.trim_supports.get(&(segment, end)) else {
        return Some(vec![0]);
    };

    match &entry.value {
        IntegerSupport::Enumerated(values) => {
            let values: Vec<u32> = values
                .iter()
                .copied()
                .filter(|value| (0..=u16::MAX as i64).contains(value))
                .map(|value| value as u32)
                .collect();
            if values.is_empty() {
                None
            } else {
                Some(values)
            }
        }
        IntegerSupport::Unavailable => {
            errors.push(plan_scoped_error(contract_precondition(
                "productive_junction_frame",
                format!(
                    "{} support is not enumerable, so productive frame preconditions cannot be validated",
                    trim_address(segment, end)
                ),
            )));
            None
        }
        IntegerSupport::Invalid(_) => None,
    }
}

fn required_np_length_residues(
    segment: Segment,
    facts: &CompileFactIndex,
    errors: &mut Vec<CompileError>,
) -> Option<HashSet<u8>> {
    let Some(entry) = facts.np_length_supports.get(&segment) else {
        errors.push(plan_scoped_error(contract_precondition(
            "productive_junction_frame",
            format!(
                "{} support is required for productive frame validation",
                np_length_address(segment)
            ),
        )));
        return None;
    };
    np_length_residues_from_entry(segment, entry, true, errors)
}

fn optional_np_length_residues(
    segment: Segment,
    facts: &CompileFactIndex,
    errors: &mut Vec<CompileError>,
) -> Option<HashSet<u8>> {
    let Some(entry) = facts.np_length_supports.get(&segment) else {
        return Some(HashSet::from([0]));
    };
    np_length_residues_from_entry(segment, entry, false, errors)
}

fn np_length_residues_from_entry(
    segment: Segment,
    entry: &Located<IntegerSupport>,
    required: bool,
    errors: &mut Vec<CompileError>,
) -> Option<HashSet<u8>> {
    match &entry.value {
        IntegerSupport::Enumerated(values) => {
            let residues: HashSet<u8> = values
                .iter()
                .copied()
                .filter(|value| (0..=u32::MAX as i64).contains(value))
                .map(|value| (value % 3) as u8)
                .collect();
            if residues.is_empty() {
                if required {
                    errors.push(plan_scoped_error(contract_precondition(
                        "productive_junction_frame",
                        format!(
                            "{} has no valid non-negative support",
                            np_length_address(segment)
                        ),
                    )));
                }
                None
            } else {
                Some(residues)
            }
        }
        IntegerSupport::Unavailable => {
            errors.push(plan_scoped_error(contract_precondition(
                "productive_junction_frame",
                format!(
                    "{} support is not enumerable, so productive frame preconditions cannot be validated",
                    np_length_address(segment)
                ),
            )));
            None
        }
        IntegerSupport::Invalid(_) => None,
    }
}

fn anchor_tail_residues(
    segment: Segment,
    refdata: &RefDataConfig,
    facts: &CompileFactIndex,
    errors: &mut Vec<CompileError>,
) -> Option<HashSet<u8>> {
    let ids = allele_ids_for_segment(segment, refdata, facts);
    let anchored_count = ids
        .iter()
        .filter(|id| {
            refdata
                .get(segment, **id)
                .is_some_and(|allele| valid_anchor(allele.anchor, allele.len()))
        })
        .count();
    if anchored_count == 0 {
        errors.push(plan_scoped_error(contract_precondition(
            "productive_junction_frame",
            format!(
                "{:?} sample support has no alleles with valid anchors",
                segment
            ),
        )));
        return None;
    }

    let trim_5 = trim_values(segment, TrimEnd::Five, facts, errors)?;
    let trim_3 = trim_values(segment, TrimEnd::Three, facts, errors)?;
    let mut residues = HashSet::new();
    for id in ids {
        let Some(allele) = refdata.get(segment, id) else {
            continue;
        };
        let Some(anchor) = allele.anchor.map(|anchor| anchor as u32) else {
            continue;
        };
        if !valid_anchor(Some(anchor as u16), allele.len()) {
            continue;
        }
        for five in &trim_5 {
            for three in &trim_3 {
                let retained_end = allele.len().saturating_sub(*three);
                if *five <= anchor && anchor + 3 <= retained_end {
                    residues.insert(((retained_end - anchor) % 3) as u8);
                }
            }
        }
    }

    if residues.is_empty() {
        errors.push(plan_scoped_error(contract_precondition(
            "productive_junction_frame",
            format!("{:?} trim support removes every valid anchor", segment),
        )));
        None
    } else {
        Some(residues)
    }
}

fn anchor_head_residues(
    segment: Segment,
    refdata: &RefDataConfig,
    facts: &CompileFactIndex,
    errors: &mut Vec<CompileError>,
) -> Option<HashSet<u8>> {
    let ids = allele_ids_for_segment(segment, refdata, facts);
    let anchored_count = ids
        .iter()
        .filter(|id| {
            refdata
                .get(segment, **id)
                .is_some_and(|allele| valid_anchor(allele.anchor, allele.len()))
        })
        .count();
    if anchored_count == 0 {
        errors.push(plan_scoped_error(contract_precondition(
            "productive_junction_frame",
            format!(
                "{:?} sample support has no alleles with valid anchors",
                segment
            ),
        )));
        return None;
    }

    let trim_5 = trim_values(segment, TrimEnd::Five, facts, errors)?;
    let trim_3 = trim_values(segment, TrimEnd::Three, facts, errors)?;
    let mut residues = HashSet::new();
    for id in ids {
        let Some(allele) = refdata.get(segment, id) else {
            continue;
        };
        let Some(anchor) = allele.anchor.map(|anchor| anchor as u32) else {
            continue;
        };
        if !valid_anchor(Some(anchor as u16), allele.len()) {
            continue;
        }
        for five in &trim_5 {
            for three in &trim_3 {
                let retained_end = allele.len().saturating_sub(*three);
                if *five <= anchor && anchor + 3 <= retained_end {
                    residues.insert(((anchor + 3 - *five) % 3) as u8);
                }
            }
        }
    }

    if residues.is_empty() {
        errors.push(plan_scoped_error(contract_precondition(
            "productive_junction_frame",
            format!("{:?} trim support removes every valid anchor", segment),
        )));
        None
    } else {
        Some(residues)
    }
}

fn d_length_residues(
    refdata: &RefDataConfig,
    facts: &CompileFactIndex,
    errors: &mut Vec<CompileError>,
) -> Option<HashSet<u8>> {
    let ids = allele_ids_for_segment(Segment::D, refdata, facts);
    if ids.is_empty() {
        errors.push(plan_scoped_error(contract_precondition(
            "productive_junction_frame",
            "D sample support is empty for a D-containing productive plan",
        )));
        return None;
    }

    let trim_5 = trim_values(Segment::D, TrimEnd::Five, facts, errors)?;
    let trim_3 = trim_values(Segment::D, TrimEnd::Three, facts, errors)?;
    let mut residues = HashSet::new();
    for id in ids {
        let Some(allele) = refdata.get(Segment::D, id) else {
            continue;
        };
        for five in &trim_5 {
            for three in &trim_3 {
                if *five + *three <= allele.len() {
                    residues.insert(((allele.len() - *five - *three) % 3) as u8);
                }
            }
        }
    }

    if residues.is_empty() {
        errors.push(plan_scoped_error(contract_precondition(
            "productive_junction_frame",
            "D trim support leaves no valid retained D lengths",
        )));
        None
    } else {
        Some(residues)
    }
}

fn valid_anchor(anchor: Option<u16>, allele_len: u32) -> bool {
    anchor.is_some_and(|anchor| (anchor as u32) + 3 <= allele_len)
}

fn residues_admit_frame(sets: &[&HashSet<u8>]) -> bool {
    fn walk(sets: &[&HashSet<u8>], index: usize, residue: u8) -> bool {
        if index == sets.len() {
            return residue % 3 == 0;
        }
        sets[index]
            .iter()
            .any(|next| walk(sets, index + 1, (residue + *next) % 3))
    }
    walk(sets, 0, 0)
}

fn build_feasibility_context(
    contracts: Option<&ContractSet>,
    refdata: Option<&RefDataConfig>,
    facts: &CompileFactIndex,
) -> Option<FeasibilityContext> {
    let contracts = contracts?;
    if !contracts.contains_kind(ContractKind::ProductiveJunctionFrame) {
        return None;
    }
    if facts.assigned_segments.contains(&Segment::D) {
        return None;
    }

    let refdata = refdata?;
    let vj_productive = build_vj_productive_feasibility(refdata, facts)?;
    let context = FeasibilityContext::empty().with_vj_productive(vj_productive);
    if context.is_empty() {
        None
    } else {
        Some(context)
    }
}

fn build_vj_productive_feasibility(
    refdata: &RefDataConfig,
    facts: &CompileFactIndex,
) -> Option<VjProductiveFeasibility> {
    Some(VjProductiveFeasibility::new(
        allele_domain_for_feasibility(Segment::V, refdata, facts)?,
        allele_domain_for_feasibility(Segment::J, refdata, facts)?,
        trim_domain_for_feasibility(Segment::V, TrimEnd::Five, facts)?,
        trim_domain_for_feasibility(Segment::V, TrimEnd::Three, facts)?,
        trim_domain_for_feasibility(Segment::J, TrimEnd::Five, facts)?,
        trim_domain_for_feasibility(Segment::J, TrimEnd::Three, facts)?,
        np_length_domain_for_feasibility(Segment::Np1, facts)?,
    ))
}

fn allele_domain_for_feasibility(
    segment: Segment,
    refdata: &RefDataConfig,
    facts: &CompileFactIndex,
) -> Option<ChoiceDomain<AlleleId>> {
    let support = facts.allele_supports.get(&segment)?;
    let values = unique_allele_ids(allele_ids_for_segment(segment, refdata, facts));
    if values.is_empty() {
        None
    } else {
        Some(ChoiceDomain::new(Some(support.pass_index), values))
    }
}

fn trim_domain_for_feasibility(
    segment: Segment,
    end: TrimEnd,
    facts: &CompileFactIndex,
) -> Option<ChoiceDomain<u32>> {
    let Some(entry) = facts.trim_supports.get(&(segment, end)) else {
        return Some(ChoiceDomain::new(None, vec![0]));
    };

    let IntegerSupport::Enumerated(values) = &entry.value else {
        return None;
    };

    let values = unique_u32_values(
        values
            .iter()
            .copied()
            .filter(|value| (0..=u16::MAX as i64).contains(value))
            .map(|value| value as u32),
    );
    if values.is_empty() {
        None
    } else {
        Some(ChoiceDomain::new(Some(entry.pass_index), values))
    }
}

fn np_length_domain_for_feasibility(
    segment: Segment,
    facts: &CompileFactIndex,
) -> Option<Vec<u32>> {
    let entry = facts.np_length_supports.get(&segment)?;
    let IntegerSupport::Enumerated(values) = &entry.value else {
        return None;
    };

    let lengths = unique_u32_values(
        values
            .iter()
            .copied()
            .filter(|value| (0..=u32::MAX as i64).contains(value))
            .map(|value| value as u32),
    );

    if lengths.is_empty() {
        None
    } else {
        Some(lengths)
    }
}

fn unique_allele_ids(values: Vec<AlleleId>) -> Vec<AlleleId> {
    let mut seen = HashSet::new();
    let mut unique = Vec::new();
    for value in values {
        if seen.insert(value) {
            unique.push(value);
        }
    }
    unique
}

fn unique_u32_values(values: impl IntoIterator<Item = u32>) -> Vec<u32> {
    let mut seen = HashSet::new();
    let mut unique = Vec::new();
    for value in values {
        if seen.insert(value) {
            unique.push(value);
        }
    }
    unique
}

fn analyze_plan(
    plan: &PassPlan,
    refdata: Option<&RefDataConfig>,
    contracts: Option<&ContractSet>,
    policy: ExecutionPolicy,
) -> (CompileReport, Vec<CompileError>, Option<FeasibilityContext>) {
    let mut assigned: HashSet<Segment> = HashSet::new();
    let mut pass_summaries = Vec::with_capacity(plan.len());
    let mut declared_choices = Vec::new();
    let mut errors = Vec::new();
    let mut fact_index = CompileFactIndex::default();
    let mut assembled: HashSet<Segment> = HashSet::new();

    for (index, pass) in plan.passes().iter().enumerate() {
        let name = pass.name().to_string();
        let choices = pass.declared_choices();
        let requirements = pass.requirements();
        let effects = pass.effects();
        let compile_facts = pass.compile_facts();

        for req in &requirements {
            match req {
                PassRequirement::RefData if refdata.is_none() => {
                    errors.push(pass_scoped_error(
                        index,
                        &name,
                        CompileErrorKind::MissingRefData,
                    ));
                }
                PassRequirement::AlleleAssignment(segment) if !assigned.contains(segment) => {
                    errors.push(pass_scoped_error(
                        index,
                        &name,
                        CompileErrorKind::MissingAssignment { segment: *segment },
                    ));
                }
                _ => {}
            }
        }

        for fact in compile_facts {
            fact_index.record_fact(fact, index, &name, refdata, &mut errors);
        }

        for address in &choices {
            declared_choices.push(DeclaredChoice {
                pass_index: index,
                pass_name: name.clone(),
                address: address.clone(),
            });
        }

        for effect in &effects {
            match effect {
                PassEffect::AssignAllele(segment) => {
                    assigned.insert(*segment);
                    fact_index.assigned_segments.insert(*segment);
                }
                PassEffect::TrimAllele(segment) if assembled.contains(segment) => {
                    errors.push(pass_scoped_error(
                        index,
                        &name,
                        invalid_pass_order(format!("trim_after_assembly.{:?}", segment)),
                    ));
                }
                PassEffect::AssembleSegment(segment) => {
                    assembled.insert(*segment);
                }
                _ => {}
            }
        }

        pass_summaries.push(PassSummary {
            index,
            name,
            declared_choices: choices,
            requirements,
            effects,
        });
    }

    let active_contracts = contracts
        .map(|set| {
            set.iter()
                .map(|contract| contract.name().to_string())
                .collect()
        })
        .unwrap_or_default();

    validate_contract_preconditions(contracts, refdata, &fact_index, &mut errors);
    let feasibility = if errors.is_empty() {
        build_feasibility_context(contracts, refdata, &fact_index)
    } else {
        None
    };

    (
        CompileReport {
            policy,
            pass_summaries,
            declared_choices,
            active_contracts,
            warnings: Vec::new(),
        },
        errors,
        feasibility,
    )
}


#[cfg(test)]
#[path = "compiled_tests.rs"]
mod tests;
