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
mod tests {
    use super::*;
    use crate::airr_record::build_airr_record;
    use crate::assignment::TrimEnd;
    use crate::contract::{Contract, ContractViolation};
    use crate::dist::{AllelePoolDist, EmpiricalLengthDist, UniformBase};
    use crate::ir::{NucFlags, Nucleotide, Segment};
    use crate::pass::{Pass, PassRuntime};
    use crate::passes::{AssembleSegmentPass, GenerateNPPass, SampleAllelePass, TrimPass};
    use crate::refdata::{Allele, ChainType};
    use crate::trace::ChoiceValue;

    fn vj_refdata() -> RefDataConfig {
        let mut cfg = RefDataConfig::empty(ChainType::Vj);
        let _ = cfg.v_pool.push(Allele {
            name: "v*01".into(),
            gene: "v".into(),
            seq: b"AAACCCGGG".to_vec(),
            segment: Segment::V,
            anchor: Some(6),
        });
        let _ = cfg.j_pool.push(Allele {
            name: "j*01".into(),
            gene: "j".into(),
            seq: b"TTTAAA".to_vec(),
            segment: Segment::J,
            anchor: Some(0),
        });
        cfg
    }

    fn vj_plan(cfg: &RefDataConfig) -> PassPlan {
        vj_plan_with_np_lengths(cfg, vec![(0, 1.0)])
    }

    fn vj_plan_with_np_lengths(cfg: &RefDataConfig, np_lengths: Vec<(i64, f64)>) -> PassPlan {
        let mut plan = PassPlan::new();
        plan.push(Box::new(SampleAllelePass::new(
            Segment::V,
            Box::new(AllelePoolDist::uniform(&cfg.v_pool)),
        )));
        plan.push(Box::new(SampleAllelePass::new(
            Segment::J,
            Box::new(AllelePoolDist::uniform(&cfg.j_pool)),
        )));
        plan.push(Box::new(AssembleSegmentPass::new(Segment::V)));
        plan.push(Box::new(GenerateNPPass::new(
            Segment::Np1,
            Box::new(EmpiricalLengthDist::from_pairs(np_lengths)),
            Box::new(UniformBase),
        )));
        plan.push(Box::new(AssembleSegmentPass::new(Segment::J)));
        plan
    }

    fn vj_refdata_with_runtime_j_residue() -> (RefDataConfig, AlleleId, AlleleId) {
        let mut cfg = RefDataConfig::empty(ChainType::Vj);
        let _ = cfg.v_pool.push(Allele {
            name: "v*01".into(),
            gene: "v".into(),
            seq: b"AAACCCGGG".to_vec(),
            segment: Segment::V,
            anchor: Some(6),
        });
        let j_good = cfg.j_pool.push(Allele {
            name: "j_good*01".into(),
            gene: "j_good".into(),
            seq: b"TTTAAA".to_vec(),
            segment: Segment::J,
            anchor: Some(0),
        });
        let j_bad = cfg.j_pool.push(Allele {
            name: "j_bad*01".into(),
            gene: "j_bad".into(),
            seq: b"TTTAAA".to_vec(),
            segment: Segment::J,
            anchor: Some(1),
        });
        (cfg, j_good, j_bad)
    }

    struct AppendBasePass {
        name: &'static str,
        base: u8,
    }

    impl AppendBasePass {
        fn new(name: &'static str, base: u8) -> Self {
            Self { name, base }
        }
    }

    impl Pass for AppendBasePass {
        fn name(&self) -> &str {
            self.name
        }

        fn execute(&self, sim: &Simulation, _ctx: &mut PassContext) -> Simulation {
            let (next, _handle) = sim.with_nucleotide_pushed(Nucleotide::synthetic(
                self.base,
                Segment::Np1,
                NucFlags::empty(),
            ));
            next
        }

        fn effects(&self) -> Vec<PassEffect> {
            vec![PassEffect::AppendNucleotides]
        }
    }

    struct TraceProbePass {
        name: &'static str,
        address: &'static str,
        base: u8,
    }

    impl TraceProbePass {
        fn new(name: &'static str, address: &'static str, base: u8) -> Self {
            Self {
                name,
                address,
                base,
            }
        }
    }

    impl Pass for TraceProbePass {
        fn name(&self) -> &str {
            self.name
        }

        fn execute(&self, sim: &Simulation, ctx: &mut PassContext) -> Simulation {
            assert!(
                ctx.trace.is_empty(),
                "{} saw previously committed trace records",
                self.name
            );
            ctx.trace.record(self.address, ChoiceValue::Base(self.base));
            let (next, _handle) = sim.with_nucleotide_pushed(Nucleotide::synthetic(
                self.base,
                Segment::Np1,
                NucFlags::empty(),
            ));
            next
        }

        fn declared_choices(&self) -> Vec<String> {
            vec![self.address.to_string()]
        }

        fn effects(&self) -> Vec<PassEffect> {
            vec![PassEffect::AppendNucleotides]
        }
    }

    struct MaxPoolLen {
        max: usize,
    }

    impl MaxPoolLen {
        fn new(max: usize) -> Self {
            Self { max }
        }
    }

    impl Contract for MaxPoolLen {
        fn name(&self) -> &str {
            "test.max_pool_len"
        }

        fn verify(
            &self,
            sim: &Simulation,
            _refdata: Option<&RefDataConfig>,
        ) -> Result<(), ContractViolation> {
            if sim.pool.len() <= self.max {
                Ok(())
            } else {
                Err(ContractViolation::new(
                    self.name(),
                    format!("pool length {} exceeds {}", sim.pool.len(), self.max),
                ))
            }
        }
    }

    #[test]
    fn compiled_simulator_collects_report_metadata() {
        let cfg = vj_refdata();
        let plan = vj_plan(&cfg);
        let contracts = crate::contract::productive();
        let compiled = CompiledSimulator::compile(
            &plan,
            Some(&cfg),
            Some(&contracts),
            ExecutionPolicy::Strict,
        )
        .expect("plan should compile");

        assert_eq!(compiled.policy(), ExecutionPolicy::Strict);
        assert_eq!(
            compiled.report().pass_names(),
            vec![
                "sample_allele.v",
                "sample_allele.j",
                "assemble.v",
                "generate_np.np1",
                "assemble.j",
            ]
        );
        assert!(compiled
            .report()
            .declared_choices
            .iter()
            .any(|choice| choice.address == "np.np1.length"));
        assert_eq!(
            compiled.report().active_contracts,
            vec![
                "productive_junction_frame",
                "no_stop_codon_in_junction",
                "anchor_preserved.v",
                "anchor_preserved.j",
            ]
        );
    }

    #[test]
    fn compiled_simulators_own_reference_match_index_when_refdata_exists() {
        let cfg = vj_refdata();
        let plan = vj_plan(&cfg);
        let compiled =
            CompiledSimulator::compile(&plan, Some(&cfg), None, ExecutionPolicy::Permissive)
                .expect("plan should compile");

        let index = compiled
            .reference_index()
            .expect("refdata-backed compile should build a reference index");
        assert_eq!(index.v.allele_count(), cfg.v_pool.len());
        assert_eq!(index.j.allele_count(), cfg.j_pool.len());
        assert_eq!(index.d.allele_count(), 0);
        assert_eq!(
            index.v.kmer_hits(b"AAACCCG").unwrap(),
            &[crate::live_call::KmerHit {
                allele_id: AlleleId::new(0),
                ref_pos: 0,
            }]
        );

        let owned =
            OwnedCompiledSimulator::compile(plan, Some(cfg), None, ExecutionPolicy::Permissive)
                .expect("owned compile should also build the index");
        assert!(owned.reference_index().is_some());

        let empty_plan = PassPlan::new();
        let no_refdata =
            CompiledSimulator::compile(&empty_plan, None, None, ExecutionPolicy::Permissive)
                .expect("empty plan without refdata should compile");
        assert!(no_refdata.reference_index().is_none());
    }

    #[test]
    fn compiled_assembly_initializes_exact_live_call_from_current_region() {
        let mut cfg = RefDataConfig::empty(ChainType::Vdj);
        let v0 = cfg.v_pool.push(Allele {
            name: "v1*01".into(),
            gene: "v1".into(),
            seq: b"GGAAACCC".to_vec(),
            segment: Segment::V,
            anchor: None,
        });
        let v1 = cfg.v_pool.push(Allele {
            name: "v2*01".into(),
            gene: "v2".into(),
            seq: b"TTAAACCC".to_vec(),
            segment: Segment::V,
            anchor: None,
        });

        let mut plan = PassPlan::new();
        plan.push(Box::new(SampleAllelePass::new(
            Segment::V,
            Box::new(AllelePoolDist::restricted_uniform(&cfg.v_pool, vec![v0])),
        )));
        plan.push(Box::new(TrimPass::new(
            Segment::V,
            TrimEnd::Five,
            Box::new(EmpiricalLengthDist::from_pairs(vec![(2, 1.0)])),
        )));
        plan.push(Box::new(AssembleSegmentPass::new(Segment::V)));

        let compiled =
            CompiledSimulator::compile(&plan, Some(&cfg), None, ExecutionPolicy::Permissive)
                .expect("plan should compile");
        let outcome = compiled.run_one(0).expect("run should succeed");
        let final_sim = outcome.final_simulation();
        let live = final_sim
            .live_calls
            .as_ref()
            .expect("compiled assembly should initialize live calls");
        let v_call = live.get(Segment::V).expect("V live call should exist");

        assert_eq!(final_sim.assignments.get(Segment::V).unwrap().allele_id, v0);
        assert_eq!(v_call.allele_call.to_ids(), vec![v0, v1]);
        assert_eq!(
            v_call.confidence,
            crate::live_call::LiveCallConfidence::ExactAmbiguous
        );
        assert_eq!(
            v_call.boundary_summary.ref_start,
            crate::live_call::BoundaryValue::Single(2)
        );
        assert_eq!(
            v_call.boundary_summary.ref_end,
            crate::live_call::BoundaryValue::Single(8)
        );
        assert_eq!(
            v_call.boundary_summary.seq_start,
            crate::live_call::BoundaryValue::Single(0)
        );
        assert_eq!(
            v_call.boundary_summary.seq_end,
            crate::live_call::BoundaryValue::Single(6)
        );
    }

    #[test]
    fn compiled_v_three_prime_trim_widens_live_and_airr_call() {
        let mut cfg = RefDataConfig::empty(ChainType::Vdj);
        let shared_prefix = b"AAACCCGGGTTT";
        let suffixes: [&[u8]; 5] = [
            b"AAAAAAAAAA",
            b"CCCCCCCCCC",
            b"GGGGGGGGGG",
            b"TTTTTTTTTT",
            b"ACGTACGTAC",
        ];
        let mut ids = Vec::new();
        for (index, suffix) in suffixes.iter().enumerate() {
            let mut seq = shared_prefix.to_vec();
            seq.extend_from_slice(suffix);
            ids.push(cfg.v_pool.push(Allele {
                name: format!("IGHVtrim*0{}", index + 1),
                gene: "IGHVtrim".into(),
                seq,
                segment: Segment::V,
                anchor: None,
            }));
        }

        let mut plan = PassPlan::new();
        plan.push(Box::new(SampleAllelePass::new(
            Segment::V,
            Box::new(AllelePoolDist::restricted_uniform(
                &cfg.v_pool,
                vec![ids[0]],
            )),
        )));
        plan.push(Box::new(TrimPass::new(
            Segment::V,
            TrimEnd::Three,
            Box::new(EmpiricalLengthDist::from_pairs(vec![(10, 1.0)])),
        )));
        plan.push(Box::new(AssembleSegmentPass::new(Segment::V)));

        let compiled =
            CompiledSimulator::compile(&plan, Some(&cfg), None, ExecutionPolicy::Permissive)
                .expect("trim-before-assembly plan should compile");
        let outcome = compiled.run_one(0).expect("run should succeed");
        let final_sim = outcome.final_simulation();
        let live = final_sim.live_calls.as_ref().expect("live calls exist");
        let v_call = live.get(Segment::V).expect("V live call exists");

        assert_eq!(
            final_sim.assignments.get(Segment::V).unwrap().allele_id,
            ids[0]
        );
        assert_eq!(v_call.allele_call.to_ids(), ids);
        assert_eq!(
            v_call.confidence,
            crate::live_call::LiveCallConfidence::ExactAmbiguous
        );
        assert_eq!(
            v_call.boundary_summary.ref_end,
            crate::live_call::BoundaryValue::Single(shared_prefix.len() as u32)
        );

        let rec = build_airr_record(&outcome, &cfg, "trim-ambiguity");
        assert_eq!(
            rec.v_call,
            "IGHVtrim*01,IGHVtrim*02,IGHVtrim*03,IGHVtrim*04,IGHVtrim*05"
        );
    }

    #[test]
    fn compiled_simulator_rejects_trim_after_assembly() {
        let cfg = vj_refdata();
        let mut plan = PassPlan::new();
        plan.push(Box::new(SampleAllelePass::new(
            Segment::V,
            Box::new(AllelePoolDist::uniform(&cfg.v_pool)),
        )));
        plan.push(Box::new(AssembleSegmentPass::new(Segment::V)));
        plan.push(Box::new(TrimPass::new(
            Segment::V,
            TrimEnd::Three,
            Box::new(EmpiricalLengthDist::from_pairs(vec![(1, 1.0)])),
        )));

        let err = match CompiledSimulator::compile(
            &plan,
            Some(&cfg),
            None,
            ExecutionPolicy::Permissive,
        ) {
            Ok(_) => panic!("trim after assembly should fail compilation"),
            Err(err) => err,
        };

        assert!(err.errors.iter().any(|error| {
            matches!(
                &error.kind,
                CompileErrorKind::InvalidPassOrder { reason }
                    if reason == "trim_after_assembly.V"
            )
        }));
    }

    #[test]
    fn compiled_simulator_rejects_assemble_without_refdata() {
        let cfg = vj_refdata();
        let plan = vj_plan(&cfg);
        let err = match CompiledSimulator::compile(&plan, None, None, ExecutionPolicy::Permissive) {
            Ok(_) => panic!("assembly without refdata should fail at compile time"),
            Err(err) => err,
        };

        assert!(err
            .errors
            .iter()
            .any(|e| matches!(e.kind, CompileErrorKind::MissingRefData)));
    }

    #[test]
    fn compiled_simulator_rejects_trim_before_sample_allele() {
        let mut plan = PassPlan::new();
        plan.push(Box::new(TrimPass::new(
            Segment::V,
            TrimEnd::Three,
            Box::new(EmpiricalLengthDist::from_pairs(vec![(1, 1.0)])),
        )));

        let err = match CompiledSimulator::compile(&plan, None, None, ExecutionPolicy::Permissive) {
            Ok(_) => panic!("trim before sample allele should fail"),
            Err(err) => err,
        };
        assert_eq!(
            err.errors[0].kind,
            CompileErrorKind::MissingAssignment {
                segment: Segment::V
            }
        );
    }

    #[test]
    fn compiled_simulator_rejects_negative_np_length_support() {
        let cfg = vj_refdata();
        let plan = vj_plan_with_np_lengths(&cfg, vec![(-1, 1.0), (0, 1.0)]);

        let err = match CompiledSimulator::compile(
            &plan,
            Some(&cfg),
            None,
            ExecutionPolicy::Permissive,
        ) {
            Ok(_) => panic!("negative NP length support should fail at compile time"),
            Err(err) => err,
        };

        assert!(err.errors.iter().any(|e| matches!(
            &e.kind,
            CompileErrorKind::InvalidParameterSupport { address, reason }
                if address == "np.np1.length" && reason == "negative_length"
        )));
    }

    #[test]
    fn compiled_simulator_rejects_productive_when_v_support_has_no_anchor() {
        let mut cfg = RefDataConfig::empty(ChainType::Vj);
        let _ = cfg.v_pool.push(Allele {
            name: "v_anchorless*01".into(),
            gene: "v_anchorless".into(),
            seq: b"AAACCCGGG".to_vec(),
            segment: Segment::V,
            anchor: None,
        });
        let _ = cfg.j_pool.push(Allele {
            name: "j*01".into(),
            gene: "j".into(),
            seq: b"TTTAAA".to_vec(),
            segment: Segment::J,
            anchor: Some(0),
        });
        let plan = vj_plan_with_np_lengths(&cfg, vec![(0, 1.0), (3, 1.0)]);
        let contracts = crate::contract::productive();

        let err = match CompiledSimulator::compile(
            &plan,
            Some(&cfg),
            Some(&contracts),
            ExecutionPolicy::Permissive,
        ) {
            Ok(_) => panic!("productive compile should reject anchorless V support"),
            Err(err) => err,
        };

        assert!(err.errors.iter().any(|e| matches!(
            &e.kind,
            CompileErrorKind::ContractPrecondition {
                contract_name,
                reason,
            } if contract_name == "productive_junction_frame"
                && reason.contains("V sample support has no alleles with valid anchors")
        )));
    }

    #[test]
    fn compiled_simulator_rejects_productive_without_in_frame_np_mass() {
        let cfg = vj_refdata();
        let plan = vj_plan_with_np_lengths(&cfg, vec![(1, 1.0), (2, 1.0)]);
        let contracts = crate::contract::productive();

        let err = match CompiledSimulator::compile(
            &plan,
            Some(&cfg),
            Some(&contracts),
            ExecutionPolicy::Permissive,
        ) {
            Ok(_) => panic!("productive compile should reject NP support with no in-frame mass"),
            Err(err) => err,
        };

        assert!(err.errors.iter().any(|e| matches!(
            &e.kind,
            CompileErrorKind::ContractPrecondition {
                contract_name,
                reason,
            } if contract_name == "productive_junction_frame"
                && reason.contains("NP1 length support has no in-frame mass")
        )));
    }

    #[test]
    fn productive_runtime_filters_anchorless_v_allele_candidates() {
        let mut cfg = RefDataConfig::empty(ChainType::Vj);
        let _ = cfg.v_pool.push(Allele {
            name: "v_anchorless*01".into(),
            gene: "v_anchorless".into(),
            seq: b"AAACCCGGG".to_vec(),
            segment: Segment::V,
            anchor: None,
        });
        let anchored_v = cfg.v_pool.push(Allele {
            name: "v_anchored*01".into(),
            gene: "v_anchored".into(),
            seq: b"AAACCCGGG".to_vec(),
            segment: Segment::V,
            anchor: Some(6),
        });
        let _ = cfg.j_pool.push(Allele {
            name: "j*01".into(),
            gene: "j".into(),
            seq: b"TTTAAA".to_vec(),
            segment: Segment::J,
            anchor: Some(0),
        });
        let plan = vj_plan(&cfg);
        let contracts = crate::contract::productive();
        let compiled = CompiledSimulator::compile(
            &plan,
            Some(&cfg),
            Some(&contracts),
            ExecutionPolicy::Strict,
        )
        .expect("mixed anchored/anchorless support should compile");

        let outcome = compiled.run_one(0).expect("run should succeed");

        assert_eq!(
            outcome.final_simulation().assignments.v.unwrap().allele_id,
            anchored_v
        );
    }

    #[test]
    fn vj_productive_feasibility_filters_j_before_np_length_sampling() {
        let (cfg, j_good, j_bad) = vj_refdata_with_runtime_j_residue();
        let plan = vj_plan_with_np_lengths(&cfg, vec![(0, 1.0)]);
        let contracts = crate::contract::productive();
        let compiled = CompiledSimulator::compile(
            &plan,
            Some(&cfg),
            Some(&contracts),
            ExecutionPolicy::Strict,
        )
        .expect("mixed J support has one feasible completion");

        for seed in 0..32 {
            let outcome = compiled.run_one(seed).expect("run should stay feasible");
            let final_sim = outcome.final_simulation();
            assert_eq!(final_sim.assignments.j.unwrap().allele_id, j_good);
            assert_ne!(final_sim.assignments.j.unwrap().allele_id, j_bad);
            assert!(
                contracts.verify(final_sim, Some(&cfg)).is_ok(),
                "seed {seed} should satisfy productive()"
            );
        }
    }

    #[test]
    fn vj_productive_feasibility_respects_future_trim_support() {
        let mut cfg = RefDataConfig::empty(ChainType::Vj);
        let _ = cfg.v_pool.push(Allele {
            name: "v*01".into(),
            gene: "v".into(),
            seq: b"AAACCCGGGAAA".to_vec(),
            segment: Segment::V,
            anchor: Some(6),
        });
        let j_good = cfg.j_pool.push(Allele {
            name: "j_good*01".into(),
            gene: "j_good".into(),
            seq: b"TTTAAA".to_vec(),
            segment: Segment::J,
            anchor: Some(0),
        });
        let j_needs_future_v_trim = cfg.j_pool.push(Allele {
            name: "j_future_trim*01".into(),
            gene: "j_future_trim".into(),
            seq: b"TTTAAA".to_vec(),
            segment: Segment::J,
            anchor: Some(1),
        });

        let mut plan = PassPlan::new();
        plan.push(Box::new(SampleAllelePass::new(
            Segment::V,
            Box::new(AllelePoolDist::uniform(&cfg.v_pool)),
        )));
        plan.push(Box::new(SampleAllelePass::new(
            Segment::J,
            Box::new(AllelePoolDist::uniform(&cfg.j_pool)),
        )));
        plan.push(Box::new(TrimPass::new(
            Segment::V,
            TrimEnd::Three,
            Box::new(EmpiricalLengthDist::from_pairs(vec![(1, 1.0)])),
        )));
        plan.push(Box::new(AssembleSegmentPass::new(Segment::V)));
        plan.push(Box::new(GenerateNPPass::new(
            Segment::Np1,
            Box::new(EmpiricalLengthDist::from_pairs(vec![(0, 1.0)])),
            Box::new(UniformBase),
        )));
        plan.push(Box::new(AssembleSegmentPass::new(Segment::J)));

        let contracts = crate::contract::productive();
        let compiled = CompiledSimulator::compile(
            &plan,
            Some(&cfg),
            Some(&contracts),
            ExecutionPolicy::Strict,
        )
        .expect("future V trim support makes one J allele feasible");

        let outcome = compiled.run_one(0).expect("run should succeed");
        let final_sim = outcome.final_simulation();

        assert_eq!(
            final_sim.assignments.j.unwrap().allele_id,
            j_needs_future_v_trim
        );
        assert_ne!(final_sim.assignments.j.unwrap().allele_id, j_good);
        assert_eq!(final_sim.assignments.v.unwrap().trim_3, 1);
        contracts
            .verify(final_sim, Some(&cfg))
            .expect("final simulation should satisfy productive()");
    }

    #[test]
    fn vj_productive_feasibility_filters_upstream_known_stop() {
        let mut cfg = RefDataConfig::empty(ChainType::Vj);
        let _v_stop = cfg.v_pool.push(Allele {
            name: "v_stop*01".into(),
            gene: "v_stop".into(),
            seq: b"AAATAAGGG".to_vec(),
            segment: Segment::V,
            anchor: Some(3),
        });
        let v_productive = cfg.v_pool.push(Allele {
            name: "v_productive*01".into(),
            gene: "v_productive".into(),
            seq: b"AAATGTGGG".to_vec(),
            segment: Segment::V,
            anchor: Some(3),
        });
        let _ = cfg.j_pool.push(Allele {
            name: "j*01".into(),
            gene: "j".into(),
            seq: b"TGGAAA".to_vec(),
            segment: Segment::J,
            anchor: Some(0),
        });
        let plan = vj_plan_with_np_lengths(&cfg, vec![(0, 1.0)]);
        let contracts = crate::contract::productive();
        let compiled = CompiledSimulator::compile(
            &plan,
            Some(&cfg),
            Some(&contracts),
            ExecutionPolicy::Strict,
        )
        .expect("mixed V support has one stop-free productive completion");

        for seed in 0..32 {
            let outcome = compiled.run_one(seed).expect("run should stay feasible");
            let final_sim = outcome.final_simulation();
            assert_eq!(final_sim.assignments.v.unwrap().allele_id, v_productive);
            contracts
                .verify(final_sim, Some(&cfg))
                .unwrap_or_else(|_| panic!("seed {seed} should satisfy productive()"));
        }
    }

    #[test]
    fn compiled_run_one_matches_direct_runtime_for_valid_plan() {
        let cfg = vj_refdata();
        let plan = vj_plan(&cfg);
        let compiled =
            CompiledSimulator::compile(&plan, Some(&cfg), None, ExecutionPolicy::Permissive)
                .expect("plan should compile");

        let direct = PassRuntime::execute_with_refdata(&plan, Simulation::new(), 42, &cfg);
        let via_compiled = compiled.run_one(42).expect("run should succeed");

        assert_eq!(direct.trace.choices(), via_compiled.trace.choices());
        assert!(direct.events.is_empty());
        assert_eq!(via_compiled.events.len(), plan.len());
        assert_eq!(
            direct.final_simulation().pool.as_slice(),
            via_compiled.final_simulation().pool.as_slice()
        );
    }

    #[test]
    fn owned_compiled_simulator_owns_inputs_and_runs_repeatedly() {
        let cfg = vj_refdata();
        let plan = vj_plan(&cfg);
        let contracts = ContractSet::new().with(Box::new(MaxPoolLen::new(64)));
        let borrowed_report = {
            let borrowed = CompiledSimulator::compile(
                &plan,
                Some(&cfg),
                Some(&contracts),
                ExecutionPolicy::Strict,
            )
            .expect("borrowed simulator should compile");
            borrowed.report().clone()
        };

        let compiled = OwnedCompiledSimulator::compile(
            plan,
            Some(cfg.clone()),
            Some(contracts.clone()),
            ExecutionPolicy::Strict,
        )
        .expect("owned simulator should compile");

        let first = compiled.run_one(7).expect("first run should succeed");
        let second = compiled.run_one(7).expect("second run should succeed");

        assert_eq!(compiled.report(), &borrowed_report);
        assert_eq!(compiled.policy(), ExecutionPolicy::Strict);
        assert_eq!(compiled.plan().len(), 5);
        assert_eq!(compiled.refdata().unwrap().v_pool.len(), 1);
        assert_eq!(compiled.contracts().unwrap().len(), 1);
        assert_eq!(first.trace.choices(), second.trace.choices());
        assert_eq!(first.events.len(), compiled.plan().len());
        assert_eq!(second.events.len(), compiled.plan().len());
    }

    #[test]
    fn compiled_execution_uses_per_pass_trace_deltas() {
        for policy in [ExecutionPolicy::Permissive, ExecutionPolicy::Strict] {
            let mut plan = PassPlan::new();
            plan.push(Box::new(TraceProbePass::new(
                "trace_probe_one",
                "trace.probe.one",
                b'A',
            )));
            plan.push(Box::new(TraceProbePass::new(
                "trace_probe_two",
                "trace.probe.two",
                b'C',
            )));
            let compiled =
                CompiledSimulator::compile(&plan, None, None, policy).expect("plan should compile");

            let outcome = compiled.run_one(0).expect("run should succeed");

            assert_eq!(outcome.trace.len(), 2);
            assert_eq!(outcome.trace.choices()[0].address, "trace.probe.one");
            assert_eq!(outcome.trace.choices()[1].address, "trace.probe.two");
            assert_eq!(outcome.events.len(), 2);
            assert_eq!(outcome.events[0].pass_index, 0);
            assert_eq!(outcome.events[0].pass_name, "trace_probe_one");
            assert_eq!(outcome.events[0].trace_span, TraceSpan::new(0, 1));
            assert_eq!(outcome.events[0].pre.pool_len, 0);
            assert_eq!(outcome.events[0].post.pool_len, 1);
            assert_eq!(outcome.events[1].pass_index, 1);
            assert_eq!(outcome.events[1].pass_name, "trace_probe_two");
            assert_eq!(outcome.events[1].trace_span, TraceSpan::new(1, 2));
            assert_eq!(outcome.events[1].pre.pool_len, 1);
            assert_eq!(outcome.events[1].post.pool_len, 2);
            assert_eq!(
                outcome.pass_names,
                vec!["trace_probe_one", "trace_probe_two"]
            );
            assert_eq!(outcome.final_simulation().pool.len(), 2);
        }
    }

    #[test]
    fn strict_compiled_run_verifies_contracts_after_each_pass() {
        let mut plan = PassPlan::new();
        plan.push(Box::new(AppendBasePass::new("append_one", b'A')));
        plan.push(Box::new(AppendBasePass::new("append_two", b'C')));
        let contracts = ContractSet::new().with(Box::new(MaxPoolLen::new(1)));
        let compiled =
            CompiledSimulator::compile(&plan, None, Some(&contracts), ExecutionPolicy::Strict)
                .expect("plan should compile");

        let err = compiled
            .run_one(0)
            .expect_err("second pass should violate the strict contract fence");

        match err {
            PassError::ContractViolation {
                pass_name,
                violations,
            } => {
                assert_eq!(pass_name, "append_two");
                assert_eq!(violations.len(), 1);
                assert_eq!(violations[0].contract_name, "test.max_pool_len");
            }
            other => panic!("expected contract violation, got {other:?}"),
        }
    }

    #[test]
    fn owned_compiled_simulator_allows_runtime_policy_override() {
        let mut plan = PassPlan::new();
        plan.push(Box::new(AppendBasePass::new("append_one", b'A')));
        plan.push(Box::new(AppendBasePass::new("append_two", b'C')));
        let contracts = ContractSet::new().with(Box::new(MaxPoolLen::new(1)));
        let compiled = OwnedCompiledSimulator::compile(
            plan,
            None,
            Some(contracts),
            ExecutionPolicy::Permissive,
        )
        .expect("plan should compile");

        assert!(compiled.run_one(0).is_ok());
        let err = compiled
            .run_one_with_policy(0, ExecutionPolicy::Strict)
            .expect_err("strict override should enforce contract fence");

        assert!(matches!(err, PassError::ContractViolation { .. }));
    }

    #[test]
    fn strict_abort_does_not_commit_rejected_pass_trace_or_state() {
        let mut plan = PassPlan::new();
        plan.push(Box::new(TraceProbePass::new(
            "append_one",
            "trace.append.one",
            b'A',
        )));
        plan.push(Box::new(TraceProbePass::new(
            "append_two",
            "trace.append.two",
            b'C',
        )));
        let contracts = ContractSet::new().with(Box::new(MaxPoolLen::new(1)));
        let compiled =
            CompiledSimulator::compile(&plan, None, Some(&contracts), ExecutionPolicy::Strict)
                .expect("plan should compile");

        let abort = compiled
            .execute_transactional(Simulation::new(), 0)
            .expect_err("second pass should abort before commit");

        match abort.error {
            PassError::ContractViolation {
                pass_name,
                violations,
            } => {
                assert_eq!(pass_name, "append_two");
                assert_eq!(violations[0].contract_name, "test.max_pool_len");
            }
            other => panic!("expected contract violation, got {other:?}"),
        }

        assert_eq!(abort.committed.pass_names, vec!["append_one"]);
        assert_eq!(abort.committed.revisions.len(), 2);
        assert_eq!(abort.committed.final_simulation().pool.len(), 1);
        assert_eq!(abort.committed.trace.len(), 1);
        assert_eq!(abort.committed.events.len(), 1);
        assert_eq!(abort.committed.events[0].pass_name, "append_one");
        assert_eq!(abort.committed.events[0].trace_span, TraceSpan::new(0, 1));
        assert!(abort.committed.trace.find("trace.append.one").is_some());
        assert!(abort.committed.trace.find("trace.append.two").is_none());
    }

    // ──────────────────────────────────────────────────────────────
    // Base edit live-call refresh fixtures.
    //
    // Each fixture builds a tiny VDJ refdata, samples a chosen V
    // allele, assembles V, then runs a deterministic test-only base
    // edit pass to surface one of the four shrink / widen / switch /
    // unsupported behaviours the design doc calls out for base
    // edits.
    // ──────────────────────────────────────────────────────────────

    /// Test-only deterministic base distribution. Always samples the
    /// configured byte; reports it as the only point in `support()`.
    /// Used by mutation- and NP-pass fixtures to produce known base
    /// sequences without seed-pinning.
    #[derive(Clone, Debug)]
    struct ConstBaseDist(u8);

    impl crate::dist::Distribution for ConstBaseDist {
        type Output = u8;
        fn sample(&self, _rng: &mut crate::rng::Rng) -> u8 {
            self.0
        }
        fn support(&self) -> Option<Vec<(u8, f64)>> {
            Some(vec![(self.0, 1.0)])
        }
    }

    /// Test-only deterministic Pass that edits a single pool position
    /// to a chosen base and reports `PassEffect::EditBases`. Drives
    /// the live-call refresh path without any RNG involvement so the
    /// fixtures stay focused on the post-edit live-call shape.
    #[derive(Clone, Debug)]
    struct EditBaseAtPass {
        handle: crate::ir::NucHandle,
        new_base: u8,
    }

    impl EditBaseAtPass {
        fn new(pool_index: u32, new_base: u8) -> Self {
            Self {
                handle: crate::ir::NucHandle::new(pool_index),
                new_base,
            }
        }
    }

    impl Pass for EditBaseAtPass {
        fn name(&self) -> &str {
            "test.edit_base_at"
        }

        fn execute(
            &self,
            sim: &Simulation,
            _ctx: &mut crate::pass::PassContext,
        ) -> Simulation {
            sim.with_base_changed(self.handle, self.new_base)
        }

        fn effects(&self) -> Vec<PassEffect> {
            vec![PassEffect::EditBases]
        }
    }

    /// V-only refdata where every allele has a unique 9-base sequence
    /// (no shared prefix). Used by the widen / switch / unsupported
    /// fixtures where each ref position cleanly distinguishes alleles.
    fn distinct_v_refdata(seqs: &[(&str, &[u8])]) -> RefDataConfig {
        let mut cfg = RefDataConfig::empty(ChainType::Vdj);
        for (name, seq) in seqs {
            let _ = cfg.v_pool.push(Allele {
                name: (*name).to_string(),
                gene: name.split('*').next().unwrap_or(name).to_string(),
                seq: seq.to_vec(),
                segment: Segment::V,
                anchor: None,
            });
        }
        cfg
    }

    /// Drive a fixture: sample-allele(V) → assemble(V) → edit(pool, base).
    /// Returns the final `Outcome` and the V live call extracted from the
    /// final revision.
    fn run_edit(
        cfg: &RefDataConfig,
        sampled_id: AlleleId,
        edits: Vec<(u32, u8)>,
    ) -> (
        crate::pass::Outcome,
        crate::live_call::SegmentLiveCall,
    ) {
        let mut plan = PassPlan::new();
        plan.push(Box::new(SampleAllelePass::new(
            Segment::V,
            Box::new(AllelePoolDist::restricted_uniform(
                &cfg.v_pool,
                vec![sampled_id],
            )),
        )));
        plan.push(Box::new(AssembleSegmentPass::new(Segment::V)));
        for (pos, new_base) in edits {
            plan.push(Box::new(EditBaseAtPass::new(pos, new_base)));
        }

        let compiled =
            CompiledSimulator::compile(&plan, Some(cfg), None, ExecutionPolicy::Permissive)
                .expect("fixture plan should compile");
        let outcome = compiled.run_one(0).expect("plan should run");
        let final_sim = outcome.final_simulation().clone();
        let v_call = final_sim
            .live_calls
            .as_ref()
            .expect("live calls populated after assembly")
            .get(Segment::V)
            .cloned()
            .expect("V live call exists after assembly");
        (outcome, v_call)
    }

    #[test]
    fn mutation_widens_live_call_when_distinguishing_base_is_lost() {
        // Three alleles whose only difference is at position 3:
        //   A1: A A A G A A A A A
        //   A2: A A A C A A A A A
        //   A3: A A A C A A A A A   (identical to A2)
        //
        // Sampling A1 produces an assembled "AAAGAAAAA" → live call {A1}.
        // Mutating position 3 from G→C makes the assembled bases match
        // both A2 and A3 exactly → live call widens to {A2, A3}.
        let cfg = distinct_v_refdata(&[
            ("V*01", b"AAAGAAAAA"),
            ("V*02", b"AAACAAAAA"),
            ("V*03", b"AAACAAAAA"),
        ]);
        let v01 = AlleleId::new(0);
        let v02 = AlleleId::new(1);
        let v03 = AlleleId::new(2);

        let (_outcome, v_call) = run_edit(&cfg, v01, vec![(3, b'C')]);
        let mut expected = vec![v02, v03];
        expected.sort_by_key(|id| id.index());
        let mut actual = v_call.allele_call.to_ids();
        actual.sort_by_key(|id| id.index());
        assert_eq!(actual, expected, "expected widening to {{A2, A3}}");
        assert!(!v_call.allele_call.contains(v01), "A1 must drop out");
    }

    #[test]
    fn mutation_switches_live_call_to_a_different_singleton() {
        // Two alleles differing at position 3 only:
        //   A1: A A A G A A A A A   (sampled)
        //   A2: A A A C A A A A A
        //
        // Pre-edit live call = {A1}. After editing position 3 G→C, the
        // assembled bases match A2 exactly → live call = {A2}. Same
        // size, but the membership has switched.
        let cfg = distinct_v_refdata(&[
            ("V*01", b"AAAGAAAAA"),
            ("V*02", b"AAACAAAAA"),
        ]);
        let v01 = AlleleId::new(0);
        let v02 = AlleleId::new(1);

        let (_outcome, v_call) = run_edit(&cfg, v01, vec![(3, b'C')]);
        assert_eq!(v_call.allele_call.to_ids(), vec![v02]);
    }

    #[test]
    fn mutation_shrinks_live_call_from_three_to_one() {
        // Four alleles. The first three are identical (perfect ambiguity),
        // the fourth differs at position 2:
        //   A1: A A A C C C   (sampled — identical to A2 / A3)
        //   A2: A A A C C C
        //   A3: A A A C C C
        //   A4: A A T C C C
        //
        // Pre-edit assembled "AAACCC" matches A1, A2, A3 but not A4 →
        // live call {A1, A2, A3} (size 3). After editing position 2
        // A→T the assembled becomes "AATCCC" which only matches A4 →
        // live call shrinks to {A4} (size 1).
        let cfg = distinct_v_refdata(&[
            ("V*01", b"AAACCC"),
            ("V*02", b"AAACCC"),
            ("V*03", b"AAACCC"),
            ("V*04", b"AATCCC"),
        ]);
        let v01 = AlleleId::new(0);
        let v04 = AlleleId::new(3);

        let (_outcome, v_call) = run_edit(&cfg, v01, vec![(2, b'T')]);
        assert_eq!(
            v_call.allele_call.len(),
            1,
            "expected post-edit set size 1, got {}",
            v_call.allele_call.len()
        );
        assert_eq!(v_call.allele_call.to_ids(), vec![v04]);
    }

    #[test]
    fn mutation_to_orphan_base_keeps_truth_in_tie_set() {
        // Single allele with a fully-determined sequence. Mutating one
        // assembled base to a value no allele has at that ref position
        // simply leaves that position non-informative (its evidence
        // contributes zero to every allele's score). The remaining
        // five positions still score the truth allele uniquely highest,
        // so the tie-set keeps V*01 — exactly the evidence-resilient
        // behavior the score-and-tie caller is designed to deliver.
        let cfg = distinct_v_refdata(&[("V*01", b"AAAAAA")]);
        let v01 = AlleleId::new(0);

        let (_outcome, v_call) = run_edit(&cfg, v01, vec![(2, b'T')]);
        assert_eq!(
            v_call.allele_call.to_ids(),
            vec![v01],
            "truth allele should remain at max score (5/6 positions still match)"
        );
    }

    #[test]
    fn pre_edit_live_call_visible_in_intermediate_revision() {
        // The compiled runtime stores a revision after every committed
        // pass. The pre-edit live call should be reachable via
        // revisions[2] (after assemble) while the final revision shows
        // the post-edit (widened) call. This guards against accidentally
        // rewriting earlier revisions when refreshing live calls.
        let cfg = distinct_v_refdata(&[
            ("V*01", b"AAAGAAAAA"),
            ("V*02", b"AAACAAAAA"),
            ("V*03", b"AAACAAAAA"),
        ]);
        let v01 = AlleleId::new(0);

        let (outcome, post_edit_call) = run_edit(&cfg, v01, vec![(3, b'C')]);

        // revisions: [initial, post-sample, post-assemble, post-edit].
        assert_eq!(outcome.revisions.len(), 4);
        let after_assemble = &outcome.revisions[2];
        let pre_edit_call = after_assemble
            .live_calls
            .as_ref()
            .expect("live calls populated after assembly")
            .get(Segment::V)
            .cloned()
            .expect("V live call exists after assembly");

        assert_eq!(pre_edit_call.allele_call.to_ids(), vec![v01]);
        assert_eq!(post_edit_call.allele_call.len(), 2);
        assert!(!post_edit_call.allele_call.contains(v01));
    }

    #[test]
    fn real_uniform_mutation_pass_increments_live_call_version() {
        // Integration-style coverage: confirm the live-call refresh
        // path also fires when the real `UniformMutationPass` runs,
        // not just our synthetic `EditBaseAtPass`. UniformMutationPass
        // samples positions with replacement so the post-state isn't
        // deterministic — but the refresh hook IS, and we can assert
        // that against the live-call `evidence_version` counter:
        // every successful refresh bumps it, so the post-mutation
        // version must be strictly greater than the post-assemble
        // version.
        use crate::passes::UniformMutationPass;

        let cfg = distinct_v_refdata(&[("V*01", b"AAAAAAAAA")]);
        let v01 = AlleleId::new(0);

        let mut plan = PassPlan::new();
        plan.push(Box::new(SampleAllelePass::new(
            Segment::V,
            Box::new(AllelePoolDist::restricted_uniform(
                &cfg.v_pool,
                vec![v01],
            )),
        )));
        plan.push(Box::new(AssembleSegmentPass::new(Segment::V)));
        plan.push(Box::new(UniformMutationPass::new(
            Box::new(EmpiricalLengthDist::from_pairs(vec![(3, 1.0)])),
            Box::new(ConstBaseDist(b'A')),
        )));

        let compiled =
            CompiledSimulator::compile(&plan, Some(&cfg), None, ExecutionPolicy::Permissive)
                .expect("integration plan should compile");
        let outcome = compiled.run_one(0).expect("plan should run");

        // revisions: [initial, post-sample, post-assemble, post-mutate].
        let post_assemble_version = outcome.revisions[2]
            .live_calls
            .as_ref()
            .expect("live calls populated after assembly")
            .version;
        let post_mutate_version = outcome.revisions[3]
            .live_calls
            .as_ref()
            .expect("live calls populated after mutation")
            .version;
        assert!(
            post_mutate_version > post_assemble_version,
            "EditBases pass must bump live-call evidence_version \
             ({post_mutate_version} should be > {post_assemble_version})"
        );
    }

    #[test]
    fn edit_outside_assembled_segment_leaves_live_call_unchanged() {
        // Sanity: editing a pool position that is NOT inside V's
        // assembled region must not invent extra hypotheses or break
        // the existing one. The test allele has length 6; we run two
        // edit passes — the first edits position 0 (inside V) so we
        // know recompute is happening, the second edits position 5
        // (still inside V's region for this fixture, so live call
        // stays consistent with the new bases at every step).
        //
        // We use this fixture to confirm that two consecutive
        // EditBases passes both trigger refresh and the final state
        // matches a manual recomputation from the post-edit bases.
        let cfg = distinct_v_refdata(&[
            ("V*01", b"AAAAAA"),
            ("V*02", b"AAAAAA"),
        ]);
        let v01 = AlleleId::new(0);
        let v02 = AlleleId::new(1);

        // Both alleles are identical; the call remains {V*01, V*02}
        // through every base edit because every position has the same
        // base across alleles.
        let (_outcome, v_call) = run_edit(
            &cfg,
            v01,
            vec![(0, b'A'), (5, b'A')], // no-ops in terms of base change
        );
        let mut expected = vec![v01, v02];
        expected.sort_by_key(|id| id.index());
        let mut actual = v_call.allele_call.to_ids();
        actual.sort_by_key(|id| id.index());
        assert_eq!(actual, expected);
    }

    // ──────────────────────────────────────────────────────────────
    // V right-boundary extension into NP1.
    //
    // The acceptance criterion from the design doc: a trimmed V suffix
    // can be recreated by NP1 bases, and when that happens the live
    // `v_call` should shrink to the allele(s) that the recreated suffix
    // exactly extends.
    // ──────────────────────────────────────────────────────────────

    /// Build a VDJ refdata holding two V alleles that share a 9-base
    /// prefix and differ in their 3-base suffix; the J allele is a
    /// minimal stub (the fixture stops after NP1 generation, so D / J
    /// pools are unused at runtime). Returns the refdata and the two
    /// V allele ids in declaration order.
    fn v_extension_refdata() -> (RefDataConfig, AlleleId, AlleleId) {
        let mut cfg = RefDataConfig::empty(ChainType::Vdj);
        let v01 = cfg.v_pool.push(Allele {
            name: "V*01".into(),
            gene: "V".into(),
            seq: b"AAACCCGGGTTT".to_vec(),
            segment: Segment::V,
            anchor: None,
        });
        let v02 = cfg.v_pool.push(Allele {
            name: "V*02".into(),
            gene: "V".into(),
            seq: b"AAACCCGGGAAA".to_vec(),
            segment: Segment::V,
            anchor: None,
        });
        // Minimal J entry — required by ChainType::Vdj construction
        // (refdata builder does not enforce this) but never executed
        // because the plan stops after NP1 generation.
        let _ = cfg.j_pool.push(Allele {
            name: "J*01".into(),
            gene: "J".into(),
            seq: b"AAA".to_vec(),
            segment: Segment::J,
            anchor: None,
        });
        (cfg, v01, v02)
    }

    /// Build the standard V-extension plan:
    ///   sample(V) → trim(V_3, by) → assemble(V) → generate(NP1, len, base)
    /// `np_len` and `np_base` control the NP1 region's content
    /// deterministically.
    fn v_extension_plan(
        cfg: &RefDataConfig,
        sampled_v: AlleleId,
        v_trim_3: i64,
        np_len: i64,
        np_base: u8,
    ) -> PassPlan {
        let mut plan = PassPlan::new();
        plan.push(Box::new(SampleAllelePass::new(
            Segment::V,
            Box::new(AllelePoolDist::restricted_uniform(
                &cfg.v_pool,
                vec![sampled_v],
            )),
        )));
        plan.push(Box::new(TrimPass::new(
            Segment::V,
            TrimEnd::Three,
            Box::new(EmpiricalLengthDist::from_pairs(vec![(v_trim_3, 1.0)])),
        )));
        plan.push(Box::new(AssembleSegmentPass::new(Segment::V)));
        plan.push(Box::new(GenerateNPPass::new(
            Segment::Np1,
            Box::new(EmpiricalLengthDist::from_pairs(vec![(np_len, 1.0)])),
            Box::new(ConstBaseDist(np_base)),
        )));
        plan
    }

    #[test]
    fn v_call_shrinks_when_np1_recreates_trimmed_suffix() {
        // V*01 = AAACCCGGG TTT (suffix TTT distinguishes it).
        // V*02 = AAACCCGGG AAA (suffix AAA distinguishes it).
        // Sample V*01, trim 3' by 3 → assembled V is AAACCCGGG.
        // Both alleles match the assembled bases at every position →
        // post-assemble v_call = {V*01, V*02}.
        // Then GenerateNP(NP1, length=3, 'T') puts TTT right after V.
        // V right-extension into NP1 walks T → matches V*01 ref pos 9
        // (V*01[9]='T'); V*02[9]='A' → V*02 drops out. After 3 NP1
        // bases the call has shrunk back to {V*01}.
        let (cfg, v01, _v02) = v_extension_refdata();

        let plan = v_extension_plan(&cfg, v01, 3, 3, b'T');
        let compiled =
            CompiledSimulator::compile(&plan, Some(&cfg), None, ExecutionPolicy::Permissive)
                .expect("plan should compile");
        let outcome = compiled.run_one(0).expect("plan should run");

        // revisions: [initial, post-sample, post-trim, post-assemble,
        //             post-np1].
        assert_eq!(outcome.revisions.len(), 5);

        // After AssembleSegment(V) but before NP1: both alleles match.
        let post_assemble_call = outcome.revisions[3]
            .live_calls
            .as_ref()
            .expect("live calls populated after assembly")
            .get(Segment::V)
            .cloned()
            .expect("V live call exists after assembly");
        let mut post_assemble_ids = post_assemble_call.allele_call.to_ids();
        post_assemble_ids.sort_by_key(|id| id.index());
        assert_eq!(
            post_assemble_ids,
            vec![v01, AlleleId::new(1)],
            "post-assemble v_call should hold both indistinguishable V alleles"
        );

        // After NP1 generation: the right-extension walk has
        // narrowed to V*01 only.
        let post_np1_call = outcome.revisions[4]
            .live_calls
            .as_ref()
            .expect("live calls still populated after NP1")
            .get(Segment::V)
            .cloned()
            .expect("V live call exists after NP1");
        assert_eq!(post_np1_call.allele_call.to_ids(), vec![v01]);

        // The hypothesis must record the elastic boundary so downstream
        // tooling can tell that the right end was extended past the
        // structural V region.
        let elastic_hypothesis = post_np1_call
            .hypotheses
            .iter()
            .find(|h| {
                h.flags
                    .contains(crate::live_call::HypothesisFlags::BOUNDARY_ELASTIC)
            })
            .expect("at least one hypothesis should be flagged BOUNDARY_ELASTIC");
        // ref_end advanced past the (post-trim) V region's ref_end of 9
        // by exactly 3 (TTT extension covers ref positions 9..12).
        assert_eq!(
            elastic_hypothesis.ref_end, 12,
            "ref_end should reach allele length 12 after extending into NP1"
        );
        assert_eq!(
            elastic_hypothesis.seq_end,
            outcome
                .final_simulation()
                .pool
                .len() as u32,
            "seq_end should reach the end of the pool (V region + 3 NP1 bases)"
        );
    }

    #[test]
    fn v_call_stays_widened_when_np1_does_not_match_any_allele() {
        // Same fixture as above but NP1 emits 'C', which is NOT the
        // next ref base for either allele (V*01[9]='T', V*02[9]='A').
        // The extension walk halts immediately — the live call should
        // remain the post-assemble widened {V*01, V*02}.
        let (cfg, v01, v02) = v_extension_refdata();

        let plan = v_extension_plan(&cfg, v01, 3, 3, b'C');
        let compiled =
            CompiledSimulator::compile(&plan, Some(&cfg), None, ExecutionPolicy::Permissive)
                .expect("plan should compile");
        let outcome = compiled.run_one(0).expect("plan should run");

        let post_np1_call = outcome
            .final_simulation()
            .live_calls
            .as_ref()
            .expect("live calls populated")
            .get(Segment::V)
            .cloned()
            .expect("V live call exists");
        let mut ids = post_np1_call.allele_call.to_ids();
        ids.sort_by_key(|id| id.index());
        assert_eq!(ids, vec![v01, v02]);

        // No extension → no BOUNDARY_ELASTIC flag.
        for h in &post_np1_call.hypotheses {
            assert!(
                !h.flags
                    .contains(crate::live_call::HypothesisFlags::BOUNDARY_ELASTIC),
                "no NP1 base extended; BOUNDARY_ELASTIC must NOT be set"
            );
        }
    }

    #[test]
    fn v_call_partially_extends_when_np1_matches_only_a_prefix() {
        // V*01 = AAACCCGGG TTT, V*02 = AAACCCGGG AAA. Trim V_3 by 3 →
        // assembled V is AAACCCGGG. NP1 emits 'T' for length 5 → bases
        // T-T-T-?-?. The first 3 NP1 bases match V*01's suffix exactly
        // (positions 9, 10, 11 are all 'T'). Position 12 doesn't exist
        // in V*01 (ref length is 12) so the extension halts at
        // ref_pos=12 / seq_end = (V end) + 3 = 12. The remaining NP1
        // bases stay outside the V hypothesis.
        let (cfg, v01, _v02) = v_extension_refdata();

        let plan = v_extension_plan(&cfg, v01, 3, 5, b'T');
        let compiled =
            CompiledSimulator::compile(&plan, Some(&cfg), None, ExecutionPolicy::Permissive)
                .expect("plan should compile");
        let outcome = compiled.run_one(0).expect("plan should run");

        let post_np1_call = outcome
            .final_simulation()
            .live_calls
            .as_ref()
            .expect("live calls populated")
            .get(Segment::V)
            .cloned()
            .expect("V live call exists");
        assert_eq!(post_np1_call.allele_call.to_ids(), vec![v01]);

        let elastic_hypothesis = post_np1_call
            .hypotheses
            .iter()
            .find(|h| {
                h.flags
                    .contains(crate::live_call::HypothesisFlags::BOUNDARY_ELASTIC)
            })
            .expect("hypothesis should be flagged BOUNDARY_ELASTIC");
        assert_eq!(
            elastic_hypothesis.ref_end, 12,
            "extension should stop when ref position reaches V*01 allele length 12"
        );
        // seq_end is V_region.end (9) + 3 extension bases = 12. The
        // remaining 2 NP1 bases (seq positions 12 and 13) are outside
        // the V hypothesis.
        assert_eq!(elastic_hypothesis.seq_end, 12);
    }

    #[test]
    fn v_call_extension_no_op_when_no_trim() {
        // Sanity: with no V_3 trim, V's region already covers the
        // whole allele. There is no ref position past V's end, so the
        // extension walk has nothing to do — the call should remain a
        // singleton {sampled} and BOUNDARY_ELASTIC should NOT be set.
        let (cfg, v01, _v02) = v_extension_refdata();

        let plan = v_extension_plan(&cfg, v01, 0, 3, b'T');
        let compiled =
            CompiledSimulator::compile(&plan, Some(&cfg), None, ExecutionPolicy::Permissive)
                .expect("plan should compile");
        let outcome = compiled.run_one(0).expect("plan should run");

        let post_np1_call = outcome
            .final_simulation()
            .live_calls
            .as_ref()
            .expect("live calls populated")
            .get(Segment::V)
            .cloned()
            .expect("V live call exists");
        assert_eq!(post_np1_call.allele_call.to_ids(), vec![v01]);
        for h in &post_np1_call.hypotheses {
            assert!(
                !h.flags
                    .contains(crate::live_call::HypothesisFlags::BOUNDARY_ELASTIC),
                "no extension should occur when V already covers full allele"
            );
        }
    }

    #[test]
    fn append_region_np1_bumps_live_call_version() {
        // Plumbing check: every successful refresh bumps the live-call
        // evidence_version. AppendRegion(Np1) must trigger a V refresh
        // and therefore advance the version past the post-assemble
        // value. (Used to guard against accidentally removing the
        // hook from `apply_live_call_updates`.)
        let (cfg, v01, _v02) = v_extension_refdata();
        let plan = v_extension_plan(&cfg, v01, 3, 3, b'T');

        let compiled =
            CompiledSimulator::compile(&plan, Some(&cfg), None, ExecutionPolicy::Permissive)
                .expect("plan should compile");
        let outcome = compiled.run_one(0).expect("plan should run");

        let post_assemble_version = outcome.revisions[3]
            .live_calls
            .as_ref()
            .expect("post-assemble live calls present")
            .version;
        let post_np1_version = outcome.revisions[4]
            .live_calls
            .as_ref()
            .expect("post-np1 live calls present")
            .version;
        assert!(
            post_np1_version > post_assemble_version,
            "AppendRegion(Np1) must bump live-call version: \
             {post_np1_version} should be > {post_assemble_version}"
        );
    }

    // ──────────────────────────────────────────────────────────────
    // J left-boundary extension into the immediately-
    // preceding NP region (NP1 in VJ chains, NP2 in VDJ chains).
    //
    // Mirrors the V right-boundary case but for the left side of J:
    // the chosen NP bases must extend backward into the J allele's
    // reference prefix.
    // Acceptance: trim J 5' so live j_call widens, then have NP bases
    // happen to recreate the trimmed J prefix → j_call shrinks back.
    // ──────────────────────────────────────────────────────────────

    /// Build a VJ-chain refdata with two J alleles that share a 9-base
    /// suffix and differ in their 3-base 5' prefix. V is a minimal
    /// stub used only so the J assembly pass has a chain context.
    fn j_extension_vj_refdata() -> (RefDataConfig, AlleleId, AlleleId) {
        let mut cfg = RefDataConfig::empty(ChainType::Vj);
        // Stub V — never sampled in the fixture plan.
        let _ = cfg.v_pool.push(Allele {
            name: "V*01".into(),
            gene: "V".into(),
            seq: b"AAA".to_vec(),
            segment: Segment::V,
            anchor: None,
        });
        // J*01 has prefix "TTT" + shared "ACGTACGTA". J*02 has
        // prefix "GGG" + the same shared suffix. With J_5 trim = 3,
        // the assembled J region exposes only the shared "ACGTACGTA"
        // → both alleles match equally → live j_call widens to {J*01, J*02}.
        let j01 = cfg.j_pool.push(Allele {
            name: "J*01".into(),
            gene: "J".into(),
            seq: b"TTTACGTACGTA".to_vec(),
            segment: Segment::J,
            anchor: None,
        });
        let j02 = cfg.j_pool.push(Allele {
            name: "J*02".into(),
            gene: "J".into(),
            seq: b"GGGACGTACGTA".to_vec(),
            segment: Segment::J,
            anchor: None,
        });
        (cfg, j01, j02)
    }

    /// Plan that:
    ///   sample(V, single allele) → assemble(V) →
    ///     generate(NP1, len, base) → sample(J, single allele) →
    ///     trim(J_5, by) → assemble(J)
    ///
    /// The V is a 3-base stub so it doesn't dominate test reasoning.
    /// `np_len` and `np_base` control NP1 contents deterministically;
    /// `j_trim_5` controls how many J prefix bases are stripped before
    /// J is assembled.
    fn j_extension_plan(
        cfg: &RefDataConfig,
        sampled_j: AlleleId,
        j_trim_5: i64,
        np_len: i64,
        np_base: u8,
    ) -> PassPlan {
        let v_id = AlleleId::new(0);
        let mut plan = PassPlan::new();
        plan.push(Box::new(SampleAllelePass::new(
            Segment::V,
            Box::new(AllelePoolDist::restricted_uniform(
                &cfg.v_pool,
                vec![v_id],
            )),
        )));
        plan.push(Box::new(AssembleSegmentPass::new(Segment::V)));
        plan.push(Box::new(GenerateNPPass::new(
            Segment::Np1,
            Box::new(EmpiricalLengthDist::from_pairs(vec![(np_len, 1.0)])),
            Box::new(ConstBaseDist(np_base)),
        )));
        plan.push(Box::new(SampleAllelePass::new(
            Segment::J,
            Box::new(AllelePoolDist::restricted_uniform(
                &cfg.j_pool,
                vec![sampled_j],
            )),
        )));
        plan.push(Box::new(TrimPass::new(
            Segment::J,
            TrimEnd::Five,
            Box::new(EmpiricalLengthDist::from_pairs(vec![(j_trim_5, 1.0)])),
        )));
        plan.push(Box::new(AssembleSegmentPass::new(Segment::J)));
        plan
    }

    #[test]
    fn j_call_shrinks_when_np1_recreates_trimmed_prefix_vj() {
        // J*01 = TTT ACGTACGTA, J*02 = GGG ACGTACGTA.
        // Sample J*01, trim 5' by 3 → assembled J ref window starts at
        // pos 3 (covering ACGTACGTA).
        // Both alleles match the assembled J bases → post-assemble
        // j_call = {J*01, J*02}.
        // NP1 bases = TTT (length 3, all 'T'). The walker checks NP1's
        // RIGHTMOST base first against ref pos 2 (J*01[2]='T'). Then
        // pos 1 (J*01[1]='T'). Then pos 0 (J*01[0]='T'). All match
        // J*01 only → j_call shrinks back to {J*01}.
        let (cfg, j01, _j02) = j_extension_vj_refdata();

        let plan = j_extension_plan(&cfg, j01, 3, 3, b'T');
        let compiled =
            CompiledSimulator::compile(&plan, Some(&cfg), None, ExecutionPolicy::Permissive)
                .expect("plan should compile");
        let outcome = compiled.run_one(0).expect("plan should run");

        // J refresh fires on AssembleSegment(J) which is the last pass.
        // The final revision should hold the post-extension j_call.
        let final_call = outcome
            .final_simulation()
            .live_calls
            .as_ref()
            .expect("live calls populated")
            .get(Segment::J)
            .cloned()
            .expect("J live call exists after assembly");
        assert_eq!(final_call.allele_call.to_ids(), vec![j01]);

        let elastic = final_call
            .hypotheses
            .iter()
            .find(|h| {
                h.flags
                    .contains(crate::live_call::HypothesisFlags::BOUNDARY_ELASTIC)
            })
            .expect("J hypothesis should be flagged BOUNDARY_ELASTIC");
        // The extension brought ref_start from 3 (post-trim) back to 0
        // (matching all 3 NP1 bases against J*01[0..3]).
        assert_eq!(elastic.ref_start, 0);
    }

    #[test]
    fn j_call_stays_widened_when_np1_does_not_recreate_prefix_vj() {
        // Same fixture as above, but NP1 emits 'C' which is not the
        // expected base for either allele's prefix at any position.
        // The extension halts immediately; j_call stays widened.
        let (cfg, j01, j02) = j_extension_vj_refdata();

        let plan = j_extension_plan(&cfg, j01, 3, 3, b'C');
        let compiled =
            CompiledSimulator::compile(&plan, Some(&cfg), None, ExecutionPolicy::Permissive)
                .expect("plan should compile");
        let outcome = compiled.run_one(0).expect("plan should run");

        let final_call = outcome
            .final_simulation()
            .live_calls
            .as_ref()
            .expect("live calls populated")
            .get(Segment::J)
            .cloned()
            .expect("J live call exists");
        let mut ids = final_call.allele_call.to_ids();
        ids.sort_by_key(|id| id.index());
        assert_eq!(ids, vec![j01, j02]);
        for h in &final_call.hypotheses {
            assert!(
                !h.flags
                    .contains(crate::live_call::HypothesisFlags::BOUNDARY_ELASTIC),
                "no extension occurred → BOUNDARY_ELASTIC must not be set"
            );
        }
    }

    #[test]
    fn j_call_partially_extends_when_np1_matches_only_a_suffix_of_prefix() {
        // J*01 = TTT ACGTACGTA, J*02 = GGG ACGTACGTA.
        // Trim J_5 by 3, NP1 length 5, base 'T'. NP1 = T-T-T-T-T.
        // Walker checks rightmost NP1 byte vs ref_start - 1 = 2 first.
        //   NP1[4]='T' vs J*01[2]='T' ✓ → candidates={J*01}.
        //   NP1[3]='T' vs J*01[1]='T' ✓.
        //   NP1[2]='T' vs J*01[0]='T' ✓.
        //   NP1[1]: ref_start now 0, can't go further → halt.
        // Two of the five NP1 bases stay outside the J hypothesis.
        let (cfg, j01, _j02) = j_extension_vj_refdata();

        let plan = j_extension_plan(&cfg, j01, 3, 5, b'T');
        let compiled =
            CompiledSimulator::compile(&plan, Some(&cfg), None, ExecutionPolicy::Permissive)
                .expect("plan should compile");
        let outcome = compiled.run_one(0).expect("plan should run");

        let final_call = outcome
            .final_simulation()
            .live_calls
            .as_ref()
            .expect("live calls populated")
            .get(Segment::J)
            .cloned()
            .expect("J live call exists");
        assert_eq!(final_call.allele_call.to_ids(), vec![j01]);

        let elastic = final_call
            .hypotheses
            .iter()
            .find(|h| {
                h.flags
                    .contains(crate::live_call::HypothesisFlags::BOUNDARY_ELASTIC)
            })
            .expect("J hypothesis should be flagged BOUNDARY_ELASTIC");
        assert_eq!(
            elastic.ref_start, 0,
            "extension reached ref_start = 0 (the start of the allele)"
        );
        // J's seq_start moved left by 3 (3 NP1 bases consumed). The
        // remaining 2 NP1 bases (positions 0 and 1 in the pool's NP1
        // span) sit outside the J hypothesis.
        // Pool layout: [V (3 bases) 0..3] [NP1 (5 bases) 3..8] [J (9
        // bases) 8..17]. Post-extension J seq_start = 8 - 3 = 5.
        assert_eq!(elastic.seq_start, 5);
        assert_eq!(elastic.seq_end, 17);
    }

    #[test]
    fn j_call_extension_no_op_when_no_trim() {
        // No J_5 trim → J's ref_start is already 0 → there is no ref
        // position to the left to extend into. The walker halts on
        // the very first iteration because `extended_ref_start == 0`.
        let (cfg, j01, _j02) = j_extension_vj_refdata();

        let plan = j_extension_plan(&cfg, j01, 0, 3, b'T');
        let compiled =
            CompiledSimulator::compile(&plan, Some(&cfg), None, ExecutionPolicy::Permissive)
                .expect("plan should compile");
        let outcome = compiled.run_one(0).expect("plan should run");

        let final_call = outcome
            .final_simulation()
            .live_calls
            .as_ref()
            .expect("live calls populated")
            .get(Segment::J)
            .cloned()
            .expect("J live call exists");
        assert_eq!(final_call.allele_call.to_ids(), vec![j01]);
        for h in &final_call.hypotheses {
            assert!(
                !h.flags
                    .contains(crate::live_call::HypothesisFlags::BOUNDARY_ELASTIC),
                "no extension when J's ref_start is already 0"
            );
        }
    }

    #[test]
    fn j_left_extension_works_for_vdj_chain_via_np2() {
        // Verify the chain-agnostic neighbour lookup: in a VDJ chain,
        // J's left neighbour is NP2 (not NP1). Build a minimal VDJ
        // refdata, plan V→NP1→D→NP2→trim(J_5)→J, and check that NP2
        // bases drive the J left-extension.
        let mut cfg = RefDataConfig::empty(ChainType::Vdj);
        let v_id = cfg.v_pool.push(Allele {
            name: "V*01".into(),
            gene: "V".into(),
            seq: b"AAA".to_vec(),
            segment: Segment::V,
            anchor: None,
        });
        let d_id = cfg.d_pool.push(Allele {
            name: "D*01".into(),
            gene: "D".into(),
            seq: b"CCC".to_vec(),
            segment: Segment::D,
            anchor: None,
        });
        let j01 = cfg.j_pool.push(Allele {
            name: "J*01".into(),
            gene: "J".into(),
            seq: b"TTTACGTACGTA".to_vec(),
            segment: Segment::J,
            anchor: None,
        });
        let _j02 = cfg.j_pool.push(Allele {
            name: "J*02".into(),
            gene: "J".into(),
            seq: b"GGGACGTACGTA".to_vec(),
            segment: Segment::J,
            anchor: None,
        });

        let mut plan = PassPlan::new();
        plan.push(Box::new(SampleAllelePass::new(
            Segment::V,
            Box::new(AllelePoolDist::restricted_uniform(
                &cfg.v_pool,
                vec![v_id],
            )),
        )));
        plan.push(Box::new(SampleAllelePass::new(
            Segment::D,
            Box::new(AllelePoolDist::restricted_uniform(
                &cfg.d_pool,
                vec![d_id],
            )),
        )));
        plan.push(Box::new(AssembleSegmentPass::new(Segment::V)));
        plan.push(Box::new(GenerateNPPass::new(
            Segment::Np1,
            Box::new(EmpiricalLengthDist::from_pairs(vec![(0, 1.0)])),
            Box::new(ConstBaseDist(b'A')),
        )));
        plan.push(Box::new(AssembleSegmentPass::new(Segment::D)));
        plan.push(Box::new(GenerateNPPass::new(
            Segment::Np2,
            Box::new(EmpiricalLengthDist::from_pairs(vec![(3, 1.0)])),
            Box::new(ConstBaseDist(b'T')),
        )));
        plan.push(Box::new(SampleAllelePass::new(
            Segment::J,
            Box::new(AllelePoolDist::restricted_uniform(
                &cfg.j_pool,
                vec![j01],
            )),
        )));
        plan.push(Box::new(TrimPass::new(
            Segment::J,
            TrimEnd::Five,
            Box::new(EmpiricalLengthDist::from_pairs(vec![(3, 1.0)])),
        )));
        plan.push(Box::new(AssembleSegmentPass::new(Segment::J)));

        let compiled =
            CompiledSimulator::compile(&plan, Some(&cfg), None, ExecutionPolicy::Permissive)
                .expect("VDJ plan should compile");
        let outcome = compiled.run_one(0).expect("plan should run");
        let final_call = outcome
            .final_simulation()
            .live_calls
            .as_ref()
            .expect("live calls populated")
            .get(Segment::J)
            .cloned()
            .expect("J live call exists");
        assert_eq!(
            final_call.allele_call.to_ids(),
            vec![j01],
            "VDJ chain: NP2 bases TTT should narrow J back to J*01"
        );
        let elastic = final_call
            .hypotheses
            .iter()
            .find(|h| {
                h.flags
                    .contains(crate::live_call::HypothesisFlags::BOUNDARY_ELASTIC)
            })
            .expect("J hypothesis should be flagged BOUNDARY_ELASTIC");
        assert_eq!(elastic.ref_start, 0);
    }

    // ──────────────────────────────────────────────────────────────
    // D both-side extension into NP1 / NP2.
    //
    // D is the "hardest segment" per the design doc — it has elastic
    // boundaries on both sides. We test both directions independently
    // and then together. Pipeline order in VDJ is:
    //
    //   sample(V) → sample(D) → assemble(V) → np1 → assemble(D) → np2 → ...
    //
    // - D LEFT extension fires automatically on AssembleSegment(D)
    //   because NP1 already exists by then.
    // - D RIGHT extension needs the new AppendRegion(Np2) hook
    //   (NP2 is generated AFTER D is assembled).
    // ──────────────────────────────────────────────────────────────

    /// Build a VDJ refdata with two D alleles that share a 5-base core
    /// "CCCCC" and differ in their 3-base 5' prefix and 3-base 3'
    /// suffix.
    ///
    ///   D*01 = TTT CCCCC AAA  (unique prefix TTT, unique suffix AAA)
    ///   D*02 = GGG CCCCC TGT  (different prefix GGG, different suffix TGT)
    ///
    /// V and J pools are minimal stubs — needed so the chain can
    /// assemble end-to-end but irrelevant to the D-call assertions.
    fn d_extension_refdata() -> (RefDataConfig, AlleleId, AlleleId) {
        let mut cfg = RefDataConfig::empty(ChainType::Vdj);
        let _ = cfg.v_pool.push(Allele {
            name: "V*01".into(),
            gene: "V".into(),
            seq: b"AAA".to_vec(),
            segment: Segment::V,
            anchor: None,
        });
        let d01 = cfg.d_pool.push(Allele {
            name: "D*01".into(),
            gene: "D".into(),
            seq: b"TTTCCCCCAAA".to_vec(),
            segment: Segment::D,
            anchor: None,
        });
        let d02 = cfg.d_pool.push(Allele {
            name: "D*02".into(),
            gene: "D".into(),
            seq: b"GGGCCCCCTGT".to_vec(),
            segment: Segment::D,
            anchor: None,
        });
        let _ = cfg.j_pool.push(Allele {
            name: "J*01".into(),
            gene: "J".into(),
            seq: b"AAA".to_vec(),
            segment: Segment::J,
            anchor: None,
        });
        (cfg, d01, d02)
    }

    /// Build the standard D-extension plan:
    ///
    ///   sample(V, stub) → sample(D, sampled_d) → assemble(V) →
    ///   np1(np1_len, np1_base) → trim(D_5, by) → trim(D_3, by) →
    ///   assemble(D) → np2(np2_len, np2_base) → sample(J, stub) →
    ///   assemble(J)
    ///
    /// V and J are stubs; the assertions only inspect the D live call.
    fn d_extension_plan(
        cfg: &RefDataConfig,
        sampled_d: AlleleId,
        d_trim_5: i64,
        d_trim_3: i64,
        np1_len: i64,
        np1_base: u8,
        np2_len: i64,
        np2_base: u8,
    ) -> PassPlan {
        let v_id = AlleleId::new(0);
        let j_id = AlleleId::new(0);
        let mut plan = PassPlan::new();
        plan.push(Box::new(SampleAllelePass::new(
            Segment::V,
            Box::new(AllelePoolDist::restricted_uniform(
                &cfg.v_pool,
                vec![v_id],
            )),
        )));
        plan.push(Box::new(SampleAllelePass::new(
            Segment::D,
            Box::new(AllelePoolDist::restricted_uniform(
                &cfg.d_pool,
                vec![sampled_d],
            )),
        )));
        plan.push(Box::new(TrimPass::new(
            Segment::D,
            TrimEnd::Five,
            Box::new(EmpiricalLengthDist::from_pairs(vec![(d_trim_5, 1.0)])),
        )));
        plan.push(Box::new(TrimPass::new(
            Segment::D,
            TrimEnd::Three,
            Box::new(EmpiricalLengthDist::from_pairs(vec![(d_trim_3, 1.0)])),
        )));
        plan.push(Box::new(AssembleSegmentPass::new(Segment::V)));
        plan.push(Box::new(GenerateNPPass::new(
            Segment::Np1,
            Box::new(EmpiricalLengthDist::from_pairs(vec![(np1_len, 1.0)])),
            Box::new(ConstBaseDist(np1_base)),
        )));
        plan.push(Box::new(AssembleSegmentPass::new(Segment::D)));
        plan.push(Box::new(GenerateNPPass::new(
            Segment::Np2,
            Box::new(EmpiricalLengthDist::from_pairs(vec![(np2_len, 1.0)])),
            Box::new(ConstBaseDist(np2_base)),
        )));
        plan.push(Box::new(SampleAllelePass::new(
            Segment::J,
            Box::new(AllelePoolDist::restricted_uniform(
                &cfg.j_pool,
                vec![j_id],
            )),
        )));
        plan.push(Box::new(AssembleSegmentPass::new(Segment::J)));
        plan
    }

    #[test]
    fn d_call_shrinks_when_np1_recreates_trimmed_prefix() {
        // Sample D*01 (TTTCCCCCAAA), trim D_5 by 3, D_3 by 0.
        // Assembled D ref window starts at pos 3, covers 3..11 (CCCCCAAA).
        // BOTH D alleles match the assembled bases at those positions →
        // post-assemble d_call = {D*01, D*02}.
        // NP1 = TTT (length 3, all 'T'). The walker checks NP1's
        // RIGHTMOST base first against ref pos 2 (D*01[2]='T',
        // D*02[2]='G'). D*02 drops out. Two more 'T's match D*01[1]
        // and D*01[0] respectively → d_call shrinks to {D*01}.
        let (cfg, d01, _d02) = d_extension_refdata();

        let plan = d_extension_plan(&cfg, d01, 3, 0, 3, b'T', 0, b'A');
        let compiled =
            CompiledSimulator::compile(&plan, Some(&cfg), None, ExecutionPolicy::Permissive)
                .expect("plan should compile");
        let outcome = compiled.run_one(0).expect("plan should run");

        let final_call = outcome
            .final_simulation()
            .live_calls
            .as_ref()
            .expect("live calls populated")
            .get(Segment::D)
            .cloned()
            .expect("D live call exists");
        assert_eq!(final_call.allele_call.to_ids(), vec![d01]);

        let elastic = final_call
            .hypotheses
            .iter()
            .find(|h| {
                h.flags
                    .contains(crate::live_call::HypothesisFlags::BOUNDARY_ELASTIC)
            })
            .expect("D hypothesis should be flagged BOUNDARY_ELASTIC");
        // ref_start should be 0 — the extension recovered all 3
        // trimmed prefix bases.
        assert_eq!(elastic.ref_start, 0);
    }

    #[test]
    fn d_call_shrinks_when_np2_recreates_trimmed_suffix() {
        // Sample D*01 (TTTCCCCCAAA), trim D_5 by 0, D_3 by 3.
        // Assembled D covers 0..8 (TTTCCCCC). Both alleles match the
        // assembled bases at those positions only after first checking
        // pos 0 (D*01[0]='T', D*02[0]='G' → D*02 drops out at primary
        // walk!). Hmm — that means d_call is already {D*01} at the
        // primary walk. So this test really checks RIGHT extension on
        // an already-singleton call: NP2 should still narrow ref_end
        // and set BOUNDARY_ELASTIC.
        //
        // To get widening via D_3 trim alone, we'd need alleles that
        // share the prefix. Let's adjust:
        let mut cfg = RefDataConfig::empty(ChainType::Vdj);
        let _ = cfg.v_pool.push(Allele {
            name: "V*01".into(),
            gene: "V".into(),
            seq: b"AAA".to_vec(),
            segment: Segment::V,
            anchor: None,
        });
        let d01 = cfg.d_pool.push(Allele {
            name: "D*01".into(),
            gene: "D".into(),
            seq: b"AAACCCCCTTT".to_vec(),
            segment: Segment::D,
            anchor: None,
        });
        let _d02 = cfg.d_pool.push(Allele {
            name: "D*02".into(),
            gene: "D".into(),
            seq: b"AAACCCCCGGG".to_vec(),
            segment: Segment::D,
            anchor: None,
        });
        let _ = cfg.j_pool.push(Allele {
            name: "J*01".into(),
            gene: "J".into(),
            seq: b"AAA".to_vec(),
            segment: Segment::J,
            anchor: None,
        });
        // Sample D*01, trim D_3 by 3, D_5 by 0. Assembled D = AAACCCCC
        // (ref pos 0..8). Both alleles match → live call = {D*01, D*02}.
        // NP2 = TTT (length 3). The walker extends right: NP2[0]='T'
        // vs ref pos 8 (D*01[8]='T' ✓, D*02[8]='G' ✗) → D*02 out.
        // Continue: NP2[1] vs pos 9 (D*01[9]='T' ✓), NP2[2] vs pos 10
        // (D*01[10]='T' ✓). Final d_call = {D*01}.
        let plan = d_extension_plan(&cfg, d01, 0, 3, 0, b'A', 3, b'T');
        let compiled =
            CompiledSimulator::compile(&plan, Some(&cfg), None, ExecutionPolicy::Permissive)
                .expect("plan should compile");
        let outcome = compiled.run_one(0).expect("plan should run");

        let final_call = outcome
            .final_simulation()
            .live_calls
            .as_ref()
            .expect("live calls populated")
            .get(Segment::D)
            .cloned()
            .expect("D live call exists");
        assert_eq!(final_call.allele_call.to_ids(), vec![d01]);

        let elastic = final_call
            .hypotheses
            .iter()
            .find(|h| {
                h.flags
                    .contains(crate::live_call::HypothesisFlags::BOUNDARY_ELASTIC)
            })
            .expect("D hypothesis should be flagged BOUNDARY_ELASTIC");
        // ref_end should now be 11 — extension covered the full
        // post-trim allele length.
        assert_eq!(elastic.ref_end, 11);
    }

    #[test]
    fn d_call_shrinks_via_both_sides_simultaneously() {
        // D*01 = TTT CCCCC AAA, D*02 = GGG CCCCC TGT. Trim D_5 = 3,
        // D_3 = 3. Assembled D = CCCCC (ref pos 3..8). Both match
        // → d_call = {D*01, D*02}.
        //
        // NP1 = TTT (length 3) — narrows D left to {D*01} on
        // AssembleSegment(D).
        // NP2 = AAA (length 3) — D right-extension into NP2 is
        // already only checking against {D*01}; ref_end advances to
        // 11 (full allele length) after AppendRegion(Np2).
        //
        // Final state: d_call = {D*01}, BOUNDARY_ELASTIC flag set,
        // ref_start = 0, ref_end = 11.
        let (cfg, d01, _d02) = d_extension_refdata();
        let plan = d_extension_plan(&cfg, d01, 3, 3, 3, b'T', 3, b'A');

        let compiled =
            CompiledSimulator::compile(&plan, Some(&cfg), None, ExecutionPolicy::Permissive)
                .expect("plan should compile");
        let outcome = compiled.run_one(0).expect("plan should run");

        let final_call = outcome
            .final_simulation()
            .live_calls
            .as_ref()
            .expect("live calls populated")
            .get(Segment::D)
            .cloned()
            .expect("D live call exists");
        assert_eq!(final_call.allele_call.to_ids(), vec![d01]);

        let elastic = final_call
            .hypotheses
            .iter()
            .find(|h| {
                h.flags
                    .contains(crate::live_call::HypothesisFlags::BOUNDARY_ELASTIC)
            })
            .expect("D hypothesis should be flagged BOUNDARY_ELASTIC");
        assert_eq!(elastic.ref_start, 0);
        assert_eq!(elastic.ref_end, 11);
    }

    #[test]
    fn d_call_stays_widened_when_neither_np_matches() {
        // NP1 emits 'C' which is neither D*01[2]='T' nor D*02[2]='G'.
        // NP2 emits 'C' which is neither D*01[8]='A' nor D*02[8]='T'.
        // Both extension walks halt immediately; d_call stays widened.
        let (cfg, d01, d02) = d_extension_refdata();
        let plan = d_extension_plan(&cfg, d01, 3, 3, 3, b'C', 3, b'C');

        let compiled =
            CompiledSimulator::compile(&plan, Some(&cfg), None, ExecutionPolicy::Permissive)
                .expect("plan should compile");
        let outcome = compiled.run_one(0).expect("plan should run");

        let final_call = outcome
            .final_simulation()
            .live_calls
            .as_ref()
            .expect("live calls populated")
            .get(Segment::D)
            .cloned()
            .expect("D live call exists");
        let mut ids = final_call.allele_call.to_ids();
        ids.sort_by_key(|id| id.index());
        assert_eq!(ids, vec![d01, d02]);
        for h in &final_call.hypotheses {
            assert!(
                !h.flags
                    .contains(crate::live_call::HypothesisFlags::BOUNDARY_ELASTIC),
                "no extension occurred → BOUNDARY_ELASTIC must not be set"
            );
        }
    }

    #[test]
    fn append_region_np2_bumps_d_live_call_version() {
        // Plumbing check: AppendRegion(Np2) must trigger D refresh,
        // bumping the live-call evidence_version. Catches accidental
        // removal of the hook from `apply_live_call_updates`.
        let (cfg, d01, _d02) = d_extension_refdata();
        let plan = d_extension_plan(&cfg, d01, 3, 3, 3, b'T', 3, b'A');

        let compiled =
            CompiledSimulator::compile(&plan, Some(&cfg), None, ExecutionPolicy::Permissive)
                .expect("plan should compile");
        let outcome = compiled.run_one(0).expect("plan should run");

        // Find indices of the AssembleSegment(D) and AppendRegion(Np2)
        // pass commits by walking pass_names. revisions are
        // [initial, post-pass-0, post-pass-1, ...].
        let pass_names = &outcome.pass_names;
        let assemble_d_idx = pass_names
            .iter()
            .position(|n| n == "assemble.d")
            .expect("plan must include assemble.d");
        let np2_idx = pass_names
            .iter()
            .position(|n| n == "generate_np.np2")
            .expect("plan must include generate_np.np2");
        assert!(
            np2_idx > assemble_d_idx,
            "np2 generation should happen after D assembly"
        );

        let post_assemble_d_version = outcome.revisions[assemble_d_idx + 1]
            .live_calls
            .as_ref()
            .expect("post-assemble-D live calls present")
            .version;
        let post_np2_version = outcome.revisions[np2_idx + 1]
            .live_calls
            .as_ref()
            .expect("post-NP2 live calls present")
            .version;
        assert!(
            post_np2_version > post_assemble_d_version,
            "AppendRegion(Np2) must bump live-call version: \
             {post_np2_version} should be > {post_assemble_d_version}"
        );
    }

    // ──────────────────────────────────────────────────────────────
    // Structural indel local recomputation.
    //
    // Indels (insertions / deletions) shift the pool layout under
    // V/D/J regions. The walker now tolerates:
    //   - synthetic (indel-inserted) nucleotides inside V/D/J
    //     regions (skip without failing),
    //   - forward jumps in `germline_pos` caused by deletions
    //     between observed positions.
    // The `apply_live_call_updates` hook on `PassEffect::StructuralIndel`
    // refreshes V/D/J live calls so post-indel evidence drives the
    // call set.
    // ──────────────────────────────────────────────────────────────

    /// Test-only deterministic Pass that deletes one nucleotide at a
    /// fixed pool position and reports `PassEffect::StructuralIndel`.
    /// Mirrors the production `IndelPass`'s deletion path but without
    /// RNG involvement.
    #[derive(Clone, Debug)]
    struct DeleteAtPass {
        at: u32,
    }

    impl DeleteAtPass {
        fn new(at: u32) -> Self {
            Self { at }
        }
    }

    impl Pass for DeleteAtPass {
        fn name(&self) -> &str {
            "test.delete_at"
        }
        fn execute(
            &self,
            sim: &Simulation,
            _ctx: &mut crate::pass::PassContext,
        ) -> Simulation {
            sim.with_indel_deleted(self.at)
        }
        fn effects(&self) -> Vec<PassEffect> {
            vec![PassEffect::StructuralIndel]
        }
    }

    /// Test-only deterministic Pass that inserts one nucleotide at a
    /// fixed pool position with the chosen segment / base, and reports
    /// `PassEffect::StructuralIndel`. Mirrors `IndelPass`'s insertion
    /// path but without RNG involvement.
    #[derive(Clone, Debug)]
    struct InsertAtPass {
        at: u32,
        base: u8,
        segment: Segment,
    }

    impl InsertAtPass {
        fn new(at: u32, base: u8, segment: Segment) -> Self {
            Self { at, base, segment }
        }
    }

    impl Pass for InsertAtPass {
        fn name(&self) -> &str {
            "test.insert_at"
        }
        fn execute(
            &self,
            sim: &Simulation,
            _ctx: &mut crate::pass::PassContext,
        ) -> Simulation {
            let nuc = crate::ir::Nucleotide::synthetic(
                self.base,
                self.segment,
                crate::ir::flag::INDEL_INSERTED,
            );
            sim.with_indel_inserted(self.at, nuc)
        }
        fn effects(&self) -> Vec<PassEffect> {
            vec![PassEffect::StructuralIndel]
        }
    }

    /// V-only refdata holding alleles where the differing position is
    /// known and predictable.
    fn v_refdata(seqs: &[(&str, &[u8])]) -> RefDataConfig {
        let mut cfg = RefDataConfig::empty(ChainType::Vdj);
        for (name, seq) in seqs {
            let _ = cfg.v_pool.push(Allele {
                name: (*name).to_string(),
                gene: name.split('*').next().unwrap_or(name).to_string(),
                seq: seq.to_vec(),
                segment: Segment::V,
                anchor: None,
            });
        }
        cfg
    }

    /// Build a sample(V) → assemble(V) → indel plan and run it. Returns
    /// the final V live call.
    fn run_indel_plan(
        cfg: &RefDataConfig,
        sampled: AlleleId,
        indels: Vec<Box<dyn Pass>>,
    ) -> (
        crate::pass::Outcome,
        crate::live_call::SegmentLiveCall,
    ) {
        let mut plan = PassPlan::new();
        plan.push(Box::new(SampleAllelePass::new(
            Segment::V,
            Box::new(AllelePoolDist::restricted_uniform(
                &cfg.v_pool,
                vec![sampled],
            )),
        )));
        plan.push(Box::new(AssembleSegmentPass::new(Segment::V)));
        for indel in indels {
            plan.push(indel);
        }
        let compiled =
            CompiledSimulator::compile(&plan, Some(cfg), None, ExecutionPolicy::Permissive)
                .expect("fixture plan should compile");
        let outcome = compiled.run_one(0).expect("plan should run");
        let v_call = outcome
            .final_simulation()
            .live_calls
            .as_ref()
            .expect("live calls populated after assembly")
            .get(Segment::V)
            .cloned()
            .expect("V live call exists after assembly");
        (outcome, v_call)
    }

    #[test]
    fn deletion_inside_v_widens_live_call_when_distinguishing_base_removed() {
        // V*01 = AAAGAAAA (position 3 = G distinguishes it).
        // V*02 = AAACAAAA (position 3 = C distinguishes it).
        // Sample V*01 → assembled "AAAGAAAA" → live call {V*01} (the
        // G at pos 3 disambiguates).
        // Delete position 3 → assembled becomes "AAAAAAA" (7 bases),
        // covering germline positions 0,1,2,4,5,6,7. With pos 3 absent
        // the alleles are pairwise indistinguishable at the remaining
        // positions → live call widens to {V*01, V*02}.
        let cfg = v_refdata(&[
            ("V*01", b"AAAGAAAA"),
            ("V*02", b"AAACAAAA"),
        ]);
        let v01 = AlleleId::new(0);
        let v02 = AlleleId::new(1);

        let (_outcome, v_call) = run_indel_plan(
            &cfg,
            v01,
            vec![Box::new(DeleteAtPass::new(3))],
        );
        let mut ids = v_call.allele_call.to_ids();
        ids.sort_by_key(|id| id.index());
        assert_eq!(ids, vec![v01, v02]);
    }

    #[test]
    fn deletion_preserves_call_when_distinguishing_base_remains() {
        // Same fixture as above, but delete a SHARED base (position 0,
        // both alleles have 'A'). The distinguishing position 3 is
        // still there → live call stays singleton {V*01}.
        let cfg = v_refdata(&[
            ("V*01", b"AAAGAAAA"),
            ("V*02", b"AAACAAAA"),
        ]);
        let v01 = AlleleId::new(0);

        let (_outcome, v_call) = run_indel_plan(
            &cfg,
            v01,
            vec![Box::new(DeleteAtPass::new(0))],
        );
        assert_eq!(v_call.allele_call.to_ids(), vec![v01]);
    }

    #[test]
    fn insertion_inside_v_does_not_fail_the_call() {
        // Sample V*01 → assembled "AAAGAAAA" → live call {V*01}.
        // Insert a synthetic 'C' at pool position 3 (inside V's
        // region). The walker should SKIP the inserted nucleotide
        // (it has NO_GERMLINE_POS) without failing the call. The
        // remaining positions still uniquely identify V*01.
        let cfg = v_refdata(&[
            ("V*01", b"AAAGAAAA"),
            ("V*02", b"AAACAAAA"),
        ]);
        let v01 = AlleleId::new(0);

        let (_outcome, v_call) = run_indel_plan(
            &cfg,
            v01,
            vec![Box::new(InsertAtPass::new(3, b'C', Segment::V))],
        );
        assert_eq!(
            v_call.allele_call.to_ids(),
            vec![v01],
            "insertion should not collapse the live call to Unsupported"
        );
        assert!(
            !matches!(
                v_call.confidence,
                crate::live_call::LiveCallConfidence::Unsupported
            ),
            "live call must remain supported after a single in-region insertion"
        );
    }

    #[test]
    fn combined_insertion_and_deletion_recomputes_correctly() {
        // Sample V*01 → live call {V*01} via the distinguishing G at
        // position 3. Then:
        //   1. Insert a synthetic 'C' at pool position 5 (a benign
        //      mid-region insertion the walker skips).
        //   2. Delete the original distinguishing base at germline
        //      position 3 (which is now at pool position 3 still —
        //      the insertion at 5 didn't shift pos 0..4).
        // Expected post-state: live call widens to {V*01, V*02}
        // because the distinguishing germline position 3 is gone.
        let cfg = v_refdata(&[
            ("V*01", b"AAAGAAAA"),
            ("V*02", b"AAACAAAA"),
        ]);
        let v01 = AlleleId::new(0);
        let v02 = AlleleId::new(1);

        let (_outcome, v_call) = run_indel_plan(
            &cfg,
            v01,
            vec![
                Box::new(InsertAtPass::new(5, b'C', Segment::V)),
                Box::new(DeleteAtPass::new(3)),
            ],
        );
        let mut ids = v_call.allele_call.to_ids();
        ids.sort_by_key(|id| id.index());
        assert_eq!(ids, vec![v01, v02]);
    }

    #[test]
    fn structural_indel_bumps_live_call_version() {
        // Plumbing check: every refresh fired by `PassEffect::StructuralIndel`
        // bumps `evidence_version`. Catches accidental removal of the
        // hook from `apply_live_call_updates`.
        let cfg = v_refdata(&[
            ("V*01", b"AAAGAAAA"),
            ("V*02", b"AAACAAAA"),
        ]);
        let v01 = AlleleId::new(0);

        let (outcome, _final_call) = run_indel_plan(
            &cfg,
            v01,
            vec![Box::new(DeleteAtPass::new(3))],
        );
        let pass_names = &outcome.pass_names;
        let assemble_idx = pass_names
            .iter()
            .position(|n| n == "assemble.v")
            .expect("plan must include assemble.v");
        let indel_idx = pass_names
            .iter()
            .position(|n| n == "test.delete_at")
            .expect("plan must include test.delete_at");
        let post_assemble_version = outcome.revisions[assemble_idx + 1]
            .live_calls
            .as_ref()
            .expect("post-assemble live calls present")
            .version;
        let post_indel_version = outcome.revisions[indel_idx + 1]
            .live_calls
            .as_ref()
            .expect("post-indel live calls present")
            .version;
        assert!(
            post_indel_version > post_assemble_version,
            "StructuralIndel must bump live-call version: \
             {post_indel_version} should be > {post_assemble_version}"
        );
    }

    // ──────────────────────────────────────────────────────────────
    // Overlap hypotheses.
    //
    // When V's allele suffix happens to match D's leading bases (or
    // D's allele suffix matches J's leading bases, or symmetric
    // left-side cases), the live graph must retain BOTH placements
    // internally. We surface this via:
    // - the upstream segment's hypothesis growing seq_end past the
    //   downstream segment's seq_start,
    // - the upstream hypothesis carrying the OVERLAPS_OTHER_SEGMENT
    //   flag.
    // ──────────────────────────────────────────────────────────────

    /// Build a VDJ refdata where V*01's distinguishing 3' suffix
    /// "TTT" coincides with D*01's leading 3 bases. With V trimmed
    /// 3' by 3 and an empty NP1, V's right-extension can reach into
    /// D's region for an exact-equivalent overlap placement.
    fn v_d_overlap_refdata() -> (RefDataConfig, AlleleId, AlleleId, AlleleId) {
        let mut cfg = RefDataConfig::empty(ChainType::Vdj);
        // Two V alleles sharing AAACCCGGG and differing in their 3'
        // trinucleotide. After V_3 trim of 3, both match the
        // structural V region equally → primary live call is widened.
        let v01 = cfg.v_pool.push(Allele {
            name: "V*01".into(),
            gene: "V".into(),
            seq: b"AAACCCGGGTTT".to_vec(),
            segment: Segment::V,
            anchor: None,
        });
        let v02 = cfg.v_pool.push(Allele {
            name: "V*02".into(),
            gene: "V".into(),
            seq: b"AAACCCGGGAAA".to_vec(),
            segment: Segment::V,
            anchor: None,
        });
        // D*01 starts with "TTT" — that's exactly the trimmed-off
        // V*01 suffix, so V right-extension into D produces overlap.
        let d01 = cfg.d_pool.push(Allele {
            name: "D*01".into(),
            gene: "D".into(),
            seq: b"TTTACGTAC".to_vec(),
            segment: Segment::D,
            anchor: None,
        });
        // Stub J — never executed in fixtures that stop after
        // AssembleSegment(D).
        let _ = cfg.j_pool.push(Allele {
            name: "J*01".into(),
            gene: "J".into(),
            seq: b"AAA".to_vec(),
            segment: Segment::J,
            anchor: None,
        });
        (cfg, v01, v02, d01)
    }

    /// Plan that ends right after assembling D, so the V live call
    /// reflects the post-D state including any overlap.
    fn v_overlap_plan(
        cfg: &RefDataConfig,
        sampled_v: AlleleId,
        sampled_d: AlleleId,
        v_trim_3: i64,
        np1_len: i64,
    ) -> PassPlan {
        let mut plan = PassPlan::new();
        plan.push(Box::new(SampleAllelePass::new(
            Segment::V,
            Box::new(AllelePoolDist::restricted_uniform(
                &cfg.v_pool,
                vec![sampled_v],
            )),
        )));
        plan.push(Box::new(SampleAllelePass::new(
            Segment::D,
            Box::new(AllelePoolDist::restricted_uniform(
                &cfg.d_pool,
                vec![sampled_d],
            )),
        )));
        plan.push(Box::new(TrimPass::new(
            Segment::V,
            TrimEnd::Three,
            Box::new(EmpiricalLengthDist::from_pairs(vec![(v_trim_3, 1.0)])),
        )));
        plan.push(Box::new(AssembleSegmentPass::new(Segment::V)));
        plan.push(Box::new(GenerateNPPass::new(
            Segment::Np1,
            Box::new(EmpiricalLengthDist::from_pairs(vec![(np1_len, 1.0)])),
            Box::new(ConstBaseDist(b'C')),
        )));
        plan.push(Box::new(AssembleSegmentPass::new(Segment::D)));
        plan
    }

    #[test]
    fn v_right_overlaps_d_when_d_starts_with_v_suffix() {
        // V*01 trimmed 3' by 3 → V region = AAACCCGGG (ref 0..9).
        // V live call after AssembleSegment(V) = {V*01, V*02}.
        // After GenerateNP(NP1, length=0), nothing changes (NP1 empty).
        // After AssembleSegment(D): the cross-segment hook retriggers V.
        //   V right-extension walker enters with NP1 [9, 9), empty,
        //   loops `9..pool_len`. seq_pos = 9 is D's first base 'T'.
        //   compatible_alleles_at(ref_pos=9, base=T) → V*01[9]='T' ✓,
        //   V*02[9]='A' ✗ → V*02 drops, candidates = {V*01}.
        //   Continue 10, 11 (V*01 ref pos 10, 11 = T, T ✓).
        //   At seq_pos=12 (D's 4th base 'A'), V*01 ref pos 12 = out of
        //   range → halt.
        //
        // Final V hypothesis: candidates={V*01}, seq_end=12,
        // ref_end=12, BOUNDARY_ELASTIC + OVERLAPS_OTHER_SEGMENT set.
        // D's region.start = 9. V hypothesis seq_end (12) > D start
        // (9) → 3 bases of overlap.
        let (cfg, v01, _v02, d01) = v_d_overlap_refdata();
        let plan = v_overlap_plan(&cfg, v01, d01, 3, 0);

        let compiled =
            CompiledSimulator::compile(&plan, Some(&cfg), None, ExecutionPolicy::Permissive)
                .expect("plan should compile");
        let outcome = compiled.run_one(0).expect("plan should run");

        let final_sim = outcome.final_simulation();
        let v_call = final_sim
            .live_calls
            .as_ref()
            .expect("live calls populated")
            .get(Segment::V)
            .cloned()
            .expect("V live call exists");
        let d_call = final_sim
            .live_calls
            .as_ref()
            .expect("live calls populated")
            .get(Segment::D)
            .cloned()
            .expect("D live call exists");

        // V's call has shrunk back to {V*01} — the overlap walk
        // recovered the trimmed suffix from D's leading bases.
        assert_eq!(v_call.allele_call.to_ids(), vec![v01]);

        let v_hypothesis = v_call
            .hypotheses
            .iter()
            .find(|h| {
                h.flags
                    .contains(crate::live_call::HypothesisFlags::BOUNDARY_ELASTIC)
            })
            .expect("V hypothesis should be flagged BOUNDARY_ELASTIC");
        assert!(
            v_hypothesis
                .flags
                .contains(crate::live_call::HypothesisFlags::OVERLAPS_OTHER_SEGMENT),
            "V hypothesis must carry OVERLAPS_OTHER_SEGMENT after \
             extending past NP1 into D's region"
        );
        assert_eq!(v_hypothesis.seq_end, 12);
        assert_eq!(v_hypothesis.ref_end, 12);

        // D's hypothesis is computed independently from D's region;
        // its seq_start should still match D.region.start = 9.
        let d_hypothesis = &d_call.hypotheses[0];
        assert_eq!(
            d_hypothesis.seq_start, 9,
            "D hypothesis is derived from D's region; seq_start unchanged"
        );

        // Internal-state property the design doc calls out:
        // "V end may be greater than D start internally when evidence
        //  supports both."
        assert!(
            v_hypothesis.seq_end > d_hypothesis.seq_start,
            "internal overlap: V seq_end ({}) should exceed D seq_start ({})",
            v_hypothesis.seq_end,
            d_hypothesis.seq_start,
        );
    }

    #[test]
    fn v_right_does_not_overlap_when_d_does_not_match_v_suffix() {
        // Same V refdata, but D*01 starts with non-matching bases:
        let mut cfg = RefDataConfig::empty(ChainType::Vdj);
        let v01 = cfg.v_pool.push(Allele {
            name: "V*01".into(),
            gene: "V".into(),
            seq: b"AAACCCGGGTTT".to_vec(),
            segment: Segment::V,
            anchor: None,
        });
        let v02 = cfg.v_pool.push(Allele {
            name: "V*02".into(),
            gene: "V".into(),
            seq: b"AAACCCGGGAAA".to_vec(),
            segment: Segment::V,
            anchor: None,
        });
        // D starts with C — neither V allele has C at ref pos 9
        // (V*01[9]=T, V*02[9]=A) → walker halts on first attempt.
        let d01 = cfg.d_pool.push(Allele {
            name: "D*01".into(),
            gene: "D".into(),
            seq: b"CCCACGTAC".to_vec(),
            segment: Segment::D,
            anchor: None,
        });
        let _ = cfg.j_pool.push(Allele {
            name: "J*01".into(),
            gene: "J".into(),
            seq: b"AAA".to_vec(),
            segment: Segment::J,
            anchor: None,
        });
        let plan = v_overlap_plan(&cfg, v01, d01, 3, 0);

        let compiled =
            CompiledSimulator::compile(&plan, Some(&cfg), None, ExecutionPolicy::Permissive)
                .expect("plan should compile");
        let outcome = compiled.run_one(0).expect("plan should run");
        let v_call = outcome
            .final_simulation()
            .live_calls
            .as_ref()
            .expect("live calls populated")
            .get(Segment::V)
            .cloned()
            .expect("V live call exists");

        // V live call stays at the post-trim widened {V*01, V*02} —
        // no allele's ref pos 9 was 'C'.
        let mut ids = v_call.allele_call.to_ids();
        ids.sort_by_key(|id| id.index());
        assert_eq!(ids, vec![v01, v02]);
        for h in &v_call.hypotheses {
            assert!(
                !h.flags
                    .contains(crate::live_call::HypothesisFlags::OVERLAPS_OTHER_SEGMENT),
                "non-matching D bases must NOT produce overlap"
            );
        }
    }

    #[test]
    fn d_left_overlaps_v_when_v_ends_with_d_prefix() {
        // Mirror of the V→D overlap, but for D's left-extension
        // walker reaching backward into V's region.
        //
        // D*01 has prefix "GGG" — that's V*01's last 3 bases.
        // Trim D_5 by 3 → D's structural region drops the GGG prefix
        // and starts at ref pos 3. After AssembleSegment(D), D's
        // left-extension walks backward from D.start through the
        // empty NP1 into V's region. V's last 3 bases are GGG, which
        // happen to match D*01's trimmed prefix → overlap.
        let mut cfg = RefDataConfig::empty(ChainType::Vdj);
        let v01 = cfg.v_pool.push(Allele {
            name: "V*01".into(),
            gene: "V".into(),
            seq: b"AAACCCGGG".to_vec(),
            segment: Segment::V,
            anchor: None,
        });
        let _ = cfg.d_pool.push(Allele {
            name: "D*01".into(),
            gene: "D".into(),
            seq: b"GGGACGTAC".to_vec(),
            segment: Segment::D,
            anchor: None,
        });
        let _ = cfg.d_pool.push(Allele {
            name: "D*02".into(),
            gene: "D".into(),
            seq: b"TTTACGTAC".to_vec(),
            segment: Segment::D,
            anchor: None,
        });
        let _ = cfg.j_pool.push(Allele {
            name: "J*01".into(),
            gene: "J".into(),
            seq: b"AAA".to_vec(),
            segment: Segment::J,
            anchor: None,
        });
        let d01 = AlleleId::new(0);

        let mut plan = PassPlan::new();
        plan.push(Box::new(SampleAllelePass::new(
            Segment::V,
            Box::new(AllelePoolDist::restricted_uniform(
                &cfg.v_pool,
                vec![v01],
            )),
        )));
        plan.push(Box::new(SampleAllelePass::new(
            Segment::D,
            Box::new(AllelePoolDist::restricted_uniform(
                &cfg.d_pool,
                vec![d01],
            )),
        )));
        plan.push(Box::new(TrimPass::new(
            Segment::D,
            TrimEnd::Five,
            Box::new(EmpiricalLengthDist::from_pairs(vec![(3, 1.0)])),
        )));
        plan.push(Box::new(AssembleSegmentPass::new(Segment::V)));
        plan.push(Box::new(GenerateNPPass::new(
            Segment::Np1,
            Box::new(EmpiricalLengthDist::from_pairs(vec![(0, 1.0)])),
            Box::new(ConstBaseDist(b'C')),
        )));
        plan.push(Box::new(AssembleSegmentPass::new(Segment::D)));

        let compiled =
            CompiledSimulator::compile(&plan, Some(&cfg), None, ExecutionPolicy::Permissive)
                .expect("plan should compile");
        let outcome = compiled.run_one(0).expect("plan should run");

        let final_sim = outcome.final_simulation();
        let v_call = final_sim
            .live_calls
            .as_ref()
            .expect("live calls populated")
            .get(Segment::V)
            .cloned()
            .expect("V live call exists");
        let d_call = final_sim
            .live_calls
            .as_ref()
            .expect("live calls populated")
            .get(Segment::D)
            .cloned()
            .expect("D live call exists");

        // D's call should have BOUNDARY_ELASTIC + OVERLAPS_OTHER_SEGMENT
        // because the left-extension walker reached back through NP1
        // (empty) into V's last 3 bases (GGG) and they exactly
        // matched D*01's trimmed prefix.
        let d_hypothesis = d_call
            .hypotheses
            .iter()
            .find(|h| {
                h.flags
                    .contains(crate::live_call::HypothesisFlags::BOUNDARY_ELASTIC)
            })
            .expect("D hypothesis should be flagged BOUNDARY_ELASTIC");
        assert!(
            d_hypothesis
                .flags
                .contains(crate::live_call::HypothesisFlags::OVERLAPS_OTHER_SEGMENT),
            "D's left-extension overlap into V must set OVERLAPS_OTHER_SEGMENT"
        );

        // D hypothesis seq_start should now be 6 (V's last 3 bases
        // started at pool position 6); ref_start should be 0.
        let v_region_end = v_call.hypotheses[0].seq_end;
        assert!(
            d_hypothesis.seq_start < v_region_end,
            "D's left boundary should reach into V's region (D start {} < V end {})",
            d_hypothesis.seq_start,
            v_region_end,
        );
        assert_eq!(d_hypothesis.ref_start, 0);
    }

    #[test]
    fn overlap_walker_halts_at_pool_end() {
        // Sanity: even if D's bases continue to match V's allele's
        // continuation, the walker must halt at pool_len (no
        // out-of-bounds reads). Construct a fixture where V's allele
        // continuation matches arbitrarily many bases — by giving V
        // an allele whose tail matches D's bases — and verify the
        // walker stops cleanly.
        let mut cfg = RefDataConfig::empty(ChainType::Vdj);
        let v01 = cfg.v_pool.push(Allele {
            name: "V*01".into(),
            gene: "V".into(),
            // V's allele runs WAY past V's region: 100 bases long.
            seq: vec![b'A'; 100],
            segment: Segment::V,
            anchor: None,
        });
        let _ = cfg.d_pool.push(Allele {
            name: "D*01".into(),
            // D is also all 'A' bases, so V's continuation matches
            // arbitrarily into D — extension walker must halt at
            // pool_len, not run forever.
            gene: "D".into(),
            seq: vec![b'A'; 5],
            segment: Segment::D,
            anchor: None,
        });
        let _ = cfg.j_pool.push(Allele {
            name: "J*01".into(),
            gene: "J".into(),
            seq: b"AAA".to_vec(),
            segment: Segment::J,
            anchor: None,
        });
        let d01 = AlleleId::new(0);
        // Trim V_3 by 95 so V's structural region is just the first
        // 5 bases. V right-extension can extend up to 95 ref positions
        // — but D has only 5 bases of pool. Walker should halt at
        // pool_len after extending into all 5 D bases.
        let plan = v_overlap_plan(&cfg, v01, d01, 95, 0);

        let compiled =
            CompiledSimulator::compile(&plan, Some(&cfg), None, ExecutionPolicy::Permissive)
                .expect("plan should compile");
        let outcome = compiled.run_one(0).expect("plan should run");
        let v_call = outcome
            .final_simulation()
            .live_calls
            .as_ref()
            .expect("live calls populated")
            .get(Segment::V)
            .cloned()
            .expect("V live call exists");

        let v_hypothesis = v_call
            .hypotheses
            .iter()
            .find(|h| {
                h.flags
                    .contains(crate::live_call::HypothesisFlags::BOUNDARY_ELASTIC)
            })
            .expect("V hypothesis should be elastic");
        // Pool layout: V_region [0, 5) (after V_3 trim of 95), NP1
        // [5, 5), D_region [5, 10). pool_len = 10.
        assert_eq!(
            v_hypothesis.seq_end, 10,
            "walker should extend exactly to pool_len = 10, not beyond"
        );
        assert!(
            v_hypothesis
                .flags
                .contains(crate::live_call::HypothesisFlags::OVERLAPS_OTHER_SEGMENT),
            "extending past NP1 into D's region must set OVERLAPS_OTHER_SEGMENT"
        );
    }

    #[test]
    fn v_overlap_into_d_bumps_v_live_call_version() {
        // Plumbing check: AssembleSegment(D) must trigger a V
        // refresh under the cross-segment hook.
        let (cfg, v01, _v02, d01) = v_d_overlap_refdata();
        let plan = v_overlap_plan(&cfg, v01, d01, 3, 0);

        let compiled =
            CompiledSimulator::compile(&plan, Some(&cfg), None, ExecutionPolicy::Permissive)
                .expect("plan should compile");
        let outcome = compiled.run_one(0).expect("plan should run");

        let pass_names = &outcome.pass_names;
        let assemble_v_idx = pass_names
            .iter()
            .position(|n| n == "assemble.v")
            .expect("plan must include assemble.v");
        let assemble_d_idx = pass_names
            .iter()
            .position(|n| n == "assemble.d")
            .expect("plan must include assemble.d");
        let post_assemble_v_version = outcome.revisions[assemble_v_idx + 1]
            .live_calls
            .as_ref()
            .expect("post-assemble-V live calls present")
            .version;
        let post_assemble_d_version = outcome.revisions[assemble_d_idx + 1]
            .live_calls
            .as_ref()
            .expect("post-assemble-D live calls present")
            .version;
        assert!(
            post_assemble_d_version > post_assemble_v_version,
            "AssembleSegment(D) must bump the version (refreshing V) for \
             overlap detection: {post_assemble_d_version} should be > \
             {post_assemble_v_version}"
        );
    }

    // ──────────────────────────────────────────────────────────────
    // Curated end-to-end allele evaluation suite.
    //
    // A small, hand-crafted V/D/J refdata where every base is
    // predictable. Each allele has:
    //   - V (12bp, anchor at 6): shared 9bp prefix `AAACCCTGT` (with
    //     the conserved Cys "TGT" at 6-8) plus a 3bp distinguishing
    //     suffix (`AAA` / `CCC` / `GGG`),
    //   - D (12bp): 3bp distinguishing prefix + shared `GGGCCC` core
    //     + 3bp distinguishing suffix,
    //   - J (9bp, anchor at 3): 3bp distinguishing prefix + shared
    //     `TGGACG` (W codon TGG at 3-5).
    //
    // Sampling V1+D1+J1 with no trim and zero NP gives a 33bp
    // assembly whose junction (pool 6..30) is in-frame and codes
    // `CKKGPFKW` with no stop codons → productive.
    //
    // The shared cores let trim widen the live call to all three
    // alleles per segment; the distinguishing edges let NP bases
    // narrow the call back. Mutations / indels / corruption can
    // then be applied at known positions and the resulting AIRR
    // metadata is hand-checkable.
    // ──────────────────────────────────────────────────────────────

    /// Build the curated VDJ refdata. Allele ids are guaranteed in
    /// declaration order (V1=0, V2=1, V3=2; same for D / J pools).
    fn curated_v_d_j_refdata() -> RefDataConfig {
        let mut cfg = RefDataConfig::empty(ChainType::Vdj);
        // V pool — anchor=6 marks the conserved Cys (TGT) codon.
        let _ = cfg.v_pool.push(Allele {
            name: "V1*01".into(),
            gene: "V1".into(),
            seq: b"AAACCCTGTAAA".to_vec(),
            segment: Segment::V,
            anchor: Some(6),
        });
        let _ = cfg.v_pool.push(Allele {
            name: "V2*01".into(),
            gene: "V2".into(),
            seq: b"AAACCCTGTCCC".to_vec(),
            segment: Segment::V,
            anchor: Some(6),
        });
        let _ = cfg.v_pool.push(Allele {
            name: "V3*01".into(),
            gene: "V3".into(),
            seq: b"AAACCCTGTGGG".to_vec(),
            segment: Segment::V,
            anchor: Some(6),
        });
        // D pool — distinguishing 3bp prefix (always T-starting, so it
        // can never extend a V allele's distinguishing A/C/G suffix)
        // + 6bp shared core (`GGGCCC`) + distinguishing 3bp suffix
        // (always G-starting, so it can never extend backward into a
        // J allele's distinguishing A/C/T prefix). These guards keep
        // the boundary-overlap walker from accidentally narrowing
        // calls during simple trim-widening tests.
        let _ = cfg.d_pool.push(Allele {
            name: "D1*01".into(),
            gene: "D1".into(),
            seq: b"TTCGGGCCCGAG".to_vec(),
            segment: Segment::D,
            anchor: None,
        });
        let _ = cfg.d_pool.push(Allele {
            name: "D2*01".into(),
            gene: "D2".into(),
            seq: b"TATGGGCCCGCG".to_vec(),
            segment: Segment::D,
            anchor: None,
        });
        let _ = cfg.d_pool.push(Allele {
            name: "D3*01".into(),
            gene: "D3".into(),
            seq: b"TCGGGGCCCGTG".to_vec(),
            segment: Segment::D,
            anchor: None,
        });
        // J pool — anchor=3 marks the conserved W (TGG) codon.
        let _ = cfg.j_pool.push(Allele {
            name: "J1*01".into(),
            gene: "J1".into(),
            seq: b"AAATGGACG".to_vec(),
            segment: Segment::J,
            anchor: Some(3),
        });
        let _ = cfg.j_pool.push(Allele {
            name: "J2*01".into(),
            gene: "J2".into(),
            seq: b"CCCTGGACG".to_vec(),
            segment: Segment::J,
            anchor: Some(3),
        });
        let _ = cfg.j_pool.push(Allele {
            name: "J3*01".into(),
            gene: "J3".into(),
            seq: b"TTTTGGACG".to_vec(),
            segment: Segment::J,
            anchor: Some(3),
        });
        cfg
    }

    /// Build a curated VDJ plan parameterised by the events to apply.
    /// All defaults are no-op (no trim, zero-length NP), so a test
    /// only specifies what it needs to vary.
    #[allow(clippy::too_many_arguments)]
    fn curated_plan(
        cfg: &RefDataConfig,
        v_id: AlleleId,
        d_id: AlleleId,
        j_id: AlleleId,
        v_trim_3: i64,
        d_trim_5: i64,
        d_trim_3: i64,
        j_trim_5: i64,
        np1_len: i64,
        np1_base: u8,
        np2_len: i64,
        np2_base: u8,
    ) -> PassPlan {
        let mut plan = PassPlan::new();
        plan.push(Box::new(SampleAllelePass::new(
            Segment::V,
            Box::new(AllelePoolDist::restricted_uniform(&cfg.v_pool, vec![v_id])),
        )));
        plan.push(Box::new(SampleAllelePass::new(
            Segment::D,
            Box::new(AllelePoolDist::restricted_uniform(&cfg.d_pool, vec![d_id])),
        )));
        plan.push(Box::new(SampleAllelePass::new(
            Segment::J,
            Box::new(AllelePoolDist::restricted_uniform(&cfg.j_pool, vec![j_id])),
        )));
        if v_trim_3 > 0 {
            plan.push(Box::new(TrimPass::new(
                Segment::V,
                TrimEnd::Three,
                Box::new(EmpiricalLengthDist::from_pairs(vec![(v_trim_3, 1.0)])),
            )));
        }
        if d_trim_5 > 0 {
            plan.push(Box::new(TrimPass::new(
                Segment::D,
                TrimEnd::Five,
                Box::new(EmpiricalLengthDist::from_pairs(vec![(d_trim_5, 1.0)])),
            )));
        }
        if d_trim_3 > 0 {
            plan.push(Box::new(TrimPass::new(
                Segment::D,
                TrimEnd::Three,
                Box::new(EmpiricalLengthDist::from_pairs(vec![(d_trim_3, 1.0)])),
            )));
        }
        if j_trim_5 > 0 {
            plan.push(Box::new(TrimPass::new(
                Segment::J,
                TrimEnd::Five,
                Box::new(EmpiricalLengthDist::from_pairs(vec![(j_trim_5, 1.0)])),
            )));
        }
        plan.push(Box::new(AssembleSegmentPass::new(Segment::V)));
        plan.push(Box::new(GenerateNPPass::new(
            Segment::Np1,
            Box::new(EmpiricalLengthDist::from_pairs(vec![(np1_len, 1.0)])),
            Box::new(ConstBaseDist(np1_base)),
        )));
        plan.push(Box::new(AssembleSegmentPass::new(Segment::D)));
        plan.push(Box::new(GenerateNPPass::new(
            Segment::Np2,
            Box::new(EmpiricalLengthDist::from_pairs(vec![(np2_len, 1.0)])),
            Box::new(ConstBaseDist(np2_base)),
        )));
        plan.push(Box::new(AssembleSegmentPass::new(Segment::J)));
        plan
    }

    /// Run a curated plan, build the AIRR record, return both.
    fn curated_run(
        cfg: &RefDataConfig,
        plan: PassPlan,
    ) -> (crate::pass::Outcome, crate::airr_record::AirrRecord) {
        let compiled =
            CompiledSimulator::compile(&plan, Some(cfg), None, ExecutionPolicy::Permissive)
                .expect("curated plan should compile");
        let outcome = compiled.run_one(0).expect("plan should run");
        let rec = build_airr_record(&outcome, cfg, "curated");
        (outcome, rec)
    }

    // V1 / V2 / V3 / D1 / D2 / D3 / J1 / J2 / J3 ids for legibility.
    fn curated_ids() -> [AlleleId; 9] {
        [
            AlleleId::new(0), // V1
            AlleleId::new(1), // V2
            AlleleId::new(2), // V3
            AlleleId::new(0), // D1
            AlleleId::new(1), // D2
            AlleleId::new(2), // D3
            AlleleId::new(0), // J1
            AlleleId::new(1), // J2
            AlleleId::new(2), // J3
        ]
    }

    #[test]
    fn curated_baseline_productive_no_corruption() {
        // Sample V1+D1+J1, no trim, zero NP. The assembled sequence
        // is exactly V1 + D1 + J1 = 33bp with the conserved Cys at
        // pool 6 and W at pool 27. Junction = pool[6..30] (24bp)
        // codes `CKKGPFKW` with no stops → productive.
        let cfg = curated_v_d_j_refdata();
        let [v1, _, _, d1, _, _, j1, _, _] = curated_ids();
        let plan = curated_plan(&cfg, v1, d1, j1, 0, 0, 0, 0, 0, b'A', 0, b'A');
        let (_outcome, rec) = curated_run(&cfg, plan);

        // Sequence is exact concatenation of V1 + D1 + J1.
        assert_eq!(rec.sequence, "AAACCCTGTAAATTCGGGCCCGAGAAATGGACG");
        assert_eq!(rec.sequence_length, 33);

        // Live calls narrow to the sampled allele (distinguishing
        // bases at every segment edge → exactly one allele matches).
        assert_eq!(rec.v_call, "V1*01");
        assert_eq!(rec.d_call, "D1*01");
        assert_eq!(rec.j_call, "J1*01");

        // Junction = pool[6..30] = TGT + AAA + TTCGGGCCCGAG + AAA + TGG.
        // Codons: TGT AAA TTC GGG CCC GAG AAA TGG = C K F G P E K W.
        assert_eq!(rec.junction, "TGTAAATTCGGGCCCGAGAAATGG");
        assert_eq!(rec.junction_length, Some(24));
        assert_eq!(rec.junction_aa, "CKFGPEKW");
        assert_eq!(rec.productive, Some(true));
        assert_eq!(rec.vj_in_frame, Some(true));
        assert_eq!(rec.stop_codon, Some(false));

        // Pure recombination → no mutations, indels, or errors.
        assert_eq!(rec.n_mutations, 0);
        assert_eq!(rec.n_indels, 0);
        assert_eq!(rec.n_pcr_errors, 0);
        assert!(!rec.is_contaminant);

        // CIGARs are pure M with each segment's full length.
        assert_eq!(rec.v_cigar, "12M");
        assert_eq!(rec.d_cigar, "12M");
        assert_eq!(rec.j_cigar, "9M");

        // Identity is 1.0 (every base matches the source allele).
        assert_eq!(rec.v_identity, Some(1.0));
        assert_eq!(rec.d_identity, Some(1.0));
        assert_eq!(rec.j_identity, Some(1.0));

        // Coordinate self-consistency.
        assert_eq!(rec.v_sequence_start, Some(0));
        assert_eq!(rec.v_sequence_end, Some(12));
        assert_eq!(rec.d_sequence_start, Some(12));
        assert_eq!(rec.d_sequence_end, Some(24));
        assert_eq!(rec.j_sequence_start, Some(24));
        assert_eq!(rec.j_sequence_end, Some(33));
        assert_eq!(rec.v_germline_start, Some(0));
        assert_eq!(rec.v_germline_end, Some(12));
        assert_eq!(rec.j_germline_end, Some(9));

        // Locus is derived from V's "V1*01" prefix → empty (our
        // synthetic alleles don't start with IGH/IGK/etc).
        assert_eq!(rec.locus, "");
    }

    #[test]
    fn curated_v_trim_widens_v_call() {
        // Trim V_3 by 3 → V coding region drops the distinguishing
        // suffix → all three V alleles match the assembled V bases.
        // The live `v_call` must reflect that ambiguity.
        let cfg = curated_v_d_j_refdata();
        let [v1, _, _, d1, _, _, j1, _, _] = curated_ids();
        let plan = curated_plan(&cfg, v1, d1, j1, 3, 0, 0, 0, 0, b'A', 0, b'A');
        let (_outcome, rec) = curated_run(&cfg, plan);

        // V coding region is now 9 bases (12 - 3).
        assert_eq!(rec.v_sequence_end, Some(9));
        assert_eq!(rec.v_germline_end, Some(9));
        assert_eq!(rec.v_trim_3, 3);

        // Live v_call lists ALL THREE V alleles (declaration order).
        assert_eq!(rec.v_call, "V1*01,V2*01,V3*01");

        // D and J calls are unchanged singletons.
        assert_eq!(rec.d_call, "D1*01");
        assert_eq!(rec.j_call, "J1*01");

        // Sequence has dropped the trimmed-off V suffix; total
        // length is 33 - 3 = 30.
        assert_eq!(rec.sequence_length, 30);
        assert_eq!(rec.v_cigar, "9M");
    }

    #[test]
    fn curated_np1_recreates_v_suffix_narrows_v_call_back() {
        // Trim V_3 by 3 (live call widens to {V1,V2,V3}) AND
        // generate NP1 = "AAA" (V1's distinguishing suffix). The
        // V right-extension walker reaches into NP1 and narrows
        // the call back to {V1}.
        //
        // AIRR `v_germline_end` reads from the live-call hypothesis
        // bounds (12, the full V allele length after NP1 extension
        // claimed 3 bases), and `v_cigar` covers the extended span
        // (12M = 9 structural + 3 NP1 columns claimed by V's right
        // extension).
        let cfg = curated_v_d_j_refdata();
        let [v1, _, _, d1, _, _, j1, _, _] = curated_ids();
        let plan = curated_plan(&cfg, v1, d1, j1, 3, 0, 0, 0, 3, b'A', 0, b'A');
        let (_outcome, rec) = curated_run(&cfg, plan);

        assert_eq!(rec.v_call, "V1*01");
        assert_eq!(rec.v_germline_end, Some(12));
        // CIGAR runs 12M after column-relabelling claims the NP bases.
        assert_eq!(rec.v_cigar, "12M");
    }

    #[test]
    fn curated_d_trim_both_sides_widens_d_call() {
        // Trim D_5 by 3 AND D_3 by 3 → D coding region drops both
        // distinguishing edges, exposing only the 6bp shared core.
        // All three D alleles match.
        let cfg = curated_v_d_j_refdata();
        let [v1, _, _, d1, _, _, j1, _, _] = curated_ids();
        let plan = curated_plan(&cfg, v1, d1, j1, 0, 3, 3, 0, 0, b'A', 0, b'A');
        let (_outcome, rec) = curated_run(&cfg, plan);

        assert_eq!(rec.d_call, "D1*01,D2*01,D3*01");
        assert_eq!(rec.d_germline_start, Some(3));
        assert_eq!(rec.d_germline_end, Some(9));
        assert_eq!(rec.d_cigar, "6M");
    }

    #[test]
    fn curated_np_bases_narrow_d_call_from_both_sides() {
        // Trim D_5+D_3 widens d_call to all three. NP1 = D1's
        // distinguishing prefix bases narrow D's left back; NP2 = D1's
        // distinguishing suffix bases narrow D's right back.
        //
        // Use single-base ConstBaseDist for NPs, so we have to pick
        // NPs as homogeneous runs. D1 prefix = "TTC", suffix = "GAT" —
        // not single-base. The fixture's `curated_plan` builder uses
        // a single-base NP, so we drive this test by trimming only
        // the side we narrow per-NP. Here: D_3 trim only, NP2 = "G"
        // ×3 (matches the first base of D's suffix; D2/D3 also have
        // 'G' as suffix's first base by design, so this still narrows
        // to {D1} only when paired with a base hit at pos 10 / 11).
        //
        // For a fully narrow-from-both-sides exercise that matches
        // the homogeneous-NP constraint, we test D_5 trim with NP1
        // narrowing to a *subset* of D candidates — D1 prefix starts
        // 'T', as do D2 and D3. Pos 0 'T' alone keeps all three. So
        // single-base NP narrowing for D is best demonstrated via
        // *one side*; the multi-base mixed-content version requires
        // a per-position base distribution and is left to a richer
        // helper if needed.
        //
        // For simplicity, this test verifies that NP1 = "A" doesn't
        // narrow D (no D allele has 'A' at prefix pos 2 = the position
        // the left-extension walker checks first). The d_call stays
        // widened to all three after both trims expose only the
        // shared core.
        //
        // Pos-2 bases of the D prefixes: D1='C', D2='T', D3='G'.
        // 'A' matches none, so the walker halts immediately.
        let cfg = curated_v_d_j_refdata();
        let [v1, _, _, d1, _, _, j1, _, _] = curated_ids();
        let plan = curated_plan(&cfg, v1, d1, j1, 0, 3, 3, 0, 1, b'A', 0, b'A');
        let (_outcome, rec) = curated_run(&cfg, plan);

        let mut ids: Vec<&str> = rec.d_call.split(',').collect();
        ids.sort();
        assert_eq!(ids, vec!["D1*01", "D2*01", "D3*01"]);
    }

    #[test]
    fn curated_j_trim_widens_j_call() {
        // Trim J_5 by 3 → J coding region drops the distinguishing
        // 3bp prefix → all three J alleles match the remaining
        // shared `TGGACG`.
        let cfg = curated_v_d_j_refdata();
        let [v1, _, _, d1, _, _, j1, _, _] = curated_ids();
        let plan = curated_plan(&cfg, v1, d1, j1, 0, 0, 0, 3, 0, b'A', 0, b'A');
        let (_outcome, rec) = curated_run(&cfg, plan);

        assert_eq!(rec.j_call, "J1*01,J2*01,J3*01");
        assert_eq!(rec.j_germline_start, Some(3));
        assert_eq!(rec.j_cigar, "6M");
    }

    #[test]
    fn curated_np2_recreates_j_prefix_narrows_j_call_back() {
        // Trim J_5 by 3 (widens) AND NP2 = "AAA" (matches J1's
        // distinguishing prefix). J left-extension reaches backward
        // into NP2 and narrows j_call back to {J1}.
        //
        // `j_germline_start` reads from the live-call hypothesis (0 —
        // extension into NP2 dragged ref_start from 3 back to 0), and
        // `j_cigar` is 9M (3 claimed NP2 columns + 6 structural J
        // columns).
        let cfg = curated_v_d_j_refdata();
        let [v1, _, _, d1, _, _, j1, _, _] = curated_ids();
        let plan = curated_plan(&cfg, v1, d1, j1, 0, 0, 0, 3, 0, b'A', 3, b'A');
        let (_outcome, rec) = curated_run(&cfg, plan);

        assert_eq!(rec.j_call, "J1*01");
        assert_eq!(rec.j_germline_start, Some(0));
        // CIGAR runs 9M after column-relabelling claims the NP bases.
        assert_eq!(rec.j_cigar, "9M");
    }

    #[test]
    fn curated_mutation_switches_v_call_to_a_different_allele() {
        // Sample V1 (distinguishing suffix AAA at pool 9-11). Then
        // edit positions 9, 10, 11 to 'C' → assembled V matches V2's
        // suffix exactly. Live v_call should switch to {V2}.
        let cfg = curated_v_d_j_refdata();
        let [v1, _, _, d1, _, _, j1, _, _] = curated_ids();
        let mut plan = curated_plan(&cfg, v1, d1, j1, 0, 0, 0, 0, 0, b'A', 0, b'A');
        plan.push(Box::new(EditBaseAtPass::new(9, b'C')));
        plan.push(Box::new(EditBaseAtPass::new(10, b'C')));
        plan.push(Box::new(EditBaseAtPass::new(11, b'C')));
        let (_outcome, rec) = curated_run(&cfg, plan);

        // Provenance allele is still V1 (origin unchanged); live
        // call switched to V2.
        assert_eq!(rec.v_call, "V2*01");
    }

    #[test]
    fn curated_indel_deletion_widens_v_call() {
        // Sample V1, no trim. Delete the three distinguishing
        // V suffix bases (positions 9, 10, 11) one at a time.
        // Note: each deletion shrinks the pool, so subsequent
        // deletions need shifted positions. Easier: delete the
        // same position 9 three times (each time the new "9" is
        // the next surviving distinguishing base because the prior
        // deletion shifted things left by one).
        let cfg = curated_v_d_j_refdata();
        let [v1, _, _, d1, _, _, j1, _, _] = curated_ids();
        let mut plan = curated_plan(&cfg, v1, d1, j1, 0, 0, 0, 0, 0, b'A', 0, b'A');
        plan.push(Box::new(DeleteAtPass::new(9)));
        plan.push(Box::new(DeleteAtPass::new(9)));
        plan.push(Box::new(DeleteAtPass::new(9)));
        let (_outcome, rec) = curated_run(&cfg, plan);

        // V's distinguishing suffix is fully deleted → v_call
        // widens to all three V alleles.
        assert_eq!(rec.v_call, "V1*01,V2*01,V3*01");
        // Three deletion ops in V's CIGAR. (We don't assert
        // `rec.n_indels` because the test-only `DeleteAtPass` doesn't
        // write to `corrupt.indel.count` — that counter is owned by
        // the production indel pass.)
        let v_d_count: u32 = rec
            .v_cigar
            .split_terminator(|c: char| c.is_ascii_alphabetic())
            .filter_map(|s| s.parse::<u32>().ok())
            .zip(rec.v_cigar.matches(|c: char| c.is_ascii_alphabetic()))
            .filter(|(_, op)| *op == "D")
            .map(|(n, _)| n)
            .sum();
        assert_eq!(v_d_count, 3);
    }

    #[test]
    fn curated_indel_insertion_inside_v_does_not_break_call() {
        // Sample V1, insert a synthetic 'C' at pool position 3
        // (inside V's coding region). Walker tolerates the synthetic
        // base — call stays {V1}.
        let cfg = curated_v_d_j_refdata();
        let [v1, _, _, d1, _, _, j1, _, _] = curated_ids();
        let mut plan = curated_plan(&cfg, v1, d1, j1, 0, 0, 0, 0, 0, b'A', 0, b'A');
        plan.push(Box::new(InsertAtPass::new(3, b'C', Segment::V)));
        let (_outcome, rec) = curated_run(&cfg, plan);

        // Insertion is absorbed without collapsing the call.
        assert_eq!(rec.v_call, "V1*01");
        // Sequence grew by one base.
        assert_eq!(rec.sequence_length, 34);
        // V CIGAR contains exactly one I op (the insertion).
        assert!(
            rec.v_cigar.contains("1I"),
            "expected an I op in v_cigar, got {}",
            rec.v_cigar
        );
    }

    #[test]
    fn v_germline_end_reflects_np_extension() {
        // Invariant: `v_germline_end` reads from the live-call
        // hypothesis instead of the trim-derived structural range.
        // Trim V_3 by 3 (structural end = 9), then NP1 = "AAA"
        // recovers V1's allele bases at ref pos 9, 10, 11 → live
        // ref_end advances to 12. AIRR `v_germline_end` should be 12.
        let cfg = curated_v_d_j_refdata();
        let [v1, _, _, d1, _, _, j1, _, _] = curated_ids();
        let plan = curated_plan(&cfg, v1, d1, j1, 3, 0, 0, 0, 3, b'A', 0, b'A');
        let (_outcome, rec) = curated_run(&cfg, plan);

        // Trim says structural V end is 9; live extension says 12.
        // The live value wins.
        assert_eq!(rec.v_trim_3, 3);
        assert_eq!(rec.v_germline_end, Some(12));
    }

    #[test]
    fn d_germline_bounds_reflect_np_extension() {
        // Trim D_5 = 3, NP1 = "C" (single base) → D's left-extension
        // walker checks ref pos 2; D1[2]='C' matches → ref_start
        // moves from 3 → 2. AIRR `d_germline_start` should be 2.
        let cfg = curated_v_d_j_refdata();
        let [v1, _, _, d1, _, _, j1, _, _] = curated_ids();
        let plan = curated_plan(&cfg, v1, d1, j1, 0, 3, 0, 0, 1, b'C', 0, b'A');
        let (_outcome, rec) = curated_run(&cfg, plan);

        // Live d_call narrows to D1 (only D1 has C at pos 2).
        assert_eq!(rec.d_call, "D1*01");
        // Live ref_start = 2 (post-extension); structural would be 3.
        assert_eq!(rec.d_trim_5, 3);
        assert_eq!(rec.d_germline_start, Some(2));
    }

    #[test]
    fn v_sequence_end_reflects_np_extension() {
        // Invariant: `v_sequence_end` reads from the live-call
        // hypothesis seq_end, so when NP1 bases extend V's right
        // boundary the AIRR sequence coord grows past the structural
        // V-region end. Concrete: V_3 trim 3 + NP1="AAA" (V1 suffix)
        // → live seq_end = 12 (3 NP1 bases claimed).
        let cfg = curated_v_d_j_refdata();
        let [v1, _, _, d1, _, _, j1, _, _] = curated_ids();
        let plan = curated_plan(&cfg, v1, d1, j1, 3, 0, 0, 0, 3, b'A', 0, b'A');
        let (_outcome, rec) = curated_run(&cfg, plan);

        // Structural V region ended at pool position 9; live
        // hypothesis grew to 12 → AIRR record reports 12.
        assert_eq!(rec.v_sequence_end, Some(12));
    }

    #[test]
    fn j_sequence_start_reflects_np_extension() {
        // Symmetric: J_5 trim 3 + NP2="AAA" (J1 prefix) → J's
        // live seq_start moves left into NP2.
        //
        // Pool layout: V(12) + NP1(0) + D(12) + NP2(3) + J(6) = 33.
        //   J structural region at [27, 33).
        //   NP2 at [24, 27).
        // After J left-extension into NP2, live seq_start = 24.
        let cfg = curated_v_d_j_refdata();
        let [v1, _, _, d1, _, _, j1, _, _] = curated_ids();
        let plan = curated_plan(&cfg, v1, d1, j1, 0, 0, 0, 3, 0, b'A', 3, b'A');
        let (_outcome, rec) = curated_run(&cfg, plan);

        assert_eq!(rec.j_sequence_start, Some(24));
        // Structural J end stays where it was; only seq_start moves.
        assert_eq!(rec.j_sequence_end, Some(33));
    }

    #[test]
    fn no_extension_preserves_structural_germline_bounds() {
        // Sanity: when NO extension fires (no NP, or NP doesn't
        // match), AIRR germline bounds match the structural trim
        // range. (Live-call hypothesis ref_start/ref_end match the
        // structural values in this case.)
        let cfg = curated_v_d_j_refdata();
        let [v1, _, _, d1, _, _, j1, _, _] = curated_ids();
        // Trim V_3 = 3, no NP1 → V's right-extension halts
        // immediately at the V/D boundary (D1's first base 'T'
        // doesn't match any V allele's pos 9). v_germline_end stays
        // at 9 = structural post-trim end.
        let plan = curated_plan(&cfg, v1, d1, j1, 3, 0, 0, 0, 0, b'A', 0, b'A');
        let (_outcome, rec) = curated_run(&cfg, plan);

        assert_eq!(rec.v_trim_3, 3);
        assert_eq!(rec.v_germline_end, Some(9));
    }

    #[test]
    fn v_cigar_extends_into_claimed_np1_columns() {
        // V's CIGAR includes claimed NP1 columns. With
        // V_3 trim 3 and NP1="AAA" (V1's distinguishing suffix), V
        // right-extends 3 bases into NP1 → CIGAR runs 12 M ops total.
        let cfg = curated_v_d_j_refdata();
        let [v1, _, _, d1, _, _, j1, _, _] = curated_ids();
        let plan = curated_plan(&cfg, v1, d1, j1, 3, 0, 0, 0, 3, b'A', 0, b'A');
        let (_outcome, rec) = curated_run(&cfg, plan);

        assert_eq!(rec.v_cigar, "12M");
    }

    #[test]
    fn j_cigar_extends_into_claimed_np2_columns() {
        // Mirror case: J left-extends 3 bases into NP2.
        // Structural J = 6M; with extension J's CIGAR = 9M.
        let cfg = curated_v_d_j_refdata();
        let [v1, _, _, d1, _, _, j1, _, _] = curated_ids();
        let plan = curated_plan(&cfg, v1, d1, j1, 0, 0, 0, 3, 0, b'A', 3, b'A');
        let (_outcome, rec) = curated_run(&cfg, plan);

        assert_eq!(rec.j_cigar, "9M");
    }

    #[test]
    fn np1_string_drops_columns_claimed_by_v() {
        // when V's right extension reabsorbs NP1 bases,
        // those bases are no longer "non-templated" — np1 / np1_length
        // must drop them. With V_3 trim 3 and NP1="AAA" all three NP1
        // bases are claimed by V → np1 is empty.
        let cfg = curated_v_d_j_refdata();
        let [v1, _, _, d1, _, _, j1, _, _] = curated_ids();
        let plan = curated_plan(&cfg, v1, d1, j1, 3, 0, 0, 0, 3, b'A', 0, b'A');
        let (_outcome, rec) = curated_run(&cfg, plan);

        assert_eq!(rec.np1, "");
        assert_eq!(rec.np1_length, 0);
    }

    #[test]
    fn np2_string_drops_columns_claimed_by_j() {
        // Mirror of 11.4 for NP2. With J_5 trim 3 and NP2="AAA"
        // all three NP2 bases get claimed by J's left extension.
        let cfg = curated_v_d_j_refdata();
        let [v1, _, _, d1, _, _, j1, _, _] = curated_ids();
        let plan = curated_plan(&cfg, v1, d1, j1, 0, 0, 0, 3, 0, b'A', 3, b'A');
        let (_outcome, rec) = curated_run(&cfg, plan);

        assert_eq!(rec.np2, "");
        assert_eq!(rec.np2_length, 0);
    }

    #[test]
    fn v_alignment_end_extends_with_np_claim() {
        // `v_alignment_end` covers the claimed NP columns
        // too. Without extension v_alignment_end = 12 (full V region);
        // with V_3 trim 3 + NP1="AAA" the structural V region is 9
        // columns long but V claims 3 NP1 columns → v_alignment_end
        // moves out to 12 (still 12 columns but now 9 structural +
        // 3 NP).
        let cfg = curated_v_d_j_refdata();
        let [v1, _, _, d1, _, _, j1, _, _] = curated_ids();
        let plan = curated_plan(&cfg, v1, d1, j1, 3, 0, 0, 0, 3, b'A', 0, b'A');
        let (_outcome, rec) = curated_run(&cfg, plan);

        assert_eq!(rec.v_alignment_start, Some(0));
        assert_eq!(rec.v_alignment_end, Some(12));
    }

    #[test]
    fn j_alignment_start_extends_with_np_claim() {
        // Mirror: when J left-extends into NP2, j_alignment_start
        // moves leftward. Pool layout is V(12)+NP1(0)+D(12)+NP2(3)+J(6).
        // Structural J spans columns [24, 30); with NP2 claimed,
        // j_alignment_start = 24 - 3 = 21 (note: column = pool index
        // here because there are no synthetic insertions / deletions).
        let cfg = curated_v_d_j_refdata();
        let [v1, _, _, d1, _, _, j1, _, _] = curated_ids();
        let plan = curated_plan(&cfg, v1, d1, j1, 0, 0, 0, 3, 0, b'A', 3, b'A');
        let (_outcome, rec) = curated_run(&cfg, plan);

        assert_eq!(rec.j_alignment_start, Some(24));
        assert_eq!(rec.j_alignment_end, Some(33));
    }

    #[test]
    fn v_identity_counts_extended_columns() {
        // V's identity reflects matches over the full
        // claimed span. With V_3 trim 3 + NP1="AAA" (which match V1
        // exactly at pos 9-11), every column in V's extended span is
        // a match → identity = 1.0.
        let cfg = curated_v_d_j_refdata();
        let [v1, _, _, d1, _, _, j1, _, _] = curated_ids();
        let plan = curated_plan(&cfg, v1, d1, j1, 3, 0, 0, 0, 3, b'A', 0, b'A');
        let (_outcome, rec) = curated_run(&cfg, plan);

        assert_eq!(rec.v_identity, Some(1.0));
    }

    #[test]
    fn junction_locates_anchors_via_germline_pos() {
        // Baseline: the curated default plan places V's Cys anchor
        // at pool position 6 (V allele anchor=6, V at pool[0..12])
        // and J's W anchor at pool position 27
        // (J allele anchor=3, J at pool[24..33]). Junction =
        // pool[6..30], length 24 = the 3bp Cys + (V tail 3) +
        // (D 12) + (J head 3) + 3bp W.
        let cfg = curated_v_d_j_refdata();
        let [v1, _, _, d1, _, _, j1, _, _] = curated_ids();
        let plan = curated_plan(&cfg, v1, d1, j1, 0, 0, 0, 0, 0, b'A', 0, b'A');
        let (_outcome, rec) = curated_run(&cfg, plan);

        assert_eq!(rec.junction_start, Some(6));
        assert_eq!(rec.junction_end, Some(30));
        assert_eq!(rec.junction_length, Some(24));
        // The anchor codons frame the junction.
        assert!(
            rec.junction.starts_with("TGT"),
            "junction should start with V Cys codon, got {:?}",
            rec.junction,
        );
        assert!(
            rec.junction.ends_with("TGG"),
            "junction should end with J W codon, got {:?}",
            rec.junction,
        );
    }

    #[test]
    fn junction_shifts_with_v_insertion_before_anchor() {
        // Insert a synthetic 'C' at pool position 3 (inside V, BEFORE
        // the anchor codon at pool 6). The junction must follow the
        // anchor's actual pool position — anchor germline position 6
        // now resides at pool index 7, so junction_start shifts from
        // 6 → 7.
        let cfg = curated_v_d_j_refdata();
        let [v1, _, _, d1, _, _, j1, _, _] = curated_ids();
        let mut plan = curated_plan(&cfg, v1, d1, j1, 0, 0, 0, 0, 0, b'A', 0, b'A');
        plan.push(Box::new(InsertAtPass::new(3, b'C', Segment::V)));
        let (_outcome, rec) = curated_run(&cfg, plan);

        assert_eq!(rec.junction_start, Some(7));
        assert!(
            rec.junction.starts_with("TGT"),
            "junction should still frame the V Cys codon after V insertion, got {:?}",
            rec.junction,
        );
    }

    #[test]
    fn junction_shifts_with_v_deletion_before_anchor() {
        // Symmetric: delete at pool position 3 (inside V, before the
        // anchor codon). Anchor germline_pos=6 now sits at pool 5
        // (one earlier). junction_start should follow.
        let cfg = curated_v_d_j_refdata();
        let [v1, _, _, d1, _, _, j1, _, _] = curated_ids();
        let mut plan = curated_plan(&cfg, v1, d1, j1, 0, 0, 0, 0, 0, b'A', 0, b'A');
        plan.push(Box::new(DeleteAtPass::new(3)));
        let (_outcome, rec) = curated_run(&cfg, plan);

        assert_eq!(rec.junction_start, Some(5));
        assert!(
            rec.junction.starts_with("TGT"),
            "junction should still frame the V Cys codon after V deletion, got {:?}",
            rec.junction,
        );
    }

    #[test]
    fn junction_shifts_with_j_insertion_before_anchor() {
        // J side mirror: inserting a base inside J (before its W
        // anchor at allele pos 3) pushes J's anchor pool position
        // rightward, so junction_end grows accordingly.
        //
        // Pool layout: V(12) + NP1(0) + D(12) + NP2(0) + J(9) = 33.
        //   J at pool[24..33]; anchor germline_pos=3 sits at pool 27.
        //   junction_end = 27 + 3 = 30 baseline.
        // After inserting at pool 25 (within J, before anchor), the
        // anchor germline_pos=3 node shifts to pool 28, so
        // junction_end = 28 + 3 = 31.
        let cfg = curated_v_d_j_refdata();
        let [v1, _, _, d1, _, _, j1, _, _] = curated_ids();
        let mut plan = curated_plan(&cfg, v1, d1, j1, 0, 0, 0, 0, 0, b'A', 0, b'A');
        plan.push(Box::new(InsertAtPass::new(25, b'C', Segment::J)));
        let (_outcome, rec) = curated_run(&cfg, plan);

        assert_eq!(rec.junction_end, Some(31));
        assert!(
            rec.junction.ends_with("TGG"),
            "junction should still frame the J W codon after J insertion, got {:?}",
            rec.junction,
        );
    }

    #[test]
    fn productive_uses_germline_pos_anchor_after_indels() {
        // Insert a synthetic 'C' at pool position 3 (inside V, before
        // the anchor). Both V and J anchors shift right by 1, so
        // junction_start = 7, junction_end = 31, junction_length = 24
        // — same length as baseline, still in-frame.
        //
        // The biological invariant: the anchor codon is still
        // recognised as Cys (germline pos 6,7,8 still map to T,G,T
        // in the IR). A structural-offset reader would have read
        // pool[6..9] (now AAC after shift) and rejected the codon,
        // marking productive=false even though the underlying allele
        // is intact.
        let cfg = curated_v_d_j_refdata();
        let [v1, _, _, d1, _, _, j1, _, _] = curated_ids();
        let mut plan = curated_plan(&cfg, v1, d1, j1, 0, 0, 0, 0, 0, b'A', 0, b'A');
        plan.push(Box::new(InsertAtPass::new(3, b'C', Segment::V)));
        let (_outcome, rec) = curated_run(&cfg, plan);

        assert_eq!(rec.junction_start, Some(7));
        assert_eq!(rec.junction_end, Some(31));
        assert_eq!(rec.junction_length, Some(24));
        assert_eq!(rec.vj_in_frame, Some(true));
        assert_eq!(rec.productive, Some(true));
    }

    #[test]
    fn germline_span_equals_m_plus_d() {
        // Invariant: AIRR `*_germline_start/_end` come from the
        // column walker's `ref_ranges` (the union of ref positions
        // consumed by `M` and `D` ops), so the strict identity
        // `germline_span == M + D` holds by construction.
        //
        // The curated default plan is no-op trims, no NP, no
        // corruption — V/D/J each contribute a clean structural
        // CIGAR (12M / 12M / 9M) over their full allele span, and
        // germline coords reflect that 1:1.
        let cfg = curated_v_d_j_refdata();
        let [v1, _, _, d1, _, _, j1, _, _] = curated_ids();
        let plan = curated_plan(&cfg, v1, d1, j1, 0, 0, 0, 0, 0, b'A', 0, b'A');
        let (_outcome, rec) = curated_run(&cfg, plan);

        for (cig, g_s, g_e) in [
            (rec.v_cigar.as_str(), rec.v_germline_start, rec.v_germline_end),
            (rec.d_cigar.as_str(), rec.d_germline_start, rec.d_germline_end),
            (rec.j_cigar.as_str(), rec.j_germline_start, rec.j_germline_end),
        ] {
            let m_count: u32 = cig
                .split_terminator(|c: char| c.is_ascii_alphabetic())
                .filter_map(|s| s.parse::<u32>().ok())
                .zip(cig.matches(|c: char| c.is_ascii_alphabetic()))
                .filter(|(_, op)| *op == "M" || *op == "D")
                .map(|(n, _)| n)
                .sum();
            let span = g_e.unwrap() - g_s.unwrap();
            assert_eq!(
                span, m_count as i64,
                "germline_span ({span}) != M+D ({m_count}) for cigar {cig:?}",
            );
        }
    }

    #[test]
    fn germline_span_under_v_indel_deletion() {
        // V_3 trim 3 + NP1="AAA" extends V to ref 12 via NP1 claim,
        // then a structural deletion at pool 1 removes one V base
        // (germline_pos=1). CIGAR walks the structural region and
        // emits a D op for the missing ref=1 — `ref_ranges` covers
        // [0, 12) (NP1 claim at ref=9..12, structural at ref=0..9
        // including the D-fill at ref=1). M+D == 12 == germline_span.
        let cfg = curated_v_d_j_refdata();
        let [v1, _, _, d1, _, _, j1, _, _] = curated_ids();
        let mut plan = curated_plan(&cfg, v1, d1, j1, 3, 0, 0, 0, 3, b'A', 0, b'A');
        plan.push(Box::new(DeleteAtPass::new(1)));
        let (_outcome, rec) = curated_run(&cfg, plan);

        let v_span = rec.v_germline_end.unwrap() - rec.v_germline_start.unwrap();
        let cig = rec.v_cigar.as_str();
        let m_count: u32 = cig
            .split_terminator(|c: char| c.is_ascii_alphabetic())
            .filter_map(|s| s.parse::<u32>().ok())
            .zip(cig.matches(|c: char| c.is_ascii_alphabetic()))
            .filter(|(_, op)| *op == "M" || *op == "D")
            .map(|(n, _)| n)
            .sum();
        assert_eq!(
            v_span, m_count as i64,
            "germline_span {v_span} != M+D {m_count} for v_cigar {cig:?}",
        );
    }

    #[test]
    fn curated_full_corruption_stack_keeps_metadata_self_consistent() {
        // Stack a mutation, an insertion, and a deletion across V's
        // region. All metadata (sequence, alignments, CIGAR, identity)
        // must remain mutually consistent.
        let cfg = curated_v_d_j_refdata();
        let [v1, _, _, d1, _, _, j1, _, _] = curated_ids();
        let mut plan = curated_plan(&cfg, v1, d1, j1, 0, 0, 0, 0, 0, b'A', 0, b'A');
        // Mutate position 3 (inside V's shared prefix).
        plan.push(Box::new(EditBaseAtPass::new(3, b'T')));
        // Insert a synthetic 'G' at position 6 (inside V).
        plan.push(Box::new(InsertAtPass::new(6, b'G', Segment::V)));
        // Delete position 1 (also inside V).
        plan.push(Box::new(DeleteAtPass::new(1)));
        let (_outcome, rec) = curated_run(&cfg, plan);

        // Mutual consistency: alignment columns = (sequence chars) +
        // (gap columns where sequence is "-").
        let n_seq_gaps = rec.sequence_alignment.matches('-').count();
        assert_eq!(
            rec.sequence_alignment.len(),
            rec.sequence.len() + n_seq_gaps,
        );

        // Removing gaps from sequence_alignment yields the
        // upper-cased sequence.
        let sa_no_gaps: String = rec
            .sequence_alignment
            .chars()
            .filter(|c| *c != '-')
            .collect();
        assert_eq!(sa_no_gaps, rec.sequence.to_uppercase());

        // sequence_alignment / germline_alignment / d_mask all share
        // the same length.
        assert_eq!(rec.sequence_alignment.len(), rec.germline_alignment.len());
        assert_eq!(
            rec.sequence_alignment.len(),
            rec.germline_alignment_d_mask.len()
        );

        // (We don't assert `rec.n_indels` because the test-only
        // `InsertAtPass` / `DeleteAtPass` don't write to
        // `corrupt.indel.count`. The structural-consistency checks
        // above are what this fixture is validating.)
    }
}
