//! Compiled simulator boundary.
//!
//! `PassPlan` is the engine's intermediate representation: an ordered
//! sequence of concrete simulation steps. It should not be the final
//! execution artifact. A `CompiledSimulator` binds that IR to the
//! reference data, active contracts, execution policy, and compile-time
//! analysis report needed to run safely and introspectably.

use crate::assignment::TrimEnd;
use crate::contract::ContractSet;
use crate::feasibility::FeasibilityContext;
use crate::ir::{Segment, Simulation};
use crate::live_call::ReferenceMatchIndex;
use crate::pass::{Outcome, PassError, PassPlan};
use crate::refdata::RefDataConfig;

// Test-only re-imports — the `compiled_tests.rs` submodule pulls these
// names through `use super::*`. Gated under `cfg(test)` so non-test
// builds don't carry unused imports.
#[cfg(test)]
use crate::event::TraceSpan;
#[cfg(test)]
use crate::pass::{PassContext, PassEffect};
#[cfg(test)]
use crate::refdata::AlleleId;

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

mod report;
pub use report::{CompileReport, CompileWarning, DeclaredChoice, PassSummary};

mod error;
pub use error::{CompileError, CompileErrorKind, CompileErrors};

mod execute;
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

mod analyze;
use analyze::analyze_plan;

mod feasibility_builder;

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

#[cfg(test)]
#[path = "../compiled_tests.rs"]
mod tests;
