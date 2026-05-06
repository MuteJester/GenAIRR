//! Pass manager — the orchestrator that walks a fixed sequence of
//! passes, applying each one to the simulation IR and recording any
//! random choices to the trace.
//!
//! ## Architectural commitments reflected here
//!
//! - **Persistent IR (D1):** `Pass::execute` takes `&Simulation` and
//!   returns a new `Simulation`. The receiver is never mutated.
//! - **Side-channel trace (B.1):** the `Trace` lives on the runtime,
//!   not on `Simulation`. Passes record choices into `ctx.trace`
//!   while constructing their next IR revision.
//! - **Fixed user-specified order (§4 of design doc):** the runtime
//!   walks the plan exactly as given. No reordering, no automatic
//!   parallelization within a plan. Biology has a natural temporal
//!   order and we follow it.
//! - **Per-action passes (D2):** every IR revision in the history
//!   corresponds to exactly one `Pass::execute` call. The history
//!   length equals the plan length plus one (for the initial IR).
//!
//! ## Phase B.2 scope
//!
//! Just the orchestration plumbing. No real biology yet — the only
//! pass implementations come in B.4 (`EchoPass`) and B.5
//! (`SampleBasePass`). This step proves the runtime correctly
//! threads `&Simulation → Simulation` through arbitrarily many
//! passes while collecting trace + history.

use crate::contract::{ContractSet, ContractViolation};
use crate::dist::FilteredSampleError;
use crate::event::EventRecord;
use crate::ir::{Segment, Simulation};
use crate::refdata::RefDataConfig;
use crate::rng::Rng;
use crate::trace::Trace;
use std::fmt;

// ──────────────────────────────────────────────────────────────────
// Pass metadata — compile-time requirements/effects
// ──────────────────────────────────────────────────────────────────

/// Static requirement a pass has before it can execute safely.
///
/// The compiled-simulator boundary consumes this metadata to reject
/// invalid plans before runtime. It is intentionally small and typed:
/// build-time analysis should not infer semantics from pass names or
/// trace-address strings.
#[derive(Copy, Clone, Debug, Eq, PartialEq)]
pub enum PassRequirement {
    /// Pass needs reference data in `PassContext::refdata`.
    RefData,
    /// Pass needs an allele assignment for the given segment to have
    /// been produced by an earlier pass.
    AlleleAssignment(Segment),
}

/// Static effect a pass has on the simulation state.
///
/// This is not a full dataflow language. It is the durable surface the
/// compiler layer needs for basic ordering validation and introspection.
/// More precise event/effect metadata can be added without changing the
/// existing pass execution contract.
#[derive(Copy, Clone, Debug, Eq, PartialEq)]
pub enum PassEffect {
    /// Assigns an allele id to a recombinable segment.
    AssignAllele(Segment),
    /// Updates trim metadata on an assigned allele instance.
    TrimAllele(Segment),
    /// Assembles a germline segment into the nucleotide pool/sequence.
    AssembleSegment(Segment),
    /// Appends a region to the sequence.
    AppendRegion(Segment),
    /// Appends one or more nucleotides without creating a region.
    AppendNucleotides,
    /// Edits existing nucleotide bases without changing pool length.
    EditBases,
    /// Changes pool length and shifts downstream handles/ranges.
    StructuralIndel,
}

// ──────────────────────────────────────────────────────────────────
// PassContext — what every pass sees in addition to the IR
// ──────────────────────────────────────────────────────────────────

/// Context handed to every `Pass::execute` invocation.
///
/// Carries the trace (where choices get recorded), the RNG (where
/// random draws come from), an optional `RefDataConfig` reference
/// (for passes that need allele sequences — C.8 `AssembleSegmentPass`
/// and beyond), and an optional `ContractSet` (D.6 — sampling
/// passes consult this for filtering candidates).
///
/// **`refdata: Option<&RefDataConfig>`** is opt-in:
/// - Passes that don't need reference data (sampling, trim, NP
///   generation, the Echo / SampleBase reference passes) ignore
///   it. They never panic on `None`.
/// - Passes that need it (assembly, future SHM mutation against
///   germline) call `ctx.refdata.expect("...")` with a clear
///   message describing the requirement.
/// - Tests that don't exercise reference-data-consuming passes
///   continue to work with `PassRuntime::execute`, which threads
///   `None`. Production paths use `execute_with_refdata` or
///   `execute_with_context`.
///
/// **`contracts: Option<&ContractSet>`** is opt-in:
/// - Constraint-aware sampling passes consult the contract set
///   via `dist::sample_filtered` to draw only admissible
///   candidates. `None` means no contracts active (MIXED mode).
///
/// **Forward-compatibility (`#[non_exhaustive]`):** Phase E will
/// almost certainly add fields (S5F kernel state, SHM mutation
/// rate, indel position buffers, …). The non-exhaustive marker
/// means *external* construction of `PassContext` is impossible,
/// so future field additions are non-breaking for any consumer
/// outside this crate. Internal construction (in `PassRuntime`)
/// still uses struct-literal syntax and will need to acknowledge
/// each new field — which has audit value.
#[non_exhaustive]
pub struct PassContext<'a> {
    pub trace: &'a mut Trace,
    pub rng: &'a mut Rng,
    pub refdata: Option<&'a RefDataConfig>,
    /// Active contract set (D.6). Sampling passes consult this
    /// via `dist::sample_filtered` to draw only admissible
    /// candidates. `None` means no contracts are active —
    /// equivalent to MIXED mode.
    pub contracts: Option<&'a ContractSet>,
}

// ──────────────────────────────────────────────────────────────────
// PassError — structured fallible execution errors
// ──────────────────────────────────────────────────────────────────

/// Structured failure from a fallible pass execution.
///
/// Existing `execute*` runtime entry points remain permissive and
/// return `Outcome` directly. Strict entry points return this error
/// instead of silently falling back when a declared contract cannot be
/// satisfied at sampling time.
#[derive(Clone, Debug, Eq, PartialEq)]
pub enum PassError {
    /// A contract-aware sampling pass could not produce an admissible
    /// candidate at `address`.
    ConstraintSampling {
        pass_name: String,
        address: String,
        reason: FilteredSampleError,
    },
    /// A pass that requires reference data was run without it.
    MissingRefData { pass_name: String },
    /// A pass requires an allele assignment that is absent in the
    /// current IR revision. Usually indicates invalid pass ordering.
    MissingAssignment { pass_name: String, segment: Segment },
    /// The IR references an allele id that is absent from the supplied
    /// reference data.
    MissingAllele {
        pass_name: String,
        segment: Segment,
        allele_id: u32,
    },
    /// A distribution emitted a value outside the pass's valid domain.
    InvalidDistributionOutput {
        pass_name: String,
        address: String,
        value: i64,
        reason: String,
    },
    /// The current IR and plan state are inconsistent for this pass.
    InvalidPlanState { pass_name: String, reason: String },
    /// A pass produced a post-state that violates active contracts.
    ContractViolation {
        pass_name: String,
        violations: Vec<ContractViolation>,
    },
}

impl PassError {
    pub fn constraint_sampling(
        pass_name: impl Into<String>,
        address: impl Into<String>,
        reason: FilteredSampleError,
    ) -> Self {
        Self::ConstraintSampling {
            pass_name: pass_name.into(),
            address: address.into(),
            reason,
        }
    }

    pub fn missing_refdata(pass_name: impl Into<String>) -> Self {
        Self::MissingRefData {
            pass_name: pass_name.into(),
        }
    }

    pub fn missing_assignment(pass_name: impl Into<String>, segment: Segment) -> Self {
        Self::MissingAssignment {
            pass_name: pass_name.into(),
            segment,
        }
    }

    pub fn missing_allele(pass_name: impl Into<String>, segment: Segment, allele_id: u32) -> Self {
        Self::MissingAllele {
            pass_name: pass_name.into(),
            segment,
            allele_id,
        }
    }

    pub fn invalid_distribution_output(
        pass_name: impl Into<String>,
        address: impl Into<String>,
        value: i64,
        reason: impl Into<String>,
    ) -> Self {
        Self::InvalidDistributionOutput {
            pass_name: pass_name.into(),
            address: address.into(),
            value,
            reason: reason.into(),
        }
    }

    pub fn invalid_plan_state(pass_name: impl Into<String>, reason: impl Into<String>) -> Self {
        Self::InvalidPlanState {
            pass_name: pass_name.into(),
            reason: reason.into(),
        }
    }

    pub fn contract_violation(
        pass_name: impl Into<String>,
        violations: Vec<ContractViolation>,
    ) -> Self {
        Self::ContractViolation {
            pass_name: pass_name.into(),
            violations,
        }
    }

    pub fn pass_name(&self) -> &str {
        match self {
            Self::ConstraintSampling { pass_name, .. } => pass_name,
            Self::MissingRefData { pass_name } => pass_name,
            Self::MissingAssignment { pass_name, .. } => pass_name,
            Self::MissingAllele { pass_name, .. } => pass_name,
            Self::InvalidDistributionOutput { pass_name, .. } => pass_name,
            Self::InvalidPlanState { pass_name, .. } => pass_name,
            Self::ContractViolation { pass_name, .. } => pass_name,
        }
    }

    pub fn address(&self) -> &str {
        match self {
            Self::ConstraintSampling { address, .. } => address,
            Self::InvalidDistributionOutput { address, .. } => address,
            Self::MissingRefData { .. }
            | Self::MissingAssignment { .. }
            | Self::MissingAllele { .. }
            | Self::InvalidPlanState { .. }
            | Self::ContractViolation { .. } => "",
        }
    }

    pub fn constraint_reason(&self) -> Option<FilteredSampleError> {
        match self {
            Self::ConstraintSampling { reason, .. } => Some(*reason),
            _ => None,
        }
    }
}

impl fmt::Display for PassError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::ConstraintSampling {
                pass_name,
                address,
                reason,
            } => {
                let detail = match reason {
                    FilteredSampleError::SupportUnavailable => {
                        "distribution support is not enumerable"
                    }
                    FilteredSampleError::EmptyAdmissibleSupport => {
                        "no admissible candidates remain after contract filtering"
                    }
                    FilteredSampleError::InvalidFilteredSupport => {
                        "filtered support has invalid weights"
                    }
                };
                write!(
                    f,
                    "{} failed strict constrained sampling at {}: {}",
                    pass_name, address, detail
                )
            }
            Self::MissingRefData { pass_name } => {
                write!(
                    f,
                    "{} requires reference data in strict execution",
                    pass_name
                )
            }
            Self::MissingAssignment { pass_name, segment } => write!(
                f,
                "{} requires a {:?} allele assignment before execution",
                pass_name, segment
            ),
            Self::MissingAllele {
                pass_name,
                segment,
                allele_id,
            } => write!(
                f,
                "{} could not resolve {:?} allele id {} in reference data",
                pass_name, segment, allele_id
            ),
            Self::InvalidDistributionOutput {
                pass_name,
                address,
                value,
                reason,
            } => write!(
                f,
                "{} received invalid distribution output at {}: {} ({})",
                pass_name, address, value, reason
            ),
            Self::InvalidPlanState { pass_name, reason } => {
                write!(
                    f,
                    "{} cannot execute in current plan state: {}",
                    pass_name, reason
                )
            }
            Self::ContractViolation {
                pass_name,
                violations,
            } => {
                write!(f, "{} produced a contract-invalid state", pass_name)?;
                for violation in violations {
                    write!(f, "; {}: {}", violation.contract_name, violation.reason)?;
                }
                Ok(())
            }
        }
    }
}

impl std::error::Error for PassError {}

// ──────────────────────────────────────────────────────────────────
// Pass — the trait every step in the pipeline implements
// ──────────────────────────────────────────────────────────────────

/// A single step in the simulation pipeline.
///
/// Every concrete biology operation (sample-allele, trim, generate-NP,
/// mutate, corrupt, …) is a `Pass`. The runtime walks a `PassPlan`,
/// invoking each pass's `execute` once.
///
/// **Contract:** `execute` must produce a new `Simulation` value
/// derived from `sim` via the persistent `with_*` API. It must not
/// mutate `sim` (the type system enforces this — `&Simulation` is
/// shared-immutable). Random choices must go through `ctx.rng` and
/// be recorded via `ctx.trace.record(...)` so the trace stays
/// faithful.
///
/// **Forward-compatibility (D6 + D7):** the trait will grow as
/// contracts arrive in Phase D. To minimise churn for existing
/// implementations, all future additions land as **defaulted
/// methods** — concrete passes only override what's relevant. The
/// current defaulted method is `declared_choices()`; expect
/// contract-related additions later.
pub trait Pass {
    /// Stable, human-readable identifier for this pass. Used in the
    /// `History` API for `revision_after(name)` queries (D9) and in
    /// debug output. Should be a short, lowercase-with-underscores
    /// label tracking the address-namespace convention from D3.
    fn name(&self) -> &str;

    /// Apply this pass to `sim`, returning the next IR revision.
    /// Records any random choices via `ctx.trace`.
    fn execute(&self, sim: &Simulation, ctx: &mut PassContext) -> Simulation;

    /// Fallible execution hook used by strict runtime entry points.
    ///
    /// Most passes are deterministic transforms or unconstrained
    /// samplers and can inherit this default. Passes that perform
    /// contract-aware sampling override it to return a structured
    /// `PassError` when the active contract set makes a draw
    /// impossible.
    fn execute_checked(
        &self,
        sim: &Simulation,
        ctx: &mut PassContext,
    ) -> Result<Simulation, PassError> {
        Ok(self.execute(sim, ctx))
    }

    /// Address specs for the random choices this pass intends to make
    /// during a typical execution.
    ///
    /// Used by:
    /// - **Phase D — contract upstream-bound propagation:** a
    ///   contract that needs to know what addresses will be drawn at
    ///   in the future of the plan can call this on each remaining
    ///   pass to gather the address universe.
    /// - **Phase D — build-time validator (D7 Phase 1):** the
    ///   validator walks the plan, asks each pass what addresses it
    ///   declares, and checks every contract against that universe
    ///   to detect statically-impossible configurations.
    /// - **Future tooling:** trace-completeness audits, replay
    ///   shape verification.
    ///
    /// Default: `vec![]`. Transform passes (which make no random
    /// choices) keep the default. Sampling passes override and
    /// return the addresses they will record at.
    ///
    /// Address grammar matches D3: hierarchical dotted strings,
    /// e.g. `"trim.v_3"`, `"np.np1.length"`,
    /// `"mutate.s5f.site[0..n]"`. Indexed addresses use the literal
    /// `[0..n]` form to indicate variable-length expansion at
    /// execution time.
    fn declared_choices(&self) -> Vec<String> {
        Vec::new()
    }

    /// Static requirements this pass has before it can execute safely.
    ///
    /// The compiled-simulator layer uses this to reject invalid plans
    /// early (for example, assembly without refdata or without a prior
    /// sample-allele pass). Defaults to no requirements.
    fn requirements(&self) -> Vec<PassRequirement> {
        Vec::new()
    }

    /// Static effects this pass has on the simulation state.
    ///
    /// The compiled-simulator layer uses this as a typed plan summary
    /// and for coarse ordering validation. Defaults to no declared
    /// effects, which is appropriate for pure analysis/no-op passes.
    fn effects(&self) -> Vec<PassEffect> {
        Vec::new()
    }
}

// ──────────────────────────────────────────────────────────────────
// PassPlan — the ordered list of passes to execute
// ──────────────────────────────────────────────────────────────────

/// An ordered sequence of passes to run as a single simulation.
///
/// `Box<dyn Pass>` because passes are heterogeneous (different
/// concrete types, often with different internal state) but share
/// the trait. The DSL compiler builds a `PassPlan` from the user's
/// fluent builder calls; the runtime executes it.
pub struct PassPlan {
    passes: Vec<Box<dyn Pass>>,
}

impl PassPlan {
    /// Empty plan — runs no passes, leaves the initial IR unchanged.
    pub fn new() -> Self {
        Self { passes: Vec::new() }
    }

    /// Append a pass to the end of the plan.
    pub fn push(&mut self, pass: Box<dyn Pass>) -> &mut Self {
        self.passes.push(pass);
        self
    }

    /// Number of passes in the plan.
    pub fn len(&self) -> usize {
        self.passes.len()
    }

    /// Whether the plan contains zero passes.
    pub fn is_empty(&self) -> bool {
        self.passes.is_empty()
    }

    /// Read-only view of the underlying pass vector. Used by the
    /// runtime; `History` queries from D9 will also walk this.
    pub fn passes(&self) -> &[Box<dyn Pass>] {
        &self.passes
    }
}

impl Default for PassPlan {
    fn default() -> Self {
        Self::new()
    }
}

// ──────────────────────────────────────────────────────────────────
// Outcome — what the runtime returns after executing a plan
// ──────────────────────────────────────────────────────────────────

/// The result of executing a `PassPlan` on an initial `Simulation`.
///
/// `revisions` is the IR history: `revisions[0]` is the initial IR
/// (the input to the first pass), `revisions[i+1]` is the output of
/// the i-th pass. So `revisions.len() == plan.len() + 1`.
///
/// `pass_names` tracks what produced each revision *after* index 0:
/// `pass_names[i]` is the name of the pass that produced
/// `revisions[i + 1]`. So `pass_names.len() == plan.len()`.
///
/// `trace` is the addressed-choice record made during the run —
/// every choice from every pass, in chronological order.
///
/// `events` is the committed event ledger. The compiled simulator
/// emits one event for each committed pass transition. Low-level
/// `PassRuntime` entry points leave this empty because they do not
/// own the compiled transaction boundary.
#[derive(Clone, Debug)]
pub struct Outcome {
    pub revisions: Vec<Simulation>,
    pub pass_names: Vec<String>,
    pub trace: Trace,
    pub events: Vec<EventRecord>,
}

impl Outcome {
    /// The final IR revision after every pass has run. Equivalent
    /// to `revisions.last().unwrap()` but with a clearer name.
    pub fn final_simulation(&self) -> &Simulation {
        self.revisions
            .last()
            .expect("Outcome must always contain at least the initial revision")
    }

    /// First revision produced by the pass with the given name.
    /// Returns `None` if no pass with that name was in the plan.
    /// Future `History::revision_after` (D9) wraps this.
    pub fn revision_after(&self, name: &str) -> Option<&Simulation> {
        self.pass_names
            .iter()
            .position(|n| n == name)
            .map(|i| &self.revisions[i + 1])
    }

    /// The committed event ledger for this run.
    pub fn events(&self) -> &[EventRecord] {
        &self.events
    }
}

// ──────────────────────────────────────────────────────────────────
// PassRuntime — the executor
// ──────────────────────────────────────────────────────────────────

/// Executes a `PassPlan` against an initial `Simulation`, threading
/// the trace and RNG through every pass.
pub struct PassRuntime;

impl PassRuntime {
    /// Run `plan` starting from `initial`, seeded with `seed`.
    /// Reference data is **not** provided to passes — use this
    /// entry point for plans that don't consume `RefDataConfig`
    /// (sampling, trim, NP generation, the Echo / SampleBase
    /// reference passes, all tests up through C.7).
    ///
    /// For plans that include `AssembleSegmentPass` or any other
    /// reference-data-consuming pass, use `execute_with_refdata`.
    pub fn execute(plan: &PassPlan, initial: Simulation, seed: u64) -> Outcome {
        Self::execute_inner(plan, initial, seed, None, None)
    }

    /// Run `plan` with reference data threaded into every
    /// `PassContext`. Required for any plan containing assembly,
    /// SHM-against-germline, or other passes that need allele
    /// sequences from the `RefDataConfig`.
    pub fn execute_with_refdata(
        plan: &PassPlan,
        initial: Simulation,
        seed: u64,
        refdata: &RefDataConfig,
    ) -> Outcome {
        Self::execute_inner(plan, initial, seed, Some(refdata), None)
    }

    /// Run `plan` with both reference data and an active contract
    /// set. The most expressive entry point — required for any
    /// plan that uses constraint-aware sampling (D.6).
    pub fn execute_with_context(
        plan: &PassPlan,
        initial: Simulation,
        seed: u64,
        refdata: Option<&RefDataConfig>,
        contracts: Option<&ContractSet>,
    ) -> Outcome {
        Self::execute_inner(plan, initial, seed, refdata, contracts)
    }

    /// Strict counterpart to `execute_with_context`.
    ///
    /// Contract-aware passes return `Err(PassError)` instead of
    /// falling back to unconstrained sampling when a contract filter
    /// has no admissible support or cannot enumerate support.
    pub fn execute_strict_with_context(
        plan: &PassPlan,
        initial: Simulation,
        seed: u64,
        refdata: Option<&RefDataConfig>,
        contracts: Option<&ContractSet>,
    ) -> Result<Outcome, PassError> {
        Self::execute_inner_checked(plan, initial, seed, refdata, contracts)
    }

    /// Returns an `Outcome` containing the full revision history,
    /// the per-revision pass names, and the chronological trace
    /// of every random choice.
    fn execute_inner(
        plan: &PassPlan,
        initial: Simulation,
        seed: u64,
        refdata: Option<&RefDataConfig>,
        contracts: Option<&ContractSet>,
    ) -> Outcome {
        let mut trace = Trace::new();
        let mut rng = Rng::new(seed);

        let mut revisions: Vec<Simulation> = Vec::with_capacity(plan.len() + 1);
        let mut pass_names: Vec<String> = Vec::with_capacity(plan.len());

        revisions.push(initial);

        for pass in plan.passes() {
            let prev = revisions.last().expect("non-empty by construction");
            let mut ctx = PassContext {
                trace: &mut trace,
                rng: &mut rng,
                refdata,
                contracts,
            };
            let next = pass.execute(prev, &mut ctx);
            pass_names.push(pass.name().to_string());
            revisions.push(next);
        }

        Outcome {
            revisions,
            pass_names,
            trace,
            events: Vec::new(),
        }
    }

    /// Fallible runtime core used by strict execution entry points.
    fn execute_inner_checked(
        plan: &PassPlan,
        initial: Simulation,
        seed: u64,
        refdata: Option<&RefDataConfig>,
        contracts: Option<&ContractSet>,
    ) -> Result<Outcome, PassError> {
        let mut trace = Trace::new();
        let mut rng = Rng::new(seed);

        let mut revisions: Vec<Simulation> = Vec::with_capacity(plan.len() + 1);
        let mut pass_names: Vec<String> = Vec::with_capacity(plan.len());

        revisions.push(initial);

        for pass in plan.passes() {
            let prev = revisions.last().expect("non-empty by construction");
            let mut ctx = PassContext {
                trace: &mut trace,
                rng: &mut rng,
                refdata,
                contracts,
            };
            let next = pass.execute_checked(prev, &mut ctx)?;
            pass_names.push(pass.name().to_string());
            revisions.push(next);
        }

        Ok(Outcome {
            revisions,
            pass_names,
            trace,
            events: Vec::new(),
        })
    }
}

// ──────────────────────────────────────────────────────────────────
// Tests
// ──────────────────────────────────────────────────────────────────

#[cfg(test)]
mod tests {
    use super::*;
    use crate::ir::{NucHandle, Nucleotide, Segment, Simulation};
    use crate::trace::ChoiceValue;

    /// Test-only pass: pushes a single nucleotide whose base is
    /// fixed at construction time. Records its name to the trace
    /// at a fixed address so we can verify trace integration.
    struct TestEchoPass {
        base: u8,
        germline_pos: u16,
    }

    impl Pass for TestEchoPass {
        fn name(&self) -> &str {
            "test_echo"
        }
        fn execute(&self, sim: &Simulation, ctx: &mut PassContext) -> Simulation {
            // Fixed deterministic mutation for testing — push one nuc.
            let (next, _h) = sim.with_nucleotide_pushed(Nucleotide::germline(
                self.base,
                self.germline_pos,
                Segment::V,
            ));
            // Record the action to the trace so we can verify
            // trace plumbing works.
            ctx.trace
                .record("test_echo.base", ChoiceValue::Base(self.base));
            next
        }
    }

    /// Test-only pass: draws one base byte from the RNG and pushes
    /// a synthetic NP1 nucleotide with that base. Used to verify
    /// RNG plumbing and replay determinism.
    struct TestRandomBasePass;

    impl Pass for TestRandomBasePass {
        fn name(&self) -> &str {
            "test_random_base"
        }
        fn execute(&self, sim: &Simulation, ctx: &mut PassContext) -> Simulation {
            const BASES: [u8; 4] = [b'A', b'C', b'G', b'T'];
            let idx = ctx.rng.range_u32(4) as usize;
            let base = BASES[idx];
            let (next, _h) = sim.with_nucleotide_pushed(Nucleotide::synthetic(
                base,
                Segment::Np1,
                crate::ir::NucFlags::empty(),
            ));
            ctx.trace
                .record("test_random_base.base", ChoiceValue::Base(base));
            next
        }
    }

    #[test]
    fn empty_plan_returns_initial_unchanged() {
        let initial = Simulation::new();
        let plan = PassPlan::new();
        let outcome = PassRuntime::execute(&plan, initial.clone(), 42);

        assert_eq!(outcome.revisions.len(), 1);
        assert!(outcome.pass_names.is_empty());
        assert!(outcome.trace.is_empty());
        assert_eq!(outcome.final_simulation().pool.len(), 0);
    }

    #[test]
    fn single_pass_plan_produces_two_revisions() {
        let initial = Simulation::new();
        let mut plan = PassPlan::new();
        plan.push(Box::new(TestEchoPass {
            base: b'A',
            germline_pos: 0,
        }));

        let outcome = PassRuntime::execute(&plan, initial.clone(), 0);

        // Initial + one pass output = 2 revisions.
        assert_eq!(outcome.revisions.len(), 2);
        assert_eq!(outcome.pass_names, vec!["test_echo".to_string()]);

        // Initial revision unchanged.
        assert_eq!(outcome.revisions[0].pool.len(), 0);
        // Post-pass revision has the pushed nucleotide.
        assert_eq!(outcome.revisions[1].pool.len(), 1);
        assert_eq!(
            outcome.revisions[1]
                .pool
                .get(NucHandle::new(0))
                .unwrap()
                .base,
            b'A'
        );

        // Trace recorded the choice.
        assert_eq!(outcome.trace.len(), 1);
        assert_eq!(
            outcome.trace.find("test_echo.base").unwrap().value,
            ChoiceValue::Base(b'A')
        );
    }

    #[test]
    fn multi_pass_plan_applies_in_order_and_history_is_consistent() {
        let initial = Simulation::new();
        let mut plan = PassPlan::new();
        plan.push(Box::new(TestEchoPass {
            base: b'A',
            germline_pos: 0,
        }));
        plan.push(Box::new(TestEchoPass {
            base: b'C',
            germline_pos: 1,
        }));
        plan.push(Box::new(TestEchoPass {
            base: b'G',
            germline_pos: 2,
        }));

        let outcome = PassRuntime::execute(&plan, initial, 0);

        // 4 revisions: initial + 3 pass outputs.
        assert_eq!(outcome.revisions.len(), 4);
        assert_eq!(outcome.pass_names.len(), 3);

        // Each successive revision has one more nucleotide.
        assert_eq!(outcome.revisions[0].pool.len(), 0);
        assert_eq!(outcome.revisions[1].pool.len(), 1);
        assert_eq!(outcome.revisions[2].pool.len(), 2);
        assert_eq!(outcome.revisions[3].pool.len(), 3);

        // Bases pushed in order: A, C, G.
        let pool = &outcome.revisions[3].pool;
        assert_eq!(pool.get(NucHandle::new(0)).unwrap().base, b'A');
        assert_eq!(pool.get(NucHandle::new(1)).unwrap().base, b'C');
        assert_eq!(pool.get(NucHandle::new(2)).unwrap().base, b'G');

        // Earlier revisions retain their state (persistent IR).
        assert_eq!(outcome.revisions[1].pool.len(), 1);
        assert!(outcome.revisions[1].pool.get(NucHandle::new(1)).is_none());
    }

    #[test]
    fn runtime_threads_rng_through_passes_deterministically() {
        // Two separate runs with the same seed should produce the
        // exact same trace. This is the D11 reproducibility contract
        // at the runtime level. (B.5 expands the test surface to
        // assert it across the full sampling-pass machinery.)
        let mut plan_a = PassPlan::new();
        for _ in 0..10 {
            plan_a.push(Box::new(TestRandomBasePass));
        }
        let mut plan_b = PassPlan::new();
        for _ in 0..10 {
            plan_b.push(Box::new(TestRandomBasePass));
        }

        let oa = PassRuntime::execute(&plan_a, Simulation::new(), 0xc0ff_ee);
        let ob = PassRuntime::execute(&plan_b, Simulation::new(), 0xc0ff_ee);

        assert_eq!(oa.trace.len(), 10);
        assert_eq!(ob.trace.len(), 10);
        for (a, b) in oa.trace.choices().iter().zip(ob.trace.choices().iter()) {
            assert_eq!(a.address, b.address);
            assert_eq!(a.value, b.value);
        }
        // Final IRs must also agree byte for byte.
        let pa = &oa.final_simulation().pool;
        let pb = &ob.final_simulation().pool;
        assert_eq!(pa.len(), pb.len());
        for i in 0..pa.len() {
            let h = NucHandle::new(i as u32);
            assert_eq!(pa.get(h).unwrap().base, pb.get(h).unwrap().base);
        }
    }

    #[test]
    fn runtime_with_different_seeds_produces_different_traces() {
        let mut plan = PassPlan::new();
        for _ in 0..10 {
            plan.push(Box::new(TestRandomBasePass));
        }

        let o42 = PassRuntime::execute(&plan, Simulation::new(), 42);
        let o43 = PassRuntime::execute(&plan, Simulation::new(), 43);

        // With overwhelming probability at least one of the 10 random
        // bases differs between the two seeds. Asserting "any" is
        // robust against a chance collision on a single byte.
        let any_diff = o42
            .trace
            .choices()
            .iter()
            .zip(o43.trace.choices().iter())
            .any(|(a, b)| a.value != b.value);
        assert!(
            any_diff,
            "different seeds produced identical traces — RNG plumbing broken"
        );
    }

    #[test]
    fn outcome_revision_after_returns_correct_revision() {
        let mut plan = PassPlan::new();
        plan.push(Box::new(TestEchoPass {
            base: b'A',
            germline_pos: 0,
        }));
        plan.push(Box::new(TestEchoPass {
            base: b'C',
            germline_pos: 1,
        }));

        let outcome = PassRuntime::execute(&plan, Simulation::new(), 0);

        // Both passes have the same name "test_echo"; revision_after
        // returns the FIRST match, which is revisions[1] (after the
        // first echo, has 1 nucleotide).
        let after_first_echo = outcome.revision_after("test_echo").unwrap();
        assert_eq!(after_first_echo.pool.len(), 1);

        // No pass has this name.
        assert!(outcome.revision_after("does_not_exist").is_none());
    }

    #[test]
    fn pass_plan_construction_helpers() {
        let mut plan = PassPlan::new();
        assert!(plan.is_empty());
        assert_eq!(plan.len(), 0);

        plan.push(Box::new(TestEchoPass {
            base: b'A',
            germline_pos: 0,
        }));
        assert!(!plan.is_empty());
        assert_eq!(plan.len(), 1);
        assert_eq!(plan.passes()[0].name(), "test_echo");
    }
}
