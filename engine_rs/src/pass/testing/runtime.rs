//! `PassRuntime` — test-only execution helper that shares the
//! production transactional loop in `compiled/execute.rs`.
//!
//! Pre-refactor `PassRuntime` had its own hand-rolled execute loop,
//! which over time drifted from `CompiledSimulator`'s loop (different
//! semantics for the live-call refresh, no effect hooks, no schedule
//! topo-sort). The audit flagged this as "two runtimes" — bugs the
//! production loop fixed could persist undetected through `PassRuntime`
//! test fixtures.
//!
//! Today `PassRuntime` is unambiguously test-only and delegates to the
//! production `execute_transactional` loop in [`crate::compiled::execute`].
//! Tests get the same per-pass commit transaction, the same `Outcome`
//! shape, and any future loop-level fixes apply to both paths.
//!
//! What `PassRuntime` does NOT do (deliberately, for test isolation):
//!
//! - Skip [`crate::compiled::analyze::analyze_plan`] — no compile-time
//!   validation of distribution supports, contract preconditions,
//!   etc. Tests can build "minimal" plans that wouldn't pass the
//!   compile gate but should still execute. Production code paths
//!   always go through `CompiledSimulator::compile`.
//! - Skip the topo-sorted [`crate::pass::Schedule::compile`] —
//!   `PassRuntime` honors the caller's insertion order verbatim.
//! - Build a [`crate::live_call::ReferenceMatchIndex`] — tests that
//!   need the index pass it explicitly or go through `CompiledSimulator`.
//! - Register [`crate::pass::EffectHook`]s — there are no
//!   post-pass derived-state refreshes. Tests asserting on live-call
//!   sidecars should use `CompiledSimulator`.

use crate::compiled::execute::{execute_transactional, ExecutionInputs};
use crate::compiled::ExecutionPolicy;
use crate::contract::ContractSet;
use crate::ir::Simulation;
use crate::refdata::RefDataConfig;

use crate::pass::{NodeId, Outcome, PassError, PassPlan};

/// Test-only execution entry point. Iterates the plan in insertion
/// order through the production transactional loop. See module docs.
pub struct PassRuntime;

impl PassRuntime {
    /// Run `plan` starting from `initial`, seeded with `seed`. No
    /// refdata, no contracts, permissive policy.
    pub fn execute(plan: &PassPlan, initial: Simulation, seed: u64) -> Outcome {
        Self::run(plan, initial, seed, None, None, ExecutionPolicy::Permissive)
            .expect("PassRuntime::execute: permissive run should not return PassError")
    }

    /// Run `plan` with reference data threaded into every
    /// `PassContext`. Permissive policy.
    pub fn execute_with_refdata(
        plan: &PassPlan,
        initial: Simulation,
        seed: u64,
        refdata: &RefDataConfig,
    ) -> Outcome {
        Self::run(
            plan,
            initial,
            seed,
            Some(refdata),
            None,
            ExecutionPolicy::Permissive,
        )
        .expect("PassRuntime::execute_with_refdata: permissive run should not return PassError")
    }

    /// Run `plan` with both reference data and an active contract
    /// set. Permissive policy.
    pub fn execute_with_context(
        plan: &PassPlan,
        initial: Simulation,
        seed: u64,
        refdata: Option<&RefDataConfig>,
        contracts: Option<&ContractSet>,
    ) -> Outcome {
        Self::run(
            plan,
            initial,
            seed,
            refdata,
            contracts,
            ExecutionPolicy::Permissive,
        )
        .expect("PassRuntime::execute_with_context: permissive run should not return PassError")
    }

    /// Strict counterpart to [`Self::execute_with_context`]. Returns
    /// the first `PassError` instead of panicking.
    pub fn execute_strict_with_context(
        plan: &PassPlan,
        initial: Simulation,
        seed: u64,
        refdata: Option<&RefDataConfig>,
        contracts: Option<&ContractSet>,
    ) -> Result<Outcome, PassError> {
        Self::run(
            plan,
            initial,
            seed,
            refdata,
            contracts,
            ExecutionPolicy::Strict,
        )
    }

    fn run(
        plan: &PassPlan,
        initial: Simulation,
        seed: u64,
        refdata: Option<&RefDataConfig>,
        contracts: Option<&ContractSet>,
        policy: ExecutionPolicy,
    ) -> Result<Outcome, PassError> {
        // Insertion-order schedule — no topo-sort. Tests that
        // intentionally test schedule-aware behaviour should go
        // through `CompiledSimulator`.
        let sorted_order: Vec<NodeId> = (0..plan.len()).map(|i| NodeId(i as u32)).collect();
        let inputs = ExecutionInputs {
            plan,
            sorted_order: &sorted_order,
            hooks: &[],
            refdata,
            contracts,
            feasibility: None,
            reference_index: None,
            policy,
            replay_records: None,
        };
        execute_transactional(inputs, initial, seed).map_err(|abort| abort.error)
    }
}
