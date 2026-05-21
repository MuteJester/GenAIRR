use crate::contract::ContractSet;
use crate::ir::Simulation;
use crate::refdata::RefDataConfig;
use crate::rng::Rng;
use crate::trace::Trace;

use super::{Outcome, PassContext, PassError, PassPlan};

/// Executes a `PassPlan` against an initial `Simulation`, threading
/// the trace and RNG through every pass.
pub struct PassRuntime;

impl PassRuntime {
    /// Run `plan` starting from `initial`, seeded with `seed`.
    pub fn execute(plan: &PassPlan, initial: Simulation, seed: u64) -> Outcome {
        Self::execute_inner(plan, initial, seed, None, None)
    }

    /// Run `plan` with reference data threaded into every `PassContext`.
    pub fn execute_with_refdata(
        plan: &PassPlan,
        initial: Simulation,
        seed: u64,
        refdata: &RefDataConfig,
    ) -> Outcome {
        Self::execute_inner(plan, initial, seed, Some(refdata), None)
    }

    /// Run `plan` with both reference data and an active contract set.
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
    pub fn execute_strict_with_context(
        plan: &PassPlan,
        initial: Simulation,
        seed: u64,
        refdata: Option<&RefDataConfig>,
        contracts: Option<&ContractSet>,
    ) -> Result<Outcome, PassError> {
        Self::execute_inner_checked(plan, initial, seed, refdata, contracts)
    }

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

        for (pass_index, pass) in plan.passes().iter().enumerate() {
            let prev = revisions.last().expect("non-empty by construction");
            let mut ctx = PassContext {
                trace: &mut trace,
                rng: &mut rng,
                pass_index,
                refdata,
                contracts,
                feasibility: None,
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

        for (pass_index, pass) in plan.passes().iter().enumerate() {
            let prev = revisions.last().expect("non-empty by construction");
            let mut ctx = PassContext {
                trace: &mut trace,
                rng: &mut rng,
                pass_index,
                refdata,
                contracts,
                feasibility: None,
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
