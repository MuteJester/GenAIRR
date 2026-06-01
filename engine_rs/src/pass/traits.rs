use crate::address::ChoiceAddressPattern;
use crate::ir::Simulation;

use super::{PassCompileFact, PassContext, PassEffect, PassError, PassRequirement};

/// A single step in the simulation pipeline.
pub trait Pass {
    /// Stable, human-readable identifier for this pass.
    fn name(&self) -> &str;

    /// Deterministic, replay-safety-relevant compile-time parameter
    /// digest for this pass.
    ///
    /// The string is folded into [`crate::trace_file::pass_plan_signature`]
    /// (Slice A — "Pass Parameter Signature") so that a trace
    /// recorded with one set of pass parameters (rates,
    /// distributions, probabilities, kernel identity, …) refuses to
    /// replay against a plan with mismatched parameters. Without
    /// this hook the plan signature is name-only and a parameter
    /// change silently produces different output at the same
    /// recorded addresses.
    ///
    /// Implementation rules:
    /// - Return value must be **deterministic** across runs (no
    ///   addresses, no timestamps, no `HashMap` iteration order).
    /// - Include **only compile-time parameters that affect
    ///   proposal support or the per-address output distribution**.
    ///   Names, plan position, and runtime state are NOT here.
    /// - Behaviourally-equivalent inputs must produce equal
    ///   strings — e.g. the default rate vector and an explicit
    ///   all-ones rate vector both serialise to `""` for an
    ///   `is_default()` short-circuit, OR to the same canonical
    ///   string. See [`crate::passes::paramsig`] for shared
    ///   formatting helpers.
    /// - Default implementation returns `""` (the pass has no
    ///   compile-time parameters or chooses not to participate in
    ///   the signature). Passes with parameters MUST override.
    fn parameter_signature(&self) -> String {
        String::new()
    }

    /// Apply this pass to `sim`, returning the next IR revision.
    fn execute(&self, sim: &Simulation, ctx: &mut PassContext) -> Simulation;

    /// Fallible execution hook used by strict runtime entry points.
    fn execute_checked(
        &self,
        sim: &Simulation,
        ctx: &mut PassContext,
    ) -> Result<Simulation, PassError> {
        Ok(self.execute(sim, ctx))
    }

    /// Address specs for the random choices this pass intends to make.
    fn declared_choices(&self) -> Vec<String> {
        self.declared_choice_patterns()
            .into_iter()
            .map(String::from)
            .collect()
    }

    /// Typed address specs for the random choices this pass intends
    /// to make.
    ///
    /// Built-in stochastic passes override this method and the
    /// compile/report path projects these typed patterns to the
    /// stable public string surface. Custom-address passes that do
    /// not use the built-in address vocabulary may still override
    /// [`Self::declared_choices`] directly, but they will not appear
    /// in typed compile-report metadata until they also expose a
    /// typed pattern.
    fn declared_choice_patterns(&self) -> Vec<ChoiceAddressPattern> {
        Vec::new()
    }

    /// Static requirements this pass has before it can execute safely.
    fn requirements(&self) -> Vec<PassRequirement> {
        Vec::new()
    }

    /// Static effects this pass has on the simulation state.
    fn effects(&self) -> Vec<PassEffect> {
        Vec::new()
    }

    /// Typed semantic facts available to the compile-time validator.
    fn compile_facts(&self) -> Vec<PassCompileFact> {
        Vec::new()
    }
}
