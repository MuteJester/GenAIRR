use crate::address::ChoiceAddressPattern;
use crate::ir::Simulation;

use super::{PassCompileFact, PassContext, PassEffect, PassError, PassRequirement};

/// A single step in the simulation pipeline.
pub trait Pass {
    /// Stable, human-readable identifier for this pass.
    fn name(&self) -> &str;

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
