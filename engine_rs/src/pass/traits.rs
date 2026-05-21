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
