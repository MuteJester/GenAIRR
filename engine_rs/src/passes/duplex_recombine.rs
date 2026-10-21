//! One atomic source draw from a finite, completely preflighted repair table.
use crate::address::{ChoiceAddress, ChoiceAddressPattern};
use crate::duplex_repair::model::{DuplexContextKey, FiniteDuplexModel, PreparedDuplexHistory};
use crate::ir::{Segment, Simulation, SimulationBuilder};
use crate::pass::{Pass, PassCompileFact, PassContext, PassEffect, PassError, PassRequirement};
use crate::trace::ChoiceValue;
use std::sync::Arc;

pub struct DuplexRecombinePass {
    model: Arc<FiniteDuplexModel>,
    digest: String,
}

impl DuplexRecombinePass {
    pub fn new(model: Arc<FiniteDuplexModel>) -> Result<Self, String> {
        model.validate_schema()?;
        let digest = model.digest();
        Ok(Self { model, digest })
    }

    pub fn model(&self) -> &Arc<FiniteDuplexModel> {
        &self.model
    }

    fn apply(&self, sim: &Simulation, ctx: &mut PassContext) -> Result<Simulation, PassError> {
        let fail = |reason: &str| PassError::invalid_plan_state(self.name(), reason);
        let refdata = ctx
            .refdata
            .ok_or_else(|| PassError::missing_refdata(self.name()))?;
        let context = DuplexContextKey::from_simulation(sim)
            .ok_or_else(|| fail("duplex source requires a complete sampled VDJ context"))?;
        let rows = ctx
            .feasibility
            .and_then(|f| f.duplex_history_rows())
            .ok_or_else(|| fail("duplex source requires its compiled finite support artifact"))?;
        let eligible: Vec<&PreparedDuplexHistory> = rows
            .iter()
            .filter(|prepared| prepared.context == context)
            .filter(|prepared| {
                self.model
                    .rows
                    .get(prepared.row_index)
                    .is_some_and(|row| row.weight > 0.0)
            })
            .collect();
        if eligible.is_empty() {
            return Err(fail(
                "sampled context has no admissible positive duplex history",
            ));
        }
        self.model
            .validate_weight_support(&eligible.iter().map(|row| row.row_index).collect::<Vec<_>>())
            .map_err(|reason| fail(&reason))?;
        let address = ChoiceAddress::DuplexRepairHistory;
        let replaying = ctx.replay_cursor.is_some();
        let selected = if let Some(cursor) = ctx.replay_cursor.as_deref_mut() {
            let value = cursor
                .expect_int(address)
                .map_err(|reason| PassError::replay(self.name(), reason))?;
            let index = usize::try_from(value)
                .map_err(|_| fail("replayed duplex history index is negative or out of range"))?;
            *eligible
                .iter()
                .find(|row| row.row_index == index)
                .ok_or_else(|| {
                    fail("replayed duplex history is outside the positive admissible context")
                })?
        } else {
            // Normalize only the admissible support, retaining authored row
            // identity and original local weights. Finite huge weights cannot
            // overflow the total; no history mass changes upstream draws.
            let max = eligible
                .iter()
                .map(|r| self.model.rows[r.row_index].weight)
                .fold(0.0f64, f64::max);
            let total: f64 = eligible
                .iter()
                .map(|r| self.model.rows[r.row_index].weight / max)
                .sum();
            let target = ctx.rng.next_f64() * total;
            let mut cumulative = 0.0;
            let mut selected = None;
            for row in &eligible {
                cumulative += self.model.rows[row.row_index].weight / max;
                if target < cumulative {
                    selected = Some(*row);
                    break;
                }
            }
            selected
                .ok_or_else(|| fail("finite duplex CDF failed to represent its sampled mass"))?
        };
        let selection = selected
            .seed
            .duplex_provenance
            .as_ref()
            .and_then(|p| p.founder.selection.as_ref())
            .ok_or_else(|| fail("compiled duplex history lacks selection evidence"))?;
        if selection.model_digest != self.digest
            || selection.row_index != selected.row_index
            || selection.history_id != self.model.rows[selected.row_index].id
        {
            return Err(fail(
                "compiled duplex history belongs to a different authored model",
            ));
        }
        // Replay executes the actual recorded operations through the same
        // kernel. It never substitutes another viable row or consumes RNG.
        let replayed;
        let seed = if replaying {
            replayed = self
                .model
                .prepare_row_with_digest(selected.row_index, refdata, &self.digest)
                .map_err(|reason| fail(&reason))?;
            if replayed.context != selected.context
                || replayed.seed.duplex_provenance != selected.seed.duplex_provenance
            {
                return Err(fail(
                    "replayed duplex history differs from compiled evidence",
                ));
            }
            &replayed.seed
        } else {
            &selected.seed
        };
        if let Some(contracts) = ctx.contracts {
            contracts
                .verify(seed, Some(refdata))
                .map_err(|violations| PassError::contract_violation(self.name(), violations))?;
            if contracts.requires_productive_emission()
                && crate::productivity::assess(seed, refdata).productive != Some(true)
            {
                return Err(fail(
                    "selected duplex history lacks a positive productivity assessment",
                ));
            }
        }
        let mut builder = SimulationBuilder::from_simulation(sim.clone());
        if ctx.event_log_sink.is_some() {
            builder.attach_event_log_observer();
        }
        builder
            .install_duplex_seed(seed)
            .map_err(|reason| fail(&reason))?;
        ctx.trace
            .record_choice(address, ChoiceValue::Int(selected.row_index as i64));
        if let Some(sink) = ctx.event_log_sink.as_deref_mut() {
            sink.extend(builder.seal_event_log_observer());
        }
        Ok(builder.seal())
    }
}

#[cfg(test)]
mod tests;

impl Pass for DuplexRecombinePass {
    fn family_replay_spec(&self) -> Option<crate::family_replay::plan::ReplayPassSpecV1> {
        #[allow(unused_imports)]
        use crate::family_replay::plan::{
            bits, rate_bits, subregion_bits, CountSpec, KernelSpecV1, ReplayPassSpecV1 as P,
        };
        Some(P::DuplexSource {
            model: crate::family_replay::plan::DuplexModelSpecV1::capture(&self.model)?,
        })
    }
    fn name(&self) -> &str {
        "duplex.recombine"
    }
    fn parameter_signature(&self) -> String {
        format!(
            "source=finite_authored_duplex_history_v1;law={};model_digest={}",
            crate::duplex_repair::model::FINITE_DUPLEX_LAW,
            self.digest
        )
    }
    fn execute(&self, sim: &Simulation, ctx: &mut PassContext) -> Simulation {
        self.apply(sim, ctx)
            .expect("finite duplex source requires a compiled valid context")
    }
    fn execute_checked(
        &self,
        sim: &Simulation,
        ctx: &mut PassContext,
    ) -> Result<Simulation, PassError> {
        self.apply(sim, ctx)
    }
    fn execute_with_policy(
        &self,
        sim: &Simulation,
        ctx: &mut PassContext,
        _strict: bool,
    ) -> Result<Simulation, PassError> {
        self.apply(sim, ctx)
    }
    fn requirements(&self) -> Vec<PassRequirement> {
        vec![
            PassRequirement::RefData,
            PassRequirement::AlleleAssignment(Segment::V),
            PassRequirement::AlleleAssignment(Segment::D),
            PassRequirement::AlleleAssignment(Segment::J),
        ]
    }
    fn effects(&self) -> Vec<PassEffect> {
        vec![
            PassEffect::AssembleSegment(Segment::V),
            PassEffect::AssembleSegment(Segment::D),
            PassEffect::AssembleSegment(Segment::J),
            PassEffect::AppendRegion(Segment::Np1),
            PassEffect::AppendRegion(Segment::Np2),
        ]
    }
    fn declared_choice_patterns(&self) -> Vec<ChoiceAddressPattern> {
        vec![ChoiceAddressPattern::DuplexRepairHistory]
    }
    fn compile_facts(&self) -> Vec<PassCompileFact> {
        vec![PassCompileFact::FiniteDuplexSource {
            model: self.model.clone(),
        }]
    }
}
