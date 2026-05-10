//! `RevCompPass` — flag-recording reverse-complement pass (Phase 12.D).
//!
//! Models the fact that ~50% of immune-seq reads come from the
//! antisense strand. With probability `apply_prob` this pass records
//! a Boolean flag in the trace; downstream record-builders (the AIRR
//! Rearrangement projection in `airr_record.rs`) then flip the
//! `sequence` and per-segment sequence-coord fields, leaving
//! alignment / germline / CIGAR fields in forward orientation per
//! the AIRR Rearrangement spec.
//!
//! The pass itself does NOT mutate the IR. Subsequent passes (if
//! any) continue to see the forward sequence, so any later
//! observation-stage corruption still operates on the forward strand.
//! This is a deliberate choice — flipping the IR mid-pipeline would
//! break every downstream pass's pool/region invariants.
//!
//! Trace addresses (D3):
//! - `corrupt.rev_comp.applied` — `Bool(true)` if the read was
//!   marked antisense, `Bool(false)` otherwise.

use crate::ir::Simulation;
use crate::pass::{Pass, PassContext, PassEffect, PassError};
use crate::trace::ChoiceValue;

pub struct RevCompPass {
    apply_prob: f64,
}

impl RevCompPass {
    /// Construct a reverse-complement flag pass.
    ///
    /// Panics if `apply_prob` is not in `[0.0, 1.0]` or is non-finite.
    pub fn new(apply_prob: f64) -> Self {
        assert!(
            apply_prob.is_finite() && (0.0..=1.0).contains(&apply_prob),
            "RevCompPass: apply_prob must be in [0.0, 1.0], got {}",
            apply_prob
        );
        Self { apply_prob }
    }

    pub fn apply_prob(&self) -> f64 {
        self.apply_prob
    }
}

impl Pass for RevCompPass {
    fn name(&self) -> &str {
        "corrupt.rev_comp"
    }

    fn execute(&self, sim: &Simulation, ctx: &mut PassContext) -> Simulation {
        let coin = ctx.rng.next_f64();
        let applied = coin < self.apply_prob;
        ctx.trace
            .record("corrupt.rev_comp.applied", ChoiceValue::Bool(applied));
        // No IR mutation — the AIRR record builder reads the trace
        // flag and post-flips the sequence at projection time.
        sim.clone()
    }

    fn execute_checked(
        &self,
        sim: &Simulation,
        ctx: &mut PassContext,
    ) -> Result<Simulation, PassError> {
        Ok(self.execute(sim, ctx))
    }

    fn declared_choices(&self) -> Vec<String> {
        vec!["corrupt.rev_comp.applied".to_string()]
    }

    fn effects(&self) -> Vec<PassEffect> {
        // The pass doesn't modify base content or region structure —
        // it's a metadata-only flag. No live-call refresh is needed.
        Vec::new()
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::ir::{Nucleotide, Segment};
    use crate::pass::{PassPlan, PassRuntime};
    use crate::trace::Trace;

    fn rev_comp_test_sim() -> Simulation {
        let mut sim = Simulation::new();
        for (i, b) in b"AAACCCGGG".iter().enumerate() {
            let (next, _) =
                sim.with_nucleotide_pushed(Nucleotide::germline(*b, i as u16, Segment::V));
            sim = next;
        }
        sim
    }

    fn run_pass(prob: f64, seed: u64) -> Trace {
        let initial = rev_comp_test_sim();
        let mut plan = PassPlan::new();
        plan.push(Box::new(RevCompPass::new(prob)));
        let outcome = PassRuntime::execute(&plan, initial, seed);
        outcome.trace
    }

    #[test]
    fn rev_comp_records_applied_flag() {
        // prob=1.0 → always applied.
        let trace = run_pass(1.0, 0);
        let rec = trace
            .find("corrupt.rev_comp.applied")
            .expect("flag must be recorded");
        assert_eq!(rec.value, ChoiceValue::Bool(true));
    }

    #[test]
    fn rev_comp_records_not_applied_flag_when_prob_zero() {
        let trace = run_pass(0.0, 0);
        let rec = trace
            .find("corrupt.rev_comp.applied")
            .expect("flag must be recorded even when not applied");
        assert_eq!(rec.value, ChoiceValue::Bool(false));
    }

    #[test]
    fn rev_comp_does_not_mutate_pool() {
        // Run with prob=1.0; the IR pool stays intact.
        let initial = rev_comp_test_sim();
        let mut plan = PassPlan::new();
        plan.push(Box::new(RevCompPass::new(1.0)));
        let outcome = PassRuntime::execute(&plan, initial.clone(), 0);
        let final_sim = outcome.final_simulation();
        let initial_bases: Vec<u8> = initial.pool.as_slice().iter().map(|n| n.base).collect();
        let final_bases: Vec<u8> = final_sim
            .pool
            .as_slice()
            .iter()
            .map(|n| n.base)
            .collect();
        assert_eq!(initial_bases, final_bases);
    }

    #[test]
    #[should_panic(expected = "apply_prob must be in [0.0, 1.0]")]
    fn rev_comp_rejects_out_of_range_prob() {
        let _ = RevCompPass::new(1.5);
    }
}
