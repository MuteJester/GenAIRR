//! `RevCompPass` — flag-recording reverse-complement pass.
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

use crate::address;
use crate::ir::{Simulation, SimulationBuilder};
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

impl RevCompPass {
    fn execute_with_sampling_mode(
        &self,
        sim: &Simulation,
        ctx: &mut PassContext,
    ) -> Result<Simulation, PassError> {
        // Trace-injected replay: consume the recorded Bool from the
        // cursor instead of flipping a coin. The downstream
        // application path is identical — record the flag, return a
        // clone of `sim` (this pass never mutates the IR).
        let applied = if let Some(cursor) = ctx.replay_cursor.as_deref_mut() {
            cursor
                .expect_bool(address::ChoiceAddress::CorruptRevCompApplied)
                .map_err(|reason| PassError::replay(self.name(), reason))?
        } else {
            ctx.rng.next_f64() < self.apply_prob
        };
        ctx.trace.record_choice(
            address::ChoiceAddress::CorruptRevCompApplied,
            ChoiceValue::Bool(applied),
        );
        // No IR mutation — the AIRR record builder reads the trace
        // flag and post-flips the sequence at projection time. But
        // we still route through a `SimulationBuilder` so the
        // `ReverseComplementFlagRecorded` event flows to any
        // attached sink, and forward into `ctx.event_log_sink`
        // when supplied. The builder's `seal()` returns the
        // unchanged simulation; the only side effect is the event
        // broadcast.
        let mut builder = SimulationBuilder::from_simulation(sim.clone());
        if ctx.event_log_sink.is_some() {
            builder.attach_event_log_observer();
        }
        builder.record_reverse_complement_flag(applied);
        if let Some(sink) = ctx.event_log_sink.as_deref_mut() {
            sink.extend(builder.seal_event_log_observer());
        }
        Ok(builder.seal())
    }
}

impl Pass for RevCompPass {
    fn name(&self) -> &str {
        address::CORRUPT_REV_COMP
    }

    fn execute(&self, sim: &Simulation, ctx: &mut PassContext) -> Simulation {
        self.execute_with_sampling_mode(sim, ctx)
            .expect("RevCompPass permissive execution must not return PassError")
    }

    fn execute_checked(
        &self,
        sim: &Simulation,
        ctx: &mut PassContext,
    ) -> Result<Simulation, PassError> {
        self.execute_with_sampling_mode(sim, ctx)
    }

    fn declared_choice_patterns(&self) -> Vec<address::ChoiceAddressPattern> {
        vec![address::ChoiceAddressPattern::CorruptRevCompApplied]
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
    use crate::pass::testing::PassRuntime;
    use crate::pass::PassPlan;
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
        let final_bases: Vec<u8> = final_sim.pool.as_slice().iter().map(|n| n.base).collect();
        assert_eq!(initial_bases, final_bases);
    }

    #[test]
    #[should_panic(expected = "apply_prob must be in [0.0, 1.0]")]
    fn rev_comp_rejects_out_of_range_prob() {
        let _ = RevCompPass::new(1.5);
    }

    // ── Trace-injected replay path (Option B Tier 1) ──────────────
    //
    // The full-plan e2e for `corrupt.rev_comp.applied` is blocked
    // by unmigrated NP base records sitting between `np.np1.length`
    // and the rev_comp record (strict-positional cursor sees the
    // NP base address when it tries to consume the Bool). These
    // direct tests construct a `PassContext` with a cursor in
    // isolation so the consume path can be validated without the
    // intermediate unmigrated records.

    use crate::pass::PassError;
    use crate::replay::TraceCursor;
    use crate::rng::Rng;

    fn run_with_replay(
        pass: &RevCompPass,
        cursor: &mut TraceCursor,
    ) -> Result<(Trace, Simulation), PassError> {
        let initial = rev_comp_test_sim();
        let mut trace = Trace::new();
        let mut rng = Rng::new(0);
        let mut ctx = crate::pass::PassContext {
            trace: &mut trace,
            rng: &mut rng,
            pass_index: 0,
            refdata: None,
            contracts: None,
            feasibility: None,
            reference_index: None,
            replay_cursor: Some(cursor),
            event_log_sink: None,
        };
        let next = pass.execute_checked(&initial, &mut ctx)?;
        Ok((trace, next))
    }

    #[test]
    fn rev_comp_replay_consumes_recorded_bool_true() {
        let pass = RevCompPass::new(0.0); // Would always be false under RNG.
        let mut input = Trace::new();
        input.record_choice(
            address::ChoiceAddress::CorruptRevCompApplied,
            ChoiceValue::Bool(true),
        );
        let mut cursor = TraceCursor::from_owned(input.choices().to_vec());

        let (trace, _next) = run_with_replay(&pass, &mut cursor).unwrap();
        // The output trace must reflect the cursor's value, not what
        // the RNG path would have produced.
        let rec = trace.find("corrupt.rev_comp.applied").unwrap();
        assert_eq!(rec.value, ChoiceValue::Bool(true));
        assert!(cursor.is_drained());
    }

    #[test]
    fn rev_comp_replay_consumes_recorded_bool_false() {
        let pass = RevCompPass::new(1.0); // Would always be true under RNG.
        let mut input = Trace::new();
        input.record_choice(
            address::ChoiceAddress::CorruptRevCompApplied,
            ChoiceValue::Bool(false),
        );
        let mut cursor = TraceCursor::from_owned(input.choices().to_vec());

        let (trace, _next) = run_with_replay(&pass, &mut cursor).unwrap();
        let rec = trace.find("corrupt.rev_comp.applied").unwrap();
        assert_eq!(rec.value, ChoiceValue::Bool(false));
    }

    #[test]
    fn rev_comp_replay_value_kind_mismatch_surfaces_replay_error() {
        let pass = RevCompPass::new(0.5);
        let mut input = Trace::new();
        // Wrong kind: Int instead of Bool.
        input.record_choice(
            address::ChoiceAddress::CorruptRevCompApplied,
            ChoiceValue::Int(1),
        );
        let mut cursor = TraceCursor::from_owned(input.choices().to_vec());

        let err = run_with_replay(&pass, &mut cursor).unwrap_err();
        match err {
            PassError::Replay { pass_name, reason } => {
                assert_eq!(pass_name, "corrupt.rev_comp");
                let msg = format!("{reason}");
                assert!(msg.contains("value-kind mismatch"));
                assert!(msg.contains("Bool"));
                assert!(msg.contains("Int"));
            }
            other => panic!("expected PassError::Replay, got {other:?}"),
        }
    }

    #[test]
    fn rev_comp_replay_exhausted_cursor_surfaces_replay_error() {
        let pass = RevCompPass::new(0.5);
        let mut cursor = TraceCursor::from_owned(Vec::new());
        let err = run_with_replay(&pass, &mut cursor).unwrap_err();
        match err {
            PassError::Replay { pass_name, reason } => {
                assert_eq!(pass_name, "corrupt.rev_comp");
                let msg = format!("{reason}");
                assert!(msg.contains("exhausted"));
                assert!(msg.contains("corrupt.rev_comp.applied"));
            }
            other => panic!("expected PassError::Replay, got {other:?}"),
        }
    }
}
