//! `InvertDPass` — recombination-stage commit of a D-segment
//! orientation decision.
//!
//! Models V(D)J inversion: with probability `prob` the sampled D
//! allele is committed in [`SegmentOrientation::ReverseComplement`]
//! orientation instead of [`SegmentOrientation::Forward`]. The
//! decision is recorded under
//! [`address::ChoiceAddress::SampleAlleleDInverted`] so trace
//! replay can re-fire it deterministically. Assembly emission
//! (Slice B, [`crate::passes::AssembleSegmentPass`]) consumes the
//! orientation and emits reverse-complemented bytes for the D slice
//! when [`SegmentOrientation::ReverseComplement`] is set.
//!
//! Pipeline position: **between** `SampleAllelePass(Segment::D)`
//! and `AssembleSegmentPass(Segment::D)`. The pass declares the
//! [`PassRequirement::AlleleAssignment(Segment::D)`] dependency so
//! the schedule analyzer's auto-derivation does the rest.
//!
//! Trace addresses (Slice C):
//! - `sample_allele.d.inverted` — `Bool(true)` when the D
//!   assignment is committed reverse-complemented, `Bool(false)`
//!   otherwise. One record per simulation.
//!
//! Slice C ships the pass surface only — no public DSL method
//! yet, no AIRR `d_inverted` field yet. Those land in Slice D / E.

use crate::address;
use crate::assignment::SegmentOrientation;
use crate::ir::{Segment, Simulation, SimulationBuilder};
use crate::pass::{Pass, PassContext, PassEffect, PassError, PassRequirement};
use crate::trace::ChoiceValue;

pub struct InvertDPass {
    prob: f64,
}

impl InvertDPass {
    /// Construct a D-inversion sampling pass.
    ///
    /// Panics if `prob` is not finite or falls outside `[0.0, 1.0]`.
    /// Validation mirrors [`crate::passes::RevCompPass`] — a
    /// Bool-record pass with the same probability invariant.
    pub fn new(prob: f64) -> Self {
        assert!(
            prob.is_finite() && (0.0..=1.0).contains(&prob),
            "InvertDPass: prob must be in [0.0, 1.0], got {}",
            prob,
        );
        Self { prob }
    }

    /// Read-only accessor for tests / report builders.
    pub fn prob(&self) -> f64 {
        self.prob
    }

    fn execute_with_sampling_mode(
        &self,
        sim: &Simulation,
        ctx: &mut PassContext,
    ) -> Result<Simulation, PassError> {
        // Pre-condition: D allele must be assigned. The schedule
        // analyzer enforces this through PassRequirement, but the
        // pass also guards the checked path so a hand-built plan
        // surfaces a structured `InvalidPlanState` rather than
        // panicking inside the builder.
        if !sim.assignments.has(Segment::D) {
            return Err(PassError::missing_assignment(self.name(), Segment::D));
        }

        // Replay path consumes the recorded Bool from the cursor.
        // Mirrors `RevCompPass`: when a cursor is attached, the
        // RNG is bypassed entirely (the run is deterministic
        // against the trace).
        let inverted = if let Some(cursor) = ctx.replay_cursor.as_deref_mut() {
            cursor
                .expect_bool(address::ChoiceAddress::SampleAlleleDInverted)
                .map_err(|reason| PassError::replay(self.name(), reason))?
        } else {
            ctx.rng.next_f64() < self.prob
        };
        ctx.trace.record_choice(
            address::ChoiceAddress::SampleAlleleDInverted,
            ChoiceValue::Bool(inverted),
        );

        // Commit the orientation through a builder so any attached
        // event-log sink sees the `OrientationChanged` event. The
        // builder method asserts assignment presence — already
        // guarded above.
        let new_orientation = if inverted {
            SegmentOrientation::ReverseComplement
        } else {
            SegmentOrientation::Forward
        };

        let mut builder = SimulationBuilder::from_simulation(sim.clone());
        if ctx.event_log_sink.is_some() {
            builder.attach_event_log_observer();
        }
        builder.update_allele_orientation(Segment::D, new_orientation);
        if let Some(sink) = ctx.event_log_sink.as_deref_mut() {
            sink.extend(builder.seal_event_log_observer());
        }
        Ok(builder.seal())
    }
}

impl Pass for InvertDPass {
    fn name(&self) -> &str {
        address::INVERT_D
    }

    fn parameter_signature(&self) -> String {
        crate::passes::paramsig::fmt_prob("prob", self.prob)
    }

    fn execute(&self, sim: &Simulation, ctx: &mut PassContext) -> Simulation {
        self.execute_with_sampling_mode(sim, ctx)
            .expect("InvertDPass permissive execution must not return PassError")
    }

    fn execute_checked(
        &self,
        sim: &Simulation,
        ctx: &mut PassContext,
    ) -> Result<Simulation, PassError> {
        self.execute_with_sampling_mode(sim, ctx)
    }

    fn declared_choice_patterns(&self) -> Vec<address::ChoiceAddressPattern> {
        vec![address::ChoiceAddressPattern::SampleAlleleDInverted]
    }

    fn requirements(&self) -> Vec<PassRequirement> {
        // Requires the D allele to be assigned first. The schedule
        // analyzer auto-derives the `SampleAllelePass(D) → InvertDPass`
        // edge from this declaration.
        vec![PassRequirement::AlleleAssignment(Segment::D)]
    }

    fn effects(&self) -> Vec<PassEffect> {
        // The pass mutates `AlleleInstance.orientation` — sidecar
        // metadata, not pool bytes. No existing `PassCompileEffect`
        // variant fits, and the live-call refresh hook doesn't
        // need to react (orientation alone moves zero pool bytes
        // until AssembleSegment fires). Empty vec is consistent
        // with `RevCompPass`, which is also a Bool-record pass with
        // no compile-effect declaration.
        //
        // The execution-order edge to AssembleSegment(D) is
        // expressed structurally — the schedule analyzer already
        // orders AssembleSegment(D) after every pass that has
        // `AlleleAssignment(D)` as a requirement, so adding
        // InvertDPass between sample and assemble drops in
        // naturally without a new effect type.
        Vec::new()
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::assignment::AlleleInstance;
    use crate::ir::SimulationEvent;
    use crate::pass::{PassError, PassPlan};
    use crate::pass::testing::PassRuntime;
    use crate::refdata::AlleleId;
    use crate::replay::TraceCursor;
    use crate::rng::Rng;
    use crate::trace::Trace;

    /// Pre-stage a Simulation with D allele 0 assigned (default
    /// Forward orientation). InvertDPass is the only pass that runs
    /// in these tests — sample and assemble live elsewhere.
    fn sim_with_d_assigned() -> Simulation {
        Simulation::new()
            .with_allele_assigned(Segment::D, AlleleInstance::new(AlleleId::new(0)))
    }

    fn run_pass(prob: f64, seed: u64, sim: Simulation) -> (Trace, Simulation) {
        let mut plan = PassPlan::new();
        plan.push(Box::new(InvertDPass::new(prob)));
        let outcome = PassRuntime::execute(&plan, sim, seed);
        let final_sim = outcome.final_simulation().clone();
        (outcome.trace, final_sim)
    }

    /// Direct-context execution for replay / event tests that need
    /// to drive cursor + event sink in isolation. The builder is
    /// the same shape `PassRuntime` constructs internally; staying
    /// out of PassRuntime keeps the test focused on InvertDPass.
    fn run_with_ctx(
        pass: &InvertDPass,
        initial: Simulation,
        cursor: Option<&mut TraceCursor>,
        event_sink: Option<&mut Vec<SimulationEvent>>,
    ) -> Result<(Trace, Simulation), PassError> {
        let mut trace = Trace::new();
        let mut rng = Rng::new(0xc0ff_ee);
        let mut ctx = PassContext {
            trace: &mut trace,
            rng: &mut rng,
            pass_index: 0,
            refdata: None,
            contracts: None,
            feasibility: None,
            reference_index: None,
            replay_cursor: cursor,
            event_log_sink: event_sink,
        };
        let next = pass.execute_checked(&initial, &mut ctx)?;
        Ok((trace, next))
    }

    // ── Construction guards ──────────────────────────────────────

    #[test]
    #[should_panic(expected = "prob must be in [0.0, 1.0]")]
    fn invert_d_rejects_out_of_range_prob() {
        let _ = InvertDPass::new(1.5);
    }

    #[test]
    #[should_panic(expected = "prob must be in [0.0, 1.0]")]
    fn invert_d_rejects_negative_prob() {
        let _ = InvertDPass::new(-0.1);
    }

    #[test]
    #[should_panic(expected = "prob must be in [0.0, 1.0]")]
    fn invert_d_rejects_nan_prob() {
        let _ = InvertDPass::new(f64::NAN);
    }

    // ── Fresh-sampling Bool record ───────────────────────────────

    #[test]
    fn invert_d_prob_zero_records_false_and_keeps_forward_orientation() {
        let (trace, after) = run_pass(0.0, 0, sim_with_d_assigned());
        let rec = trace
            .find("sample_allele.d.inverted")
            .expect("Bool must be recorded even when prob = 0");
        assert_eq!(rec.value, ChoiceValue::Bool(false));
        assert_eq!(
            after.assignments.get(Segment::D).unwrap().orientation,
            SegmentOrientation::Forward,
        );
    }

    #[test]
    fn invert_d_prob_one_records_true_and_flips_to_reverse_complement() {
        let (trace, after) = run_pass(1.0, 0, sim_with_d_assigned());
        let rec = trace
            .find("sample_allele.d.inverted")
            .expect("Bool must be recorded when prob = 1");
        assert_eq!(rec.value, ChoiceValue::Bool(true));
        assert_eq!(
            after.assignments.get(Segment::D).unwrap().orientation,
            SegmentOrientation::ReverseComplement,
        );
    }

    #[test]
    fn invert_d_preserves_allele_id_and_trims() {
        // The pass touches `orientation` only. allele_id and trim
        // values must round-trip identically — Slice B's assembly
        // path reads all three and a regression here would silently
        // corrupt the post-assembly pool.
        let trimmed = AlleleInstance::new(AlleleId::new(0))
            .with_trim_5(2)
            .with_trim_3(3);
        let sim = Simulation::new().with_allele_assigned(Segment::D, trimmed);
        let (_, after) = run_pass(1.0, 0, sim);
        let inst = after.assignments.get(Segment::D).unwrap();
        assert_eq!(inst.allele_id, AlleleId::new(0));
        assert_eq!(inst.trim_5, 2);
        assert_eq!(inst.trim_3, 3);
        assert_eq!(inst.orientation, SegmentOrientation::ReverseComplement);
    }

    // ── Replay ─────────────────────────────────────────────────

    #[test]
    fn invert_d_replay_true_sets_reverse_complement_without_rng() {
        // prob=0.0 would always sample false under the RNG, but a
        // replay cursor with Bool(true) must override and commit
        // ReverseComplement. Pin the "trace is the source of truth"
        // contract.
        let pass = InvertDPass::new(0.0);
        let mut input = Trace::new();
        input.record_choice(
            address::ChoiceAddress::SampleAlleleDInverted,
            ChoiceValue::Bool(true),
        );
        let mut cursor = TraceCursor::from_owned(input.choices().to_vec());

        let (trace, after) =
            run_with_ctx(&pass, sim_with_d_assigned(), Some(&mut cursor), None).unwrap();

        assert_eq!(
            trace.find("sample_allele.d.inverted").unwrap().value,
            ChoiceValue::Bool(true),
        );
        assert_eq!(
            after.assignments.get(Segment::D).unwrap().orientation,
            SegmentOrientation::ReverseComplement,
        );
        assert!(cursor.is_drained(), "cursor should consume exactly the one Bool");
    }

    #[test]
    fn invert_d_replay_false_sets_forward_without_rng() {
        // Symmetric: prob=1.0 would always sample true under the
        // RNG, but Bool(false) in the trace must commit Forward.
        let pass = InvertDPass::new(1.0);
        let mut input = Trace::new();
        input.record_choice(
            address::ChoiceAddress::SampleAlleleDInverted,
            ChoiceValue::Bool(false),
        );
        let mut cursor = TraceCursor::from_owned(input.choices().to_vec());

        let (trace, after) =
            run_with_ctx(&pass, sim_with_d_assigned(), Some(&mut cursor), None).unwrap();

        assert_eq!(
            trace.find("sample_allele.d.inverted").unwrap().value,
            ChoiceValue::Bool(false),
        );
        assert_eq!(
            after.assignments.get(Segment::D).unwrap().orientation,
            SegmentOrientation::Forward,
        );
    }

    #[test]
    fn invert_d_replay_wrong_value_kind_surfaces_replay_error() {
        // The cursor returns an Int where Bool was expected.
        let pass = InvertDPass::new(0.5);
        let mut input = Trace::new();
        input.record_choice(
            address::ChoiceAddress::SampleAlleleDInverted,
            ChoiceValue::Int(1),
        );
        let mut cursor = TraceCursor::from_owned(input.choices().to_vec());

        let err = run_with_ctx(&pass, sim_with_d_assigned(), Some(&mut cursor), None).unwrap_err();
        match err {
            PassError::Replay { pass_name, reason } => {
                assert_eq!(pass_name, "invert_d");
                let msg = format!("{reason}");
                assert!(msg.contains("value-kind mismatch"));
                assert!(msg.contains("Bool"));
                assert!(msg.contains("Int"));
            }
            other => panic!("expected PassError::Replay, got {other:?}"),
        }
    }

    #[test]
    fn invert_d_replay_wrong_address_surfaces_replay_error() {
        // A Bool at a different address — the cursor sees a stale
        // address-string and refuses to consume it for the
        // InvertD step.
        let pass = InvertDPass::new(0.5);
        let mut input = Trace::new();
        input.record_choice(
            address::ChoiceAddress::CorruptRevCompApplied,
            ChoiceValue::Bool(true),
        );
        let mut cursor = TraceCursor::from_owned(input.choices().to_vec());

        let err = run_with_ctx(&pass, sim_with_d_assigned(), Some(&mut cursor), None).unwrap_err();
        match err {
            PassError::Replay { pass_name, reason } => {
                assert_eq!(pass_name, "invert_d");
                let msg = format!("{reason}");
                assert!(msg.contains("sample_allele.d.inverted"));
            }
            other => panic!("expected PassError::Replay, got {other:?}"),
        }
    }

    #[test]
    fn invert_d_replay_exhausted_cursor_surfaces_replay_error() {
        let pass = InvertDPass::new(0.5);
        let mut cursor = TraceCursor::from_owned(Vec::new());

        let err = run_with_ctx(&pass, sim_with_d_assigned(), Some(&mut cursor), None).unwrap_err();
        match err {
            PassError::Replay { pass_name, reason } => {
                assert_eq!(pass_name, "invert_d");
                let msg = format!("{reason}");
                assert!(msg.contains("exhausted"));
                assert!(msg.contains("sample_allele.d.inverted"));
            }
            other => panic!("expected PassError::Replay, got {other:?}"),
        }
    }

    // ── Plan-state guards (checked path) ─────────────────────────

    #[test]
    fn invert_d_returns_invalid_plan_state_when_d_unassigned_checked() {
        // No D assignment, no cursor, checked path → InvalidPlanState
        // (NOT a panic). This is the canonical "trace replay /
        // schedule analyzer would have caught this earlier" surface.
        let pass = InvertDPass::new(0.5);
        let err = run_with_ctx(&pass, Simulation::new(), None, None).unwrap_err();
        match err {
            PassError::MissingAssignment { pass_name, segment, .. } => {
                assert_eq!(pass_name, "invert_d");
                assert_eq!(segment, Segment::D);
            }
            other => panic!(
                "expected PassError::MissingAssignment, got {other:?}",
            ),
        }
    }

    #[test]
    #[should_panic(expected = "permissive execution must not return PassError")]
    fn invert_d_permissive_run_with_no_d_panics_documenting_caller_bug() {
        // Permissive `execute` propagates `InvalidPlanState` through
        // `.expect(...)`. The pin is: this is a caller bug (plan
        // built without SampleAllele(D)) and the error message has
        // to be loud. Matches the existing AssembleSegmentPass
        // convention.
        let mut plan = PassPlan::new();
        plan.push(Box::new(InvertDPass::new(0.0)));
        let _ = PassRuntime::execute(&plan, Simulation::new(), 0);
    }

    // ── Event log ────────────────────────────────────────────────

    #[test]
    fn invert_d_true_emits_orientation_changed_forward_to_reverse() {
        let pass = InvertDPass::new(1.0);
        let mut events: Vec<SimulationEvent> = Vec::new();
        let (_, _) = run_with_ctx(
            &pass,
            sim_with_d_assigned(),
            None,
            Some(&mut events),
        )
        .unwrap();
        let orient: Vec<_> = events
            .iter()
            .filter_map(|e| match e {
                SimulationEvent::OrientationChanged { segment, old, new } => {
                    Some((*segment, *old, *new))
                }
                _ => None,
            })
            .collect();
        assert_eq!(orient.len(), 1, "exactly one OrientationChanged event");
        assert_eq!(
            orient[0],
            (
                Segment::D,
                SegmentOrientation::Forward,
                SegmentOrientation::ReverseComplement,
            ),
        );
    }

    #[test]
    fn invert_d_false_from_already_rc_emits_reverse_to_forward() {
        // Pre-stage the Simulation with the D assignment already in
        // ReverseComplement (a future replay scenario or a chained
        // pipeline). prob = 0.0 commits false → Forward. The event
        // payload must surface (old=ReverseComplement, new=Forward).
        let pass = InvertDPass::new(0.0);
        let pre = Simulation::new()
            .with_allele_assigned(Segment::D, AlleleInstance::new(AlleleId::new(0)))
            .with_allele_orientation(Segment::D, SegmentOrientation::ReverseComplement);

        let mut events: Vec<SimulationEvent> = Vec::new();
        let (_, after) = run_with_ctx(&pass, pre, None, Some(&mut events)).unwrap();

        // Orientation actually flipped back.
        assert_eq!(
            after.assignments.get(Segment::D).unwrap().orientation,
            SegmentOrientation::Forward,
        );

        let orient: Vec<_> = events
            .iter()
            .filter_map(|e| match e {
                SimulationEvent::OrientationChanged { segment, old, new } => {
                    Some((*segment, *old, *new))
                }
                _ => None,
            })
            .collect();
        assert_eq!(orient.len(), 1);
        assert_eq!(
            orient[0],
            (
                Segment::D,
                SegmentOrientation::ReverseComplement,
                SegmentOrientation::Forward,
            ),
        );
    }

    #[test]
    fn invert_d_false_from_forward_still_emits_event() {
        // Even when orientation doesn't change (Forward → Forward),
        // the event fires. Sinks that count commits must see one
        // entry per pass execution; suppressing the event on a
        // no-op would require sinks to special-case the pass.
        let pass = InvertDPass::new(0.0);
        let mut events: Vec<SimulationEvent> = Vec::new();
        let _ = run_with_ctx(
            &pass,
            sim_with_d_assigned(),
            None,
            Some(&mut events),
        )
        .unwrap();
        let orient = events
            .iter()
            .filter(|e| matches!(e, SimulationEvent::OrientationChanged { .. }))
            .count();
        assert_eq!(orient, 1, "OrientationChanged fires even for no-op commits");
    }

    // ── Pass metadata ────────────────────────────────────────────

    #[test]
    fn invert_d_declares_requirement_on_d_assignment() {
        let pass = InvertDPass::new(0.5);
        assert_eq!(
            pass.requirements(),
            vec![PassRequirement::AlleleAssignment(Segment::D)],
        );
    }

    #[test]
    fn invert_d_declares_choice_pattern() {
        let pass = InvertDPass::new(0.5);
        assert_eq!(
            pass.declared_choice_patterns(),
            vec![address::ChoiceAddressPattern::SampleAlleleDInverted],
        );
    }

    #[test]
    fn invert_d_name_is_stable() {
        // The frozen pass-name string is referenced by
        // `pass_plan_signature` — any rename here breaks every
        // existing trace file. Pin it.
        let pass = InvertDPass::new(0.5);
        assert_eq!(pass.name(), "invert_d");
    }

    #[test]
    fn invert_d_declares_no_compile_effects() {
        // The orientation field is sidecar metadata — assembly
        // reads it but no live-call observer needs to react. Pin
        // the empty declaration so a future Slice E (AIRR field
        // wiring) that wants to add a compile effect surfaces here
        // first.
        let pass = InvertDPass::new(0.5);
        assert!(pass.effects().is_empty());
    }
}
