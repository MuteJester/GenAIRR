//! `PairedEndSamplingPass` — internal trace-only sampling pass
//! that records the three paired-end geometry choices.
//!
//! Slice C of the paired-end roadmap. The pass is the bridge
//! between a configured [`PairedEndLayoutSpec`] (three
//! `Distribution<Output = i64>` knobs) and the
//! [`crate::airr_record::sequence::project_paired_end_layout`]
//! kernel landed in Slice B: it samples (or, under replay,
//! consumes from the cursor) `r1_length` / `r2_length` /
//! `insert_size` and records them on the trace at:
//!
//! - `paired_end.r1_length` — `Int`
//! - `paired_end.r2_length` — `Int`
//! - `paired_end.insert_size` — `Int`
//!
//! The AIRR builder (`crate::airr_record::builder::build_airr_record`)
//! reads the three trace records when assembling a record; if all
//! three are present, the projection kernel runs and populates the
//! eight paired-end AIRR fields. Absent records (no
//! `PairedEndSamplingPass` in the plan) leave the fields at their
//! defaults — the Slice A no-layout invariant continues to hold.
//!
//! ## Architectural boundary (per the audit §5.2)
//!
//! - **No IR mutation.** The pass returns `sim` unchanged.
//! - **No `SimulationEvent` emission.** No
//!   `SimulationBuilder` is involved.
//! - **No live-call refresh.** Paired-end remains a projection-
//!   only view of the final molecule.
//!
//! Slice D will add the public `Experiment.paired_end(…)` DSL
//! and the lowering that pushes a `PairedEndSamplingPass` with
//! distributions sourced from the user's call. Until then the
//! pass is reachable only via the `pass::testing::PassRuntime`
//! harness.

use crate::address;
use crate::dist::Distribution;
use crate::ir::Simulation;
use crate::pass::{Pass, PassContext, PassError};
use crate::trace::ChoiceValue;

/// Internal layout spec. Holds the three integer distributions the
/// pass samples from to produce the per-simulation paired-end
/// geometry. Slice D wires this from `Experiment.paired_end(…)`.
pub struct PairedEndLayoutSpec {
    pub r1_length: Box<dyn Distribution<Output = i64>>,
    pub r2_length: Box<dyn Distribution<Output = i64>>,
    pub insert_size: Box<dyn Distribution<Output = i64>>,
}

impl PairedEndLayoutSpec {
    pub fn new(
        r1_length: Box<dyn Distribution<Output = i64>>,
        r2_length: Box<dyn Distribution<Output = i64>>,
        insert_size: Box<dyn Distribution<Output = i64>>,
    ) -> Self {
        Self {
            r1_length,
            r2_length,
            insert_size,
        }
    }
}

/// Trace-only sampling pass that records the three paired-end
/// geometry choices for one simulation. Returns the input
/// `Simulation` unchanged.
pub struct PairedEndSamplingPass {
    spec: PairedEndLayoutSpec,
}

impl PairedEndSamplingPass {
    pub fn new(spec: PairedEndLayoutSpec) -> Self {
        Self { spec }
    }

    fn execute_with_sampling_mode(
        &self,
        sim: &Simulation,
        ctx: &mut PassContext,
    ) -> Result<Simulation, PassError> {
        // Sample (or consume from cursor) the three values in the
        // documented order: r1_length first, r2_length next,
        // insert_size last. Replay errors at each step surface
        // immediately so the failing address is visible in the
        // diagnostic.
        let r1_length = self.sample_or_replay(
            ctx,
            address::ChoiceAddress::PairedEndR1Length,
            self.spec.r1_length.as_ref(),
        )?;
        let r2_length = self.sample_or_replay(
            ctx,
            address::ChoiceAddress::PairedEndR2Length,
            self.spec.r2_length.as_ref(),
        )?;
        let insert_size = self.sample_or_replay(
            ctx,
            address::ChoiceAddress::PairedEndInsertSize,
            self.spec.insert_size.as_ref(),
        )?;

        // Validate the inter-value relationships the projection
        // kernel will need to honour. The remaining kernel
        // constraint (`insert_size <= sequence_length`) requires
        // the record's `sequence` length and is enforced by the
        // AIRR-build path via the projection error → record
        // default → validator-out-of-bounds chain.
        self.validate_relationships(r1_length, r2_length, insert_size)?;

        // Record after validation so a failed replay leaves the
        // trace consistent with the input.
        ctx.trace.record_choice(
            address::ChoiceAddress::PairedEndR1Length,
            ChoiceValue::Int(r1_length),
        );
        ctx.trace.record_choice(
            address::ChoiceAddress::PairedEndR2Length,
            ChoiceValue::Int(r2_length),
        );
        ctx.trace.record_choice(
            address::ChoiceAddress::PairedEndInsertSize,
            ChoiceValue::Int(insert_size),
        );

        Ok(sim.clone())
    }

    fn sample_or_replay(
        &self,
        ctx: &mut PassContext,
        address: address::ChoiceAddress,
        dist: &dyn Distribution<Output = i64>,
    ) -> Result<i64, PassError> {
        if let Some(cursor) = ctx.replay_cursor.as_deref_mut() {
            cursor
                .expect_int(address)
                .map_err(|reason| PassError::replay(self.name(), reason))
        } else {
            Ok(dist.sample(ctx.rng))
        }
    }

    fn validate_relationships(
        &self,
        r1_length: i64,
        r2_length: i64,
        insert_size: i64,
    ) -> Result<(), PassError> {
        if r1_length <= 0 {
            return Err(PassError::invalid_distribution_output(
                self.name(),
                address::ChoiceAddress::PairedEndR1Length.to_string(),
                r1_length,
                "r1_length_must_be_positive",
            ));
        }
        if r2_length <= 0 {
            return Err(PassError::invalid_distribution_output(
                self.name(),
                address::ChoiceAddress::PairedEndR2Length.to_string(),
                r2_length,
                "r2_length_must_be_positive",
            ));
        }
        if insert_size < 0 {
            return Err(PassError::invalid_distribution_output(
                self.name(),
                address::ChoiceAddress::PairedEndInsertSize.to_string(),
                insert_size,
                "insert_size_must_be_non_negative",
            ));
        }
        if r1_length > insert_size {
            return Err(PassError::invalid_distribution_output(
                self.name(),
                address::ChoiceAddress::PairedEndR1Length.to_string(),
                r1_length,
                "r1_length_exceeds_insert_size",
            ));
        }
        if r2_length > insert_size {
            return Err(PassError::invalid_distribution_output(
                self.name(),
                address::ChoiceAddress::PairedEndR2Length.to_string(),
                r2_length,
                "r2_length_exceeds_insert_size",
            ));
        }
        Ok(())
    }
}

impl Pass for PairedEndSamplingPass {
    fn name(&self) -> &str {
        address::PAIRED_END
    }

    fn parameter_signature(&self) -> String {
        use crate::passes::paramsig::{fmt_int_dist, join_parts};
        join_parts([
            format!("r1={}", fmt_int_dist(self.spec.r1_length.as_ref())),
            format!("r2={}", fmt_int_dist(self.spec.r2_length.as_ref())),
            format!(
                "insert={}",
                fmt_int_dist(self.spec.insert_size.as_ref())
            ),
        ])
    }

    fn execute(&self, sim: &Simulation, ctx: &mut PassContext) -> Simulation {
        self.execute_with_sampling_mode(sim, ctx)
            .expect("PairedEndSamplingPass permissive execution must not return PassError")
    }

    fn execute_checked(
        &self,
        sim: &Simulation,
        ctx: &mut PassContext,
    ) -> Result<Simulation, PassError> {
        self.execute_with_sampling_mode(sim, ctx)
    }

    fn declared_choice_patterns(&self) -> Vec<address::ChoiceAddressPattern> {
        vec![
            address::ChoiceAddressPattern::PairedEndR1Length,
            address::ChoiceAddressPattern::PairedEndR2Length,
            address::ChoiceAddressPattern::PairedEndInsertSize,
        ]
    }

    fn effects(&self) -> Vec<crate::pass::PassEffect> {
        // Trace-only pass. The AIRR builder reads the recorded
        // values at projection time; no IR side-effect.
        Vec::new()
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::dist::EmpiricalLengthDist;
    use crate::ir::{Nucleotide, Segment};
    use crate::pass::testing::PassRuntime;
    use crate::pass::PassPlan;
    use crate::pass::PassError;
    use crate::replay::TraceCursor;
    use crate::rng::Rng;
    use crate::trace::Trace;

    fn fixed_dist(value: i64) -> Box<dyn Distribution<Output = i64>> {
        Box::new(EmpiricalLengthDist::from_pairs(vec![(value, 1.0)]))
    }

    fn empty_sim() -> Simulation {
        // Tiny pool so `Simulation::new()` doesn't trip downstream
        // length assumptions; this pass doesn't read pool bytes.
        let mut sim = Simulation::new();
        for (i, b) in b"AAAA".iter().enumerate() {
            let (next, _) = sim.with_nucleotide_pushed(Nucleotide::germline(
                *b,
                i as u16,
                Segment::V,
            ));
            sim = next;
        }
        sim
    }

    fn pass_with_constants(r1: i64, r2: i64, insert: i64) -> PairedEndSamplingPass {
        PairedEndSamplingPass::new(PairedEndLayoutSpec::new(
            fixed_dist(r1),
            fixed_dist(r2),
            fixed_dist(insert),
        ))
    }

    fn run_with_ctx(
        pass: &PairedEndSamplingPass,
        cursor: Option<&mut TraceCursor>,
    ) -> Result<Trace, PassError> {
        let mut trace = Trace::new();
        let mut rng = Rng::new(0);
        let mut ctx = PassContext {
            trace: &mut trace,
            rng: &mut rng,
            pass_index: 0,
            refdata: None,
            contracts: None,
            feasibility: None,
            reference_index: None,
            replay_cursor: cursor,
            event_log_sink: None,
        };
        pass.execute_checked(&empty_sim(), &mut ctx)?;
        Ok(trace)
    }

    // ── Fresh sampling ──────────────────────────────────────────

    #[test]
    fn fresh_run_records_three_addresses_in_documented_order() {
        let pass = pass_with_constants(150, 140, 300);
        let trace = run_with_ctx(&pass, None).unwrap();
        let recorded: Vec<&str> = trace
            .choices()
            .iter()
            .map(|c| c.address.as_str())
            .collect();
        assert_eq!(
            recorded,
            vec![
                "paired_end.r1_length",
                "paired_end.r2_length",
                "paired_end.insert_size",
            ],
            "trace records must appear in r1 → r2 → insert order"
        );
    }

    #[test]
    fn fresh_run_records_distribution_values() {
        let pass = pass_with_constants(150, 140, 300);
        let trace = run_with_ctx(&pass, None).unwrap();
        assert_eq!(
            trace.find("paired_end.r1_length").unwrap().value,
            ChoiceValue::Int(150)
        );
        assert_eq!(
            trace.find("paired_end.r2_length").unwrap().value,
            ChoiceValue::Int(140)
        );
        assert_eq!(
            trace.find("paired_end.insert_size").unwrap().value,
            ChoiceValue::Int(300)
        );
    }

    #[test]
    fn pass_returns_simulation_unchanged() {
        let pass = pass_with_constants(150, 140, 300);
        let initial = empty_sim();
        let mut plan = PassPlan::new();
        plan.push(Box::new(pass_with_constants(150, 140, 300)));
        let outcome = PassRuntime::execute(&plan, initial.clone(), 0);
        let final_sim = outcome.final_simulation();
        let initial_bases: Vec<u8> =
            initial.pool.as_slice().iter().map(|n| n.base).collect();
        let final_bases: Vec<u8> =
            final_sim.pool.as_slice().iter().map(|n| n.base).collect();
        assert_eq!(initial_bases, final_bases);
        let _ = pass;
    }

    // ── Replay ───────────────────────────────────────────────────

    #[test]
    fn replay_consumes_recorded_values_in_order() {
        let pass = pass_with_constants(1, 1, 1); // any fresh defaults
        let mut input = Trace::new();
        input.record_choice(
            address::ChoiceAddress::PairedEndR1Length,
            ChoiceValue::Int(150),
        );
        input.record_choice(
            address::ChoiceAddress::PairedEndR2Length,
            ChoiceValue::Int(140),
        );
        input.record_choice(
            address::ChoiceAddress::PairedEndInsertSize,
            ChoiceValue::Int(300),
        );
        let mut cursor = TraceCursor::from_owned(input.choices().to_vec());
        let trace = run_with_ctx(&pass, Some(&mut cursor)).unwrap();
        assert_eq!(
            trace.find("paired_end.r1_length").unwrap().value,
            ChoiceValue::Int(150),
            "replay must override the pass's fresh-sampling distribution"
        );
        assert_eq!(
            trace.find("paired_end.r2_length").unwrap().value,
            ChoiceValue::Int(140)
        );
        assert_eq!(
            trace.find("paired_end.insert_size").unwrap().value,
            ChoiceValue::Int(300)
        );
        assert!(cursor.is_drained());
    }

    #[test]
    fn replay_missing_r1_address_surfaces_replay_error() {
        let pass = pass_with_constants(150, 140, 300);
        // Cursor has a wrong-address record where r1_length should be.
        let mut input = Trace::new();
        input.record_choice(
            address::ChoiceAddress::ReceptorRevisionApplied,
            ChoiceValue::Bool(false),
        );
        let mut cursor = TraceCursor::from_owned(input.choices().to_vec());
        let err = run_with_ctx(&pass, Some(&mut cursor)).unwrap_err();
        match err {
            PassError::Replay { pass_name, reason } => {
                assert_eq!(pass_name, "paired_end");
                let msg = format!("{reason}");
                assert!(
                    msg.contains("paired_end.r1_length"),
                    "replay error must name the failing address; got {msg}",
                );
            }
            other => panic!("expected PassError::Replay, got {other:?}"),
        }
    }

    #[test]
    fn replay_missing_insert_size_at_end_surfaces_replay_error() {
        let pass = pass_with_constants(150, 140, 300);
        let mut input = Trace::new();
        input.record_choice(
            address::ChoiceAddress::PairedEndR1Length,
            ChoiceValue::Int(150),
        );
        input.record_choice(
            address::ChoiceAddress::PairedEndR2Length,
            ChoiceValue::Int(140),
        );
        // No insert_size record.
        let mut cursor = TraceCursor::from_owned(input.choices().to_vec());
        let err = run_with_ctx(&pass, Some(&mut cursor)).unwrap_err();
        match err {
            PassError::Replay { pass_name, reason } => {
                assert_eq!(pass_name, "paired_end");
                let msg = format!("{reason}");
                assert!(
                    msg.contains("paired_end.insert_size"),
                    "replay error must name the failing address; got {msg}",
                );
                assert!(msg.contains("exhausted"));
            }
            other => panic!("expected PassError::Replay, got {other:?}"),
        }
    }

    #[test]
    fn replay_wrong_value_kind_surfaces_replay_error() {
        let pass = pass_with_constants(150, 140, 300);
        let mut input = Trace::new();
        input.record_choice(
            address::ChoiceAddress::PairedEndR1Length,
            ChoiceValue::Bool(true), // wrong kind
        );
        let mut cursor = TraceCursor::from_owned(input.choices().to_vec());
        let err = run_with_ctx(&pass, Some(&mut cursor)).unwrap_err();
        match err {
            PassError::Replay { pass_name, reason } => {
                assert_eq!(pass_name, "paired_end");
                let msg = format!("{reason}");
                assert!(msg.contains("value-kind mismatch"));
                assert!(msg.contains("Int"));
                assert!(msg.contains("Bool"));
            }
            other => panic!("expected PassError::Replay, got {other:?}"),
        }
    }

    // ── Geometry validation (fresh + replay) ─────────────────────

    #[test]
    fn replay_with_insert_smaller_than_r1_length_fails_clearly() {
        let pass = pass_with_constants(150, 140, 300);
        let mut input = Trace::new();
        input.record_choice(
            address::ChoiceAddress::PairedEndR1Length,
            ChoiceValue::Int(150),
        );
        input.record_choice(
            address::ChoiceAddress::PairedEndR2Length,
            ChoiceValue::Int(140),
        );
        // insert < r1 — caught by the relationship check.
        input.record_choice(
            address::ChoiceAddress::PairedEndInsertSize,
            ChoiceValue::Int(100),
        );
        let mut cursor = TraceCursor::from_owned(input.choices().to_vec());
        let err = run_with_ctx(&pass, Some(&mut cursor)).unwrap_err();
        match err {
            PassError::InvalidDistributionOutput {
                pass_name,
                address,
                value,
                reason,
            } => {
                assert_eq!(pass_name, "paired_end");
                assert_eq!(address, "paired_end.r1_length");
                assert_eq!(value, 150);
                assert_eq!(reason, "r1_length_exceeds_insert_size");
            }
            other => {
                panic!("expected InvalidDistributionOutput for insert<r1; got {other:?}")
            }
        }
    }

    #[test]
    fn replay_with_zero_r2_length_fails_with_documented_reason() {
        let pass = pass_with_constants(150, 140, 300);
        let mut input = Trace::new();
        input.record_choice(
            address::ChoiceAddress::PairedEndR1Length,
            ChoiceValue::Int(150),
        );
        input.record_choice(
            address::ChoiceAddress::PairedEndR2Length,
            ChoiceValue::Int(0),
        );
        input.record_choice(
            address::ChoiceAddress::PairedEndInsertSize,
            ChoiceValue::Int(300),
        );
        let mut cursor = TraceCursor::from_owned(input.choices().to_vec());
        let err = run_with_ctx(&pass, Some(&mut cursor)).unwrap_err();
        match err {
            PassError::InvalidDistributionOutput {
                address,
                reason,
                ..
            } => {
                assert_eq!(address, "paired_end.r2_length");
                assert_eq!(reason, "r2_length_must_be_positive");
            }
            other => panic!("expected InvalidDistributionOutput; got {other:?}"),
        }
    }

    // ── Pass metadata ────────────────────────────────────────────

    #[test]
    fn declares_three_choice_patterns_in_documented_order() {
        let pass = pass_with_constants(150, 140, 300);
        assert_eq!(
            pass.declared_choice_patterns(),
            vec![
                address::ChoiceAddressPattern::PairedEndR1Length,
                address::ChoiceAddressPattern::PairedEndR2Length,
                address::ChoiceAddressPattern::PairedEndInsertSize,
            ]
        );
    }

    #[test]
    fn pass_name_is_paired_end() {
        let pass = pass_with_constants(150, 140, 300);
        assert_eq!(pass.name(), "paired_end");
    }

    #[test]
    fn declares_no_compile_effects() {
        // Trace-only pass; no IR effect to surface.
        let pass = pass_with_constants(150, 140, 300);
        assert!(pass.effects().is_empty());
    }

    // ── Baseline (no pass = no records) ──────────────────────────

    #[test]
    fn baseline_run_without_pass_emits_no_paired_end_records() {
        let plan = PassPlan::new();
        // Empty plan; no PairedEndSamplingPass.
        let outcome = PassRuntime::execute(&plan, empty_sim(), 0);
        for choice in outcome.trace.choices() {
            assert!(
                !choice.address.starts_with("paired_end."),
                "baseline run emitted paired_end.* address {:?}",
                choice.address,
            );
        }
    }
}
