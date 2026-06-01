//! `EndLossPass` — strip bases from the 5' or 3' end.
//!
//! Models primer-region trimming and read-end degradation, both
//! standard sources of length variation in real immune-seq data.
//! With a length sampled from `length_dist`, this pass removes that
//! many bases from either the 5' (pool start) or 3' (pool end) of
//! the assembled sequence. Region ranges and codon rails are
//! adjusted via the IR's `with_indel_deleted` primitive, so the
//! resulting `Simulation` stays well-formed and downstream live-call
//! refreshes pick up the structural change correctly.
//!
//! **Architectural shape:**
//! - Length-driven (vs the count-driven PCR / quality / indel
//!   passes, or the probability-driven contaminant / rev-comp).
//! - Modifies the IR — the loss is permanent for downstream passes.
//!   Subsequent corruption operates on the shorter pool.
//!
//! **End-loss vs recombination trim (v3.0 semantics):**
//! End-loss is the *observation-stage* primer / read-degradation
//! artifact. It physically deletes pool bytes via
//! [`MutationTransaction::delete_base_admitting`], which routes
//! through the bundle's `admits_post_event`. It does NOT update
//! allele `trim_5` / `trim_3` metadata — that's the
//! recombination-stage `TrimPass`'s job. Under `productive()`
//! this means: if the 5'-loss tries to bite into V (or 3'-loss
//! into J), [`AnchorPreserved`] sees the live anchor codon's
//! bytes shifted relative to the reference allele and rejects.
//! Permissive mode stops trimming at that boundary (records the
//! realized count); strict mode surfaces a structured error. End-
//! loss outside V/J regions (e.g. trimming pre-V padding) is
//! unaffected. This is the architecturally honest behavior:
//! conflating end-loss with recombination trim metadata would
//! corrupt distinct biological provenance.
//!
//! Trace addresses (D3):
//! - `corrupt.end_loss.5` — `Int(actual)`: bases removed from the
//!   5' end (clamped to the pool length, so a length sample > pool
//!   length records the actual removal, not the requested one).
//! - `corrupt.end_loss.3` — `Int(actual)`: same on the 3' end.

use crate::address;
use crate::dist::Distribution;
use crate::ir::Simulation;
use crate::pass::{Pass, PassContext, PassEffect, PassError};
use crate::passes::mutation_transaction::MutationTransaction;
use crate::trace::ChoiceValue;

/// Which end the loss is applied to.
#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub enum LossEnd {
    /// 5' end — pool position 0, walk forward.
    Five,
    /// 3' end — pool position pool_len-1, walk backward.
    Three,
}

pub struct EndLossPass {
    end: LossEnd,
    length_dist: Box<dyn Distribution<Output = i64>>,
}

impl EndLossPass {
    /// Construct an end-loss pass.
    pub fn new(end: LossEnd, length_dist: Box<dyn Distribution<Output = i64>>) -> Self {
        Self { end, length_dist }
    }

    pub fn end(&self) -> LossEnd {
        self.end
    }

    fn address(&self) -> &'static str {
        match self.end {
            LossEnd::Five => address::CORRUPT_END_LOSS_5,
            LossEnd::Three => address::CORRUPT_END_LOSS_3,
        }
    }

    fn prime_end(&self) -> address::PrimeEnd {
        match self.end {
            LossEnd::Five => address::PrimeEnd::Five,
            LossEnd::Three => address::PrimeEnd::Three,
        }
    }

    fn choice_address(&self) -> address::ChoiceAddress {
        address::ChoiceAddress::CorruptEndLoss(self.prime_end())
    }

    fn execute_inner(
        &self,
        sim: &Simulation,
        ctx: &mut PassContext,
        strict: bool,
    ) -> Result<Simulation, PassError> {
        // Trace-injected replay (Tier 3): consume the recorded
        // *realized* count (the fresh-RNG path records `applied`,
        // not the raw distribution draw). Replay then re-applies
        // exactly that many deletions through the same
        // `delete_base_admitting` path. If any deletion is rejected
        // mid-loop by the current contract bundle, the realized
        // count would diverge from the recorded count — surface
        // that as a structured error rather than silently apply
        // fewer deletions.
        let replay_count = if let Some(cursor) = ctx.replay_cursor.as_deref_mut() {
            let n = cursor
                .expect_int(self.choice_address())
                .map_err(|reason| PassError::replay(self.name(), reason))?;
            if n < 0 {
                return Err(PassError::invalid_plan_state(
                    self.name().to_string(),
                    format!(
                        "end_loss: replayed count {} is negative; trace records realized counts",
                        n
                    ),
                ));
            }
            Some(n as u32)
        } else {
            None
        };

        let pool_len = sim.pool.len() as u32;
        let target = match replay_count {
            Some(n) => {
                if n > pool_len {
                    return Err(PassError::invalid_plan_state(
                        self.name().to_string(),
                        format!(
                            "end_loss: replayed count {} exceeds pool length {}",
                            n, pool_len
                        ),
                    ));
                }
                n
            }
            None => {
                let raw = self.length_dist.sample(ctx.rng);
                let requested = raw.max(0) as u32;
                requested.min(pool_len)
            }
        };

        if target == 0 {
            ctx.trace
                .record_choice(self.choice_address(), ChoiceValue::Int(0));
            return Ok(sim.clone());
        }

        // Route deletes through MutationTransaction so each removal
        // emits `on_indel_deleted` to attached observers and the
        // post-pass live-call refresh hook absorbs the structural
        // change correctly. Use the contract-aware `delete_base_admitting`
        // path so an active `productive_only` contract can stop the
        // trim before it crosses a V/J anchor — without this, primer
        // trim that ate into the anchor codon would silently produce
        // non-productive output downstream.
        let mut tx = MutationTransaction::open(sim, ctx, self.name(), strict);

        let address = self.choice_address();
        let mut applied: u32 = 0;
        let in_replay = replay_count.is_some();
        match self.end {
            LossEnd::Five => {
                // Delete pool position 0 repeatedly. The contract
                // reduces the effective trim count when a candidate
                // deletion would cross the V anchor.
                for _ in 0..target {
                    let did_delete = tx.delete_base_admitting(0, address)?;
                    if !did_delete {
                        if in_replay {
                            return Err(PassError::constraint_sampling(
                                self.name(),
                                self.address(),
                                crate::dist::FilteredSampleError::EmptyAdmissibleSupport,
                            ));
                        }
                        break;
                    }
                    applied += 1;
                }
            }
            LossEnd::Three => {
                // Delete the last position repeatedly. After each
                // call the new pool length is one less, so the next
                // "last" index decreases.
                for _ in 0..target {
                    let cur_len = tx.peek().pool.len() as u32;
                    if cur_len == 0 {
                        if in_replay {
                            return Err(PassError::invalid_plan_state(
                                self.name().to_string(),
                                "end_loss: replayed count exceeds available pool mid-loop"
                                    .to_string(),
                            ));
                        }
                        break;
                    }
                    let last = cur_len - 1;
                    let did_delete = tx.delete_base_admitting(last, address)?;
                    if !did_delete {
                        if in_replay {
                            return Err(PassError::constraint_sampling(
                                self.name(),
                                self.address(),
                                crate::dist::FilteredSampleError::EmptyAdmissibleSupport,
                            ));
                        }
                        break;
                    }
                    applied += 1;
                }
            }
        }

        // Trace records the *realized* trim count. Under
        // `productive_only`, the realized count clamps to the
        // contract-safe maximum. Under replay, `applied == target`
        // is guaranteed by the loop body (any premature stop has
        // already errored above), so the recorded count round-trips
        // verbatim.
        tx.trace()
            .record_choice(self.choice_address(), ChoiceValue::Int(applied as i64));

        tx.commit()
    }
}

impl Pass for EndLossPass {
    fn name(&self) -> &str {
        self.address()
    }

    fn parameter_signature(&self) -> String {
        // The end (5' vs 3') is encoded in `name()`
        // (`corrupt.end_loss.5` / `.3`); pin only the length
        // distribution.
        format!(
            "length={}",
            crate::passes::paramsig::fmt_int_dist(self.length_dist.as_ref())
        )
    }

    fn execute(&self, sim: &Simulation, ctx: &mut PassContext) -> Simulation {
        self.execute_inner(sim, ctx, false)
            .expect("EndLossPass permissive execution must not return PassError")
    }

    fn execute_checked(
        &self,
        sim: &Simulation,
        ctx: &mut PassContext,
    ) -> Result<Simulation, PassError> {
        self.execute_inner(sim, ctx, true)
    }

    fn declared_choice_patterns(&self) -> Vec<address::ChoiceAddressPattern> {
        vec![address::ChoiceAddressPattern::CorruptEndLoss(
            self.prime_end(),
        )]
    }

    fn effects(&self) -> Vec<PassEffect> {
        // The loss permanently changes pool length and shifts every
        // surviving region — same shape as a corruption-stage indel
        // pass. Declaring `StructuralIndel` triggers the live-call
        // refresh hook, so V/D/J calls reflect the post-loss bounds.
        vec![PassEffect::StructuralIndel]
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::dist::EmpiricalLengthDist;
    use crate::ir::{Nucleotide, Region, Segment};
    use crate::pass::testing::PassRuntime;
    use crate::pass::PassPlan;

    fn end_loss_test_sim() -> Simulation {
        let mut sim = Simulation::new();
        for (i, b) in b"AAACCCGGGTTT".iter().enumerate() {
            let (next, _) =
                sim.with_nucleotide_pushed(Nucleotide::germline(*b, i as u16, Segment::V));
            sim = next;
        }
        let region = Region::new(
            Segment::V,
            crate::ir::NucHandle::new(0),
            crate::ir::NucHandle::new(12),
        );
        sim.with_region_added(region)
    }

    fn run_loss(end: LossEnd, length: i64) -> Simulation {
        let initial = end_loss_test_sim();
        let mut plan = PassPlan::new();
        plan.push(Box::new(EndLossPass::new(
            end,
            Box::new(EmpiricalLengthDist::from_pairs(vec![(length, 1.0)])),
        )));
        let outcome = PassRuntime::execute(&plan, initial, 0);
        outcome.final_simulation().clone()
    }

    fn pool_string(sim: &Simulation) -> String {
        sim.pool.as_slice().iter().map(|n| n.base as char).collect()
    }

    #[test]
    fn five_prime_loss_removes_leading_bases() {
        let result = run_loss(LossEnd::Five, 3);
        assert_eq!(pool_string(&result), "CCCGGGTTT");
    }

    #[test]
    fn three_prime_loss_removes_trailing_bases() {
        let result = run_loss(LossEnd::Three, 4);
        assert_eq!(pool_string(&result), "AAACCCGG");
    }

    #[test]
    fn loss_clamps_to_pool_length() {
        // Asking for 100 bases of loss on a 12-base pool drops the
        // entire pool. The trace records the clamped value, not 100.
        let initial = end_loss_test_sim();
        let mut plan = PassPlan::new();
        plan.push(Box::new(EndLossPass::new(
            LossEnd::Five,
            Box::new(EmpiricalLengthDist::from_pairs(vec![(100, 1.0)])),
        )));
        let outcome = PassRuntime::execute(&plan, initial, 0);
        let final_sim = outcome.final_simulation();
        assert_eq!(pool_string(final_sim), "");
        let rec = outcome
            .trace
            .find("corrupt.end_loss.5")
            .expect("trace must record the clamped length");
        assert_eq!(rec.value, ChoiceValue::Int(12));
    }

    #[test]
    fn loss_zero_is_no_op() {
        let result = run_loss(LossEnd::Five, 0);
        assert_eq!(pool_string(&result), "AAACCCGGGTTT");
    }

    #[test]
    fn loss_records_actual_in_trace() {
        let initial = end_loss_test_sim();
        let mut plan = PassPlan::new();
        plan.push(Box::new(EndLossPass::new(
            LossEnd::Three,
            Box::new(EmpiricalLengthDist::from_pairs(vec![(5, 1.0)])),
        )));
        let outcome = PassRuntime::execute(&plan, initial, 0);
        let rec = outcome
            .trace
            .find("corrupt.end_loss.3")
            .expect("trace must record the actual length");
        assert_eq!(rec.value, ChoiceValue::Int(5));
    }

    #[test]
    fn region_shrinks_with_five_prime_loss() {
        let result = run_loss(LossEnd::Five, 5);
        let region = result
            .sequence
            .regions
            .iter()
            .find(|r| r.segment == Segment::V)
            .expect("V region must survive partial 5' loss");
        // 12 - 5 = 7 bases left; region should reflect the new pool.
        assert_eq!(region.start.index(), 0);
        assert_eq!(region.end.index(), 7);
    }

    // ── Trace-injected replay (Tier 3) ─────────────────────────────

    use crate::address::ChoiceAddress;
    use crate::pass::PassContext;
    use crate::replay::{ReplayError, TraceCursor};
    use crate::rng::Rng;
    use crate::trace::{ChoiceRecord, Trace};

    fn rec(addr: ChoiceAddress, v: ChoiceValue) -> ChoiceRecord {
        ChoiceRecord::new(addr.to_string(), v)
    }

    fn run_end_loss_replay(
        end: LossEnd,
        records: Vec<ChoiceRecord>,
    ) -> (Result<Simulation, PassError>, Trace, u64) {
        let initial = end_loss_test_sim();
        let pass = EndLossPass::new(
            end,
            Box::new(EmpiricalLengthDist::from_pairs(vec![(7, 1.0)])),
        );
        let mut cursor = TraceCursor::from_owned(records);
        let mut trace = Trace::new();
        let mut rng = Rng::new(0xc0ff_ee);
        let result = {
            let mut ctx = PassContext {
                trace: &mut trace,
                rng: &mut rng,
                pass_index: 0,
                refdata: None,
                contracts: None,
                feasibility: None,
                reference_index: None,
                replay_cursor: Some(&mut cursor),
                event_log_sink: None,
            };
            pass.execute_checked(&initial, &mut ctx)
        };
        let words = rng.words_consumed();
        (result, trace, words)
    }

    #[test]
    fn end_loss_replay_consumes_recorded_count_and_applies_deletions() {
        // Recorded realized count = 3. Replay deletes 3 from the 5'
        // end and writes Int(3) to the output trace.
        let (result, trace, rng_words) = run_end_loss_replay(
            LossEnd::Five,
            vec![rec(
                ChoiceAddress::CorruptEndLoss(address::PrimeEnd::Five),
                ChoiceValue::Int(3),
            )],
        );
        let next = result.unwrap();
        assert_eq!(pool_string(&next), "CCCGGGTTT");
        assert_eq!(
            trace.find("corrupt.end_loss.5").unwrap().value,
            ChoiceValue::Int(3),
        );
        // RNG bypassed entirely.
        assert_eq!(rng_words, 0);
    }

    #[test]
    fn end_loss_replay_zero_count_is_noop() {
        let (result, trace, _) = run_end_loss_replay(
            LossEnd::Three,
            vec![rec(
                ChoiceAddress::CorruptEndLoss(address::PrimeEnd::Three),
                ChoiceValue::Int(0),
            )],
        );
        let next = result.unwrap();
        assert_eq!(pool_string(&next), "AAACCCGGGTTT");
        assert_eq!(
            trace.find("corrupt.end_loss.3").unwrap().value,
            ChoiceValue::Int(0),
        );
    }

    #[test]
    fn end_loss_replay_negative_count_rejects() {
        let (result, _, _) = run_end_loss_replay(
            LossEnd::Five,
            vec![rec(
                ChoiceAddress::CorruptEndLoss(address::PrimeEnd::Five),
                ChoiceValue::Int(-1),
            )],
        );
        match result.unwrap_err() {
            PassError::InvalidPlanState { reason, .. } => {
                assert!(reason.contains("negative"));
            }
            other => panic!("expected InvalidPlanState, got {other:?}"),
        }
    }

    #[test]
    fn end_loss_replay_count_exceeds_pool_rejects() {
        let (result, _, _) = run_end_loss_replay(
            LossEnd::Five,
            vec![rec(
                ChoiceAddress::CorruptEndLoss(address::PrimeEnd::Five),
                ChoiceValue::Int(99),
            )],
        );
        match result.unwrap_err() {
            PassError::InvalidPlanState { reason, .. } => {
                assert!(reason.contains("99"));
                assert!(reason.contains("exceeds pool length"));
            }
            other => panic!("expected InvalidPlanState, got {other:?}"),
        }
    }

    #[test]
    fn end_loss_replay_wrong_value_kind_surfaces_replay_error() {
        let (result, _, _) = run_end_loss_replay(
            LossEnd::Five,
            vec![rec(
                ChoiceAddress::CorruptEndLoss(address::PrimeEnd::Five),
                ChoiceValue::Base(b'A'),
            )],
        );
        match result.unwrap_err() {
            PassError::Replay { reason, .. } => {
                assert!(matches!(reason, ReplayError::ValueKindMismatch { .. }));
            }
            other => panic!("expected Replay::ValueKindMismatch, got {other:?}"),
        }
    }
}
