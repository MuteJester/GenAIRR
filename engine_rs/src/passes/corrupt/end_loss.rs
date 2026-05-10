//! `EndLossPass` — strip bases from the 5' or 3' end (Phase 12.D).
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
//! Trace addresses (D3):
//! - `corrupt.end_loss.5` — `Int(actual)`: bases removed from the
//!   5' end (clamped to the pool length, so a length sample > pool
//!   length records the actual removal, not the requested one).
//! - `corrupt.end_loss.3` — `Int(actual)`: same on the 3' end.

use crate::dist::Distribution;
use crate::ir::Simulation;
use crate::pass::{Pass, PassContext, PassEffect, PassError};
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
            LossEnd::Five => "corrupt.end_loss.5",
            LossEnd::Three => "corrupt.end_loss.3",
        }
    }
}

impl Pass for EndLossPass {
    fn name(&self) -> &str {
        match self.end {
            LossEnd::Five => "corrupt.end_loss.5",
            LossEnd::Three => "corrupt.end_loss.3",
        }
    }

    fn execute(&self, sim: &Simulation, ctx: &mut PassContext) -> Simulation {
        let raw = self.length_dist.sample(ctx.rng);
        let requested = raw.max(0) as u32;
        let pool_len = sim.pool.len() as u32;
        let actual = requested.min(pool_len);
        ctx.trace
            .record(self.address(), ChoiceValue::Int(actual as i64));
        if actual == 0 {
            return sim.clone();
        }
        let mut current = sim.clone();
        match self.end {
            LossEnd::Five => {
                // Delete from pool position 0, `actual` times. Each
                // call shifts every later position left by 1, so the
                // next position 0 is the original position 1, etc.
                for _ in 0..actual {
                    current = current.with_indel_deleted(0);
                }
            }
            LossEnd::Three => {
                // Delete the last position, `actual` times. After
                // each call the new pool length is one less, so the
                // next "last" index decreases.
                for _ in 0..actual {
                    let last = current.pool.len() as u32 - 1;
                    current = current.with_indel_deleted(last);
                }
            }
        }
        current
    }

    fn execute_checked(
        &self,
        sim: &Simulation,
        ctx: &mut PassContext,
    ) -> Result<Simulation, PassError> {
        Ok(self.execute(sim, ctx))
    }

    fn declared_choices(&self) -> Vec<String> {
        vec![self.address().to_string()]
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
    use crate::pass::{PassPlan, PassRuntime};

    fn end_loss_test_sim() -> Simulation {
        let mut sim = Simulation::new();
        for (i, b) in b"AAACCCGGGTTT".iter().enumerate() {
            let (next, _) = sim.with_nucleotide_pushed(Nucleotide::germline(
                *b,
                i as u16,
                Segment::V,
            ));
            sim = next;
        }
        let region = Region::new(
            Segment::V,
            crate::ir::NucHandle::new(0),
            crate::ir::NucHandle::new(12),
        )
        .with_codon_rail_recomputed(&sim.pool);
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
        sim.pool
            .as_slice()
            .iter()
            .map(|n| n.base as char)
            .collect()
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
}
