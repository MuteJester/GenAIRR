//! `NCorruptionPass` — sprinkle ambiguous `N` bases.
//!
//! Real sequencing reads contain `N` nucleotides at low-quality
//! positions where the base caller could not commit to a definite
//! call. Aligners and IUPAC-aware tooling treat `N` specially —
//! benchmarking either against GenAIRR data needs `N` to actually
//! appear in the output.
//!
//! Mechanically identical to `QualityErrorPass` minus the base-
//! sampling distribution: a count is sampled, that many random
//! pool positions are picked, and each is overwritten with `b'N'`.
//!
//! Trace addresses (D3):
//! - `corrupt.ns.count` — number of `N` substitutions applied.
//! - `corrupt.ns.site[i]` — pool position of the i-th substitution.

use crate::address;
use crate::dist::Distribution;
use crate::ir::{NucHandle, Simulation};
use crate::pass::{Pass, PassContext, PassEffect, PassError};
use crate::passes::count_source::{sample_validated_count, CountSource};
use crate::passes::mutation_transaction::MutationTransaction;
use crate::trace::ChoiceValue;

pub struct NCorruptionPass {
    count_source: CountSource,
}

impl NCorruptionPass {
    pub fn new(count_dist: Box<dyn Distribution<Output = i64>>) -> Self {
        Self {
            count_source: CountSource::Distribution(count_dist),
        }
    }

    fn execute_with_sampling_mode(
        &self,
        sim: &Simulation,
        ctx: &mut PassContext,
        strict: bool,
    ) -> Result<Simulation, PassError> {
        let pool_len = sim.pool.len() as u32;
        let count = sample_validated_count(
            &self.count_source,
            ctx,
            pool_len,
            self.name(),
            address::ChoiceAddress::CorruptNsCount,
            strict,
        )?;

        if pool_len == 0 || count == 0 {
            return Ok(sim.clone());
        }

        // N-corruption pins the candidate value to `b'N'` but lets
        // the active contract bundle arbitrate per-site. Contracts
        // that reject a pinned `N` write (for example an anchor
        // codon preservation check) cause a permissive no-op at
        // that sampled site; strict mode surfaces the rejection.
        //
        // Trace-injected replay (Tier 2 sub-slice 2): the per-step
        // site Int is recorded by this pass (the pinned `N` base
        // is not a separate trace record, so only the site needs
        // consuming). Branching here keeps `substitute_base_admitting`
        // unchanged — the helper still receives a known site and
        // the pinned byte, and runs the same contract arbitration.
        let mut tx = MutationTransaction::open(sim, ctx, self.name(), strict);
        let pass_name = self.name();

        for i in 0..count {
            let site_address = address::ChoiceAddress::CorruptNsSite(i);
            let site = if let Some(cursor) = tx.replay_cursor() {
                let s = cursor
                    .expect_int(site_address)
                    .map_err(|reason| PassError::replay(pass_name, reason))?;
                if s < 0 || s >= pool_len as i64 {
                    return Err(PassError::invalid_plan_state(
                        pass_name.to_string(),
                        format!(
                            "N-corrupt: replayed site {} out of pool range [0, {})",
                            s, pool_len
                        ),
                    ));
                }
                s as u32
            } else {
                tx.rng().range_u32(pool_len)
            };
            let site_handle = NucHandle::new(site);
            tx.trace()
                .record_choice(site_address, ChoiceValue::Int(site as i64));
            tx.substitute_base_admitting(site_handle, b'N', site_address)?;
        }

        tx.commit()
    }
}

impl Pass for NCorruptionPass {
    fn name(&self) -> &str {
        address::CORRUPT_NS
    }

    fn execute(&self, sim: &Simulation, ctx: &mut PassContext) -> Simulation {
        self.execute_with_sampling_mode(sim, ctx, false)
            .expect("NCorruptionPass permissive execution must not return PassError")
    }

    fn execute_checked(
        &self,
        sim: &Simulation,
        ctx: &mut PassContext,
    ) -> Result<Simulation, PassError> {
        self.execute_with_sampling_mode(sim, ctx, true)
    }

    fn declared_choice_patterns(&self) -> Vec<address::ChoiceAddressPattern> {
        vec![
            address::ChoiceAddressPattern::CorruptNsCount,
            address::ChoiceAddressPattern::CorruptNsSite,
        ]
    }

    fn effects(&self) -> Vec<PassEffect> {
        vec![PassEffect::EditBases]
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::dist::EmpiricalLengthDist;
    use crate::ir::{Nucleotide, Segment};
    use crate::pass::testing::PassRuntime;
    use crate::pass::PassPlan;

    fn ns_test_sim() -> Simulation {
        let mut sim = Simulation::new();
        for (i, b) in b"AAACCCGGGTTT".iter().enumerate() {
            let (next, _) =
                sim.with_nucleotide_pushed(Nucleotide::germline(*b, i as u16, Segment::V));
            sim = next;
        }
        sim
    }

    fn run_ns(count: i64) -> Simulation {
        let initial = ns_test_sim();
        let mut plan = PassPlan::new();
        plan.push(Box::new(NCorruptionPass::new(Box::new(
            EmpiricalLengthDist::from_pairs(vec![(count, 1.0)]),
        ))));
        let outcome = PassRuntime::execute(&plan, initial, 0);
        outcome.final_simulation().clone()
    }

    fn n_count(sim: &Simulation) -> usize {
        sim.pool
            .as_slice()
            .iter()
            .filter(|n| n.base == b'N')
            .count()
    }

    #[test]
    fn ncorrupt_zero_count_is_no_op() {
        let result = run_ns(0);
        assert_eq!(n_count(&result), 0);
    }

    #[test]
    fn ncorrupt_count_5_writes_5_ns() {
        let result = run_ns(5);
        assert_eq!(n_count(&result), 5);
    }

    #[test]
    fn ncorrupt_records_count_in_trace() {
        let initial = ns_test_sim();
        let mut plan = PassPlan::new();
        plan.push(Box::new(NCorruptionPass::new(Box::new(
            EmpiricalLengthDist::from_pairs(vec![(3, 1.0)]),
        ))));
        let outcome = PassRuntime::execute(&plan, initial, 0);
        let rec = outcome
            .trace
            .find("corrupt.ns.count")
            .expect("trace must record count");
        assert_eq!(rec.value, ChoiceValue::Int(3));
    }

    #[test]
    fn ncorrupt_count_can_collide_within_count() {
        // count=12 on a 12-base pool may produce repeated sites
        // (uniform sampling with replacement). The number of
        // resulting Ns is therefore <= 12, not exactly 12.
        let result = run_ns(12);
        assert!(n_count(&result) <= 12);
    }
}
