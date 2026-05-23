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
use crate::trace::ChoiceValue;

pub struct NCorruptionPass {
    count_dist: Box<dyn Distribution<Output = i64>>,
}

impl NCorruptionPass {
    pub fn new(count_dist: Box<dyn Distribution<Output = i64>>) -> Self {
        Self { count_dist }
    }

    fn execute_with_sampling_mode(
        &self,
        sim: &Simulation,
        ctx: &mut PassContext,
        strict: bool,
    ) -> Result<Simulation, PassError> {
        let count_raw = self.count_dist.sample(ctx.rng);
        if strict && count_raw < 0 {
            return Err(PassError::invalid_distribution_output(
                self.name(),
                address::CORRUPT_NS_COUNT,
                count_raw,
                "negative_count",
            ));
        }
        if strict && count_raw > u32::MAX as i64 {
            return Err(PassError::invalid_distribution_output(
                self.name(),
                address::CORRUPT_NS_COUNT,
                count_raw,
                "count_exceeds_u32",
            ));
        }
        assert!(
            count_raw >= 0,
            "NCorruptionPass: count distribution returned negative {}",
            count_raw
        );
        assert!(
            count_raw <= u32::MAX as i64,
            "NCorruptionPass: count distribution returned {} > u32::MAX",
            count_raw
        );
        let count = count_raw as u32;
        ctx.trace
            .record(address::CORRUPT_NS_COUNT, ChoiceValue::Int(count_raw));

        let pool_len = sim.pool.len() as u32;
        if pool_len == 0 || count == 0 {
            return Ok(sim.clone());
        }

        // Phase 8: builder-pattern port (see quality.rs / s5f.rs).
        // Phase 11: also attach walker observers when ref_index is
        // available so the post-pass walker refresh is suppressed.
        let mut builder = crate::ir::SimulationBuilder::from_simulation(sim.clone());
        builder.attach_standard_mutation_observers(ctx.reference_index);

        for i in 0..count {
            let site = ctx.rng.range_u32(pool_len);
            let site_handle = NucHandle::new(site);
            ctx.trace
                .record(address::corrupt_ns_site(i), ChoiceValue::Int(site as i64));
            builder.change_base(site_handle, b'N');
        }

        Ok(if let Some(ref_index) = ctx.reference_index {
            builder.seal_with_committed_live_calls(ref_index)
        } else {
            builder.seal_with_committed_codon_rails()
        })
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

    fn declared_choices(&self) -> Vec<String> {
        vec![
            address::CORRUPT_NS_COUNT.to_string(),
            address::CORRUPT_NS_SITE_PATTERN.to_string(),
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
    use crate::pass::{PassPlan, PassRuntime};

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
