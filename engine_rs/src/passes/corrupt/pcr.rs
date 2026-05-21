//! `PCRErrorPass` — observation-stage PCR amplification errors (E.4).

use crate::address;
use crate::dist::Distribution;
use crate::ir::{NucHandle, Simulation};
use crate::pass::{Pass, PassContext, PassEffect, PassError};
use crate::passes::constrained::{sample_targeted_base, TargetedBaseChoice};
use crate::trace::ChoiceValue;

/// Models PCR amplification errors as a small number of random
/// base substitutions across the sequence. Conceptually applied at
/// observation stage (after biology) but mechanically identical to
/// `UniformMutationPass` — a count + per-mutation (site, base) draw
/// applied through `with_base_changed`.
///
/// The biological distinction (PCR errors vs SHM) is preserved
/// through:
/// - Pass `name()` is `"corrupt.pcr"` (vs `"mutate.uniform"`)
/// - Trace addresses use `corrupt.pcr.*` prefix
/// - The count distribution typically encodes a much lower rate
///   (~10⁻⁵ per base per cycle in real PCR, but the rate is
///   user-supplied so it could be anything)
///
/// Like all substitution passes, every `with_base_changed` call
/// triggers automatic codon-rail refresh (post-D audit fix).
/// When contracts are active, replacement bases are sampled from
/// the admissible subset for the chosen target site, so observation
/// errors cannot silently break enforced contracts when a safe
/// replacement exists.
///
/// Trace addresses (D3):
/// - `corrupt.pcr.count` — total PCR errors applied
/// - `corrupt.pcr.error_site[i]` — pool position of the i-th error
/// - `corrupt.pcr.error_base[i]` — replacement base at the i-th error
pub struct PCRErrorPass {
    count_dist: Box<dyn Distribution<Output = i64>>,
    base_dist: Box<dyn Distribution<Output = u8>>,
}

impl PCRErrorPass {
    pub fn new(
        count_dist: Box<dyn Distribution<Output = i64>>,
        base_dist: Box<dyn Distribution<Output = u8>>,
    ) -> Self {
        Self {
            count_dist,
            base_dist,
        }
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
                address::CORRUPT_PCR_COUNT,
                count_raw,
                "negative_count",
            ));
        }
        if strict && count_raw > u32::MAX as i64 {
            return Err(PassError::invalid_distribution_output(
                self.name(),
                address::CORRUPT_PCR_COUNT,
                count_raw,
                "count_exceeds_u32",
            ));
        }
        assert!(
            count_raw >= 0,
            "PCRErrorPass: count distribution returned negative {}",
            count_raw
        );
        assert!(
            count_raw <= u32::MAX as i64,
            "PCRErrorPass: count distribution returned {} > u32::MAX",
            count_raw
        );
        let count = count_raw as u32;
        ctx.trace
            .record(address::CORRUPT_PCR_COUNT, ChoiceValue::Int(count_raw));

        let pool_len = sim.pool.len() as u32;
        if pool_len == 0 || count == 0 {
            return Ok(sim.clone());
        }

        let mut current = sim.clone();
        for i in 0..count {
            let site = ctx.rng.range_u32(pool_len);
            let site_handle = NucHandle::new(site);
            let base_address = address::corrupt_pcr_base(i);

            ctx.trace
                .record(address::corrupt_pcr_site(i), ChoiceValue::Int(site as i64));
            let new_base = sample_targeted_base(
                &current,
                ctx,
                self.base_dist.as_ref(),
                TargetedBaseChoice::new(self.name(), &base_address, i, count, site_handle, strict),
            )?;
            ctx.trace.record(base_address, ChoiceValue::Base(new_base));

            current = current.with_base_changed(site_handle, new_base);
        }

        Ok(current)
    }
}

impl Pass for PCRErrorPass {
    fn name(&self) -> &str {
        address::CORRUPT_PCR
    }

    fn execute(&self, sim: &Simulation, ctx: &mut PassContext) -> Simulation {
        self.execute_with_sampling_mode(sim, ctx, false)
            .expect("PCRErrorPass permissive execution must not return PassError")
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
            address::CORRUPT_PCR_COUNT.to_string(),
            address::CORRUPT_PCR_SITE_PATTERN.to_string(),
            address::CORRUPT_PCR_BASE_PATTERN.to_string(),
        ]
    }

    fn effects(&self) -> Vec<PassEffect> {
        vec![PassEffect::EditBases]
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::contract::productive;
    use crate::dist::{EmpiricalLengthDist, FilteredSampleError, UniformBase};
    use crate::ir::{Nucleotide, Region, Segment};
    use crate::pass::{PassPlan, PassRuntime};
    use crate::passes::test_support::{
        make_substitution_productive_vj_fixture, StopOnlyMutationBaseDist,
        StopThenSafeMutationBaseDist,
    };

    fn pcr_test_sim() -> Simulation {
        let mut sim = Simulation::new();
        for (i, b) in b"AAACCCGGGTTTAAACCCGGG".iter().enumerate() {
            let (next, _) =
                sim.with_nucleotide_pushed(Nucleotide::germline(*b, i as u16, Segment::V));
            sim = next;
        }
        let region = Region::new(Segment::V, NucHandle::new(0), NucHandle::new(21))
            .with_codon_rail_recomputed(&sim.pool);
        sim.with_region_added(region)
    }

    fn pcr_single_error_plan(base_dist: Box<dyn Distribution<Output = u8>>) -> PassPlan {
        let mut plan = PassPlan::new();
        plan.push(Box::new(PCRErrorPass::new(
            Box::new(EmpiricalLengthDist::from_pairs(vec![(1, 1.0)])),
            base_dist,
        )));
        plan
    }

    fn pcr_error_site(seed: u64, sim: Simulation) -> u32 {
        let outcome = PassRuntime::execute(
            &pcr_single_error_plan(Box::new(StopOnlyMutationBaseDist)),
            sim,
            seed,
        );
        match outcome
            .trace
            .find("corrupt.pcr.error_site[0]")
            .unwrap()
            .value
        {
            ChoiceValue::Int(site) => site as u32,
            _ => panic!("wrong variant"),
        }
    }

    fn find_seed_for_pcr_error_site(sim: &Simulation, target_site: u32) -> u64 {
        for seed in 0..512u64 {
            if pcr_error_site(seed, sim.clone()) == target_site {
                return seed;
            }
        }
        panic!("no seed in search range targeted PCR site {}", target_site);
    }

    #[test]
    fn pcr_error_pass_zero_count_is_noop() {
        let mut plan = PassPlan::new();
        plan.push(Box::new(PCRErrorPass::new(
            Box::new(EmpiricalLengthDist::from_pairs(vec![(0, 1.0)])),
            Box::new(UniformBase),
        )));
        let outcome = PassRuntime::execute(&plan, pcr_test_sim(), 0);
        assert_eq!(outcome.trace.len(), 1);
        assert_eq!(
            outcome.trace.find("corrupt.pcr.count").unwrap().value,
            ChoiceValue::Int(0)
        );
    }

    #[test]
    fn pcr_error_pass_applies_n_errors_at_canonical_addresses() {
        let mut plan = PassPlan::new();
        plan.push(Box::new(PCRErrorPass::new(
            Box::new(EmpiricalLengthDist::from_pairs(vec![(4, 1.0)])),
            Box::new(UniformBase),
        )));
        let outcome = PassRuntime::execute(&plan, pcr_test_sim(), 7);

        // 1 count + 4 error_site + 4 error_base = 9 records.
        assert_eq!(outcome.trace.len(), 9);
        for i in 0..4 {
            assert!(outcome
                .trace
                .find(&format!("corrupt.pcr.error_site[{}]", i))
                .is_some());
            assert!(outcome
                .trace
                .find(&format!("corrupt.pcr.error_base[{}]", i))
                .is_some());
        }
    }

    #[test]
    fn pcr_error_pass_pool_reflects_recorded_errors() {
        let mut plan = PassPlan::new();
        plan.push(Box::new(PCRErrorPass::new(
            Box::new(EmpiricalLengthDist::from_pairs(vec![(3, 1.0)])),
            Box::new(UniformBase),
        )));
        let outcome = PassRuntime::execute(&plan, pcr_test_sim(), 99);
        let final_sim = outcome.final_simulation();

        let mut last_at_site: std::collections::HashMap<u32, u8> = std::collections::HashMap::new();
        for i in 0..3 {
            let s = match outcome
                .trace
                .find(&format!("corrupt.pcr.error_site[{}]", i))
                .unwrap()
                .value
            {
                ChoiceValue::Int(n) => n as u32,
                _ => unreachable!(),
            };
            let b = match outcome
                .trace
                .find(&format!("corrupt.pcr.error_base[{}]", i))
                .unwrap()
                .value
            {
                ChoiceValue::Base(b) => b,
                _ => unreachable!(),
            };
            last_at_site.insert(s, b);
        }
        for (&site, &expected) in last_at_site.iter() {
            assert_eq!(
                final_sim.pool.get(NucHandle::new(site)).unwrap().base,
                expected
            );
        }
    }

    #[test]
    fn pcr_error_pass_codon_rail_stays_consistent() {
        // Same post-D fix invariant: stored amino_acids matches
        // a fresh recompute against the post-error pool.
        let mut plan = PassPlan::new();
        plan.push(Box::new(PCRErrorPass::new(
            Box::new(EmpiricalLengthDist::from_pairs(vec![(5, 1.0)])),
            Box::new(UniformBase),
        )));
        let outcome = PassRuntime::execute(&plan, pcr_test_sim(), 1);
        let final_sim = outcome.final_simulation();

        let stored = &final_sim.sequence.regions[0].amino_acids;
        let fresh = final_sim.sequence.regions[0].with_codon_rail_recomputed(&final_sim.pool);
        assert_eq!(stored, &fresh.amino_acids);
    }

    #[test]
    fn pcr_error_pass_is_deterministic() {
        let plan = || {
            let mut p = PassPlan::new();
            p.push(Box::new(PCRErrorPass::new(
                Box::new(EmpiricalLengthDist::from_pairs(vec![(3, 1.0)])),
                Box::new(UniformBase),
            )));
            p
        };
        let oa = PassRuntime::execute(&plan(), pcr_test_sim(), 0xc0ff_ee);
        let ob = PassRuntime::execute(&plan(), pcr_test_sim(), 0xc0ff_ee);
        assert_eq!(oa.trace.choices(), ob.trace.choices());
    }

    #[test]
    fn pcr_error_pass_declared_choices() {
        let pass = PCRErrorPass::new(
            Box::new(EmpiricalLengthDist::from_pairs(vec![(0, 1.0)])),
            Box::new(UniformBase),
        );
        let declared = pass.declared_choices();
        assert!(declared.contains(&"corrupt.pcr.count".to_string()));
        assert!(declared.contains(&"corrupt.pcr.error_site[0..n]".to_string()));
        assert!(declared.contains(&"corrupt.pcr.error_base[0..n]".to_string()));
    }

    #[test]
    fn pcr_error_productive_filters_base_that_would_create_stop() {
        let (cfg, sim) = make_substitution_productive_vj_fixture();
        let seed = find_seed_for_pcr_error_site(&sim, 2);
        let contracts = productive();

        let constrained = PassRuntime::execute_with_context(
            &pcr_single_error_plan(Box::new(StopThenSafeMutationBaseDist)),
            sim.clone(),
            seed,
            Some(&cfg),
            Some(&contracts),
        );

        assert_eq!(
            constrained
                .trace
                .find("corrupt.pcr.error_site[0]")
                .unwrap()
                .value,
            ChoiceValue::Int(2)
        );
        assert_eq!(
            constrained
                .trace
                .find("corrupt.pcr.error_base[0]")
                .unwrap()
                .value,
            ChoiceValue::Base(b'C')
        );
        assert!(contracts
            .verify(constrained.final_simulation(), Some(&cfg))
            .is_ok());

        let unconstrained = PassRuntime::execute_with_context(
            &pcr_single_error_plan(Box::new(StopThenSafeMutationBaseDist)),
            sim,
            seed,
            Some(&cfg),
            None,
        );
        assert_eq!(
            unconstrained
                .trace
                .find("corrupt.pcr.error_base[0]")
                .unwrap()
                .value,
            ChoiceValue::Base(b'A')
        );
        let violations = productive()
            .verify(unconstrained.final_simulation(), Some(&cfg))
            .unwrap_err();
        assert!(violations
            .iter()
            .any(|v| v.contract_name == "no_stop_codon_in_junction"));
    }

    #[test]
    fn pcr_error_strict_errors_when_base_filter_empty() {
        let (cfg, sim) = make_substitution_productive_vj_fixture();
        let seed = find_seed_for_pcr_error_site(&sim, 2);
        let contracts = productive();

        let err = PassRuntime::execute_strict_with_context(
            &pcr_single_error_plan(Box::new(StopOnlyMutationBaseDist)),
            sim,
            seed,
            Some(&cfg),
            Some(&contracts),
        )
        .unwrap_err();

        assert_eq!(err.pass_name(), "corrupt.pcr");
        assert_eq!(err.address(), "corrupt.pcr.error_base[0]");
        assert_eq!(
            err.constraint_reason(),
            Some(FilteredSampleError::EmptyAdmissibleSupport)
        );
    }
}
