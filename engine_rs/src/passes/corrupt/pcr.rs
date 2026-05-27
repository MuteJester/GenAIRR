//! `PCRErrorPass` — observation-stage PCR amplification errors (E.4).

use crate::address;
use crate::dist::Distribution;
#[cfg(test)]
use crate::ir::NucHandle;
use crate::ir::Simulation;
use crate::pass::{Pass, PassContext, PassEffect, PassError};
use crate::passes::count_source::sample_validated_count;
use crate::passes::mutation_transaction::MutationTransaction;
#[cfg(test)]
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
/// Like all substitution passes, the pool is the authoritative
/// source — codon-rail data is not stored on `Region`. When
/// contracts are active, replacement bases are sampled from the
/// admissible subset for the chosen target site, so observation
/// errors cannot silently break enforced contracts when a safe
/// replacement exists.
///
/// Trace addresses (D3):
/// - `corrupt.pcr.count` — total PCR errors applied
/// - `corrupt.pcr.error_site[i]` — pool position of the i-th error
/// - `corrupt.pcr.error_base[i]` — replacement base at the i-th error
pub struct PCRErrorPass {
    count_source: crate::passes::count_source::CountSource,
    base_dist: Box<dyn Distribution<Output = u8>>,
}

impl PCRErrorPass {
    /// Construct from an explicit count distribution. Per pass
    /// execution, the count is sampled once independently of pool
    /// length — matching v1 semantics.
    pub fn new(
        count_dist: Box<dyn Distribution<Output = i64>>,
        base_dist: Box<dyn Distribution<Output = u8>>,
    ) -> Self {
        Self {
            count_source: crate::passes::count_source::CountSource::Distribution(count_dist),
            base_dist,
        }
    }

    /// Construct from a per-base error rate. Per pass execution, the
    /// count is drawn from `Poisson(rate * pool_len)` against the
    /// current pool length — matching how PCR error is universally
    /// reported in the literature (per-cycle, per-base rate).
    pub fn new_rate(rate: f64, base_dist: Box<dyn Distribution<Output = u8>>) -> Self {
        assert!(
            rate.is_finite() && (0.0..=1.0).contains(&rate),
            "PCRErrorPass: rate must be in [0.0, 1.0], got {}",
            rate
        );
        Self {
            count_source: crate::passes::count_source::CountSource::Rate(rate),
            base_dist,
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
            address::ChoiceAddress::CorruptPcrCount,
            strict,
        )?;

        if pool_len == 0 || count == 0 {
            return Ok(sim.clone());
        }

        let mut tx = MutationTransaction::open(sim, ctx, self.name(), strict);

        for i in 0..count {
            tx.substitute_position_constrained(
                self.base_dist.as_ref(),
                address::ChoiceAddress::CorruptPcrSite(i),
                address::ChoiceAddress::CorruptPcrBase(i),
                None,
            )?;
        }

        tx.commit()
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

    fn declared_choice_patterns(&self) -> Vec<address::ChoiceAddressPattern> {
        vec![
            address::ChoiceAddressPattern::CorruptPcrCount,
            address::ChoiceAddressPattern::CorruptPcrSite,
            address::ChoiceAddressPattern::CorruptPcrBase,
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
    use crate::pass::testing::PassRuntime;
    use crate::pass::PassPlan;
    use crate::passes::test_support::{
        make_fully_locked_vj_fixture, make_substitution_productive_vj_fixture,
        StopOnlyMutationBaseDist, StopThenSafeMutationBaseDist,
    };
    use crate::rng::Rng;
    use crate::trace::Trace;

    fn pcr_test_sim() -> Simulation {
        let mut sim = Simulation::new();
        for (i, b) in b"AAACCCGGGTTTAAACCCGGG".iter().enumerate() {
            let (next, _) =
                sim.with_nucleotide_pushed(Nucleotide::germline(*b, i as u16, Segment::V));
            sim = next;
        }
        let region = Region::new(Segment::V, NucHandle::new(0), NucHandle::new(21));
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

    /// Constrained-path seed finder for the v3.0 site-weighted
    /// PCR path: find a seed under which the contract-aware sample
    /// lands at `target_site` and the unconstrained sample at the
    /// same seed produces a stop violation.
    fn find_seed_for_constrained_pcr_error_site(
        cfg: &crate::refdata::RefDataConfig,
        contracts: &crate::contract::ContractSet,
        sim: &Simulation,
        target_site: u32,
    ) -> u64 {
        for seed in 0..2048u64 {
            let outcome = PassRuntime::execute_with_context(
                &pcr_single_error_plan(Box::new(StopThenSafeMutationBaseDist)),
                sim.clone(),
                seed,
                Some(cfg),
                Some(contracts),
            );
            let site = match outcome.trace.find("corrupt.pcr.error_site[0]") {
                Some(rec) => match rec.value {
                    ChoiceValue::Int(s) => s as u32,
                    _ => continue,
                },
                None => continue,
            };
            if site != target_site {
                continue;
            }
            let unconstrained = PassRuntime::execute_with_context(
                &pcr_single_error_plan(Box::new(StopThenSafeMutationBaseDist)),
                sim.clone(),
                seed,
                Some(cfg),
                None,
            );
            if contracts
                .verify(unconstrained.final_simulation(), Some(cfg))
                .is_err()
            {
                return seed;
            }
        }
        panic!(
            "no seed in search range produced constrained PCR site {} with unconstrained stop violation",
            target_site
        );
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
    fn pcr_error_no_contracts_rng_and_trace_shape_are_characterized() {
        let mut sim = Simulation::new();
        for (i, b) in b"AAAACCCCGGGGTTTT".iter().enumerate() {
            let (next, _) =
                sim.with_nucleotide_pushed(Nucleotide::germline(*b, i as u16, Segment::V));
            sim = next;
        }
        sim = sim.with_region_added(Region::new(
            Segment::V,
            NucHandle::new(0),
            NucHandle::new(16),
        ));

        let pass = PCRErrorPass::new(
            Box::new(EmpiricalLengthDist::from_pairs(vec![(2, 1.0)])),
            Box::new(UniformBase),
        );
        let mut trace = Trace::new();
        let mut rng = Rng::new(0xbeef);
        {
            let mut ctx = PassContext {
                replay_cursor: None,
                trace: &mut trace,
                rng: &mut rng,
                pass_index: 0,
                refdata: None,
                contracts: None,
                feasibility: None,
                reference_index: None,
                event_log_sink: None,
            };
            let _ = pass.execute(&sim, &mut ctx);
        }

        // Fixed count distribution consumes one f64. Each PCR error
        // consumes one word for the site and one word for the base.
        assert_eq!(rng.words_consumed(), 1 + 2 * 2);
        assert_eq!(trace.len(), 1 + 2 * 2);
        assert_eq!(
            trace.find("corrupt.pcr.count").unwrap().value,
            ChoiceValue::Int(2)
        );
        for i in 0..2 {
            assert!(trace
                .find(&format!("corrupt.pcr.error_site[{}]", i))
                .is_some());
            assert!(trace
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
        // Verify that on-demand recompute against the post-error pool
        // is deterministic.
        let mut plan = PassPlan::new();
        plan.push(Box::new(PCRErrorPass::new(
            Box::new(EmpiricalLengthDist::from_pairs(vec![(5, 1.0)])),
            Box::new(UniformBase),
        )));
        let outcome = PassRuntime::execute(&plan, pcr_test_sim(), 1);
        let final_sim = outcome.final_simulation();

        let a = crate::ir::compute_codon_rail(&final_sim.sequence.regions[0], &final_sim.pool);
        let b = crate::ir::compute_codon_rail(&final_sim.sequence.regions[0], &final_sim.pool);
        assert_eq!(a.amino_acids, b.amino_acids);
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
        let contracts = productive();
        // v3.0 constrain-before-propose weights site selection by
        // admissible mass — the seed → site map differs from raw
        // uniform sampling, so the search runs through the
        // contract-aware path.
        let seed = find_seed_for_constrained_pcr_error_site(&cfg, &contracts, &sim, 2);

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
        // Under v3.0 constrain-before-propose, strict-mode
        // `EmptyAdmissibleSupport` fires when *no* (site, base)
        // combination across the pool admits the draw. Use the
        // fully-locked V=TGG / J=TGG fixture so every site rejects
        // `{A}`.
        let (cfg, sim) = make_fully_locked_vj_fixture();
        let contracts = productive();

        let err = PassRuntime::execute_strict_with_context(
            &pcr_single_error_plan(Box::new(StopOnlyMutationBaseDist)),
            sim,
            0,
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
