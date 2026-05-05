//! `ContaminantPass` — wholesale sequence replacement (E.6).

use crate::contract::ChoiceContext;
use crate::dist::{sample_filtered_result, Distribution, FilteredSampleError};
use crate::ir::{NucHandle, Simulation};
use crate::pass::{Pass, PassContext, PassError};
use crate::trace::ChoiceValue;

/// Models read contamination: with probability `apply_prob` the
/// entire assembled pool is overwritten with bases drawn from a
/// contaminant distribution. Used to simulate primer dimers,
/// bacterial DNA, or any non-receptor sequence that ends up in a
/// receptor-sequencing library.
///
/// **Architectural shape vs other corruption passes:**
/// - PCR / quality errors are *count-driven* — sample N positions,
///   substitute each.
/// - Contaminant is *probability-driven at the read level* — one
///   coin flip decides whether the *entire* read is wiped, then if
///   yes, every base gets a contaminant draw.
///
/// This means the trace begins with a single Boolean choice:
/// `corrupt.contaminant.applied`. When that's `Bool(true)`, the
/// trace continues with one base entry per pool position; when it's
/// `Bool(false)`, no further records are emitted (the pool is
/// returned unchanged).
///
/// Codon rail consistency is preserved per-base by
/// `with_base_changed`. After a contamination event, the affected
/// region's `amino_acids` reflects the contaminant content (not the
/// original germline) — which is the desired post-contamination
/// semantics. When contracts are active, each replacement base is
/// filtered against the current intermediate IR and the target site,
/// so a contamination event cannot transiently violate enforced
/// contracts when an admissible replacement exists.
///
/// Trace addresses (D3):
/// - `corrupt.contaminant.applied` — `Bool(true)` if contamination
///   was applied, `Bool(false)` otherwise.
/// - `corrupt.contaminant.bases[i]` — i-th replacement base, only
///   present when `applied = true`.
pub struct ContaminantPass {
    apply_prob: f64,
    base_dist: Box<dyn Distribution<Output = u8>>,
}

impl ContaminantPass {
    /// Construct a contaminant pass.
    ///
    /// Panics if `apply_prob` is not in `[0.0, 1.0]` or is non-finite.
    pub fn new(apply_prob: f64, base_dist: Box<dyn Distribution<Output = u8>>) -> Self {
        assert!(
            apply_prob.is_finite() && (0.0..=1.0).contains(&apply_prob),
            "ContaminantPass: apply_prob must be in [0.0, 1.0], got {}",
            apply_prob
        );
        Self {
            apply_prob,
            base_dist,
        }
    }

    pub fn apply_prob(&self) -> f64 {
        self.apply_prob
    }

    fn constraint_sampling_error(
        &self,
        address: impl Into<String>,
        reason: FilteredSampleError,
    ) -> PassError {
        PassError::constraint_sampling(self.name(), address, reason)
    }

    fn sample_base(
        &self,
        sim: &Simulation,
        ctx: &mut PassContext,
        address: &str,
        index: u32,
        count: u32,
        site: NucHandle,
        strict: bool,
    ) -> Result<u8, PassError> {
        let refdata = ctx.refdata;
        let contracts = ctx.contracts;

        if let Some(contracts) = contracts {
            let context = ChoiceContext::indexed_target(index, count, site);
            match sample_filtered_result(ctx.rng, self.base_dist.as_ref(), |candidate: &u8| {
                contracts
                    .admits_with_context(
                        sim,
                        refdata,
                        address,
                        &ChoiceValue::Base(*candidate),
                        context,
                    )
                    .is_ok()
            }) {
                Ok(value) => return Ok(value),
                Err(reason) if strict => {
                    return Err(self.constraint_sampling_error(address, reason));
                }
                Err(_) => {
                    // Permissive legacy path: if filtering cannot
                    // produce a value, preserve unconstrained
                    // contaminant replacement.
                }
            }
        }

        Ok(self.base_dist.sample(ctx.rng))
    }

    fn execute_with_sampling_mode(
        &self,
        sim: &Simulation,
        ctx: &mut PassContext,
        strict: bool,
    ) -> Result<Simulation, PassError> {
        // 1. Coin flip: is this read contaminated?
        let coin = ctx.rng.next_f64();
        let applied = coin < self.apply_prob;
        ctx.trace
            .record("corrupt.contaminant.applied", ChoiceValue::Bool(applied));

        if !applied {
            return Ok(sim.clone());
        }

        let pool_len = sim.pool.len() as u32;
        if pool_len == 0 {
            return Ok(sim.clone());
        }

        // 2. Replace every base in the pool with a contaminant draw.
        let mut current = sim.clone();
        for i in 0..pool_len {
            let site = NucHandle::new(i);
            let address = format!("corrupt.contaminant.bases[{}]", i);
            let new_base = self.sample_base(&current, ctx, &address, i, pool_len, site, strict)?;
            ctx.trace.record(address, ChoiceValue::Base(new_base));
            current = current.with_base_changed(site, new_base);
        }
        Ok(current)
    }
}

impl Pass for ContaminantPass {
    fn name(&self) -> &str {
        "corrupt.contaminant"
    }

    fn execute(&self, sim: &Simulation, ctx: &mut PassContext) -> Simulation {
        self.execute_with_sampling_mode(sim, ctx, false)
            .expect("ContaminantPass permissive execution must not return PassError")
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
            "corrupt.contaminant.applied".to_string(),
            "corrupt.contaminant.bases[0..n]".to_string(),
        ]
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::assignment::AlleleInstance;
    use crate::contract::productive;
    use crate::dist::UniformBase;
    use crate::ir::{Nucleotide, Region, Segment};
    use crate::pass::{PassPlan, PassRuntime};
    use crate::refdata::{Allele, AlleleId, ChainType, RefDataConfig};

    fn contaminant_test_sim() -> Simulation {
        let mut sim = Simulation::new();
        for (i, b) in b"AAACCCGGG".iter().enumerate() {
            let (next, _) =
                sim.with_nucleotide_pushed(Nucleotide::germline(*b, i as u16, Segment::V));
            sim = next;
        }
        sim
    }

    #[derive(Clone, Debug)]
    struct StopOnlyContaminantBaseDist;

    impl Distribution for StopOnlyContaminantBaseDist {
        type Output = u8;

        fn sample(&self, _rng: &mut crate::rng::Rng) -> u8 {
            b'T'
        }

        fn support(&self) -> Option<Vec<(u8, f64)>> {
            Some(vec![(b'T', 1.0)])
        }
    }

    fn contaminant_plan(base_dist: Box<dyn Distribution<Output = u8>>) -> PassPlan {
        let mut plan = PassPlan::new();
        plan.push(Box::new(ContaminantPass::new(1.0, base_dist)));
        plan
    }

    fn make_contaminant_productive_vj_fixture() -> (RefDataConfig, Simulation) {
        let mut cfg = RefDataConfig::empty(ChainType::Vj);
        let _ = cfg.v_pool.push(Allele {
            name: "v_contam*01".into(),
            gene: "v_contam".into(),
            seq: b"AAA".to_vec(),
            segment: Segment::V,
            anchor: Some(0),
        });
        let _ = cfg.j_pool.push(Allele {
            name: "j_contam*01".into(),
            gene: "j_contam".into(),
            seq: b"TGG".to_vec(),
            segment: Segment::J,
            anchor: Some(0),
        });

        let mut sim = Simulation::new();
        for (i, &b) in b"AAA".iter().enumerate() {
            let (next, _) =
                sim.with_nucleotide_pushed(Nucleotide::germline(b, i as u16, Segment::V));
            sim = next;
        }
        let v_region = Region::new(Segment::V, NucHandle::new(0), NucHandle::new(3))
            .with_codon_rail_recomputed(&sim.pool);
        sim = sim.with_region_added(v_region);

        for (i, &b) in b"TGG".iter().enumerate() {
            let (next, _) =
                sim.with_nucleotide_pushed(Nucleotide::germline(b, i as u16, Segment::J));
            sim = next;
        }
        let j_region = Region::new(Segment::J, NucHandle::new(3), NucHandle::new(6))
            .with_codon_rail_recomputed(&sim.pool);
        sim = sim.with_region_added(j_region);

        sim = sim
            .with_allele_assigned(Segment::V, AlleleInstance::new(AlleleId::new(0)))
            .with_allele_assigned(Segment::J, AlleleInstance::new(AlleleId::new(0)));

        (cfg, sim)
    }

    fn find_seed_for_unconstrained_contaminant_prefix(sim: &Simulation, expected: &[u8]) -> u64 {
        for seed in 0..4096u64 {
            let outcome =
                PassRuntime::execute(&contaminant_plan(Box::new(UniformBase)), sim.clone(), seed);
            let matches = expected.iter().enumerate().all(|(i, &base)| {
                outcome
                    .trace
                    .find(&format!("corrupt.contaminant.bases[{}]", i))
                    .map(|rec| rec.value == ChoiceValue::Base(base))
                    .unwrap_or(false)
            });
            if matches {
                return seed;
            }
        }
        panic!(
            "no seed in search range produced contaminant prefix {:?}",
            expected
        );
    }

    #[test]
    #[should_panic(expected = "apply_prob must be in [0.0, 1.0]")]
    fn contaminant_pass_rejects_negative_probability() {
        let _ = ContaminantPass::new(-0.1, Box::new(UniformBase));
    }

    #[test]
    #[should_panic(expected = "apply_prob must be in [0.0, 1.0]")]
    fn contaminant_pass_rejects_probability_above_one() {
        let _ = ContaminantPass::new(1.5, Box::new(UniformBase));
    }

    #[test]
    #[should_panic(expected = "apply_prob must be in [0.0, 1.0]")]
    fn contaminant_pass_rejects_nan_probability() {
        let _ = ContaminantPass::new(f64::NAN, Box::new(UniformBase));
    }

    #[test]
    fn contaminant_pass_zero_probability_never_applies() {
        let mut plan = PassPlan::new();
        plan.push(Box::new(ContaminantPass::new(0.0, Box::new(UniformBase))));

        for seed in 0..20u64 {
            let outcome = PassRuntime::execute(&plan, contaminant_test_sim(), seed);
            // Trace records `applied: Bool(false)` and nothing else.
            assert_eq!(outcome.trace.len(), 1);
            assert_eq!(
                outcome
                    .trace
                    .find("corrupt.contaminant.applied")
                    .unwrap()
                    .value,
                ChoiceValue::Bool(false)
            );
            // Pool unchanged.
            for i in 0..9 {
                let b = outcome
                    .final_simulation()
                    .pool
                    .get(NucHandle::new(i as u32))
                    .unwrap()
                    .base;
                assert!(matches!(b, b'A' | b'C' | b'G'));
            }
        }
    }

    #[test]
    fn contaminant_pass_one_probability_always_applies() {
        let mut plan = PassPlan::new();
        plan.push(Box::new(ContaminantPass::new(1.0, Box::new(UniformBase))));

        for seed in 0..20u64 {
            let outcome = PassRuntime::execute(&plan, contaminant_test_sim(), seed);
            // Trace records `applied: Bool(true)` followed by 9 base records.
            assert_eq!(outcome.trace.len(), 10);
            assert_eq!(
                outcome
                    .trace
                    .find("corrupt.contaminant.applied")
                    .unwrap()
                    .value,
                ChoiceValue::Bool(true)
            );
            for i in 0..9 {
                assert!(outcome
                    .trace
                    .find(&format!("corrupt.contaminant.bases[{}]", i))
                    .is_some());
            }
        }
    }

    #[test]
    fn contaminant_pass_one_probability_pool_reflects_recorded_bases() {
        let mut plan = PassPlan::new();
        plan.push(Box::new(ContaminantPass::new(1.0, Box::new(UniformBase))));
        let outcome = PassRuntime::execute(&plan, contaminant_test_sim(), 42);
        let final_sim = outcome.final_simulation();

        for i in 0..9u32 {
            let recorded = match outcome
                .trace
                .find(&format!("corrupt.contaminant.bases[{}]", i))
                .unwrap()
                .value
            {
                ChoiceValue::Base(b) => b,
                _ => unreachable!(),
            };
            let pool_base = final_sim.pool.get(NucHandle::new(i)).unwrap().base;
            assert_eq!(recorded, pool_base, "trace lies at position {}", i);
        }
    }

    #[test]
    fn contaminant_pass_half_probability_mixed_outcomes() {
        // Across many seeds, p=0.5 should produce both applied and
        // not-applied outcomes. Negative-control proof that the
        // coin flip is actually firing.
        let mut plan = PassPlan::new();
        plan.push(Box::new(ContaminantPass::new(0.5, Box::new(UniformBase))));

        let mut applied_count = 0;
        let mut not_applied_count = 0;
        for seed in 0..50u64 {
            let outcome = PassRuntime::execute(&plan, contaminant_test_sim(), seed);
            match outcome
                .trace
                .find("corrupt.contaminant.applied")
                .unwrap()
                .value
            {
                ChoiceValue::Bool(true) => applied_count += 1,
                ChoiceValue::Bool(false) => not_applied_count += 1,
                _ => unreachable!(),
            }
        }
        assert!(applied_count > 0, "expected at least one applied outcome");
        assert!(
            not_applied_count > 0,
            "expected at least one not-applied outcome"
        );
    }

    #[test]
    fn contaminant_pass_codon_rail_refresh_after_application() {
        // When contamination applies, the assembled region's codon
        // rail must reflect the new (contaminant) bases.
        let mut sim = Simulation::new();
        for (i, b) in b"ATGGGGGGG".iter().enumerate() {
            let (next, _) =
                sim.with_nucleotide_pushed(Nucleotide::germline(*b, i as u16, Segment::V));
            sim = next;
        }
        let region = Region::new(Segment::V, NucHandle::new(0), NucHandle::new(9))
            .with_codon_rail_recomputed(&sim.pool);
        sim = sim.with_region_added(region);
        // Before contamination: M G G.
        assert_eq!(sim.sequence.regions[0].amino_acids, b"MGG");

        let mut plan = PassPlan::new();
        plan.push(Box::new(ContaminantPass::new(1.0, Box::new(UniformBase))));
        let outcome = PassRuntime::execute(&plan, sim, 0);
        let final_sim = outcome.final_simulation();

        // After contamination: codon rail recomputed, may differ.
        let stored = &final_sim.sequence.regions[0].amino_acids;
        let fresh = final_sim.sequence.regions[0].with_codon_rail_recomputed(&final_sim.pool);
        assert_eq!(stored, &fresh.amino_acids);
    }

    #[test]
    fn contaminant_pass_is_deterministic() {
        let plan = || {
            let mut p = PassPlan::new();
            p.push(Box::new(ContaminantPass::new(0.5, Box::new(UniformBase))));
            p
        };
        let oa = PassRuntime::execute(&plan(), contaminant_test_sim(), 0xc0ff_ee);
        let ob = PassRuntime::execute(&plan(), contaminant_test_sim(), 0xc0ff_ee);
        assert_eq!(oa.trace.choices(), ob.trace.choices());
        for i in 0..9u32 {
            assert_eq!(
                oa.final_simulation()
                    .pool
                    .get(NucHandle::new(i))
                    .unwrap()
                    .base,
                ob.final_simulation()
                    .pool
                    .get(NucHandle::new(i))
                    .unwrap()
                    .base
            );
        }
    }

    #[test]
    fn contaminant_pass_empty_pool_skips_replacement() {
        let mut plan = PassPlan::new();
        plan.push(Box::new(ContaminantPass::new(1.0, Box::new(UniformBase))));
        let outcome = PassRuntime::execute(&plan, Simulation::new(), 0);
        // Coin flip happened (Bool recorded) but no bases written.
        assert_eq!(outcome.trace.len(), 1);
        assert_eq!(
            outcome
                .trace
                .find("corrupt.contaminant.applied")
                .unwrap()
                .value,
            ChoiceValue::Bool(true)
        );
    }

    #[test]
    fn contaminant_pass_declared_choices() {
        let pass = ContaminantPass::new(0.1, Box::new(UniformBase));
        let declared = pass.declared_choices();
        assert!(declared.contains(&"corrupt.contaminant.applied".to_string()));
        assert!(declared.contains(&"corrupt.contaminant.bases[0..n]".to_string()));
    }

    #[test]
    fn contaminant_productive_filters_base_that_would_create_stop() {
        let (cfg, sim) = make_contaminant_productive_vj_fixture();
        let seed = find_seed_for_unconstrained_contaminant_prefix(&sim, b"TAA");
        let contracts = productive();

        let constrained = PassRuntime::execute_with_context(
            &contaminant_plan(Box::new(UniformBase)),
            sim.clone(),
            seed,
            Some(&cfg),
            Some(&contracts),
        );

        assert_eq!(
            constrained
                .trace
                .find("corrupt.contaminant.applied")
                .unwrap()
                .value,
            ChoiceValue::Bool(true)
        );
        assert_ne!(
            constrained
                .trace
                .find("corrupt.contaminant.bases[0]")
                .unwrap()
                .value,
            ChoiceValue::Base(b'T')
        );
        assert!(contracts
            .verify(constrained.final_simulation(), Some(&cfg))
            .is_ok());

        let unconstrained = PassRuntime::execute_with_context(
            &contaminant_plan(Box::new(UniformBase)),
            sim,
            seed,
            Some(&cfg),
            None,
        );
        let first_three: Vec<u8> = (0..3)
            .map(|i| {
                match unconstrained
                    .trace
                    .find(&format!("corrupt.contaminant.bases[{}]", i))
                    .unwrap()
                    .value
                {
                    ChoiceValue::Base(b) => b,
                    _ => unreachable!(),
                }
            })
            .collect();
        assert_eq!(first_three, b"TAA");
        let violations = productive()
            .verify(unconstrained.final_simulation(), Some(&cfg))
            .unwrap_err();
        assert!(violations
            .iter()
            .any(|v| v.contract_name == "no_stop_codon_in_junction"));
    }

    #[test]
    fn contaminant_strict_errors_when_base_filter_empty() {
        let (cfg, sim) = make_contaminant_productive_vj_fixture();
        let contracts = productive();

        let err = PassRuntime::execute_strict_with_context(
            &contaminant_plan(Box::new(StopOnlyContaminantBaseDist)),
            sim,
            0,
            Some(&cfg),
            Some(&contracts),
        )
        .unwrap_err();

        assert_eq!(err.pass_name(), "corrupt.contaminant");
        assert_eq!(err.address(), "corrupt.contaminant.bases[0]");
        assert_eq!(
            err.constraint_reason(),
            Some(FilteredSampleError::EmptyAdmissibleSupport)
        );
    }
}
