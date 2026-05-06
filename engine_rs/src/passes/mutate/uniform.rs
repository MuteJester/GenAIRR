//! `UniformMutationPass` — simplest SHM model (Phase E.1).

use crate::dist::Distribution;
use crate::ir::{NucHandle, Simulation};
use crate::pass::{Pass, PassContext, PassEffect, PassError};
use crate::passes::constrained::{sample_targeted_base, TargetedBaseChoice};
use crate::trace::ChoiceValue;

/// The simplest mutation pass: pick `N` positions uniformly across
/// the assembled pool and substitute each with a base drawn from
/// `base_dist`.
///
/// Models position-independent point mutations. Real biology uses
/// the context-dependent S5F model (Phase E.3); this pass is the
/// architectural reference for any SHM-like pass — establishes the
/// trace-address shape, the per-mutation IR-revision flow, and the
/// integration with constraint-aware sampling.
///
/// **Determinism:** consumes RNG words deterministically — one for
/// the count, then two per mutation (site index + base). Same seed
/// → same mutations.
///
/// **Codon-rail consistency:** every `with_base_changed` call
/// auto-refreshes the affected region's codon rail (per the
/// post-D audit fix). Stop codons that get introduced will be
/// visible in `Region.amino_acids` immediately.
///
/// **Constraint awareness:** every candidate base is filtered
/// against the active `ContractSet` with the chosen target site in
/// [`ChoiceContext`]. This lets contracts such as
/// `NoStopCodonInJunction` reject substitutions that would leave the
/// junction non-productive before the mutation is committed. The
/// mutation count and site choice remain unconstrained for now; the
/// base draw is the first biologically meaningful hardening point.
///
/// Trace addresses (D3 hierarchical strings):
/// - `mutate.uniform.count` — total mutations applied
/// - `mutate.uniform.site[i]` — position of the i-th mutation
/// - `mutate.uniform.base[i]` — new base at the i-th mutation
pub struct UniformMutationPass {
    count_dist: Box<dyn Distribution<Output = i64>>,
    base_dist: Box<dyn Distribution<Output = u8>>,
}

impl UniformMutationPass {
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
        // 1. Sample the number of mutations to apply.
        let count_raw = self.count_dist.sample(ctx.rng);
        if strict && count_raw < 0 {
            return Err(PassError::invalid_distribution_output(
                self.name(),
                "mutate.uniform.count",
                count_raw,
                "negative_count",
            ));
        }
        if strict && count_raw > u32::MAX as i64 {
            return Err(PassError::invalid_distribution_output(
                self.name(),
                "mutate.uniform.count",
                count_raw,
                "count_exceeds_u32",
            ));
        }
        assert!(
            count_raw >= 0,
            "UniformMutationPass: count distribution returned negative {}",
            count_raw
        );
        assert!(
            count_raw <= u32::MAX as i64,
            "UniformMutationPass: count distribution returned {} > u32::MAX",
            count_raw
        );
        let count = count_raw as u32;
        ctx.trace
            .record("mutate.uniform.count", ChoiceValue::Int(count_raw));

        // No-op if the pool is empty — nothing to mutate.
        let pool_len = sim.pool.len() as u32;
        if pool_len == 0 || count == 0 {
            return Ok(sim.clone());
        }

        // 2. Apply `count` mutations sequentially. Each mutation
        //    samples a site (uniform in [0, pool_len)) and then draws
        //    an admissible base for that target when contracts are
        //    active. The IR evolves one mutation at a time through
        //    the persistent API.
        let mut current = sim.clone();
        for i in 0..count {
            let site = ctx.rng.range_u32(pool_len);
            let site_handle = NucHandle::new(site);
            let base_address = format!("mutate.uniform.base[{}]", i);

            ctx.trace.record(
                format!("mutate.uniform.site[{}]", i),
                ChoiceValue::Int(site as i64),
            );
            let new_base = sample_targeted_base(
                &current,
                ctx,
                self.base_dist.as_ref(),
                TargetedBaseChoice::new(self.name(), &base_address, i, count, site_handle, strict),
            )?;
            ctx.trace.record(base_address, ChoiceValue::Base(new_base));

            // `with_base_changed` auto-refreshes the codon rail of
            // any region containing `site` (post-D audit fix).
            current = current.with_base_changed(site_handle, new_base);
        }

        Ok(current)
    }
}

impl Pass for UniformMutationPass {
    fn name(&self) -> &str {
        "mutate.uniform"
    }

    fn execute(&self, sim: &Simulation, ctx: &mut PassContext) -> Simulation {
        self.execute_with_sampling_mode(sim, ctx, false)
            .expect("UniformMutationPass permissive execution must not return PassError")
    }

    fn execute_checked(
        &self,
        sim: &Simulation,
        ctx: &mut PassContext,
    ) -> Result<Simulation, PassError> {
        self.execute_with_sampling_mode(sim, ctx, true)
    }

    fn declared_choices(&self) -> Vec<String> {
        // Count is fixed-address; site and base are variable per
        // mutation. Use the [0..n] expansion convention from D3.
        vec![
            "mutate.uniform.count".to_string(),
            "mutate.uniform.site[0..n]".to_string(),
            "mutate.uniform.base[0..n]".to_string(),
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

    fn uniform_mutation_plan(base_dist: Box<dyn Distribution<Output = u8>>) -> PassPlan {
        let mut plan = PassPlan::new();
        plan.push(Box::new(UniformMutationPass::new(
            Box::new(EmpiricalLengthDist::from_pairs(vec![(1, 1.0)])),
            base_dist,
        )));
        plan
    }

    fn uniform_mutation_site(seed: u64, sim: Simulation) -> u32 {
        let outcome = PassRuntime::execute(
            &uniform_mutation_plan(Box::new(StopOnlyMutationBaseDist)),
            sim,
            seed,
        );
        match outcome.trace.find("mutate.uniform.site[0]").unwrap().value {
            ChoiceValue::Int(site) => site as u32,
            _ => panic!("wrong variant"),
        }
    }

    fn find_seed_for_uniform_mutation_site(sim: &Simulation, target_site: u32) -> u64 {
        for seed in 0..512u64 {
            if uniform_mutation_site(seed, sim.clone()) == target_site {
                return seed;
            }
        }
        panic!("no seed in search range targeted site {}", target_site);
    }

    #[test]
    fn uniform_mutation_pass_zero_count_is_noop() {
        // Build a sim with a few nucleotides + region, run uniform
        // mutation with count=0, verify nothing changed.
        let mut sim = Simulation::new();
        for (i, b) in b"AAACCCGGG".iter().enumerate() {
            let (next, _) =
                sim.with_nucleotide_pushed(Nucleotide::germline(*b, i as u16, Segment::V));
            sim = next;
        }
        let region = Region::new(Segment::V, NucHandle::new(0), NucHandle::new(9))
            .with_codon_rail_recomputed(&sim.pool);
        sim = sim.with_region_added(region);

        let mut plan = PassPlan::new();
        plan.push(Box::new(UniformMutationPass::new(
            Box::new(EmpiricalLengthDist::from_pairs(vec![(0, 1.0)])),
            Box::new(UniformBase),
        )));

        let outcome = PassRuntime::execute(&plan, sim.clone(), 42);
        let final_sim = outcome.final_simulation();

        // Pool unchanged.
        for i in 0..9 {
            assert_eq!(
                final_sim.pool.get(NucHandle::new(i)).unwrap().base,
                sim.pool.get(NucHandle::new(i)).unwrap().base
            );
        }
        // Codon rail unchanged.
        assert_eq!(
            final_sim.sequence.regions[0].amino_acids,
            sim.sequence.regions[0].amino_acids
        );
        // Trace recorded count=0, no site/base entries.
        assert_eq!(outcome.trace.len(), 1);
        assert_eq!(
            outcome.trace.find("mutate.uniform.count").unwrap().value,
            ChoiceValue::Int(0)
        );
    }

    #[test]
    fn uniform_mutation_pass_applies_n_mutations_with_traced_addresses() {
        let mut sim = Simulation::new();
        for (i, b) in b"AAACCCGGGTTT".iter().enumerate() {
            let (next, _) =
                sim.with_nucleotide_pushed(Nucleotide::germline(*b, i as u16, Segment::V));
            sim = next;
        }
        let region = Region::new(Segment::V, NucHandle::new(0), NucHandle::new(12))
            .with_codon_rail_recomputed(&sim.pool);
        sim = sim.with_region_added(region);

        let mut plan = PassPlan::new();
        plan.push(Box::new(UniformMutationPass::new(
            Box::new(EmpiricalLengthDist::from_pairs(vec![(5, 1.0)])),
            Box::new(UniformBase),
        )));

        let outcome = PassRuntime::execute(&plan, sim, 7);

        // Trace: 1 count + 5 site + 5 base = 11 records.
        assert_eq!(outcome.trace.len(), 11);
        assert_eq!(
            outcome.trace.find("mutate.uniform.count").unwrap().value,
            ChoiceValue::Int(5)
        );
        for i in 0..5 {
            let site_addr = format!("mutate.uniform.site[{}]", i);
            let base_addr = format!("mutate.uniform.base[{}]", i);
            assert!(outcome.trace.find(&site_addr).is_some());
            assert!(outcome.trace.find(&base_addr).is_some());

            // Site is in [0, 12).
            match outcome.trace.find(&site_addr).unwrap().value {
                ChoiceValue::Int(s) => assert!(s >= 0 && s < 12),
                _ => panic!("wrong variant"),
            }
            // Base is one of A/C/G/T.
            match outcome.trace.find(&base_addr).unwrap().value {
                ChoiceValue::Base(b) => assert!(matches!(b, b'A' | b'C' | b'G' | b'T')),
                _ => panic!("wrong variant"),
            }
        }
    }

    #[test]
    fn uniform_mutation_pass_pool_reflects_recorded_mutations() {
        // Faithfulness: the pool's base at each recorded site equals
        // the recorded new base (trace is honest).
        let mut sim = Simulation::new();
        for (i, b) in b"AAAAAAAAAA".iter().enumerate() {
            let (next, _) =
                sim.with_nucleotide_pushed(Nucleotide::germline(*b, i as u16, Segment::V));
            sim = next;
        }

        let mut plan = PassPlan::new();
        plan.push(Box::new(UniformMutationPass::new(
            Box::new(EmpiricalLengthDist::from_pairs(vec![(7, 1.0)])),
            Box::new(UniformBase),
        )));
        let outcome = PassRuntime::execute(&plan, sim, 0xc0ff_ee);
        let final_sim = outcome.final_simulation();

        // Walk the trace: at each (site[i], base[i]) the pool should
        // hold base[i]. (NOTE: later mutations could overwrite earlier
        // ones at the same site, so we can only check the LAST
        // mutation per site.)
        let mut last_at_site = std::collections::HashMap::new();
        for i in 0..7 {
            let s = match outcome
                .trace
                .find(&format!("mutate.uniform.site[{}]", i))
                .unwrap()
                .value
            {
                ChoiceValue::Int(n) => n as u32,
                _ => unreachable!(),
            };
            let b = match outcome
                .trace
                .find(&format!("mutate.uniform.base[{}]", i))
                .unwrap()
                .value
            {
                ChoiceValue::Base(b) => b,
                _ => unreachable!(),
            };
            last_at_site.insert(s, b);
        }
        for (&site, &expected_base) in last_at_site.iter() {
            let actual = final_sim.pool.get(NucHandle::new(site)).unwrap().base;
            assert_eq!(
                actual, expected_base,
                "trace says site {} got base {}, but pool has {}",
                site, expected_base as char, actual as char
            );
        }
    }

    #[test]
    fn uniform_mutation_pass_refreshes_codon_rail_through_persistent_api() {
        // The post-D audit fix in action: every `with_base_changed`
        // call inside the pass updates Region.amino_acids
        // automatically. After a mutation pass runs, the codon rail
        // should still match a fresh recomputation against the pool.
        let mut sim = Simulation::new();
        for (i, b) in b"ATGGGGAAATTT".iter().enumerate() {
            let (next, _) =
                sim.with_nucleotide_pushed(Nucleotide::germline(*b, i as u16, Segment::V));
            sim = next;
        }
        let region = Region::new(Segment::V, NucHandle::new(0), NucHandle::new(12))
            .with_codon_rail_recomputed(&sim.pool);
        sim = sim.with_region_added(region);

        let mut plan = PassPlan::new();
        plan.push(Box::new(UniformMutationPass::new(
            Box::new(EmpiricalLengthDist::from_pairs(vec![(8, 1.0)])),
            Box::new(UniformBase),
        )));

        let outcome = PassRuntime::execute(&plan, sim, 99);
        let final_sim = outcome.final_simulation();

        // Stored codon rail equals a fresh recomputation against
        // the post-mutation pool. This is the staleness invariant.
        let stored_aa = &final_sim.sequence.regions[0].amino_acids;
        let fresh_region =
            final_sim.sequence.regions[0].with_codon_rail_recomputed(&final_sim.pool);
        assert_eq!(stored_aa, &fresh_region.amino_acids);
        assert_eq!(
            final_sim.sequence.regions[0].stop_codon_positions,
            fresh_region.stop_codon_positions
        );
    }

    #[test]
    fn uniform_mutation_pass_is_deterministic_under_same_seed() {
        let plan = || {
            let mut p = PassPlan::new();
            p.push(Box::new(UniformMutationPass::new(
                Box::new(EmpiricalLengthDist::from_pairs(vec![
                    (3, 1.0),
                    (5, 2.0),
                    (7, 1.0),
                ])),
                Box::new(UniformBase),
            )));
            p
        };
        let build_sim = || {
            let mut s = Simulation::new();
            for (i, b) in b"AAACCCGGGTTTAAA".iter().enumerate() {
                let (next, _) =
                    s.with_nucleotide_pushed(Nucleotide::germline(*b, i as u16, Segment::V));
                s = next;
            }
            s
        };

        let oa = PassRuntime::execute(&plan(), build_sim(), 0xfeed);
        let ob = PassRuntime::execute(&plan(), build_sim(), 0xfeed);

        assert_eq!(oa.trace.choices(), ob.trace.choices());
        // Pool byte-identical.
        for i in 0..15 {
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
    fn uniform_mutation_pass_preserves_persistent_ir() {
        // Original sim must be untouched after the mutation pass
        // runs against a clone of it.
        let mut sim = Simulation::new();
        for (i, b) in b"AAACCC".iter().enumerate() {
            let (next, _) =
                sim.with_nucleotide_pushed(Nucleotide::germline(*b, i as u16, Segment::V));
            sim = next;
        }
        let pre_pool: Vec<u8> = (0..6)
            .map(|i| sim.pool.get(NucHandle::new(i)).unwrap().base)
            .collect();

        let mut plan = PassPlan::new();
        plan.push(Box::new(UniformMutationPass::new(
            Box::new(EmpiricalLengthDist::from_pairs(vec![(20, 1.0)])),
            Box::new(UniformBase),
        )));
        let _ = PassRuntime::execute(&plan, sim.clone(), 5);

        // Verify original sim unchanged.
        let post_pool: Vec<u8> = (0..6)
            .map(|i| sim.pool.get(NucHandle::new(i)).unwrap().base)
            .collect();
        assert_eq!(pre_pool, post_pool);
    }

    #[test]
    fn uniform_mutation_pass_declared_choices_uses_indexed_pattern() {
        let pass = UniformMutationPass::new(
            Box::new(EmpiricalLengthDist::from_pairs(vec![(3, 1.0)])),
            Box::new(UniformBase),
        );
        let declared = pass.declared_choices();
        assert_eq!(declared.len(), 3);
        assert!(declared.contains(&"mutate.uniform.count".to_string()));
        assert!(declared.contains(&"mutate.uniform.site[0..n]".to_string()));
        assert!(declared.contains(&"mutate.uniform.base[0..n]".to_string()));
    }

    #[test]
    fn uniform_mutation_productive_filters_base_that_would_create_stop() {
        let (cfg, sim) = make_substitution_productive_vj_fixture();
        let seed = find_seed_for_uniform_mutation_site(&sim, 2);
        let contracts = productive();

        let constrained = PassRuntime::execute_with_context(
            &uniform_mutation_plan(Box::new(StopThenSafeMutationBaseDist)),
            sim.clone(),
            seed,
            Some(&cfg),
            Some(&contracts),
        );

        assert_eq!(
            constrained
                .trace
                .find("mutate.uniform.site[0]")
                .unwrap()
                .value,
            ChoiceValue::Int(2)
        );
        assert_eq!(
            constrained
                .trace
                .find("mutate.uniform.base[0]")
                .unwrap()
                .value,
            ChoiceValue::Base(b'C')
        );
        assert!(contracts
            .verify(constrained.final_simulation(), Some(&cfg))
            .is_ok());

        let unconstrained = PassRuntime::execute_with_context(
            &uniform_mutation_plan(Box::new(StopThenSafeMutationBaseDist)),
            sim,
            seed,
            Some(&cfg),
            None,
        );
        assert_eq!(
            unconstrained
                .trace
                .find("mutate.uniform.base[0]")
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
    fn uniform_mutation_strict_errors_when_base_filter_empty() {
        let (cfg, sim) = make_substitution_productive_vj_fixture();
        let seed = find_seed_for_uniform_mutation_site(&sim, 2);
        let contracts = productive();

        let err = PassRuntime::execute_strict_with_context(
            &uniform_mutation_plan(Box::new(StopOnlyMutationBaseDist)),
            sim,
            seed,
            Some(&cfg),
            Some(&contracts),
        )
        .unwrap_err();

        assert_eq!(err.pass_name(), "mutate.uniform");
        assert_eq!(err.address(), "mutate.uniform.base[0]");
        assert_eq!(
            err.constraint_reason(),
            Some(FilteredSampleError::EmptyAdmissibleSupport)
        );
    }

    #[test]
    fn uniform_mutation_pass_works_on_empty_pool() {
        // Edge case: no nucleotides yet → mutation is a no-op
        // regardless of the count distribution.
        let mut plan = PassPlan::new();
        plan.push(Box::new(UniformMutationPass::new(
            Box::new(EmpiricalLengthDist::from_pairs(vec![(5, 1.0)])),
            Box::new(UniformBase),
        )));

        let outcome = PassRuntime::execute(&plan, Simulation::new(), 0);
        assert_eq!(outcome.final_simulation().pool.len(), 0);
        // Trace still records the count (sample happened) but no
        // site/base entries because the loop didn't iterate.
        assert_eq!(outcome.trace.len(), 1);
        assert_eq!(
            outcome.trace.find("mutate.uniform.count").unwrap().value,
            ChoiceValue::Int(5)
        );
    }
}
