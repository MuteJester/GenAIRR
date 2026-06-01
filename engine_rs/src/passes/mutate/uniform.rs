//! `UniformMutationPass` — simplest SHM model.

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

/// The simplest mutation pass: pick `N` positions uniformly across
/// the assembled pool and substitute each with a base drawn from
/// `base_dist`.
///
/// Models position-independent point mutations. Real biology uses
/// the context-dependent S5F model; this pass is the
/// architectural reference for any SHM-like pass — establishes the
/// trace-address shape, the per-mutation IR-revision flow, and the
/// integration with constraint-aware sampling.
///
/// **Determinism:** consumes RNG words deterministically — one for
/// the count, then two per mutation (site index + base). Same seed
/// → same mutations.
///
/// **Codon-rail consistency:** codon-rail data is not stored on
/// `Region` — the pool is the authoritative source.
/// `NoStopCodonInJunction` reads pool bytes directly; callers that
/// need amino-acid translation call
/// [`crate::ir::compute_codon_rail`] against the post-mutation pool.
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
    count_source: super::CountSource,
    base_dist: Box<dyn Distribution<Output = u8>>,
    /// Per-biological-segment SHM rate scalars (V / D / J / NP).
    /// Default is flat (all 1.0) — produces byte-identical output
    /// to the pre-slice engine. Non-default values weight site
    /// selection by segment; zero-rate sites drop from support
    /// before contract admissibility.
    segment_rates: super::SegmentRateWeights,
    /// Per-V-subregion SHM rate scalars (Slice B). Default is
    /// flat (all 1.0) — byte-identical to the pre-slice engine.
    /// Composes multiplicatively with `segment_rates` for V
    /// sites; non-V sites are unaffected. V sites on alleles
    /// without subregion annotations receive factor `1.0`.
    v_subregion_rates: super::VSubregionRateWeights,
}

impl UniformMutationPass {
    /// Construct from an explicit count distribution (e.g. fixed `8`
    /// or empirical). The count is sampled once per pass execution
    /// independently of the pool length.
    pub fn new(
        count_dist: Box<dyn Distribution<Output = i64>>,
        base_dist: Box<dyn Distribution<Output = u8>>,
    ) -> Self {
        Self {
            count_source: super::CountSource::Distribution(count_dist),
            base_dist,
            segment_rates: super::SegmentRateWeights::default(),
            v_subregion_rates: super::VSubregionRateWeights::default(),
        }
    }

    /// Construct from a per-base mutation rate (e.g. `0.03` for 3 %
    /// SHM). The count is drawn from `Poisson(rate * pool_len)` per
    /// pass execution — matching how immunologists report SHM in
    /// the literature. See [`super::CountSource`].
    pub fn new_rate(rate: f64, base_dist: Box<dyn Distribution<Output = u8>>) -> Self {
        assert!(
            rate.is_finite() && (0.0..=1.0).contains(&rate),
            "UniformMutationPass: rate must be in [0.0, 1.0], got {}",
            rate
        );
        Self {
            count_source: super::CountSource::Rate(rate),
            base_dist,
            segment_rates: super::SegmentRateWeights::default(),
            v_subregion_rates: super::VSubregionRateWeights::default(),
        }
    }

    /// Persistent setter for per-segment SHM rate scalars. The DSL
    /// boundary in Python validates the four-bucket dict; the
    /// engine binding constructs the [`super::SegmentRateWeights`]
    /// value and installs it via this method. Default (flat 1.0
    /// across V/D/J/NP) is byte-identical to the pre-slice engine.
    #[must_use]
    pub fn with_segment_rates(mut self, rates: super::SegmentRateWeights) -> Self {
        self.segment_rates = rates;
        self
    }

    /// Persistent setter for per-V-subregion SHM rate scalars
    /// (Slice B — `docs/v_subregion_shm_rate_design.md`). The
    /// Python DSL boundary
    /// (`experiment._validate_v_subregion_rates`) validates the
    /// five-label dict with `FWR` / `CDR` alias expansion and
    /// override semantics; the engine binding constructs the
    /// [`super::VSubregionRateWeights`] tuple and installs it via
    /// this method. Default (flat 1.0 across FWR1 / CDR1 / FWR2 /
    /// CDR2 / FWR3) is byte-identical to the pre-slice engine.
    #[must_use]
    pub fn with_v_subregion_rates(mut self, rates: super::VSubregionRateWeights) -> Self {
        self.v_subregion_rates = rates;
        self
    }

    fn execute_with_sampling_mode(
        &self,
        sim: &Simulation,
        ctx: &mut PassContext,
        strict: bool,
    ) -> Result<Simulation, PassError> {
        // 1. Sample the number of mutations to apply. Rate-mode
        //    consults the current pool length; distribution-mode
        //    ignores it. Validation + trace recording shared with
        //    PCR / quality / N-corrupt via `sample_validated_count`.
        let pool_len = sim.pool.len() as u32;
        let count = sample_validated_count(
            &self.count_source,
            ctx,
            pool_len,
            self.name(),
            address::ChoiceAddress::MutateUniformCount,
            strict,
        )?;

        // No-op if the pool is empty — nothing to mutate.
        if pool_len == 0 || count == 0 {
            return Ok(sim.clone());
        }

        // 2. Apply `count` mutations sequentially through a scoped
        //    MutationTransaction. Under v3.0 constrain-before-propose,
        //    `substitute_position_constrained` weights the per-site
        //    selection by admissible mass so sites with narrower
        //    masks (e.g. only 1 base admissible vs the full 4) draw
        //    proportionally less of the per-step probability. When
        //    no contracts are active the helper takes the original
        //    uniform-site + base-sample fast path, preserving the
        //    legacy RNG consumption shape.
        let mut tx = MutationTransaction::open(sim, ctx, self.name(), strict);
        let mut applied: u32 = 0;

        // Pass the configured segment-rate + V-subregion-rate
        // vectors into each per-site sampling call. Flat-default
        // vectors are detected inside the transaction's fast path
        // and skip the per-position lookup. Slice B (the
        // v_subregion_rates kwarg) composes multiplicatively with
        // segment_rates on V sites.
        let seg_rates_ref = Some(&self.segment_rates);
        let v_sub_rates_ref = Some(&self.v_subregion_rates);
        for i in 0..count {
            let wrote = tx.substitute_position_constrained(
                self.base_dist.as_ref(),
                address::ChoiceAddress::MutateUniformSite(i),
                address::ChoiceAddress::MutateUniformBase(i),
                None,
                seg_rates_ref,
                v_sub_rates_ref,
            )?;
            if wrote {
                applied += 1;
            }
        }

        // Realized-count semantics: bump n_mutations by the number
        // of mutations actually written, not the count drawn. When
        // contracts narrow the admissible support to empty across
        // every site, permissive mode skips the slot silently and
        // strict mode surfaces a ConstraintSampling error above.
        if applied > 0 {
            tx.add_to_mutation_count(applied);
        }
        tx.commit()
    }
}

impl Pass for UniformMutationPass {
    fn name(&self) -> &str {
        address::MUTATE_UNIFORM
    }

    fn parameter_signature(&self) -> String {
        use crate::passes::paramsig::{
            fmt_byte_dist, fmt_count_source, fmt_segment_rates,
            fmt_v_subregion_rates, join_parts,
        };
        join_parts([
            fmt_count_source(&self.count_source),
            format!("base={}", fmt_byte_dist(self.base_dist.as_ref())),
            fmt_segment_rates(&self.segment_rates),
            fmt_v_subregion_rates(&self.v_subregion_rates),
        ])
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

    fn declared_choice_patterns(&self) -> Vec<address::ChoiceAddressPattern> {
        vec![
            address::ChoiceAddressPattern::MutateUniformCount,
            address::ChoiceAddressPattern::MutateUniformSite,
            address::ChoiceAddressPattern::MutateUniformBase,
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
    use crate::ir::{compute_codon_rail, Nucleotide, Region, Segment};
    use crate::pass::testing::PassRuntime;
    use crate::pass::PassPlan;
    use crate::passes::test_support::{
        make_fully_locked_vj_fixture, make_substitution_productive_vj_fixture,
        StopOnlyMutationBaseDist, StopThenSafeMutationBaseDist,
    };
    use crate::rng::Rng;
    use crate::trace::Trace;

    fn uniform_mutation_plan(base_dist: Box<dyn Distribution<Output = u8>>) -> PassPlan {
        let mut plan = PassPlan::new();
        plan.push(Box::new(UniformMutationPass::new(
            Box::new(EmpiricalLengthDist::from_pairs(vec![(1, 1.0)])),
            base_dist,
        )));
        plan
    }

    /// Constrained-path seed finder: search for a seed under which the
    /// v3.0 contract-aware uniform pass lands at `target_site`.
    /// Unlike the unconstrained finder, the constrained path weights
    /// site selection by admissible mass, so its seed→site mapping
    /// differs from raw uniform sampling.
    fn find_seed_for_constrained_uniform_mutation_site(
        cfg: &crate::refdata::RefDataConfig,
        contracts: &crate::contract::ContractSet,
        sim: &Simulation,
        target_site: u32,
    ) -> u64 {
        for seed in 0..2048u64 {
            let outcome = PassRuntime::execute_with_context(
                &uniform_mutation_plan(Box::new(StopThenSafeMutationBaseDist)),
                sim.clone(),
                seed,
                Some(cfg),
                Some(contracts),
            );
            let site = match outcome.trace.find("mutate.uniform.site[0]") {
                Some(rec) => match rec.value {
                    ChoiceValue::Int(s) => s as u32,
                    _ => continue,
                },
                None => continue,
            };
            if site != target_site {
                continue;
            }
            // Also verify the unconstrained run at the same seed
            // produces a stop violation — so the comparison test can
            // demonstrate the contract is what diverts the outcome.
            let unconstrained = PassRuntime::execute_with_context(
                &uniform_mutation_plan(Box::new(StopThenSafeMutationBaseDist)),
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
            "no seed in search range produced constrained site {} with unconstrained stop violation",
            target_site
        );
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
        let region = Region::new(Segment::V, NucHandle::new(0), NucHandle::new(9));
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
            compute_codon_rail(&final_sim.sequence.regions[0], &final_sim.pool).amino_acids,
            compute_codon_rail(&sim.sequence.regions[0], &sim.pool).amino_acids
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
        let region = Region::new(Segment::V, NucHandle::new(0), NucHandle::new(12));
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
    fn uniform_mutation_no_contracts_rng_and_trace_shape_are_characterized() {
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

        let pass = UniformMutationPass::new(
            Box::new(EmpiricalLengthDist::from_pairs(vec![(3, 1.0)])),
            Box::new(UniformBase),
        );
        let mut trace = Trace::new();
        let mut rng = Rng::new(0x1234);
        let final_sim = {
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
            pass.execute(&sim, &mut ctx)
        };

        // Fixed count distribution consumes one f64. Each mutation
        // consumes one word for the site and one word for the base.
        assert_eq!(rng.words_consumed(), 1 + 3 * 2);
        assert_eq!(trace.len(), 1 + 3 * 2);
        assert_eq!(
            trace.find("mutate.uniform.count").unwrap().value,
            ChoiceValue::Int(3)
        );
        for i in 0..3 {
            assert!(trace.find(&format!("mutate.uniform.site[{}]", i)).is_some());
            assert!(trace.find(&format!("mutate.uniform.base[{}]", i)).is_some());
        }
        assert_eq!(final_sim.mutation_count, 3);
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
    fn uniform_mutation_pass_codon_rail_recompute_is_deterministic() {
        // Codon-rail data isn't cached on `Region` — every reader
        // calls `compute_codon_rail` against the current pool. This
        // test verifies the pure function is deterministic: two
        // consecutive computes against the post-mutation pool produce
        // identical results.
        let mut sim = Simulation::new();
        for (i, b) in b"ATGGGGAAATTT".iter().enumerate() {
            let (next, _) =
                sim.with_nucleotide_pushed(Nucleotide::germline(*b, i as u16, Segment::V));
            sim = next;
        }
        let region = Region::new(Segment::V, NucHandle::new(0), NucHandle::new(12));
        sim = sim.with_region_added(region);

        let mut plan = PassPlan::new();
        plan.push(Box::new(UniformMutationPass::new(
            Box::new(EmpiricalLengthDist::from_pairs(vec![(8, 1.0)])),
            Box::new(UniformBase),
        )));

        let outcome = PassRuntime::execute(&plan, sim, 99);
        let final_sim = outcome.final_simulation();

        // Verify that an on-demand recompute against the post-
        // mutation pool is deterministic.
        let a = compute_codon_rail(&final_sim.sequence.regions[0], &final_sim.pool);
        let b = compute_codon_rail(&final_sim.sequence.regions[0], &final_sim.pool);
        assert_eq!(a.amino_acids, b.amino_acids);
        assert_eq!(a.stop_codon_positions, b.stop_codon_positions);
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

        assert_eq!(
            pass.declared_choice_patterns(),
            vec![
                address::ChoiceAddressPattern::MutateUniformCount,
                address::ChoiceAddressPattern::MutateUniformSite,
                address::ChoiceAddressPattern::MutateUniformBase,
            ]
        );
    }

    #[test]
    fn uniform_mutation_productive_filters_base_that_would_create_stop() {
        let (cfg, sim) = make_substitution_productive_vj_fixture();
        let contracts = productive();
        // Under v3.0 constrain-before-propose, site selection is
        // weighted by per-site admissible mass — the seed → site
        // map differs from raw uniform sampling, so the search has
        // to run through the constrained path.
        let seed = find_seed_for_constrained_uniform_mutation_site(&cfg, &contracts, &sim, 2);

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

    /// Custom base distribution that draws uniformly but does NOT
    /// expose `support()`. Used to exercise the
    /// constrain-before-propose helper's behavior when the
    /// distribution can't be enumerated — under v3.0 that path must
    /// either skip (permissive) or error (strict), never fall
    /// through to the reject-after-propose loop.
    #[derive(Clone, Debug)]
    struct NoSupportBaseDist;

    impl Distribution for NoSupportBaseDist {
        type Output = u8;

        fn sample(&self, rng: &mut crate::rng::Rng) -> u8 {
            const BASES: [u8; 4] = [b'A', b'C', b'G', b'T'];
            BASES[rng.range_u32(4) as usize]
        }

        // Deliberately no `support()` override — returns `None` via
        // the trait default. This is the architectural escape hatch
        // we are pinning down.
    }

    #[test]
    fn uniform_mutation_support_unavailable_under_contracts_permissive_skips_slot() {
        let (cfg, sim) = make_substitution_productive_vj_fixture();
        let contracts = productive();
        let outcome = PassRuntime::execute_with_context(
            &uniform_mutation_plan(Box::new(NoSupportBaseDist)),
            sim.clone(),
            0,
            Some(&cfg),
            Some(&contracts),
        );
        // No per-event trace entries: the slot was skipped, not
        // sampled via a contract-blind fallback.
        assert!(outcome.trace.find("mutate.uniform.site[0]").is_none());
        assert!(outcome.trace.find("mutate.uniform.base[0]").is_none());
        // Realized-count semantics: zero mutations applied.
        assert_eq!(outcome.final_simulation().mutation_count, 0);
        // Pool unchanged.
        for i in 0..sim.pool.len() {
            let h = NucHandle::new(i as u32);
            assert_eq!(
                outcome.final_simulation().pool.get(h).unwrap().base,
                sim.pool.get(h).unwrap().base,
            );
        }
    }

    #[test]
    fn uniform_mutation_support_unavailable_under_contracts_strict_errors() {
        let (cfg, sim) = make_substitution_productive_vj_fixture();
        let contracts = productive();
        let err = PassRuntime::execute_strict_with_context(
            &uniform_mutation_plan(Box::new(NoSupportBaseDist)),
            sim,
            0,
            Some(&cfg),
            Some(&contracts),
        )
        .unwrap_err();
        assert_eq!(err.pass_name(), "mutate.uniform");
        assert_eq!(err.address(), "mutate.uniform.base[0]");
        assert_eq!(
            err.constraint_reason(),
            Some(FilteredSampleError::SupportUnavailable)
        );
    }

    #[test]
    fn uniform_mutation_strict_errors_when_base_filter_empty() {
        // v3.0 constrain-before-propose: strict-mode
        // `EmptyAdmissibleSupport` fires when *no* (site, base)
        // combination across the pool admits a draw — not just the
        // one site the seed-uniform path happened to pick. The
        // locked-down V=TGG / J=TGG fixture has every junction site
        // restricted to `{T}` or `{G}`, so a dist whose support is
        // `{A}` is rejected everywhere.
        let (cfg, sim) = make_fully_locked_vj_fixture();
        let contracts = productive();

        let err = PassRuntime::execute_strict_with_context(
            &uniform_mutation_plan(Box::new(StopOnlyMutationBaseDist)),
            sim,
            0,
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
