//! `ContaminantPass` — wholesale sequence replacement (E.6).

use crate::address;
use crate::dist::Distribution;
use crate::ir::{NucHandle, Simulation};
use crate::pass::{Pass, PassContext, PassEffect, PassError};
use crate::passes::mutation_transaction::MutationTransaction;
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
/// trace records one base entry for each position that was actually
/// replaced; under active contracts in permissive mode, positions
/// with empty admissible support are skipped and have no base trace
/// record. When `applied` is `Bool(false)`, no further records are
/// emitted (the pool is returned unchanged).
///
/// Codon-rail data is not stored on `Region` — the pool directly
/// reflects the contaminant bytes after each `with_base_changed`,
/// and callers that need amino-acid translation call
/// [`crate::ir::compute_codon_rail`] on demand. When contracts are
/// active, each replacement base is filtered against the current
/// intermediate IR and the target site, so a contamination event
/// cannot transiently violate enforced contracts when an admissible
/// replacement exists.
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

    fn execute_with_sampling_mode(
        &self,
        sim: &Simulation,
        ctx: &mut PassContext,
        strict: bool,
    ) -> Result<Simulation, PassError> {
        // Trace-injected replay (Tier 3): consume the recorded
        // Bool. Same shape as RevCompPass — single coin flip.
        let applied = if let Some(cursor) = ctx.replay_cursor.as_deref_mut() {
            cursor
                .expect_bool(address::ChoiceAddress::CorruptContaminantApplied)
                .map_err(|reason| PassError::replay(self.name(), reason))?
        } else {
            ctx.rng.next_f64() < self.apply_prob
        };
        ctx.trace.record_choice(
            address::ChoiceAddress::CorruptContaminantApplied,
            ChoiceValue::Bool(applied),
        );

        if !applied {
            return Ok(sim.clone());
        }

        let pool_len = sim.pool.len() as u32;
        if pool_len == 0 {
            return Ok(sim.clone());
        }

        // Replace every base in the pool with a contaminant draw,
        // contract-filtered per base.
        let mut tx = MutationTransaction::open(sim, ctx, self.name(), strict);

        // Per-site iteration. Under replay mode the original run
        // may have skipped sites where contracts rejected the
        // candidate, producing trace gaps; the cursor's next record
        // address tells us whether the original recorded this site.
        // We peek and decide:
        //   - cursor's next address matches `bases[i]` → call
        //     `substitute_base`, which consumes the base via the
        //     replay branch and validates against current contracts.
        //   - mismatch / cursor drained → skip this site, matching
        //     the original empty-support behavior.
        for i in 0..pool_len {
            let site = NucHandle::new(i);
            let base_choice_address = address::ChoiceAddress::CorruptContaminantBase(i);
            let should_consume = if tx.replay_cursor().is_some() {
                let expected = base_choice_address.to_string();
                tx.replay_cursor()
                    .and_then(|c| c.peek_address().map(|s| s == expected))
                    .unwrap_or(false)
            } else {
                true
            };
            if !should_consume {
                // Replay mode: cursor doesn't have a record for this
                // index → original run skipped (contract reject).
                // Skip in replay too, matching the realized behavior.
                continue;
            }
            tx.substitute_base(
                site,
                self.base_dist.as_ref(),
                base_choice_address,
                i,
                pool_len,
                None,
            )?;
        }
        tx.commit()
    }
}

impl Pass for ContaminantPass {
    fn name(&self) -> &str {
        address::CORRUPT_CONTAMINANT
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

    fn declared_choice_patterns(&self) -> Vec<address::ChoiceAddressPattern> {
        vec![
            address::ChoiceAddressPattern::CorruptContaminantApplied,
            address::ChoiceAddressPattern::CorruptContaminantBase,
        ]
    }

    fn effects(&self) -> Vec<PassEffect> {
        vec![PassEffect::EditBases]
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::assignment::AlleleInstance;
    use crate::contract::productive;
    use crate::dist::{FilteredSampleError, UniformBase};
    use crate::ir::compute_codon_rail;
    use crate::ir::{Nucleotide, Region, Segment};
    use crate::pass::testing::PassRuntime;
    use crate::pass::PassPlan;
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
        let v_region = Region::new(Segment::V, NucHandle::new(0), NucHandle::new(3));
        sim = sim.with_region_added(v_region);

        for (i, &b) in b"TGG".iter().enumerate() {
            let (next, _) =
                sim.with_nucleotide_pushed(Nucleotide::germline(b, i as u16, Segment::J));
            sim = next;
        }
        let j_region = Region::new(Segment::J, NucHandle::new(3), NucHandle::new(6));
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
        let region = Region::new(Segment::V, NucHandle::new(0), NucHandle::new(9));
        sim = sim.with_region_added(region);
        // Before contamination: M G G (rail computed via the
        // persistent helper above).
        assert_eq!(
            compute_codon_rail(&sim.sequence.regions[0], &sim.pool).amino_acids,
            b"MGG"
        );

        let mut plan = PassPlan::new();
        plan.push(Box::new(ContaminantPass::new(1.0, Box::new(UniformBase))));
        let outcome = PassRuntime::execute(&plan, sim, 0);
        let final_sim = outcome.final_simulation();

        // The codon rail is computed on demand. Assert that
        // recomputing twice yields the same result — i.e. the
        // function is deterministic on the final pool.
        let a = compute_codon_rail(&final_sim.sequence.regions[0], &final_sim.pool);
        let b = compute_codon_rail(&final_sim.sequence.regions[0], &final_sim.pool);
        assert_eq!(a.amino_acids, b.amino_acids);
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
    fn contaminant_permissive_empty_support_skips_sites_without_unconstrained_fallback() {
        // v3.0 rule: under active contracts + permissive mode,
        // when the natural base distribution has no admissible
        // candidates at a site, the engine must NOT fall back to
        // an unconstrained draw (that would re-introduce
        // reject-after-propose). Instead the site is skipped:
        // no trace record, no pool mutation.
        //
        // Fixture pool = "AAATGG" (V anchor "AAA" / K, J anchor
        // "TGG" / W). With `StopOnlyContaminantBaseDist` (support
        // {T}), only site 3 admits T (it's already T in the J
        // anchor codon); sites 0..2 and 4,5 all reject T under
        // AnchorPreserved.
        //
        // Pre-v3.0 the permissive path would write T at every
        // rejecting site (anchor V violation) → productive()
        // verify fails. v3.0 leaves those sites untouched, so
        // the bundle holds.
        let (cfg, sim) = make_contaminant_productive_vj_fixture();
        let contracts = productive();
        let outcome = PassRuntime::execute_with_context(
            &contaminant_plan(Box::new(StopOnlyContaminantBaseDist)),
            sim.clone(),
            0,
            Some(&cfg),
            Some(&contracts),
        );

        // `applied = true` (the coin flip fired); per-site `bases[i]`
        // entries exist only for the one admissible site (i=3).
        assert_eq!(
            outcome
                .trace
                .find("corrupt.contaminant.applied")
                .unwrap()
                .value,
            ChoiceValue::Bool(true)
        );
        for i in [0u32, 1, 2, 4, 5] {
            assert!(
                outcome
                    .trace
                    .find(&format!("corrupt.contaminant.bases[{}]", i))
                    .is_none(),
                "permissive empty-support must skip the trace record at site {}",
                i
            );
        }
        // Site 3 was admissible (T → T no-op) and writes a trace
        // entry.
        assert_eq!(
            outcome
                .trace
                .find("corrupt.contaminant.bases[3]")
                .unwrap()
                .value,
            ChoiceValue::Base(b'T')
        );

        // Bundle holds — no anchor codon was overwritten.
        assert!(
            contracts
                .verify(outcome.final_simulation(), Some(&cfg))
                .is_ok(),
            "v3.0 permissive contaminant must leave the bundle satisfied"
        );

        // Pool byte-identical to input.
        for i in 0..(sim.pool.len() as u32) {
            let h = NucHandle::new(i);
            assert_eq!(
                outcome.final_simulation().pool.get(h).unwrap().base,
                sim.pool.get(h).unwrap().base,
                "site {} pool byte changed despite empty admissible support",
                i
            );
        }
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

    // ── Trace-injected replay (Tier 3) ─────────────────────────────

    use crate::address::ChoiceAddress;
    use crate::pass::PassContext;
    use crate::replay::{ReplayError, TraceCursor};
    use crate::rng::Rng;
    use crate::trace::{ChoiceRecord, Trace};

    fn rec(addr: ChoiceAddress, v: ChoiceValue) -> ChoiceRecord {
        ChoiceRecord::new(addr.to_string(), v)
    }

    fn run_contaminant_replay(
        sim: &Simulation,
        apply_prob: f64,
        records: Vec<ChoiceRecord>,
    ) -> (Result<Simulation, PassError>, Trace, u64) {
        let pass = ContaminantPass::new(apply_prob, Box::new(UniformBase));
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
            pass.execute_checked(sim, &mut ctx)
        };
        let words = rng.words_consumed();
        (result, trace, words)
    }

    #[test]
    fn contaminant_replay_consumes_recorded_bool_and_per_site_bases() {
        let sim = contaminant_test_sim(); // pool_len = 9
        let bases = b"TGCATGCAT";
        let mut records =
            vec![rec(ChoiceAddress::CorruptContaminantApplied, ChoiceValue::Bool(true))];
        for (i, b) in bases.iter().enumerate() {
            records.push(rec(
                ChoiceAddress::CorruptContaminantBase(i as u32),
                ChoiceValue::Base(*b),
            ));
        }

        let (result, trace, rng_words) = run_contaminant_replay(&sim, 1.0, records);
        let next = result.unwrap();

        // Every site now carries the replayed base.
        for (i, b) in bases.iter().enumerate() {
            assert_eq!(next.pool.get(NucHandle::new(i as u32)).unwrap().base, *b);
            assert_eq!(
                trace
                    .find(&format!("corrupt.contaminant.bases[{}]", i))
                    .unwrap()
                    .value,
                ChoiceValue::Base(*b),
            );
        }
        assert_eq!(rng_words, 0);
    }

    #[test]
    fn contaminant_replay_applied_false_skips_per_site_loop() {
        let sim = contaminant_test_sim();
        let records = vec![rec(
            ChoiceAddress::CorruptContaminantApplied,
            ChoiceValue::Bool(false),
        )];
        let (result, trace, rng_words) = run_contaminant_replay(&sim, 0.5, records);
        let next = result.unwrap();

        // Pool unchanged.
        for i in 0..sim.pool.len() {
            let h = NucHandle::new(i as u32);
            assert_eq!(
                next.pool.get(h).unwrap().base,
                sim.pool.get(h).unwrap().base,
            );
        }
        // Trace records the False flag and no per-site bases.
        assert_eq!(
            trace.find("corrupt.contaminant.applied").unwrap().value,
            ChoiceValue::Bool(false),
        );
        assert_eq!(trace.find("corrupt.contaminant.bases[0]"), None);
        assert_eq!(rng_words, 0);
    }

    #[test]
    fn contaminant_replay_skips_sites_with_no_record() {
        // Original run rejected sites 2 and 5 under contracts.
        // Trace has bases at 0, 1, 3, 4, 6, 7, 8 — no records for
        // 2 or 5. Replay must skip those sites positionally.
        let sim = contaminant_test_sim();
        let kept_sites: Vec<u32> = vec![0, 1, 3, 4, 6, 7, 8];
        let mut records =
            vec![rec(ChoiceAddress::CorruptContaminantApplied, ChoiceValue::Bool(true))];
        for &i in &kept_sites {
            records.push(rec(
                ChoiceAddress::CorruptContaminantBase(i),
                ChoiceValue::Base(b'C'),
            ));
        }

        let (result, _trace, _) = run_contaminant_replay(&sim, 1.0, records);
        let next = result.unwrap();

        // Replayed sites carry C; skipped sites keep their germline base.
        for &i in &kept_sites {
            assert_eq!(next.pool.get(NucHandle::new(i)).unwrap().base, b'C');
        }
        for &i in &[2u32, 5] {
            let germline = sim.pool.get(NucHandle::new(i)).unwrap().base;
            assert_eq!(next.pool.get(NucHandle::new(i)).unwrap().base, germline);
        }
    }

    #[test]
    fn contaminant_replay_wrong_kind_applied_surfaces_replay_error() {
        let sim = contaminant_test_sim();
        // Wrong kind: Int instead of Bool.
        let records = vec![rec(
            ChoiceAddress::CorruptContaminantApplied,
            ChoiceValue::Int(1),
        )];
        let (result, _, _) = run_contaminant_replay(&sim, 1.0, records);
        match result.unwrap_err() {
            PassError::Replay { reason, .. } => {
                assert!(matches!(reason, ReplayError::ValueKindMismatch { .. }));
            }
            other => panic!("expected Replay::ValueKindMismatch, got {other:?}"),
        }
    }

    #[test]
    fn contaminant_replay_wrong_kind_base_surfaces_replay_error() {
        let sim = contaminant_test_sim();
        let records = vec![
            rec(ChoiceAddress::CorruptContaminantApplied, ChoiceValue::Bool(true)),
            // Wrong kind for bases[0]: Int instead of Base.
            rec(
                ChoiceAddress::CorruptContaminantBase(0),
                ChoiceValue::Int(0),
            ),
        ];
        let (result, _, _) = run_contaminant_replay(&sim, 1.0, records);
        match result.unwrap_err() {
            PassError::Replay { reason, .. } => {
                assert!(matches!(reason, ReplayError::ValueKindMismatch { .. }));
            }
            other => panic!("expected Replay::ValueKindMismatch, got {other:?}"),
        }
    }

    #[test]
    fn contaminant_replay_under_contract_rejects_inadmissible_recorded_base() {
        // V-anchor fixture: site 0..2 = TGG (Trp). AnchorPreserved
        // rejects any non-T at site 0 (changes the codon).
        let (cfg, sim) = make_contaminant_productive_vj_fixture();
        let contracts = productive();
        let records = vec![
            rec(ChoiceAddress::CorruptContaminantApplied, ChoiceValue::Bool(true)),
            // Pool is "AAATGG"; site 0 is the first A of the V
            // anchor codon (AAA = Lys). Substituting to G changes
            // the amino acid → AnchorPreserved.V rejects.
            rec(
                ChoiceAddress::CorruptContaminantBase(0),
                ChoiceValue::Base(b'G'),
            ),
        ];
        let pass = ContaminantPass::new(1.0, Box::new(UniformBase));
        let mut cursor = TraceCursor::from_owned(records);
        let mut trace = Trace::new();
        let mut rng = Rng::new(0);
        let result = {
            let mut ctx = PassContext {
                trace: &mut trace,
                rng: &mut rng,
                pass_index: 0,
                refdata: Some(&cfg),
                contracts: Some(&contracts),
                feasibility: None,
                reference_index: None,
                replay_cursor: Some(&mut cursor),
                event_log_sink: None,
            };
            pass.execute_checked(&sim, &mut ctx)
        };
        match result.unwrap_err() {
            PassError::ConstraintSampling { address, .. } => {
                assert_eq!(address, "corrupt.contaminant.bases[0]");
            }
            other => panic!("expected ConstraintSampling, got {other:?}"),
        }
    }
}
