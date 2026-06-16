//! Concrete pass implementations.
//!
//! Each biology operation lives in its own submodule; this module
//! re-exports the pass types so external code keeps the flat
//! `crate::passes::FooPass` import surface.
//!
//! Submodule layout:
//! - [`echo`] — `EchoPass`, the deterministic transform reference.
//! - [`sample_base`] — `SampleBasePass`, the sampling reference.
//! - [`sample_allele`] — V/D/J allele sampling.
//! - [`trim`] — recombination trim sampling.
//! - [`assemble_segment`] — copy a germline allele slice into
//!   the pool.
//! - [`generate_np`] — TdT-like N-nucleotide region generation
//!  .
//! - [`mutate`] — SHM passes (uniform + S5F).
//! - [`corrupt`] — observation-stage perturbations (PCR error,
//!   quality error, contamination, indels).

pub mod assemble_segment;
pub mod corrupt;
pub(crate) mod count_source;
pub mod echo;
pub mod generate_np;
pub mod invert_d;
pub mod mutate;
pub(crate) mod mutation_transaction;
pub mod p_addition;
pub mod paired_end;
pub(crate) mod paramsig;
pub mod receptor_revision;
pub mod sample_allele;
pub mod sample_base;
pub mod sample_haplotype;
pub mod trim;

#[cfg(test)]
pub(crate) mod test_support;

pub use assemble_segment::AssembleSegmentPass;
pub use corrupt::{
    ContaminantPass, EndLossPass, IndelPass, LossEnd, NCorruptionPass, PCRErrorPass,
    QualityErrorPass, RevCompPass,
};
pub use echo::EchoPass;
pub use generate_np::GenerateNPPass;
pub use invert_d::InvertDPass;
pub use mutate::{S5FMutationPass, UniformMutationPass};
pub use p_addition::PAdditionPass;
pub use paired_end::{PairedEndLayoutSpec, PairedEndSamplingPass};
pub use receptor_revision::ReceptorRevisionPass;
pub use sample_allele::SampleAllelePass;
pub use sample_base::SampleBasePass;
pub use trim::TrimPass;

// ──────────────────────────────────────────────────────────────────
// Cross-cutting integration tests (D.6 productive bundle)
// ──────────────────────────────────────────────────────────────────

#[cfg(test)]
mod tests {
    use super::*;
    use crate::assignment::TrimEnd;
    use crate::contract::{
        productive, ChoiceContext, Contract, ContractSet, ContractViolation,
        ProductiveJunctionFrame,
    };
    use crate::dist::{
        AllelePoolDist, Distribution, EmpiricalLengthDist, FilteredSampleError, UniformBase,
    };
    use crate::ir::{NucHandle, Segment, Simulation};
    use crate::pass::testing::PassRuntime;
    use crate::pass::{Pass, PassPlan};
    use crate::passes::sample_allele::test_support::make_test_pool;
    use crate::s5f::{S5FKernel, S5F_NUM_CONTEXTS, S5F_SUBSTITUTION_LEN};
    use crate::trace::ChoiceValue;

    fn fixed_count_dist() -> Box<dyn Distribution<Output = i64>> {
        Box::new(EmpiricalLengthDist::from_pairs(vec![(1, 1.0)]))
    }

    fn test_s5f_kernel() -> S5FKernel {
        S5FKernel::new(
            vec![1.0; S5F_NUM_CONTEXTS],
            vec![0.25; S5F_SUBSTITUTION_LEN],
        )
    }

    fn assert_typed_patterns_project_to_declared_choices(pass: &dyn Pass) {
        let projected: Vec<String> = pass
            .declared_choice_patterns()
            .into_iter()
            .map(String::from)
            .collect();
        assert_eq!(
            projected,
            pass.declared_choices(),
            "{} typed declared-choice patterns must project to the stable string surface",
            pass.name()
        );
    }

    #[test]
    fn built_in_pass_declared_choice_patterns_match_string_surface() {
        let v_pool = make_test_pool(2, Segment::V);
        let d_pool = make_test_pool(2, Segment::D);
        let j_pool = make_test_pool(2, Segment::J);

        let passes: Vec<Box<dyn Pass>> = vec![
            Box::new(EchoPass::new(b'A', 0, Segment::V)),
            Box::new(AssembleSegmentPass::new(Segment::V)),
            Box::new(SampleAllelePass::new(
                Segment::V,
                Box::new(AllelePoolDist::uniform(&v_pool)),
            )),
            Box::new(SampleAllelePass::new(
                Segment::D,
                Box::new(AllelePoolDist::uniform(&d_pool)),
            )),
            Box::new(SampleAllelePass::new(
                Segment::J,
                Box::new(AllelePoolDist::uniform(&j_pool)),
            )),
            Box::new(TrimPass::new(Segment::V, TrimEnd::Five, fixed_count_dist())),
            Box::new(TrimPass::new(
                Segment::D,
                TrimEnd::Three,
                fixed_count_dist(),
            )),
            Box::new(TrimPass::new(Segment::J, TrimEnd::Five, fixed_count_dist())),
            Box::new(GenerateNPPass::new(
                Segment::Np1,
                fixed_count_dist(),
                Box::new(UniformBase),
            )),
            Box::new(GenerateNPPass::new(
                Segment::Np2,
                fixed_count_dist(),
                Box::new(UniformBase),
            )),
            Box::new(UniformMutationPass::new(
                fixed_count_dist(),
                Box::new(UniformBase),
            )),
            Box::new(S5FMutationPass::new(test_s5f_kernel(), fixed_count_dist())),
            Box::new(PCRErrorPass::new(fixed_count_dist(), Box::new(UniformBase))),
            Box::new(QualityErrorPass::new(
                fixed_count_dist(),
                Box::new(UniformBase),
            )),
            Box::new(ContaminantPass::new(0.5, Box::new(UniformBase))),
            Box::new(IndelPass::new(
                fixed_count_dist(),
                0.5,
                Box::new(UniformBase),
            )),
            Box::new(NCorruptionPass::new(fixed_count_dist())),
            Box::new(EndLossPass::new(LossEnd::Five, fixed_count_dist())),
            Box::new(EndLossPass::new(LossEnd::Three, fixed_count_dist())),
            Box::new(RevCompPass::new(0.5)),
        ];

        for pass in passes {
            assert_typed_patterns_project_to_declared_choices(pass.as_ref());
        }
    }

    /// Build a synthetic VJ refdata: V "AAACCCGGG" (9bp, anchor 6),
    /// J "TTTAAA" (6bp, anchor 0). Junction = V_anchor_to_end (3bp)
    /// + NP1 + J_anchor_to_W3 (3bp). For productive frame:
    /// (3 + NP1 + 3) % 3 == 0 → NP1 % 3 == 0.
    fn make_vj_refdata_for_filter() -> crate::refdata::RefDataConfig {
        let mut cfg = crate::refdata::RefDataConfig::empty(crate::refdata::ChainType::Vj);
        let _ = cfg.v_pool.push(crate::refdata::Allele {
            name: "v_test*01".into(),
            gene: "v_test".into(),
            seq: b"AAACCCGGG".to_vec(),
            segment: Segment::V,
            anchor: Some(6),
            functional_status: None,
            subregions: Vec::new(),
        });
        let _ = cfg.j_pool.push(crate::refdata::Allele {
            name: "j_test*01".into(),
            gene: "j_test".into(),
            seq: b"TTTAAA".to_vec(),
            segment: Segment::J,
            anchor: Some(0),
            functional_status: None,
            subregions: Vec::new(),
        });
        cfg
    }

    #[test]
    fn productive_admits_filters_np1_for_vj_chain() {
        let cfg = make_vj_refdata_for_filter();
        let dist = EmpiricalLengthDist::from_pairs((0..7).map(|i| (i, 1.0)).collect::<Vec<_>>());

        let mut plan = PassPlan::new();
        plan.push(Box::new(SampleAllelePass::new(
            Segment::V,
            Box::new(AllelePoolDist::uniform(&cfg.v_pool)),
        )));
        plan.push(Box::new(SampleAllelePass::new(
            Segment::J,
            Box::new(AllelePoolDist::uniform(&cfg.j_pool)),
        )));
        plan.push(Box::new(AssembleSegmentPass::new(Segment::V)));
        plan.push(Box::new(GenerateNPPass::new(
            Segment::Np1,
            Box::new(dist),
            Box::new(UniformBase),
        )));

        let contracts = ContractSet::new().with(Box::new(ProductiveJunctionFrame::new()));

        for seed in 0..20u64 {
            let outcome = PassRuntime::execute_with_context(
                &plan,
                Simulation::new(),
                seed,
                Some(&cfg),
                Some(&contracts),
            );
            let np1_len = match outcome.trace.find("np.np1.length").unwrap().value {
                ChoiceValue::Int(n) => n,
                _ => panic!("wrong variant"),
            };
            assert!(
                np1_len % 3 == 0,
                "seed {} produced NP1 length {} (not divisible by 3)",
                seed,
                np1_len
            );
        }
    }

    #[test]
    fn productive_admits_unconstrained_without_contracts() {
        let cfg = make_vj_refdata_for_filter();
        let dist = EmpiricalLengthDist::from_pairs((0..7).map(|i| (i, 1.0)).collect::<Vec<_>>());

        let mut plan = PassPlan::new();
        plan.push(Box::new(SampleAllelePass::new(
            Segment::V,
            Box::new(AllelePoolDist::uniform(&cfg.v_pool)),
        )));
        plan.push(Box::new(SampleAllelePass::new(
            Segment::J,
            Box::new(AllelePoolDist::uniform(&cfg.j_pool)),
        )));
        plan.push(Box::new(AssembleSegmentPass::new(Segment::V)));
        plan.push(Box::new(GenerateNPPass::new(
            Segment::Np1,
            Box::new(dist),
            Box::new(UniformBase),
        )));

        // Without contracts: any value in [0, 6] is allowed.
        let mut seen_out_of_frame = false;
        for seed in 0..50u64 {
            let outcome = PassRuntime::execute_with_refdata(&plan, Simulation::new(), seed, &cfg);
            let np1_len = match outcome.trace.find("np.np1.length").unwrap().value {
                ChoiceValue::Int(n) => n,
                _ => panic!("wrong variant"),
            };
            assert!(np1_len >= 0 && np1_len <= 6);
            if np1_len % 3 != 0 {
                seen_out_of_frame = true;
            }
        }
        assert!(
            seen_out_of_frame,
            "Without contracts, expected at least one out-of-frame NP1 sample"
        );
    }

    #[test]
    fn productive_admits_empty_filter_returns_explicit_zero_length() {
        // v3.0 rule: when natural ∩ admissible is empty under
        // active contracts, permissive mode must NOT fall back
        // to an unconstrained draw from the natural support
        // (that would re-introduce reject-after-propose). The
        // architectural no-op for a length sampler is 0 — same
        // shape as `TrimPass::sample_trim`'s empty-support
        // fallback. (Pre-v3.0 the pass fell through to
        // `length_dist.sample(rng)` and produced one of
        // {1,2,4,5}, all of which violate frame in this fixture.)
        let cfg = make_vj_refdata_for_filter();
        let dist = EmpiricalLengthDist::from_pairs(vec![(1, 1.0), (2, 1.0), (4, 1.0), (5, 1.0)]);

        let mut plan = PassPlan::new();
        plan.push(Box::new(SampleAllelePass::new(
            Segment::V,
            Box::new(AllelePoolDist::uniform(&cfg.v_pool)),
        )));
        plan.push(Box::new(SampleAllelePass::new(
            Segment::J,
            Box::new(AllelePoolDist::uniform(&cfg.j_pool)),
        )));
        plan.push(Box::new(AssembleSegmentPass::new(Segment::V)));
        plan.push(Box::new(GenerateNPPass::new(
            Segment::Np1,
            Box::new(dist),
            Box::new(UniformBase),
        )));

        let contracts = ContractSet::new().with(Box::new(ProductiveJunctionFrame::new()));

        let outcome = PassRuntime::execute_with_context(
            &plan,
            Simulation::new(),
            0,
            Some(&cfg),
            Some(&contracts),
        );
        let np1_len = match outcome.trace.find("np.np1.length").unwrap().value {
            ChoiceValue::Int(n) => n,
            _ => panic!("wrong variant"),
        };
        assert_eq!(np1_len, 0);
    }

    #[derive(Clone, Debug)]
    struct UnenumerableLengthDist;

    impl Distribution for UnenumerableLengthDist {
        type Output = i64;

        fn sample(&self, _rng: &mut crate::rng::Rng) -> i64 {
            0
        }
    }

    #[test]
    fn productive_strict_errors_when_length_filter_empty() {
        let cfg = make_vj_refdata_for_filter();
        let dist = EmpiricalLengthDist::from_pairs(vec![(1, 1.0), (2, 1.0), (4, 1.0), (5, 1.0)]);

        let mut plan = PassPlan::new();
        plan.push(Box::new(SampleAllelePass::new(
            Segment::V,
            Box::new(AllelePoolDist::uniform(&cfg.v_pool)),
        )));
        plan.push(Box::new(SampleAllelePass::new(
            Segment::J,
            Box::new(AllelePoolDist::uniform(&cfg.j_pool)),
        )));
        plan.push(Box::new(AssembleSegmentPass::new(Segment::V)));
        plan.push(Box::new(GenerateNPPass::new(
            Segment::Np1,
            Box::new(dist),
            Box::new(UniformBase),
        )));

        let contracts = ContractSet::new().with(Box::new(ProductiveJunctionFrame::new()));
        let err = PassRuntime::execute_strict_with_context(
            &plan,
            Simulation::new(),
            0,
            Some(&cfg),
            Some(&contracts),
        )
        .unwrap_err();

        assert_eq!(err.pass_name(), "generate_np.np1");
        assert_eq!(err.address(), "np.np1.length");
        assert_eq!(
            err.constraint_reason(),
            Some(FilteredSampleError::EmptyAdmissibleSupport)
        );
        assert!(err.to_string().contains("no admissible candidates"));
    }

    #[test]
    fn productive_strict_errors_when_length_support_unavailable() {
        let cfg = make_vj_refdata_for_filter();

        let mut plan = PassPlan::new();
        plan.push(Box::new(SampleAllelePass::new(
            Segment::V,
            Box::new(AllelePoolDist::uniform(&cfg.v_pool)),
        )));
        plan.push(Box::new(SampleAllelePass::new(
            Segment::J,
            Box::new(AllelePoolDist::uniform(&cfg.j_pool)),
        )));
        plan.push(Box::new(AssembleSegmentPass::new(Segment::V)));
        plan.push(Box::new(GenerateNPPass::new(
            Segment::Np1,
            Box::new(UnenumerableLengthDist),
            Box::new(UniformBase),
        )));

        let contracts = ContractSet::new().with(Box::new(ProductiveJunctionFrame::new()));
        let err = PassRuntime::execute_strict_with_context(
            &plan,
            Simulation::new(),
            0,
            Some(&cfg),
            Some(&contracts),
        )
        .unwrap_err();

        assert_eq!(err.pass_name(), "generate_np.np1");
        assert_eq!(err.address(), "np.np1.length");
        assert_eq!(
            err.constraint_reason(),
            Some(FilteredSampleError::SupportUnavailable)
        );
    }

    struct RejectNpBases;

    impl Contract for RejectNpBases {
        fn name(&self) -> &str {
            "reject_np_bases"
        }

        fn verify(
            &self,
            _sim: &Simulation,
            _refdata: Option<&crate::refdata::RefDataConfig>,
        ) -> Result<(), ContractViolation> {
            Ok(())
        }

        // Post-flip the trait's primary surface is `admits_typed`,
        // so this test contract dispatches on `ChoiceAddress::NpBase`
        // directly rather than prefix-matching the legacy string.
        fn admits_typed(
            &self,
            _sim: &Simulation,
            _refdata: Option<&crate::refdata::RefDataConfig>,
            context: ChoiceContext<'_>,
            _candidate: &ChoiceValue,
        ) -> Result<(), ContractViolation> {
            if matches!(
                context.address,
                Some(crate::address::ChoiceAddress::NpBase {
                    segment: crate::address::NpSegment::Np1,
                    ..
                })
            ) {
                return Err(ContractViolation::new(self.name(), "rejected by test"));
            }
            Ok(())
        }
    }

    #[test]
    fn productive_admits_empty_base_filter_writes_n_sentinel_in_permissive() {
        // v3.0 rule (mirror of `sample_length`'s explicit-zero
        // fallback): when the bundle leaves no admissible base
        // for an NP slot, permissive mode writes the IUPAC `N`
        // sentinel rather than falling back to an unconstrained
        // draw. `N` translates to amino acid `X`, which is not
        // a stop and not bound by anchor preservation (NP slots
        // are between V/J anchor codons), so it satisfies the
        // productive bundle's full admit check.
        //
        // The matching strict test below (`productive_strict_errors_…`)
        // pins that strict mode surfaces the same empty support
        // as `EmptyAdmissibleSupport` instead.
        let mut plan = PassPlan::new();
        plan.push(Box::new(GenerateNPPass::new(
            Segment::Np1,
            Box::new(EmpiricalLengthDist::from_pairs(vec![(1, 1.0)])),
            Box::new(UniformBase),
        )));
        let contracts = ContractSet::new().with(Box::new(RejectNpBases));

        let outcome =
            PassRuntime::execute_with_context(&plan, Simulation::new(), 0, None, Some(&contracts));

        // The trace's recorded base for NP slot 0 is `b'N'`.
        assert_eq!(
            outcome.trace.find("np.np1.bases[0]").unwrap().value,
            ChoiceValue::Base(b'N')
        );
        // The pool reflects the same sentinel.
        let final_sim = outcome.final_simulation();
        assert_eq!(final_sim.pool.len(), 1);
        assert_eq!(final_sim.pool.get(NucHandle::new(0)).unwrap().base, b'N');
    }

    #[test]
    fn productive_strict_errors_when_base_filter_empty() {
        let mut plan = PassPlan::new();
        plan.push(Box::new(GenerateNPPass::new(
            Segment::Np1,
            Box::new(EmpiricalLengthDist::from_pairs(vec![(1, 1.0)])),
            Box::new(UniformBase),
        )));
        let contracts = ContractSet::new().with(Box::new(RejectNpBases));

        let err = PassRuntime::execute_strict_with_context(
            &plan,
            Simulation::new(),
            0,
            None,
            Some(&contracts),
        )
        .unwrap_err();

        assert_eq!(err.pass_name(), "generate_np.np1");
        assert_eq!(err.address(), "np.np1.bases[0]");
        assert_eq!(
            err.constraint_reason(),
            Some(FilteredSampleError::EmptyAdmissibleSupport)
        );
    }

    #[test]
    fn productive_strict_succeeds_when_admissible_length_exists() {
        let cfg = make_vj_refdata_for_filter();
        let dist = EmpiricalLengthDist::from_pairs((0..7).map(|i| (i, 1.0)).collect::<Vec<_>>());

        let mut plan = PassPlan::new();
        plan.push(Box::new(SampleAllelePass::new(
            Segment::V,
            Box::new(AllelePoolDist::uniform(&cfg.v_pool)),
        )));
        plan.push(Box::new(SampleAllelePass::new(
            Segment::J,
            Box::new(AllelePoolDist::uniform(&cfg.j_pool)),
        )));
        plan.push(Box::new(AssembleSegmentPass::new(Segment::V)));
        plan.push(Box::new(GenerateNPPass::new(
            Segment::Np1,
            Box::new(dist),
            Box::new(UniformBase),
        )));

        let contracts = ContractSet::new().with(Box::new(ProductiveJunctionFrame::new()));
        let outcome = PassRuntime::execute_strict_with_context(
            &plan,
            Simulation::new(),
            0,
            Some(&cfg),
            Some(&contracts),
        )
        .unwrap();
        let np1_len = match outcome.trace.find("np.np1.length").unwrap().value {
            ChoiceValue::Int(n) => n,
            _ => panic!("wrong variant"),
        };

        assert_eq!(np1_len % 3, 0);
    }

    #[test]
    fn productive_admits_makes_full_pipeline_in_frame_for_vj() {
        use crate::junction::compute_junction;

        let cfg = make_vj_refdata_for_filter();
        let dist = EmpiricalLengthDist::from_pairs((0..7).map(|i| (i, 1.0)).collect::<Vec<_>>());

        let mut plan = PassPlan::new();
        plan.push(Box::new(SampleAllelePass::new(
            Segment::V,
            Box::new(AllelePoolDist::uniform(&cfg.v_pool)),
        )));
        plan.push(Box::new(SampleAllelePass::new(
            Segment::J,
            Box::new(AllelePoolDist::uniform(&cfg.j_pool)),
        )));
        plan.push(Box::new(AssembleSegmentPass::new(Segment::V)));
        plan.push(Box::new(GenerateNPPass::new(
            Segment::Np1,
            Box::new(dist),
            Box::new(UniformBase),
        )));
        plan.push(Box::new(AssembleSegmentPass::new(Segment::J)));

        let contracts = ContractSet::new().with(Box::new(ProductiveJunctionFrame::new()));

        for seed in 0..30u64 {
            let outcome = PassRuntime::execute_with_context(
                &plan,
                Simulation::new(),
                seed,
                Some(&cfg),
                Some(&contracts),
            );
            let junction = compute_junction(outcome.final_simulation(), &cfg)
                .expect("junction should be defined");
            assert!(
                junction.is_in_frame(),
                "seed {} produced out-of-frame junction (length {})",
                seed,
                junction.length
            );
        }
    }

    #[test]
    fn productive_full_bundle_in_frame_and_admits_dispatch_works() {
        use crate::junction::compute_junction;

        let cfg = make_vj_refdata_for_filter();
        let dist = EmpiricalLengthDist::from_pairs((0..10).map(|i| (i, 1.0)).collect::<Vec<_>>());

        let mut plan = PassPlan::new();
        plan.push(Box::new(SampleAllelePass::new(
            Segment::V,
            Box::new(AllelePoolDist::uniform(&cfg.v_pool)),
        )));
        plan.push(Box::new(SampleAllelePass::new(
            Segment::J,
            Box::new(AllelePoolDist::uniform(&cfg.j_pool)),
        )));
        plan.push(Box::new(AssembleSegmentPass::new(Segment::V)));
        plan.push(Box::new(GenerateNPPass::new(
            Segment::Np1,
            Box::new(dist),
            Box::new(UniformBase),
        )));
        plan.push(Box::new(AssembleSegmentPass::new(Segment::J)));

        let contracts = productive();

        for seed in 0..20u64 {
            let outcome = PassRuntime::execute_with_context(
                &plan,
                Simulation::new(),
                seed,
                Some(&cfg),
                Some(&contracts),
            );
            let junction = compute_junction(outcome.final_simulation(), &cfg)
                .expect("junction should be defined");
            assert!(junction.is_in_frame());
        }
    }

    #[test]
    fn sample_allele_pass_replay_in_mixed_plan_with_echo_passes() {
        // Mixed plan: SampleAllele + Echo (transform) interleaved.
        // The trace records only the sampling choices; the IR
        // accumulates both the assignment and the echoed nucleotides.
        let pool = make_test_pool(3, Segment::V);
        let mut plan = PassPlan::new();
        plan.push(Box::new(EchoPass::new(b'A', 0, Segment::V)));
        plan.push(Box::new(SampleAllelePass::new(
            Segment::V,
            Box::new(AllelePoolDist::uniform(&pool)),
        )));
        plan.push(Box::new(EchoPass::new(b'C', 1, Segment::V)));

        let oa = PassRuntime::execute(&plan, Simulation::new(), 0xfeed);
        let ob = PassRuntime::execute(&plan, Simulation::new(), 0xfeed);

        assert_eq!(
            oa.final_simulation().pool.len(),
            ob.final_simulation().pool.len()
        );
        assert_eq!(
            oa.final_simulation()
                .assignments
                .get(Segment::V)
                .copied()
                .unwrap()
                .allele_id,
            ob.final_simulation()
                .assignments
                .get(Segment::V)
                .copied()
                .unwrap()
                .allele_id
        );
        // Trace contains exactly one entry — only the sampling pass records.
        assert_eq!(oa.trace.len(), 1);
        assert_eq!(oa.trace.choices()[0].address, "sample_allele.v");
    }

    /// A representative mixed plan: alternating EchoPass (no RNG)
    /// and SampleBasePass (RNG-consuming). Used by the cross-cutting
    /// replay-determinism test below.
    fn mixed_plan() -> PassPlan {
        use crate::ir::flag;
        let mut plan = PassPlan::new();
        for i in 0..8 {
            plan.push(Box::new(EchoPass::new(b'A', i as u16, Segment::V)));
            plan.push(Box::new(SampleBasePass::new(
                format!("np.np1.bases[{}]", i),
                Box::new(UniformBase),
                Segment::Np1,
                flag::N_NUC,
            )));
        }
        plan
    }

    #[test]
    fn replay_determinism_holds_under_repeated_runs() {
        // Run the same plan with the same seed five times. Every
        // trace and every final IR must be identical to the first.
        let baseline = PassRuntime::execute(&mixed_plan(), Simulation::new(), 0xfeed_face);
        for _ in 0..5 {
            let other = PassRuntime::execute(&mixed_plan(), Simulation::new(), 0xfeed_face);
            assert_eq!(other.trace.len(), baseline.trace.len());
            for (a, b) in baseline
                .trace
                .choices()
                .iter()
                .zip(other.trace.choices().iter())
            {
                assert_eq!(a.address, b.address);
                assert_eq!(a.value, b.value);
            }
            assert_eq!(
                other.final_simulation().pool.len(),
                baseline.final_simulation().pool.len()
            );
        }
    }

    // ──────────────────────────────────────────────────────────────
    // Per-pass `EventRecord.simulation_events` outcome contract
    //
    // Every committed pass now carries two complementary surfaces:
    //   - `trace_span` → choices the pass consumed/emitted
    //   - `simulation_events` → state consequences of those choices
    //
    // These tests pin which events each pass produces, indexed off
    // the committed `outcome.events` ledger — i.e. observability
    // from the *outcome* surface, not the test-only PassContext
    // hook.
    // ──────────────────────────────────────────────────────────────

    /// Find exactly-one event matching `pat` inside `events`. Fails
    /// the test if zero or more-than-one match. Returns the match.
    fn find_one<F>(
        events: &[crate::ir::SimulationEvent],
        pat: F,
    ) -> &crate::ir::SimulationEvent
    where
        F: Fn(&crate::ir::SimulationEvent) -> bool,
    {
        let matches: Vec<&crate::ir::SimulationEvent> =
            events.iter().filter(|e| pat(e)).collect();
        assert_eq!(
            matches.len(),
            1,
            "expected exactly one matching event, got {} from {:?}",
            matches.len(),
            events
        );
        matches[0]
    }

    #[test]
    fn outcome_event_records_carry_assignment_trim_and_region_events() {
        // Plan: sample V → trim V 5' → assemble V → generate NP1 →
        // assemble J. Each committed pass's `EventRecord` should
        // carry the consequence event it actually fired.
        let cfg = make_vj_refdata_for_filter();
        let trim_dist = || -> Box<dyn Distribution<Output = i64>> {
            Box::new(EmpiricalLengthDist::from_pairs(vec![(0, 1.0)]))
        };
        let np_dist = || -> Box<dyn Distribution<Output = i64>> {
            Box::new(EmpiricalLengthDist::from_pairs(vec![(3, 1.0)]))
        };

        let mut plan = PassPlan::new();
        plan.push(Box::new(SampleAllelePass::new(
            Segment::V,
            Box::new(AllelePoolDist::uniform(&cfg.v_pool)),
        )));
        plan.push(Box::new(SampleAllelePass::new(
            Segment::J,
            Box::new(AllelePoolDist::uniform(&cfg.j_pool)),
        )));
        plan.push(Box::new(TrimPass::new(
            Segment::V,
            TrimEnd::Five,
            trim_dist(),
        )));
        plan.push(Box::new(AssembleSegmentPass::new(Segment::V)));
        plan.push(Box::new(GenerateNPPass::new(
            Segment::Np1,
            np_dist(),
            Box::new(UniformBase),
        )));
        plan.push(Box::new(AssembleSegmentPass::new(Segment::J)));

        let outcome = PassRuntime::execute_with_refdata(&plan, Simulation::new(), 42, &cfg);

        // Sanity: pass count matches plan length.
        assert_eq!(outcome.events.len(), plan.len());

        // Pass 0 — SampleAllele V → one AssignmentChanged for V.
        let ev0 = &outcome.events[0];
        assert_eq!(ev0.pass_name, "sample_allele.v");
        let assign = find_one(&ev0.simulation_events, |e| {
            matches!(
                e,
                crate::ir::SimulationEvent::AssignmentChanged {
                    segment: Segment::V,
                    ..
                }
            )
        });
        match assign {
            crate::ir::SimulationEvent::AssignmentChanged { segment, old, .. } => {
                assert_eq!(*segment, Segment::V);
                assert!(old.is_none(), "first-install AssignmentChanged has old=None");
            }
            _ => unreachable!(),
        }

        // Pass 2 — TrimPass V 5' → one TrimChanged for (V, Five).
        let ev2 = &outcome.events[2];
        assert_eq!(ev2.pass_name, "trim.v_5");
        let trim_event = find_one(&ev2.simulation_events, |e| {
            matches!(
                e,
                crate::ir::SimulationEvent::TrimChanged {
                    segment: Segment::V,
                    end: TrimEnd::Five,
                    ..
                }
            )
        });
        match trim_event {
            crate::ir::SimulationEvent::TrimChanged { new, .. } => {
                assert_eq!(*new, 0, "trim distribution fixed at 0");
            }
            _ => unreachable!(),
        }

        // Pass 3 — AssembleSegment V → one RegionAdded for V.
        let ev3 = &outcome.events[3];
        assert_eq!(ev3.pass_name, "assemble.v");
        let v_region = find_one(&ev3.simulation_events, |e| {
            matches!(
                e,
                crate::ir::SimulationEvent::RegionAdded { region } if region.segment == Segment::V
            )
        });
        match v_region {
            crate::ir::SimulationEvent::RegionAdded { region } => {
                // V allele is 9 bases long with trim_5=0 → region [0, 9).
                assert_eq!(region.start, NucHandle::new(0));
                assert_eq!(region.end, NucHandle::new(9));
            }
            _ => unreachable!(),
        }

        // Pass 4 — GenerateNP Np1 → one RegionAdded for Np1.
        let ev4 = &outcome.events[4];
        assert_eq!(ev4.pass_name, "generate_np.np1");
        let np_region = find_one(&ev4.simulation_events, |e| {
            matches!(
                e,
                crate::ir::SimulationEvent::RegionAdded { region } if region.segment == Segment::Np1
            )
        });
        match np_region {
            crate::ir::SimulationEvent::RegionAdded { region } => {
                // NP1 starts where V ended (handle 9), length 3.
                assert_eq!(region.start, NucHandle::new(9));
                assert_eq!(region.end, NucHandle::new(12));
            }
            _ => unreachable!(),
        }

        // Pass 5 — AssembleSegment J → one RegionAdded for J.
        let ev5 = &outcome.events[5];
        assert_eq!(ev5.pass_name, "assemble.j");
        let _ = find_one(&ev5.simulation_events, |e| {
            matches!(
                e,
                crate::ir::SimulationEvent::RegionAdded { region } if region.segment == Segment::J
            )
        });
    }

    #[test]
    fn outcome_event_records_carry_base_changed_and_mutation_count_for_mutation_pass() {
        use crate::passes::UniformMutationPass;
        // Plan: sample V → assemble V → uniform-mutation 1.
        // The mutation pass's EventRecord must carry exactly one
        // `BaseChanged` (the single applied substitution) and
        // exactly one `MutationCountChanged` (the commit-time count
        // bump).
        let cfg = make_vj_refdata_for_filter();
        let mutation_count_dist = || -> Box<dyn Distribution<Output = i64>> {
            Box::new(EmpiricalLengthDist::from_pairs(vec![(1, 1.0)]))
        };

        let mut plan = PassPlan::new();
        plan.push(Box::new(SampleAllelePass::new(
            Segment::V,
            Box::new(AllelePoolDist::uniform(&cfg.v_pool)),
        )));
        plan.push(Box::new(AssembleSegmentPass::new(Segment::V)));
        plan.push(Box::new(UniformMutationPass::new(
            mutation_count_dist(),
            Box::new(UniformBase),
        )));

        let outcome = PassRuntime::execute_with_refdata(&plan, Simulation::new(), 42, &cfg);

        // The mutation pass is the third committed pass.
        let mut_event = outcome
            .events
            .iter()
            .find(|e| e.pass_name == "mutate.uniform")
            .expect("uniform mutation pass committed an event record");

        let _ = find_one(&mut_event.simulation_events, |e| {
            matches!(e, crate::ir::SimulationEvent::BaseChanged { .. })
        });
        let mc = find_one(&mut_event.simulation_events, |e| {
            matches!(e, crate::ir::SimulationEvent::MutationCountChanged { .. })
        });
        match mc {
            crate::ir::SimulationEvent::MutationCountChanged {
                old,
                new,
                delta,
            } => {
                assert_eq!(*old, 0);
                assert_eq!(*new, 1);
                assert_eq!(*delta, 1);
            }
            _ => unreachable!(),
        }

        // Final mutation_count carried on the sealed sim matches the
        // event's `new` field — pin the outcome-ledger / sim
        // consistency contract.
        assert_eq!(outcome.final_simulation().mutation_count, 1);
    }

    // ──────────────────────────────────────────────────────────────
    // Event-coverage audit + invariant tests
    //
    // Three coordinated tests that lock down the consequence-event
    // ledger before any derived-state lifecycle refactor:
    //
    //   1. Cross-pass coverage — every meaningful pass emits the
    //      consequence variants it owns.
    //   2. Event/state consistency — outcome counts reconstructed
    //      from `simulation_events` match the final `Simulation`.
    //   3. Trace/event alignment — `TraceSpan`s are contiguous,
    //      cover the full trace, and replay equality holds.
    // ──────────────────────────────────────────────────────────────

    /// Find every event matching `pat` in a slice. Used in the
    /// audit tests where we care about "at least one of X" and
    /// counts, not exact-one-of-X.
    fn count_matching<F>(events: &[crate::ir::SimulationEvent], pat: F) -> usize
    where
        F: Fn(&crate::ir::SimulationEvent) -> bool,
    {
        events.iter().filter(|e| pat(e)).count()
    }

    /// Build the heavy-stack VJ plan used by the coverage audit.
    /// Every category the user listed (sample/trim/assemble/np/
    /// mutation/pcr/quality/ncorrupt/indel/rev-comp) is represented
    /// by exactly one pass. Distributions are pinned to fixed
    /// counts so the audit is deterministic across RNG paths.
    fn build_coverage_audit_plan(
        cfg: &crate::refdata::RefDataConfig,
    ) -> PassPlan {
        use crate::passes::{
            IndelPass, NCorruptionPass, PCRErrorPass, QualityErrorPass, RevCompPass,
            UniformMutationPass,
        };
        let count_one = || -> Box<dyn Distribution<Output = i64>> {
            Box::new(EmpiricalLengthDist::from_pairs(vec![(1, 1.0)]))
        };
        let count_zero = || -> Box<dyn Distribution<Output = i64>> {
            Box::new(EmpiricalLengthDist::from_pairs(vec![(0, 1.0)]))
        };
        let np_three = || -> Box<dyn Distribution<Output = i64>> {
            Box::new(EmpiricalLengthDist::from_pairs(vec![(3, 1.0)]))
        };

        let mut plan = PassPlan::new();
        plan.push(Box::new(SampleAllelePass::new(
            Segment::V,
            Box::new(AllelePoolDist::uniform(&cfg.v_pool)),
        )));
        plan.push(Box::new(SampleAllelePass::new(
            Segment::J,
            Box::new(AllelePoolDist::uniform(&cfg.j_pool)),
        )));
        plan.push(Box::new(TrimPass::new(
            Segment::V,
            TrimEnd::Five,
            count_zero(),
        )));
        plan.push(Box::new(AssembleSegmentPass::new(Segment::V)));
        plan.push(Box::new(GenerateNPPass::new(
            Segment::Np1,
            np_three(),
            Box::new(UniformBase),
        )));
        plan.push(Box::new(AssembleSegmentPass::new(Segment::J)));
        plan.push(Box::new(UniformMutationPass::new(
            count_one(),
            Box::new(UniformBase),
        )));
        plan.push(Box::new(PCRErrorPass::new(count_one(), Box::new(UniformBase))));
        plan.push(Box::new(QualityErrorPass::new(
            count_one(),
            Box::new(UniformBase),
        )));
        plan.push(Box::new(NCorruptionPass::new(count_one())));
        plan.push(Box::new(IndelPass::new(
            count_one(),
            0.5,
            Box::new(UniformBase),
        )));
        plan.push(Box::new(RevCompPass::new(1.0)));
        plan
    }

    #[test]
    fn cross_pass_event_coverage_emits_expected_consequence_types() {
        // Drive a heavy-stack VJ plan covering every event-emitting
        // pass category. Assert each pass's `EventRecord.simulation_events`
        // carries at least one event of the expected variant — this
        // is the load-bearing "no holes in the ledger" check before
        // we lean on it for lifecycle decisions.
        use crate::ir::SimulationEvent;

        let cfg = make_vj_refdata_for_filter();
        let plan = build_coverage_audit_plan(&cfg);
        let outcome = PassRuntime::execute_with_refdata(&plan, Simulation::new(), 42, &cfg);

        // Convenience: locate a committed pass's event record by
        // its `pass_name`. Panics if missing because the plan is
        // hard-coded — any missing pass is a regression.
        let find = |name: &str| {
            outcome
                .events
                .iter()
                .find(|e| e.pass_name == name)
                .unwrap_or_else(|| panic!("missing event record for {name}"))
        };

        // Sample allele: AssignmentChanged for V and J.
        assert!(
            count_matching(&find("sample_allele.v").simulation_events, |e| matches!(
                e,
                SimulationEvent::AssignmentChanged { segment: Segment::V, .. }
            )) >= 1
        );
        assert!(
            count_matching(&find("sample_allele.j").simulation_events, |e| matches!(
                e,
                SimulationEvent::AssignmentChanged { segment: Segment::J, .. }
            )) >= 1
        );

        // Trim: TrimChanged.
        assert!(
            count_matching(&find("trim.v_5").simulation_events, |e| matches!(
                e,
                SimulationEvent::TrimChanged { segment: Segment::V, .. }
            )) >= 1
        );

        // Assemble: many BasePushed + one RegionAdded per segment.
        for (name, seg) in [("assemble.v", Segment::V), ("assemble.j", Segment::J)] {
            let ev = find(name);
            assert!(
                count_matching(&ev.simulation_events, |e| matches!(e, SimulationEvent::BasePushed { .. })) > 0,
                "{name} should emit BasePushed events"
            );
            assert!(
                count_matching(&ev.simulation_events, |e| matches!(
                    e,
                    SimulationEvent::RegionAdded { region } if region.segment == seg
                )) == 1,
                "{name} should emit exactly one RegionAdded for its segment"
            );
        }

        // NP: BasePushed (length 3) + one RegionAdded for Np1.
        let np = find("generate_np.np1");
        assert!(
            count_matching(&np.simulation_events, |e| matches!(e, SimulationEvent::BasePushed { .. })) >= 1
        );
        assert!(
            count_matching(&np.simulation_events, |e| matches!(
                e,
                SimulationEvent::RegionAdded { region } if region.segment == Segment::Np1
            )) == 1
        );

        // Uniform mutation: BaseChanged + MutationCountChanged.
        let mu = find("mutate.uniform");
        assert!(
            count_matching(&mu.simulation_events, |e| matches!(e, SimulationEvent::BaseChanged { .. })) >= 1
        );
        assert!(
            count_matching(&mu.simulation_events, |e| matches!(
                e,
                SimulationEvent::MutationCountChanged { .. }
            )) == 1
        );

        // PCR / quality / ncorrupt: BaseChanged (no mutation-count
        // bump — these are sequencing artifacts, not biological
        // mutations, and don't tag `add_to_mutation_count`).
        for name in ["corrupt.pcr", "corrupt.quality", "corrupt.ns"] {
            let ev = find(name);
            assert!(
                count_matching(&ev.simulation_events, |e| matches!(e, SimulationEvent::BaseChanged { .. })) >= 1,
                "{name} should emit at least one BaseChanged"
            );
            assert!(
                count_matching(&ev.simulation_events, |e| matches!(
                    e,
                    SimulationEvent::MutationCountChanged { .. }
                )) == 0,
                "{name} is a sequencing-artifact pass and must not bump n_mutations"
            );
        }

        // Indel: at least one of IndelInserted / IndelDeleted.
        let indel = find("corrupt.indel");
        let indels = count_matching(&indel.simulation_events, |e| {
            matches!(
                e,
                SimulationEvent::IndelInserted { .. } | SimulationEvent::IndelDeleted { .. }
            )
        });
        assert!(indels >= 1, "corrupt.indel should emit at least one indel event");

        // Rev-comp: one ReverseComplementFlagRecorded with applied=true
        // (apply_prob = 1.0 in the audit plan).
        let rc = find("corrupt.rev_comp");
        assert_eq!(
            count_matching(&rc.simulation_events, |e| matches!(
                e,
                SimulationEvent::ReverseComplementFlagRecorded { applied: true }
            )),
            1
        );
    }

    #[test]
    fn outcome_events_are_consistent_with_final_simulation_state() {
        // Reconstruct three coarse counts purely from the
        // `simulation_events` stream and assert they match the
        // sealed final `Simulation`. This is the load-bearing
        // self-consistency proof: a future derived-state observer
        // can trust the ledger as the source of truth.
        use crate::ir::SimulationEvent;

        let cfg = make_vj_refdata_for_filter();
        let count_one = || -> Box<dyn Distribution<Output = i64>> {
            Box::new(EmpiricalLengthDist::from_pairs(vec![(1, 1.0)]))
        };
        let np_three = || -> Box<dyn Distribution<Output = i64>> {
            Box::new(EmpiricalLengthDist::from_pairs(vec![(3, 1.0)]))
        };

        let mut plan = PassPlan::new();
        plan.push(Box::new(SampleAllelePass::new(
            Segment::V,
            Box::new(AllelePoolDist::uniform(&cfg.v_pool)),
        )));
        plan.push(Box::new(SampleAllelePass::new(
            Segment::J,
            Box::new(AllelePoolDist::uniform(&cfg.j_pool)),
        )));
        plan.push(Box::new(AssembleSegmentPass::new(Segment::V)));
        plan.push(Box::new(GenerateNPPass::new(
            Segment::Np1,
            np_three(),
            Box::new(UniformBase),
        )));
        plan.push(Box::new(AssembleSegmentPass::new(Segment::J)));
        plan.push(Box::new(crate::passes::UniformMutationPass::new(
            count_one(),
            Box::new(UniformBase),
        )));

        let outcome = PassRuntime::execute_with_refdata(&plan, Simulation::new(), 0, &cfg);
        let final_sim = outcome.final_simulation();

        // ── Invariant 1: regions ────────────────────────────────
        // Indel-free plan → final region count equals the total
        // `RegionAdded` events across all passes.
        let region_adds: usize = outcome
            .events
            .iter()
            .flat_map(|e| &e.simulation_events)
            .filter(|e| matches!(e, SimulationEvent::RegionAdded { .. }))
            .count();
        assert_eq!(
            region_adds,
            final_sim.sequence.regions.len(),
            "sum of RegionAdded events must equal final region count for indel-free plan"
        );

        // ── Invariant 2: assignments ────────────────────────────
        // The last `AssignmentChanged` per segment must equal the
        // segment's final `AlleleInstance` in `final_sim.assignments`.
        use crate::assignment::AlleleInstance;
        let mut last_by_segment: std::collections::HashMap<Segment, AlleleInstance> =
            std::collections::HashMap::new();
        for ev in outcome.events.iter().flat_map(|e| &e.simulation_events) {
            if let SimulationEvent::AssignmentChanged { segment, new, .. } = ev {
                last_by_segment.insert(*segment, *new);
            }
        }
        for (seg, last_assignment) in &last_by_segment {
            assert_eq!(
                final_sim.assignments.get(*seg).copied(),
                Some(*last_assignment),
                "last AssignmentChanged for {:?} must equal final assignment",
                seg
            );
        }
        // V and J should both have been assigned.
        assert!(last_by_segment.contains_key(&Segment::V));
        assert!(last_by_segment.contains_key(&Segment::J));

        // ── Invariant 3: mutation_count ─────────────────────────
        // The last `MutationCountChanged.new` across the run must
        // equal `final_sim.mutation_count`. (Plan has exactly one
        // count-bumping pass; this still holds with multiple.)
        let final_mc = outcome
            .events
            .iter()
            .flat_map(|e| &e.simulation_events)
            .filter_map(|e| match e {
                SimulationEvent::MutationCountChanged { new, .. } => Some(*new),
                _ => None,
            })
            .last()
            .expect("plan has one mutation pass");
        assert_eq!(final_mc, final_sim.mutation_count);
    }

    #[test]
    fn trace_spans_partition_outcome_trace_and_replay_preserves_trace() {
        // Two-part alignment proof:
        //
        //   (a) Every committed pass yields one `EventRecord` whose
        //       `trace_span` is contiguous with its neighbours and
        //       whose union covers the entire outcome trace. No
        //       gaps, no overlaps.
        //
        //   (b) Replaying the captured trace through `CompiledSimulator`
        //       produces an outcome whose trace matches the original
        //       record-for-record — proving event capture is
        //       transparent w.r.t. the trace, the only persisted
        //       artifact.
        use crate::compiled::ExecutionPolicy;

        let cfg = make_vj_refdata_for_filter();
        let np_three = || -> Box<dyn Distribution<Output = i64>> {
            Box::new(EmpiricalLengthDist::from_pairs(vec![(3, 1.0)]))
        };
        let mut plan = PassPlan::new();
        plan.push(Box::new(SampleAllelePass::new(
            Segment::V,
            Box::new(AllelePoolDist::uniform(&cfg.v_pool)),
        )));
        plan.push(Box::new(SampleAllelePass::new(
            Segment::J,
            Box::new(AllelePoolDist::uniform(&cfg.j_pool)),
        )));
        plan.push(Box::new(AssembleSegmentPass::new(Segment::V)));
        plan.push(Box::new(GenerateNPPass::new(
            Segment::Np1,
            np_three(),
            Box::new(UniformBase),
        )));
        plan.push(Box::new(AssembleSegmentPass::new(Segment::J)));

        // ── (a) Trace-span partition ─────────────────────────────
        let outcome = PassRuntime::execute_with_refdata(&plan, Simulation::new(), 7, &cfg);
        assert_eq!(
            outcome.events.len(),
            plan.len(),
            "every committed pass must contribute exactly one EventRecord"
        );
        assert_eq!(outcome.events[0].trace_span.start, 0);
        for i in 1..outcome.events.len() {
            assert_eq!(
                outcome.events[i - 1].trace_span.end,
                outcome.events[i].trace_span.start,
                "trace spans must be contiguous (gap or overlap at boundary {})",
                i
            );
        }
        assert_eq!(
            outcome
                .events
                .last()
                .map(|e| e.trace_span.end)
                .unwrap_or(0),
            outcome.trace.len(),
            "the union of all trace spans must cover the full outcome trace"
        );

        // ── (b) Replay equality ──────────────────────────────────
        // Compile the same plan + refdata into an owned simulator,
        // replay the captured trace, and assert the replayed
        // trace matches the original record-for-record. Event
        // capture happens on both runs (the compiled path always
        // supplies a sink) but neither side persists events — so a
        // matching trace is sufficient proof that capture is
        // transparent.
        //
        // `OwnedCompiledSimulator::compile` takes ownership of the
        // plan + refdata, so build a fresh equivalent plan for the
        // replay leg. (The borrowing `CompiledSimulator` doesn't
        // expose `replay_from_trace_records`.)
        let mut replay_plan = PassPlan::new();
        replay_plan.push(Box::new(SampleAllelePass::new(
            Segment::V,
            Box::new(AllelePoolDist::uniform(&cfg.v_pool)),
        )));
        replay_plan.push(Box::new(SampleAllelePass::new(
            Segment::J,
            Box::new(AllelePoolDist::uniform(&cfg.j_pool)),
        )));
        replay_plan.push(Box::new(AssembleSegmentPass::new(Segment::V)));
        replay_plan.push(Box::new(GenerateNPPass::new(
            Segment::Np1,
            np_three(),
            Box::new(UniformBase),
        )));
        replay_plan.push(Box::new(AssembleSegmentPass::new(Segment::J)));
        let compiled = crate::compiled::OwnedCompiledSimulator::compile_with_options(
            replay_plan,
            Some(cfg.clone()),
            None,
            ExecutionPolicy::Permissive,
            crate::compiled::CompileOptions::skip_refdata_validation(),
        )
        .expect("audit plan should compile");
        let original_records: Vec<_> = outcome.trace.choices().to_vec();
        let replayed = compiled
            .replay_from_trace_records(
                &original_records,
                /* seed unused for migrated sites */ 0,
                ExecutionPolicy::Strict,
            )
            .expect("replay should succeed against the same plan");
        let replayed_records: Vec<_> = replayed.trace.choices().to_vec();
        assert_eq!(
            original_records.len(),
            replayed_records.len(),
            "replay must produce a trace of the same length"
        );
        for (i, (orig, repl)) in original_records.iter().zip(replayed_records.iter()).enumerate() {
            assert_eq!(orig.address, repl.address, "address divergence at record {}", i);
            assert_eq!(orig.value, repl.value, "value divergence at record {} (addr={})", i, orig.address);
        }
    }

    // ──────────────────────────────────────────────────────────────
    // Compile-effect / event-emission policy conformance
    //
    // Built-in mutating passes that declare a `PassCompileEffect`
    // are expected to emit at least one corresponding
    // `SimulationEvent` when they actually mutate state. This test
    // exercises the heavy-stack VJ plan and asserts every
    // `EventRecord` passes the policy check in
    // [`crate::event::check_event_emission_consistency`].
    //
    // It is the regression net for the "new pass mutates Simulation
    // directly via sim.with_* and forgets to emit events" anti-
    // pattern. The runtime live-call refresh trusts events; a
    // silent zero-event mutating pass would skip the refresh
    // (this is also what the divergence tests in
    // `compiled::tests::live_call_edits` demonstrate by construction).
    // ──────────────────────────────────────────────────────────────

    #[test]
    fn builtin_passes_emit_events_consistent_with_declared_compile_effects() {
        let cfg = make_vj_refdata_for_filter();
        let plan = build_coverage_audit_plan(&cfg);
        let outcome = PassRuntime::execute_with_refdata(&plan, Simulation::new(), 42, &cfg);

        // Apply the policy to every committed pass. The plan's
        // distributions are pinned to deterministic non-zero
        // counts so the EditBases / StructuralIndel zero-
        // exemption never fires here — meaningful coverage.
        let mut violations: Vec<String> = Vec::new();
        for record in &outcome.events {
            if let Err(err) = crate::event::check_event_emission_consistency(record) {
                violations.push(format!(
                    "{}: declared {:?} but {}",
                    err.pass_name, err.compile_effect, err.reason
                ));
            }
        }
        assert!(
            violations.is_empty(),
            "built-in passes must emit events matching their compile-effect declarations.\n\
             Violations:\n  {}\n\n\
             Likely cause: the pass mutates `Simulation` via `sim.with_*` directly \
             instead of routing through `SimulationBuilder` — runtime derived-state \
             refresh trusts events, so a silent mutator gets no refresh.",
            violations.join("\n  ")
        );
    }
}
