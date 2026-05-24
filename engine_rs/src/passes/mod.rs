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
pub(crate) mod constrained;
pub mod corrupt;
pub mod echo;
pub mod generate_np;
pub mod mutate;
pub mod sample_allele;
pub mod sample_base;
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
pub use mutate::{S5FMutationPass, UniformMutationPass};
pub use sample_allele::SampleAllelePass;
pub use sample_base::SampleBasePass;
pub use trim::TrimPass;

// ──────────────────────────────────────────────────────────────────
// Cross-cutting integration tests (D.6 productive bundle)
// ──────────────────────────────────────────────────────────────────

#[cfg(test)]
mod tests {
    use super::*;
    use crate::contract::{
        productive, Contract, ContractSet, ContractViolation, ProductiveJunctionFrame,
    };
    use crate::dist::{
        AllelePoolDist, Distribution, EmpiricalLengthDist, FilteredSampleError, UniformBase,
    };
    use crate::ir::{Segment, Simulation};
    use crate::pass::{PassPlan, PassRuntime};
    use crate::passes::sample_allele::test_support::make_test_pool;
    use crate::trace::ChoiceValue;

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
        });
        let _ = cfg.j_pool.push(crate::refdata::Allele {
            name: "j_test*01".into(),
            gene: "j_test".into(),
            seq: b"TTTAAA".to_vec(),
            segment: Segment::J,
            anchor: Some(0),
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
    fn productive_admits_falls_back_when_filter_empty() {
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
        assert!(matches!(np1_len, 1 | 2 | 4 | 5));
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

        fn admits(
            &self,
            _sim: &Simulation,
            _refdata: Option<&crate::refdata::RefDataConfig>,
            address: &str,
            _candidate: &ChoiceValue,
        ) -> Result<(), ContractViolation> {
            if address.starts_with("np.np1.bases[") {
                return Err(ContractViolation::new(self.name(), "rejected by test"));
            }
            Ok(())
        }
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
            oa.final_simulation().assignments.get(Segment::V).copied().unwrap().allele_id,
            ob.final_simulation().assignments.get(Segment::V).copied().unwrap().allele_id
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
}
