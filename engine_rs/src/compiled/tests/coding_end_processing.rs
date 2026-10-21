use super::*;
use crate::address::PEnd;
use crate::end_processing::EndOutcome;
use crate::passes::{EmitCodingEndPPass, ProcessCodingEndsPass};

#[test]
fn edits_must_follow_the_original_join_even_without_productive_contracts() {
    let (cfg, _) = fixture(&[6], -1, 0, true);
    let mut plan = PassPlan::new();
    plan.push(Box::new(SampleAllelePass::new(
        Segment::D,
        Box::new(AllelePoolDist::uniform(&cfg.d_pool)),
    )));
    plan.push(Box::new(
        ProcessCodingEndsPass::new(
            Segment::D,
            vec![EndOutcome {
                five: -1,
                three: 0,
                weight: 1.0,
            }],
        )
        .unwrap(),
    ));
    plan.push(Box::new(EmitCodingEndPPass::new(PEnd::D5)));
    plan.push(Box::new(crate::passes::UniformMutationPass::new(
        Box::new(EmpiricalLengthDist::from_pairs(vec![(0, 1.0)])),
        Box::new(UniformBase),
    )));
    plan.push(Box::new(AssembleSegmentPass::new(Segment::D)));
    plan.push(Box::new(EmitCodingEndPPass::new(PEnd::D3)));
    assert!(compile_test_fixture(&plan, Some(&cfg), None, ExecutionPolicy::Strict).is_err());
    let (cfg, mut plan) = fixture(&[6], -1, 0, true);
    plan.push(Box::new(crate::passes::UniformMutationPass::new(
        Box::new(EmpiricalLengthDist::from_pairs(vec![(0, 1.0)])),
        Box::new(UniformBase),
    )));
    assert!(compile_test_fixture(&plan, Some(&cfg), None, ExecutionPolicy::Strict).is_ok());
}

fn fixture(lengths: &[usize], five: i64, three: i64, emit: bool) -> (RefDataConfig, PassPlan) {
    let mut cfg = RefDataConfig::empty(ChainType::Vdj);
    for (index, length) in lengths.iter().enumerate() {
        let _ = cfg.d_pool.push(Allele {
            name: format!("D{index}*01"),
            gene: format!("D{index}"),
            seq: vec![b'C'; *length],
            segment: Segment::D,
            anchor: None,
            functional_status: None,
            subregions: vec![],
        });
    }
    let mut plan = PassPlan::new();
    plan.push(Box::new(SampleAllelePass::new(
        Segment::D,
        Box::new(AllelePoolDist::uniform(&cfg.d_pool)),
    )));
    plan.push(Box::new(
        ProcessCodingEndsPass::new(
            Segment::D,
            vec![EndOutcome {
                five,
                three,
                weight: 1.0,
            }],
        )
        .unwrap(),
    ));
    if emit {
        plan.push(Box::new(EmitCodingEndPPass::new(PEnd::D5)));
    }
    plan.push(Box::new(AssembleSegmentPass::new(Segment::D)));
    if emit {
        plan.push(Box::new(EmitCodingEndPPass::new(PEnd::D3)));
    }
    (cfg, plan)
}

#[test]
fn short_alleles_are_filtered_before_sampling_without_productivity_contracts() {
    let (cfg, plan) = fixture(&[1, 6], -3, 2, true);
    for policy in [ExecutionPolicy::Strict, ExecutionPolicy::Permissive] {
        let compiled = compile_test_fixture(&plan, Some(&cfg), None, policy).unwrap();
        for seed in 0..64 {
            let outcome = compiled.run_one(seed).unwrap();
            let sim = outcome.final_simulation();
            assert_eq!(
                sim.assignments.get(Segment::D).unwrap().allele_id.index(),
                1
            );
            assert_eq!(sim.pool.len(), 7);
        }
    }
}

#[test]
fn incomplete_emission_and_legacy_trim_conflicts_are_compile_errors() {
    let (cfg, plan) = fixture(&[6], 0, 0, false);
    assert!(compile_test_fixture(&plan, Some(&cfg), None, ExecutionPolicy::Strict).is_err());
    let (cfg, mut plan) = fixture(&[6], 0, 0, true);
    plan.push(Box::new(TrimPass::new(
        Segment::D,
        TrimEnd::Five,
        Box::new(EmpiricalLengthDist::from_pairs(vec![(0, 1.0)])),
    )));
    assert!(compile_test_fixture(&plan, Some(&cfg), None, ExecutionPolicy::Strict).is_err());
    let (mut cfg, mut plan) = fixture(&[6], 0, 0, true);
    let _ = cfg.j_pool.push(Allele {
        name: "J*01".into(),
        gene: "J".into(),
        seq: b"TTTGGG".to_vec(),
        segment: Segment::J,
        anchor: Some(0),
        functional_status: None,
        subregions: vec![],
    });
    plan.push(Box::new(SampleAllelePass::new(
        Segment::J,
        Box::new(AllelePoolDist::uniform(&cfg.j_pool)),
    )));
    // A legacy P pass on another assigned segment has no new-model ancestry.
    plan.push(Box::new(crate::passes::PAdditionPass::new(
        PEnd::J5,
        Box::new(EmpiricalLengthDist::from_pairs(vec![(1, 1.0)])),
    )));
    plan.push(Box::new(AssembleSegmentPass::new(Segment::J)));
    let error = match compile_test_fixture(&plan, Some(&cfg), None, ExecutionPolicy::Strict) {
        Err(error) => error,
        Ok(_) => panic!("mixed P models must be rejected"),
    };
    assert!(error.errors.iter().any(|error| matches!(&error.kind,
        CompileErrorKind::InvalidParameterSupport { reason, .. }
            if reason.contains("cannot be mixed with legacy P-addition")
    )));
}

#[test]
fn complete_d_deletion_is_allowed_and_crossing_deletions_fail_at_compile() {
    let (cfg, plan) = fixture(&[3], 1, 2, true);
    let compiled = compile_test_fixture(&plan, Some(&cfg), None, ExecutionPolicy::Strict).unwrap();
    assert_eq!(
        compiled.run_one(0).unwrap().final_simulation().pool.len(),
        0
    );
    let (cfg, plan) = fixture(&[3], 2, 2, true);
    assert!(compile_test_fixture(&plan, Some(&cfg), None, ExecutionPolicy::Permissive).is_err());
}

#[test]
fn d5_emission_before_inversion_is_rejected_without_contracts() {
    let (cfg, _) = fixture(&[6], -2, 0, true);
    let plan_for = |late_inversion: bool| {
        let mut plan = PassPlan::new();
        plan.push(Box::new(SampleAllelePass::new(
            Segment::D,
            Box::new(AllelePoolDist::uniform(&cfg.d_pool)),
        )));
        plan.push(Box::new(
            ProcessCodingEndsPass::new(
                Segment::D,
                vec![EndOutcome {
                    five: -2,
                    three: 0,
                    weight: 1.0,
                }],
            )
            .unwrap(),
        ));
        if !late_inversion {
            plan.push(Box::new(crate::passes::InvertDPass::new(1.0)));
        }
        plan.push(Box::new(EmitCodingEndPPass::new(PEnd::D5)));
        if late_inversion {
            plan.push(Box::new(crate::passes::InvertDPass::new(1.0)));
        }
        plan.push(Box::new(AssembleSegmentPass::new(Segment::D)));
        plan.push(Box::new(EmitCodingEndPPass::new(PEnd::D3)));
        plan
    };
    let valid = plan_for(false);
    let invalid = plan_for(true);
    for policy in [ExecutionPolicy::Strict, ExecutionPolicy::Permissive] {
        assert!(compile_test_fixture(&valid, Some(&cfg), None, policy).is_ok());
        let errors = match compile_test_fixture(&invalid, Some(&cfg), None, policy) {
            Ok(_) => panic!("D5 emission must see the committed D orientation"),
            Err(errors) => errors,
        };
        assert!(
            errors.errors.iter().any(|error| matches!(
                &error.kind,
                CompileErrorKind::InvalidParameterSupport { reason, .. }
                    if reason.contains("coding-order")
            )),
            "{errors:?}"
        );
    }
}

#[test]
fn v3_emission_after_n1_is_rejected_without_contracts() {
    let cfg = vj_refdata();
    let plan_for = |late_emission: bool| {
        let mut plan = PassPlan::new();
        for (segment, pool) in [(Segment::V, &cfg.v_pool), (Segment::J, &cfg.j_pool)] {
            plan.push(Box::new(SampleAllelePass::new(
                segment,
                Box::new(AllelePoolDist::uniform(pool)),
            )));
        }
        plan.push(Box::new(
            ProcessCodingEndsPass::new(
                Segment::V,
                vec![EndOutcome {
                    five: 0,
                    three: -2,
                    weight: 1.0,
                }],
            )
            .unwrap(),
        ));
        plan.push(Box::new(AssembleSegmentPass::new(Segment::V)));
        if !late_emission {
            plan.push(Box::new(EmitCodingEndPPass::new(PEnd::V3)));
        }
        plan.push(Box::new(GenerateNPPass::new(
            Segment::Np1,
            Box::new(EmpiricalLengthDist::from_pairs(vec![(2, 1.0)])),
            Box::new(UniformBase),
        )));
        if late_emission {
            plan.push(Box::new(EmitCodingEndPPass::new(PEnd::V3)));
        }
        plan.push(Box::new(AssembleSegmentPass::new(Segment::J)));
        plan
    };
    let valid = plan_for(false);
    let invalid = plan_for(true);
    for policy in [ExecutionPolicy::Strict, ExecutionPolicy::Permissive] {
        assert!(compile_test_fixture(&valid, Some(&cfg), None, policy).is_ok());
        let errors = match compile_test_fixture(&invalid, Some(&cfg), None, policy) {
            Ok(_) => panic!("V3 P must be emitted before N1"),
            Err(errors) => errors,
        };
        assert!(
            errors.errors.iter().any(|error| matches!(
                &error.kind,
                CompileErrorKind::InvalidParameterSupport { reason, .. }
                    if reason.contains("coding-order")
            )),
            "{errors:?}"
        );
    }
}
