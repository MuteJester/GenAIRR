use super::*;
use crate::ir::Segment;
use crate::refdata::{Allele, ChainType};

fn refs() -> RefDataConfig {
    let mut r = RefDataConfig::empty(ChainType::Vdj);
    for (segment, pool, name, seq, anchor) in [
        (
            Segment::V,
            &mut r.v_pool,
            "V*01",
            b"AAATGT".as_slice(),
            Some(3),
        ),
        (Segment::D, &mut r.d_pool, "D*01", b"ACT".as_slice(), None),
        (
            Segment::J,
            &mut r.j_pool,
            "J*01",
            b"TTTAAA".as_slice(),
            Some(0),
        ),
    ] {
        let _ = pool.push(Allele {
            name: name.into(),
            gene: name.into(),
            seq: seq.to_vec(),
            segment,
            anchor,
            functional_status: None,
            subregions: vec![],
        });
    }
    r
}
fn integer() -> Box<dyn Distribution<Output = i64>> {
    Box::new(EmpiricalLengthDist::from_pairs([
        (0, 0.1),
        (2, 0.2),
        (2, 0.7),
    ]))
}

#[test]
fn all_post_read_native_descriptors_preserve_parameters_and_rng() {
    let r = refs();
    let k = S5FKernel::new(vec![0.5; 1024], vec![0.25; 4096]);
    let mut plan = PassPlan::new();
    plan.push(Box::new(UniformMutationPass::new(
        integer(),
        Box::new(UniformBase),
    )));
    plan.push(Box::new(UniformMutationPass::new_rate(
        0.05,
        Box::new(CategoricalBase::from_weights([0.0, 1.0, 2.0, 3.0])),
    )));
    plan.push(Box::new(
        S5FMutationPass::new(k.clone(), integer()).with_kernel_name("same-name"),
    ));
    plan.push(Box::new(S5FMutationPass::new_rate(k, 0.025)));
    plan.push(Box::new(PCRErrorPass::new(
        integer(),
        Box::new(UniformBase),
    )));
    plan.push(Box::new(PCRErrorPass::new_rate(0.1, Box::new(UniformBase))));
    plan.push(Box::new(QualityErrorPass::new(
        integer(),
        Box::new(UniformBase),
    )));
    plan.push(Box::new(QualityErrorPass::new_rate(
        0.2,
        Box::new(UniformBase),
    )));
    plan.push(Box::new(NCorruptionPass::new(integer())));
    plan.push(Box::new(
        IndelPass::new(integer(), 0.75, Box::new(UniformBase)).with_sampling_options(
            IndelSamplingOptions {
                rejection_attempts: 7,
                exact_node_budget: 12345,
                legacy_replay_sentinels: false,
            },
        ),
    ));
    plan.push(Box::new(EndLossPass::new(
        LossEnd::Five,
        Box::new(UniformInt::new(0, 7)),
    )));
    plan.push(Box::new(EndLossPass::new(LossEnd::Three, integer())));
    plan.push(Box::new(RevCompPass::new(0.125)));
    let spec = ReplayPlanSpecV1::capture(&plan, &r, None, false, false).unwrap();
    let wire = serde_json::to_string(&spec).unwrap();
    let decoded: ReplayPlanSpecV1 = serde_json::from_str(&wire).unwrap();
    let (restored, rr, c) = decoded.reconstruct().unwrap();
    assert!(c.is_none());
    assert_eq!(
        spec,
        ReplayPlanSpecV1::capture(&restored, &rr, None, false, false).unwrap()
    );
    for (a, b) in plan.passes().iter().zip(restored.passes()) {
        assert_eq!(a.name(), b.name());
        assert_eq!(a.sampling_law_signature(), b.sampling_law_signature());
    }
}

#[test]
fn exact_cdf_round_trips_without_reconstructing_from_subtracted_weights() {
    let d = EmpiricalLengthDist::from_pairs([(0, 0.1), (7, 0.2), (2, 0.7)]);
    let spec = d.family_replay_spec().unwrap();
    let rebuilt = spec.integer().unwrap();
    assert_eq!(spec, rebuilt.family_replay_spec().unwrap());
    for seed in 0..128 {
        let mut a = crate::rng::Rng::new(seed);
        let mut b = crate::rng::Rng::new(seed);
        for _ in 0..64 {
            assert_eq!(d.sample(&mut a), rebuilt.sample(&mut b));
        }
        assert_eq!(a.words_consumed(), b.words_consumed());
    }
}

#[test]
fn aliasing_legacy_identity_messages_have_distinct_typed_commitments() {
    let mut a = refs();
    let mut b = refs();
    a.identity.species = Some("a\nidentity.locus=b".into());
    a.identity.locus = Some("c".into());
    b.identity.species = Some("a".into());
    b.identity.locus = Some("b\nidentity.locus=c".into());
    assert_eq!(
        crate::trace_file::refdata_content_hash(&a),
        crate::trace_file::refdata_content_hash(&b)
    );
    let sa = ReferenceSpecV1::capture(&a).unwrap();
    let sb = ReferenceSpecV1::capture(&b).unwrap();
    assert_ne!(sa.digest().unwrap(), sb.digest().unwrap());
    assert_eq!(
        sa,
        ReferenceSpecV1::capture(&sa.reconstruct().unwrap()).unwrap()
    );
}

#[test]
fn malformed_numeric_and_distribution_parameters_return_errors_without_panics() {
    let bad = DistributionSpec::EmpiricalInt {
        values: vec![1],
        cumulative: vec![],
        total: F64Bits::new(1.0).unwrap(),
    };
    assert!(bad.integer().is_err());
    assert!(DistributionSpec::UniformInt {
        min: i64::MIN,
        max: i64::MAX
    }
    .integer()
    .is_err());
    let mut kernel =
        KernelSpecV1::capture(&S5FKernel::new(vec![1.0; 1024], vec![0.25; 4096])).unwrap();
    kernel.mutability[0] = F64Bits::new(-1.0).unwrap();
    assert!(kernel.reconstruct().is_err());
}

#[test]
fn absent_and_present_empty_contract_sets_stay_distinct() {
    let p = PassPlan::new();
    let r = refs();
    let c = crate::contract::ContractSet::new();
    let a = ReplayPlanSpecV1::capture(&p, &r, None, false, false).unwrap();
    let b = ReplayPlanSpecV1::capture(&p, &r, Some(&c), false, false).unwrap();
    assert_ne!(a, b);
    assert!(a.reconstruct().unwrap().2.is_none());
    assert!(b.reconstruct().unwrap().2.is_some());
}

#[test]
fn same_name_changed_numeric_kernel_changes_strong_plan_identity() {
    let r = refs();
    let mut a = PassPlan::new();
    let mut b = PassPlan::new();
    a.push(Box::new(
        S5FMutationPass::new_rate(S5FKernel::new(vec![1.0; 1024], vec![0.25; 4096]), 0.1)
            .with_kernel_name("name"),
    ));
    b.push(Box::new(
        S5FMutationPass::new_rate(S5FKernel::new(vec![0.5; 1024], vec![0.25; 4096]), 0.1)
            .with_kernel_name("name"),
    ));
    assert_ne!(
        ReplayPlanSpecV1::capture(&a, &r, None, false, false)
            .unwrap()
            .digest()
            .unwrap(),
        ReplayPlanSpecV1::capture(&b, &r, None, false, false)
            .unwrap()
            .digest()
            .unwrap()
    );
}

fn replay_pass(
    pass: &dyn Pass,
    sim: &crate::ir::Simulation,
    records: Vec<crate::trace::ChoiceRecord>,
    strict: bool,
    reference: Option<&RefDataConfig>,
    contracts: Option<&crate::contract::ContractSet>,
    family: bool,
) -> Result<crate::ir::Simulation, String> {
    let run = || {
        let mut trace = crate::trace::Trace::new();
        let mut rng = crate::rng::Rng::new(42);
        let mut cursor = crate::replay::TraceCursor::from_owned(records);
        let mut ctx = crate::pass::PassContext {
            trace: &mut trace,
            rng: &mut rng,
            pass_index: 0,
            refdata: reference,
            contracts,
            feasibility: None,
            reference_index: None,
            event_log_sink: None,
            replay_cursor: Some(&mut cursor),
        };
        pass.execute_with_policy(sim, &mut ctx, strict)
            .map_err(|e| e.to_string())
    };
    if family {
        crate::family_replay::with_replay_scope(run)
    } else {
        run()
    }
}
fn rec(address: &str, value: crate::trace::ChoiceValue) -> crate::trace::ChoiceRecord {
    crate::trace::ChoiceRecord::new(address, value)
}
fn fixed(n: i64) -> Box<dyn Distribution<Output = i64>> {
    Box::new(EmpiricalLengthDist::from_pairs([(n, 1.0)]))
}

#[test]
fn sampler_kinds_allele_id_order_and_zero_mass_entries_round_trip_exactly() {
    let mut r = refs();
    let mut extra = r.v_pool.iter().next().unwrap().1.clone();
    extra.name = "V*02".into();
    let _ = r.v_pool.push(extra);
    let samplers = [
        AllelePoolDist::uniform(&r.v_pool),
        AllelePoolDist::from_weights(&r.v_pool, vec![0.1, 0.2]),
        AllelePoolDist::restricted_uniform(&r.v_pool, vec![AlleleId::new(1), AlleleId::new(0)]),
    ];
    for original in samplers {
        let spec = original.family_replay_spec().unwrap();
        let rebuilt = spec.allele(&r.v_pool).unwrap();
        assert_eq!(spec, rebuilt.family_replay_spec().unwrap());
        for seed in 0..32 {
            let mut a = crate::rng::Rng::new(seed);
            let mut b = crate::rng::Rng::new(seed);
            for _ in 0..32 {
                assert_eq!(original.sample(&mut a), rebuilt.sample(&mut b));
            }
            assert_eq!(a.words_consumed(), b.words_consumed());
        }
    }
    let zero = DistributionSpec::AllelePool {
        ids: None,
        cumulative: vec![bits(0.0).unwrap(), bits(1.0).unwrap()],
        total: bits(1.0).unwrap(),
    };
    let rebuilt = zero.allele(&r.v_pool).unwrap();
    assert_eq!(zero, rebuilt.family_replay_spec().unwrap());
    for seed in 0..32 {
        assert_eq!(
            rebuilt.sample(&mut crate::rng::Rng::new(seed)),
            AlleleId::new(1)
        );
    }
    let zero_int = DistributionSpec::EmpiricalInt {
        values: vec![99, 7],
        cumulative: vec![bits(0.0).unwrap(), bits(1.0).unwrap()],
        total: bits(1.0).unwrap(),
    };
    assert!(!zero_int.supports_integer(99).unwrap());
    assert!(zero_int.supports_integer(7).unwrap());
    assert_eq!(
        zero_int,
        zero_int.integer().unwrap().family_replay_spec().unwrap()
    );
    assert_ne!(
        UniformBase.family_replay_spec(),
        CategoricalBase::from_weights([1.0; 4]).family_replay_spec()
    );
}

#[test]
fn orphan_reference_evidence_cannot_be_silently_dropped() {
    use crate::refdata::{AlleleEvidence, AnchorEvidenceStatus, CodingBoundaryStatus};
    let mut r = refs();
    r.evidence.insert(
        (Segment::V, "unknown*01".into()),
        AlleleEvidence {
            schema_version: 1,
            sequence_sha256: "sha256:orphan".into(),
            coding_phase: None,
            anchor_status: AnchorEvidenceStatus::Unknown,
            boundary_5: CodingBoundaryStatus::Unknown,
            boundary_3: CodingBoundaryStatus::Unknown,
            source: None,
            method: "test".into(),
            notes: vec![],
        },
    );
    assert!(ReferenceSpecV1::capture(&r).unwrap_err().contains("orphan"));
}

#[test]
fn family_count_replay_requires_positive_natural_support_in_both_policies() {
    use crate::trace::ChoiceValue::Int;
    let sim = crate::ir::Simulation::new();
    let pass = PCRErrorPass::new(fixed(0), Box::new(UniformBase));
    for strict in [false, true] {
        for count in [-1, 1, u32::MAX as i64 + 1] {
            let result = replay_pass(
                &pass,
                &sim,
                vec![rec("corrupt.pcr.count", Int(count))],
                strict,
                None,
                None,
                true,
            );
            assert!(result.is_err(), "{strict} {count}");
        }
        assert!(replay_pass(
            &pass,
            &sim,
            vec![rec("corrupt.pcr.count", Int(0))],
            strict,
            None,
            None,
            true
        )
        .is_ok());
    }
    // Legacy traces intentionally retain their historical range-only contract.
    assert!(replay_pass(
        &pass,
        &sim,
        vec![rec("corrupt.pcr.count", Int(1))],
        true,
        None,
        None,
        false
    )
    .is_ok());
    for rate in [0.0, 0.5] {
        let pass = PCRErrorPass::new_rate(rate, Box::new(UniformBase));
        for strict in [false, true] {
            assert!(replay_pass(
                &pass,
                &sim,
                vec![rec("corrupt.pcr.count", Int(1))],
                strict,
                None,
                None,
                true
            )
            .is_err());
        }
    }
}

#[test]
fn family_substitution_support_uses_the_actual_transformed_byte() {
    use crate::trace::ChoiceValue::{Base, Int};
    let (reference, sim) = crate::passes::test_support::make_substitution_productive_vj_fixture();
    let pass = QualityErrorPass::new(
        fixed(1),
        Box::new(CategoricalBase::from_weights([1.0, 0.0, 0.0, 0.0])),
    );
    for strict in [false, true] {
        for (base, valid) in [(b'a', true), (b'A', false), (b'c', false)] {
            let records = vec![
                rec("corrupt.quality.count", Int(1)),
                rec("corrupt.quality.error_site[0]", Int(0)),
                rec("corrupt.quality.error_base[0]", Base(base)),
            ];
            let result = replay_pass(&pass, &sim, records, strict, Some(&reference), None, true);
            assert_eq!(
                result.is_ok(),
                valid,
                "{strict} {} {result:?}",
                base as char
            );
        }
    }
}

#[test]
fn family_orientation_replay_rejects_zero_probability_without_panicking() {
    use crate::trace::ChoiceValue::Bool;
    for strict in [false, true] {
        for (p, value) in [(0.0, true), (1.0, false)] {
            let pass = RevCompPass::new(p);
            let records = vec![rec("corrupt.rev_comp.applied", Bool(value))];
            assert!(replay_pass(
                &pass,
                &crate::ir::Simulation::new(),
                records.clone(),
                strict,
                None,
                None,
                true
            )
            .is_err());
            assert!(replay_pass(
                &pass,
                &crate::ir::Simulation::new(),
                records,
                strict,
                None,
                None,
                false
            )
            .is_ok());
        }
    }
}

#[test]
fn end_loss_preimages_include_clamps_wrapping_and_zero_mass() {
    let fixed_spec = |raw| DistributionSpec::EmpiricalInt {
        values: vec![raw],
        cumulative: vec![bits(1.0).unwrap()],
        total: bits(1.0).unwrap(),
    };
    assert!(fixed_spec(-4).end_loss_preimage(0, 6, false).unwrap());
    assert!(fixed_spec(40).end_loss_preimage(6, 6, false).unwrap());
    assert!(!fixed_spec(2).end_loss_preimage(1, 6, false).unwrap());
    let wrap = DistributionSpec::UniformInt {
        min: (1i64 << 32) - 1,
        max: (1i64 << 32) + 2,
    };
    for n in [0, 1, 6] {
        assert!(wrap.end_loss_preimage(n, 6, false).unwrap());
    }
    assert!(!wrap.end_loss_preimage(2, 6, false).unwrap());
    assert!(wrap.end_loss_preimage(1, 6, true).unwrap());
    let zero = DistributionSpec::EmpiricalInt {
        values: vec![1, 3],
        cumulative: vec![bits(0.0).unwrap(), bits(1.0).unwrap()],
        total: bits(1.0).unwrap(),
    };
    assert!(!zero.end_loss_preimage(1, 6, false).unwrap());
}

#[test]
fn family_end_loss_validates_realized_count_and_permissive_contract_stop() {
    use crate::trace::ChoiceValue::Int;
    let (reference, sim) = crate::passes::test_support::make_substitution_productive_vj_fixture();
    let pass = EndLossPass::new(LossEnd::Five, fixed(2));
    for strict in [false, true] {
        for count in [0, 1, u32::MAX as i64 + 1] {
            assert!(replay_pass(
                &pass,
                &sim,
                vec![rec("corrupt.end_loss.5", Int(count))],
                strict,
                Some(&reference),
                None,
                true
            )
            .is_err());
        }
        assert!(replay_pass(
            &pass,
            &sim,
            vec![rec("corrupt.end_loss.5", Int(2))],
            strict,
            Some(&reference),
            None,
            true
        )
        .is_ok());
    }
    let contracts = crate::contract::productive();
    let records = vec![rec("corrupt.end_loss.5", Int(0))];
    assert!(replay_pass(
        &pass,
        &sim,
        records.clone(),
        false,
        Some(&reference),
        Some(&contracts),
        true
    )
    .is_ok());
    assert!(replay_pass(
        &pass,
        &sim,
        records,
        true,
        Some(&reference),
        Some(&contracts),
        true
    )
    .is_err());
}

#[test]
fn family_rate_support_respects_native_poisson_cap_and_numerical_zero() {
    use crate::trace::ChoiceValue::Int;
    let (reference, sim) = crate::passes::test_support::make_substitution_productive_vj_fixture();
    for strict in [false, true] {
        for (rate, count) in [(f64::MIN_POSITIVE, 1), (0.5, 1_000_002)] {
            let pass = PCRErrorPass::new_rate(rate, Box::new(UniformBase));
            let error = replay_pass(
                &pass,
                &sim,
                vec![rec("corrupt.pcr.count", Int(count))],
                strict,
                Some(&reference),
                None,
                true,
            )
            .unwrap_err();
            assert!(error.contains("positive natural support"), "{error}");
        }
    }
}

#[test]
fn count_codecs_reject_invalid_positive_support_but_preserve_end_loss_clamps() {
    for value in [-1, u32::MAX as i64 + 1] {
        let spec = fixed(value).family_replay_spec().unwrap();
        assert!(spec.count_integer().is_err());
        assert!(spec.integer().is_ok());
        let pass = ReplayPassSpecV1::Pcr {
            count: CountSpec::Distribution { distribution: spec },
            base: DistributionSpec::UniformBase,
        };
        assert!(pass.reconstruct(&refs()).is_err());
    }
    for (min, max) in [(-1, 1), (u32::MAX as i64, u32::MAX as i64 + 2)] {
        assert!(DistributionSpec::UniformInt { min, max }
            .count_integer()
            .is_err());
    }
    assert!(DistributionSpec::UniformInt {
        min: u32::MAX as i64,
        max: u32::MAX as i64 + 1
    }
    .count_integer()
    .is_ok());
    let disabled = DistributionSpec::EmpiricalInt {
        values: vec![-1, 2],
        cumulative: vec![bits(0.0).unwrap(), bits(1.0).unwrap()],
        total: bits(1.0).unwrap(),
    };
    assert!(disabled.count_integer().is_ok());
}
