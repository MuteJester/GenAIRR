use super::*;
use crate::duplex_repair::model::FiniteDuplexModel;
use crate::duplex_repair::{RepairOperation, Side, Strand};
use crate::genotype::{GeneCopy, Genotype, Haplotype};
use crate::passes::sample_genotype::SampleGenotypePass;
use crate::passes::{DuplexRecombinePass, InvertDPass};
use crate::refdata::GeneId;
use std::sync::Arc;

fn fixture(d_sequences: &[&str]) -> (RefDataConfig, FiniteDuplexModel) {
    let (mut cfg, mut model) = crate::duplex_repair::model::tests::fixture();
    // Keep allele 0 from the shared fixture and replace its sequence by
    // constructing the D pool afresh, so each test names the exact substrate.
    cfg.d_pool = Default::default();
    let template = model.rows[0].clone();
    model.rows.clear();
    for (index, sequence) in d_sequences.iter().enumerate() {
        let name = format!("D{index}*01");
        let _ = cfg.d_pool.push(Allele {
            name: name.clone(),
            gene: format!("D{index}"),
            seq: sequence.as_bytes().to_vec(),
            segment: Segment::D,
            anchor: None,
            functional_status: None,
            subregions: vec![],
        });
        let mut row = template.clone();
        row.id = format!("history{index}");
        row.d_call = name;
        row.script.d = (*sequence).into();
        row.script.joins[1].operations[2] = RepairOperation::Synapse {
            right_start: sequence.len() as i64,
        };
        model.rows.push(row);
    }
    model.reference_content_hash = crate::trace_file::refdata_content_hash(&cfg);
    (cfg, model)
}

fn flat_plan(
    cfg: &RefDataConfig,
    model: FiniteDuplexModel,
    d_ids: Option<Vec<AlleleId>>,
    inversion: Option<f64>,
) -> PassPlan {
    let mut plan = PassPlan::new();
    for segment in [Segment::V, Segment::D, Segment::J] {
        let pool = cfg.pool_for(segment).unwrap();
        let dist = match (&d_ids, segment) {
            (Some(ids), Segment::D) => AllelePoolDist::restricted_uniform(pool, ids.clone()),
            _ => AllelePoolDist::uniform(pool),
        };
        plan.push(Box::new(SampleAllelePass::new(segment, Box::new(dist))));
    }
    if let Some(prob) = inversion {
        plan.push(Box::new(InvertDPass::new(prob)));
    }
    plan.push(Box::new(DuplexRecombinePass::new(Arc::new(model)).unwrap()));
    plan
}

fn failure(plan: &PassPlan, cfg: &RefDataConfig, contracts: Option<&ContractSet>) -> CompileErrors {
    match compile_test_fixture(plan, Some(cfg), contracts, ExecutionPolicy::Permissive) {
        Err(error) => error,
        Ok(_) => panic!("invalid finite support must not compile"),
    }
}

#[test]
fn every_natural_context_is_required_before_productivity_masking_but_locks_narrow_it() {
    let (cfg, mut model) = fixture(&["TAA", "CCC"]);
    model.rows.remove(0);
    let contracts = crate::contract::productive();
    let plan = flat_plan(&cfg, model.clone(), None, None);
    let error = failure(&plan, &cfg, Some(&contracts));
    assert!(error.to_string().contains("lacks positive history support"));
    let plan = flat_plan(&cfg, model, Some(vec![AlleleId::new(1)]), None);
    let compiled =
        compile_test_fixture(&plan, Some(&cfg), Some(&contracts), ExecutionPolicy::Strict).unwrap();
    assert_eq!(
        compiled
            .run_one(1)
            .unwrap()
            .final_simulation()
            .assignments
            .get(Segment::D)
            .unwrap()
            .allele_id,
        AlleleId::new(1)
    );
}

#[test]
fn invalid_active_history_cannot_be_hidden_by_an_unproductive_source() {
    let (cfg, mut model) = fixture(&["TAA", "CCC"]);
    model.rows[0].script.joins[0].operations[0] = RepairOperation::Open {
        side: Side::Left,
        offset: 100,
    };
    let contracts = crate::contract::productive();
    let plan = flat_plan(&cfg, model.clone(), None, None);
    let error = failure(&plan, &cfg, Some(&contracts));
    assert!(!error
        .to_string()
        .contains("no positive admissible completion"));
    // A positive but unreachable row is still bound by the digest; it need
    // not be executable for a deliberately locked, different context.
    let plan = flat_plan(&cfg, model, Some(vec![AlleleId::new(1)]), None);
    compile_test_fixture(&plan, Some(&cfg), Some(&contracts), ExecutionPolicy::Strict).unwrap();
}

#[test]
fn productive_masks_restrict_alleles_in_both_execution_modes_and_refuse_bad_replay() {
    let (cfg, model) = fixture(&["TAA", "CCC"]);
    for policy in [ExecutionPolicy::Strict, ExecutionPolicy::Permissive] {
        let compiled = compile_test_fixture_owned(
            flat_plan(&cfg, model.clone(), None, None),
            Some(cfg.clone()),
            Some(crate::contract::productive()),
            policy,
        )
        .unwrap();
        for seed in 0..32 {
            let outcome = compiled.run_one(seed).unwrap();
            let sim = outcome.final_simulation();
            assert_eq!(
                sim.assignments.get(Segment::D).unwrap().allele_id.index(),
                1
            );
            assert_eq!(
                crate::productivity::assess(sim, &cfg).productive,
                Some(true)
            );
        }
        let original = compiled.run_one(12).unwrap();
        let mut choices = original.trace.choices().to_vec();
        choices
            .iter_mut()
            .find(|c| c.address == "sample_allele.d")
            .unwrap()
            .value = ChoiceValue::AlleleId(0);
        assert!(compiled
            .replay_from_trace_records(&choices, 99, policy)
            .is_err());
    }
    let (cfg, model) = fixture(&["TAA"]);
    let plan = flat_plan(&cfg, model, None, None);
    let contracts = crate::contract::productive();
    assert!(failure(&plan, &cfg, Some(&contracts))
        .to_string()
        .contains("no positive admissible completion"));
}

#[test]
fn orientation_is_a_future_local_choice_and_only_positive_orientations_need_rows() {
    let (cfg, mut model) = fixture(&["TAA"]);
    let mut reverse = model.rows[0].clone();
    reverse.id = "reverse".into();
    reverse.script.d_inverted = true;
    model.rows.push(reverse);
    let contracts = crate::contract::productive();
    for policy in [ExecutionPolicy::Strict, ExecutionPolicy::Permissive] {
        let plan = flat_plan(&cfg, model.clone(), None, Some(0.25));
        let compiled = compile_test_fixture(&plan, Some(&cfg), Some(&contracts), policy).unwrap();
        for seed in 0..32 {
            let outcome = compiled.run_one(seed).unwrap();
            assert!(outcome
                .final_simulation()
                .assignments
                .get(Segment::D)
                .unwrap()
                .orientation
                .is_reverse());
        }
    }
    model.rows.pop();
    let plan = flat_plan(&cfg, model.clone(), None, Some(1.0));
    assert!(failure(&plan, &cfg, None)
        .to_string()
        .contains("lacks positive history support"));
    let plan = flat_plan(&cfg, model, None, Some(0.0));
    compile_test_fixture(&plan, Some(&cfg), None, ExecutionPolicy::Strict).unwrap();
}

#[test]
fn original_stop_can_be_removed_by_authored_repair_before_coding_assessment() {
    let (mut cfg, mut model) = fixture(&["CCC"]);
    cfg.v_pool = Default::default();
    let _ = cfg.v_pool.push(Allele {
        name: "V*01".into(),
        gene: "V".into(),
        seq: b"AAATGTTAACCC".to_vec(),
        segment: Segment::V,
        anchor: Some(3),
        functional_status: None,
        subregions: vec![],
    });
    model.rows[0].script.v = "AAATGTTAACCC".into();
    model.rows[0].script.joins[0].operations.splice(
        2..2,
        [
            RepairOperation::Resect {
                side: Side::Left,
                strand: Strand::Top,
                count: 6,
            },
            RepairOperation::Resect {
                side: Side::Left,
                strand: Strand::Bottom,
                count: 6,
            },
        ],
    );
    model.reference_content_hash = crate::trace_file::refdata_content_hash(&cfg);
    let plan = flat_plan(&cfg, model, None, None);
    let contracts = crate::contract::productive();
    for policy in [ExecutionPolicy::Strict, ExecutionPolicy::Permissive] {
        let compiled = compile_test_fixture(&plan, Some(&cfg), Some(&contracts), policy).unwrap();
        let outcome = compiled.run_one(1).unwrap();
        let sim = outcome.final_simulation();
        assert_eq!(sim.assignments.get(Segment::V).unwrap().trim_3, 6);
        assert_eq!(
            crate::productivity::assess(sim, &cfg).productive,
            Some(true)
        );
    }
}

#[test]
fn reservation_and_operation_limits_are_explicit_errors_not_biological_zero_support() {
    let (cfg, model) = fixture(&["TAA", "CCC"]);
    let contracts = crate::contract::productive();
    for bases in [false, true] {
        let mut limited = model.clone();
        if bases {
            limited.max_preflight_bases = 1;
        } else {
            limited.max_preflight_operations = 1;
        }
        let plan = flat_plan(&cfg, limited, None, None);
        let error = failure(&plan, &cfg, Some(&contracts));
        assert!(error
            .errors
            .iter()
            .any(|e| matches!(e.kind, CompileErrorKind::ResourceLimit { .. })));
        assert!(!error
            .to_string()
            .contains("no positive admissible completion"));
    }
}

#[test]
fn positive_weight_underflow_is_rejected_before_productivity_filtering() {
    let (cfg, mut model) = fixture(&["TAA", "CCC"]);
    model.rows[0].weight = f64::from_bits(1);
    let mut duplicate = model.rows[0].clone();
    duplicate.id = "large_same_context".into();
    duplicate.weight = 1e308;
    model.rows.push(duplicate);
    let plan = flat_plan(&cfg, model, None, None);
    let contracts = crate::contract::productive();
    assert!(failure(&plan, &cfg, Some(&contracts))
        .to_string()
        .contains("numerical"));
}

fn genotype_plan(
    cfg: &RefDataConfig,
    model: FiniteDuplexModel,
    chromosome_weights: [f32; 2],
    usage_d: Vec<(GeneId, f64)>,
) -> PassPlan {
    let mut haplotypes = [Haplotype::new(), Haplotype::new()];
    for (c, haplotype) in haplotypes.iter_mut().enumerate() {
        for segment in [Segment::V, Segment::D, Segment::J] {
            let id = if segment == Segment::D { c as u32 } else { 0 };
            haplotype.set(
                segment,
                GeneId::new(id),
                vec![GeneCopy {
                    allele: AlleleId::new(id),
                    copies: 1,
                    weight: 1.0,
                }],
            );
        }
    }
    let genotype = Genotype::new(
        haplotypes,
        chromosome_weights,
        Some("donor".into()),
        crate::trace_file::refdata_content_hash(cfg),
    );
    let mut plan = PassPlan::new();
    plan.push(Box::new(SampleGenotypePass::new(
        Arc::new(genotype),
        true,
        vec![],
        usage_d,
        vec![],
    )));
    plan.push(Box::new(DuplexRecombinePass::new(Arc::new(model)).unwrap()));
    plan
}

#[test]
fn genotype_future_choices_stay_on_selected_haplotype_and_keep_donor_support() {
    let (cfg, model) = fixture(&["TAA", "CCC"]);
    for policy in [ExecutionPolicy::Strict, ExecutionPolicy::Permissive] {
        let plan = genotype_plan(&cfg, model.clone(), [9.0, 1.0], vec![]);
        let contracts = crate::contract::productive();
        let compiled = compile_test_fixture(&plan, Some(&cfg), Some(&contracts), policy).unwrap();
        for seed in 0..32 {
            let outcome = compiled.run_one(seed).unwrap();
            let sim = outcome.final_simulation();
            for segment in [Segment::V, Segment::D, Segment::J] {
                assert_eq!(sim.assignments.get(segment).unwrap().haplotype, Some(1));
            }
            assert_eq!(
                sim.assignments.get(Segment::D).unwrap().allele_id.index(),
                1
            );
            assert!(sim.live_call_support.is_some());
        }
    }
}

#[test]
fn zero_weight_haplotype_or_gene_does_not_require_inactive_context() {
    let (cfg, mut model) = fixture(&["TAA", "CCC"]);
    model.rows.remove(0);
    for (weights, usage) in [
        ([0.0, 1.0], vec![]),
        ([1.0, 1.0], vec![(GeneId::new(0), 0.0)]),
    ] {
        let plan = genotype_plan(&cfg, model.clone(), weights, usage);
        let compiled =
            compile_test_fixture(&plan, Some(&cfg), None, ExecutionPolicy::Permissive).unwrap();
        for seed in 0..16 {
            assert_eq!(
                compiled
                    .run_one(seed)
                    .unwrap()
                    .final_simulation()
                    .assignments
                    .get(Segment::D)
                    .unwrap()
                    .haplotype,
                Some(1)
            );
        }
    }
}

#[test]
fn unsupported_suffix_or_legacy_assembly_fails_before_generation() {
    let (cfg, model) = fixture(&["CCC"]);
    let mut plan = flat_plan(&cfg, model.clone(), None, None);
    plan.push(Box::new(AssembleSegmentPass::new(Segment::V)));
    assert!(failure(&plan, &cfg, None)
        .to_string()
        .contains("atomic duplex repair"));
    let mut plan = flat_plan(&cfg, model, None, None);
    plan.push(Box::new(TraceProbePass::new("unknown", "unknown", b'A')));
    assert!(failure(&plan, &cfg, None)
        .to_string()
        .contains("unsupported"));
}

#[test]
fn finite_native_genotype_rejects_nonfinite_and_unrepresentable_weights() {
    let (cfg, model) = fixture(&["CCC", "ACT"]);
    for weight in [f64::NAN, f64::INFINITY, -1.0] {
        let plan = genotype_plan(
            &cfg,
            model.clone(),
            [1.0, 1.0],
            vec![(GeneId::new(0), weight)],
        );
        assert!(failure(&plan, &cfg, None)
            .to_string()
            .contains("finite and nonnegative"));
    }
    let plan = genotype_plan(&cfg, model, [1.0, 1e-30], vec![]);
    assert!(failure(&plan, &cfg, None)
        .to_string()
        .contains("numerical support"));
}

#[test]
fn genotype_contexts_are_joint_not_the_cross_haplotype_cartesian_product() {
    let (mut cfg, mut model) = fixture(&["CCC", "ACT"]);
    for segment in [Segment::V, Segment::J] {
        let pool = if segment == Segment::V {
            &mut cfg.v_pool
        } else {
            &mut cfg.j_pool
        };
        let mut second = pool.get(AlleleId::new(0)).unwrap().clone();
        second.name = format!("{segment:?}2*01");
        second.gene = format!("{segment:?}2");
        let _ = pool.push(second);
    }
    model.rows[1].v_call = "V2*01".into();
    model.rows[1].j_call = "J2*01".into();
    model.reference_content_hash = crate::trace_file::refdata_content_hash(&cfg);
    let mut haplotypes = [Haplotype::new(), Haplotype::new()];
    for (c, haplotype) in haplotypes.iter_mut().enumerate() {
        for segment in [Segment::V, Segment::D, Segment::J] {
            haplotype.set(
                segment,
                GeneId::new(c as u32),
                vec![GeneCopy {
                    allele: AlleleId::new(c as u32),
                    copies: 1,
                    weight: 1.0,
                }],
            );
        }
    }
    let genotype = Genotype::new(
        haplotypes,
        [1.0, 1.0],
        None,
        crate::trace_file::refdata_content_hash(&cfg),
    );
    let mut plan = PassPlan::new();
    plan.push(Box::new(SampleGenotypePass::new(
        Arc::new(genotype),
        true,
        vec![],
        vec![],
        vec![],
    )));
    plan.push(Box::new(DuplexRecombinePass::new(Arc::new(model)).unwrap()));
    let contracts = crate::contract::productive();
    let compiled =
        compile_test_fixture(&plan, Some(&cfg), Some(&contracts), ExecutionPolicy::Strict).unwrap();
    let mut seen = [false; 2];
    for seed in 0..64 {
        let outcome = compiled.run_one(seed).unwrap();
        let sim = outcome.final_simulation();
        let haplotype = sim.assignments.get(Segment::V).unwrap().haplotype.unwrap();
        seen[usize::from(haplotype)] = true;
        for segment in [Segment::V, Segment::D, Segment::J] {
            assert_eq!(
                sim.assignments.get(segment).unwrap().allele_id.index(),
                u32::from(haplotype)
            );
        }
    }
    assert_eq!(seen, [true, true]);
}

#[test]
fn history_count_and_mass_do_not_reweight_upstream_allele_choices() {
    let (cfg, mut model) = fixture(&["CCC", "ACT"]);
    model.rows[0].weight = 100.0;
    model.rows[1].weight = 0.01;
    for index in 0..19 {
        let mut extra = model.rows[0].clone();
        extra.id = format!("extra{index}");
        model.rows.push(extra);
    }
    model.max_preflight_bases = 10000;
    let mut plan = PassPlan::new();
    for segment in [Segment::V, Segment::D, Segment::J] {
        let pool = cfg.pool_for(segment).unwrap();
        let dist = if segment == Segment::D {
            AllelePoolDist::from_weights(pool, vec![1.0, 3.0])
        } else {
            AllelePoolDist::uniform(pool)
        };
        plan.push(Box::new(SampleAllelePass::new(segment, Box::new(dist))));
    }
    plan.push(Box::new(DuplexRecombinePass::new(Arc::new(model)).unwrap()));
    let contracts = crate::contract::productive();
    let compiled =
        compile_test_fixture(&plan, Some(&cfg), Some(&contracts), ExecutionPolicy::Strict).unwrap();
    let n = 2048;
    let mut d1 = 0;
    for seed in 0..n {
        let outcome = compiled.run_one(seed).unwrap();
        d1 += usize::from(
            outcome
                .final_simulation()
                .assignments
                .get(Segment::D)
                .unwrap()
                .allele_id
                .index()
                == 1,
        );
    }
    let expected = n as f64 * 0.75;
    let sigma = (n as f64 * 0.75 * 0.25).sqrt();
    assert!((d1 as f64 - expected).abs() < 5.0 * sigma, "{d1}/{n}");
}

#[test]
fn zero_mass_gene_replay_is_rejected_even_if_same_allele_is_positive_elsewhere() {
    let (cfg, model) = fixture(&["CCC"]);
    let mut hap = Haplotype::new();
    for segment in [Segment::V, Segment::D, Segment::J] {
        hap.set(
            segment,
            GeneId::new(0),
            vec![GeneCopy {
                allele: AlleleId::new(0),
                copies: 1,
                weight: 1.0,
            }],
        );
    }
    hap.set(
        Segment::D,
        GeneId::new(1),
        vec![GeneCopy {
            allele: AlleleId::new(0),
            copies: 1,
            weight: 1.0,
        }],
    );
    let genotype = Genotype::new(
        [hap, Haplotype::new()],
        [1.0, 0.0],
        None,
        crate::trace_file::refdata_content_hash(&cfg),
    );
    let mut plan = PassPlan::new();
    plan.push(Box::new(SampleGenotypePass::new(
        Arc::new(genotype),
        true,
        vec![],
        vec![(GeneId::new(1), 0.0)],
        vec![],
    )));
    plan.push(Box::new(DuplexRecombinePass::new(Arc::new(model)).unwrap()));
    let compiled =
        compile_test_fixture_owned(plan, Some(cfg), None, ExecutionPolicy::Permissive).unwrap();
    let original = compiled.run_one(0).unwrap();
    let mut choices = original.trace.choices().to_vec();
    choices
        .iter_mut()
        .find(|c| c.address == "sample_gene.d")
        .unwrap()
        .value = ChoiceValue::GeneId(1);
    assert!(compiled
        .replay_from_trace_records(&choices, 0, ExecutionPolicy::Permissive)
        .is_err());
}

#[test]
fn preservation_only_contract_is_an_explicit_capability_error() {
    let (cfg, model) = fixture(&["CCC"]);
    let plan = flat_plan(&cfg, model, None, None);
    let contracts =
        ContractSet::new().with(Box::new(crate::contract::AnchorPreserved::new(Segment::V)));
    assert!(failure(&plan, &cfg, Some(&contracts))
        .errors
        .iter()
        .any(|error| matches!(
            error.kind,
            CompileErrorKind::UnsupportedContractCombo { .. }
        )));
}

#[test]
fn empty_productive_support_reports_unknown_separately_from_false() {
    let (cfg, model) = fixture(&["TAA"]);
    let contracts = crate::contract::productive();
    let plan = flat_plan(&cfg, model, None, None);
    let error = failure(&plan, &cfg, Some(&contracts)).to_string();
    assert!(error.contains("proven_false=1, unresolved=0"), "{error}");

    let (mut cfg, mut model) = fixture(&["CCC"]);
    let mut unannotated = cfg.v_pool.get(AlleleId::new(0)).unwrap().clone();
    unannotated.anchor = None;
    cfg.v_pool = Default::default();
    let _ = cfg.v_pool.push(unannotated);
    model.reference_content_hash = crate::trace_file::refdata_content_hash(&cfg);
    let plan = flat_plan(&cfg, model, None, None);
    let error = failure(&plan, &cfg, Some(&contracts)).to_string();
    assert!(error.contains("proven_false=0, unresolved=1"), "{error}");
}

#[test]
fn physically_valid_ambiguous_body_projection_cannot_be_productivity_filtered_away() {
    let (cfg, mut model) = fixture(&["TAA", "CCC"]);
    let script = &mut model.rows[0].script;
    script.joins[0].operations = vec![
        RepairOperation::Open {
            side: Side::Left,
            offset: 0,
        },
        RepairOperation::Open {
            side: Side::Right,
            offset: 0,
        },
        RepairOperation::Resect {
            side: Side::Left,
            strand: Strand::Bottom,
            count: 1,
        },
        RepairOperation::Resect {
            side: Side::Right,
            strand: Strand::Top,
            count: 1,
        },
        RepairOperation::Synapse { right_start: 5 },
        RepairOperation::Ligate {},
    ];
    // The V-terminal T and D-initial T occupy one valid paired column;
    // two body claims are meaningful, but not representable in this IR.
    crate::duplex_repair::repair_coupled_annotated(
        &script.v,
        &script.d,
        &script.j,
        false,
        &script.joins,
        script.max_operations,
        script.max_bases,
    )
    .unwrap();
    let contracts = crate::contract::productive();
    let plan = flat_plan(&cfg, model, None, None);
    let error = failure(&plan, &cfg, Some(&contracts)).to_string();
    assert!(error.contains("history history0"), "{error}");
    assert!(!error.contains("no positive admissible completion"));
}

#[test]
fn finite_sampling_signature_binds_genotype_copy_order_and_usage() {
    let make = |reverse: bool, weight: f64| {
        let mut hap = Haplotype::new();
        let mut copies = vec![
            GeneCopy {
                allele: AlleleId::new(0),
                copies: 1,
                weight: 1.0,
            },
            GeneCopy {
                allele: AlleleId::new(1),
                copies: 1,
                weight: 2.0,
            },
        ];
        if reverse {
            copies.reverse();
        }
        hap.set(Segment::V, GeneId::new(0), copies);
        SampleGenotypePass::new(
            Arc::new(Genotype::new(
                [hap, Haplotype::new()],
                [1.0, 0.0],
                None,
                "source".into(),
            )),
            true,
            vec![(GeneId::new(0), weight)],
            vec![],
            vec![],
        )
    };
    let original = make(false, 1.0);
    for changed in [make(true, 1.0), make(false, 2.0)] {
        assert_eq!(
            original.parameter_signature(),
            changed.parameter_signature()
        );
        assert_ne!(
            original.sampling_law_signature(),
            changed.sampling_law_signature()
        );
    }
}

struct ExactWeights(Vec<(AlleleId, f64)>);

impl crate::dist::Distribution for ExactWeights {
    type Output = AlleleId;
    fn sample(&self, _rng: &mut crate::rng::Rng) -> AlleleId {
        panic!("finite source must use its explicit support");
    }
    fn support(&self) -> Option<Vec<(AlleleId, f64)>> {
        Some(self.0.clone())
    }
}

fn exact_weight_plan(cfg: &RefDataConfig, model: FiniteDuplexModel, weights: &[f64]) -> PassPlan {
    let mut plan = PassPlan::new();
    for segment in [Segment::V, Segment::D, Segment::J] {
        let dist: Box<dyn crate::dist::Distribution<Output = AlleleId>> = if segment == Segment::D {
            Box::new(ExactWeights(
                weights
                    .iter()
                    .enumerate()
                    .map(|(id, w)| (AlleleId::new(id as u32), *w))
                    .collect(),
            ))
        } else {
            Box::new(AllelePoolDist::uniform(cfg.pool_for(segment).unwrap()))
        };
        plan.push(Box::new(SampleAllelePass::new(segment, dist)));
    }
    plan.push(Box::new(DuplexRecombinePass::new(Arc::new(model)).unwrap()));
    plan
}

#[test]
fn native_advertised_allele_weights_are_checked_for_overflow_and_cdf_loss() {
    let (cfg, model) = fixture(&["CCC", "ACT"]);
    for weights in [[1e308, 1e308], [1.0, 1e-20]] {
        let plan = exact_weight_plan(&cfg, model.clone(), &weights);
        let error = failure(&plan, &cfg, None).to_string();
        assert!(error.contains("numerical"), "{error}");
    }
}

#[test]
fn masked_half_ulp_allele_mass_is_rejected_at_compile_even_when_natural_cdf_advances() {
    let weights = [1.0, f64::EPSILON, f64::EPSILON / 2.0];
    crate::passes::sample_allele::validate_finite_weights(&weights).unwrap();
    assert!(
        crate::passes::sample_allele::validate_finite_weights(&[weights[0], weights[2]]).is_err()
    );
    let (cfg, model) = fixture(&["CCC", "TAA", "ACT"]);
    let plan = exact_weight_plan(&cfg, model, &weights);
    let contracts = crate::contract::productive();
    let error = failure(&plan, &cfg, Some(&contracts)).to_string();
    assert!(error.contains("after completion masking"), "{error}");
}

#[test]
fn numeric_support_preflight_is_charged_to_the_declared_operation_budget() {
    let (cfg, mut model) = fixture(&["CCC"]);
    // One context + eight repair operations + three natural source atoms
    // + three distinct masked source atoms = fifteen operations.
    model.max_preflight_operations = 14;
    let plan = flat_plan(&cfg, model.clone(), None, None);
    assert!(failure(&plan, &cfg, None).errors.iter().any(|e| matches!(
        e.kind,
        CompileErrorKind::ResourceLimit {
            required: 15,
            limit: 14,
            ..
        }
    )));
    model.max_preflight_operations = 15;
    let plan = flat_plan(&cfg, model, None, None);
    compile_test_fixture(&plan, Some(&cfg), None, ExecutionPolicy::Strict).unwrap();
}

#[test]
fn genotype_masked_half_ulp_mass_errors_before_any_local_draw_in_both_modes() {
    let (cfg, model) = fixture(&["CCC", "TAA", "ACT"]);
    let mut hap = Haplotype::new();
    for segment in [Segment::V, Segment::J] {
        hap.set(
            segment,
            GeneId::new(0),
            vec![GeneCopy {
                allele: AlleleId::new(0),
                copies: 1,
                weight: 1.0,
            }],
        );
    }
    for id in 0..3 {
        hap.set(
            Segment::D,
            GeneId::new(id),
            vec![GeneCopy {
                allele: AlleleId::new(id),
                copies: 1,
                weight: 1.0,
            }],
        );
    }
    let genotype = Arc::new(Genotype::new(
        [hap, Haplotype::new()],
        [1.0, 0.0],
        None,
        crate::trace_file::refdata_content_hash(&cfg),
    ));
    for policy in [ExecutionPolicy::Strict, ExecutionPolicy::Permissive] {
        let mut plan = PassPlan::new();
        plan.push(Box::new(SampleGenotypePass::new(
            genotype.clone(),
            true,
            vec![],
            vec![
                (GeneId::new(0), 1.0),
                (GeneId::new(1), f64::EPSILON),
                (GeneId::new(2), f64::EPSILON / 2.0),
            ],
            vec![],
        )));
        plan.push(Box::new(
            DuplexRecombinePass::new(Arc::new(model.clone())).unwrap(),
        ));
        let contracts = crate::contract::productive();
        let compiled = compile_test_fixture(&plan, Some(&cfg), Some(&contracts), policy).unwrap();
        let error = compiled.run_one(0).unwrap_err().to_string();
        assert!(error.contains("numerical"), "{error}");
    }
}

#[test]
fn repeated_weight_atoms_share_one_biological_context_without_changing_the_cdf() {
    let (cfg, mut model) = fixture(&["CCC"]);
    let atoms = 64;
    model.max_preflight_operations = 6 * atoms + 9;
    let mut plan = PassPlan::new();
    for segment in [Segment::V, Segment::D, Segment::J] {
        plan.push(Box::new(SampleAllelePass::new(
            segment,
            Box::new(ExactWeights(vec![(AlleleId::new(0), 1.0); atoms as usize])),
        )));
    }
    plan.push(Box::new(DuplexRecombinePass::new(Arc::new(model)).unwrap()));
    let compiled = compile_test_fixture(&plan, Some(&cfg), None, ExecutionPolicy::Strict).unwrap();
    assert!(compiled
        .run_one(3)
        .unwrap()
        .final_simulation()
        .duplex_provenance
        .is_some());
}

#[test]
fn genotype_keeps_natural_gene_dosage_when_only_some_copies_have_productive_completion() {
    let (mut cfg, mut model) = fixture(&["CCC", "TAA", "ACT"]);
    let old_d = cfg.d_pool.clone();
    cfg.d_pool = Default::default();
    for (id, (name, gene)) in [("DA*01", "DA"), ("DA*02", "DA"), ("DB*01", "DB")]
        .into_iter()
        .enumerate()
    {
        let mut allele = old_d.get(AlleleId::new(id as u32)).unwrap().clone();
        allele.name = name.into();
        allele.gene = gene.into();
        let _ = cfg.d_pool.push(allele);
        model.rows[id].d_call = name.into();
    }
    model.reference_content_hash = crate::trace_file::refdata_content_hash(&cfg);
    let mut hap = Haplotype::new();
    for segment in [Segment::V, Segment::J] {
        hap.set(
            segment,
            GeneId::new(0),
            vec![GeneCopy {
                allele: AlleleId::new(0),
                copies: 1,
                weight: 1.0,
            }],
        );
    }
    hap.set(
        Segment::D,
        GeneId::new(0),
        vec![
            GeneCopy {
                allele: AlleleId::new(0),
                copies: 1,
                weight: 1.0,
            },
            GeneCopy {
                allele: AlleleId::new(1),
                copies: 1,
                weight: 9.0,
            },
        ],
    );
    hap.set(
        Segment::D,
        GeneId::new(1),
        vec![GeneCopy {
            allele: AlleleId::new(2),
            copies: 1,
            weight: 1.0,
        }],
    );
    let genotype = Arc::new(Genotype::new(
        [hap, Haplotype::new()],
        [1.0, 0.0],
        None,
        crate::trace_file::refdata_content_hash(&cfg),
    ));
    // Independent SplitMix64 word-four oracle. The preceding choices are
    // haplotype, V gene and V copy; all are singletons but each consumes a word.
    let d_gene_uniform = |seed: u64| {
        let mut word = seed.wrapping_add(0x9e37_79b9_7f4a_7c15u64.wrapping_mul(4));
        word = (word ^ (word >> 30)).wrapping_mul(0xbf58_476d_1ce4_e5b9);
        word = (word ^ (word >> 27)).wrapping_mul(0x94d0_49bb_1331_11eb);
        word ^= word >> 31;
        ((word >> 11) as f64) / ((1u64 << 53) as f64)
    };
    for policy in [ExecutionPolicy::Strict, ExecutionPolicy::Permissive] {
        let mut plan = PassPlan::new();
        plan.push(Box::new(SampleGenotypePass::new(
            genotype.clone(),
            true,
            vec![],
            vec![],
            vec![],
        )));
        plan.push(Box::new(
            DuplexRecombinePass::new(Arc::new(model.clone())).unwrap(),
        ));
        let compiled = compile_test_fixture_owned(
            plan,
            Some(cfg.clone()),
            Some(crate::contract::productive()),
            policy,
        )
        .unwrap();
        let mut counts = [0; 2];
        let mut differs_from_filtered_dosage = 0;
        for seed in 1..=256 {
            let u = d_gene_uniform(seed);
            let expected_gene = u32::from(u * 11.0 >= 10.0);
            differs_from_filtered_dosage += usize::from(expected_gene != u32::from(u * 2.0 >= 1.0));
            let outcome = compiled.run_one(seed).unwrap();
            let chosen = outcome
                .trace
                .choices()
                .iter()
                .find(|choice| choice.address == "sample_gene.d")
                .unwrap();
            assert_eq!(
                chosen.value,
                ChoiceValue::GeneId(expected_gene),
                "seed {seed}"
            );
            counts[expected_gene as usize] += 1;
            let expected_allele = if expected_gene == 0 { 0 } else { 2 };
            assert_eq!(
                outcome
                    .final_simulation()
                    .assignments
                    .get(Segment::D)
                    .unwrap()
                    .allele_id
                    .index(),
                expected_allele
            );
            let replay = compiled
                .replay_from_trace_records(outcome.trace.choices(), 9876, policy)
                .unwrap();
            assert_eq!(outcome.trace.choices(), replay.trace.choices());
            assert_eq!(
                outcome.final_simulation().duplex_provenance,
                replay.final_simulation().duplex_provenance
            );
            if expected_gene == 0 && seed < 16 {
                let mut bad = outcome.trace.choices().to_vec();
                for choice in &mut bad {
                    if choice.address == "sample_allele.d"
                        || choice.address == "sample_allele_in_slot.d"
                    {
                        choice.value = ChoiceValue::AlleleId(1);
                    }
                }
                assert!(compiled
                    .replay_from_trace_records(&bad, seed, policy)
                    .is_err());
            }
        }
        assert!(counts.iter().all(|count| *count > 0));
        assert!(differs_from_filtered_dosage > 64);
    }
}
