use super::*;
use crate::passes::corrupt::indel::edit_view::EditIndex;
use crate::passes::corrupt::indel::event::IndelEvent;
use crate::passes::corrupt::indel::tuple_sampler::apply_event_hypothetical;
use crate::productivity::SequenceView;

fn events(length: usize) -> Vec<IndelEvent> {
    let mut result = Vec::new();
    for site in 0..=length {
        for &base in b"ACGT" {
            result.push(IndelEvent::Insertion {
                site: site as u32,
                base,
            });
        }
        if site < length {
            result.push(IndelEvent::Deletion {
                site: Some(site as u32),
            });
        }
    }
    if length == 0 {
        result.push(IndelEvent::Deletion { site: None });
    }
    result
}

#[test]
fn two_edit_view_matches_materialized_ir_and_codon_evidence_exhaustively() {
    let (cfg, original) = balanced_dp_distribution_fixture();
    // Ambiguous coding bases exercise definite-stop versus possible-stop summaries.
    for ambiguous in [false, true] {
        let initial = if ambiguous {
            original.with_base_changed(NucHandle::new(4), b'N')
        } else {
            original.clone()
        };
        let index = EditIndex::new(&initial);
        for first in events(initial.pool.len()) {
            let one = apply_event_hypothetical(&initial, first);
            let one_view = index.view().edited(first);
            for second in events(one.pool.len()) {
                let materialized = apply_event_hypothetical(&one, second);
                let view = one_view.edited(second);
                assert_eq!(view.len(), materialized.pool.len());
                for p in 0..view.len() {
                    assert_eq!(
                        view.nucleotide(p),
                        materialized.nucleotide(p),
                        "{first:?}, {second:?}, {p}"
                    );
                }
                for segment in [Segment::V, Segment::J] {
                    assert_eq!(view.region(segment), materialized.region(segment));
                    assert_eq!(
                        view.last_germline(segment),
                        materialized.last_germline(segment)
                    );
                    for g in 0..6 {
                        assert_eq!(
                            view.germline_position(segment, g),
                            materialized.germline_position(segment, g)
                        );
                    }
                }
                for phase in 0..3 {
                    assert_eq!(
                        view.coding_stop(phase, view.len()),
                        materialized.coding_stop(phase, view.len()),
                        "{first:?}, {second:?}, phase={phase}"
                    );
                }
                let (seen, built) = (
                    crate::productivity::assess(&view, &cfg),
                    crate::productivity::assess(&materialized, &cfg),
                );
                assert_eq!(seen.productive, built.productive);
                // The junction, including its end after an edit inside the J codon.
                assert_eq!(seen.junction, built.junction, "{first:?}, {second:?}");
                for segment in [Segment::V, Segment::J] {
                    for g in 0..6 {
                        assert_eq!(
                            view.anchor_codon_end(segment, g),
                            materialized.anchor_codon_end(segment, g)
                        );
                    }
                }
            }
        }
    }
}

#[test]
fn extended_edit_view_matches_ir_through_sixty_four_ordered_edits() {
    let (cfg, original) = balanced_dp_distribution_fixture();
    let initial = original.with_base_changed(NucHandle::new(4), b'N');
    let index = EditIndex::new(&initial);
    for seed in 0..32 {
        let mut rng = crate::rng::Rng::new(seed);
        let mut view = index.view_with_capacity::<129>();
        let mut materialized = initial.clone();
        for _ in 0..64 {
            let options = events(materialized.pool.len());
            let event = options[rng.range_u32(options.len() as u32) as usize];
            view = view.edited(event);
            materialized = apply_event_hypothetical(&materialized, event);
            assert_eq!(view.len(), materialized.pool.len());
            for p in 0..view.len() {
                assert_eq!(view.nucleotide(p), materialized.nucleotide(p));
            }
            for segment in [Segment::V, Segment::J] {
                assert_eq!(view.region(segment), materialized.region(segment));
                assert_eq!(
                    view.last_germline(segment),
                    materialized.last_germline(segment)
                );
                for g in 0..6 {
                    assert_eq!(
                        view.germline_position(segment, g),
                        materialized.germline_position(segment, g)
                    );
                }
            }
            for phase in 0..3 {
                assert_eq!(
                    view.coding_stop(phase, view.len()),
                    materialized.coding_stop(phase, view.len())
                );
            }
            let (seen, built) = (
                crate::productivity::assess(&view, &cfg),
                crate::productivity::assess(&materialized, &cfg),
            );
            assert_eq!(seen.productive, built.productive);
            assert_eq!(seen.junction, built.junction);
        }
    }
}

#[test]
fn natural_sampler_can_cancel_new_sites_and_temporarily_interrupted_anchors() {
    let (cfg, sim) = make_substitution_productive_vj_fixture();
    let contracts = productive();
    let mut seen_junction_cancellation = false;
    let mut seen_anchor_repair = false;
    let mut plan = PassPlan::new();
    plan.push(Box::new(IndelPass::new(
        fixed_count(2),
        0.5,
        Box::new(UniformBase),
    )));
    for seed in 0..256 {
        let out = PassRuntime::execute_strict_with_context(
            &plan,
            sim.clone(),
            seed,
            Some(&cfg),
            Some(&contracts),
        )
        .unwrap();
        let value = |name| &out.trace.find(name).unwrap().value;
        if value("corrupt.indel.kind[0]") == &ChoiceValue::Bool(true)
            && value("corrupt.indel.kind[1]") == &ChoiceValue::Bool(false)
            && value("corrupt.indel.site[0]") == value("corrupt.indel.site[1]")
        {
            if value("corrupt.indel.site[0]") == &ChoiceValue::Int(3) {
                seen_junction_cancellation = true;
            }
            if matches!(
                value("corrupt.indel.site[0]"),
                ChoiceValue::Int(1 | 2 | 4 | 5)
            ) {
                seen_anchor_repair = true;
            }
        }
    }
    assert!(
        seen_junction_cancellation,
        "newly created deletion support must be reachable"
    );
    assert!(
        seen_anchor_repair,
        "temporary interruption is not irreversible loss"
    );
}

#[test]
fn known_nonempty_support_with_zero_work_budget_is_not_reported_empty() {
    let (cfg, sim) = balanced_dp_distribution_fixture();
    let contracts = productive();
    for strict in [false, true] {
        let pass = IndelPass::new(fixed_count(2), 0.5, Box::new(UniformBase))
            .with_sampling_options(IndelSamplingOptions {
                rejection_attempts: 0,
                exact_node_budget: 0,
                ..Default::default()
            });
        let mut trace = crate::trace::Trace::new();
        let mut rng = crate::rng::Rng::new(5);
        let mut ctx = crate::pass::PassContext {
            trace: &mut trace,
            rng: &mut rng,
            pass_index: 0,
            refdata: Some(&cfg),
            contracts: Some(&contracts),
            feasibility: None,
            reference_index: None,
            replay_cursor: None,
            event_log_sink: None,
        };
        let error = pass
            .execute_with_policy(&sim, &mut ctx, strict)
            .unwrap_err();
        assert!(matches!(
            error,
            PassError::SamplingBudgetExceeded {
                requested_count: 2,
                explored_states: 0,
                ..
            }
        ));
        assert_eq!(sim.pool.len(), 12);
        assert!(!trace
            .choices()
            .iter()
            .any(|r| r.address.starts_with("corrupt.indel.site")));
    }
}

#[test]
fn exact_empty_support_has_a_completion_diagnostic_in_permissive_mode() {
    let (cfg, sim) = make_substitution_productive_vj_fixture();
    let contracts = productive();
    let mut plan = PassPlan::new();
    plan.push(Box::new(IndelPass::new(
        fixed_count(1),
        0.0,
        Box::new(UniformBase),
    )));
    let out = PassRuntime::execute_with_context(&plan, sim, 1, Some(&cfg), Some(&contracts));
    let report = &out.trace.diagnostics()[0];
    assert_eq!(report.status, "proven_empty_support");
    assert_eq!(report.applied_count, 0);
    assert_eq!(report.requested_count, 1);
    assert_eq!(report.model, "natural_sequential_v1");
    assert!(report.explored_states > 0);
}

#[test]
fn higher_count_exhaustion_never_returns_a_partial_exact_reservoir() {
    let (cfg, sim) = balanced_dp_distribution_fixture();
    let contracts = productive();
    let pass = IndelPass::new(fixed_count(3), 0.5, Box::new(UniformBase)).with_sampling_options(
        IndelSamplingOptions {
            rejection_attempts: 0,
            exact_node_budget: 500,
            ..Default::default()
        },
    );
    let mut trace = crate::trace::Trace::new();
    let mut rng = crate::rng::Rng::new(2);
    let mut ctx = crate::pass::PassContext {
        trace: &mut trace,
        rng: &mut rng,
        pass_index: 0,
        refdata: Some(&cfg),
        contracts: Some(&contracts),
        feasibility: None,
        reference_index: None,
        replay_cursor: None,
        event_log_sink: None,
    };
    assert!(matches!(
        pass.execute_checked(&sim, &mut ctx),
        Err(PassError::SamplingBudgetExceeded {
            explored_states: 500,
            ..
        })
    ));
    assert!(!trace
        .choices()
        .iter()
        .any(|r| r.address.starts_with("corrupt.indel.site")));
}

#[test]
fn sampler_version_and_work_limits_change_plan_identity() {
    let standard = IndelPass::new(fixed_count(2), 0.5, Box::new(UniformBase));
    let limited = IndelPass::new(fixed_count(2), 0.5, Box::new(UniformBase)).with_sampling_options(
        IndelSamplingOptions {
            exact_node_budget: 1,
            ..Default::default()
        },
    );
    assert!(standard
        .parameter_signature()
        .contains("sampler=natural_sequential_v1"));
    assert_ne!(
        standard.parameter_signature(),
        limited.parameter_signature()
    );
}

#[test]
fn finite_extreme_weights_preserve_positive_support_and_rejection_distribution() {
    use crate::dist::CategoricalBase;
    let tiny = IndelPass::new(
        fixed_count(1),
        1.0,
        Box::new(CategoricalBase::from_weights([1e-300, 1e300, 0.0, 0.0])),
    );
    let support = tiny.base_log_support().unwrap();
    assert_eq!(support.len(), 2);
    assert!(support.iter().all(|(_, w)| w.is_finite()));
    assert!(
        support[0].1 < -1000.0,
        "positive support must survive probability underflow"
    );

    let (cfg, sim) = make_substitution_productive_vj_fixture();
    let mut plan = PassPlan::new();
    plan.push(Box::new(IndelPass::new(
        fixed_count(1),
        1.0,
        Box::new(CategoricalBase::from_weights([1e308, 1e308, 1e308, 1e308])),
    )));
    let mut counts = [0u32; 4];
    for seed in 0..2048 {
        let out = PassRuntime::execute_strict_with_context(
            &plan,
            sim.clone(),
            seed,
            Some(&cfg),
            Some(&productive()),
        )
        .unwrap();
        let ChoiceValue::Base(base) = out.trace.find("corrupt.indel.base[0]").unwrap().value else {
            panic!()
        };
        counts[b"ACGT".iter().position(|b| *b == base).unwrap()] += 1;
    }
    for count in counts {
        let sigma = (2048.0_f64 * 0.25 * 0.75).sqrt();
        assert!((f64::from(count) - 512.0).abs() < 5.0 * sigma);
    }
}

#[test]
fn exact_fallback_and_rejection_preserve_nonuniform_conditional_weights() {
    use crate::dist::CategoricalBase;
    let (cfg, sim) = make_substitution_productive_vj_fixture();
    // This molecule admits only insertion after J, independent of its base.
    // Therefore the accepted A:C marginal must remain exactly 1:3.
    for rejection_attempts in [0, 256] {
        let mut plan = PassPlan::new();
        plan.push(Box::new(
            IndelPass::new(
                fixed_count(1),
                1.0,
                Box::new(CategoricalBase::from_weights([1.0, 3.0, 0.0, 0.0])),
            )
            .with_sampling_options(IndelSamplingOptions {
                rejection_attempts,
                ..Default::default()
            }),
        ));
        let mut adenines = 0;
        for seed in 0..2048 {
            let out = PassRuntime::execute_strict_with_context(
                &plan,
                sim.clone(),
                seed,
                Some(&cfg),
                Some(&productive()),
            )
            .unwrap();
            assert_eq!(
                out.trace.find("corrupt.indel.site[0]").unwrap().value,
                ChoiceValue::Int(6)
            );
            match &out.trace.find("corrupt.indel.base[0]").unwrap().value {
                ChoiceValue::Base(b'A') => adenines += 1,
                ChoiceValue::Base(b'C') => (),
                other => panic!("unsupported base: {other:?}"),
            }
        }
        let sigma = (2048.0_f64 * 0.25 * 0.75).sqrt();
        assert!(
            (f64::from(adenines) - 512.0).abs() < 5.0 * sigma,
            "attempts={rejection_attempts}, A={adenines}"
        );
    }
}

#[test]
#[ignore = "release timing probe for the exact two-edit operating range"]
fn benchmark_exact_two_edits_at_receptor_lengths() {
    use crate::assignment::AlleleInstance;
    use crate::refdata::{Allele, ChainType, RefDataConfig};
    for length in [369usize, 501, 651] {
        let mut cfg = RefDataConfig::empty(ChainType::Vj);
        let mut sim = Simulation::new();
        let mut v = vec![b'A'; length - 6];
        let anchor = v.len() - 3;
        v[anchor..].copy_from_slice(b"TGT");
        for (segment, seq, anchor) in [(Segment::V, v, anchor), (Segment::J, b"TGGAAA".to_vec(), 0)]
        {
            let allele = Allele {
                name: format!("{segment:?}*01"),
                gene: format!("{segment:?}"),
                seq: seq.clone(),
                segment,
                anchor: Some(anchor as u16),
                functional_status: None,
                subregions: Vec::new(),
            };
            let id = if segment == Segment::V {
                cfg.v_pool.push(allele)
            } else {
                cfg.j_pool.push(allele)
            };
            sim = sim.with_allele_assigned(segment, AlleleInstance::new(id));
            let start = sim.pool.len() as u32;
            for (g, b) in seq.iter().enumerate() {
                sim = sim
                    .with_nucleotide_pushed(Nucleotide::germline(*b, g as u16, segment))
                    .0;
            }
            sim = sim.with_region_added(Region::new(
                segment,
                NucHandle::new(start),
                NucHandle::new(start + seq.len() as u32),
            ));
        }
        let mut plan = PassPlan::new();
        plan.push(Box::new(
            IndelPass::new(fixed_count(2), 0.5, Box::new(UniformBase)).with_sampling_options(
                IndelSamplingOptions {
                    rejection_attempts: 0,
                    exact_node_budget: 20_000_000,
                    ..Default::default()
                },
            ),
        ));
        let start = std::time::Instant::now();
        let out = PassRuntime::execute_strict_with_context(
            &plan,
            sim,
            20261028,
            Some(&cfg),
            Some(&productive()),
        )
        .unwrap();
        eprintln!(
            "EXACT_BENCH length={length} elapsed_ms={} explored={}",
            start.elapsed().as_millis(),
            out.trace.diagnostics()[0].explored_states
        );
    }
}
