use super::*;
use crate::airr_record::{build_airr_record, validate_airr_record, RecordValidationIssue};
use crate::dist::{AllelePoolDist, EmpiricalLengthDist, UniformBase};
use crate::pass::{testing::PassRuntime, PassPlan};
use crate::passes::{
    AssembleSegmentPass, EndLossPass, GenerateNPPass, LossEnd, SampleAllelePass,
    UniformMutationPass,
};

fn fixture() -> (RefDataConfig, PassPlan) {
    let mut cfg = RefDataConfig::empty(ChainType::Vj);
    for (segment, sequence, anchor) in [
        (Segment::V, b"AAACCCTGT".as_slice(), 6),
        (Segment::J, b"TTTGGG".as_slice(), 0),
    ] {
        let allele = Allele {
            name: format!("{segment:?}*01"),
            gene: format!("{segment:?}"),
            segment,
            seq: sequence.to_vec(),
            anchor: Some(anchor),
            functional_status: None,
            subregions: Vec::new(),
        };
        match segment {
            Segment::V => {
                let _ = cfg.v_pool.push(allele);
            }
            Segment::J => {
                let _ = cfg.j_pool.push(allele);
            }
            _ => unreachable!(),
        }
    }
    let mut plan = PassPlan::new();
    for (segment, pool) in [(Segment::V, &cfg.v_pool), (Segment::J, &cfg.j_pool)] {
        plan.push(Box::new(SampleAllelePass::new(
            segment,
            Box::new(AllelePoolDist::uniform(pool)),
        )));
    }
    plan.push(Box::new(AssembleSegmentPass::new(Segment::V)));
    plan.push(Box::new(GenerateNPPass::new(
        Segment::Np1,
        Box::new(EmpiricalLengthDist::from_pairs(vec![(3, 1.0)])),
        Box::new(UniformBase),
    )));
    plan.push(Box::new(AssembleSegmentPass::new(Segment::J)));
    (cfg, plan)
}

fn mutate(plan: &mut PassPlan, count: i64) {
    plan.push(Box::new(UniformMutationPass::new(
        Box::new(EmpiricalLengthDist::from_pairs(vec![(count, 1.0)])),
        Box::new(UniformBase),
    )));
}

fn lose(plan: &mut PassPlan, count: i64) {
    plan.push(Box::new(EndLossPass::new(
        LossEnd::Three,
        Box::new(EmpiricalLengthDist::from_pairs(vec![(count, 1.0)])),
    )));
}

#[test]
fn mutation_burden_preserves_pre_observation_denominator_even_for_empty_reads() {
    for loss in [0, 6, 99] {
        let (cfg, mut plan) = fixture();
        mutate(&mut plan, 3);
        lose(&mut plan, loss);
        let outcome = PassRuntime::execute_with_refdata(&plan, Simulation::new(), 91, &cfg);
        let rec = build_airr_record(&outcome, &cfg, "loss");
        assert_eq!(rec.n_mutations, 3);
        assert_eq!(rec.sequence_length, 18 - loss.min(18));
        assert_eq!(rec.mutation_rate, 3.0 / 18.0);
        assert!(!validate_airr_record(&rec, &outcome, &cfg)
            .iter()
            .any(|issue| matches!(issue, RecordValidationIssue::MutationRateMismatch { .. })));
    }
}

#[test]
fn mutation_burden_sums_stage_specific_lengths_and_ignores_zero_count() {
    let (cfg, mut plan) = fixture();
    mutate(&mut plan, 2);
    lose(&mut plan, 6);
    mutate(&mut plan, 1);
    mutate(&mut plan, 0);
    let outcome = PassRuntime::execute_with_refdata(&plan, Simulation::new(), 91, &cfg);
    let mut rec = build_airr_record(&outcome, &cfg, "stages");
    assert_eq!(rec.n_mutations, 3);
    assert_eq!(rec.mutation_rate, 2.0 / 18.0 + 1.0 / 12.0);
    assert!(!validate_airr_record(&rec, &outcome, &cfg)
        .iter()
        .any(|issue| matches!(issue, RecordValidationIssue::MutationRateMismatch { .. })));
    rec.mutation_rate = 3.0 / 12.0;
    assert!(validate_airr_record(&rec, &outcome, &cfg)
        .iter()
        .any(|issue| matches!(issue, RecordValidationIssue::MutationRateMismatch { .. })));
}

#[test]
fn continuation_preserves_ancestor_burden_counters_and_rebases_trace_spans() {
    let (cfg, mut plan) = fixture();
    mutate(&mut plan, 2);
    let parent = PassRuntime::execute_with_refdata(&plan, Simulation::new(), 90, &cfg);
    let mut child_plan = PassPlan::new();
    lose(&mut child_plan, 6);
    mutate(&mut child_plan, 1);
    let child =
        PassRuntime::execute_with_refdata(&child_plan, parent.final_simulation().clone(), 91, &cfg);
    let merged = parent.with_continuation(&child).unwrap();
    let rec = build_airr_record(&merged, &cfg, "child");
    assert_eq!(rec.n_mutations, 3);
    assert_eq!(
        rec.n_v_mutations + rec.n_d_mutations + rec.n_j_mutations + rec.n_np_mutations,
        3
    );
    assert_eq!(rec.mutation_rate, 2.0 / 18.0 + 1.0 / 12.0);
    for (index, event) in merged.events.iter().enumerate() {
        assert_eq!(event.pass_index, index);
        assert_eq!(event.pass_name, merged.pass_names[index]);
        assert!(event.trace_span.end <= merged.trace.len());
        if index > 0 {
            assert_eq!(
                event.trace_span.start,
                merged.events[index - 1].trace_span.end
            );
        }
    }
    let issues = validate_airr_record(&rec, &merged, &cfg);
    assert!(issues.is_empty(), "{issues:?}");
    // This direct PassRuntime continuation did not refresh V inference.
    // Generating truth remains available without inventing a re-call.
    assert_eq!(rec.v_call, "V*01");
    assert!(rec.recalled_v_call.is_empty());
    let mut invented_recall = rec.clone();
    invented_recall.recalled_v_call = rec.v_call.clone();
    assert!(validate_airr_record(&invented_recall, &merged, &cfg)
        .iter()
        .any(|issue| matches!(
            issue,
            RecordValidationIssue::RecalledAlleleCallMismatch {
                segment: Segment::V,
                ..
            }
        )));
    let mut grand_plan = PassPlan::new();
    mutate(&mut grand_plan, 1);
    let grand =
        PassRuntime::execute_with_refdata(&grand_plan, merged.final_simulation().clone(), 92, &cfg);
    let grand = merged.with_continuation(&grand).unwrap();
    let rec = build_airr_record(&grand, &cfg, "grand");
    assert_eq!(rec.n_mutations, 4);
    assert!((rec.mutation_rate - (2.0 / 18.0 + 2.0 / 12.0)).abs() < 1e-12);
    assert!(validate_airr_record(&rec, &grand, &cfg).is_empty());
    let empty = PassRuntime::execute_with_refdata(
        &PassPlan::new(),
        grand.final_simulation().clone(),
        93,
        &cfg,
    );
    let unchanged = grand.with_continuation(&empty).unwrap();
    assert_eq!(unchanged.events, grand.events);
    assert_eq!(unchanged.revisions.len(), grand.revisions.len());
    assert!(parent.with_continuation(&empty).is_err());
}

#[test]
fn zero_mutations_have_zero_burden_on_nonempty_and_empty_molecules() {
    for length in [0, 3] {
        let mut initial = Simulation::new();
        for _ in 0..length {
            initial = initial
                .with_nucleotide_pushed(Nucleotide::synthetic(
                    b'A',
                    Segment::Np1,
                    crate::ir::NucFlags::default(),
                ))
                .0;
        }
        let mut plan = PassPlan::new();
        mutate(&mut plan, 0);
        let outcome = PassRuntime::execute(&plan, initial, 0);
        assert_eq!(outcome.final_simulation().mutation_count, 0);
        assert_eq!(outcome.final_simulation().mutation_burden, 0.0);
    }
}
