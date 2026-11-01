//! The record validator is the oracle. This module proves it rejects.
//!
//! Every other module here asks "is the metadata true?". This one asks the
//! question behind that one: when the metadata is *false*, does anything
//! notice? A validator nobody has watched reject a wrong value is
//! decoration, and a test that only asserts a clean record stays quiet
//! cannot tell a working check from a missing one.
//!
//! The gap this closes was concrete. `validate.rs` defines 47 issue kinds,
//! and the C2 counter family -- `NMutationsMismatch`, `NPcrErrorsMismatch`,
//! `NQualityErrorsMismatch`, `NAmbiguousBasesMismatch` -- re-derives its
//! expected value from the event ledger specifically so a tampered or
//! hand-built record surfaces. Not one of those four had a test that
//! tampered with anything.
//!
//! Three issue kinds *were* already proven, in-crate and by the same means:
//! `validator_flags_tampered_n_cdr1_mutations`,
//! `validator_flags_p_length_mismatch_on_tampered_field` and
//! `validator_flags_a_tampered_observable_mutation_count`, all in
//! `src/airr_record/tests/projection.rs`. This module follows that pattern
//! and extends it to the counters it had not reached. It differs in one way
//! on purpose: those tests inject a hand-built fixture outcome with a known
//! event, while this one drives the real compiled passes, so it also
//! exercises what those passes actually produce -- which is how the
//! realized-versus-requested finding recorded below surfaced.
//!
//! What cannot serve as this proof is the Python side.
//! `test_validator_flags_tampered_n_v_mutations` does not tamper despite
//! its name: `validate_records` re-projects the record from the outcome, so
//! a record-dict edit is discarded before the validator sees it. Tampering
//! is only possible against `validate_airr_record` directly.
//!
//! The counters are checked as a family rather than one at a time because
//! they share a failure mode: each is a count the builder derives from the
//! ledger, and a check that silently stopped firing would leave the field
//! exported, plausible and unguarded.

use super::common;

use genairr_engine::airr_record::{
    build_airr_record, validate_airr_record, AirrRecord, RecordValidationIssue,
};
use genairr_engine::compiled::{CompiledSimulator, ExecutionPolicy};
use genairr_engine::dist::{AllelePoolDist, EmpiricalLengthDist, UniformBase};
use genairr_engine::ir::Segment;
use genairr_engine::pass::{Outcome, PassPlan};
use genairr_engine::passes::{
    AssembleSegmentPass, GenerateNPPass, NCorruptionPass, PCRErrorPass, QualityErrorPass,
    SampleAllelePass, UniformMutationPass,
};
use genairr_engine::refdata::RefDataConfig;

// ──────────────────────────────────────────────────────────────────
// Local helpers
// ──────────────────────────────────────────────────────────────────

/// A plan that gives every counter under test a non-zero truth. No indels
/// and no rev-comp -- this module is about the counters, and an indel moves
/// coordinates that other issue kinds would then report on, which would
/// blur what the assertions below are pinning.
///
/// The counts are 4 rather than 1 for a reason worth recording: a requested
/// count is not a realized count. Every substitution pass draws its
/// replacement base from all four bases *including the base already there*,
/// and an identity write is correctly not counted -- so `count = 1` yields
/// a counter of 0 about a quarter of the time. The first version of this
/// plan asked for one PCR error, got an identity write on seed 17, and
/// `n_pcr_errors` was 0. Measured over 300 human IGH records, the realized
/// fraction is 0.757 for PCR errors, 0.757 for quality errors and 0.78 for
/// uniform SHM; `ambiguous_base_calls` is ~1.0 because it pins the write to
/// `N` rather than sampling. Asking for 4 keeps every counter reliably
/// above zero without depending on which seed this runs under.
fn counted_plan(cfg: &RefDataConfig) -> PassPlan {
    let v_id = common::allele_id_by_name(cfg, Segment::V, "v01*01");
    let d_id = common::allele_id_by_name(cfg, Segment::D, "d01*01");
    let j_id = common::allele_id_by_name(cfg, Segment::J, "j01*01");

    let mut plan = PassPlan::new();

    plan.push(Box::new(SampleAllelePass::new(
        Segment::V,
        Box::new(AllelePoolDist::restricted_uniform(&cfg.v_pool, vec![v_id])),
    )));
    plan.push(Box::new(SampleAllelePass::new(
        Segment::D,
        Box::new(AllelePoolDist::restricted_uniform(&cfg.d_pool, vec![d_id])),
    )));
    plan.push(Box::new(SampleAllelePass::new(
        Segment::J,
        Box::new(AllelePoolDist::restricted_uniform(&cfg.j_pool, vec![j_id])),
    )));

    plan.push(Box::new(AssembleSegmentPass::new(Segment::V)));
    plan.push(Box::new(GenerateNPPass::new(
        Segment::Np1,
        Box::new(EmpiricalLengthDist::from_pairs(vec![(3, 1.0)])),
        Box::new(UniformBase),
    )));
    plan.push(Box::new(AssembleSegmentPass::new(Segment::D)));
    plan.push(Box::new(GenerateNPPass::new(
        Segment::Np2,
        Box::new(EmpiricalLengthDist::from_pairs(vec![(3, 1.0)])),
        Box::new(UniformBase),
    )));
    plan.push(Box::new(AssembleSegmentPass::new(Segment::J)));

    plan.push(Box::new(UniformMutationPass::new(
        Box::new(EmpiricalLengthDist::from_pairs(vec![(4, 1.0)])),
        Box::new(UniformBase),
    )));
    plan.push(Box::new(PCRErrorPass::new(
        Box::new(EmpiricalLengthDist::from_pairs(vec![(4, 1.0)])),
        Box::new(UniformBase),
    )));
    plan.push(Box::new(QualityErrorPass::new(
        Box::new(EmpiricalLengthDist::from_pairs(vec![(4, 1.0)])),
        Box::new(UniformBase),
    )));
    plan.push(Box::new(NCorruptionPass::new(Box::new(
        EmpiricalLengthDist::from_pairs(vec![(4, 1.0)]),
    ))));

    plan
}

fn run_counted(cfg: &RefDataConfig, seed: u64) -> (Outcome, AirrRecord) {
    let plan = counted_plan(cfg);
    let compiled = CompiledSimulator::compile(&plan, Some(cfg), None, ExecutionPolicy::Permissive)
        .expect("counted plan compiles");
    let outcome = compiled.run_one(seed).expect("counted plan runs");
    let rec = build_airr_record(&outcome, cfg, "tamper");
    (outcome, rec)
}

/// Every C2 counter issue kind this module tampers with, as a predicate
/// over the issue enum. Written as a match rather than a string compare so
/// that renaming a variant is a compile error here, not a test that
/// quietly stops asserting anything.
fn is_counter_mismatch(issue: &RecordValidationIssue) -> bool {
    matches!(
        issue,
        RecordValidationIssue::NMutationsMismatch { .. }
            | RecordValidationIssue::NPcrErrorsMismatch { .. }
            | RecordValidationIssue::NQualityErrorsMismatch { .. }
            | RecordValidationIssue::NAmbiguousBasesMismatch { .. }
    )
}

// ──────────────────────────────────────────────────────────────────
// Tests
// ──────────────────────────────────────────────────────────────────

#[test]
fn an_honest_record_trips_no_counter_check() {
    // The floor under every tampering test below: if a clean record
    // already tripped these checks, a firing check would prove nothing
    // about the tamper.
    let cfg = common::vdj_ambiguous_refdata();
    let (outcome, rec) = run_counted(&cfg, 17);

    let issues = validate_airr_record(&rec, &outcome, &cfg);
    let counter_issues: Vec<_> = issues.iter().filter(|i| is_counter_mismatch(i)).collect();
    assert!(
        counter_issues.is_empty(),
        "clean record tripped counter checks: {:?}",
        counter_issues,
    );
}

#[test]
fn the_counters_under_test_are_actually_exercised() {
    // A tamper against a truth of 0 is a weaker test than one against a
    // real count: it cannot catch a check that compares only to zero. Pin
    // that the plan above really does produce each kind of base change.
    let cfg = common::vdj_ambiguous_refdata();
    let (_outcome, rec) = run_counted(&cfg, 17);

    for (field, value) in [
        ("n_mutations", rec.n_mutations),
        ("n_pcr_errors", rec.n_pcr_errors),
        ("n_quality_errors", rec.n_quality_errors),
        ("n_ambiguous_bases", rec.n_ambiguous_bases),
    ] {
        assert!(
            value > 0,
            "{} is {} under the counted plan, so tampering with it \
             would only test the zero case",
            field,
            value,
        );
    }
}

#[test]
fn a_tampered_ambiguous_base_count_is_rejected() {
    // The check this module was written for. `n_ambiguous_bases` is
    // derived from the `corrupt.ns` ledger, and the validator re-walks
    // that ledger itself -- asserting more than the builder does, namely
    // that the pass wrote `N` and nothing else.
    let cfg = common::vdj_ambiguous_refdata();
    let (outcome, mut rec) = run_counted(&cfg, 17);
    let truth = rec.n_ambiguous_bases;

    rec.n_ambiguous_bases = truth + 7;

    let issues = validate_airr_record(&rec, &outcome, &cfg);
    let found = issues.iter().find_map(|i| match i {
        RecordValidationIssue::NAmbiguousBasesMismatch {
            reported,
            event_count,
        } => Some((*reported, *event_count)),
        _ => None,
    });
    assert_eq!(
        found,
        Some((truth + 7, truth)),
        "validator did not reject n_ambiguous_bases={} against a ledger \
         of {}; issues were {:?}",
        truth + 7,
        truth,
        issues,
    );
}

#[test]
fn a_tampered_ambiguous_base_count_is_rejected_downward_too() {
    // Under-reporting is the direction a bug takes: the duplex path
    // dropped this counter entirely and reported 0. An over-count test
    // alone would pass against a check that only looked for excess.
    let cfg = common::vdj_ambiguous_refdata();
    let (outcome, mut rec) = run_counted(&cfg, 17);
    let truth = rec.n_ambiguous_bases;

    rec.n_ambiguous_bases = 0;

    let issues = validate_airr_record(&rec, &outcome, &cfg);
    assert!(
        issues.iter().any(|i| matches!(
            i,
            RecordValidationIssue::NAmbiguousBasesMismatch {
                reported: 0,
                event_count,
            } if *event_count == truth
        )),
        "validator accepted n_ambiguous_bases=0 against a ledger of {}; \
         issues were {:?}",
        truth,
        issues,
    );
}

#[test]
fn a_tampered_mutation_count_is_rejected() {
    let cfg = common::vdj_ambiguous_refdata();
    let (outcome, mut rec) = run_counted(&cfg, 17);
    rec.n_mutations += 7;

    let issues = validate_airr_record(&rec, &outcome, &cfg);
    assert!(
        issues
            .iter()
            .any(|i| matches!(i, RecordValidationIssue::NMutationsMismatch { .. })),
        "validator accepted an inflated n_mutations; issues were {:?}",
        issues,
    );
}

#[test]
fn a_tampered_pcr_error_count_is_rejected() {
    let cfg = common::vdj_ambiguous_refdata();
    let (outcome, mut rec) = run_counted(&cfg, 17);
    rec.n_pcr_errors += 7;

    let issues = validate_airr_record(&rec, &outcome, &cfg);
    assert!(
        issues
            .iter()
            .any(|i| matches!(i, RecordValidationIssue::NPcrErrorsMismatch { .. })),
        "validator accepted an inflated n_pcr_errors; issues were {:?}",
        issues,
    );
}

#[test]
fn a_tampered_quality_error_count_is_rejected() {
    let cfg = common::vdj_ambiguous_refdata();
    let (outcome, mut rec) = run_counted(&cfg, 17);
    rec.n_quality_errors += 7;

    let issues = validate_airr_record(&rec, &outcome, &cfg);
    assert!(
        issues
            .iter()
            .any(|i| matches!(i, RecordValidationIssue::NQualityErrorsMismatch { .. })),
        "validator accepted an inflated n_quality_errors; issues were {:?}",
        issues,
    );
}

#[test]
fn each_counter_check_rejects_only_its_own_field() {
    // Four checks that all compare a count to a ledger walk are easy to
    // cross-wire. Tamper with one field at a time and assert exactly one
    // counter kind fires -- which is also what makes the issue useful to
    // whoever reads the report.
    let cfg = common::vdj_ambiguous_refdata();

    for field in ["ambiguous", "mutations", "pcr", "quality"] {
        let (outcome, mut rec) = run_counted(&cfg, 17);
        match field {
            "ambiguous" => rec.n_ambiguous_bases += 7,
            "mutations" => rec.n_mutations += 7,
            "pcr" => rec.n_pcr_errors += 7,
            "quality" => rec.n_quality_errors += 7,
            _ => unreachable!(),
        }

        let fired: Vec<_> = validate_airr_record(&rec, &outcome, &cfg)
            .into_iter()
            .filter(is_counter_mismatch)
            .collect();
        assert_eq!(
            fired.len(),
            1,
            "tampering with {} fired {} counter checks, expected exactly \
             one: {:?}",
            field,
            fired.len(),
            fired,
        );
    }
}

#[test]
fn an_honest_record_trips_no_np_partition_check() {
    // The floor under the two below, and a real assertion in its own right:
    // the np origin partition is re-derived on every record, so a clean
    // record tripping it would mean the walk and the builder disagree about
    // which bucket owns a base.
    let cfg = common::vdj_ambiguous_refdata();
    let (outcome, rec) = run_counted(&cfg, 17);

    let np_issues: Vec<_> = validate_airr_record(&rec, &outcome, &cfg)
        .into_iter()
        .filter(|i| {
            matches!(
                i,
                RecordValidationIssue::NpPartitionMismatch { .. }
                    | RecordValidationIssue::NpPartitionOriginAmbiguous { .. }
            )
        })
        .collect();
    assert!(
        np_issues.is_empty(),
        "clean record tripped the np partition check: {:?}",
        np_issues,
    );
}

#[test]
fn a_tampered_np_origin_count_is_rejected() {
    // The np interval is partitioned by what wrote each base, and the
    // validator walks the pool flags itself rather than calling the helper
    // the builder uses -- a helper shared with the thing it checks can only
    // catch a tampered record, never a wrong projection.
    let cfg = common::vdj_ambiguous_refdata();
    let (outcome, mut rec) = run_counted(&cfg, 17);
    let truth = rec.np1_n_tdt;

    rec.np1_n_tdt = truth + 3;

    let issues = validate_airr_record(&rec, &outcome, &cfg);
    let found = issues.iter().find_map(|i| match i {
        RecordValidationIssue::NpPartitionMismatch {
            region: "np1",
            bucket: "tdt",
            reported,
            flag_count,
        } => Some((*reported, *flag_count)),
        _ => None,
    });
    assert_eq!(
        found,
        Some((truth + 3, truth)),
        "validator did not reject np1_n_tdt={} against {} flagged positions; \
         issues were {:?}",
        truth + 3,
        truth,
        issues,
    );
}

#[test]
fn moving_a_base_between_np_buckets_is_rejected_though_the_sum_still_holds() {
    // The case the sum identity cannot see, and the reason the validator
    // compares each bucket rather than only their total. Take one base from
    // the TdT count and give it to the P count: the three still add to
    // np1_length, and every base is still accounted for, but two buckets now
    // lie about which pass wrote it.
    //
    // Demonstrated rather than assumed: replacing the builder's projection
    // with "every flagged base is TdT" left
    // test_the_np_partition_sums_to_the_np_length passing on both of its
    // parameters.
    let cfg = common::vdj_ambiguous_refdata();
    let (outcome, mut rec) = run_counted(&cfg, 17);
    assert!(
        rec.np1_n_tdt > 0,
        "np1 carries no TdT addition under the counted plan, so there is \
         nothing to move out of that bucket",
    );
    let sum_before = rec.np1_n_tdt + rec.np1_n_p + rec.np1_n_inserted;

    rec.np1_n_tdt -= 1;
    rec.np1_n_p += 1;
    assert_eq!(
        rec.np1_n_tdt + rec.np1_n_p + rec.np1_n_inserted,
        sum_before,
        "the tamper was supposed to preserve the sum",
    );

    let fired: Vec<_> = validate_airr_record(&rec, &outcome, &cfg)
        .into_iter()
        .filter(|i| matches!(i, RecordValidationIssue::NpPartitionMismatch { .. }))
        .collect();
    assert_eq!(
        fired.len(),
        2,
        "moving a base between buckets should be reported against both \
         buckets, got {:?}",
        fired,
    );
}
