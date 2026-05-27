//! Trace-injected replay tests for `IndelPass` (Tier 3 final).
//!
//! These tests build a `PassContext` with a cursor in isolation
//! so the per-tuple validation chain (shape, bounds, frame,
//! post-event) can be exercised without depending on the
//! full-plan staged-migration replay surface.

use super::indel_test_sim;
use super::*;
use crate::address::ChoiceAddress;
use crate::pass::PassContext;
use crate::replay::{ReplayError, TraceCursor};
use crate::rng::Rng;
use crate::trace::{ChoiceRecord, Trace};

fn rec(addr: ChoiceAddress, v: ChoiceValue) -> ChoiceRecord {
    ChoiceRecord::new(addr.to_string(), v)
}

fn run_indel_replay(
    pass: &IndelPass,
    sim: Simulation,
    records: Vec<ChoiceRecord>,
    contracts: Option<&crate::contract::ContractSet>,
    refdata: Option<&crate::refdata::RefDataConfig>,
) -> (Result<Simulation, PassError>, Trace, u64) {
    let mut cursor = TraceCursor::from_owned(records);
    let mut trace = Trace::new();
    let mut rng = Rng::new(0xc0ff_ee);
    let result = {
        let mut ctx = PassContext {
            trace: &mut trace,
            rng: &mut rng,
            pass_index: 0,
            refdata,
            contracts,
            feasibility: None,
            reference_index: None,
            replay_cursor: Some(&mut cursor),
            event_log_sink: None,
        };
        pass.execute_checked(&sim, &mut ctx)
    };
    let words = rng.words_consumed();
    (result, trace, words)
}

fn count_rec(count: i64) -> ChoiceRecord {
    rec(ChoiceAddress::CorruptIndelCount, ChoiceValue::Int(count))
}

// ─────────────────────────────────────────────────────────────────
// 1. Valid insertion + deletion tuple replays with zero RNG
// ─────────────────────────────────────────────────────────────────

#[test]
fn indel_replay_consumes_insertion_with_zero_rng() {
    // Single insertion at site 5, base = 'A'. Pool was "AAACCCGGGTTT"
    // (12 bytes); after insertion the pool has 13 bytes with 'A' at
    // index 5.
    let pass = IndelPass::new(fixed_count(1), 1.0, Box::new(UniformBase));
    let records = vec![
        count_rec(1),
        rec(ChoiceAddress::CorruptIndelKind(0), ChoiceValue::Bool(true)),
        rec(ChoiceAddress::CorruptIndelSite(0), ChoiceValue::Int(5)),
        rec(ChoiceAddress::CorruptIndelBase(0), ChoiceValue::Base(b'A')),
    ];
    let (result, trace, rng_words) =
        run_indel_replay(&pass, indel_test_sim(), records, None, None);
    let next = result.unwrap();

    assert_eq!(next.pool.len(), 13);
    assert_eq!(next.pool.get(NucHandle::new(5)).unwrap().base, b'A');
    // Trace echoes the replayed values.
    assert_eq!(
        trace.find("corrupt.indel.kind[0]").unwrap().value,
        ChoiceValue::Bool(true),
    );
    assert_eq!(
        trace.find("corrupt.indel.site[0]").unwrap().value,
        ChoiceValue::Int(5),
    );
    assert_eq!(
        trace.find("corrupt.indel.base[0]").unwrap().value,
        ChoiceValue::Base(b'A'),
    );
    // No RNG consumed.
    assert_eq!(rng_words, 0);
}

#[test]
fn indel_replay_consumes_deletion_with_zero_rng() {
    // Single deletion at site 3. Pool was 12 bytes; after deletion 11.
    let pass = IndelPass::new(fixed_count(1), 0.0, Box::new(UniformBase));
    let initial = indel_test_sim();
    let bases_before: Vec<u8> = initial.pool.as_slice().iter().map(|n| n.base).collect();
    let records = vec![
        count_rec(1),
        rec(ChoiceAddress::CorruptIndelKind(0), ChoiceValue::Bool(false)),
        rec(ChoiceAddress::CorruptIndelSite(0), ChoiceValue::Int(3)),
    ];
    let (result, _trace, rng_words) =
        run_indel_replay(&pass, initial, records, None, None);
    let next = result.unwrap();

    assert_eq!(next.pool.len(), 11);
    // Bases shifted left at site 3 onward.
    for i in 0..3usize {
        assert_eq!(
            next.pool.get(NucHandle::new(i as u32)).unwrap().base,
            bases_before[i],
        );
    }
    for i in 3..11usize {
        assert_eq!(
            next.pool.get(NucHandle::new(i as u32)).unwrap().base,
            bases_before[i + 1],
        );
    }
    assert_eq!(rng_words, 0);
}

#[test]
fn indel_replay_consumes_mixed_insertion_and_deletion_tuple() {
    // 2-event tuple: deletion at site 5, then insertion at site 5
    // (against the shrunken pool). Net effect is a 1-byte swap at
    // position 5.
    let pass = IndelPass::new(fixed_count(2), 0.5, Box::new(UniformBase));
    let initial = indel_test_sim();
    let records = vec![
        count_rec(2),
        rec(ChoiceAddress::CorruptIndelKind(0), ChoiceValue::Bool(false)),
        rec(ChoiceAddress::CorruptIndelSite(0), ChoiceValue::Int(5)),
        rec(ChoiceAddress::CorruptIndelKind(1), ChoiceValue::Bool(true)),
        rec(ChoiceAddress::CorruptIndelSite(1), ChoiceValue::Int(5)),
        rec(ChoiceAddress::CorruptIndelBase(1), ChoiceValue::Base(b'T')),
    ];
    let (result, _, _) = run_indel_replay(&pass, initial, records, None, None);
    let next = result.unwrap();
    assert_eq!(next.pool.len(), 12); // deletion -1, insertion +1, net 0.
    assert_eq!(next.pool.get(NucHandle::new(5)).unwrap().base, b'T');
}

// ─────────────────────────────────────────────────────────────────
// 2. No-op sentinel tuple replays without mutation
// ─────────────────────────────────────────────────────────────────

#[test]
fn indel_replay_no_op_sentinel_applies_without_pool_mutation() {
    // count=2, every slot is kind=false / site=-1 (the no-op
    // sentinel). Pool stays unchanged; trace mirrors the sentinel.
    let pass = IndelPass::new(fixed_count(2), 0.5, Box::new(UniformBase));
    let initial = indel_test_sim();
    let initial_bases: Vec<u8> = initial.pool.as_slice().iter().map(|n| n.base).collect();
    let records = vec![
        count_rec(2),
        rec(ChoiceAddress::CorruptIndelKind(0), ChoiceValue::Bool(false)),
        rec(ChoiceAddress::CorruptIndelSite(0), ChoiceValue::Int(-1)),
        rec(ChoiceAddress::CorruptIndelKind(1), ChoiceValue::Bool(false)),
        rec(ChoiceAddress::CorruptIndelSite(1), ChoiceValue::Int(-1)),
    ];
    let (result, trace, _) = run_indel_replay(&pass, initial, records, None, None);
    let next = result.unwrap();

    assert_eq!(next.pool.len(), initial_bases.len());
    for (i, b) in initial_bases.iter().enumerate() {
        assert_eq!(next.pool.get(NucHandle::new(i as u32)).unwrap().base, *b);
    }
    assert_eq!(
        trace.find("corrupt.indel.site[0]").unwrap().value,
        ChoiceValue::Int(-1),
    );
    assert_eq!(
        trace.find("corrupt.indel.site[1]").unwrap().value,
        ChoiceValue::Int(-1),
    );
}

// ─────────────────────────────────────────────────────────────────
// 3. Malformed trace shape → PassError::Replay
// ─────────────────────────────────────────────────────────────────

#[test]
fn indel_replay_wrong_kind_record_surfaces_replay_error() {
    let pass = IndelPass::new(fixed_count(1), 1.0, Box::new(UniformBase));
    let records = vec![
        count_rec(1),
        // Wrong kind: site comes first, missing kind record.
        rec(ChoiceAddress::CorruptIndelSite(0), ChoiceValue::Int(5)),
        rec(ChoiceAddress::CorruptIndelBase(0), ChoiceValue::Base(b'A')),
    ];
    let (result, _, _) = run_indel_replay(&pass, indel_test_sim(), records, None, None);
    match result.unwrap_err() {
        PassError::Replay { reason, .. } => {
            assert!(matches!(reason, ReplayError::AddressMismatch { .. }));
        }
        other => panic!("expected Replay::AddressMismatch, got {other:?}"),
    }
}

#[test]
fn indel_replay_wrong_site_kind_surfaces_replay_error() {
    let pass = IndelPass::new(fixed_count(1), 1.0, Box::new(UniformBase));
    let records = vec![
        count_rec(1),
        rec(ChoiceAddress::CorruptIndelKind(0), ChoiceValue::Bool(true)),
        // Wrong kind: Bool instead of Int.
        rec(ChoiceAddress::CorruptIndelSite(0), ChoiceValue::Bool(true)),
        rec(ChoiceAddress::CorruptIndelBase(0), ChoiceValue::Base(b'A')),
    ];
    let (result, _, _) = run_indel_replay(&pass, indel_test_sim(), records, None, None);
    match result.unwrap_err() {
        PassError::Replay { reason, .. } => {
            assert!(matches!(reason, ReplayError::ValueKindMismatch { .. }));
        }
        other => panic!("expected Replay::ValueKindMismatch, got {other:?}"),
    }
}

#[test]
fn indel_replay_missing_insertion_base_surfaces_replay_error() {
    // Insertion slot but no base record at the cursor.
    let pass = IndelPass::new(fixed_count(1), 1.0, Box::new(UniformBase));
    let records = vec![
        count_rec(1),
        rec(ChoiceAddress::CorruptIndelKind(0), ChoiceValue::Bool(true)),
        rec(ChoiceAddress::CorruptIndelSite(0), ChoiceValue::Int(0)),
        // No CorruptIndelBase(0). Cursor exhausts.
    ];
    let (result, _, _) = run_indel_replay(&pass, indel_test_sim(), records, None, None);
    match result.unwrap_err() {
        PassError::Replay { reason, .. } => {
            assert!(matches!(reason, ReplayError::Exhausted { .. }));
        }
        other => panic!("expected Replay::Exhausted, got {other:?}"),
    }
}

// ─────────────────────────────────────────────────────────────────
// 4. Bounds validation
// ─────────────────────────────────────────────────────────────────

#[test]
fn indel_replay_insertion_with_negative_site_rejects() {
    let pass = IndelPass::new(fixed_count(1), 1.0, Box::new(UniformBase));
    let records = vec![
        count_rec(1),
        rec(ChoiceAddress::CorruptIndelKind(0), ChoiceValue::Bool(true)),
        // Insertion at site -1 is not a valid sentinel.
        rec(ChoiceAddress::CorruptIndelSite(0), ChoiceValue::Int(-1)),
        rec(ChoiceAddress::CorruptIndelBase(0), ChoiceValue::Base(b'A')),
    ];
    let (result, _, _) = run_indel_replay(&pass, indel_test_sim(), records, None, None);
    match result.unwrap_err() {
        PassError::InvalidPlanState { reason, .. } => {
            assert!(reason.contains("must be non-negative"));
        }
        other => panic!("expected InvalidPlanState, got {other:?}"),
    }
}

#[test]
fn indel_replay_insertion_site_out_of_range_rejects() {
    // Pool has 12 bytes; insertion at site 99 is out of bounds.
    let pass = IndelPass::new(fixed_count(1), 1.0, Box::new(UniformBase));
    let records = vec![
        count_rec(1),
        rec(ChoiceAddress::CorruptIndelKind(0), ChoiceValue::Bool(true)),
        rec(ChoiceAddress::CorruptIndelSite(0), ChoiceValue::Int(99)),
        rec(ChoiceAddress::CorruptIndelBase(0), ChoiceValue::Base(b'A')),
    ];
    let (result, _, _) = run_indel_replay(&pass, indel_test_sim(), records, None, None);
    match result.unwrap_err() {
        PassError::InvalidPlanState { reason, .. } => {
            assert!(reason.contains("99"));
            assert!(reason.contains("pool length"));
        }
        other => panic!("expected InvalidPlanState, got {other:?}"),
    }
}

#[test]
fn indel_replay_deletion_site_out_of_range_rejects() {
    let pass = IndelPass::new(fixed_count(1), 0.0, Box::new(UniformBase));
    let records = vec![
        count_rec(1),
        rec(ChoiceAddress::CorruptIndelKind(0), ChoiceValue::Bool(false)),
        rec(ChoiceAddress::CorruptIndelSite(0), ChoiceValue::Int(99)),
    ];
    let (result, _, _) = run_indel_replay(&pass, indel_test_sim(), records, None, None);
    match result.unwrap_err() {
        PassError::InvalidPlanState { reason, .. } => {
            assert!(reason.contains("99"));
            assert!(reason.contains("pool length"));
        }
        other => panic!("expected InvalidPlanState, got {other:?}"),
    }
}

// ─────────────────────────────────────────────────────────────────
// 5. Insertion base not in distribution support
// ─────────────────────────────────────────────────────────────────

#[test]
fn indel_replay_insertion_base_not_in_support_rejects() {
    // UniformBase support = {A,C,G,T}. Recorded base = 'N' is not
    // in the support → the live sampler could not have produced it,
    // so replay rejects rather than force-applies.
    let pass = IndelPass::new(fixed_count(1), 1.0, Box::new(UniformBase));
    let records = vec![
        count_rec(1),
        rec(ChoiceAddress::CorruptIndelKind(0), ChoiceValue::Bool(true)),
        rec(ChoiceAddress::CorruptIndelSite(0), ChoiceValue::Int(5)),
        rec(ChoiceAddress::CorruptIndelBase(0), ChoiceValue::Base(b'N')),
    ];
    let (result, _, _) = run_indel_replay(&pass, indel_test_sim(), records, None, None);
    match result.unwrap_err() {
        PassError::ConstraintSampling { address, .. } => {
            assert_eq!(address, "corrupt.indel.base[0]");
        }
        other => panic!("expected ConstraintSampling, got {other:?}"),
    }
}

// ─────────────────────────────────────────────────────────────────
// 6. Frame-invalid productive tuple rejects
// ─────────────────────────────────────────────────────────────────

#[test]
fn indel_replay_frame_invalid_tuple_rejects_under_productive() {
    // Productive bundle requires the cumulative frame delta to be
    // 0 mod 3. A single +1 insertion in the productive sample
    // space (where the FrameDelta site is reachable) violates the
    // junction-frame invariant.
    let (cfg, sim) = balanced_dp_distribution_fixture();
    let contracts = productive();
    let pass = IndelPass::new(fixed_count(1), 1.0, Box::new(UniformBase));

    // Single insertion at site 3 in V (FrameDelta(+1) territory).
    // The tuple has a single +1 frame contribution; cumulative
    // delta = 1, not 0 mod 3 → frame check rejects.
    let records = vec![
        count_rec(1),
        rec(ChoiceAddress::CorruptIndelKind(0), ChoiceValue::Bool(true)),
        rec(ChoiceAddress::CorruptIndelSite(0), ChoiceValue::Int(3)),
        rec(ChoiceAddress::CorruptIndelBase(0), ChoiceValue::Base(b'A')),
    ];
    let (result, _, _) =
        run_indel_replay(&pass, sim, records, Some(&contracts), Some(&cfg));
    match result.unwrap_err() {
        PassError::ConstraintSampling { .. } => { /* expected */ }
        other => panic!("expected ConstraintSampling, got {other:?}"),
    }
}

#[test]
fn indel_replay_forbidden_site_rejects_under_productive() {
    // Sites inside the V anchor codon are Forbidden under
    // AnchorPreserved.V. A replayed insertion there must reject.
    let (cfg, sim) = balanced_dp_distribution_fixture();
    let contracts = productive();
    let pass = IndelPass::new(fixed_count(1), 1.0, Box::new(UniformBase));
    let records = vec![
        count_rec(1),
        rec(ChoiceAddress::CorruptIndelKind(0), ChoiceValue::Bool(true)),
        // Site 0 = first base of V anchor codon (TGT). Forbidden.
        rec(ChoiceAddress::CorruptIndelSite(0), ChoiceValue::Int(0)),
        rec(ChoiceAddress::CorruptIndelBase(0), ChoiceValue::Base(b'A')),
    ];
    let (result, _, _) =
        run_indel_replay(&pass, sim, records, Some(&contracts), Some(&cfg));
    match result.unwrap_err() {
        PassError::ConstraintSampling { .. } => { /* expected */ }
        other => panic!("expected ConstraintSampling, got {other:?}"),
    }
}

// ─────────────────────────────────────────────────────────────────
// 7. Post-event validator rejection
// ─────────────────────────────────────────────────────────────────
//
// A custom contract that classifies every event as FrameNeutral
// (so the frame check passes) but rejects the assembled
// post-state via `verify`. Replay must surface the post-event
// failure as `ConstraintSampling`.

struct RejectAssembledPostState;
impl crate::contract::Contract for RejectAssembledPostState {
    fn name(&self) -> &str {
        "reject_assembled_post_state"
    }
    fn verify(
        &self,
        _sim: &Simulation,
        _refdata: Option<&crate::refdata::RefDataConfig>,
    ) -> Result<(), crate::contract::ContractViolation> {
        Err(crate::contract::ContractViolation::new(
            self.name(),
            "always rejected at post-event",
        ))
    }
    fn admissible_indel_class_at(
        &self,
        _sim: &Simulation,
        _refdata: Option<&crate::refdata::RefDataConfig>,
        _site: u32,
        _kind: crate::contract::IndelKindHint,
    ) -> crate::contract::IndelEventClass {
        crate::contract::IndelEventClass::FrameNeutral
    }
}

#[test]
fn indel_replay_post_event_rejection_surfaces_constraint_sampling() {
    use crate::contract::ContractSet;
    let sim = indel_test_sim();
    let contracts = ContractSet::new().with(Box::new(RejectAssembledPostState));
    let pass = IndelPass::new(fixed_count(1), 1.0, Box::new(UniformBase));

    let records = vec![
        count_rec(1),
        rec(ChoiceAddress::CorruptIndelKind(0), ChoiceValue::Bool(true)),
        rec(ChoiceAddress::CorruptIndelSite(0), ChoiceValue::Int(5)),
        rec(ChoiceAddress::CorruptIndelBase(0), ChoiceValue::Base(b'A')),
    ];
    let (result, _, _) = run_indel_replay(&pass, sim, records, Some(&contracts), None);
    match result.unwrap_err() {
        PassError::ConstraintSampling { .. } => { /* expected */ }
        other => panic!("expected ConstraintSampling, got {other:?}"),
    }
    let _ = FilteredSampleError::EmptyAdmissibleSupport; // keep import load-bearing
}
