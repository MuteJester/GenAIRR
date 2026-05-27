//! Trace-injected replay tests for `S5FMutationPass` (Tier 2, S5F).
//!
//! The replay branch lives at the top of S5F's per-step loop. These
//! tests construct a `PassContext` with a cursor directly so they
//! can observe RNG consumption and exercise the validation chain
//! (S5F naturality + contract admissibility) without depending on
//! full-plan execution.
//!
//! Per the migration invariant: a recorded value that violates the
//! current contract bundle or has zero S5F natural weight is a
//! **structured error** — not a silent skip — in replay mode.

use super::{s5f_test_sim, s5f_uniform_kernel, s5f_zero_kernel};
use crate::address::ChoiceAddress;
use crate::contract::{AnchorPreserved, ContractSet};
use crate::dist::EmpiricalLengthDist;
use crate::ir::{Region, Segment, Simulation};
use crate::pass::{Pass, PassContext, PassError};
use crate::passes::S5FMutationPass;
use crate::refdata::{Allele, AlleleId, ChainType, RefDataConfig};
use crate::replay::{ReplayError, TraceCursor};
use crate::rng::Rng;
use crate::trace::{ChoiceRecord, ChoiceValue, Trace};

fn rec(addr: ChoiceAddress, v: ChoiceValue) -> ChoiceRecord {
    ChoiceRecord::new(addr.to_string(), v)
}

fn make_pass(kernel: crate::s5f::S5FKernel, count: i64) -> S5FMutationPass {
    S5FMutationPass::new(
        kernel,
        Box::new(EmpiricalLengthDist::from_pairs(vec![(count, 1.0)])),
    )
}

/// Drive S5F once against a caller-supplied cursor + (optional)
/// contracts. Returns `(result, output trace, total RNG words
/// consumed by the run)`.
fn run_s5f_replay(
    pass: &S5FMutationPass,
    sim: Simulation,
    records: Vec<ChoiceRecord>,
    contracts: Option<&ContractSet>,
    refdata: Option<&RefDataConfig>,
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

// ── Helpers for finding a valid (site, base) pair ────────────────

/// Sample once under fresh RNG to discover a known-good
/// `(site, base)` the kernel + pool agree on. The discovered pair
/// is what replay must reproduce.
fn discover_valid_event(
    pass: &S5FMutationPass,
    sim: &Simulation,
) -> (i64, u8) {
    use crate::pass::PassPlan;
    use crate::pass::testing::PassRuntime;
    let mut plan = PassPlan::new();
    plan.push(Box::new(S5FMutationPass::new(
        s5f_uniform_kernel(),
        Box::new(EmpiricalLengthDist::from_pairs(vec![(1, 1.0)])),
    )));
    let _ = pass; // future-proof if the pass-kernel decoupling drifts.
    let outcome = PassRuntime::execute(&plan, sim.clone(), 42);
    let site = match outcome.trace.find("mutate.s5f.site[0]").unwrap().value {
        ChoiceValue::Int(s) => s,
        _ => unreachable!(),
    };
    let base = match outcome.trace.find("mutate.s5f.base[0]").unwrap().value {
        ChoiceValue::Base(b) => b,
        _ => unreachable!(),
    };
    (site, base)
}

// ─────────────────────────────────────────────────────────────────
// 1. Valid replay: recorded site/base applies, RNG untouched
// ─────────────────────────────────────────────────────────────────

#[test]
fn s5f_replay_consumes_recorded_event_with_zero_rng_consumption() {
    let pass = make_pass(s5f_uniform_kernel(), 1);
    let sim = s5f_test_sim();
    let (site, base) = discover_valid_event(&pass, &sim);

    let records = vec![
        // Count = 1 (consumed by sample_validated_count).
        rec(ChoiceAddress::MutateS5fCount, ChoiceValue::Int(1)),
        // The one (site, base) draw.
        rec(ChoiceAddress::MutateS5fSite(0), ChoiceValue::Int(site)),
        rec(ChoiceAddress::MutateS5fBase(0), ChoiceValue::Base(base)),
    ];
    let (result, trace, rng_words) = run_s5f_replay(&pass, sim.clone(), records, None, None);
    let next = result.unwrap();

    // The pool now carries the replayed base at the replayed site.
    assert_eq!(
        next.pool.get(crate::ir::NucHandle::new(site as u32)).unwrap().base,
        base,
        "replay must write the recorded base at the recorded site",
    );
    // The output trace echoes the replayed (site, base) at the
    // canonical S5F addresses.
    assert_eq!(
        trace.find("mutate.s5f.site[0]").unwrap().value,
        ChoiceValue::Int(site),
    );
    assert_eq!(
        trace.find("mutate.s5f.base[0]").unwrap().value,
        ChoiceValue::Base(base),
    );
    // **Zero RNG consumption for the (site, base) event** — the
    // count helper bumps the cursor and not the rng either, so the
    // total expected is 0 words consumed.
    assert_eq!(rng_words, 0);
}

// ─────────────────────────────────────────────────────────────────
// 2. Wrong address / wrong kind errors
// ─────────────────────────────────────────────────────────────────

#[test]
fn s5f_replay_wrong_site_kind_surfaces_replay_error() {
    let pass = make_pass(s5f_uniform_kernel(), 1);
    let sim = s5f_test_sim();
    let records = vec![
        rec(ChoiceAddress::MutateS5fCount, ChoiceValue::Int(1)),
        // Should be Int; we put Bool.
        rec(ChoiceAddress::MutateS5fSite(0), ChoiceValue::Bool(true)),
        rec(ChoiceAddress::MutateS5fBase(0), ChoiceValue::Base(b'A')),
    ];
    let (result, _, _) = run_s5f_replay(&pass, sim, records, None, None);
    match result.unwrap_err() {
        PassError::Replay { reason, .. } => {
            assert!(matches!(reason, ReplayError::ValueKindMismatch { .. }));
        }
        other => panic!("expected Replay::ValueKindMismatch, got {other:?}"),
    }
}

#[test]
fn s5f_replay_wrong_base_kind_surfaces_replay_error() {
    let pass = make_pass(s5f_uniform_kernel(), 1);
    let sim = s5f_test_sim();
    let records = vec![
        rec(ChoiceAddress::MutateS5fCount, ChoiceValue::Int(1)),
        rec(ChoiceAddress::MutateS5fSite(0), ChoiceValue::Int(5)),
        // Should be Base; we put Int.
        rec(ChoiceAddress::MutateS5fBase(0), ChoiceValue::Int(0)),
    ];
    let (result, _, _) = run_s5f_replay(&pass, sim, records, None, None);
    match result.unwrap_err() {
        PassError::Replay { reason, .. } => {
            assert!(matches!(reason, ReplayError::ValueKindMismatch { .. }));
        }
        other => panic!("expected Replay::ValueKindMismatch, got {other:?}"),
    }
}

#[test]
fn s5f_replay_wrong_site_address_surfaces_replay_error() {
    let pass = make_pass(s5f_uniform_kernel(), 1);
    let sim = s5f_test_sim();
    let records = vec![
        rec(ChoiceAddress::MutateS5fCount, ChoiceValue::Int(1)),
        // Wrong address (uniform's site, not S5F's).
        rec(ChoiceAddress::MutateUniformSite(0), ChoiceValue::Int(5)),
        rec(ChoiceAddress::MutateS5fBase(0), ChoiceValue::Base(b'A')),
    ];
    let (result, _, _) = run_s5f_replay(&pass, sim, records, None, None);
    match result.unwrap_err() {
        PassError::Replay { reason, .. } => {
            assert!(matches!(reason, ReplayError::AddressMismatch { .. }));
        }
        other => panic!("expected Replay::AddressMismatch, got {other:?}"),
    }
}

// ─────────────────────────────────────────────────────────────────
// 3. S5F natural weight is zero → structured error
// ─────────────────────────────────────────────────────────────────

#[test]
fn s5f_replay_site_with_no_5mer_context_rejects() {
    // The first/last positions of the pool have no 5-mer context
    // (S5F needs `[pos-2, pos+2]`). A recorded site at pool boundary
    // can't have any positive natural weight.
    let pass = make_pass(s5f_uniform_kernel(), 1);
    let sim = s5f_test_sim();
    let records = vec![
        rec(ChoiceAddress::MutateS5fCount, ChoiceValue::Int(1)),
        // Site 0: no 5-mer context (pos < 2).
        rec(ChoiceAddress::MutateS5fSite(0), ChoiceValue::Int(0)),
        rec(ChoiceAddress::MutateS5fBase(0), ChoiceValue::Base(b'A')),
    ];
    let (result, _, _) = run_s5f_replay(&pass, sim, records, None, None);
    match result.unwrap_err() {
        PassError::ConstraintSampling { address, .. } => {
            assert_eq!(address, "mutate.s5f.base[0]");
        }
        other => panic!("expected ConstraintSampling, got {other:?}"),
    }
}

#[test]
fn s5f_replay_zero_natural_weight_rejects() {
    // Build a sim whose 5-mer at site 5 is well-defined, then use
    // a kernel with mutability=0 everywhere so every base at every
    // site has natural weight 0.
    let pass = make_pass(s5f_zero_kernel(), 1);
    let sim = s5f_test_sim();
    let records = vec![
        rec(ChoiceAddress::MutateS5fCount, ChoiceValue::Int(1)),
        rec(ChoiceAddress::MutateS5fSite(0), ChoiceValue::Int(5)),
        rec(ChoiceAddress::MutateS5fBase(0), ChoiceValue::Base(b'A')),
    ];
    let (result, _, _) = run_s5f_replay(&pass, sim, records, None, None);
    match result.unwrap_err() {
        PassError::ConstraintSampling { address, .. } => {
            assert_eq!(address, "mutate.s5f.base[0]");
        }
        other => panic!("expected ConstraintSampling, got {other:?}"),
    }
}

#[test]
fn s5f_replay_non_canonical_base_rejects() {
    // S5F's row is indexed over A/C/G/T only. A recorded `b'N'`
    // can't be a valid S5F event and must be rejected before any
    // pool mutation.
    let pass = make_pass(s5f_uniform_kernel(), 1);
    let sim = s5f_test_sim();
    let records = vec![
        rec(ChoiceAddress::MutateS5fCount, ChoiceValue::Int(1)),
        rec(ChoiceAddress::MutateS5fSite(0), ChoiceValue::Int(5)),
        rec(ChoiceAddress::MutateS5fBase(0), ChoiceValue::Base(b'N')),
    ];
    let (result, trace, _) = run_s5f_replay(&pass, sim.clone(), records, None, None);
    match result.unwrap_err() {
        PassError::ConstraintSampling { address, .. } => {
            assert_eq!(address, "mutate.s5f.base[0]");
        }
        other => panic!("expected ConstraintSampling, got {other:?}"),
    }
    // Pool untouched.
    assert_eq!(
        trace.find("mutate.s5f.site[0]"),
        None,
        "rejected replay must not write to the output trace",
    );
}

// ─────────────────────────────────────────────────────────────────
// 4. Contract rejects replayed candidate
// ─────────────────────────────────────────────────────────────────

/// Build a VJ fixture where the V anchor codon is `TGG` (W) so
/// `AnchorPreserved::V` constrains the first three pool positions
/// (the W codon). Any non-synonymous substitution at positions 0/1/2
/// changes the anchor amino acid and is rejected by the contract.
///
/// Pool layout: V[0..3] = TGG, V[3..8] = padding. Long enough so
/// position ≥ 2 has a 5-mer context.
fn anchor_locked_v_sim() -> (RefDataConfig, Simulation) {
    let mut cfg = RefDataConfig::empty(ChainType::Vj);
    let _ = cfg.v_pool.push(Allele {
        name: "v_lock*01".into(),
        gene: "v_lock".into(),
        seq: b"TGGAAAAA".to_vec(),
        segment: Segment::V,
        anchor: Some(0),
    });
    let _ = cfg.j_pool.push(Allele {
        name: "j_lock*01".into(),
        gene: "j_lock".into(),
        seq: b"TGG".to_vec(),
        segment: Segment::J,
        anchor: Some(0),
    });
    let mut sim = Simulation::new();
    for (i, &b) in b"TGGAAAAA".iter().enumerate() {
        let (next, _) = sim.with_nucleotide_pushed(crate::ir::Nucleotide::germline(
            b,
            i as u16,
            Segment::V,
        ));
        sim = next;
    }
    sim = sim.with_region_added(Region::new(
        Segment::V,
        crate::ir::NucHandle::new(0),
        crate::ir::NucHandle::new(8),
    ));
    sim = sim.with_allele_assigned(
        Segment::V,
        crate::assignment::AlleleInstance::new(AlleleId::new(0)),
    );
    (cfg, sim)
}

#[test]
fn s5f_replay_contract_rejects_replayed_candidate() {
    // Position 2 sits inside the V anchor codon TGG (the third G).
    // Substituting G → A makes the codon TGA — a stop, definitely
    // not Trp. `AnchorPreserved::V` rejects.
    let (cfg, sim) = anchor_locked_v_sim();
    let contracts = ContractSet::new().with(Box::new(AnchorPreserved::new(Segment::V)));
    let pass = make_pass(s5f_uniform_kernel(), 1);

    let records = vec![
        rec(ChoiceAddress::MutateS5fCount, ChoiceValue::Int(1)),
        rec(ChoiceAddress::MutateS5fSite(0), ChoiceValue::Int(2)),
        rec(ChoiceAddress::MutateS5fBase(0), ChoiceValue::Base(b'A')),
    ];
    let (result, _, _) =
        run_s5f_replay(&pass, sim, records, Some(&contracts), Some(&cfg));
    match result.unwrap_err() {
        PassError::ConstraintSampling { address, .. } => {
            assert_eq!(address, "mutate.s5f.base[0]");
        }
        other => panic!("expected ConstraintSampling, got {other:?}"),
    }
}

#[test]
fn s5f_replay_contract_accepts_admissible_replayed_candidate() {
    // Position 3 sits OUTSIDE the V anchor codon (TGG at [0..3));
    // position 3 is part of the padding AAAAA. Any substitution
    // there is admissible by AnchorPreserved.V. We pick base A
    // (same as current) so it's always natural.
    let (cfg, sim) = anchor_locked_v_sim();
    let contracts = ContractSet::new().with(Box::new(AnchorPreserved::new(Segment::V)));
    let pass = make_pass(s5f_uniform_kernel(), 1);

    let records = vec![
        rec(ChoiceAddress::MutateS5fCount, ChoiceValue::Int(1)),
        rec(ChoiceAddress::MutateS5fSite(0), ChoiceValue::Int(3)),
        rec(ChoiceAddress::MutateS5fBase(0), ChoiceValue::Base(b'A')),
    ];
    let (result, trace, rng_words) =
        run_s5f_replay(&pass, sim, records, Some(&contracts), Some(&cfg));
    assert!(result.is_ok());
    assert_eq!(
        trace.find("mutate.s5f.base[0]").unwrap().value,
        ChoiceValue::Base(b'A'),
    );
    // RNG untouched even with contracts active.
    assert_eq!(rng_words, 0);
}

// ─────────────────────────────────────────────────────────────────
// 5. Out-of-range site
// ─────────────────────────────────────────────────────────────────

#[test]
fn s5f_replay_site_out_of_range_surfaces_invalid_plan_state() {
    let pass = make_pass(s5f_uniform_kernel(), 1);
    let sim = s5f_test_sim();
    let records = vec![
        rec(ChoiceAddress::MutateS5fCount, ChoiceValue::Int(1)),
        rec(ChoiceAddress::MutateS5fSite(0), ChoiceValue::Int(9999)),
        rec(ChoiceAddress::MutateS5fBase(0), ChoiceValue::Base(b'A')),
    ];
    let (result, _, _) = run_s5f_replay(&pass, sim, records, None, None);
    match result.unwrap_err() {
        PassError::InvalidPlanState { reason, .. } => {
            assert!(reason.contains("9999"));
            assert!(reason.contains("out of pool range"));
        }
        other => panic!("expected InvalidPlanState, got {other:?}"),
    }
}
