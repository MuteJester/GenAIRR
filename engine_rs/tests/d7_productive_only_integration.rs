//! End-to-end PRODUCTIVE_ONLY integration test.
//!
//! Under `respect=[productive()]`, the simulation produces 100%
//! in-frame junctions **by construction** — no retries, no
//! after-the-fact filtering. Every random choice the engine makes
//! is constrained at sampling time by the active `ContractSet`,
//! so the output stream contains only admissible sequences.
//!
//! This integration test exercises the full architecture stack:
//! - persistent IR (D1) carrying every revision
//! - addressed traces (D3) recording every choice
//! - the four-pattern architecture (passes + queries + traces +
//!   contracts)
//! - constraint-aware sampling via `sample_filtered` (D.6)
//! - the canonical `productive()` bundle (D.5)
//!
//! Both VJ (light chain) and VDJ (heavy chain) cases run, since
//! the joint NP1+NP2 constraint in VDJ (with NP2 compensating) is
//! the more sophisticated case.

use genairr_engine::contract::productive;
use genairr_engine::dist::{AllelePoolDist, Distribution, EmpiricalLengthDist, UniformBase};
use genairr_engine::ir::Segment;
use genairr_engine::junction::compute_junction;
use genairr_engine::pass::testing::PassRuntime;
use genairr_engine::pass::PassPlan;
use genairr_engine::passes::{AssembleSegmentPass, GenerateNPPass, SampleAllelePass};
use genairr_engine::refdata::{Allele, ChainType, RefDataConfig};
use genairr_engine::rng::Rng;
use genairr_engine::trace::ChoiceValue;

// ──────────────────────────────────────────────────────────────────
// Test fixtures: synthetic VJ and VDJ refdata
// ──────────────────────────────────────────────────────────────────

/// Synthetic VJ chain: V "AAACCCGGG" (anchor 6), J "TTTAAA" (anchor 0).
/// V_anchor_to_end = 3bp, J_anchor_to_W3 = 3bp → fixed junction = 6.
/// NP1 must satisfy (6 + NP1) % 3 == 0 → NP1 % 3 == 0.
fn make_vj_refdata() -> RefDataConfig {
    let mut cfg = RefDataConfig::empty(ChainType::Vj);
    let _ = cfg.v_pool.push(Allele {
        name: "v_test*01".into(),
        gene: "v_test".into(),
        seq: b"AAACCCGGG".to_vec(),
        segment: Segment::V,
        anchor: Some(6),
        functional_status: None,
        subregions: Vec::new(),
    });
    let _ = cfg.j_pool.push(Allele {
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

/// VJ fixture where the first NP1 base can complete either a stop codon
/// or a safe codon:
/// - V = "GGGTA", anchor 0. Positions 3..4 are "TA".
/// - NP1 length = 1.
/// - Candidate NP1 base 'A' creates "TAA" (stop); candidate 'C'
///   creates "TAC" (safe).
/// - J = "TTTAAA", anchor 0. With NP1 length 1, junction length is 9.
fn make_vj_np_stop_filter_refdata() -> RefDataConfig {
    let mut cfg = RefDataConfig::empty(ChainType::Vj);
    let _ = cfg.v_pool.push(Allele {
        name: "v_stop_filter*01".into(),
        gene: "v_stop_filter".into(),
        seq: b"GGGTA".to_vec(),
        segment: Segment::V,
        anchor: Some(0),
        functional_status: None,
        subregions: Vec::new(),
    });
    let _ = cfg.j_pool.push(Allele {
        name: "j_stop_filter*01".into(),
        gene: "j_stop_filter".into(),
        seq: b"TTTAAA".to_vec(),
        segment: Segment::J,
        anchor: Some(0),
        functional_status: None,
        subregions: Vec::new(),
    });
    cfg
}

/// VDJ fixture where the NP1 base can force a stop with the future D
/// prefix:
/// - V = "GGG", anchor 0.
/// - NP1 length = 1.
/// - D starts "AA...".
/// - Candidate NP1 base 'T' creates "TAA" using the first two D bases.
/// - Candidate 'C' creates "CAA" (safe).
/// - NP2 length = 2 so the final junction is in-frame.
fn make_vdj_np1_future_d_stop_refdata() -> RefDataConfig {
    let mut cfg = RefDataConfig::empty(ChainType::Vdj);
    let _ = cfg.v_pool.push(Allele {
        name: "v_future_d_stop*01".into(),
        gene: "v_future_d_stop".into(),
        seq: b"GGG".to_vec(),
        segment: Segment::V,
        anchor: Some(0),
        functional_status: None,
        subregions: Vec::new(),
    });
    let _ = cfg.d_pool.push(Allele {
        name: "d_future_d_stop*01".into(),
        gene: "d_future_d_stop".into(),
        seq: b"AAC".to_vec(),
        segment: Segment::D,
        anchor: None,
        functional_status: None,
        subregions: Vec::new(),
    });
    let _ = cfg.j_pool.push(Allele {
        name: "j_future_d_stop*01".into(),
        gene: "j_future_d_stop".into(),
        seq: b"TTTAAA".to_vec(),
        segment: Segment::J,
        anchor: Some(0),
        functional_status: None,
        subregions: Vec::new(),
    });
    cfg
}

#[derive(Clone, Debug)]
struct StopOrSafeBaseDist;

impl Distribution for StopOrSafeBaseDist {
    type Output = u8;

    fn sample(&self, _rng: &mut Rng) -> u8 {
        b'A'
    }

    fn support(&self) -> Option<Vec<(u8, f64)>> {
        Some(vec![(b'A', 1.0), (b'C', 1.0)])
    }
}

#[derive(Clone, Debug)]
struct FutureDStopOrSafeBaseDist;

impl Distribution for FutureDStopOrSafeBaseDist {
    type Output = u8;

    fn sample(&self, _rng: &mut Rng) -> u8 {
        b'T'
    }

    fn support(&self) -> Option<Vec<(u8, f64)>> {
        Some(vec![(b'T', 1.0), (b'C', 1.0)])
    }
}

/// Synthetic VDJ chain: V "AAACCCGGG" (anchor 6), D "TTTTTT"
/// (anchorless, 6bp), J "TTTAAA" (anchor 0).
/// At NP2 sampling time: V (3 from anchor to end) + NP1 + D (6) +
/// NP2 + J (3 from start to W+3) = junction. NP1 is free; NP2
/// compensates to make the total divisible by 3.
fn make_vdj_refdata() -> RefDataConfig {
    let mut cfg = RefDataConfig::empty(ChainType::Vdj);
    let _ = cfg.v_pool.push(Allele {
        name: "v_test*01".into(),
        gene: "v_test".into(),
        seq: b"AAACCCGGG".to_vec(),
        segment: Segment::V,
        anchor: Some(6),
        functional_status: None,
        subregions: Vec::new(),
    });
    let _ = cfg.d_pool.push(Allele {
        name: "d_test*01".into(),
        gene: "d_test".into(),
        seq: b"TTTTTT".to_vec(),
        segment: Segment::D,
        anchor: None,
        functional_status: None,
        subregions: Vec::new(),
    });
    let _ = cfg.j_pool.push(Allele {
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

fn build_vj_plan(refdata: &RefDataConfig) -> PassPlan {
    let length_dist = || {
        Box::new(EmpiricalLengthDist::from_pairs(
            (0..10).map(|i| (i, 1.0)).collect::<Vec<_>>(),
        ))
    };
    let mut plan = PassPlan::new();
    plan.push(Box::new(SampleAllelePass::new(
        Segment::V,
        Box::new(AllelePoolDist::uniform(&refdata.v_pool)),
    )));
    plan.push(Box::new(SampleAllelePass::new(
        Segment::J,
        Box::new(AllelePoolDist::uniform(&refdata.j_pool)),
    )));
    plan.push(Box::new(AssembleSegmentPass::new(Segment::V)));
    plan.push(Box::new(GenerateNPPass::new(
        Segment::Np1,
        length_dist(),
        Box::new(UniformBase),
    )));
    plan.push(Box::new(AssembleSegmentPass::new(Segment::J)));
    plan
}

fn build_vj_np_stop_filter_plan(refdata: &RefDataConfig) -> PassPlan {
    let mut plan = PassPlan::new();
    plan.push(Box::new(SampleAllelePass::new(
        Segment::V,
        Box::new(AllelePoolDist::uniform(&refdata.v_pool)),
    )));
    plan.push(Box::new(SampleAllelePass::new(
        Segment::J,
        Box::new(AllelePoolDist::uniform(&refdata.j_pool)),
    )));
    plan.push(Box::new(AssembleSegmentPass::new(Segment::V)));
    plan.push(Box::new(GenerateNPPass::new(
        Segment::Np1,
        Box::new(EmpiricalLengthDist::from_pairs(vec![(1, 1.0)])),
        Box::new(StopOrSafeBaseDist),
    )));
    plan.push(Box::new(AssembleSegmentPass::new(Segment::J)));
    plan
}

fn build_vdj_np1_future_d_stop_filter_plan(refdata: &RefDataConfig) -> PassPlan {
    let mut plan = PassPlan::new();
    plan.push(Box::new(SampleAllelePass::new(
        Segment::V,
        Box::new(AllelePoolDist::uniform(&refdata.v_pool)),
    )));
    plan.push(Box::new(SampleAllelePass::new(
        Segment::D,
        Box::new(AllelePoolDist::uniform(&refdata.d_pool)),
    )));
    plan.push(Box::new(SampleAllelePass::new(
        Segment::J,
        Box::new(AllelePoolDist::uniform(&refdata.j_pool)),
    )));
    plan.push(Box::new(AssembleSegmentPass::new(Segment::V)));
    plan.push(Box::new(GenerateNPPass::new(
        Segment::Np1,
        Box::new(EmpiricalLengthDist::from_pairs(vec![(1, 1.0)])),
        Box::new(FutureDStopOrSafeBaseDist),
    )));
    plan.push(Box::new(AssembleSegmentPass::new(Segment::D)));
    plan.push(Box::new(GenerateNPPass::new(
        Segment::Np2,
        Box::new(EmpiricalLengthDist::from_pairs(vec![(2, 1.0)])),
        Box::new(UniformBase),
    )));
    plan.push(Box::new(AssembleSegmentPass::new(Segment::J)));
    plan
}

fn build_vdj_plan(refdata: &RefDataConfig) -> PassPlan {
    let length_dist = || {
        Box::new(EmpiricalLengthDist::from_pairs(
            (0..10).map(|i| (i, 1.0)).collect::<Vec<_>>(),
        ))
    };
    let mut plan = PassPlan::new();
    plan.push(Box::new(SampleAllelePass::new(
        Segment::V,
        Box::new(AllelePoolDist::uniform(&refdata.v_pool)),
    )));
    plan.push(Box::new(SampleAllelePass::new(
        Segment::D,
        Box::new(AllelePoolDist::uniform(&refdata.d_pool)),
    )));
    plan.push(Box::new(SampleAllelePass::new(
        Segment::J,
        Box::new(AllelePoolDist::uniform(&refdata.j_pool)),
    )));
    plan.push(Box::new(AssembleSegmentPass::new(Segment::V)));
    plan.push(Box::new(GenerateNPPass::new(
        Segment::Np1,
        length_dist(),
        Box::new(UniformBase),
    )));
    plan.push(Box::new(AssembleSegmentPass::new(Segment::D)));
    plan.push(Box::new(GenerateNPPass::new(
        Segment::Np2,
        length_dist(),
        Box::new(UniformBase),
    )));
    plan.push(Box::new(AssembleSegmentPass::new(Segment::J)));
    plan
}

// ──────────────────────────────────────────────────────────────────
// Tests
// ──────────────────────────────────────────────────────────────────

const SEED_RANGE: u64 = 50;

#[test]
fn vj_productive_produces_in_frame_junctions_across_many_seeds() {
    let refdata = make_vj_refdata();
    let plan = build_vj_plan(&refdata);
    let contracts = productive();

    for seed in 0..SEED_RANGE {
        let outcome = PassRuntime::execute_with_context(
            &plan,
            Default::default(),
            seed,
            Some(&refdata),
            Some(&contracts),
        );
        let junction = compute_junction(outcome.final_simulation(), &refdata)
            .expect("junction should be defined for VJ with assembled regions");
        assert!(
            junction.is_in_frame(),
            "VJ seed {} produced out-of-frame junction (length {})",
            seed,
            junction.length
        );
    }
}

#[test]
fn vdj_productive_produces_in_frame_junctions_across_many_seeds() {
    let refdata = make_vdj_refdata();
    let plan = build_vdj_plan(&refdata);
    let contracts = productive();

    for seed in 0..SEED_RANGE {
        let outcome = PassRuntime::execute_with_context(
            &plan,
            Default::default(),
            seed,
            Some(&refdata),
            Some(&contracts),
        );
        let junction = compute_junction(outcome.final_simulation(), &refdata)
            .expect("junction should be defined for VDJ with assembled regions");
        assert!(
            junction.is_in_frame(),
            "VDJ seed {} produced out-of-frame junction (length {})",
            seed,
            junction.length
        );
    }
}

#[test]
fn vj_productive_full_bundle_verifies_post_pass() {
    // Not just "in frame" — the full productive() bundle (frame +
    // anchors + no-stops) should verify after each run. This is
    // the end-to-end soundness check.
    let refdata = make_vj_refdata();
    let plan = build_vj_plan(&refdata);
    let contracts = productive();

    for seed in 0..SEED_RANGE {
        let outcome = PassRuntime::execute_with_context(
            &plan,
            Default::default(),
            seed,
            Some(&refdata),
            Some(&contracts),
        );
        assert!(
            contracts
                .verify(outcome.final_simulation(), Some(&refdata))
                .is_ok(),
            "VJ seed {} produced a productive() violation",
            seed
        );
    }
}

#[test]
fn vj_productive_filters_np_base_that_would_complete_stop_codon() {
    let refdata = make_vj_np_stop_filter_refdata();
    let plan = build_vj_np_stop_filter_plan(&refdata);
    let contracts = productive();

    let outcome = PassRuntime::execute_with_context(
        &plan,
        Default::default(),
        0,
        Some(&refdata),
        Some(&contracts),
    );

    let sampled_base = match outcome.trace.find("np.np1.bases[0]").unwrap().value {
        ChoiceValue::Base(b) => b,
        _ => panic!("wrong variant"),
    };
    assert_eq!(
        sampled_base, b'C',
        "productive() must reject the stop-forming candidate base A"
    );
    assert!(contracts
        .verify(outcome.final_simulation(), Some(&refdata))
        .is_ok());
}

#[test]
fn vj_unconstrained_stop_filter_fixture_exposes_stop_codon() {
    let refdata = make_vj_np_stop_filter_refdata();
    let plan = build_vj_np_stop_filter_plan(&refdata);
    let contracts = productive();

    let outcome = PassRuntime::execute_with_refdata(&plan, Default::default(), 0, &refdata);

    let sampled_base = match outcome.trace.find("np.np1.bases[0]").unwrap().value {
        ChoiceValue::Base(b) => b,
        _ => panic!("wrong variant"),
    };
    assert_eq!(sampled_base, b'A');
    let violations = contracts
        .verify(outcome.final_simulation(), Some(&refdata))
        .unwrap_err();
    assert!(violations
        .iter()
        .any(|v| v.contract_name == "no_stop_codon_in_junction"));
}

#[test]
fn vdj_productive_filters_np1_base_that_would_force_future_d_stop_codon() {
    let refdata = make_vdj_np1_future_d_stop_refdata();
    let plan = build_vdj_np1_future_d_stop_filter_plan(&refdata);
    let contracts = productive();

    let outcome = PassRuntime::execute_with_context(
        &plan,
        Default::default(),
        0,
        Some(&refdata),
        Some(&contracts),
    );

    let sampled_base = match outcome.trace.find("np.np1.bases[0]").unwrap().value {
        ChoiceValue::Base(b) => b,
        _ => panic!("wrong variant"),
    };
    assert_eq!(
        sampled_base, b'C',
        "productive() must reject NP1 base T because future D would make TAA"
    );
    assert!(contracts
        .verify(outcome.final_simulation(), Some(&refdata))
        .is_ok());
}

#[test]
fn vdj_unconstrained_future_d_stop_fixture_exposes_stop_codon() {
    let refdata = make_vdj_np1_future_d_stop_refdata();
    let plan = build_vdj_np1_future_d_stop_filter_plan(&refdata);
    let contracts = productive();

    let outcome = PassRuntime::execute_with_refdata(&plan, Default::default(), 0, &refdata);

    let sampled_base = match outcome.trace.find("np.np1.bases[0]").unwrap().value {
        ChoiceValue::Base(b) => b,
        _ => panic!("wrong variant"),
    };
    assert_eq!(sampled_base, b'T');
    let violations = contracts
        .verify(outcome.final_simulation(), Some(&refdata))
        .unwrap_err();
    assert!(violations
        .iter()
        .any(|v| v.contract_name == "no_stop_codon_in_junction"));
}

#[test]
fn vj_unconstrained_produces_out_of_frame_junctions_too() {
    // Sanity check: WITHOUT the productive contract, the same plan
    // produces out-of-frame junctions sometimes (uniform NP1 length
    // distribution → 6/10 lengths give out-of-frame). This is the
    // "negative control" — proves the constraint actually does work.
    let refdata = make_vj_refdata();
    let plan = build_vj_plan(&refdata);

    let mut in_frame = 0u32;
    let mut out_of_frame = 0u32;
    for seed in 0..SEED_RANGE {
        let outcome = PassRuntime::execute_with_refdata(&plan, Default::default(), seed, &refdata);
        let junction =
            compute_junction(outcome.final_simulation(), &refdata).expect("junction defined");
        if junction.is_in_frame() {
            in_frame += 1;
        } else {
            out_of_frame += 1;
        }
    }

    // Without contracts, expect ~33% in-frame (every third length).
    // Both buckets should be non-trivially populated — proves the
    // constraint *isn't* in effect.
    assert!(
        out_of_frame > 0,
        "Expected at least one out-of-frame junction without contracts"
    );
    assert!(
        in_frame > 0,
        "Expected at least one in-frame junction without contracts"
    );
}

#[test]
fn vj_productive_is_deterministic_under_same_seed() {
    let refdata = make_vj_refdata();
    let plan = build_vj_plan(&refdata);
    let contracts = productive();

    let oa = PassRuntime::execute_with_context(
        &plan,
        Default::default(),
        0xcafe,
        Some(&refdata),
        Some(&contracts),
    );
    let ob = PassRuntime::execute_with_context(
        &plan,
        Default::default(),
        0xcafe,
        Some(&refdata),
        Some(&contracts),
    );

    // Trace identical.
    assert_eq!(oa.trace.choices(), ob.trace.choices());
    // Pool identical.
    let pa = &oa.final_simulation().pool;
    let pb = &ob.final_simulation().pool;
    assert_eq!(pa.len(), pb.len());
}

#[test]
fn vj_productive_trace_only_contains_in_frame_np1_lengths() {
    // Direct evidence that the constraint kicks in at sampling time:
    // every NP1 length recorded to the trace satisfies NP1 % 3 == 0
    // (the constraint logic for VJ with V_anchor_to_end = 3 +
    // J_anchor_to_W3 = 3 → fixed = 6 → NP1 must be 0 mod 3).
    let refdata = make_vj_refdata();
    let plan = build_vj_plan(&refdata);
    let contracts = productive();

    for seed in 0..SEED_RANGE {
        let outcome = PassRuntime::execute_with_context(
            &plan,
            Default::default(),
            seed,
            Some(&refdata),
            Some(&contracts),
        );
        let np1_len = match outcome.trace.find("np.np1.length").unwrap().value {
            ChoiceValue::Int(n) => n,
            _ => panic!("wrong variant"),
        };
        assert_eq!(
            np1_len % 3,
            0,
            "VJ seed {} sampled NP1 length {} (not divisible by 3)",
            seed,
            np1_len
        );
    }
}

#[test]
fn vdj_productive_trace_constrains_np2_for_frame() {
    // For VDJ, the constraint is on NP2 (NP1 is free). Verify the
    // recorded NP2 length always makes the junction in-frame.
    let refdata = make_vdj_refdata();
    let plan = build_vdj_plan(&refdata);
    let contracts = productive();

    for seed in 0..SEED_RANGE {
        let outcome = PassRuntime::execute_with_context(
            &plan,
            Default::default(),
            seed,
            Some(&refdata),
            Some(&contracts),
        );
        // The proof is via the resulting junction, since NP2's
        // valid value depends on the specific NP1 + D + V/J state.
        let junction =
            compute_junction(outcome.final_simulation(), &refdata).expect("junction defined");
        assert!(
            junction.is_in_frame(),
            "VDJ seed {} produced out-of-frame junction (length {})",
            seed,
            junction.length
        );
    }
}
