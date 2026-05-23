//! Cross-segment overlap: V right edge into D, D left edge into V,
//! D right edge into J, J left edge into D, etc.
//!
//! When one segment's allele suffix happens to match the next
//! segment's allele prefix, the live-call walker can extend its
//! own range *past* the structural boundary into the neighbor's
//! bytes. This produces:
//!
//! - A non-zero NP claim by the extending segment.
//! - `OVERLAPS_OTHER_SEGMENT` flag on the live-call hypothesis.
//! - Shrunk NP region string in AIRR projection.
//! - Adjusted `*_sequence_end` / `*_sequence_start` for the
//!   extending segment.
//!
//! ## Invariants this module should test
//!
//! - **V right extends into D when D shares V's suffix**.
//! - **D left extends into V when V ends with D's prefix**.
//! - **D right extends into J when J starts with D's suffix**.
//! - **J left extends into D when D ends with J's prefix**.
//! - **Walker halts at pool end** (don't index past).
//! - **Overlap does NOT happen when no shared bytes exist** —
//!   safety check that we don't extend on noise.
//! - **Cross-segment overlap reflects in CIGAR** — the M-run
//!   extends past the structural region boundary in the AIRR
//!   alignment.
//! - **Overlap interacts correctly with NP extension** (V right
//!   first extends into NP1, then can extend further into D if
//!   D's prefix matches; the two extensions don't double-claim
//!   bytes).
//! - **Conflicting overlap is resolved**: if V wants to extend right
//!   and D wants to extend left into the same NP byte, exactly one
//!   wins per byte (the walker's claim policy in
//!   `compiled/walk/`).
//!
//! ## Adjacent in-tree coverage
//!
//! - `engine_rs/src/compiled/tests/overlaps.rs`

use super::common;
use genairr_engine::airr_record::build_airr_record;
use genairr_engine::assignment::TrimEnd;
use genairr_engine::compiled::{CompiledSimulator, ExecutionPolicy};
use genairr_engine::dist::{AllelePoolDist, Distribution, EmpiricalLengthDist};
use genairr_engine::ir::Segment;
use genairr_engine::live_call::HypothesisFlags;
use genairr_engine::pass::PassPlan;
use genairr_engine::passes::{AssembleSegmentPass, GenerateNPPass, SampleAllelePass, TrimPass};
use genairr_engine::refdata::{Allele, AlleleId, ChainType, RefDataConfig};

#[derive(Clone, Debug)]
struct ConstBaseDist(u8);

impl Distribution for ConstBaseDist {
    type Output = u8;
    fn sample(&self, _rng: &mut genairr_engine::rng::Rng) -> u8 {
        self.0
    }
    fn support(&self) -> Option<Vec<(u8, f64)>> {
        Some(vec![(self.0, 1.0)])
    }
}

/// Build a VDJ refdata where V*01's distinguishing 3' suffix
/// "TTT" coincides with D*01's leading 3 bases. With V trimmed
/// 3' by 3 and an empty NP1, V's right-extension can reach into
/// D's region for an exact-equivalent overlap placement.
fn v_d_overlap_refdata() -> (RefDataConfig, AlleleId, AlleleId, AlleleId, AlleleId) {
    let mut cfg = RefDataConfig::empty(ChainType::Vdj);
    let v01 = cfg.v_pool.push(Allele {
        name: "V*01".into(),
        gene: "V".into(),
        seq: b"AAACCCGGGTTT".to_vec(),
        segment: Segment::V,
        anchor: None,
    });
    let v02 = cfg.v_pool.push(Allele {
        name: "V*02".into(),
        gene: "V".into(),
        seq: b"AAACCCGGGAAA".to_vec(),
        segment: Segment::V,
        anchor: None,
    });
    // D*01 starts with TTT — exactly V*01's distinguishing suffix.
    let d01 = cfg.d_pool.push(Allele {
        name: "D*01".into(),
        gene: "D".into(),
        seq: b"TTTACGTAC".to_vec(),
        segment: Segment::D,
        anchor: None,
    });
    let j01 = cfg.j_pool.push(Allele {
        name: "J*01".into(),
        gene: "J".into(),
        seq: b"AAA".to_vec(),
        segment: Segment::J,
        anchor: None,
    });
    (cfg, v01, v02, d01, j01)
}

/// Plan that ends right after assembling D, so the V live call
/// reflects the post-D state including any overlap.
fn v_overlap_plan(
    cfg: &RefDataConfig,
    v_id: AlleleId,
    d_id: AlleleId,
    j_id: AlleleId,
    v_trim_3: i64,
    np1_len: i64,
    np1_base: u8,
) -> PassPlan {
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
    if v_trim_3 > 0 {
        plan.push(Box::new(TrimPass::new(
            Segment::V,
            TrimEnd::Three,
            Box::new(EmpiricalLengthDist::from_pairs(vec![(v_trim_3, 1.0)])),
        )));
    }
    plan.push(Box::new(AssembleSegmentPass::new(Segment::V)));
    plan.push(Box::new(GenerateNPPass::new(
        Segment::Np1,
        Box::new(EmpiricalLengthDist::from_pairs(vec![(np1_len, 1.0)])),
        Box::new(ConstBaseDist(np1_base)),
    )));
    plan.push(Box::new(AssembleSegmentPass::new(Segment::D)));
    plan.push(Box::new(GenerateNPPass::new(
        Segment::Np2,
        Box::new(EmpiricalLengthDist::from_pairs(vec![(0, 1.0)])),
        Box::new(ConstBaseDist(b'A')),
    )));
    plan.push(Box::new(AssembleSegmentPass::new(Segment::J)));
    plan
}

fn compile_and_run<'a>(
    cfg: &'a RefDataConfig,
    plan: &'a PassPlan,
) -> genairr_engine::pass::Outcome {
    let compiled = CompiledSimulator::compile(plan, Some(cfg), None, ExecutionPolicy::Permissive)
        .expect("plan should compile");
    compiled
        .run_one(common::DEFAULT_SEED)
        .expect("plan should run")
}

fn v_hypothesis_has_overlap_flag(sim: &genairr_engine::ir::Simulation) -> bool {
    sim.live_calls
        .as_ref()
        .and_then(|lc| lc.get(Segment::V))
        .map(|sc| {
            sc.hypotheses
                .iter()
                .any(|h| h.flags.contains(HypothesisFlags::OVERLAPS_OTHER_SEGMENT))
        })
        .unwrap_or(false)
}

fn v_hypothesis_seq_end(sim: &genairr_engine::ir::Simulation) -> Option<u32> {
    sim.live_calls
        .as_ref()
        .and_then(|lc| lc.get(Segment::V))
        .and_then(|sc| sc.hypotheses.iter().map(|h| h.seq_end).max())
}

// ──────────────────────────────────────────────────────────────
// Tests
// ──────────────────────────────────────────────────────────────

#[test]
fn v_right_overlaps_d_when_d_starts_with_v_suffix() {
    // V*01 trimmed 3' by 3; D*01 starts with TTT (V*01's exact
    // suffix). After AssembleSegment(D), the V live call should
    // have narrowed back to {V*01} and the V hypothesis should
    // carry OVERLAPS_OTHER_SEGMENT.
    let (cfg, v01, _v02, d01, j01) = v_d_overlap_refdata();
    let plan = v_overlap_plan(&cfg, v01, d01, j01, 3, 0, b'A');
    let outcome = compile_and_run(&cfg, &plan);
    let sim = outcome.final_simulation();

    let names = common::v_call_names(sim, &cfg);
    assert_eq!(
        names,
        vec!["V*01".to_string()],
        "V's overlap into D's matching prefix should narrow v_call to V*01",
    );

    assert!(
        v_hypothesis_has_overlap_flag(sim),
        "V hypothesis must carry OVERLAPS_OTHER_SEGMENT after extending into D",
    );
}

#[test]
fn v_right_does_not_overlap_when_d_prefix_does_not_match() {
    // V*01's suffix is TTT; rebuild D so its prefix is "CCC"
    // (matches neither V*01 nor V*02 at ref pos 9). V's right-
    // extension walker should halt immediately; v_call remains
    // the post-trim widened {V*01, V*02}.
    let mut cfg = RefDataConfig::empty(ChainType::Vdj);
    let v01 = cfg.v_pool.push(Allele {
        name: "V*01".into(),
        gene: "V".into(),
        seq: b"AAACCCGGGTTT".to_vec(),
        segment: Segment::V,
        anchor: None,
    });
    let _v02 = cfg.v_pool.push(Allele {
        name: "V*02".into(),
        gene: "V".into(),
        seq: b"AAACCCGGGAAA".to_vec(),
        segment: Segment::V,
        anchor: None,
    });
    let d01 = cfg.d_pool.push(Allele {
        name: "D*01".into(),
        gene: "D".into(),
        seq: b"CCCACGTAC".to_vec(),
        segment: Segment::D,
        anchor: None,
    });
    let j01 = cfg.j_pool.push(Allele {
        name: "J*01".into(),
        gene: "J".into(),
        seq: b"AAA".to_vec(),
        segment: Segment::J,
        anchor: None,
    });
    let plan = v_overlap_plan(&cfg, v01, d01, j01, 3, 0, b'A');
    let outcome = compile_and_run(&cfg, &plan);
    let sim = outcome.final_simulation();

    let names = common::v_call_names(sim, &cfg);
    assert_eq!(
        names,
        vec!["V*01".to_string(), "V*02".to_string()],
        "non-matching D prefix; v_call stays widened",
    );

    assert!(
        !v_hypothesis_has_overlap_flag(sim),
        "no shared bytes → OVERLAPS_OTHER_SEGMENT must not be set",
    );
}

#[test]
fn d_left_overlaps_v_when_v_ends_with_d_prefix() {
    // V*01 = AAACCCGGG (9 bytes). D*01 = GGGACGTAC (prefix GGG
    // matches V*01's tail). Trim D_5 by 3 — D's structural
    // region exposes only ACGTAC, but D's left-extension walker
    // reads backward through the empty NP1 into V's last 3 bases
    // (GGG) which match D*01's trimmed prefix.
    //
    // D*02 = TTTACGTAC (prefix TTT) → drops once the walker
    // encounters the first backward byte.
    let mut cfg = RefDataConfig::empty(ChainType::Vdj);
    let v01 = cfg.v_pool.push(Allele {
        name: "V*01".into(),
        gene: "V".into(),
        seq: b"AAACCCGGG".to_vec(),
        segment: Segment::V,
        anchor: None,
    });
    let d01 = cfg.d_pool.push(Allele {
        name: "D*01".into(),
        gene: "D".into(),
        seq: b"GGGACGTAC".to_vec(),
        segment: Segment::D,
        anchor: None,
    });
    let _d02 = cfg.d_pool.push(Allele {
        name: "D*02".into(),
        gene: "D".into(),
        seq: b"TTTACGTAC".to_vec(),
        segment: Segment::D,
        anchor: None,
    });
    let j01 = cfg.j_pool.push(Allele {
        name: "J*01".into(),
        gene: "J".into(),
        seq: b"AAA".to_vec(),
        segment: Segment::J,
        anchor: None,
    });

    let mut plan = PassPlan::new();
    plan.push(Box::new(SampleAllelePass::new(
        Segment::V,
        Box::new(AllelePoolDist::restricted_uniform(&cfg.v_pool, vec![v01])),
    )));
    plan.push(Box::new(SampleAllelePass::new(
        Segment::D,
        Box::new(AllelePoolDist::restricted_uniform(&cfg.d_pool, vec![d01])),
    )));
    plan.push(Box::new(SampleAllelePass::new(
        Segment::J,
        Box::new(AllelePoolDist::restricted_uniform(&cfg.j_pool, vec![j01])),
    )));
    plan.push(Box::new(TrimPass::new(
        Segment::D,
        TrimEnd::Five,
        Box::new(EmpiricalLengthDist::from_pairs(vec![(3, 1.0)])),
    )));
    plan.push(Box::new(AssembleSegmentPass::new(Segment::V)));
    plan.push(Box::new(GenerateNPPass::new(
        Segment::Np1,
        Box::new(EmpiricalLengthDist::from_pairs(vec![(0, 1.0)])),
        Box::new(ConstBaseDist(b'A')),
    )));
    plan.push(Box::new(AssembleSegmentPass::new(Segment::D)));
    plan.push(Box::new(GenerateNPPass::new(
        Segment::Np2,
        Box::new(EmpiricalLengthDist::from_pairs(vec![(0, 1.0)])),
        Box::new(ConstBaseDist(b'A')),
    )));
    plan.push(Box::new(AssembleSegmentPass::new(Segment::J)));

    let outcome = compile_and_run(&cfg, &plan);
    let sim = outcome.final_simulation();
    let d_call = sim
        .live_calls
        .as_ref()
        .expect("live calls populated")
        .get(Segment::D)
        .cloned()
        .expect("D live call exists");
    assert_eq!(
        d_call.allele_call.to_ids(),
        vec![d01],
        "D's left-extension into V's GGG tail should narrow d_call to D*01",
    );
    let elastic = d_call
        .hypotheses
        .iter()
        .find(|h| h.flags.contains(HypothesisFlags::BOUNDARY_ELASTIC))
        .expect("D hypothesis should be flagged BOUNDARY_ELASTIC");
    assert!(
        elastic.flags.contains(HypothesisFlags::OVERLAPS_OTHER_SEGMENT),
        "D's left-extension must set OVERLAPS_OTHER_SEGMENT when reaching into V",
    );
    // Under conservative extension, only the byte that narrows
    // the call is claimed: ref_start moves from 3 to 2, not 0.
    assert_eq!(
        elastic.ref_start, 2,
        "conservative extension claims only the narrowing byte"
    );
}

#[test]
fn overlap_walker_halts_at_pool_end() {
    // Under the conservative extension policy (Phase 20) the
    // walker halts as soon as the call set collapses to a
    // singleton — the old "could run forever" worry is bounded
    // by the call narrowing, well before pool_len. This test
    // exercises the conservative invariant: even when the
    // downstream segment's bytes all match the upstream allele's
    // continuation, the walker stops the moment the call is fully
    // resolved.
    //
    // Two V alleles share a 5-base prefix and diverge only at
    // ref pos 5 (V*01='T', V*02='A'). D = all-T's. The first D
    // byte narrows {V01,V02}→{V01}; subsequent bytes can't narrow
    // further. seq_end is bounded by the resolution point, not
    // by pool_len.
    let mut cfg = RefDataConfig::empty(ChainType::Vdj);
    let v01 = cfg.v_pool.push(Allele {
        name: "V*01".into(),
        gene: "V".into(),
        seq: {
            let mut s = vec![b'A'; 12];
            s[5] = b'T';
            s
        },
        segment: Segment::V,
        anchor: None,
    });
    let _v02 = cfg.v_pool.push(Allele {
        name: "V*02".into(),
        gene: "V".into(),
        seq: vec![b'A'; 12],
        segment: Segment::V,
        anchor: None,
    });
    let d01 = cfg.d_pool.push(Allele {
        name: "D*01".into(),
        gene: "D".into(),
        seq: vec![b'T'; 5],
        segment: Segment::D,
        anchor: None,
    });
    let j01 = cfg.j_pool.push(Allele {
        name: "J*01".into(),
        gene: "J".into(),
        seq: b"AAA".to_vec(),
        segment: Segment::J,
        anchor: None,
    });
    // Trim V_3 by 7 → V structural region = first 5 'A's.
    // Post-assemble tie set = {V*01, V*02} (both match AAAAA).
    let plan = v_overlap_plan(&cfg, v01, d01, j01, 7, 0, b'A');
    let outcome = compile_and_run(&cfg, &plan);
    let sim = outcome.final_simulation();
    let seq_end = v_hypothesis_seq_end(sim).expect("V live call should exist");
    // Pool layout: V[0..5] + NP1[5..5] + D[5..10]. The first D
    // byte (pool pos 5, 'T') narrows {V01,V02}→{V01}; walker
    // halts. seq_end = 6 (well inside the pool).
    assert_eq!(
        seq_end, 6,
        "conservative walker halts at the byte that fully resolves the call",
    );
    assert!(seq_end < 10, "seq_end stays within pool_len");
    assert!(
        v_hypothesis_has_overlap_flag(sim),
        "extending into D must set OVERLAPS_OTHER_SEGMENT",
    );
}

#[test]
fn overlap_into_d_keeps_airr_v_sequence_end_at_structural_boundary() {
    // When V's right-extension reaches into D (overlap), the live
    // hypothesis records seq_end past V's structural end — but the
    // AIRR projection's column walker keeps V's `v_sequence_end`
    // at the structural V-region end (9). This is the documented
    // "internal overlap stays internal" boundary: NP-claim
    // extensions propagate to AIRR coordinates, but cross-segment
    // overlap does not.
    //
    // Under conservative extension (Phase 20), the live seq_end
    // advances by exactly the single byte that narrows
    // {V01,V02}→{V01} → seq_end = 10 internally; the AIRR
    // projection still reports 9.
    let (cfg, v01, _v02, d01, j01) = v_d_overlap_refdata();
    let plan = v_overlap_plan(&cfg, v01, d01, j01, 3, 0, b'A');
    let outcome = compile_and_run(&cfg, &plan);
    let sim = outcome.final_simulation();

    // Internal evidence: live hypothesis seq_end = 10
    // (1 D byte claimed under conservative).
    let seq_end = v_hypothesis_seq_end(sim).expect("V live call should exist");
    assert_eq!(seq_end, 10, "live hypothesis records overlap into D (1 byte)");

    // External (AIRR) projection: v_sequence_end stays at the
    // structural V end (9), NP1 is empty (no claim because NP1
    // was zero-length to start with).
    let rec = build_airr_record(&outcome, &cfg, "overlap-airr-bound");
    assert_eq!(rec.v_sequence_end, Some(9));
    assert_eq!(rec.np1, "");
    assert_eq!(rec.np1_length, 0);
}

#[test]
fn overlap_into_d_keeps_airr_v_cigar_at_structural_length() {
    // Companion to `overlap_into_d_keeps_airr_v_sequence_end...`:
    // even when V's live call has narrowed back to a singleton
    // via overlap into D, v_cigar reports the structural V length
    // (9M) — not 12M. This is by design (the AIRR walker keeps
    // each segment's M-run within its structural region).
    let (cfg, v01, _v02, d01, j01) = v_d_overlap_refdata();
    let plan = v_overlap_plan(&cfg, v01, d01, j01, 3, 0, b'A');
    let outcome = compile_and_run(&cfg, &plan);
    let rec = build_airr_record(&outcome, &cfg, "overlap-cigar-bound");
    // Live call has narrowed (verified by other tests).
    let names = common::v_call_names(outcome.final_simulation(), &cfg);
    assert_eq!(names, vec!["V*01".to_string()]);

    // But AIRR v_cigar still reports the structural span.
    assert_eq!(rec.v_cigar, "9M");
    // And v_germline_end stays at the structural post-trim end too
    // (NOT extended to 12 via overlap, unlike the NP-claim path).
    assert_eq!(rec.v_germline_end, Some(9));
}

#[test]
fn no_overlap_when_ample_bytes_present_but_none_match() {
    // V*01 trimmed 3' by 3; D has plenty of bytes but starts with
    // 'X' that's neither V*01 (T) nor V*02 (A) at ref pos 9.
    // Confirms the walker doesn't extend on noise even when there's
    // a long sequence to read.
    let mut cfg = RefDataConfig::empty(ChainType::Vdj);
    let v01 = cfg.v_pool.push(Allele {
        name: "V*01".into(),
        gene: "V".into(),
        seq: b"AAACCCGGGTTT".to_vec(),
        segment: Segment::V,
        anchor: None,
    });
    let _v02 = cfg.v_pool.push(Allele {
        name: "V*02".into(),
        gene: "V".into(),
        seq: b"AAACCCGGGAAA".to_vec(),
        segment: Segment::V,
        anchor: None,
    });
    // D's first base is 'G' — matches neither V*01[9]=T nor
    // V*02[9]=A. Then has 9 more bases the walker would otherwise
    // extend through.
    let d01 = cfg.d_pool.push(Allele {
        name: "D*01".into(),
        gene: "D".into(),
        seq: b"GACGTACGT".to_vec(),
        segment: Segment::D,
        anchor: None,
    });
    let j01 = cfg.j_pool.push(Allele {
        name: "J*01".into(),
        gene: "J".into(),
        seq: b"AAA".to_vec(),
        segment: Segment::J,
        anchor: None,
    });
    let plan = v_overlap_plan(&cfg, v01, d01, j01, 3, 0, b'A');
    let outcome = compile_and_run(&cfg, &plan);
    let sim = outcome.final_simulation();
    let names = common::v_call_names(sim, &cfg);
    assert_eq!(
        names,
        vec!["V*01".to_string(), "V*02".to_string()],
        "no shared bytes → v_call stays widened",
    );
    assert!(
        !v_hypothesis_has_overlap_flag(sim),
        "no shared bytes → OVERLAPS_OTHER_SEGMENT must not be set",
    );

    // And the AIRR projection reflects no extension.
    let rec = build_airr_record(&outcome, &cfg, "no-overlap");
    assert_eq!(rec.v_sequence_end, Some(9), "no overlap; v_sequence_end at structural V end");
    assert_eq!(rec.v_germline_end, Some(9));
}

// ──────────────────────────────────────────────────────────────
// VJ-chain V↔J overlap tests (no D between).
// ──────────────────────────────────────────────────────────────

#[test]
fn v_to_j_overlap_fires_when_v_suffix_matches_j_prefix_in_vj_chain() {
    // VJ chain: V*01's distinguishing 3' suffix is TTT; J*01 starts
    // with TTT. With V trimmed 3' by 3 and an empty NP1, V's right-
    // extension walker can reach into J's region (the bytes match
    // V*01's continuation). Result: v_call narrows back to {V*01},
    // V hypothesis carries OVERLAPS_OTHER_SEGMENT.
    let mut cfg = RefDataConfig::empty(ChainType::Vj);
    let v01 = cfg.v_pool.push(Allele {
        name: "V*01".into(),
        gene: "V".into(),
        seq: b"AAACCCGGGTTT".to_vec(),
        segment: Segment::V,
        anchor: None,
    });
    let _v02 = cfg.v_pool.push(Allele {
        name: "V*02".into(),
        gene: "V".into(),
        seq: b"AAACCCGGGAAA".to_vec(),
        segment: Segment::V,
        anchor: None,
    });
    // J*01 starts with TTT — exactly V*01's distinguishing suffix.
    let j01 = cfg.j_pool.push(Allele {
        name: "J*01".into(),
        gene: "J".into(),
        seq: b"TTTAAACCC".to_vec(),
        segment: Segment::J,
        anchor: None,
    });

    let mut plan = PassPlan::new();
    plan.push(Box::new(SampleAllelePass::new(
        Segment::V,
        Box::new(AllelePoolDist::restricted_uniform(&cfg.v_pool, vec![v01])),
    )));
    plan.push(Box::new(SampleAllelePass::new(
        Segment::J,
        Box::new(AllelePoolDist::restricted_uniform(&cfg.j_pool, vec![j01])),
    )));
    plan.push(Box::new(TrimPass::new(
        Segment::V,
        TrimEnd::Three,
        Box::new(EmpiricalLengthDist::from_pairs(vec![(3, 1.0)])),
    )));
    plan.push(Box::new(AssembleSegmentPass::new(Segment::V)));
    plan.push(Box::new(GenerateNPPass::new(
        Segment::Np1,
        Box::new(EmpiricalLengthDist::from_pairs(vec![(0, 1.0)])),
        Box::new(ConstBaseDist(b'A')),
    )));
    plan.push(Box::new(AssembleSegmentPass::new(Segment::J)));

    let outcome = compile_and_run(&cfg, &plan);
    let sim = outcome.final_simulation();

    let names = common::v_call_names(sim, &cfg);
    assert_eq!(
        names,
        vec!["V*01".to_string()],
        "V's overlap into J's matching prefix should narrow v_call to V*01 in a VJ chain",
    );
    assert!(
        v_hypothesis_has_overlap_flag(sim),
        "V hypothesis must carry OVERLAPS_OTHER_SEGMENT after extending into J in a VJ chain",
    );
}

#[test]
fn v_to_j_overlap_does_not_fire_when_np1_separates_them() {
    // Same V*01 / V*02 pool as above. J*01 starts with TTT. With a
    // single NP1 byte 'A' between V and J (NP1 first byte at ref pos
    // 9 should be 'T' to match V*01), the NP byte 'A' does NOT match
    // V*01[9]='T' OR V*02[9]='A' uniquely:
    //   V*01[9] = 'T' — NP byte 'A' → mismatch
    //   V*02[9] = 'A' — NP byte 'A' → MATCH
    // So the walker actually narrows to V*02 via the NP byte. To
    // demonstrate "no overlap fires", use an NP byte that doesn't
    // match either: 'C' (V*01[9]='T', V*02[9]='A' — 'C' matches
    // neither). The walker halts immediately; no extension.
    let mut cfg = RefDataConfig::empty(ChainType::Vj);
    let v01 = cfg.v_pool.push(Allele {
        name: "V*01".into(),
        gene: "V".into(),
        seq: b"AAACCCGGGTTT".to_vec(),
        segment: Segment::V,
        anchor: None,
    });
    let _v02 = cfg.v_pool.push(Allele {
        name: "V*02".into(),
        gene: "V".into(),
        seq: b"AAACCCGGGAAA".to_vec(),
        segment: Segment::V,
        anchor: None,
    });
    let j01 = cfg.j_pool.push(Allele {
        name: "J*01".into(),
        gene: "J".into(),
        seq: b"TTTAAACCC".to_vec(),
        segment: Segment::J,
        anchor: None,
    });

    let mut plan = PassPlan::new();
    plan.push(Box::new(SampleAllelePass::new(
        Segment::V,
        Box::new(AllelePoolDist::restricted_uniform(&cfg.v_pool, vec![v01])),
    )));
    plan.push(Box::new(SampleAllelePass::new(
        Segment::J,
        Box::new(AllelePoolDist::restricted_uniform(&cfg.j_pool, vec![j01])),
    )));
    plan.push(Box::new(TrimPass::new(
        Segment::V,
        TrimEnd::Three,
        Box::new(EmpiricalLengthDist::from_pairs(vec![(3, 1.0)])),
    )));
    plan.push(Box::new(AssembleSegmentPass::new(Segment::V)));
    plan.push(Box::new(GenerateNPPass::new(
        Segment::Np1,
        Box::new(EmpiricalLengthDist::from_pairs(vec![(1, 1.0)])),
        Box::new(ConstBaseDist(b'C')),
    )));
    plan.push(Box::new(AssembleSegmentPass::new(Segment::J)));

    let outcome = compile_and_run(&cfg, &plan);
    let sim = outcome.final_simulation();

    // 'C' matches neither V*01[9]='T' nor V*02[9]='A' → no V right-
    // extension, no overlap. v_call stays widened.
    let names = common::v_call_names(sim, &cfg);
    assert_eq!(
        names,
        vec!["V*01".to_string(), "V*02".to_string()],
        "NP byte separates V and J — v_call must stay widened; got {:?}",
        names,
    );
    assert!(
        !v_hypothesis_has_overlap_flag(sim),
        "NP byte 'C' matches no V allele — V must NOT carry OVERLAPS_OTHER_SEGMENT",
    );
}

#[test]
fn np1_bytes_take_precedence_over_overlap_when_both_could_extend() {
    // V*01 suffix TTT. NP1 emits 3 'T' bytes (matches V*01).
    // D's prefix is also TTT (also matches). Under conservative
    // extension the walker proceeds one byte at a time until the
    // call is fully resolved. The first NP1 byte at ref pos 9
    // narrows {V*01,V*02}→{V*01}; the walker halts there because
    // subsequent bytes can no longer narrow the singleton tie set.
    //
    // The V hypothesis reaches ref_end=10 and seq_end=10 — only
    // the first NP1 byte is claimed; the remaining 2 NP1 bytes
    // (and all of D) stay outside V's hypothesis.
    let mut cfg = RefDataConfig::empty(ChainType::Vdj);
    let v01 = cfg.v_pool.push(Allele {
        name: "V*01".into(),
        gene: "V".into(),
        seq: b"AAACCCGGGTTT".to_vec(),
        segment: Segment::V,
        anchor: None,
    });
    let _v02 = cfg.v_pool.push(Allele {
        name: "V*02".into(),
        gene: "V".into(),
        seq: b"AAACCCGGGAAA".to_vec(),
        segment: Segment::V,
        anchor: None,
    });
    let d01 = cfg.d_pool.push(Allele {
        name: "D*01".into(),
        gene: "D".into(),
        seq: b"TTTACGTAC".to_vec(),
        segment: Segment::D,
        anchor: None,
    });
    let j01 = cfg.j_pool.push(Allele {
        name: "J*01".into(),
        gene: "J".into(),
        seq: b"AAA".to_vec(),
        segment: Segment::J,
        anchor: None,
    });
    let plan = v_overlap_plan(&cfg, v01, d01, j01, 3, 3, b'T');
    let outcome = compile_and_run(&cfg, &plan);
    let sim = outcome.final_simulation();
    let seq_end = v_hypothesis_seq_end(sim).expect("V live call should exist");
    // Pool layout: V[0..9] + NP1[9..12] + D[12..21].
    // V claims only the first NP1 column (the byte that narrows
    // V01,V02→V01); seq_end = 10.
    assert_eq!(
        seq_end, 10,
        "conservative walker stops after the single narrowing NP1 byte",
    );
    let names = common::v_call_names(sim, &cfg);
    assert_eq!(names, vec!["V*01".to_string()]);
}
