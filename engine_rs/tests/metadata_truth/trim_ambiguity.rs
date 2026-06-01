//! Trim-induced allele-call ambiguity.
//!
//! When a V/D/J allele is trimmed at its 5' or 3' end, the surviving
//! bases may no longer be sufficient to distinguish the truth allele
//! from one or more other alleles in the same gene's pool. The engine
//! must surface ALL indistinguishable alleles in the corresponding
//! AIRR `v_call` / `d_call` / `j_call` field — no aligner / biologist /
//! oracle given only the bases could favor one of them over the
//! others, so the ground truth call set must reflect that.
//!
//! ## Invariants this module should test
//!
//! - **Aliases**: when two alleles in the pool share the same
//!   nucleotide sequence under different names, both names appear in
//!   the call set whenever either is the sampled truth.
//! - **Symmetric ambiguity**: when 5'/3' trim length is large enough
//!   that the surviving bytes match N>1 alleles equally well, all N
//!   alleles appear in the call set, regardless of which one was the
//!   sampled truth.
//! - **Boundary correctness**: `v_sequence_end` / `d_sequence_start`
//!   / `d_sequence_end` / `j_sequence_start` reflect the live-call
//!   hypothesis range, not the pre-trim allele bounds.
//! - **Germline coords reflect trim**: `v_germline_end` /
//!   `d_germline_*` / `j_germline_start` shift inward when the
//!   corresponding edge is trimmed (i.e. they index into the
//!   allele reference at the surviving range).
//! - **Truth allele always retained**: the call set ALWAYS contains
//!   the originally-sampled allele unless mutation / corruption
//!   actively excluded it.
//! - **Score-based narrowing**: an allele whose germline at the
//!   surviving positions matches strictly more bases than any other
//!   is the sole call (no spurious ambiguity).
//!
//! ## Adjacent in-tree coverage (not duplicated here)
//!
//! - `engine_rs/src/compiled/tests/v_extensions.rs::v_call_shrinks_when_np1_recreates_trimmed_suffix`
//! - `engine_rs/src/compiled/tests/curated.rs::curated_v_trim_widens_v_call`
//! - `engine_rs/src/live_call/tests/caller.rs::assembled_segment_live_call_keeps_trim_induced_ambiguity`
//! - `engine_rs/src/live_call/tests/walker_observer.rs::observer_matches_oracle_trim_induced_ambiguity`

use super::common;
use genairr_engine::airr_record::build_airr_record;
use genairr_engine::assignment::TrimEnd;
use genairr_engine::compiled::{CompiledSimulator, ExecutionPolicy};
use genairr_engine::dist::{AllelePoolDist, Distribution, EmpiricalLengthDist};
use genairr_engine::ir::Segment;
use genairr_engine::pass::PassPlan;
use genairr_engine::passes::{AssembleSegmentPass, GenerateNPPass, SampleAllelePass, TrimPass};
use genairr_engine::refdata::{AlleleId, RefDataConfig};

/// Test-only deterministic base distribution. Always samples the
/// configured byte. Redeclared per file because the compiled/tests
/// version is not `pub`.
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

// ──────────────────────────────────────────────────────────────
// Helpers: build VJ / VDJ plans driving the shared fixtures from
// `common.rs`. Pass `None` to disable a trim segment.
// ──────────────────────────────────────────────────────────────

#[allow(clippy::too_many_arguments)]
fn vj_plan(
    cfg: &RefDataConfig,
    v_id: AlleleId,
    j_id: AlleleId,
    v_trim_3: i64,
    j_trim_5: i64,
    np1_len: i64,
    np1_base: u8,
) -> PassPlan {
    let mut plan = PassPlan::new();
    plan.push(Box::new(SampleAllelePass::new(
        Segment::V,
        Box::new(AllelePoolDist::restricted_uniform(&cfg.v_pool, vec![v_id])),
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
    if j_trim_5 > 0 {
        plan.push(Box::new(TrimPass::new(
            Segment::J,
            TrimEnd::Five,
            Box::new(EmpiricalLengthDist::from_pairs(vec![(j_trim_5, 1.0)])),
        )));
    }
    plan.push(Box::new(AssembleSegmentPass::new(Segment::V)));
    plan.push(Box::new(GenerateNPPass::new(
        Segment::Np1,
        Box::new(EmpiricalLengthDist::from_pairs(vec![(np1_len, 1.0)])),
        Box::new(ConstBaseDist(np1_base)),
    )));
    plan.push(Box::new(AssembleSegmentPass::new(Segment::J)));
    plan
}

#[allow(clippy::too_many_arguments)]
fn vdj_plan(
    cfg: &RefDataConfig,
    v_id: AlleleId,
    d_id: AlleleId,
    j_id: AlleleId,
    v_trim_3: i64,
    d_trim_5: i64,
    d_trim_3: i64,
    j_trim_5: i64,
    np1_len: i64,
    np1_base: u8,
    np2_len: i64,
    np2_base: u8,
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
    if d_trim_5 > 0 {
        plan.push(Box::new(TrimPass::new(
            Segment::D,
            TrimEnd::Five,
            Box::new(EmpiricalLengthDist::from_pairs(vec![(d_trim_5, 1.0)])),
        )));
    }
    if d_trim_3 > 0 {
        plan.push(Box::new(TrimPass::new(
            Segment::D,
            TrimEnd::Three,
            Box::new(EmpiricalLengthDist::from_pairs(vec![(d_trim_3, 1.0)])),
        )));
    }
    if j_trim_5 > 0 {
        plan.push(Box::new(TrimPass::new(
            Segment::J,
            TrimEnd::Five,
            Box::new(EmpiricalLengthDist::from_pairs(vec![(j_trim_5, 1.0)])),
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
        Box::new(EmpiricalLengthDist::from_pairs(vec![(np2_len, 1.0)])),
        Box::new(ConstBaseDist(np2_base)),
    )));
    plan.push(Box::new(AssembleSegmentPass::new(Segment::J)));
    plan
}

fn run_vj(cfg: &RefDataConfig, plan: PassPlan) -> genairr_engine::pass::Outcome {
    let compiled = CompiledSimulator::compile(&plan, Some(cfg), None, ExecutionPolicy::Permissive)
        .expect("plan should compile");
    compiled
        .run_one(common::DEFAULT_SEED)
        .expect("plan should run")
}

// ──────────────────────────────────────────────────────────────
// Tests
// ──────────────────────────────────────────────────────────────

#[test]
fn aliases_both_appear_in_v_call_when_v01_is_truth() {
    // v01*01 and v01*02 are identical 12-base sequences under
    // different names. Sampling v01*01 with no trim should always
    // surface BOTH names in v_call — the bases give no way to
    // distinguish them.
    let cfg = common::vj_ambiguous_refdata();
    let v01 = common::allele_id_by_name(&cfg, Segment::V, "v01*01");
    let j01 = common::allele_id_by_name(&cfg, Segment::J, "j01*01");
    let plan = vj_plan(&cfg, v01, j01, 0, 0, 0, b'A');
    let outcome = run_vj(&cfg, plan);
    let sim = outcome.final_simulation();
    let names = common::v_call_names(sim, &cfg);
    assert_eq!(
        names,
        vec!["v01*01".to_string(), "v01*02".to_string()],
        "aliased V01 alleles must both appear in v_call",
    );
}

#[test]
fn aliases_both_appear_in_v_call_when_alias_is_truth() {
    // Symmetric — sample v01*02 (the alias). v_call still
    // contains both names, regardless of which one was sampled.
    let cfg = common::vj_ambiguous_refdata();
    let alias = common::allele_id_by_name(&cfg, Segment::V, "v01*02");
    let j01 = common::allele_id_by_name(&cfg, Segment::J, "j01*01");
    let plan = vj_plan(&cfg, alias, j01, 0, 0, 0, b'A');
    let outcome = run_vj(&cfg, plan);
    let sim = outcome.final_simulation();
    let names = common::v_call_names(sim, &cfg);
    assert_eq!(
        names,
        vec!["v01*01".to_string(), "v01*02".to_string()],
        "alias truth must still surface both names in v_call",
    );
}

#[test]
fn trim_widens_v_call_to_all_three_real_alleles() {
    // All four V entries share the 9-byte prefix AAACCCGGG. A
    // 3-byte 3' trim leaves only that prefix, which matches v01
    // (both aliases), v02, and v03 equally. Live v_call must hold
    // ALL FOUR.
    //
    // We pair with j02*01 (prefix GGG) rather than j01*01 (prefix
    // TTT): J's first byte is the first evidence position past V's
    // trim boundary, and wired V's right-extension to walk
    // into J's bytes in VJ chains. j02's G matches NONE of v01
    // (ref pos 9 = T) / v02 (= A) / v03 (= C) — extension halts
    // immediately, the widened call stays at 4. With j01 (T) v01
    // would narrow back to {v01*01, v01*02} via cross-segment
    // overlap, which is correct behavior but defeats THIS test's
    // intent of pinning the "trim widens to all 4" invariant.
    let cfg = common::vj_ambiguous_refdata();
    let v01 = common::allele_id_by_name(&cfg, Segment::V, "v01*01");
    let j02 = common::allele_id_by_name(&cfg, Segment::J, "j02*01");
    let plan = vj_plan(&cfg, v01, j02, 3, 0, 0, b'A');
    let outcome = run_vj(&cfg, plan);
    let sim = outcome.final_simulation();
    let names = common::v_call_names(sim, &cfg);
    assert_eq!(
        names,
        vec![
            "v01*01".to_string(),
            "v01*02".to_string(),
            "v02*01".to_string(),
            "v03*01".to_string(),
        ],
        "trimmed V leaves only the shared 9bp prefix → all 4 alleles match",
    );
}

#[test]
fn truth_allele_retained_in_widened_set_regardless_of_sampled_truth() {
    // The widened call set must contain the truth allele
    // whichever V was sampled (v01, v02, or v03).
    // Uses j02*01 (GGG prefix) so V's right-extension into J in the
    // VJ chain can't narrow the call back to the v01-alias group
    // (see `trim_widens_v_call_to_all_three_real_alleles` for the
    // rationale).
    let cfg = common::vj_ambiguous_refdata();
    let j02 = common::allele_id_by_name(&cfg, Segment::J, "j02*01");
    for truth_name in ["v01*01", "v02*01", "v03*01"] {
        let truth = common::allele_id_by_name(&cfg, Segment::V, truth_name);
        let plan = vj_plan(&cfg, truth, j02, 3, 0, 0, b'A');
        let outcome = run_vj(&cfg, plan);
        let sim = outcome.final_simulation();
        let names = common::v_call_names(sim, &cfg);
        assert!(
            names.iter().any(|n| n == truth_name),
            "truth allele {truth_name} must always be in v_call; got {names:?}",
        );
        // And the widened set is the full 4-allele set (because all
        // surviving bases are in the shared prefix).
        assert_eq!(
            names.len(),
            4,
            "trim-widened v_call should hold every pool entry sharing the prefix",
        );
    }
}

#[test]
fn score_based_narrowing_when_np_recreates_unique_suffix() {
    // Sample v02*01 (suffix AAA), trim V_3 by 3 → v_call widens
    // to all 4. NP1 = "AAA" then exactly recreates v02's
    // distinguishing suffix. Walker narrows v_call to {v02*01}
    // (v02 is the ONLY allele whose pos 9..12 is AAA).
    let cfg = common::vj_ambiguous_refdata();
    let v02 = common::allele_id_by_name(&cfg, Segment::V, "v02*01");
    let j01 = common::allele_id_by_name(&cfg, Segment::J, "j01*01");
    let plan = vj_plan(&cfg, v02, j01, 3, 0, 3, b'A');
    let outcome = run_vj(&cfg, plan);
    let sim = outcome.final_simulation();
    let names = common::v_call_names(sim, &cfg);
    assert_eq!(
        names,
        vec!["v02*01".to_string()],
        "NP1 = AAA matches v02's unique suffix; call should narrow",
    );
}

#[test]
fn v_sequence_end_reflects_live_hypothesis_after_trim() {
    // After V_3 trim 3 with no NP recovery, the V live hypothesis
    // ends at seq pos 12 (the post-trim structural V end: V length
    // 15 minus 3 = 12). AIRR `v_sequence_end` should reflect that,
    // not the pre-trim value of 15.
    let cfg = common::vj_ambiguous_refdata();
    let v01 = common::allele_id_by_name(&cfg, Segment::V, "v01*01");
    let j01 = common::allele_id_by_name(&cfg, Segment::J, "j01*01");
    let plan = vj_plan(&cfg, v01, j01, 3, 0, 0, b'A');
    let compiled = CompiledSimulator::compile(&plan, Some(&cfg), None, ExecutionPolicy::Permissive)
        .expect("plan should compile");
    let outcome = compiled
        .run_one(common::DEFAULT_SEED)
        .expect("plan should run");
    let rec = build_airr_record(&outcome, &cfg, "trim-test");
    assert_eq!(rec.v_trim_3, 3, "trim should record 3");
    assert_eq!(
        rec.v_sequence_end,
        Some(12),
        "v_sequence_end should track post-trim live hypothesis",
    );
}

#[test]
fn v_germline_end_shifts_inward_with_trim() {
    // The germline range tracks the surviving allele bases. A
    // 3-byte 3' trim with no NP recovery → v_germline_end is 12,
    // not 15 (V length minus the 3-byte trim).
    let cfg = common::vj_ambiguous_refdata();
    let v01 = common::allele_id_by_name(&cfg, Segment::V, "v01*01");
    let j01 = common::allele_id_by_name(&cfg, Segment::J, "j01*01");
    let plan = vj_plan(&cfg, v01, j01, 3, 0, 0, b'A');
    let compiled = CompiledSimulator::compile(&plan, Some(&cfg), None, ExecutionPolicy::Permissive)
        .expect("plan should compile");
    let outcome = compiled
        .run_one(common::DEFAULT_SEED)
        .expect("plan should run");
    let rec = build_airr_record(&outcome, &cfg, "trim-test");
    assert_eq!(
        rec.v_germline_end,
        Some(12),
        "trimmed V germline_end should shift inward to the surviving allele length",
    );
    assert_eq!(rec.v_germline_start, Some(0));
}

#[test]
fn trimmed_5_prime_j_widens_j_call_vj_chain() {
    // J pool: j01 = TTTAAACCC, j02 = GGGAAACCC. Trim J_5 by 3 →
    // surviving J bases are the shared AAACCC and both alleles
    // match equally → j_call = {j01*01, j02*01}.
    let cfg = common::vj_ambiguous_refdata();
    let v01 = common::allele_id_by_name(&cfg, Segment::V, "v01*01");
    let j01 = common::allele_id_by_name(&cfg, Segment::J, "j01*01");
    // NP1 = 'C' (single base) — NEITHER J prefix has 'C' at the
    // recoverable extension position, so no narrowing happens.
    // Pos checked first is ref_start - 1 = 2 (J*01[2]='T',
    // J*02[2]='G'). 'C' matches neither → no extension.
    let plan = vj_plan(&cfg, v01, j01, 0, 3, 1, b'C');
    let compiled = CompiledSimulator::compile(&plan, Some(&cfg), None, ExecutionPolicy::Permissive)
        .expect("plan should compile");
    let outcome = compiled
        .run_one(common::DEFAULT_SEED)
        .expect("plan should run");
    let sim = outcome.final_simulation();
    let names = common::j_call_names(sim, &cfg);
    assert_eq!(
        names,
        vec!["j01*01".to_string(), "j02*01".to_string()],
        "5' J trim should widen j_call to both J alleles",
    );

    // Also verify the AIRR projection records the trim and the
    // shifted germline_start.
    let rec = build_airr_record(&outcome, &cfg, "j-trim-test");
    assert_eq!(rec.j_trim_5, 3, "j_trim_5 should be recorded");
    assert_eq!(
        rec.j_germline_start,
        Some(3),
        "j_germline_start should shift inward to position 3 after 5' trim",
    );
}

#[test]
fn zero_trim_call_includes_aliases_but_excludes_distinguishable_alleles() {
    // Corrected from the original audit framing: zero-trim does NOT
    // widen the call to "the full pool" — it narrows to "the alias
    // group of truth". v01*01 and v01*02 are byte-identical, so both
    // appear; v02*01 (suffix AAA) and v03*01 (suffix CCC) are
    // distinguishable from v01 (suffix TTT) in the surviving bytes and
    // must NOT appear.
    let cfg = common::vj_ambiguous_refdata();
    let v01 = common::allele_id_by_name(&cfg, Segment::V, "v01*01");
    let j01 = common::allele_id_by_name(&cfg, Segment::J, "j01*01");
    let plan = vj_plan(&cfg, v01, j01, 0, 0, 0, b'A');
    let outcome = run_vj(&cfg, plan);
    let sim = outcome.final_simulation();
    let names = common::v_call_names(sim, &cfg);
    assert_eq!(
        names,
        vec!["v01*01".to_string(), "v01*02".to_string()],
        "zero trim narrows to the alias group of truth, NOT the full pool",
    );
    assert!(
        !names.contains(&"v02*01".to_string()),
        "distinguishable v02*01 must NOT be in call set under zero trim; got {:?}",
        names,
    );
    assert!(
        !names.contains(&"v03*01".to_string()),
        "distinguishable v03*01 must NOT be in call set under zero trim; got {:?}",
        names,
    );
}

#[test]
fn three_identical_aliases_all_appear_in_call_under_zero_trim() {
    // The triplet alias fixture has three byte-identical V alleles
    // (v_triplet*01..*03) plus a distinguishable v_unique*01. Sampling
    // any one of the triplet under zero trim must surface ALL THREE
    // triplet names in v_call.
    let cfg = common::vj_triplet_alias_refdata();
    let j_id = common::allele_id_by_name(&cfg, Segment::J, "j_single*01");
    for truth_name in ["v_triplet*01", "v_triplet*02", "v_triplet*03"] {
        let v_id = common::allele_id_by_name(&cfg, Segment::V, truth_name);
        let plan = vj_plan(&cfg, v_id, j_id, 0, 0, 0, b'A');
        let outcome = run_vj(&cfg, plan);
        let sim = outcome.final_simulation();
        let names = common::v_call_names(sim, &cfg);
        assert_eq!(
            names,
            vec![
                "v_triplet*01".to_string(),
                "v_triplet*02".to_string(),
                "v_triplet*03".to_string(),
            ],
            "triplet alleles must all appear when {} is truth; got {:?}",
            truth_name,
            names,
        );
        assert!(
            !names.contains(&"v_unique*01".to_string()),
            "v_unique*01 must drop out (distinguishable in tail); got {:?}",
            names,
        );
    }
}

#[test]
fn triplet_alias_set_narrows_correctly_when_distinguishable_allele_present() {
    // Sample v_triplet*01, no trim, in the pool that also contains
    // the distinguishable v_unique*01. The call must include all
    // three triplet aliases but NOT v_unique*01.
    let cfg = common::vj_triplet_alias_refdata();
    let v_id = common::allele_id_by_name(&cfg, Segment::V, "v_triplet*01");
    let j_id = common::allele_id_by_name(&cfg, Segment::J, "j_single*01");
    let plan = vj_plan(&cfg, v_id, j_id, 0, 0, 0, b'A');
    let outcome = run_vj(&cfg, plan);
    let sim = outcome.final_simulation();
    let names = common::v_call_names(sim, &cfg);
    assert_eq!(
        names,
        vec![
            "v_triplet*01".to_string(),
            "v_triplet*02".to_string(),
            "v_triplet*03".to_string(),
        ],
        "triplet aliases must all appear; v_unique must drop; got {:?}",
        names,
    );
}

#[test]
fn d_trim_both_sides_widens_d_call_vdj_chain() {
    // D pool: d01 = AAATTTGGG, d02 = AAATTTCCC (shares 6-byte
    // prefix with d01), d03 = CCCTTTGGG (shares 6-byte suffix
    // with d01). Trim D_5 by 3 (drops prefix) AND D_3 by 3 (drops
    // suffix) → surviving D bases are just TTT (3 bytes) and all
    // 3 D alleles match.
    //
    // NP1/NP2 = 'X' that matches nothing — single base 'X' is not
    // a valid base, use 'N'? We use 'A' but verify D pool's
    // pos 2 (left-extension first check) is C/G/A. d01[2]='A',
    // d02[2]='A', d03[2]='C'. NP1='A' matches d01 + d02 (drops d03).
    // To avoid that, use 'N'? Wait, simpler: use 'X' would crash.
    // Instead, let's use the empty NP (length 0) so no extension fires.
    let cfg = common::vdj_ambiguous_refdata();
    let v01 = common::allele_id_by_name(&cfg, Segment::V, "v01*01");
    let d01 = common::allele_id_by_name(&cfg, Segment::D, "d01*01");
    let j01 = common::allele_id_by_name(&cfg, Segment::J, "j01*01");
    let plan = vdj_plan(&cfg, v01, d01, j01, 0, 3, 3, 0, 0, b'A', 0, b'A');
    let compiled = CompiledSimulator::compile(&plan, Some(&cfg), None, ExecutionPolicy::Permissive)
        .expect("plan should compile");
    let outcome = compiled
        .run_one(common::DEFAULT_SEED)
        .expect("plan should run");
    let sim = outcome.final_simulation();
    let names = common::d_call_names(sim, &cfg);
    assert_eq!(
        names,
        vec![
            "d01*01".to_string(),
            "d02*01".to_string(),
            "d03*01".to_string(),
        ],
        "both-sides D trim leaves only the shared TTT → all 3 D alleles match",
    );

    let rec = build_airr_record(&outcome, &cfg, "d-trim-test");
    assert_eq!(rec.d_trim_5, 3);
    assert_eq!(rec.d_trim_3, 3);
    // d_germline_start/end shift to the surviving range [3, 6).
    assert_eq!(rec.d_germline_start, Some(3));
    assert_eq!(rec.d_germline_end, Some(6));
}
