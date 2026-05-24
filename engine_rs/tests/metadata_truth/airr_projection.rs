//! AIRR projection: live-call → AIRR record surfaces correctly.
//!
//! The live-call walker produces a `SegmentLiveCall` per V/D/J
//! segment containing the candidate allele set, hypothesis ranges,
//! and evidence score. The AIRR projection layer reads from this
//! and the structural region geometry to produce 69 named fields.
//! This module checks the projection path is faithful.
//!
//! ## Invariants this module should test
//!
//! - **`v_call` / `d_call` / `j_call` is the comma-joined sorted
//!   set of allele names from the live call** — verbatim, no
//!   reordering, no dedup beyond what the AlleleBitSet enforces.
//! - **AIRR call falls back to the origin assignment when no live
//!   call is staged** (test path / pre-compile path).
//! - **AIRR call falls back to truth when live call is
//!   `Unsupported`** (malformed IR or zero-score case).
//! - **Empty NP regions produce empty `np1` / `np2` strings**, not
//!   `null`. `np1_length == 0`, `np1_aa == ""`.
//! - **NP-claimed columns DON'T appear in `np1` / `np2`** strings
//!   (V/D claim takes precedence).
//! - **`junction_length == junction.len()` post-string-slice**,
//!   accounting for synthetic indel-inserted bytes.
//! - **`stop_codon` is true iff junction contains an in-frame
//!   TAA/TAG/TGA**, regardless of whether the record is productive.
//! - **`vj_in_frame` is true iff junction length is divisible by
//!   3**, regardless of stop codons.
//! - **`productive` is true iff vj_in_frame && !stop_codon &&
//!   anchors_preserved**.
//! - **`sequence_aa` translation is in junction frame** (not in
//!   V-region's pool frame — these are different when trim is not
//!   a multiple of 3; see the analysis).
//! - **`mutation_rate` == `n_mutations / sequence_length`** exactly
//!   (no float rounding drift past 6 digits).
//!
//! ## Adjacent in-tree coverage
//!
//! - `engine_rs/src/airr_record/tests/projection.rs`
//! - `engine_rs/src/airr_record/tests/anchors.rs`
//! - Golden TSV regressions: `tests/golden/seed42_n10_*.tsv`,
//!   `tests/golden/seed42_n20_human_igh_mixed_features.tsv`

use super::common;

use genairr_engine::airr_record::build_airr_record;
use genairr_engine::compiled::{CompiledSimulator, ExecutionPolicy};
use genairr_engine::dist::{AllelePoolDist, Distribution, EmpiricalLengthDist, UniformBase};
use genairr_engine::ir::Segment;
use genairr_engine::pass::PassPlan;
use genairr_engine::passes::{
    AssembleSegmentPass, GenerateNPPass, SampleAllelePass, UniformMutationPass,
};
use genairr_engine::refdata::{AlleleId, RefDataConfig};
use genairr_engine::rng::Rng;

/// Local clone of `compiled::tests::ConstBaseDist`: always samples
/// the configured byte. The in-tree version is private; this one is
/// local to the metadata-truth suite.
#[derive(Clone, Debug)]
struct ConstBaseDist(u8);

impl Distribution for ConstBaseDist {
    type Output = u8;
    fn sample(&self, _rng: &mut Rng) -> u8 {
        self.0
    }
    fn support(&self) -> Option<Vec<(u8, f64)>> {
        Some(vec![(self.0, 1.0)])
    }
}

// ──────────────────────────────────────────────────────────────────
// Local helpers
// ──────────────────────────────────────────────────────────────────

/// Compile + run + project to AIRR. Permissive policy throughout.
fn run_compiled(
    cfg: &RefDataConfig,
    plan: PassPlan,
    seed: u64,
) -> genairr_engine::airr_record::AirrRecord {
    let compiled =
        CompiledSimulator::compile(&plan, Some(cfg), None, ExecutionPolicy::Permissive)
            .expect("plan should compile");
    let outcome = compiled.run_one(seed).expect("plan should run");
    build_airr_record(&outcome, cfg, "proj")
}

/// Plain VJ plan that builds a zero-NP recombination with the
/// specified V/J truth alleles.
fn vj_recombine_plan(
    cfg: &RefDataConfig,
    v_id: AlleleId,
    j_id: AlleleId,
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
    plan.push(Box::new(AssembleSegmentPass::new(Segment::V)));
    plan.push(Box::new(GenerateNPPass::new(
        Segment::Np1,
        Box::new(EmpiricalLengthDist::from_pairs(vec![(np1_len, 1.0)])),
        Box::new(ConstBaseDist(np1_base)),
    )));
    plan.push(Box::new(AssembleSegmentPass::new(Segment::J)));
    plan
}

/// VDJ plan: sample given truth alleles, zero NP.
fn vdj_recombine_plan(
    cfg: &RefDataConfig,
    v_id: AlleleId,
    d_id: AlleleId,
    j_id: AlleleId,
    np1_len: i64,
    np2_len: i64,
    np_base: u8,
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
    plan.push(Box::new(AssembleSegmentPass::new(Segment::V)));
    plan.push(Box::new(GenerateNPPass::new(
        Segment::Np1,
        Box::new(EmpiricalLengthDist::from_pairs(vec![(np1_len, 1.0)])),
        Box::new(ConstBaseDist(np_base)),
    )));
    plan.push(Box::new(AssembleSegmentPass::new(Segment::D)));
    plan.push(Box::new(GenerateNPPass::new(
        Segment::Np2,
        Box::new(EmpiricalLengthDist::from_pairs(vec![(np2_len, 1.0)])),
        Box::new(ConstBaseDist(np_base)),
    )));
    plan.push(Box::new(AssembleSegmentPass::new(Segment::J)));
    plan
}

// ──────────────────────────────────────────────────────────────────
// Tests
// ──────────────────────────────────────────────────────────────────

#[test]
fn empty_np1_region_produces_empty_string_and_zero_length() {
    // Zero-length NP1: AIRR `np1` must be "" (not null / not missing)
    // and `np1_length == 0`, `np1_aa == ""`.
    let cfg = common::vdj_ambiguous_refdata();
    let v_id = common::allele_id_by_name(&cfg, Segment::V, "v01*01");
    let d_id = common::allele_id_by_name(&cfg, Segment::D, "d01*01");
    let j_id = common::allele_id_by_name(&cfg, Segment::J, "j01*01");

    // NP1 = 0, NP2 = 0.
    let plan = vdj_recombine_plan(&cfg, v_id, d_id, j_id, 0, 0, b'A');
    let rec = run_compiled(&cfg, plan, 0);

    assert_eq!(rec.np1, "");
    assert_eq!(rec.np1_length, 0);
    assert_eq!(rec.np1_aa, "");
    assert_eq!(rec.np2, "");
    assert_eq!(rec.np2_length, 0);
    assert_eq!(rec.np2_aa, "");
}

#[test]
fn np1_string_is_the_unclaimed_span_after_segment_extension() {
    // V1 = AAACCCGGGTTT, J1 = TTTAAACCC. Trim V3=0; no V suffix
    // claim; NP1 has 3 bases starting with 'A'. With NP1=3*A and no
    // V trim, V's right-extension walker compares NP1[0..3]="AAA"
    // against ref bytes at allele pos 12.. — V has only 12 bytes,
    // so the walker halts immediately (no extension). NP1 stays as
    // a 3-base unclaimed span.
    //
    // BUT: J's left-extension walker would also try to reach
    // backward into NP1. J1 starts with TTT, so the closest NP1 byte
    // (pos 11 = NP1's last 'A') would have to match J1[0-1]='T' to
    // extend. It doesn't, so no claim.
    let cfg = common::vj_ambiguous_refdata();
    let v_id = common::allele_id_by_name(&cfg, Segment::V, "v01*01");
    let j_id = common::allele_id_by_name(&cfg, Segment::J, "j01*01");

    // No trim, 3-base NP1 of 'A'.
    let plan = vj_recombine_plan(&cfg, v_id, j_id, 3, b'A');
    let rec = run_compiled(&cfg, plan, 0);

    // np1 is the full NP1 region (no segment claim).
    assert_eq!(rec.np1, "AAA");
    assert_eq!(rec.np1_length, 3);
}

#[test]
fn np1_string_drops_columns_claimed_by_v_extension() {
    // Trim V3=3 widens v_call to {v01*01, v01*02, v02*01, v03*01}.
    // NP1 = 'T' x 3 — only v01 (alleles 0 and 1) has T at ref
    // pos 9. The first NP1 byte narrows
    // {v01,v01*02,v02,v03}→{v01,v01*02}; the remaining 'T' bytes
    // match BOTH surviving aliases identically, so they cannot
    // narrow the tie set further. Under conservative extension
    // the walker stops after the single narrowing byte.
    //
    // 1 of 3 NP1 columns claimed → np1 = "TT" (2 chars).
    let cfg = common::vj_ambiguous_refdata();
    let v_id = common::allele_id_by_name(&cfg, Segment::V, "v01*01");
    let j_id = common::allele_id_by_name(&cfg, Segment::J, "j01*01");

    let mut plan = PassPlan::new();
    plan.push(Box::new(SampleAllelePass::new(
        Segment::V,
        Box::new(AllelePoolDist::restricted_uniform(&cfg.v_pool, vec![v_id])),
    )));
    plan.push(Box::new(SampleAllelePass::new(
        Segment::J,
        Box::new(AllelePoolDist::restricted_uniform(&cfg.j_pool, vec![j_id])),
    )));
    // Trim V_3 by 3.
    plan.push(Box::new(genairr_engine::passes::TrimPass::new(
        Segment::V,
        genairr_engine::assignment::TrimEnd::Three,
        Box::new(EmpiricalLengthDist::from_pairs(vec![(3, 1.0)])),
    )));
    plan.push(Box::new(AssembleSegmentPass::new(Segment::V)));
    plan.push(Box::new(GenerateNPPass::new(
        Segment::Np1,
        Box::new(EmpiricalLengthDist::from_pairs(vec![(3, 1.0)])),
        Box::new(ConstBaseDist(b'T')),
    )));
    plan.push(Box::new(AssembleSegmentPass::new(Segment::J)));

    let rec = run_compiled(&cfg, plan, 0);

    // V claimed the single narrowing NP1 column; 2 remain.
    assert_eq!(rec.np1, "TT");
    assert_eq!(rec.np1_length, 2);
    // v_call narrowed to the v01 aliases.
    let mut v_calls: Vec<&str> = rec.v_call.split(',').collect();
    v_calls.sort();
    assert_eq!(v_calls, vec!["v01*01", "v01*02"]);
}

#[test]
fn vj_in_frame_reflects_junction_length_mod_three() {
    // The junction is V[anchor..end] + NP1 + J[start..anchor+3].
    // For the VJ fixture: V anchor=9, V len=12 → V tail 3 bytes.
    // J anchor=0, J len=9 → J head 3 bytes.
    // Total without NP = 6 bytes (in frame).
    // Add NP1=1 → 7 bytes (out of frame).
    // Add NP1=3 → 9 bytes (in frame).
    let cfg = common::vj_ambiguous_refdata();
    let v_id = common::allele_id_by_name(&cfg, Segment::V, "v01*01");
    let j_id = common::allele_id_by_name(&cfg, Segment::J, "j01*01");

    // 0-base NP → 6-byte junction → in frame.
    let plan = vj_recombine_plan(&cfg, v_id, j_id, 0, b'A');
    let rec_a = run_compiled(&cfg, plan, 0);
    assert_eq!(rec_a.vj_in_frame, Some(true));
    assert_eq!(rec_a.junction_length.unwrap() % 3, 0);

    // 1-base NP → 7-byte junction → out of frame.
    let plan = vj_recombine_plan(&cfg, v_id, j_id, 1, b'A');
    let rec_b = run_compiled(&cfg, plan, 0);
    assert_eq!(rec_b.vj_in_frame, Some(false));
    assert!(rec_b.junction_length.unwrap() % 3 != 0);
    // Out-of-frame → productive must be false.
    assert_eq!(rec_b.productive, Some(false));
}

#[test]
fn junction_length_equals_junction_string_byte_count() {
    // `junction_length` is set to `junction_nt.len() as i64` after
    // slicing, so `junction_length == rec.junction.len()` by
    // construction. Verify across a few configurations.
    let cfg = common::vj_ambiguous_refdata();
    let v_id = common::allele_id_by_name(&cfg, Segment::V, "v01*01");
    let j_id = common::allele_id_by_name(&cfg, Segment::J, "j01*01");

    for np1_len in 0..6i64 {
        let plan = vj_recombine_plan(&cfg, v_id, j_id, np1_len, b'A');
        let rec = run_compiled(&cfg, plan, 0);
        assert_eq!(
            rec.junction_length.unwrap() as usize,
            rec.junction.len(),
            "np1_len={} : junction_length / junction.len() mismatch",
            np1_len,
        );
    }
}

#[test]
fn stop_codon_flag_matches_actual_codon_scan_of_junction() {
    // The `stop_codon` flag must be true iff the junction nt
    // contains an in-frame TAA / TAG / TGA. Drive an NP that forces
    // a stop: V1 tail = "TGT" + NP1 + J1 head = "TTT" + remaining.
    // With NP1 = "AAT" (3 bases): junction = "TGT" + "AAT" + "TTT"
    // = TGT AAT TTT. No stop. To get a stop, NP1 needs to make a
    // TAA/TAG/TGA — set NP1 = "AA" (1 base — out of frame though)
    // or align it so a codon falls on a stop.
    //
    // Concrete trick: J1 starts with "TTT" but to get a TAA codon,
    // we need an NP that places TAA at a codon boundary inside the
    // junction. Junction starts at V anchor (TGT). Codons are
    // TGT|...|...|TTT (or extended). We can force a stop by making
    // NP1 = "TAA" (3 bytes): codons are TGT|TAA|TTT → TAA at codon
    // 2 → has_stop = true.
    let cfg = common::vj_ambiguous_refdata();
    let v_id = common::allele_id_by_name(&cfg, Segment::V, "v01*01");
    let j_id = common::allele_id_by_name(&cfg, Segment::J, "j01*01");

    // We have to build NP1 with mixed bases — use a per-base sample
    // pass with `ConstBaseDist` isn't possible (we use one base per
    // NP). Easier: NP1 length 6 of "AAA" → junction codons are
    // V_tail="TGT" + NP1=AAAAAA + J_head="TTT" → "TGT AAA AAA TTT" →
    // K, K, F — no stop. Try NP1 length 6 with mix? We can only emit
    // one base per ConstBaseDist. So pick NP1=6 of 'A': codons all
    // K/F. No stop. For a stop we need TAA/TAG/TGA, which requires a
    // T plus AA/AG/GA. Easier path: junction always contains "TGT"
    // (V anchor) — that's just C. And "TGG" (J anchor for VDJ
    // chains; for VJ J1 we use J1's start = "TTT"). So TGT|...|TTT.
    // The middle (NP1) is fixed-base. With NP1=AA (2 bases) →
    // out-of-frame. With NP1=3*A → TGT|AAA|TTT → no stop.
    //
    // Since the fixture restricts us to a single base per NP, we
    // cover the no-stop case cleanly and verify a separate
    // contradiction: a deliberately constructed junction with NP1
    // = 3*'T' → TGT|TTT|TTT codons → F F F → no stop.
    //
    // For an actual stop, we must use a DIFFERENT route: trim J5
    // so J head doesn't start with TTT but with whatever survives.
    // J1 = TTTAAACCC, anchor=0; trim J5=0 keeps J head as the first
    // 3 bytes = TTT. We can't easily force a stop with just one
    // ConstBaseDist for NP.
    //
    // INSTEAD: verify the contrapositive. If the junction has NO
    // in-frame stop codon, `stop_codon` must be Some(false). Build
    // the baseline (NP1=3*'A') and check.
    let plan = vj_recombine_plan(&cfg, v_id, j_id, 3, b'A');
    let rec = run_compiled(&cfg, plan, 0);
    assert_eq!(rec.vj_in_frame, Some(true));
    assert_eq!(rec.stop_codon, Some(false));
    assert!(!rec.junction.is_empty());
    // Verify manually by scanning the junction.
    let junction_bytes = rec.junction.as_bytes();
    let mut found_stop = false;
    for chunk in junction_bytes.chunks_exact(3) {
        let upper = [
            chunk[0].to_ascii_uppercase(),
            chunk[1].to_ascii_uppercase(),
            chunk[2].to_ascii_uppercase(),
        ];
        if upper == *b"TAA" || upper == *b"TAG" || upper == *b"TGA" {
            found_stop = true;
            break;
        }
    }
    assert_eq!(rec.stop_codon, Some(found_stop));
}

#[test]
fn productive_equals_in_frame_and_no_stop_and_anchors_preserved() {
    // For a clean baseline (no mutations, no trim, no NP that
    // destroys anchors): productive == (vj_in_frame && !stop_codon).
    let cfg = common::vj_ambiguous_refdata();
    let v_id = common::allele_id_by_name(&cfg, Segment::V, "v01*01");
    let j_id = common::allele_id_by_name(&cfg, Segment::J, "j01*01");

    // Productive in-frame baseline (NP1=3*A).
    let plan = vj_recombine_plan(&cfg, v_id, j_id, 3, b'A');
    let rec = run_compiled(&cfg, plan, 0);
    assert_eq!(rec.vj_in_frame, Some(true));
    assert_eq!(rec.stop_codon, Some(false));
    assert_eq!(rec.productive, Some(true));

    // Out-of-frame: productive forced to false.
    let plan = vj_recombine_plan(&cfg, v_id, j_id, 1, b'A');
    let rec = run_compiled(&cfg, plan, 0);
    assert_eq!(rec.vj_in_frame, Some(false));
    assert_eq!(rec.productive, Some(false));
}

#[test]
fn mutation_rate_equals_n_mutations_over_sequence_length_exactly() {
    // mutation_rate is computed as `rec.n_mutations as f64 /
    // rec.sequence_length as f64`. Verify the value is exactly the
    // float division of those two integers (no rounding past f64).
    let cfg = common::vdj_ambiguous_refdata();
    let v_id = common::allele_id_by_name(&cfg, Segment::V, "v01*01");
    let d_id = common::allele_id_by_name(&cfg, Segment::D, "d01*01");
    let j_id = common::allele_id_by_name(&cfg, Segment::J, "j01*01");

    let mut plan = vdj_recombine_plan(&cfg, v_id, d_id, j_id, 0, 0, b'A');
    plan.push(Box::new(UniformMutationPass::new(
        Box::new(EmpiricalLengthDist::from_pairs(vec![(5, 1.0)])),
        Box::new(UniformBase),
    )));

    let rec = run_compiled(&cfg, plan, 11);
    let expected = rec.n_mutations as f64 / rec.sequence_length as f64;
    assert_eq!(rec.mutation_rate, expected);
    // And sanity check that the value is in (0, 1].
    assert!(rec.mutation_rate > 0.0 && rec.mutation_rate <= 1.0);
}

#[test]
fn n_mutations_reflects_live_call_state_mutation_count() {
    // AIRR `n_mutations` reads from `LiveCallState.mutation_count`
    //. Verify the AIRR field matches
    // the IR sidecar after a clean recombine + mutate run.
    let cfg = common::vdj_ambiguous_refdata();
    let v_id = common::allele_id_by_name(&cfg, Segment::V, "v01*01");
    let d_id = common::allele_id_by_name(&cfg, Segment::D, "d01*01");
    let j_id = common::allele_id_by_name(&cfg, Segment::J, "j01*01");

    let mut plan = vdj_recombine_plan(&cfg, v_id, d_id, j_id, 0, 0, b'A');
    plan.push(Box::new(UniformMutationPass::new(
        Box::new(EmpiricalLengthDist::from_pairs(vec![(4, 1.0)])),
        Box::new(UniformBase),
    )));

    let compiled =
        CompiledSimulator::compile(&plan, Some(&cfg), None, ExecutionPolicy::Permissive)
            .expect("plan should compile");
    let outcome = compiled.run_one(7).expect("plan should run");
    let rec = build_airr_record(&outcome, &cfg, "n_mut");

    let live_count = outcome.final_simulation().mutation_count as i64;

    assert_eq!(rec.n_mutations, live_count);
    assert_eq!(rec.n_mutations, 4);
}

#[test]
fn v_call_falls_back_to_origin_when_no_live_call_is_present() {
    // When the pipeline runs via `PassRuntime::execute_with_refdata`
    // (no compiled flow → no walker observers attached → no live
    // calls populated), the AIRR `v_call` must fall back to the
    // origin assignment's allele name.
    //
    // Here we use PassRuntime::execute (not compiled) so live_calls
    // stays None.
    let cfg = common::vj_ambiguous_refdata();
    let v_id = common::allele_id_by_name(&cfg, Segment::V, "v01*01");
    let j_id = common::allele_id_by_name(&cfg, Segment::J, "j01*01");

    let plan = vj_recombine_plan(&cfg, v_id, j_id, 0, b'A');
    let outcome = genairr_engine::pass::PassRuntime::execute_with_refdata(
        &plan,
        genairr_engine::ir::Simulation::new(),
        0,
        &cfg,
    );
    let rec = build_airr_record(&outcome, &cfg, "fallback");

    // segment_calls is empty (no walker observer ran) — the v_call
    // must come from the origin.
    assert_eq!(
        outcome.final_simulation().segment_calls.version,
        0,
        "PassRuntime::execute_with_refdata should NOT populate segment_calls",
    );
    assert_eq!(rec.v_call, "v01*01");
    assert_eq!(rec.j_call, "j01*01");
}

#[test]
fn v_call_listing_orders_truth_first_then_remaining_alleles() {
    // When the live call sets a tie-set {V01*01, V01*02, V02*01, V03*01}
    // (e.g. trim V3=3), AIRR `v_call` lists the truth allele FIRST,
    // then the other alleles in iter_ids() ascending order.
    // (See `live_call_name` in airr_record/projection.rs.)
    //
    // Uses j02*01 (GGG prefix) to keep the v_call tie set at 4 in
    // the VJ chain. wired V's right-extension to walk into
    // J's bytes; with j01 (TTT prefix) v01's continuation would
    // narrow the call back to {v01*01, v01*02} and break this
    // ordering test. GGG matches no V allele at the boundary, so
    // the widened set survives intact.
    let cfg = common::vj_ambiguous_refdata();
    let v_id = common::allele_id_by_name(&cfg, Segment::V, "v02*01");
    let j_id = common::allele_id_by_name(&cfg, Segment::J, "j02*01");

    let mut plan = PassPlan::new();
    plan.push(Box::new(SampleAllelePass::new(
        Segment::V,
        Box::new(AllelePoolDist::restricted_uniform(&cfg.v_pool, vec![v_id])),
    )));
    plan.push(Box::new(SampleAllelePass::new(
        Segment::J,
        Box::new(AllelePoolDist::restricted_uniform(&cfg.j_pool, vec![j_id])),
    )));
    plan.push(Box::new(genairr_engine::passes::TrimPass::new(
        Segment::V,
        genairr_engine::assignment::TrimEnd::Three,
        Box::new(EmpiricalLengthDist::from_pairs(vec![(3, 1.0)])),
    )));
    plan.push(Box::new(AssembleSegmentPass::new(Segment::V)));
    plan.push(Box::new(GenerateNPPass::new(
        Segment::Np1,
        Box::new(EmpiricalLengthDist::from_pairs(vec![(0, 1.0)])),
        Box::new(UniformBase),
    )));
    plan.push(Box::new(AssembleSegmentPass::new(Segment::J)));

    let rec = run_compiled(&cfg, plan, 0);

    // Truth (v02*01) is FIRST in the call list; rest follow iter_ids
    // order (which is ascending pool-index, i.e. v01*01 (id 0),
    // v01*02 (id 1), v03*01 (id 3); v02*01 already listed first).
    let v_calls: Vec<&str> = rec.v_call.split(',').collect();
    assert_eq!(v_calls.first(), Some(&"v02*01"));
    // The set must contain all 4 alleles (ambiguity after V3=3 trim).
    let mut sorted = v_calls.clone();
    sorted.sort();
    assert_eq!(sorted, vec!["v01*01", "v01*02", "v02*01", "v03*01"]);
}
