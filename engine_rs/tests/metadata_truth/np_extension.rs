//! NP region extending the V/D/J call back into trimmed bytes.
//!
//! After 5'/3' trimming widens a V/D/J call to a set of
//! indistinguishable alleles, the adjacent NP region's random bases
//! may "recreate" enough of the trimmed germline that the call
//! narrows back. The engine must propagate this NP-extension event
//! into both the call set AND the structural coordinates:
//!
//! - `v_sequence_end` extends into NP1 by the number of recreated
//!   bytes.
//! - `np1` (the AIRR string) drops those columns from the start.
//! - `v_germline_end` extends in the allele reference by the same
//!   count.
//! - `v_cigar` reflects the extended M-run.
//! - Equivalent invariants for D (NP1 left + NP2 right) and J (NP1
//!   left in VJ; NP2 left in VDJ).
//!
//! ## Invariants this module should test
//!
//! - **Exact recreation narrows call to single allele** (V01 was the
//!   truth, trim removed 3 bytes, NP1 first 3 bytes happen to be the
//!   V01 trimmed bytes verbatim — call narrows to just V01).
//! - **Partial recreation narrows partially** (NP1 recreates only
//!   the first 2 of 3 trimmed bytes — call narrows to alleles that
//!   match those 2, stays widened for the third).
//! - **No recreation leaves call widened** (NP1 bytes don't match
//!   any candidate's trimmed suffix).
//! - **NP extension crosses INTO the next region** (V extends past
//!   the entire NP1 region into D's bytes — D's bytes form a suffix
//!   of a V allele).
//! - **NP1 / NP2 strings shrink correctly** when V/D/J claim
//!   adjacent NP columns.
//! - **Bidirectional D extension** (D's left edge extends into NP1
//!   AND D's right edge extends into NP2 simultaneously).
//! - **VJ chain J extension into NP1** (no D region, NP1 is V→J
//!   spacer, J's 5' extends back into it).
//! - **Frame-phase invariance**: the extended region's frame_phase
//!   is unchanged (extension doesn't move the V anchor frame).
//!
//! ## Adjacent in-tree coverage
//!
//! - `engine_rs/src/compiled/tests/v_extensions.rs`
//! - `engine_rs/src/compiled/tests/d_extensions.rs`
//! - `engine_rs/src/compiled/tests/j_extensions.rs`
//! - `engine_rs/src/compiled/tests/curated.rs::curated_np1_recreates_v_suffix_narrows_v_call_back`
//! - `engine_rs/src/compiled/tests/curated.rs::v_germline_end_reflects_np_extension`
//! - `engine_rs/src/compiled/tests/curated.rs::v_cigar_extends_into_claimed_np1_columns`
//! - `engine_rs/src/compiled/tests/curated.rs::np1_string_drops_columns_claimed_by_v`

use super::common;
use genairr_engine::airr_record::build_airr_record;
use genairr_engine::assignment::TrimEnd;
use genairr_engine::compiled::{CompiledSimulator, ExecutionPolicy};
use genairr_engine::dist::{AllelePoolDist, Distribution, EmpiricalLengthDist};
use genairr_engine::ir::Segment;
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

// VJ plan with optional V_3 trim and configurable NP1 content.
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

// ──────────────────────────────────────────────────────────────
// Tests
// ──────────────────────────────────────────────────────────────

#[test]
fn exact_recreation_narrows_v_call_to_single_truth_allele() {
    // Sample v03*01 (distinguishing suffix CCC). Trim V_3 by 3 →
    // v_call widens to all 4. NP1 = "CCC". Under conservative
    // extension, the first NP1 byte 'C' is enough to
    // narrow {v01,v01*02,v02,v03}→{v03} (only v03 has C at pos 9).
    // The remaining two 'C' bytes cannot narrow {v03} further, so
    // the walker halts. v_call collapses to {v03*01} after a
    // single byte of extension; ref_end advances by 1 only.
    let cfg = common::vj_ambiguous_refdata();
    let v03 = common::allele_id_by_name(&cfg, Segment::V, "v03*01");
    let j01 = common::allele_id_by_name(&cfg, Segment::J, "j01*01");
    let plan = vj_plan(&cfg, v03, j01, 3, 0, 3, b'C');
    let outcome = compile_and_run(&cfg, &plan);
    let sim = outcome.final_simulation();
    let names = common::v_call_names(sim, &cfg);
    assert_eq!(
        names,
        vec!["v03*01".to_string()],
        "NP1=CCC recreates v03's unique suffix; call should narrow to v03 only",
    );

    // AIRR projection: v_germline_end extends by 1 byte (the
    // narrowing step) under conservative extension, from
    // structural 9 to 10.
    let rec = build_airr_record(&outcome, &cfg, "exact-recreation");
    assert_eq!(rec.v_germline_end, Some(10));
    assert_eq!(rec.v_sequence_end, Some(10));
    // NP1 string drops only the single claimed column; the other
    // two 'C' bytes stay non-templated.
    assert_eq!(rec.np1, "CC");
    assert_eq!(rec.np1_length, 2);
    // CIGAR covers 9 structural + 1 NP1 = 10M.
    assert_eq!(rec.v_cigar, "10M");
}

#[test]
fn partial_recreation_narrows_partway() {
    // Sample v03 (distinguishing suffix CCC). Trim V_3 by 3,
    // NP1 length 5 of 'C'. Under conservative extension, the
    // first 'C' narrows {v01,v01*02,v02,v03}→{v03}; subsequent
    // bytes cannot narrow further so the walker stops after 1
    // byte. v_call narrows to {v03}, ref_end advances 9→10,
    // and 4 of the 5 NP1 bytes stay in the np1 string.
    let cfg = common::vj_ambiguous_refdata();
    let v03 = common::allele_id_by_name(&cfg, Segment::V, "v03*01");
    let j01 = common::allele_id_by_name(&cfg, Segment::J, "j01*01");
    let plan = vj_plan(&cfg, v03, j01, 3, 0, 5, b'C');
    let outcome = compile_and_run(&cfg, &plan);
    let sim = outcome.final_simulation();
    let names = common::v_call_names(sim, &cfg);
    assert_eq!(names, vec!["v03*01".to_string()]);

    let rec = build_airr_record(&outcome, &cfg, "partial-recreation");
    assert_eq!(rec.v_germline_end, Some(10));
    // V claimed 1 of the 5 NP1 bases; remaining 4 stay in np1.
    assert_eq!(rec.np1_length, 4);
    assert_eq!(rec.np1, "CCCC");
}

#[test]
fn no_recreation_leaves_v_call_widened() {
    // Sample v01, trim V_3 by 3, NP1 = "X" — use a base that
    // matches NO V allele's pos 9. v01[9]='T', v02[9]='A',
    // v03[9]='C'. 'G' matches none → no extension; v_call stays
    // at the 4-allele widened set.
    let cfg = common::vj_ambiguous_refdata();
    let v01 = common::allele_id_by_name(&cfg, Segment::V, "v01*01");
    let j01 = common::allele_id_by_name(&cfg, Segment::J, "j01*01");
    let plan = vj_plan(&cfg, v01, j01, 3, 0, 3, b'G');
    let outcome = compile_and_run(&cfg, &plan);
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
        "NP1=G matches no V allele suffix; call must remain widened",
    );

    let rec = build_airr_record(&outcome, &cfg, "no-recreation");
    assert_eq!(
        rec.v_germline_end,
        Some(9),
        "no extension → germline_end at structural post-trim"
    );
    assert_eq!(rec.np1_length, 3, "no NP1 columns claimed");
}

#[test]
fn bidirectional_d_extension_recreates_both_sides() {
    // VDJ: trim D_5 by 3 AND D_3 by 3 → D structural region is
    // just TTT (the shared core). NP1 = "AAA" (3 bytes) and
    // NP2 = "GGG" (3 bytes).
    //
    // D pool: d01=AAATTTGGG, d02=AAATTTCCC, d03=CCCTTTGGG.
    //
    // NP1 left walk (against ref pos 2,1,0):
    //   pre tie set after assembly = {d01,d02,d03} (all match TTT)
    //   step 1 at ref pos 2: d01[2]='A' ✓, d02[2]='A' ✓,
    //     d03[2]='C' ✗ → narrows {d01,d02,d03}→{d01,d02}. Extends.
    //   step 2 at ref pos 1: tie set is {d01,d02}, both match 'A'.
    //     No "any_missed" → conservative declines. Walker stops.
    //   ref_start moves from 3 to 2 (one byte).
    //
    // NP2 right walk (against ref pos 6,7,8):
    //   pre tie set after NP1 extension = {d01,d02} (d03 dropped)
    //   step 1 at ref pos 6: d01[6]='G' ✓, d02[6]='C' ✗ → narrows
    //     {d01,d02}→{d01}. Extends.
    //   step 2: tie {d01}; conservative declines. Walker stops.
    //   ref_end moves from 6 to 7 (one byte).
    let cfg = common::vdj_ambiguous_refdata();
    let v01 = common::allele_id_by_name(&cfg, Segment::V, "v01*01");
    let d01 = common::allele_id_by_name(&cfg, Segment::D, "d01*01");
    let j01 = common::allele_id_by_name(&cfg, Segment::J, "j01*01");
    let plan = vdj_plan(&cfg, v01, d01, j01, 0, 3, 3, 0, 3, b'A', 3, b'G');
    let outcome = compile_and_run(&cfg, &plan);
    let sim = outcome.final_simulation();
    let names = common::d_call_names(sim, &cfg);
    assert_eq!(
        names,
        vec!["d01*01".to_string()],
        "bidirectional NP extension narrows d_call to {{d01}}",
    );

    let rec = build_airr_record(&outcome, &cfg, "bidirectional-d");
    // Each side recovers 1 byte under conservative extension.
    assert_eq!(rec.d_germline_start, Some(2));
    assert_eq!(rec.d_germline_end, Some(7));
    // Each NP region: 1 byte claimed, 2 remaining.
    assert_eq!(rec.np1, "AA");
    assert_eq!(rec.np2, "GG");
    assert_eq!(rec.np1_length, 2);
    assert_eq!(rec.np2_length, 2);
}

#[test]
fn v_extends_past_empty_np1_into_d_bytes() {
    // Build a small one-off refdata where V01's distinguishing
    // suffix is "TTT" and D01's first 3 bytes are also "TTT".
    // Trim V_3 by 3, NP1 length 0 — V's right extension walker
    // crosses the empty NP1 and reads D's leading bytes directly,
    // narrowing v_call to {V01}.
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
        seq: b"TTTAAGCG".to_vec(),
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
    let plan = vdj_plan(&cfg, v01, d01, j01, 3, 0, 0, 0, 0, b'A', 0, b'A');
    let outcome = compile_and_run(&cfg, &plan);
    let sim = outcome.final_simulation();
    let names = common::v_call_names(sim, &cfg);
    assert_eq!(
        names,
        vec!["V*01".to_string()],
        "V should extend across empty NP1 into D's TTT prefix",
    );
}

#[test]
fn np1_string_shrinks_when_v_claims_bytes() {
    // Sample v01, V_3 trim 3, NP1 length 6, base 'T'. v01[9]='T',
    // v01*02[9]='T', v02[9]='A', v03[9]='C'. The first NP1 byte
    // at ref pos 9 narrows {v01,v01*02,v02,v03}→{v01,v01*02}
    // (two indistinguishable aliases remain). The second byte at
    // ref pos 10 cannot narrow further (both surviving alleles
    // match T at pos 10), so under conservative extension the
    // walker stops after just 1 byte.
    let cfg = common::vj_ambiguous_refdata();
    let v01 = common::allele_id_by_name(&cfg, Segment::V, "v01*01");
    let j01 = common::allele_id_by_name(&cfg, Segment::J, "j01*01");
    let plan = vj_plan(&cfg, v01, j01, 3, 0, 6, b'T');
    let outcome = compile_and_run(&cfg, &plan);
    let rec = build_airr_record(&outcome, &cfg, "np1-shrink");
    assert_eq!(
        rec.np1_length, 5,
        "1 of the 6 NP1 bytes is claimed by V; np1_length should be 5",
    );
    assert_eq!(
        rec.np1, "TTTTT",
        "remaining NP1 string holds the un-claimed bytes",
    );
}

#[test]
fn vj_chain_j_extension_into_np1() {
    // VJ chain: J's left neighbour is NP1 (no D). Sample j01,
    // trim J_5 by 3 → j_call widens to {j01, j02}. NP1 = "TTT"
    // (j01's distinguishing prefix). Walker checks right-to-left:
    //   NP1[2]='T' vs j01[2]='T' ✓, j02[2]='G' ✗ → narrows the
    //     tie set {j01,j02}→{j01}. Extends.
    //   NP1[1]: tie {j01}; conservative declines. Walker stops.
    // j_call narrows to {j01*01} via a single byte of extension;
    // ref_start moves from 3 to 2, not all the way to 0.
    let cfg = common::vj_ambiguous_refdata();
    let v01 = common::allele_id_by_name(&cfg, Segment::V, "v01*01");
    let j01 = common::allele_id_by_name(&cfg, Segment::J, "j01*01");
    let plan = vj_plan(&cfg, v01, j01, 0, 3, 3, b'T');
    let outcome = compile_and_run(&cfg, &plan);
    let sim = outcome.final_simulation();
    let names = common::j_call_names(sim, &cfg);
    assert_eq!(
        names,
        vec!["j01*01".to_string()],
        "VJ chain: NP1=TTT recreates j01's trimmed prefix; call narrows",
    );

    let rec = build_airr_record(&outcome, &cfg, "vj-j-np1");
    assert_eq!(
        rec.j_germline_start,
        Some(2),
        "single byte of left extension"
    );
    // Only 1 NP1 column absorbed; remaining 2 stay non-templated.
    assert_eq!(rec.np1, "TT");
    assert_eq!(rec.np1_length, 2);
}

#[test]
fn vdj_chain_j_extension_into_np2() {
    // VDJ chain: J's left neighbour is NP2. Trim J_5 by 3, NP2 =
    // "TTT" — narrows j_call to {j01*01} via left-extension.
    // Under conservative extension only the first NP2 byte (which
    // narrows {j01,j02}→{j01}) is claimed; subsequent bytes stay
    // non-templated.
    let cfg = common::vdj_ambiguous_refdata();
    let v01 = common::allele_id_by_name(&cfg, Segment::V, "v01*01");
    let d01 = common::allele_id_by_name(&cfg, Segment::D, "d01*01");
    let j01 = common::allele_id_by_name(&cfg, Segment::J, "j01*01");
    let plan = vdj_plan(&cfg, v01, d01, j01, 0, 0, 0, 3, 0, b'A', 3, b'T');
    let outcome = compile_and_run(&cfg, &plan);
    let sim = outcome.final_simulation();
    let names = common::j_call_names(sim, &cfg);
    assert_eq!(
        names,
        vec!["j01*01".to_string()],
        "VDJ chain: NP2=TTT recreates j01's prefix; call narrows",
    );

    let rec = build_airr_record(&outcome, &cfg, "vdj-j-np2");
    assert_eq!(rec.j_germline_start, Some(2));
    assert_eq!(rec.np2, "TT");
}

#[test]
fn v_cigar_extends_across_claimed_np1_columns() {
    // With V_3 trim 3 and NP1=TTT against v01 (suffix TTT,
    // aliased by v01*02), the first NP1 byte at ref pos 9
    // narrows {v01,v01*02,v02,v03}→{v01,v01*02}; subsequent
    // 'T' bytes match both surviving alleles and cannot
    // narrow the tie set further. Under conservative extension
    // the walker stops after 1 byte → v_cigar = 10M (9
    // structural + 1 NP1 column claimed).
    let cfg = common::vj_ambiguous_refdata();
    let v01 = common::allele_id_by_name(&cfg, Segment::V, "v01*01");
    let j01 = common::allele_id_by_name(&cfg, Segment::J, "j01*01");
    let plan = vj_plan(&cfg, v01, j01, 3, 0, 3, b'T');
    let outcome = compile_and_run(&cfg, &plan);
    let rec = build_airr_record(&outcome, &cfg, "cigar-extend");
    assert_eq!(rec.v_cigar, "10M");
}
