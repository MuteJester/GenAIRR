//! NP evidence can narrow live re-calls without changing generation truth.
//!
//! Matching adjacent junction bases may extend an inferred placement. Standard
//! AIRR calls and alignment bounds retain the generating assignments, and all
//! generated NP bytes remain in their actual junction intervals. These tests
//! check the live inference and the independent standard projection together.

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
fn exact_recreation_narrows_recall_without_reassigning_np() {
    // NP1=CCC identifies v03 from the trimmed V tie set; all three
    // bases remain generated junction bases in the standard record.
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

    let rec = build_airr_record(&outcome, &cfg, "exact-recreation");
    assert_eq!(rec.v_germline_end, Some(12));
    assert_eq!(rec.v_sequence_end, Some(12));
    assert_eq!(rec.np1, "CCC");
    assert_eq!(rec.np1_length, 3);
    assert_eq!(rec.v_cigar, "12M");
}

#[test]
fn partial_recreation_narrows_recall_without_shortening_np() {
    // NP1 narrows the live call after one matching base. Standard
    // generation truth keeps the full five-base junction.
    let cfg = common::vj_ambiguous_refdata();
    let v03 = common::allele_id_by_name(&cfg, Segment::V, "v03*01");
    let j01 = common::allele_id_by_name(&cfg, Segment::J, "j01*01");
    let plan = vj_plan(&cfg, v03, j01, 3, 0, 5, b'C');
    let outcome = compile_and_run(&cfg, &plan);
    let sim = outcome.final_simulation();
    let names = common::v_call_names(sim, &cfg);
    assert_eq!(names, vec!["v03*01".to_string()]);

    let rec = build_airr_record(&outcome, &cfg, "partial-recreation");
    assert_eq!(rec.v_germline_end, Some(12));
    assert_eq!(rec.np1_length, 5);
    assert_eq!(rec.np1, "CCCCC");
}

#[test]
fn no_recreation_leaves_v_call_widened() {
    // Sample v01, trim V_3 by 3, NP1 = "T" — the one base that
    // matches NO V allele's first distinguishing-tail byte at
    // pos 12. v01[12]='G', v02[12]='A', v03[12]='C'. 'T' matches
    // none → no extension; v_call stays at the 4-allele widened
    // set. (J's first byte is also 'T', so J-side right-extension
    // halts immediately too.)
    let cfg = common::vj_ambiguous_refdata();
    let v01 = common::allele_id_by_name(&cfg, Segment::V, "v01*01");
    let j01 = common::allele_id_by_name(&cfg, Segment::J, "j01*01");
    let plan = vj_plan(&cfg, v01, j01, 3, 0, 3, b'T');
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
        "NP1=T matches no V allele suffix; call must remain widened",
    );

    let rec = build_airr_record(&outcome, &cfg, "no-recreation");
    assert_eq!(
        rec.v_germline_end,
        Some(12),
        "no extension → germline_end at structural post-trim"
    );
    assert_eq!(rec.np1_length, 3, "no NP1 columns claimed");
}

#[test]
fn bidirectional_d_recall_extension_preserves_both_generated_junctions() {
    // NP1 and NP2 together resolve the D re-call. Standard D bounds
    // remain the retained reference interval 3..6 and both NPs remain length 3.
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
    assert_eq!(rec.d_germline_start, Some(3));
    assert_eq!(rec.d_germline_end, Some(6));
    assert_eq!(rec.np1, "AAA");
    assert_eq!(rec.np2, "GGG");
    assert_eq!(rec.np1_length, 3);
    assert_eq!(rec.np2_length, 3);
}

#[test]
fn v_extends_past_empty_np1_into_d_bytes() {
    // Build a small one-off refdata where V01's distinguishing
    // suffix is "TTT" and D01's first 3 bytes are also "TTT".
    // Trim V_3 by 3, NP1 length 0 — V's right extension walker
    // crosses the empty NP1 and reads D's leading bytes directly,
    // narrowing v_call to {V01}.
    let mut cfg = RefDataConfig::empty(ChainType::Vdj);
    // Test uses intentionally anchorless V/J alleles; relax
    // the default required-anchor rules so the catalogue
    // passes strict validation.
    cfg.rules.v_anchor.required = false;
    cfg.rules.j_anchor.required = false;
    let v01 = cfg.v_pool.push(Allele {
        name: "V*01".into(),
        gene: "V".into(),
        seq: b"AAACCCGGGTTT".to_vec(),
        segment: Segment::V,
        anchor: None,
        functional_status: None,
        subregions: Vec::new(),
    });
    let _v02 = cfg.v_pool.push(Allele {
        name: "V*02".into(),
        gene: "V".into(),
        seq: b"AAACCCGGGAAA".to_vec(),
        segment: Segment::V,
        anchor: None,
        functional_status: None,
        subregions: Vec::new(),
    });
    let d01 = cfg.d_pool.push(Allele {
        name: "D*01".into(),
        gene: "D".into(),
        seq: b"TTTAAGCG".to_vec(),
        segment: Segment::D,
        anchor: None,
        functional_status: None,
        subregions: Vec::new(),
    });
    let j01 = cfg.j_pool.push(Allele {
        name: "J*01".into(),
        gene: "J".into(),
        seq: b"AAA".to_vec(),
        segment: Segment::J,
        anchor: None,
        functional_status: None,
        subregions: Vec::new(),
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
fn np1_string_retains_all_generated_bases_despite_v_recall_extension() {
    // Matching NP evidence narrows V to its two indistinguishable aliases,
    // while the standard NP1 field retains all six generated bases.
    let cfg = common::vj_ambiguous_refdata();
    let v01 = common::allele_id_by_name(&cfg, Segment::V, "v01*01");
    let j01 = common::allele_id_by_name(&cfg, Segment::J, "j01*01");
    let plan = vj_plan(&cfg, v01, j01, 3, 0, 6, b'G');
    let outcome = compile_and_run(&cfg, &plan);
    let rec = build_airr_record(&outcome, &cfg, "np1-shrink");
    assert_eq!(
        rec.np1_length, 6,
        "all six generated bases must remain in NP1",
    );
    assert_eq!(
        rec.np1, "GGGGGG",
        "NP1 string preserves the generated junction",
    );
}

#[test]
fn vj_j_recall_extension_does_not_move_generating_j_boundary() {
    // NP1=TTT narrows the J re-call after trimming its prefix. The
    // generating J still starts at reference 3 and NP1 retains all three bases.
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
        Some(3),
        "generating J begins at the sampled trim boundary"
    );
    assert_eq!(rec.np1, "TTT");
    assert_eq!(rec.np1_length, 3);
}

#[test]
fn vdj_j_recall_extension_does_not_move_generating_j_boundary() {
    // NP2=TTT identifies J from its trimmed prefix without changing
    // the generating J bounds or the three-base NP2 interval.
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
    assert_eq!(rec.j_germline_start, Some(3));
    assert_eq!(rec.np2, "TTT");
}

#[test]
fn v_cigar_excludes_np_bases_used_by_recall() {
    // The standard V CIGAR covers its twelve generating bases. The
    // three generated NP1 bases are independent of re-call extension evidence.
    let cfg = common::vj_ambiguous_refdata();
    let v01 = common::allele_id_by_name(&cfg, Segment::V, "v01*01");
    let j01 = common::allele_id_by_name(&cfg, Segment::J, "j01*01");
    let plan = vj_plan(&cfg, v01, j01, 3, 0, 3, b'G');
    let outcome = compile_and_run(&cfg, &plan);
    let rec = build_airr_record(&outcome, &cfg, "cigar-extend");
    assert_eq!(rec.v_cigar, "12M");
}
