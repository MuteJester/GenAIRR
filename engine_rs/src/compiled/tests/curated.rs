use super::indels::{DeleteAtPass, InsertAtPass};
use super::live_call_edits::EditBaseAtPass;
use super::*;

// ──────────────────────────────────────────────────────────────
// Curated end-to-end allele evaluation suite.
//
// A small, hand-crafted V/D/J refdata where every base is
// predictable. Each allele has:
//   - V (12bp, anchor at 6): shared 9bp prefix `AAACCCTGT` (with
//     the conserved Cys "TGT" at 6-8) plus a 3bp distinguishing
//     suffix (`AAA` / `CCC` / `GGG`),
//   - D (12bp): 3bp distinguishing prefix + shared `GGGCCC` core
//     + 3bp distinguishing suffix,
//   - J (9bp, anchor at 3): 3bp distinguishing prefix + shared
//     `TGGACG` (W codon TGG at 3-5).
//
// Sampling V1+D1+J1 with no trim and zero NP gives a 33bp
// assembly whose junction (pool 6..30) is in-frame and codes
// `CKKGPFKW` with no stop codons → productive.
//
// The shared cores let trim widen the live call to all three
// alleles per segment; the distinguishing edges let NP bases
// narrow the call back. Mutations / indels / corruption can
// then be applied at known positions and the resulting AIRR
// metadata is hand-checkable.
// ──────────────────────────────────────────────────────────────

/// Build the curated VDJ refdata. Allele ids are guaranteed in
/// declaration order (V1=0, V2=1, V3=2; same for D / J pools).
fn curated_v_d_j_refdata() -> RefDataConfig {
    let mut cfg = RefDataConfig::empty(ChainType::Vdj);
    // V pool — anchor=6 marks the conserved Cys (TGT) codon.
    let _ = cfg.v_pool.push(Allele {
        name: "V1*01".into(),
        gene: "V1".into(),
        seq: b"AAACCCTGTAAA".to_vec(),
        segment: Segment::V,
        anchor: Some(6),
        functional_status: None,
        subregions: Vec::new(),
    });
    let _ = cfg.v_pool.push(Allele {
        name: "V2*01".into(),
        gene: "V2".into(),
        seq: b"AAACCCTGTCCC".to_vec(),
        segment: Segment::V,
        anchor: Some(6),
        functional_status: None,
        subregions: Vec::new(),
    });
    let _ = cfg.v_pool.push(Allele {
        name: "V3*01".into(),
        gene: "V3".into(),
        seq: b"AAACCCTGTGGG".to_vec(),
        segment: Segment::V,
        anchor: Some(6),
        functional_status: None,
        subregions: Vec::new(),
    });
    // D pool — distinguishing 3bp prefix (always T-starting, so it
    // can never extend a V allele's distinguishing A/C/G suffix)
    // + 6bp shared core (`GGGCCC`) + distinguishing 3bp suffix
    // (always G-starting, so it can never extend backward into a
    // J allele's distinguishing A/C/T prefix). These guards keep
    // the boundary-overlap walker from accidentally narrowing
    // calls during simple trim-widening tests.
    let _ = cfg.d_pool.push(Allele {
        name: "D1*01".into(),
        gene: "D1".into(),
        seq: b"TTCGGGCCCGAG".to_vec(),
        segment: Segment::D,
        anchor: None,
        functional_status: None,
        subregions: Vec::new(),
    });
    let _ = cfg.d_pool.push(Allele {
        name: "D2*01".into(),
        gene: "D2".into(),
        seq: b"TATGGGCCCGCG".to_vec(),
        segment: Segment::D,
        anchor: None,
        functional_status: None,
        subregions: Vec::new(),
    });
    let _ = cfg.d_pool.push(Allele {
        name: "D3*01".into(),
        gene: "D3".into(),
        seq: b"TCGGGGCCCGTG".to_vec(),
        segment: Segment::D,
        anchor: None,
        functional_status: None,
        subregions: Vec::new(),
    });
    // J pool — anchor=3 marks the conserved W (TGG) codon.
    let _ = cfg.j_pool.push(Allele {
        name: "J1*01".into(),
        gene: "J1".into(),
        seq: b"AAATGGACG".to_vec(),
        segment: Segment::J,
        anchor: Some(3),
        functional_status: None,
        subregions: Vec::new(),
    });
    let _ = cfg.j_pool.push(Allele {
        name: "J2*01".into(),
        gene: "J2".into(),
        seq: b"CCCTGGACG".to_vec(),
        segment: Segment::J,
        anchor: Some(3),
        functional_status: None,
        subregions: Vec::new(),
    });
    let _ = cfg.j_pool.push(Allele {
        name: "J3*01".into(),
        gene: "J3".into(),
        seq: b"TTTTGGACG".to_vec(),
        segment: Segment::J,
        anchor: Some(3),
        functional_status: None,
        subregions: Vec::new(),
    });
    cfg
}

/// Build a curated VDJ plan parameterised by the events to apply.
/// All defaults are no-op (no trim, zero-length NP), so a test
/// only specifies what it needs to vary.
#[allow(clippy::too_many_arguments)]
fn curated_plan(
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

/// Run a curated plan, build the AIRR record, return both.
fn curated_run(
    cfg: &RefDataConfig,
    plan: PassPlan,
) -> (crate::pass::Outcome, crate::airr_record::AirrRecord) {
    let compiled = CompiledSimulator::compile(&plan, Some(cfg), None, ExecutionPolicy::Permissive)
        .expect("curated plan should compile");
    let outcome = compiled.run_one(0).expect("plan should run");
    let rec = build_airr_record(&outcome, cfg, "curated");
    (outcome, rec)
}

// V1 / V2 / V3 / D1 / D2 / D3 / J1 / J2 / J3 ids for legibility.
fn curated_ids() -> [AlleleId; 9] {
    [
        AlleleId::new(0), // V1
        AlleleId::new(1), // V2
        AlleleId::new(2), // V3
        AlleleId::new(0), // D1
        AlleleId::new(1), // D2
        AlleleId::new(2), // D3
        AlleleId::new(0), // J1
        AlleleId::new(1), // J2
        AlleleId::new(2), // J3
    ]
}

#[test]
fn curated_baseline_productive_no_corruption() {
    // Sample V1+D1+J1, no trim, zero NP. The assembled sequence
    // is exactly V1 + D1 + J1 = 33bp with the conserved Cys at
    // pool 6 and W at pool 27. Junction = pool[6..30] (24bp)
    // codes `CKKGPFKW` with no stops → productive.
    let cfg = curated_v_d_j_refdata();
    let [v1, _, _, d1, _, _, j1, _, _] = curated_ids();
    let plan = curated_plan(&cfg, v1, d1, j1, 0, 0, 0, 0, 0, b'A', 0, b'A');
    let (_outcome, rec) = curated_run(&cfg, plan);

    // Sequence is exact concatenation of V1 + D1 + J1.
    assert_eq!(rec.sequence, "AAACCCTGTAAATTCGGGCCCGAGAAATGGACG");
    assert_eq!(rec.sequence_length, 33);

    // Live calls narrow to the sampled allele (distinguishing
    // bases at every segment edge → exactly one allele matches).
    assert_eq!(rec.v_call, "V1*01");
    assert_eq!(rec.d_call, "D1*01");
    assert_eq!(rec.j_call, "J1*01");

    // Junction = pool[6..30] = TGT + AAA + TTCGGGCCCGAG + AAA + TGG.
    // Codons: TGT AAA TTC GGG CCC GAG AAA TGG = C K F G P E K W.
    assert_eq!(rec.junction, "TGTAAATTCGGGCCCGAGAAATGG");
    assert_eq!(rec.junction_length, Some(24));
    assert_eq!(rec.junction_aa, "CKFGPEKW");
    assert_eq!(rec.productive, Some(true));
    assert_eq!(rec.vj_in_frame, Some(true));
    assert_eq!(rec.stop_codon, Some(false));

    // Pure recombination → no mutations, indels, or errors.
    assert_eq!(rec.n_mutations, 0);
    assert_eq!(rec.n_indels, 0);
    assert_eq!(rec.n_pcr_errors, 0);
    assert!(!rec.is_contaminant);

    // CIGARs are pure M with each segment's full length.
    assert_eq!(rec.v_cigar, "12M");
    assert_eq!(rec.d_cigar, "12M");
    assert_eq!(rec.j_cigar, "9M");

    // Identity is 1.0 (every base matches the source allele).
    assert_eq!(rec.v_identity, Some(1.0));
    assert_eq!(rec.d_identity, Some(1.0));
    assert_eq!(rec.j_identity, Some(1.0));

    // Coordinate self-consistency.
    assert_eq!(rec.v_sequence_start, Some(0));
    assert_eq!(rec.v_sequence_end, Some(12));
    assert_eq!(rec.d_sequence_start, Some(12));
    assert_eq!(rec.d_sequence_end, Some(24));
    assert_eq!(rec.j_sequence_start, Some(24));
    assert_eq!(rec.j_sequence_end, Some(33));
    assert_eq!(rec.v_germline_start, Some(0));
    assert_eq!(rec.v_germline_end, Some(12));
    assert_eq!(rec.j_germline_end, Some(9));

    // Locus is derived from V's "V1*01" prefix → empty (our
    // synthetic alleles don't start with IGH/IGK/etc).
    assert_eq!(rec.locus, "");
}

#[test]
fn curated_v_trim_widens_recall_but_preserves_generating_call() {
    // Trimming removes evidence that distinguished the generating allele.
    // The re-call becomes ambiguous; standard fields retain generation truth.
    let cfg = curated_v_d_j_refdata();
    let [v1, _, _, d1, _, _, j1, _, _] = curated_ids();
    let plan = curated_plan(&cfg, v1, d1, j1, 3, 0, 0, 0, 0, b'A', 0, b'A');
    let (_outcome, rec) = curated_run(&cfg, plan);

    assert_eq!(rec.v_sequence_end, Some(9));
    assert_eq!(rec.v_germline_end, Some(9));
    assert_eq!(rec.v_trim_3, 3);

    assert_eq!(rec.v_call, "V1*01");
    assert_eq!(rec.recalled_v_call, "V1*01,V2*01,V3*01");

    assert_eq!(rec.d_call, "D1*01");
    assert_eq!(rec.j_call, "J1*01");

    assert_eq!(rec.sequence_length, 30);
    assert_eq!(rec.v_cigar, "9M");
}

#[test]
fn curated_np1_recreates_v_suffix_narrows_recall_without_changing_truth_bounds() {
    // Matching NP evidence narrows the V re-call. Its generated origin and
    // the generating reference interval stay independent of that inference.
    let cfg = curated_v_d_j_refdata();
    let [v1, _, _, d1, _, _, j1, _, _] = curated_ids();
    let plan = curated_plan(&cfg, v1, d1, j1, 3, 0, 0, 0, 3, b'A', 0, b'A');
    let (_outcome, rec) = curated_run(&cfg, plan);

    assert_eq!(rec.v_call, "V1*01");
    assert_eq!(rec.recalled_v_call, "V1*01");
    assert_eq!(rec.v_germline_end, Some(9));
    assert_eq!(rec.v_cigar, "9M");
}

#[test]
fn curated_d_trim_both_sides_widens_recall_but_preserves_generating_call() {
    // Removing distinguishing D ends widens only the evidence-based re-call.
    let cfg = curated_v_d_j_refdata();
    let [v1, _, _, d1, _, _, j1, _, _] = curated_ids();
    let plan = curated_plan(&cfg, v1, d1, j1, 0, 3, 3, 0, 0, b'A', 0, b'A');
    let (_outcome, rec) = curated_run(&cfg, plan);

    assert_eq!(rec.d_call, "D1*01");
    assert_eq!(rec.recalled_d_call, "D1*01,D2*01,D3*01");
    assert_eq!(rec.d_germline_start, Some(3));
    assert_eq!(rec.d_germline_end, Some(9));
    assert_eq!(rec.d_cigar, "6M");
}

#[test]
fn curated_unmatched_np_keeps_d_recall_ambiguous() {
    // After both D ends are trimmed, NP1=A matches none of the distinguishing
    // prefix bases. The re-call remains ambiguous and the generating D is known.
    let cfg = curated_v_d_j_refdata();
    let [v1, _, _, d1, _, _, j1, _, _] = curated_ids();
    let plan = curated_plan(&cfg, v1, d1, j1, 0, 3, 3, 0, 1, b'A', 0, b'A');
    let (_outcome, rec) = curated_run(&cfg, plan);

    assert_eq!(rec.d_call, "D1*01");
    let mut ids: Vec<&str> = rec.recalled_d_call.split(',').collect();
    ids.sort();
    assert_eq!(ids, vec!["D1*01", "D2*01", "D3*01"]);
}

#[test]
fn curated_j_trim_widens_recall_but_preserves_generating_call() {
    // Removing the distinguishing J prefix widens only the re-call.
    let cfg = curated_v_d_j_refdata();
    let [v1, _, _, d1, _, _, j1, _, _] = curated_ids();
    let plan = curated_plan(&cfg, v1, d1, j1, 0, 0, 0, 3, 0, b'A', 0, b'A');
    let (_outcome, rec) = curated_run(&cfg, plan);

    assert_eq!(rec.j_call, "J1*01");
    assert_eq!(rec.recalled_j_call, "J1*01,J2*01,J3*01");
    assert_eq!(rec.j_germline_start, Some(3));
    assert_eq!(rec.j_cigar, "6M");
}

#[test]
fn curated_np2_recreates_j_prefix_narrows_recall_without_changing_truth_bounds() {
    // NP2 matching the removed J prefix narrows its re-call. The generating
    // J still begins at reference coordinate 3 and contains six bases.
    let cfg = curated_v_d_j_refdata();
    let [v1, _, _, d1, _, _, j1, _, _] = curated_ids();
    let plan = curated_plan(&cfg, v1, d1, j1, 0, 0, 0, 3, 0, b'A', 3, b'A');
    let (_outcome, rec) = curated_run(&cfg, plan);

    assert_eq!(rec.j_call, "J1*01");
    assert_eq!(rec.recalled_j_call, "J1*01");
    assert_eq!(rec.j_germline_start, Some(3));
    assert_eq!(rec.j_cigar, "6M");
}

#[test]
fn curated_mutation_switches_recall_without_replacing_generating_reference() {
    // Mutation changes the three V1 suffix bases to match V2. The re-call
    // changes while the standard call and mismatch reference remain V1.
    let cfg = curated_v_d_j_refdata();
    let [v1, _, _, d1, _, _, j1, _, _] = curated_ids();
    let mut plan = curated_plan(&cfg, v1, d1, j1, 0, 0, 0, 0, 0, b'A', 0, b'A');
    plan.push(Box::new(EditBaseAtPass::new(9, b'C')));
    plan.push(Box::new(EditBaseAtPass::new(10, b'C')));
    plan.push(Box::new(EditBaseAtPass::new(11, b'C')));
    let (_outcome, rec) = curated_run(&cfg, plan);

    assert_eq!(rec.v_call, "V1*01");
    assert_eq!(rec.recalled_v_call, "V2*01");
    assert_eq!(&rec.germline_alignment[..12], "AAACCCTGTAAA");
    assert_eq!(rec.v_identity, Some(0.75));
}

#[test]
fn curated_indel_deletion_widens_recall_but_preserves_source_deletion() {
    // Deleting all three distinguishing suffix bases widens the re-call.
    // Their generating-reference deletion columns remain explicit.
    let cfg = curated_v_d_j_refdata();
    let [v1, _, _, d1, _, _, j1, _, _] = curated_ids();
    let mut plan = curated_plan(&cfg, v1, d1, j1, 0, 0, 0, 0, 0, b'A', 0, b'A');
    plan.push(Box::new(DeleteAtPass::new(9)));
    plan.push(Box::new(DeleteAtPass::new(9)));
    plan.push(Box::new(DeleteAtPass::new(9)));
    let (_outcome, rec) = curated_run(&cfg, plan);

    assert_eq!(rec.v_call, "V1*01");
    assert_eq!(rec.recalled_v_call, "V1*01,V2*01,V3*01");
    let v_d_count: u32 = rec
        .v_cigar
        .split_terminator(|c: char| c.is_ascii_alphabetic())
        .filter_map(|s| s.parse::<u32>().ok())
        .zip(rec.v_cigar.matches(|c: char| c.is_ascii_alphabetic()))
        .filter(|(_, op)| *op == "D")
        .map(|(n, _)| n)
        .sum();
    assert_eq!(v_d_count, 3);
}

#[test]
fn curated_indel_insertion_inside_v_does_not_break_call() {
    // Sample V1, insert a synthetic 'C' at pool position 3
    // (inside V's coding region). Walker tolerates the synthetic
    // base — call stays {V1}.
    let cfg = curated_v_d_j_refdata();
    let [v1, _, _, d1, _, _, j1, _, _] = curated_ids();
    let mut plan = curated_plan(&cfg, v1, d1, j1, 0, 0, 0, 0, 0, b'A', 0, b'A');
    plan.push(Box::new(InsertAtPass::new(3, b'C', Segment::V)));
    let (_outcome, rec) = curated_run(&cfg, plan);

    // Insertion is absorbed without collapsing the call.
    assert_eq!(rec.v_call, "V1*01");
    // Sequence grew by one base.
    assert_eq!(rec.sequence_length, 34);
    // V CIGAR contains exactly one I op (the insertion).
    assert!(
        rec.v_cigar.contains("1I"),
        "expected an I op in v_cigar, got {}",
        rec.v_cigar
    );
}

#[test]
fn v_germline_end_retains_generating_boundary_despite_recall_extension() {
    // The live V hypothesis extends into NP, but the generating V ends at 9.
    let cfg = curated_v_d_j_refdata();
    let [v1, _, _, d1, _, _, j1, _, _] = curated_ids();
    let plan = curated_plan(&cfg, v1, d1, j1, 3, 0, 0, 0, 3, b'A', 0, b'A');
    let (_outcome, rec) = curated_run(&cfg, plan);

    assert_eq!(rec.v_trim_3, 3);
    assert_eq!(rec.v_germline_end, Some(9));
}

#[test]
fn d_germline_bounds_retain_generating_trim() {
    // The primary D evidence already resolves its re-call to D1. The
    // generating D reference starts at the realized trim boundary 3.
    let cfg = curated_v_d_j_refdata();
    let [v1, _, _, d1, _, _, j1, _, _] = curated_ids();
    let plan = curated_plan(&cfg, v1, d1, j1, 0, 3, 0, 0, 1, b'C', 0, b'A');
    let (_outcome, rec) = curated_run(&cfg, plan);

    assert_eq!(rec.d_call, "D1*01");
    assert_eq!(rec.d_trim_5, 3);
    assert_eq!(rec.d_germline_start, Some(3));
}

#[test]
fn v_sequence_end_excludes_np_used_by_recall() {
    // The V re-call can use matching NP evidence without moving the
    // generating V sequence boundary from 9.
    let cfg = curated_v_d_j_refdata();
    let [v1, _, _, d1, _, _, j1, _, _] = curated_ids();
    let plan = curated_plan(&cfg, v1, d1, j1, 3, 0, 0, 0, 3, b'A', 0, b'A');
    let (_outcome, rec) = curated_run(&cfg, plan);

    assert_eq!(rec.v_sequence_end, Some(9));
}

#[test]
fn j_sequence_start_excludes_np_used_by_recall() {
    // Layout: V(12) + D(12) + NP2(3) + J6. Generating J begins at 27;
    // a re-call extension cannot move that boundary into NP.
    let cfg = curated_v_d_j_refdata();
    let [v1, _, _, d1, _, _, j1, _, _] = curated_ids();
    let plan = curated_plan(&cfg, v1, d1, j1, 0, 0, 0, 3, 0, b'A', 3, b'A');
    let (_outcome, rec) = curated_run(&cfg, plan);

    assert_eq!(rec.j_sequence_start, Some(27));
    assert_eq!(rec.j_sequence_end, Some(33));
}

#[test]
fn no_extension_preserves_structural_germline_bounds() {
    // Sanity: when NO extension fires (no NP, or NP doesn't
    // match), AIRR germline bounds match the structural trim
    // range. (Live-call hypothesis ref_start/ref_end match the
    // structural values in this case.)
    let cfg = curated_v_d_j_refdata();
    let [v1, _, _, d1, _, _, j1, _, _] = curated_ids();
    // Trim V_3 = 3, no NP1 → V's right-extension halts
    // immediately at the V/D boundary (D1's first base 'T'
    // doesn't match any V allele's pos 9). v_germline_end stays
    // at 9 = structural post-trim end.
    let plan = curated_plan(&cfg, v1, d1, j1, 3, 0, 0, 0, 0, b'A', 0, b'A');
    let (_outcome, rec) = curated_run(&cfg, plan);

    assert_eq!(rec.v_trim_3, 3);
    assert_eq!(rec.v_germline_end, Some(9));
}

#[test]
fn v_cigar_excludes_np1_used_by_recall() {
    // The generating V contains nine bases even though NP1 narrows its re-call.
    let cfg = curated_v_d_j_refdata();
    let [v1, _, _, d1, _, _, j1, _, _] = curated_ids();
    let plan = curated_plan(&cfg, v1, d1, j1, 3, 0, 0, 0, 3, b'A', 0, b'A');
    let (_outcome, rec) = curated_run(&cfg, plan);

    assert_eq!(rec.v_cigar, "9M");
}

#[test]
fn j_cigar_excludes_np2_used_by_recall() {
    // The generating J contains six bases; matching NP2 bases remain junction bases.
    let cfg = curated_v_d_j_refdata();
    let [v1, _, _, d1, _, _, j1, _, _] = curated_ids();
    let plan = curated_plan(&cfg, v1, d1, j1, 0, 0, 0, 3, 0, b'A', 3, b'A');
    let (_outcome, rec) = curated_run(&cfg, plan);

    assert_eq!(rec.j_cigar, "6M");
}

#[test]
fn np1_string_retains_generated_bases_used_by_v_recall() {
    // All three generated NP1 bases remain NP, including the base that
    // resolves the V re-call ambiguity.
    let cfg = curated_v_d_j_refdata();
    let [v1, _, _, d1, _, _, j1, _, _] = curated_ids();
    let plan = curated_plan(&cfg, v1, d1, j1, 3, 0, 0, 0, 3, b'A', 0, b'A');
    let (_outcome, rec) = curated_run(&cfg, plan);

    assert_eq!(rec.np1, "AAA");
    assert_eq!(rec.np1_length, 3);
}

#[test]
fn np2_string_retains_generated_bases_used_by_j_recall() {
    // All three generated NP2 bases remain NP, including the base that
    // resolves the J re-call ambiguity.
    let cfg = curated_v_d_j_refdata();
    let [v1, _, _, d1, _, _, j1, _, _] = curated_ids();
    let plan = curated_plan(&cfg, v1, d1, j1, 0, 0, 0, 3, 0, b'A', 3, b'A');
    let (_outcome, rec) = curated_run(&cfg, plan);

    assert_eq!(rec.np2, "AAA");
    assert_eq!(rec.np2_length, 3);
}

#[test]
fn v_alignment_end_excludes_np_recall_evidence() {
    // The nine generating V columns determine its alignment span.
    let cfg = curated_v_d_j_refdata();
    let [v1, _, _, d1, _, _, j1, _, _] = curated_ids();
    let plan = curated_plan(&cfg, v1, d1, j1, 3, 0, 0, 0, 3, b'A', 0, b'A');
    let (_outcome, rec) = curated_run(&cfg, plan);

    assert_eq!(rec.v_alignment_start, Some(0));
    assert_eq!(rec.v_alignment_end, Some(9));
}

#[test]
fn j_alignment_start_excludes_np_recall_evidence() {
    // Generating J starts after V(12) + D(12) + NP2(3), at column 27.
    let cfg = curated_v_d_j_refdata();
    let [v1, _, _, d1, _, _, j1, _, _] = curated_ids();
    let plan = curated_plan(&cfg, v1, d1, j1, 0, 0, 0, 3, 0, b'A', 3, b'A');
    let (_outcome, rec) = curated_run(&cfg, plan);

    assert_eq!(rec.j_alignment_start, Some(27));
    assert_eq!(rec.j_alignment_end, Some(33));
}

#[test]
fn v_identity_counts_generating_segment_columns() {
    // Identity compares the nine generating V bases to their source reference.
    let cfg = curated_v_d_j_refdata();
    let [v1, _, _, d1, _, _, j1, _, _] = curated_ids();
    let plan = curated_plan(&cfg, v1, d1, j1, 3, 0, 0, 0, 3, b'A', 0, b'A');
    let (_outcome, rec) = curated_run(&cfg, plan);

    assert_eq!(rec.v_identity, Some(1.0));
}

#[test]
fn junction_locates_anchors_via_germline_pos() {
    // Baseline: the curated default plan places V's Cys anchor
    // at pool position 6 (V allele anchor=6, V at pool[0..12])
    // and J's W anchor at pool position 27
    // (J allele anchor=3, J at pool[24..33]). Junction =
    // pool[6..30], length 24 = the 3bp Cys + (V tail 3) +
    // (D 12) + (J head 3) + 3bp W.
    let cfg = curated_v_d_j_refdata();
    let [v1, _, _, d1, _, _, j1, _, _] = curated_ids();
    let plan = curated_plan(&cfg, v1, d1, j1, 0, 0, 0, 0, 0, b'A', 0, b'A');
    let (_outcome, rec) = curated_run(&cfg, plan);

    assert_eq!(rec.junction_start, Some(6));
    assert_eq!(rec.junction_end, Some(30));
    assert_eq!(rec.junction_length, Some(24));
    // The anchor codons frame the junction.
    assert!(
        rec.junction.starts_with("TGT"),
        "junction should start with V Cys codon, got {:?}",
        rec.junction,
    );
    assert!(
        rec.junction.ends_with("TGG"),
        "junction should end with J W codon, got {:?}",
        rec.junction,
    );
}

#[test]
fn junction_shifts_with_v_insertion_before_anchor() {
    // Insert a synthetic 'C' at pool position 3 (inside V, BEFORE
    // the anchor codon at pool 6). The junction must follow the
    // anchor's actual pool position — anchor germline position 6
    // now resides at pool index 7, so junction_start shifts from
    // 6 → 7.
    let cfg = curated_v_d_j_refdata();
    let [v1, _, _, d1, _, _, j1, _, _] = curated_ids();
    let mut plan = curated_plan(&cfg, v1, d1, j1, 0, 0, 0, 0, 0, b'A', 0, b'A');
    plan.push(Box::new(InsertAtPass::new(3, b'C', Segment::V)));
    let (_outcome, rec) = curated_run(&cfg, plan);

    assert_eq!(rec.junction_start, Some(7));
    assert!(
        rec.junction.starts_with("TGT"),
        "junction should still frame the V Cys codon after V insertion, got {:?}",
        rec.junction,
    );
}

#[test]
fn junction_shifts_with_v_deletion_before_anchor() {
    // Symmetric: delete at pool position 3 (inside V, before the
    // anchor codon). Anchor germline_pos=6 now sits at pool 5
    // (one earlier). junction_start should follow.
    let cfg = curated_v_d_j_refdata();
    let [v1, _, _, d1, _, _, j1, _, _] = curated_ids();
    let mut plan = curated_plan(&cfg, v1, d1, j1, 0, 0, 0, 0, 0, b'A', 0, b'A');
    plan.push(Box::new(DeleteAtPass::new(3)));
    let (_outcome, rec) = curated_run(&cfg, plan);

    assert_eq!(rec.junction_start, Some(5));
    assert!(
        rec.junction.starts_with("TGT"),
        "junction should still frame the V Cys codon after V deletion, got {:?}",
        rec.junction,
    );
}

#[test]
fn junction_shifts_with_j_insertion_before_anchor() {
    // J side mirror: inserting a base inside J (before its W
    // anchor at allele pos 3) pushes J's anchor pool position
    // rightward, so junction_end grows accordingly.
    //
    // Pool layout: V(12) + NP1(0) + D(12) + NP2(0) + J(9) = 33.
    //   J at pool[24..33]; anchor germline_pos=3 sits at pool 27.
    //   junction_end = 27 + 3 = 30 baseline.
    // After inserting at pool 25 (within J, before anchor), the
    // anchor germline_pos=3 node shifts to pool 28, so
    // junction_end = 28 + 3 = 31.
    let cfg = curated_v_d_j_refdata();
    let [v1, _, _, d1, _, _, j1, _, _] = curated_ids();
    let mut plan = curated_plan(&cfg, v1, d1, j1, 0, 0, 0, 0, 0, b'A', 0, b'A');
    plan.push(Box::new(InsertAtPass::new(25, b'C', Segment::J)));
    let (_outcome, rec) = curated_run(&cfg, plan);

    assert_eq!(rec.junction_end, Some(31));
    assert!(
        rec.junction.ends_with("TGG"),
        "junction should still frame the J W codon after J insertion, got {:?}",
        rec.junction,
    );
}

#[test]
fn upstream_indel_preserves_anchor_coordinates_but_disrupts_coding_frame() {
    // Insert a synthetic 'C' at pool position 3 (inside V, before
    // the anchor). Both V and J anchors shift right by 1, so
    // junction_start = 7, junction_end = 31, junction_length = 24
    // — same junction length as baseline. The observed V coding origin
    // stays fixed, however, so the displaced anchor is out of frame.
    let cfg = curated_v_d_j_refdata();
    let [v1, _, _, d1, _, _, j1, _, _] = curated_ids();
    let mut plan = curated_plan(&cfg, v1, d1, j1, 0, 0, 0, 0, 0, b'A', 0, b'A');
    plan.push(Box::new(InsertAtPass::new(3, b'C', Segment::V)));
    let (_outcome, rec) = curated_run(&cfg, plan);

    assert_eq!(rec.junction_start, Some(7));
    assert_eq!(rec.junction_end, Some(31));
    assert_eq!(rec.junction_length, Some(24));
    assert_eq!(rec.vj_in_frame, Some(false));
    assert_eq!(rec.productive, Some(false));
}

#[test]
fn germline_span_equals_m_plus_d() {
    // Invariant: AIRR `*_germline_start/_end` come from the
    // column walker's `ref_ranges` (the union of ref positions
    // consumed by `M` and `D` ops), so the strict identity
    // `germline_span == M + D` holds by construction.
    //
    // The curated default plan is no-op trims, no NP, no
    // corruption — V/D/J each contribute a clean structural
    // CIGAR (12M / 12M / 9M) over their full allele span, and
    // germline coords reflect that 1:1.
    let cfg = curated_v_d_j_refdata();
    let [v1, _, _, d1, _, _, j1, _, _] = curated_ids();
    let plan = curated_plan(&cfg, v1, d1, j1, 0, 0, 0, 0, 0, b'A', 0, b'A');
    let (_outcome, rec) = curated_run(&cfg, plan);

    for (cig, g_s, g_e) in [
        (
            rec.v_cigar.as_str(),
            rec.v_germline_start,
            rec.v_germline_end,
        ),
        (
            rec.d_cigar.as_str(),
            rec.d_germline_start,
            rec.d_germline_end,
        ),
        (
            rec.j_cigar.as_str(),
            rec.j_germline_start,
            rec.j_germline_end,
        ),
    ] {
        let m_count: u32 = cig
            .split_terminator(|c: char| c.is_ascii_alphabetic())
            .filter_map(|s| s.parse::<u32>().ok())
            .zip(cig.matches(|c: char| c.is_ascii_alphabetic()))
            .filter(|(_, op)| *op == "M" || *op == "D")
            .map(|(n, _)| n)
            .sum();
        let span = g_e.unwrap() - g_s.unwrap();
        assert_eq!(
            span, m_count as i64,
            "germline_span ({span}) != M+D ({m_count}) for cigar {cig:?}",
        );
    }
}

#[test]
fn germline_span_under_v_indel_deletion() {
    // V_3 trim 3 + NP1="AAA" extends V to ref 12 via NP1 claim,
    // then a structural deletion at pool 1 removes one V base
    // (germline_pos=1). CIGAR walks the structural region and
    // emits a D op for the missing ref=1 — `ref_ranges` covers
    // [0, 12) (NP1 claim at ref=9..12, structural at ref=0..9
    // including the D-fill at ref=1). M+D == 12 == germline_span.
    let cfg = curated_v_d_j_refdata();
    let [v1, _, _, d1, _, _, j1, _, _] = curated_ids();
    let mut plan = curated_plan(&cfg, v1, d1, j1, 3, 0, 0, 0, 3, b'A', 0, b'A');
    plan.push(Box::new(DeleteAtPass::new(1)));
    let (_outcome, rec) = curated_run(&cfg, plan);

    let v_span = rec.v_germline_end.unwrap() - rec.v_germline_start.unwrap();
    let cig = rec.v_cigar.as_str();
    let m_count: u32 = cig
        .split_terminator(|c: char| c.is_ascii_alphabetic())
        .filter_map(|s| s.parse::<u32>().ok())
        .zip(cig.matches(|c: char| c.is_ascii_alphabetic()))
        .filter(|(_, op)| *op == "M" || *op == "D")
        .map(|(n, _)| n)
        .sum();
    assert_eq!(
        v_span, m_count as i64,
        "germline_span {v_span} != M+D {m_count} for v_cigar {cig:?}",
    );
}

#[test]
fn curated_full_corruption_stack_keeps_metadata_self_consistent() {
    // Stack a mutation, an insertion, and a deletion across V's
    // region. All metadata (sequence, alignments, CIGAR, identity)
    // must remain mutually consistent.
    let cfg = curated_v_d_j_refdata();
    let [v1, _, _, d1, _, _, j1, _, _] = curated_ids();
    let mut plan = curated_plan(&cfg, v1, d1, j1, 0, 0, 0, 0, 0, b'A', 0, b'A');
    // Mutate position 3 (inside V's shared prefix).
    plan.push(Box::new(EditBaseAtPass::new(3, b'T')));
    // Insert a synthetic 'G' at position 6 (inside V).
    plan.push(Box::new(InsertAtPass::new(6, b'G', Segment::V)));
    // Delete position 1 (also inside V).
    plan.push(Box::new(DeleteAtPass::new(1)));
    let (_outcome, rec) = curated_run(&cfg, plan);

    // Mutual consistency: alignment columns = (sequence chars) +
    // (gap columns where sequence is "-").
    let n_seq_gaps = rec.sequence_alignment.matches('-').count();
    assert_eq!(
        rec.sequence_alignment.len(),
        rec.sequence.len() + n_seq_gaps,
    );

    // Removing gaps from sequence_alignment yields the
    // upper-cased sequence.
    let sa_no_gaps: String = rec
        .sequence_alignment
        .chars()
        .filter(|c| *c != '-')
        .collect();
    assert_eq!(sa_no_gaps, rec.sequence.to_uppercase());

    // sequence_alignment / germline_alignment / d_mask all share
    // the same length.
    assert_eq!(rec.sequence_alignment.len(), rec.germline_alignment.len());
    assert_eq!(
        rec.sequence_alignment.len(),
        rec.germline_alignment_d_mask.len()
    );

    // (We don't assert `rec.n_indels` because the test-only
    // `InsertAtPass` / `DeleteAtPass` don't write to
    // `corrupt.indel.count`. The structural-consistency checks
    // above are what this fixture is validating.)
}
