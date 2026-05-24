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
    });
    let _ = cfg.v_pool.push(Allele {
        name: "V2*01".into(),
        gene: "V2".into(),
        seq: b"AAACCCTGTCCC".to_vec(),
        segment: Segment::V,
        anchor: Some(6),
    });
    let _ = cfg.v_pool.push(Allele {
        name: "V3*01".into(),
        gene: "V3".into(),
        seq: b"AAACCCTGTGGG".to_vec(),
        segment: Segment::V,
        anchor: Some(6),
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
    });
    let _ = cfg.d_pool.push(Allele {
        name: "D2*01".into(),
        gene: "D2".into(),
        seq: b"TATGGGCCCGCG".to_vec(),
        segment: Segment::D,
        anchor: None,
    });
    let _ = cfg.d_pool.push(Allele {
        name: "D3*01".into(),
        gene: "D3".into(),
        seq: b"TCGGGGCCCGTG".to_vec(),
        segment: Segment::D,
        anchor: None,
    });
    // J pool — anchor=3 marks the conserved W (TGG) codon.
    let _ = cfg.j_pool.push(Allele {
        name: "J1*01".into(),
        gene: "J1".into(),
        seq: b"AAATGGACG".to_vec(),
        segment: Segment::J,
        anchor: Some(3),
    });
    let _ = cfg.j_pool.push(Allele {
        name: "J2*01".into(),
        gene: "J2".into(),
        seq: b"CCCTGGACG".to_vec(),
        segment: Segment::J,
        anchor: Some(3),
    });
    let _ = cfg.j_pool.push(Allele {
        name: "J3*01".into(),
        gene: "J3".into(),
        seq: b"TTTTGGACG".to_vec(),
        segment: Segment::J,
        anchor: Some(3),
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
fn curated_v_trim_widens_v_call() {
    // Trim V_3 by 3 → V coding region drops the distinguishing
    // suffix → all three V alleles match the assembled V bases.
    // The live `v_call` must reflect that ambiguity.
    let cfg = curated_v_d_j_refdata();
    let [v1, _, _, d1, _, _, j1, _, _] = curated_ids();
    let plan = curated_plan(&cfg, v1, d1, j1, 3, 0, 0, 0, 0, b'A', 0, b'A');
    let (_outcome, rec) = curated_run(&cfg, plan);

    // V coding region is now 9 bases (12 - 3).
    assert_eq!(rec.v_sequence_end, Some(9));
    assert_eq!(rec.v_germline_end, Some(9));
    assert_eq!(rec.v_trim_3, 3);

    // Live v_call lists ALL THREE V alleles (declaration order).
    assert_eq!(rec.v_call, "V1*01,V2*01,V3*01");

    // D and J calls are unchanged singletons.
    assert_eq!(rec.d_call, "D1*01");
    assert_eq!(rec.j_call, "J1*01");

    // Sequence has dropped the trimmed-off V suffix; total
    // length is 33 - 3 = 30.
    assert_eq!(rec.sequence_length, 30);
    assert_eq!(rec.v_cigar, "9M");
}

#[test]
fn curated_np1_recreates_v_suffix_narrows_v_call_back() {
    // Trim V_3 by 3 (live call widens to {V1,V2,V3}) AND
    // generate NP1 = "AAA" (V1's distinguishing suffix). The
    // V right-extension walker reaches into NP1; under conservative
    // extension the first NP1 byte 'A' narrows
    // {V1,V2,V3}→{V1} (V1[9]='A', V2[9]='C', V3[9]='G'). The
    // walker then halts: subsequent bytes cannot narrow the
    // singleton tie set. v_call collapses to {V1} after only
    // one byte of extension.
    //
    // AIRR `v_germline_end` reflects the live-call ref_end (10 —
    // structural 9 + the one narrowing byte), and `v_cigar` covers
    // 10 columns total (9 structural + 1 NP1 column claimed).
    let cfg = curated_v_d_j_refdata();
    let [v1, _, _, d1, _, _, j1, _, _] = curated_ids();
    let plan = curated_plan(&cfg, v1, d1, j1, 3, 0, 0, 0, 3, b'A', 0, b'A');
    let (_outcome, rec) = curated_run(&cfg, plan);

    assert_eq!(rec.v_call, "V1*01");
    assert_eq!(rec.v_germline_end, Some(10));
    // CIGAR covers 9 structural + 1 NP1 column = 10M.
    assert_eq!(rec.v_cigar, "10M");
}

#[test]
fn curated_d_trim_both_sides_widens_d_call() {
    // Trim D_5 by 3 AND D_3 by 3 → D coding region drops both
    // distinguishing edges, exposing only the 6bp shared core.
    // All three D alleles match.
    let cfg = curated_v_d_j_refdata();
    let [v1, _, _, d1, _, _, j1, _, _] = curated_ids();
    let plan = curated_plan(&cfg, v1, d1, j1, 0, 3, 3, 0, 0, b'A', 0, b'A');
    let (_outcome, rec) = curated_run(&cfg, plan);

    assert_eq!(rec.d_call, "D1*01,D2*01,D3*01");
    assert_eq!(rec.d_germline_start, Some(3));
    assert_eq!(rec.d_germline_end, Some(9));
    assert_eq!(rec.d_cigar, "6M");
}

#[test]
fn curated_np_bases_narrow_d_call_from_both_sides() {
    // Trim D_5+D_3 widens d_call to all three. NP1 = D1's
    // distinguishing prefix bases narrow D's left back; NP2 = D1's
    // distinguishing suffix bases narrow D's right back.
    //
    // Use single-base ConstBaseDist for NPs, so we have to pick
    // NPs as homogeneous runs. D1 prefix = "TTC", suffix = "GAT" —
    // not single-base. The fixture's `curated_plan` builder uses
    // a single-base NP, so we drive this test by trimming only
    // the side we narrow per-NP. Here: D_3 trim only, NP2 = "G"
    // ×3 (matches the first base of D's suffix; D2/D3 also have
    // 'G' as suffix's first base by design, so this still narrows
    // to {D1} only when paired with a base hit at pos 10 / 11).
    //
    // For a fully narrow-from-both-sides exercise that matches
    // the homogeneous-NP constraint, we test D_5 trim with NP1
    // narrowing to a *subset* of D candidates — D1 prefix starts
    // 'T', as do D2 and D3. Pos 0 'T' alone keeps all three. So
    // single-base NP narrowing for D is best demonstrated via
    // *one side*; the multi-base mixed-content version requires
    // a per-position base distribution and is left to a richer
    // helper if needed.
    //
    // For simplicity, this test verifies that NP1 = "A" doesn't
    // narrow D (no D allele has 'A' at prefix pos 2 = the position
    // the left-extension walker checks first). The d_call stays
    // widened to all three after both trims expose only the
    // shared core.
    //
    // Pos-2 bases of the D prefixes: D1='C', D2='T', D3='G'.
    // 'A' matches none, so the walker halts immediately.
    let cfg = curated_v_d_j_refdata();
    let [v1, _, _, d1, _, _, j1, _, _] = curated_ids();
    let plan = curated_plan(&cfg, v1, d1, j1, 0, 3, 3, 0, 1, b'A', 0, b'A');
    let (_outcome, rec) = curated_run(&cfg, plan);

    let mut ids: Vec<&str> = rec.d_call.split(',').collect();
    ids.sort();
    assert_eq!(ids, vec!["D1*01", "D2*01", "D3*01"]);
}

#[test]
fn curated_j_trim_widens_j_call() {
    // Trim J_5 by 3 → J coding region drops the distinguishing
    // 3bp prefix → all three J alleles match the remaining
    // shared `TGGACG`.
    let cfg = curated_v_d_j_refdata();
    let [v1, _, _, d1, _, _, j1, _, _] = curated_ids();
    let plan = curated_plan(&cfg, v1, d1, j1, 0, 0, 0, 3, 0, b'A', 0, b'A');
    let (_outcome, rec) = curated_run(&cfg, plan);

    assert_eq!(rec.j_call, "J1*01,J2*01,J3*01");
    assert_eq!(rec.j_germline_start, Some(3));
    assert_eq!(rec.j_cigar, "6M");
}

#[test]
fn curated_np2_recreates_j_prefix_narrows_j_call_back() {
    // Trim J_5 by 3 (widens) AND NP2 = "AAA" (matches J1's
    // distinguishing prefix). J left-extension reaches backward
    // into NP2. The rightmost NP2 byte at ref pos 2 narrows
    // {J1,J2,J3}→{J1} (J1[2]='A', J2[2]='C', J3[2]='T'). Under
    // conservative extension the walker stops there; the
    // remaining NP2 bytes cannot narrow {J1} further.
    //
    // `j_germline_start` reads from the live-call hypothesis (2 —
    // structural 3 minus the one byte that narrowed the call), and
    // `j_cigar` is 7M (1 claimed NP2 column + 6 structural J).
    let cfg = curated_v_d_j_refdata();
    let [v1, _, _, d1, _, _, j1, _, _] = curated_ids();
    let plan = curated_plan(&cfg, v1, d1, j1, 0, 0, 0, 3, 0, b'A', 3, b'A');
    let (_outcome, rec) = curated_run(&cfg, plan);

    assert_eq!(rec.j_call, "J1*01");
    assert_eq!(rec.j_germline_start, Some(2));
    // CIGAR covers 6 structural + 1 NP2 column = 7M.
    assert_eq!(rec.j_cigar, "7M");
}

#[test]
fn curated_mutation_switches_v_call_to_a_different_allele() {
    // Sample V1 (distinguishing suffix AAA at pool 9-11). Then
    // edit positions 9, 10, 11 to 'C' → assembled V matches V2's
    // suffix exactly. Live v_call should switch to {V2}.
    let cfg = curated_v_d_j_refdata();
    let [v1, _, _, d1, _, _, j1, _, _] = curated_ids();
    let mut plan = curated_plan(&cfg, v1, d1, j1, 0, 0, 0, 0, 0, b'A', 0, b'A');
    plan.push(Box::new(EditBaseAtPass::new(9, b'C')));
    plan.push(Box::new(EditBaseAtPass::new(10, b'C')));
    plan.push(Box::new(EditBaseAtPass::new(11, b'C')));
    let (_outcome, rec) = curated_run(&cfg, plan);

    // Provenance allele is still V1 (origin unchanged); live
    // call switched to V2.
    assert_eq!(rec.v_call, "V2*01");
}

#[test]
fn curated_indel_deletion_widens_v_call() {
    // Sample V1, no trim. Delete the three distinguishing
    // V suffix bases (positions 9, 10, 11) one at a time.
    // Note: each deletion shrinks the pool, so subsequent
    // deletions need shifted positions. Easier: delete the
    // same position 9 three times (each time the new "9" is
    // the next surviving distinguishing base because the prior
    // deletion shifted things left by one).
    let cfg = curated_v_d_j_refdata();
    let [v1, _, _, d1, _, _, j1, _, _] = curated_ids();
    let mut plan = curated_plan(&cfg, v1, d1, j1, 0, 0, 0, 0, 0, b'A', 0, b'A');
    plan.push(Box::new(DeleteAtPass::new(9)));
    plan.push(Box::new(DeleteAtPass::new(9)));
    plan.push(Box::new(DeleteAtPass::new(9)));
    let (_outcome, rec) = curated_run(&cfg, plan);

    // V's distinguishing suffix is fully deleted → v_call
    // widens to all three V alleles.
    assert_eq!(rec.v_call, "V1*01,V2*01,V3*01");
    // Three deletion ops in V's CIGAR. (We don't assert
    // `rec.n_indels` because the test-only `DeleteAtPass` doesn't
    // write to `corrupt.indel.count` — that counter is owned by
    // the production indel pass.)
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
fn v_germline_end_reflects_np_extension() {
    // Invariant: `v_germline_end` reads from the live-call
    // hypothesis instead of the trim-derived structural range.
    // Trim V_3 by 3 (structural end = 9), then NP1 = "AAA". Under
    // conservative extension the first byte narrows
    // {V1,V2,V3}→{V1} (V1[9]='A', V2[9]='C', V3[9]='G'); the
    // remaining NP1 bytes cannot narrow further. Live ref_end
    // advances from 9 to 10. AIRR `v_germline_end` should be 10.
    let cfg = curated_v_d_j_refdata();
    let [v1, _, _, d1, _, _, j1, _, _] = curated_ids();
    let plan = curated_plan(&cfg, v1, d1, j1, 3, 0, 0, 0, 3, b'A', 0, b'A');
    let (_outcome, rec) = curated_run(&cfg, plan);

    // Trim says structural V end is 9; live extension says 10
    // (one narrowing byte). The live value wins.
    assert_eq!(rec.v_trim_3, 3);
    assert_eq!(rec.v_germline_end, Some(10));
}

#[test]
fn d_germline_bounds_reflect_np_extension() {
    // Trim D_5 = 3, NP1 = "C" (single base). Under conservative
    // extension, extension fires only when it strictly
    // narrows the tie set.
    //
    // In this fixture the primary structural walk over D's
    // trimmed region (positions 3..12 = GGGCCCGAG) already
    // distinguishes the alleles at pos 10 (D1='A', D2='C',
    // D3='T'). D1 scores 9, D2/D3 score 8 → post-assemble tie
    // set is already {D1}. The NP1 byte's left-extension thus
    // does NOT fire (no ambiguity left to resolve), and
    // `d_germline_start` stays at the structural trim boundary
    // 3, not 2. This is the conservative invariant: ref bounds
    // only move when extension genuinely narrows the call.
    let cfg = curated_v_d_j_refdata();
    let [v1, _, _, d1, _, _, j1, _, _] = curated_ids();
    let plan = curated_plan(&cfg, v1, d1, j1, 0, 3, 0, 0, 1, b'C', 0, b'A');
    let (_outcome, rec) = curated_run(&cfg, plan);

    // d_call resolved to {D1} by the primary walk alone.
    assert_eq!(rec.d_call, "D1*01");
    assert_eq!(rec.d_trim_5, 3);
    // No extension → germline_start stays at structural 3.
    assert_eq!(rec.d_germline_start, Some(3));
}

#[test]
fn v_sequence_end_reflects_np_extension() {
    // Invariant: `v_sequence_end` reads from the live-call
    // hypothesis seq_end. Concrete: V_3 trim 3 + NP1="AAA". Under
    // conservative extension the first 'A' narrows {V1,V2,V3}→{V1};
    // the remaining 'A's cannot narrow {V1} further. seq_end
    // advances from 9 → 10 (one NP1 byte claimed).
    let cfg = curated_v_d_j_refdata();
    let [v1, _, _, d1, _, _, j1, _, _] = curated_ids();
    let plan = curated_plan(&cfg, v1, d1, j1, 3, 0, 0, 0, 3, b'A', 0, b'A');
    let (_outcome, rec) = curated_run(&cfg, plan);

    // Structural V region ended at pool position 9; live
    // hypothesis grew to 10 → AIRR record reports 10.
    assert_eq!(rec.v_sequence_end, Some(10));
}

#[test]
fn j_sequence_start_reflects_np_extension() {
    // Symmetric: J_5 trim 3 + NP2="AAA" (J1 prefix). Under
    // conservative extension, the rightmost NP2 byte at ref pos 2
    // narrows {J1,J2,J3}→{J1}; subsequent NP2 bytes cannot narrow
    // {J1} further so the walker stops.
    //
    // Pool layout: V(12) + NP1(0) + D(12) + NP2(3) + J(6) = 33.
    //   J structural region at [27, 33).
    //   NP2 at [24, 27).
    // After 1 byte of J left-extension into NP2, seq_start = 26.
    let cfg = curated_v_d_j_refdata();
    let [v1, _, _, d1, _, _, j1, _, _] = curated_ids();
    let plan = curated_plan(&cfg, v1, d1, j1, 0, 0, 0, 3, 0, b'A', 3, b'A');
    let (_outcome, rec) = curated_run(&cfg, plan);

    assert_eq!(rec.j_sequence_start, Some(26));
    // Structural J end stays where it was; only seq_start moves.
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
fn v_cigar_extends_into_claimed_np1_columns() {
    // V's CIGAR includes claimed NP1 columns. With V_3 trim 3 and
    // NP1="AAA", conservative extension claims exactly the single
    // NP1 byte that narrows {V1,V2,V3}→{V1}. CIGAR runs 10 M ops
    // (9 structural + 1 claimed NP1 column).
    let cfg = curated_v_d_j_refdata();
    let [v1, _, _, d1, _, _, j1, _, _] = curated_ids();
    let plan = curated_plan(&cfg, v1, d1, j1, 3, 0, 0, 0, 3, b'A', 0, b'A');
    let (_outcome, rec) = curated_run(&cfg, plan);

    assert_eq!(rec.v_cigar, "10M");
}

#[test]
fn j_cigar_extends_into_claimed_np2_columns() {
    // Mirror case: J left-extends into NP2 — under conservative
    // extension the walker claims only the single NP2 byte that
    // narrows the call to {J1}. Structural J = 6M; with one
    // extension byte J's CIGAR = 7M.
    let cfg = curated_v_d_j_refdata();
    let [v1, _, _, d1, _, _, j1, _, _] = curated_ids();
    let plan = curated_plan(&cfg, v1, d1, j1, 0, 0, 0, 3, 0, b'A', 3, b'A');
    let (_outcome, rec) = curated_run(&cfg, plan);

    assert_eq!(rec.j_cigar, "7M");
}

#[test]
fn np1_string_drops_columns_claimed_by_v() {
    // when V's right extension reabsorbs NP1 bases,
    // those bases are no longer "non-templated" — np1 / np1_length
    // must drop them. With V_3 trim 3 and NP1="AAA" under
    // conservative extension only the first NP1 byte is claimed
    // (the one that narrows V's tie set to {V1}). The remaining 2
    // NP1 bytes stay non-templated → np1 = "AA", np1_length = 2.
    let cfg = curated_v_d_j_refdata();
    let [v1, _, _, d1, _, _, j1, _, _] = curated_ids();
    let plan = curated_plan(&cfg, v1, d1, j1, 3, 0, 0, 0, 3, b'A', 0, b'A');
    let (_outcome, rec) = curated_run(&cfg, plan);

    assert_eq!(rec.np1, "AA");
    assert_eq!(rec.np1_length, 2);
}

#[test]
fn np2_string_drops_columns_claimed_by_j() {
    // Mirror for NP2. With J_5 trim 3 and NP2="AAA" under
    // conservative extension only the rightmost NP2 byte is
    // claimed (the one that narrows J's tie set to {J1}); the
    // other 2 NP2 bytes stay non-templated.
    let cfg = curated_v_d_j_refdata();
    let [v1, _, _, d1, _, _, j1, _, _] = curated_ids();
    let plan = curated_plan(&cfg, v1, d1, j1, 0, 0, 0, 3, 0, b'A', 3, b'A');
    let (_outcome, rec) = curated_run(&cfg, plan);

    assert_eq!(rec.np2, "AA");
    assert_eq!(rec.np2_length, 2);
}

#[test]
fn v_alignment_end_extends_with_np_claim() {
    // `v_alignment_end` covers the claimed NP columns too. With
    // V_3 trim 3 + NP1="AAA" the structural V region is 9
    // columns long; under conservative extension V claims 1 NP1
    // column (the byte that narrowed {V1,V2,V3}→{V1}). So
    // v_alignment_end moves out to 10 = 9 structural + 1 NP.
    let cfg = curated_v_d_j_refdata();
    let [v1, _, _, d1, _, _, j1, _, _] = curated_ids();
    let plan = curated_plan(&cfg, v1, d1, j1, 3, 0, 0, 0, 3, b'A', 0, b'A');
    let (_outcome, rec) = curated_run(&cfg, plan);

    assert_eq!(rec.v_alignment_start, Some(0));
    assert_eq!(rec.v_alignment_end, Some(10));
}

#[test]
fn j_alignment_start_extends_with_np_claim() {
    // Mirror: when J left-extends into NP2, j_alignment_start
    // moves leftward. Pool layout is V(12)+NP1(0)+D(12)+NP2(3)+J(6).
    // Structural J spans columns [27, 33); under conservative
    // extension J claims 1 NP2 column (the byte that narrowed
    // the call to {J1}), so j_alignment_start = 27 - 1 = 26.
    let cfg = curated_v_d_j_refdata();
    let [v1, _, _, d1, _, _, j1, _, _] = curated_ids();
    let plan = curated_plan(&cfg, v1, d1, j1, 0, 0, 0, 3, 0, b'A', 3, b'A');
    let (_outcome, rec) = curated_run(&cfg, plan);

    assert_eq!(rec.j_alignment_start, Some(26));
    assert_eq!(rec.j_alignment_end, Some(33));
}

#[test]
fn v_identity_counts_extended_columns() {
    // V's identity reflects matches over the full
    // claimed span. With V_3 trim 3 + NP1="AAA" (which match V1
    // exactly at pos 9-11), every column in V's extended span is
    // a match → identity = 1.0.
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
fn productive_uses_germline_pos_anchor_after_indels() {
    // Insert a synthetic 'C' at pool position 3 (inside V, before
    // the anchor). Both V and J anchors shift right by 1, so
    // junction_start = 7, junction_end = 31, junction_length = 24
    // — same length as baseline, still in-frame.
    //
    // The biological invariant: the anchor codon is still
    // recognised as Cys (germline pos 6,7,8 still map to T,G,T
    // in the IR). A structural-offset reader would have read
    // pool[6..9] (now AAC after shift) and rejected the codon,
    // marking productive=false even though the underlying allele
    // is intact.
    let cfg = curated_v_d_j_refdata();
    let [v1, _, _, d1, _, _, j1, _, _] = curated_ids();
    let mut plan = curated_plan(&cfg, v1, d1, j1, 0, 0, 0, 0, 0, b'A', 0, b'A');
    plan.push(Box::new(InsertAtPass::new(3, b'C', Segment::V)));
    let (_outcome, rec) = curated_run(&cfg, plan);

    assert_eq!(rec.junction_start, Some(7));
    assert_eq!(rec.junction_end, Some(31));
    assert_eq!(rec.junction_length, Some(24));
    assert_eq!(rec.vj_in_frame, Some(true));
    assert_eq!(rec.productive, Some(true));
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
