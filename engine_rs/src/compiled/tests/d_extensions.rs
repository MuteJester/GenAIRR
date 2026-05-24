use super::*;

//
// D is the "hardest segment" per the design doc — it has elastic
// boundaries on both sides. We test both directions independently
// and then together. Pipeline order in VDJ is:
//
//   sample(V) → sample(D) → assemble(V) → np1 → assemble(D) → np2 → ...
//
// - D LEFT extension fires automatically on AssembleSegment(D)
//   because NP1 already exists by then.
// - D RIGHT extension needs the new AppendRegion(Np2) hook
//   (NP2 is generated AFTER D is assembled).
//
// extension into NP is *conservative*. It proceeds one
// byte at a time only while the byte strictly narrows the current
// max-score tie set. If the primary structural walk has already
// resolved the call to a single allele, conservative extension
// declines to walk further; ref_start / ref_end stay at the
// structural boundary.
// ──────────────────────────────────────────────────────────────

/// Build a VDJ refdata with two D alleles that share a 5-base core
/// "CCCCC" and differ in their 3-base 5' prefix and 3-base 3'
/// suffix.
///
///   D*01 = TTT CCCCC AAA  (unique prefix TTT, unique suffix AAA)
///   D*02 = GGG CCCCC TGT  (different prefix GGG, different suffix TGT)
///
/// V and J pools are minimal stubs — needed so the chain can
/// assemble end-to-end but irrelevant to the D-call assertions.
fn d_extension_refdata() -> (RefDataConfig, AlleleId, AlleleId) {
    let mut cfg = RefDataConfig::empty(ChainType::Vdj);
    let _ = cfg.v_pool.push(Allele {
        name: "V*01".into(),
        gene: "V".into(),
        seq: b"AAA".to_vec(),
        segment: Segment::V,
        anchor: None,
    });
    let d01 = cfg.d_pool.push(Allele {
        name: "D*01".into(),
        gene: "D".into(),
        seq: b"TTTCCCCCAAA".to_vec(),
        segment: Segment::D,
        anchor: None,
    });
    let d02 = cfg.d_pool.push(Allele {
        name: "D*02".into(),
        gene: "D".into(),
        seq: b"GGGCCCCCTGT".to_vec(),
        segment: Segment::D,
        anchor: None,
    });
    let _ = cfg.j_pool.push(Allele {
        name: "J*01".into(),
        gene: "J".into(),
        seq: b"AAA".to_vec(),
        segment: Segment::J,
        anchor: None,
    });
    (cfg, d01, d02)
}

/// Build the standard D-extension plan:
///
///   sample(V, stub) → sample(D, sampled_d) → assemble(V) →
///   np1(np1_len, np1_base) → trim(D_5, by) → trim(D_3, by) →
///   assemble(D) → np2(np2_len, np2_base) → sample(J, stub) →
///   assemble(J)
///
/// V and J are stubs; the assertions only inspect the D live call.
fn d_extension_plan(
    cfg: &RefDataConfig,
    sampled_d: AlleleId,
    d_trim_5: i64,
    d_trim_3: i64,
    np1_len: i64,
    np1_base: u8,
    np2_len: i64,
    np2_base: u8,
) -> PassPlan {
    let v_id = AlleleId::new(0);
    let j_id = AlleleId::new(0);
    let mut plan = PassPlan::new();
    plan.push(Box::new(SampleAllelePass::new(
        Segment::V,
        Box::new(AllelePoolDist::restricted_uniform(&cfg.v_pool, vec![v_id])),
    )));
    plan.push(Box::new(SampleAllelePass::new(
        Segment::D,
        Box::new(AllelePoolDist::restricted_uniform(
            &cfg.d_pool,
            vec![sampled_d],
        )),
    )));
    plan.push(Box::new(TrimPass::new(
        Segment::D,
        TrimEnd::Five,
        Box::new(EmpiricalLengthDist::from_pairs(vec![(d_trim_5, 1.0)])),
    )));
    plan.push(Box::new(TrimPass::new(
        Segment::D,
        TrimEnd::Three,
        Box::new(EmpiricalLengthDist::from_pairs(vec![(d_trim_3, 1.0)])),
    )));
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
    plan.push(Box::new(SampleAllelePass::new(
        Segment::J,
        Box::new(AllelePoolDist::restricted_uniform(&cfg.j_pool, vec![j_id])),
    )));
    plan.push(Box::new(AssembleSegmentPass::new(Segment::J)));
    plan
}

#[test]
fn d_call_shrinks_when_np1_recreates_trimmed_prefix() {
    // Sample D*01 (TTTCCCCCAAA), trim D_5 by 3, D_3 by 3.
    // Assembled D ref window starts at pos 3, covers 3..8 (CCCCC).
    // BOTH D alleles match the assembled bases at those positions →
    // post-assemble d_call = {D*01, D*02}.
    // NP1 = TTT (length 3, all 'T'). The walker checks NP1's
    // RIGHTMOST base first against ref pos 2 (D*01[2]='T',
    // D*02[2]='G'). D*02 drops out; the tie set narrows from
    // {D01,D02} to {D01}.
    //
    // Under conservative extension, this is the only narrowing
    // step: at ref pos 1 the tie set is already {D01}, so no
    // further byte can narrow it. ref_start advances from 3 to 2
    // (one byte recovered), not to 0.
    let (cfg, d01, _d02) = d_extension_refdata();

    // Trim both ends by 3 so the post-assemble tie set is genuinely
    // {D*01, D*02} (both match the CCCCC core).
    let plan = d_extension_plan(&cfg, d01, 3, 3, 3, b'T', 0, b'A');
    let compiled = CompiledSimulator::compile(&plan, Some(&cfg), None, ExecutionPolicy::Permissive)
        .expect("plan should compile");
    let outcome = compiled.run_one(0).expect("plan should run");

    let final_call = outcome
        .final_simulation()
        .segment_calls
        .get(Segment::D)
        .cloned()
        .expect("D live call exists");
    assert_eq!(final_call.allele_call.to_ids(), vec![d01]);

    let elastic = final_call
        .hypotheses
        .iter()
        .find(|h| {
            h.flags
                .contains(crate::live_call::HypothesisFlags::BOUNDARY_ELASTIC)
        })
        .expect("D hypothesis should be flagged BOUNDARY_ELASTIC");
    // Under conservative extension, ref_start moves from 3 (the
    // post-trim structural start) to 2 — the single NP1 byte that
    // narrowed the call from {D01,D02} to {D01}. It does NOT reach
    // 0 (would require greedy extension into bytes that no longer
    // narrow anything).
    assert_eq!(elastic.ref_start, 2);
}

#[test]
fn d_call_shrinks_when_np2_recreates_trimmed_suffix() {
    // Sample D*01 (TTTCCCCCAAA), trim D_5 by 0, D_3 by 3.
    // Assembled D covers 0..8 (TTTCCCCC). Both alleles match the
    // assembled bases at those positions only after first checking
    // pos 0 (D*01[0]='T', D*02[0]='G' → D*02 drops out at primary
    // walk!). Hmm — that means d_call is already {D*01} at the
    // primary walk. So this test really checks RIGHT extension on
    // an already-singleton call: NP2 should still narrow ref_end
    // and set BOUNDARY_ELASTIC.
    //
    // To get widening via D_3 trim alone, we'd need alleles that
    // share the prefix. Let's adjust:
    let mut cfg = RefDataConfig::empty(ChainType::Vdj);
    let _ = cfg.v_pool.push(Allele {
        name: "V*01".into(),
        gene: "V".into(),
        seq: b"AAA".to_vec(),
        segment: Segment::V,
        anchor: None,
    });
    let d01 = cfg.d_pool.push(Allele {
        name: "D*01".into(),
        gene: "D".into(),
        seq: b"AAACCCCCTTT".to_vec(),
        segment: Segment::D,
        anchor: None,
    });
    let _d02 = cfg.d_pool.push(Allele {
        name: "D*02".into(),
        gene: "D".into(),
        seq: b"AAACCCCCGGG".to_vec(),
        segment: Segment::D,
        anchor: None,
    });
    let _ = cfg.j_pool.push(Allele {
        name: "J*01".into(),
        gene: "J".into(),
        seq: b"AAA".to_vec(),
        segment: Segment::J,
        anchor: None,
    });
    // Sample D*01, trim D_3 by 3, D_5 by 0. Assembled D = AAACCCCC
    // (ref pos 0..8). Both alleles match → live call = {D*01, D*02}.
    // NP2 = TTT (length 3). The walker extends right: NP2[0]='T'
    // vs ref pos 8 (D*01[8]='T' ✓, D*02[8]='G' ✗) → D*02 out, tie
    // set narrows {D01,D02}→{D01}. Under conservative extension
    // this is the only narrowing step; positions 9 and 10 cannot
    // narrow the singleton further so the walker stops. Final
    // d_call = {D*01}, ref_end = 9 (not 11).
    let plan = d_extension_plan(&cfg, d01, 0, 3, 0, b'A', 3, b'T');
    let compiled = CompiledSimulator::compile(&plan, Some(&cfg), None, ExecutionPolicy::Permissive)
        .expect("plan should compile");
    let outcome = compiled.run_one(0).expect("plan should run");

    let final_call = outcome
        .final_simulation()
        .segment_calls
        .get(Segment::D)
        .cloned()
        .expect("D live call exists");
    assert_eq!(final_call.allele_call.to_ids(), vec![d01]);

    let elastic = final_call
        .hypotheses
        .iter()
        .find(|h| {
            h.flags
                .contains(crate::live_call::HypothesisFlags::BOUNDARY_ELASTIC)
        })
        .expect("D hypothesis should be flagged BOUNDARY_ELASTIC");
    // Under conservative extension, ref_end advances by exactly
    // the single NP2 byte that narrowed {D01,D02}→{D01}. It does
    // not reach the full allele length (11).
    assert_eq!(elastic.ref_end, 9);
}

#[test]
fn d_call_shrinks_via_both_sides_simultaneously() {
    // D*01 = TTT CCCCC AAA, D*02 = GGG CCCCC TGT. Trim D_5 = 3,
    // D_3 = 3. Assembled D = CCCCC (ref pos 3..8). Both match
    // → d_call = {D*01, D*02}.
    //
    // NP1 = TTT (length 3) — first byte at ref pos 2 narrows
    // {D01,D02}→{D01} (D*01[2]='T' vs D*02[2]='G'). Conservative
    // extension stops there: ref_start advances from 3 to 2.
    // NP2 = AAA (length 3) — the right walk starts with a
    // singleton tie set {D01}, so no NP2 byte can narrow it
    // further. ref_end stays at the structural boundary 8.
    //
    // Under conservative the call still collapses to {D*01} (the
    // single NP1 byte was enough), but the recovered boundaries
    // only advance by the bytes that actually narrowed the tie
    // set: ref_start = 2 (not 0), ref_end = 8 (not 11).
    let (cfg, d01, _d02) = d_extension_refdata();
    let plan = d_extension_plan(&cfg, d01, 3, 3, 3, b'T', 3, b'A');

    let compiled = CompiledSimulator::compile(&plan, Some(&cfg), None, ExecutionPolicy::Permissive)
        .expect("plan should compile");
    let outcome = compiled.run_one(0).expect("plan should run");

    let final_call = outcome
        .final_simulation()
        .segment_calls
        .get(Segment::D)
        .cloned()
        .expect("D live call exists");
    assert_eq!(final_call.allele_call.to_ids(), vec![d01]);

    let elastic = final_call
        .hypotheses
        .iter()
        .find(|h| {
            h.flags
                .contains(crate::live_call::HypothesisFlags::BOUNDARY_ELASTIC)
        })
        .expect("D hypothesis should be flagged BOUNDARY_ELASTIC");
    assert_eq!(elastic.ref_start, 2);
    assert_eq!(elastic.ref_end, 8);
}

#[test]
fn d_call_stays_widened_when_neither_np_matches() {
    // NP1 emits 'C' which is neither D*01[2]='T' nor D*02[2]='G'.
    // NP2 emits 'C' which is neither D*01[8]='A' nor D*02[8]='T'.
    // Both extension walks halt immediately; d_call stays widened.
    let (cfg, d01, d02) = d_extension_refdata();
    let plan = d_extension_plan(&cfg, d01, 3, 3, 3, b'C', 3, b'C');

    let compiled = CompiledSimulator::compile(&plan, Some(&cfg), None, ExecutionPolicy::Permissive)
        .expect("plan should compile");
    let outcome = compiled.run_one(0).expect("plan should run");

    let final_call = outcome
        .final_simulation()
        .segment_calls
        .get(Segment::D)
        .cloned()
        .expect("D live call exists");
    let mut ids = final_call.allele_call.to_ids();
    ids.sort_by_key(|id| id.index());
    assert_eq!(ids, vec![d01, d02]);
    for h in &final_call.hypotheses {
        assert!(
            !h.flags
                .contains(crate::live_call::HypothesisFlags::BOUNDARY_ELASTIC),
            "no extension occurred → BOUNDARY_ELASTIC must not be set"
        );
    }
}

#[test]
fn append_region_np2_bumps_d_live_call_version() {
    // Plumbing check: AppendRegion(Np2) must trigger D refresh,
    // bumping the live-call evidence_version. Catches accidental
    // removal of the hook from `apply_live_call_updates`.
    let (cfg, d01, _d02) = d_extension_refdata();
    let plan = d_extension_plan(&cfg, d01, 3, 3, 3, b'T', 3, b'A');

    let compiled = CompiledSimulator::compile(&plan, Some(&cfg), None, ExecutionPolicy::Permissive)
        .expect("plan should compile");
    let outcome = compiled.run_one(0).expect("plan should run");

    // Find indices of the AssembleSegment(D) and AppendRegion(Np2)
    // pass commits by walking pass_names. revisions are
    // [initial, post-pass-0, post-pass-1, ...].
    let pass_names = &outcome.pass_names;
    let assemble_d_idx = pass_names
        .iter()
        .position(|n| n == "assemble.d")
        .expect("plan must include assemble.d");
    let np2_idx = pass_names
        .iter()
        .position(|n| n == "generate_np.np2")
        .expect("plan must include generate_np.np2");
    assert!(
        np2_idx > assemble_d_idx,
        "np2 generation should happen after D assembly"
    );

    let post_assemble_d_version = outcome.revisions[assemble_d_idx + 1]
        .segment_calls
        .version;
    let post_np2_version = outcome.revisions[np2_idx + 1]
        .segment_calls
        .version;
    assert!(
        post_np2_version > post_assemble_d_version,
        "AppendRegion(Np2) must bump live-call version: \
             {post_np2_version} should be > {post_assemble_d_version}"
    );
}
