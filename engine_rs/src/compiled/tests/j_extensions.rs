use super::*;

// ──────────────────────────────────────────────────────────────
// J left-boundary extension into the immediately-
// preceding NP region (NP1 in VJ chains, NP2 in VDJ chains).
//
// Mirrors the V right-boundary case but for the left side of J:
// the chosen NP bases must extend backward into the J allele's
// reference prefix.
// Acceptance: trim J 5' so live j_call widens, then have NP bases
// happen to recreate the trimmed J prefix → j_call shrinks back.
//
// extension is *conservative*. The walker takes one byte
// at a time and stops as soon as the byte fails to strictly narrow
// the current max-score tie set. In practice this means: once a
// single NP byte has collapsed the call to one allele, the walker
// halts; subsequent matching bytes do not advance ref_start.
// ──────────────────────────────────────────────────────────────

/// Build a VJ-chain refdata with two J alleles that share a 9-base
/// suffix and differ in their 3-base 5' prefix. V is a minimal
/// stub used only so the J assembly pass has a chain context.
fn j_extension_vj_refdata() -> (RefDataConfig, AlleleId, AlleleId) {
    let mut cfg = RefDataConfig::empty(ChainType::Vj);
    // Stub V — never sampled in the fixture plan.
    let _ = cfg.v_pool.push(Allele {
        name: "V*01".into(),
        gene: "V".into(),
        seq: b"AAA".to_vec(),
        segment: Segment::V,
        anchor: None,
    });
    // J*01 has prefix "TTT" + shared "ACGTACGTA". J*02 has
    // prefix "GGG" + the same shared suffix. With J_5 trim = 3,
    // the assembled J region exposes only the shared "ACGTACGTA"
    // → both alleles match equally → live j_call widens to {J*01, J*02}.
    let j01 = cfg.j_pool.push(Allele {
        name: "J*01".into(),
        gene: "J".into(),
        seq: b"TTTACGTACGTA".to_vec(),
        segment: Segment::J,
        anchor: None,
    });
    let j02 = cfg.j_pool.push(Allele {
        name: "J*02".into(),
        gene: "J".into(),
        seq: b"GGGACGTACGTA".to_vec(),
        segment: Segment::J,
        anchor: None,
    });
    (cfg, j01, j02)
}

/// Plan that:
///   sample(V, single allele) → assemble(V) →
///     generate(NP1, len, base) → sample(J, single allele) →
///     trim(J_5, by) → assemble(J)
///
/// The V is a 3-base stub so it doesn't dominate test reasoning.
/// `np_len` and `np_base` control NP1 contents deterministically;
/// `j_trim_5` controls how many J prefix bases are stripped before
/// J is assembled.
fn j_extension_plan(
    cfg: &RefDataConfig,
    sampled_j: AlleleId,
    j_trim_5: i64,
    np_len: i64,
    np_base: u8,
) -> PassPlan {
    let v_id = AlleleId::new(0);
    let mut plan = PassPlan::new();
    plan.push(Box::new(SampleAllelePass::new(
        Segment::V,
        Box::new(AllelePoolDist::restricted_uniform(&cfg.v_pool, vec![v_id])),
    )));
    plan.push(Box::new(AssembleSegmentPass::new(Segment::V)));
    plan.push(Box::new(GenerateNPPass::new(
        Segment::Np1,
        Box::new(EmpiricalLengthDist::from_pairs(vec![(np_len, 1.0)])),
        Box::new(ConstBaseDist(np_base)),
    )));
    plan.push(Box::new(SampleAllelePass::new(
        Segment::J,
        Box::new(AllelePoolDist::restricted_uniform(
            &cfg.j_pool,
            vec![sampled_j],
        )),
    )));
    plan.push(Box::new(TrimPass::new(
        Segment::J,
        TrimEnd::Five,
        Box::new(EmpiricalLengthDist::from_pairs(vec![(j_trim_5, 1.0)])),
    )));
    plan.push(Box::new(AssembleSegmentPass::new(Segment::J)));
    plan
}

#[test]
fn j_call_shrinks_when_np1_recreates_trimmed_prefix_vj() {
    // J*01 = TTT ACGTACGTA, J*02 = GGG ACGTACGTA.
    // Sample J*01, trim 5' by 3 → assembled J ref window starts at
    // pos 3 (covering ACGTACGTA).
    // Both alleles match the assembled J bases → post-assemble
    // j_call = {J*01, J*02}.
    // NP1 bases = TTT (length 3, all 'T'). The walker checks NP1's
    // RIGHTMOST base first against ref pos 2 (J*01[2]='T' ✓,
    // J*02[2]='G' ✗) → tie set {J01,J02}→{J01}. Under conservative
    // extension, that single narrowing step is enough; bytes at
    // ref pos 1 and 0 cannot narrow the singleton {J01} further
    // and the walker halts. ref_start advances from 3 to 2 only.
    let (cfg, j01, _j02) = j_extension_vj_refdata();

    let plan = j_extension_plan(&cfg, j01, 3, 3, b'T');
    let compiled = CompiledSimulator::compile(&plan, Some(&cfg), None, ExecutionPolicy::Permissive)
        .expect("plan should compile");
    let outcome = compiled.run_one(0).expect("plan should run");

    // J refresh fires on AssembleSegment(J) which is the last pass.
    // The final revision should hold the post-extension j_call.
    let final_call = outcome
        .final_simulation()
        .live_calls
        .as_ref()
        .expect("live calls populated")
        .get(Segment::J)
        .cloned()
        .expect("J live call exists after assembly");
    assert_eq!(final_call.allele_call.to_ids(), vec![j01]);

    let elastic = final_call
        .hypotheses
        .iter()
        .find(|h| {
            h.flags
                .contains(crate::live_call::HypothesisFlags::BOUNDARY_ELASTIC)
        })
        .expect("J hypothesis should be flagged BOUNDARY_ELASTIC");
    // Under conservative extension, ref_start advances from 3 to
    // 2 (one narrowing byte) — not all the way back to 0.
    assert_eq!(elastic.ref_start, 2);
}

#[test]
fn j_call_stays_widened_when_np1_does_not_recreate_prefix_vj() {
    // Same fixture as above, but NP1 emits 'C' which is not the
    // expected base for either allele's prefix at any position.
    // The extension halts immediately; j_call stays widened.
    let (cfg, j01, j02) = j_extension_vj_refdata();

    let plan = j_extension_plan(&cfg, j01, 3, 3, b'C');
    let compiled = CompiledSimulator::compile(&plan, Some(&cfg), None, ExecutionPolicy::Permissive)
        .expect("plan should compile");
    let outcome = compiled.run_one(0).expect("plan should run");

    let final_call = outcome
        .final_simulation()
        .live_calls
        .as_ref()
        .expect("live calls populated")
        .get(Segment::J)
        .cloned()
        .expect("J live call exists");
    let mut ids = final_call.allele_call.to_ids();
    ids.sort_by_key(|id| id.index());
    assert_eq!(ids, vec![j01, j02]);
    for h in &final_call.hypotheses {
        assert!(
            !h.flags
                .contains(crate::live_call::HypothesisFlags::BOUNDARY_ELASTIC),
            "no extension occurred → BOUNDARY_ELASTIC must not be set"
        );
    }
}

#[test]
fn j_call_partially_extends_when_np1_matches_only_a_suffix_of_prefix() {
    // J*01 = TTT ACGTACGTA, J*02 = GGG ACGTACGTA.
    // Trim J_5 by 3, NP1 length 5, base 'T'. NP1 = T-T-T-T-T.
    // Walker checks rightmost NP1 byte vs ref_start - 1 = 2 first.
    //   NP1[4]='T' vs J*01[2]='T' ✓, J*02[2]='G' ✗ → tie set
    //     {J01,J02}→{J01}.
    //   NP1[3]: tie set is now {J01}; conservative extension
    //     declines (no allele in the tie set could MISS the byte
    //     to narrow it further). Walker halts.
    // Four of the five NP1 bases stay outside the J hypothesis.
    let (cfg, j01, _j02) = j_extension_vj_refdata();

    let plan = j_extension_plan(&cfg, j01, 3, 5, b'T');
    let compiled = CompiledSimulator::compile(&plan, Some(&cfg), None, ExecutionPolicy::Permissive)
        .expect("plan should compile");
    let outcome = compiled.run_one(0).expect("plan should run");

    let final_call = outcome
        .final_simulation()
        .live_calls
        .as_ref()
        .expect("live calls populated")
        .get(Segment::J)
        .cloned()
        .expect("J live call exists");
    assert_eq!(final_call.allele_call.to_ids(), vec![j01]);

    let elastic = final_call
        .hypotheses
        .iter()
        .find(|h| {
            h.flags
                .contains(crate::live_call::HypothesisFlags::BOUNDARY_ELASTIC)
        })
        .expect("J hypothesis should be flagged BOUNDARY_ELASTIC");
    assert_eq!(
        elastic.ref_start, 2,
        "conservative extension consumes only the single byte that narrows the tie set"
    );
    // J's seq_start moved left by 1 (1 NP1 base consumed). The
    // remaining 4 NP1 bases stay outside the J hypothesis.
    // Pool layout: [V (3 bases) 0..3] [NP1 (5 bases) 3..8] [J (9
    // bases) 8..17]. Post-extension J seq_start = 8 - 1 = 7.
    assert_eq!(elastic.seq_start, 7);
    assert_eq!(elastic.seq_end, 17);
}

#[test]
fn j_call_extension_no_op_when_no_trim() {
    // No J_5 trim → J's ref_start is already 0 → there is no ref
    // position to the left to extend into. The walker halts on
    // the very first iteration because `extended_ref_start == 0`.
    let (cfg, j01, _j02) = j_extension_vj_refdata();

    let plan = j_extension_plan(&cfg, j01, 0, 3, b'T');
    let compiled = CompiledSimulator::compile(&plan, Some(&cfg), None, ExecutionPolicy::Permissive)
        .expect("plan should compile");
    let outcome = compiled.run_one(0).expect("plan should run");

    let final_call = outcome
        .final_simulation()
        .live_calls
        .as_ref()
        .expect("live calls populated")
        .get(Segment::J)
        .cloned()
        .expect("J live call exists");
    assert_eq!(final_call.allele_call.to_ids(), vec![j01]);
    for h in &final_call.hypotheses {
        assert!(
            !h.flags
                .contains(crate::live_call::HypothesisFlags::BOUNDARY_ELASTIC),
            "no extension when J's ref_start is already 0"
        );
    }
}

#[test]
fn j_left_extension_works_for_vdj_chain_via_np2() {
    // Verify the chain-agnostic neighbour lookup: in a VDJ chain,
    // J's left neighbour is NP2 (not NP1). Build a minimal VDJ
    // refdata, plan V→NP1→D→NP2→trim(J_5)→J, and check that NP2
    // bases drive the J left-extension.
    let mut cfg = RefDataConfig::empty(ChainType::Vdj);
    let v_id = cfg.v_pool.push(Allele {
        name: "V*01".into(),
        gene: "V".into(),
        seq: b"AAA".to_vec(),
        segment: Segment::V,
        anchor: None,
    });
    let d_id = cfg.d_pool.push(Allele {
        name: "D*01".into(),
        gene: "D".into(),
        seq: b"CCC".to_vec(),
        segment: Segment::D,
        anchor: None,
    });
    let j01 = cfg.j_pool.push(Allele {
        name: "J*01".into(),
        gene: "J".into(),
        seq: b"TTTACGTACGTA".to_vec(),
        segment: Segment::J,
        anchor: None,
    });
    let _j02 = cfg.j_pool.push(Allele {
        name: "J*02".into(),
        gene: "J".into(),
        seq: b"GGGACGTACGTA".to_vec(),
        segment: Segment::J,
        anchor: None,
    });

    let mut plan = PassPlan::new();
    plan.push(Box::new(SampleAllelePass::new(
        Segment::V,
        Box::new(AllelePoolDist::restricted_uniform(&cfg.v_pool, vec![v_id])),
    )));
    plan.push(Box::new(SampleAllelePass::new(
        Segment::D,
        Box::new(AllelePoolDist::restricted_uniform(&cfg.d_pool, vec![d_id])),
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
        Box::new(EmpiricalLengthDist::from_pairs(vec![(3, 1.0)])),
        Box::new(ConstBaseDist(b'T')),
    )));
    plan.push(Box::new(SampleAllelePass::new(
        Segment::J,
        Box::new(AllelePoolDist::restricted_uniform(&cfg.j_pool, vec![j01])),
    )));
    plan.push(Box::new(TrimPass::new(
        Segment::J,
        TrimEnd::Five,
        Box::new(EmpiricalLengthDist::from_pairs(vec![(3, 1.0)])),
    )));
    plan.push(Box::new(AssembleSegmentPass::new(Segment::J)));

    let compiled = CompiledSimulator::compile(&plan, Some(&cfg), None, ExecutionPolicy::Permissive)
        .expect("VDJ plan should compile");
    let outcome = compiled.run_one(0).expect("plan should run");
    let final_call = outcome
        .final_simulation()
        .live_calls
        .as_ref()
        .expect("live calls populated")
        .get(Segment::J)
        .cloned()
        .expect("J live call exists");
    assert_eq!(
        final_call.allele_call.to_ids(),
        vec![j01],
        "VDJ chain: NP2 bases TTT should narrow J back to J*01"
    );
    let elastic = final_call
        .hypotheses
        .iter()
        .find(|h| {
            h.flags
                .contains(crate::live_call::HypothesisFlags::BOUNDARY_ELASTIC)
        })
        .expect("J hypothesis should be flagged BOUNDARY_ELASTIC");
    // Conservative extension: the first NP2 byte narrows
    // {J01,J02}→{J01}. The remaining bytes cannot narrow the
    // singleton further, so ref_start stops at 2 (not 0).
    assert_eq!(elastic.ref_start, 2);
}

// ──────────────────────────────────────────────────────────────
