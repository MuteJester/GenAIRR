use super::*;

// ──────────────────────────────────────────────────────────────
// V right-boundary extension into NP1.
//
// The acceptance criterion from the design doc: a trimmed V suffix
// can be recreated by NP1 bases, and when that happens the live
// `v_call` should narrow to the allele(s) that the recreated bytes
// distinguish. Under the conservative extension policy,
// extension proceeds one byte at a time only while the byte strictly
// narrows the current max-score tie set. Once the call set has
// collapsed to a single allele, further matching bytes do NOT extend
// the boundary — there is no ambiguity left to resolve.
// ──────────────────────────────────────────────────────────────

/// Build a VDJ refdata holding two V alleles that share a 9-base
/// prefix and differ in their 3-base suffix; the J allele is a
/// minimal stub (the fixture stops after NP1 generation, so D / J
/// pools are unused at runtime). Returns the refdata and the two
/// V allele ids in declaration order.
fn v_extension_refdata() -> (RefDataConfig, AlleleId, AlleleId) {
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
    // Minimal J entry — required by ChainType::Vdj construction
    // (refdata builder does not enforce this) but never executed
    // because the plan stops after NP1 generation.
    let _ = cfg.j_pool.push(Allele {
        name: "J*01".into(),
        gene: "J".into(),
        seq: b"AAA".to_vec(),
        segment: Segment::J,
        anchor: None,
    });
    (cfg, v01, v02)
}

/// Build the standard V-extension plan:
///   sample(V) → trim(V_3, by) → assemble(V) → generate(NP1, len, base)
/// `np_len` and `np_base` control the NP1 region's content
/// deterministically.
fn v_extension_plan(
    cfg: &RefDataConfig,
    sampled_v: AlleleId,
    v_trim_3: i64,
    np_len: i64,
    np_base: u8,
) -> PassPlan {
    let mut plan = PassPlan::new();
    plan.push(Box::new(SampleAllelePass::new(
        Segment::V,
        Box::new(AllelePoolDist::restricted_uniform(
            &cfg.v_pool,
            vec![sampled_v],
        )),
    )));
    plan.push(Box::new(TrimPass::new(
        Segment::V,
        TrimEnd::Three,
        Box::new(EmpiricalLengthDist::from_pairs(vec![(v_trim_3, 1.0)])),
    )));
    plan.push(Box::new(AssembleSegmentPass::new(Segment::V)));
    plan.push(Box::new(GenerateNPPass::new(
        Segment::Np1,
        Box::new(EmpiricalLengthDist::from_pairs(vec![(np_len, 1.0)])),
        Box::new(ConstBaseDist(np_base)),
    )));
    plan
}

#[test]
fn v_call_shrinks_when_np1_recreates_trimmed_suffix() {
    // V*01 = AAACCCGGG TTT (suffix TTT distinguishes it).
    // V*02 = AAACCCGGG AAA (suffix AAA distinguishes it).
    // Sample V*01, trim 3' by 3 → assembled V is AAACCCGGG.
    // Both alleles match the assembled bases at every position →
    // post-assemble v_call = {V*01, V*02}.
    // Then GenerateNP(NP1, length=3, 'T') puts TTT right after V.
    // V right-extension into NP1 walks T → matches V*01 ref pos 9
    // (V*01[9]='T'); V*02[9]='A' → V*02 drops out, narrowing
    // {V01,V02}→{V01}. Under conservative extension this is the
    // last narrowing step: position 10 has only V*01 in the tie
    // set, so the next byte cannot narrow further and the walker
    // stops. ref_end therefore reaches 10 (not 12), and seq_end
    // advances by exactly 1 NP1 base.
    let (cfg, v01, _v02) = v_extension_refdata();

    let plan = v_extension_plan(&cfg, v01, 3, 3, b'T');
    let compiled = CompiledSimulator::compile(&plan, Some(&cfg), None, ExecutionPolicy::Permissive)
        .expect("plan should compile");
    let outcome = compiled.run_one(0).expect("plan should run");

    // revisions: [initial, post-sample, post-trim, post-assemble,
    //             post-np1].
    assert_eq!(outcome.revisions.len(), 5);

    // After AssembleSegment(V) but before NP1: both alleles match.
    let post_assemble_call = outcome.revisions[3]
        .live_calls
        .as_ref()
        .expect("live calls populated after assembly")
        .get(Segment::V)
        .cloned()
        .expect("V live call exists after assembly");
    let mut post_assemble_ids = post_assemble_call.allele_call.to_ids();
    post_assemble_ids.sort_by_key(|id| id.index());
    assert_eq!(
        post_assemble_ids,
        vec![v01, AlleleId::new(1)],
        "post-assemble v_call should hold both indistinguishable V alleles"
    );

    // After NP1 generation: the right-extension walk has
    // narrowed to V*01 only.
    let post_np1_call = outcome.revisions[4]
        .live_calls
        .as_ref()
        .expect("live calls still populated after NP1")
        .get(Segment::V)
        .cloned()
        .expect("V live call exists after NP1");
    assert_eq!(post_np1_call.allele_call.to_ids(), vec![v01]);

    // The hypothesis must record the elastic boundary so downstream
    // tooling can tell that the right end was extended past the
    // structural V region.
    let elastic_hypothesis = post_np1_call
        .hypotheses
        .iter()
        .find(|h| {
            h.flags
                .contains(crate::live_call::HypothesisFlags::BOUNDARY_ELASTIC)
        })
        .expect("at least one hypothesis should be flagged BOUNDARY_ELASTIC");
    // Under conservative extension, ref_end advances past the
    // (post-trim) V region's ref_end of 9 by exactly 1: the first
    // NP1 byte 'T' narrows {V01,V02}→{V01}; the next byte cannot
    // narrow further (tie set is now {V01} only) so the walker
    // stops. ref_end therefore reaches 10, not the allele end 12.
    assert_eq!(
        elastic_hypothesis.ref_end, 10,
        "conservative extension stops once the call set narrows to a single allele"
    );
    // seq_end advanced from V end (9) by 1 NP1 byte = 10. The
    // remaining 2 NP1 bytes (seq positions 10, 11) stay outside
    // the V hypothesis.
    assert_eq!(
        elastic_hypothesis.seq_end, 10,
        "seq_end advances by the single NP1 byte that narrowed the call"
    );
}

#[test]
fn v_call_stays_widened_when_np1_does_not_match_any_allele() {
    // Same fixture as above but NP1 emits 'C', which is NOT the
    // next ref base for either allele (V*01[9]='T', V*02[9]='A').
    // The extension walk halts immediately — the live call should
    // remain the post-assemble widened {V*01, V*02}.
    let (cfg, v01, v02) = v_extension_refdata();

    let plan = v_extension_plan(&cfg, v01, 3, 3, b'C');
    let compiled = CompiledSimulator::compile(&plan, Some(&cfg), None, ExecutionPolicy::Permissive)
        .expect("plan should compile");
    let outcome = compiled.run_one(0).expect("plan should run");

    let post_np1_call = outcome
        .final_simulation()
        .live_calls
        .as_ref()
        .expect("live calls populated")
        .get(Segment::V)
        .cloned()
        .expect("V live call exists");
    let mut ids = post_np1_call.allele_call.to_ids();
    ids.sort_by_key(|id| id.index());
    assert_eq!(ids, vec![v01, v02]);

    // No extension → no BOUNDARY_ELASTIC flag.
    for h in &post_np1_call.hypotheses {
        assert!(
            !h.flags
                .contains(crate::live_call::HypothesisFlags::BOUNDARY_ELASTIC),
            "no NP1 base extended; BOUNDARY_ELASTIC must NOT be set"
        );
    }
}

#[test]
fn v_call_partially_extends_when_np1_matches_only_a_prefix() {
    // V*01 = AAACCCGGG TTT, V*02 = AAACCCGGG AAA. Trim V_3 by 3 →
    // assembled V is AAACCCGGG. NP1 emits 'T' for length 5 → bases
    // T-T-T-?-?. Under conservative extension, the first NP1 byte
    // narrows {V01,V02}→{V01}. The second byte cannot narrow the
    // tie set further (the set is already {V01}), so the walker
    // stops at ref_pos=10 / seq_end=10. The remaining 4 NP1 bytes
    // stay outside the V hypothesis.
    let (cfg, v01, _v02) = v_extension_refdata();

    let plan = v_extension_plan(&cfg, v01, 3, 5, b'T');
    let compiled = CompiledSimulator::compile(&plan, Some(&cfg), None, ExecutionPolicy::Permissive)
        .expect("plan should compile");
    let outcome = compiled.run_one(0).expect("plan should run");

    let post_np1_call = outcome
        .final_simulation()
        .live_calls
        .as_ref()
        .expect("live calls populated")
        .get(Segment::V)
        .cloned()
        .expect("V live call exists");
    assert_eq!(post_np1_call.allele_call.to_ids(), vec![v01]);

    let elastic_hypothesis = post_np1_call
        .hypotheses
        .iter()
        .find(|h| {
            h.flags
                .contains(crate::live_call::HypothesisFlags::BOUNDARY_ELASTIC)
        })
        .expect("hypothesis should be flagged BOUNDARY_ELASTIC");
    assert_eq!(
        elastic_hypothesis.ref_end, 10,
        "conservative extension stops at the first byte where the call set is fully resolved"
    );
    // seq_end advances from V_region.end (9) by the single NP1
    // byte that narrowed the call = 10. The remaining 4 NP1 bytes
    // (seq positions 10..14) stay outside the V hypothesis.
    assert_eq!(elastic_hypothesis.seq_end, 10);
}

#[test]
fn v_call_extension_no_op_when_no_trim() {
    // Sanity: with no V_3 trim, V's region already covers the
    // whole allele. There is no ref position past V's end, so the
    // extension walk has nothing to do — the call should remain a
    // singleton {sampled} and BOUNDARY_ELASTIC should NOT be set.
    let (cfg, v01, _v02) = v_extension_refdata();

    let plan = v_extension_plan(&cfg, v01, 0, 3, b'T');
    let compiled = CompiledSimulator::compile(&plan, Some(&cfg), None, ExecutionPolicy::Permissive)
        .expect("plan should compile");
    let outcome = compiled.run_one(0).expect("plan should run");

    let post_np1_call = outcome
        .final_simulation()
        .live_calls
        .as_ref()
        .expect("live calls populated")
        .get(Segment::V)
        .cloned()
        .expect("V live call exists");
    assert_eq!(post_np1_call.allele_call.to_ids(), vec![v01]);
    for h in &post_np1_call.hypotheses {
        assert!(
            !h.flags
                .contains(crate::live_call::HypothesisFlags::BOUNDARY_ELASTIC),
            "no extension should occur when V already covers full allele"
        );
    }
}

#[test]
fn append_region_np1_bumps_live_call_version() {
    // Plumbing check: every successful refresh bumps the live-call
    // evidence_version. AppendRegion(Np1) must trigger a V refresh
    // and therefore advance the version past the post-assemble
    // value. (Used to guard against accidentally removing the
    // hook from `apply_live_call_updates`.)
    let (cfg, v01, _v02) = v_extension_refdata();
    let plan = v_extension_plan(&cfg, v01, 3, 3, b'T');

    let compiled = CompiledSimulator::compile(&plan, Some(&cfg), None, ExecutionPolicy::Permissive)
        .expect("plan should compile");
    let outcome = compiled.run_one(0).expect("plan should run");

    let post_assemble_version = outcome.revisions[3]
        .live_calls
        .as_ref()
        .expect("post-assemble live calls present")
        .version;
    let post_np1_version = outcome.revisions[4]
        .live_calls
        .as_ref()
        .expect("post-np1 live calls present")
        .version;
    assert!(
        post_np1_version > post_assemble_version,
        "AppendRegion(Np1) must bump live-call version: \
             {post_np1_version} should be > {post_assemble_version}"
    );
}
