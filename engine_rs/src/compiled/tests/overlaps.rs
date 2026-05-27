use super::*;

// ──────────────────────────────────────────────────────────────
// Overlap hypotheses.
//
// When V's allele suffix happens to match D's leading bases (or
// D's allele suffix matches J's leading bases, or symmetric
// left-side cases), the live graph must retain BOTH placements
// internally. We surface this via:
// - the upstream segment's hypothesis growing seq_end past the
//   downstream segment's seq_start,
// - the upstream hypothesis carrying the OVERLAPS_OTHER_SEGMENT
//   flag.
//
// extension is *conservative*. Even when D's leading
// bases continue to match the upstream allele's continuation, the
// walker only consumes bytes that strictly narrow the upstream
// call's tie set. Once the call has resolved to a singleton, the
// overlap walk halts — even mid-D — and seq_end reflects just the
// bytes that contributed evidence.
// ──────────────────────────────────────────────────────────────

/// Build a VDJ refdata where V*01's distinguishing 3' suffix
/// "TTT" coincides with D*01's leading 3 bases. With V trimmed
/// 3' by 3 and an empty NP1, V's right-extension can reach into
/// D's region for an exact-equivalent overlap placement.
fn v_d_overlap_refdata() -> (RefDataConfig, AlleleId, AlleleId, AlleleId) {
    let mut cfg = RefDataConfig::empty(ChainType::Vdj);
    // Two V alleles sharing AAACCCGGG and differing in their 3'
    // trinucleotide. After V_3 trim of 3, both match the
    // structural V region equally → primary live call is widened.
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
    // D*01 starts with "TTT" — that's exactly the trimmed-off
    // V*01 suffix, so V right-extension into D produces overlap.
    let d01 = cfg.d_pool.push(Allele {
        name: "D*01".into(),
        gene: "D".into(),
        seq: b"TTTACGTAC".to_vec(),
        segment: Segment::D,
        anchor: None,
    });
    // Stub J — never executed in fixtures that stop after
    // AssembleSegment(D).
    let _ = cfg.j_pool.push(Allele {
        name: "J*01".into(),
        gene: "J".into(),
        seq: b"AAA".to_vec(),
        segment: Segment::J,
        anchor: None,
    });
    (cfg, v01, v02, d01)
}

/// Plan that ends right after assembling D, so the V live call
/// reflects the post-D state including any overlap.
fn v_overlap_plan(
    cfg: &RefDataConfig,
    sampled_v: AlleleId,
    sampled_d: AlleleId,
    v_trim_3: i64,
    np1_len: i64,
) -> PassPlan {
    let mut plan = PassPlan::new();
    plan.push(Box::new(SampleAllelePass::new(
        Segment::V,
        Box::new(AllelePoolDist::restricted_uniform(
            &cfg.v_pool,
            vec![sampled_v],
        )),
    )));
    plan.push(Box::new(SampleAllelePass::new(
        Segment::D,
        Box::new(AllelePoolDist::restricted_uniform(
            &cfg.d_pool,
            vec![sampled_d],
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
        Box::new(EmpiricalLengthDist::from_pairs(vec![(np1_len, 1.0)])),
        Box::new(ConstBaseDist(b'C')),
    )));
    plan.push(Box::new(AssembleSegmentPass::new(Segment::D)));
    plan
}

#[test]
fn v_right_overlaps_d_when_d_starts_with_v_suffix() {
    // V*01 trimmed 3' by 3 → V region = AAACCCGGG (ref 0..9).
    // V live call after AssembleSegment(V) = {V*01, V*02}.
    // After GenerateNP(NP1, length=0), nothing changes (NP1 empty).
    // After AssembleSegment(D): the cross-segment hook retriggers V.
    //   V right-extension walker enters with NP1 [9, 9), empty,
    //   loops `9..pool_len`. seq_pos = 9 is D's first base 'T'.
    //   compatible_alleles_at(ref_pos=9, base=T) → V*01[9]='T' ✓,
    //   V*02[9]='A' ✗ → tie set {V01,V02}→{V01}. Extends.
    //   At seq_pos=10 the tie set is already {V01}; conservative
    //   declines. Walker halts.
    //
    // Final V hypothesis: candidates={V*01}, seq_end=10,
    // ref_end=10, BOUNDARY_ELASTIC + OVERLAPS_OTHER_SEGMENT set.
    // D's region.start = 9. V hypothesis seq_end (10) > D start
    // (9) → 1 base of overlap.
    let (cfg, v01, _v02, d01) = v_d_overlap_refdata();
    let plan = v_overlap_plan(&cfg, v01, d01, 3, 0);

    let compiled = CompiledSimulator::compile(&plan, Some(&cfg), None, ExecutionPolicy::Permissive)
        .expect("plan should compile");
    let outcome = compiled.run_one(0).expect("plan should run");

    let final_sim = outcome.final_simulation();
    let v_call = final_sim
        .segment_calls
        .get(Segment::V)
        .cloned()
        .expect("V live call exists");
    let d_call = final_sim
        .segment_calls
        .get(Segment::D)
        .cloned()
        .expect("D live call exists");

    // V's call has shrunk back to {V*01} — the overlap walk
    // recovered the trimmed suffix from D's leading bases.
    assert_eq!(v_call.allele_call.to_ids(), vec![v01]);

    let v_hypothesis = v_call
        .hypotheses
        .iter()
        .find(|h| {
            h.flags
                .contains(crate::live_call::HypothesisFlags::BOUNDARY_ELASTIC)
        })
        .expect("V hypothesis should be flagged BOUNDARY_ELASTIC");
    assert!(
        v_hypothesis
            .flags
            .contains(crate::live_call::HypothesisFlags::OVERLAPS_OTHER_SEGMENT),
        "V hypothesis must carry OVERLAPS_OTHER_SEGMENT after \
             extending past NP1 into D's region"
    );
    assert_eq!(v_hypothesis.seq_end, 10);
    assert_eq!(v_hypothesis.ref_end, 10);

    // D's hypothesis is computed independently from D's region;
    // its seq_start should still match D.region.start = 9.
    let d_hypothesis = &d_call.hypotheses[0];
    assert_eq!(
        d_hypothesis.seq_start, 9,
        "D hypothesis is derived from D's region; seq_start unchanged"
    );

    // Internal-state property the design doc calls out:
    // "V end may be greater than D start internally when evidence
    //  supports both."
    assert!(
        v_hypothesis.seq_end > d_hypothesis.seq_start,
        "internal overlap: V seq_end ({}) should exceed D seq_start ({})",
        v_hypothesis.seq_end,
        d_hypothesis.seq_start,
    );
}

#[test]
fn v_right_does_not_overlap_when_d_does_not_match_v_suffix() {
    // Same V refdata, but D*01 starts with non-matching bases:
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
    // D starts with C — neither V allele has C at ref pos 9
    // (V*01[9]=T, V*02[9]=A) → walker halts on first attempt.
    let d01 = cfg.d_pool.push(Allele {
        name: "D*01".into(),
        gene: "D".into(),
        seq: b"CCCACGTAC".to_vec(),
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
    let plan = v_overlap_plan(&cfg, v01, d01, 3, 0);

    let compiled = CompiledSimulator::compile(&plan, Some(&cfg), None, ExecutionPolicy::Permissive)
        .expect("plan should compile");
    let outcome = compiled.run_one(0).expect("plan should run");
    let v_call = outcome
        .final_simulation()
        .segment_calls
        .get(Segment::V)
        .cloned()
        .expect("V live call exists");

    // V live call stays at the post-trim widened {V*01, V*02} —
    // no allele's ref pos 9 was 'C'.
    let mut ids = v_call.allele_call.to_ids();
    ids.sort_by_key(|id| id.index());
    assert_eq!(ids, vec![v01, v02]);
    for h in &v_call.hypotheses {
        assert!(
            !h.flags
                .contains(crate::live_call::HypothesisFlags::OVERLAPS_OTHER_SEGMENT),
            "non-matching D bases must NOT produce overlap"
        );
    }
}

#[test]
fn d_left_overlaps_v_when_v_ends_with_d_prefix() {
    // Mirror of the V→D overlap, but for D's left-extension
    // walker reaching backward into V's region.
    //
    // D*01 has prefix "GGG" — that's V*01's last 3 bases.
    // D*02 has prefix "TTT" (mismatches V's tail).
    // Trim D_5 by 3 → D's structural region drops the GGG prefix
    // and starts at ref pos 3. After AssembleSegment(D), D's
    // post-assemble live call is {D*01, D*02} (both match the
    // shared ACGTAC tail). D's left-extension walks backward
    // through the empty NP1 into V's region. V's last base 'G'
    // matches D*01[2]='G' but not D*02[2]='T', narrowing the tie
    // set to {D*01} — and the walker halts at that single byte
    // under conservative extension (ref_start = 2, not 0).
    let mut cfg = RefDataConfig::empty(ChainType::Vdj);
    let v01 = cfg.v_pool.push(Allele {
        name: "V*01".into(),
        gene: "V".into(),
        seq: b"AAACCCGGG".to_vec(),
        segment: Segment::V,
        anchor: None,
    });
    let _ = cfg.d_pool.push(Allele {
        name: "D*01".into(),
        gene: "D".into(),
        seq: b"GGGACGTAC".to_vec(),
        segment: Segment::D,
        anchor: None,
    });
    let _ = cfg.d_pool.push(Allele {
        name: "D*02".into(),
        gene: "D".into(),
        seq: b"TTTACGTAC".to_vec(),
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
    let d01 = AlleleId::new(0);

    let mut plan = PassPlan::new();
    plan.push(Box::new(SampleAllelePass::new(
        Segment::V,
        Box::new(AllelePoolDist::restricted_uniform(&cfg.v_pool, vec![v01])),
    )));
    plan.push(Box::new(SampleAllelePass::new(
        Segment::D,
        Box::new(AllelePoolDist::restricted_uniform(&cfg.d_pool, vec![d01])),
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
        Box::new(ConstBaseDist(b'C')),
    )));
    plan.push(Box::new(AssembleSegmentPass::new(Segment::D)));

    let compiled = CompiledSimulator::compile(&plan, Some(&cfg), None, ExecutionPolicy::Permissive)
        .expect("plan should compile");
    let outcome = compiled.run_one(0).expect("plan should run");

    let final_sim = outcome.final_simulation();
    let v_call = final_sim
        .segment_calls
        .get(Segment::V)
        .cloned()
        .expect("V live call exists");
    let d_call = final_sim
        .segment_calls
        .get(Segment::D)
        .cloned()
        .expect("D live call exists");

    // D's call should have BOUNDARY_ELASTIC + OVERLAPS_OTHER_SEGMENT
    // because the left-extension walker reached back through NP1
    // (empty) into V's last 3 bases (GGG) and they exactly
    // matched D*01's trimmed prefix.
    let d_hypothesis = d_call
        .hypotheses
        .iter()
        .find(|h| {
            h.flags
                .contains(crate::live_call::HypothesisFlags::BOUNDARY_ELASTIC)
        })
        .expect("D hypothesis should be flagged BOUNDARY_ELASTIC");
    assert!(
        d_hypothesis
            .flags
            .contains(crate::live_call::HypothesisFlags::OVERLAPS_OTHER_SEGMENT),
        "D's left-extension overlap into V must set OVERLAPS_OTHER_SEGMENT"
    );

    // Under conservative extension, D hypothesis seq_start moves
    // from the structural D start (9) back by exactly the single
    // narrowing byte → seq_start = 8 (V's last 'G' position).
    // ref_start advances from 3 to 2, not all the way to 0.
    let v_region_end = v_call.hypotheses[0].seq_end;
    assert!(
        d_hypothesis.seq_start < v_region_end,
        "D's left boundary should reach into V's region (D start {} < V end {})",
        d_hypothesis.seq_start,
        v_region_end,
    );
    assert_eq!(d_hypothesis.ref_start, 2);
}

#[test]
fn overlap_walker_halts_at_pool_end() {
    // Conservative-extension counterpart of the old "pool_len
    // upper-bound" safety check. Under the conservative policy
    //, the walker halts as soon as the call set
    // collapses to a singleton — long before it could possibly
    // run off the end of the pool. This test exercises that
    // halt-on-resolution behaviour:
    //
    // - Two V alleles with a shared prefix and differing on a
    //   single byte at ref pos 5 (V*01='T', V*02='A').
    // - V_3 trim of 7 → V's structural region is just the first
    //   5 bases (the shared prefix). Post-assemble tie set =
    //   {V*01, V*02}.
    // - D is long all-'T' bases. The first byte of D narrows the
    //   tie set to {V*01}; the subsequent matching 'T' bytes in
    //   D no longer narrow anything, so the walker stops mid-D.
    //
    // seq_end must therefore reflect "one byte into D", not the
    // entire D length — and certainly never beyond pool_len.
    let mut cfg = RefDataConfig::empty(ChainType::Vdj);
    let v01 = cfg.v_pool.push(Allele {
        name: "V*01".into(),
        gene: "V".into(),
        // V*01: AAAAA T AAAAA... (12 bases). First 5 'A's shared,
        // pos 5 = 'T'.
        seq: {
            let mut s = vec![b'A'; 12];
            s[5] = b'T';
            s
        },
        segment: Segment::V,
        anchor: None,
    });
    let _ = cfg.v_pool.push(Allele {
        name: "V*02".into(),
        gene: "V".into(),
        // V*02: AAAAA A AAAAA... (12 bases of 'A'). Diverges from
        // V*01 only at pos 5.
        seq: vec![b'A'; 12],
        segment: Segment::V,
        anchor: None,
    });
    let _ = cfg.d_pool.push(Allele {
        name: "D*01".into(),
        gene: "D".into(),
        seq: vec![b'T'; 5],
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
    let d01 = AlleleId::new(0);
    // Trim V_3 by 7 so V's structural region is just the first
    // 5 'A's (shared between both alleles → tie set has 2).
    let plan = v_overlap_plan(&cfg, v01, d01, 7, 0);

    let compiled = CompiledSimulator::compile(&plan, Some(&cfg), None, ExecutionPolicy::Permissive)
        .expect("plan should compile");
    let outcome = compiled.run_one(0).expect("plan should run");
    let v_call = outcome
        .final_simulation()
        .segment_calls
        .get(Segment::V)
        .cloned()
        .expect("V live call exists");

    let v_hypothesis = v_call
        .hypotheses
        .iter()
        .find(|h| {
            h.flags
                .contains(crate::live_call::HypothesisFlags::BOUNDARY_ELASTIC)
        })
        .expect("V hypothesis should be elastic");
    // Pool layout: V_region [0, 5), NP1 [5, 5), D_region [5, 10).
    // pool_len = 10. The first D byte (pool pos 5 = 'T') narrows
    // {V01,V02}→{V01} (V*01[5]='T', V*02[5]='A'). seq_end = 6.
    // The remaining D bytes can't narrow {V01} further; walker
    // halts safely without running to pool_len.
    assert_eq!(
        v_hypothesis.seq_end, 6,
        "conservative walker halts at the first byte that fully resolves the call"
    );
    assert!(v_hypothesis.seq_end < 10, "seq_end stays within pool_len");
    assert!(
        v_hypothesis
            .flags
            .contains(crate::live_call::HypothesisFlags::OVERLAPS_OTHER_SEGMENT),
        "extending past NP1 into D's region must set OVERLAPS_OTHER_SEGMENT"
    );
}

#[test]
fn v_overlap_into_d_bumps_v_live_call_version() {
    // Plumbing check: AssembleSegment(D) must trigger a V
    // refresh under the cross-segment hook.
    let (cfg, v01, _v02, d01) = v_d_overlap_refdata();
    let plan = v_overlap_plan(&cfg, v01, d01, 3, 0);

    let compiled = CompiledSimulator::compile(&plan, Some(&cfg), None, ExecutionPolicy::Permissive)
        .expect("plan should compile");
    let outcome = compiled.run_one(0).expect("plan should run");

    let pass_names = &outcome.pass_names;
    let assemble_v_idx = pass_names
        .iter()
        .position(|n| n == "assemble.v")
        .expect("plan must include assemble.v");
    let assemble_d_idx = pass_names
        .iter()
        .position(|n| n == "assemble.d")
        .expect("plan must include assemble.d");
    let post_assemble_v_version = outcome.revisions[assemble_v_idx + 1].segment_calls.version;
    let post_assemble_d_version = outcome.revisions[assemble_d_idx + 1].segment_calls.version;
    assert!(
        post_assemble_d_version > post_assemble_v_version,
        "AssembleSegment(D) must bump the version (refreshing V) for \
             overlap detection: {post_assemble_d_version} should be > \
             {post_assemble_v_version}"
    );
}
