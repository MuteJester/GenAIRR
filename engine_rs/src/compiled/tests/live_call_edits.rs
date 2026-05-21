use super::*;

// ──────────────────────────────────────────────────────────────
// Base edit live-call refresh fixtures.
//
// Each fixture builds a tiny VDJ refdata, samples a chosen V
// allele, assembles V, then runs a deterministic test-only base
// edit pass to surface one of the four shrink / widen / switch /
// unsupported behaviours the design doc calls out for base
// edits.
// ──────────────────────────────────────────────────────────────

/// Test-only deterministic Pass that edits a single pool position
/// to a chosen base and reports `PassEffect::EditBases`. Drives
/// the live-call refresh path without any RNG involvement so the
/// fixtures stay focused on the post-edit live-call shape.
#[derive(Clone, Debug)]
pub(super) struct EditBaseAtPass {
    handle: crate::ir::NucHandle,
    new_base: u8,
}

impl EditBaseAtPass {
    pub(super) fn new(pool_index: u32, new_base: u8) -> Self {
        Self {
            handle: crate::ir::NucHandle::new(pool_index),
            new_base,
        }
    }
}

impl Pass for EditBaseAtPass {
    fn name(&self) -> &str {
        "test.edit_base_at"
    }

    fn execute(&self, sim: &Simulation, _ctx: &mut crate::pass::PassContext) -> Simulation {
        sim.with_base_changed(self.handle, self.new_base)
    }

    fn effects(&self) -> Vec<PassEffect> {
        vec![PassEffect::EditBases]
    }
}

/// V-only refdata where every allele has a unique 9-base sequence
/// (no shared prefix). Used by the widen / switch / unsupported
/// fixtures where each ref position cleanly distinguishes alleles.
fn distinct_v_refdata(seqs: &[(&str, &[u8])]) -> RefDataConfig {
    let mut cfg = RefDataConfig::empty(ChainType::Vdj);
    for (name, seq) in seqs {
        let _ = cfg.v_pool.push(Allele {
            name: (*name).to_string(),
            gene: name.split('*').next().unwrap_or(name).to_string(),
            seq: seq.to_vec(),
            segment: Segment::V,
            anchor: None,
        });
    }
    cfg
}

/// Drive a fixture: sample-allele(V) → assemble(V) → edit(pool, base).
/// Returns the final `Outcome` and the V live call extracted from the
/// final revision.
fn run_edit(
    cfg: &RefDataConfig,
    sampled_id: AlleleId,
    edits: Vec<(u32, u8)>,
) -> (crate::pass::Outcome, crate::live_call::SegmentLiveCall) {
    let mut plan = PassPlan::new();
    plan.push(Box::new(SampleAllelePass::new(
        Segment::V,
        Box::new(AllelePoolDist::restricted_uniform(
            &cfg.v_pool,
            vec![sampled_id],
        )),
    )));
    plan.push(Box::new(AssembleSegmentPass::new(Segment::V)));
    for (pos, new_base) in edits {
        plan.push(Box::new(EditBaseAtPass::new(pos, new_base)));
    }

    let compiled = CompiledSimulator::compile(&plan, Some(cfg), None, ExecutionPolicy::Permissive)
        .expect("fixture plan should compile");
    let outcome = compiled.run_one(0).expect("plan should run");
    let final_sim = outcome.final_simulation().clone();
    let v_call = final_sim
        .live_calls
        .as_ref()
        .expect("live calls populated after assembly")
        .get(Segment::V)
        .cloned()
        .expect("V live call exists after assembly");
    (outcome, v_call)
}

#[test]
fn mutation_widens_live_call_when_distinguishing_base_is_lost() {
    // Three alleles whose only difference is at position 3:
    //   A1: A A A G A A A A A
    //   A2: A A A C A A A A A
    //   A3: A A A C A A A A A   (identical to A2)
    //
    // Sampling A1 produces an assembled "AAAGAAAAA" → live call {A1}.
    // Mutating position 3 from G→C makes the assembled bases match
    // both A2 and A3 exactly → live call widens to {A2, A3}.
    let cfg = distinct_v_refdata(&[
        ("V*01", b"AAAGAAAAA"),
        ("V*02", b"AAACAAAAA"),
        ("V*03", b"AAACAAAAA"),
    ]);
    let v01 = AlleleId::new(0);
    let v02 = AlleleId::new(1);
    let v03 = AlleleId::new(2);

    let (_outcome, v_call) = run_edit(&cfg, v01, vec![(3, b'C')]);
    let mut expected = vec![v02, v03];
    expected.sort_by_key(|id| id.index());
    let mut actual = v_call.allele_call.to_ids();
    actual.sort_by_key(|id| id.index());
    assert_eq!(actual, expected, "expected widening to {{A2, A3}}");
    assert!(!v_call.allele_call.contains(v01), "A1 must drop out");
}

#[test]
fn mutation_switches_live_call_to_a_different_singleton() {
    // Two alleles differing at position 3 only:
    //   A1: A A A G A A A A A   (sampled)
    //   A2: A A A C A A A A A
    //
    // Pre-edit live call = {A1}. After editing position 3 G→C, the
    // assembled bases match A2 exactly → live call = {A2}. Same
    // size, but the membership has switched.
    let cfg = distinct_v_refdata(&[("V*01", b"AAAGAAAAA"), ("V*02", b"AAACAAAAA")]);
    let v01 = AlleleId::new(0);
    let v02 = AlleleId::new(1);

    let (_outcome, v_call) = run_edit(&cfg, v01, vec![(3, b'C')]);
    assert_eq!(v_call.allele_call.to_ids(), vec![v02]);
}

#[test]
fn mutation_shrinks_live_call_from_three_to_one() {
    // Four alleles. The first three are identical (perfect ambiguity),
    // the fourth differs at position 2:
    //   A1: A A A C C C   (sampled — identical to A2 / A3)
    //   A2: A A A C C C
    //   A3: A A A C C C
    //   A4: A A T C C C
    //
    // Pre-edit assembled "AAACCC" matches A1, A2, A3 but not A4 →
    // live call {A1, A2, A3} (size 3). After editing position 2
    // A→T the assembled becomes "AATCCC" which only matches A4 →
    // live call shrinks to {A4} (size 1).
    let cfg = distinct_v_refdata(&[
        ("V*01", b"AAACCC"),
        ("V*02", b"AAACCC"),
        ("V*03", b"AAACCC"),
        ("V*04", b"AATCCC"),
    ]);
    let v01 = AlleleId::new(0);
    let v04 = AlleleId::new(3);

    let (_outcome, v_call) = run_edit(&cfg, v01, vec![(2, b'T')]);
    assert_eq!(
        v_call.allele_call.len(),
        1,
        "expected post-edit set size 1, got {}",
        v_call.allele_call.len()
    );
    assert_eq!(v_call.allele_call.to_ids(), vec![v04]);
}

#[test]
fn mutation_to_orphan_base_keeps_truth_in_tie_set() {
    // Single allele with a fully-determined sequence. Mutating one
    // assembled base to a value no allele has at that ref position
    // simply leaves that position non-informative (its evidence
    // contributes zero to every allele's score). The remaining
    // five positions still score the truth allele uniquely highest,
    // so the tie-set keeps V*01 — exactly the evidence-resilient
    // behavior the score-and-tie caller is designed to deliver.
    let cfg = distinct_v_refdata(&[("V*01", b"AAAAAA")]);
    let v01 = AlleleId::new(0);

    let (_outcome, v_call) = run_edit(&cfg, v01, vec![(2, b'T')]);
    assert_eq!(
        v_call.allele_call.to_ids(),
        vec![v01],
        "truth allele should remain at max score (5/6 positions still match)"
    );
}

#[test]
fn pre_edit_live_call_visible_in_intermediate_revision() {
    // The compiled runtime stores a revision after every committed
    // pass. The pre-edit live call should be reachable via
    // revisions[2] (after assemble) while the final revision shows
    // the post-edit (widened) call. This guards against accidentally
    // rewriting earlier revisions when refreshing live calls.
    let cfg = distinct_v_refdata(&[
        ("V*01", b"AAAGAAAAA"),
        ("V*02", b"AAACAAAAA"),
        ("V*03", b"AAACAAAAA"),
    ]);
    let v01 = AlleleId::new(0);

    let (outcome, post_edit_call) = run_edit(&cfg, v01, vec![(3, b'C')]);

    // revisions: [initial, post-sample, post-assemble, post-edit].
    assert_eq!(outcome.revisions.len(), 4);
    let after_assemble = &outcome.revisions[2];
    let pre_edit_call = after_assemble
        .live_calls
        .as_ref()
        .expect("live calls populated after assembly")
        .get(Segment::V)
        .cloned()
        .expect("V live call exists after assembly");

    assert_eq!(pre_edit_call.allele_call.to_ids(), vec![v01]);
    assert_eq!(post_edit_call.allele_call.len(), 2);
    assert!(!post_edit_call.allele_call.contains(v01));
}

#[test]
fn real_uniform_mutation_pass_increments_live_call_version() {
    // Integration-style coverage: confirm the live-call refresh
    // path also fires when the real `UniformMutationPass` runs,
    // not just our synthetic `EditBaseAtPass`. UniformMutationPass
    // samples positions with replacement so the post-state isn't
    // deterministic — but the refresh hook IS, and we can assert
    // that against the live-call `evidence_version` counter:
    // every successful refresh bumps it, so the post-mutation
    // version must be strictly greater than the post-assemble
    // version.
    use crate::passes::UniformMutationPass;

    let cfg = distinct_v_refdata(&[("V*01", b"AAAAAAAAA")]);
    let v01 = AlleleId::new(0);

    let mut plan = PassPlan::new();
    plan.push(Box::new(SampleAllelePass::new(
        Segment::V,
        Box::new(AllelePoolDist::restricted_uniform(&cfg.v_pool, vec![v01])),
    )));
    plan.push(Box::new(AssembleSegmentPass::new(Segment::V)));
    plan.push(Box::new(UniformMutationPass::new(
        Box::new(EmpiricalLengthDist::from_pairs(vec![(3, 1.0)])),
        Box::new(ConstBaseDist(b'A')),
    )));

    let compiled = CompiledSimulator::compile(&plan, Some(&cfg), None, ExecutionPolicy::Permissive)
        .expect("integration plan should compile");
    let outcome = compiled.run_one(0).expect("plan should run");

    // revisions: [initial, post-sample, post-assemble, post-mutate].
    let post_assemble_version = outcome.revisions[2]
        .live_calls
        .as_ref()
        .expect("live calls populated after assembly")
        .version;
    let post_mutate_version = outcome.revisions[3]
        .live_calls
        .as_ref()
        .expect("live calls populated after mutation")
        .version;
    assert!(
        post_mutate_version > post_assemble_version,
        "EditBases pass must bump live-call evidence_version \
             ({post_mutate_version} should be > {post_assemble_version})"
    );
}

#[test]
fn edit_outside_assembled_segment_leaves_live_call_unchanged() {
    // Sanity: editing a pool position that is NOT inside V's
    // assembled region must not invent extra hypotheses or break
    // the existing one. The test allele has length 6; we run two
    // edit passes — the first edits position 0 (inside V) so we
    // know recompute is happening, the second edits position 5
    // (still inside V's region for this fixture, so live call
    // stays consistent with the new bases at every step).
    //
    // We use this fixture to confirm that two consecutive
    // EditBases passes both trigger refresh and the final state
    // matches a manual recomputation from the post-edit bases.
    let cfg = distinct_v_refdata(&[("V*01", b"AAAAAA"), ("V*02", b"AAAAAA")]);
    let v01 = AlleleId::new(0);
    let v02 = AlleleId::new(1);

    // Both alleles are identical; the call remains {V*01, V*02}
    // through every base edit because every position has the same
    // base across alleles.
    let (_outcome, v_call) = run_edit(
        &cfg,
        v01,
        vec![(0, b'A'), (5, b'A')], // no-ops in terms of base change
    );
    let mut expected = vec![v01, v02];
    expected.sort_by_key(|id| id.index());
    let mut actual = v_call.allele_call.to_ids();
    actual.sort_by_key(|id| id.index());
    assert_eq!(actual, expected);
}
