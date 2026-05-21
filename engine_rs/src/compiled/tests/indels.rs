use super::*;

// ──────────────────────────────────────────────────────────────
// Structural indel local recomputation.
//
// Indels (insertions / deletions) shift the pool layout under
// V/D/J regions. The walker now tolerates:
//   - synthetic (indel-inserted) nucleotides inside V/D/J
//     regions (skip without failing),
//   - forward jumps in `germline_pos` caused by deletions
//     between observed positions.
// The `apply_live_call_updates` hook on `PassEffect::StructuralIndel`
// refreshes V/D/J live calls so post-indel evidence drives the
// call set.
// ──────────────────────────────────────────────────────────────

/// Test-only deterministic Pass that deletes one nucleotide at a
/// fixed pool position and reports `PassEffect::StructuralIndel`.
/// Mirrors the production `IndelPass`'s deletion path but without
/// RNG involvement.
#[derive(Clone, Debug)]
pub(super) struct DeleteAtPass {
    at: u32,
}

impl DeleteAtPass {
    pub(super) fn new(at: u32) -> Self {
        Self { at }
    }
}

impl Pass for DeleteAtPass {
    fn name(&self) -> &str {
        "test.delete_at"
    }
    fn execute(&self, sim: &Simulation, _ctx: &mut crate::pass::PassContext) -> Simulation {
        sim.with_indel_deleted(self.at)
    }
    fn effects(&self) -> Vec<PassEffect> {
        vec![PassEffect::StructuralIndel]
    }
}

/// Test-only deterministic Pass that inserts one nucleotide at a
/// fixed pool position with the chosen segment / base, and reports
/// `PassEffect::StructuralIndel`. Mirrors `IndelPass`'s insertion
/// path but without RNG involvement.
#[derive(Clone, Debug)]
pub(super) struct InsertAtPass {
    at: u32,
    base: u8,
    segment: Segment,
}

impl InsertAtPass {
    pub(super) fn new(at: u32, base: u8, segment: Segment) -> Self {
        Self { at, base, segment }
    }
}

impl Pass for InsertAtPass {
    fn name(&self) -> &str {
        "test.insert_at"
    }
    fn execute(&self, sim: &Simulation, _ctx: &mut crate::pass::PassContext) -> Simulation {
        let nuc = crate::ir::Nucleotide::synthetic(
            self.base,
            self.segment,
            crate::ir::flag::INDEL_INSERTED,
        );
        sim.with_indel_inserted(self.at, nuc)
    }
    fn effects(&self) -> Vec<PassEffect> {
        vec![PassEffect::StructuralIndel]
    }
}

/// V-only refdata holding alleles where the differing position is
/// known and predictable.
fn v_refdata(seqs: &[(&str, &[u8])]) -> RefDataConfig {
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

/// Build a sample(V) → assemble(V) → indel plan and run it. Returns
/// the final V live call.
fn run_indel_plan(
    cfg: &RefDataConfig,
    sampled: AlleleId,
    indels: Vec<Box<dyn Pass>>,
) -> (crate::pass::Outcome, crate::live_call::SegmentLiveCall) {
    let mut plan = PassPlan::new();
    plan.push(Box::new(SampleAllelePass::new(
        Segment::V,
        Box::new(AllelePoolDist::restricted_uniform(
            &cfg.v_pool,
            vec![sampled],
        )),
    )));
    plan.push(Box::new(AssembleSegmentPass::new(Segment::V)));
    for indel in indels {
        plan.push(indel);
    }
    let compiled = CompiledSimulator::compile(&plan, Some(cfg), None, ExecutionPolicy::Permissive)
        .expect("fixture plan should compile");
    let outcome = compiled.run_one(0).expect("plan should run");
    let v_call = outcome
        .final_simulation()
        .live_calls
        .as_ref()
        .expect("live calls populated after assembly")
        .get(Segment::V)
        .cloned()
        .expect("V live call exists after assembly");
    (outcome, v_call)
}

#[test]
fn deletion_inside_v_widens_live_call_when_distinguishing_base_removed() {
    // V*01 = AAAGAAAA (position 3 = G distinguishes it).
    // V*02 = AAACAAAA (position 3 = C distinguishes it).
    // Sample V*01 → assembled "AAAGAAAA" → live call {V*01} (the
    // G at pos 3 disambiguates).
    // Delete position 3 → assembled becomes "AAAAAAA" (7 bases),
    // covering germline positions 0,1,2,4,5,6,7. With pos 3 absent
    // the alleles are pairwise indistinguishable at the remaining
    // positions → live call widens to {V*01, V*02}.
    let cfg = v_refdata(&[("V*01", b"AAAGAAAA"), ("V*02", b"AAACAAAA")]);
    let v01 = AlleleId::new(0);
    let v02 = AlleleId::new(1);

    let (_outcome, v_call) = run_indel_plan(&cfg, v01, vec![Box::new(DeleteAtPass::new(3))]);
    let mut ids = v_call.allele_call.to_ids();
    ids.sort_by_key(|id| id.index());
    assert_eq!(ids, vec![v01, v02]);
}

#[test]
fn deletion_preserves_call_when_distinguishing_base_remains() {
    // Same fixture as above, but delete a SHARED base (position 0,
    // both alleles have 'A'). The distinguishing position 3 is
    // still there → live call stays singleton {V*01}.
    let cfg = v_refdata(&[("V*01", b"AAAGAAAA"), ("V*02", b"AAACAAAA")]);
    let v01 = AlleleId::new(0);

    let (_outcome, v_call) = run_indel_plan(&cfg, v01, vec![Box::new(DeleteAtPass::new(0))]);
    assert_eq!(v_call.allele_call.to_ids(), vec![v01]);
}

#[test]
fn insertion_inside_v_does_not_fail_the_call() {
    // Sample V*01 → assembled "AAAGAAAA" → live call {V*01}.
    // Insert a synthetic 'C' at pool position 3 (inside V's
    // region). The walker should SKIP the inserted nucleotide
    // (it has GermlinePos::NONE) without failing the call. The
    // remaining positions still uniquely identify V*01.
    let cfg = v_refdata(&[("V*01", b"AAAGAAAA"), ("V*02", b"AAACAAAA")]);
    let v01 = AlleleId::new(0);

    let (_outcome, v_call) = run_indel_plan(
        &cfg,
        v01,
        vec![Box::new(InsertAtPass::new(3, b'C', Segment::V))],
    );
    assert_eq!(
        v_call.allele_call.to_ids(),
        vec![v01],
        "insertion should not collapse the live call to Unsupported"
    );
    assert!(
        !matches!(
            v_call.confidence,
            crate::live_call::LiveCallConfidence::Unsupported
        ),
        "live call must remain supported after a single in-region insertion"
    );
}

#[test]
fn combined_insertion_and_deletion_recomputes_correctly() {
    // Sample V*01 → live call {V*01} via the distinguishing G at
    // position 3. Then:
    //   1. Insert a synthetic 'C' at pool position 5 (a benign
    //      mid-region insertion the walker skips).
    //   2. Delete the original distinguishing base at germline
    //      position 3 (which is now at pool position 3 still —
    //      the insertion at 5 didn't shift pos 0..4).
    // Expected post-state: live call widens to {V*01, V*02}
    // because the distinguishing germline position 3 is gone.
    let cfg = v_refdata(&[("V*01", b"AAAGAAAA"), ("V*02", b"AAACAAAA")]);
    let v01 = AlleleId::new(0);
    let v02 = AlleleId::new(1);

    let (_outcome, v_call) = run_indel_plan(
        &cfg,
        v01,
        vec![
            Box::new(InsertAtPass::new(5, b'C', Segment::V)),
            Box::new(DeleteAtPass::new(3)),
        ],
    );
    let mut ids = v_call.allele_call.to_ids();
    ids.sort_by_key(|id| id.index());
    assert_eq!(ids, vec![v01, v02]);
}

#[test]
fn structural_indel_bumps_live_call_version() {
    // Plumbing check: every refresh fired by `PassEffect::StructuralIndel`
    // bumps `evidence_version`. Catches accidental removal of the
    // hook from `apply_live_call_updates`.
    let cfg = v_refdata(&[("V*01", b"AAAGAAAA"), ("V*02", b"AAACAAAA")]);
    let v01 = AlleleId::new(0);

    let (outcome, _final_call) = run_indel_plan(&cfg, v01, vec![Box::new(DeleteAtPass::new(3))]);
    let pass_names = &outcome.pass_names;
    let assemble_idx = pass_names
        .iter()
        .position(|n| n == "assemble.v")
        .expect("plan must include assemble.v");
    let indel_idx = pass_names
        .iter()
        .position(|n| n == "test.delete_at")
        .expect("plan must include test.delete_at");
    let post_assemble_version = outcome.revisions[assemble_idx + 1]
        .live_calls
        .as_ref()
        .expect("post-assemble live calls present")
        .version;
    let post_indel_version = outcome.revisions[indel_idx + 1]
        .live_calls
        .as_ref()
        .expect("post-indel live calls present")
        .version;
    assert!(
        post_indel_version > post_assemble_version,
        "StructuralIndel must bump live-call version: \
             {post_indel_version} should be > {post_assemble_version}"
    );
}
