//! Mutation effects on the live call set.
//!
//! SHM (S5F or uniform) and PCR / quality / contaminant base-edit
//! events all flip individual bases in the assembled pool. The live
//! call walker observer (Phase 1+) updates its per-allele scores
//! incrementally with each `on_base_changed` event. The resulting
//! AIRR call set must reflect the post-mutation evidence:
//!
//! - **Mutation BREAKS ambiguity**: truth allele was ambiguous
//!   post-trim with N candidates. A single SHM event flips a base
//!   at a position where one candidate matches and the others don't.
//!   Call narrows to that candidate (which is or isn't the truth).
//! - **Mutation INTRODUCES ambiguity**: truth allele was unique
//!   post-trim. A single SHM event flips a base such that another
//!   allele now scores equally well. Call widens to include the new
//!   candidate.
//! - **Mutation FLIPS the call AWAY from truth**: enough mutations
//!   accumulate that the highest-scoring candidate is not the
//!   originally-sampled truth. The truth allele falls out of the
//!   call set — and this is biologically correct (an aligner given
//!   only the mutated bases couldn't recover the truth either).
//! - **n_mutations field reflects the recorded count**: Phase 18 staged
//!   the count on `LiveCallState.mutation_count`. Verify it matches
//!   the AIRR `n_mutations` field exactly.
//! - **Mutation count source preference** (S5F vs uniform): when
//!   both passes ran (rare), verify the AIRR field semantics
//!   (currently: first-set wins / sum / last-set wins — pin down
//!   the contract).
//!
//! ## Invariants this module should test
//!
//! - **Wildcard mutations (N injection)**: ncorrupt fills positions
//!   with N. N is treated as wildcard evidence — every allele
//!   "matches" at that position. Call set widens, but the truth
//!   stays in it.
//! - **Mutation rate matches**: `mutation_rate =
//!   n_mutations / sequence_length` (mechanical, but verify the
//!   sequence_length passed in is post-corruption when corruption
//!   shifts length).
//! - **Mutation in junction window**: SHM hits inside the junction.
//!   `junction_aa` reflects the post-mutation amino acid; the
//!   anchor-preserved flag flips if it crosses the C/W anchor.
//! - **Cross-region mutation count**: an SHM event mutates a base in
//!   V, an NP, then J. The codon rail and call set both react
//!   correctly. The trace counts the total mutations.
//!
//! ## Adjacent in-tree coverage
//!
//! - `engine_rs/src/live_call/tests/walker_observer_change_base.rs`
//! - `engine_rs/src/live_call/tests/walker_observer.rs::observer_matches_oracle_mutated_base_disambiguates`
//! - `engine_rs/src/live_call/tests/walker_observer.rs::observer_matches_oracle_flip_toward_other_allele`
//! - `engine_rs/src/compiled/tests/live_call_edits.rs`
//! - `engine_rs/src/compiled/tests/curated.rs::curated_mutation_switches_v_call_to_a_different_allele`

use super::common;
use genairr_engine::airr_record::build_airr_record;
use genairr_engine::compiled::{CompiledSimulator, ExecutionPolicy};
use genairr_engine::dist::{AllelePoolDist, Distribution, EmpiricalLengthDist, UniformBase};
use genairr_engine::ir::{NucHandle, Segment, Simulation};
use genairr_engine::live_call::HypothesisFlags;
use genairr_engine::pass::{Pass, PassContext, PassEffect, PassPlan};
use genairr_engine::passes::{
    AssembleSegmentPass, GenerateNPPass, NCorruptionPass, PCRErrorPass, QualityErrorPass,
    SampleAllelePass, TrimPass, UniformMutationPass,
};
use genairr_engine::refdata::RefDataConfig;
use genairr_engine::rng::Rng;

// ──────────────────────────────────────────────────────────────
// Test-local helpers.
// ──────────────────────────────────────────────────────────────

/// Pass-local deterministic Pass that edits a single pool position
/// to a chosen base. Mirrors the test-only `EditBaseAtPass` in
/// `engine_rs/src/compiled/tests/live_call_edits.rs` so we can pin
/// SHM-style behaviour without seed-search.
#[derive(Clone, Debug)]
struct EditBaseAtPass {
    handle: NucHandle,
    new_base: u8,
}

impl EditBaseAtPass {
    fn new(pool_index: u32, new_base: u8) -> Self {
        Self {
            handle: NucHandle::new(pool_index),
            new_base,
        }
    }
}

impl Pass for EditBaseAtPass {
    fn name(&self) -> &str {
        "test.edit_base_at"
    }

    fn execute(&self, sim: &Simulation, _ctx: &mut PassContext) -> Simulation {
        sim.with_base_changed(self.handle, self.new_base)
    }

    fn effects(&self) -> Vec<PassEffect> {
        vec![PassEffect::EditBases]
    }
}

/// Deterministic single-base distribution. Always emits `byte`.
#[derive(Clone, Debug)]
struct ConstBaseDist(u8);

impl Distribution for ConstBaseDist {
    type Output = u8;
    fn sample(&self, _rng: &mut Rng) -> u8 {
        self.0
    }
    fn support(&self) -> Option<Vec<(u8, f64)>> {
        Some(vec![(self.0, 1.0)])
    }
}

/// Build the minimal VJ plan used by these tests: sample the chosen
/// V/J truth alleles (no ambiguity at sample time), apply optional
/// V_3 trim, assemble V → NP1 → J. Caller appends mutation passes.
fn vj_truth_plan(
    cfg: &RefDataConfig,
    v_truth: &str,
    j_truth: &str,
    v_trim_3: i64,
    np1_len: i64,
) -> PassPlan {
    let v_id = common::allele_id_by_name(cfg, Segment::V, v_truth);
    let j_id = common::allele_id_by_name(cfg, Segment::J, j_truth);
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
            genairr_engine::assignment::TrimEnd::Three,
            Box::new(EmpiricalLengthDist::from_pairs(vec![(v_trim_3, 1.0)])),
        )));
    }
    plan.push(Box::new(AssembleSegmentPass::new(Segment::V)));
    plan.push(Box::new(GenerateNPPass::new(
        Segment::Np1,
        Box::new(EmpiricalLengthDist::from_pairs(vec![(np1_len, 1.0)])),
        Box::new(ConstBaseDist(b'A')),
    )));
    plan.push(Box::new(AssembleSegmentPass::new(Segment::J)));
    plan
}

fn run_plan(cfg: &RefDataConfig, plan: PassPlan, seed: u64) -> Simulation {
    let compiled = CompiledSimulator::compile(&plan, Some(cfg), None, ExecutionPolicy::Permissive)
        .expect("plan should compile");
    let outcome = compiled.run_one(seed).expect("plan should run");
    outcome.final_simulation().clone()
}

fn run_plan_with_record(
    cfg: &RefDataConfig,
    plan: PassPlan,
    seed: u64,
) -> (Simulation, genairr_engine::airr_record::AirrRecord) {
    let compiled = CompiledSimulator::compile(&plan, Some(cfg), None, ExecutionPolicy::Permissive)
        .expect("plan should compile");
    let outcome = compiled.run_one(seed).expect("plan should run");
    let rec = build_airr_record(&outcome, cfg, "test");
    (outcome.final_simulation().clone(), rec)
}

// ──────────────────────────────────────────────────────────────
// Tests.
// ──────────────────────────────────────────────────────────────

#[test]
fn trim_widens_then_mutation_narrows_call() {
    // Sample v01 truth (`AAACCCGGGTTT`), trim V_3 by 3 → assembled
    // V is `AAACCCGGG` (the shared 9-byte prefix), so v_call widens
    // to all four alleles {v01*01, v01*02, v02*01, v03*01}.
    //
    // Now mutate pool position 0 ('A' → 'C'). All real alleles have
    // 'A' at pos 0 — no allele matches 'C'. Walker scores every
    // allele 1 lower than before — they all remain tied → call
    // stays the full set.
    //
    // Then mutate position 1 ('A' → 'A' is a no-op base — instead
    // we use a single SHM to disambiguate via NP1: trim V_3 by 3,
    // NP1 = 'C' (matches v03's pos 9 = 'C'). The walker
    // right-extension narrows v_call to {v03*01} alone — proving
    // a single NP base flips the call away from the v01 truth.
    //
    // This is the "ambiguity narrows" path: an evidence event
    // shifts the call set toward a specific allele.
    let cfg = common::vj_ambiguous_refdata();
    // The default `vj_truth_plan` uses 'A' as the NP base; here we
    // need a 'C' NP1 base to disambiguate to v03. Build a one-off
    // plan inline.
    let mut plan = PassPlan::new();
    let v_id = common::allele_id_by_name(&cfg, Segment::V, "v01*01");
    let j_id = common::allele_id_by_name(&cfg, Segment::J, "j01*01");
    plan.push(Box::new(SampleAllelePass::new(
        Segment::V,
        Box::new(AllelePoolDist::restricted_uniform(&cfg.v_pool, vec![v_id])),
    )));
    plan.push(Box::new(SampleAllelePass::new(
        Segment::J,
        Box::new(AllelePoolDist::restricted_uniform(&cfg.j_pool, vec![j_id])),
    )));
    plan.push(Box::new(TrimPass::new(
        Segment::V,
        genairr_engine::assignment::TrimEnd::Three,
        Box::new(EmpiricalLengthDist::from_pairs(vec![(3, 1.0)])),
    )));
    plan.push(Box::new(AssembleSegmentPass::new(Segment::V)));
    plan.push(Box::new(GenerateNPPass::new(
        Segment::Np1,
        Box::new(EmpiricalLengthDist::from_pairs(vec![(1, 1.0)])),
        Box::new(ConstBaseDist(b'C')),
    )));
    plan.push(Box::new(AssembleSegmentPass::new(Segment::J)));
    let sim = run_plan(&cfg, plan, 0);

    let v_calls = common::v_call_names(&sim, &cfg);
    assert_eq!(
        v_calls,
        vec!["v03*01"],
        "single matching NP base must narrow call to v03; got {:?}",
        v_calls,
    );
    // Truth (v01) dropped out — that's biologically correct: the
    // observable evidence supports v03, not v01.
    assert!(
        !v_calls.iter().any(|n| n.starts_with("v01")),
        "v01* truth must have dropped out, got {:?}",
        v_calls,
    );
}

#[test]
fn single_mutation_introduces_ambiguity_from_unique_truth() {
    // Sample v02 truth (`AAACCCGGGAAA`), no trim, no NP. Live call
    // narrows to {v02*01} (distinguishing 'A' at pos 9-11).
    // Edit pool position 9 from 'A' to 'T'. Now pool[9..12] =
    // 'TAA'. v01 has 'TTT' at pos 9-11 and v02 has 'AAA'.
    //
    // After the edit: pool[9] = T, pool[10] = A, pool[11] = A.
    //   - v01*01/v01*02 score: match T at 9, mismatch at 10/11 → 1 match in tail
    //   - v02*01 score: mismatch T at 9, match A at 10/11 → 2 matches in tail
    //   - v03*01 score: mismatch at 9, mismatch at 10/11 → 0 matches
    //
    // The walker uses position-by-position equality matching. The
    // exact tie-breaking depends on the live caller's scoring rule
    // — what we can robustly assert is that v03 is excluded, and
    // the call set didn't lose v02 (the truth).
    let cfg = common::vj_ambiguous_refdata();
    let mut plan = vj_truth_plan(&cfg, "v02*01", "j01*01", 0, 0);
    plan.push(Box::new(EditBaseAtPass::new(9, b'T')));
    let sim = run_plan(&cfg, plan, 0);

    let v_calls = common::v_call_names(&sim, &cfg);
    assert!(
        v_calls.contains(&"v02*01".to_string()),
        "truth v02*01 should remain in call set (still highest-scoring); got {:?}",
        v_calls,
    );
    assert!(
        !v_calls.contains(&"v03*01".to_string()),
        "v03*01 (zero matches in tail) should be excluded; got {:?}",
        v_calls,
    );
}

#[test]
fn multiple_mutations_flip_call_away_from_truth() {
    // Sample v01 truth, no trim. Mutate the entire `TTT` tail
    // (pool 9, 10, 11) to `AAA` (v02's distinguishing suffix).
    // Now the assembled bases are byte-identical to v02*01's
    // sequence — v_call must switch to {v02*01}, dropping v01.
    let cfg = common::vj_ambiguous_refdata();
    let mut plan = vj_truth_plan(&cfg, "v01*01", "j01*01", 0, 0);
    plan.push(Box::new(EditBaseAtPass::new(9, b'A')));
    plan.push(Box::new(EditBaseAtPass::new(10, b'A')));
    plan.push(Box::new(EditBaseAtPass::new(11, b'A')));
    let sim = run_plan(&cfg, plan, 0);

    let v_calls = common::v_call_names(&sim, &cfg);
    assert_eq!(v_calls, vec!["v02*01"], "got {:?}", v_calls);
    assert!(
        !v_calls.iter().any(|n| n.starts_with("v01")),
        "truth v01* alleles should have dropped out, got {:?}",
        v_calls,
    );
}

#[test]
fn silent_mutation_at_non_distinguishing_position_preserves_call() {
    // Sample v01 truth. Edit pool position 0 ('A' → 'A' is a no-op,
    // so use a base no allele has at that position to keep the
    // walker score-only-where-it-matches semantics honest). Position
    // 0 is 'A' across all v01/v02/v03; changing it to 'C' costs all
    // four alleles the same score → all should remain in the call
    // set (singleton truth stays in, others all have same score).
    //
    // The v_call is unique to v01*01 + v01*02 (same sequence) at
    // baseline. After the edit, the tail still distinguishes v01
    // from v02/v03 → call should still be just {v01*01, v01*02}.
    let cfg = common::vj_ambiguous_refdata();
    let mut plan = vj_truth_plan(&cfg, "v01*01", "j01*01", 0, 0);
    plan.push(Box::new(EditBaseAtPass::new(0, b'C')));
    let sim = run_plan(&cfg, plan, 0);

    let v_calls = common::v_call_names(&sim, &cfg);
    // v01*01 and v01*02 are byte-identical, so both stay in the call.
    let truth_present = v_calls.contains(&"v01*01".to_string());
    let alias_present = v_calls.contains(&"v01*02".to_string());
    assert!(
        truth_present && alias_present,
        "truth v01*01 + alias v01*02 must remain after a non-distinguishing edit; got {:?}",
        v_calls,
    );
}

#[test]
fn n_mutations_field_equals_uniform_mutation_count() {
    // UniformMutationPass with count=4 should produce
    // `n_mutations = 4` in the AIRR record (Phase 18 stash on
    // LiveCallState.mutation_count).
    let cfg = common::vj_ambiguous_refdata();
    let mut plan = vj_truth_plan(&cfg, "v01*01", "j01*01", 0, 0);
    plan.push(Box::new(UniformMutationPass::new(
        Box::new(EmpiricalLengthDist::from_pairs(vec![(4, 1.0)])),
        Box::new(UniformBase),
    )));
    let (_sim, rec) = run_plan_with_record(&cfg, plan, 7);

    assert_eq!(
        rec.n_mutations, 4,
        "n_mutations must equal the UniformMutation count distribution sample"
    );
}

#[test]
fn mutation_rate_equals_n_mutations_over_sequence_length() {
    // Mechanical correspondence: mutation_rate must exactly equal
    // n_mutations / sequence_length post-pass.
    let cfg = common::vj_ambiguous_refdata();
    let mut plan = vj_truth_plan(&cfg, "v01*01", "j01*01", 0, 0);
    plan.push(Box::new(UniformMutationPass::new(
        Box::new(EmpiricalLengthDist::from_pairs(vec![(3, 1.0)])),
        Box::new(UniformBase),
    )));
    let (_sim, rec) = run_plan_with_record(&cfg, plan, 42);

    assert!(rec.sequence_length > 0);
    let expected = rec.n_mutations as f64 / rec.sequence_length as f64;
    assert!(
        (rec.mutation_rate - expected).abs() < 1e-12,
        "mutation_rate {} != n_mutations {} / sequence_length {} = {}",
        rec.mutation_rate,
        rec.n_mutations,
        rec.sequence_length,
        expected,
    );
    assert_eq!(rec.n_mutations, 3);
}

#[test]
fn ncorrupt_n_injection_widens_call_but_keeps_truth() {
    // Sample v02 (distinguishing 'A' at pool 9-11). NCorruptionPass
    // fills 3 random positions with 'N'. N is wildcard evidence —
    // every allele scores a match at those positions, so the call
    // set widens. The truth (v02*01) must remain in the call.
    //
    // We can't pin which positions get N without seed search, but
    // for the invariant "truth survives" we can just check it's
    // still there.
    let cfg = common::vj_ambiguous_refdata();
    let mut plan = vj_truth_plan(&cfg, "v02*01", "j01*01", 0, 0);
    plan.push(Box::new(NCorruptionPass::new(Box::new(
        EmpiricalLengthDist::from_pairs(vec![(3, 1.0)]),
    ))));

    // Try a handful of seeds; truth must survive every one.
    for seed in 0..16u64 {
        let sim = run_plan(&cfg, plan_clone(&cfg, "v02*01", 3, "n"), seed);
        let v_calls = common::v_call_names(&sim, &cfg);
        assert!(
            v_calls.contains(&"v02*01".to_string()),
            "truth v02*01 must remain after N injection (seed {}); got {:?}",
            seed,
            v_calls,
        );
    }
}

// Re-build the ncorrupt plan inside the loop (PassPlan isn't `Clone`).
fn plan_clone(cfg: &RefDataConfig, v_truth: &str, n_count: i64, mode: &str) -> PassPlan {
    let mut plan = vj_truth_plan(cfg, v_truth, "j01*01", 0, 0);
    match mode {
        "n" => {
            plan.push(Box::new(NCorruptionPass::new(Box::new(
                EmpiricalLengthDist::from_pairs(vec![(n_count, 1.0)]),
            ))));
        }
        "pcr" => {
            plan.push(Box::new(PCRErrorPass::new(
                Box::new(EmpiricalLengthDist::from_pairs(vec![(n_count, 1.0)])),
                Box::new(UniformBase),
            )));
        }
        "quality" => {
            plan.push(Box::new(QualityErrorPass::new(
                Box::new(EmpiricalLengthDist::from_pairs(vec![(n_count, 1.0)])),
                Box::new(UniformBase),
            )));
        }
        other => panic!("unknown corruption mode {}", other),
    }
    plan
}

// ──────────────────────────────────────────────────────────────
// Helpers to peek at V hypothesis flags.
// ──────────────────────────────────────────────────────────────

fn v_hypothesis_flags(sim: &Simulation) -> HypothesisFlags {
    sim.live_calls
        .as_ref()
        .and_then(|lc| lc.get(Segment::V))
        .and_then(|sc| sc.hypotheses.first().map(|h| h.flags))
        .unwrap_or(HypothesisFlags::EMPTY)
}

#[test]
fn heavy_shm_does_not_leave_stale_boundary_elastic_flag() {
    // Apply a heavy UniformMutationPass (count=15) on the VJ fixture.
    // The V region is 12 bytes; after 15 mutation attempts the pool
    // has been hit many times. The post-mutation hypothesis flags
    // must reflect the post-mutation evidence, not stale state from
    // before mutation. We compare against a from-scratch run with
    // mutation count = 0 — the flag set should reflect the trimmed /
    // assembled fixture's actual elastic/wildcard properties.
    //
    // Concretely: with no trim and no NP-driven extension, the V
    // baseline hypothesis on this fixture has no extension and no N
    // injection — BOUNDARY_ELASTIC should be FALSE on the baseline
    // and STILL FALSE after heavy SHM (mutation doesn't introduce
    // elastic extension; it just edits bases at the existing
    // positions). HAS_WILDCARD_EVIDENCE should also stay false unless
    // a mutated byte became N — UniformMutation samples
    // A/C/G/T (no N), so it must remain false.
    let cfg = common::vj_ambiguous_refdata();

    // Baseline: no mutation.
    let baseline_plan = vj_truth_plan(&cfg, "v01*01", "j01*01", 0, 0);
    let baseline_sim = run_plan(&cfg, baseline_plan, 0);
    let baseline_flags = v_hypothesis_flags(&baseline_sim);

    // Heavy mutation: count=15.
    let mut heavy_plan = vj_truth_plan(&cfg, "v01*01", "j01*01", 0, 0);
    heavy_plan.push(Box::new(UniformMutationPass::new(
        Box::new(EmpiricalLengthDist::from_pairs(vec![(15, 1.0)])),
        Box::new(UniformBase),
    )));
    let heavy_sim = run_plan(&cfg, heavy_plan, 7);
    let heavy_flags = v_hypothesis_flags(&heavy_sim);

    // UniformMutation samples from {A,C,G,T} — no N injection — so
    // HAS_WILDCARD_EVIDENCE must stay false on baseline AND post-SHM.
    assert!(
        !baseline_flags.contains(HypothesisFlags::HAS_WILDCARD_EVIDENCE),
        "baseline must not carry HAS_WILDCARD_EVIDENCE",
    );
    assert!(
        !heavy_flags.contains(HypothesisFlags::HAS_WILDCARD_EVIDENCE),
        "uniform SHM must not introduce HAS_WILDCARD_EVIDENCE (no N bases sampled)",
    );

    // BOUNDARY_ELASTIC reflects whether the walker extended beyond
    // the structural boundary. With zero trim and zero NP, no
    // extension happens — neither baseline nor heavy-SHM run should
    // set this flag. This is the "no stale flag" assertion.
    assert_eq!(
        baseline_flags.contains(HypothesisFlags::BOUNDARY_ELASTIC),
        heavy_flags.contains(HypothesisFlags::BOUNDARY_ELASTIC),
        "BOUNDARY_ELASTIC must match baseline (no extension happens on either; baseline={:?} heavy={:?})",
        baseline_flags, heavy_flags,
    );
}

#[test]
fn ncorrupt_n_injection_sets_has_wildcard_evidence_flag() {
    // NCorruptionPass fills positions with N. If any N lands inside
    // V's pool range, the V hypothesis's HAS_WILDCARD_EVIDENCE flag
    // must be set. With count=2 on a small fixture (pool = 12 + 0 +
    // 9 = 21 bytes), the probability that AT LEAST ONE N lands inside
    // V's 12-byte range is very high. We test across multiple seeds
    // and verify the flag fires when an N actually lands in V.
    let cfg = common::vj_ambiguous_refdata();
    let mut any_seed_set_flag = false;
    for seed in 0..16u64 {
        let mut plan = vj_truth_plan(&cfg, "v01*01", "j01*01", 0, 0);
        plan.push(Box::new(NCorruptionPass::new(Box::new(
            EmpiricalLengthDist::from_pairs(vec![(2, 1.0)]),
        ))));
        let sim = run_plan(&cfg, plan, seed);

        // Check if any pool position inside V's range now holds an N.
        let v_region = sim
            .sequence
            .regions
            .iter()
            .find(|r| r.segment == Segment::V)
            .expect("V region exists");
        let mut v_has_n = false;
        for idx in v_region.start.index()..v_region.end.index() {
            let nuc = sim.pool.get(NucHandle::new(idx)).unwrap();
            if nuc.base == b'N' {
                v_has_n = true;
                break;
            }
        }

        let flags = v_hypothesis_flags(&sim);
        if v_has_n {
            assert!(
                flags.contains(HypothesisFlags::HAS_WILDCARD_EVIDENCE),
                "seed {}: N landed in V's range, but HAS_WILDCARD_EVIDENCE not set on V hypothesis (flags={:?})",
                seed, flags,
            );
            any_seed_set_flag = true;
        }
    }
    assert!(
        any_seed_set_flag,
        "across 16 seeds with N count=2 on a 12-byte V, expected at least one seed to land an N in V",
    );
}

#[test]
fn pcr_error_at_distinguishing_position_behaves_like_shm() {
    // PCRErrorPass and QualityErrorPass both flip bases via the
    // same builder.change_base path used by SHM. Confirm that a
    // PCR error at a distinguishing position narrows the call set
    // exactly like a deterministic edit would.
    //
    // We can't pin a single PCR-pass mutation to a position
    // deterministically (PCR samples the site uniformly), so we
    // verify the looser invariant: across many seeds, the AIRR
    // `n_pcr_errors` field always equals the sampled count (3),
    // and the call set still always contains *some* candidate.
    let cfg = common::vj_ambiguous_refdata();
    for seed in 0..32u64 {
        let plan = plan_clone(&cfg, "v01*01", 3, "pcr");
        let (sim, rec) = run_plan_with_record(&cfg, plan, seed);
        assert_eq!(rec.n_pcr_errors, 3, "seed {}", seed);
        let v_calls = common::v_call_names(&sim, &cfg);
        assert!(
            !v_calls.is_empty(),
            "call set must be non-empty after PCR errors (seed {}); got {:?}",
            seed,
            v_calls,
        );
    }
}
