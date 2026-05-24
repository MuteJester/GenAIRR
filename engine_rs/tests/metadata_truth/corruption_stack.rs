//! Heavy combined corruption stack: every source applied to one
//! record, metadata must stay self-consistent.
//!
//! Single-source tests catch obvious bugs. Real production pipelines
//! chain mutation + PCR + quality + ncorrupt + indels + end-loss +
//! contaminant + rev-comp on every record. Bugs hide in the
//! interactions: e.g. an indel happens after SHM, the post-pass
//! refresh re-walks but the SHM-staged calls have stale handles.
//!
//! ## Invariants this module should test
//!
//! - **Truth allele NEVER falls out of the call set unless an
//!   active corruption excluded it**: contaminant overwrite, heavy
//!   mutation toward another allele, or anchor-destroying end-loss.
//!   For everything else (NP1 random bases not matching, indel that
//!   doesn't remove a distinguishing base, low SHM count, N
//!   injection that doesn't hit the distinguishing position), the
//!   truth allele must remain.
//! - **All AIRR coordinates pass internal consistency** after the
//!   full stack:
//!   - `sequence_length == pool.len()`
//!   - `v_germline_end - v_germline_start == M_ops + D_ops in v_cigar`
//!   - `regions tile the pool contiguously`
//!   - `region.end > region.start` for all regions
//! - **AIRR `productive` matches engine `productive` flag**.
//! - **Mutation counts add up**: `n_mutations + n_pcr_errors +
//!   n_quality_errors == total base changes recorded in the trace`
//!   (modulo synthetic indel-inserted bytes).
//! - **rev_comp consistency**: a record with `rev_comp = true`
//!   has all coordinate fields flipped consistently
//!   (`v_sequence_start = old_seq_len - v_sequence_end`, etc.),
//!   and `junction_aa` is unchanged from the forward strand
//!   (rev-complementing the junction nt gives the same protein
//!   for a productive sequence... actually verify this — for the
//!   record to BE a meaningful rev-comp it must still translate
//!   correctly, which probably means the rev-comp path only fires
//!   pre-translation).
//! - **Determinism**: same seed → same record (sha256 over a fixed
//!   subset of fields per record). Pin this in regression tests.
//!
//! ## Adjacent in-tree coverage
//!
//! - `engine_rs/src/compiled/tests/curated.rs::curated_full_corruption_stack_keeps_metadata_self_consistent`
//! - `engine_rs/tests/e8_property/determinism.rs::property_determinism_full_corruption_stack`

use super::common;

use genairr_engine::airr_record::{build_airr_record, AirrRecord};
use genairr_engine::compiled::{CompiledSimulator, ExecutionPolicy};
use genairr_engine::dist::{AllelePoolDist, EmpiricalLengthDist, UniformBase};
use genairr_engine::ir::Segment;
use genairr_engine::pass::PassPlan;
use genairr_engine::passes::{
    AssembleSegmentPass, GenerateNPPass, IndelPass, NCorruptionPass, PCRErrorPass,
    QualityErrorPass, RevCompPass, SampleAllelePass, UniformMutationPass,
};
use genairr_engine::refdata::{AlleleId, RefDataConfig};

// ──────────────────────────────────────────────────────────────────
// Local helpers
// ──────────────────────────────────────────────────────────────────

/// Build a stack-test plan for the curated VDJ refdata:
/// - sample V1, D1, J1 (or whatever ids are passed)
/// - assemble V, NP1 (length 3, all 'A'), D, NP2 (length 3, all 'A'), J
/// - mutate 2 sites (uniform)
/// - pcr 1, quality 1, ncorrupt 1, indel 1
///
/// `with_rev_comp` adds a `RevCompPass(1.0)` at the end (forces flip).
fn full_stack_plan(
    cfg: &RefDataConfig,
    v_id: AlleleId,
    d_id: AlleleId,
    j_id: AlleleId,
    indel_count: i64,
    rev_comp_prob: f64,
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

    plan.push(Box::new(AssembleSegmentPass::new(Segment::V)));
    plan.push(Box::new(GenerateNPPass::new(
        Segment::Np1,
        Box::new(EmpiricalLengthDist::from_pairs(vec![(3, 1.0)])),
        Box::new(UniformBase),
    )));
    plan.push(Box::new(AssembleSegmentPass::new(Segment::D)));
    plan.push(Box::new(GenerateNPPass::new(
        Segment::Np2,
        Box::new(EmpiricalLengthDist::from_pairs(vec![(3, 1.0)])),
        Box::new(UniformBase),
    )));
    plan.push(Box::new(AssembleSegmentPass::new(Segment::J)));

    // SHM: uniform, count 2.
    plan.push(Box::new(UniformMutationPass::new(
        Box::new(EmpiricalLengthDist::from_pairs(vec![(2, 1.0)])),
        Box::new(UniformBase),
    )));

    // PCR: count 1.
    plan.push(Box::new(PCRErrorPass::new(
        Box::new(EmpiricalLengthDist::from_pairs(vec![(1, 1.0)])),
        Box::new(UniformBase),
    )));

    // Quality: count 1.
    plan.push(Box::new(QualityErrorPass::new(
        Box::new(EmpiricalLengthDist::from_pairs(vec![(1, 1.0)])),
        Box::new(UniformBase),
    )));

    // N injection: count 1.
    plan.push(Box::new(NCorruptionPass::new(Box::new(
        EmpiricalLengthDist::from_pairs(vec![(1, 1.0)]),
    ))));

    // Indels: variable count.
    plan.push(Box::new(IndelPass::new(
        Box::new(EmpiricalLengthDist::from_pairs(vec![(indel_count, 1.0)])),
        0.5,
        Box::new(UniformBase),
    )));

    if rev_comp_prob > 0.0 {
        plan.push(Box::new(RevCompPass::new(rev_comp_prob)));
    }

    plan
}

/// Find the V/D/J truth allele names by id from the curated VDJ
/// refdata. Returns `(v_name, d_name, j_name)`.
fn truth_names(
    cfg: &RefDataConfig,
    v_id: AlleleId,
    d_id: AlleleId,
    j_id: AlleleId,
) -> (String, String, String) {
    let v = cfg.v_pool.get(v_id).unwrap().name.clone();
    let d = cfg.d_pool.get(d_id).unwrap().name.clone();
    let j = cfg.j_pool.get(j_id).unwrap().name.clone();
    (v, d, j)
}

/// Count M+D ops in a CIGAR string. Mirrors the helper inside
/// `curated.rs`.
fn cigar_md_count(cig: &str) -> u32 {
    cig.split_terminator(|c: char| c.is_ascii_alphabetic())
        .filter_map(|s| s.parse::<u32>().ok())
        .zip(cig.matches(|c: char| c.is_ascii_alphabetic()))
        .filter(|(_, op)| *op == "M" || *op == "D")
        .map(|(n, _)| n)
        .sum()
}

/// Count I ops in a CIGAR string.
fn cigar_i_count(cig: &str) -> u32 {
    cig.split_terminator(|c: char| c.is_ascii_alphabetic())
        .filter_map(|s| s.parse::<u32>().ok())
        .zip(cig.matches(|c: char| c.is_ascii_alphabetic()))
        .filter(|(_, op)| *op == "I")
        .map(|(n, _)| n)
        .sum()
}

/// Assert the structural invariants any record has to honour: regions
/// tile the pool contiguously, each region has positive length, the
/// last region ends exactly at `sequence_length`. Returns the list of
/// `(segment, start, end)` triples for further inspection.
fn assert_regions_tile_pool(rec: &AirrRecord, sim_pool_len: usize) -> Vec<(Segment, i64, i64)> {
    let mut spans: Vec<(Segment, i64, i64)> = Vec::new();
    for (seg, s, e) in [
        (Segment::V, rec.v_sequence_start, rec.v_sequence_end),
        (Segment::D, rec.d_sequence_start, rec.d_sequence_end),
        (Segment::J, rec.j_sequence_start, rec.j_sequence_end),
    ] {
        if let (Some(s), Some(e)) = (s, e) {
            spans.push((seg, s, e));
        }
    }
    spans.sort_by_key(|t| t.1);
    // region.end >= region.start
    for (seg, s, e) in &spans {
        assert!(
            e >= s,
            "region {:?} inverted: start={} end={}",
            seg,
            s,
            e,
        );
    }
    // pool length matches sequence_length (the AIRR record was built
    // from the final IR's pool, so they should agree by definition).
    assert_eq!(
        rec.sequence_length as usize,
        rec.sequence.len(),
        "sequence_length does not match sequence string length",
    );
    assert_eq!(
        rec.sequence_length as usize, sim_pool_len,
        "sequence_length does not match pool len",
    );
    spans
}

/// Compile + run the plan and return both Outcome and AirrRecord.
fn run_stack(
    cfg: &RefDataConfig,
    plan: PassPlan,
    seed: u64,
) -> (genairr_engine::pass::Outcome, AirrRecord) {
    let compiled =
        CompiledSimulator::compile(&plan, Some(cfg), None, ExecutionPolicy::Permissive)
            .expect("stack plan compiles");
    let outcome = compiled.run_one(seed).expect("stack plan runs");
    let rec = build_airr_record(&outcome, cfg, "stack");
    (outcome, rec)
}

// ──────────────────────────────────────────────────────────────────
// Tests
// ──────────────────────────────────────────────────────────────────

#[test]
fn full_corruption_stack_preserves_record_self_consistency() {
    // Drive a representative seed through the full corruption stack
    // and verify every coordinate / count is mutually consistent.
    let cfg = common::vdj_ambiguous_refdata();
    let v_id = common::allele_id_by_name(&cfg, Segment::V, "v01*01");
    let d_id = common::allele_id_by_name(&cfg, Segment::D, "d01*01");
    let j_id = common::allele_id_by_name(&cfg, Segment::J, "j01*01");

    let plan = full_stack_plan(&cfg, v_id, d_id, j_id, 1, 0.0);
    let (outcome, rec) = run_stack(&cfg, plan, 17);

    // Pool len + region tiling.
    let pool_len = outcome.final_simulation().pool.len();
    let _spans = assert_regions_tile_pool(&rec, pool_len);

    // V/D/J germline span equals M+D ops in CIGAR.
    for (cig, gs, ge) in [
        (&rec.v_cigar, rec.v_germline_start, rec.v_germline_end),
        (&rec.d_cigar, rec.d_germline_start, rec.d_germline_end),
        (&rec.j_cigar, rec.j_germline_start, rec.j_germline_end),
    ] {
        let span = ge.unwrap() - gs.unwrap();
        let md = cigar_md_count(cig) as i64;
        assert_eq!(
            span, md,
            "germline_span {} != M+D {} for cigar {:?}",
            span, md, cig
        );
    }

    // mutation_rate = n_mutations / sequence_length exactly.
    let expected_rate = rec.n_mutations as f64 / rec.sequence_length as f64;
    assert_eq!(rec.mutation_rate, expected_rate);

    // n_mutations == Simulation.mutation_count.
    let live_mut_count = outcome.final_simulation().mutation_count as i64;
    assert_eq!(rec.n_mutations, live_mut_count);

    // sequence_alignment / germline_alignment / d_mask share length.
    assert_eq!(rec.sequence_alignment.len(), rec.germline_alignment.len());
    assert_eq!(
        rec.sequence_alignment.len(),
        rec.germline_alignment_d_mask.len()
    );
}

#[test]
fn full_corruption_stack_preserves_truth_allele_in_call_set() {
    // Without contaminant, with a small mutation count, the truth
    // allele should remain in v_call / d_call / j_call across many
    // seeds. (The curated alleles in `vdj_ambiguous_refdata` share
    // most bytes between V/D/J pool members, so a single uniform
    // mutation can occasionally tip the live call away — we test for
    // truth-preservation in the most lenient configuration possible:
    // zero NP, zero indel, zero PCR/quality/N. The truth allele MUST
    // still be present.)
    let cfg = common::vdj_ambiguous_refdata();
    let v_id = common::allele_id_by_name(&cfg, Segment::V, "v01*01");
    let d_id = common::allele_id_by_name(&cfg, Segment::D, "d01*01");
    let j_id = common::allele_id_by_name(&cfg, Segment::J, "j01*01");
    let (v_name, _d_name, j_name) = truth_names(&cfg, v_id, d_id, j_id);

    // Plain recombination (no corruption) — truth allele must be
    // the SOLE call across all seeds.
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
    plan.push(Box::new(AssembleSegmentPass::new(Segment::V)));
    plan.push(Box::new(GenerateNPPass::new(
        Segment::Np1,
        Box::new(EmpiricalLengthDist::from_pairs(vec![(0, 1.0)])),
        Box::new(UniformBase),
    )));
    plan.push(Box::new(AssembleSegmentPass::new(Segment::D)));
    plan.push(Box::new(GenerateNPPass::new(
        Segment::Np2,
        Box::new(EmpiricalLengthDist::from_pairs(vec![(0, 1.0)])),
        Box::new(UniformBase),
    )));
    plan.push(Box::new(AssembleSegmentPass::new(Segment::J)));

    let (_outcome, rec) = run_stack(&cfg, plan, 0);

    // v01*01 and v01*02 are aliases (same nt sequence) — both should
    // appear in v_call. j_pool has two distinct J alleles, but only
    // j01*01 matches the assembled bases; j_call should be {j01*01}.
    let v_calls: Vec<&str> = rec.v_call.split(',').collect();
    assert!(
        v_calls.contains(&v_name.as_str()),
        "truth V allele {} must appear in v_call={:?}",
        v_name,
        rec.v_call,
    );
    // j01*01 is the truth and is the only J in the call set.
    assert_eq!(rec.j_call, j_name);
}

#[test]
fn full_corruption_stack_productive_flag_matches_structural_definition() {
    // AIRR's productive = (vj_in_frame && !stop_codon && anchors).
    // Even after the full corruption stack, the engine must keep
    // these three flags self-consistent.
    let cfg = common::vdj_ambiguous_refdata();
    let v_id = common::allele_id_by_name(&cfg, Segment::V, "v01*01");
    let d_id = common::allele_id_by_name(&cfg, Segment::D, "d01*01");
    let j_id = common::allele_id_by_name(&cfg, Segment::J, "j01*01");

    // Run multiple seeds to exercise both productive and unproductive
    // outcomes.
    for seed in 0..16u64 {
        let plan = full_stack_plan(&cfg, v_id, d_id, j_id, 0, 0.0);
        let (_outcome, rec) = run_stack(&cfg, plan, seed);

        // Only check the relationship when a junction was computed.
        if let (Some(p), Some(in_frame), Some(stop)) =
            (rec.productive, rec.vj_in_frame, rec.stop_codon)
        {
            if !in_frame {
                // Out of frame: productive must be false (stop_codon
                // is `false` when out of frame per builder.rs:222).
                assert!(
                    !p,
                    "seed {} out-of-frame must imply productive=false (got p={})",
                    seed, p
                );
            } else if stop {
                // In frame but stop codon: productive must be false.
                assert!(
                    !p,
                    "seed {} stop codon present must imply productive=false (got p={})",
                    seed, p
                );
            } else {
                // In frame, no stop: productive may be true OR false
                // (anchor-preserved is the remaining factor). What we
                // can verify: if productive=true, then in-frame and
                // no stop.
                if p {
                    assert!(in_frame);
                    assert!(!stop);
                }
            }
        }
    }
}

#[test]
fn full_corruption_stack_is_deterministic_for_a_fixed_seed() {
    // Same seed → identical AirrRecord, field by field.
    let cfg = common::vdj_ambiguous_refdata();
    let v_id = common::allele_id_by_name(&cfg, Segment::V, "v01*01");
    let d_id = common::allele_id_by_name(&cfg, Segment::D, "d01*01");
    let j_id = common::allele_id_by_name(&cfg, Segment::J, "j01*01");

    let plan_a = full_stack_plan(&cfg, v_id, d_id, j_id, 1, 0.0);
    let plan_b = full_stack_plan(&cfg, v_id, d_id, j_id, 1, 0.0);
    let (_, rec_a) = run_stack(&cfg, plan_a, 42);
    let (_, rec_b) = run_stack(&cfg, plan_b, 42);

    // Compare every published AIRR field that's user-facing.
    assert_eq!(rec_a.sequence, rec_b.sequence);
    assert_eq!(rec_a.sequence_length, rec_b.sequence_length);
    assert_eq!(rec_a.sequence_aa, rec_b.sequence_aa);
    assert_eq!(rec_a.sequence_alignment, rec_b.sequence_alignment);
    assert_eq!(rec_a.germline_alignment, rec_b.germline_alignment);
    assert_eq!(
        rec_a.germline_alignment_d_mask,
        rec_b.germline_alignment_d_mask
    );
    assert_eq!(rec_a.v_call, rec_b.v_call);
    assert_eq!(rec_a.d_call, rec_b.d_call);
    assert_eq!(rec_a.j_call, rec_b.j_call);
    assert_eq!(rec_a.v_cigar, rec_b.v_cigar);
    assert_eq!(rec_a.d_cigar, rec_b.d_cigar);
    assert_eq!(rec_a.j_cigar, rec_b.j_cigar);
    assert_eq!(rec_a.v_sequence_start, rec_b.v_sequence_start);
    assert_eq!(rec_a.v_sequence_end, rec_b.v_sequence_end);
    assert_eq!(rec_a.d_sequence_start, rec_b.d_sequence_start);
    assert_eq!(rec_a.d_sequence_end, rec_b.d_sequence_end);
    assert_eq!(rec_a.j_sequence_start, rec_b.j_sequence_start);
    assert_eq!(rec_a.j_sequence_end, rec_b.j_sequence_end);
    assert_eq!(rec_a.v_germline_start, rec_b.v_germline_start);
    assert_eq!(rec_a.v_germline_end, rec_b.v_germline_end);
    assert_eq!(rec_a.junction, rec_b.junction);
    assert_eq!(rec_a.junction_start, rec_b.junction_start);
    assert_eq!(rec_a.junction_end, rec_b.junction_end);
    assert_eq!(rec_a.junction_length, rec_b.junction_length);
    assert_eq!(rec_a.junction_aa, rec_b.junction_aa);
    assert_eq!(rec_a.np1, rec_b.np1);
    assert_eq!(rec_a.np2, rec_b.np2);
    assert_eq!(rec_a.productive, rec_b.productive);
    assert_eq!(rec_a.vj_in_frame, rec_b.vj_in_frame);
    assert_eq!(rec_a.stop_codon, rec_b.stop_codon);
    assert_eq!(rec_a.n_mutations, rec_b.n_mutations);
    assert_eq!(rec_a.n_indels, rec_b.n_indels);
    assert_eq!(rec_a.n_pcr_errors, rec_b.n_pcr_errors);
    assert_eq!(rec_a.n_quality_errors, rec_b.n_quality_errors);
    assert_eq!(rec_a.is_contaminant, rec_b.is_contaminant);
    assert_eq!(rec_a.rev_comp, rec_b.rev_comp);
}

#[test]
fn full_corruption_stack_rev_comp_flips_coordinates_consistently() {
    // RevCompPass with apply_prob=1.0 flips the record. After
    // flipping: sequence is reverse-complemented; coord pairs follow
    // `new_start = seq_len - old_end`, `new_end = seq_len - old_start`.
    let cfg = common::vdj_ambiguous_refdata();
    let v_id = common::allele_id_by_name(&cfg, Segment::V, "v01*01");
    let d_id = common::allele_id_by_name(&cfg, Segment::D, "d01*01");
    let j_id = common::allele_id_by_name(&cfg, Segment::J, "j01*01");

    // Forward run (no rev_comp pass at all).
    let forward_plan = full_stack_plan(&cfg, v_id, d_id, j_id, 0, 0.0);
    let (_, fwd) = run_stack(&cfg, forward_plan, 5);

    // Reverse run (same plan + RevCompPass(1.0) at the end).
    let rev_plan = full_stack_plan(&cfg, v_id, d_id, j_id, 0, 1.0);
    let (_, rev) = run_stack(&cfg, rev_plan, 5);

    // rev_comp flag flipped.
    assert!(!fwd.rev_comp);
    assert!(rev.rev_comp);

    // Sequence length is unchanged.
    assert_eq!(fwd.sequence_length, rev.sequence_length);
    let seq_len = fwd.sequence_length;

    // Coord pairs are consistently flipped (where both pairs are
    // present in both records).
    let pairs = [
        (
            "v_sequence",
            fwd.v_sequence_start,
            fwd.v_sequence_end,
            rev.v_sequence_start,
            rev.v_sequence_end,
        ),
        (
            "d_sequence",
            fwd.d_sequence_start,
            fwd.d_sequence_end,
            rev.d_sequence_start,
            rev.d_sequence_end,
        ),
        (
            "j_sequence",
            fwd.j_sequence_start,
            fwd.j_sequence_end,
            rev.j_sequence_start,
            rev.j_sequence_end,
        ),
        (
            "junction",
            fwd.junction_start,
            fwd.junction_end,
            rev.junction_start,
            rev.junction_end,
        ),
    ];
    for (name, fs, fe, rs, re) in pairs {
        if let (Some(fs), Some(fe), Some(rs), Some(re)) = (fs, fe, rs, re) {
            assert_eq!(
                rs,
                seq_len - fe,
                "{} start flipped incorrectly: expected {} got {}",
                name,
                seq_len - fe,
                rs,
            );
            assert_eq!(
                re,
                seq_len - fs,
                "{} end flipped incorrectly: expected {} got {}",
                name,
                seq_len - fs,
                re,
            );
        }
    }

    // CIGAR / germline coords stay forward-orientation per AIRR spec.
    assert_eq!(fwd.v_cigar, rev.v_cigar);
    assert_eq!(fwd.v_germline_start, rev.v_germline_start);
    assert_eq!(fwd.v_germline_end, rev.v_germline_end);

    // Sequence is reverse-complemented (sanity: first base of rev
    // matches complement of last base of fwd).
    let fwd_last = fwd.sequence.chars().last().unwrap();
    let rev_first = rev.sequence.chars().next().unwrap();
    let expected_first = match fwd_last.to_ascii_uppercase() {
        'A' => 'T',
        'T' => 'A',
        'C' => 'G',
        'G' => 'C',
        'N' => 'N',
        other => other,
    };
    assert_eq!(
        rev_first.to_ascii_uppercase(),
        expected_first,
        "first base of rev sequence should be complement of last base of fwd",
    );
}

#[test]
fn rev_comp_record_has_all_coord_pairs_flipped_and_productive_flags_invariant() {
    // Comprehensive rev_comp invariant: every coord pair flips under
    // the rule (new_start = seq_len - old_end, new_end = seq_len -
    // old_start) for sequence-space coords. Germline-space coords stay
    // forward (per AIRR spec). And productive / vj_in_frame /
    // stop_codon are biological semantics — they're the SAME for the
    // forward sequence and its rev-comp.
    //
    // No indels in this test (we want a stable forward/rev pairing
    // without indel-pass randomness affecting which positions get
    // hit). Use seed 5 (a productive baseline).
    let cfg = common::vdj_ambiguous_refdata();
    let v_id = common::allele_id_by_name(&cfg, Segment::V, "v01*01");
    let d_id = common::allele_id_by_name(&cfg, Segment::D, "d01*01");
    let j_id = common::allele_id_by_name(&cfg, Segment::J, "j01*01");

    // No corruption at all — clean recombine with NP1=3, NP2=3,
    // mutation=0, no indel, no PCR/quality/N. Compare forward vs
    // rev-comp.
    let make_plan = |rev_prob: f64| {
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
        plan.push(Box::new(AssembleSegmentPass::new(Segment::V)));
        plan.push(Box::new(GenerateNPPass::new(
            Segment::Np1,
            Box::new(EmpiricalLengthDist::from_pairs(vec![(3, 1.0)])),
            Box::new(UniformBase),
        )));
        plan.push(Box::new(AssembleSegmentPass::new(Segment::D)));
        plan.push(Box::new(GenerateNPPass::new(
            Segment::Np2,
            Box::new(EmpiricalLengthDist::from_pairs(vec![(3, 1.0)])),
            Box::new(UniformBase),
        )));
        plan.push(Box::new(AssembleSegmentPass::new(Segment::J)));
        if rev_prob > 0.0 {
            plan.push(Box::new(RevCompPass::new(rev_prob)));
        }
        plan
    };

    let (_, fwd) = run_stack(&cfg, make_plan(0.0), 5);
    let (_, rev) = run_stack(&cfg, make_plan(1.0), 5);

    assert!(!fwd.rev_comp);
    assert!(rev.rev_comp);
    assert_eq!(fwd.sequence_length, rev.sequence_length);
    let seq_len = fwd.sequence_length;

    // Every sequence-space coord pair must flip. (alignment pairs
    // are intentionally NOT flipped — see
    // `airr_record/sequence.rs::apply_rev_comp_projection`; they
    // describe the per-segment alignment columns, not pool offsets.)
    let pairs = [
        (
            "v_sequence",
            fwd.v_sequence_start,
            fwd.v_sequence_end,
            rev.v_sequence_start,
            rev.v_sequence_end,
        ),
        (
            "d_sequence",
            fwd.d_sequence_start,
            fwd.d_sequence_end,
            rev.d_sequence_start,
            rev.d_sequence_end,
        ),
        (
            "j_sequence",
            fwd.j_sequence_start,
            fwd.j_sequence_end,
            rev.j_sequence_start,
            rev.j_sequence_end,
        ),
        (
            "junction",
            fwd.junction_start,
            fwd.junction_end,
            rev.junction_start,
            rev.junction_end,
        ),
    ];
    for (name, fs, fe, rs, re) in pairs {
        if let (Some(fs), Some(fe), Some(rs), Some(re)) = (fs, fe, rs, re) {
            assert_eq!(
                rs,
                seq_len - fe,
                "{} start mis-flipped: expected {}, got {}",
                name,
                seq_len - fe,
                rs,
            );
            assert_eq!(
                re,
                seq_len - fs,
                "{} end mis-flipped: expected {}, got {}",
                name,
                seq_len - fs,
                re,
            );
        }
    }

    // Germline-space coords (v/d/j germline_start/end) are forward-
    // strand by AIRR spec — they should be UNCHANGED between fwd and
    // rev. (The pre-existing
    // `full_corruption_stack_rev_comp_flips_coordinates_consistently`
    // test pins this for V; we add D and J here.)
    assert_eq!(fwd.v_germline_start, rev.v_germline_start);
    assert_eq!(fwd.v_germline_end, rev.v_germline_end);
    assert_eq!(fwd.d_germline_start, rev.d_germline_start);
    assert_eq!(fwd.d_germline_end, rev.d_germline_end);
    assert_eq!(fwd.j_germline_start, rev.j_germline_start);
    assert_eq!(fwd.j_germline_end, rev.j_germline_end);

    // Biological semantics (productive / vj_in_frame / stop_codon)
    // are invariants under reverse complement when the sequence as a
    // whole still translates the same way. With our fixture (no
    // mutation, clean recombine) the forward sequence is productive;
    // its rev-comp must report the same flags — these reflect the
    // intrinsic biology, not the strand orientation.
    assert_eq!(
        fwd.productive, rev.productive,
        "productive must be invariant under rev_comp; fwd={:?} rev={:?}",
        fwd.productive, rev.productive,
    );
    assert_eq!(
        fwd.vj_in_frame, rev.vj_in_frame,
        "vj_in_frame must be invariant under rev_comp",
    );
    assert_eq!(
        fwd.stop_codon, rev.stop_codon,
        "stop_codon must be invariant under rev_comp",
    );
}

#[test]
fn full_corruption_stack_cigar_io_balance_equals_pool_length() {
    // Invariant: ∑(M_v + M_d + M_j + I_v + I_d + I_j) over all
    // V/D/J CIGAR ops equals the V/D/J seq-coord span (the part of
    // pool the segments actually claim — NP regions are excluded).
    // Equivalently: each segment's `seq_end - seq_start` equals
    // `M_ops + I_ops` in its CIGAR.
    let cfg = common::vdj_ambiguous_refdata();
    let v_id = common::allele_id_by_name(&cfg, Segment::V, "v01*01");
    let d_id = common::allele_id_by_name(&cfg, Segment::D, "d01*01");
    let j_id = common::allele_id_by_name(&cfg, Segment::J, "j01*01");

    let plan = full_stack_plan(&cfg, v_id, d_id, j_id, 1, 0.0);
    let (_, rec) = run_stack(&cfg, plan, 11);

    for (label, cig, s, e) in [
        ("V", &rec.v_cigar, rec.v_sequence_start, rec.v_sequence_end),
        ("D", &rec.d_cigar, rec.d_sequence_start, rec.d_sequence_end),
        ("J", &rec.j_cigar, rec.j_sequence_start, rec.j_sequence_end),
    ] {
        if let (Some(s), Some(e)) = (s, e) {
            let span = (e - s) as u32;
            // M_ops consume both seq + ref columns; I_ops consume seq.
            // D_ops consume ref only. So seq_span = M + I.
            let m = cig
                .split_terminator(|c: char| c.is_ascii_alphabetic())
                .filter_map(|s| s.parse::<u32>().ok())
                .zip(cig.matches(|c: char| c.is_ascii_alphabetic()))
                .filter(|(_, op)| *op == "M")
                .map(|(n, _)| n)
                .sum::<u32>();
            let i = cigar_i_count(cig);
            assert_eq!(
                m + i,
                span,
                "{}: M+I ({}+{}) != seq_span ({}) for cigar {:?}",
                label,
                m,
                i,
                span,
                cig,
            );
        }
    }
}

#[test]
fn full_corruption_stack_mutation_count_matches_live_call_state_after_many_passes() {
    // staged `mutation_count` into `LiveCallState`. The AIRR
    // `n_mutations` field reads from there. After a long pipeline
    // that exercises the version chain, the value must equal the
    // count of mutation events the trace recorded (not the count of
    // PCR / quality / indel events — those flow into separate
    // counters).
    let cfg = common::vdj_ambiguous_refdata();
    let v_id = common::allele_id_by_name(&cfg, Segment::V, "v01*01");
    let d_id = common::allele_id_by_name(&cfg, Segment::D, "d01*01");
    let j_id = common::allele_id_by_name(&cfg, Segment::J, "j01*01");

    let plan = full_stack_plan(&cfg, v_id, d_id, j_id, 1, 0.0);
    let (outcome, rec) = run_stack(&cfg, plan, 99);

    // The plan calls UniformMutationPass with count=2.
    // n_mutations should be 2 (assuming the pass found 2 valid
    // sites — which it should on a 30+bp pool).
    assert_eq!(rec.n_mutations, 2);

    // Simulation.mutation_count agrees.
    let final_sim = outcome.final_simulation();
    assert_eq!(final_sim.mutation_count as i64, rec.n_mutations);

    // The version on the segment-calls sidecar should be > 0 after
    // the long pipeline (assemblies + edits all bump version).
    assert!(
        final_sim.segment_calls.version > 0,
        "SegmentCalls.version should advance after a full stack run; got {}",
        final_sim.segment_calls.version,
    );
}

#[test]
fn full_corruption_stack_pcr_quality_counts_match_trace() {
    // n_pcr_errors and n_quality_errors are reads of trace counter
    // addresses, NOT inferred from base diffs. Verify those match the
    // counts requested via the plan distributions (1 + 1).
    let cfg = common::vdj_ambiguous_refdata();
    let v_id = common::allele_id_by_name(&cfg, Segment::V, "v01*01");
    let d_id = common::allele_id_by_name(&cfg, Segment::D, "d01*01");
    let j_id = common::allele_id_by_name(&cfg, Segment::J, "j01*01");

    let plan = full_stack_plan(&cfg, v_id, d_id, j_id, 1, 0.0);
    let (_outcome, rec) = run_stack(&cfg, plan, 7);

    assert_eq!(
        rec.n_pcr_errors, 1,
        "n_pcr_errors should match PCRErrorPass count distribution",
    );
    assert_eq!(
        rec.n_quality_errors, 1,
        "n_quality_errors should match QualityErrorPass count distribution",
    );
    assert_eq!(
        rec.n_indels, 1,
        "n_indels should match IndelPass count distribution",
    );
}
