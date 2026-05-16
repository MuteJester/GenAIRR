//! Property-based tests over key engine invariants.
//!
//! These tests sweep many seeds (typically 100+) through realistic
//! pipelines — full V(D)J recombination plus corruption / mutation
//! passes — and assert structural invariants that must hold across
//! **every** run. They are the architectural soundness gate:
//! if any pass silently breaks codon-rail consistency, persistent IR,
//! trace faithfulness, or determinism, a property test fires.
//!
//! No external proptest crate is used. The seed sweep drives the
//! existing deterministic SplitMix64 pipeline; with 100–200 seeds
//! and the synthetic VJ / VDJ refdata fixtures from D.7, every test
//! finishes in milliseconds.
//!
//! Invariants under test:
//!
//! 1. **Persistent IR** — every pass leaves its input simulation
//!    byte-identical to its pre-pass state, regardless of seed.
//! 2. **Codon-rail consistency** — for every region in every
//!    post-pass IR, `region.amino_acids` equals a fresh recompute
//!    against the post-pass pool.
//! 3. **Pool-length math after IndelPass** — final pool length is
//!    `initial + insertions − deletions` (modulo the empty-pool
//!    deletion clamp).
//! 4. **INDEL_INSERTED bookkeeping** — count of flagged nucleotides
//!    in the post-pass pool equals count of `kind[i]=Bool(true)`
//!    in the trace.
//! 5. **Trace faithfulness** — for `UniformMutationPass` the
//!    last-recorded base at each site appears at that pool index
//!    in the post-pass IR.
//! 6. **Region-range validity** — every region satisfies
//!    `0 ≤ start < end ≤ pool.len()` after every pass.
//! 7. **Determinism** — same seed + same plan ⇒ identical trace
//!    *and* identical pool bytes.
//! 8. **Productive admits ⇒ in-frame junction** — under
//!    `respect=[productive()]`, the assembled junction is in-frame
//!    on every seed.
//! 9. **Frame-phase consistency** — every region's `frame_phase`
//!    equals cumulative upstream region length modulo 3.

use genairr_engine::contract::productive;
use genairr_engine::dist::{AllelePoolDist, Distribution, EmpiricalLengthDist, UniformBase};
use genairr_engine::ir::{flag, NucHandle, Nucleotide, Region, Segment, Simulation};
use genairr_engine::junction::compute_junction;
use genairr_engine::pass::{PassPlan, PassRuntime};
use genairr_engine::passes::{
    AssembleSegmentPass, ContaminantPass, GenerateNPPass, IndelPass, PCRErrorPass,
    QualityErrorPass, S5FMutationPass, SampleAllelePass, UniformMutationPass,
};
use genairr_engine::refdata::{Allele, ChainType, RefDataConfig};
use genairr_engine::s5f::{S5FKernel, S5F_NUM_CONTEXTS, S5F_SUBSTITUTION_LEN};
use genairr_engine::trace::ChoiceValue;

// ──────────────────────────────────────────────────────────────────
// Fixtures — refdata + base sims + uniform S5F kernel
// ──────────────────────────────────────────────────────────────────

const SEED_RANGE: u64 = 100;

/// Synthetic VJ refdata. Mirrors the D.7 fixture so the productive
/// invariant test exercises the same constraint surface.
fn vj_refdata() -> RefDataConfig {
    let mut cfg = RefDataConfig::empty(ChainType::Vj);
    let _ = cfg.v_pool.push(Allele {
        name: "v_test*01".into(),
        gene: "v_test".into(),
        seq: b"AAACCCGGG".to_vec(),
        segment: Segment::V,
        anchor: Some(6),
    });
    let _ = cfg.j_pool.push(Allele {
        name: "j_test*01".into(),
        gene: "j_test".into(),
        seq: b"TTTAAA".to_vec(),
        segment: Segment::J,
        anchor: Some(0),
    });
    cfg
}

/// Build a VJ recombination plan. Same shape as the D.7 plan.
fn vj_plan(refdata: &RefDataConfig) -> PassPlan {
    let length_dist = || {
        Box::new(EmpiricalLengthDist::from_pairs(
            (0..10).map(|i| (i, 1.0)).collect::<Vec<_>>(),
        ))
    };
    let mut plan = PassPlan::new();
    plan.push(Box::new(SampleAllelePass::new(
        Segment::V,
        Box::new(AllelePoolDist::uniform(&refdata.v_pool)),
    )));
    plan.push(Box::new(SampleAllelePass::new(
        Segment::J,
        Box::new(AllelePoolDist::uniform(&refdata.j_pool)),
    )));
    plan.push(Box::new(AssembleSegmentPass::new(Segment::V)));
    plan.push(Box::new(GenerateNPPass::new(
        Segment::Np1,
        length_dist(),
        Box::new(UniformBase),
    )));
    plan.push(Box::new(AssembleSegmentPass::new(Segment::J)));
    plan
}

/// Helper: a 12bp germline V region. Used as the canonical input
/// to single-pass invariant sweeps so we don't pay the recombination
/// cost on every seed.
fn assembled_v_sim() -> Simulation {
    let mut sim = Simulation::new();
    for (i, b) in b"AAACCCGGGTTT".iter().enumerate() {
        let (next, _) = sim.with_nucleotide_pushed(Nucleotide::germline(*b, i as u16, Segment::V));
        sim = next;
    }
    let region = Region::new(Segment::V, NucHandle::new(0), NucHandle::new(12))
        .with_codon_rail_recomputed(&sim.pool);
    sim.with_region_added(region)
}

/// Uniform S5F kernel: every context has equal mutability and equal
/// substitution probability. Used so the property tests don't depend
/// on a specific empirical kernel; mechanical invariants are what's
/// under test.
fn uniform_s5f() -> S5FKernel {
    S5FKernel::new(
        vec![1.0; S5F_NUM_CONTEXTS],
        vec![0.25; S5F_SUBSTITUTION_LEN],
    )
}

/// Helper: assert every region's stored amino_acids matches a fresh
/// recompute against the given pool. The single "is the codon rail
/// in sync with the pool?" check used by every property test.
fn assert_codon_rails_consistent(sim: &Simulation, label: &str) {
    for (i, region) in sim.sequence.regions.iter().enumerate() {
        let fresh = region.with_codon_rail_recomputed(&sim.pool);
        assert_eq!(
            region.amino_acids, fresh.amino_acids,
            "{}: region[{}] codon rail stale (stored {:?}, fresh {:?})",
            label, i, region.amino_acids, fresh.amino_acids
        );
    }
}

/// Helper: assert every region has `0 ≤ start ≤ end ≤ pool.len()`.
/// Empty regions (`start == end`) are valid — e.g., a zero-length
/// NP1 is biologically meaningful and `Region::is_empty()` exists
/// to express it. Inverted ranges (`start > end`) are not.
fn assert_region_ranges_valid(sim: &Simulation, label: &str) {
    let pool_len = sim.pool.len() as u32;
    for (i, region) in sim.sequence.regions.iter().enumerate() {
        let s = region.start.index();
        let e = region.end.index();
        assert!(
            s <= e,
            "{}: region[{}] has inverted range start={} end={}",
            label,
            i,
            s,
            e
        );
        assert!(
            e <= pool_len,
            "{}: region[{}] end {} exceeds pool len {}",
            label,
            i,
            e,
            pool_len
        );
    }
}

/// Helper: assert each region's frame phase matches the sequence-level
/// codon frame implied by all preceding region lengths.
fn assert_frame_phases_consistent(sim: &Simulation, label: &str) {
    let mut cumulative_len = 0u64;
    for (i, region) in sim.sequence.regions.iter().enumerate() {
        let expected = (cumulative_len % 3) as u8;
        assert_eq!(
            region.frame_phase, expected,
            "{}: region[{}] has stale frame_phase {} (expected {})",
            label, i, region.frame_phase, expected
        );
        cumulative_len = cumulative_len.saturating_add(region.len() as u64);
    }
}

// ──────────────────────────────────────────────────────────────────
// Property 1: Persistent IR — input simulation untouched after pass
// ──────────────────────────────────────────────────────────────────

#[test]
fn property_persistent_ir_uniform_mutation() {
    // Every corruption / mutation pass must respect the persistent
    // IR contract: running the pass on a clone leaves the original
    // byte-identical.
    for seed in 0..SEED_RANGE {
        let sim = assembled_v_sim();
        let pre_len = sim.pool.len();
        let pre_region_end = sim.sequence.regions[0].end.index();
        let pre_amino = sim.sequence.regions[0].amino_acids.clone();

        let mut plan = PassPlan::new();
        plan.push(Box::new(UniformMutationPass::new(
            Box::new(EmpiricalLengthDist::from_pairs(vec![(3, 1.0)])),
            Box::new(UniformBase),
        )));
        let _ = PassRuntime::execute(&plan, sim.clone(), seed);

        // sim still byte-identical.
        assert_eq!(sim.pool.len(), pre_len);
        assert_eq!(sim.sequence.regions[0].end.index(), pre_region_end);
        assert_eq!(sim.sequence.regions[0].amino_acids, pre_amino);
    }
}

#[test]
fn property_persistent_ir_indel_pass() {
    // Same persistence invariant for the pool-length-changing pass.
    for seed in 0..SEED_RANGE {
        let sim = assembled_v_sim();
        let pre_len = sim.pool.len();
        let pre_region_end = sim.sequence.regions[0].end.index();
        let pre_amino = sim.sequence.regions[0].amino_acids.clone();

        let mut plan = PassPlan::new();
        plan.push(Box::new(IndelPass::new(
            Box::new(EmpiricalLengthDist::from_pairs(vec![(4, 1.0)])),
            0.5,
            Box::new(UniformBase),
        )));
        let _ = PassRuntime::execute(&plan, sim.clone(), seed);

        assert_eq!(sim.pool.len(), pre_len);
        assert_eq!(sim.sequence.regions[0].end.index(), pre_region_end);
        assert_eq!(sim.sequence.regions[0].amino_acids, pre_amino);
    }
}

// ──────────────────────────────────────────────────────────────────
// Property 2: Codon-rail consistency across every corruption pass
// ──────────────────────────────────────────────────────────────────

#[test]
fn property_codon_rail_consistency_uniform_mutation() {
    let mut plan = PassPlan::new();
    plan.push(Box::new(UniformMutationPass::new(
        Box::new(EmpiricalLengthDist::from_pairs(vec![(3, 1.0)])),
        Box::new(UniformBase),
    )));

    for seed in 0..SEED_RANGE {
        let outcome = PassRuntime::execute(&plan, assembled_v_sim(), seed);
        assert_codon_rails_consistent(
            outcome.final_simulation(),
            &format!("UniformMutation seed={}", seed),
        );
    }
}

#[test]
fn property_codon_rail_consistency_s5f() {
    let mut plan = PassPlan::new();
    plan.push(Box::new(S5FMutationPass::new(
        uniform_s5f(),
        Box::new(EmpiricalLengthDist::from_pairs(vec![(3, 1.0)])),
    )));

    for seed in 0..SEED_RANGE {
        let outcome = PassRuntime::execute(&plan, assembled_v_sim(), seed);
        assert_codon_rails_consistent(
            outcome.final_simulation(),
            &format!("S5FMutation seed={}", seed),
        );
    }
}

#[test]
fn property_codon_rail_consistency_pcr() {
    let mut plan = PassPlan::new();
    plan.push(Box::new(PCRErrorPass::new(
        Box::new(EmpiricalLengthDist::from_pairs(vec![(3, 1.0)])),
        Box::new(UniformBase),
    )));

    for seed in 0..SEED_RANGE {
        let outcome = PassRuntime::execute(&plan, assembled_v_sim(), seed);
        assert_codon_rails_consistent(
            outcome.final_simulation(),
            &format!("PCRError seed={}", seed),
        );
    }
}

#[test]
fn property_codon_rail_consistency_quality() {
    let mut plan = PassPlan::new();
    plan.push(Box::new(QualityErrorPass::new(
        Box::new(EmpiricalLengthDist::from_pairs(vec![(3, 1.0)])),
        Box::new(UniformBase),
    )));

    for seed in 0..SEED_RANGE {
        let outcome = PassRuntime::execute(&plan, assembled_v_sim(), seed);
        assert_codon_rails_consistent(
            outcome.final_simulation(),
            &format!("QualityError seed={}", seed),
        );
    }
}

#[test]
fn property_codon_rail_consistency_contaminant() {
    let mut plan = PassPlan::new();
    plan.push(Box::new(ContaminantPass::new(0.5, Box::new(UniformBase))));

    for seed in 0..SEED_RANGE {
        let outcome = PassRuntime::execute(&plan, assembled_v_sim(), seed);
        assert_codon_rails_consistent(
            outcome.final_simulation(),
            &format!("Contaminant seed={}", seed),
        );
    }
}

#[test]
fn property_codon_rail_consistency_indel() {
    let mut plan = PassPlan::new();
    plan.push(Box::new(IndelPass::new(
        Box::new(EmpiricalLengthDist::from_pairs(vec![(3, 1.0)])),
        0.5,
        Box::new(UniformBase),
    )));

    for seed in 0..SEED_RANGE {
        let outcome = PassRuntime::execute(&plan, assembled_v_sim(), seed);
        assert_codon_rails_consistent(outcome.final_simulation(), &format!("Indel seed={}", seed));
    }
}

#[test]
fn property_codon_rail_consistency_full_corruption_stack() {
    // Compound pass: SHM + PCR + quality + indel + contaminant. The
    // compound is the realistic corruption pipeline shape, and every
    // intermediate IR must keep the rail consistent so the next pass
    // reads accurate amino_acids.
    let mut plan = PassPlan::new();
    plan.push(Box::new(S5FMutationPass::new(
        uniform_s5f(),
        Box::new(EmpiricalLengthDist::from_pairs(vec![(2, 1.0)])),
    )));
    plan.push(Box::new(PCRErrorPass::new(
        Box::new(EmpiricalLengthDist::from_pairs(vec![(1, 1.0)])),
        Box::new(UniformBase),
    )));
    plan.push(Box::new(IndelPass::new(
        Box::new(EmpiricalLengthDist::from_pairs(vec![(2, 1.0)])),
        0.5,
        Box::new(UniformBase),
    )));
    plan.push(Box::new(QualityErrorPass::new(
        Box::new(EmpiricalLengthDist::from_pairs(vec![(1, 1.0)])),
        Box::new(UniformBase),
    )));
    plan.push(Box::new(ContaminantPass::new(0.1, Box::new(UniformBase))));

    for seed in 0..SEED_RANGE {
        let outcome = PassRuntime::execute(&plan, assembled_v_sim(), seed);
        assert_codon_rails_consistent(
            outcome.final_simulation(),
            &format!("full-corruption seed={}", seed),
        );
        assert_region_ranges_valid(
            outcome.final_simulation(),
            &format!("full-corruption seed={}", seed),
        );
    }
}

// ──────────────────────────────────────────────────────────────────
// Property 3: IndelPass pool-length math
// ──────────────────────────────────────────────────────────────────

#[test]
fn property_indel_pool_length_matches_trace_arithmetic() {
    // For every seed: count insertions / deletions in the trace,
    // verify final_len = initial_len + insertions − deletions.
    // (The empty-pool clamp records site=−1 and is a no-op, so we
    // exclude clamped entries from the deletion count.)
    let mut plan = PassPlan::new();
    plan.push(Box::new(IndelPass::new(
        Box::new(EmpiricalLengthDist::from_pairs(vec![(5, 1.0)])),
        0.5,
        Box::new(UniformBase),
    )));

    for seed in 0..SEED_RANGE {
        let initial_len = assembled_v_sim().pool.len() as i64;
        let outcome = PassRuntime::execute(&plan, assembled_v_sim(), seed);

        let count = match outcome.trace.find("corrupt.indel.count").unwrap().value {
            ChoiceValue::Int(n) => n,
            _ => unreachable!(),
        };
        assert_eq!(count, 5, "seed {} had unexpected count", seed);

        let mut insertions = 0i64;
        let mut deletions = 0i64;
        for i in 0..5 {
            let kind = matches!(
                outcome
                    .trace
                    .find(&format!("corrupt.indel.kind[{}]", i))
                    .unwrap()
                    .value,
                ChoiceValue::Bool(true)
            );
            let site = match outcome
                .trace
                .find(&format!("corrupt.indel.site[{}]", i))
                .unwrap()
                .value
            {
                ChoiceValue::Int(s) => s,
                _ => unreachable!(),
            };
            if kind {
                insertions += 1;
            } else if site >= 0 {
                deletions += 1;
            }
            // site == -1 means empty-pool deletion clamp; not counted.
        }

        let expected_len = initial_len + insertions - deletions;
        assert_eq!(
            outcome.final_simulation().pool.len() as i64,
            expected_len,
            "seed {}: expected pool len {} (initial {} + {} ins − {} del), got {}",
            seed,
            expected_len,
            initial_len,
            insertions,
            deletions,
            outcome.final_simulation().pool.len()
        );
    }
}

// ──────────────────────────────────────────────────────────────────
// Property 4: INDEL_INSERTED bookkeeping
// ──────────────────────────────────────────────────────────────────

#[test]
fn property_indel_inserted_flag_count_at_most_traced_insertions() {
    // Every flagged nucleotide came from an insertion event, but a
    // *later* deletion can remove an earlier insertion. So the
    // post-pass flagged count is bounded above by the trace's
    // insertion count — equality only holds when no event undoes
    // an earlier one. We assert the bound and additionally check
    // every flagged nucleotide carries `GermlinePos::NONE` (the other
    // half of the insertion provenance contract).
    let mut plan = PassPlan::new();
    plan.push(Box::new(IndelPass::new(
        Box::new(EmpiricalLengthDist::from_pairs(vec![(6, 1.0)])),
        0.5,
        Box::new(UniformBase),
    )));

    for seed in 0..SEED_RANGE {
        let outcome = PassRuntime::execute(&plan, assembled_v_sim(), seed);
        let final_sim = outcome.final_simulation();

        let mut flagged = 0u32;
        for i in 0..final_sim.pool.len() as u32 {
            let n = final_sim.pool.get(NucHandle::new(i)).unwrap();
            if n.flags.contains(flag::INDEL_INSERTED) {
                flagged += 1;
                assert!(
                    n.germline_pos.is_none(),
                    "seed {}: flagged nucleotide at {} has germline pos {:?}",
                    seed,
                    i,
                    n.germline_pos.get()
                );
            }
        }

        let traced_insertions = (0..6)
            .filter(|i| {
                matches!(
                    outcome
                        .trace
                        .find(&format!("corrupt.indel.kind[{}]", i))
                        .unwrap()
                        .value,
                    ChoiceValue::Bool(true)
                )
            })
            .count() as u32;

        assert!(
            flagged <= traced_insertions,
            "seed {}: flagged count {} exceeds traced insertions {}",
            seed,
            flagged,
            traced_insertions
        );
    }
}

#[test]
fn property_indel_inserted_flag_count_equals_insertions_when_no_deletions() {
    // The strong form of the bookkeeping property: with p_ins=1.0
    // every event is an insertion, so flagged count exactly equals
    // count.
    let mut plan = PassPlan::new();
    plan.push(Box::new(IndelPass::new(
        Box::new(EmpiricalLengthDist::from_pairs(vec![(5, 1.0)])),
        1.0,
        Box::new(UniformBase),
    )));

    for seed in 0..SEED_RANGE {
        let outcome = PassRuntime::execute(&plan, assembled_v_sim(), seed);
        let final_sim = outcome.final_simulation();

        let flagged = (0..final_sim.pool.len() as u32)
            .filter(|i| {
                final_sim
                    .pool
                    .get(NucHandle::new(*i))
                    .unwrap()
                    .flags
                    .contains(flag::INDEL_INSERTED)
            })
            .count();

        assert_eq!(
            flagged, 5,
            "seed {}: expected 5 flagged nucleotides, got {}",
            seed, flagged
        );
    }
}

// ──────────────────────────────────────────────────────────────────
// Property 5: Trace faithfulness — UniformMutationPass
// ──────────────────────────────────────────────────────────────────

#[test]
fn property_uniform_mutation_trace_faithfulness() {
    // For each mutation event the trace records (site, base). The
    // pool's base at that site, after the pass, must match the
    // last-recorded base for that site (later mutations to the
    // same site overwrite earlier ones).
    let count = 5usize;
    let mut plan = PassPlan::new();
    plan.push(Box::new(UniformMutationPass::new(
        Box::new(EmpiricalLengthDist::from_pairs(vec![(count as i64, 1.0)])),
        Box::new(UniformBase),
    )));

    for seed in 0..SEED_RANGE {
        let outcome = PassRuntime::execute(&plan, assembled_v_sim(), seed);
        let final_sim = outcome.final_simulation();

        // Walk events in order; later events at the same site
        // overwrite earlier ones, so the post-pass pool reflects
        // the *last* recorded base per site.
        let mut last: std::collections::HashMap<u32, u8> = std::collections::HashMap::new();
        for i in 0..count {
            let site = match outcome
                .trace
                .find(&format!("mutate.uniform.site[{}]", i))
                .unwrap()
                .value
            {
                ChoiceValue::Int(s) => s as u32,
                _ => unreachable!(),
            };
            let base = match outcome
                .trace
                .find(&format!("mutate.uniform.base[{}]", i))
                .unwrap()
                .value
            {
                ChoiceValue::Base(b) => b,
                _ => unreachable!(),
            };
            last.insert(site, base);
        }

        for (&site, &expected) in last.iter() {
            assert_eq!(
                final_sim.pool.get(NucHandle::new(site)).unwrap().base,
                expected,
                "seed {}: pool base at site {} doesn't match traced {} ",
                seed,
                site,
                expected as char
            );
        }
    }
}

// ──────────────────────────────────────────────────────────────────
// Property 6: Region range validity after every pass
// ──────────────────────────────────────────────────────────────────

#[test]
fn property_region_ranges_valid_after_indel() {
    let mut plan = PassPlan::new();
    plan.push(Box::new(IndelPass::new(
        Box::new(EmpiricalLengthDist::from_pairs(vec![(5, 1.0)])),
        0.5,
        Box::new(UniformBase),
    )));

    for seed in 0..SEED_RANGE {
        let outcome = PassRuntime::execute(&plan, assembled_v_sim(), seed);
        assert_region_ranges_valid(outcome.final_simulation(), &format!("indel seed={}", seed));
    }
}

#[test]
fn property_region_ranges_valid_after_full_vj_pipeline() {
    // After a full VJ recombination + corruption stack, every
    // region must still be in-bounds.
    let refdata = vj_refdata();
    let recomb = vj_plan(&refdata);

    for seed in 0..SEED_RANGE {
        // Recombination first.
        let recomb_outcome =
            PassRuntime::execute_with_refdata(&recomb, Simulation::default(), seed, &refdata);
        let assembled = recomb_outcome.final_simulation().clone();

        // Then corruption.
        let mut corrupt = PassPlan::new();
        corrupt.push(Box::new(S5FMutationPass::new(
            uniform_s5f(),
            Box::new(EmpiricalLengthDist::from_pairs(vec![(2, 1.0)])),
        )));
        corrupt.push(Box::new(IndelPass::new(
            Box::new(EmpiricalLengthDist::from_pairs(vec![(2, 1.0)])),
            0.5,
            Box::new(UniformBase),
        )));
        let corr_outcome = PassRuntime::execute(&corrupt, assembled, seed);

        assert_region_ranges_valid(
            corr_outcome.final_simulation(),
            &format!("vj+corruption seed={}", seed),
        );
        assert_frame_phases_consistent(
            corr_outcome.final_simulation(),
            &format!("vj+corruption seed={}", seed),
        );
        assert_codon_rails_consistent(
            corr_outcome.final_simulation(),
            &format!("vj+corruption seed={}", seed),
        );
    }
}

#[test]
fn property_frame_phases_valid_after_indel_in_vj_pipeline() {
    // Single-region indel tests cannot expose stale downstream phases.
    // Build a real VJ product first, then apply indels across the
    // multi-region sequence.
    let refdata = vj_refdata();
    let recomb = vj_plan(&refdata);

    let mut indel = PassPlan::new();
    indel.push(Box::new(IndelPass::new(
        Box::new(EmpiricalLengthDist::from_pairs(vec![(5, 1.0)])),
        0.5,
        Box::new(UniformBase),
    )));

    for seed in 0..SEED_RANGE {
        let recomb_outcome =
            PassRuntime::execute_with_refdata(&recomb, Simulation::default(), seed, &refdata);
        let outcome = PassRuntime::execute(
            &indel,
            recomb_outcome.final_simulation().clone(),
            seed ^ 0x9e37_79b9,
        );

        assert_frame_phases_consistent(
            outcome.final_simulation(),
            &format!("vj+indel seed={}", seed),
        );
        assert_codon_rails_consistent(
            outcome.final_simulation(),
            &format!("vj+indel seed={}", seed),
        );
    }
}

// ──────────────────────────────────────────────────────────────────
// Property 7: Determinism — same seed ⇒ identical trace + IR
// ──────────────────────────────────────────────────────────────────

#[test]
fn property_determinism_full_corruption_stack() {
    let plan = || {
        let mut p = PassPlan::new();
        p.push(Box::new(S5FMutationPass::new(
            uniform_s5f(),
            Box::new(EmpiricalLengthDist::from_pairs(vec![(3, 1.0)])),
        )));
        p.push(Box::new(PCRErrorPass::new(
            Box::new(EmpiricalLengthDist::from_pairs(vec![(2, 1.0)])),
            Box::new(UniformBase),
        )));
        p.push(Box::new(IndelPass::new(
            Box::new(EmpiricalLengthDist::from_pairs(vec![(2, 1.0)])),
            0.5,
            Box::new(UniformBase),
        )));
        p.push(Box::new(QualityErrorPass::new(
            Box::new(EmpiricalLengthDist::from_pairs(vec![(1, 1.0)])),
            Box::new(UniformBase),
        )));
        p.push(Box::new(ContaminantPass::new(0.2, Box::new(UniformBase))));
        p
    };

    for seed in 0..SEED_RANGE {
        let oa = PassRuntime::execute(&plan(), assembled_v_sim(), seed);
        let ob = PassRuntime::execute(&plan(), assembled_v_sim(), seed);

        // Trace identical.
        assert_eq!(
            oa.trace.choices(),
            ob.trace.choices(),
            "seed {} non-deterministic trace",
            seed
        );

        // Pool byte-identical.
        let pa = &oa.final_simulation().pool;
        let pb = &ob.final_simulation().pool;
        assert_eq!(pa.len(), pb.len(), "seed {} non-deterministic len", seed);
        for i in 0..pa.len() as u32 {
            let na = pa.get(NucHandle::new(i)).unwrap();
            let nb = pb.get(NucHandle::new(i)).unwrap();
            assert_eq!(
                na.base, nb.base,
                "seed {} non-deterministic base at {}",
                seed, i
            );
            assert_eq!(
                na.flags.bits(),
                nb.flags.bits(),
                "seed {} non-deterministic flags at {}",
                seed,
                i
            );
        }
    }
}

// ──────────────────────────────────────────────────────────────────
// Property 8: Productive admits ⇒ in-frame junction
// ──────────────────────────────────────────────────────────────────

#[test]
fn property_productive_implies_in_frame_vj() {
    // Under respect=[productive()], the assembled junction must be
    // in-frame on every seed. The D.7 test covers 50 seeds; this
    // is the same invariant at 100-seed sweep with a single
    // assertion failure surface.
    let refdata = vj_refdata();
    let plan = vj_plan(&refdata);
    let contracts = productive();

    for seed in 0..SEED_RANGE {
        let outcome = PassRuntime::execute_with_context(
            &plan,
            Simulation::default(),
            seed,
            Some(&refdata),
            Some(&contracts),
        );
        let junction = compute_junction(outcome.final_simulation(), &refdata)
            .expect("VJ junction should be defined");
        assert!(
            junction.is_in_frame(),
            "seed {} produced out-of-frame junction (length {})",
            seed,
            junction.length
        );
    }
}

// ──────────────────────────────────────────────────────────────────
// Property 9: Distribution support cardinality matches sampled values
// ──────────────────────────────────────────────────────────────────

#[test]
fn property_uniform_int_samples_lie_in_declared_support() {
    // A general distribution-trait property: for `EmpiricalLengthDist`
    // and `UniformBase`, every sampled value lies in the declared
    // support. This is the contract that `sample_filtered` /
    // constraint-aware sampling depends on.
    use genairr_engine::dist::UniformInt;
    use genairr_engine::rng::Rng;

    let dist = UniformInt::new(-3, 7); // [-3, 7)
    let support = dist.support().expect("UniformInt has finite support");
    let support_set: std::collections::HashSet<i64> = support.iter().map(|(v, _)| *v).collect();

    let mut rng = Rng::new(0xdead_beef);
    for _ in 0..2000 {
        let v = dist.sample(&mut rng);
        assert!(
            support_set.contains(&v),
            "sampled value {} outside declared support",
            v
        );
    }
}

#[test]
fn property_empirical_length_dist_samples_lie_in_declared_support() {
    use genairr_engine::rng::Rng;

    let dist = EmpiricalLengthDist::from_pairs(vec![(2, 1.0), (5, 2.0), (11, 0.5)]);
    let support = dist
        .support()
        .expect("EmpiricalLengthDist has finite support");
    let support_set: std::collections::HashSet<i64> = support.iter().map(|(v, _)| *v).collect();

    let mut rng = Rng::new(0xfeed_face);
    for _ in 0..2000 {
        let v = dist.sample(&mut rng);
        assert!(
            support_set.contains(&v),
            "sampled value {} outside declared support",
            v
        );
    }
}
