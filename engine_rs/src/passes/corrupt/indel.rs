//! `IndelPass` — insertions + deletions at the observation stage (E.7).

use crate::dist::Distribution;
use crate::ir::{flag, NucHandle, Nucleotide, Segment, Simulation};
use crate::pass::{Pass, PassContext};
use crate::trace::ChoiceValue;

/// Models PCR / sequencing indels: a small number of insertions
/// and deletions sprinkled across the assembled sequence.
///
/// **Architectural distinction from substitution passes:** indels
/// change pool length and shift downstream handle positions. The
/// underlying primitives (`Simulation::with_indel_inserted` /
/// `with_indel_deleted`) handle the position-shifting and
/// region-range adjustments. Codon rails refresh automatically
/// for any region whose range changed.
///
/// **Per-indel semantics:**
/// - Each indel is independently chosen to be an insertion or
///   deletion (50/50 by default — `insertion_prob` field gives
///   user control).
/// - The position is sampled uniformly within the current pool
///   (which means later indels may target positions that didn't
///   exist before earlier indels — that's fine, the pass walks
///   the pool's *current* state at each step).
/// - For insertions: a base is sampled from `base_dist`. Inserted
///   nucleotides are tagged with `flag::INDEL_INSERTED` and
///   carry no germline provenance (`germline_pos = NO_GERMLINE_POS`).
///
/// **Trace addresses (D3):**
/// - `corrupt.indel.count` — total indel events
/// - `corrupt.indel.kind[i]` — `Bool(true)` for insertion,
///   `Bool(false)` for deletion
/// - `corrupt.indel.site[i]` — pool position
/// - `corrupt.indel.base[i]` — only for insertions; the inserted base
///
/// **Edge cases handled:**
/// - Empty pool: deletion is a no-op for that step (insertion can
///   still happen at position 0).
/// - Deletion when pool length = 0: the kind is recorded but the
///   IR isn't modified (no-op).
pub struct IndelPass {
    count_dist: Box<dyn Distribution<Output = i64>>,
    insertion_prob: f64,
    base_dist: Box<dyn Distribution<Output = u8>>,
}

impl IndelPass {
    /// Construct an indel pass.
    ///
    /// Panics if `insertion_prob` is not in `[0.0, 1.0]` or is
    /// non-finite.
    pub fn new(
        count_dist: Box<dyn Distribution<Output = i64>>,
        insertion_prob: f64,
        base_dist: Box<dyn Distribution<Output = u8>>,
    ) -> Self {
        assert!(
            insertion_prob.is_finite() && (0.0..=1.0).contains(&insertion_prob),
            "IndelPass: insertion_prob must be in [0.0, 1.0], got {}",
            insertion_prob
        );
        Self {
            count_dist,
            insertion_prob,
            base_dist,
        }
    }

    /// Best-effort segment provenance for an inserted observation base.
    ///
    /// Insertions are synthetic, but downstream metadata still needs the
    /// nucleotide's segment to agree with the sequence context it lands in.
    /// This mirrors `Sequence::with_indel_adjusted`: if the insertion point
    /// is inside a region, or exactly at the start of a following region,
    /// the inserted nucleotide belongs to that region. Outside all regions,
    /// fall back to the nearest nucleotide's segment.
    fn insertion_segment(sim: &Simulation, at: u32) -> Segment {
        for region in &sim.sequence.regions {
            let start = region.start.index();
            let end = region.end.index();
            if start <= at && at < end {
                return region.segment;
            }
        }

        if at < sim.pool.len() as u32 {
            return sim
                .pool
                .get(NucHandle::new(at))
                .map(|n| n.segment)
                .unwrap_or(Segment::V);
        }

        if at > 0 {
            return sim
                .pool
                .get(NucHandle::new(at - 1))
                .map(|n| n.segment)
                .unwrap_or(Segment::V);
        }

        Segment::V
    }
}

impl Pass for IndelPass {
    fn name(&self) -> &str {
        "corrupt.indel"
    }

    fn execute(&self, sim: &Simulation, ctx: &mut PassContext) -> Simulation {
        let count_raw = self.count_dist.sample(ctx.rng);
        assert!(
            count_raw >= 0,
            "IndelPass: count distribution returned negative {}",
            count_raw
        );
        let count = count_raw as u32;
        ctx.trace
            .record("corrupt.indel.count", ChoiceValue::Int(count_raw));

        if count == 0 {
            return sim.clone();
        }

        let mut current = sim.clone();
        for i in 0..count {
            // Decide insertion vs deletion.
            let insertion = ctx.rng.next_f64() < self.insertion_prob;
            ctx.trace.record(
                format!("corrupt.indel.kind[{}]", i),
                ChoiceValue::Bool(insertion),
            );

            let pool_len = current.pool.len() as u32;

            if insertion {
                // Pool position can be in [0, pool_len] inclusive
                // (insertion at end is allowed).
                let site = if pool_len == 0 {
                    0
                } else {
                    ctx.rng.range_u32(pool_len + 1)
                };
                ctx.trace.record(
                    format!("corrupt.indel.site[{}]", i),
                    ChoiceValue::Int(site as i64),
                );

                let base = self.base_dist.sample(ctx.rng);
                ctx.trace.record(
                    format!("corrupt.indel.base[{}]", i),
                    ChoiceValue::Base(base),
                );

                let segment = Self::insertion_segment(&current, site);
                let new_nuc = Nucleotide::synthetic(base, segment, flag::INDEL_INSERTED);
                current = current.with_indel_inserted(site, new_nuc);
            } else {
                // Deletion: skip if pool empty.
                if pool_len == 0 {
                    // Record a sentinel position so the trace stays
                    // structurally consistent; downstream queries
                    // can detect the no-op via this value.
                    ctx.trace
                        .record(format!("corrupt.indel.site[{}]", i), ChoiceValue::Int(-1));
                    continue;
                }
                let site = ctx.rng.range_u32(pool_len);
                ctx.trace.record(
                    format!("corrupt.indel.site[{}]", i),
                    ChoiceValue::Int(site as i64),
                );
                current = current.with_indel_deleted(site);
            }
        }

        current
    }

    fn declared_choices(&self) -> Vec<String> {
        vec![
            "corrupt.indel.count".to_string(),
            "corrupt.indel.kind[0..n]".to_string(),
            "corrupt.indel.site[0..n]".to_string(),
            "corrupt.indel.base[0..n]".to_string(),
        ]
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::dist::{EmpiricalLengthDist, UniformBase};
    use crate::ir::{NucHandle, Region};
    use crate::pass::{PassPlan, PassRuntime};

    /// Helper: build a sim with a 12-base germline V region. Used as
    /// the canonical input for indel-pass tests so length/region
    /// invariants are easy to assert against the post-pass state.
    fn indel_test_sim() -> Simulation {
        let mut sim = Simulation::new();
        for (i, b) in b"AAACCCGGGTTT".iter().enumerate() {
            let (next, _) =
                sim.with_nucleotide_pushed(Nucleotide::germline(*b, i as u16, Segment::V));
            sim = next;
        }
        let region = Region::new(Segment::V, NucHandle::new(0), NucHandle::new(12))
            .with_codon_rail_recomputed(&sim.pool);
        sim.with_region_added(region)
    }

    /// Helper: build a contiguous multi-segment simulation so indel
    /// insertion provenance can be checked outside V.
    fn multi_segment_indel_context_sim() -> Simulation {
        let layout = [
            (Segment::V, b"AAA".as_slice()),
            (Segment::Np1, b"CC".as_slice()),
            (Segment::J, b"GGG".as_slice()),
        ];

        let mut sim = Simulation::new();
        let mut regions = Vec::new();
        for (segment, bases) in layout {
            let start = sim.pool.len() as u32;
            for (i, b) in bases.iter().enumerate() {
                let (next, _) =
                    sim.with_nucleotide_pushed(Nucleotide::germline(*b, i as u16, segment));
                sim = next;
            }
            let end = sim.pool.len() as u32;
            regions.push(Region::new(
                segment,
                NucHandle::new(start),
                NucHandle::new(end),
            ));
        }

        for region in regions {
            sim = sim.with_region_added(region.with_codon_rail_recomputed(&sim.pool));
        }
        sim
    }

    /// Helper: a count distribution that always returns the given
    /// value. Built on `EmpiricalLengthDist` for ergonomic test setup.
    fn fixed_count(n: i64) -> Box<dyn Distribution<Output = i64>> {
        Box::new(EmpiricalLengthDist::from_pairs(vec![(n, 1.0)]))
    }

    #[test]
    #[should_panic(expected = "insertion_prob must be in [0.0, 1.0]")]
    fn indel_pass_rejects_negative_insertion_prob() {
        let _ = IndelPass::new(fixed_count(0), -0.5, Box::new(UniformBase));
    }

    #[test]
    #[should_panic(expected = "insertion_prob must be in [0.0, 1.0]")]
    fn indel_pass_rejects_insertion_prob_above_one() {
        let _ = IndelPass::new(fixed_count(0), 1.5, Box::new(UniformBase));
    }

    #[test]
    #[should_panic(expected = "insertion_prob must be in [0.0, 1.0]")]
    fn indel_pass_rejects_nan_insertion_prob() {
        let _ = IndelPass::new(fixed_count(0), f64::NAN, Box::new(UniformBase));
    }

    #[test]
    fn indel_pass_zero_count_is_noop() {
        let mut plan = PassPlan::new();
        plan.push(Box::new(IndelPass::new(
            fixed_count(0),
            0.5,
            Box::new(UniformBase),
        )));
        let outcome = PassRuntime::execute(&plan, indel_test_sim(), 0);

        // Only the count is recorded; no per-indel addresses.
        assert_eq!(outcome.trace.len(), 1);
        assert_eq!(
            outcome.trace.find("corrupt.indel.count").unwrap().value,
            ChoiceValue::Int(0)
        );

        // Pool length unchanged.
        let final_sim = outcome.final_simulation();
        assert_eq!(final_sim.pool.len(), 12);
    }

    #[test]
    #[should_panic(expected = "count distribution returned negative")]
    fn indel_pass_negative_count_panics() {
        // crate::dist::UniformInt::new(-3, -2) always emits -3 —
        // exercises the negative-count guard.
        use crate::dist::UniformInt;
        let mut plan = PassPlan::new();
        plan.push(Box::new(IndelPass::new(
            Box::new(UniformInt::new(-3, -2)),
            0.5,
            Box::new(UniformBase),
        )));
        let _ = PassRuntime::execute(&plan, indel_test_sim(), 0);
    }

    #[test]
    fn indel_pass_insertion_prob_one_grows_pool() {
        // With p_ins = 1.0, every event is an insertion → pool grows
        // by exactly count.
        let mut plan = PassPlan::new();
        plan.push(Box::new(IndelPass::new(
            fixed_count(4),
            1.0,
            Box::new(UniformBase),
        )));
        let outcome = PassRuntime::execute(&plan, indel_test_sim(), 7);
        assert_eq!(outcome.final_simulation().pool.len(), 12 + 4);

        // Every kind[i] = Bool(true).
        for i in 0..4 {
            assert_eq!(
                outcome
                    .trace
                    .find(&format!("corrupt.indel.kind[{}]", i))
                    .unwrap()
                    .value,
                ChoiceValue::Bool(true)
            );
        }
    }

    #[test]
    fn indel_pass_insertion_prob_zero_shrinks_pool() {
        // With p_ins = 0.0, every event is a deletion → pool shrinks
        // by exactly count (provided pool stays non-empty).
        let mut plan = PassPlan::new();
        plan.push(Box::new(IndelPass::new(
            fixed_count(4),
            0.0,
            Box::new(UniformBase),
        )));
        let outcome = PassRuntime::execute(&plan, indel_test_sim(), 11);
        assert_eq!(outcome.final_simulation().pool.len(), 12 - 4);

        for i in 0..4 {
            assert_eq!(
                outcome
                    .trace
                    .find(&format!("corrupt.indel.kind[{}]", i))
                    .unwrap()
                    .value,
                ChoiceValue::Bool(false)
            );
        }
    }

    #[test]
    fn indel_pass_inserted_nucleotides_carry_indel_flag() {
        // Every nucleotide that landed via insertion must have the
        // INDEL_INSERTED flag set; existing germline nucleotides
        // must not.
        let mut plan = PassPlan::new();
        plan.push(Box::new(IndelPass::new(
            fixed_count(3),
            1.0,
            Box::new(UniformBase),
        )));
        let outcome = PassRuntime::execute(&plan, indel_test_sim(), 19);
        let final_sim = outcome.final_simulation();

        let mut inserted = 0;
        let mut germline = 0;
        for i in 0..final_sim.pool.len() as u32 {
            let n = final_sim.pool.get(NucHandle::new(i)).unwrap();
            if n.flags.contains(flag::INDEL_INSERTED) {
                inserted += 1;
                assert_eq!(n.germline_pos, Nucleotide::NO_GERMLINE_POS);
            } else {
                germline += 1;
            }
        }
        assert_eq!(inserted, 3);
        assert_eq!(germline, 12);
    }

    #[test]
    fn indel_pass_insertion_segment_uses_sequence_context() {
        let sim = multi_segment_indel_context_sim();

        assert_eq!(IndelPass::insertion_segment(&sim, 1), Segment::V);
        assert_eq!(IndelPass::insertion_segment(&sim, 3), Segment::Np1);
        assert_eq!(IndelPass::insertion_segment(&sim, 5), Segment::J);
        assert_eq!(IndelPass::insertion_segment(&sim, 8), Segment::J);
        assert_eq!(
            IndelPass::insertion_segment(&Simulation::new(), 0),
            Segment::V
        );
    }

    #[test]
    fn indel_pass_inserted_nucleotide_inherits_non_v_segment() {
        let sim = multi_segment_indel_context_sim();
        let mut plan = PassPlan::new();
        plan.push(Box::new(IndelPass::new(
            fixed_count(1),
            1.0,
            Box::new(UniformBase),
        )));

        for seed in 0..200u64 {
            let outcome = PassRuntime::execute(&plan, sim.clone(), seed);
            let site = match outcome.trace.find("corrupt.indel.site[0]").unwrap().value {
                ChoiceValue::Int(s) => s as u32,
                _ => unreachable!(),
            };

            if site >= 5 {
                let inserted = outcome
                    .final_simulation()
                    .pool
                    .get(NucHandle::new(site))
                    .unwrap();
                assert_eq!(inserted.segment, Segment::J);
                assert!(inserted.flags.contains(flag::INDEL_INSERTED));
                return;
            }
        }

        panic!("test setup did not produce an insertion in J context");
    }

    #[test]
    fn indel_pass_records_canonical_addresses() {
        // 1 count + 3 kind[i] + 3 site[i] + 3 base[i] = 10 records
        // (since p_ins=1.0 every event is an insertion and base[i]
        // is always recorded).
        let mut plan = PassPlan::new();
        plan.push(Box::new(IndelPass::new(
            fixed_count(3),
            1.0,
            Box::new(UniformBase),
        )));
        let outcome = PassRuntime::execute(&plan, indel_test_sim(), 5);
        assert_eq!(outcome.trace.len(), 10);
        for i in 0..3 {
            assert!(outcome
                .trace
                .find(&format!("corrupt.indel.kind[{}]", i))
                .is_some());
            assert!(outcome
                .trace
                .find(&format!("corrupt.indel.site[{}]", i))
                .is_some());
            assert!(outcome
                .trace
                .find(&format!("corrupt.indel.base[{}]", i))
                .is_some());
        }
    }

    #[test]
    fn indel_pass_deletion_path_does_not_record_base() {
        // p_ins = 0.0 → every event is a deletion → no base[i]
        // records are emitted.
        let mut plan = PassPlan::new();
        plan.push(Box::new(IndelPass::new(
            fixed_count(3),
            0.0,
            Box::new(UniformBase),
        )));
        let outcome = PassRuntime::execute(&plan, indel_test_sim(), 5);
        // 1 count + 3 kind + 3 site = 7. No base entries.
        assert_eq!(outcome.trace.len(), 7);
        for i in 0..3 {
            assert!(outcome
                .trace
                .find(&format!("corrupt.indel.base[{}]", i))
                .is_none());
        }
    }

    #[test]
    fn indel_pass_inserted_base_in_pool_matches_traced_base() {
        // After insertion, the freshly-inserted nucleotide at the
        // recorded site carries the recorded base. Subsequent indels
        // can shift earlier positions, so we replay the events
        // forward and track each insert's final position.
        let mut plan = PassPlan::new();
        plan.push(Box::new(IndelPass::new(
            fixed_count(2),
            1.0,
            Box::new(UniformBase),
        )));
        let outcome = PassRuntime::execute(&plan, indel_test_sim(), 42);
        let final_sim = outcome.final_simulation();

        // Replay the recorded events to know each insert's final
        // resting position after subsequent shifts.
        let mut sites: Vec<u32> = (0..2)
            .map(|i| {
                match outcome
                    .trace
                    .find(&format!("corrupt.indel.site[{}]", i))
                    .unwrap()
                    .value
                {
                    ChoiceValue::Int(s) => s as u32,
                    _ => unreachable!(),
                }
            })
            .collect();
        let bases: Vec<u8> = (0..2)
            .map(|i| {
                match outcome
                    .trace
                    .find(&format!("corrupt.indel.base[{}]", i))
                    .unwrap()
                    .value
                {
                    ChoiceValue::Base(b) => b,
                    _ => unreachable!(),
                }
            })
            .collect();

        // Each later insertion at a position ≤ an earlier site shifts
        // the earlier site up by 1. Walk the events left to right,
        // updating the resting positions of all previously-inserted
        // nucleotides.
        let mut resting: Vec<u32> = Vec::new();
        for i in 0..2 {
            let s = sites[i];
            for r in resting.iter_mut() {
                if *r >= s {
                    *r += 1;
                }
            }
            resting.push(s);
        }
        sites = resting;

        for (s, expected) in sites.iter().zip(bases.iter()) {
            let pool_base = final_sim.pool.get(NucHandle::new(*s)).unwrap().base;
            assert_eq!(pool_base, *expected, "trace lies at site {}", s);
        }
    }

    #[test]
    fn indel_pass_insertion_grows_region_only_when_spanning() {
        // The pre-pass region is [0, 12). An insertion at site < 12
        // spans the region and grows `end` by 1; an insertion at
        // site == pool_len lands strictly after the region and
        // leaves it unchanged. Replay the recorded events to compute
        // the expected final end and assert against it.
        let mut plan = PassPlan::new();
        plan.push(Box::new(IndelPass::new(
            fixed_count(3),
            1.0,
            Box::new(UniformBase),
        )));
        let outcome = PassRuntime::execute(&plan, indel_test_sim(), 0);
        let final_sim = outcome.final_simulation();

        // Replay: track the moving region end and shift it whenever
        // an insert site falls strictly inside the current region.
        let mut region_end: u32 = 12;
        for i in 0..3 {
            let site = match outcome
                .trace
                .find(&format!("corrupt.indel.site[{}]", i))
                .unwrap()
                .value
            {
                ChoiceValue::Int(s) => s as u32,
                _ => unreachable!(),
            };
            // Spanning rule: region.start (0) ≤ site < region.end.
            if site < region_end {
                region_end += 1;
            }
            // site == region_end is the "entirely before" case for
            // the region (region.end ≤ pos), so no shift.
        }

        assert_eq!(final_sim.sequence.regions.len(), 1);
        assert_eq!(final_sim.sequence.regions[0].start.index(), 0);
        assert_eq!(final_sim.sequence.regions[0].end.index(), region_end);
        assert_eq!(final_sim.pool.len(), 12 + 3);
    }

    #[test]
    fn indel_pass_codon_rail_refresh_after_insertion() {
        // Stored amino_acids on the spanning region must equal a
        // fresh recompute against the post-pass pool. Same invariant
        // proved for substitution passes; indel-aware version exercises
        // the rail-refresh path through `with_indel_inserted`.
        let mut plan = PassPlan::new();
        plan.push(Box::new(IndelPass::new(
            fixed_count(2),
            1.0,
            Box::new(UniformBase),
        )));
        let outcome = PassRuntime::execute(&plan, indel_test_sim(), 1);
        let final_sim = outcome.final_simulation();

        let stored = &final_sim.sequence.regions[0].amino_acids;
        let fresh = final_sim.sequence.regions[0].with_codon_rail_recomputed(&final_sim.pool);
        assert_eq!(stored, &fresh.amino_acids);
    }

    #[test]
    fn indel_pass_codon_rail_refresh_after_deletion() {
        // Same invariant for the deletion path.
        let mut plan = PassPlan::new();
        plan.push(Box::new(IndelPass::new(
            fixed_count(2),
            0.0,
            Box::new(UniformBase),
        )));
        let outcome = PassRuntime::execute(&plan, indel_test_sim(), 1);
        let final_sim = outcome.final_simulation();

        let stored = &final_sim.sequence.regions[0].amino_acids;
        let fresh = final_sim.sequence.regions[0].with_codon_rail_recomputed(&final_sim.pool);
        assert_eq!(stored, &fresh.amino_acids);
    }

    #[test]
    fn indel_pass_is_deterministic_under_same_seed() {
        let plan = || {
            let mut p = PassPlan::new();
            p.push(Box::new(IndelPass::new(
                fixed_count(3),
                0.5,
                Box::new(UniformBase),
            )));
            p
        };
        let oa = PassRuntime::execute(&plan(), indel_test_sim(), 0xc0ff_ee);
        let ob = PassRuntime::execute(&plan(), indel_test_sim(), 0xc0ff_ee);
        assert_eq!(oa.trace.choices(), ob.trace.choices());
        assert_eq!(
            oa.final_simulation().pool.len(),
            ob.final_simulation().pool.len()
        );
        for i in 0..oa.final_simulation().pool.len() as u32 {
            let h = NucHandle::new(i);
            assert_eq!(
                oa.final_simulation().pool.get(h).unwrap().base,
                ob.final_simulation().pool.get(h).unwrap().base
            );
        }
    }

    #[test]
    fn indel_pass_empty_pool_deletion_records_minus_one_site() {
        // Deletion against an empty pool: no IR change, but the trace
        // records site = -1 to keep address layout structurally
        // consistent with the non-empty case.
        let mut plan = PassPlan::new();
        plan.push(Box::new(IndelPass::new(
            fixed_count(2),
            0.0,
            Box::new(UniformBase),
        )));
        let outcome = PassRuntime::execute(&plan, Simulation::new(), 0);

        // Pool stays empty.
        assert_eq!(outcome.final_simulation().pool.len(), 0);

        // Both site[i] entries are -1.
        for i in 0..2 {
            assert_eq!(
                outcome
                    .trace
                    .find(&format!("corrupt.indel.site[{}]", i))
                    .unwrap()
                    .value,
                ChoiceValue::Int(-1)
            );
        }
    }

    #[test]
    fn indel_pass_half_probability_produces_mixed_kinds() {
        // p_ins = 0.5 should produce both insertions and deletions
        // across many seeds. Negative-control proof that the kind
        // coin actually fires.
        let plan = || {
            let mut p = PassPlan::new();
            p.push(Box::new(IndelPass::new(
                fixed_count(4),
                0.5,
                Box::new(UniformBase),
            )));
            p
        };

        let mut saw_insertion = false;
        let mut saw_deletion = false;
        for seed in 0..50u64 {
            let outcome = PassRuntime::execute(&plan(), indel_test_sim(), seed);
            for i in 0..4 {
                match outcome
                    .trace
                    .find(&format!("corrupt.indel.kind[{}]", i))
                    .unwrap()
                    .value
                {
                    ChoiceValue::Bool(true) => saw_insertion = true,
                    ChoiceValue::Bool(false) => saw_deletion = true,
                    _ => unreachable!(),
                }
            }
            if saw_insertion && saw_deletion {
                break;
            }
        }
        assert!(saw_insertion, "expected at least one insertion");
        assert!(saw_deletion, "expected at least one deletion");
    }

    #[test]
    fn indel_pass_persistent_ir_preserves_input() {
        // Persistent IR contract: input simulation isn't mutated.
        let sim = indel_test_sim();
        let original_len = sim.pool.len();
        let original_region_end = sim.sequence.regions[0].end.index();

        let mut plan = PassPlan::new();
        plan.push(Box::new(IndelPass::new(
            fixed_count(3),
            1.0,
            Box::new(UniformBase),
        )));
        let _outcome = PassRuntime::execute(&plan, sim.clone(), 99);

        assert_eq!(sim.pool.len(), original_len);
        assert_eq!(sim.sequence.regions[0].end.index(), original_region_end);
    }

    #[test]
    fn indel_pass_declared_choices() {
        let pass = IndelPass::new(fixed_count(0), 0.5, Box::new(UniformBase));
        let declared = pass.declared_choices();
        assert!(declared.contains(&"corrupt.indel.count".to_string()));
        assert!(declared.contains(&"corrupt.indel.kind[0..n]".to_string()));
        assert!(declared.contains(&"corrupt.indel.site[0..n]".to_string()));
        assert!(declared.contains(&"corrupt.indel.base[0..n]".to_string()));
    }
}
