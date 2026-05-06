//! `IndelPass` — insertions + deletions at the observation stage (E.7).

use crate::contract::ChoiceContext;
use crate::dist::{Distribution, FilteredSampleError};
use crate::ir::{flag, NucHandle, Nucleotide, Segment, Simulation};
use crate::pass::{Pass, PassContext, PassEffect, PassError};
use crate::passes::constrained::{sample_contract_verified_event, PostEventCandidate};
use crate::trace::ChoiceValue;

#[derive(Clone, Copy, Debug, Eq, PartialEq)]
enum IndelEvent {
    Insertion { site: u32, base: u8 },
    Deletion { site: Option<u32> },
}

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
///
/// **Constraint awareness:** with active contracts, kind/site/base are
/// treated as one finite structural event. The pass builds each
/// candidate's hypothetical post-event `Simulation` and samples only
/// from candidates whose post-state verifies against the active
/// contracts. No-contract execution stays on the legacy RNG path.
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

    fn base_support(&self) -> Result<Vec<(u8, f64)>, FilteredSampleError> {
        let support = self
            .base_dist
            .support()
            .ok_or(FilteredSampleError::SupportUnavailable)?;
        let total: f64 = support.iter().map(|(_, weight)| *weight).sum();
        if support.is_empty() || !total.is_finite() || total <= 0.0 {
            return Err(FilteredSampleError::InvalidFilteredSupport);
        }
        if support
            .iter()
            .any(|(_, weight)| !weight.is_finite() || *weight <= 0.0)
        {
            return Err(FilteredSampleError::InvalidFilteredSupport);
        }
        Ok(support
            .into_iter()
            .map(|(base, weight)| (base, weight / total))
            .collect())
    }

    fn constrained_candidates(
        &self,
        sim: &Simulation,
        index: u32,
        count: u32,
    ) -> Result<Vec<PostEventCandidate<IndelEvent>>, FilteredSampleError> {
        let pool_len = sim.pool.len() as u32;
        let mut candidates = Vec::new();

        if self.insertion_prob > 0.0 {
            let base_support = self.base_support()?;
            let site_weight = self.insertion_prob / (pool_len + 1) as f64;

            for site in 0..=pool_len {
                let segment = Self::insertion_segment(sim, site);
                for &(base, base_weight) in &base_support {
                    let nuc = Nucleotide::synthetic(base, segment, flag::INDEL_INSERTED);
                    let post_sim = sim.with_indel_inserted(site, nuc);
                    candidates.push(PostEventCandidate::new(
                        IndelEvent::Insertion { site, base },
                        site_weight * base_weight,
                        post_sim,
                        ChoiceContext::indel_insertion(index, count, NucHandle::new(site)),
                    ));
                }
            }
        }

        let deletion_prob = 1.0 - self.insertion_prob;
        if deletion_prob > 0.0 {
            if pool_len == 0 {
                candidates.push(PostEventCandidate::new(
                    IndelEvent::Deletion { site: None },
                    deletion_prob,
                    sim.clone(),
                    ChoiceContext::indel_deletion_noop(index, count),
                ));
            } else {
                let site_weight = deletion_prob / pool_len as f64;
                for site in 0..pool_len {
                    let post_sim = sim.with_indel_deleted(site);
                    candidates.push(PostEventCandidate::new(
                        IndelEvent::Deletion { site: Some(site) },
                        site_weight,
                        post_sim,
                        ChoiceContext::indel_deletion(index, count, NucHandle::new(site)),
                    ));
                }
            }
        }

        Ok(candidates)
    }

    fn sample_legacy_event(&self, sim: &Simulation, ctx: &mut PassContext) -> IndelEvent {
        let insertion = ctx.rng.next_f64() < self.insertion_prob;
        let pool_len = sim.pool.len() as u32;

        if insertion {
            let site = if pool_len == 0 {
                0
            } else {
                ctx.rng.range_u32(pool_len + 1)
            };
            let base = self.base_dist.sample(ctx.rng);
            IndelEvent::Insertion { site, base }
        } else if pool_len == 0 {
            IndelEvent::Deletion { site: None }
        } else {
            IndelEvent::Deletion {
                site: Some(ctx.rng.range_u32(pool_len)),
            }
        }
    }

    fn sample_event(
        &self,
        sim: &Simulation,
        ctx: &mut PassContext,
        index: u32,
        count: u32,
        strict: bool,
    ) -> Result<IndelEvent, PassError> {
        if ctx.contracts.is_some() {
            let event_address = format!("corrupt.indel.site[{}]", index);
            let candidates = self.constrained_candidates(sim, index, count);
            if let Some(event) = sample_contract_verified_event(
                sim,
                ctx,
                self.name(),
                &event_address,
                strict,
                candidates,
            )? {
                return Ok(event);
            }
        }

        Ok(self.sample_legacy_event(sim, ctx))
    }

    fn record_event(ctx: &mut PassContext, index: u32, event: IndelEvent) {
        match event {
            IndelEvent::Insertion { site, base } => {
                ctx.trace.record(
                    format!("corrupt.indel.kind[{}]", index),
                    ChoiceValue::Bool(true),
                );
                ctx.trace.record(
                    format!("corrupt.indel.site[{}]", index),
                    ChoiceValue::Int(site as i64),
                );
                ctx.trace.record(
                    format!("corrupt.indel.base[{}]", index),
                    ChoiceValue::Base(base),
                );
            }
            IndelEvent::Deletion { site } => {
                ctx.trace.record(
                    format!("corrupt.indel.kind[{}]", index),
                    ChoiceValue::Bool(false),
                );
                ctx.trace.record(
                    format!("corrupt.indel.site[{}]", index),
                    ChoiceValue::Int(site.map(i64::from).unwrap_or(-1)),
                );
            }
        }
    }

    fn apply_event(sim: &Simulation, event: IndelEvent) -> Simulation {
        match event {
            IndelEvent::Insertion { site, base } => {
                let segment = Self::insertion_segment(sim, site);
                let new_nuc = Nucleotide::synthetic(base, segment, flag::INDEL_INSERTED);
                sim.with_indel_inserted(site, new_nuc)
            }
            IndelEvent::Deletion { site: Some(site) } => sim.with_indel_deleted(site),
            IndelEvent::Deletion { site: None } => sim.clone(),
        }
    }

    fn execute_with_sampling_mode(
        &self,
        sim: &Simulation,
        ctx: &mut PassContext,
        strict: bool,
    ) -> Result<Simulation, PassError> {
        let count_raw = self.count_dist.sample(ctx.rng);
        if strict && count_raw < 0 {
            return Err(PassError::invalid_distribution_output(
                self.name(),
                "corrupt.indel.count",
                count_raw,
                "negative_count",
            ));
        }
        if strict && count_raw > u32::MAX as i64 {
            return Err(PassError::invalid_distribution_output(
                self.name(),
                "corrupt.indel.count",
                count_raw,
                "count_exceeds_u32",
            ));
        }
        assert!(
            count_raw >= 0,
            "IndelPass: count distribution returned negative {}",
            count_raw
        );
        assert!(
            count_raw <= u32::MAX as i64,
            "IndelPass: count distribution returned {} > u32::MAX",
            count_raw
        );
        let count = count_raw as u32;
        ctx.trace
            .record("corrupt.indel.count", ChoiceValue::Int(count_raw));

        if count == 0 {
            return Ok(sim.clone());
        }

        let mut current = sim.clone();
        for i in 0..count {
            let event = self.sample_event(&current, ctx, i, count, strict)?;
            Self::record_event(ctx, i, event);
            current = Self::apply_event(&current, event);
        }

        Ok(current)
    }
}

impl Pass for IndelPass {
    fn name(&self) -> &str {
        "corrupt.indel"
    }

    fn execute(&self, sim: &Simulation, ctx: &mut PassContext) -> Simulation {
        self.execute_with_sampling_mode(sim, ctx, false)
            .expect("IndelPass permissive execution must not return PassError")
    }

    fn execute_checked(
        &self,
        sim: &Simulation,
        ctx: &mut PassContext,
    ) -> Result<Simulation, PassError> {
        self.execute_with_sampling_mode(sim, ctx, true)
    }

    fn declared_choices(&self) -> Vec<String> {
        vec![
            "corrupt.indel.count".to_string(),
            "corrupt.indel.kind[0..n]".to_string(),
            "corrupt.indel.site[0..n]".to_string(),
            "corrupt.indel.base[0..n]".to_string(),
        ]
    }

    fn effects(&self) -> Vec<PassEffect> {
        vec![PassEffect::StructuralIndel]
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::contract::productive;
    use crate::dist::{EmpiricalLengthDist, FilteredSampleError, UniformBase};
    use crate::ir::{NucHandle, Region};
    use crate::pass::{PassError, PassPlan, PassRuntime};
    use crate::passes::test_support::make_substitution_productive_vj_fixture;

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
    fn indel_pass_strict_errors_on_negative_count() {
        use crate::dist::UniformInt;
        let mut plan = PassPlan::new();
        plan.push(Box::new(IndelPass::new(
            Box::new(UniformInt::new(-3, -2)),
            0.5,
            Box::new(UniformBase),
        )));

        let err = PassRuntime::execute_strict_with_context(&plan, indel_test_sim(), 0, None, None)
            .unwrap_err();

        assert_eq!(err.pass_name(), "corrupt.indel");
        assert_eq!(err.address(), "corrupt.indel.count");
        assert!(matches!(
            err,
            PassError::InvalidDistributionOutput { value: -3, .. }
        ));
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

    #[test]
    fn indel_pass_productive_filters_insertion_to_structurally_safe_site() {
        let (cfg, sim) = make_substitution_productive_vj_fixture();
        let contracts = productive();

        let mut plan = PassPlan::new();
        plan.push(Box::new(IndelPass::new(
            fixed_count(1),
            1.0,
            Box::new(UniformBase),
        )));

        let outcome =
            PassRuntime::execute_with_context(&plan, sim, 0, Some(&cfg), Some(&contracts));

        assert_eq!(
            outcome.trace.find("corrupt.indel.kind[0]").unwrap().value,
            ChoiceValue::Bool(true)
        );
        assert_eq!(
            outcome.trace.find("corrupt.indel.site[0]").unwrap().value,
            ChoiceValue::Int(6)
        );
        assert!(contracts
            .verify(outcome.final_simulation(), Some(&cfg))
            .is_ok());
    }

    #[test]
    fn indel_pass_permissive_falls_back_when_deletions_cannot_satisfy_contracts() {
        let (cfg, sim) = make_substitution_productive_vj_fixture();
        let contracts = productive();

        let mut plan = PassPlan::new();
        plan.push(Box::new(IndelPass::new(
            fixed_count(1),
            0.0,
            Box::new(UniformBase),
        )));

        let outcome =
            PassRuntime::execute_with_context(&plan, sim, 0, Some(&cfg), Some(&contracts));

        assert_eq!(
            outcome.trace.find("corrupt.indel.kind[0]").unwrap().value,
            ChoiceValue::Bool(false)
        );
        let violations = contracts
            .verify(outcome.final_simulation(), Some(&cfg))
            .unwrap_err();
        assert!(
            violations.iter().any(|v| {
                v.contract_name == "productive_junction_frame"
                    || v.contract_name == "anchor_preserved.v"
                    || v.contract_name == "anchor_preserved.j"
            }),
            "expected productive-bundle structural violation, got {:?}",
            violations
        );
    }

    #[test]
    fn indel_pass_strict_errors_when_deletion_filter_empty() {
        let (cfg, sim) = make_substitution_productive_vj_fixture();
        let contracts = productive();

        let mut plan = PassPlan::new();
        plan.push(Box::new(IndelPass::new(
            fixed_count(1),
            0.0,
            Box::new(UniformBase),
        )));

        let err =
            PassRuntime::execute_strict_with_context(&plan, sim, 0, Some(&cfg), Some(&contracts))
                .unwrap_err();

        assert_eq!(err.pass_name(), "corrupt.indel");
        assert_eq!(err.address(), "corrupt.indel.site[0]");
        assert_eq!(
            err.constraint_reason(),
            Some(FilteredSampleError::EmptyAdmissibleSupport)
        );
    }
}
