//! `IndelPass` — insertions + deletions at the observation stage (E.7).
//!
//! Models PCR / sequencing indels: a small number of insertions
//! and deletions sprinkled across the assembled sequence.
//!
//! **Architectural distinction from substitution passes:** indels
//! change pool length and shift downstream handle positions. The
//! underlying primitives (`Simulation::with_indel_inserted` /
//! `with_indel_deleted`) handle the position-shifting and
//! region-range adjustments. Codon-rail data is not stored on
//! `Region`; callers that need amino-acid translation call
//! `compute_codon_rail(region, pool)` against the post-indel pool.
//!
//! **Per-indel semantics:**
//! - Each indel is independently chosen to be an insertion or
//!   deletion (50/50 by default — `insertion_prob` field gives
//!   user control).
//! - The position is sampled uniformly within the current pool
//!   (which means later indels may target positions that didn't
//!   exist before earlier indels — that's fine, the pass walks
//!   the pool's *current* state at each step).
//! - For insertions: a base is sampled from `base_dist`. Inserted
//!   nucleotides are tagged with `flag::INDEL_INSERTED` and
//!   carry no germline provenance (`germline_pos = GermlinePos::NONE`).
//!
//! **Trace addresses (D3):**
//! - `corrupt.indel.count` — total indel events
//! - `corrupt.indel.kind[i]` — `Bool(true)` for insertion,
//!   `Bool(false)` for deletion
//! - `corrupt.indel.site[i]` — pool position
//! - `corrupt.indel.base[i]` — only for insertions; the inserted base
//!
//! **Edge cases handled:**
//! - Empty pool: deletion is a no-op for that step (insertion can
//!   still happen at position 0).
//! - Deletion when pool length = 0: the kind is recorded but the
//!   IR isn't modified (no-op).
//!
//! **Constraint awareness:** with active contracts, kind/site/base are
//! treated as one finite structural event. The pass builds each
//! candidate's hypothetical post-event `Simulation` and samples only
//! from candidates whose post-state verifies against the active
//! contracts. No-contract execution stays on the legacy RNG path.

use crate::address;
use crate::dist::Distribution;
use crate::ir::Simulation;
use crate::pass::{Pass, PassContext, PassEffect, PassError};
use crate::passes::count_source::{sample_validated_count, CountSource};
use crate::passes::mutation_transaction::MutationTransaction;

mod event;
mod replay;
mod sampling;
mod tuple_sampler;

pub struct IndelPass {
    count_source: CountSource,
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
            count_source: CountSource::Distribution(count_dist),
            insertion_prob,
            base_dist,
        }
    }

    fn execute_with_sampling_mode(
        &self,
        sim: &Simulation,
        ctx: &mut PassContext,
        strict: bool,
    ) -> Result<Simulation, PassError> {
        // Count validation + trace recording shared with the other
        // four count-driven passes via `sample_validated_count`. The
        // helper carries trace-injected replay for the count slot
        // automatically; the per-event tuple replay is handled
        // below.
        let count = sample_validated_count(
            &self.count_source,
            ctx,
            sim.pool.len() as u32,
            self.name(),
            address::ChoiceAddress::CorruptIndelCount,
            strict,
        )?;

        if count == 0 {
            return Ok(sim.clone());
        }

        // Route indel mutations through `MutationTransaction` so each
        // insertion/deletion emits an `on_indel_*` event to attached
        // observers, and the seal step picks the reference-index-aware
        // path automatically.
        //
        // Codon-rail and walker observers implement `on_indel_*`
        // properly. External indels (before a region) are absorbed as
        // `seq_start` shifts; internal indels mark the observer
        // `needs_rebuild` so the seal helpers replace it with a fresh
        // `from_existing_region` rebuild against the post-mutation
        // pool. The staged live calls then absorb via the Phase-1
        // fast path in the post-pass `PassEffect::StructuralIndel`
        // dispatch, replacing the prior unconditional from-scratch
        // `call_from_region` re-walk per V/D/J segment.
        // v3.0 architecture: under active contracts, single indels
        // around the junction break frame and would all be rejected
        // by the post-event check. The
        // `sample_admissible_tuple` path picks a frame-balanced
        // tuple in one shot using a mod-3 reachability DP +
        // continuation weighting, then validates the assembled
        // post-state against the full contract bundle. On
        // permissive-mode rejection (DP exhausted, or final
        // validator fails) every slot records as NoOp; in strict
        // mode the failure surfaces as `PassError::ConstraintSampling`
        // / `ContractViolation`. When no contracts are active we
        // stay on the per-step legacy path so RNG behavior matches
        // the pre-v3 indel pass exactly.
        let mut tx = MutationTransaction::open(sim, ctx, self.name(), strict);

        // Trace-injected replay (Tier 3 — final): consume the
        // pre-recorded tuple from the cursor, validate each event
        // against current state + contracts, and apply through the
        // same `apply_event_via_tx` path as fresh sampling. No DP
        // sampling, no RNG draws.
        if tx.replay_cursor().is_some() {
            self.replay_tuple(&mut tx, count)?;
            return tx.commit();
        }

        if tx.ctx().contracts.is_some() {
            let tuple_opt = {
                let (sim_ref, ctx_ref) = tx.split_borrows();
                self.sample_admissible_tuple(sim_ref, ctx_ref, count, strict)?
            };
            match tuple_opt {
                Some(events) => {
                    for (i, event) in events.into_iter().enumerate() {
                        let i = i as u32;
                        Self::record_event(tx.ctx(), i, event);
                        Self::apply_event_via_tx(&mut tx, event)?;
                    }
                }
                None => {
                    // Permissive: no admissible tuple. Record NoOp
                    // for each slot so trace consumers can audit
                    // the reduction; no pool mutation.
                    for i in 0..count {
                        Self::record_event(tx.ctx(), i, event::IndelEvent::NoOp);
                    }
                }
            }
            return tx.commit();
        }

        // Legacy path (no contracts): per-step independent
        // sampling via the unconstrained natural distribution.
        // When contracts are active the path above (tuple sampler)
        // runs instead, so this branch is reached only with
        // `ctx.contracts == None`.
        let _ = strict;
        for i in 0..count {
            let event = {
                let (sim_ref, ctx_ref) = tx.split_borrows();
                self.sample_legacy_event(sim_ref, ctx_ref)
            };
            Self::record_event(tx.ctx(), i, event);
            Self::apply_event_via_tx(&mut tx, event)?;
        }

        tx.commit()
    }
}

impl Pass for IndelPass {
    fn name(&self) -> &str {
        address::CORRUPT_INDEL
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

    fn declared_choice_patterns(&self) -> Vec<address::ChoiceAddressPattern> {
        vec![
            address::ChoiceAddressPattern::CorruptIndelCount,
            address::ChoiceAddressPattern::CorruptIndelKind,
            address::ChoiceAddressPattern::CorruptIndelSite,
            address::ChoiceAddressPattern::CorruptIndelBase,
        ]
    }

    fn effects(&self) -> Vec<PassEffect> {
        vec![PassEffect::StructuralIndel]
    }
}

#[cfg(test)]
mod tests;
