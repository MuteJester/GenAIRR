//! `IndelPass` — insertions + deletions at the observation stage (E.7).
//!
//! Models PCR / sequencing indels: a small number of insertions
//! and deletions sprinkled across the assembled sequence.
//!
//! **Architectural distinction from substitution passes:** indels
//! change pool length and shift downstream handle positions. The
//! underlying primitives (`Simulation::with_indel_inserted` /
//! `with_indel_deleted`) handle the position-shifting and
//! region-range adjustments. Codon rails refresh automatically
//! for any region whose range changed.
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
use crate::passes::mutation_transaction::MutationTransaction;
use crate::trace::ChoiceValue;

mod event;
mod sampling;

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
                address::CORRUPT_INDEL_COUNT,
                count_raw,
                "negative_count",
            ));
        }
        if strict && count_raw > u32::MAX as i64 {
            return Err(PassError::invalid_distribution_output(
                self.name(),
                address::CORRUPT_INDEL_COUNT,
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
            .record(address::CORRUPT_INDEL_COUNT, ChoiceValue::Int(count_raw));

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
        let mut tx = MutationTransaction::open(sim, ctx, self.name(), strict);

        for i in 0..count {
            let event = {
                let (sim_ref, ctx_ref) = tx.split_borrows();
                self.sample_event(sim_ref, ctx_ref, i, count, strict)?
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

    fn declared_choices(&self) -> Vec<String> {
        vec![
            address::CORRUPT_INDEL_COUNT.to_string(),
            address::CORRUPT_INDEL_KIND_PATTERN.to_string(),
            address::CORRUPT_INDEL_SITE_PATTERN.to_string(),
            address::CORRUPT_INDEL_BASE_PATTERN.to_string(),
        ]
    }

    fn effects(&self) -> Vec<PassEffect> {
        vec![PassEffect::StructuralIndel]
    }
}

#[cfg(test)]
mod tests;
