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

mod edit_view;
mod event;
mod replay;
mod sampling;
mod tuple_sampler;
pub use tuple_sampler::IndelSamplingOptions;

pub struct IndelPass {
    count_source: CountSource,
    insertion_prob: f64,
    base_dist: Box<dyn Distribution<Output = u8>>,
    sampling_options: IndelSamplingOptions,
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
            sampling_options: IndelSamplingOptions::default(),
        }
    }

    #[must_use]
    pub fn with_sampling_options(mut self, options: IndelSamplingOptions) -> Self {
        self.sampling_options = options;
        self
    }

    fn record_sampling_report(
        ctx: &mut PassContext,
        count: u32,
        applied: u32,
        status: &str,
        report: &tuple_sampler::SamplingReport,
    ) {
        ctx.trace
            .record_sampling_diagnostic(crate::trace::SamplingDiagnostic {
                pass_name: address::CORRUPT_INDEL.to_string(),
                model: tuple_sampler::SAMPLER_MODEL.to_string(),
                strategy: report.strategy.to_string(),
                status: status.to_string(),
                requested_count: count,
                applied_count: applied,
                attempts: report.attempts,
                explored_states: report.explored_states,
                choice_offset: 0,
            });
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
            true,
        )?;
        if ctx.replay_cursor.is_some() {
            let supported =
                self.count_source
                    .supports_value(i64::from(count))
                    .map_err(|reason| {
                        PassError::constraint_sampling(
                            self.name(),
                            address::ChoiceAddress::CorruptIndelCount.to_string(),
                            reason,
                        )
                    })?;
            if !supported {
                return Err(PassError::invalid_distribution_output(
                    self.name(),
                    address::ChoiceAddress::CorruptIndelCount.to_string(),
                    i64::from(count),
                    "replay_count_outside_positive_support",
                ));
            }
        }

        if count == 0 {
            Self::record_sampling_report(
                ctx,
                0,
                0,
                "completed",
                &tuple_sampler::SamplingReport {
                    strategy: if ctx.replay_cursor.is_some() {
                        "replay"
                    } else {
                        "zero_count"
                    },
                    attempts: 0,
                    explored_states: 0,
                },
            );
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
        // Active contracts condition the natural sequential history on the
        // complete final state. Temporary, repairable changes are allowed.
        // Only a completed empty-support proof permits a permissive no-op;
        // work exhaustion remains an error in both execution modes.
        let mut tx = MutationTransaction::open(sim, ctx, self.name(), strict);

        // Trace-injected replay (Tier 3 — final): consume the
        // pre-recorded tuple from the cursor, validate each event
        // against current state + contracts, and apply through the
        // same `apply_event_via_tx` path as fresh sampling. No proposal
        // sampling, no RNG draws.
        if tx.replay_cursor().is_some() {
            self.replay_tuple(&mut tx, count)?;
            return tx.commit();
        }

        if tx.ctx().contracts.is_some() {
            let sampled = {
                let (sim_ref, ctx_ref) = tx.split_borrows();
                self.sample_admissible_tuple(sim_ref, ctx_ref, count)?
            };
            match sampled.events {
                Some(events) => {
                    let applied = events
                        .iter()
                        .filter(|event| {
                            matches!(
                                event,
                                event::IndelEvent::Insertion { .. }
                                    | event::IndelEvent::Deletion { site: Some(_) }
                            )
                        })
                        .count() as u32;
                    for (i, event) in events.into_iter().enumerate() {
                        let i = i as u32;
                        Self::record_event(tx.ctx(), i, event);
                        Self::apply_event_via_tx(&mut tx, event)?;
                    }
                    Self::record_sampling_report(
                        tx.ctx(),
                        count,
                        applied,
                        "completed",
                        &sampled.report,
                    );
                }
                None => {
                    if strict {
                        return Err(PassError::constraint_sampling(
                            self.name(),
                            address::ChoiceAddress::CorruptIndelCount.to_string(),
                            crate::dist::FilteredSampleError::EmptyAdmissibleSupport,
                        ));
                    }
                    // Permissive: no admissible tuple. Record NoOp
                    // for each slot so trace consumers can audit
                    // the reduction; no pool mutation.
                    for i in 0..count {
                        Self::record_event(tx.ctx(), i, event::IndelEvent::NoOp);
                    }
                    Self::record_sampling_report(
                        tx.ctx(),
                        count,
                        0,
                        "proven_empty_support",
                        &sampled.report,
                    );
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
        let mut applied = 0;
        for i in 0..count {
            let event = {
                let (sim_ref, ctx_ref) = tx.split_borrows();
                self.sample_legacy_event(sim_ref, ctx_ref)
            };
            if matches!(
                event,
                event::IndelEvent::Insertion { .. } | event::IndelEvent::Deletion { site: Some(_) }
            ) {
                applied += 1;
            }
            Self::record_event(tx.ctx(), i, event);
            Self::apply_event_via_tx(&mut tx, event)?;
        }
        Self::record_sampling_report(
            tx.ctx(),
            count,
            applied,
            "completed",
            &tuple_sampler::SamplingReport {
                strategy: "unconstrained_sequential",
                attempts: 1,
                explored_states: 0,
            },
        );

        tx.commit()
    }
}

impl Pass for IndelPass {
    fn name(&self) -> &str {
        address::CORRUPT_INDEL
    }

    fn parameter_signature(&self) -> String {
        use crate::passes::paramsig::{fmt_byte_dist, fmt_count_source, fmt_prob, join_parts};
        join_parts([
            format!(
                "sampler={},search=bounded_interleaved_v1,rejection_attempts={},exact_node_budget={},legacy_replay_sentinels={}",
                tuple_sampler::SAMPLER_MODEL,
                self.sampling_options.rejection_attempts,
                self.sampling_options.exact_node_budget,
                self.sampling_options.legacy_replay_sentinels
            ),
            fmt_count_source(&self.count_source),
            fmt_prob("insertion_prob", self.insertion_prob),
            format!("base={}", fmt_byte_dist(self.base_dist.as_ref())),
        ])
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

    fn execute_with_policy(
        &self,
        sim: &Simulation,
        ctx: &mut PassContext,
        strict: bool,
    ) -> Result<Simulation, PassError> {
        self.execute_with_sampling_mode(sim, ctx, strict)
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
