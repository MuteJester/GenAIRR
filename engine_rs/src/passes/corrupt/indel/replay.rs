//! Trace-injected replay for the indel pass (Tier 3 final).
//!
//! Replay does **not** run the mod-3 DP sampler. The cursor's
//! per-slot `kind[i]` / `site[i]` / `base[i]` records are the
//! proposal source. Each replayed event still passes:
//!
//! 1. **Shape validation** — wrong value-kind / missing
//!    insertion-base / address mismatch surfaces as
//!    `PassError::Replay`.
//! 2. **Bounds validation** — insertion site in
//!    `[0, current_pool_len]`, deletion site in
//!    `[0, current_pool_len)` (or the `-1` no-op sentinel),
//!    insertion base in `base_dist.support()` when enumerable.
//!    Mismatches surface as `InvalidPlanState`.
//! 3. **Frame compatibility** — under active contracts every
//!    event is classified via `ContractSet::admissible_indel_class_at`;
//!    `Forbidden` classes reject, and the cumulative mod-3 frame
//!    delta over all events must be 0 at end of tuple (the
//!    invariant the original DP enforced).
//! 4. **Final post-event check** — under active contracts the
//!    bundle's `admits_post_event` runs on the assembled
//!    hypothetical post-state. Failure surfaces as
//!    `ConstraintSampling`.
//!
//! The `kind=false, site=-1` slot is the only no-op shape replay
//! tolerates. It covers two original-run sources: an empty-pool
//! deletion (legacy no-contracts path) and the all-NoOp permissive
//! sentinel (contracts rejected every tuple). Both apply
//! identically — no pool mutation. Anywhere else, validation
//! failure errors rather than silently no-op'ing.

use crate::address::ChoiceAddress;
use crate::contract::{ChoiceContext, IndelEventClass, IndelKindHint};
use crate::dist::FilteredSampleError;
use crate::ir::NucHandle;
use crate::pass::{Pass, PassError};
use crate::passes::mutation_transaction::MutationTransaction;

use super::event::IndelEvent;
use super::tuple_sampler::apply_event_hypothetical;
use super::IndelPass;

/// Sentinel value the indel trace uses for both empty-pool
/// deletions and the all-NoOp permissive marker.
const NO_OP_SITE: i64 = -1;

impl IndelPass {
    /// Replay-mode tuple consumer. Called from
    /// `execute_with_sampling_mode` after the count is consumed and
    /// the transaction is open. Validates each event, runs the
    /// final post-event check under active contracts, then applies
    /// every event through the existing TX path.
    pub(super) fn replay_tuple(
        &self,
        tx: &mut MutationTransaction<'_, '_>,
        count: u32,
    ) -> Result<(), PassError> {
        // Step 1: consume all events from the cursor into an owned
        // Vec. The cursor borrow ends here; we hold the events for
        // the rest of the function.
        let mut events: Vec<IndelEvent> = Vec::with_capacity(count as usize);
        let pass_name = self.name();
        let support = self.base_dist.support();

        for i in 0..count {
            let kind_addr = ChoiceAddress::CorruptIndelKind(i);
            let site_addr = ChoiceAddress::CorruptIndelSite(i);
            let base_addr = ChoiceAddress::CorruptIndelBase(i);

            // Consume kind + site (always present).
            let (is_insertion, site_i64) = {
                let cursor = tx
                    .replay_cursor()
                    .expect("replay_cursor presence checked at top level");
                let kind = cursor
                    .expect_bool(kind_addr)
                    .map_err(|reason| PassError::replay(pass_name, reason))?;
                let site = cursor
                    .expect_int(site_addr)
                    .map_err(|reason| PassError::replay(pass_name, reason))?;
                (kind, site)
            };

            // Validate site value: deletions allow the `-1` no-op
            // sentinel; insertions never do.
            if is_insertion && site_i64 < 0 {
                return Err(PassError::invalid_plan_state(
                    pass_name.to_string(),
                    format!(
                        "indel: replayed insertion at slot {} has site {} (must be non-negative)",
                        i, site_i64
                    ),
                ));
            }

            if is_insertion {
                let base = {
                    let cursor = tx
                        .replay_cursor()
                        .expect("replay_cursor presence checked at top level");
                    cursor
                        .expect_base(base_addr)
                        .map_err(|reason| PassError::replay(pass_name, reason))?
                };

                // Insertion base must be in base_dist.support() if
                // the distribution exposes its support — otherwise
                // the live sampler could not have produced this
                // byte. Non-enumerable support is best-effort
                // accept.
                if let Some(ref supp) = support {
                    if !supp.iter().any(|(b, _)| *b == base) {
                        return Err(PassError::constraint_sampling(
                            pass_name,
                            base_addr.to_string(),
                            FilteredSampleError::EmptyAdmissibleSupport,
                        ));
                    }
                }

                events.push(IndelEvent::Insertion {
                    site: site_i64 as u32,
                    base,
                });
            } else if site_i64 == NO_OP_SITE {
                // Empty-pool deletion OR all-NoOp permissive sentinel.
                events.push(IndelEvent::Deletion { site: None });
            } else {
                events.push(IndelEvent::Deletion {
                    site: Some(site_i64 as u32),
                });
            }
        }

        // Step 2: validate against the in-progress sim. We walk a
        // hypothetical chain so each event sees the post-state of
        // its predecessors — the same shape `apply_event_hypothetical`
        // builds inside the DP sampler.
        let mut hypothetical = tx.peek().clone();
        let initial = hypothetical.clone();
        let mut frame_delta_total: i32 = 0;
        let contracts = tx.ctx().contracts;
        let refdata = tx.ctx().refdata;

        for (i, &event) in events.iter().enumerate() {
            // Per-event bounds against the hypothetical pool.
            let pool_len = hypothetical.pool.len() as u32;
            match event {
                IndelEvent::Insertion { site, .. } => {
                    if site > pool_len {
                        return Err(PassError::invalid_plan_state(
                            pass_name.to_string(),
                            format!(
                                "indel: replayed insertion at slot {} has site {} > pool length {}",
                                i, site, pool_len
                            ),
                        ));
                    }
                }
                IndelEvent::Deletion { site: Some(site) } => {
                    if site >= pool_len {
                        return Err(PassError::invalid_plan_state(
                            pass_name.to_string(),
                            format!(
                                "indel: replayed deletion at slot {} has site {} >= pool length {}",
                                i, site, pool_len
                            ),
                        ));
                    }
                }
                IndelEvent::Deletion { site: None } => {}
                IndelEvent::NoOp => {}
            }

            // Per-event frame classification under contracts.
            if let Some(contracts) = contracts {
                if let Some((site, kind_hint)) = indel_class_args(event) {
                    let class = contracts.admissible_indel_class_at(
                        &hypothetical,
                        refdata,
                        site,
                        kind_hint,
                    );
                    match class {
                        IndelEventClass::Forbidden => {
                            return Err(PassError::constraint_sampling(
                                pass_name,
                                ChoiceAddress::CorruptIndelSite(i as u32).to_string(),
                                FilteredSampleError::EmptyAdmissibleSupport,
                            ));
                        }
                        IndelEventClass::FrameDelta(d) => {
                            frame_delta_total += d as i32;
                        }
                        IndelEventClass::FrameNeutral => {}
                    }
                }
            }

            hypothetical = apply_event_hypothetical(&hypothetical, event);
        }

        if contracts.is_some() {
            // Tuple must preserve junction frame mod 3 — the DP
            // sampler in fresh mode only emits frame-balanced
            // tuples. Replay enforces the same invariant.
            if frame_delta_total.rem_euclid(3) != 0 {
                return Err(PassError::constraint_sampling(
                    pass_name,
                    ChoiceAddress::CorruptIndelCount.to_string(),
                    FilteredSampleError::EmptyAdmissibleSupport,
                ));
            }
        }

        // Step 3: final `admits_post_event` under active contracts.
        // The DP sampler runs this same check after building the
        // tuple; replay reproduces it on the assembled hypothetical
        // state.
        if let Some(contracts) = contracts {
            let final_index = count.saturating_sub(1);
            let validator_context =
                ChoiceContext::indel_insertion(final_index, count, NucHandle::new(0));
            if let Err(_violation) =
                contracts.admits_post_event(&initial, &hypothetical, refdata, validator_context)
            {
                return Err(PassError::constraint_sampling(
                    pass_name,
                    ChoiceAddress::CorruptIndelCount.to_string(),
                    FilteredSampleError::EmptyAdmissibleSupport,
                ));
            }
        }

        // Step 4: validation passed. Record + apply through the TX
        // exactly as fresh sampling does.
        for (i, event) in events.into_iter().enumerate() {
            let i = i as u32;
            Self::record_event(tx.ctx(), i, event);
            Self::apply_event_via_tx(tx, event)?;
        }
        Ok(())
    }
}

/// Return the `(site, kind_hint)` pair `admissible_indel_class_at`
/// queries against. Returns `None` for events with no site (no-op
/// shapes), which carry no class contribution.
fn indel_class_args(event: IndelEvent) -> Option<(u32, IndelKindHint)> {
    match event {
        IndelEvent::Insertion { site, .. } => Some((site, IndelKindHint::Insertion)),
        IndelEvent::Deletion { site: Some(site) } => Some((site, IndelKindHint::Deletion)),
        IndelEvent::Deletion { site: None } | IndelEvent::NoOp => None,
    }
}
