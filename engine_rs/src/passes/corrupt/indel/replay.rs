//! Replay uses the same positive natural support and final tuple predicate as
//! fresh sampling. Intermediate frameshifts and anchor insertions may recover;
//! replay never draws replacement proposals.
use crate::address::ChoiceAddress;
use crate::dist::FilteredSampleError;
use crate::pass::{Pass, PassError};
use crate::passes::mutation_transaction::MutationTransaction;

use super::event::IndelEvent;
use super::tuple_sampler::{apply_event_hypothetical, tuple_admitted, SamplingReport};
use super::IndelPass;

impl IndelPass {
    pub(super) fn replay_tuple(
        &self,
        tx: &mut MutationTransaction<'_, '_>,
        count: u32,
    ) -> Result<(), PassError> {
        let mut events = Vec::with_capacity((count as usize).min(4096));
        let pass = self.name();
        // Validate positive finite support once, before any tuple can commit.
        let bases = if self.insertion_prob > 0.0 {
            Some(self.base_log_support().map_err(|reason| {
                PassError::constraint_sampling(
                    pass,
                    ChoiceAddress::CorruptIndelCount.to_string(),
                    reason,
                )
            })?)
        } else {
            None
        };
        for i in 0..count {
            let cursor = tx.replay_cursor().expect("replay cursor checked by caller");
            let insertion = cursor
                .expect_bool(ChoiceAddress::CorruptIndelKind(i))
                .map_err(|e| PassError::replay(pass, e))?;
            let site = cursor
                .expect_int(ChoiceAddress::CorruptIndelSite(i))
                .map_err(|e| PassError::replay(pass, e))?;
            if insertion {
                let base = cursor
                    .expect_base(ChoiceAddress::CorruptIndelBase(i))
                    .map_err(|e| PassError::replay(pass, e))?;
                let site =
                    u32::try_from(site).map_err(|_| {
                        PassError::invalid_plan_state(pass,
                    format!("replayed insertion site must be non-negative and fit u32, got {site}"))
                    })?;
                if self.insertion_prob == 0.0 {
                    return Err(PassError::invalid_plan_state(
                        pass,
                        "replay insertion has zero natural probability",
                    ));
                }
                if !bases.as_ref().unwrap().iter().any(|(b, _)| *b == base) {
                    return Err(PassError::constraint_sampling(
                        pass,
                        ChoiceAddress::CorruptIndelBase(i).to_string(),
                        FilteredSampleError::EmptyAdmissibleSupport,
                    ));
                }
                events.push(IndelEvent::Insertion { site, base });
            } else if site == -1 {
                events.push(IndelEvent::Deletion { site: None });
            } else {
                let site = u32::try_from(site).map_err(|_| {
                    PassError::invalid_plan_state(
                        pass,
                        format!("replayed deletion site must be -1 or fit u32, got {site}"),
                    )
                })?;
                if self.insertion_prob == 1.0 {
                    return Err(PassError::invalid_plan_state(
                        pass,
                        "replay deletion has zero natural probability",
                    ));
                }
                events.push(IndelEvent::Deletion { site: Some(site) });
            }
        }
        let initial = tx.peek().clone();
        let all_noop = events
            .iter()
            .all(|e| matches!(e, IndelEvent::Deletion { site: None }));
        let natural_empty_deletions = initial.pool.is_empty() && self.insertion_prob < 1.0;
        let mut status = "completed";
        let mut work = 0;
        let compatibility_skip = if all_noop && !natural_empty_deletions {
            if self.sampling_options.legacy_replay_sentinels {
                status = "legacy_skipped_batch";
                true
            } else if let Some(contracts) = tx.ctx().contracts {
                let (empty, states) =
                    self.prove_empty_support(&initial, contracts, tx.ctx().refdata, count)?;
                work = states;
                if !empty {
                    return Err(PassError::invalid_plan_state(pass,
                        "all-no-op replay is neither an empty-pool deletion nor a proven-empty batch"));
                }
                status = "proven_empty_support";
                true
            } else {
                return Err(PassError::invalid_plan_state(
                    pass,
                    "no-op replay requires explicit legacy compatibility",
                ));
            }
        } else {
            false
        };

        let mut hypothetical = initial.clone();
        if !compatibility_skip {
            for (i, &event) in events.iter().enumerate() {
                let len = hypothetical.pool.len() as u32;
                match event {
                    IndelEvent::Insertion { site, .. } if site > len => return Err(PassError::invalid_plan_state(pass,
                        format!("replayed insertion at slot {i} has site {site} > pool length {len}"))),
                    IndelEvent::Deletion { site: Some(site) } if site >= len => return Err(PassError::invalid_plan_state(pass,
                        format!("replayed deletion at slot {i} has site {site} >= pool length {len}"))),
                    IndelEvent::Deletion { site: None } if len != 0 || self.insertion_prob >= 1.0 => {
                        return Err(PassError::invalid_plan_state(pass,
                            "mixed no-op replay event has no positive natural probability"));
                    }
                    _ => {}
                }
                hypothetical = apply_event_hypothetical(&hypothetical, event);
            }
            if let Some(contracts) = tx.ctx().contracts {
                if !tuple_admitted(&initial, &hypothetical, contracts, tx.ctx().refdata, count) {
                    return Err(PassError::constraint_sampling(
                        pass,
                        ChoiceAddress::CorruptIndelCount.to_string(),
                        FilteredSampleError::EmptyAdmissibleSupport,
                    ));
                }
            }
        }
        let applied = events
            .iter()
            .filter(|event| {
                matches!(
                    event,
                    IndelEvent::Insertion { .. } | IndelEvent::Deletion { site: Some(_) }
                )
            })
            .count() as u32;
        for (i, event) in events.into_iter().enumerate() {
            Self::record_event(tx.ctx(), i as u32, event);
            Self::apply_event_via_tx(tx, event)?;
        }
        Self::record_sampling_report(
            tx.ctx(),
            count,
            if compatibility_skip { 0 } else { applied },
            status,
            &SamplingReport {
                strategy: "replay",
                attempts: 0,
                explored_states: work,
            },
        );
        Ok(())
    }
}
