//! Structural-edit API for [`MutationTransaction`].
//!
//! Exposes the four indel entry points:
//! - [`MutationTransaction::insert_base`] — raw insertion that
//!   bypasses the active contract bundle.
//! - [`MutationTransaction::insert_base_admitting`] — insertion
//!   with contract arbitration via `admits_post_event`.
//! - [`MutationTransaction::delete_base`] — raw deletion that
//!   bypasses contracts.
//! - [`MutationTransaction::delete_base_admitting`] — deletion
//!   with contract arbitration.

use crate::address::ChoiceAddress;
use crate::contract::ChoiceContext;
use crate::ir::{NucHandle, Nucleotide, Segment};
use crate::pass::PassError;

use super::MutationTransaction;

#[allow(dead_code)]
impl<'a, 'idx> MutationTransaction<'a, 'idx> {
    /// Insert a nucleotide at pool position `at`. Position must be
    /// in `[0, pool_len]` (insertion at the end is allowed). Returns
    /// `Err` for out-of-range positions.
    ///
    /// **Bypasses the active contract bundle**. See
    /// [`Self::insert_base_admitting`] for the contract-aware
    /// variant.
    pub(crate) fn insert_base(&mut self, at: u32, nucleotide: Nucleotide) -> Result<(), PassError> {
        let pool_len = self.builder.peek().pool.len() as u32;
        if at > pool_len {
            return Err(PassError::invalid_plan_state(
                self.pass_name.to_string(),
                format!(
                    "insert_base: position {} out of pool range [0, {}]",
                    at, pool_len
                ),
            ));
        }
        self.builder.insert_indel(at, nucleotide);
        Ok(())
    }

    /// Insert a nucleotide at pool position `at`, contract-aware.
    /// Builds the hypothetical post-insertion IR and asks the
    /// contract bundle to arbitrate. Skips in permissive mode if
    /// rejected; returns `Err` in strict mode.
    ///
    /// Returns:
    /// - `Ok(true)` — insertion applied
    /// - `Ok(false)` — contract rejected, skipped (permissive)
    /// - `Err(..)` — out-of-range position, or strict-mode rejection
    pub(crate) fn insert_base_admitting(
        &mut self,
        at: u32,
        nucleotide: Nucleotide,
        address: ChoiceAddress,
    ) -> Result<bool, PassError> {
        let pool_len = self.builder.peek().pool.len() as u32;
        if at > pool_len {
            return Err(PassError::invalid_plan_state(
                self.pass_name.to_string(),
                format!(
                    "insert_base_admitting: position {} out of pool range [0, {}]",
                    at, pool_len
                ),
            ));
        }
        if let Some(contracts) = self.ctx.contracts {
            let post_sim = self.builder.peek().with_indel_inserted(at, nucleotide);
            let context =
                ChoiceContext::indel_insertion(0, 1, NucHandle::new(at)).with_address(address);
            if contracts
                .admits_post_event(self.builder.peek(), &post_sim, self.ctx.refdata, context)
                .is_err()
            {
                if self.strict {
                    return Err(PassError::constraint_sampling(
                        self.pass_name.clone(),
                        address.to_string(),
                        crate::dist::FilteredSampleError::EmptyAdmissibleSupport,
                    ));
                }
                return Ok(false);
            }
        }
        self.builder.insert_indel(at, nucleotide);
        Ok(true)
    }

    /// Delete the nucleotide at pool position `at`. Returns the
    /// removed nucleotide, or `Err` if `at >= pool_len` or the pool
    /// is empty.
    ///
    /// **Bypasses the active contract bundle**. See
    /// [`Self::delete_base_admitting`] for the contract-aware
    /// variant.
    pub(crate) fn delete_base(&mut self, at: u32) -> Result<Nucleotide, PassError> {
        let pool_len = self.builder.peek().pool.len() as u32;
        if at >= pool_len {
            return Err(PassError::invalid_plan_state(
                self.pass_name.to_string(),
                format!(
                    "delete_base: position {} out of pool range [0, {})",
                    at, pool_len
                ),
            ));
        }
        self.builder.delete_indel(at).ok_or_else(|| {
            PassError::invalid_plan_state(
                self.pass_name.to_string(),
                format!("delete_base: pool returned None for position {}", at),
            )
        })
    }

    /// Delete the nucleotide at pool position `at`, contract-aware.
    /// Builds the hypothetical post-deletion IR and asks the
    /// contract bundle to arbitrate. Skips in permissive mode if
    /// rejected; returns `Err` in strict mode.
    ///
    /// Returns:
    /// - `Ok(true)` — deletion applied
    /// - `Ok(false)` — contract rejected, skipped (permissive)
    /// - `Err(..)` — out-of-range position, or strict-mode rejection
    pub(crate) fn delete_base_admitting(
        &mut self,
        at: u32,
        address: ChoiceAddress,
    ) -> Result<bool, PassError> {
        let pool_len = self.builder.peek().pool.len() as u32;
        if at >= pool_len {
            return Err(PassError::invalid_plan_state(
                self.pass_name.to_string(),
                format!(
                    "delete_base_admitting: position {} out of pool range [0, {})",
                    at, pool_len
                ),
            ));
        }
        // End loss changes which portion of the original coding frame
        // is observed. Carry that distinction in persistent state so
        // the first-base structural deletion cannot masquerade as loss
        // of read coverage, including after a clonal fork.
        let removed = &self.builder.peek().pool.as_slice()[at as usize];
        let observed_loss_end = match address {
            ChoiceAddress::CorruptEndLoss(end) => Some(end),
            _ => None,
        };
        let mut observed_loss: Vec<_> = observed_loss_end
            .and_then(|end| {
                let segment = removed.segment;
                self.builder
                    .peek()
                    .assignments
                    .get(segment)
                    .copied()
                    .map(|mut inst| {
                        use crate::address::PrimeEnd;
                        if segment == Segment::V && end == PrimeEnd::Five {
                            inst.observation_loss_5 = inst.observation_loss_5.saturating_add(1);
                        } else if segment == Segment::J && end == PrimeEnd::Three {
                            inst.observation_loss_3 = inst.observation_loss_3.saturating_add(1);
                        }
                        // Reference coverage follows the removed germline
                        // coordinate, not the number of query bases removed.
                        // An inserted base consumes no reference coverage.
                        if let Some(position) = removed.germline_pos.get().map(u32::from) {
                            let clips_reference_start =
                                (end == PrimeEnd::Five) != inst.orientation.is_reverse();
                            if clips_reference_start {
                                inst.observed_ref_start =
                                    Some(inst.observed_ref_start.unwrap_or(0).max(position + 1));
                            } else {
                                inst.observed_ref_end =
                                    Some(inst.observed_ref_end.unwrap_or(u32::MAX).min(position));
                            }
                        }
                        // Removing the last observed query base makes this
                        // segment wholly unobserved, including any old terminal
                        // deletion gaps that no longer have read coverage.
                        if self.builder.peek().sequence.regions.iter().any(|region| {
                            region.segment == segment
                                && region.start.index() == at
                                && region.end.index() == at + 1
                        }) {
                            inst.observed_ref_start = Some(0);
                            inst.observed_ref_end = Some(0);
                        }
                        (segment, inst)
                    })
            })
            .into_iter()
            .collect();
        if let Some(end) = observed_loss_end {
            // A segment can already be query-empty after structural
            // deletions. Crossing its zero-width boundary removes coverage
            // even though no deleted nucleotide carries that segment tag.
            let crossed_boundary = match end {
                crate::address::PrimeEnd::Five => at,
                crate::address::PrimeEnd::Three => at + 1,
            };
            for region in &self.builder.peek().sequence.regions {
                if region.start.index() == crossed_boundary
                    && region.end.index() == crossed_boundary
                {
                    if let Some(mut instance) =
                        self.builder.peek().assignments.get(region.segment).copied()
                    {
                        instance.observed_ref_start = Some(0);
                        instance.observed_ref_end = Some(0);
                        observed_loss.push((region.segment, instance));
                    }
                }
            }
        }
        if let Some(contracts) = self.ctx.contracts {
            let mut post_sim = self.builder.peek().with_indel_deleted(at);
            for &(segment, inst) in &observed_loss {
                post_sim = post_sim.with_allele_assigned(segment, inst);
            }
            let context =
                ChoiceContext::indel_deletion(0, 1, NucHandle::new(at)).with_address(address);
            if contracts
                .admits_post_event(self.builder.peek(), &post_sim, self.ctx.refdata, context)
                .is_err()
            {
                if self.strict {
                    return Err(PassError::constraint_sampling(
                        self.pass_name.clone(),
                        address.to_string(),
                        crate::dist::FilteredSampleError::EmptyAdmissibleSupport,
                    ));
                }
                return Ok(false);
            }
        }
        self.builder.delete_indel(at).ok_or_else(|| {
            PassError::invalid_plan_state(
                self.pass_name.to_string(),
                format!(
                    "delete_base_admitting: pool returned None for position {}",
                    at
                ),
            )
        })?;
        for (segment, inst) in observed_loss {
            self.builder.assign_allele(segment, inst);
        }
        Ok(true)
    }
}
