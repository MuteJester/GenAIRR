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
use crate::ir::{NucHandle, Nucleotide};
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
        if let Some(contracts) = self.ctx.contracts {
            let post_sim = self.builder.peek().with_indel_deleted(at);
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
        Ok(true)
    }
}
