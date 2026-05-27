use crate::dist::FilteredSampleError;
use crate::ir::Simulation;
use crate::pass::PassContext;

use super::event::IndelEvent;
use super::IndelPass;

impl IndelPass {
    /// Normalize the configured insertion-base distribution's
    /// support into a `(base, weight)` list with weights summing
    /// to 1. Used by [`Self::sample_admissible_tuple`] when
    /// computing per-event natural weights for the mod-3 DP.
    pub(super) fn base_support(&self) -> Result<Vec<(u8, f64)>, FilteredSampleError> {
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

    /// Unconstrained per-step indel draw, used by the no-contracts
    /// path in [`IndelPass::execute_with_sampling_mode`]. Under
    /// active contracts the indel pass runs
    /// [`Self::sample_admissible_tuple`] instead.
    pub(super) fn sample_legacy_event(
        &self,
        sim: &Simulation,
        ctx: &mut PassContext,
    ) -> IndelEvent {
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
}
