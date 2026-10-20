use crate::dist::FilteredSampleError;
use crate::ir::Simulation;
use crate::pass::PassContext;

use super::event::IndelEvent;
use super::IndelPass;

impl IndelPass {
    /// Positive support in log probability space. Zero-weight entries are not
    /// possible events; finite extreme weights must not disappear by underflow.
    pub(super) fn base_log_support(&self) -> Result<Vec<(u8, f64)>, FilteredSampleError> {
        let support = self
            .base_dist
            .support()
            .ok_or(FilteredSampleError::SupportUnavailable)?;
        if support.iter().any(|(_, w)| !w.is_finite() || *w < 0.0) {
            return Err(FilteredSampleError::InvalidFilteredSupport);
        }
        let mut combined = std::collections::BTreeMap::<u8, Vec<f64>>::new();
        for (base, weight) in support {
            if weight > 0.0 {
                combined.entry(base).or_default().push(weight.ln());
            }
        }
        if combined.is_empty() {
            return Err(FilteredSampleError::InvalidFilteredSupport);
        }
        let sum_logs = |values: &[f64]| {
            let max = values.iter().copied().fold(f64::NEG_INFINITY, f64::max);
            max + values.iter().map(|x| (*x - max).exp()).sum::<f64>().ln()
        };
        let mut support: Vec<_> = combined
            .into_iter()
            .map(|(base, weights)| (base, sum_logs(&weights)))
            .collect();
        let total = sum_logs(&support.iter().map(|(_, w)| *w).collect::<Vec<_>>());
        for (_, log_weight) in &mut support {
            *log_weight -= total;
        }
        Ok(support)
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
