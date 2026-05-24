use crate::address;
use crate::contract::ChoiceContext;
use crate::dist::FilteredSampleError;
use crate::ir::{flag, NucHandle, Nucleotide, Simulation};
use crate::pass::{Pass, PassContext, PassError};
use crate::passes::constrained::{sample_contract_verified_event, PostEventCandidate};

use super::event::IndelEvent;
use super::IndelPass;

impl IndelPass {
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

    pub(super) fn constrained_candidates(
        &self,
        sim: &Simulation,
        index: u32,
        count: u32,
    ) -> Result<Vec<PostEventCandidate<'static, IndelEvent>>, FilteredSampleError> {
        let pool_len = sim.pool.len() as u32;
        let mut candidates = Vec::new();

        if self.insertion_prob > 0.0 {
            let base_support = self.base_support()?;
            let site_weight = self.insertion_prob / (pool_len + 1) as f64;

            for site in 0..=pool_len {
                let segment = Self::insertion_segment(sim, site);
                for &(base, base_weight) in &base_support {
                    let nuc = Nucleotide::synthetic(base, segment, flag::INDEL_INSERTED);
                    // Skip the codon-rail recompute on the candidate
                    // post-sim: contracts evaluate against `sim.pool`
                    // directly (no contract reads
                    // `region.amino_acids` / `stop_codon_positions`),
                    // and the candidate sim is discarded once the
                    // weight is computed.
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

    pub(super) fn sample_event(
        &self,
        sim: &Simulation,
        ctx: &mut PassContext,
        index: u32,
        count: u32,
        strict: bool,
    ) -> Result<IndelEvent, PassError> {
        if ctx.contracts.is_some() {
            let event_address = address::corrupt_indel_site(index);
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
}
