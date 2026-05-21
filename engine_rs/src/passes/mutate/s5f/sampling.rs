use crate::contract::ChoiceContext;
use crate::dist::FilteredSampleError;
use crate::ir::{NucHandle, NucleotidePool, Simulation};

use super::event::{S5FMutationEvent, WeightedS5FMutationEvent};
use super::S5FMutationPass;

impl S5FMutationPass {
    /// Build a 5-mer at `pos` from the current pool, encoded as a
    /// context index. Returns `None` if `pos` is too close to the
    /// pool boundary or if any base in the 5-mer is non-A/C/G/T.
    pub(super) fn context_at(pool: &NucleotidePool, pos: u32) -> Option<u16> {
        if pos < 2 || pos + 2 >= pool.len() as u32 {
            return None;
        }
        let b1 = pool.get(NucHandle::new(pos - 2))?.base;
        let b2 = pool.get(NucHandle::new(pos - 1))?.base;
        let b3 = pool.get(NucHandle::new(pos))?.base;
        let b4 = pool.get(NucHandle::new(pos + 1))?.base;
        let b5 = pool.get(NucHandle::new(pos + 2))?.base;
        crate::s5f::S5FKernel::encode_context(b1, b2, b3, b4, b5)
    }

    /// Walk the pool and build the per-position mutability profile.
    pub(super) fn build_profile(&self, pool: &NucleotidePool) -> Vec<(u32, f64)> {
        let n = pool.len() as u32;
        if n < 5 {
            return Vec::new();
        }
        let mut profile = Vec::with_capacity((n - 4) as usize);
        for pos in 2..n - 2 {
            if let Some(ctx) = Self::context_at(pool, pos) {
                let mu = self.kernel.mutability(ctx);
                if mu > 0.0 {
                    profile.push((pos, mu));
                }
            }
        }
        profile
    }

    pub(super) fn positive_row_candidates(row: [f64; 4]) -> Vec<(u8, f64)> {
        const BASES: [u8; 4] = [b'A', b'C', b'G', b'T'];
        BASES
            .iter()
            .copied()
            .zip(row)
            .filter(|(_, weight)| *weight > 0.0)
            .collect()
    }

    pub(super) fn sample_weighted_base(
        rng: &mut crate::rng::Rng,
        candidates: &[(u8, f64)],
    ) -> Result<Option<u8>, FilteredSampleError> {
        if candidates.is_empty() {
            return Ok(None);
        }

        let total: f64 = candidates.iter().map(|(_, weight)| weight).sum();
        if !total.is_finite() || total <= 0.0 {
            return Err(FilteredSampleError::InvalidFilteredSupport);
        }

        let r = rng.next_f64() * total;
        let mut cum = 0.0;
        for &(base, weight) in candidates {
            cum += weight;
            if r < cum {
                return Ok(Some(base));
            }
        }

        Ok(Some(candidates.last().expect("non-empty candidates").0))
    }

    pub(super) fn event_candidates(
        &self,
        sim: &Simulation,
        profile: &[(u32, f64)],
        index: u32,
        count: u32,
    ) -> Result<Vec<WeightedS5FMutationEvent>, FilteredSampleError> {
        let mut events = Vec::new();

        for &(site, mutability) in profile {
            if !mutability.is_finite() || mutability <= 0.0 {
                continue;
            }

            let Some(context) = Self::context_at(&sim.pool, site) else {
                continue;
            };
            let row = self.kernel.substitution_row(context);
            let row_candidates = Self::positive_row_candidates(row);
            if row_candidates.is_empty() {
                continue;
            }
            let row_total: f64 = row_candidates.iter().map(|(_, weight)| *weight).sum();
            if !row_total.is_finite() || row_total <= 0.0 {
                return Err(FilteredSampleError::InvalidFilteredSupport);
            }

            for (base, weight) in row_candidates {
                events.push(WeightedS5FMutationEvent {
                    event: S5FMutationEvent { site, base },
                    weight: mutability * (weight / row_total),
                    context: ChoiceContext::targeted_base_substitution(
                        index,
                        count,
                        NucHandle::new(site),
                    ),
                });
            }
        }

        Ok(events)
    }

    pub(super) fn sample_weighted_event(
        rng: &mut crate::rng::Rng,
        candidates: &[WeightedS5FMutationEvent],
    ) -> Result<Option<S5FMutationEvent>, FilteredSampleError> {
        if candidates.is_empty() {
            return Ok(None);
        }

        let total: f64 = candidates.iter().map(|candidate| candidate.weight).sum();
        if !total.is_finite() || total <= 0.0 {
            return Err(FilteredSampleError::InvalidFilteredSupport);
        }

        let r = rng.next_f64() * total;
        let mut cum = 0.0;
        for candidate in candidates {
            cum += candidate.weight;
            if r < cum {
                return Ok(Some(candidate.event));
            }
        }

        Ok(Some(
            candidates
                .last()
                .expect("non-empty candidates checked above")
                .event,
        ))
    }
}
