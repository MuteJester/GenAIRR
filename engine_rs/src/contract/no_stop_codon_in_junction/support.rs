use crate::contract::BaseMask;
use crate::ir::{translate_codon, NucHandle, Simulation, AMINO_STOP};
use crate::junction::compute_junction;
use crate::refdata::RefDataConfig;

use super::NoStopCodonInJunction;

impl NoStopCodonInJunction {
    /// v3.0 constrain-before-propose: per-site admissible-base mask.
    pub(super) fn admissible_bases_at_impl(
        &self,
        sim: &Simulation,
        refdata: Option<&RefDataConfig>,
        site: NucHandle,
    ) -> BaseMask {
        let Some(refdata) = refdata else {
            return BaseMask::UNCONSTRAINED;
        };
        let Some(junction) = compute_junction(sim, refdata) else {
            return BaseMask::UNCONSTRAINED;
        };
        if !junction.is_in_frame() {
            return BaseMask::UNCONSTRAINED;
        }
        let start = junction.start.index();
        let end = junction.end.index();
        let target_idx = site.index();
        if target_idx < start || target_idx >= end {
            return BaseMask::UNCONSTRAINED;
        }

        let mut mask: u8 = 0;
        for bit in 0u8..4 {
            let candidate_base = BaseMask::base_for_bit(bit);
            let mut admissible = true;
            let mut pos = start;
            while pos + 3 <= end {
                let b1 = pool_base_with_candidate(sim, NucHandle::new(pos), site, candidate_base);
                let b2 =
                    pool_base_with_candidate(sim, NucHandle::new(pos + 1), site, candidate_base);
                let b3 =
                    pool_base_with_candidate(sim, NucHandle::new(pos + 2), site, candidate_base);
                if let (Some(b1), Some(b2), Some(b3)) = (b1, b2, b3) {
                    if translate_codon(b1, b2, b3) == AMINO_STOP {
                        admissible = false;
                        break;
                    }
                }
                pos += 3;
            }
            if admissible {
                mask |= 1 << bit;
            }
        }
        BaseMask(mask)
    }

    /// Pinned non-canonical write: walk junction codons with `byte`
    /// substituted at `site`, reject if any resulting codon
    /// translates to a stop.
    pub(super) fn admits_fixed_base_at_impl(
        &self,
        sim: &Simulation,
        refdata: Option<&RefDataConfig>,
        site: NucHandle,
        byte: u8,
    ) -> bool {
        let Some(refdata) = refdata else {
            return true;
        };
        let Some(junction) = compute_junction(sim, refdata) else {
            return true;
        };
        if !junction.is_in_frame() {
            return true;
        }
        let start = junction.start.index();
        let end = junction.end.index();
        let target_idx = site.index();
        if target_idx < start || target_idx >= end {
            return true;
        }
        let mut pos = start;
        while pos + 3 <= end {
            let b1 = pool_base_with_candidate(sim, NucHandle::new(pos), site, byte);
            let b2 = pool_base_with_candidate(sim, NucHandle::new(pos + 1), site, byte);
            let b3 = pool_base_with_candidate(sim, NucHandle::new(pos + 2), site, byte);
            if let (Some(b1), Some(b2), Some(b3)) = (b1, b2, b3) {
                if translate_codon(b1, b2, b3) == AMINO_STOP {
                    return false;
                }
            }
            pos += 3;
        }
        true
    }
}

/// Pool base at `handle`, substituted with `candidate_base` if
/// `handle == target`.
fn pool_base_with_candidate(
    sim: &Simulation,
    handle: NucHandle,
    target: NucHandle,
    candidate_base: u8,
) -> Option<u8> {
    if handle == target {
        Some(candidate_base)
    } else {
        sim.pool.get(handle).map(|n| n.base)
    }
}
