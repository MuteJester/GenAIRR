use crate::ir::{translate_codon, NucHandle, Segment, Simulation, AMINO_STOP};

use super::layout::np_segment_label;
use super::model::{SlotSource, CONTRACT_NAME};
use super::{ContractViolation, JunctionStopState};

impl JunctionStopState {
    /// Test whether `candidate` at NP `(np_segment, np_index)` would
    /// force a stop codon in the junction. Mirrors the slow path's
    /// `reject_known_stop` decision but reads at most one codon.
    pub fn admits_np_candidate(
        &self,
        sim: &Simulation,
        np_segment: Segment,
        np_index: u32,
        candidate: u8,
    ) -> Result<(), ContractViolation> {
        if let Some(violation) = self.static_violation {
            return Err(ContractViolation::new(
                CONTRACT_NAME,
                format!(
                    "candidate at np.{}.bases[{}] would force stop codon '{}{}{}' at junction offset {}",
                    np_segment_label(np_segment),
                    np_index,
                    violation.codon[0] as char,
                    violation.codon[1] as char,
                    violation.codon[2] as char,
                    violation.offset
                ),
            ));
        }

        let Some(cand_offset) = self.find_candidate_offset(np_segment, np_index) else {
            return Ok(());
        };

        let codon_start = (cand_offset / 3) * 3;
        if codon_start + 3 > self.slots.len() as u32 {
            return Ok(());
        }

        let mut codon = [0u8; 3];
        for k in 0..3u32 {
            let off = codon_start + k;
            if off == cand_offset {
                codon[k as usize] = candidate;
                continue;
            }
            match self.read_slot(off, sim, np_segment, np_index, candidate) {
                Some(b) => codon[k as usize] = b,
                None => return Ok(()),
            }
        }

        if translate_codon(codon[0], codon[1], codon[2]) == AMINO_STOP {
            return Err(ContractViolation::new(
                CONTRACT_NAME,
                format!(
                    "candidate at np.{}.bases[{}] would force stop codon '{}{}{}' at junction offset {}",
                    np_segment_label(np_segment),
                    np_index,
                    codon[0] as char,
                    codon[1] as char,
                    codon[2] as char,
                    codon_start
                ),
            ));
        }

        Ok(())
    }

    fn find_candidate_offset(&self, np_segment: Segment, np_index: u32) -> Option<u32> {
        for (i, slot) in self.slots.iter().enumerate() {
            if let SlotSource::Np {
                segment,
                np_index: ni,
                ..
            } = slot.source
            {
                if segment == np_segment && ni == np_index {
                    return Some(i as u32);
                }
            }
        }
        None
    }

    fn read_slot(
        &self,
        offset: u32,
        sim: &Simulation,
        cand_segment: Segment,
        cand_index: u32,
        candidate: u8,
    ) -> Option<u8> {
        let slot = self.slots.get(offset as usize)?;
        match slot.source {
            SlotSource::Fixed { fixed_index } => {
                self.fixed_bytes.get(fixed_index as usize).copied()
            }
            SlotSource::Np {
                segment,
                np_index,
                pool_pos,
            } => {
                if segment == cand_segment && np_index == cand_index {
                    return Some(candidate);
                }
                if (pool_pos as usize) < sim.pool.len() {
                    sim.pool.get(NucHandle::new(pool_pos)).map(|n| n.base)
                } else {
                    None
                }
            }
        }
    }
}
