use crate::ir::{GermlinePos, NucHandle, Simulation};
use crate::refdata::Allele;

use super::super::{Contract, ContractViolation};
use super::AnchorPreserved;

impl AnchorPreserved {
    pub(super) fn anchor_pool_start(
        &self,
        sim: &Simulation,
        allele: &Allele,
        anchor: u32,
        trim_5: u16,
    ) -> Result<Option<u32>, ContractViolation> {
        let Some(region) = sim
            .sequence
            .regions
            .iter()
            .find(|r| r.segment == self.segment)
        else {
            return Ok(None);
        };

        let trim_5 = trim_5 as u32;
        if anchor < trim_5 {
            return Err(ContractViolation::new(
                self.name(),
                format!(
                    "anchor at allele position {} but trim_5 = {} \
                     (anchor 5'-trimmed away from {} '{}')",
                    anchor, trim_5, allele.name, allele.gene
                ),
            ));
        }

        let anchor_pool_start = region.start.index() + (anchor - trim_5);
        if anchor_pool_start + 3 > region.end.index() {
            return Err(ContractViolation::new(
                self.name(),
                format!(
                    "anchor codon maps to pool range [{}, {}) but assembled {:?} \
                     region ends at {}",
                    anchor_pool_start,
                    anchor_pool_start + 3,
                    self.segment,
                    region.end.index()
                ),
            ));
        }

        Ok(Some(anchor_pool_start))
    }

    pub(super) fn live_anchor_codon(
        &self,
        sim: &Simulation,
        anchor_pool_start: u32,
    ) -> Result<[u8; 3], ContractViolation> {
        let mut codon = [b'N'; 3];
        for offset in 0..3 {
            let pool_pos = anchor_pool_start + offset;
            let Some(nuc) = sim.pool.get(NucHandle::new(pool_pos)) else {
                return Err(ContractViolation::new(
                    self.name(),
                    format!(
                        "anchor codon maps to missing pool position {} in {:?} region",
                        pool_pos, self.segment
                    ),
                ));
            };
            codon[offset as usize] = nuc.base;
        }
        Ok(codon)
    }

    pub(super) fn verify_live_anchor(
        &self,
        sim: &Simulation,
        allele: &Allele,
        anchor: u32,
        trim_5: u16,
    ) -> Result<(), ContractViolation> {
        let Some(anchor_pool_start) = self.anchor_pool_start(sim, allele, anchor, trim_5)? else {
            return Ok(());
        };

        let mut live_codon = [b'N'; 3];
        for offset in 0..3 {
            let pool_pos = anchor_pool_start + offset;
            let expected_germline_pos = GermlinePos::pos((anchor + offset) as u16);
            let Some(nuc) = sim.pool.get(NucHandle::new(pool_pos)) else {
                return Err(ContractViolation::new(
                    self.name(),
                    format!(
                        "anchor codon maps to missing pool position {} in {:?} region",
                        pool_pos, self.segment
                    ),
                ));
            };

            if nuc.segment != self.segment || nuc.germline_pos != expected_germline_pos {
                return Err(ContractViolation::new(
                    self.name(),
                    format!(
                        "anchor codon provenance mismatch at pool position {}: \
                         expected {:?} germline position {:?}, got {:?} germline position {:?}",
                        pool_pos,
                        self.segment,
                        expected_germline_pos.get(),
                        nuc.segment,
                        nuc.germline_pos.get()
                    ),
                ));
            }
            live_codon[offset as usize] = nuc.base;
        }

        let reference_codon = Self::anchor_codon(allele, anchor);
        self.require_anchor_amino_acid_preserved(allele, reference_codon, live_codon)
    }
}
