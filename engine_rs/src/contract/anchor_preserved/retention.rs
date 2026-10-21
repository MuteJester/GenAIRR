use crate::ir::translate_codon;
use crate::ir::Simulation;
use crate::refdata::{Allele, RefDataConfig};

use super::super::{Contract, ContractViolation};
use super::AnchorPreserved;

impl AnchorPreserved {
    /// A duplex continuation may carry a Verified positional landmark whose
    /// source residue was noncanonical but whose current codon is accepted.
    /// Its live mutation support depends on that codon, not the source class.
    /// Legacy candidate filtering retains its existing source-residue rule.
    pub(super) fn require_live_anchor_retained(
        &self,
        sim: &Simulation,
        refdata: &RefDataConfig,
        allele: &Allele,
        trim_5: u16,
        trim_3: u16,
    ) -> Result<(), ContractViolation> {
        if self.canonical_class && sim.duplex_provenance.is_some() {
            let ready = sim.assignments.get(self.segment).is_some_and(|inst| {
                crate::productivity::reference_ready(refdata, self.segment, inst.allele_id)
            });
            if !ready {
                return Err(ContractViolation::new(
                    self.name(),
                    "reference_anchor_unresolved",
                ));
            }
            return Self::new(self.segment)
                .require_anchor_retained(allele, trim_5, trim_3)
                .map_err(|violation| ContractViolation::new(self.name(), violation.reason));
        }
        self.require_anchor_retained(allele, trim_5, trim_3)
    }

    pub(super) fn require_anchor_retained(
        &self,
        allele: &Allele,
        trim_5: u16,
        trim_3: u16,
    ) -> Result<(), ContractViolation> {
        if self.canonical_class && !crate::productivity::reference_anchor_ready(allele) {
            return Err(ContractViolation::new(
                self.name(),
                "reference_anchor_unresolved",
            ));
        }
        let Some(anchor) = allele.anchor.map(|anchor| anchor as u32) else {
            return Err(ContractViolation::new(
                self.name(),
                format!("{} '{}' has no anchor", allele.name, allele.gene),
            ));
        };
        if anchor + 3 > allele.len() {
            return Err(ContractViolation::new(
                self.name(),
                format!(
                    "anchor codon at [{}, {}) exceeds allele length {} for {} '{}'",
                    anchor,
                    anchor + 3,
                    allele.len(),
                    allele.name,
                    allele.gene
                ),
            ));
        }

        let trim_5 = trim_5 as u32;
        let trim_3 = trim_3 as u32;
        let retained_end = allele.len().saturating_sub(trim_3);
        if anchor < trim_5 {
            return Err(ContractViolation::new(
                self.name(),
                format!(
                    "candidate would 5'-trim anchor at allele position {} from {} '{}'",
                    anchor, allele.name, allele.gene
                ),
            ));
        }
        if anchor + 3 > retained_end {
            return Err(ContractViolation::new(
                self.name(),
                format!(
                    "candidate would 3'-trim anchor codon [{}, {}) from {} '{}'",
                    anchor,
                    anchor + 3,
                    allele.name,
                    allele.gene
                ),
            ));
        }
        Ok(())
    }

    pub(super) fn anchor_codon(allele: &Allele, anchor: u32) -> [u8; 3] {
        [
            allele.seq[anchor as usize],
            allele.seq[anchor as usize + 1],
            allele.seq[anchor as usize + 2],
        ]
    }

    pub(super) fn codon_amino_acid(codon: [u8; 3]) -> u8 {
        translate_codon(codon[0], codon[1], codon[2])
    }

    pub(super) fn require_anchor_amino_acid_preserved(
        &self,
        allele: &Allele,
        reference_codon: [u8; 3],
        live_codon: [u8; 3],
    ) -> Result<(), ContractViolation> {
        if self.canonical_class {
            return if crate::productivity::anchor_class(self.segment, live_codon) == Some(true) {
                Ok(())
            } else {
                Err(ContractViolation::new(
                    self.name(),
                    "live anchor is not definitely in the accepted residue class",
                ))
            };
        }
        let reference_aa = Self::codon_amino_acid(reference_codon);
        let live_aa = Self::codon_amino_acid(live_codon);
        if live_aa == reference_aa {
            return Ok(());
        }

        Err(ContractViolation::new(
            self.name(),
            format!(
                "anchor codon amino acid changed for {} '{}': reference {}{}{} -> {}, live {}{}{} -> {}",
                allele.name,
                allele.gene,
                reference_codon[0] as char,
                reference_codon[1] as char,
                reference_codon[2] as char,
                reference_aa as char,
                live_codon[0] as char,
                live_codon[1] as char,
                live_codon[2] as char,
                live_aa as char
            ),
        ))
    }
}
