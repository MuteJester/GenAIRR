use crate::ir::Simulation;
use crate::refdata::RefDataConfig;
use crate::trace::ChoiceValue;

use super::super::{ChoiceContext, ChoiceKind, Contract, ContractViolation};
use super::AnchorPreserved;

impl AnchorPreserved {
    pub(super) fn admits_targeted_anchor_substitution(
        &self,
        sim: &Simulation,
        refdata: Option<&RefDataConfig>,
        address: &str,
        candidate: &ChoiceValue,
        context: ChoiceContext<'_>,
    ) -> Result<(), ContractViolation> {
        if context.kind != ChoiceKind::TargetedBaseSubstitution {
            return Ok(());
        }

        let candidate_base = match candidate {
            ChoiceValue::Base(base) => *base,
            _ => return Ok(()),
        };
        let target = match context.target {
            Some(target) => target,
            None => return Ok(()),
        };
        let refdata = match refdata {
            Some(refdata) => refdata,
            None => return Ok(()),
        };
        let inst = match sim.assignments.get(self.segment) {
            Some(inst) => inst,
            None => return Ok(()),
        };
        let allele = match refdata.get(self.segment, inst.allele_id) {
            Some(allele) => allele,
            None => return Ok(()),
        };
        let anchor = match allele.anchor {
            Some(anchor) => anchor as u32,
            None => return Ok(()),
        };

        self.require_anchor_retained(allele, inst.trim_5, inst.trim_3)?;
        let Some(anchor_pool_start) = self.anchor_pool_start(sim, allele, anchor, inst.trim_5)?
        else {
            return Ok(());
        };

        let target_idx = target.index();
        if target_idx < anchor_pool_start || target_idx >= anchor_pool_start + 3 {
            return Ok(());
        }

        let reference_codon = Self::anchor_codon(allele, anchor);
        let mut live_codon = self.live_anchor_codon(sim, anchor_pool_start)?;
        live_codon[(target_idx - anchor_pool_start) as usize] = candidate_base;

        self.require_anchor_amino_acid_preserved(allele, reference_codon, live_codon)
            .map_err(|_| {
                ContractViolation::new(
                    self.name(),
                    format!(
                        "candidate at {} for target {} would change {:?} anchor amino acid",
                        address, target_idx, self.segment
                    ),
                )
            })
    }
}
