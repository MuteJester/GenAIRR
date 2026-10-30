//! Whole observed coding-sequence protection for the productive bundle.
//! Junction generation still uses its dedicated incremental NP contract.

use crate::coding::{anchor_position, coding_frame, j_tail_in_frame};
use crate::ir::{NucHandle, Segment, Simulation};
use crate::productivity::SequenceView;
use crate::refdata::RefDataConfig;
use crate::trace::ChoiceValue;

use super::{BaseMask, ChoiceContext, ChoiceKind, Contract, ContractKind, ContractViolation};

pub struct CodingSequenceIntegrity;

impl CodingSequenceIntegrity {
    fn allows_base(
        &self,
        sim: &Simulation,
        refdata: Option<&RefDataConfig>,
        site: NucHandle,
        base: u8,
    ) -> bool {
        let Some(frame) = refdata.and_then(|r| coding_frame(sim, r)) else {
            return true;
        };
        let p = site.index() as usize;
        if p < frame.start || p >= frame.end {
            return true;
        }
        let codon_start = p - (p - frame.start) % 3;
        if codon_start + 3 > frame.end {
            return true;
        }
        // The compiled allele/trim support and structural post-event checks
        // establish the stop-free invariant. A substitution can change only
        // this one codon; scanning every other codon for every site/base would
        // make S5F support construction quadratic in sequence length.
        let mut codon = [b'N'; 3];
        for (offset, byte) in codon.iter_mut().enumerate() {
            let position = codon_start + offset;
            *byte = if position == p {
                base
            } else {
                sim.pool
                    .get(NucHandle::new(position as u32))
                    .map(|n| n.base)
                    .unwrap_or(b'N')
            };
        }
        crate::productivity::codon_stop(codon) == Some(false)
    }
}

impl Contract for CodingSequenceIntegrity {
    fn name(&self) -> &str {
        "coding_sequence_integrity"
    }

    fn kind(&self) -> ContractKind {
        ContractKind::CodingSequenceIntegrity
    }

    fn verify(
        &self,
        sim: &Simulation,
        refdata: Option<&RefDataConfig>,
    ) -> Result<(), ContractViolation> {
        let Some(refdata) = refdata else {
            return Ok(());
        };
        let Some(frame) = coding_frame(sim, refdata) else {
            return Ok(());
        };
        if sim.coding_stop(frame.start, frame.end) != Some(false) {
            return Err(ContractViolation::new(
                self.name(),
                "coding sequence is not definitely stop-free",
            ));
        }
        for segment in [Segment::V, Segment::J] {
            let Some(instance) = sim.assignments.get(segment) else {
                continue;
            };
            let Some(allele) = refdata.get(segment, instance.allele_id) else {
                continue;
            };
            let Some(anchor) = allele.anchor else {
                continue;
            };
            let Some(region) = sim.sequence.regions.iter().find(|r| r.segment == segment) else {
                continue;
            };
            if let Some(position) = anchor_position(sim, region, u32::from(anchor)) {
                if position as usize % 3 != frame.phase {
                    return Err(ContractViolation::new(
                        self.name(),
                        format!("{segment:?} anchor is out of the retained V coding frame"),
                    ));
                }
            }
        }
        if !j_tail_in_frame(sim, refdata, frame) {
            return Err(ContractViolation::new(
                self.name(),
                "J coding tail has a frameshift",
            ));
        }
        Ok(())
    }

    fn admits_typed(
        &self,
        sim: &Simulation,
        refdata: Option<&RefDataConfig>,
        context: ChoiceContext<'_>,
        candidate: &ChoiceValue,
    ) -> Result<(), ContractViolation> {
        if context.kind == ChoiceKind::TargetedBaseSubstitution {
            if let (Some(site), ChoiceValue::Base(base)) = (context.target, candidate) {
                if !self.allows_base(sim, refdata, site, *base) {
                    return Err(ContractViolation::new(
                        self.name(),
                        "candidate creates a coding stop",
                    ));
                }
            }
        }
        Ok(())
    }

    fn admissible_bases_at(
        &self,
        sim: &Simulation,
        refdata: Option<&RefDataConfig>,
        site: NucHandle,
    ) -> BaseMask {
        // Everything here except the candidate byte is the same for all four
        // bases, so it is computed once rather than four times. This hook runs
        // per site while constrained sampling builds its support, and this
        // contract joined the productive bundle after the recorded budget
        // baselines were taken -- which is most of why the two productive
        // budgets drifted while the other four did not.
        //
        // The repeated work was the coding frame, which walks the simulation
        // through `coding_span`, plus the two codon bytes that are not this
        // site. `allows_base` still exists for `admits_fixed_base_at`, which
        // answers about one byte and cannot amortize anything.
        const ALL: BaseMask = BaseMask(0b1111);
        let Some(frame) = refdata.and_then(|r| coding_frame(sim, r)) else {
            return ALL;
        };
        let p = site.index() as usize;
        if p < frame.start || p >= frame.end {
            return ALL;
        }
        let codon_start = p - (p - frame.start) % 3;
        if codon_start + 3 > frame.end {
            return ALL;
        }
        let mut codon = [b'N'; 3];
        for (offset, byte) in codon.iter_mut().enumerate() {
            let position = codon_start + offset;
            if position != p {
                *byte = sim
                    .pool
                    .get(NucHandle::new(position as u32))
                    .map(|n| n.base)
                    .unwrap_or(b'N');
            }
        }
        let slot = p - codon_start;
        let mut mask = 0u8;
        for bit in 0..4 {
            codon[slot] = BaseMask::base_for_bit(bit);
            if crate::productivity::codon_stop(codon) == Some(false) {
                mask |= 1 << bit;
            }
        }
        BaseMask(mask)
    }

    fn admits_fixed_base_at(
        &self,
        sim: &Simulation,
        refdata: Option<&RefDataConfig>,
        site: NucHandle,
        byte: u8,
    ) -> bool {
        self.allows_base(sim, refdata, site, byte)
    }
}
