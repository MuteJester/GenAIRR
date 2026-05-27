use crate::contract::{
    BaseMask, IndelEventClass, IndelKindHint, LengthSupport, TrimEnd, TrimTarget,
};
use crate::ir::{NucHandle, Simulation};
use crate::refdata::RefDataConfig;

use super::AnchorPreserved;

impl AnchorPreserved {
    /// v3.0 constrain-before-propose: per-site admissible-base mask
    /// for the anchor codon.
    pub(super) fn admissible_bases_at_impl(
        &self,
        sim: &Simulation,
        refdata: Option<&RefDataConfig>,
        site: NucHandle,
    ) -> BaseMask {
        let Some(refdata) = refdata else {
            return BaseMask::UNCONSTRAINED;
        };
        let Some(inst) = sim.assignments.get(self.segment).copied() else {
            return BaseMask::UNCONSTRAINED;
        };
        let Some(allele) = refdata.get(self.segment, inst.allele_id) else {
            return BaseMask::UNCONSTRAINED;
        };
        let Some(anchor) = allele.anchor.map(|a| a as u32) else {
            return BaseMask::UNCONSTRAINED;
        };

        // If the anchor is already not retained, there is no live
        // anchor codon to classify locally. Verification / trim
        // support report the broken invariant; base masks stay
        // unconstrained.
        if self
            .require_anchor_retained(allele, inst.trim_5, inst.trim_3)
            .is_err()
        {
            return BaseMask::UNCONSTRAINED;
        }
        let anchor_pool_start = match self.anchor_pool_start(sim, allele, anchor, inst.trim_5) {
            Ok(Some(start)) => start,
            Ok(None) | Err(_) => return BaseMask::UNCONSTRAINED,
        };

        let target_idx = site.index();
        if target_idx < anchor_pool_start || target_idx >= anchor_pool_start + 3 {
            return BaseMask::UNCONSTRAINED;
        }

        let reference_codon = Self::anchor_codon(allele, anchor);
        let live_codon = match self.live_anchor_codon(sim, anchor_pool_start) {
            Ok(codon) => codon,
            Err(_) => return BaseMask::UNCONSTRAINED,
        };
        let codon_offset = (target_idx - anchor_pool_start) as usize;

        let mut mask: u8 = 0;
        for bit in 0u8..4 {
            let candidate_base = BaseMask::base_for_bit(bit);
            let mut hypothetical = live_codon;
            hypothetical[codon_offset] = candidate_base;
            if self
                .require_anchor_amino_acid_preserved(allele, reference_codon, hypothetical)
                .is_ok()
            {
                mask |= 1 << bit;
            }
        }
        BaseMask(mask)
    }

    /// Classify an indel candidate's effect on this segment's
    /// anchor codon.
    pub(super) fn admissible_indel_class_at_impl(
        &self,
        sim: &Simulation,
        refdata: Option<&RefDataConfig>,
        site: u32,
        kind: IndelKindHint,
    ) -> IndelEventClass {
        let Some(refdata) = refdata else {
            return IndelEventClass::FrameNeutral;
        };
        let Some(inst) = sim.assignments.get(self.segment).copied() else {
            return IndelEventClass::FrameNeutral;
        };
        let Some(allele) = refdata.get(self.segment, inst.allele_id) else {
            return IndelEventClass::FrameNeutral;
        };
        let Some(anchor) = allele.anchor.map(|a| a as u32) else {
            return IndelEventClass::FrameNeutral;
        };
        if self
            .require_anchor_retained(allele, inst.trim_5, inst.trim_3)
            .is_err()
        {
            return IndelEventClass::FrameNeutral;
        }
        let region = match sim
            .sequence
            .regions
            .iter()
            .find(|r| r.segment == self.segment)
        {
            Some(r) => r,
            None => return IndelEventClass::FrameNeutral,
        };
        let region_start = region.start.index();
        let trim_5 = inst.trim_5 as u32;
        let anchor_pool = region_start + (anchor - trim_5);
        let anchor_end = anchor_pool + 3;
        if site >= region_start && site < anchor_end {
            IndelEventClass::Forbidden
        } else {
            let _ = kind;
            IndelEventClass::FrameNeutral
        }
    }

    /// Closed-form admissible trim-length support for this
    /// segment's anchor.
    pub(super) fn admissible_trim_lengths_impl(
        &self,
        sim: &Simulation,
        refdata: Option<&RefDataConfig>,
        target: TrimTarget,
        requested_max: u32,
    ) -> LengthSupport {
        if target.segment != self.segment {
            return LengthSupport::Full(requested_max);
        }
        let Some(refdata) = refdata else {
            return LengthSupport::Full(requested_max);
        };
        let Some(inst) = sim.assignments.get(self.segment).copied() else {
            return LengthSupport::Full(requested_max);
        };
        let Some(allele) = refdata.get(self.segment, inst.allele_id) else {
            return LengthSupport::Full(requested_max);
        };
        let Some(anchor) = allele.anchor.map(|a| a as u32) else {
            return LengthSupport::Full(requested_max);
        };
        let allele_len = allele.len();
        let current_5 = inst.trim_5 as u32;
        let current_3 = inst.trim_3 as u32;

        let max_admit = match target.end {
            TrimEnd::Five => {
                if anchor + 3 > allele_len.saturating_sub(current_3) {
                    return LengthSupport::Empty;
                }
                anchor
            }
            TrimEnd::Three => {
                if anchor < current_5 {
                    return LengthSupport::Empty;
                }
                allele_len.saturating_sub(anchor + 3)
            }
        };
        LengthSupport::Full(requested_max.min(max_admit))
    }

    /// Pinned non-canonical write: substitute `byte` at `site`,
    /// translate the resulting anchor codon, admit iff the amino
    /// acid is preserved.
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
        let Some(inst) = sim.assignments.get(self.segment).copied() else {
            return true;
        };
        let Some(allele) = refdata.get(self.segment, inst.allele_id) else {
            return true;
        };
        let Some(anchor) = allele.anchor.map(|a| a as u32) else {
            return true;
        };
        if self
            .require_anchor_retained(allele, inst.trim_5, inst.trim_3)
            .is_err()
        {
            return true;
        }
        let anchor_pool_start = match self.anchor_pool_start(sim, allele, anchor, inst.trim_5) {
            Ok(Some(start)) => start,
            Ok(None) | Err(_) => return true,
        };
        let target_idx = site.index();
        if target_idx < anchor_pool_start || target_idx >= anchor_pool_start + 3 {
            return true;
        }
        let reference_codon = Self::anchor_codon(allele, anchor);
        let live_codon = match self.live_anchor_codon(sim, anchor_pool_start) {
            Ok(codon) => codon,
            Err(_) => return true,
        };
        let codon_offset = (target_idx - anchor_pool_start) as usize;
        let mut hypothetical = live_codon;
        hypothetical[codon_offset] = byte;
        self.require_anchor_amino_acid_preserved(allele, reference_codon, hypothetical)
            .is_ok()
    }
}
