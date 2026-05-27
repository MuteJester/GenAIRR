//! `AnchorPreserved` — anchor codon must remain in the retained slice.

use crate::address::{ChoiceAddress, VdjSegment};
use crate::assignment::TrimEnd;
use crate::ir::{NucHandle, Segment, Simulation};
use crate::refdata::{AlleleId, RefDataConfig};
use crate::trace::ChoiceValue;

use super::{ChoiceContext, Contract, ContractKind, ContractViolation};

mod admission;
mod live;
mod retention;
mod support;

/// Verifies that the assigned allele's anchor codon (3 bases
/// starting at `allele.anchor`) sits entirely within the post-trim
/// retained slice for a given segment and, when that segment has
/// already been assembled, that the live pool still carries those
/// three anchor-provenance nucleotides at the mapped anchor position,
/// and that the live anchor codon still translates to the same amino
/// acid as the assigned reference allele's anchor codon.
///
/// Vacuously satisfied (returns `Ok`) when:
/// - no allele is assigned to `segment` yet,
/// - no `RefDataConfig` is provided to look the allele up in,
/// - the allele has no anchor (`Allele::anchor == None` —
///   pseudogenes / partial alleles).
///
/// Returns a `ContractViolation` when:
/// - `trim_5 > allele.anchor` (anchor 5'-trimmed away), or
/// - `allele.anchor + 3 > allele.len() - trim_3` (anchor
///   3'-trimmed away), or
/// - a post-assembly structural edit has removed/shifted/replaced the
///   live anchor codon provenance, or
/// - a targeted substitution changes the anchor codon's amino acid.
pub struct AnchorPreserved {
    segment: Segment,
}

impl AnchorPreserved {
    /// Construct the contract for a given segment.
    /// Panics if `segment` is not V, D, or J.
    pub fn new(segment: Segment) -> Self {
        match segment {
            Segment::V | Segment::D | Segment::J => {}
            _ => panic!(
                "AnchorPreserved: segment must be V, D, or J — got {:?}",
                segment
            ),
        }
        Self { segment }
    }

    pub fn segment(&self) -> Segment {
        self.segment
    }

    fn name_for(segment: Segment) -> &'static str {
        match segment {
            Segment::V => "anchor_preserved.v",
            Segment::D => "anchor_preserved.d",
            Segment::J => "anchor_preserved.j",
            _ => unreachable!("AnchorPreserved with non-V/D/J segment"),
        }
    }

    fn typed_segment_matches(&self, segment: VdjSegment) -> bool {
        Segment::from(segment) == self.segment
    }

    fn admits_sample_allele_candidate(
        &self,
        refdata: &RefDataConfig,
        candidate: &ChoiceValue,
    ) -> Result<(), ContractViolation> {
        let ChoiceValue::AlleleId(raw_id) = candidate else {
            return Ok(());
        };
        let allele = refdata
            .get(self.segment, AlleleId::new(*raw_id))
            .ok_or_else(|| {
                ContractViolation::new(
                    self.name(),
                    format!(
                        "candidate allele id {} is absent from {:?} pool",
                        raw_id, self.segment
                    ),
                )
            })?;
        self.require_anchor_retained(allele, 0, 0)
    }

    fn admits_trim_candidate(
        &self,
        sim: &Simulation,
        refdata: &RefDataConfig,
        end: TrimEnd,
        candidate: &ChoiceValue,
    ) -> Result<(), ContractViolation> {
        let value = match candidate {
            ChoiceValue::Int(value) if (0..=u16::MAX as i64).contains(value) => *value as u16,
            _ => return Ok(()),
        };
        let Some(inst) = sim.assignments.get(self.segment).copied() else {
            return Ok(());
        };
        let Some(allele) = refdata.get(self.segment, inst.allele_id) else {
            return Ok(());
        };

        let (trim_5, trim_3) = match end {
            TrimEnd::Five => (value, inst.trim_3),
            TrimEnd::Three => (inst.trim_5, value),
        };
        self.require_anchor_retained(allele, trim_5, trim_3)
    }
}

impl Contract for AnchorPreserved {
    fn name(&self) -> &str {
        Self::name_for(self.segment)
    }

    fn kind(&self) -> ContractKind {
        ContractKind::AnchorPreserved {
            segment: self.segment,
        }
    }

    fn verify(
        &self,
        sim: &Simulation,
        refdata: Option<&RefDataConfig>,
    ) -> Result<(), ContractViolation> {
        // No allele yet → vacuously satisfied.
        let inst = match sim.assignments.get(self.segment) {
            None => return Ok(()),
            Some(i) => i,
        };

        // No refdata → can't verify; treat as satisfied.
        let refdata = match refdata {
            None => return Ok(()),
            Some(r) => r,
        };

        // Allele not in pool → can't verify; treat as satisfied.
        let allele = match refdata.get(self.segment, inst.allele_id) {
            None => return Ok(()),
            Some(a) => a,
        };

        // Anchorless allele → can't violate.
        let anchor = match allele.anchor {
            None => return Ok(()),
            Some(a) => a as u32,
        };

        let trim_5 = inst.trim_5 as u32;
        let trim_3 = inst.trim_3 as u32;
        let allele_len = allele.len();
        let retained_end = allele_len.saturating_sub(trim_3);

        // Anchor codon spans [anchor, anchor + 3). Must be fully
        // within the retained slice [trim_5, retained_end).
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
        if anchor + 3 > retained_end {
            return Err(ContractViolation::new(
                self.name(),
                format!(
                    "anchor codon at allele positions [{}, {}) but retained \
                     slice ends at {} (anchor 3'-trimmed from {} '{}')",
                    anchor,
                    anchor + 3,
                    retained_end,
                    allele.name,
                    allele.gene
                ),
            ));
        }

        self.verify_live_anchor(sim, allele, anchor, inst.trim_5)?;

        Ok(())
    }

    fn admits_typed(
        &self,
        sim: &Simulation,
        refdata: Option<&RefDataConfig>,
        context: ChoiceContext<'_>,
        candidate: &ChoiceValue,
    ) -> Result<(), ContractViolation> {
        let address = context.address_string().unwrap_or_default();
        if let Some(refdata) = refdata {
            match context.address {
                Some(ChoiceAddress::SampleAllele(segment))
                    if self.typed_segment_matches(segment) =>
                {
                    self.admits_sample_allele_candidate(refdata, candidate)?;
                }
                Some(ChoiceAddress::Trim { segment, end })
                    if self.typed_segment_matches(segment) =>
                {
                    self.admits_trim_candidate(sim, refdata, end, candidate)?;
                }
                _ => {}
            }
        }

        self.admits_targeted_anchor_substitution(sim, refdata, &address, candidate, context)
    }

    /// v3.0 constrain-before-propose: per-site admissible-base mask
    /// for the anchor codon.
    ///
    /// Mirrors `admits_targeted_anchor_substitution` semantics
    /// exactly:
    /// - No refdata / no allele assignment / no anchor on the allele
    ///   → all four bases admissible (this contract has nothing to
    ///   guard at this site).
    /// - Anchor not retained at all (trim already removes it) →
    ///   the existing predicate already returns an error; in the
    ///   mask context that means no admissible base at the anchor
    ///   site, so we return [`BaseMask::EMPTY`] when the query is
    ///   *inside* the (no-longer-present) anchor codon — defensive,
    ///   the pass should then skip. Sites outside the anchor codon
    ///   stay unconstrained.
    /// - Site outside the live anchor codon → unconstrained.
    /// - Site inside the live anchor codon → admit only those
    ///   candidate bases that preserve the reference anchor amino
    ///   acid.
    fn admissible_bases_at(
        &self,
        sim: &Simulation,
        refdata: Option<&RefDataConfig>,
        site: NucHandle,
    ) -> super::BaseMask {
        self.admissible_bases_at_impl(sim, refdata, site)
    }

    /// Classify an indel candidate's effect on this segment's
    /// anchor codon.
    ///
    /// Anchor preservation cares about whether the indel disturbs
    /// the **bytes** at `[anchor_pool, anchor_pool + 3)`. Under
    /// `with_indel_adjusted(pos, ±1)`, `region.start` shifts iff
    /// `region.start > pos` (strict). When the region doesn't
    /// shift but the indel lands inside the region, the new byte
    /// (or the deletion gap) corrupts the anchor codon bytes.
    ///
    /// So the anchor codon bytes are preserved iff:
    /// - Insertion at `s < region.start`: region (and anchor)
    ///   shift uniformly right → bytes unchanged → `FrameNeutral`.
    /// - Insertion at `s ≥ anchor_pool + 3`: anchor bytes unchanged
    ///   (the insertion lands strictly past the codon) →
    ///   `FrameNeutral`.
    /// - Insertion at `s ∈ [region.start, anchor_pool + 3)`:
    ///   region.start stays but a byte slips into the region
    ///   (and shifts the bytes within the anchor codon) → anchor
    ///   BROKEN → `Forbidden`.
    /// - Deletion at `s < region.start`: region (and anchor) shift
    ///   uniformly left → bytes unchanged → `FrameNeutral`.
    /// - Deletion at `s ≥ anchor_pool + 3`: anchor bytes unchanged
    ///   → `FrameNeutral`.
    /// - Deletion at `s ∈ [region.start, anchor_pool + 3)`:
    ///   region.start stays but the deletion removes a byte from
    ///   the anchor codon (or shifts subsequent bytes into it) →
    ///   `Forbidden`.
    ///
    /// Vacuous (`FrameNeutral`) when this segment isn't yet
    /// assembled (no region), the allele has no anchor, or the
    /// anchor has been trimmed away. Note: the contract only
    /// describes this segment's anchor; the bundle composes
    /// per-segment classifications via [`IndelEventClass::compose`].
    fn admissible_indel_class_at(
        &self,
        sim: &Simulation,
        refdata: Option<&RefDataConfig>,
        site: u32,
        kind: super::IndelKindHint,
    ) -> super::IndelEventClass {
        self.admissible_indel_class_at_impl(sim, refdata, site, kind)
    }

    /// Closed-form admissible trim-length support for this
    /// segment's anchor. A trim of length `L` at `target.end`
    /// replaces (not adds to) the current trim metadata, so the
    /// admissible bound is:
    ///   - `TrimEnd::Five`: `L ≤ anchor_offset` (the 5' trim
    ///     must not bite past the anchor base).
    ///   - `TrimEnd::Three`: `L ≤ allele_len − anchor_offset − 3`
    ///     (the 3' trim must not bite past the third anchor
    ///     codon byte).
    /// In addition, the *other* end's current trim must already
    /// satisfy anchor retention — otherwise the anchor is gone
    /// regardless of this candidate, and the support is empty.
    ///
    /// Vacuous (`Full(requested_max)`) when this segment isn't
    /// the trim target, no refdata is bound, no allele has been
    /// assigned, or the allele has no anchor — same shape as the
    /// other admit hooks.
    fn admissible_trim_lengths(
        &self,
        sim: &Simulation,
        refdata: Option<&RefDataConfig>,
        target: super::TrimTarget,
        requested_max: u32,
    ) -> super::LengthSupport {
        self.admissible_trim_lengths_impl(sim, refdata, target, requested_max)
    }

    /// Pinned non-canonical write: substitute `byte` at `site`,
    /// translate the resulting anchor codon, admit iff the amino
    /// acid is preserved. `translate_codon` is case-insensitive, so
    /// lowercase quality writes are evaluated identically to
    /// uppercase; non-canonical bytes (`N`, IUPAC ambiguities)
    /// translate to `X`, which by definition does not match the
    /// reference amino acid — so they are rejected at anchor sites.
    fn admits_fixed_base_at(
        &self,
        sim: &Simulation,
        refdata: Option<&RefDataConfig>,
        site: NucHandle,
        byte: u8,
    ) -> bool {
        self.admits_fixed_base_at_impl(sim, refdata, site, byte)
    }
}

#[cfg(test)]
mod tests;
