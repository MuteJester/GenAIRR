//! `AnchorPreserved` — anchor codon must remain in the retained slice.

use crate::assignment::TrimEnd;
use crate::ir::{translate_codon, GermlinePos, NucHandle, Segment, Simulation};
use crate::refdata::{Allele, AlleleId, RefDataConfig};
use crate::trace::ChoiceValue;

use super::{ChoiceContext, ChoiceKind, Contract, ContractKind, ContractViolation};

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

    fn sample_address_for(segment: Segment) -> &'static str {
        match segment {
            Segment::V => "sample_allele.v",
            Segment::D => "sample_allele.d",
            Segment::J => "sample_allele.j",
            _ => unreachable!("AnchorPreserved with non-V/D/J segment"),
        }
    }

    fn trim_end_for_address(segment: Segment, address: &str) -> Option<TrimEnd> {
        match (segment, address) {
            (Segment::V, "trim.v_5") | (Segment::D, "trim.d_5") | (Segment::J, "trim.j_5") => {
                Some(TrimEnd::Five)
            }
            (Segment::V, "trim.v_3") | (Segment::D, "trim.d_3") | (Segment::J, "trim.j_3") => {
                Some(TrimEnd::Three)
            }
            _ => None,
        }
    }

    fn require_anchor_retained(
        &self,
        allele: &Allele,
        trim_5: u16,
        trim_3: u16,
    ) -> Result<(), ContractViolation> {
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

    fn anchor_codon(allele: &Allele, anchor: u32) -> [u8; 3] {
        [
            allele.seq[anchor as usize],
            allele.seq[anchor as usize + 1],
            allele.seq[anchor as usize + 2],
        ]
    }

    fn codon_amino_acid(codon: [u8; 3]) -> u8 {
        translate_codon(codon[0], codon[1], codon[2])
    }

    fn anchor_pool_start(
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

    fn live_anchor_codon(
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

    fn require_anchor_amino_acid_preserved(
        &self,
        allele: &Allele,
        reference_codon: [u8; 3],
        live_codon: [u8; 3],
    ) -> Result<(), ContractViolation> {
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

    fn admits_targeted_anchor_substitution(
        &self,
        sim: &Simulation,
        refdata: Option<&RefDataConfig>,
        address: &str,
        candidate: &ChoiceValue,
        context: ChoiceContext,
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

        let Some(anchor_pool_start) = self.anchor_pool_start(sim, allele, anchor, inst.trim_5)?
        else {
            return Ok(());
        };

        let mut live_codon = [b'N'; 3];
        for offset in 0..3 {
            let pool_pos = anchor_pool_start + offset;
            let expected_germline_pos = GermlinePos::pos((anchor + offset) as u16);
            let Some(nuc) = sim.pool.get(crate::ir::NucHandle::new(pool_pos)) else {
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
        self.require_anchor_amino_acid_preserved(allele, reference_codon, live_codon)?;

        Ok(())
    }

    fn admits(
        &self,
        sim: &Simulation,
        refdata: Option<&RefDataConfig>,
        address: &str,
        candidate: &ChoiceValue,
    ) -> Result<(), ContractViolation> {
        let Some(refdata) = refdata else {
            return Ok(());
        };

        if address == Self::sample_address_for(self.segment) {
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
            return self.require_anchor_retained(allele, 0, 0);
        }

        let Some(end) = Self::trim_end_for_address(self.segment, address) else {
            return Ok(());
        };
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

    fn admits_with_context(
        &self,
        sim: &Simulation,
        refdata: Option<&RefDataConfig>,
        address: &str,
        candidate: &ChoiceValue,
        context: ChoiceContext,
    ) -> Result<(), ContractViolation> {
        self.admits(sim, refdata, address, candidate)?;
        self.admits_targeted_anchor_substitution(sim, refdata, address, candidate, context)
    }
}

#[cfg(test)]
mod tests {
    use super::super::test_support::{
        make_assembled_sim_from_refdata, make_v_anchor_at, make_vj_with_anchor_codons,
    };
    use super::*;
    use crate::assignment::AlleleInstance;
    use crate::ir::{flag, Nucleotide};
    use crate::refdata::{Allele, AlleleId, ChainType, RefDataConfig};

    #[test]
    #[should_panic(expected = "segment must be V, D, or J")]
    fn anchor_preserved_rejects_np_segment() {
        let _ = AnchorPreserved::new(Segment::Np1);
    }

    #[test]
    fn anchor_preserved_name_per_segment() {
        assert_eq!(
            AnchorPreserved::new(Segment::V).name(),
            "anchor_preserved.v"
        );
        assert_eq!(
            AnchorPreserved::new(Segment::D).name(),
            "anchor_preserved.d"
        );
        assert_eq!(
            AnchorPreserved::new(Segment::J).name(),
            "anchor_preserved.j"
        );
    }

    #[test]
    fn anchor_preserved_no_allele_assigned_passes_vacuously() {
        let cfg = make_v_anchor_at(10, Some(5));
        let contract = AnchorPreserved::new(Segment::V);
        let sim = Simulation::new();

        // No assignment → vacuously satisfied.
        assert!(contract.verify(&sim, Some(&cfg)).is_ok());
    }

    #[test]
    fn anchor_preserved_no_refdata_passes_vacuously() {
        let cfg = make_v_anchor_at(10, Some(5));
        let contract = AnchorPreserved::new(Segment::V);
        let sim = Simulation::new()
            .with_allele_assigned(Segment::V, AlleleInstance::new(AlleleId::new(0)));

        // Refdata None → can't verify; treat as satisfied.
        assert!(contract.verify(&sim, None).is_ok());
        // Sanity: WITH refdata it would also pass (no trims).
        assert!(contract.verify(&sim, Some(&cfg)).is_ok());
    }

    #[test]
    fn anchor_preserved_anchorless_allele_passes_vacuously() {
        let cfg = make_v_anchor_at(10, None); // anchorless
        let contract = AnchorPreserved::new(Segment::V);
        let sim = Simulation::new()
            .with_allele_assigned(Segment::V, AlleleInstance::new(AlleleId::new(0)));

        assert!(contract.verify(&sim, Some(&cfg)).is_ok());
    }

    #[test]
    fn anchor_preserved_zero_trims_and_anchor_in_range_passes() {
        // 10-base allele, anchor at 3. No trims. Anchor codon spans
        // [3, 6) which fits entirely in [0, 10).
        let cfg = make_v_anchor_at(10, Some(3));
        let contract = AnchorPreserved::new(Segment::V);
        let sim = Simulation::new()
            .with_allele_assigned(Segment::V, AlleleInstance::new(AlleleId::new(0)));

        assert!(contract.verify(&sim, Some(&cfg)).is_ok());
    }

    #[test]
    fn anchor_preserved_rejects_live_anchor_provenance_shift_after_indel() {
        let cfg = make_vj_with_anchor_codons(b"GGG", b"TTT");
        let sim = make_assembled_sim_from_refdata(&cfg);
        let shifted = sim.with_indel_inserted(
            6,
            Nucleotide::synthetic(b'A', Segment::V, flag::INDEL_INSERTED),
        );

        let err = AnchorPreserved::new(Segment::V)
            .verify(&shifted, Some(&cfg))
            .unwrap_err();

        assert_eq!(err.contract_name, "anchor_preserved.v");
        assert!(
            err.reason.contains("provenance mismatch"),
            "reason should mention live provenance mismatch, got: {}",
            err.reason
        );
    }

    #[test]
    fn anchor_preserved_rejects_live_anchor_amino_acid_change() {
        let cfg = make_vj_with_anchor_codons(b"TGT", b"TGG");
        let sim = make_assembled_sim_from_refdata(&cfg);
        let changed = sim.with_base_changed(NucHandle::new(6), b'A');

        let err = AnchorPreserved::new(Segment::V)
            .verify(&changed, Some(&cfg))
            .unwrap_err();

        assert_eq!(err.contract_name, "anchor_preserved.v");
        assert!(
            err.reason.contains("anchor codon amino acid changed"),
            "reason should mention anchor amino-acid change, got: {}",
            err.reason
        );
    }

    #[test]
    fn anchor_preserved_allows_synonymous_anchor_base_change() {
        let cfg = make_vj_with_anchor_codons(b"TGT", b"TGG");
        let sim = make_assembled_sim_from_refdata(&cfg);
        let synonymous = sim.with_base_changed(NucHandle::new(2), b'C');

        assert!(AnchorPreserved::new(Segment::V)
            .verify(&synonymous, Some(&cfg))
            .is_ok());
    }

    #[test]
    fn anchor_preserved_rejects_targeted_anchor_substitution_candidate() {
        let cfg = make_vj_with_anchor_codons(b"TGT", b"TGG");
        let sim = make_assembled_sim_from_refdata(&cfg);
        let contract = AnchorPreserved::new(Segment::V);

        let err = contract
            .admits_with_context(
                &sim,
                Some(&cfg),
                "mutate.uniform.base[0]",
                &ChoiceValue::Base(b'A'),
                ChoiceContext::targeted_base_substitution(0, 1, NucHandle::new(6)),
            )
            .unwrap_err();

        assert_eq!(err.contract_name, "anchor_preserved.v");
        assert!(
            err.reason.contains("would change V anchor amino acid"),
            "reason should mention anchor amino-acid filtering, got: {}",
            err.reason
        );
    }

    #[test]
    fn anchor_preserved_violates_when_5prime_trim_eats_anchor() {
        // Anchor at 3. trim_5 = 4 → anchor 5'-trimmed away.
        use crate::assignment::TrimEnd;

        let cfg = make_v_anchor_at(10, Some(3));
        let contract = AnchorPreserved::new(Segment::V);
        let sim = Simulation::new()
            .with_allele_assigned(Segment::V, AlleleInstance::new(AlleleId::new(0)))
            .with_trim(Segment::V, TrimEnd::Five, 4);

        let err = contract.verify(&sim, Some(&cfg)).unwrap_err();
        assert_eq!(err.contract_name, "anchor_preserved.v");
        assert!(
            err.reason.contains("trim_5"),
            "reason should mention trim_5, got: {}",
            err.reason
        );
    }

    #[test]
    fn anchor_preserved_violates_when_3prime_trim_eats_anchor() {
        // 10-base allele, anchor at 6. Anchor codon spans [6, 9).
        // trim_3 = 2 → retained_end = 8, so anchor codon's last
        // position (8) is at the boundary — anchor + 3 = 9 > 8 → fail.
        use crate::assignment::TrimEnd;

        let cfg = make_v_anchor_at(10, Some(6));
        let contract = AnchorPreserved::new(Segment::V);
        let sim = Simulation::new()
            .with_allele_assigned(Segment::V, AlleleInstance::new(AlleleId::new(0)))
            .with_trim(Segment::V, TrimEnd::Three, 2);

        let err = contract.verify(&sim, Some(&cfg)).unwrap_err();
        assert_eq!(err.contract_name, "anchor_preserved.v");
        assert!(
            err.reason.contains("retained slice ends at"),
            "reason should describe retained slice, got: {}",
            err.reason
        );
    }

    #[test]
    fn anchor_preserved_anchor_at_5prime_boundary_passes() {
        // anchor at 0 → trim_5 = 0 leaves anchor at the very start.
        let cfg = make_v_anchor_at(10, Some(0));
        let contract = AnchorPreserved::new(Segment::V);
        let sim = Simulation::new()
            .with_allele_assigned(Segment::V, AlleleInstance::new(AlleleId::new(0)));

        assert!(contract.verify(&sim, Some(&cfg)).is_ok());
    }

    #[test]
    fn anchor_preserved_anchor_at_3prime_boundary_passes() {
        // 10-base allele, anchor at 7. Anchor codon spans [7, 10).
        // trim_3 = 0 → retained_end = 10. Boundary case — passes.
        let cfg = make_v_anchor_at(10, Some(7));
        let contract = AnchorPreserved::new(Segment::V);
        let sim = Simulation::new()
            .with_allele_assigned(Segment::V, AlleleInstance::new(AlleleId::new(0)));

        assert!(contract.verify(&sim, Some(&cfg)).is_ok());
    }

    #[test]
    fn anchor_preserved_one_more_trim_just_past_boundary_violates() {
        // Same as above but trim_3 = 1 → retained_end = 9, anchor
        // codon end (10) > 9 → fail.
        use crate::assignment::TrimEnd;

        let cfg = make_v_anchor_at(10, Some(7));
        let contract = AnchorPreserved::new(Segment::V);
        let sim = Simulation::new()
            .with_allele_assigned(Segment::V, AlleleInstance::new(AlleleId::new(0)))
            .with_trim(Segment::V, TrimEnd::Three, 1);

        assert!(contract.verify(&sim, Some(&cfg)).is_err());
    }

    #[test]
    fn anchor_preserved_works_for_j_segment_too() {
        // Build a J pool with an anchor at position 0 (typical J
        // anchor is near the 5' end). trim_5 = 0 → anchor preserved.
        // trim_5 = 1 → anchor 5'-trimmed → fail.
        use crate::assignment::TrimEnd;

        let mut cfg = RefDataConfig::empty(ChainType::Vdj);
        let _ = cfg.j_pool.push(Allele {
            name: "j_test*01".to_string(),
            gene: "j_test".to_string(),
            seq: vec![b'A'; 12],
            segment: Segment::J,
            anchor: Some(0),
        });
        let contract = AnchorPreserved::new(Segment::J);

        let ok = Simulation::new()
            .with_allele_assigned(Segment::J, AlleleInstance::new(AlleleId::new(0)));
        assert!(contract.verify(&ok, Some(&cfg)).is_ok());

        let bad = ok.with_trim(Segment::J, TrimEnd::Five, 1);
        assert!(contract.verify(&bad, Some(&cfg)).is_err());
    }

    #[test]
    fn contract_works_through_box_dyn() {
        // Trait-object usage is dyn-compatible.
        let cfg = make_v_anchor_at(10, Some(3));
        let contract: Box<dyn Contract> = Box::new(AnchorPreserved::new(Segment::V));
        let sim = Simulation::new()
            .with_allele_assigned(Segment::V, AlleleInstance::new(AlleleId::new(0)));

        assert!(contract.verify(&sim, Some(&cfg)).is_ok());
        assert_eq!(contract.name(), "anchor_preserved.v");
    }
}
