//! `AnchorPreserved` — anchor codon must remain in the retained slice.

use crate::ir::{Segment, Simulation};
use crate::refdata::RefDataConfig;

use super::{Contract, ContractViolation};

/// Verifies that the assigned allele's anchor codon (3 bases
/// starting at `allele.anchor`) sits entirely within the post-trim
/// retained slice for a given segment.
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
///   3'-trimmed away).
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
}

impl Contract for AnchorPreserved {
    fn name(&self) -> &str {
        Self::name_for(self.segment)
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

        Ok(())
    }
}

#[cfg(test)]
mod tests {
    use super::super::test_support::make_v_anchor_at;
    use super::*;
    use crate::assignment::AlleleInstance;
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
