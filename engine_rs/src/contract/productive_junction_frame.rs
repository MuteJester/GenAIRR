//! `ProductiveJunctionFrame` — junction length divisible by 3.

use crate::address::ChoiceAddress;
use crate::ir::{Segment, Simulation};
use crate::junction::compute_junction;
use crate::refdata::RefDataConfig;
use crate::trace::ChoiceValue;

use super::{Contract, ContractKind, ContractViolation};

/// Verifies that the junction (V Cys → J W/F + 3) has a length
/// divisible by 3 — i.e., the codon frame closes cleanly across
/// the V/NP/D/NP/J junction.
///
/// Vacuously satisfied (returns `Ok`) when:
/// - no `RefDataConfig` is provided to look anchors up in,
/// - the junction is undefined (no V or J assignment, anchorless
///   allele, V or J region not yet assembled, etc. — see
///   `compute_junction`).
///
/// Returns a `ContractViolation` when the junction is materialized
/// and its length is not divisible by 3.
///
/// **Note:** this contract checks frame *only*. Stop codons inside
/// an in-frame junction are the `NoStopCodonInJunction` contract's
/// concern (D.3). The `productive()` bundle (D.5) composes both.
pub struct ProductiveJunctionFrame;

impl ProductiveJunctionFrame {
    pub fn new() -> Self {
        Self
    }

    fn admits_np_length_candidate(
        &self,
        sim: &Simulation,
        refdata: Option<&RefDataConfig>,
        np_segment: Segment,
        address: &str,
        length: u32,
    ) -> Result<(), ContractViolation> {
        let refdata = match refdata {
            None => return Ok(()),
            Some(r) => r,
        };

        let is_vj = sim.assignments.get(Segment::D).is_none();
        let applicable = match np_segment {
            Segment::Np1 => is_vj,
            Segment::Np2 => !is_vj,
            _ => false,
        };
        if !applicable {
            return Ok(());
        }

        let v_inst = match sim.assignments.get(Segment::V).copied() {
            None => return Ok(()),
            Some(v) => v,
        };
        let j_inst = match sim.assignments.get(Segment::J).copied() {
            None => return Ok(()),
            Some(j) => j,
        };
        let v_allele = match refdata.get(Segment::V, v_inst.allele_id) {
            None => return Ok(()),
            Some(a) => a,
        };
        let j_allele = match refdata.get(Segment::J, j_inst.allele_id) {
            None => return Ok(()),
            Some(a) => a,
        };
        let v_anchor = match v_allele.anchor {
            None => return Ok(()),
            Some(a) => a as u32,
        };
        let j_anchor = match j_allele.anchor {
            None => return Ok(()),
            Some(a) => a as u32,
        };

        let v_trim_5 = v_inst.trim_5 as u32;
        let j_trim_5 = j_inst.trim_5 as u32;
        if v_trim_5 > v_anchor || j_trim_5 > j_anchor {
            return Ok(());
        }

        let v_region = match sim
            .sequence
            .regions
            .iter()
            .find(|r| r.segment == Segment::V)
        {
            None => return Ok(()),
            Some(r) => r,
        };
        let v_anchor_pool = v_region.start.index() + (v_anchor - v_trim_5);

        debug_assert!(
            sim.sequence.regions.iter().all(|r| r.segment != Segment::J),
            "ProductiveJunctionFrame::admits at {}: J region is already \
             assembled in sim.sequence.regions; the hypothetical-J-start \
             math below assumes J has not been added to the pool yet. \
             This indicates a pass-ordering bug — NP generation must run \
             before J assembly.",
            address
        );

        let hypothetical_j_start = sim.pool.len() as u32 + length;
        let hypothetical_j_anchor_pool = hypothetical_j_start + (j_anchor - j_trim_5);
        if hypothetical_j_anchor_pool + 3 <= v_anchor_pool {
            return Ok(());
        }
        let junction_length = (hypothetical_j_anchor_pool + 3) - v_anchor_pool;

        if junction_length % 3 == 0 {
            Ok(())
        } else {
            Err(ContractViolation::new(
                self.name(),
                format!(
                    "NP length {} at {} would produce out-of-frame \
                     junction (hypothetical length {})",
                    length, address, junction_length
                ),
            ))
        }
    }
}

impl Default for ProductiveJunctionFrame {
    fn default() -> Self {
        Self::new()
    }
}

impl Contract for ProductiveJunctionFrame {
    fn name(&self) -> &str {
        "productive_junction_frame"
    }

    fn kind(&self) -> ContractKind {
        ContractKind::ProductiveJunctionFrame
    }

    fn verify(
        &self,
        sim: &Simulation,
        refdata: Option<&RefDataConfig>,
    ) -> Result<(), ContractViolation> {
        let refdata = match refdata {
            None => return Ok(()),
            Some(r) => r,
        };

        let junction = match compute_junction(sim, refdata) {
            None => return Ok(()),
            Some(j) => j,
        };

        if junction.is_in_frame() {
            return Ok(());
        }

        Err(ContractViolation::new(
            self.name(),
            format!(
                "junction length {} is not divisible by 3 (start={}, end={})",
                junction.length,
                junction.start.index(),
                junction.end.index()
            ),
        ))
    }

    // `admits` and `admits_with_context` are intentionally NOT
    // overridden: the post-bridge-flip trait defaults parse the
    // legacy string address into the typed `ChoiceContext` and
    // route through `admits_typed` below, which dispatches on
    // `ChoiceAddress::NpLength(...)` directly. The string-prefix
    // match against `address::NP1_LENGTH` / `NP2_LENGTH` that the
    // old `admits` override carried is therefore dead — every
    // built-in caller (and any legacy string caller, via the
    // trait default) reaches the typed path.

    fn admits_typed(
        &self,
        sim: &Simulation,
        refdata: Option<&RefDataConfig>,
        context: super::ChoiceContext<'_>,
        candidate: &ChoiceValue,
    ) -> Result<(), ContractViolation> {
        let length = match candidate {
            ChoiceValue::Int(n) if *n >= 0 => *n as u32,
            _ => return Ok(()),
        };
        let Some(ChoiceAddress::NpLength(segment)) = context.address else {
            return Ok(());
        };
        let np_segment: Segment = segment.into();
        let address = context.address_string().unwrap_or_default();
        self.admits_np_length_candidate(sim, refdata, np_segment, &address, length)
    }

    /// Classify an indel candidate's effect on the junction frame.
    ///
    /// `compute_junction` derives the junction window as
    /// `[V_region.start + V_anchor_offset, J_region.start +
    /// J_anchor_offset + 3)`, so the junction *length* equals
    /// `(J_region.start - V_region.start) + (J_anchor_offset −
    /// V_anchor_offset) + 3`. The anchor offsets are allele-static —
    /// only the region starts shift under indels. Therefore an
    /// indel changes the length iff it shifts the V and J region
    /// starts asymmetrically.
    ///
    /// `Sequence::with_indel_adjusted(pos, ±1)` shifts `region.start`
    /// iff `region.start > pos` (strict — `region.start == pos`
    /// does NOT shift). So an indel at site `s`:
    /// - `s < V_region.start`: both V and J start shift → length
    ///   unchanged → `FrameNeutral`.
    /// - `V_region.start ≤ s < J_region.start`: V stays, J shifts
    ///   → length ± 1 → `FrameDelta(±1)` (sign follows kind).
    /// - `s ≥ J_region.start`: neither shifts → length unchanged →
    ///   `FrameNeutral`.
    ///
    /// Vacuous (`FrameNeutral`) whenever V or J isn't yet assembled
    /// (no region to read region.start from). The mod-3 DP in the
    /// indel pass aggregates this classification across the
    /// candidate tuple; final stop-codon / anchor-codon checks
    /// happen via `admits_post_event` on the assembled tuple.
    fn admissible_indel_class_at(
        &self,
        sim: &Simulation,
        _refdata: Option<&RefDataConfig>,
        site: u32,
        kind: super::IndelKindHint,
    ) -> super::IndelEventClass {
        use super::IndelEventClass;
        let v_start = match sim
            .sequence
            .regions
            .iter()
            .find(|r| r.segment == Segment::V)
        {
            Some(r) => r.start.index(),
            None => return IndelEventClass::FrameNeutral,
        };
        let j_start = match sim
            .sequence
            .regions
            .iter()
            .find(|r| r.segment == Segment::J)
        {
            Some(r) => r.start.index(),
            None => return IndelEventClass::FrameNeutral,
        };
        if j_start <= v_start {
            // Degenerate / not-yet-assembled layout; no junction
            // length to defend.
            return IndelEventClass::FrameNeutral;
        }
        if site >= v_start && site < j_start {
            match kind {
                super::IndelKindHint::Insertion => IndelEventClass::FrameDelta(1),
                super::IndelKindHint::Deletion => IndelEventClass::FrameDelta(-1),
            }
        } else {
            IndelEventClass::FrameNeutral
        }
    }
}

#[cfg(test)]
mod tests {
    use super::super::test_support::{make_assembled_sim, make_vj_for_frame_test};
    use super::*;
    use crate::address::{ChoiceAddress, NpSegment};
    use crate::assignment::AlleleInstance;
    use crate::contract::ChoiceContext;
    use crate::ir::{NucHandle, Nucleotide, Region};
    use crate::refdata::AlleleId;

    #[test]
    fn productive_junction_frame_name_is_canonical() {
        let c = ProductiveJunctionFrame::new();
        assert_eq!(c.name(), "productive_junction_frame");
    }

    #[test]
    fn productive_junction_frame_no_refdata_passes_vacuously() {
        let c = ProductiveJunctionFrame::new();
        let sim = Simulation::new();
        assert!(c.verify(&sim, None).is_ok());
    }

    #[test]
    fn productive_junction_frame_no_v_assignment_passes_vacuously() {
        let cfg = make_vj_for_frame_test(Some(6), Some(0));
        let c = ProductiveJunctionFrame::new();
        let sim = Simulation::new()
            .with_allele_assigned(Segment::J, AlleleInstance::new(AlleleId::new(0)));
        assert!(c.verify(&sim, Some(&cfg)).is_ok());
    }

    #[test]
    fn productive_junction_frame_anchorless_v_passes_vacuously() {
        // Junction undefined → contract has nothing to check.
        let cfg = make_vj_for_frame_test(None, Some(0));
        let c = ProductiveJunctionFrame::new();
        let v_inst = AlleleInstance::new(AlleleId::new(0));
        let j_inst = AlleleInstance::new(AlleleId::new(0));
        let sim = make_assembled_sim(0, 9, 9, 6, v_inst, j_inst);
        assert!(c.verify(&sim, Some(&cfg)).is_ok());
    }

    #[test]
    fn productive_junction_frame_in_frame_junction_passes() {
        // V anchor 6 in 9bp V → V anchor pool = 6.
        // J anchor 0 in 6bp J → J anchor pool = 9.
        // Junction = [6, 12). Length = 6 (divisible by 3).
        let cfg = make_vj_for_frame_test(Some(6), Some(0));
        let c = ProductiveJunctionFrame::new();
        let v_inst = AlleleInstance::new(AlleleId::new(0));
        let j_inst = AlleleInstance::new(AlleleId::new(0));
        let sim = make_assembled_sim(0, 9, 9, 6, v_inst, j_inst);
        assert!(c.verify(&sim, Some(&cfg)).is_ok());
    }

    #[test]
    fn productive_junction_frame_out_of_frame_junction_violates() {
        // V anchor 6, J anchor 2, J trim_5=1 → J anchor in trimmed J = 1.
        // J at pool [9, 14), J anchor pool = 9 + 1 = 10. End = 13.
        // V anchor pool = 6. Junction = [6, 13). Length 7 → not divisible by 3.
        let cfg = make_vj_for_frame_test(Some(6), Some(2));
        let c = ProductiveJunctionFrame::new();
        let v_inst = AlleleInstance::new(AlleleId::new(0));
        let j_inst = AlleleInstance::new(AlleleId::new(0)).with_trim_5(1);
        let sim = make_assembled_sim(0, 9, 9, 5, v_inst, j_inst);

        let err = c.verify(&sim, Some(&cfg)).unwrap_err();
        assert_eq!(err.contract_name, "productive_junction_frame");
        assert!(
            err.reason.contains("not divisible by 3"),
            "reason should mention frame: {}",
            err.reason
        );
        assert!(
            err.reason.contains("7"),
            "reason should mention length 7: {}",
            err.reason
        );
    }

    #[test]
    fn productive_junction_frame_in_frame_with_np_padding_passes() {
        // V at [0, 9). NP1 padding 3bp at [9, 12). J at [12, 18).
        // V anchor pool = 6. J anchor pool = 12. Junction = [6, 15). Length 9.
        let cfg = make_vj_for_frame_test(Some(6), Some(0));
        let c = ProductiveJunctionFrame::new();
        let v_inst = AlleleInstance::new(AlleleId::new(0));
        let j_inst = AlleleInstance::new(AlleleId::new(0));
        let sim = make_assembled_sim(0, 9, 12, 6, v_inst, j_inst);
        assert!(c.verify(&sim, Some(&cfg)).is_ok());
    }

    #[test]
    fn productive_junction_frame_typed_np_length_filters_by_frame() {
        // Pre-J assembly state: V is materialized, J is assigned but
        // not yet pushed. This is the exact state NP1 length sampling
        // filters in the VJ recombination plan.
        let cfg = make_vj_for_frame_test(Some(6), Some(0));
        let c = ProductiveJunctionFrame::new();

        let mut sim = Simulation::new();
        for i in 0..9u32 {
            let (next, _) =
                sim.with_nucleotide_pushed(Nucleotide::germline(b'A', i as u16, Segment::V));
            sim = next;
        }
        sim = sim
            .with_region_added(Region::new(
                Segment::V,
                NucHandle::new(0),
                NucHandle::new(9),
            ))
            .with_allele_assigned(Segment::V, AlleleInstance::new(AlleleId::new(0)))
            .with_allele_assigned(Segment::J, AlleleInstance::new(AlleleId::new(0)));

        // Junction-frame divisibility check is the only filter the
        // contract applies at NP1-length time: a length is admitted
        // iff (V_anchor_to_end + length + J_to_anchor) is divisible
        // by 3. For this fixture that means length % 3 == 0.
        for length in 0..6i64 {
            let candidate = ChoiceValue::Int(length);
            let admitted = c
                .admits_typed(
                    &sim,
                    Some(&cfg),
                    ChoiceContext::none().with_address(ChoiceAddress::NpLength(NpSegment::Np1)),
                    &candidate,
                )
                .is_ok();
            assert_eq!(admitted, length % 3 == 0, "length {length}");
        }
    }

    #[test]
    fn productive_junction_frame_typed_np2_length_filters_by_frame() {
        // VDJ branch: NP1 has D + NP2 + J between it and the
        // junction end, so the contract's frame filter applies to
        // NP2 (the LAST NP before J). The NP1 filter is a no-op
        // in this configuration.
        //
        // Set up the exact state NP2 length sampling sees: V
        // assembled, NP1 already padded into the pool, D assembled,
        // J assigned but not yet pushed.
        use crate::refdata::{Allele, ChainType, RefDataConfig};

        let mut cfg = RefDataConfig::empty(ChainType::Vdj);
        let _ = cfg.v_pool.push(Allele {
            name: "v_test*01".into(),
            gene: "v_test".into(),
            seq: b"AAACCCGGG".to_vec(),
            segment: Segment::V,
            anchor: Some(6),
            functional_status: None,
            subregions: Vec::new(),
        });
        let _ = cfg.d_pool.push(Allele {
            name: "d_test*01".into(),
            gene: "d_test".into(),
            seq: b"CCC".to_vec(),
            segment: Segment::D,
            anchor: None,
            functional_status: None,
            subregions: Vec::new(),
        });
        let _ = cfg.j_pool.push(Allele {
            name: "j_test*01".into(),
            gene: "j_test".into(),
            seq: b"TTTAAA".to_vec(),
            segment: Segment::J,
            anchor: Some(0),
            functional_status: None,
            subregions: Vec::new(),
        });
        let c = ProductiveJunctionFrame::new();

        // V[0..9) + NP1[9..12) + D[12..15); J not yet assembled.
        let mut sim = Simulation::new();
        for i in 0..9u32 {
            let (next, _) =
                sim.with_nucleotide_pushed(Nucleotide::germline(b'A', i as u16, Segment::V));
            sim = next;
        }
        for i in 0..3u32 {
            let (next, _) =
                sim.with_nucleotide_pushed(Nucleotide::germline(b'A', i as u16, Segment::Np1));
            sim = next;
        }
        for i in 0..3u32 {
            let (next, _) =
                sim.with_nucleotide_pushed(Nucleotide::germline(b'C', i as u16, Segment::D));
            sim = next;
        }
        sim = sim
            .with_region_added(Region::new(
                Segment::V,
                NucHandle::new(0),
                NucHandle::new(9),
            ))
            .with_region_added(Region::new(
                Segment::Np1,
                NucHandle::new(9),
                NucHandle::new(12),
            ))
            .with_region_added(Region::new(
                Segment::D,
                NucHandle::new(12),
                NucHandle::new(15),
            ))
            .with_allele_assigned(Segment::V, AlleleInstance::new(AlleleId::new(0)))
            .with_allele_assigned(Segment::D, AlleleInstance::new(AlleleId::new(0)))
            .with_allele_assigned(Segment::J, AlleleInstance::new(AlleleId::new(0)));

        // NP2 frame filter is the only meaningful constraint at NP2-
        // length time in VDJ. With V_anchor-to-end = 3, NP1 = 3,
        // D = 3, J_to_anchor = 0, total mod 3 is 0 + length mod 3;
        // so length % 3 == 0 is admissible.
        for length in 0..6i64 {
            let candidate = ChoiceValue::Int(length);
            let admitted = c
                .admits_typed(
                    &sim,
                    Some(&cfg),
                    ChoiceContext::none().with_address(ChoiceAddress::NpLength(NpSegment::Np2)),
                    &candidate,
                )
                .is_ok();
            assert_eq!(admitted, length % 3 == 0, "NP2 length {length}");
        }

        // The NP1 filter is a no-op in this VDJ configuration — the
        // contract intentionally defers to NP2 in
        // `admits_np_length_candidate`'s `applicable` check. Every
        // NP1 length is admitted vacuously.
        for length in 0..6i64 {
            let candidate = ChoiceValue::Int(length);
            assert!(
                c.admits_typed(
                    &sim,
                    Some(&cfg),
                    ChoiceContext::none().with_address(ChoiceAddress::NpLength(NpSegment::Np1)),
                    &candidate,
                )
                .is_ok(),
                "NP1 filter must be vacuous in VDJ at length {length}",
            );
        }
    }

    #[test]
    fn productive_junction_frame_works_through_box_dyn() {
        let cfg = make_vj_for_frame_test(Some(6), Some(0));
        let v_inst = AlleleInstance::new(AlleleId::new(0));
        let j_inst = AlleleInstance::new(AlleleId::new(0));
        let sim = make_assembled_sim(0, 9, 9, 6, v_inst, j_inst);

        let c: Box<dyn Contract> = Box::new(ProductiveJunctionFrame::new());
        assert!(c.verify(&sim, Some(&cfg)).is_ok());
        assert_eq!(c.name(), "productive_junction_frame");
    }

    mod admissible_indel_class_at {
        use super::*;
        use crate::contract::{IndelEventClass, IndelKindHint};

        fn fixture() -> (RefDataConfig, Simulation) {
            // V.region [0, 9), J.region [9, 15). Anchors at V[6] and
            // J[0] so junction = [6, 12), length 6.
            let cfg = make_vj_for_frame_test(Some(6), Some(0));
            let sim = make_assembled_sim(
                0,
                9,
                9,
                6,
                AlleleInstance::new(AlleleId::new(0)),
                AlleleInstance::new(AlleleId::new(0)),
            );
            (cfg, sim)
        }

        #[test]
        fn insertion_inside_v_region_is_frame_delta_plus_one() {
            // Any insertion at site in [V.start=0, J.start=9) shifts
            // J.region.start by +1 but leaves V.region.start fixed
            // → junction length += 1.
            let (cfg, sim) = fixture();
            let c = ProductiveJunctionFrame::new();
            for site in 0..9u32 {
                assert_eq!(
                    c.admissible_indel_class_at(&sim, Some(&cfg), site, IndelKindHint::Insertion,),
                    IndelEventClass::FrameDelta(1),
                    "site {} should be FrameDelta(+1) for insertion",
                    site
                );
            }
        }

        #[test]
        fn deletion_inside_v_region_is_frame_delta_minus_one() {
            let (cfg, sim) = fixture();
            let c = ProductiveJunctionFrame::new();
            for site in 0..9u32 {
                assert_eq!(
                    c.admissible_indel_class_at(&sim, Some(&cfg), site, IndelKindHint::Deletion),
                    IndelEventClass::FrameDelta(-1),
                    "site {} should be FrameDelta(-1) for deletion",
                    site
                );
            }
        }

        #[test]
        fn insertion_inside_j_region_is_frame_neutral() {
            // Site s ≥ J.start: neither V.start nor J.start shifts
            // → length unchanged.
            let (cfg, sim) = fixture();
            let c = ProductiveJunctionFrame::new();
            for site in 9..=15u32 {
                assert_eq!(
                    c.admissible_indel_class_at(&sim, Some(&cfg), site, IndelKindHint::Insertion,),
                    IndelEventClass::FrameNeutral,
                    "site {} should be FrameNeutral for insertion",
                    site
                );
            }
        }

        #[test]
        fn deletion_inside_j_region_is_frame_neutral() {
            let (cfg, sim) = fixture();
            let c = ProductiveJunctionFrame::new();
            for site in 9..15u32 {
                assert_eq!(
                    c.admissible_indel_class_at(&sim, Some(&cfg), site, IndelKindHint::Deletion),
                    IndelEventClass::FrameNeutral,
                    "site {} should be FrameNeutral for deletion",
                    site
                );
            }
        }

        #[test]
        fn unassembled_layout_is_frame_neutral() {
            let cfg = make_vj_for_frame_test(Some(6), Some(0));
            let c = ProductiveJunctionFrame::new();
            let sim = Simulation::new();
            // No V/J regions present → classification has nothing
            // to defend.
            for site in 0..10u32 {
                assert_eq!(
                    c.admissible_indel_class_at(&sim, Some(&cfg), site, IndelKindHint::Insertion,),
                    IndelEventClass::FrameNeutral
                );
            }
        }
    }
}
