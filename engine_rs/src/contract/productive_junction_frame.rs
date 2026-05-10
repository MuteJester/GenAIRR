//! `ProductiveJunctionFrame` — junction length divisible by 3.

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

    fn admits(
        &self,
        sim: &Simulation,
        refdata: Option<&RefDataConfig>,
        address: &str,
        candidate: &ChoiceValue,
    ) -> Result<(), ContractViolation> {
        // Filter NP length samples to values that produce an in-frame
        // junction. The architectural pivot from D.4 → D.6: instead
        // of verifying after-the-fact, prune the candidate distribution
        // before sampling.
        //
        // We can compute the hypothetical junction length when:
        // - the NP being sampled is the LAST event before J assembly,
        // - we have V assembled (so V's pool position is known),
        // - we have V and J anchors known.
        //
        // For VJ chains: NP1 is the last event before J → filter np.np1.length.
        // For VDJ chains: NP2 is the last event before J → filter np.np2.length.
        // VDJ NP1 has D + NP2 + J between it and the junction end — NP2 will
        // compensate, so we don't filter NP1 in the VDJ case.
        //
        // **Pre-condition (asserted in debug builds):** the J region
        // must NOT be assembled yet. The hypothetical-J-start math
        // below assumes `sim.pool.len() + length` is where J will
        // start; that's only true if J hasn't already been pushed
        // into the pool. Standard plans honor this by always running
        // NP generation before J assembly. A future plan (e.g., one
        // that pushes contaminant or adapter bases into the pool
        // before J assembly) that violates this invariant would
        // silently produce wrong frame predictions.
        // The debug assertion catches that drift in dev builds; in
        // release builds the contract may produce an over-eager
        // rejection (false-negative admits) but never a hidden
        // soundness violation.
        let refdata = match refdata {
            None => return Ok(()),
            Some(r) => r,
        };

        let length = match candidate {
            ChoiceValue::Int(n) if *n >= 0 => *n as u32,
            _ => return Ok(()),
        };

        let is_vj = sim.assignments.d.is_none();
        let applicable = match address {
            "np.np1.length" => is_vj,
            "np.np2.length" => !is_vj,
            _ => false,
        };
        if !applicable {
            return Ok(());
        }

        // Need V/J alleles + anchors + V region in pool.
        let v_inst = match sim.assignments.v {
            None => return Ok(()),
            Some(v) => v,
        };
        let j_inst = match sim.assignments.j {
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

        // J must not yet be assembled — see method-level doc above.
        debug_assert!(
            sim.sequence.regions.iter().all(|r| r.segment != Segment::J),
            "ProductiveJunctionFrame::admits at {}: J region is already \
             assembled in sim.sequence.regions; the hypothetical-J-start \
             math below assumes J has not been added to the pool yet. \
             This indicates a pass-ordering bug — NP generation must run \
             before J assembly.",
            address
        );

        // Hypothetical J start: pool grows by `length` more bases
        // when the NP gets generated, then J starts immediately
        // after.
        let hypothetical_j_start = sim.pool.len() as u32 + length;
        let hypothetical_j_anchor_pool = hypothetical_j_start + (j_anchor - j_trim_5);
        // Defensive: junction must be a positive window.
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

#[cfg(test)]
mod tests {
    use super::super::test_support::{make_assembled_sim, make_vj_for_frame_test};
    use super::*;
    use crate::assignment::AlleleInstance;
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
    fn productive_junction_frame_works_through_box_dyn() {
        let cfg = make_vj_for_frame_test(Some(6), Some(0));
        let v_inst = AlleleInstance::new(AlleleId::new(0));
        let j_inst = AlleleInstance::new(AlleleId::new(0));
        let sim = make_assembled_sim(0, 9, 9, 6, v_inst, j_inst);

        let c: Box<dyn Contract> = Box::new(ProductiveJunctionFrame::new());
        assert!(c.verify(&sim, Some(&cfg)).is_ok());
        assert_eq!(c.name(), "productive_junction_frame");
    }
}
