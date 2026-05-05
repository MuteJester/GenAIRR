//! Junction — the virtual entity from V Cys to J W/F + 3.
//!
//! ## What this is
//!
//! Per design doc §10.1, the **junction** is the canonical biological
//! window from the V conserved Cys codon (inclusive) through the J
//! conserved W/F codon (inclusive). It crosses the V/NP1/D/NP2/J
//! region boundaries, so it is not a `Region` of its own — it's
//! computed on demand from the assigned alleles + their assemblies.
//!
//! Productivity contracts in D.2 / D.3 read the junction. Future
//! AIRR projection (Phase F) will read it to populate the
//! `junction_*` fields of the output record.
//!
//! ## D.1 scope
//!
//! Just the data structure plus a pure function that materializes
//! the junction from a `Simulation` and a `RefDataConfig`. No
//! contracts yet — those land in D.2 / D.3.
//!
//! ## Why a function, not a method on `Simulation`
//!
//! Computing the junction requires both the IR (to know where
//! regions are in the pool) and the reference data (to know the
//! anchor positions on the original alleles). `Simulation` does not
//! own the reference data — it lives separately, threaded through
//! `PassContext` (D.8). A free function takes both and reads each
//! immutably.

use crate::ir::{NucHandle, Segment, Simulation};
use crate::refdata::RefDataConfig;

// ──────────────────────────────────────────────────────────────────
// Junction — the materialized junction window
// ──────────────────────────────────────────────────────────────────

/// The junction window inside a `Simulation`'s nucleotide pool.
///
/// `start` is the pool position of the V Cys's first base (the
/// anchor base in the V allele, mapped through the V region's
/// pool placement). `end` is `start + length` — i.e., one past the
/// last base of the J W/F codon. `length = end - start` is the
/// junction nucleotide count, the number whose divisibility by 3
/// determines productive frame.
#[derive(Copy, Clone, Eq, PartialEq, Debug)]
pub struct Junction {
    pub start: NucHandle,
    pub end: NucHandle,
    pub length: u32,
}

impl Junction {
    /// Whether the junction length is divisible by 3 (i.e., the
    /// productive frame condition is satisfied for this metric).
    /// Note: this does *not* check stop codons — that's the
    /// `NoStopCodonInJunction` contract's job (D.3).
    pub fn is_in_frame(&self) -> bool {
        self.length % 3 == 0
    }
}

// ──────────────────────────────────────────────────────────────────
// compute_junction — materialize the junction or return None
// ──────────────────────────────────────────────────────────────────

/// Compute the junction window from the simulation IR and reference
/// data, or return `None` when the junction is undefined.
///
/// Returns `None` when:
/// - V or J assignments are missing (pre-recombination),
/// - the assigned V or J allele is anchorless,
/// - the V or J allele cannot be looked up in `refdata` (wrong
///   refdata for this simulation — caller bug),
/// - V trim_5 has eaten past V's anchor (anchor 5'-trimmed),
/// - J trim_5 has eaten past J's anchor (J anchor 5'-trimmed),
/// - the V region (or J region) has not been assembled yet (no
///   `Region` for that segment in `sim.sequence.regions`).
///
/// A `None` result is informational, not an error — it just means
/// "no junction at this state." Contracts handle `None` according
/// to their own semantics (typically: vacuously satisfied, since
/// there's nothing to verify yet).
pub fn compute_junction(sim: &Simulation, refdata: &RefDataConfig) -> Option<Junction> {
    let v_inst = sim.assignments.v?;
    let j_inst = sim.assignments.j?;

    let v_allele = refdata.get(Segment::V, v_inst.allele_id)?;
    let j_allele = refdata.get(Segment::J, j_inst.allele_id)?;

    let v_anchor = v_allele.anchor? as u32;
    let j_anchor = j_allele.anchor? as u32;

    let v_trim_5 = v_inst.trim_5 as u32;
    let j_trim_5 = j_inst.trim_5 as u32;

    // Anchor must lie within the retained slice — i.e., trim_5 must
    // not have crossed the anchor base. (The full anchor codon
    // position is checked by `AnchorPreserved`; here we just check
    // that the start of the codon is reachable.)
    if v_trim_5 > v_anchor {
        return None;
    }
    if j_trim_5 > j_anchor {
        return None;
    }

    // Find the V and J regions in the sequence's region list.
    // Regions are appended in assembly order; either may be absent
    // if the corresponding `AssembleSegmentPass` hasn't run.
    let v_region = sim
        .sequence
        .regions
        .iter()
        .find(|r| r.segment == Segment::V)?;
    let j_region = sim
        .sequence
        .regions
        .iter()
        .find(|r| r.segment == Segment::J)?;

    // Map allele-coordinate anchors into pool coordinates.
    //
    //   v_anchor_in_trimmed_v = v_anchor - v_trim_5
    //   v_anchor_in_pool       = v_region.start + above
    //
    //   j_anchor_in_trimmed_j  = j_anchor - j_trim_5
    //   j_anchor_in_pool       = j_region.start + above
    //
    // The j_region.start uses the assembled J's pool position, so
    // the chain of intervening regions (NP1, D, NP2) is implicitly
    // accounted for.
    let v_anchor_in_pool = v_region.start.index() + (v_anchor - v_trim_5);
    let j_anchor_in_pool = j_region.start.index() + (j_anchor - j_trim_5);

    let start = NucHandle::new(v_anchor_in_pool);
    let end = NucHandle::new(j_anchor_in_pool + 3);

    // Defensive: end must be strictly greater than start. If the
    // J anchor maps to a position before the V anchor, the simulation
    // is malformed.
    if end.index() <= start.index() {
        return None;
    }

    let length = end.index() - start.index();

    Some(Junction { start, end, length })
}

// ──────────────────────────────────────────────────────────────────
// Tests
// ──────────────────────────────────────────────────────────────────

#[cfg(test)]
mod tests {
    use super::*;
    use crate::assignment::AlleleInstance;
    use crate::ir::{Nucleotide, Region};
    use crate::refdata::{Allele, AlleleId, ChainType};

    /// Build a minimal V+J refdata with controllable anchor positions.
    /// V allele "AAACCCGGG" (9bp) anchor at v_anchor.
    /// J allele "TTTAAA" (6bp) anchor at j_anchor.
    fn make_vj_refdata(v_anchor: Option<u16>, j_anchor: Option<u16>) -> RefDataConfig {
        let mut cfg = RefDataConfig::empty(ChainType::Vdj);
        let _ = cfg.v_pool.push(Allele {
            name: "v_test*01".into(),
            gene: "v_test".into(),
            seq: b"AAACCCGGG".to_vec(),
            segment: Segment::V,
            anchor: v_anchor,
        });
        let _ = cfg.j_pool.push(Allele {
            name: "j_test*01".into(),
            gene: "j_test".into(),
            seq: b"TTTAAA".to_vec(),
            segment: Segment::J,
            anchor: j_anchor,
        });
        cfg
    }

    /// Build a simulation that has V and J assigned + V/J regions
    /// in the sequence at given pool offsets. No NP1/D/NP2 — just
    /// V then J. Returns the simulation.
    fn make_sim_with_v_and_j_assemblies(
        v_pool_start: u32,
        v_len: u32,
        j_pool_start: u32,
        j_len: u32,
        v_inst: AlleleInstance,
        j_inst: AlleleInstance,
    ) -> Simulation {
        let mut sim = Simulation::new();

        // Pad pool with placeholder bases up to v_pool_start.
        for i in 0..v_pool_start {
            let (next, _) =
                sim.with_nucleotide_pushed(Nucleotide::germline(b'X', i as u16, Segment::V));
            sim = next;
        }
        for i in 0..v_len {
            let (next, _) =
                sim.with_nucleotide_pushed(Nucleotide::germline(b'A', i as u16, Segment::V));
            sim = next;
        }

        // J pool padding (between V end and J start).
        for _ in v_pool_start + v_len..j_pool_start {
            let (next, _) = sim.with_nucleotide_pushed(Nucleotide::synthetic(
                b'a',
                Segment::Np1,
                crate::ir::flag::N_NUC,
            ));
            sim = next;
        }
        for i in 0..j_len {
            let (next, _) =
                sim.with_nucleotide_pushed(Nucleotide::germline(b'T', i as u16, Segment::J));
            sim = next;
        }

        // V region.
        sim = sim.with_region_added(Region::new(
            Segment::V,
            NucHandle::new(v_pool_start),
            NucHandle::new(v_pool_start + v_len),
        ));
        // J region.
        sim = sim.with_region_added(Region::new(
            Segment::J,
            NucHandle::new(j_pool_start),
            NucHandle::new(j_pool_start + j_len),
        ));

        sim.with_allele_assigned(Segment::V, v_inst)
            .with_allele_assigned(Segment::J, j_inst)
    }

    #[test]
    fn junction_in_frame_predicate() {
        let j = Junction {
            start: NucHandle::new(0),
            end: NucHandle::new(9),
            length: 9,
        };
        assert!(j.is_in_frame());

        let j = Junction {
            start: NucHandle::new(0),
            end: NucHandle::new(10),
            length: 10,
        };
        assert!(!j.is_in_frame());
    }

    #[test]
    fn compute_junction_returns_none_without_v_assignment() {
        let cfg = make_vj_refdata(Some(6), Some(0));
        let sim = Simulation::new()
            .with_allele_assigned(Segment::J, AlleleInstance::new(AlleleId::new(0)));
        assert!(compute_junction(&sim, &cfg).is_none());
    }

    #[test]
    fn compute_junction_returns_none_without_j_assignment() {
        let cfg = make_vj_refdata(Some(6), Some(0));
        let sim = Simulation::new()
            .with_allele_assigned(Segment::V, AlleleInstance::new(AlleleId::new(0)));
        assert!(compute_junction(&sim, &cfg).is_none());
    }

    #[test]
    fn compute_junction_returns_none_for_anchorless_v() {
        let cfg = make_vj_refdata(None, Some(0));
        let v_inst = AlleleInstance::new(AlleleId::new(0));
        let j_inst = AlleleInstance::new(AlleleId::new(0));
        let sim = make_sim_with_v_and_j_assemblies(0, 9, 9, 6, v_inst, j_inst);
        assert!(compute_junction(&sim, &cfg).is_none());
    }

    #[test]
    fn compute_junction_returns_none_for_anchorless_j() {
        let cfg = make_vj_refdata(Some(6), None);
        let v_inst = AlleleInstance::new(AlleleId::new(0));
        let j_inst = AlleleInstance::new(AlleleId::new(0));
        let sim = make_sim_with_v_and_j_assemblies(0, 9, 9, 6, v_inst, j_inst);
        assert!(compute_junction(&sim, &cfg).is_none());
    }

    #[test]
    fn compute_junction_returns_none_when_v_trim_5_past_anchor() {
        let cfg = make_vj_refdata(Some(6), Some(0));
        let v_inst = AlleleInstance::new(AlleleId::new(0)).with_trim_5(7);
        let j_inst = AlleleInstance::new(AlleleId::new(0));
        let sim = make_sim_with_v_and_j_assemblies(0, 9, 9, 6, v_inst, j_inst);
        assert!(compute_junction(&sim, &cfg).is_none());
    }

    #[test]
    fn compute_junction_returns_none_when_j_region_missing() {
        let cfg = make_vj_refdata(Some(6), Some(0));
        let v_inst = AlleleInstance::new(AlleleId::new(0));
        let j_inst = AlleleInstance::new(AlleleId::new(0));
        let mut sim = Simulation::new();
        // V region exists, J doesn't.
        sim = sim.with_region_added(Region::new(
            Segment::V,
            NucHandle::new(0),
            NucHandle::new(9),
        ));
        sim = sim
            .with_allele_assigned(Segment::V, v_inst)
            .with_allele_assigned(Segment::J, j_inst);

        assert!(compute_junction(&sim, &cfg).is_none());
    }

    #[test]
    fn compute_junction_basic_case_no_trims_no_np() {
        // V len 9, J len 6. V at [0, 9), J at [9, 15).
        // V anchor 6 in V allele → V anchor in pool = 0 + 6 = 6.
        // J anchor 0 in J allele → J anchor in pool = 9 + 0 = 9.
        // Junction = [6, 9 + 3) = [6, 12). Length = 6.
        let cfg = make_vj_refdata(Some(6), Some(0));
        let v_inst = AlleleInstance::new(AlleleId::new(0));
        let j_inst = AlleleInstance::new(AlleleId::new(0));
        let sim = make_sim_with_v_and_j_assemblies(0, 9, 9, 6, v_inst, j_inst);

        let j = compute_junction(&sim, &cfg).expect("junction defined");
        assert_eq!(j.start.index(), 6);
        assert_eq!(j.end.index(), 12);
        assert_eq!(j.length, 6);
        assert!(j.is_in_frame());
    }

    #[test]
    fn compute_junction_accounts_for_np1_pool_padding() {
        // V at [0, 9), 3-base NP1 padding at [9, 12), J at [12, 18).
        // V anchor 6 → pool 6. J anchor 0 → pool 12. Junction len = 12+3-6 = 9.
        let cfg = make_vj_refdata(Some(6), Some(0));
        let v_inst = AlleleInstance::new(AlleleId::new(0));
        let j_inst = AlleleInstance::new(AlleleId::new(0));
        let sim = make_sim_with_v_and_j_assemblies(0, 9, 12, 6, v_inst, j_inst);

        let j = compute_junction(&sim, &cfg).expect("junction defined");
        assert_eq!(j.start.index(), 6);
        assert_eq!(j.end.index(), 15);
        assert_eq!(j.length, 9);
        assert!(j.is_in_frame());
    }

    #[test]
    fn compute_junction_with_v_trim_5() {
        // V trim_5 = 2 → V anchor 6 maps to anchor-trim_5 = 4 in
        // trimmed V. V region starts at pool 0, so V anchor pool = 4.
        // J anchor 0, no trim. J at pool [7, 13). J anchor pool = 7.
        // Junction = [4, 7+3) = [4, 10). Length = 6.
        let cfg = make_vj_refdata(Some(6), Some(0));
        let v_inst = AlleleInstance::new(AlleleId::new(0)).with_trim_5(2);
        let j_inst = AlleleInstance::new(AlleleId::new(0));
        // V len after trim_5=2 is 7.
        let sim = make_sim_with_v_and_j_assemblies(0, 7, 7, 6, v_inst, j_inst);

        let j = compute_junction(&sim, &cfg).expect("junction defined");
        assert_eq!(j.start.index(), 4);
        assert_eq!(j.end.index(), 10);
        assert_eq!(j.length, 6);
    }

    #[test]
    fn compute_junction_with_j_trim_5() {
        // J anchor 2, J trim_5 = 1 → anchor-trim_5 = 1 in trimmed J.
        // J at pool [9, 14) (J len 5 after trim_5=1).
        // J anchor pool = 9 + 1 = 10. Junction end = 10 + 3 = 13.
        // V anchor pool = 6, no trim. Junction = [6, 13). Length = 7.
        let cfg = make_vj_refdata(Some(6), Some(2));
        let v_inst = AlleleInstance::new(AlleleId::new(0));
        let j_inst = AlleleInstance::new(AlleleId::new(0)).with_trim_5(1);
        let sim = make_sim_with_v_and_j_assemblies(0, 9, 9, 5, v_inst, j_inst);

        let j = compute_junction(&sim, &cfg).expect("junction defined");
        assert_eq!(j.start.index(), 6);
        assert_eq!(j.end.index(), 13);
        assert_eq!(j.length, 7);
        assert!(!j.is_in_frame());
    }

    #[test]
    fn compute_junction_in_frame_predicate_works_per_test_case() {
        // length 6 → in frame
        let cfg = make_vj_refdata(Some(6), Some(0));
        let v_inst = AlleleInstance::new(AlleleId::new(0));
        let j_inst = AlleleInstance::new(AlleleId::new(0));
        let sim = make_sim_with_v_and_j_assemblies(0, 9, 9, 6, v_inst, j_inst);
        assert!(compute_junction(&sim, &cfg).unwrap().is_in_frame());

        // length 7 → out of frame
        let cfg = make_vj_refdata(Some(6), Some(2));
        let v_inst = AlleleInstance::new(AlleleId::new(0));
        let j_inst = AlleleInstance::new(AlleleId::new(0)).with_trim_5(1);
        let sim = make_sim_with_v_and_j_assemblies(0, 9, 9, 5, v_inst, j_inst);
        assert!(!compute_junction(&sim, &cfg).unwrap().is_in_frame());
    }

    #[test]
    fn compute_junction_extent_stable_across_pool_revisions() {
        // Audit finding: the junction extent depends only on V/J
        // anchors + trims + region boundaries, NOT on individual
        // base content. Therefore mutating a base inside the
        // junction window must NOT change `compute_junction(...)`
        // — the junction is a coordinate window, not a content
        // window.
        //
        // This invariant is what makes constraint-aware sampling
        // stable through SHM in Phase E: a mutation pass that
        // changes bases inside the junction shouldn't shift the
        // junction's position.
        let cfg = make_vj_refdata(Some(6), Some(0));
        let v_inst = AlleleInstance::new(AlleleId::new(0));
        let j_inst = AlleleInstance::new(AlleleId::new(0));
        let sim_a = make_sim_with_v_and_j_assemblies(0, 9, 9, 6, v_inst, j_inst);
        let junction_a = compute_junction(&sim_a, &cfg).expect("junction defined for sim_a");

        // Mutate position 7 (inside the junction window [6, 12)).
        let sim_b = sim_a.with_base_changed(NucHandle::new(7), b'T');
        let junction_b = compute_junction(&sim_b, &cfg).expect("junction defined for sim_b");

        assert_eq!(junction_a.start, junction_b.start);
        assert_eq!(junction_a.end, junction_b.end);
        assert_eq!(junction_a.length, junction_b.length);

        // Mutate position 4 (V region, OUTSIDE the junction at [6, 12)).
        let sim_c = sim_a.with_base_changed(NucHandle::new(4), b'G');
        let junction_c = compute_junction(&sim_c, &cfg).expect("junction defined for sim_c");
        assert_eq!(junction_a.start, junction_c.start);
        assert_eq!(junction_a.end, junction_c.end);
    }

    #[test]
    fn compute_junction_v_anchor_at_zero() {
        // V anchor at 0. With trim_5 = 0, anchor pool = 0.
        // J anchor at 0 → J anchor pool = j_pool_start.
        let cfg = make_vj_refdata(Some(0), Some(0));
        let v_inst = AlleleInstance::new(AlleleId::new(0));
        let j_inst = AlleleInstance::new(AlleleId::new(0));
        let sim = make_sim_with_v_and_j_assemblies(0, 9, 9, 6, v_inst, j_inst);

        let j = compute_junction(&sim, &cfg).expect("junction defined");
        assert_eq!(j.start.index(), 0);
        assert_eq!(j.end.index(), 12);
    }
}
