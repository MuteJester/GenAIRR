use super::*;
use crate::contract::productive;
use crate::dist::{Distribution, EmpiricalLengthDist, FilteredSampleError, UniformBase};
use crate::ir::{NucHandle, Nucleotide, Region, Segment, Simulation};
use crate::pass::testing::PassRuntime;
use crate::pass::{PassError, PassPlan};
use crate::passes::test_support::make_substitution_productive_vj_fixture;
use crate::trace::ChoiceValue;

mod constraints;
mod core;
mod provenance;
mod replay;

/// Helper: build a sim with a 12-base germline V region. Used as
/// the canonical input for indel-pass tests so length/region
/// invariants are easy to assert against the post-pass state.
fn indel_test_sim() -> Simulation {
    let mut sim = Simulation::new();
    for (i, b) in b"AAACCCGGGTTT".iter().enumerate() {
        let (next, _) = sim.with_nucleotide_pushed(Nucleotide::germline(*b, i as u16, Segment::V));
        sim = next;
    }
    let region = Region::new(Segment::V, NucHandle::new(0), NucHandle::new(12));
    sim.with_region_added(region)
}

/// Helper: build a contiguous multi-segment simulation so indel
/// insertion provenance can be checked outside V.
fn multi_segment_indel_context_sim() -> Simulation {
    let layout = [
        (Segment::V, b"AAA".as_slice()),
        (Segment::Np1, b"CC".as_slice()),
        (Segment::J, b"GGG".as_slice()),
    ];

    let mut sim = Simulation::new();
    let mut regions = Vec::new();
    for (segment, bases) in layout {
        let start = sim.pool.len() as u32;
        for (i, b) in bases.iter().enumerate() {
            let (next, _) = sim.with_nucleotide_pushed(Nucleotide::germline(*b, i as u16, segment));
            sim = next;
        }
        let end = sim.pool.len() as u32;
        regions.push(Region::new(
            segment,
            NucHandle::new(start),
            NucHandle::new(end),
        ));
    }

    for region in regions {
        sim = sim.with_region_added(region);
    }
    sim
}

/// Helper: a count distribution that always returns the given
/// value. Built on `EmpiricalLengthDist` for ergonomic test setup.
fn fixed_count(n: i64) -> Box<dyn Distribution<Output = i64>> {
    Box::new(EmpiricalLengthDist::from_pairs(vec![(n, 1.0)]))
}

/// Smaller fixture sized for hand-enumerated distribution math:
/// V allele "TGTGGG" (6 bytes, anchor at 0 → V codon TGT = Cys),
/// J allele "TGGAAA" (6 bytes, anchor at 0 → J codon TGG = Trp).
/// Pool = V[0..6) + J[6..12); junction = [0..9), length 9.
///
/// Under `productive()` the per-site indel-class composition is:
///   - Sites 0,1,2 (V anchor codon): Forbidden.
///   - Sites 3,4,5 (V past anchor): FrameDelta(+1 ins / −1 del).
///   - Sites 6,7,8 (J anchor codon): Forbidden.
///   - Sites 9,10,11 (J past anchor): FrameNeutral.
///   - Site 12 (post-pool insertion): FrameNeutral.
///
/// This shape is the smallest fixture where the mod-3 DP has a
/// non-trivial mix of FrameNeutral and FrameDelta candidates,
/// suitable for exact-enumeration probability calculations.
pub(super) fn balanced_dp_distribution_fixture(
) -> (crate::refdata::RefDataConfig, crate::ir::Simulation) {
    use crate::assignment::AlleleInstance;
    use crate::ir::{Region, Segment};
    use crate::refdata::{Allele, AlleleId, ChainType, RefDataConfig};

    let mut cfg = RefDataConfig::empty(ChainType::Vdj);
    let _ = cfg.v_pool.push(Allele {
        name: "v_dp*01".into(),
        gene: "v_dp".into(),
        seq: b"TGTGGG".to_vec(),
        segment: Segment::V,
        anchor: Some(0),
        functional_status: None,
        subregions: Vec::new(),
    });
    let _ = cfg.j_pool.push(Allele {
        name: "j_dp*01".into(),
        gene: "j_dp".into(),
        seq: b"TGGAAA".to_vec(),
        segment: Segment::J,
        anchor: Some(0),
        functional_status: None,
        subregions: Vec::new(),
    });

    let mut sim = Simulation::new();
    for (i, &b) in b"TGTGGG".iter().enumerate() {
        let (next, _) = sim.with_nucleotide_pushed(Nucleotide::germline(b, i as u16, Segment::V));
        sim = next;
    }
    sim = sim.with_region_added(Region::new(
        Segment::V,
        NucHandle::new(0),
        NucHandle::new(6),
    ));
    for (i, &b) in b"TGGAAA".iter().enumerate() {
        let (next, _) = sim.with_nucleotide_pushed(Nucleotide::germline(b, i as u16, Segment::J));
        sim = next;
    }
    sim = sim.with_region_added(Region::new(
        Segment::J,
        NucHandle::new(6),
        NucHandle::new(12),
    ));
    sim = sim
        .with_allele_assigned(Segment::V, AlleleInstance::new(AlleleId::new(0)))
        .with_allele_assigned(Segment::J, AlleleInstance::new(AlleleId::new(0)));
    (cfg, sim)
}

/// Productive-compatible VJ fixture with enough sample-space for
/// the mod-3 DP to balance multi-event indel tuples.
///
/// Layout (using `make_vj_with_anchor_codons` so allele bytes
/// match pool bytes — AnchorPreserved verifies cleanly):
///   - V seq = "AAACCC" + V_anchor_codon, 9 bytes, anchor at 6.
///   - J seq = J_anchor_codon + "AAA", 6 bytes, anchor at 0.
///   - Pool = V[0..9) + J[9..15) — junction = [6, 12) = anchor codons.
///
/// V anchor = `TGT` (Cys), J anchor = `TGG` (Trp) — the canonical
/// IGH Cys-Trp bracket. Both have synonyms (TGC for Cys; W is
/// unique to TGG) so anchor preservation is satisfiable, and
/// neither codon is a stop.
///
/// FrameDelta sample space (under productive() bundle):
///   - Insertion FrameDelta(+1): sites s ∈ [V.start=0, J.start=9)
///     minus V anchor Forbidden range [0, anchor_end_V=9).
///     Net: empty — V[0..9) is entirely inside the anchor's
///     "Forbidden via shifted bytes" range. So FrameDelta only
///     viable through J: sites s ∈ ...wait, J.region.start=9 and
///     the J anchor codon is [9, 12), so J[12..15) is FrameNeutral.
///
/// Actually with V_anchor_offset=6, V's Forbidden range is
/// [V.start=0, V.start + anchor_offset + 3) = [0, 9) = entire V.
/// So FrameDelta sites within V are all Forbidden by the V anchor
/// contract. The only FrameDelta opportunity is the impossible
/// case (no overlap with V available without anchor Forbidden).
///
/// For the test to exercise +1/-1 pairing we'd need a V anchor at
/// offset 0 instead of 6 — then V[3..9) becomes FrameDelta(±1)
/// without anchor V Forbidden. Use V_anchor=0 + a J anchor codon
/// that still preserves the Cys-W junction shape: V allele =
/// `TGT` + "AAACCC" (anchor at 0), J = `TGG` + "AAA" (anchor at 0).
pub(super) fn make_productive_indel_balance_fixture(
) -> (crate::refdata::RefDataConfig, crate::ir::Simulation) {
    use crate::assignment::AlleleInstance;
    use crate::ir::{Region, Segment};
    use crate::refdata::{Allele, AlleleId, ChainType, RefDataConfig};

    let mut cfg = RefDataConfig::empty(ChainType::Vdj);
    let _ = cfg.v_pool.push(Allele {
        name: "v_balance*01".into(),
        gene: "v_balance".into(),
        seq: b"TGTAAACCC".to_vec(), // V_anchor=0 → codon TGT (Cys)
        segment: Segment::V,
        anchor: Some(0),
        functional_status: None,
        subregions: Vec::new(),
    });
    let _ = cfg.j_pool.push(Allele {
        name: "j_balance*01".into(),
        gene: "j_balance".into(),
        seq: b"TGGAAA".to_vec(), // J_anchor=0 → codon TGG (Trp)
        segment: Segment::J,
        anchor: Some(0),
        functional_status: None,
        subregions: Vec::new(),
    });

    let mut sim = Simulation::new();
    for (i, &b) in b"TGTAAACCC".iter().enumerate() {
        let (next, _) = sim.with_nucleotide_pushed(Nucleotide::germline(b, i as u16, Segment::V));
        sim = next;
    }
    sim = sim.with_region_added(Region::new(
        Segment::V,
        NucHandle::new(0),
        NucHandle::new(9),
    ));
    for (i, &b) in b"TGGAAA".iter().enumerate() {
        let (next, _) = sim.with_nucleotide_pushed(Nucleotide::germline(b, i as u16, Segment::J));
        sim = next;
    }
    sim = sim.with_region_added(Region::new(
        Segment::J,
        NucHandle::new(9),
        NucHandle::new(15),
    ));
    sim = sim
        .with_allele_assigned(Segment::V, AlleleInstance::new(AlleleId::new(0)))
        .with_allele_assigned(Segment::J, AlleleInstance::new(AlleleId::new(0)));

    (cfg, sim)
}
