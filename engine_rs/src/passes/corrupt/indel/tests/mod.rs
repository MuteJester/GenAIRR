use super::*;
use crate::contract::productive;
use crate::dist::{Distribution, EmpiricalLengthDist, FilteredSampleError, UniformBase};
use crate::ir::{NucHandle, Nucleotide, Region, Segment, Simulation};
use crate::pass::testing::PassRuntime;
use crate::pass::{PassError, PassPlan};
use crate::trace::ChoiceValue;

mod constraints;
mod core;
mod natural;
mod provenance;
mod replay;

fn make_substitution_productive_vj_fixture() -> (crate::refdata::RefDataConfig, Simulation) {
    use crate::assignment::AlleleInstance;
    use crate::refdata::{Allele, ChainType, RefDataConfig};
    let mut cfg = RefDataConfig::empty(ChainType::Vj);
    let mut sim = Simulation::new();
    for (segment, seq) in [(Segment::V, b"TGT"), (Segment::J, b"TGG")] {
        let allele = Allele {
            name: format!("{segment:?}*01"),
            gene: format!("{segment:?}"),
            seq: seq.to_vec(),
            segment,
            anchor: Some(0),
            functional_status: None,
            subregions: Vec::new(),
        };
        let id = if segment == Segment::V {
            cfg.v_pool.push(allele)
        } else {
            cfg.j_pool.push(allele)
        };
        sim = sim.with_allele_assigned(segment, AlleleInstance::new(id));
        let start = sim.pool.len() as u32;
        for (g, b) in seq.iter().enumerate() {
            sim = sim
                .with_nucleotide_pushed(Nucleotide::germline(*b, g as u16, segment))
                .0;
        }
        sim = sim.with_region_added(Region::new(
            segment,
            NucHandle::new(start),
            NucHandle::new(start + 3),
        ));
    }
    (cfg, sim)
}

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
/// Full-tuple acceptance checks original anchor provenance, observed coding
/// phase, coding stops, and the retained J endpoint. Temporary insertions
/// inside anchors may be repaired later. The short sequence makes independent
/// enumeration of every two-event history practical.
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

/// Productive-compatible VJ fixture with enough support for compensating
/// indel tuples: V=TGTAAACCC at [0..9), J=TGGAAA at [9..15).
/// The V Cys anchor is [0..3), the J Trp anchor is [9..12), and
/// sites between them admit junction-frame deltas of ±1. The complete
/// productive bundle additionally checks the coding frame and stops.
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
