use super::*;
use crate::assignment::AlleleInstance;
use crate::ir::{flag, NucHandle, Nucleotide, Region, Segment, Simulation};
use crate::refdata::{Allele, AlleleId, AllelePool, ChainType, RefDataConfig};

mod np1_frames;
mod preconditions;
mod static_stops;
mod vdj;

/// VJ chain refdata with a V allele "AAACCCxxx" (anchor at 6)
/// and J allele "yyyAAA" (anchor at 0).
fn make_vj_refdata(v_anchor_codon: &[u8; 3], j_anchor_codon: &[u8; 3]) -> RefDataConfig {
    let mut v_seq = b"AAACCC".to_vec();
    v_seq.extend_from_slice(v_anchor_codon);
    let mut j_seq = j_anchor_codon.to_vec();
    j_seq.extend_from_slice(b"AAA");

    let mut cfg = RefDataConfig::empty(ChainType::Vj);
    let _ = cfg.v_pool.push(Allele {
        name: "v_test*01".into(),
        gene: "v_test".into(),
        seq: v_seq,
        segment: Segment::V,
        anchor: Some(6),
        functional_status: None,
        subregions: Vec::new(),
    });
    let _ = cfg.j_pool.push(Allele {
        name: "j_test*01".into(),
        gene: "j_test".into(),
        seq: j_seq,
        segment: Segment::J,
        anchor: Some(0),
        functional_status: None,
        subregions: Vec::new(),
    });
    cfg
}

/// VDJ chain refdata. V "AAACCCxxx" (anchor 6), D "ddd", J
/// "yyyAAA" (anchor 0).
fn make_vdj_refdata(
    v_anchor_codon: &[u8; 3],
    d_body: &[u8],
    j_anchor_codon: &[u8; 3],
) -> RefDataConfig {
    let mut v_seq = b"AAACCC".to_vec();
    v_seq.extend_from_slice(v_anchor_codon);
    let mut j_seq = j_anchor_codon.to_vec();
    j_seq.extend_from_slice(b"AAA");

    let mut cfg = RefDataConfig::empty(ChainType::Vdj);
    let _ = cfg.v_pool.push(Allele {
        name: "v_test*01".into(),
        gene: "v_test".into(),
        seq: v_seq,
        segment: Segment::V,
        anchor: Some(6),
        functional_status: None,
        subregions: Vec::new(),
    });
    let _ = cfg.d_pool.push(Allele {
        name: "d_test*01".into(),
        gene: "d_test".into(),
        seq: d_body.to_vec(),
        segment: Segment::D,
        anchor: None,
        functional_status: None,
        subregions: Vec::new(),
    });
    let _ = cfg.j_pool.push(Allele {
        name: "j_test*01".into(),
        gene: "j_test".into(),
        seq: j_seq,
        segment: Segment::J,
        anchor: Some(0),
        functional_status: None,
        subregions: Vec::new(),
    });
    cfg
}

/// Build a sim with V already assembled in the pool, and the
/// V/J/(D) alleles assigned. NP not yet drawn.
fn make_sim_v_assembled(cfg: &RefDataConfig, has_d: bool) -> Simulation {
    let v_allele = cfg.v_pool.get(AlleleId::new(0)).unwrap();
    let mut sim = Simulation::new();
    for (i, &b) in v_allele.seq.iter().enumerate() {
        let (next, _) = sim.with_nucleotide_pushed(Nucleotide::germline(b, i as u16, Segment::V));
        sim = next;
    }
    sim = sim.with_region_added(Region::new(
        Segment::V,
        NucHandle::new(0),
        NucHandle::new(v_allele.seq.len() as u32),
    ));
    sim = sim.with_allele_assigned(Segment::V, AlleleInstance::new(AlleleId::new(0)));
    sim = sim.with_allele_assigned(Segment::J, AlleleInstance::new(AlleleId::new(0)));
    if has_d {
        sim = sim.with_allele_assigned(Segment::D, AlleleInstance::new(AlleleId::new(0)));
    }
    sim
}
