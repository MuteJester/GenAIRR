use super::*;
use crate::address::ChoiceAddress;
use crate::assignment::AlleleInstance;
use crate::contract::test_support::{make_assembled_sim_from_refdata, make_vj_with_anchor_codons};
use crate::ir::{NucHandle, Nucleotide, Region, Segment, Simulation};
use crate::refdata::{Allele, AlleleId, AllelePool, ChainType, RefDataConfig};

mod admits_fixed;
mod np_candidates;
mod targeted;
mod verify;

fn make_partial_np_stop_filter_case() -> (RefDataConfig, Simulation) {
    let mut cfg = RefDataConfig::empty(ChainType::Vj);
    let _ = cfg.v_pool.push(Allele {
        name: "v_test*01".into(),
        gene: "v_test".into(),
        seq: b"GGGTA".to_vec(),
        segment: Segment::V,
        anchor: Some(0),
    });
    let _ = cfg.j_pool.push(Allele {
        name: "j_test*01".into(),
        gene: "j_test".into(),
        seq: b"TTTAAA".to_vec(),
        segment: Segment::J,
        anchor: Some(0),
    });

    let mut sim = Simulation::new();
    for (i, &b) in b"GGGTA".iter().enumerate() {
        let (next, _) = sim.with_nucleotide_pushed(Nucleotide::germline(b, i as u16, Segment::V));
        sim = next;
    }
    sim = sim.with_region_added(Region::new(
        Segment::V,
        NucHandle::new(0),
        NucHandle::new(5),
    ));
    sim = sim
        .with_allele_assigned(Segment::V, AlleleInstance::new(AlleleId::new(0)))
        .with_allele_assigned(Segment::J, AlleleInstance::new(AlleleId::new(0)));

    (cfg, sim)
}

fn make_partial_np1_future_d_stop_case(v_seq: &[u8]) -> (RefDataConfig, Simulation) {
    let mut cfg = RefDataConfig::empty(ChainType::Vdj);
    let _ = cfg.v_pool.push(Allele {
        name: "v_test*01".into(),
        gene: "v_test".into(),
        seq: v_seq.to_vec(),
        segment: Segment::V,
        anchor: Some(0),
    });
    let _ = cfg.d_pool.push(Allele {
        name: "d_test*01".into(),
        gene: "d_test".into(),
        seq: b"AAC".to_vec(),
        segment: Segment::D,
        anchor: None,
    });
    let _ = cfg.j_pool.push(Allele {
        name: "j_test*01".into(),
        gene: "j_test".into(),
        seq: b"TTTAAA".to_vec(),
        segment: Segment::J,
        anchor: Some(0),
    });

    let mut sim = Simulation::new();
    for (i, &b) in v_seq.iter().enumerate() {
        let (next, _) = sim.with_nucleotide_pushed(Nucleotide::germline(b, i as u16, Segment::V));
        sim = next;
    }
    sim = sim.with_region_added(Region::new(
        Segment::V,
        NucHandle::new(0),
        NucHandle::new(v_seq.len() as u32),
    ));
    sim = sim
        .with_allele_assigned(Segment::V, AlleleInstance::new(AlleleId::new(0)))
        .with_allele_assigned(Segment::D, AlleleInstance::new(AlleleId::new(0)))
        .with_allele_assigned(Segment::J, AlleleInstance::new(AlleleId::new(0)));

    (cfg, sim)
}
