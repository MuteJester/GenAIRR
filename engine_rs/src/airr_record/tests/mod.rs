use crate::assignment::AlleleInstance;
use crate::ir::{NucHandle, Nucleotide, Region, Segment, Simulation};
use crate::pass::Outcome;
use crate::refdata::{Allele, AlleleId, ChainType, RefDataConfig};
use crate::trace::Trace;

mod anchors;
mod projection;
mod utilities;

fn anchor_record_fixture() -> (RefDataConfig, Simulation) {
    let mut cfg = RefDataConfig::empty(ChainType::Vj);
    let _ = cfg.v_pool.push(Allele {
        name: "v_airr*01".into(),
        gene: "v_airr".into(),
        seq: b"TGT".to_vec(),
        segment: Segment::V,
        anchor: Some(0),
    });

    let mut sim = Simulation::new();
    for (i, &b) in b"TGT".iter().enumerate() {
        let (next, _) = sim.with_nucleotide_pushed(Nucleotide::germline(b, i as u16, Segment::V));
        sim = next;
    }
    let region = Region::new(Segment::V, NucHandle::new(0), NucHandle::new(3));
    sim = sim.with_region_added(region).with_allele_assigned(
        Segment::V,
        crate::assignment::AlleleInstance::new(AlleleId::new(0)),
    );

    (cfg, sim)
}

fn outcome_from_sim(sim: Simulation) -> Outcome {
    Outcome {
        revisions: vec![sim],
        pass_names: Vec::new(),
        trace: Trace::new(),
        events: Vec::new(),
    }
}

fn call_projection_fixture() -> (RefDataConfig, Simulation, AlleleId, AlleleId) {
    let mut cfg = RefDataConfig::empty(ChainType::Vj);
    let v0 = cfg.v_pool.push(Allele {
        name: "IGHV1-1*01".into(),
        gene: "IGHV1-1".into(),
        seq: b"AAACCC".to_vec(),
        segment: Segment::V,
        anchor: None,
    });
    let v1 = cfg.v_pool.push(Allele {
        name: "IGHV1-1*02".into(),
        gene: "IGHV1-1".into(),
        seq: b"AAACCC".to_vec(),
        segment: Segment::V,
        anchor: None,
    });

    let mut sim = Simulation::new();
    for (i, &base) in b"AAACCC".iter().enumerate() {
        let (next, _) =
            sim.with_nucleotide_pushed(Nucleotide::germline(base, i as u16, Segment::V));
        sim = next;
    }
    let region = Region::new(Segment::V, NucHandle::new(0), NucHandle::new(6));
    sim = sim
        .with_region_added(region)
        .with_allele_assigned(Segment::V, AlleleInstance::new(v0));

    (cfg, sim, v0, v1)
}
