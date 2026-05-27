use crate::ir::{NucHandle, Nucleotide, Region, Segment, Simulation};
use crate::refdata::{Allele, AlleleId};

mod bitset;
mod caller;
mod reference_index;
mod scoring;
mod state;
mod walker_observer_change_base;
mod walker_observer_from_existing_region;

fn id(index: u32) -> AlleleId {
    AlleleId::new(index)
}

fn allele(segment: Segment, name: &str, seq: &[u8]) -> Allele {
    Allele {
        name: name.to_string(),
        gene: name.split('*').next().unwrap_or(name).to_string(),
        seq: seq.to_vec(),
        segment,
        anchor: None,
    }
}

fn simulation_with_region(segment: Segment, bases: &[u8], ref_start: u16) -> Simulation {
    let mut sim = Simulation::new();
    for (offset, base) in bases.iter().enumerate() {
        let (next, _) = sim.with_nucleotide_pushed(Nucleotide::germline(
            *base,
            ref_start + offset as u16,
            segment,
        ));
        sim = next;
    }
    sim.with_region_added(Region::new(
        segment,
        NucHandle::new(0),
        NucHandle::new(bases.len() as u32),
    ))
}
