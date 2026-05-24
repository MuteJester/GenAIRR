use super::*;
use crate::contract::productive;
use crate::dist::{Distribution, EmpiricalLengthDist, FilteredSampleError, UniformBase};
use crate::ir::{NucHandle, Nucleotide, Region, Segment, Simulation};
use crate::pass::{PassError, PassPlan};
    use crate::pass::testing::PassRuntime;
use crate::passes::test_support::make_substitution_productive_vj_fixture;
use crate::trace::ChoiceValue;

mod constraints;
mod core;
mod provenance;

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
