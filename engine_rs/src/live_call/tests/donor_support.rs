use super::*;
use crate::genotype::{GeneCopy, Genotype, Haplotype};
use crate::ir::SimulationBuilder;
use crate::live_call::{
    assembled_segment_live_call, check_segment_calls_parity, LiveCallSupport, ReferenceMatchIndex,
};
use crate::refdata::{ChainType, GeneId, RefDataConfig};
use std::sync::Arc;

fn support(hap0: &[u32], hap1: &[u32]) -> Arc<LiveCallSupport> {
    let hap = |ids: &[u32]| {
        let mut h = Haplotype::new();
        for &index in ids {
            h.set(
                Segment::V,
                GeneId::new(index),
                vec![GeneCopy {
                    allele: id(index),
                    copies: 1,
                    weight: 1.0,
                }],
            );
        }
        h
    };
    Arc::new(LiveCallSupport::from_genotype(&Genotype::new(
        [hap(hap0), hap(hap1)],
        [1.0, 0.0],
        None,
        "fixture".into(),
    )))
}

#[test]
fn donor_support_ranks_both_haplotypes_without_renumbering_alleles() {
    let mut cfg = RefDataConfig::empty(ChainType::Vj);
    for (name, seq) in [("V0", b"AAAA"), ("V1", b"CCCC"), ("V2", b"GGGG")] {
        let _ = cfg.v_pool.push(allele(Segment::V, name, seq));
    }
    let index = ReferenceMatchIndex::build(&cfg);
    let sim = simulation_with_region(Segment::V, b"CCGA", 0);
    let unrestricted = assembled_segment_live_call(&sim, &index, Segment::V, 1).unwrap();
    assert_eq!(unrestricted.allele_call.to_ids(), vec![id(1)]);
    let sim = sim.with_live_call_support(support(&[0], &[2]));
    let restricted = assembled_segment_live_call(&sim, &index, Segment::V, 1).unwrap();
    assert_eq!(restricted.allele_call.to_ids(), vec![id(0), id(2)]);
    // Haplotype 1 has zero generative weight, but its allele is still carried.
    assert_eq!(sim.live_call_alleles(Segment::V).unwrap(), &[id(0), id(2)]);

    let mut builder = SimulationBuilder::from_simulation(sim.clone());
    builder.attach_standard_mutation_observers(Some(&index));
    builder.change_base(NucHandle::new(0), b'G');
    let updated = builder.seal_with_committed_live_calls(&index);
    assert_eq!(
        updated
            .segment_calls
            .get(Segment::V)
            .unwrap()
            .allele_call
            .to_ids(),
        vec![id(2)]
    );
    for parity in check_segment_calls_parity(&updated, &cfg) {
        assert!(parity.tie_set_matches);
        assert_ne!(parity.hypothesis_bounds_match, Some(false));
    }
}

#[test]
fn donor_support_zero_score_fallback_excludes_uncarried_alleles() {
    let mut cfg = RefDataConfig::empty(ChainType::Vj);
    for (name, seq) in [("V0", b"AAAA"), ("V1", b"CCCC"), ("V2", b"GGGG")] {
        let _ = cfg.v_pool.push(allele(Segment::V, name, seq));
    }
    let index = ReferenceMatchIndex::build(&cfg);
    let sim =
        simulation_with_region(Segment::V, b"TTTT", 0).with_live_call_support(support(&[0], &[2]));
    let call = assembled_segment_live_call(&sim, &index, Segment::V, 1).unwrap();
    assert_eq!(call.allele_call.to_ids(), vec![id(0), id(2)]);
}

#[test]
fn donor_support_prevents_extensions_that_only_resolve_uncarried_alleles() {
    let mut cfg = RefDataConfig::empty(ChainType::Vj);
    let _ = cfg.v_pool.push(allele(Segment::V, "V0", b"AAAC"));
    let _ = cfg.v_pool.push(allele(Segment::V, "V1", b"AAAG"));
    let index = ReferenceMatchIndex::build(&cfg);
    let mut sim = simulation_with_region(Segment::V, b"AAA", 0).with_allele_assigned(
        Segment::V,
        crate::assignment::AlleleInstance::new(id(0)).with_trim_3(1),
    );
    let (next, _) = sim.with_nucleotide_pushed(Nucleotide::synthetic(
        b'G',
        Segment::Np1,
        crate::ir::NucFlags::default(),
    ));
    sim = next.with_region_added(Region::new(
        Segment::Np1,
        NucHandle::new(3),
        NucHandle::new(4),
    ));
    assert_eq!(
        assembled_segment_live_call(&sim, &index, Segment::V, 1)
            .unwrap()
            .hypotheses[0]
            .seq_end,
        4
    );
    let sim = sim.with_live_call_support(support(&[0], &[0]));
    let call = assembled_segment_live_call(&sim, &index, Segment::V, 1).unwrap();
    assert_eq!(call.allele_call.to_ids(), vec![id(0)]);
    assert_eq!(call.hypotheses[0].seq_end, 3);
    let scores = crate::live_call::scoring::score_alleles_with_extensions(
        &sim,
        Segment::V,
        cfg.v_pool.as_slice(),
        &sim.sequence.regions[0],
        0,
        1,
        crate::assignment::SegmentOrientation::Forward,
    );
    assert_eq!(
        scores,
        vec![3, 3],
        "independent oracle must also decline the extension"
    );
}
