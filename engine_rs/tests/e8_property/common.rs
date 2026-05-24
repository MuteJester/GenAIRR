use genairr_engine::dist::{AllelePoolDist, EmpiricalLengthDist, UniformBase};
use genairr_engine::ir::{NucHandle, Nucleotide, Region, Segment, Simulation};
use genairr_engine::pass::PassPlan;
use genairr_engine::passes::{AssembleSegmentPass, GenerateNPPass, SampleAllelePass};
use genairr_engine::refdata::{Allele, ChainType, RefDataConfig};
use genairr_engine::s5f::{S5FKernel, S5F_NUM_CONTEXTS, S5F_SUBSTITUTION_LEN};

pub const SEED_RANGE: u64 = 100;

pub fn vj_refdata() -> RefDataConfig {
    let mut cfg = RefDataConfig::empty(ChainType::Vj);
    let _ = cfg.v_pool.push(Allele {
        name: "v_test*01".into(),
        gene: "v_test".into(),
        seq: b"AAACCCGGG".to_vec(),
        segment: Segment::V,
        anchor: Some(6),
    });
    let _ = cfg.j_pool.push(Allele {
        name: "j_test*01".into(),
        gene: "j_test".into(),
        seq: b"TTTAAA".to_vec(),
        segment: Segment::J,
        anchor: Some(0),
    });
    cfg
}

pub fn vj_plan(refdata: &RefDataConfig) -> PassPlan {
    let length_dist = || {
        Box::new(EmpiricalLengthDist::from_pairs(
            (0..10).map(|i| (i, 1.0)).collect::<Vec<_>>(),
        ))
    };
    let mut plan = PassPlan::new();
    plan.push(Box::new(SampleAllelePass::new(
        Segment::V,
        Box::new(AllelePoolDist::uniform(&refdata.v_pool)),
    )));
    plan.push(Box::new(SampleAllelePass::new(
        Segment::J,
        Box::new(AllelePoolDist::uniform(&refdata.j_pool)),
    )));
    plan.push(Box::new(AssembleSegmentPass::new(Segment::V)));
    plan.push(Box::new(GenerateNPPass::new(
        Segment::Np1,
        length_dist(),
        Box::new(UniformBase),
    )));
    plan.push(Box::new(AssembleSegmentPass::new(Segment::J)));
    plan
}

pub fn assembled_v_sim() -> Simulation {
    let mut sim = Simulation::new();
    for (i, b) in b"AAACCCGGGTTT".iter().enumerate() {
        let (next, _) = sim.with_nucleotide_pushed(Nucleotide::germline(*b, i as u16, Segment::V));
        sim = next;
    }
    let region = Region::new(Segment::V, NucHandle::new(0), NucHandle::new(12));
    sim.with_region_added(region)
}

pub fn uniform_s5f() -> S5FKernel {
    S5FKernel::new(
        vec![1.0; S5F_NUM_CONTEXTS],
        vec![0.25; S5F_SUBSTITUTION_LEN],
    )
}

pub fn assert_codon_rails_consistent(sim: &Simulation, label: &str) {
    // the per-region rail is no longer maintained in the
    // hot path. Region.amino_acids stays empty after assembly /
    // mutation / indel passes; consumers compute the rail on demand
    // via `with_codon_rail_recomputed`. So "stored == fresh" no
    // longer means anything (stored is always empty). The invariant
    // we DO want to verify is that a fresh recompute is deterministic
    // and produces a rail that's internally consistent with the
    // post-pipeline pool — i.e., every stop_codon_position handle
    // falls within the region range and the amino_acids count
    // matches the codon count derivable from `(region.len(), skip)`.
    for (i, region) in sim.sequence.regions.iter().enumerate() {
        let fresh = genairr_engine::ir::compute_codon_rail(region, &sim.pool);
        let skip = (3 - (region.frame_phase as u32)) % 3;
        let coding_bytes = region.len().saturating_sub(skip);
        let expected_codons = (coding_bytes / 3) as usize;
        assert_eq!(
            fresh.amino_acids.len(),
            expected_codons,
            "{}: region[{}] codon count mismatch — bytes={} skip={} expected={} got={}",
            label,
            i,
            region.len(),
            skip,
            expected_codons,
            fresh.amino_acids.len(),
        );
        for h in &fresh.stop_codon_positions {
            assert!(
                h.index() >= region.start.index() && h.index() < region.end.index(),
                "{}: region[{}] stop_codon at handle {} outside region [{}, {})",
                label,
                i,
                h.index(),
                region.start.index(),
                region.end.index(),
            );
        }
    }
}

pub fn assert_region_ranges_valid(sim: &Simulation, label: &str) {
    let pool_len = sim.pool.len() as u32;
    for (i, region) in sim.sequence.regions.iter().enumerate() {
        let s = region.start.index();
        let e = region.end.index();
        assert!(
            s <= e,
            "{}: region[{}] has inverted range start={} end={}",
            label,
            i,
            s,
            e
        );
        assert!(
            e <= pool_len,
            "{}: region[{}] end {} exceeds pool len {}",
            label,
            i,
            e,
            pool_len
        );
    }
}

pub fn assert_frame_phases_consistent(sim: &Simulation, label: &str) {
    let mut cumulative_len = 0u64;
    for (i, region) in sim.sequence.regions.iter().enumerate() {
        let expected = (cumulative_len % 3) as u8;
        assert_eq!(
            region.frame_phase, expected,
            "{}: region[{}] has stale frame_phase {} (expected {})",
            label, i, region.frame_phase, expected
        );
        cumulative_len = cumulative_len.saturating_add(region.len() as u64);
    }
}
