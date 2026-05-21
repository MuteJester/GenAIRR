use super::super::{
    assembled_segment_live_call, with_assembled_segment_live_call, BoundaryValue, EvidenceScore,
    HypothesisFlags, LiveCallConfidence, ReferenceMatchIndex,
};
use super::{allele, simulation_with_region};
use crate::ir::{NucFlags, NucHandle, Nucleotide, Region, Segment, Simulation};
use crate::refdata::{ChainType, RefDataConfig};

#[test]
fn assembled_segment_live_call_keeps_trim_induced_ambiguity() {
    let mut cfg = RefDataConfig::empty(ChainType::Vdj);
    let id0 = cfg.v_pool.push(allele(Segment::V, "v1*01", b"GGAAACCC"));
    let id1 = cfg.v_pool.push(allele(Segment::V, "v2*01", b"TTAAACCC"));
    let _id2 = cfg.v_pool.push(allele(Segment::V, "v3*01", b"GGTATCCC"));
    let index = ReferenceMatchIndex::build(&cfg);
    let sim = simulation_with_region(Segment::V, b"AAACCC", 2);

    let call = assembled_segment_live_call(&sim, &index, Segment::V, 1)
        .expect("V region should produce a call");

    assert_eq!(call.allele_call.to_ids(), vec![id0, id1]);
    assert_eq!(call.confidence, LiveCallConfidence::ExactAmbiguous);
    assert_eq!(call.boundary_summary.seq_start, BoundaryValue::Single(0));
    assert_eq!(call.boundary_summary.seq_end, BoundaryValue::Single(6));
    assert_eq!(call.boundary_summary.ref_start, BoundaryValue::Single(2));
    assert_eq!(call.boundary_summary.ref_end, BoundaryValue::Single(8));
    assert_eq!(call.hypotheses[0].score, EvidenceScore::exact(6, 0));
}

#[test]
fn assembled_segment_live_call_uses_current_observed_bases() {
    let mut cfg = RefDataConfig::empty(ChainType::Vdj);
    let _id0 = cfg.v_pool.push(allele(Segment::V, "v1*01", b"AACT"));
    let id1 = cfg.v_pool.push(allele(Segment::V, "v1*02", b"AAGT"));
    let index = ReferenceMatchIndex::build(&cfg);
    let sim = simulation_with_region(Segment::V, b"AACT", 0);
    let sim = sim.with_base_changed(NucHandle::new(2), b'G');

    let call = assembled_segment_live_call(&sim, &index, Segment::V, 1)
        .expect("V region should produce a call");

    assert_eq!(call.allele_call.to_ids(), vec![id1]);
    assert_eq!(call.confidence, LiveCallConfidence::ExactSingle);
    assert_eq!(call.hypotheses[0].score, EvidenceScore::exact(4, 0));
}

#[test]
fn y6_truth_allele_retained_under_single_position_mutation() {
    let mut cfg = RefDataConfig::empty(ChainType::Vdj);
    let v0 = cfg.v_pool.push(allele(Segment::V, "v*01", b"ACGTAC"));
    let v1 = cfg.v_pool.push(allele(Segment::V, "v*02", b"ACGTAG"));
    let v2 = cfg.v_pool.push(allele(Segment::V, "v*03", b"AGGTAC"));
    let index = ReferenceMatchIndex::build(&cfg);

    let sim = simulation_with_region(Segment::V, b"ACGTAC", 0);
    let sim = sim.with_base_changed(NucHandle::new(2), b'T');

    let call = assembled_segment_live_call(&sim, &index, Segment::V, 1)
        .expect("V region should produce a call");

    assert_eq!(call.allele_call.to_ids(), vec![v0]);
    let _unused = (v1, v2);
}

#[test]
fn y6_v_call_diverges_from_truth_when_mutations_flip_toward_another_allele() {
    let mut cfg = RefDataConfig::empty(ChainType::Vdj);
    let v0 = cfg.v_pool.push(allele(Segment::V, "v*01", b"ACGTAC"));
    let v1 = cfg.v_pool.push(allele(Segment::V, "v*02", b"AGGTGA"));
    let index = ReferenceMatchIndex::build(&cfg);

    let sim = simulation_with_region(Segment::V, b"ACGTAC", 0);
    let sim = sim.with_base_changed(NucHandle::new(1), b'G');
    let sim = sim.with_base_changed(NucHandle::new(4), b'G');
    let sim = sim.with_base_changed(NucHandle::new(5), b'A');

    let call = assembled_segment_live_call(&sim, &index, Segment::V, 1)
        .expect("V region should produce a call");

    assert_eq!(call.allele_call.to_ids(), vec![v1]);
    let _unused = v0;
}

#[test]
fn y6_trim_ambiguity_narrows_via_np_extension() {
    let mut cfg = RefDataConfig::empty(ChainType::Vdj);
    let v0 = cfg.v_pool.push(allele(Segment::V, "v*01", b"AACCAA"));
    let v1 = cfg.v_pool.push(allele(Segment::V, "v*02", b"AACCTG"));
    let v2 = cfg.v_pool.push(allele(Segment::V, "v*03", b"AACCGG"));
    let index = ReferenceMatchIndex::build(&cfg);

    let mut sim = Simulation::new();
    for (i, &b) in b"AACC".iter().enumerate() {
        let (next, _) = sim.with_nucleotide_pushed(Nucleotide::germline(b, i as u16, Segment::V));
        sim = next;
    }
    sim = sim.with_region_added(Region::new(
        Segment::V,
        NucHandle::new(0),
        NucHandle::new(4),
    ));

    for &b in &[b'T', b'G'] {
        let (next, _) =
            sim.with_nucleotide_pushed(Nucleotide::synthetic(b, Segment::Np1, NucFlags::empty()));
        sim = next;
    }
    sim = sim.with_region_added(Region::new(
        Segment::Np1,
        NucHandle::new(4),
        NucHandle::new(6),
    ));

    let call = assembled_segment_live_call(&sim, &index, Segment::V, 1)
        .expect("V region should produce a call");

    assert_eq!(call.allele_call.to_ids(), vec![v1]);
    let _unused = (v0, v2);
}

#[test]
fn y6_full_pool_returned_when_no_allele_matches_any_position() {
    let mut cfg = RefDataConfig::empty(ChainType::Vdj);
    let v0 = cfg.v_pool.push(allele(Segment::V, "v*01", b"AACT"));
    let v1 = cfg.v_pool.push(allele(Segment::V, "v*02", b"GGCA"));
    let index = ReferenceMatchIndex::build(&cfg);

    let sim = simulation_with_region(Segment::V, b"AACT", 0);
    let sim = sim.with_base_changed(NucHandle::new(0), b'T');
    let sim = sim.with_base_changed(NucHandle::new(1), b'T');
    let sim = sim.with_base_changed(NucHandle::new(2), b'A');
    let sim = sim.with_base_changed(NucHandle::new(3), b'G');

    let call = assembled_segment_live_call(&sim, &index, Segment::V, 1)
        .expect("V region should produce a call");

    let mut ids = call.allele_call.to_ids();
    ids.sort_by_key(|i| i.index());
    assert_eq!(ids, vec![v0, v1]);
}

#[test]
fn assembled_segment_live_call_treats_n_as_non_informative() {
    let mut cfg = RefDataConfig::empty(ChainType::Vdj);
    let id0 = cfg.v_pool.push(allele(Segment::V, "v1*01", b"AACT"));
    let id1 = cfg.v_pool.push(allele(Segment::V, "v1*02", b"AAGT"));
    let index = ReferenceMatchIndex::build(&cfg);
    let sim = simulation_with_region(Segment::V, b"AACT", 0);
    let sim = sim.with_base_changed(NucHandle::new(2), b'N');

    let call = assembled_segment_live_call(&sim, &index, Segment::V, 1)
        .expect("V region should produce a call");

    assert_eq!(call.allele_call.to_ids(), vec![id0, id1]);
    assert_eq!(call.confidence, LiveCallConfidence::ExactAmbiguous);
    assert_eq!(call.hypotheses[0].score, EvidenceScore::exact(3, 1));
    assert!(call.hypotheses[0]
        .flags
        .contains(HypothesisFlags::HAS_WILDCARD_EVIDENCE));
}

#[test]
fn with_assembled_segment_live_call_persists_state_update() {
    let mut cfg = RefDataConfig::empty(ChainType::Vdj);
    let id0 = cfg.v_pool.push(allele(Segment::V, "v1*01", b"AACT"));
    let index = ReferenceMatchIndex::build(&cfg);
    let sim = simulation_with_region(Segment::V, b"AACT", 0);

    let next = with_assembled_segment_live_call(&sim, &index, Segment::V);

    assert!(sim.live_calls.is_none());
    let live = next
        .live_calls
        .expect("live call state should be initialized");
    let v = live.get(Segment::V).expect("V call should be present");
    assert_eq!(live.version, 1);
    assert_eq!(v.evidence_version, 1);
    assert_eq!(v.allele_call.to_ids(), vec![id0]);
}
