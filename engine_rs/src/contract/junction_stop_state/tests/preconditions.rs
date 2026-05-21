use super::*;

#[test]
fn build_returns_none_without_v_assignment() {
    let cfg = make_vj_refdata(b"GGG", b"TTT");
    let sim = Simulation::new();
    let state = JunctionStopState::build(&sim, &cfg, Segment::Np1, 3);
    assert!(state.is_none());
}

#[test]
fn build_returns_none_without_j_assignment() {
    let cfg = make_vj_refdata(b"GGG", b"TTT");
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

    let state = JunctionStopState::build(&sim, &cfg, Segment::Np1, 3);
    assert!(state.is_none());
}

#[test]
fn build_returns_none_with_anchorless_v() {
    let mut cfg = make_vj_refdata(b"GGG", b"TTT");
    let v = cfg.v_pool.get(AlleleId::new(0)).unwrap().clone();
    cfg.v_pool = AllelePool::from_vec(vec![Allele { anchor: None, ..v }]);

    let sim = make_sim_v_assembled(&cfg, false);
    let state = JunctionStopState::build(&sim, &cfg, Segment::Np1, 3);
    assert!(state.is_none());
}

#[test]
fn build_returns_some_on_valid_sim() {
    let cfg = make_vj_refdata(b"GGG", b"TTT");
    let sim = make_sim_v_assembled(&cfg, false);
    let state = JunctionStopState::build(&sim, &cfg, Segment::Np1, 0);
    assert!(state.is_some());
}
