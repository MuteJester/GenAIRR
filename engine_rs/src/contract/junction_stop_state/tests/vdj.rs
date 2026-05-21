use super::*;

#[test]
fn admits_np_candidate_vdj_spans_np1_to_d_rejects_stop() {
    let cfg = make_vdj_refdata(b"GGG", b"AAC", b"GGG");
    let sim = make_sim_v_assembled(&cfg, true);
    let state = JunctionStopState::build(&sim, &cfg, Segment::Np1, 1).unwrap();

    let err = state
        .admits_np_candidate(&sim, Segment::Np1, 0, b'T')
        .unwrap_err();
    assert!(err.reason.contains("TAA"), "{}", err.reason);

    assert!(state
        .admits_np_candidate(&sim, Segment::Np1, 0, b'C')
        .is_ok());
}

#[test]
fn admits_np_candidate_vdj_np2_rejects_completing_stop() {
    let cfg = make_vdj_refdata(b"GGG", b"TA", b"GGG");
    let sim = make_sim_v_assembled(&cfg, true);

    let np1_start = NucHandle::new(sim.pool.len() as u32);
    let sim = sim.with_region_added(Region::new(Segment::Np1, np1_start, np1_start));

    let state = JunctionStopState::build(&sim, &cfg, Segment::Np2, 1).unwrap();
    let err = state
        .admits_np_candidate(&sim, Segment::Np2, 0, b'A')
        .unwrap_err();
    assert!(err.reason.contains("TAA"), "{}", err.reason);

    assert!(state
        .admits_np_candidate(&sim, Segment::Np2, 0, b'C')
        .is_ok());
}
