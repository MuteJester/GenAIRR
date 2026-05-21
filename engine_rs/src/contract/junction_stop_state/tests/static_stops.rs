use super::*;

#[test]
fn build_detects_static_stop_in_v_anchor_codon() {
    let cfg = make_vj_refdata(b"TAA", b"GGG");
    let sim = make_sim_v_assembled(&cfg, false);
    let state = JunctionStopState::build(&sim, &cfg, Segment::Np1, 3).unwrap();
    let err = state
        .admits_np_candidate(&sim, Segment::Np1, 0, b'A')
        .unwrap_err();
    assert!(err.reason.contains("TAA"), "{}", err.reason);
}

#[test]
fn build_detects_static_stop_in_j_head() {
    let cfg = make_vj_refdata(b"GGG", b"TAA");
    let sim = make_sim_v_assembled(&cfg, false);
    let state = JunctionStopState::build(&sim, &cfg, Segment::Np1, 0).unwrap();
    let err = state
        .admits_np_candidate(&sim, Segment::Np1, 0, b'C')
        .unwrap_err();
    assert!(err.reason.contains("TAA"), "{}", err.reason);
}
