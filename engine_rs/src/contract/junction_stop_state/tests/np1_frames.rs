use super::*;

#[test]
fn admits_np_candidate_frame_phase_2_rejects_completing_stop() {
    let cfg = make_vj_refdata(b"GGG", b"GGG");
    let build_sim = make_sim_v_assembled(&cfg, false);
    let state = JunctionStopState::build(&build_sim, &cfg, Segment::Np1, 3).unwrap();

    let mut sim = build_sim;
    let (next, _) =
        sim.with_nucleotide_pushed(Nucleotide::synthetic(b'T', Segment::Np1, flag::N_NUC));
    sim = next;
    let (next, _) =
        sim.with_nucleotide_pushed(Nucleotide::synthetic(b'A', Segment::Np1, flag::N_NUC));
    sim = next;

    let err = state
        .admits_np_candidate(&sim, Segment::Np1, 2, b'A')
        .unwrap_err();
    assert!(err.reason.contains("TAA"), "{}", err.reason);

    assert!(state
        .admits_np_candidate(&sim, Segment::Np1, 2, b'C')
        .is_ok());
}

#[test]
fn admits_np_candidate_frame_phase_0_passes_incomplete_codon() {
    let cfg = make_vj_refdata(b"GGG", b"GGG");
    let sim = make_sim_v_assembled(&cfg, false);
    let state = JunctionStopState::build(&sim, &cfg, Segment::Np1, 3).unwrap();
    assert!(state
        .admits_np_candidate(&sim, Segment::Np1, 0, b'T')
        .is_ok());
    assert!(state
        .admits_np_candidate(&sim, Segment::Np1, 0, b'A')
        .is_ok());
}

#[test]
fn admits_np_candidate_frame_phase_1_passes_incomplete_codon() {
    let cfg = make_vj_refdata(b"GGG", b"GGG");
    let build_sim = make_sim_v_assembled(&cfg, false);
    let state = JunctionStopState::build(&build_sim, &cfg, Segment::Np1, 3).unwrap();

    let mut sim = build_sim;
    let (next, _) =
        sim.with_nucleotide_pushed(Nucleotide::synthetic(b'T', Segment::Np1, flag::N_NUC));
    sim = next;

    assert!(state
        .admits_np_candidate(&sim, Segment::Np1, 1, b'A')
        .is_ok());
}
