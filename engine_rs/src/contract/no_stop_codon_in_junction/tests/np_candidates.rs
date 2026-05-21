use super::*;

#[test]
fn no_stop_codon_admits_rejects_np_base_that_completes_stop() {
    let (cfg, sim) = make_partial_np_stop_filter_case();
    let c = NoStopCodonInJunction::new();

    let err = c
        .admits(
            &sim,
            Some(&cfg),
            "np.np1.bases[0]",
            &ChoiceValue::Base(b'A'),
        )
        .unwrap_err();
    assert_eq!(err.contract_name, "no_stop_codon_in_junction");
    assert!(err.reason.contains("TAA"), "{}", err.reason);

    assert!(c
        .admits(
            &sim,
            Some(&cfg),
            "np.np1.bases[0]",
            &ChoiceValue::Base(b'C')
        )
        .is_ok());
}

#[test]
fn no_stop_codon_admits_ignores_np_base_until_codon_complete() {
    let (cfg, sim) = make_partial_np_stop_filter_case();
    let c = NoStopCodonInJunction::new();
    let (sim, _) = sim.with_nucleotide_pushed(Nucleotide::synthetic(
        b'T',
        Segment::Np1,
        crate::ir::flag::N_NUC,
    ));

    assert!(c
        .admits(
            &sim,
            Some(&cfg),
            "np.np1.bases[1]",
            &ChoiceValue::Base(b'A')
        )
        .is_ok());
}

#[test]
fn no_stop_codon_admits_rejects_np1_base_that_forces_future_d_stop() {
    let (cfg, sim) = make_partial_np1_future_d_stop_case(b"GGG");
    let c = NoStopCodonInJunction::new();

    let err = c
        .admits_with_context(
            &sim,
            Some(&cfg),
            "np.np1.bases[0]",
            &ChoiceValue::Base(b'T'),
            ChoiceContext::indexed(0, 1),
        )
        .unwrap_err();
    assert_eq!(err.contract_name, "no_stop_codon_in_junction");
    assert!(err.reason.contains("TAA"), "{}", err.reason);

    assert!(c
        .admits_with_context(
            &sim,
            Some(&cfg),
            "np.np1.bases[0]",
            &ChoiceValue::Base(b'C'),
            ChoiceContext::indexed(0, 1),
        )
        .is_ok());
}

#[test]
fn no_stop_codon_admits_rejects_zero_np1_length_that_forces_future_d_stop() {
    let (cfg, sim) = make_partial_np1_future_d_stop_case(b"GGGT");
    let c = NoStopCodonInJunction::new();

    let err = c
        .admits(&sim, Some(&cfg), "np.np1.length", &ChoiceValue::Int(0))
        .unwrap_err();
    assert_eq!(err.contract_name, "no_stop_codon_in_junction");
    assert!(err.reason.contains("TAA"), "{}", err.reason);

    assert!(c
        .admits(&sim, Some(&cfg), "np.np1.length", &ChoiceValue::Int(1))
        .is_ok());
}
