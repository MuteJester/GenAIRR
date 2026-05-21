use super::*;

#[test]
fn no_stop_codon_name_is_canonical() {
    assert_eq!(
        NoStopCodonInJunction::new().name(),
        "no_stop_codon_in_junction"
    );
}

#[test]
fn no_stop_codon_no_refdata_passes_vacuously() {
    let c = NoStopCodonInJunction::new();
    assert!(c.verify(&Simulation::new(), None).is_ok());
}

#[test]
fn no_stop_codon_undefined_junction_passes_vacuously() {
    let cfg = make_vj_with_anchor_codons(b"GGG", b"TTT");
    let c = NoStopCodonInJunction::new();
    assert!(c.verify(&Simulation::new(), Some(&cfg)).is_ok());
}

#[test]
fn no_stop_codon_in_frame_no_stops_passes() {
    let cfg = make_vj_with_anchor_codons(b"GGG", b"TTT");
    let c = NoStopCodonInJunction::new();
    let sim = make_assembled_sim_from_refdata(&cfg);

    assert!(c.verify(&sim, Some(&cfg)).is_ok());
}

#[test]
fn no_stop_codon_v_anchor_taa_violates() {
    let cfg = make_vj_with_anchor_codons(b"TAA", b"GGG");
    let c = NoStopCodonInJunction::new();
    let sim = make_assembled_sim_from_refdata(&cfg);

    let err = c.verify(&sim, Some(&cfg)).unwrap_err();
    assert_eq!(err.contract_name, "no_stop_codon_in_junction");
    assert!(err.reason.contains("stop"), "{}", err.reason);
    assert!(err.reason.contains("TAA"), "{}", err.reason);
}

#[test]
fn no_stop_codon_v_anchor_tag_violates() {
    let cfg = make_vj_with_anchor_codons(b"TAG", b"GGG");
    let c = NoStopCodonInJunction::new();
    let sim = make_assembled_sim_from_refdata(&cfg);

    let err = c.verify(&sim, Some(&cfg)).unwrap_err();
    assert!(err.reason.contains("TAG"), "{}", err.reason);
}

#[test]
fn no_stop_codon_v_anchor_tga_violates() {
    let cfg = make_vj_with_anchor_codons(b"TGA", b"GGG");
    let c = NoStopCodonInJunction::new();
    let sim = make_assembled_sim_from_refdata(&cfg);

    let err = c.verify(&sim, Some(&cfg)).unwrap_err();
    assert!(err.reason.contains("TGA"), "{}", err.reason);
}

#[test]
fn no_stop_codon_j_anchor_stop_violates() {
    let cfg = make_vj_with_anchor_codons(b"GGG", b"TAA");
    let c = NoStopCodonInJunction::new();
    let sim = make_assembled_sim_from_refdata(&cfg);

    let err = c.verify(&sim, Some(&cfg)).unwrap_err();
    assert!(err.reason.contains("TAA"), "{}", err.reason);
}

#[test]
fn no_stop_codon_first_stop_reported() {
    let cfg = make_vj_with_anchor_codons(b"TAA", b"TAG");
    let c = NoStopCodonInJunction::new();
    let sim = make_assembled_sim_from_refdata(&cfg);

    let err = c.verify(&sim, Some(&cfg)).unwrap_err();
    assert!(err.reason.contains("TAA"), "{}", err.reason);
    assert!(err.reason.contains("position 6"), "{}", err.reason);
}

#[test]
fn no_stop_codon_out_of_frame_junction_passes_vacuously() {
    let mut cfg = make_vj_with_anchor_codons(b"GGG", b"TTT");
    let j_allele = cfg.j_pool.get(AlleleId::new(0)).unwrap().clone();
    cfg.j_pool = AllelePool::from_vec(vec![Allele {
        anchor: Some(2),
        ..j_allele
    }]);

    let c = NoStopCodonInJunction::new();
    let sim = make_assembled_sim_from_refdata(&cfg);
    assert!(c.verify(&sim, Some(&cfg)).is_ok());
}

#[test]
fn no_stop_codon_works_through_box_dyn() {
    let cfg = make_vj_with_anchor_codons(b"GGG", b"TTT");
    let sim = make_assembled_sim_from_refdata(&cfg);
    let c: Box<dyn Contract> = Box::new(NoStopCodonInJunction::new());
    assert!(c.verify(&sim, Some(&cfg)).is_ok());
    assert_eq!(c.name(), "no_stop_codon_in_junction");
}
