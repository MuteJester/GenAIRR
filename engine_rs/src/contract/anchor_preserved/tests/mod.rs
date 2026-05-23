use super::*;
use crate::assignment::AlleleInstance;
use crate::contract::test_support::{
    make_assembled_sim_from_refdata, make_v_anchor_at, make_vj_with_anchor_codons,
};
use crate::ir::{flag, NucHandle, Nucleotide};
use crate::refdata::{Allele, AlleleId, ChainType, RefDataConfig};

#[test]
#[should_panic(expected = "segment must be V, D, or J")]
fn anchor_preserved_rejects_np_segment() {
    let _ = AnchorPreserved::new(Segment::Np1);
}

#[test]
fn anchor_preserved_name_per_segment() {
    assert_eq!(
        AnchorPreserved::new(Segment::V).name(),
        "anchor_preserved.v"
    );
    assert_eq!(
        AnchorPreserved::new(Segment::D).name(),
        "anchor_preserved.d"
    );
    assert_eq!(
        AnchorPreserved::new(Segment::J).name(),
        "anchor_preserved.j"
    );
}

#[test]
fn anchor_preserved_no_allele_assigned_passes_vacuously() {
    let cfg = make_v_anchor_at(10, Some(5));
    let contract = AnchorPreserved::new(Segment::V);
    let sim = Simulation::new();

    // No assignment -> vacuously satisfied.
    assert!(contract.verify(&sim, Some(&cfg)).is_ok());
}

#[test]
fn anchor_preserved_no_refdata_passes_vacuously() {
    let cfg = make_v_anchor_at(10, Some(5));
    let contract = AnchorPreserved::new(Segment::V);
    let sim =
        Simulation::new().with_allele_assigned(Segment::V, AlleleInstance::new(AlleleId::new(0)));

    // Refdata None -> can't verify; treat as satisfied.
    assert!(contract.verify(&sim, None).is_ok());
    // Sanity: WITH refdata it would also pass (no trims).
    assert!(contract.verify(&sim, Some(&cfg)).is_ok());
}

#[test]
fn anchor_preserved_anchorless_allele_passes_vacuously() {
    let cfg = make_v_anchor_at(10, None);
    let contract = AnchorPreserved::new(Segment::V);
    let sim =
        Simulation::new().with_allele_assigned(Segment::V, AlleleInstance::new(AlleleId::new(0)));

    assert!(contract.verify(&sim, Some(&cfg)).is_ok());
}

#[test]
fn anchor_preserved_zero_trims_and_anchor_in_range_passes() {
    let cfg = make_v_anchor_at(10, Some(3));
    let contract = AnchorPreserved::new(Segment::V);
    let sim =
        Simulation::new().with_allele_assigned(Segment::V, AlleleInstance::new(AlleleId::new(0)));

    assert!(contract.verify(&sim, Some(&cfg)).is_ok());
}

#[test]
fn anchor_preserved_rejects_live_anchor_provenance_shift_after_indel() {
    let cfg = make_vj_with_anchor_codons(b"GGG", b"TTT");
    let sim = make_assembled_sim_from_refdata(&cfg);
    let shifted = sim.with_indel_inserted(
        6,
        Nucleotide::synthetic(b'A', Segment::V, flag::INDEL_INSERTED),
    );

    let err = AnchorPreserved::new(Segment::V)
        .verify(&shifted, Some(&cfg))
        .unwrap_err();

    assert_eq!(err.contract_name, "anchor_preserved.v");
    assert!(
        err.reason.contains("provenance mismatch"),
        "reason should mention live provenance mismatch, got: {}",
        err.reason
    );
}

#[test]
fn anchor_preserved_rejects_live_anchor_amino_acid_change() {
    let cfg = make_vj_with_anchor_codons(b"TGT", b"TGG");
    let sim = make_assembled_sim_from_refdata(&cfg);
    let changed = sim.with_base_changed(NucHandle::new(6), b'A');

    let err = AnchorPreserved::new(Segment::V)
        .verify(&changed, Some(&cfg))
        .unwrap_err();

    assert_eq!(err.contract_name, "anchor_preserved.v");
    assert!(
        err.reason.contains("anchor codon amino acid changed"),
        "reason should mention anchor amino-acid change, got: {}",
        err.reason
    );
}

#[test]
fn anchor_preserved_allows_synonymous_anchor_base_change() {
    let cfg = make_vj_with_anchor_codons(b"TGT", b"TGG");
    let sim = make_assembled_sim_from_refdata(&cfg);
    let synonymous = sim.with_base_changed(NucHandle::new(2), b'C');

    assert!(AnchorPreserved::new(Segment::V)
        .verify(&synonymous, Some(&cfg))
        .is_ok());
}

#[test]
fn anchor_preserved_rejects_targeted_anchor_substitution_candidate() {
    let cfg = make_vj_with_anchor_codons(b"TGT", b"TGG");
    let sim = make_assembled_sim_from_refdata(&cfg);
    let contract = AnchorPreserved::new(Segment::V);

    let err = contract
        .admits_with_context(
            &sim,
            Some(&cfg),
            "mutate.uniform.base[0]",
            &ChoiceValue::Base(b'A'),
            ChoiceContext::targeted_base_substitution(0, 1, NucHandle::new(6)),
        )
        .unwrap_err();

    assert_eq!(err.contract_name, "anchor_preserved.v");
    assert!(
        err.reason.contains("would change V anchor amino acid"),
        "reason should mention anchor amino-acid filtering, got: {}",
        err.reason
    );
}

#[test]
fn anchor_preserved_violates_when_5prime_trim_eats_anchor() {
    use crate::assignment::TrimEnd;

    let cfg = make_v_anchor_at(10, Some(3));
    let contract = AnchorPreserved::new(Segment::V);
    let sim = Simulation::new()
        .with_allele_assigned(Segment::V, AlleleInstance::new(AlleleId::new(0)))
        .with_trim(Segment::V, TrimEnd::Five, 4);

    let err = contract.verify(&sim, Some(&cfg)).unwrap_err();
    assert_eq!(err.contract_name, "anchor_preserved.v");
    assert!(
        err.reason.contains("trim_5"),
        "reason should mention trim_5, got: {}",
        err.reason
    );
}

#[test]
fn anchor_preserved_violates_when_3prime_trim_eats_anchor() {
    use crate::assignment::TrimEnd;

    let cfg = make_v_anchor_at(10, Some(6));
    let contract = AnchorPreserved::new(Segment::V);
    let sim = Simulation::new()
        .with_allele_assigned(Segment::V, AlleleInstance::new(AlleleId::new(0)))
        .with_trim(Segment::V, TrimEnd::Three, 2);

    let err = contract.verify(&sim, Some(&cfg)).unwrap_err();
    assert_eq!(err.contract_name, "anchor_preserved.v");
    assert!(
        err.reason.contains("retained slice ends at"),
        "reason should describe retained slice, got: {}",
        err.reason
    );
}

#[test]
fn anchor_preserved_anchor_at_5prime_boundary_passes() {
    let cfg = make_v_anchor_at(10, Some(0));
    let contract = AnchorPreserved::new(Segment::V);
    let sim =
        Simulation::new().with_allele_assigned(Segment::V, AlleleInstance::new(AlleleId::new(0)));

    assert!(contract.verify(&sim, Some(&cfg)).is_ok());
}

#[test]
fn anchor_preserved_anchor_at_3prime_boundary_passes() {
    let cfg = make_v_anchor_at(10, Some(7));
    let contract = AnchorPreserved::new(Segment::V);
    let sim =
        Simulation::new().with_allele_assigned(Segment::V, AlleleInstance::new(AlleleId::new(0)));

    assert!(contract.verify(&sim, Some(&cfg)).is_ok());
}

#[test]
fn anchor_preserved_one_more_trim_just_past_boundary_violates() {
    use crate::assignment::TrimEnd;

    let cfg = make_v_anchor_at(10, Some(7));
    let contract = AnchorPreserved::new(Segment::V);
    let sim = Simulation::new()
        .with_allele_assigned(Segment::V, AlleleInstance::new(AlleleId::new(0)))
        .with_trim(Segment::V, TrimEnd::Three, 1);

    assert!(contract.verify(&sim, Some(&cfg)).is_err());
}

#[test]
fn anchor_preserved_works_for_j_segment_too() {
    use crate::assignment::TrimEnd;

    let mut cfg = RefDataConfig::empty(ChainType::Vdj);
    let _ = cfg.j_pool.push(Allele {
        name: "j_test*01".to_string(),
        gene: "j_test".to_string(),
        seq: vec![b'A'; 12],
        segment: Segment::J,
        anchor: Some(0),
    });
    let contract = AnchorPreserved::new(Segment::J);

    let ok =
        Simulation::new().with_allele_assigned(Segment::J, AlleleInstance::new(AlleleId::new(0)));
    assert!(contract.verify(&ok, Some(&cfg)).is_ok());

    let bad = ok.with_trim(Segment::J, TrimEnd::Five, 1);
    assert!(contract.verify(&bad, Some(&cfg)).is_err());
}

#[test]
fn contract_works_through_box_dyn() {
    let cfg = make_v_anchor_at(10, Some(3));
    let contract: Box<dyn Contract> = Box::new(AnchorPreserved::new(Segment::V));
    let sim =
        Simulation::new().with_allele_assigned(Segment::V, AlleleInstance::new(AlleleId::new(0)));

    assert!(contract.verify(&sim, Some(&cfg)).is_ok());
    assert_eq!(contract.name(), "anchor_preserved.v");
}
