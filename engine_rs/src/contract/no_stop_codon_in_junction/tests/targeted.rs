use super::*;

fn assert_rejects_targeted_address(address: &str, target: NucHandle) {
    let cfg = make_vj_with_anchor_codons(b"TAC", b"GGG");
    let c = NoStopCodonInJunction::new();
    let sim = make_assembled_sim_from_refdata(&cfg);

    let context = ChoiceContext::targeted_base_substitution(0, 1, target)
        .with_address_if_missing(ChoiceAddress::parse(address));
    let err = c
        .admits_typed(&sim, Some(&cfg), context, &ChoiceValue::Base(b'A'))
        .unwrap_err();

    assert_eq!(err.contract_name, "no_stop_codon_in_junction");
    assert!(err.reason.contains("TAA"), "{}", err.reason);
}

#[test]
fn no_stop_codon_admits_rejects_uniform_mutation_base_that_creates_stop() {
    assert_rejects_targeted_address("mutate.uniform.base[0]", NucHandle::new(8));

    let cfg = make_vj_with_anchor_codons(b"TAC", b"GGG");
    let c = NoStopCodonInJunction::new();
    let sim = make_assembled_sim_from_refdata(&cfg);
    let context = ChoiceContext::targeted_base_substitution(0, 1, NucHandle::new(8))
        .with_address_if_missing(ChoiceAddress::parse("mutate.uniform.base[0]"));
    assert!(c
        .admits_typed(&sim, Some(&cfg), context, &ChoiceValue::Base(b'C'))
        .is_ok());
}

#[test]
fn no_stop_codon_admits_rejects_s5f_mutation_base_that_creates_stop() {
    assert_rejects_targeted_address("mutate.s5f.base[0]", NucHandle::new(8));
}

#[test]
fn no_stop_codon_admits_rejects_pcr_error_base_that_creates_stop() {
    assert_rejects_targeted_address("corrupt.pcr.error_base[0]", NucHandle::new(8));
}

#[test]
fn no_stop_codon_admits_rejects_contaminant_base_that_creates_stop() {
    let cfg = make_vj_with_anchor_codons(b"AAA", b"GGG");
    let c = NoStopCodonInJunction::new();
    let sim = make_assembled_sim_from_refdata(&cfg);

    let context = ChoiceContext::targeted_base_substitution(6, 12, NucHandle::new(6))
        .with_address_if_missing(ChoiceAddress::parse("corrupt.contaminant.bases[6]"));
    let err = c
        .admits_typed(&sim, Some(&cfg), context, &ChoiceValue::Base(b'T'))
        .unwrap_err();

    assert_eq!(err.contract_name, "no_stop_codon_in_junction");
    assert!(err.reason.contains("TAA"), "{}", err.reason);
}
