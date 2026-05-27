use super::*;
use crate::address::{ChoiceAddress, NpSegment};
use crate::contract::ChoiceContext;

fn np_base_addr(index: u32) -> ChoiceAddress {
    ChoiceAddress::NpBase {
        segment: NpSegment::Np1,
        index,
    }
}

const NP1_LENGTH_ADDR: ChoiceAddress = ChoiceAddress::NpLength(NpSegment::Np1);

#[test]
fn no_stop_codon_admits_rejects_np_base_that_completes_stop() {
    let (cfg, sim) = make_partial_np_stop_filter_case();
    let c = NoStopCodonInJunction::new();

    let context = ChoiceContext::none().with_address(np_base_addr(0));
    let err = c
        .admits_typed(&sim, Some(&cfg), context, &ChoiceValue::Base(b'A'))
        .unwrap_err();
    assert_eq!(err.contract_name, "no_stop_codon_in_junction");
    assert!(err.reason.contains("TAA"), "{}", err.reason);

    assert!(c
        .admits_typed(&sim, Some(&cfg), context, &ChoiceValue::Base(b'C'))
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

    let context = ChoiceContext::none().with_address(np_base_addr(1));
    assert!(c
        .admits_typed(&sim, Some(&cfg), context, &ChoiceValue::Base(b'A'))
        .is_ok());
}

#[test]
fn no_stop_codon_admits_rejects_np1_base_that_forces_future_d_stop() {
    let (cfg, sim) = make_partial_np1_future_d_stop_case(b"GGG");
    let c = NoStopCodonInJunction::new();

    let context = ChoiceContext::indexed(0, 1).with_address(np_base_addr(0));
    let err = c
        .admits_typed(&sim, Some(&cfg), context, &ChoiceValue::Base(b'T'))
        .unwrap_err();
    assert_eq!(err.contract_name, "no_stop_codon_in_junction");
    assert!(err.reason.contains("TAA"), "{}", err.reason);

    assert!(c
        .admits_typed(&sim, Some(&cfg), context, &ChoiceValue::Base(b'C'))
        .is_ok());
}

#[test]
fn no_stop_codon_admits_rejects_zero_np1_length_that_forces_future_d_stop() {
    let (cfg, sim) = make_partial_np1_future_d_stop_case(b"GGGT");
    let c = NoStopCodonInJunction::new();

    let context = ChoiceContext::none().with_address(NP1_LENGTH_ADDR);
    let err = c
        .admits_typed(&sim, Some(&cfg), context, &ChoiceValue::Int(0))
        .unwrap_err();
    assert_eq!(err.contract_name, "no_stop_codon_in_junction");
    assert!(err.reason.contains("TAA"), "{}", err.reason);

    assert!(c
        .admits_typed(&sim, Some(&cfg), context, &ChoiceValue::Int(1))
        .is_ok());
}
