use super::super::junction::anchor_amino_acid_preserved;
use super::*;

#[test]
fn anchor_amino_acid_preserved_allows_synonymous_change() {
    let (cfg, sim) = anchor_record_fixture();
    let changed = sim.with_base_changed(NucHandle::new(2), b'C');
    let region = changed
        .sequence
        .regions
        .iter()
        .find(|r| r.segment == Segment::V)
        .unwrap();

    assert!(anchor_amino_acid_preserved(
        &changed,
        &cfg,
        Segment::V,
        region,
        Some(AlleleId::new(0)),
        0,
    ));
}

#[test]
fn anchor_amino_acid_preserved_rejects_nonsynonymous_change() {
    let (cfg, sim) = anchor_record_fixture();
    let changed = sim.with_base_changed(NucHandle::new(0), b'A');
    let region = changed
        .sequence
        .regions
        .iter()
        .find(|r| r.segment == Segment::V)
        .unwrap();

    assert!(!anchor_amino_acid_preserved(
        &changed,
        &cfg,
        Segment::V,
        region,
        Some(AlleleId::new(0)),
        0,
    ));
}
