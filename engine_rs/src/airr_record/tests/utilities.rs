use super::super::junction::{aa_slice_for_region, junction_has_stop};
use super::super::locus::{derive_locus, locus_from_refdata};
use super::super::walk;
use super::*;
use crate::codon::{translate_codon_slice, translate_seq};

#[test]
fn translate_codon_basic() {
    assert_eq!(translate_codon_slice(b"ATG"), 'M');
    assert_eq!(translate_codon_slice(b"taa"), '*');
    assert_eq!(translate_codon_slice(b"NNN"), 'X');
}

#[test]
fn translate_seq_truncates_partial_codon() {
    assert_eq!(translate_seq("ATGCCC"), "MP");
    assert_eq!(translate_seq("ATGCCCA"), "MP");
}

#[test]
fn derive_locus_picks_first_known_prefix() {
    assert_eq!(derive_locus("IGHV1-2*01", "", ""), "IGH");
    assert_eq!(derive_locus("", "IGKJ4*01", ""), "IGK");
    assert_eq!(derive_locus("ighv1-2*01", "", ""), "IGH");
    assert_eq!(derive_locus("", "", ""), "");
    assert_eq!(derive_locus("XYZ1*01", "", ""), "");
}

#[test]
fn locus_from_refdata_falls_back_to_pool_allele_names() {
    let mut cfg = RefDataConfig::empty(ChainType::Vdj);
    let _ = cfg.v_pool.push(Allele {
        name: "IGHV1-2*01".into(),
        gene: "IGHV1-2".into(),
        seq: b"AAA".to_vec(),
        segment: Segment::V,
        anchor: None,
        functional_status: None,
        subregions: Vec::new(),
    });
    assert_eq!(locus_from_refdata(&cfg), "IGH");

    let mut alien = RefDataConfig::empty(ChainType::Vdj);
    let _ = alien.v_pool.push(Allele {
        name: "XYZV1*01".into(),
        gene: "XYZV1".into(),
        seq: b"AAA".to_vec(),
        segment: Segment::V,
        anchor: None,
        functional_status: None,
        subregions: Vec::new(),
    });
    assert_eq!(locus_from_refdata(&alien), "");

    let empty = RefDataConfig::empty(ChainType::Vj);
    assert_eq!(locus_from_refdata(&empty), "");
}

#[test]
fn runlength_collapses_repeated_ops() {
    let runs = vec![(5, b'M'), (2, b'I'), (3, b'M'), (1, b'D')];
    assert_eq!(walk::helpers::runlength_to_string(&runs), "5M2I3M1D");
}

#[test]
fn runlength_empty_is_empty_string() {
    assert_eq!(walk::helpers::runlength_to_string(&[]), "");
}

#[test]
fn aa_slice_for_region_basic() {
    assert_eq!(aa_slice_for_region(3, 12, 0, "ABCDE"), "BCD");
    assert_eq!(aa_slice_for_region(4, 12, 0, "ABCDE"), "CD");
    assert_eq!(aa_slice_for_region(5, 5, 0, "ABCDE"), "");
    assert_eq!(aa_slice_for_region(4, 5, 0, "ABCDE"), "");
}

#[test]
fn junction_has_stop_detects_stops() {
    assert!(junction_has_stop("TAA"));
    assert!(junction_has_stop("ATGTAA"));
    assert!(junction_has_stop("ATGTAG"));
    assert!(junction_has_stop("ATGTGA"));
    assert!(!junction_has_stop("ATG"));
    assert!(!junction_has_stop("ATGCCC"));
}
