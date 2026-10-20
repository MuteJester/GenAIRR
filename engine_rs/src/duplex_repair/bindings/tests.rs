use super::*;
use crate::refdata::{Allele, AnchorEvidenceStatus, ChainType};

fn reference() -> RefDataConfig {
    let mut data = RefDataConfig::empty(ChainType::Vdj);
    for (segment, name, sequence, anchor) in [
        (Segment::V, "V*01", "AAATGT", Some(3)),
        (Segment::D, "D*01", "CAA", None),
        (Segment::J, "J*01", "TTTAAA", Some(0)),
    ] {
        let allele = Allele {
            name: name.into(),
            gene: name.split('*').next().unwrap().into(),
            seq: sequence.as_bytes().to_vec(),
            segment,
            anchor,
            functional_status: None,
            subregions: vec![],
        };
        let _ = match segment {
            Segment::V => data.v_pool.push(allele),
            Segment::D => data.d_pool.push(allele),
            Segment::J => data.j_pool.push(allele),
            _ => unreachable!(),
        };
    }
    data
}

fn bind(data: &RefDataConfig, asserted: bool) -> Result<DuplexReferenceBindings, String> {
    DuplexReferenceBindings::new(
        data,
        ["V*01", "D*01", "J*01"],
        ["AAATGT", "CAA", "TTTAAA"],
        false,
        asserted,
    )
}

fn set_boundary(
    data: &mut RefDataConfig,
    segment: Segment,
    five: CodingBoundaryStatus,
    three: CodingBoundaryStatus,
) {
    let allele = data.get(segment, AlleleId::new(0)).unwrap();
    let evidence = AlleleEvidence {
        schema_version: 1,
        sequence_sha256: sequence_sha256(&allele.seq),
        coding_phase: None,
        anchor_status: if segment == Segment::D {
            AnchorEvidenceStatus::Unknown
        } else {
            AnchorEvidenceStatus::Inferred
        },
        boundary_5: five,
        boundary_3: three,
        source: None,
        method: "authored_test_boundary".into(),
        notes: vec![],
    };
    data.set_allele_evidence(segment, AlleleId::new(0), evidence)
        .unwrap();
}

#[test]
fn unknown_ends_require_assertion_and_preserve_unknown_identity() {
    let data = reference();
    assert!(bind(&data, false)
        .unwrap_err()
        .contains("unknown coding end"));
    let bound = bind(&data, true).unwrap();
    assert_eq!(bound.v.boundary_contract, "caller_asserted");
    assert_eq!(bound.v.boundary_3, CodingBoundaryStatus::Unknown);
    assert!(bound.reference_identity["species"].is_null());
    assert_eq!(bound.v.anchor, Some(3));
    assert_eq!(bound.v.coding_phase, Some(0));
    assert!(bound.v.landmark_ready);
}

#[test]
fn known_partial_facing_end_cannot_be_overridden_but_nonfacing_coverage_is_retained() {
    let mut data = reference();
    set_boundary(
        &mut data,
        Segment::V,
        CodingBoundaryStatus::Partial,
        CodingBoundaryStatus::Complete,
    );
    assert_eq!(
        bind(&data, true).unwrap().v.boundary_5,
        CodingBoundaryStatus::Partial
    );
    for segment in [Segment::V, Segment::D, Segment::J] {
        let mut invalid = reference();
        let (five, three) = if segment == Segment::V {
            (CodingBoundaryStatus::Unknown, CodingBoundaryStatus::Partial)
        } else {
            (CodingBoundaryStatus::Partial, CodingBoundaryStatus::Unknown)
        };
        set_boundary(&mut invalid, segment, five, three);
        assert!(bind(&invalid, true).unwrap_err().contains("known partial"));
    }
}

#[test]
fn declared_complete_ends_do_not_require_caller_assertion() {
    let mut data = reference();
    for segment in [Segment::V, Segment::D, Segment::J] {
        set_boundary(
            &mut data,
            segment,
            CodingBoundaryStatus::Complete,
            CodingBoundaryStatus::Complete,
        );
    }
    let bound = bind(&data, false).unwrap();
    assert_eq!(bound.d.boundary_contract, "declared_complete");
    assert!(!bound.reference_ends_are_coding_ends);
}

#[test]
fn bindings_preserve_distinct_alleles_with_identical_sequences() {
    let mut data = reference();
    let mut duplicate_sequence = data.v_pool.get(AlleleId::new(0)).unwrap().clone();
    duplicate_sequence.name = "V*02".into();
    let _ = data.v_pool.push(duplicate_sequence);
    let bound = DuplexReferenceBindings::new(
        &data,
        ["V*02", "D*01", "J*01"],
        ["AAATGT", "CAA", "TTTAAA"],
        false,
        true,
    )
    .unwrap();
    assert_eq!(bound.allele_id(Substrate::V), AlleleId::new(1));
    assert_eq!(bound.v.allele_name, "V*02");
    assert!(DuplexReferenceBindings::new(
        &data,
        ["invented", "D*01", "J*01"],
        ["AAATGT", "CAA", "TTTAAA"],
        false,
        true,
    )
    .unwrap_err()
    .contains("unknown generating"));
    assert!(DuplexReferenceBindings::new(
        &data,
        ["V*01", "D*01", "J*01"],
        ["TGT", "CAA", "TTTAAA"],
        false,
        true,
    )
    .unwrap_err()
    .contains("entire canonical"));
}

#[test]
fn validation_rejects_changed_evidence_binding_or_product_orientation() {
    use crate::duplex_repair::{repair_coupled, JoinName, JoinScript, RepairOperation, Side};
    let data = reference();
    let bound = bind(&data, true).unwrap();
    let scripts: Vec<_> = [(JoinName::VD, 6), (JoinName::DJ, 3)]
        .into_iter()
        .map(|(join, right_start)| JoinScript {
            join,
            operations: vec![
                RepairOperation::Open {
                    side: Side::Left,
                    offset: 0,
                },
                RepairOperation::Open {
                    side: Side::Right,
                    offset: 0,
                },
                RepairOperation::Synapse { right_start },
                RepairOperation::Ligate {},
            ],
        })
        .collect();
    let mut product =
        repair_coupled("AAATGT", "CAA", "TTTAAA", false, &scripts, 256, 100_000).unwrap();
    bound.validate(&product, &data).unwrap();
    let mut edited = bound.clone();
    edited.v.coding_phase = Some(1);
    assert!(edited.validate(&product, &data).is_err());
    let mut changed = data.clone();
    changed.identity.source = Some("changed reference declaration".into());
    assert!(bound.validate(&product, &changed).is_err());
    product.d_inverted = true;
    assert!(bound
        .validate(&product, &data)
        .unwrap_err()
        .contains("orientation"));
}
