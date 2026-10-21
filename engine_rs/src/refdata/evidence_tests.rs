use super::evidence::{resolve_anchor, sequence_sha256};
use super::*;
use crate::trace_file::refdata_content_hash;

fn fixture() -> RefDataConfig {
    let mut r = RefDataConfig::empty(ChainType::Vj);
    let _ = r.v_pool.push(Allele {
        name: "V1*01".into(),
        gene: "V1".into(),
        seq: b"AAATGTAAA".to_vec(),
        segment: Segment::V,
        anchor: Some(3),
        functional_status: None,
        subregions: vec![],
    });
    let _ = r.j_pool.push(Allele {
        name: "J1*01".into(),
        gene: "J1".into(),
        seq: b"TTTGGGAAAGGG".to_vec(),
        segment: Segment::J,
        anchor: Some(0),
        functional_status: None,
        subregions: vec![],
    });
    r
}
fn evidence(seq: &[u8]) -> AlleleEvidence {
    AlleleEvidence {
        schema_version: 1,
        sequence_sha256: sequence_sha256(seq),
        coding_phase: Some(0),
        anchor_status: AnchorEvidenceStatus::Verified,
        boundary_5: CodingBoundaryStatus::Unknown,
        boundary_3: CodingBoundaryStatus::Complete,
        source: Some("synthetic_fixture".into()),
        method: "explicit_test_annotation".into(),
        notes: vec![],
    }
}
#[test]
fn reference_evidence_roundtrip_hash_and_stale_rejection() {
    let mut r = fixture();
    let baseline = refdata_content_hash(&r);
    let id = AlleleId::new(0);
    let e = evidence(&r.v_pool.get(id).unwrap().seq);
    r.set_allele_evidence(Segment::V, id, e.clone()).unwrap();
    assert_ne!(baseline, refdata_content_hash(&r));
    assert_eq!(r.allele_evidence(Segment::V, id), Some(&e));
    assert!(r.validate().is_empty());
    let mut stale = e.clone();
    stale.sequence_sha256 = "sha256:stale".into();
    assert!(r.set_allele_evidence(Segment::V, id, stale).is_err());
    let mut bad = e;
    bad.coding_phase = Some(1);
    assert!(r.set_allele_evidence(Segment::V, id, bad).is_err());
    r.evidence.clear();
    assert_eq!(baseline, refdata_content_hash(&r));
}
#[test]
fn reference_evidence_hash_is_order_independent() {
    let mut a = fixture();
    let mut b = fixture();
    for s in [Segment::V, Segment::J] {
        let e = evidence(&a.get(s, AlleleId::new(0)).unwrap().seq);
        a.set_allele_evidence(s, AlleleId::new(0), e).unwrap();
    }
    for s in [Segment::J, Segment::V] {
        let e = evidence(&b.get(s, AlleleId::new(0)).unwrap().seq);
        b.set_allele_evidence(s, AlleleId::new(0), e).unwrap();
    }
    assert_eq!(refdata_content_hash(&a), refdata_content_hash(&b));
}
#[test]
fn native_anchor_resolver_does_not_guess_v_terminal_cys() {
    let r = resolve_anchor(b"AAATGTAAA", Segment::V, "IGH", None, false, None);
    assert_eq!(r.position, None);
    assert_eq!(r.status, AnchorEvidenceStatus::Unknown);
}
#[test]
fn native_anchor_resolver_maps_declared_imgt_gaps_and_rejects_wrong_codon() {
    let mut g = vec![b'A'; 315];
    g[3..6].copy_from_slice(b"...");
    g[309..312].copy_from_slice(b"TGT");
    let seq: Vec<_> = g.iter().copied().filter(|b| *b != b'.').collect();
    let r = resolve_anchor(&seq, Segment::V, "IGH", Some(&g), true, None);
    assert_eq!(r.position, Some(306));
    assert_eq!(r.status, AnchorEvidenceStatus::Inferred);
    g[309..312].copy_from_slice(b"AAA");
    let seq: Vec<_> = g.iter().copied().filter(|b| *b != b'.').collect();
    assert_eq!(
        resolve_anchor(&seq, Segment::V, "IGH", Some(&g), true, None).status,
        AnchorEvidenceStatus::Invalid
    );
}
#[test]
fn native_j_resolver_uses_locus_and_does_not_choose_ambiguous_motifs() {
    assert_eq!(
        resolve_anchor(b"TTTGGGAAAGGG", Segment::J, "TRB", None, false, None).position,
        Some(0)
    );
    assert_eq!(
        resolve_anchor(b"TTTGGGAAAGGG", Segment::J, "IGH", None, false, None).position,
        None
    );
    let r = resolve_anchor(
        b"TTTGGGAAAGGGTTTGGGAAAGGG",
        Segment::J,
        "TRB",
        None,
        false,
        None,
    );
    assert_eq!(r.position, None);
    assert_eq!(r.candidates, vec![0, 12]);
}
#[test]
fn explicit_unresolved_evidence_is_reported_even_if_a_legacy_codon_matches() {
    let mut r = fixture();
    let mut e = evidence(&r.v_pool.get(AlleleId::new(0)).unwrap().seq);
    e.anchor_status = AnchorEvidenceStatus::Unknown;
    e.coding_phase = None;
    r.set_allele_evidence(Segment::V, AlleleId::new(0), e)
        .unwrap();
    assert!(r
        .validate()
        .iter()
        .any(|i| matches!(i, RefDataValidationIssue::ReferenceEvidenceProblem { .. })));
    assert!(r
        .validate_with_mode(RefDataValidationMode::AllowCuratable)
        .is_ok());
}

#[test]
fn known_partial_coding_end_is_checked_independently_of_landmark_and_phase() {
    let mut r = fixture();
    let id = AlleleId::new(0);
    let mut e = evidence(&r.v_pool.get(id).unwrap().seq);
    e.boundary_5 = CodingBoundaryStatus::Partial;
    r.set_allele_evidence(Segment::V, id, e.clone()).unwrap();
    assert!(r.coding_ends_are_not_known_partial(Segment::V, id));
    e.boundary_3 = CodingBoundaryStatus::Partial;
    r.set_allele_evidence(Segment::V, id, e).unwrap();
    assert!(!r.coding_ends_are_not_known_partial(Segment::V, id));
    assert!(r.validate().is_empty()); // useful partial references remain loadable
}

#[test]
fn direct_coding_end_passes_reject_partial_tips_before_trace_or_pool_commit() {
    use crate::assignment::AlleleInstance;
    use crate::end_processing::EndOutcome;
    use crate::ir::Simulation;
    use crate::passes::test_support::run_pass_capturing_events;
    use crate::passes::{EmitCodingEndPPass, ProcessCodingEndsPass};
    let mut r = fixture();
    let id = AlleleId::new(0);
    let mut e = evidence(&r.v_pool.get(id).unwrap().seq);
    e.boundary_3 = CodingBoundaryStatus::Partial;
    r.set_allele_evidence(Segment::V, id, e).unwrap();
    let sim = Simulation::new().with_allele_assigned(Segment::V, AlleleInstance::new(id));
    let process = ProcessCodingEndsPass::new(
        Segment::V,
        vec![EndOutcome {
            five: 0,
            three: -1,
            weight: 1.0,
        }],
    )
    .unwrap();
    assert!(run_pass_capturing_events(&process, &sim, Some(&r), None).is_err());
    let forged = sim.with_coding_end_choice(Segment::V, id, 0, -1, 0);
    let emit = EmitCodingEndPPass::new(crate::address::PEnd::V3);
    assert!(run_pass_capturing_events(&emit, &forged, Some(&r), None).is_err());
}
