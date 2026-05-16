    use super::*;
    use crate::assignment::AlleleInstance;
    use crate::live_call::{
        AlleleBitSet, EvidenceScore, HypothesisFlags, LiveCallState, PlacementHypothesis,
        SegmentLiveCall,
    };
    use crate::refdata::{Allele, AlleleId, ChainType};
    use crate::trace::Trace;

    #[test]
    fn translate_codon_basic() {
        assert_eq!(translate_codon_slice(b"ATG"), 'M');
        assert_eq!(translate_codon_slice(b"taa"), '*');
        assert_eq!(translate_codon_slice(b"NNN"), 'X');
    }

    #[test]
    fn translate_seq_truncates_partial_codon() {
        assert_eq!(translate_seq("ATGCCC"), "MP");
        assert_eq!(translate_seq("ATGCCCA"), "MP"); // partial codon dropped
    }

    #[test]
    fn derive_locus_picks_first_known_prefix() {
        assert_eq!(derive_locus("IGHV1-2*01", "", ""), "IGH");
        assert_eq!(derive_locus("", "IGKJ4*01", ""), "IGK");
        assert_eq!(derive_locus("ighv1-2*01", "", ""), "IGH"); // case-insensitive
        assert_eq!(derive_locus("", "", ""), "");
        assert_eq!(derive_locus("XYZ1*01", "", ""), "");
    }

    #[test]
    fn locus_from_refdata_falls_back_to_pool_allele_names() {
        // when every live call has been wiped, the
        // refdata's pool allele names still tell us the locus.
        let mut cfg = RefDataConfig::empty(ChainType::Vdj);
        let _ = cfg.v_pool.push(Allele {
            name: "IGHV1-2*01".into(),
            gene: "IGHV1-2".into(),
            seq: b"AAA".to_vec(),
            segment: Segment::V,
            anchor: None,
        });
        assert_eq!(locus_from_refdata(&cfg), "IGH");

        // Non-AIRR-prefixed names should yield empty (the helper
        // requires a recognisable AIRR locus).
        let mut alien = RefDataConfig::empty(ChainType::Vdj);
        let _ = alien.v_pool.push(Allele {
            name: "XYZV1*01".into(),
            gene: "XYZV1".into(),
            seq: b"AAA".to_vec(),
            segment: Segment::V,
            anchor: None,
        });
        assert_eq!(locus_from_refdata(&alien), "");

        // Empty pool → empty locus.
        let empty = RefDataConfig::empty(ChainType::Vj);
        assert_eq!(locus_from_refdata(&empty), "");
    }

    #[test]
    fn runlength_collapses_repeated_ops() {
        let runs = vec![(5, b'M'), (2, b'I'), (3, b'M'), (1, b'D')];
        assert_eq!(runlength_to_string(&runs), "5M2I3M1D");
    }

    #[test]
    fn runlength_empty_is_empty_string() {
        assert_eq!(runlength_to_string(&[]), "");
    }

    #[test]
    fn aa_slice_for_region_basic() {
        // sequence_aa "ABCDE", frame_offset 0, region [3, 12) → "BCD".
        assert_eq!(aa_slice_for_region(3, 12, 0, "ABCDE"), "BCD");
        // Region [4, 12) starts mid-codon → first complete codon at 6.
        assert_eq!(aa_slice_for_region(4, 12, 0, "ABCDE"), "CD");
        // Empty region.
        assert_eq!(aa_slice_for_region(5, 5, 0, "ABCDE"), "");
        // Region too short for a complete codon.
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

    fn anchor_record_fixture() -> (RefDataConfig, Simulation) {
        let mut cfg = RefDataConfig::empty(ChainType::Vj);
        let _ = cfg.v_pool.push(Allele {
            name: "v_airr*01".into(),
            gene: "v_airr".into(),
            seq: b"TGT".to_vec(),
            segment: Segment::V,
            anchor: Some(0),
        });

        let mut sim = Simulation::new();
        for (i, &b) in b"TGT".iter().enumerate() {
            let (next, _) =
                sim.with_nucleotide_pushed(Nucleotide::germline(b, i as u16, Segment::V));
            sim = next;
        }
        let region = Region::new(Segment::V, NucHandle::new(0), NucHandle::new(3))
            .with_codon_rail_recomputed(&sim.pool);
        sim = sim.with_region_added(region).with_allele_assigned(
            Segment::V,
            crate::assignment::AlleleInstance::new(AlleleId::new(0)),
        );

        (cfg, sim)
    }

    fn outcome_from_sim(sim: Simulation) -> Outcome {
        Outcome {
            revisions: vec![sim],
            pass_names: Vec::new(),
            trace: Trace::new(),
            events: Vec::new(),
        }
    }

    fn call_projection_fixture() -> (RefDataConfig, Simulation, AlleleId, AlleleId) {
        let mut cfg = RefDataConfig::empty(ChainType::Vj);
        let v0 = cfg.v_pool.push(Allele {
            name: "IGHV1-1*01".into(),
            gene: "IGHV1-1".into(),
            seq: b"AAACCC".to_vec(),
            segment: Segment::V,
            anchor: None,
        });
        let v1 = cfg.v_pool.push(Allele {
            name: "IGHV1-1*02".into(),
            gene: "IGHV1-1".into(),
            seq: b"AAACCC".to_vec(),
            segment: Segment::V,
            anchor: None,
        });

        let mut sim = Simulation::new();
        for (i, &base) in b"AAACCC".iter().enumerate() {
            let (next, _) =
                sim.with_nucleotide_pushed(Nucleotide::germline(base, i as u16, Segment::V));
            sim = next;
        }
        let region = Region::new(Segment::V, NucHandle::new(0), NucHandle::new(6))
            .with_codon_rail_recomputed(&sim.pool);
        sim = sim
            .with_region_added(region)
            .with_allele_assigned(Segment::V, AlleleInstance::new(v0));

        (cfg, sim, v0, v1)
    }

    #[test]
    fn airr_call_projection_falls_back_to_origin_assignment_without_live_call() {
        let (cfg, sim, _v0, _v1) = call_projection_fixture();
        let outcome = outcome_from_sim(sim);
        let rec = build_airr_record(&outcome, &cfg, "fallback");

        assert_eq!(rec.v_call, "IGHV1-1*01");
        assert_eq!(rec.locus, "IGH");
    }

    #[test]
    fn airr_call_projection_prefers_live_call_allele_set() {
        let (cfg, sim, v0, v1) = call_projection_fixture();
        let hypothesis = PlacementHypothesis::new(
            Segment::V,
            0,
            6,
            0,
            6,
            AlleleBitSet::from_ids(cfg.v_pool.len(), [v0, v1]),
            EvidenceScore::exact(6, 0),
            HypothesisFlags::EMPTY,
        );
        let call =
            SegmentLiveCall::from_hypotheses(Segment::V, cfg.v_pool.len(), vec![hypothesis], 1);
        let live = LiveCallState::empty().with_segment_call(call);
        let outcome = outcome_from_sim(sim.with_live_calls(live));

        let rec = build_airr_record(&outcome, &cfg, "live");

        assert_eq!(rec.v_call, "IGHV1-1*01,IGHV1-1*02");
        assert_eq!(rec.locus, "IGH");
    }

    #[test]
    fn airr_call_projection_falls_back_to_truth_when_live_call_is_unsupported() {
        // An "unsupported" live call carries zero hypotheses and so
        // produces an empty `allele_call`. Under the score-and-tie
        // caller this state cannot arise from a real run (max=0
        // returns the full pool); it can only be constructed manually
        // as in this test. When it does arise, `projected_call_name`
        // falls through to the originally-sampled allele name, keeping
        // `v_call` consistent with `projected_allele_id`'s existing
        // fallback and ensuring the field is never silently empty when
        // a truth assignment exists.
        let (cfg, sim, _v0, _v1) = call_projection_fixture();
        let call = SegmentLiveCall::from_hypotheses(Segment::V, cfg.v_pool.len(), Vec::new(), 1);
        let live = LiveCallState::empty().with_segment_call(call);
        let outcome = outcome_from_sim(sim.with_live_calls(live));

        let rec = build_airr_record(&outcome, &cfg, "unsupported");

        assert_eq!(rec.v_call, "IGHV1-1*01");
    }

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
