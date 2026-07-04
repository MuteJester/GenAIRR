    use super::*;
    use std::mem::size_of;

    /// Helper: a tiny synthetic V allele for tests.
    fn make_v(name: &str, gene: &str, seq: &[u8], anchor: Option<u16>) -> Allele {
        Allele {
            name: name.to_string(),
            gene: gene.to_string(),
            seq: seq.to_vec(),
            segment: Segment::V,
            anchor,
            functional_status: None,
            subregions: Vec::new(),
        }
    }

    #[test]
    fn allele_id_is_zero_cost_newtype() {
        assert_eq!(size_of::<AlleleId>(), size_of::<u32>());
    }

    #[test]
    fn allele_id_round_trip() {
        let id = AlleleId::new(7);
        assert_eq!(id.index(), 7);
        assert_eq!(id.as_usize(), 7);
    }

    #[test]
    fn chain_type_has_d_distinction() {
        assert!(!ChainType::Vj.has_d());
        assert!(ChainType::Vdj.has_d());
    }

    #[test]
    fn allele_basic_accessors() {
        let a = make_v("IGHV1-2*01", "IGHV1-2", b"ACGTACGT", Some(3));
        assert_eq!(a.len(), 8);
        assert!(a.has_anchor());
        assert_eq!(a.anchor, Some(3));
        assert_eq!(a.segment, Segment::V);
    }

    #[test]
    fn allele_anchorless_round_trip() {
        let a = make_v("IGHV-pseudo*01", "IGHV-pseudo", b"ACGT", None);
        assert!(!a.has_anchor());
        assert_eq!(a.anchor, None);
    }

    #[test]
    fn allele_pool_starts_empty() {
        let p = AllelePool::new();
        assert_eq!(p.len(), 0);
        assert!(p.is_empty());
        assert!(p.get(AlleleId::new(0)).is_none());
        assert!(p.find_by_name("nonexistent").is_none());
    }

    #[test]
    fn allele_pool_push_returns_sequential_ids() {
        let mut p = AllelePool::new();
        let id0 = p.push(make_v("a*01", "a", b"AA", Some(0)));
        let id1 = p.push(make_v("b*01", "b", b"CC", Some(0)));
        let id2 = p.push(make_v("c*01", "c", b"GG", None));

        assert_eq!(id0.index(), 0);
        assert_eq!(id1.index(), 1);
        assert_eq!(id2.index(), 2);
        assert_eq!(p.len(), 3);
    }

    #[test]
    fn allele_pool_get_returns_stored_allele() {
        let mut p = AllelePool::new();
        let id = p.push(make_v("x*01", "x", b"ATGC", Some(2)));
        let got = p.get(id).expect("just pushed");
        assert_eq!(got.name, "x*01");
        assert_eq!(got.seq, b"ATGC");
        assert_eq!(got.anchor, Some(2));
    }

    #[test]
    fn allele_pool_get_out_of_bounds_returns_none() {
        let p = AllelePool::new();
        assert!(p.get(AlleleId::new(99)).is_none());
    }

    #[test]
    fn allele_pool_iter_yields_id_allele_pairs_in_order() {
        let mut p = AllelePool::new();
        let _ = p.push(make_v("a*01", "a", b"AA", None));
        let _ = p.push(make_v("b*01", "b", b"CC", None));
        let _ = p.push(make_v("c*01", "c", b"GG", None));

        let collected: Vec<(u32, String)> = p
            .iter()
            .map(|(id, a)| (id.index(), a.name.clone()))
            .collect();
        assert_eq!(
            collected,
            vec![
                (0, "a*01".to_string()),
                (1, "b*01".to_string()),
                (2, "c*01".to_string()),
            ]
        );
    }

    #[test]
    fn allele_pool_find_by_name_locates_exact_match() {
        let mut p = AllelePool::new();
        let _ = p.push(make_v("IGHV1-2*01", "IGHV1-2", b"AA", None));
        let target_id = p.push(make_v("IGHV1-2*02", "IGHV1-2", b"AC", None));
        let _ = p.push(make_v("IGHV3-23*01", "IGHV3-23", b"GG", None));

        let (id, allele) = p.find_by_name("IGHV1-2*02").expect("name should exist");
        assert_eq!(id, target_id);
        assert_eq!(allele.seq, b"AC");

        // Partial / similar names should not match.
        assert!(p.find_by_name("IGHV1-2").is_none());
        assert!(p.find_by_name("IGHV1-2*03").is_none());
    }

    #[test]
    fn ref_data_config_empty_for_chain_type() {
        let cfg = RefDataConfig::empty(ChainType::Vdj);
        assert_eq!(cfg.chain_type, ChainType::Vdj);
        assert!(cfg.v_pool.is_empty());
        assert!(cfg.d_pool.is_empty());
        assert!(cfg.j_pool.is_empty());
        assert!(cfg.c_pool.is_empty());
    }

    #[test]
    fn ref_data_config_pool_for_segment_routes_correctly() {
        let mut cfg = RefDataConfig::empty(ChainType::Vdj);
        let _ = cfg.v_pool.push(make_v("v*01", "v", b"AA", None));
        let _ = cfg.d_pool.push(Allele {
            name: "d*01".into(),
            gene: "d".into(),
            seq: b"GG".to_vec(),
            segment: Segment::D,
            anchor: None,
            functional_status: None,
            subregions: Vec::new(),
        });
        let _ = cfg.j_pool.push(Allele {
            name: "j*01".into(),
            gene: "j".into(),
            seq: b"TT".to_vec(),
            segment: Segment::J,
            anchor: Some(0),
            functional_status: None,
            subregions: Vec::new(),
        });

        assert_eq!(cfg.pool_for(Segment::V).unwrap().len(), 1);
        assert_eq!(cfg.pool_for(Segment::D).unwrap().len(), 1);
        assert_eq!(cfg.pool_for(Segment::J).unwrap().len(), 1);
        assert!(cfg.pool_for(Segment::Np1).is_none());
        assert!(cfg.pool_for(Segment::Np2).is_none());
    }

    #[test]
    fn ref_data_config_get_resolves_segment_id_pair() {
        let mut cfg = RefDataConfig::empty(ChainType::Vdj);
        let v_id = cfg.v_pool.push(make_v("v*01", "v", b"AAAT", Some(1)));

        let v = cfg.get(Segment::V, v_id).expect("v*01 should resolve");
        assert_eq!(v.name, "v*01");
        assert_eq!(v.anchor, Some(1));

        // Wrong segment -> None.
        assert!(cfg.get(Segment::J, v_id).is_none());
        // NP segment -> None defensively.
        assert!(cfg.get(Segment::Np1, v_id).is_none());
    }

    // ── Curation policy ───────────────────────────────────────────

    fn allele_with_anchor(name: &str, seq: &[u8], anchor: Option<u16>, segment: Segment) -> Allele {
        Allele {
            name: name.to_string(),
            gene: name.split('*').next().unwrap_or(name).to_string(),
            seq: seq.to_vec(),
            segment,
            anchor,
            functional_status: None,
            subregions: Vec::new(),
        }
    }

    fn vj_cfg_for_curation() -> RefDataConfig {
        // Two V alleles: one Cys-anchored, one Gly-anchored.
        // Two J alleles: one W-anchored, one anchorless.
        let mut cfg = RefDataConfig::empty(ChainType::Vj);
        let _ = cfg
            .v_pool
            .push(allele_with_anchor("MYV-good*01", b"TGTAAACCC", Some(0), Segment::V));
        let _ = cfg
            .v_pool
            .push(allele_with_anchor("MYV-gly*01", b"GGGAAACCC", Some(0), Segment::V));
        let _ = cfg
            .j_pool
            .push(allele_with_anchor("MYJ-good*01", b"TGGAAACCC", Some(0), Segment::J));
        let _ = cfg
            .j_pool
            .push(allele_with_anchor("MYJ-orphan*01", b"GGG", None, Segment::J));
        cfg
    }

    #[test]
    fn curation_raw_is_identity() {
        let cfg = vj_cfg_for_curation();
        let curated = cfg.curated(RefDataCurationPolicy::Raw);
        assert_eq!(curated.v_pool.len(), 2);
        assert_eq!(curated.j_pool.len(), 2);
        // Raw policy does NOT touch identity.
        assert_eq!(curated.identity.source, None);
    }

    #[test]
    fn curation_functional_anchors_only_filters_v_and_j() {
        let cfg = vj_cfg_for_curation();
        let curated = cfg.curated(RefDataCurationPolicy::FunctionalAnchorsOnly);
        assert_eq!(curated.v_pool.len(), 1);
        assert_eq!(curated.v_pool.iter().next().unwrap().1.name, "MYV-good*01");
        assert_eq!(curated.j_pool.len(), 1);
        assert_eq!(curated.j_pool.iter().next().unwrap().1.name, "MYJ-good*01");
    }

    #[test]
    fn curation_tags_identity_source() {
        let cfg = vj_cfg_for_curation();
        let curated = cfg.curated(RefDataCurationPolicy::FunctionalAnchorsOnly);
        assert_eq!(
            curated.identity.source.as_deref(),
            Some("curated:functional_anchors_only"),
        );
    }

    #[test]
    fn curation_extends_existing_identity_source() {
        let mut cfg = vj_cfg_for_curation();
        cfg.identity.source = Some("DataConfig".to_string());
        let curated = cfg.curated(RefDataCurationPolicy::FunctionalAnchorsOnly);
        assert_eq!(
            curated.identity.source.as_deref(),
            Some("DataConfig|curated:functional_anchors_only"),
        );
    }

    #[test]
    fn curation_preserves_other_identity_fields() {
        let mut cfg = vj_cfg_for_curation();
        cfg.identity.species = Some("Human".into());
        cfg.identity.locus = Some("IGK".into());
        let curated = cfg.curated(RefDataCurationPolicy::FunctionalAnchorsOnly);
        assert_eq!(curated.identity.species.as_deref(), Some("Human"));
        assert_eq!(curated.identity.locus.as_deref(), Some("IGK"));
    }

    #[test]
    fn curation_does_not_silently_fix_duplicate_names() {
        // Add a duplicate V allele AFTER the good one. Both are
        // Cys-anchored — they survive curation, then validation
        // surfaces the duplicate as a Fatal issue. Curation never
        // removes a structurally-corrupt entry to make a catalogue
        // look healthier.
        let mut cfg = vj_cfg_for_curation();
        let _ = cfg.v_pool.push(allele_with_anchor(
            "MYV-good*01",
            b"TGTAAACCC",
            Some(0),
            Segment::V,
        ));
        let curated = cfg.curated(RefDataCurationPolicy::FunctionalAnchorsOnly);
        let issues = curated.validate();
        assert!(issues.iter().any(|i| matches!(
            i,
            RefDataValidationIssue::DuplicateAlleleName { .. }
        )));
    }

    #[test]
    fn curation_does_not_silently_fix_invalid_bytes() {
        // Add an allele with a gap byte ('.') and a valid Cys anchor.
        // The invalid byte is Fatal and must still surface post-curation.
        let mut cfg = vj_cfg_for_curation();
        let _ = cfg.v_pool.push(allele_with_anchor(
            "MYV-badbyte*01",
            b"TGT.AAACC",
            Some(0),
            Segment::V,
        ));
        let curated = cfg.curated(RefDataCurationPolicy::FunctionalAnchorsOnly);
        assert!(curated.validate().iter().any(|i| matches!(
            i,
            RefDataValidationIssue::InvalidAlleleByte { .. }
        )));
    }

    #[test]
    fn curation_to_empty_pool_surfaces_empty_required_pool() {
        // Every V allele is non-Cys — curation empties V.
        let mut cfg = RefDataConfig::empty(ChainType::Vj);
        let _ = cfg
            .v_pool
            .push(allele_with_anchor("MYV-gly*01", b"GGGAAACCC", Some(0), Segment::V));
        let _ = cfg
            .j_pool
            .push(allele_with_anchor("MYJ-good*01", b"TGGAAACCC", Some(0), Segment::J));
        let curated = cfg.curated(RefDataCurationPolicy::FunctionalAnchorsOnly);
        assert_eq!(curated.v_pool.len(), 0);
        let issues = curated.validate();
        assert!(issues.iter().any(|i| matches!(
            i,
            RefDataValidationIssue::EmptyRequiredPool { segment: Segment::V }
        )));
    }

    #[test]
    fn curation_d_pool_passes_through_unchanged() {
        let mut cfg = RefDataConfig::empty(ChainType::Vdj);
        let _ = cfg
            .v_pool
            .push(allele_with_anchor("MYV*01", b"TGTAAACCC", Some(0), Segment::V));
        let _ = cfg.d_pool.push(allele_with_anchor(
            "MYD*01",
            b"GGGCCCAAA",
            None,
            Segment::D,
        ));
        let _ = cfg
            .j_pool
            .push(allele_with_anchor("MYJ*01", b"TGGAAACCC", Some(0), Segment::J));
        let curated = cfg.curated(RefDataCurationPolicy::FunctionalAnchorsOnly);
        // Anchorless D survives — D pool is never filtered.
        assert_eq!(curated.d_pool.len(), 1);
    }

    #[test]
    fn curation_respects_custom_j_anchor_rule() {
        // Rule expects Y at J anchor. TGG (W) and TTC (F) are dropped;
        // TAT (Y) is kept.
        let mut cfg = RefDataConfig::empty(ChainType::Vj);
        cfg.rules.j_anchor.expected_amino_acids = vec!['Y'];
        let _ = cfg
            .v_pool
            .push(allele_with_anchor("MYV*01", b"TGTAAACCC", Some(0), Segment::V));
        let _ = cfg
            .j_pool
            .push(allele_with_anchor("MYJ-w*01", b"TGGAAACCC", Some(0), Segment::J));
        let _ = cfg
            .j_pool
            .push(allele_with_anchor("MYJ-y*01", b"TATAAACCC", Some(0), Segment::J));
        let curated = cfg.curated(RefDataCurationPolicy::FunctionalAnchorsOnly);
        assert_eq!(curated.j_pool.len(), 1);
        assert_eq!(curated.j_pool.iter().next().unwrap().1.name, "MYJ-y*01");
    }

    #[test]
    fn curation_anchor_required_false_keeps_anchorless_alleles() {
        // If the rule says anchor is optional, anchorless alleles
        // are kept under FunctionalAnchorsOnly.
        let mut cfg = RefDataConfig::empty(ChainType::Vj);
        cfg.rules.j_anchor.required = false;
        let _ = cfg
            .v_pool
            .push(allele_with_anchor("MYV*01", b"TGTAAACCC", Some(0), Segment::V));
        let _ = cfg
            .j_pool
            .push(allele_with_anchor("MYJ-orphan*01", b"GGG", None, Segment::J));
        let curated = cfg.curated(RefDataCurationPolicy::FunctionalAnchorsOnly);
        assert_eq!(curated.j_pool.len(), 1);
    }

    // ── Functional-status curation ────────────────────────────────

    fn allele_with_status(
        name: &str,
        seq: &[u8],
        segment: Segment,
        anchor: Option<u16>,
        status: Option<FunctionalStatus>,
    ) -> Allele {
        Allele {
            name: name.to_string(),
            gene: name.split('*').next().unwrap_or(name).to_string(),
            seq: seq.to_vec(),
            segment,
            anchor,
            functional_status: status,
            subregions: Vec::new(),
        }
    }

    fn mixed_status_cfg() -> RefDataConfig {
        let mut cfg = RefDataConfig::empty(ChainType::Vdj);
        let _ = cfg.v_pool.push(allele_with_status(
            "v-f*01", b"TGTAAACCC", Segment::V, Some(0),
            Some(FunctionalStatus::Functional),
        ));
        let _ = cfg.v_pool.push(allele_with_status(
            "v-o*01", b"TGTAAACCC", Segment::V, Some(0),
            Some(FunctionalStatus::Orf),
        ));
        let _ = cfg.v_pool.push(allele_with_status(
            "v-p*01", b"TGTAAACCC", Segment::V, Some(0),
            Some(FunctionalStatus::Pseudogene),
        ));
        let _ = cfg.v_pool.push(allele_with_status(
            "v-na*01", b"TGTAAACCC", Segment::V, Some(0), None,
        ));
        let _ = cfg.d_pool.push(allele_with_status(
            "d-f*01", b"GGG", Segment::D, None,
            Some(FunctionalStatus::Functional),
        ));
        let _ = cfg.d_pool.push(allele_with_status(
            "d-p*01", b"GGG", Segment::D, None,
            Some(FunctionalStatus::Pseudogene),
        ));
        let _ = cfg.j_pool.push(allele_with_status(
            "j-f*01", b"TGGAAACCC", Segment::J, Some(0),
            Some(FunctionalStatus::Functional),
        ));
        cfg
    }

    #[test]
    fn curation_functional_status_filters_v_d_j() {
        let cfg = mixed_status_cfg();
        let curated = cfg.curated(RefDataCurationPolicy::FunctionalStatus {
            allowed: vec![FunctionalStatus::Functional],
            keep_unannotated: false,
        });
        let v_names: Vec<&str> =
            curated.v_pool.iter().map(|(_, a)| a.name.as_str()).collect();
        let d_names: Vec<&str> =
            curated.d_pool.iter().map(|(_, a)| a.name.as_str()).collect();
        let j_names: Vec<&str> =
            curated.j_pool.iter().map(|(_, a)| a.name.as_str()).collect();
        assert_eq!(v_names, vec!["v-f*01"]);
        assert_eq!(d_names, vec!["d-f*01"]);
        assert_eq!(j_names, vec!["j-f*01"]);
    }

    #[test]
    fn curation_functional_status_keeps_unannotated_when_flagged() {
        let cfg = mixed_status_cfg();
        let curated = cfg.curated(RefDataCurationPolicy::FunctionalStatus {
            allowed: vec![FunctionalStatus::Functional],
            keep_unannotated: true,
        });
        let v_names: Vec<&str> =
            curated.v_pool.iter().map(|(_, a)| a.name.as_str()).collect();
        assert_eq!(v_names, vec!["v-f*01", "v-na*01"]);
    }

    #[test]
    fn curation_functional_status_tag_is_canonical() {
        let p = RefDataCurationPolicy::FunctionalStatus {
            allowed: vec![FunctionalStatus::Orf, FunctionalStatus::Functional],
            keep_unannotated: false,
        };
        // allowed list is sorted into canonical lexical order; the
        // policy tag is the source of identity-source provenance, so
        // two policies producing identical curated catalogues must
        // produce identical tags regardless of input order.
        assert_eq!(
            p.tag(),
            "functional_status:functional,orf|keep_unannotated=false",
        );
    }

    #[test]
    fn curation_functional_status_dedupes_allowed_in_tag() {
        let p = RefDataCurationPolicy::FunctionalStatus {
            allowed: vec![
                FunctionalStatus::Functional,
                FunctionalStatus::Functional,
            ],
            keep_unannotated: true,
        };
        assert_eq!(p.tag(), "functional_status:functional|keep_unannotated=true");
    }

    #[test]
    fn curation_functional_status_empty_v_surfaces_empty_required_pool() {
        let mut cfg = RefDataConfig::empty(ChainType::Vj);
        let _ = cfg.v_pool.push(allele_with_status(
            "v-p*01", b"TGTAAACCC", Segment::V, Some(0),
            Some(FunctionalStatus::Pseudogene),
        ));
        let _ = cfg.j_pool.push(allele_with_status(
            "j-f*01", b"TGGAAACCC", Segment::J, Some(0),
            Some(FunctionalStatus::Functional),
        ));
        let curated = cfg.curated(RefDataCurationPolicy::FunctionalStatus {
            allowed: vec![FunctionalStatus::Functional],
            keep_unannotated: false,
        });
        assert!(curated.v_pool.is_empty());
        let issues = curated.validate();
        assert!(issues.iter().any(|i| matches!(
            i,
            RefDataValidationIssue::EmptyRequiredPool { segment: Segment::V }
        )));
    }

    #[test]
    fn curation_functional_status_tags_identity_source() {
        let cfg = mixed_status_cfg();
        let curated = cfg.curated(RefDataCurationPolicy::FunctionalStatus {
            allowed: vec![FunctionalStatus::Functional],
            keep_unannotated: true,
        });
        let src = curated.identity.source.as_deref().unwrap_or("");
        assert!(src.contains("curated:functional_status:functional|keep_unannotated=true"));
    }

    #[test]
    fn ref_data_config_supports_vj_chain_with_empty_d_pool() {
        // VJ chains: d_pool is conventionally empty. Construction
        // does not enforce this; assembly (C.8) handles VJ vs VDJ
        // explicitly via chain_type.has_d().
        let mut cfg = RefDataConfig::empty(ChainType::Vj);
        let _ = cfg.v_pool.push(make_v("v*01", "v", b"AA", Some(0)));
        let _ = cfg.j_pool.push(Allele {
            name: "j*01".into(),
            gene: "j".into(),
            seq: b"TT".to_vec(),
            segment: Segment::J,
            anchor: Some(0),
            functional_status: None,
            subregions: Vec::new(),
        });

        assert!(!cfg.chain_type.has_d());
        assert!(cfg.d_pool.is_empty());
        assert_eq!(cfg.v_pool.len(), 1);
        assert_eq!(cfg.j_pool.len(), 1);
    }
