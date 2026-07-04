    use super::*;
    use crate::refdata::{Allele, AlleleId, ChainType, RefDataConfig};

    fn v_allele(name: &str, seq: &[u8], anchor: Option<u16>) -> Allele {
        Allele {
            name: name.to_string(),
            gene: name.split('*').next().unwrap_or(name).to_string(),
            seq: seq.to_vec(),
            segment: Segment::V,
            anchor,
            functional_status: None,
            subregions: Vec::new(),
        }
    }
    fn d_allele(name: &str, seq: &[u8]) -> Allele {
        Allele {
            name: name.to_string(),
            gene: name.split('*').next().unwrap_or(name).to_string(),
            seq: seq.to_vec(),
            segment: Segment::D,
            anchor: None,
            functional_status: None,
            subregions: Vec::new(),
        }
    }
    fn j_allele(name: &str, seq: &[u8], anchor: Option<u16>) -> Allele {
        Allele {
            name: name.to_string(),
            gene: name.split('*').next().unwrap_or(name).to_string(),
            seq: seq.to_vec(),
            segment: Segment::J,
            anchor,
            functional_status: None,
            subregions: Vec::new(),
        }
    }

    /// Build a minimal valid VDJ refdata so each test only varies
    /// the one field under examination.
    fn minimal_vdj() -> RefDataConfig {
        let mut cfg = RefDataConfig::empty(ChainType::Vdj);
        // V allele: anchor codon TGT (Cys) at position 0.
        let _ = cfg.v_pool.push(v_allele("IGHV1-1*01", b"TGTAAACCC", Some(0)));
        // D allele: no anchor.
        let _ = cfg.d_pool.push(d_allele("IGHD1-1*01", b"GGGCCCAAA"));
        // J allele: anchor codon TGG (Trp = W) at position 0.
        let _ = cfg.j_pool.push(j_allele("IGHJ1*01", b"TGGAAACCC", Some(0)));
        cfg
    }

    fn minimal_vj() -> RefDataConfig {
        let mut cfg = RefDataConfig::empty(ChainType::Vj);
        let _ = cfg.v_pool.push(v_allele("IGKV1-1*01", b"TGTAAACCC", Some(0)));
        // IGK J anchor codon TTC (Phe = F).
        let _ = cfg.j_pool.push(j_allele("IGKJ1*01", b"TTCAAACCC", Some(0)));
        cfg
    }

    #[test]
    fn minimal_vdj_validates_clean() {
        assert!(minimal_vdj().validate().is_empty());
        assert!(minimal_vdj().validate_strict().is_ok());
    }

    #[test]
    fn minimal_vj_validates_clean() {
        assert!(minimal_vj().validate().is_empty());
    }

    #[test]
    fn empty_v_pool_fails_vj() {
        let mut cfg = minimal_vj();
        cfg.v_pool = super::AllelePool::new();
        let issues = cfg.validate();
        assert!(issues.contains(&RefDataValidationIssue::EmptyRequiredPool {
            segment: Segment::V,
        }));
    }

    #[test]
    fn empty_v_pool_fails_vdj() {
        let mut cfg = minimal_vdj();
        cfg.v_pool = super::AllelePool::new();
        let issues = cfg.validate();
        assert!(issues.contains(&RefDataValidationIssue::EmptyRequiredPool {
            segment: Segment::V,
        }));
    }

    #[test]
    fn empty_j_pool_fails_vj() {
        let mut cfg = minimal_vj();
        cfg.j_pool = super::AllelePool::new();
        let issues = cfg.validate();
        assert!(issues.contains(&RefDataValidationIssue::EmptyRequiredPool {
            segment: Segment::J,
        }));
    }

    #[test]
    fn empty_d_pool_allowed_vj() {
        // VJ chains have empty D pools by construction. This must
        // not flag an EmptyRequiredPool for D.
        let cfg = minimal_vj();
        let issues = cfg.validate();
        assert!(!issues
            .iter()
            .any(|i| matches!(i, RefDataValidationIssue::EmptyRequiredPool { segment } if *segment == Segment::D)));
    }

    #[test]
    fn empty_d_pool_rejected_vdj() {
        let mut cfg = minimal_vdj();
        cfg.d_pool = super::AllelePool::new();
        let issues = cfg.validate();
        assert!(issues.contains(&RefDataValidationIssue::EmptyRequiredPool {
            segment: Segment::D,
        }));
    }

    #[test]
    fn duplicate_v_allele_names_fail() {
        let mut cfg = minimal_vdj();
        let _ = cfg.v_pool.push(v_allele("IGHV1-1*01", b"TGTAAA", Some(0)));
        let issues = cfg.validate();
        assert!(issues.contains(&RefDataValidationIssue::DuplicateAlleleName {
            segment: Segment::V,
            name: "IGHV1-1*01".to_string(),
        }));
    }

    #[test]
    fn duplicate_j_allele_names_fail() {
        let mut cfg = minimal_vdj();
        let _ = cfg.j_pool.push(j_allele("IGHJ1*01", b"TGGAAA", Some(0)));
        let issues = cfg.validate();
        assert!(issues.contains(&RefDataValidationIssue::DuplicateAlleleName {
            segment: Segment::J,
            name: "IGHJ1*01".to_string(),
        }));
    }

    #[test]
    fn invalid_byte_in_v_seq_fails_with_exact_position() {
        let mut cfg = minimal_vdj();
        let _ = cfg.v_pool.push(v_allele("IGHV2*01", b"TGT.AAA", Some(0)));
        let issues = cfg.validate();
        let bad_v_id = AlleleId::new(1);
        let expected = RefDataValidationIssue::InvalidAlleleByte {
            segment: Segment::V,
            allele_id: bad_v_id,
            pos: 3,
            byte: b'.',
        };
        assert!(
            issues.contains(&expected),
            "expected invalid-byte issue at pos 3, got: {issues:?}"
        );
    }

    #[test]
    fn lowercase_bases_accepted() {
        let mut cfg = minimal_vdj();
        let _ = cfg.v_pool.push(v_allele("IGHV3*01", b"tgtaaaccc", Some(0)));
        // Lowercase A/C/G/T/N is allowed. The anchor codon `tgt` still
        // translates to Cys, so V anchor check also passes.
        let issues = cfg.validate();
        assert!(
            issues.is_empty(),
            "lowercase canonical bases should validate, got: {issues:?}"
        );
    }

    #[test]
    fn iupac_ambiguity_other_than_n_fails() {
        let mut cfg = minimal_vdj();
        let _ = cfg.v_pool.push(v_allele("IGHV4*01", b"TGTRAACCC", Some(0)));
        let issues = cfg.validate();
        let v_id = AlleleId::new(1);
        assert!(issues.contains(&RefDataValidationIssue::InvalidAlleleByte {
            segment: Segment::V,
            allele_id: v_id,
            pos: 3,
            byte: b'R',
        }));
    }

    #[test]
    fn v_anchor_out_of_bounds_fails() {
        let mut cfg = minimal_vdj();
        let _ = cfg.v_pool.push(v_allele("IGHV5*01", b"TGT", Some(5)));
        let issues = cfg.validate();
        let v_id = AlleleId::new(1);
        assert!(issues.contains(&RefDataValidationIssue::AnchorOutOfBounds {
            segment: Segment::V,
            allele_id: v_id,
            anchor: 5,
            len: 3,
        }));
    }

    #[test]
    fn v_anchor_at_boundary_with_no_full_codon_fails() {
        // seq.len() == 5, anchor == 3 → codon would span [3, 6) but
        // len is 5 → out of bounds.
        let mut cfg = minimal_vdj();
        let _ = cfg.v_pool.push(v_allele("IGHV6*01", b"AAATG", Some(3)));
        let issues = cfg.validate();
        let v_id = AlleleId::new(1);
        assert!(issues.iter().any(|i| matches!(
            i,
            RefDataValidationIssue::AnchorOutOfBounds { segment, allele_id, .. }
                if *segment == Segment::V && *allele_id == v_id
        )));
    }

    #[test]
    fn v_anchor_non_cys_fails() {
        let mut cfg = minimal_vdj();
        // GGG = Gly, not Cys.
        let _ = cfg.v_pool.push(v_allele("IGHV7*01", b"GGGAAACCC", Some(0)));
        let issues = cfg.validate();
        let v_id = AlleleId::new(1);
        assert!(issues.contains(&RefDataValidationIssue::VAnchorNotCys {
            allele_id: v_id,
            codon: [b'G', b'G', b'G'],
            aa: 'G',
            severity: RefDataIssueSeverity::Curatable,
        }));
    }

    #[test]
    fn v_anchor_tgc_is_cys_and_passes() {
        // TGC is also Cys (degenerate; both TGT and TGC code for C).
        let mut cfg = minimal_vdj();
        let _ = cfg.v_pool.push(v_allele("IGHV8*01", b"TGCAAACCC", Some(0)));
        let issues = cfg.validate();
        assert!(
            issues.is_empty(),
            "TGC anchor codon should pass V Cys check, got: {issues:?}"
        );
    }

    #[test]
    fn j_anchor_unexpected_aa_under_igh_rule_fails() {
        // The catalogue claims IGH (rules.j_anchor.expected = ['W']);
        // a J allele with anchor codon TTC (Phe = F) is flagged.
        // Allele names are no longer used as the source of truth for
        // expected J anchors — the rule on the config is.
        let mut cfg = minimal_vdj();
        cfg.rules = ReferenceRules::for_locus("IGH");
        let _ = cfg.j_pool.push(j_allele("IGHJ2*01", b"TTCAAACCC", Some(0)));
        let issues = cfg.validate();
        let j_id = AlleleId::new(1);
        assert!(issues.iter().any(|i| matches!(
            i,
            RefDataValidationIssue::JAnchorUnexpectedAa { allele_id, aa, .. }
                if *allele_id == j_id && *aa == 'F'
        )));
    }

    #[test]
    fn j_anchor_accepts_f_under_igl_rule() {
        let mut cfg = RefDataConfig::empty(ChainType::Vj);
        cfg.rules = ReferenceRules::for_locus("IGL");
        let _ = cfg.v_pool.push(v_allele("IGLV1*01", b"TGTAAA", Some(0)));
        let _ = cfg.j_pool.push(j_allele("IGLJ1*01", b"TTCAAACCC", Some(0)));
        assert!(cfg.validate().is_empty());
    }

    #[test]
    fn j_anchor_w_under_igl_rule_is_flagged() {
        let mut cfg = RefDataConfig::empty(ChainType::Vj);
        cfg.rules = ReferenceRules::for_locus("IGL");
        let _ = cfg.v_pool.push(v_allele("IGLV1*01", b"TGTAAA", Some(0)));
        // IGL rule expects F — a W-coded J anchor is flagged.
        let _ = cfg.j_pool.push(j_allele("IGLJ2*01", b"TGGAAACCC", Some(0)));
        let issues = cfg.validate();
        assert!(issues.iter().any(|i| matches!(
            i,
            RefDataValidationIssue::JAnchorUnexpectedAa { aa, .. } if *aa == 'W'
        )));
    }

    #[test]
    fn default_rule_accepts_w_or_f_for_unconfigured_locus() {
        // `RefDataConfig::empty` initialises with the lenient default
        // (J expects W or F). Any allele-name pattern that the loader
        // doesn't recognise lands here.
        let mut cfg = RefDataConfig::empty(ChainType::Vj);
        let _ = cfg.v_pool.push(v_allele("MYV1*01", b"TGTAAA", Some(0)));
        let _ = cfg.j_pool.push(j_allele("MYJ1*01", b"TGGAAACCC", Some(0)));
        assert!(
            cfg.validate().is_empty(),
            "default-locus J anchor 'W' should be accepted"
        );

        let mut cfg2 = RefDataConfig::empty(ChainType::Vj);
        let _ = cfg2.v_pool.push(v_allele("MYV1*01", b"TGTAAA", Some(0)));
        let _ = cfg2.j_pool.push(j_allele("MYJ2*01", b"TTCAAACCC", Some(0)));
        assert!(
            cfg2.validate().is_empty(),
            "default-locus J anchor 'F' should be accepted"
        );
    }

    #[test]
    fn missing_anchor_on_v_is_reported() {
        let mut cfg = minimal_vdj();
        let _ = cfg.v_pool.push(v_allele("IGHV9*01", b"AAACCCGGG", None));
        let issues = cfg.validate();
        let v_id = AlleleId::new(1);
        assert!(issues.contains(&RefDataValidationIssue::MissingAnchor {
            segment: Segment::V,
            allele_id: v_id,
            severity: RefDataIssueSeverity::Curatable,
        }));
    }

    #[test]
    fn missing_anchor_on_j_is_reported() {
        let mut cfg = minimal_vdj();
        let _ = cfg.j_pool.push(j_allele("IGHJ-orphan*01", b"GGGAAA", None));
        let issues = cfg.validate();
        let j_id = AlleleId::new(1);
        assert!(issues.contains(&RefDataValidationIssue::MissingAnchor {
            segment: Segment::J,
            allele_id: j_id,
            severity: RefDataIssueSeverity::Curatable,
        }));
    }

    #[test]
    fn d_anchor_absence_is_not_reported() {
        // The default minimal_vdj D allele has no anchor; this must
        // not flag a MissingAnchor issue (D anchors are not part of
        // the contract).
        let cfg = minimal_vdj();
        let issues = cfg.validate();
        assert!(!issues
            .iter()
            .any(|i| matches!(i, RefDataValidationIssue::MissingAnchor { segment, .. } if *segment == Segment::D)));
    }

    #[test]
    fn validate_strict_packages_errors() {
        let mut cfg = minimal_vdj();
        let _ = cfg.v_pool.push(v_allele("IGHV1-1*01", b"TGTAAA", Some(0)));
        let err = cfg.validate_strict().expect_err("duplicate must fail strict");
        assert_eq!(err.issues.len(), 1);
        assert!(format!("{err}").contains("duplicate allele name"));
    }

    #[test]
    fn validate_strict_collects_all_issues_not_just_first() {
        let mut cfg = RefDataConfig::empty(ChainType::Vdj);
        // Multiple problems at once: empty V, empty J, empty D.
        let err = cfg.validate_strict().expect_err("all empty should fail");
        let segs: Vec<_> = err
            .issues
            .iter()
            .filter_map(|i| match i {
                RefDataValidationIssue::EmptyRequiredPool { segment } => Some(*segment),
                _ => None,
            })
            .collect();
        assert!(segs.contains(&Segment::V));
        assert!(segs.contains(&Segment::J));
        assert!(segs.contains(&Segment::D));

        // Plug V, leave D and J empty: still has two issues.
        let _ = cfg.v_pool.push(v_allele("IGHV1*01", b"TGTAAA", Some(0)));
        let err2 = cfg.validate_strict().expect_err("two still missing");
        assert_eq!(err2.issues.len(), 2);
    }

    // ── Severity classification + mode-aware gating ───────────────

    /// VJ refdata with only curatable issues (Gly V anchor + no J
    /// anchor). Useful to verify mode toggling without bleed-over
    /// from Fatal issues.
    fn vj_curatable_only() -> RefDataConfig {
        let mut cfg = RefDataConfig::empty(ChainType::Vj);
        let _ = cfg.v_pool.push(v_allele("IGKV1-1*01", b"AAACCCGGG", Some(6)));
        let _ = cfg.j_pool.push(Allele {
            name: "IGKJ-orphan*01".into(),
            gene: "IGKJ-orphan".into(),
            seq: b"TTCAAACCC".to_vec(),
            segment: Segment::J,
            anchor: None,
            functional_status: None,
            subregions: Vec::new(),
        });
        cfg
    }

    #[test]
    fn validate_with_strict_mode_rejects_curatable_issues() {
        let cfg = vj_curatable_only();
        let err = cfg
            .validate_with_mode(RefDataValidationMode::Strict)
            .expect_err("strict must reject curatable issues");
        let (fatal, curatable) = err.severity_counts();
        assert_eq!(fatal, 0);
        assert!(curatable >= 2);
        assert_eq!(err.mode, RefDataValidationMode::Strict);
    }

    #[test]
    fn validate_with_allow_curatable_accepts_curatable_only() {
        let cfg = vj_curatable_only();
        cfg.validate_with_mode(RefDataValidationMode::AllowCuratable)
            .expect("allow_curatable must accept curatable-only fixtures");
    }

    #[test]
    fn validate_with_allow_curatable_still_rejects_fatal() {
        let mut cfg = vj_curatable_only();
        // Add an invalid byte allele (Fatal).
        let _ = cfg.v_pool.push(v_allele("IGKV2*01", b"TGT.AAACC", Some(0)));
        let err = cfg
            .validate_with_mode(RefDataValidationMode::AllowCuratable)
            .expect_err("allow_curatable must still reject Fatal issues");
        let (fatal, curatable) = err.severity_counts();
        assert!(fatal >= 1, "expected ≥1 fatal, got: {:?}", err.issues);
        assert!(curatable >= 1, "expected ≥1 curatable preserved, got: {:?}", err.issues);
        assert_eq!(err.mode, RefDataValidationMode::AllowCuratable);
    }

    #[test]
    fn display_remediation_hint_when_only_curatable_issues_remain() {
        let cfg = vj_curatable_only();
        let err = cfg.validate_strict().expect_err("must fail");
        let msg = format!("{err}");
        assert!(
            msg.contains("allow_curatable_refdata"),
            "msg must suggest allow_curatable_refdata when only curatable: {msg}"
        );
        assert!(
            msg.contains("pseudogene/ORF"),
            "msg must mention pseudogene/ORF: {msg}"
        );
    }

    #[test]
    fn display_omits_remediation_hint_when_fatal_present() {
        let mut cfg = vj_curatable_only();
        let _ = cfg.v_pool.push(v_allele("IGKV2*01", b"TGT.AAACC", Some(0)));
        let err = cfg.validate_strict().expect_err("must fail");
        let msg = format!("{err}");
        // Fatal issues exist — opting into allow_curatable_refdata
        // wouldn't help, so the hint must NOT appear.
        assert!(
            !msg.contains("allow_curatable_refdata"),
            "msg must NOT suggest allow_curatable_refdata when Fatal issues remain: {msg}"
        );
    }

    #[test]
    fn each_issue_variant_carries_documented_severity() {
        use RefDataIssueSeverity::*;
        let aid = AlleleId::new(0);
        // Fatal set.
        assert_eq!(
            RefDataValidationIssue::EmptyRequiredPool { segment: Segment::V }.severity(),
            Fatal
        );
        assert_eq!(
            RefDataValidationIssue::DuplicateAlleleName {
                segment: Segment::J,
                name: "x".into()
            }
            .severity(),
            Fatal
        );
        assert_eq!(
            RefDataValidationIssue::InvalidAlleleByte {
                segment: Segment::V,
                allele_id: aid,
                pos: 0,
                byte: b'.'
            }
            .severity(),
            Fatal
        );
        assert_eq!(
            RefDataValidationIssue::AnchorOutOfBounds {
                segment: Segment::V,
                allele_id: aid,
                anchor: 99,
                len: 5
            }
            .severity(),
            Fatal
        );
        // Rule-controlled set — severity is whatever the construction
        // site stamped on the issue (`severity()` just returns it).
        assert_eq!(
            RefDataValidationIssue::VAnchorNotCys {
                allele_id: aid,
                codon: [b'G', b'G', b'G'],
                aa: 'G',
                severity: Curatable,
            }
            .severity(),
            Curatable
        );
        assert_eq!(
            RefDataValidationIssue::JAnchorUnexpectedAa {
                allele_id: aid,
                codon: [b'T', b'T', b'A'],
                aa: 'L',
                expected: vec!['W'],
                severity: Curatable,
            }
            .severity(),
            Curatable
        );
        assert_eq!(
            RefDataValidationIssue::MissingAnchor {
                segment: Segment::V,
                allele_id: aid,
                severity: Curatable,
            }
            .severity(),
            Curatable
        );
        // The rule-controlled variants honour the severity carried
        // on the issue itself — a Fatal anchor mismatch (custom rule)
        // reports as Fatal.
        assert_eq!(
            RefDataValidationIssue::MissingAnchor {
                segment: Segment::V,
                allele_id: aid,
                severity: Fatal,
            }
            .severity(),
            Fatal
        );
    }

    // ── ReferenceRules v1 — configurable anchor + alphabet ──────

    #[test]
    fn default_rules_match_legacy_behavior() {
        // V → ['C'], J → ['W', 'F'], alphabet ACGTN, all Curatable.
        let rules = ReferenceRules::default();
        assert_eq!(rules.v_anchor.expected_amino_acids, vec!['C']);
        assert_eq!(rules.j_anchor.expected_amino_acids, vec!['W', 'F']);
        assert!(rules.v_anchor.required);
        assert!(rules.j_anchor.required);
        for &b in b"acgtnACGTN" {
            assert!(rules.alphabet.is_allowed(b), "should allow {}", b as char);
        }
        assert!(!rules.alphabet.is_allowed(b'.'));
        assert!(!rules.alphabet.is_allowed(b'R'));
    }

    #[test]
    fn for_locus_returns_locus_appropriate_j_rule() {
        assert_eq!(
            ReferenceRules::for_locus("IGH").j_anchor.expected_amino_acids,
            vec!['W']
        );
        assert_eq!(
            ReferenceRules::for_locus("IGK").j_anchor.expected_amino_acids,
            vec!['F']
        );
        assert_eq!(
            ReferenceRules::for_locus("IGL").j_anchor.expected_amino_acids,
            vec!['F']
        );
        for tr in ["TRA", "TRB", "TRG", "TRD"] {
            assert_eq!(
                ReferenceRules::for_locus(tr).j_anchor.expected_amino_acids,
                vec!['F'],
                "{tr} J should expect F"
            );
        }
        // Unknown prefix → default lenient set.
        assert_eq!(
            ReferenceRules::for_locus("XXX").j_anchor.expected_amino_acids,
            vec!['W', 'F']
        );
    }

    #[test]
    fn custom_j_rule_y_accepts_tat_and_tac() {
        let mut cfg = RefDataConfig::empty(ChainType::Vj);
        cfg.rules.j_anchor.expected_amino_acids = vec!['Y'];
        let _ = cfg.v_pool.push(v_allele("MYV1*01", b"TGTAAA", Some(0)));
        // TAT → Y.
        let _ = cfg.j_pool.push(j_allele("MYJ1*01", b"TATAAACCC", Some(0)));
        assert!(cfg.validate().is_empty());
        // TAC → Y.
        let mut cfg2 = cfg.clone();
        cfg2.j_pool = crate::refdata::AllelePool::new();
        let _ = cfg2.j_pool.push(j_allele("MYJ2*01", b"TACAAACCC", Some(0)));
        assert!(cfg2.validate().is_empty());
    }

    #[test]
    fn custom_j_rule_y_rejects_tgg() {
        let mut cfg = RefDataConfig::empty(ChainType::Vj);
        cfg.rules.j_anchor.expected_amino_acids = vec!['Y'];
        let _ = cfg.v_pool.push(v_allele("MYV1*01", b"TGTAAA", Some(0)));
        // TGG → W. Under a Y-only J rule, this is flagged.
        let _ = cfg.j_pool.push(j_allele("MYJ-bad*01", b"TGGAAACCC", Some(0)));
        let issues = cfg.validate();
        assert!(issues.iter().any(|i| matches!(
            i,
            RefDataValidationIssue::JAnchorUnexpectedAa { aa, .. } if *aa == 'W'
        )));
    }

    #[test]
    fn anchor_required_false_suppresses_missing_anchor_issue() {
        // Build a config where the J rule says anchor is optional.
        let mut cfg = RefDataConfig::empty(ChainType::Vj);
        cfg.rules.j_anchor.required = false;
        let _ = cfg.v_pool.push(v_allele("MYV1*01", b"TGTAAA", Some(0)));
        // J allele with no anchor → would normally emit MissingAnchor.
        let _ = cfg.j_pool.push(j_allele("MYJ-orphan*01", b"GGG", None));
        assert!(
            cfg.validate().is_empty(),
            "rule.required=false must suppress MissingAnchor"
        );
    }

    #[test]
    fn missing_severity_follows_rule() {
        // Bump missing_severity to Fatal — the resulting issue's
        // severity() reflects that, and AllowCuratable mode still
        // rejects.
        let mut cfg = RefDataConfig::empty(ChainType::Vj);
        cfg.rules.j_anchor.missing_severity = RefDataIssueSeverity::Fatal;
        let _ = cfg.v_pool.push(v_allele("MYV1*01", b"TGTAAA", Some(0)));
        let _ = cfg.j_pool.push(j_allele("MYJ*01", b"TGGAAA", None));
        let issues = cfg.validate();
        let missing = issues
            .iter()
            .find(|i| matches!(i, RefDataValidationIssue::MissingAnchor { .. }))
            .expect("must emit MissingAnchor");
        assert_eq!(missing.severity(), RefDataIssueSeverity::Fatal);

        // AllowCuratable still rejects.
        cfg.validate_with_mode(RefDataValidationMode::AllowCuratable)
            .expect_err("Fatal missing-anchor must not be opt-outable");
    }

    #[test]
    fn mismatch_severity_follows_rule() {
        let mut cfg = RefDataConfig::empty(ChainType::Vj);
        cfg.rules.v_anchor.mismatch_severity = RefDataIssueSeverity::Fatal;
        // GGG (Gly) at V anchor.
        let _ = cfg.v_pool.push(v_allele("MYV-bad*01", b"GGGAAACCC", Some(0)));
        let _ = cfg.j_pool.push(j_allele("MYJ*01", b"TGGAAA", Some(0)));
        let issues = cfg.validate();
        let mismatch = issues
            .iter()
            .find(|i| matches!(i, RefDataValidationIssue::VAnchorNotCys { .. }))
            .expect("must emit VAnchorNotCys");
        assert_eq!(mismatch.severity(), RefDataIssueSeverity::Fatal);

        // AllowCuratable still rejects this catalogue.
        cfg.validate_with_mode(RefDataValidationMode::AllowCuratable)
            .expect_err("Fatal V anchor mismatch must not be opt-outable");
    }

    #[test]
    fn invalid_byte_remains_fatal_regardless_of_rule() {
        // Even with an extremely permissive anchor rule, an invalid
        // sequence byte stays Fatal — structural problems aren't
        // rule-controlled.
        let mut cfg = RefDataConfig::empty(ChainType::Vj);
        cfg.rules.v_anchor.expected_amino_acids = vec!['A', 'C', 'D', 'E', 'F'];
        cfg.rules.j_anchor.expected_amino_acids = vec!['A', 'C', 'D', 'E', 'F'];
        let _ = cfg.v_pool.push(v_allele("MYV-bad*01", b"TGT.AAACC", Some(0)));
        let _ = cfg.j_pool.push(j_allele("MYJ*01", b"TGGAAA", Some(0)));
        let issues = cfg.validate();
        let bad = issues
            .iter()
            .find(|i| matches!(i, RefDataValidationIssue::InvalidAlleleByte { .. }))
            .expect("must emit InvalidAlleleByte");
        assert_eq!(bad.severity(), RefDataIssueSeverity::Fatal);
    }

    #[test]
    fn custom_alphabet_extends_allowed_set() {
        let mut cfg = RefDataConfig::empty(ChainType::Vj);
        cfg.rules.alphabet = ReferenceAlphabet {
            allowed: vec![b'A', b'C', b'G', b'T', b'N', b'R'], // R = puRine ambig
        };
        let _ = cfg.v_pool.push(v_allele("MYV*01", b"TGTRAA", Some(0)));
        let _ = cfg.j_pool.push(j_allele("MYJ*01", b"TGGAAA", Some(0)));
        let issues = cfg.validate();
        // R is now in-alphabet → no InvalidAlleleByte issue.
        assert!(
            !issues
                .iter()
                .any(|i| matches!(i, RefDataValidationIssue::InvalidAlleleByte { .. })),
            "R should be allowed under the extended alphabet; got: {issues:?}"
        );
    }

    #[test]
    fn reference_alphabet_is_case_insensitive() {
        let alphabet = ReferenceAlphabet::default();
        assert!(alphabet.is_allowed(b'a'));
        assert!(alphabet.is_allowed(b'A'));
        assert!(alphabet.is_allowed(b'n'));
        assert!(!alphabet.is_allowed(b'.'));
    }

    // ── Identity ↔ chain-type mismatch (Reference Identity slice) ──

    fn vdj_cfg_with_identity_locus(locus: &str) -> RefDataConfig {
        let mut cfg = minimal_vdj();
        cfg.identity.locus = Some(locus.to_string());
        cfg
    }

    fn vj_cfg_with_identity_locus(locus: &str) -> RefDataConfig {
        let mut cfg = RefDataConfig::empty(ChainType::Vj);
        let _ = cfg.v_pool.push(v_allele("IGKV1-1*01", b"TGTAAACCC", Some(0)));
        let _ = cfg.j_pool.push(j_allele("IGKJ1*01", b"TTCAAACCC", Some(0)));
        cfg.identity.locus = Some(locus.to_string());
        cfg
    }

    #[test]
    fn vj_cartridge_with_igh_locus_is_flagged() {
        // VJ topology declaring an IGH (VDJ) locus is structurally
        // wrong — recombination shape and locus disagree.
        let cfg = vj_cfg_with_identity_locus("IGH");
        let issues = cfg.validate();
        assert!(issues.iter().any(|i| matches!(
            i,
            RefDataValidationIssue::LocusChainTypeMismatch { locus, .. } if locus == "IGH"
        )));
    }

    #[test]
    fn vj_cartridge_with_igk_locus_passes() {
        // Locus matches topology — no LocusChainTypeMismatch.
        let cfg = vj_cfg_with_identity_locus("IGK");
        let issues = cfg.validate();
        assert!(!issues.iter().any(|i| matches!(
            i,
            RefDataValidationIssue::LocusChainTypeMismatch { .. }
        )));
    }

    #[test]
    fn vdj_cartridge_with_igk_locus_is_flagged() {
        let cfg = vdj_cfg_with_identity_locus("IGK");
        assert!(cfg.validate().iter().any(|i| matches!(
            i,
            RefDataValidationIssue::LocusChainTypeMismatch { locus, .. } if locus == "IGK"
        )));
    }

    #[test]
    fn locus_chain_type_mismatch_is_fatal() {
        let cfg = vj_cfg_with_identity_locus("IGH");
        let issue = cfg
            .validate()
            .into_iter()
            .find(|i| matches!(i, RefDataValidationIssue::LocusChainTypeMismatch { .. }))
            .expect("must emit LocusChainTypeMismatch");
        assert_eq!(issue.severity(), RefDataIssueSeverity::Fatal);
        // AllowCuratable mode cannot opt this out.
        cfg.validate_with_mode(RefDataValidationMode::AllowCuratable)
            .expect_err("Fatal mismatch must not be opt-outable");
    }

    #[test]
    fn unknown_locus_does_not_fail_validation() {
        let cfg = vj_cfg_with_identity_locus("XYZ");
        assert!(!cfg.validate().iter().any(|i| matches!(
            i,
            RefDataValidationIssue::LocusChainTypeMismatch { .. }
        )));
    }

    #[test]
    fn empty_identity_does_not_fail_validation() {
        // No locus declared → no mismatch issue regardless of chain.
        let cfg_vj = RefDataConfig::empty(ChainType::Vj);
        assert!(!cfg_vj.validate().iter().any(|i| matches!(
            i,
            RefDataValidationIssue::LocusChainTypeMismatch { .. }
        )));
        let cfg_vdj = RefDataConfig::empty(ChainType::Vdj);
        assert!(!cfg_vdj.validate().iter().any(|i| matches!(
            i,
            RefDataValidationIssue::LocusChainTypeMismatch { .. }
        )));
    }

    #[test]
    fn locus_is_case_insensitive() {
        // Lowercase / mixed-case locus still flags correctly.
        let cfg = vj_cfg_with_identity_locus("igh");
        let issues = cfg.validate();
        assert!(issues.iter().any(|i| matches!(
            i,
            RefDataValidationIssue::LocusChainTypeMismatch { locus, .. } if locus == "IGH"
        )));
    }
