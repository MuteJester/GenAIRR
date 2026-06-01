use super::super::test_support::{
    make_assembled_sim_from_refdata, make_v_anchor_at, make_vj_with_anchor_codons,
};
use super::*;
use crate::address::ChoiceAddress;
use crate::assignment::AlleleInstance;
use crate::contract::ChoiceContext;
use crate::ir::{flag, Nucleotide};
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

    // No assignment → vacuously satisfied.
    assert!(contract.verify(&sim, Some(&cfg)).is_ok());
}

#[test]
fn anchor_preserved_no_refdata_passes_vacuously() {
    let cfg = make_v_anchor_at(10, Some(5));
    let contract = AnchorPreserved::new(Segment::V);
    let sim =
        Simulation::new().with_allele_assigned(Segment::V, AlleleInstance::new(AlleleId::new(0)));

    // Refdata None → can't verify; treat as satisfied.
    assert!(contract.verify(&sim, None).is_ok());
    // Sanity: WITH refdata it would also pass (no trims).
    assert!(contract.verify(&sim, Some(&cfg)).is_ok());
}

#[test]
fn anchor_preserved_anchorless_allele_passes_vacuously() {
    let cfg = make_v_anchor_at(10, None); // anchorless
    let contract = AnchorPreserved::new(Segment::V);
    let sim =
        Simulation::new().with_allele_assigned(Segment::V, AlleleInstance::new(AlleleId::new(0)));

    assert!(contract.verify(&sim, Some(&cfg)).is_ok());
}

#[test]
fn anchor_preserved_zero_trims_and_anchor_in_range_passes() {
    // 10-base allele, anchor at 3. No trims. Anchor codon spans
    // [3, 6) which fits entirely in [0, 10).
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

    let context = ChoiceContext::targeted_base_substitution(0, 1, NucHandle::new(6))
        .with_address_if_missing(ChoiceAddress::parse("mutate.uniform.base[0]"));
    let err = contract
        .admits_typed(&sim, Some(&cfg), context, &ChoiceValue::Base(b'A'))
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
    // Anchor at 3. trim_5 = 4 → anchor 5'-trimmed away.
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
    // 10-base allele, anchor at 6. Anchor codon spans [6, 9).
    // trim_3 = 2 → retained_end = 8, so anchor codon's last
    // position (8) is at the boundary — anchor + 3 = 9 > 8 → fail.
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
    // anchor at 0 → trim_5 = 0 leaves anchor at the very start.
    let cfg = make_v_anchor_at(10, Some(0));
    let contract = AnchorPreserved::new(Segment::V);
    let sim =
        Simulation::new().with_allele_assigned(Segment::V, AlleleInstance::new(AlleleId::new(0)));

    assert!(contract.verify(&sim, Some(&cfg)).is_ok());
}

#[test]
fn anchor_preserved_anchor_at_3prime_boundary_passes() {
    // 10-base allele, anchor at 7. Anchor codon spans [7, 10).
    // trim_3 = 0 → retained_end = 10. Boundary case — passes.
    let cfg = make_v_anchor_at(10, Some(7));
    let contract = AnchorPreserved::new(Segment::V);
    let sim =
        Simulation::new().with_allele_assigned(Segment::V, AlleleInstance::new(AlleleId::new(0)));

    assert!(contract.verify(&sim, Some(&cfg)).is_ok());
}

#[test]
fn anchor_preserved_one_more_trim_just_past_boundary_violates() {
    // Same as above but trim_3 = 1 → retained_end = 9, anchor
    // codon end (10) > 9 → fail.
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
    // Build a J pool with an anchor at position 0 (typical J
    // anchor is near the 5' end). trim_5 = 0 → anchor preserved.
    // trim_5 = 1 → anchor 5'-trimmed → fail.
    use crate::assignment::TrimEnd;

    let mut cfg = RefDataConfig::empty(ChainType::Vdj);
    let _ = cfg.j_pool.push(Allele {
        name: "j_test*01".to_string(),
        gene: "j_test".to_string(),
        seq: vec![b'A'; 12],
        segment: Segment::J,
        anchor: Some(0),
        functional_status: None,
        subregions: Vec::new(),
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
    // Trait-object usage is dyn-compatible.
    let cfg = make_v_anchor_at(10, Some(3));
    let contract: Box<dyn Contract> = Box::new(AnchorPreserved::new(Segment::V));
    let sim =
        Simulation::new().with_allele_assigned(Segment::V, AlleleInstance::new(AlleleId::new(0)));

    assert!(contract.verify(&sim, Some(&cfg)).is_ok());
    assert_eq!(contract.name(), "anchor_preserved.v");
}

// v3.0 review fix: pin the `admits_fixed_base_at` override against
// drift. The default trait impl returns `true` for every site/byte;
// the override substitutes the byte at the site, translates the
// resulting anchor codon, and admits only if the reference amino
// acid is preserved.
mod admits_fixed {
    use super::super::super::test_support::{
        make_assembled_sim_from_refdata, make_vj_with_anchor_codons,
    };
    use super::*;
    use crate::contract::{IndelEventClass, IndelKindHint};

    /// Anchor codon `TAC` (Y) at V[6..9]. Y has two synonyms —
    /// TAC and TAT — so at codon position 2 (pool position 8)
    /// only `C` and `T` preserve Y; `A` (→ TAA = stop) and
    /// `G` (→ TAG = stop) do not.
    fn fixture() -> (RefDataConfig, Simulation) {
        let cfg = make_vj_with_anchor_codons(b"TAC", b"GGG");
        let sim = make_assembled_sim_from_refdata(&cfg);
        (cfg, sim)
    }

    #[test]
    fn admits_fixed_base_at_synonym_codon_byte_is_admitted() {
        let (cfg, sim) = fixture();
        let c = AnchorPreserved::new(Segment::V);
        assert!(c.admits_fixed_base_at(&sim, Some(&cfg), NucHandle::new(8), b'C'));
        assert!(c.admits_fixed_base_at(&sim, Some(&cfg), NucHandle::new(8), b'T'));
    }

    #[test]
    fn admits_fixed_base_at_non_synonym_byte_is_rejected() {
        let (cfg, sim) = fixture();
        let c = AnchorPreserved::new(Segment::V);
        assert!(!c.admits_fixed_base_at(&sim, Some(&cfg), NucHandle::new(8), b'A'));
        assert!(!c.admits_fixed_base_at(&sim, Some(&cfg), NucHandle::new(8), b'G'));
    }

    #[test]
    fn admits_fixed_base_at_lowercase_matches_uppercase_verdict() {
        let (cfg, sim) = fixture();
        let c = AnchorPreserved::new(Segment::V);
        assert!(c.admits_fixed_base_at(&sim, Some(&cfg), NucHandle::new(8), b'c'));
        assert!(!c.admits_fixed_base_at(&sim, Some(&cfg), NucHandle::new(8), b'a'));
    }

    #[test]
    fn admits_fixed_base_at_non_canonical_byte_at_anchor_is_rejected() {
        // `N` translates to `X`, not the reference amino acid →
        // anchor not preserved → rejected.
        let (cfg, sim) = fixture();
        let c = AnchorPreserved::new(Segment::V);
        assert!(!c.admits_fixed_base_at(&sim, Some(&cfg), NucHandle::new(8), b'N'));
    }

    #[test]
    fn admits_fixed_base_at_out_of_anchor_codon_is_unconstrained() {
        let (cfg, sim) = fixture();
        let c = AnchorPreserved::new(Segment::V);
        assert!(c.admits_fixed_base_at(&sim, Some(&cfg), NucHandle::new(0), b'A'));
        assert!(c.admits_fixed_base_at(&sim, Some(&cfg), NucHandle::new(9), b'A'));
    }

    #[test]
    fn admits_fixed_base_at_no_refdata_is_unconstrained() {
        let (_cfg, sim) = fixture();
        let c = AnchorPreserved::new(Segment::V);
        assert!(c.admits_fixed_base_at(&sim, None, NucHandle::new(8), b'A'));
    }

    #[test]
    fn admits_fixed_base_at_agrees_with_mask_for_canonical_bases() {
        // Pin the v3.0 invariant: mask and fixed-base API
        // agree on canonical A/C/G/T at every site.
        let (cfg, sim) = fixture();
        let c = AnchorPreserved::new(Segment::V);
        for site in 0..(sim.pool.len() as u32) {
            let handle = NucHandle::new(site);
            let mask = c.admissible_bases_at(&sim, Some(&cfg), handle);
            for &b in &[b'A', b'C', b'G', b'T'] {
                assert_eq!(
                    mask.admits(b),
                    c.admits_fixed_base_at(&sim, Some(&cfg), handle, b),
                    "mismatch at site {} base {}",
                    site,
                    b as char
                );
            }
        }
    }

    #[test]
    fn admissible_indel_class_inside_anchor_codon_is_forbidden() {
        let (cfg, sim) = fixture();
        let c = AnchorPreserved::new(Segment::V);
        for site in 6..9u32 {
            assert_eq!(
                c.admissible_indel_class_at(&sim, Some(&cfg), site, IndelKindHint::Insertion,),
                IndelEventClass::Forbidden,
                "insertion at site {} should be Forbidden",
                site
            );
            assert_eq!(
                c.admissible_indel_class_at(&sim, Some(&cfg), site, IndelKindHint::Deletion,),
                IndelEventClass::Forbidden,
                "deletion at site {} should be Forbidden",
                site
            );
        }
    }

    #[test]
    fn admissible_indel_class_inside_region_before_anchor_is_forbidden() {
        // Indels at sites in [region.start, anchor_codon_end)
        // shift the anchor codon's bytes → Forbidden.
        let (cfg, sim) = fixture();
        let c = AnchorPreserved::new(Segment::V);
        for site in 0..6u32 {
            assert_eq!(
                c.admissible_indel_class_at(&sim, Some(&cfg), site, IndelKindHint::Insertion,),
                IndelEventClass::Forbidden,
            );
            assert_eq!(
                c.admissible_indel_class_at(&sim, Some(&cfg), site, IndelKindHint::Deletion,),
                IndelEventClass::Forbidden,
            );
        }
    }

    #[test]
    fn admissible_indel_class_past_anchor_end_is_frame_neutral() {
        let (cfg, sim) = fixture();
        let c = AnchorPreserved::new(Segment::V);
        for site in 9..=15u32 {
            assert_eq!(
                c.admissible_indel_class_at(&sim, Some(&cfg), site, IndelKindHint::Insertion,),
                IndelEventClass::FrameNeutral,
            );
            assert_eq!(
                c.admissible_indel_class_at(&sim, Some(&cfg), site, IndelKindHint::Deletion,),
                IndelEventClass::FrameNeutral,
            );
        }
    }

    #[test]
    fn admissible_indel_class_unassembled_segment_is_frame_neutral() {
        let cfg = make_vj_with_anchor_codons(b"TAC", b"GGG");
        let c = AnchorPreserved::new(Segment::V);
        let sim = Simulation::new();
        assert_eq!(
            c.admissible_indel_class_at(&sim, Some(&cfg), 5, IndelKindHint::Insertion,),
            IndelEventClass::FrameNeutral
        );
    }

    // v3.0 Phase E: trim-length support overrides.
    //
    // Fixture: V allele "AAACCC" + "TAC" = "AAACCCTAC" (9
    // bytes), anchor at offset 6 → anchor codon "TAC" (Y).
    // Max admissible 5' trim = anchor = 6. Max admissible 3'
    // trim = allele_len − anchor − 3 = 9 − 6 − 3 = 0.
    mod admissible_trim_lengths {
        use super::*;
        use crate::contract::{LengthSupport, TrimTarget};

        fn fixture() -> (RefDataConfig, Simulation) {
            let cfg = make_vj_with_anchor_codons(b"TAC", b"GGG");
            // The default fixture only has alleles assigned, no
            // assembled regions — admissible_trim_lengths reads
            // assignment metadata only, so we don't need
            // assembled regions for the override.
            let sim = Simulation::new()
                .with_allele_assigned(Segment::V, AlleleInstance::new(AlleleId::new(0)))
                .with_allele_assigned(Segment::J, AlleleInstance::new(AlleleId::new(0)));
            (cfg, sim)
        }

        #[test]
        fn five_prime_trim_max_is_anchor_offset() {
            let (cfg, sim) = fixture();
            let c = AnchorPreserved::new(Segment::V);
            // Anchor at offset 6 → max trim_5 = 6.
            let support = c.admissible_trim_lengths(
                &sim,
                Some(&cfg),
                TrimTarget {
                    segment: Segment::V,
                    end: TrimEnd::Five,
                },
                100,
            );
            assert_eq!(support, LengthSupport::Full(6));
        }

        #[test]
        fn three_prime_trim_max_is_allele_len_minus_anchor_minus_three() {
            let (cfg, sim) = fixture();
            let c = AnchorPreserved::new(Segment::V);
            // V allele "AAACCCTAC" len=9, anchor=6
            // → max trim_3 = 9 − 6 − 3 = 0.
            let support = c.admissible_trim_lengths(
                &sim,
                Some(&cfg),
                TrimTarget {
                    segment: Segment::V,
                    end: TrimEnd::Three,
                },
                100,
            );
            assert_eq!(support, LengthSupport::Full(0));
        }

        #[test]
        fn other_segment_target_is_unconstrained() {
            let (cfg, sim) = fixture();
            // V-anchor contract has no opinion on J trim.
            let c = AnchorPreserved::new(Segment::V);
            let support = c.admissible_trim_lengths(
                &sim,
                Some(&cfg),
                TrimTarget {
                    segment: Segment::J,
                    end: TrimEnd::Five,
                },
                100,
            );
            assert_eq!(support, LengthSupport::Full(100));
        }

        #[test]
        fn no_refdata_is_unconstrained() {
            let (_cfg, sim) = fixture();
            let c = AnchorPreserved::new(Segment::V);
            let support = c.admissible_trim_lengths(
                &sim,
                None,
                TrimTarget {
                    segment: Segment::V,
                    end: TrimEnd::Five,
                },
                100,
            );
            assert_eq!(support, LengthSupport::Full(100));
        }

        #[test]
        fn no_assignment_is_unconstrained() {
            let (cfg, _sim) = fixture();
            let c = AnchorPreserved::new(Segment::V);
            let support = c.admissible_trim_lengths(
                &Simulation::new(),
                Some(&cfg),
                TrimTarget {
                    segment: Segment::V,
                    end: TrimEnd::Five,
                },
                100,
            );
            assert_eq!(support, LengthSupport::Full(100));
        }

        #[test]
        fn requested_max_caps_admissible_when_lower() {
            let (cfg, sim) = fixture();
            let c = AnchorPreserved::new(Segment::V);
            // anchor=6 but request max=3 → support capped at 3.
            let support = c.admissible_trim_lengths(
                &sim,
                Some(&cfg),
                TrimTarget {
                    segment: Segment::V,
                    end: TrimEnd::Five,
                },
                3,
            );
            assert_eq!(support, LengthSupport::Full(3));
        }

        #[test]
        fn five_prime_support_empty_when_three_already_destroyed() {
            // Set up: trim_3 large enough to destroy the
            // anchor → 5' support is Empty (defensive — the
            // current trim_3 already eliminated the anchor
            // codon regardless of the candidate L).
            let cfg = make_vj_with_anchor_codons(b"TAC", b"GGG");
            let sim = Simulation::new()
                .with_allele_assigned(Segment::V, AlleleInstance::new(AlleleId::new(0)))
                // V allele len=9, anchor=6, anchor_end=9.
                // trim_3=7 → retained_end=2 → anchor (6)+3=9 > 2.
                .with_trim(Segment::V, TrimEnd::Three, 7)
                .with_allele_assigned(Segment::J, AlleleInstance::new(AlleleId::new(0)));
            let c = AnchorPreserved::new(Segment::V);
            let support = c.admissible_trim_lengths(
                &sim,
                Some(&cfg),
                TrimTarget {
                    segment: Segment::V,
                    end: TrimEnd::Five,
                },
                100,
            );
            assert_eq!(support, LengthSupport::Empty);
        }
    }
}
