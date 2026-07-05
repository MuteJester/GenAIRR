    use super::*;

    #[test]
    fn segment_helpers_emit_existing_recombination_addresses() {
        assert_eq!(sample_allele_vdj(Segment::V), "sample_allele.v");
        assert_eq!(sample_allele_vdj(Segment::D), "sample_allele.d");
        assert_eq!(sample_allele_vdj(Segment::J), "sample_allele.j");
        assert_eq!(trim_vdj(Segment::V, TrimEnd::Five), "trim.v_5");
        assert_eq!(trim_vdj(Segment::J, TrimEnd::Three), "trim.j_3");
        assert_eq!(assemble_vdj(Segment::D), "assemble.d");
        assert_eq!(np_length_region(Segment::Np1), "np.np1.length");
    }

    #[test]
    fn parsers_accept_existing_indexed_addresses() {
        assert_eq!(
            parse_indexed("mutate.s5f.base[2]", MUTATE_S5F_BASE_PREFIX),
            Some(2)
        );
        assert_eq!(
            ChoiceAddress::parse("np.np1.bases[17]"),
            Some(ChoiceAddress::NpBase {
                segment: NpSegment::Np1,
                index: 17,
            })
        );
        assert_eq!(
            ChoiceAddress::parse("np.np2.bases[4]"),
            Some(ChoiceAddress::NpBase {
                segment: NpSegment::Np2,
                index: 4,
            })
        );
        assert_eq!(ChoiceAddress::parse("np.np2.bases[x]"), None);
        assert_eq!(
            ChoiceAddress::parse("np.np1.length"),
            Some(ChoiceAddress::NpLength(NpSegment::Np1))
        );
    }

    #[test]
    fn typed_choice_addresses_round_trip_persisted_strings() {
        let cases = [
            (ChoiceAddress::SampleAllele(VdjSegment::V), SAMPLE_ALLELE_V),
            (ChoiceAddress::SampleAllele(VdjSegment::D), SAMPLE_ALLELE_D),
            (ChoiceAddress::SampleAllele(VdjSegment::J), SAMPLE_ALLELE_J),
            (
                ChoiceAddress::Trim {
                    segment: VdjSegment::V,
                    end: TrimEnd::Three,
                },
                TRIM_V_3,
            ),
            (
                ChoiceAddress::Trim {
                    segment: VdjSegment::D,
                    end: TrimEnd::Five,
                },
                TRIM_D_5,
            ),
            (
                ChoiceAddress::Trim {
                    segment: VdjSegment::J,
                    end: TrimEnd::Five,
                },
                TRIM_J_5,
            ),
            (ChoiceAddress::NpLength(NpSegment::Np1), NP1_LENGTH),
            (
                ChoiceAddress::NpBase {
                    segment: NpSegment::Np2,
                    index: 17,
                },
                "np.np2.bases[17]",
            ),
            (ChoiceAddress::MutateUniformCount, MUTATE_UNIFORM_COUNT),
            (
                ChoiceAddress::MutateUniformSite(2),
                "mutate.uniform.site[2]",
            ),
            (
                ChoiceAddress::MutateUniformBase(3),
                "mutate.uniform.base[3]",
            ),
            (ChoiceAddress::MutateS5fCount, MUTATE_S5F_COUNT),
            (ChoiceAddress::MutateS5fSite(4), "mutate.s5f.site[4]"),
            (ChoiceAddress::MutateS5fBase(5), "mutate.s5f.base[5]"),
            (ChoiceAddress::CorruptPcrCount, CORRUPT_PCR_COUNT),
            (
                ChoiceAddress::CorruptPcrSite(6),
                "corrupt.pcr.error_site[6]",
            ),
            (
                ChoiceAddress::CorruptPcrBase(7),
                "corrupt.pcr.error_base[7]",
            ),
            (ChoiceAddress::CorruptQualityCount, CORRUPT_QUALITY_COUNT),
            (
                ChoiceAddress::CorruptQualitySite(8),
                "corrupt.quality.error_site[8]",
            ),
            (
                ChoiceAddress::CorruptQualityBase(9),
                "corrupt.quality.error_base[9]",
            ),
            (
                ChoiceAddress::CorruptContaminantApplied,
                CORRUPT_CONTAMINANT_APPLIED,
            ),
            (
                ChoiceAddress::CorruptContaminantBase(10),
                "corrupt.contaminant.bases[10]",
            ),
            (ChoiceAddress::CorruptIndelCount, CORRUPT_INDEL_COUNT),
            (
                ChoiceAddress::CorruptIndelKind(11),
                "corrupt.indel.kind[11]",
            ),
            (
                ChoiceAddress::CorruptIndelSite(12),
                "corrupt.indel.site[12]",
            ),
            (
                ChoiceAddress::CorruptIndelBase(13),
                "corrupt.indel.base[13]",
            ),
            (ChoiceAddress::CorruptNsCount, CORRUPT_NS_COUNT),
            (ChoiceAddress::CorruptNsSite(14), "corrupt.ns.site[14]"),
            (
                ChoiceAddress::CorruptEndLoss(PrimeEnd::Five),
                CORRUPT_END_LOSS_5,
            ),
            (
                ChoiceAddress::CorruptEndLoss(PrimeEnd::Three),
                CORRUPT_END_LOSS_3,
            ),
            (
                ChoiceAddress::CorruptRevCompApplied,
                CORRUPT_REV_COMP_APPLIED,
            ),
            (
                ChoiceAddress::SampleAlleleDInverted,
                SAMPLE_ALLELE_D_INVERTED,
            ),
            (
                ChoiceAddress::ReceptorRevisionApplied,
                RECEPTOR_REVISION_APPLIED,
            ),
            (
                ChoiceAddress::ReceptorRevisionVAllele,
                RECEPTOR_REVISION_V_ALLELE,
            ),
            (
                ChoiceAddress::ReceptorRevisionVTrim3,
                RECEPTOR_REVISION_V_TRIM_3,
            ),
            (ChoiceAddress::PairedEndR1Length, PAIRED_END_R1_LENGTH),
            (ChoiceAddress::PairedEndR2Length, PAIRED_END_R2_LENGTH),
            (
                ChoiceAddress::PairedEndInsertSize,
                PAIRED_END_INSERT_SIZE,
            ),
        ];

        for (typed, raw) in cases {
            assert_eq!(typed.to_string(), raw);
            assert_eq!(ChoiceAddress::parse(raw), Some(typed));
            assert_eq!(raw.parse::<ChoiceAddress>().unwrap(), typed);
        }
    }

    #[test]
    fn typed_choice_address_rejects_patterns_and_unknown_strings() {
        for raw in [
            "mutate.s5f.site[0..n]",
            "np.np1.bases[x]",
            "assemble.v",
            "sample_allele.np1",
            "trim.np1_3",
            "custom.choice",
            "",
        ] {
            assert_eq!(ChoiceAddress::parse(raw), None, "raw={raw:?}");
            assert!(raw.parse::<ChoiceAddress>().is_err(), "raw={raw:?}");
        }
    }

    #[test]
    fn typed_segment_conversions_reject_wrong_segment_family() {
        assert_eq!(VdjSegment::try_from(Segment::V), Ok(VdjSegment::V));
        assert_eq!(VdjSegment::try_from(Segment::Np1), Err(()));
        assert_eq!(NpSegment::try_from(Segment::Np2), Ok(NpSegment::Np2));
        assert_eq!(NpSegment::try_from(Segment::J), Err(()));
    }

    #[test]
    fn typed_choice_address_patterns_round_trip_report_strings() {
        let cases = [
            (
                ChoiceAddressPattern::SampleAllele(VdjSegment::V),
                SAMPLE_ALLELE_V,
            ),
            (
                ChoiceAddressPattern::SampleAllele(VdjSegment::D),
                SAMPLE_ALLELE_D,
            ),
            (
                ChoiceAddressPattern::SampleAllele(VdjSegment::J),
                SAMPLE_ALLELE_J,
            ),
            (
                ChoiceAddressPattern::Trim {
                    segment: VdjSegment::V,
                    end: TrimEnd::Five,
                },
                TRIM_V_5,
            ),
            (
                ChoiceAddressPattern::Trim {
                    segment: VdjSegment::D,
                    end: TrimEnd::Three,
                },
                TRIM_D_3,
            ),
            (
                ChoiceAddressPattern::Trim {
                    segment: VdjSegment::J,
                    end: TrimEnd::Five,
                },
                TRIM_J_5,
            ),
            (ChoiceAddressPattern::NpLength(NpSegment::Np1), NP1_LENGTH),
            (ChoiceAddressPattern::NpLength(NpSegment::Np2), NP2_LENGTH),
            (
                ChoiceAddressPattern::NpBase(NpSegment::Np1),
                NP1_BASES_PATTERN,
            ),
            (
                ChoiceAddressPattern::NpBase(NpSegment::Np2),
                NP2_BASES_PATTERN,
            ),
            (
                ChoiceAddressPattern::MutateUniformCount,
                MUTATE_UNIFORM_COUNT,
            ),
            (
                ChoiceAddressPattern::MutateUniformSite,
                MUTATE_UNIFORM_SITE_PATTERN,
            ),
            (
                ChoiceAddressPattern::MutateUniformBase,
                MUTATE_UNIFORM_BASE_PATTERN,
            ),
            (ChoiceAddressPattern::MutateS5fCount, MUTATE_S5F_COUNT),
            (ChoiceAddressPattern::MutateS5fSite, MUTATE_S5F_SITE_PATTERN),
            (ChoiceAddressPattern::MutateS5fBase, MUTATE_S5F_BASE_PATTERN),
            (ChoiceAddressPattern::CorruptPcrCount, CORRUPT_PCR_COUNT),
            (
                ChoiceAddressPattern::CorruptPcrSite,
                CORRUPT_PCR_SITE_PATTERN,
            ),
            (
                ChoiceAddressPattern::CorruptPcrBase,
                CORRUPT_PCR_BASE_PATTERN,
            ),
            (
                ChoiceAddressPattern::CorruptQualityCount,
                CORRUPT_QUALITY_COUNT,
            ),
            (
                ChoiceAddressPattern::CorruptQualitySite,
                CORRUPT_QUALITY_SITE_PATTERN,
            ),
            (
                ChoiceAddressPattern::CorruptQualityBase,
                CORRUPT_QUALITY_BASE_PATTERN,
            ),
            (
                ChoiceAddressPattern::CorruptContaminantApplied,
                CORRUPT_CONTAMINANT_APPLIED,
            ),
            (
                ChoiceAddressPattern::CorruptContaminantBase,
                CORRUPT_CONTAMINANT_BASES_PATTERN,
            ),
            (ChoiceAddressPattern::CorruptIndelCount, CORRUPT_INDEL_COUNT),
            (
                ChoiceAddressPattern::CorruptIndelKind,
                CORRUPT_INDEL_KIND_PATTERN,
            ),
            (
                ChoiceAddressPattern::CorruptIndelSite,
                CORRUPT_INDEL_SITE_PATTERN,
            ),
            (
                ChoiceAddressPattern::CorruptIndelBase,
                CORRUPT_INDEL_BASE_PATTERN,
            ),
            (ChoiceAddressPattern::CorruptNsCount, CORRUPT_NS_COUNT),
            (ChoiceAddressPattern::CorruptNsSite, CORRUPT_NS_SITE_PATTERN),
            (
                ChoiceAddressPattern::CorruptEndLoss(PrimeEnd::Five),
                CORRUPT_END_LOSS_5,
            ),
            (
                ChoiceAddressPattern::CorruptEndLoss(PrimeEnd::Three),
                CORRUPT_END_LOSS_3,
            ),
            (
                ChoiceAddressPattern::CorruptRevCompApplied,
                CORRUPT_REV_COMP_APPLIED,
            ),
            (
                ChoiceAddressPattern::SampleAlleleDInverted,
                SAMPLE_ALLELE_D_INVERTED,
            ),
            (
                ChoiceAddressPattern::ReceptorRevisionApplied,
                RECEPTOR_REVISION_APPLIED,
            ),
            (
                ChoiceAddressPattern::ReceptorRevisionVAllele,
                RECEPTOR_REVISION_V_ALLELE,
            ),
            (
                ChoiceAddressPattern::ReceptorRevisionVTrim3,
                RECEPTOR_REVISION_V_TRIM_3,
            ),
            (
                ChoiceAddressPattern::PairedEndR1Length,
                PAIRED_END_R1_LENGTH,
            ),
            (
                ChoiceAddressPattern::PairedEndR2Length,
                PAIRED_END_R2_LENGTH,
            ),
            (
                ChoiceAddressPattern::PairedEndInsertSize,
                PAIRED_END_INSERT_SIZE,
            ),
        ];

        for (typed, raw) in cases {
            assert_eq!(typed.to_string(), raw);
            assert_eq!(ChoiceAddressPattern::parse(raw), Some(typed));
            assert_eq!(raw.parse::<ChoiceAddressPattern>().unwrap(), typed);
        }
    }

    #[test]
    fn typed_choice_address_pattern_rejects_concrete_indexes_and_unknown_strings() {
        for raw in [
            "mutate.s5f.site[2]",
            "np.np1.bases[0]",
            "corrupt.indel.kind[7]",
            "mutate.s5f.site[x]",
            "assemble.v",
            "custom.choice",
            "",
        ] {
            assert_eq!(ChoiceAddressPattern::parse(raw), None, "raw={raw:?}");
            assert!(raw.parse::<ChoiceAddressPattern>().is_err(), "raw={raw:?}");
        }
    }

    // ── Frozen address vocabulary (compile-fence) ─────────────────
    //
    // Every persisted trace file carries an `address_schema_version`.
    // Bumping the constant is the explicit signal that the on-disk
    // vocabulary has changed; this test makes accidental drift loud.
    //
    // For each `ChoiceAddress` variant we pin:
    //   (a) the exact `Display` string (the on-disk spelling), and
    //   (b) that `Display → parse → Display` round-trips to the same
    //       string (the bidirectional contract).
    //
    // A code change that touches the Display or parse paths must
    // either preserve the pinned strings or bump
    // `ADDRESS_SCHEMA_VERSION` and re-pin the new spellings here.

    fn assert_pinned(addr: ChoiceAddress, expected: &str) {
        let s = addr.to_string();
        assert_eq!(
            s, expected,
            "ChoiceAddress::Display drift detected; if intentional, bump ADDRESS_SCHEMA_VERSION",
        );
        let parsed = ChoiceAddress::parse(expected)
            .unwrap_or_else(|| panic!("frozen address string fails to parse: {expected:?}"));
        assert_eq!(
            parsed, addr,
            "ChoiceAddress::parse drift; round-trip would break replay of committed traces",
        );
        let round = parsed.to_string();
        assert_eq!(round, expected, "Display ∘ parse must equal Display");
    }

    #[test]
    fn frozen_address_spellings_for_choice_address_schema_v1() {
        assert_eq!(
            ADDRESS_SCHEMA_VERSION, 1,
            "if you bumped ADDRESS_SCHEMA_VERSION, also re-pin this test",
        );

        // SampleAllele.{v,d,j}
        assert_pinned(
            ChoiceAddress::SampleAllele(VdjSegment::V),
            "sample_allele.v",
        );
        assert_pinned(
            ChoiceAddress::SampleAllele(VdjSegment::D),
            "sample_allele.d",
        );
        assert_pinned(
            ChoiceAddress::SampleAllele(VdjSegment::J),
            "sample_allele.j",
        );

        // Trim.{segment}_{end}
        for (seg, seg_str) in [
            (VdjSegment::V, "v"),
            (VdjSegment::D, "d"),
            (VdjSegment::J, "j"),
        ] {
            for (end, end_str) in [(TrimEnd::Five, "5"), (TrimEnd::Three, "3")] {
                let expected = format!("trim.{seg_str}_{end_str}");
                assert_pinned(ChoiceAddress::Trim { segment: seg, end }, &expected);
            }
        }

        // NpLength.{np1,np2}
        assert_pinned(
            ChoiceAddress::NpLength(NpSegment::Np1),
            "np.np1.length",
        );
        assert_pinned(
            ChoiceAddress::NpLength(NpSegment::Np2),
            "np.np2.length",
        );

        // NpBase[i] — pin a few representative indices.
        for i in [0u32, 1, 5, 42] {
            assert_pinned(
                ChoiceAddress::NpBase {
                    segment: NpSegment::Np1,
                    index: i,
                },
                &format!("np.np1.bases[{i}]"),
            );
            assert_pinned(
                ChoiceAddress::NpBase {
                    segment: NpSegment::Np2,
                    index: i,
                },
                &format!("np.np2.bases[{i}]"),
            );
        }

        // Mutate kernels.
        assert_pinned(ChoiceAddress::MutateUniformCount, "mutate.uniform.count");
        for i in [0u32, 7, 99] {
            assert_pinned(
                ChoiceAddress::MutateUniformSite(i),
                &format!("mutate.uniform.site[{i}]"),
            );
            assert_pinned(
                ChoiceAddress::MutateUniformBase(i),
                &format!("mutate.uniform.base[{i}]"),
            );
        }
        assert_pinned(ChoiceAddress::MutateS5fCount, "mutate.s5f.count");
        for i in [0u32, 3] {
            assert_pinned(
                ChoiceAddress::MutateS5fSite(i),
                &format!("mutate.s5f.site[{i}]"),
            );
            assert_pinned(
                ChoiceAddress::MutateS5fBase(i),
                &format!("mutate.s5f.base[{i}]"),
            );
        }

        // Corruption passes — PCR, quality, contaminant, ns, end_loss,
        // rev_comp, indel.
        assert_pinned(ChoiceAddress::CorruptPcrCount, "corrupt.pcr.count");
        for i in [0u32, 2] {
            assert_pinned(
                ChoiceAddress::CorruptPcrSite(i),
                &format!("corrupt.pcr.error_site[{i}]"),
            );
            assert_pinned(
                ChoiceAddress::CorruptPcrBase(i),
                &format!("corrupt.pcr.error_base[{i}]"),
            );
        }
        assert_pinned(ChoiceAddress::CorruptQualityCount, "corrupt.quality.count");
        for i in [0u32, 2] {
            assert_pinned(
                ChoiceAddress::CorruptQualitySite(i),
                &format!("corrupt.quality.error_site[{i}]"),
            );
            assert_pinned(
                ChoiceAddress::CorruptQualityBase(i),
                &format!("corrupt.quality.error_base[{i}]"),
            );
        }
        assert_pinned(
            ChoiceAddress::CorruptContaminantApplied,
            "corrupt.contaminant.applied",
        );
        for i in [0u32, 5] {
            assert_pinned(
                ChoiceAddress::CorruptContaminantBase(i),
                &format!("corrupt.contaminant.bases[{i}]"),
            );
        }
        assert_pinned(ChoiceAddress::CorruptIndelCount, "corrupt.indel.count");
        for i in [0u32, 4] {
            assert_pinned(
                ChoiceAddress::CorruptIndelKind(i),
                &format!("corrupt.indel.kind[{i}]"),
            );
            assert_pinned(
                ChoiceAddress::CorruptIndelSite(i),
                &format!("corrupt.indel.site[{i}]"),
            );
            assert_pinned(
                ChoiceAddress::CorruptIndelBase(i),
                &format!("corrupt.indel.base[{i}]"),
            );
        }
        assert_pinned(ChoiceAddress::CorruptNsCount, "corrupt.ns.count");
        for i in [0u32, 3] {
            assert_pinned(
                ChoiceAddress::CorruptNsSite(i),
                &format!("corrupt.ns.site[{i}]"),
            );
        }
        assert_pinned(
            ChoiceAddress::CorruptEndLoss(PrimeEnd::Five),
            "corrupt.end_loss.5",
        );
        assert_pinned(
            ChoiceAddress::CorruptEndLoss(PrimeEnd::Three),
            "corrupt.end_loss.3",
        );
        assert_pinned(
            ChoiceAddress::CorruptRevCompApplied,
            "corrupt.rev_comp.applied",
        );

        // D inversion (Slice C). The on-disk spelling sits under the
        // existing `sample_allele.d` namespace so the persisted
        // address-schema-version stays at v1 — no migration burden
        // for traces emitted by pre-Slice-C engines that simply
        // never wrote this address.
        assert_pinned(
            ChoiceAddress::SampleAlleleDInverted,
            "sample_allele.d.inverted",
        );

        // Receptor revision (Slice C of the receptor-revision roadmap).
        // New top-level `receptor_revision.*` namespace; additive
        // under the v1 vocabulary policy — old traces don't reference
        // these strings, new traces parse on engines that postdate
        // this slice. No ADDRESS_SCHEMA_VERSION bump required.
        assert_pinned(
            ChoiceAddress::ReceptorRevisionApplied,
            "receptor_revision.applied",
        );
        assert_pinned(
            ChoiceAddress::ReceptorRevisionVAllele,
            "receptor_revision.v_allele",
        );
        assert_pinned(
            ChoiceAddress::ReceptorRevisionVTrim3,
            "receptor_revision.v_trim_3",
        );

        // Paired-end / read layout (Slice C of the paired-end
        // roadmap). New top-level `paired_end.*` namespace; same
        // additive policy as receptor revision — old traces don't
        // reference these strings, no ADDRESS_SCHEMA_VERSION bump.
        assert_pinned(
            ChoiceAddress::PairedEndR1Length,
            "paired_end.r1_length",
        );
        assert_pinned(
            ChoiceAddress::PairedEndR2Length,
            "paired_end.r2_length",
        );
        assert_pinned(
            ChoiceAddress::PairedEndInsertSize,
            "paired_end.insert_size",
        );

        // P-nucleotide length per end (palindromic addition
        // slice). New top-level `p.*` namespace; same additive
        // policy as receptor revision / paired-end — old traces
        // don't reference these strings, no
        // ADDRESS_SCHEMA_VERSION bump.
        assert_pinned(ChoiceAddress::PLength { end: PEnd::V3 }, "p.v_3.length");
        assert_pinned(ChoiceAddress::PLength { end: PEnd::D5 }, "p.d_5.length");
        assert_pinned(ChoiceAddress::PLength { end: PEnd::D3 }, "p.d_3.length");
        assert_pinned(ChoiceAddress::PLength { end: PEnd::J5 }, "p.j_5.length");

        // Phased genotype (genotype-modeling PR1). New top-level
        // `sample_haplotype` + `sample_gene.*` + `sample_allele_in_slot.*`
        // namespaces; same additive policy as receptor revision /
        // paired-end — old traces don't reference these strings, no
        // ADDRESS_SCHEMA_VERSION bump.
        assert_pinned(ChoiceAddress::SampleHaplotype, "sample_haplotype");
        assert_pinned(ChoiceAddress::SampleGene(VdjSegment::V), "sample_gene.v");
        assert_pinned(ChoiceAddress::SampleGene(VdjSegment::D), "sample_gene.d");
        assert_pinned(ChoiceAddress::SampleGene(VdjSegment::J), "sample_gene.j");
        assert_pinned(
            ChoiceAddress::SampleAlleleInSlot(VdjSegment::V),
            "sample_allele_in_slot.v",
        );
        assert_pinned(
            ChoiceAddress::SampleAlleleInSlot(VdjSegment::D),
            "sample_allele_in_slot.d",
        );
        assert_pinned(
            ChoiceAddress::SampleAlleleInSlot(VdjSegment::J),
            "sample_allele_in_slot.j",
        );
    }
