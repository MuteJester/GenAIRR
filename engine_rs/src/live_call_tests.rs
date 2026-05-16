    use super::*;
    use crate::refdata::{Allele, AllelePool, ChainType, RefDataConfig};

    fn id(index: u32) -> AlleleId {
        AlleleId::new(index)
    }

    fn allele(segment: Segment, name: &str, seq: &[u8]) -> Allele {
        Allele {
            name: name.to_string(),
            gene: name.split('*').next().unwrap_or(name).to_string(),
            seq: seq.to_vec(),
            segment,
            anchor: None,
        }
    }

    fn simulation_with_region(segment: Segment, bases: &[u8], ref_start: u16) -> Simulation {
        let mut sim = Simulation::new();
        for (offset, base) in bases.iter().enumerate() {
            let (next, _) = sim.with_nucleotide_pushed(Nucleotide::germline(
                *base,
                ref_start + offset as u16,
                segment,
            ));
            sim = next;
        }
        sim.with_region_added(Region::new(
            segment,
            NucHandle::new(0),
            NucHandle::new(bases.len() as u32),
        ))
    }

    #[test]
    fn assembled_segment_live_call_keeps_trim_induced_ambiguity() {
        let mut cfg = RefDataConfig::empty(ChainType::Vdj);
        let id0 = cfg.v_pool.push(allele(Segment::V, "v1*01", b"GGAAACCC"));
        let id1 = cfg.v_pool.push(allele(Segment::V, "v2*01", b"TTAAACCC"));
        let _id2 = cfg.v_pool.push(allele(Segment::V, "v3*01", b"GGTATCCC"));
        let index = ReferenceMatchIndex::build(&cfg);
        let sim = simulation_with_region(Segment::V, b"AAACCC", 2);

        let call = assembled_segment_live_call(&sim, &index, Segment::V, 1)
            .expect("V region should produce a call");

        assert_eq!(call.allele_call.to_ids(), vec![id0, id1]);
        assert_eq!(call.confidence, LiveCallConfidence::ExactAmbiguous);
        assert_eq!(call.boundary_summary.seq_start, BoundaryValue::Single(0));
        assert_eq!(call.boundary_summary.seq_end, BoundaryValue::Single(6));
        assert_eq!(call.boundary_summary.ref_start, BoundaryValue::Single(2));
        assert_eq!(call.boundary_summary.ref_end, BoundaryValue::Single(8));
        assert_eq!(call.hypotheses[0].score, EvidenceScore::exact(6, 0));
    }

    #[test]
    fn assembled_segment_live_call_uses_current_observed_bases() {
        let mut cfg = RefDataConfig::empty(ChainType::Vdj);
        let _id0 = cfg.v_pool.push(allele(Segment::V, "v1*01", b"AACT"));
        let id1 = cfg.v_pool.push(allele(Segment::V, "v1*02", b"AAGT"));
        let index = ReferenceMatchIndex::build(&cfg);
        let sim = simulation_with_region(Segment::V, b"AACT", 0);
        let sim = sim.with_base_changed(NucHandle::new(2), b'G');

        let call = assembled_segment_live_call(&sim, &index, Segment::V, 1)
            .expect("V region should produce a call");

        assert_eq!(call.allele_call.to_ids(), vec![id1]);
        assert_eq!(call.confidence, LiveCallConfidence::ExactSingle);
        assert_eq!(call.hypotheses[0].score, EvidenceScore::exact(4, 0));
    }

    // ────────────────────────────────────────────────────────────────
    // Y6 score-and-tie behavior tests.
    //
    // These four tests pin down the score-and-tie caller's contract:
    //   1. Mutation can't silently drop the truth allele if other
    //      positions still favor it.
    //   2. A mutation that flips bases toward a different allele can
    //      legitimately divert the call away from truth — the aligner
    //      drift narrative.
    //   3. Trim-induced ambiguity gets resolved when extension into NP
    //      bases finds evidence for one of the tied alleles.
    //   4. When no allele matches any position (max_score == 0), the
    //      call is the full pool — never empty.
    // ────────────────────────────────────────────────────────────────

    #[test]
    fn y6_truth_allele_retained_under_single_position_mutation() {
        // Pool of three alleles, each distinguishable from v0 at a
        // position OTHER than the one we mutate. Sequence matches v0
        // exactly except for position 2 mutated to a base no allele has
        // at ref_pos 2 — this position becomes non-informative. The
        // other un-mutated positions still single out v0 as the unique
        // best match.
        let mut cfg = RefDataConfig::empty(ChainType::Vdj);
        // v0 = ACGTAC (truth)
        // v1 differs from v0 at pos 5 only (C vs G)
        // v2 differs from v0 at pos 1 only (C vs G)
        let v0 = cfg.v_pool.push(allele(Segment::V, "v*01", b"ACGTAC"));
        let v1 = cfg.v_pool.push(allele(Segment::V, "v*02", b"ACGTAG"));
        let v2 = cfg.v_pool.push(allele(Segment::V, "v*03", b"AGGTAC"));
        let index = ReferenceMatchIndex::build(&cfg);

        let sim = simulation_with_region(Segment::V, b"ACGTAC", 0);
        // Mutate position 2 (G→T): no allele has T at ref_pos 2.
        let sim = sim.with_base_changed(NucHandle::new(2), b'T');

        let call = assembled_segment_live_call(&sim, &index, Segment::V, 1)
            .expect("V region should produce a call");

        // Per-position scoring with pool bases ACTTAC after the edit:
        //   pos 0 (A) → v0, v1, v2 all +1
        //   pos 1 (C) → v0, v1 +1 (v2 has G at pos 1)
        //   pos 2 (T) → no allele has T at ref_pos 2 → 0 update
        //   pos 3 (T) → v0, v1, v2 all +1
        //   pos 4 (A) → v0, v1, v2 all +1
        //   pos 5 (C) → v0, v2 +1 (v1 has G at pos 5)
        // Scores: v0=5, v1=4, v2=4. Tie-set = {v0}.
        assert_eq!(call.allele_call.to_ids(), vec![v0]);
        let _unused = (v1, v2);
    }

    #[test]
    fn y6_v_call_diverges_from_truth_when_mutations_flip_toward_another_allele() {
        // Two alleles that disagree on three positions. The pool bases
        // start matching v0 exactly, then three edits flip them toward
        // v1's germline. v1 now scores strictly higher — the call
        // genuinely diverges from the originally-sampled v0.
        let mut cfg = RefDataConfig::empty(ChainType::Vdj);
        let v0 = cfg.v_pool.push(allele(Segment::V, "v*01", b"ACGTAC"));
        let v1 = cfg.v_pool.push(allele(Segment::V, "v*02", b"AGGTGA"));
        let index = ReferenceMatchIndex::build(&cfg);

        // Region built from v0's sequence...
        let sim = simulation_with_region(Segment::V, b"ACGTAC", 0);
        // ...with three positions flipped to v1's bases.
        let sim = sim.with_base_changed(NucHandle::new(1), b'G');
        let sim = sim.with_base_changed(NucHandle::new(4), b'G');
        let sim = sim.with_base_changed(NucHandle::new(5), b'A');

        let call = assembled_segment_live_call(&sim, &index, Segment::V, 1)
            .expect("V region should produce a call");

        // Current pool bases: A G G T G A → matches v1 entirely (6).
        //                                  → matches v0 at pos {0,2,3} (3).
        // Tie-set = {v1}. v_call diverges from truth_v_call (v0).
        assert_eq!(call.allele_call.to_ids(), vec![v1]);
        let _unused = v0;
    }

    #[test]
    fn y6_trim_ambiguity_narrows_via_np_extension() {
        // Three alleles share their first 4 bases (AACC) but diverge on
        // bytes 4-5. The structural region holds only the first 4
        // bases (trim ate the last 2) — all three score 4 there. An
        // adjacent NP region holds bases that match v1's continuation
        // at ref_pos 4 and 5. The extension walker scores v1 twice and
        // the tie-set narrows to {v1}.
        let mut cfg = RefDataConfig::empty(ChainType::Vdj);
        let v0 = cfg.v_pool.push(allele(Segment::V, "v*01", b"AACCAA"));
        let v1 = cfg.v_pool.push(allele(Segment::V, "v*02", b"AACCTG"));
        let v2 = cfg.v_pool.push(allele(Segment::V, "v*03", b"AACCGG"));
        let index = ReferenceMatchIndex::build(&cfg);

        // Build V region (first 4 bases), then NP1 region with bases
        // matching v1's bytes at ref_pos 4 and 5 (T, G).
        let mut sim = Simulation::new();
        for (i, &b) in b"AACC".iter().enumerate() {
            let (next, _) =
                sim.with_nucleotide_pushed(Nucleotide::germline(b, i as u16, Segment::V));
            sim = next;
        }
        sim = sim.with_region_added(Region::new(
            Segment::V,
            NucHandle::new(0),
            NucHandle::new(4),
        ));
        // NP region holding T, G. NP bases use GermlinePos::NONE so the
        // structural walker won't try to score them; the extension
        // walker uses them as evidence for v1's continuation.
        for &b in &[b'T', b'G'] {
            let (next, _) = sim.with_nucleotide_pushed(Nucleotide::synthetic(
                b,
                Segment::Np1,
                crate::ir::NucFlags::empty(),
            ));
            sim = next;
        }
        sim = sim.with_region_added(Region::new(
            Segment::Np1,
            NucHandle::new(4),
            NucHandle::new(6),
        ));

        let call = assembled_segment_live_call(&sim, &index, Segment::V, 1)
            .expect("V region should produce a call");

        // Structural: all 3 score 4. Right-extension into NP1:
        //   ref_pos 4 with 'T' → only v1 has T → v1 += 1
        //   ref_pos 5 with 'G' → v1 and v2 both have G → both += 1
        // Final scores: v0=4, v1=6, v2=5. Tie-set = {v1}.
        assert_eq!(call.allele_call.to_ids(), vec![v1]);
        let _unused = (v0, v2);
    }

    #[test]
    fn y6_full_pool_returned_when_no_allele_matches_any_position() {
        // Pool of two alleles. Region bases are mutated to a base
        // (X — non-canonical) that no allele has at any ref_pos. No
        // scores accumulate → max_score == 0 → tie-set is the full pool
        // (every allele is equally consistent with the absence of
        // evidence).
        let mut cfg = RefDataConfig::empty(ChainType::Vdj);
        let v0 = cfg.v_pool.push(allele(Segment::V, "v*01", b"AACT"));
        let v1 = cfg.v_pool.push(allele(Segment::V, "v*02", b"GGCA"));
        let index = ReferenceMatchIndex::build(&cfg);

        let sim = simulation_with_region(Segment::V, b"AACT", 0);
        // Hit every position with a base no allele has at that ref_pos.
        // v0 = AACT, v1 = GGCA — flipping every position to 'X' (any
        // canonical base that mismatches both) gives a zero-score walk.
        let sim = sim.with_base_changed(NucHandle::new(0), b'T'); // v0 has A, v1 has G
        let sim = sim.with_base_changed(NucHandle::new(1), b'T'); // v0 has A, v1 has G
        let sim = sim.with_base_changed(NucHandle::new(2), b'A'); // v0 has C, v1 has C
        let sim = sim.with_base_changed(NucHandle::new(3), b'G'); // v0 has T, v1 has A

        let call = assembled_segment_live_call(&sim, &index, Segment::V, 1)
            .expect("V region should produce a call");

        let mut ids = call.allele_call.to_ids();
        ids.sort_by_key(|i| i.index());
        // Full pool returned — v_call is honest about having no
        // distinguishing evidence rather than empty-string lying.
        assert_eq!(ids, vec![v0, v1]);
    }

    #[test]
    fn assembled_segment_live_call_treats_n_as_non_informative() {
        let mut cfg = RefDataConfig::empty(ChainType::Vdj);
        let id0 = cfg.v_pool.push(allele(Segment::V, "v1*01", b"AACT"));
        let id1 = cfg.v_pool.push(allele(Segment::V, "v1*02", b"AAGT"));
        let index = ReferenceMatchIndex::build(&cfg);
        let sim = simulation_with_region(Segment::V, b"AACT", 0);
        let sim = sim.with_base_changed(NucHandle::new(2), b'N');

        let call = assembled_segment_live_call(&sim, &index, Segment::V, 1)
            .expect("V region should produce a call");

        assert_eq!(call.allele_call.to_ids(), vec![id0, id1]);
        assert_eq!(call.confidence, LiveCallConfidence::ExactAmbiguous);
        assert_eq!(call.hypotheses[0].score, EvidenceScore::exact(3, 1));
        assert!(call.hypotheses[0]
            .flags
            .contains(HypothesisFlags::HAS_WILDCARD_EVIDENCE));
    }

    #[test]
    fn with_assembled_segment_live_call_persists_state_update() {
        let mut cfg = RefDataConfig::empty(ChainType::Vdj);
        let id0 = cfg.v_pool.push(allele(Segment::V, "v1*01", b"AACT"));
        let index = ReferenceMatchIndex::build(&cfg);
        let sim = simulation_with_region(Segment::V, b"AACT", 0);

        let next = with_assembled_segment_live_call(&sim, &index, Segment::V);

        assert!(sim.live_calls.is_none());
        let live = next
            .live_calls
            .expect("live call state should be initialized");
        let v = live.get(Segment::V).expect("V call should be present");
        assert_eq!(live.version, 1);
        assert_eq!(v.evidence_version, 1);
        assert_eq!(v.allele_call.to_ids(), vec![id0]);
    }

    #[test]
    fn segment_ref_index_maps_coordinate_base_evidence() {
        let mut pool = AllelePool::new();
        let id0 = pool.push(allele(Segment::V, "v1*01", b"AACT"));
        let id1 = pool.push(allele(Segment::V, "v1*02", b"AACG"));
        let id2 = pool.push(allele(Segment::V, "v2*01", b"AGCT"));

        let index = SegmentRefIndex::build(Segment::V, &pool, DEFAULT_REFERENCE_KMER_LEN);

        let a_at_pos1 = index.compatible_alleles_at(1, b'A').unwrap();
        assert!(a_at_pos1.informative);
        assert_eq!(a_at_pos1.allele_ids.to_ids(), vec![id0, id1]);

        let g_at_pos1 = index.compatible_alleles_at(1, b'g').unwrap();
        assert!(g_at_pos1.informative);
        assert_eq!(g_at_pos1.allele_ids.to_ids(), vec![id2]);

        let wildcard_at_pos1 = index.compatible_alleles_at(1, b'N').unwrap();
        assert!(!wildcard_at_pos1.informative);
        assert_eq!(wildcard_at_pos1.allele_ids.to_ids(), vec![id0, id1, id2]);

        assert!(index.compatible_alleles_at(1, b'X').is_none());
        assert!(index.compatible_alleles_at(99, b'A').is_none());
    }

    #[test]
    fn segment_ref_index_falls_back_to_short_kmers_for_short_segments() {
        let mut pool = AllelePool::new();
        let id0 = pool.push(allele(Segment::D, "d1*01", b"ATG"));
        let id1 = pool.push(allele(Segment::D, "d1*02", b"ATG"));
        let id2 = pool.push(allele(Segment::D, "d2*01", b"TGA"));

        let index = SegmentRefIndex::build(Segment::D, &pool, DEFAULT_REFERENCE_KMER_LEN);
        assert_eq!(index.kmer_index.kmer_len(), 3);

        let atg_hits = index.kmer_hits(b"atg").unwrap();
        assert_eq!(
            atg_hits,
            &[
                KmerHit {
                    allele_id: id0,
                    ref_pos: 0
                },
                KmerHit {
                    allele_id: id1,
                    ref_pos: 0
                },
            ]
        );

        let tga_hits = index.kmer_hits(b"TGA").unwrap();
        assert_eq!(
            tga_hits,
            &[KmerHit {
                allele_id: id2,
                ref_pos: 0
            }]
        );
        assert!(index.kmer_hits(b"AT").is_none());
        assert!(index.kmer_hits(b"ANN").is_none());
    }

    #[test]
    fn boundary_indexes_merge_indistinguishable_prefixes_and_suffixes() {
        let mut pool = AllelePool::new();
        let id0 = pool.push(allele(Segment::V, "v1*01", b"AAACCCGGG"));
        let id1 = pool.push(allele(Segment::V, "v2*01", b"TTTCCCGGG"));
        let id2 = pool.push(allele(Segment::V, "v1*02", b"AAACCCGGA"));

        let index = SegmentRefIndex::build(Segment::V, &pool, DEFAULT_REFERENCE_KMER_LEN);

        assert_eq!(
            index.prefix_index.exact_matches(b"aaa").unwrap().to_ids(),
            vec![id0, id2]
        );
        assert_eq!(
            index.suffix_index.exact_matches(b"GGG").unwrap().to_ids(),
            vec![id0, id1]
        );
        assert_eq!(
            index.suffix_index.exact_matches(b"GGA").unwrap().to_ids(),
            vec![id2]
        );
        assert!(index.prefix_index.exact_matches(b"AAN").is_none());
    }

    #[test]
    fn reference_match_index_routes_vdj_segments_and_skips_np() {
        let mut cfg = RefDataConfig::empty(ChainType::Vj);
        let _ = cfg.v_pool.push(allele(Segment::V, "v1*01", b"AAACCCGGG"));
        let _ = cfg.j_pool.push(allele(Segment::J, "j1*01", b"TTTAAA"));

        let index = ReferenceMatchIndex::build(&cfg);

        assert_eq!(index.get(Segment::V).unwrap().allele_count(), 1);
        assert_eq!(index.get(Segment::J).unwrap().allele_count(), 1);
        assert_eq!(index.get(Segment::D).unwrap().allele_count(), 0);
        assert!(index.get(Segment::Np1).is_none());
        assert!(index.d.all_alleles.is_empty());
    }

    #[test]
    fn allele_bitset_insert_remove_and_iterate_are_stable() {
        let mut set = AlleleBitSet::empty(130);
        assert!(set.is_empty());

        assert!(set.insert(id(0)));
        assert!(set.insert(id(64)));
        assert!(set.insert(id(129)));
        assert!(!set.insert(id(64)));

        assert_eq!(set.len(), 3);
        assert!(set.contains(id(0)));
        assert!(set.contains(id(64)));
        assert!(set.contains(id(129)));
        assert_eq!(set.to_ids(), vec![id(0), id(64), id(129)]);

        assert!(set.remove(id(64)));
        assert!(!set.remove(id(64)));
        assert_eq!(set.to_ids(), vec![id(0), id(129)]);
    }

    #[test]
    fn allele_bitset_union_and_intersection_require_same_universe() {
        let a = AlleleBitSet::from_ids(8, [id(1), id(2), id(5)]);
        let b = AlleleBitSet::from_ids(8, [id(2), id(3), id(5)]);

        assert_eq!(a.unioned(&b).to_ids(), vec![id(1), id(2), id(3), id(5)]);
        assert_eq!(a.intersected(&b).to_ids(), vec![id(2), id(5)]);
    }

    #[test]
    #[should_panic(expected = "outside universe length")]
    fn allele_bitset_rejects_out_of_universe_id() {
        let mut set = AlleleBitSet::empty(2);
        set.insert(id(2));
    }

    #[test]
    fn full_bitset_masks_unused_tail_bits() {
        let set = AlleleBitSet::full(65);
        assert_eq!(set.len(), 65);
        assert_eq!(set.to_ids().first(), Some(&id(0)));
        assert_eq!(set.to_ids().last(), Some(&id(64)));
    }

    #[test]
    fn segment_live_call_summarizes_hypothesis_boundaries_and_alleles() {
        let h1 = PlacementHypothesis::new(
            Segment::V,
            10,
            30,
            0,
            20,
            AlleleBitSet::from_ids(6, [id(1), id(2)]),
            EvidenceScore::exact(20, 0),
            HypothesisFlags::EMPTY,
        );
        let h2 = PlacementHypothesis::new(
            Segment::V,
            10,
            33,
            0,
            23,
            AlleleBitSet::from_ids(6, [id(2), id(4)]),
            EvidenceScore::exact(22, 1),
            HypothesisFlags::BOUNDARY_ELASTIC,
        );

        let call = SegmentLiveCall::from_hypotheses(Segment::V, 6, vec![h1, h2], 7);

        assert_eq!(call.allele_call.to_ids(), vec![id(1), id(2), id(4)]);
        assert_eq!(call.confidence, LiveCallConfidence::ExactAmbiguous);
        assert_eq!(call.boundary_summary.seq_start, BoundaryValue::Single(10));
        assert_eq!(
            call.boundary_summary.seq_end,
            BoundaryValue::Ambiguous(vec![30, 33])
        );
        assert_eq!(call.boundary_summary.ref_start, BoundaryValue::Single(0));
        assert_eq!(
            call.boundary_summary.ref_end,
            BoundaryValue::Ambiguous(vec![20, 23])
        );
        assert_eq!(call.evidence_version, 7);
    }

    #[test]
    fn live_call_state_updates_are_persistent() {
        let state = LiveCallState::empty();
        let call = SegmentLiveCall::unresolved(Segment::J, 4);
        let with_call = state.with_segment_call(call);
        let with_dirty = with_call.with_dirty_window(DirtyWindow::new(
            5,
            8,
            DirtyReason::StructuralIndel { site: 6, delta: 1 },
        ));

        assert!(state.get(Segment::J).is_none());
        assert_eq!(state.version, 0);
        assert!(with_call.get(Segment::J).is_some());
        assert_eq!(with_call.version, 1);
        assert!(with_call.dirty_windows.is_empty());
        assert_eq!(with_dirty.version, 2);
        assert_eq!(with_dirty.dirty_windows.len(), 1);
    }

    #[test]
    #[should_panic(expected = "only defined for V/D/J")]
    fn placement_hypothesis_rejects_np_segments() {
        let _ = PlacementHypothesis::new(
            Segment::Np1,
            0,
            1,
            0,
            1,
            AlleleBitSet::full(1),
            EvidenceScore::exact(1, 0),
            HypothesisFlags::EMPTY,
        );
    }
