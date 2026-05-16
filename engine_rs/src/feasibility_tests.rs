    //! Safety-net unit tests for the pure helpers in `feasibility`.
    //!
    //! Written ahead of the planned `compiled.rs → src/compile/*` split
    //! so any regression in `anchor_tail` / `anchor_head` / the stop-codon
    //! and frame checks gets caught by `cargo test` rather than by a much
    //! later integration test.

    use super::*;

    // ── anchor_tail ───────────────────────────────────────────────

    #[test]
    fn anchor_tail_no_trim_returns_anchor_onwards() {
        //  seq:     0 1 2 3 4 5 6 7 8 9
        //  anchor =                 7 (Cys codon at indices 7,8,9)
        let seq = b"AAAAACGTGCN";
        let out = anchor_tail(seq, 7, 0, 0).expect("anchor fits, no trim");
        assert_eq!(out, &seq[7..]);
    }

    #[test]
    fn anchor_tail_with_5prime_trim_under_anchor_keeps_anchor_onwards() {
        // trim_5 = 4 still leaves the anchor (at index 7) intact.
        let seq = b"AAAAACGTGCN";
        let out = anchor_tail(seq, 7, 4, 0).expect("5' trim doesn't touch anchor");
        // anchor_tail returns from `anchor` to the retained end, NOT from
        // trim_5 — verifies the function's `Some(&seq[anchor..end])` shape.
        assert_eq!(out, &seq[7..]);
    }

    #[test]
    fn anchor_tail_5prime_trim_past_anchor_is_none() {
        let seq = b"AAAAACGTGCN";
        // trim_5 = 8 chews past anchor=7 → Cys gone → None
        assert!(anchor_tail(seq, 7, 8, 0).is_none());
    }

    #[test]
    fn anchor_tail_3prime_trim_into_anchor_codon_is_none() {
        let seq = b"AAAAACGTGCN"; // 11 bases, anchor at 7 → codon spans 7..10
        // trim_3 = 2 leaves only 9 bases (0..9), anchor codon needs up to 10 → None
        assert!(anchor_tail(seq, 7, 0, 2).is_none());
    }

    #[test]
    fn anchor_tail_trim_3_larger_than_seq_is_none_not_panic() {
        let seq = b"AAAAACGTGCN";
        // trim_3 > seq.len() — guarded by checked_sub, returns None
        assert!(anchor_tail(seq, 7, 0, 999).is_none());
    }

    // ── anchor_head ───────────────────────────────────────────────

    #[test]
    fn anchor_head_no_trim_returns_start_through_anchor_codon_end() {
        // J anchor (Trp/Phe) at index 3 → codon ends at index 6.
        let seq = b"AAATGGNNN";
        let out = anchor_head(seq, 3, 0, 0).expect("anchor fits, no trim");
        assert_eq!(out, &seq[0..6]);
    }

    #[test]
    fn anchor_head_with_5prime_trim_under_anchor_starts_at_trim_5() {
        let seq = b"AAATGGNNN";
        let out = anchor_head(seq, 3, 2, 0).expect("trim_5 stays under anchor");
        assert_eq!(out, &seq[2..6]);
    }

    #[test]
    fn anchor_head_5prime_trim_past_anchor_is_none() {
        let seq = b"AAATGGNNN";
        assert!(anchor_head(seq, 3, 5, 0).is_none());
    }

    #[test]
    fn anchor_head_3prime_trim_into_anchor_codon_is_none() {
        let seq = b"AAATGGNNN"; // 9 bases, anchor codon spans 3..6
        // trim_3 = 4 leaves only 0..5, anchor codon needs up to 6 → None
        assert!(anchor_head(seq, 3, 0, 4).is_none());
    }

    // ── productive_length_is_feasible ─────────────────────────────

    #[test]
    fn productive_length_in_frame_with_no_forced_stops_is_feasible() {
        // v_tail (3) + np_len (6) + j_head (3) = 12 bases → in frame.
        // Codons spanning the junction touch NP wildcards, so no
        // forced stops can be detected.
        let v_tail = b"TGT"; // Cys
        let j_head = b"TGG"; // Trp
        assert!(productive_length_is_feasible(v_tail, 6, j_head));
    }

    #[test]
    fn productive_length_out_of_frame_is_not_feasible() {
        // 3 + 5 + 3 = 11 → not divisible by 3.
        let v_tail = b"TGT";
        let j_head = b"TGG";
        assert!(!productive_length_is_feasible(v_tail, 5, j_head));
    }

    #[test]
    fn productive_length_with_forced_stop_in_known_region_is_not_feasible() {
        // v_tail contains a TAA codon at offset 0 — fully known, no NP
        // ambiguity covers it, so it's a forced stop.
        let v_tail = b"TAA"; // stop codon at offset 0
        let j_head = b"TGG";
        assert!(!productive_length_is_feasible(v_tail, 0, j_head));
    }

    // ── no_known_stop_for_length ──────────────────────────────────

    #[test]
    fn no_known_stop_treats_np_bases_as_wildcards() {
        // The codon at offset 0 is partly v_tail (T) and partly NP wildcards
        // (None, None). Because one of the bases is unknown, the function
        // can't say it's a stop and must conservatively return true.
        let v_tail = b"T";
        let j_head = b"GG"; // not 3 bytes — only the partial codon counts
        // junction layout: T _ _ G G  (5 bases — 1 full codon scanned)
        assert!(no_known_stop_for_length(v_tail, 2, j_head));
    }

    #[test]
    fn no_known_stop_detects_stop_in_fully_known_v_tail() {
        let v_tail = b"TAA"; // stop in the first codon
        let j_head = b"GGG";
        assert!(!no_known_stop_for_length(v_tail, 0, j_head));
    }

    #[test]
    fn no_known_stop_detects_stop_in_fully_known_j_head() {
        // Junction is v_tail(3) + np(0) + j_head(3) = 6 → 2 codons,
        // first known TGT (Cys), second known TAG (stop).
        let v_tail = b"TGT";
        let j_head = b"TAG";
        assert!(!no_known_stop_for_length(v_tail, 0, j_head));
    }

    // ── is_vj_productive_choice ───────────────────────────────────

    #[test]
    fn is_vj_productive_choice_accepts_the_six_known_addresses() {
        for addr in [
            "sample_allele.v",
            "sample_allele.j",
            "trim.v_5",
            "trim.v_3",
            "trim.j_5",
            "trim.j_3",
        ] {
            assert!(is_vj_productive_choice(addr), "expected {addr:?} to count");
        }
    }

    #[test]
    fn is_vj_productive_choice_rejects_unrelated_addresses() {
        for addr in [
            "sample_allele.d",
            "trim.d_5",
            "trim.d_3",
            "np.np1.length",
            "mutate.s5f.site[0]",
            "",
            "sample_allele.v.something",
        ] {
            assert!(!is_vj_productive_choice(addr), "expected {addr:?} to NOT count");
        }
    }

    // ── trim_address / sample_allele_address helpers ──────────────

    #[test]
    fn trim_address_maps_each_supported_segment_end_pair() {
        assert_eq!(trim_address(Segment::V, TrimEnd::Five), "trim.v_5");
        assert_eq!(trim_address(Segment::V, TrimEnd::Three), "trim.v_3");
        assert_eq!(trim_address(Segment::J, TrimEnd::Five), "trim.j_5");
        assert_eq!(trim_address(Segment::J, TrimEnd::Three), "trim.j_3");
    }

    #[test]
    fn trim_address_returns_unsupported_marker_for_other_segments() {
        // D / NP1 / NP2 aren't part of the VJ-productive feasibility surface;
        // the helper returns a stable marker rather than panicking.
        assert_eq!(trim_address(Segment::D, TrimEnd::Five), "trim.<unsupported>");
        assert_eq!(trim_address(Segment::Np1, TrimEnd::Three), "trim.<unsupported>");
    }

    #[test]
    fn sample_allele_address_maps_v_and_j_only() {
        assert_eq!(sample_allele_address(Segment::V), "sample_allele.v");
        assert_eq!(sample_allele_address(Segment::J), "sample_allele.j");
        assert_eq!(
            sample_allele_address(Segment::D),
            "sample_allele.<unsupported>"
        );
    }
