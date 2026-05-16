//! Per-segment live-call walker.
//!
//! `call_from_region` is the score-and-tie evidence walker that runs
//! over one assembled structural region plus its optional NP-region
//! extensions and produces a `SegmentLiveCall`. The scoring rules,
//! elastic-boundary semantics, and the trim/structural caps live here
//! so the rest of `live_call/` stays focused on the data types.

use super::{
    AlleleBitSet, EvidenceScore, HypothesisFlags, PlacementHypothesis, SegmentLiveCall,
    SegmentRefIndex,
};
use crate::ir::{NucHandle, Nucleotide, Region, Segment, Simulation};
use crate::refdata::AlleleId;

pub(super) fn call_from_region(
    sim: &Simulation,
    segment_index: &SegmentRefIndex,
    region: &Region,
    left_extension: Option<&Region>,
    right_extension: Option<&Region>,
    evidence_version: u64,
) -> SegmentLiveCall {
    let segment = segment_index.segment;
    let allele_universe_len = segment_index.allele_count();
    if region.is_empty() {
        return SegmentLiveCall::unresolved(segment, allele_universe_len);
    }

    // Per-allele score: how many positions (structural + extensions) have a
    // current base that matches this allele's germline byte at the projected
    // reference position. The candidate `allele_call` is the strict tie-set
    // of alleles at the maximum score after all three walks complete.
    //
    // Under this scheme `v_call` is the set of alleles the current evidence
    // cannot rule out: alleles with strictly more matches win, alleles tied
    // at the maximum share the call. Mutations, PCR errors, quality drift,
    // and N corruption all contribute to the scores whichever way they
    // happen to point — the truth allele has no special standing.
    let mut scores: Vec<u32> = vec![0; allele_universe_len];
    let mut informative_matches = 0u32;
    let mut wildcard_matches = 0u32;
    let mut ref_start: Option<u32> = None;
    let mut next_ref_pos: Option<u32> = None;

    for seq_pos in region.start.index()..region.end.index() {
        let nucleotide = sim
            .pool
            .get(NucHandle::new(seq_pos))
            .expect("region range must point into the nucleotide pool");

        // an indel-inserted nucleotide ends up inside V/D/J's
        // region with `germline_pos == NO_GERMLINE_POS`. It carries
        // no allele evidence (no germline byte to compare), so we
        // skip it without failing the call.
        if nucleotide.germline_pos == Nucleotide::NO_GERMLINE_POS {
            // Sanity: only same-segment synthetic bases (e.g. indel
            // insertions placed inside V's region) belong here.
            // Cross-segment synthetic bases would be a malformed IR
            // and still warrant `unsupported_call`.
            if nucleotide.segment != segment {
                return unsupported_call(segment, allele_universe_len, evidence_version);
            }
            continue;
        }

        let Some(ref_pos) = live_ref_pos(nucleotide, segment) else {
            return unsupported_call(segment, allele_universe_len, evidence_version);
        };

        match next_ref_pos {
            // ref_pos may *jump forward* if a base was
            // deleted between this position and the previous one.
            // We allow gap-up but still reject backwards motion,
            // which would indicate a genuinely-broken IR.
            Some(expected) if ref_pos < expected => {
                return unsupported_call(segment, allele_universe_len, evidence_version);
            }
            Some(_) => {}
            None => ref_start = Some(ref_pos),
        }
        next_ref_pos = Some(ref_pos.saturating_add(1));

        // Score increment: every allele whose germline at this ref_pos
        // matches the observed (current) base picks up +1. If no allele
        // matches at this position, no scores change — we simply advance
        // and keep walking. We never abort the call on a single
        // unmatched position; only an explicitly malformed IR (above)
        // produces `unsupported_call`.
        let Some(evidence) =
            segment_index.compatible_alleles_at(ref_pos as usize, nucleotide.base)
        else {
            continue;
        };
        for id in evidence.allele_ids.iter_ids() {
            scores[id.as_usize()] = scores[id.as_usize()].saturating_add(1);
        }
        if evidence.informative {
            informative_matches = informative_matches.saturating_add(1);
        } else {
            wildcard_matches = wildcard_matches.saturating_add(1);
        }
    }

    let Some(ref_start_initial) = ref_start else {
        // Every position in the region was an indel-insert (no
        // germline_pos). Defer to the unresolved state — there is no
        // structural evidence to anchor a hypothesis.
        return SegmentLiveCall::unresolved(segment, allele_universe_len);
    };
    let mut ref_start = ref_start_initial;
    let mut ref_end = next_ref_pos.expect("non-empty region should set ref_end");
    let mut seq_start = region.start.index();
    let mut seq_end = region.end.index();
    let mut flags = HypothesisFlags::EMPTY;

    // Trim-bounded extension caps. The sampled allele had `trim_5`
    // bases removed from its 5' end and `trim_3` from its 3' end
    // during recombination. The right extension can walk at most
    // `trim_3` ref positions past `ref_end` (no further than the
    // sampled allele's original 3' edge); the left extension can
    // reach at most `trim_5` ref positions before `ref_start`.
    //
    // When the segment has no assignment (test fixtures that build a
    // simulation directly without going through recombination), the
    // cap is unset and the walker uses the pre-existing
    // halt-on-empty-evidence semantics. This keeps unit-level live-call
    // tests backwards-compatible while real pipeline runs always have
    // the trim metadata available.
    let trim_cap_3: Option<u32> = sim
        .assignments
        .get(segment)
        .map(|a| a.trim_3 as u32);
    let trim_cap_5: Option<u32> = sim
        .assignments
        .get(segment)
        .map(|a| a.trim_5 as u32);

    // Structural-boundary caps for the extension walkers.
    //
    // Overlap into the immediately-adjacent structural region is fine
    // (e.g. J's left walker overlapping into D's bases when those bases
    // happen to fit J's projected ref_pos). Crossing *past* that
    // structural region into a more-distant region's territory is not:
    // it produces score noise from coincidental matches in a region the
    // segment cannot plausibly occupy, and breaks the AIRR coordinate
    // invariant that `*_alignment_start..*_alignment_end` is a single
    // contiguous span containing the segment's own bytes only.
    //
    // For a left walker (J/D), the boundary is the structural region
    // immediately before the adjacent NP region; the walker may claim
    // positions inside that structural region, but stops at its left
    // edge. For a right walker (V/D), it is the structural region
    // immediately after the adjacent NP region; the walker may claim
    // inside it but stops at its right edge.
    let left_boundary: Option<u32> = left_extension.and_then(|np_region| {
        sim.sequence
            .regions
            .iter()
            .find(|r| {
                matches!(r.segment, Segment::V | Segment::D | Segment::J)
                    && r.end == np_region.start
            })
            .map(|r| r.start.index())
    });
    let right_boundary: Option<u32> = right_extension.and_then(|np_region| {
        sim.sequence
            .regions
            .iter()
            .find(|r| {
                matches!(r.segment, Segment::V | Segment::D | Segment::J)
                    && r.start == np_region.end
            })
            .map(|r| r.end.index())
    });

    // ── Left-side extension + overlap walk ───────────────────────────
    //
    // If `left_extension` (NP1 for J on a VJ chain, NP2 for J on a
    // VDJ chain — or NP1 for D) sits immediately to the left of
    // `region`, walk backward through it (and beyond it into
    // earlier coding territory). Each position contributes a score
    // increment to every allele whose germline at the projected
    // ref_pos matches the observed base.
    //
    // Mark `BOUNDARY_ELASTIC` whenever any extension position
    // produced an evidence record at all (whether or not it shifted
    // the tie-set). Mark `OVERLAPS_OTHER_SEGMENT` when the walker
    // crosses past `np_region.start` into the next-earlier region's
    // bases.
    //
    // Halt: when ref_pos would go below 0, when the projected ref_pos
    // has no allele covering it (compatible_alleles_at returns None),
    // or when we run out of pool positions; never fail the call.
    if let Some(np_region) = left_extension {
        if np_region.end.index() == seq_start {
            let np_start = np_region.start.index();
            let mut extended_seq_start = seq_start;
            let mut extended_ref_start = ref_start;
            let mut extended_inform = informative_matches;
            let mut extended_wildcard = wildcard_matches;
            let mut extension_added = false;
            let mut crossed_into_overlap = false;
            let mut steps_taken: u32 = 0;

            // Walk pool positions in REVERSE order from `np_region.end`
            // down to the structural-boundary cap (the start of the
            // structural region immediately before this NP region) so
            // the walker can overlap into that prior structural region
            // when evidence supports it, but never cross past it.
            let lower_bound: u32 = left_boundary.unwrap_or(0);
            for seq_pos in (lower_bound..np_region.end.index()).rev() {
                // Cap by the sampled allele's 5' trim — the walker can
                // reach at most `trim_5` ref positions before the
                // structural region's start, since that's the only
                // territory the sampled allele could have originally
                // covered.
                if let Some(cap) = trim_cap_5 {
                    if steps_taken >= cap {
                        break;
                    }
                }
                if extended_ref_start == 0 {
                    break;
                }
                let candidate_ref_pos = extended_ref_start - 1;
                let nucleotide = sim
                    .pool
                    .get(NucHandle::new(seq_pos))
                    .expect("extension range must point into the nucleotide pool");
                let Some(evidence) = segment_index
                    .compatible_alleles_at(candidate_ref_pos as usize, nucleotide.base)
                else {
                    break;
                };
                // No allele has this base at the projected ref_pos. The
                // extension has lost its anchor — halt rather than walk
                // into territory where the only thing we could find is a
                // single random later match polluting the tie-set.
                if evidence.allele_ids.is_empty() {
                    break;
                }
                for id in evidence.allele_ids.iter_ids() {
                    scores[id.as_usize()] = scores[id.as_usize()].saturating_add(1);
                }
                extended_seq_start = seq_pos;
                extended_ref_start = candidate_ref_pos;
                if evidence.informative {
                    extended_inform = extended_inform.saturating_add(1);
                } else {
                    extended_wildcard = extended_wildcard.saturating_add(1);
                }
                if seq_pos < np_start {
                    crossed_into_overlap = true;
                }
                extension_added = true;
                steps_taken = steps_taken.saturating_add(1);
            }

            if extension_added {
                seq_start = extended_seq_start;
                ref_start = extended_ref_start;
                informative_matches = extended_inform;
                wildcard_matches = extended_wildcard;
                flags.insert(HypothesisFlags::BOUNDARY_ELASTIC);
                if crossed_into_overlap {
                    flags.insert(HypothesisFlags::OVERLAPS_OTHER_SEGMENT);
                }
            }
        }
    }

    // ── Right-side extension + overlap walk ──────────────────────────
    //
    // Walk forward from `region.end`. Each position contributes a
    // score increment to every allele whose germline at the projected
    // ref_pos matches the observed base. The walker continues into the
    // next-later coding region's bases (D for V, J for D) until the
    // projected ref_pos no longer has any covering allele, marking
    // `OVERLAPS_OTHER_SEGMENT` when it crosses past `np_region.end`.
    if let Some(np_region) = right_extension {
        if np_region.start.index() == seq_end {
            let np_end = np_region.end.index();
            let pool_len = sim.pool.len() as u32;
            let mut extended_seq_end = seq_end;
            let mut extended_ref_pos = ref_end;
            let mut extended_inform = informative_matches;
            let mut extended_wildcard = wildcard_matches;
            let mut extension_added = false;
            let mut crossed_into_overlap = false;
            let mut steps_taken: u32 = 0;

            // Cap walker forward by the structural region immediately
            // after the adjacent NP region (e.g. D for V's right walker,
            // J for D's right walker). Overlap into that region is fine;
            // crossing past it into more-distant territory is not.
            let upper_bound: u32 = right_boundary.unwrap_or(pool_len).min(pool_len);
            for seq_pos in np_region.start.index()..upper_bound {
                // Cap by the sampled allele's 3' trim — the walker can
                // reach at most `trim_3` ref positions past the
                // structural region's end, since that's the only
                // territory the sampled allele could have originally
                // covered.
                if let Some(cap) = trim_cap_3 {
                    if steps_taken >= cap {
                        break;
                    }
                }
                let nucleotide = sim
                    .pool
                    .get(NucHandle::new(seq_pos))
                    .expect("extension range must point into the nucleotide pool");
                let Some(evidence) = segment_index
                    .compatible_alleles_at(extended_ref_pos as usize, nucleotide.base)
                else {
                    break;
                };
                // No allele has this base at the projected ref_pos. The
                // extension has lost its anchor — halt rather than walk
                // into territory where the only thing we could find is a
                // single random later match polluting the tie-set.
                if evidence.allele_ids.is_empty() {
                    break;
                }
                for id in evidence.allele_ids.iter_ids() {
                    scores[id.as_usize()] = scores[id.as_usize()].saturating_add(1);
                }
                extended_seq_end = seq_pos.saturating_add(1);
                extended_ref_pos = extended_ref_pos.saturating_add(1);
                if evidence.informative {
                    extended_inform = extended_inform.saturating_add(1);
                } else {
                    extended_wildcard = extended_wildcard.saturating_add(1);
                }
                if extended_seq_end > np_end {
                    crossed_into_overlap = true;
                }
                extension_added = true;
                steps_taken = steps_taken.saturating_add(1);
            }

            if extension_added {
                seq_end = extended_seq_end;
                ref_end = extended_ref_pos;
                informative_matches = extended_inform;
                wildcard_matches = extended_wildcard;
                flags.insert(HypothesisFlags::BOUNDARY_ELASTIC);
                if crossed_into_overlap {
                    flags.insert(HypothesisFlags::OVERLAPS_OTHER_SEGMENT);
                }
            }
        }
    }

    if wildcard_matches > 0 {
        flags.insert(HypothesisFlags::HAS_WILDCARD_EVIDENCE);
    }

    // ── Build the tie-set ────────────────────────────────────────────
    //
    // The candidate allele_call is the strict tie-set at max score —
    // alleles with strictly more matches win, alleles tied at the
    // maximum share the call. If no allele matched any position
    // (max_score == 0, e.g. an extreme-corruption simulation where
    // every base has been edited away from every allele's germline),
    // every allele is equally consistent with the absence of
    // evidence — return the full pool.
    let max_score = scores.iter().copied().max().unwrap_or(0);
    let allele_call = if max_score == 0 {
        segment_index.all_alleles.clone()
    } else {
        let mut bitset = AlleleBitSet::empty(allele_universe_len);
        for (idx, &score) in scores.iter().enumerate() {
            if score == max_score {
                bitset.insert(AlleleId::new(idx as u32));
            }
        }
        bitset
    };

    let hypothesis = PlacementHypothesis::new(
        segment,
        seq_start,
        seq_end,
        ref_start,
        ref_end,
        allele_call,
        EvidenceScore::exact(informative_matches, wildcard_matches),
        flags,
    );

    SegmentLiveCall::from_hypotheses(
        segment,
        allele_universe_len,
        vec![hypothesis],
        evidence_version,
    )
}

fn live_ref_pos(nucleotide: &Nucleotide, segment: Segment) -> Option<u32> {
    if nucleotide.segment != segment || nucleotide.germline_pos == Nucleotide::NO_GERMLINE_POS {
        return None;
    }
    Some(nucleotide.germline_pos as u32)
}

fn unsupported_call(
    segment: Segment,
    allele_universe_len: usize,
    evidence_version: u64,
) -> SegmentLiveCall {
    SegmentLiveCall::from_hypotheses(segment, allele_universe_len, Vec::new(), evidence_version)
}
