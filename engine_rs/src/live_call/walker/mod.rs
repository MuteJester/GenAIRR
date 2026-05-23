//! Per-segment live-call walker.
//!
//! `call_from_region` is the score-and-tie evidence walker that runs
//! over one assembled structural region plus its optional NP-region
//! extensions and produces a `SegmentLiveCall`. The scoring rules,
//! elastic-boundary semantics, and the trim/structural caps live here
//! so the rest of `live_call/` stays focused on the data types.
//!
//! The NP-region extension walks live in the [`extensions`] submodule
//! so they can be shared with the streaming
//! [`super::walker_observer::WalkerObserverState`]: both code paths
//! reuse the same `ExtensionWalkState` borrow shape and apply
//! byte-for-byte identical scoring/flag updates.

use super::{
    AlleleBitSet, EvidenceScore, HypothesisFlags, PlacementHypothesis, SegmentLiveCall,
    SegmentRefIndex,
};
use crate::ir::{NucHandle, Region, Segment, Simulation};
use crate::refdata::AlleleId;

pub(crate) mod extensions;
use extensions::{walk_left_extension, walk_right_extension, ExtensionWalkState};

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

        // Single segment check covers both branches below: cross-segment
        // nucleotides — synthetic (indel-insert) or otherwise — are a
        // malformed IR for this region and produce `unsupported_call`.
        if nucleotide.segment != segment {
            return unsupported_call(segment, allele_universe_len, evidence_version);
        }

        // an indel-inserted nucleotide ends up inside V/D/J's
        // region with `germline_pos == GermlinePos::NONE`. It carries
        // no allele evidence (no germline byte to compare), so we
        // skip it without failing the call.
        let Some(ref_pos) = nucleotide.germline_pos.get().map(|p| p as u32) else {
            continue;
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
        let Some(evidence) = segment_index.compatible_alleles_at(ref_pos as usize, nucleotide.base)
        else {
            continue;
        };
        evidence.allele_ids.for_each_id(|id| {
            let slot = &mut scores[id.as_usize()];
            *slot = slot.saturating_add(1);
        });
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
    let trim_cap_3: Option<u32> = sim.assignments.get(segment).map(|a| a.trim_3 as u32);
    let trim_cap_5: Option<u32> = sim.assignments.get(segment).map(|a| a.trim_5 as u32);

    // Delegate the left and right NP-region extension walks to the
    // shared `extensions` submodule so the streaming
    // `WalkerObserverState` can reuse the exact same code. The walks
    // mutate `scores`, `informative_matches`, `wildcard_matches`,
    // `ref_start`/`ref_end`, `seq_start`/`seq_end`, and `flags`
    // in-place via `ExtensionWalkState` borrows.
    if let Some(np_region) = left_extension {
        let mut state = ExtensionWalkState {
            scores: &mut scores,
            informative_matches: &mut informative_matches,
            wildcard_matches: &mut wildcard_matches,
            ref_start: &mut ref_start,
            ref_end: &mut ref_end,
            seq_start: &mut seq_start,
            seq_end: &mut seq_end,
            flags: &mut flags,
        };
        walk_left_extension(sim, segment_index, np_region, trim_cap_5, &mut state);
    }
    if let Some(np_region) = right_extension {
        let mut state = ExtensionWalkState {
            scores: &mut scores,
            informative_matches: &mut informative_matches,
            wildcard_matches: &mut wildcard_matches,
            ref_start: &mut ref_start,
            ref_end: &mut ref_end,
            seq_start: &mut seq_start,
            seq_end: &mut seq_end,
            flags: &mut flags,
        };
        walk_right_extension(sim, segment_index, np_region, trim_cap_3, &mut state);
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

fn unsupported_call(
    segment: Segment,
    allele_universe_len: usize,
    evidence_version: u64,
) -> SegmentLiveCall {
    SegmentLiveCall::from_hypotheses(segment, allele_universe_len, Vec::new(), evidence_version)
}
