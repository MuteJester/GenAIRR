//! Live V/D/J call evidence state.
//!
//! Data structures plus small, testable operations that back the
//! dynamic allele-call layer.

use crate::ir::{NucHandle, Nucleotide, Region, Segment, Simulation};
use crate::refdata::{Allele, AlleleId, AllelePool, RefDataConfig};
use std::collections::HashMap;

/// Default seed length for compiled reference lookups.
///
/// Segment indexes automatically fall back to a shorter value when all
/// alleles in a segment are shorter than this, which matters for D
/// alleles and tiny synthetic test fixtures.
pub const DEFAULT_REFERENCE_KMER_LEN: usize = 7;

/// Dense allele-id bitset used by live call hypotheses.
///
/// The bitset has a fixed universe size so equality is unambiguous:
/// two bitsets over different allele pools are not equal even if their
/// word storage currently contains the same bits.
#[derive(Clone, Debug, Eq, PartialEq)]
pub struct AlleleBitSet {
    universe_len: usize,
    words: Vec<u64>,
}

impl AlleleBitSet {
    pub fn empty(universe_len: usize) -> Self {
        Self {
            universe_len,
            words: vec![0; universe_len.div_ceil(64)],
        }
    }

    pub fn full(universe_len: usize) -> Self {
        let mut set = Self {
            universe_len,
            words: vec![u64::MAX; universe_len.div_ceil(64)],
        };
        set.clear_unused_tail_bits();
        set
    }

    pub fn from_ids(universe_len: usize, ids: impl IntoIterator<Item = AlleleId>) -> Self {
        let mut set = Self::empty(universe_len);
        for id in ids {
            set.insert(id);
        }
        set
    }

    pub fn universe_len(&self) -> usize {
        self.universe_len
    }

    pub fn is_empty(&self) -> bool {
        self.words.iter().all(|word| *word == 0)
    }

    pub fn len(&self) -> usize {
        self.words
            .iter()
            .map(|word| word.count_ones() as usize)
            .sum()
    }

    pub fn contains(&self, id: AlleleId) -> bool {
        let index = self.checked_index(id);
        (self.words[index / 64] & (1u64 << (index % 64))) != 0
    }

    pub fn insert(&mut self, id: AlleleId) -> bool {
        let index = self.checked_index(id);
        let word = &mut self.words[index / 64];
        let mask = 1u64 << (index % 64);
        let was_present = (*word & mask) != 0;
        *word |= mask;
        !was_present
    }

    pub fn remove(&mut self, id: AlleleId) -> bool {
        let index = self.checked_index(id);
        let word = &mut self.words[index / 64];
        let mask = 1u64 << (index % 64);
        let was_present = (*word & mask) != 0;
        *word &= !mask;
        was_present
    }

    pub fn union_with(&mut self, other: &Self) {
        self.assert_same_universe(other);
        for (lhs, rhs) in self.words.iter_mut().zip(&other.words) {
            *lhs |= *rhs;
        }
    }

    pub fn intersect_with(&mut self, other: &Self) {
        self.assert_same_universe(other);
        for (lhs, rhs) in self.words.iter_mut().zip(&other.words) {
            *lhs &= *rhs;
        }
    }

    pub fn unioned(&self, other: &Self) -> Self {
        let mut out = self.clone();
        out.union_with(other);
        out
    }

    pub fn intersected(&self, other: &Self) -> Self {
        let mut out = self.clone();
        out.intersect_with(other);
        out
    }

    pub fn iter_ids(&self) -> impl Iterator<Item = AlleleId> + '_ {
        let universe_len = self.universe_len;
        self.words
            .iter()
            .copied()
            .enumerate()
            .flat_map(move |(word_index, word)| {
                (0..64).filter_map(move |bit| {
                    let index = word_index * 64 + bit;
                    if index >= universe_len || (word & (1u64 << bit)) == 0 {
                        None
                    } else {
                        Some(AlleleId::new(index as u32))
                    }
                })
            })
    }

    pub fn to_ids(&self) -> Vec<AlleleId> {
        self.iter_ids().collect()
    }

    fn checked_index(&self, id: AlleleId) -> usize {
        let index = id.as_usize();
        assert!(
            index < self.universe_len,
            "AlleleBitSet: allele id {} outside universe length {}",
            id.index(),
            self.universe_len
        );
        index
    }

    fn assert_same_universe(&self, other: &Self) {
        assert_eq!(
            self.universe_len, other.universe_len,
            "AlleleBitSet operation requires matching universe lengths"
        );
    }

    fn clear_unused_tail_bits(&mut self) {
        let unused = self.words.len() * 64 - self.universe_len;
        if unused == 0 || self.words.is_empty() {
            return;
        }
        let used = 64 - unused;
        let mask = if used == 64 {
            u64::MAX
        } else {
            (1u64 << used) - 1
        };
        if let Some(last) = self.words.last_mut() {
            *last &= mask;
        }
    }
}

/// Immutable V/D/J reference index used by future live-call updates.
///
/// This is built once at simulator compilation from the active
/// `RefDataConfig`. It intentionally excludes C/NP data because live
/// allele calls are defined only for V/D/J.
#[derive(Clone, Debug, Eq, PartialEq)]
pub struct ReferenceMatchIndex {
    pub v: SegmentRefIndex,
    pub d: SegmentRefIndex,
    pub j: SegmentRefIndex,
}

impl ReferenceMatchIndex {
    pub fn build(refdata: &RefDataConfig) -> Self {
        Self {
            v: SegmentRefIndex::build(Segment::V, &refdata.v_pool, DEFAULT_REFERENCE_KMER_LEN),
            d: SegmentRefIndex::build(Segment::D, &refdata.d_pool, DEFAULT_REFERENCE_KMER_LEN),
            j: SegmentRefIndex::build(Segment::J, &refdata.j_pool, DEFAULT_REFERENCE_KMER_LEN),
        }
    }

    pub fn get(&self, segment: Segment) -> Option<&SegmentRefIndex> {
        match segment {
            Segment::V => Some(&self.v),
            Segment::D => Some(&self.d),
            Segment::J => Some(&self.j),
            Segment::Np1 | Segment::Np2 => None,
        }
    }
}

/// Indexed copy of one allele's stable reference information.
#[derive(Clone, Debug, Eq, PartialEq)]
pub struct IndexedAllele {
    pub id: AlleleId,
    pub name: String,
    pub gene: String,
    pub seq: Vec<u8>,
    pub len: u32,
    pub anchor: Option<u16>,
}

impl IndexedAllele {
    fn from_allele(id: AlleleId, allele: &Allele) -> Self {
        let seq: Vec<u8> = allele
            .seq
            .iter()
            .map(|base| base.to_ascii_uppercase())
            .collect();
        Self {
            id,
            name: allele.name.clone(),
            gene: allele.gene.clone(),
            len: seq.len() as u32,
            seq,
            anchor: allele.anchor,
        }
    }
}

/// Per-segment compiled reference index.
#[derive(Clone, Debug, Eq, PartialEq)]
pub struct SegmentRefIndex {
    pub segment: Segment,
    pub alleles: Vec<IndexedAllele>,
    pub all_alleles: AlleleBitSet,
    pub base_at_pos: Vec<BaseBitsets>,
    pub kmer_index: KmerIndex,
    pub prefix_index: BoundaryIndex,
    pub suffix_index: BoundaryIndex,
    pub max_len: usize,
}

impl SegmentRefIndex {
    pub fn build(segment: Segment, pool: &AllelePool, requested_kmer_len: usize) -> Self {
        assert_live_segment(segment);
        assert!(
            (1..=31).contains(&requested_kmer_len),
            "reference k-mer length must be in 1..=31"
        );

        let universe_len = pool.len();
        let max_len = pool
            .iter()
            .map(|(_, allele)| allele.seq.len())
            .max()
            .unwrap_or(0);
        let kmer_len = if max_len == 0 {
            requested_kmer_len
        } else {
            requested_kmer_len.min(max_len)
        };

        let mut alleles = Vec::with_capacity(universe_len);
        let mut all_alleles = AlleleBitSet::empty(universe_len);
        let mut base_at_pos = (0..max_len)
            .map(|_| BaseBitsets::empty(universe_len))
            .collect::<Vec<_>>();
        let mut kmer_index = KmerIndex::new(kmer_len);
        let mut prefix_index = BoundaryIndex::new(universe_len);
        let mut suffix_index = BoundaryIndex::new(universe_len);

        for (id, allele) in pool.iter() {
            debug_assert_eq!(allele.segment, segment);
            all_alleles.insert(id);

            for (pos, base) in allele.seq.iter().enumerate() {
                if let Some(base_index) = canonical_base_index(*base) {
                    base_at_pos[pos].insert_index(base_index, id);
                }
            }

            kmer_index.insert_allele(id, &allele.seq);
            prefix_index.insert_prefixes(id, &allele.seq);
            suffix_index.insert_suffixes(id, &allele.seq);
            alleles.push(IndexedAllele::from_allele(id, allele));
        }

        Self {
            segment,
            alleles,
            all_alleles,
            base_at_pos,
            kmer_index,
            prefix_index,
            suffix_index,
            max_len,
        }
    }

    pub fn allele_count(&self) -> usize {
        self.alleles.len()
    }

    pub fn compatible_alleles_at(&self, ref_pos: usize, observed_base: u8) -> Option<BaseEvidence> {
        self.base_at_pos
            .get(ref_pos)
            .and_then(|base_sets| base_sets.compatible_with_observed(observed_base))
    }

    pub fn kmer_hits(&self, kmer: &[u8]) -> Option<&[KmerHit]> {
        self.kmer_index.hits(kmer)
    }
}

/// Return a simulation with the live call for an assembled V/D/J
/// structural region refreshed from exact current-base evidence.
///
/// Discovery is restricted to the existing structural region. It
/// does not extend into NP evidence, repair indel-shifted
/// hypotheses, or alter AIRR projection.
pub fn with_assembled_segment_live_call(
    sim: &Simulation,
    reference_index: &ReferenceMatchIndex,
    segment: Segment,
) -> Simulation {
    assert_live_segment(segment);
    let base_state = sim.live_calls.clone().unwrap_or_default();
    let evidence_version = base_state.version.saturating_add(1);
    let Some(call) = assembled_segment_live_call(sim, reference_index, segment, evidence_version)
    else {
        return sim.clone();
    };
    sim.with_live_calls(base_state.with_segment_call(call))
}

/// Build the exact live call implied by the latest structural region
/// for `segment`.
pub fn assembled_segment_live_call(
    sim: &Simulation,
    reference_index: &ReferenceMatchIndex,
    segment: Segment,
    evidence_version: u64,
) -> Option<SegmentLiveCall> {
    assert_live_segment(segment);
    let segment_index = reference_index.get(segment)?;
    let region = latest_region_for_segment(sim, segment)?;
    let left_extension = left_extension_region_for(sim, segment, region);
    let right_extension = right_extension_region_for(sim, segment, region);
    Some(call_from_region(
        sim,
        segment_index,
        region,
        left_extension,
        right_extension,
        evidence_version,
    ))
}

fn latest_region_for_segment(sim: &Simulation, segment: Segment) -> Option<&Region> {
    sim.sequence
        .regions
        .iter()
        .rev()
        .find(|region| region.segment == segment)
}

/// Pick the NP region (if any) whose bases sit immediately to the right
/// of `segment`'s assembled region.
///
/// - V → NP1.
/// - D → NP2.
/// - J never has a right neighbour in either VJ or VDJ chain layouts.
fn right_extension_region_for<'a>(
    sim: &'a Simulation,
    segment: Segment,
    region: &Region,
) -> Option<&'a Region> {
    match segment {
        Segment::V => sim
            .sequence
            .regions
            .iter()
            .find(|r| r.segment == Segment::Np1 && r.start == region.end),
        Segment::D => sim
            .sequence
            .regions
            .iter()
            .find(|r| r.segment == Segment::Np2 && r.start == region.end),
        Segment::J => None,
        Segment::Np1 | Segment::Np2 => None,
    }
}

/// Pick the NP region (if any) whose bases sit immediately to the left
/// of `segment`'s assembled region.
///
/// - J → adjacent NP region: NP1 in VJ chains (V → NP1 → J), NP2 in
///   VDJ chains (… → NP2 → J). Both cases collapse to "the NP region
///   whose `end` equals `region.start`", so we don't need to know the
///   chain type.
/// - D → NP1.
/// - V never has a left neighbour in our DSL (V is always the first
///   region).
fn left_extension_region_for<'a>(
    sim: &'a Simulation,
    segment: Segment,
    region: &Region,
) -> Option<&'a Region> {
    match segment {
        Segment::J => sim.sequence.regions.iter().find(|r| {
            matches!(r.segment, Segment::Np1 | Segment::Np2) && r.end == region.start
        }),
        Segment::D => sim
            .sequence
            .regions
            .iter()
            .find(|r| r.segment == Segment::Np1 && r.end == region.start),
        Segment::V => None,
        Segment::Np1 | Segment::Np2 => None,
    }
}

fn call_from_region(
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

/// Base-partitioned allele sets for one reference coordinate.
#[derive(Clone, Debug, Eq, PartialEq)]
pub struct BaseBitsets {
    sets: [AlleleBitSet; 4],
}

impl BaseBitsets {
    pub fn empty(universe_len: usize) -> Self {
        Self {
            sets: std::array::from_fn(|_| AlleleBitSet::empty(universe_len)),
        }
    }

    pub fn for_base(&self, base: u8) -> Option<&AlleleBitSet> {
        canonical_base_index(base).map(|index| &self.sets[index])
    }

    pub fn compatible_with_observed(&self, observed_base: u8) -> Option<BaseEvidence> {
        match observed_base_kind(observed_base)? {
            ObservedBaseKind::Canonical(index) => Some(BaseEvidence {
                allele_ids: self.sets[index].clone(),
                informative: true,
            }),
            ObservedBaseKind::Wildcard => {
                let mut allele_ids = AlleleBitSet::empty(self.sets[0].universe_len());
                for set in &self.sets {
                    allele_ids.union_with(set);
                }
                Some(BaseEvidence {
                    allele_ids,
                    informative: false,
                })
            }
        }
    }

    fn insert_index(&mut self, base_index: usize, allele_id: AlleleId) {
        self.sets[base_index].insert(allele_id);
    }
}

/// Result of one observed-base compatibility lookup.
#[derive(Clone, Debug, Eq, PartialEq)]
pub struct BaseEvidence {
    pub allele_ids: AlleleBitSet,
    pub informative: bool,
}

/// One exact k-mer seed hit in a reference allele.
#[derive(Copy, Clone, Debug, Eq, PartialEq, Hash)]
pub struct KmerHit {
    pub allele_id: AlleleId,
    pub ref_pos: u32,
}

/// Exact canonical k-mer index.
#[derive(Clone, Debug, Eq, PartialEq)]
pub struct KmerIndex {
    kmer_len: usize,
    hits_by_kmer: HashMap<u64, Vec<KmerHit>>,
}

impl KmerIndex {
    pub fn new(kmer_len: usize) -> Self {
        assert!(
            (1..=31).contains(&kmer_len),
            "reference k-mer length must be in 1..=31"
        );
        Self {
            kmer_len,
            hits_by_kmer: HashMap::new(),
        }
    }

    pub fn kmer_len(&self) -> usize {
        self.kmer_len
    }

    pub fn len(&self) -> usize {
        self.hits_by_kmer.len()
    }

    pub fn is_empty(&self) -> bool {
        self.hits_by_kmer.is_empty()
    }

    pub fn hits(&self, kmer: &[u8]) -> Option<&[KmerHit]> {
        if kmer.len() != self.kmer_len {
            return None;
        }
        let encoded = encode_kmer(kmer)?;
        self.hits_by_kmer.get(&encoded).map(Vec::as_slice)
    }

    fn insert_allele(&mut self, allele_id: AlleleId, seq: &[u8]) {
        if seq.len() < self.kmer_len {
            return;
        }

        for (pos, window) in seq.windows(self.kmer_len).enumerate() {
            let Some(encoded) = encode_kmer(window) else {
                continue;
            };
            self.hits_by_kmer.entry(encoded).or_default().push(KmerHit {
                allele_id,
                ref_pos: pos as u32,
            });
        }
    }
}

/// Exact prefix/suffix boundary lookup.
#[derive(Clone, Debug, Eq, PartialEq)]
pub struct BoundaryIndex {
    universe_len: usize,
    max_len: usize,
    by_sequence: HashMap<Vec<u8>, AlleleBitSet>,
}

impl BoundaryIndex {
    pub fn new(universe_len: usize) -> Self {
        Self {
            universe_len,
            max_len: 0,
            by_sequence: HashMap::new(),
        }
    }

    pub fn len(&self) -> usize {
        self.by_sequence.len()
    }

    pub fn is_empty(&self) -> bool {
        self.by_sequence.is_empty()
    }

    pub fn max_len(&self) -> usize {
        self.max_len
    }

    pub fn exact_matches(&self, sequence: &[u8]) -> Option<&AlleleBitSet> {
        let sequence = canonicalized_sequence(sequence)?;
        self.by_sequence.get(sequence.as_slice())
    }

    fn insert_prefixes(&mut self, allele_id: AlleleId, seq: &[u8]) {
        for len in 1..=seq.len() {
            self.insert_sequence(&seq[..len], allele_id);
        }
    }

    fn insert_suffixes(&mut self, allele_id: AlleleId, seq: &[u8]) {
        for len in 1..=seq.len() {
            self.insert_sequence(&seq[seq.len() - len..], allele_id);
        }
    }

    fn insert_sequence(&mut self, sequence: &[u8], allele_id: AlleleId) {
        let Some(sequence) = canonicalized_sequence(sequence) else {
            return;
        };
        self.max_len = self.max_len.max(sequence.len());
        self.by_sequence
            .entry(sequence)
            .or_insert_with(|| AlleleBitSet::empty(self.universe_len))
            .insert(allele_id);
    }
}

#[derive(Copy, Clone, Debug, Eq, PartialEq)]
enum ObservedBaseKind {
    Canonical(usize),
    Wildcard,
}

fn observed_base_kind(base: u8) -> Option<ObservedBaseKind> {
    canonical_base_index(base)
        .map(ObservedBaseKind::Canonical)
        .or_else(|| (base == b'N' || base == b'n').then_some(ObservedBaseKind::Wildcard))
}

fn canonical_base_index(base: u8) -> Option<usize> {
    match base.to_ascii_uppercase() {
        b'A' => Some(0),
        b'C' => Some(1),
        b'G' => Some(2),
        b'T' => Some(3),
        _ => None,
    }
}

fn canonicalized_sequence(sequence: &[u8]) -> Option<Vec<u8>> {
    let mut out = Vec::with_capacity(sequence.len());
    for base in sequence {
        if canonical_base_index(*base).is_none() {
            return None;
        }
        out.push(base.to_ascii_uppercase());
    }
    Some(out)
}

fn encode_kmer(kmer: &[u8]) -> Option<u64> {
    let mut encoded = 0u64;
    for base in kmer {
        encoded <<= 2;
        encoded |= canonical_base_index(*base)? as u64;
    }
    Some(encoded)
}

/// Evidence quality for one placement hypothesis.
#[derive(Copy, Clone, Debug, Default, Eq, PartialEq)]
pub struct EvidenceScore {
    pub conflicts: u32,
    pub informative_matches: u32,
    pub wildcard_matches: u32,
    pub span_len: u32,
}

impl EvidenceScore {
    pub const fn exact(informative_matches: u32, wildcard_matches: u32) -> Self {
        Self {
            conflicts: 0,
            informative_matches,
            wildcard_matches,
            span_len: informative_matches + wildcard_matches,
        }
    }

    pub const fn is_exact(self) -> bool {
        self.conflicts == 0
    }
}

/// Flags describing how a hypothesis was discovered.
#[derive(Copy, Clone, Debug, Default, Eq, PartialEq)]
pub struct HypothesisFlags(u32);

impl HypothesisFlags {
    pub const EMPTY: Self = Self(0);
    pub const BOUNDARY_ELASTIC: Self = Self(1 << 0);
    pub const OVERLAPS_OTHER_SEGMENT: Self = Self(1 << 1);
    pub const HAS_WILDCARD_EVIDENCE: Self = Self(1 << 2);

    pub const fn contains(self, flag: Self) -> bool {
        (self.0 & flag.0) == flag.0
    }

    pub fn insert(&mut self, flag: Self) {
        self.0 |= flag.0;
    }
}

/// One exact-equivalence placement hypothesis for a V/D/J segment.
#[derive(Clone, Debug, Eq, PartialEq)]
pub struct PlacementHypothesis {
    pub segment: Segment,
    pub seq_start: u32,
    pub seq_end: u32,
    pub ref_start: u32,
    pub ref_end: u32,
    pub allele_ids: AlleleBitSet,
    pub score: EvidenceScore,
    pub flags: HypothesisFlags,
}

impl PlacementHypothesis {
    #[allow(clippy::too_many_arguments)]
    pub fn new(
        segment: Segment,
        seq_start: u32,
        seq_end: u32,
        ref_start: u32,
        ref_end: u32,
        allele_ids: AlleleBitSet,
        score: EvidenceScore,
        flags: HypothesisFlags,
    ) -> Self {
        assert_live_segment(segment);
        assert!(seq_start <= seq_end, "hypothesis seq range is inverted");
        assert!(ref_start <= ref_end, "hypothesis ref range is inverted");
        Self {
            segment,
            seq_start,
            seq_end,
            ref_start,
            ref_end,
            allele_ids,
            score,
            flags,
        }
    }
}

#[derive(Clone, Debug, Eq, PartialEq)]
pub enum BoundaryValue {
    Unresolved,
    Single(u32),
    Ambiguous(Vec<u32>),
}

impl BoundaryValue {
    pub fn from_values(values: impl IntoIterator<Item = u32>) -> Self {
        let mut values: Vec<u32> = values.into_iter().collect();
        values.sort_unstable();
        values.dedup();

        match values.len() {
            0 => Self::Unresolved,
            1 => Self::Single(values[0]),
            _ => Self::Ambiguous(values),
        }
    }
}

#[derive(Clone, Debug, Eq, PartialEq)]
pub struct BoundarySummary {
    pub seq_start: BoundaryValue,
    pub seq_end: BoundaryValue,
    pub ref_start: BoundaryValue,
    pub ref_end: BoundaryValue,
}

impl BoundarySummary {
    pub fn unresolved() -> Self {
        Self {
            seq_start: BoundaryValue::Unresolved,
            seq_end: BoundaryValue::Unresolved,
            ref_start: BoundaryValue::Unresolved,
            ref_end: BoundaryValue::Unresolved,
        }
    }

    pub fn from_hypotheses(hypotheses: &[PlacementHypothesis]) -> Self {
        Self {
            seq_start: BoundaryValue::from_values(hypotheses.iter().map(|h| h.seq_start)),
            seq_end: BoundaryValue::from_values(hypotheses.iter().map(|h| h.seq_end)),
            ref_start: BoundaryValue::from_values(hypotheses.iter().map(|h| h.ref_start)),
            ref_end: BoundaryValue::from_values(hypotheses.iter().map(|h| h.ref_end)),
        }
    }
}

#[derive(Copy, Clone, Debug, Eq, PartialEq)]
pub enum LiveCallConfidence {
    Unresolved,
    Unsupported,
    ExactSingle,
    ExactAmbiguous,
}

/// Live evidence state for one segment.
#[derive(Clone, Debug, Eq, PartialEq)]
pub struct SegmentLiveCall {
    pub segment: Segment,
    pub hypotheses: Vec<PlacementHypothesis>,
    pub allele_call: AlleleBitSet,
    pub boundary_summary: BoundarySummary,
    pub confidence: LiveCallConfidence,
    pub evidence_version: u64,
}

impl SegmentLiveCall {
    pub fn unresolved(segment: Segment, allele_universe_len: usize) -> Self {
        assert_live_segment(segment);
        Self {
            segment,
            hypotheses: Vec::new(),
            allele_call: AlleleBitSet::empty(allele_universe_len),
            boundary_summary: BoundarySummary::unresolved(),
            confidence: LiveCallConfidence::Unresolved,
            evidence_version: 0,
        }
    }

    pub fn from_hypotheses(
        segment: Segment,
        allele_universe_len: usize,
        hypotheses: Vec<PlacementHypothesis>,
        evidence_version: u64,
    ) -> Self {
        assert_live_segment(segment);
        let mut allele_call = AlleleBitSet::empty(allele_universe_len);
        for hypothesis in &hypotheses {
            assert_eq!(
                hypothesis.segment, segment,
                "SegmentLiveCall cannot contain {:?} hypothesis for {:?}",
                hypothesis.segment, segment
            );
            allele_call.union_with(&hypothesis.allele_ids);
        }

        let confidence = if hypotheses.is_empty() {
            LiveCallConfidence::Unsupported
        } else if allele_call.len() == 1 {
            LiveCallConfidence::ExactSingle
        } else {
            LiveCallConfidence::ExactAmbiguous
        };

        let boundary_summary = BoundarySummary::from_hypotheses(&hypotheses);
        Self {
            segment,
            hypotheses,
            allele_call,
            boundary_summary,
            confidence,
            evidence_version,
        }
    }
}

/// Dirty live-call recomputation interval.
#[derive(Clone, Debug, Eq, PartialEq)]
pub struct DirtyWindow {
    pub start: u32,
    pub end: u32,
    pub reason: DirtyReason,
}

impl DirtyWindow {
    pub fn new(start: u32, end: u32, reason: DirtyReason) -> Self {
        assert!(start <= end, "dirty window range is inverted");
        Self { start, end, reason }
    }
}

#[derive(Clone, Debug, Eq, PartialEq)]
pub enum DirtyReason {
    AlleleSampled(Segment),
    TrimChanged(Segment),
    SegmentAssembled(Segment),
    NpGenerated(Segment),
    BaseEdited { site: u32 },
    StructuralIndel { site: u32, delta: i32 },
    ContaminantReplacement,
}

/// Dormant top-level live-call sidecar.
#[derive(Clone, Debug, Default, Eq, PartialEq)]
pub struct LiveCallState {
    pub v: Option<SegmentLiveCall>,
    pub d: Option<SegmentLiveCall>,
    pub j: Option<SegmentLiveCall>,
    pub dirty_windows: Vec<DirtyWindow>,
    pub version: u64,
}

impl LiveCallState {
    pub fn empty() -> Self {
        Self::default()
    }

    pub fn get(&self, segment: Segment) -> Option<&SegmentLiveCall> {
        match segment {
            Segment::V => self.v.as_ref(),
            Segment::D => self.d.as_ref(),
            Segment::J => self.j.as_ref(),
            Segment::Np1 | Segment::Np2 => None,
        }
    }

    pub fn with_segment_call(&self, call: SegmentLiveCall) -> Self {
        let mut next = self.clone();
        match call.segment {
            Segment::V => next.v = Some(call),
            Segment::D => next.d = Some(call),
            Segment::J => next.j = Some(call),
            Segment::Np1 | Segment::Np2 => unreachable!("SegmentLiveCall rejects NP segments"),
        }
        next.version = next.version.saturating_add(1);
        next
    }

    pub fn with_dirty_window(&self, window: DirtyWindow) -> Self {
        let mut next = self.clone();
        next.dirty_windows.push(window);
        next.version = next.version.saturating_add(1);
        next
    }
}

fn assert_live_segment(segment: Segment) {
    assert!(
        matches!(segment, Segment::V | Segment::D | Segment::J),
        "live allele calls are only defined for V/D/J segments, got {:?}",
        segment
    );
}

#[cfg(test)]
mod tests {
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
        // NP region holding T, G. NP bases use NO_GERMLINE_POS so the
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
}
