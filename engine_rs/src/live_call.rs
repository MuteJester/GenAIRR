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

    let mut candidates = segment_index.all_alleles.clone();
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
        // skip it without failing the call. The candidate set is
        // unchanged; downstream walkers see one fewer "evidence
        // position" but otherwise behave normally.
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

        let Some(evidence) = segment_index.compatible_alleles_at(ref_pos as usize, nucleotide.base)
        else {
            return unsupported_call(segment, allele_universe_len, evidence_version);
        };

        if evidence.informative {
            informative_matches = informative_matches.saturating_add(1);
        } else {
            wildcard_matches = wildcard_matches.saturating_add(1);
        }

        candidates.intersect_with(&evidence.allele_ids);
        if candidates.is_empty() {
            return unsupported_call(segment, allele_universe_len, evidence_version);
        }
    }

    let mut ref_start = ref_start.expect("non-empty region should set ref_start");
    let mut ref_end = next_ref_pos.expect("non-empty region should set ref_end");
    let mut seq_start = region.start.index();
    let mut seq_end = region.end.index();
    let mut flags = HypothesisFlags::EMPTY;

    // ── Left-side extension + overlap walk ───────────────────────────
    //
    // If `left_extension` (NP1 for J on a VJ chain, NP2 for J on a
    // VDJ chain — or NP1 for D) sits immediately to the left of
    // `region`, walk backward through it (and beyond it into
    // earlier coding territory) as long as each base extends some
    // candidate allele's reference *prefix*.
    //
    // Mark `BOUNDARY_ELASTIC` whenever any extension occurs.
    // Additionally mark `OVERLAPS_OTHER_SEGMENT` when the walker
    // crosses past `np_region.start` into the next-earlier region's
    // bases — that's the live-graph overlap signal the design doc
    // calls out (V/D and D/J overlap when bases support both).
    //
    // Halt semantics unchanged: stop when the candidate set empties,
    // when ref_pos would go below 0, or when we run out of pool
    // positions; never fail the call.
    if let Some(np_region) = left_extension {
        if np_region.end.index() == seq_start {
            let np_start = np_region.start.index();
            let mut extended_seq_start = seq_start;
            let mut extended_ref_start = ref_start;
            let mut extended_candidates = candidates.clone();
            let mut extended_inform = informative_matches;
            let mut extended_wildcard = wildcard_matches;
            let mut extension_added = false;
            let mut crossed_into_overlap = false;

            // Walk pool positions in REVERSE order from `np_region.end`
            // down to 0, continuing past `np_region.start` so V's left
            // can reach into an earlier coding region's bases when
            // evidence supports it.
            for seq_pos in (0..np_region.end.index()).rev() {
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
                let mut next_candidates = extended_candidates.clone();
                next_candidates.intersect_with(&evidence.allele_ids);
                if next_candidates.is_empty() {
                    break;
                }
                extended_candidates = next_candidates;
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
            }

            if extension_added {
                candidates = extended_candidates;
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
    // Walk forward from `region.end` as long as each base extends
    // some allele's reference suffix. The walker continues into the
    // next-later coding region's bases (D for V, J for D) when
    // evidence supports both placements simultaneously, marking
    // `OVERLAPS_OTHER_SEGMENT` on the resulting hypothesis.
    if let Some(np_region) = right_extension {
        if np_region.start.index() == seq_end {
            let np_end = np_region.end.index();
            let pool_len = sim.pool.len() as u32;
            let mut extended_seq_end = seq_end;
            let mut extended_ref_pos = ref_end;
            let mut extended_candidates = candidates.clone();
            let mut extended_inform = informative_matches;
            let mut extended_wildcard = wildcard_matches;
            let mut extension_added = false;
            let mut crossed_into_overlap = false;

            for seq_pos in np_region.start.index()..pool_len {
                let nucleotide = sim
                    .pool
                    .get(NucHandle::new(seq_pos))
                    .expect("extension range must point into the nucleotide pool");
                let Some(evidence) = segment_index
                    .compatible_alleles_at(extended_ref_pos as usize, nucleotide.base)
                else {
                    break;
                };
                let mut next_candidates = extended_candidates.clone();
                next_candidates.intersect_with(&evidence.allele_ids);
                if next_candidates.is_empty() {
                    break;
                }
                extended_candidates = next_candidates;
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
            }

            if extension_added {
                candidates = extended_candidates;
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

    let hypothesis = PlacementHypothesis::new(
        segment,
        seq_start,
        seq_end,
        ref_start,
        ref_end,
        candidates,
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
