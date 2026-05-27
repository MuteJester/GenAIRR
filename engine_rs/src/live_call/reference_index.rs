use super::{assert_live_segment, AlleleBitSet};
use crate::ir::Segment;
use crate::refdata::{Allele, AlleleId, AllelePool, RefDataConfig};
use std::borrow::Cow;
use std::collections::HashMap;

/// Default seed length for compiled reference lookups.
///
/// Segment indexes automatically fall back to a shorter value when all
/// alleles in a segment are shorter than this, which matters for D
/// alleles and tiny synthetic test fixtures.
pub const DEFAULT_REFERENCE_KMER_LEN: usize = 7;

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

    pub fn compatible_alleles_at(
        &self,
        ref_pos: usize,
        observed_base: u8,
    ) -> Option<BaseEvidence<'_>> {
        self.base_at_pos
            .get(ref_pos)
            .and_then(|base_sets| base_sets.compatible_with_observed(observed_base))
    }

    pub fn kmer_hits(&self, kmer: &[u8]) -> Option<&[KmerHit]> {
        self.kmer_index.hits(kmer)
    }
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

    pub fn compatible_with_observed(&self, observed_base: u8) -> Option<BaseEvidence<'_>> {
        match observed_base_kind(observed_base)? {
            // Canonical-base path: hand back a borrow into the
            // precomputed per-base bitset. The walker only reads from
            // it (for_each_id / is_empty), so no clone is needed in
            // the hot path. This is the dominant case.
            ObservedBaseKind::Canonical(index) => Some(BaseEvidence {
                allele_ids: Cow::Borrowed(&self.sets[index]),
                informative: true,
            }),
            // Wildcard `N`: still owned because the union has to be
            // materialised. Hit rate on this path is low.
            ObservedBaseKind::Wildcard => {
                let mut allele_ids = AlleleBitSet::empty(self.sets[0].universe_len());
                for set in &self.sets {
                    allele_ids.union_with(set);
                }
                Some(BaseEvidence {
                    allele_ids: Cow::Owned(allele_ids),
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
///
/// `allele_ids` is a `Cow` because the canonical-base lookup (the
/// overwhelmingly common path in the walker) can hand back a borrow
/// into the precomputed per-base bitset stored on `BaseBitsets`,
/// avoiding a per-position `Vec<u64>` allocation. The wildcard `N`
/// case still owns its bitset because the union across the four
/// canonical-base sets has to be materialised.
#[derive(Clone, Debug, Eq, PartialEq)]
pub struct BaseEvidence<'a> {
    pub allele_ids: Cow<'a, AlleleBitSet>,
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
    // Route through the shared scoring kernel: scoring::classify_base
    // is the single source of truth for "what counts as canonical /
    // wildcard / invalid". The inverted-index lookup just needs the
    // 0..4 slot index for canonical bases.
    match super::scoring::classify_base(base) {
        super::scoring::BaseKind::Canonical(b) => canonical_base_index(b).map(ObservedBaseKind::Canonical),
        super::scoring::BaseKind::Wildcard => Some(ObservedBaseKind::Wildcard),
        super::scoring::BaseKind::Invalid => None,
    }
}

fn canonical_base_index(base: u8) -> Option<usize> {
    // Slot index for the inverted-index data structure (A=0..T=3).
    // The "is this canonical" question is delegated to scoring::classify_base
    // via observed_base_kind; this function is purely the slot lookup.
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
