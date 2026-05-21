use super::super::{KmerHit, ReferenceMatchIndex, SegmentRefIndex, DEFAULT_REFERENCE_KMER_LEN};
use super::allele;
use crate::ir::Segment;
use crate::refdata::{AllelePool, ChainType, RefDataConfig};

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
