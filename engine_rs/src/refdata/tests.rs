use super::*;
use crate::ir::Segment;
use std::mem::size_of;

/// Helper: a tiny synthetic V allele for tests.
fn make_v(name: &str, gene: &str, seq: &[u8], anchor: Option<u16>) -> Allele {
    Allele {
        name: name.to_string(),
        gene: gene.to_string(),
        seq: seq.to_vec(),
        segment: Segment::V,
        anchor,
    }
}

#[test]
fn allele_id_is_zero_cost_newtype() {
    assert_eq!(size_of::<AlleleId>(), size_of::<u32>());
}

#[test]
fn allele_id_round_trip() {
    let id = AlleleId::new(7);
    assert_eq!(id.index(), 7);
    assert_eq!(id.as_usize(), 7);
}

#[test]
fn chain_type_has_d_distinction() {
    assert!(!ChainType::Vj.has_d());
    assert!(ChainType::Vdj.has_d());
}

#[test]
fn allele_basic_accessors() {
    let a = make_v("IGHV1-2*01", "IGHV1-2", b"ACGTACGT", Some(3));
    assert_eq!(a.len(), 8);
    assert!(a.has_anchor());
    assert_eq!(a.anchor, Some(3));
    assert_eq!(a.segment, Segment::V);
}

#[test]
fn allele_anchorless_round_trip() {
    let a = make_v("IGHV-pseudo*01", "IGHV-pseudo", b"ACGT", None);
    assert!(!a.has_anchor());
    assert_eq!(a.anchor, None);
}

#[test]
fn allele_pool_starts_empty() {
    let p = AllelePool::new();
    assert_eq!(p.len(), 0);
    assert!(p.is_empty());
    assert!(p.get(AlleleId::new(0)).is_none());
    assert!(p.find_by_name("nonexistent").is_none());
}

#[test]
fn allele_pool_push_returns_sequential_ids() {
    let mut p = AllelePool::new();
    let id0 = p.push(make_v("a*01", "a", b"AA", Some(0)));
    let id1 = p.push(make_v("b*01", "b", b"CC", Some(0)));
    let id2 = p.push(make_v("c*01", "c", b"GG", None));

    assert_eq!(id0.index(), 0);
    assert_eq!(id1.index(), 1);
    assert_eq!(id2.index(), 2);
    assert_eq!(p.len(), 3);
}

#[test]
fn allele_pool_get_returns_stored_allele() {
    let mut p = AllelePool::new();
    let id = p.push(make_v("x*01", "x", b"ATGC", Some(2)));
    let got = p.get(id).expect("just pushed");
    assert_eq!(got.name, "x*01");
    assert_eq!(got.seq, b"ATGC");
    assert_eq!(got.anchor, Some(2));
}

#[test]
fn allele_pool_get_out_of_bounds_returns_none() {
    let p = AllelePool::new();
    assert!(p.get(AlleleId::new(99)).is_none());
}

#[test]
fn allele_pool_iter_yields_id_allele_pairs_in_order() {
    let mut p = AllelePool::new();
    let _ = p.push(make_v("a*01", "a", b"AA", None));
    let _ = p.push(make_v("b*01", "b", b"CC", None));
    let _ = p.push(make_v("c*01", "c", b"GG", None));

    let collected: Vec<(u32, String)> = p
        .iter()
        .map(|(id, a)| (id.index(), a.name.clone()))
        .collect();
    assert_eq!(
        collected,
        vec![
            (0, "a*01".to_string()),
            (1, "b*01".to_string()),
            (2, "c*01".to_string()),
        ]
    );
}

#[test]
fn allele_pool_find_by_name_locates_exact_match() {
    let mut p = AllelePool::new();
    let _ = p.push(make_v("IGHV1-2*01", "IGHV1-2", b"AA", None));
    let target_id = p.push(make_v("IGHV1-2*02", "IGHV1-2", b"AC", None));
    let _ = p.push(make_v("IGHV3-23*01", "IGHV3-23", b"GG", None));

    let (id, allele) = p.find_by_name("IGHV1-2*02").expect("name should exist");
    assert_eq!(id, target_id);
    assert_eq!(allele.seq, b"AC");

    // Partial / similar names should not match.
    assert!(p.find_by_name("IGHV1-2").is_none());
    assert!(p.find_by_name("IGHV1-2*03").is_none());
}

#[test]
fn ref_data_config_empty_for_chain_type() {
    let cfg = RefDataConfig::empty(ChainType::Vdj);
    assert_eq!(cfg.chain_type, ChainType::Vdj);
    assert!(cfg.v_pool.is_empty());
    assert!(cfg.d_pool.is_empty());
    assert!(cfg.j_pool.is_empty());
    assert!(cfg.c_pool.is_empty());
}

#[test]
fn ref_data_config_pool_for_segment_routes_correctly() {
    let mut cfg = RefDataConfig::empty(ChainType::Vdj);
    let _ = cfg.v_pool.push(make_v("v*01", "v", b"AA", None));
    let _ = cfg.d_pool.push(Allele {
        name: "d*01".into(),
        gene: "d".into(),
        seq: b"GG".to_vec(),
        segment: Segment::D,
        anchor: None,
    });
    let _ = cfg.j_pool.push(Allele {
        name: "j*01".into(),
        gene: "j".into(),
        seq: b"TT".to_vec(),
        segment: Segment::J,
        anchor: Some(0),
    });

    assert_eq!(cfg.pool_for(Segment::V).unwrap().len(), 1);
    assert_eq!(cfg.pool_for(Segment::D).unwrap().len(), 1);
    assert_eq!(cfg.pool_for(Segment::J).unwrap().len(), 1);
    assert!(cfg.pool_for(Segment::Np1).is_none());
    assert!(cfg.pool_for(Segment::Np2).is_none());
}

#[test]
fn ref_data_config_get_resolves_segment_id_pair() {
    let mut cfg = RefDataConfig::empty(ChainType::Vdj);
    let v_id = cfg.v_pool.push(make_v("v*01", "v", b"AAAT", Some(1)));

    let v = cfg.get(Segment::V, v_id).expect("v*01 should resolve");
    assert_eq!(v.name, "v*01");
    assert_eq!(v.anchor, Some(1));

    // Wrong segment -> None.
    assert!(cfg.get(Segment::J, v_id).is_none());
    // NP segment -> None defensively.
    assert!(cfg.get(Segment::Np1, v_id).is_none());
}

#[test]
fn ref_data_config_supports_vj_chain_with_empty_d_pool() {
    // VJ chains: d_pool is conventionally empty. Construction
    // does not enforce this; assembly (C.8) handles VJ vs VDJ
    // explicitly via chain_type.has_d().
    let mut cfg = RefDataConfig::empty(ChainType::Vj);
    let _ = cfg.v_pool.push(make_v("v*01", "v", b"AA", Some(0)));
    let _ = cfg.j_pool.push(Allele {
        name: "j*01".into(),
        gene: "j".into(),
        seq: b"TT".to_vec(),
        segment: Segment::J,
        anchor: Some(0),
    });

    assert!(!cfg.chain_type.has_d());
    assert!(cfg.d_pool.is_empty());
    assert_eq!(cfg.v_pool.len(), 1);
    assert_eq!(cfg.j_pool.len(), 1);
}
