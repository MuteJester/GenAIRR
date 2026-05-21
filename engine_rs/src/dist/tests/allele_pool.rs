use super::super::{AllelePoolDist, Distribution};
use super::make_pool;
use crate::refdata::{AlleleId, AllelePool};
use crate::rng::Rng;

#[test]
#[should_panic(expected = "pool must contain at least one allele")]
fn allele_pool_dist_uniform_rejects_empty_pool() {
    let p = AllelePool::new();
    let _ = AllelePoolDist::uniform(&p);
}

#[test]
#[should_panic(expected = "weights")]
fn allele_pool_dist_from_weights_rejects_size_mismatch() {
    let p = make_pool(3);
    let _ = AllelePoolDist::from_weights(&p, vec![1.0, 2.0]);
}

#[test]
#[should_panic(expected = "weight at index")]
fn allele_pool_dist_from_weights_rejects_zero_weight() {
    let p = make_pool(2);
    let _ = AllelePoolDist::from_weights(&p, vec![0.0, 1.0]);
}

#[test]
#[should_panic(expected = "weight at index")]
fn allele_pool_dist_from_weights_rejects_nan() {
    let p = make_pool(1);
    let _ = AllelePoolDist::from_weights(&p, vec![f64::NAN]);
}

#[test]
fn allele_pool_dist_uniform_construction_round_trip() {
    let p = make_pool(5);
    let dist = AllelePoolDist::uniform(&p);
    assert_eq!(dist.len(), 5);
    assert!(!dist.is_empty());
    assert!((dist.total_weight() - 5.0).abs() < 1e-12);
}

#[test]
fn allele_pool_dist_uniform_covers_all_alleles() {
    let p = make_pool(10);
    let dist = AllelePoolDist::uniform(&p);
    let mut rng = Rng::new(99);
    let mut seen = vec![false; 10];
    for _ in 0..2000 {
        let id = dist.sample(&mut rng);
        assert!(id.as_usize() < 10, "out-of-bounds AlleleId");
        seen[id.as_usize()] = true;
    }
    for (i, &was_seen) in seen.iter().enumerate() {
        assert!(was_seen, "AlleleId({}) never sampled", i);
    }
}

#[test]
fn allele_pool_dist_uniform_is_roughly_uniform() {
    let p = make_pool(4);
    let dist = AllelePoolDist::uniform(&p);
    let mut rng = Rng::new(0xfade);
    let mut counts = [0u32; 4];
    let n = 10_000;
    for _ in 0..n {
        counts[dist.sample(&mut rng).as_usize()] += 1;
    }
    for (i, &c) in counts.iter().enumerate() {
        assert!(
            (2200..=2800).contains(&c),
            "bucket {} count {} outside [2200, 2800]",
            i,
            c
        );
    }
}

#[test]
fn allele_pool_dist_weights_are_respected() {
    let p = make_pool(3);
    let dist = AllelePoolDist::from_weights(&p, vec![0.8, 0.1, 0.1]);
    let mut rng = Rng::new(0xfeed);
    let mut counts = [0u32; 3];
    let n = 10_000;
    for _ in 0..n {
        counts[dist.sample(&mut rng).as_usize()] += 1;
    }
    assert!(
        (7500..=8500).contains(&counts[0]),
        "heavy bucket count {} outside [7500, 8500]",
        counts[0]
    );
    assert!(
        (700..=1300).contains(&counts[1]) && (700..=1300).contains(&counts[2]),
        "light bucket counts {:?} outside [700, 1300]",
        &counts[1..]
    );
}

#[test]
fn allele_pool_dist_same_seed_same_stream() {
    let p = make_pool(7);
    let dist = AllelePoolDist::uniform(&p);
    let mut a = Rng::new(0xcafe);
    let mut b = Rng::new(0xcafe);
    for _ in 0..100 {
        assert_eq!(dist.sample(&mut a), dist.sample(&mut b));
    }
}

#[test]
fn allele_pool_dist_single_allele_always_returned() {
    let p = make_pool(1);
    let dist = AllelePoolDist::uniform(&p);
    let mut rng = Rng::new(1);
    for _ in 0..100 {
        assert_eq!(dist.sample(&mut rng), AlleleId::new(0));
    }
}

#[test]
fn allele_pool_dist_works_through_box_dyn() {
    let p = make_pool(3);
    let dist: Box<dyn Distribution<Output = AlleleId>> = Box::new(AllelePoolDist::uniform(&p));
    let mut rng = Rng::new(0);
    let id = dist.sample(&mut rng);
    assert!(id.as_usize() < 3);
}

#[test]
fn allele_pool_dist_sampled_ids_resolve_in_pool() {
    let p = make_pool(5);
    let dist = AllelePoolDist::uniform(&p);
    let mut rng = Rng::new(42);
    for _ in 0..50 {
        let id = dist.sample(&mut rng);
        let resolved = p.get(id);
        assert!(resolved.is_some(), "sampled id {:?} did not resolve", id);
    }
}

#[test]
fn allele_pool_dist_restricted_uniform_only_samples_allowed_ids() {
    let p = make_pool(10);
    let allowed = vec![AlleleId::new(2), AlleleId::new(5), AlleleId::new(7)];
    let dist = AllelePoolDist::restricted_uniform(&p, allowed.clone());
    assert_eq!(dist.len(), 3);
    let mut rng = Rng::new(0xbeef);
    for _ in 0..500 {
        let id = dist.sample(&mut rng);
        assert!(
            allowed.contains(&id),
            "sampled id {:?} not in allowed set",
            id
        );
    }
}

#[test]
fn allele_pool_dist_restricted_uniform_covers_all_allowed() {
    let p = make_pool(8);
    let allowed = vec![AlleleId::new(1), AlleleId::new(4), AlleleId::new(6)];
    let dist = AllelePoolDist::restricted_uniform(&p, allowed.clone());
    let mut rng = Rng::new(0xc0de);
    let mut seen = vec![false; allowed.len()];
    for _ in 0..2000 {
        let id = dist.sample(&mut rng);
        let pos = allowed.iter().position(|&a| a == id).unwrap();
        seen[pos] = true;
    }
    assert!(seen.iter().all(|&s| s), "not every allowed id was sampled");
}

#[test]
fn allele_pool_dist_restricted_uniform_single_id_always_returns_it() {
    let p = make_pool(20);
    let dist = AllelePoolDist::restricted_uniform(&p, vec![AlleleId::new(13)]);
    let mut rng = Rng::new(7);
    for _ in 0..50 {
        assert_eq!(dist.sample(&mut rng), AlleleId::new(13));
    }
}

#[test]
fn allele_pool_dist_restricted_uniform_support_round_trip() {
    let p = make_pool(6);
    let allowed = vec![AlleleId::new(0), AlleleId::new(3), AlleleId::new(5)];
    let dist = AllelePoolDist::restricted_uniform(&p, allowed.clone());
    let support = dist.support().expect("support is Some");
    let support_ids: Vec<AlleleId> = support.iter().map(|(id, _)| *id).collect();
    assert_eq!(support_ids, allowed);
    for (_, weight) in &support {
        assert!((*weight - 1.0).abs() < 1e-12);
    }
}

#[test]
#[should_panic(expected = "allowed_ids must be non-empty")]
fn allele_pool_dist_restricted_uniform_rejects_empty() {
    let p = make_pool(3);
    let _ = AllelePoolDist::restricted_uniform(&p, vec![]);
}

#[test]
#[should_panic(expected = "out of range")]
fn allele_pool_dist_restricted_uniform_rejects_out_of_range_id() {
    let p = make_pool(3);
    let _ = AllelePoolDist::restricted_uniform(&p, vec![AlleleId::new(5)]);
}

#[test]
#[should_panic(expected = "duplicate")]
fn allele_pool_dist_restricted_uniform_rejects_duplicate_ids() {
    let p = make_pool(3);
    let _ = AllelePoolDist::restricted_uniform(&p, vec![AlleleId::new(1), AlleleId::new(1)]);
}
