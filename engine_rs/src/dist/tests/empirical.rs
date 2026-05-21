use super::super::{Distribution, EmpiricalLengthDist};
use crate::rng::Rng;

#[test]
#[should_panic(expected = "must have at least one entry")]
fn empirical_rejects_empty_histogram() {
    let _ = EmpiricalLengthDist::from_pairs(std::iter::empty());
}

#[test]
#[should_panic(expected = "weight at index")]
fn empirical_rejects_zero_weight() {
    let _ = EmpiricalLengthDist::from_pairs(vec![(5, 0.0), (6, 1.0)]);
}

#[test]
#[should_panic(expected = "weight at index")]
fn empirical_rejects_negative_weight() {
    let _ = EmpiricalLengthDist::from_pairs(vec![(5, -1.0)]);
}

#[test]
#[should_panic(expected = "weight at index")]
fn empirical_rejects_nan_weight() {
    let _ = EmpiricalLengthDist::from_pairs(vec![(5, f64::NAN)]);
}

#[test]
#[should_panic(expected = "weight at index")]
fn empirical_rejects_infinite_weight() {
    let _ = EmpiricalLengthDist::from_pairs(vec![(5, f64::INFINITY)]);
}

#[test]
fn empirical_single_value_always_returned() {
    let dist = EmpiricalLengthDist::from_pairs(vec![(42, 1.0)]);
    let mut rng = Rng::new(0xc0ff_ee);
    for _ in 0..1000 {
        assert_eq!(dist.sample(&mut rng), 42);
    }
}

#[test]
fn empirical_construction_accessors_round_trip() {
    let dist = EmpiricalLengthDist::from_pairs(vec![(0, 1.0), (1, 2.0), (2, 1.0)]);
    assert_eq!(dist.len(), 3);
    assert!(!dist.is_empty());
    assert_eq!(dist.values(), &[0, 1, 2]);
    assert!((dist.total_weight() - 4.0).abs() < 1e-12);
}

#[test]
fn empirical_from_values_and_weights_matches_from_pairs() {
    let pairs = EmpiricalLengthDist::from_pairs(vec![(0, 1.0), (1, 2.0), (2, 3.0)]);
    let split = EmpiricalLengthDist::from_values_and_weights(vec![0, 1, 2], vec![1.0, 2.0, 3.0]);

    let mut a = Rng::new(7);
    let mut b = Rng::new(7);
    for _ in 0..100 {
        assert_eq!(pairs.sample(&mut a), split.sample(&mut b));
    }
}

#[test]
#[should_panic(expected = "values and weights must have equal lengths")]
fn empirical_from_values_and_weights_rejects_length_mismatch() {
    let _ = EmpiricalLengthDist::from_values_and_weights(vec![0, 1], vec![1.0]);
}

#[test]
fn empirical_stays_in_value_set() {
    let dist = EmpiricalLengthDist::from_pairs(vec![(3, 1.0), (7, 2.0), (11, 1.0), (-2, 0.5)]);
    let mut rng = Rng::new(17);
    for _ in 0..1000 {
        let v = dist.sample(&mut rng);
        assert!(matches!(v, 3 | 7 | 11 | -2));
    }
}

#[test]
fn empirical_covers_full_value_set() {
    let dist = EmpiricalLengthDist::from_pairs(vec![(0, 1.0), (1, 1.0), (2, 1.0), (3, 1.0)]);
    let mut rng = Rng::new(99);
    let mut seen = [false; 4];
    for _ in 0..1000 {
        let v = dist.sample(&mut rng);
        seen[v as usize] = true;
    }
    for (i, &was_seen) in seen.iter().enumerate() {
        assert!(was_seen, "value {} never appeared", i);
    }
}

#[test]
fn empirical_weights_are_respected() {
    let dist = EmpiricalLengthDist::from_pairs(vec![(0, 0.9), (1, 0.1)]);
    let mut rng = Rng::new(0xabba);
    let mut zero_count = 0;
    let n = 10_000;
    for _ in 0..n {
        if dist.sample(&mut rng) == 0 {
            zero_count += 1;
        }
    }
    assert!(
        (8500..=9500).contains(&zero_count),
        "expected zero_count ~9000 of 10000, got {}",
        zero_count
    );
}

#[test]
fn empirical_uniform_weights_produce_roughly_uniform_distribution() {
    let dist =
        EmpiricalLengthDist::from_pairs(vec![(0, 1.0), (1, 1.0), (2, 1.0), (3, 1.0), (4, 1.0)]);
    let mut rng = Rng::new(0xdead);
    let mut counts = [0u32; 5];
    let n = 10_000u32;
    for _ in 0..n {
        counts[dist.sample(&mut rng) as usize] += 1;
    }
    for (i, &c) in counts.iter().enumerate() {
        assert!(
            (1700..=2300).contains(&c),
            "bucket {} count {} outside [1700, 2300]",
            i,
            c
        );
    }
}

#[test]
fn empirical_same_seed_same_stream() {
    let dist = EmpiricalLengthDist::from_pairs(vec![(10, 1.0), (20, 2.0), (30, 3.0)]);
    let mut a = Rng::new(0xfeed);
    let mut b = Rng::new(0xfeed);
    for _ in 0..100 {
        assert_eq!(dist.sample(&mut a), dist.sample(&mut b));
    }
}

#[test]
fn empirical_negative_values_supported() {
    let dist = EmpiricalLengthDist::from_pairs(vec![(-5, 1.0), (0, 1.0), (5, 1.0)]);
    let mut rng = Rng::new(13);
    for _ in 0..100 {
        let v = dist.sample(&mut rng);
        assert!(matches!(v, -5 | 0 | 5));
    }
}
