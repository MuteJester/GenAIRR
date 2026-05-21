use super::super::{Distribution, UniformBase, UniformInt};
use crate::rng::Rng;

#[test]
fn uniform_base_only_emits_canonical_bases() {
    let mut rng = Rng::new(1);
    let dist = UniformBase;
    for _ in 0..1000 {
        let b = dist.sample(&mut rng);
        assert!(
            matches!(b, b'A' | b'C' | b'G' | b'T'),
            "UniformBase emitted unexpected byte 0x{:02x}",
            b
        );
    }
}

#[test]
fn uniform_base_covers_all_four_bases() {
    let mut rng = Rng::new(123);
    let dist = UniformBase;
    let mut seen = [false; 4];
    for _ in 0..1000 {
        let b = dist.sample(&mut rng);
        let idx = match b {
            b'A' => 0,
            b'C' => 1,
            b'G' => 2,
            b'T' => 3,
            _ => panic!("unexpected base"),
        };
        seen[idx] = true;
    }
    for (i, &was_seen) in seen.iter().enumerate() {
        assert!(was_seen, "base index {} never appeared in 1000 draws", i);
    }
}

#[test]
fn uniform_base_distribution_is_roughly_uniform() {
    let mut rng = Rng::new(0xfeed);
    let dist = UniformBase;
    let mut counts = [0u32; 4];
    let n = 10_000u32;
    for _ in 0..n {
        let idx = match dist.sample(&mut rng) {
            b'A' => 0,
            b'C' => 1,
            b'G' => 2,
            b'T' => 3,
            _ => unreachable!(),
        };
        counts[idx] += 1;
    }
    for &c in &counts {
        assert!(
            (2200..=2800).contains(&c),
            "UniformBase bucket count {} outside [2200, 2800]",
            c
        );
    }
}

#[test]
fn uniform_base_same_seed_same_stream() {
    let mut a = Rng::new(7);
    let mut b = Rng::new(7);
    let dist = UniformBase;
    for _ in 0..100 {
        assert_eq!(dist.sample(&mut a), dist.sample(&mut b));
    }
}

#[test]
#[should_panic(expected = "UniformInt: max")]
fn uniform_int_new_rejects_max_le_min() {
    let _ = UniformInt::new(5, 5);
}

#[test]
#[should_panic(expected = "UniformInt: max")]
fn uniform_int_new_rejects_inverted_range() {
    let _ = UniformInt::new(10, 5);
}

#[test]
#[should_panic(expected = "UniformInt: span")]
fn uniform_int_new_rejects_oversized_span() {
    let _ = UniformInt::new(0, (u32::MAX as i64) + 2);
}

#[test]
fn uniform_int_new_accepts_max_span() {
    let dist = UniformInt::new(0, u32::MAX as i64);
    assert_eq!(dist.span(), u32::MAX as u64);
}

#[test]
fn uniform_int_stays_in_bounds() {
    let mut rng = Rng::new(42);
    let dist = UniformInt::new(3, 10);
    for _ in 0..1000 {
        let v = dist.sample(&mut rng);
        assert!(v >= 3, "UniformInt(3,10) produced {} (< 3)", v);
        assert!(v < 10, "UniformInt(3,10) produced {} (>= 10)", v);
    }
}

#[test]
fn uniform_int_covers_full_range() {
    let mut rng = Rng::new(99);
    let dist = UniformInt::new(0, 5);
    let mut seen = [false; 5];
    for _ in 0..1000 {
        let v = dist.sample(&mut rng);
        seen[v as usize] = true;
    }
    for (i, &was_seen) in seen.iter().enumerate() {
        assert!(was_seen, "UniformInt(0,5) never produced {}", i);
    }
}

#[test]
fn uniform_int_negative_min_works() {
    let mut rng = Rng::new(11);
    let dist = UniformInt::new(-3, 4);
    for _ in 0..1000 {
        let v = dist.sample(&mut rng);
        assert!(v >= -3, "{}", v);
        assert!(v < 4, "{}", v);
    }
}

#[test]
fn uniform_int_span_one_returns_min() {
    let mut rng = Rng::new(1);
    let dist = UniformInt::new(7, 8);
    for _ in 0..100 {
        assert_eq!(dist.sample(&mut rng), 7);
    }
}

#[test]
fn uniform_int_same_seed_same_stream() {
    let mut a = Rng::new(0xbabe);
    let mut b = Rng::new(0xbabe);
    let dist = UniformInt::new(0, 100);
    for _ in 0..100 {
        assert_eq!(dist.sample(&mut a), dist.sample(&mut b));
    }
}

#[test]
fn uniform_int_accessors_round_trip() {
    let dist = UniformInt::new(-5, 12);
    assert_eq!(dist.min(), -5);
    assert_eq!(dist.max(), 12);
    assert_eq!(dist.span(), 17);
}
