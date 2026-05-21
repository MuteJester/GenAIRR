use genairr_engine::dist::{Distribution, EmpiricalLengthDist, UniformInt};
use genairr_engine::rng::Rng;

#[test]
fn property_uniform_int_samples_lie_in_declared_support() {
    let dist = UniformInt::new(-3, 7);
    let support = dist.support().expect("UniformInt has finite support");
    let support_set: std::collections::HashSet<i64> = support.iter().map(|(v, _)| *v).collect();

    let mut rng = Rng::new(0xdead_beef);
    for _ in 0..2000 {
        let v = dist.sample(&mut rng);
        assert!(
            support_set.contains(&v),
            "sampled value {} outside declared support",
            v
        );
    }
}

#[test]
fn property_empirical_length_dist_samples_lie_in_declared_support() {
    let dist = EmpiricalLengthDist::from_pairs(vec![(2, 1.0), (5, 2.0), (11, 0.5)]);
    let support = dist
        .support()
        .expect("EmpiricalLengthDist has finite support");
    let support_set: std::collections::HashSet<i64> = support.iter().map(|(v, _)| *v).collect();

    let mut rng = Rng::new(0xfeed_face);
    for _ in 0..2000 {
        let v = dist.sample(&mut rng);
        assert!(
            support_set.contains(&v),
            "sampled value {} outside declared support",
            v
        );
    }
}
