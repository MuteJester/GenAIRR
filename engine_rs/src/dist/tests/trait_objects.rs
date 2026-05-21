use super::super::{Distribution, EmpiricalLengthDist, UniformBase, UniformInt};
use crate::rng::Rng;

#[test]
fn box_dyn_distribution_dispatches_correctly() {
    let mut rng = Rng::new(1);

    let base_dist: Box<dyn Distribution<Output = u8>> = Box::new(UniformBase);
    let int_dist: Box<dyn Distribution<Output = i64>> = Box::new(UniformInt::new(0, 10));

    let b = base_dist.sample(&mut rng);
    let i = int_dist.sample(&mut rng);

    assert!(matches!(b, b'A' | b'C' | b'G' | b'T'));
    assert!((0..10).contains(&i));
}

#[test]
fn vec_of_homogeneous_dyn_distributions() {
    let dists: Vec<Box<dyn Distribution<Output = u8>>> = vec![
        Box::new(UniformBase),
        Box::new(UniformBase),
        Box::new(UniformBase),
    ];
    let mut rng = Rng::new(0);
    for d in &dists {
        let b = d.sample(&mut rng);
        assert!(matches!(b, b'A' | b'C' | b'G' | b'T'));
    }
}

#[test]
fn empirical_works_through_box_dyn() {
    let dist: Box<dyn Distribution<Output = i64>> =
        Box::new(EmpiricalLengthDist::from_pairs(vec![(7, 1.0)]));
    let mut rng = Rng::new(0);
    assert_eq!(dist.sample(&mut rng), 7);
}
