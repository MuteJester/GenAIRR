use super::super::{
    sample_filtered, sample_filtered_result, Distribution, FilteredSampleError, UniformInt,
};
use crate::rng::Rng;

#[derive(Clone, Debug)]
struct NoSupportDist;

impl Distribution for NoSupportDist {
    type Output = i64;

    fn sample(&self, _rng: &mut Rng) -> i64 {
        0
    }
}

#[test]
fn sample_filtered_result_reports_support_unavailable() {
    let mut rng = Rng::new(1);
    let err = sample_filtered_result(&mut rng, &NoSupportDist, |_| true).unwrap_err();
    assert_eq!(err, FilteredSampleError::SupportUnavailable);
}

#[test]
fn sample_filtered_result_reports_empty_admissible_support() {
    let mut rng = Rng::new(1);
    let dist = UniformInt::new(0, 4);
    let err = sample_filtered_result(&mut rng, &dist, |_| false).unwrap_err();
    assert_eq!(err, FilteredSampleError::EmptyAdmissibleSupport);
}

#[test]
fn sample_filtered_result_samples_from_admissible_subset() {
    let mut rng = Rng::new(1);
    let dist = UniformInt::new(0, 10);
    for _ in 0..100 {
        let value = sample_filtered_result(&mut rng, &dist, |v| *v >= 7).unwrap();
        assert!((7..10).contains(&value));
    }
}

#[test]
fn sample_filtered_permissive_collapses_filter_errors_to_none() {
    let mut rng = Rng::new(1);
    let dist = UniformInt::new(0, 4);
    assert_eq!(sample_filtered(&mut rng, &dist, |_| false), None);
    assert_eq!(sample_filtered(&mut rng, &NoSupportDist, |_| true), None);
}
