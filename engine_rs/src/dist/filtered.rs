//! Filtered sampling — draw from a distribution restricted to a predicate.
//!
//! The architectural primitive for constraint-aware sampling: a sampling
//! pass calls [`sample_filtered_result`] with a `predicate` wrapping a
//! `ContractSet::admits` check, and the pass either receives an
//! admissible value or a structured error explaining why no admissible
//! candidate exists.

use super::Distribution;
use crate::rng::Rng;

/// Why a filtered sample could not be produced.
#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub enum FilteredSampleError {
    /// The distribution does not expose enumerable support.
    SupportUnavailable,
    /// Enumerable support existed, but every candidate was rejected.
    EmptyAdmissibleSupport,
    /// Filtered support had invalid weights.
    InvalidFilteredSupport,
}

/// Draw one sample from `dist` restricted to values satisfying
/// `predicate`, preserving the reason if no filtered sample exists.
///
/// Returns `Err` when:
/// - the distribution doesn't expose a `support()` (continuous /
///   too-large categorical),
/// - the filtered support is empty (no admissible candidates), or
/// - the filtered support weights are invalid.
///
/// **The architectural primitive for constraint-aware sampling.**
/// Sampling passes (D.6) call this with `predicate` wrapping a
/// `ContractSet::admits` check. The pass receives `Ok(value)` if
/// any admissible value exists, or a structured error explaining why
/// a strict runtime cannot sample under the requested constraints.
///
/// Determinism: consumes exactly one RNG word only when a filtered
/// sample can be drawn. Failure paths consume no RNG words.
pub fn sample_filtered_result<T, D, F>(
    rng: &mut Rng,
    dist: &D,
    predicate: F,
) -> Result<T, FilteredSampleError>
where
    T: Clone,
    D: Distribution<Output = T> + ?Sized,
    F: Fn(&T) -> bool,
{
    let support = dist
        .support()
        .ok_or(FilteredSampleError::SupportUnavailable)?;
    let mut filtered: Vec<(T, f64)> = support.into_iter().filter(|(v, _)| predicate(v)).collect();

    if filtered.is_empty() {
        return Err(FilteredSampleError::EmptyAdmissibleSupport);
    }

    // Renormalize and sample by inverse CDF.
    let total: f64 = filtered.iter().map(|(_, w)| w).sum();
    if !total.is_finite() || total <= 0.0 {
        return Err(FilteredSampleError::InvalidFilteredSupport);
    }
    let r = rng.next_f64() * total;
    let mut cum = 0.0;
    for (v, w) in filtered.iter() {
        cum += w;
        if r < cum {
            return Ok(v.clone());
        }
    }
    // Defensive last-index fallback for ULP rounding.
    Ok(filtered.pop().unwrap().0)
}

/// Permissive filtered sampling helper.
///
/// Callers get a filtered value when possible and `None` when
/// filtering cannot be performed. Strict runtime paths should call
/// [`sample_filtered_result`] so they can surface structured errors
/// instead of silently falling back.
pub fn sample_filtered<T, D, F>(rng: &mut Rng, dist: &D, predicate: F) -> Option<T>
where
    T: Clone,
    D: Distribution<Output = T> + ?Sized,
    F: Fn(&T) -> bool,
{
    sample_filtered_result(rng, dist, predicate).ok()
}
