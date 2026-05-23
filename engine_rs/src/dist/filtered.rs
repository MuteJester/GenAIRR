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

/// Sample one canonical base (`A`/`C`/`G`/`T`) from `dist` restricted
/// to the admitted bases indicated by `admit_mask`.
///
/// `admit_mask` follows the convention from
/// `contract::admit_mask_observer`: bit 0 = `A`, bit 1 = `C`,
/// bit 2 = `G`, bit 3 = `T`. A base is admitted iff its bit is set.
///
/// Compared to [`sample_filtered_result`] over a predicate that
/// performs the same admit check, this function:
/// - Pulls the distribution's support exactly once (one `Vec`
///   allocation) instead of materialising a separate `filtered`
///   Vec on top.
/// - Iterates the support twice in-place (sum + inverse CDF)
///   without temporary storage.
/// - Performs a single bit-test per candidate instead of a contract
///   trait dispatch chain.
///
/// Returns `Err(FilteredSampleError)` on the same conditions
/// `sample_filtered_result` would — empty admissible support,
/// missing support, or invalid total weight.
///
/// Determinism: consumes exactly one RNG word on the success path.
/// Non-canonical bases in the support are treated as unadmitted
/// (they cannot map to any bit in the 4-bit mask) and contribute
/// zero weight; this matches the predicate path's behaviour when
/// the contract dispatch returns `Err` for non-canonical
/// candidates.
pub fn sample_base_with_admit_mask<D>(
    rng: &mut Rng,
    dist: &D,
    admit_mask: u8,
) -> Result<u8, FilteredSampleError>
where
    D: Distribution<Output = u8> + ?Sized,
{
    let support = dist
        .support()
        .ok_or(FilteredSampleError::SupportUnavailable)?;

    let is_admitted = |base: u8| -> bool {
        let idx = match base {
            b'A' | b'a' => 0,
            b'C' | b'c' => 1,
            b'G' | b'g' => 2,
            b'T' | b't' => 3,
            _ => return false,
        };
        (admit_mask >> idx) & 1 == 1
    };

    // Pass 1: sum admitted weights.
    let mut total = 0.0f64;
    for &(base, w) in support.iter() {
        if is_admitted(base) {
            total += w;
        }
    }
    if total <= 0.0 || !total.is_finite() {
        // total == 0 with admit_mask == 0 (no admissible bases) is
        // the typical case; total non-finite is a defensive guard.
        if admit_mask == 0 {
            return Err(FilteredSampleError::EmptyAdmissibleSupport);
        }
        return Err(FilteredSampleError::InvalidFilteredSupport);
    }

    // Pass 2: inverse CDF over admitted subset.
    let r = rng.next_f64() * total;
    let mut cum = 0.0;
    let mut last_admitted: Option<u8> = None;
    for &(base, w) in support.iter() {
        if is_admitted(base) {
            cum += w;
            if r < cum {
                return Ok(base);
            }
            last_admitted = Some(base);
        }
    }
    // Defensive last-admitted fallback for ULP rounding.
    last_admitted.ok_or(FilteredSampleError::EmptyAdmissibleSupport)
}
