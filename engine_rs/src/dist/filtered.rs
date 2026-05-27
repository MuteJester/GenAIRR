//! Filtered sampling — draw from a distribution restricted to a predicate.
//!
//! The architectural primitive for constraint-aware sampling: a sampling
//! pass calls [`sample_filtered_result`] with a `predicate` wrapping a
//! `ContractSet::admits` check, and the pass either receives an
//! admissible value or a structured error explaining why no admissible
//! candidate exists.
//!
//! ## Empty-support policy (v3.0)
//!
//! Strict mode always surfaces a structured `PassError` when the
//! filtered support is empty (or unenumerable / invalid). The
//! interesting design choice is *permissive* mode: the v3.0 rule
//! says the engine must NOT silently fall back to an unconstrained
//! draw — that would reintroduce reject-after-propose at the
//! per-event level. Each pass declares its policy via
//! [`EmptySupport`]:
//!
//! - [`EmptySupport::Skip`] — the slot is consumed as a no-op. No
//!   trace record, no pool mutation. Used by passes whose "skip"
//!   semantics is well-defined: per-site substitution, indel
//!   tuples that can't be balanced, contaminant slots.
//! - [`EmptySupport::Sentinel`] — a known-safe sentinel is
//!   written. Used by length samplers (NP/Trim default to `0`)
//!   and the NP base sampler (defaults to `b'N'`, the IUPAC
//!   ambiguous nucleotide).
//!
//! `sample_filtered_result` returns the raw `Result<T, _>`;
//! [`sample_filtered_with_policy`] wraps it with the strict /
//! permissive policy resolution so each pass becomes one line at
//! the call site instead of a repeated `match` block.

use super::Distribution;
use crate::pass::PassError;
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

/// Per-pass policy for what to return from a constrained sampler
/// when the filtered support is empty (in permissive mode).
///
/// **Strict mode is independent of this policy**: when the filter
/// is empty under strict, `sample_filtered_with_policy` always
/// returns `Err(PassError::ConstraintSampling)`. The policy only
/// governs permissive-mode behavior, where the v3.0 rule forbids
/// the silent unconstrained-fallback that earlier versions used.
///
/// Choose the variant that matches each pass's well-defined
/// "no-op": substitution passes use [`Skip`](Self::Skip) (no
/// trace record, no pool mutation); length samplers and the NP
/// base path use [`Sentinel`](Self::Sentinel) with a known-safe
/// stand-in value (`0` for length, `b'N'` for NP base).
#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub enum EmptySupport<T> {
    /// The call site treats permissive empty-support as a no-op.
    /// `sample_filtered_with_policy` returns `Ok(None)` and the
    /// caller skips the slot.
    Skip,
    /// The call site writes the carried sentinel value in
    /// permissive empty-support.
    /// `sample_filtered_with_policy` returns `Ok(Some(sentinel))`.
    Sentinel(T),
}

/// Filtered sample + per-pass empty-support policy resolution.
///
/// On success: `Ok(Some(value))` for an admissible draw, or
/// `Ok(Some(sentinel))` when permissive mode hit empty support
/// with [`EmptySupport::Sentinel`]. `Ok(None)` is returned for
/// permissive [`EmptySupport::Skip`] on empty support.
///
/// **Error semantics differ by failure cause:**
/// - [`FilteredSampleError::EmptyAdmissibleSupport`] and
///   [`FilteredSampleError::SupportUnavailable`]: legitimate
///   no-admissible-candidate states for the active contract
///   bundle. Strict mode surfaces
///   `PassError::ConstraintSampling`; permissive mode applies
///   the caller's [`EmptySupport`] policy.
/// - [`FilteredSampleError::InvalidFilteredSupport`] (non-finite
///   total weight, ≤ 0 total): this signals a **distribution
///   bug** — the natural distribution itself has corrupt
///   weights. Surfaces `PassError::ConstraintSampling` in
///   **both** strict and permissive modes; the policy is not
///   applied. Silently degrading to a sentinel/no-op here would
///   mask data-pipeline corruption.
///
/// Centralizes the v3.0 invariant ("engine never proposes
/// contract-violating actions on empty support") so each pass
/// declares its policy once at the call site instead of
/// reimplementing the strict/permissive match.
pub fn sample_filtered_with_policy<T, D, F>(
    rng: &mut Rng,
    dist: &D,
    predicate: F,
    strict: bool,
    pass_name: &str,
    address: &str,
    policy: EmptySupport<T>,
) -> Result<Option<T>, PassError>
where
    T: Clone,
    D: Distribution<Output = T> + ?Sized,
    F: Fn(&T) -> bool,
{
    match sample_filtered_result(rng, dist, predicate) {
        Ok(value) => Ok(Some(value)),
        Err(FilteredSampleError::InvalidFilteredSupport) => {
            // Distribution bug — strict and permissive both
            // refuse to silently mask it. The pass surfaces the
            // structured `ConstraintSampling` error with the
            // `InvalidFilteredSupport` reason so the caller can
            // diagnose the natural distribution.
            Err(PassError::constraint_sampling(
                pass_name.to_string(),
                address.to_string(),
                FilteredSampleError::InvalidFilteredSupport,
            ))
        }
        Err(reason) if strict => Err(PassError::constraint_sampling(
            pass_name.to_string(),
            address.to_string(),
            reason,
        )),
        Err(_) => match policy {
            EmptySupport::Skip => Ok(None),
            EmptySupport::Sentinel(s) => Ok(Some(s)),
        },
    }
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
