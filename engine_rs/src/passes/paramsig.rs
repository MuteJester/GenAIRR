//! Shared formatting helpers for [`Pass::parameter_signature`].
//!
//! Every pass that takes compile-time parameters formats them
//! through these helpers so the resulting string is:
//!
//! - **Deterministic** — same parameter struct in, same string
//!   out, across runs and platforms (modulo f64 formatting which
//!   is itself stable for finite normal values).
//! - **Canonical for behaviourally-equivalent inputs** — a
//!   default rate vector and an explicit all-ones rate vector
//!   produce equal strings; an empirical distribution with the
//!   same `(value, weight)` pairs in the same insertion order
//!   produces equal strings regardless of how the caller built
//!   it.
//! - **Concise** — the string is folded into the plan signature
//!   that ships in every trace file; keeping it tight matters.
//!
//! Format convention used throughout: lowercase comma-separated
//! `key=value` pairs inside the `pass_name(...)` envelope.
//! Values use a stable rendering for each domain type:
//!
//! - `f64`              → `format!("{:?}", x)` (Rust's `Debug`
//!                        rendering is stable for normal finite
//!                        values: `1.0`, `0.03`, `2.5e-5`, …).
//! - `i64`              → `format!("{}", x)`.
//! - `u32` / `u8`       → `format!("{}", x)`.
//! - distributions over `i64` (NP length, trim, end-loss length,
//!   paired-end r1/r2/insert) → `[(v:w),(v:w),…]` from
//!   `Distribution::support()`.
//! - distributions over `u8` (base draws) → `[(v:w),(v:w),…]`
//!   with `v` as the ASCII byte (the actual u8 nucleotide code).
//!
//! When a distribution returns `None` from `support()` the helper
//! emits the sentinel `"opaque"` — that path is only reachable
//! for user-defined distributions outside the built-in vocabulary
//! (none of the production passes use it).

use crate::dist::Distribution;
use crate::passes::mutate::CountSource;

/// Render a per-segment SHM rate vector for the plan signature.
///
/// The default vector (all 1.0) — whether arrived at via
/// `SegmentRateWeights::default()`, an omitted `segment_rates`
/// kwarg, or an explicit `{"V": 1.0, "D": 1.0, "J": 1.0, "NP":
/// 1.0}` dict — short-circuits to the empty string. This is
/// what makes the default-vs-explicit-all-ones equivalence test
/// pass at the signature layer.
///
/// Non-default vectors render as
/// `segment_rates=(v:Vd:Dj:Jnp:NP)` where each token uses
/// Rust's stable `Debug` formatting for f64. Example:
/// `segment_rates=(v:2.0,d:1.0,j:1.0,np:0.0)`.
pub fn fmt_segment_rates(rates: &crate::passes::mutate::SegmentRateWeights) -> String {
    if rates.is_default() {
        return String::new();
    }
    format!(
        "segment_rates=(v:{:?},d:{:?},j:{:?},np:{:?})",
        rates.v, rates.d, rates.j, rates.np
    )
}

/// Render a per-V-subregion SHM rate vector for the plan signature.
///
/// Same short-circuit discipline as [`fmt_segment_rates`]: the
/// flat default (all five labels at `1.0`) renders as the empty
/// string so omitting the kwarg or passing
/// `{"FWR1": 1.0, "CDR1": 1.0, "FWR2": 1.0, "CDR2": 1.0,
/// "FWR3": 1.0}` collide on equal signatures. Non-default
/// vectors render as
/// `v_subregion_rates=(fwr1:FWR1,cdr1:CDR1,fwr2:FWR2,cdr2:CDR2,fwr3:FWR3)`.
pub fn fmt_v_subregion_rates(
    rates: &crate::passes::mutate::VSubregionRateWeights,
) -> String {
    if rates.is_default() {
        return String::new();
    }
    format!(
        "v_subregion_rates=(fwr1:{:?},cdr1:{:?},fwr2:{:?},cdr2:{:?},fwr3:{:?})",
        rates.fwr1, rates.cdr1, rates.fwr2, rates.cdr2, rates.fwr3
    )
}

/// Render a [`CountSource`] for the plan signature.
///
/// `Rate(x)` → `count=rate:x` with `x` in Rust's f64 Debug
/// rendering. `Distribution(d)` → `count=dist:[(v:w),(v:w),…]`
/// via [`fmt_int_dist`]. `Distribution` variants over an
/// opaque (unsupported) inner type fall through to
/// `count=dist:opaque`.
pub fn fmt_count_source(src: &CountSource) -> String {
    match src {
        CountSource::Rate(r) => format!("count=rate:{:?}", r),
        CountSource::Distribution(d) => format!("count=dist:{}", fmt_int_dist(d.as_ref())),
    }
}

/// Render a boxed `Distribution<Output = i64>` (NP length, trim,
/// end-loss length, paired-end window) via its
/// [`Distribution::support`] vector.
///
/// Each `(value, weight)` pair is emitted as `(v:w)` with
/// integers formatted as decimals and weights in f64 Debug
/// rendering. The list preserves insertion order — the
/// distribution's authoring discipline must keep that
/// deterministic per construction site.
pub fn fmt_int_dist(d: &dyn Distribution<Output = i64>) -> String {
    match d.support() {
        Some(pairs) => fmt_int_pairs(&pairs),
        None => "opaque".to_string(),
    }
}

/// Render an explicit `(i64, f64)` support vector (e.g.
/// reconstructed from PyO3 `count_pairs`). Mirrors the
/// `fmt_int_dist` payload format so manually-formatted and
/// distribution-derived strings collide on equal inputs.
pub fn fmt_int_pairs(pairs: &[(i64, f64)]) -> String {
    let mut s = String::from("[");
    for (i, (v, w)) in pairs.iter().enumerate() {
        if i > 0 {
            s.push(',');
        }
        s.push_str(&format!("({}:{:?})", v, w));
    }
    s.push(']');
    s
}

/// Render a boxed `Distribution<Output = u8>` (base draws) via
/// its [`Distribution::support`] vector. The u8 is emitted as
/// its decimal integer value (the actual byte) so two
/// behaviourally-equivalent base distributions collide on equal
/// inputs.
pub fn fmt_byte_dist(d: &dyn Distribution<Output = u8>) -> String {
    match d.support() {
        Some(pairs) => {
            let mut s = String::from("[");
            for (i, (v, w)) in pairs.iter().enumerate() {
                if i > 0 {
                    s.push(',');
                }
                s.push_str(&format!("({}:{:?})", v, w));
            }
            s.push(']');
            s
        }
        None => "opaque".to_string(),
    }
}

/// Render an f64 probability scalar (Bernoulli probabilities,
/// insertion fractions, …) under the `key=value` convention.
/// Identity rendering via Rust's f64 Debug; equal values render
/// equal strings.
pub fn fmt_prob(key: &str, value: f64) -> String {
    format!("{}={:?}", key, value)
}

/// Join a list of `key=value` parts into the canonical
/// `parameter_signature` body, dropping empty entries so a
/// short-circuiting helper (e.g. [`fmt_segment_rates`] under
/// the default vector) doesn't introduce a trailing comma.
pub fn join_parts<I: IntoIterator<Item = String>>(parts: I) -> String {
    let collected: Vec<String> = parts.into_iter().filter(|p| !p.is_empty()).collect();
    collected.join(",")
}
