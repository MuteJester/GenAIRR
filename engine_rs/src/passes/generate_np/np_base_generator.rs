//! `NpBaseGenerator` — pass-level abstraction for per-position NP
//! base sampling that admits previous-base conditioning.
//!
//! The slice introduces a single trait separate from
//! [`crate::dist::Distribution`] because PCR / quality / contaminant
//! passes share the stateless `Distribution<Output = u8>` consumer
//! surface — adding an `Option<u8>` previous-base parameter to that
//! trait would force every unrelated consumer to thread a
//! meaningless value. The companion audit
//! [`docs/np_markov_base_generator_design.md`](../../../../../docs/np_markov_base_generator_design.md)
//! §1 documents the reasoning.
//!
//! Three concrete implementations ship in this module:
//!
//! - [`UniformNpGenerator`] — byte-identical wrapper around
//!   [`crate::dist::UniformBase`]; ignores `previous` / `position`.
//! - [`CategoricalNpGenerator`] — byte-identical wrapper around
//!   [`crate::dist::CategoricalBase`]; ignores `previous` /
//!   `position`.
//! - [`MarkovBaseGenerator`] — full 1-step Markov generator with
//!   per-position `first_base` row and 4 transition rows in
//!   canonical A/C/G/T order.
//!
//! Plus an internal [`DistributionNpGenerator`] adapter so existing
//! `GenerateNPPass::new(...)` call sites that pass a
//! `Box<dyn Distribution<Output = u8>>` keep compiling — the
//! adapter delegates to the inner distribution's `support()` and
//! formats `signature()` through the shared `fmt_byte_dist` helper
//! so plan signatures remain byte-identical to the pre-slice
//! behaviour.

use crate::dist::{CategoricalBase, Distribution, UniformBase};
use crate::passes::paramsig::fmt_byte_dist;

/// Pass-level abstraction for per-position NP base sampling.
///
/// Implementors return `(byte, weight)` support pairs for the
/// requested NP position; the caller materialises that support
/// once per position and runs it through the existing
/// admit-mask / `sample_filtered_*` machinery.
///
/// `previous` is the previously emitted NP base. It is `None` only
/// at position 0 (the first NP position). For positions >= 1
/// implementations that don't condition on history (uniform,
/// independent categorical) ignore it; the Markov implementation
/// dispatches on it.
///
/// `position` is the index of the NP slot being drawn (0-indexed,
/// same as the trace address `np.npN.bases[i]`). Implementations
/// that don't position-condition ignore it.
pub trait NpBaseGenerator {
    /// Enumerate the discrete `(byte, weight)` support at this
    /// position. Same shape `Distribution::support()` returns —
    /// the engine's admit-mask / filtered-sample helpers consume
    /// this verbatim.
    ///
    /// Must always return a non-empty vector with at least one
    /// strictly positive weight; the spec layer and the
    /// constructors below enforce this so the per-position
    /// vector is safe to feed into the inverse-CDF samplers.
    fn support(&self, position: usize, previous: Option<u8>) -> Vec<(u8, f64)>;

    /// Deterministic per-cartridge identity string used by
    /// [`crate::passes::generate_np::GenerateNPPass::parameter_signature`].
    ///
    /// **Hard backwards-compat constraint** for the wrapper types:
    /// `UniformNpGenerator::signature()` and
    /// `CategoricalNpGenerator::signature()` MUST return the
    /// exact string that the pre-slice `fmt_byte_dist` produced
    /// over the inner `Distribution`, so legacy trace files
    /// replay byte-identically. The `MarkovBaseGenerator`
    /// signature flattens all 5 rows into one deterministic
    /// string — see the doc on
    /// [`MarkovBaseGenerator::signature`] for the format.
    fn signature(&self) -> String;
}

// ──────────────────────────────────────────────────────────────────
// Wrapper: legacy `Distribution<Output = u8>` adapter
// ──────────────────────────────────────────────────────────────────

/// Adapter that lifts any `Box<dyn Distribution<Output = u8>>`
/// into an `NpBaseGenerator`. Used by `GenerateNPPass::new(...)`
/// to preserve the legacy constructor signature — every call
/// site that today passes `Box::new(UniformBase)` /
/// `Box::new(CategoricalBase::...)` keeps compiling.
///
/// `support()` ignores `position` / `previous` and calls
/// `inner.support()`. The inner distribution must expose
/// enumerable support; the v1 NP base distributions
/// (`UniformBase`, `CategoricalBase`) always do.
///
/// `signature()` delegates to `fmt_byte_dist` so the produced
/// string is byte-identical to the pre-slice
/// `parameter_signature` payload.
pub(crate) struct DistributionNpGenerator {
    inner: Box<dyn Distribution<Output = u8>>,
}

impl DistributionNpGenerator {
    pub(crate) fn new(inner: Box<dyn Distribution<Output = u8>>) -> Self {
        Self { inner }
    }
}

impl NpBaseGenerator for DistributionNpGenerator {
    fn support(&self, _position: usize, _previous: Option<u8>) -> Vec<(u8, f64)> {
        self.inner.support().unwrap_or_default()
    }

    fn signature(&self) -> String {
        fmt_byte_dist(self.inner.as_ref())
    }
}

// ──────────────────────────────────────────────────────────────────
// Concrete wrapper: UniformNpGenerator
// ──────────────────────────────────────────────────────────────────

/// `NpBaseGenerator` wrapper around the canonical 4-way uniform
/// distribution `{b'A', b'C', b'G', b'T'}` (each weight `1.0`).
///
/// Ignores `previous` and `position`. Produces a plan-signature
/// substring byte-identical to the pre-slice `UniformBase`
/// rendering through `fmt_byte_dist`:
/// `[(65:1.0),(67:1.0),(71:1.0),(84:1.0)]`.
pub struct UniformNpGenerator;

impl UniformNpGenerator {
    pub fn new() -> Self {
        Self
    }
}

impl Default for UniformNpGenerator {
    fn default() -> Self {
        Self::new()
    }
}

impl NpBaseGenerator for UniformNpGenerator {
    fn support(&self, _position: usize, _previous: Option<u8>) -> Vec<(u8, f64)> {
        vec![(b'A', 1.0), (b'C', 1.0), (b'G', 1.0), (b'T', 1.0)]
    }

    fn signature(&self) -> String {
        // Hard-coded to match `fmt_byte_dist(UniformBase)` for
        // byte-identical legacy plan signatures.
        fmt_byte_dist(&UniformBase as &dyn Distribution<Output = u8>)
    }
}

// ──────────────────────────────────────────────────────────────────
// Concrete wrapper: CategoricalNpGenerator
// ──────────────────────────────────────────────────────────────────

/// `NpBaseGenerator` wrapper around a position-independent
/// `CategoricalBase` (the existing `empirical_first_base`
/// surface). Ignores `previous` and `position` — every NP slot
/// samples from the same weighted A/C/G/T categorical.
///
/// Plan signature is byte-identical to the pre-slice
/// `CategoricalBase` rendering through `fmt_byte_dist`.
pub struct CategoricalNpGenerator {
    inner: CategoricalBase,
}

impl CategoricalNpGenerator {
    pub fn new(inner: CategoricalBase) -> Self {
        Self { inner }
    }

    pub fn from_pairs(pairs: Vec<(u8, f64)>) -> Self {
        Self {
            inner: CategoricalBase::from_pairs(pairs),
        }
    }
}

impl NpBaseGenerator for CategoricalNpGenerator {
    fn support(&self, _position: usize, _previous: Option<u8>) -> Vec<(u8, f64)> {
        self.inner.support().expect("CategoricalBase always exposes support")
    }

    fn signature(&self) -> String {
        fmt_byte_dist(&self.inner as &dyn Distribution<Output = u8>)
    }
}

// ──────────────────────────────────────────────────────────────────
// Concrete: MarkovBaseGenerator
// ──────────────────────────────────────────────────────────────────

/// Index in canonical A/C/G/T order. Returns `None` for any byte
/// outside the canonical alphabet — used by Markov to fall back
/// to the first-base row when `previous` is a mid-stream sentinel
/// (`b'N'`) or otherwise non-canonical.
#[inline]
fn canonical_base_index(b: u8) -> Option<usize> {
    match b {
        b'A' | b'a' => Some(0),
        b'C' | b'c' => Some(1),
        b'G' | b'g' => Some(2),
        b'T' | b't' => Some(3),
        _ => None,
    }
}

const CANONICAL_BASES: [u8; 4] = [b'A', b'C', b'G', b'T'];

/// 1-step Markov generator: position 0 samples from
/// `first_base`; positions 1+ sample from the transition row
/// keyed by the previously emitted base.
///
/// Storage shape:
///
/// - `first_base: [f64; 4]` — A/C/G/T weights for position 0.
/// - `transitions: [[f64; 4]; 4]` — row-major `from × to`
///   weights. Row `i` is the categorical for the next base when
///   the previous base was `CANONICAL_BASES[i]`.
///
/// Per-row invariants (validated at construction):
///
/// - Each weight is finite and non-negative.
/// - Each row has at least one strictly positive weight (i.e.
///   total > 0).
///
/// **Mid-stream `b'N'` fallback.** If a previously emitted base
/// is the `NP_BASE_EMPTY_SUPPORT` sentinel (`b'N'`) — emitted
/// when the admit mask × generator support intersection was
/// empty at the prior position — `support(position, previous)`
/// falls back to `first_base`. This is the recommended fallback
/// from the audit §8: mid-stream `b'N'` was already a sentinel,
/// reverting to `first_base` for the next position is the
/// cleanest non-failure path.
pub struct MarkovBaseGenerator {
    first_base: [f64; 4],
    transitions: [[f64; 4]; 4],
}

impl MarkovBaseGenerator {
    /// Construct from the two raw weight arrays in canonical
    /// A/C/G/T order. Defensive validation re-runs the
    /// finite/non-negative/row-positive checks the Python spec
    /// layer already enforces — custom Rust callers that bypass
    /// the spec layer still get a clear panic.
    pub fn from_arrays(first_base: [f64; 4], transitions: [[f64; 4]; 4]) -> Self {
        for (i, w) in first_base.iter().enumerate() {
            assert!(
                w.is_finite() && *w >= 0.0,
                "MarkovBaseGenerator.first_base[{i}] must be finite and non-negative, got {w}",
            );
        }
        let total: f64 = first_base.iter().sum();
        assert!(
            total > 0.0,
            "MarkovBaseGenerator.first_base must have at least one strictly positive weight; total was {total}",
        );
        for (i, row) in transitions.iter().enumerate() {
            for (j, w) in row.iter().enumerate() {
                assert!(
                    w.is_finite() && *w >= 0.0,
                    "MarkovBaseGenerator.transitions[{i}][{j}] must be finite and non-negative, got {w}",
                );
            }
            let row_total: f64 = row.iter().sum();
            assert!(
                row_total > 0.0,
                "MarkovBaseGenerator.transitions[{i}] must have at least one strictly positive weight; total was {row_total}",
            );
        }
        Self {
            first_base,
            transitions,
        }
    }

    /// Construct from a row-list representation suitable for the
    /// PyO3 bridge. `transitions_rows` is the four-row matrix in
    /// canonical A/C/G/T from-order; each row is a 4-element
    /// `[f64; 4]` in canonical A/C/G/T to-order.
    pub fn from_first_and_rows(
        first_base: Vec<f64>,
        transitions_rows: Vec<Vec<f64>>,
    ) -> Self {
        assert_eq!(
            first_base.len(),
            4,
            "MarkovBaseGenerator: first_base must be a 4-element vector in A/C/G/T order"
        );
        assert_eq!(
            transitions_rows.len(),
            4,
            "MarkovBaseGenerator: transitions must have 4 rows (A/C/G/T from-bases)"
        );
        let mut fb = [0.0f64; 4];
        for i in 0..4 {
            fb[i] = first_base[i];
        }
        let mut tx = [[0.0f64; 4]; 4];
        for (i, row) in transitions_rows.into_iter().enumerate() {
            assert_eq!(
                row.len(),
                4,
                "MarkovBaseGenerator: transitions row {i} must be a 4-element vector in A/C/G/T order"
            );
            for (j, w) in row.into_iter().enumerate() {
                tx[i][j] = w;
            }
        }
        Self::from_arrays(fb, tx)
    }
}

impl NpBaseGenerator for MarkovBaseGenerator {
    fn support(&self, position: usize, previous: Option<u8>) -> Vec<(u8, f64)> {
        // Position 0 OR `previous` is a non-canonical byte
        // (typically the `b'N'` mid-stream sentinel) → fall back
        // to the first-base row. The latter case is the
        // audit-§8 recommended fallback: mid-stream `b'N'` is
        // already a sentinel; reverting to `first_base` is the
        // cleanest non-failure path that keeps replay
        // deterministic without needing a Markov state machine
        // that tracks failure history.
        let row = if position == 0 {
            self.first_base
        } else {
            match previous.and_then(canonical_base_index) {
                Some(i) => self.transitions[i],
                None => self.first_base,
            }
        };
        CANONICAL_BASES
            .iter()
            .zip(row.iter())
            .map(|(&b, &w)| (b, w))
            .collect()
    }

    /// Flat deterministic encoding covering first-base + all
    /// four transition rows in canonical A/C/G/T order. Format:
    ///
    /// ```text
    /// markov:first=[(65:w_A),(67:w_C),(71:w_G),(84:w_T)]
    ///   |from=A:[(65:w_AA),(67:w_AC),(71:w_AG),(84:w_AT)]
    ///   |from=C:[…]|from=G:[…]|from=T:[…]
    /// ```
    ///
    /// Row separators are `|`; the `from=X:` prefix tags each
    /// transition row by from-base letter. Inside each
    /// `[(v:w),…]` list weights use Rust's stable f64 Debug
    /// rendering, identical to `fmt_byte_dist`'s discipline.
    fn signature(&self) -> String {
        fn fmt_row(weights: &[f64; 4]) -> String {
            let mut s = String::from("[");
            for (i, (&b, &w)) in CANONICAL_BASES.iter().zip(weights.iter()).enumerate() {
                if i > 0 {
                    s.push(',');
                }
                s.push_str(&format!("({}:{:?})", b, w));
            }
            s.push(']');
            s
        }
        let mut out = String::from("markov:first=");
        out.push_str(&fmt_row(&self.first_base));
        for (i, &from) in CANONICAL_BASES.iter().enumerate() {
            out.push_str("|from=");
            out.push(from as char);
            out.push(':');
            out.push_str(&fmt_row(&self.transitions[i]));
        }
        out
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn uniform_np_generator_support_is_canonical_quartet() {
        let g = UniformNpGenerator::new();
        let s = g.support(0, None);
        assert_eq!(s, vec![(b'A', 1.0), (b'C', 1.0), (b'G', 1.0), (b'T', 1.0)]);
        // Ignores position / previous.
        let s2 = g.support(5, Some(b'C'));
        assert_eq!(s2, s);
    }

    #[test]
    fn uniform_np_generator_signature_matches_uniform_base_substring() {
        let g = UniformNpGenerator::new();
        assert_eq!(
            g.signature(),
            "[(65:1.0),(67:1.0),(71:1.0),(84:1.0)]"
        );
    }

    #[test]
    fn categorical_np_generator_delegates_support_to_inner() {
        let cat = CategoricalBase::from_weights([1.0, 2.0, 3.0, 4.0]);
        let g = CategoricalNpGenerator::new(cat);
        let s = g.support(0, None);
        assert_eq!(s, vec![(b'A', 1.0), (b'C', 2.0), (b'G', 3.0), (b'T', 4.0)]);
    }

    #[test]
    fn markov_generator_position_zero_returns_first_base() {
        let g = MarkovBaseGenerator::from_arrays(
            [0.1, 0.2, 0.3, 0.4],
            [[1.0; 4]; 4],
        );
        let s = g.support(0, None);
        assert_eq!(s, vec![(b'A', 0.1), (b'C', 0.2), (b'G', 0.3), (b'T', 0.4)]);
    }

    #[test]
    fn markov_generator_position_one_routes_through_previous() {
        // A row says A→100% T, so support keyed on previous=A
        // returns weights (0,0,0,1).
        let g = MarkovBaseGenerator::from_arrays(
            [1.0; 4],
            [
                [0.0, 0.0, 0.0, 1.0], // from A
                [1.0, 1.0, 1.0, 1.0], // from C
                [1.0, 1.0, 1.0, 1.0], // from G
                [1.0, 1.0, 1.0, 1.0], // from T
            ],
        );
        let s = g.support(1, Some(b'A'));
        assert_eq!(s, vec![(b'A', 0.0), (b'C', 0.0), (b'G', 0.0), (b'T', 1.0)]);
    }

    #[test]
    fn markov_generator_mid_stream_sentinel_falls_back_to_first_base() {
        let g = MarkovBaseGenerator::from_arrays(
            [0.5, 0.5, 0.0, 0.0],
            [[0.0, 0.0, 0.0, 1.0]; 4], // every row says next = T
        );
        // `previous = b'N'` should NOT match a transition row;
        // fall back to first_base.
        let s = g.support(3, Some(b'N'));
        assert_eq!(s, vec![(b'A', 0.5), (b'C', 0.5), (b'G', 0.0), (b'T', 0.0)]);
    }

    #[test]
    #[should_panic(expected = "at least one strictly positive weight")]
    fn markov_generator_rejects_zero_row() {
        let _ = MarkovBaseGenerator::from_arrays(
            [1.0; 4],
            [[0.0; 4], [1.0; 4], [1.0; 4], [1.0; 4]],
        );
    }

    #[test]
    fn markov_signature_canonical_a_c_g_t_order() {
        // The format must be byte-deterministic. We pin the
        // exact substrings here so a refactor that reorders
        // rows surfaces immediately.
        let g = MarkovBaseGenerator::from_arrays(
            [1.0, 2.0, 3.0, 4.0],
            [
                [0.1, 0.2, 0.3, 0.4],
                [0.5, 0.6, 0.7, 0.8],
                [1.1, 1.2, 1.3, 1.4],
                [1.5, 1.6, 1.7, 1.8],
            ],
        );
        let sig = g.signature();
        assert!(sig.starts_with("markov:first=[(65:1.0),(67:2.0),(71:3.0),(84:4.0)]"));
        assert!(sig.contains("|from=A:[(65:0.1),(67:0.2),(71:0.3),(84:0.4)]"));
        assert!(sig.contains("|from=C:[(65:0.5),(67:0.6),(71:0.7),(84:0.8)]"));
        assert!(sig.contains("|from=G:[(65:1.1),(67:1.2),(71:1.3),(84:1.4)]"));
        assert!(sig.ends_with("|from=T:[(65:1.5),(67:1.6),(71:1.7),(84:1.8)]"));
    }

    #[test]
    fn markov_signature_differs_for_different_matrices() {
        let g1 = MarkovBaseGenerator::from_arrays([1.0; 4], [[1.0; 4]; 4]);
        let g2 = MarkovBaseGenerator::from_arrays(
            [1.0; 4],
            [
                [0.0, 0.0, 0.0, 1.0],
                [1.0, 1.0, 1.0, 1.0],
                [1.0, 1.0, 1.0, 1.0],
                [1.0, 1.0, 1.0, 1.0],
            ],
        );
        assert_ne!(g1.signature(), g2.signature());
    }
}
