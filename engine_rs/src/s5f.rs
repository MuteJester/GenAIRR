//! S5F mutation kernel — context-dependent substitution model for SHM.
//!
//! ## What this is
//!
//! The S5F (Shazam 5-mer Substitution Format) model from Yaari et
//! al. 2013 represents somatic hypermutation as two independent
//! tables indexed by 5-mer context:
//!
//! 1. **Mutability** — for each of the 4⁵ = 1024 possible 5-mer
//!    contexts, the relative likelihood that the central (3rd)
//!    base mutates.
//! 2. **Substitution** — for each context + destination base, the
//!    probability that a mutation at the central base produces
//!    that destination. 1024 × 4 = 4096 entries.
//!
//! The full model is: at each position with a defined 5-mer context
//! (i.e., not within 2 bases of the sequence end and no `N`s in the
//! context), draw a per-position Bernoulli or sample the mutability
//! distribution to decide whether to mutate, then sample a
//! destination base from the per-context substitution distribution.
//!
//! Reference: Yaari, G. et al. (2013) "Models of somatic hypermutation
//! targeting and substitution based on synonymous mutations from
//! high-throughput immunoglobulin sequencing data." Front. Immunol.
//!
//! ## Phase E.2 scope
//!
//! Just the data structure plus accessors and context encoding.
//! Phase E.3 builds the actual `S5FMutationPass` on top.

use crate::ir::encode_base;

// ──────────────────────────────────────────────────────────────────
// S5FKernel — the two empirical tables
// ──────────────────────────────────────────────────────────────────

/// Number of distinct 5-mer contexts (4 bases × 5 positions).
pub const S5F_NUM_CONTEXTS: usize = 1024;

/// Length of the substitution table (1024 contexts × 4 dest bases).
pub const S5F_SUBSTITUTION_LEN: usize = S5F_NUM_CONTEXTS * 4;

/// The S5F kernel: mutability scores + substitution probabilities
/// keyed by 5-mer context.
///
/// Context encoding: a 5-mer `b1 b2 b3 b4 b5` maps to an index
/// `b1 << 8 | b2 << 6 | b3 << 4 | b4 << 2 | b5` with A=0, C=1,
/// G=2, T=3 (matching `ir::encode_base`). Only A/C/G/T are valid;
/// any context with `N` / gap / ambiguity codes returns `None`
/// from `encode_context`.
///
/// Substitution table is laid out as a flat `Vec<f64>` of length
/// `1024 * 4`; the entry for `(context, dest_base)` is at offset
/// `context * 4 + dest_idx` where `dest_idx` follows the same
/// A=0/C=1/G=2/T=3 encoding.
#[derive(Clone, Debug)]
pub struct S5FKernel {
    mutability: Vec<f64>,
    substitution: Vec<f64>,
}

impl S5FKernel {
    /// Construct a kernel from raw mutability + substitution tables.
    ///
    /// Panics if:
    /// - `mutability.len() != 1024`,
    /// - `substitution.len() != 4096`,
    /// - any value is non-finite or negative.
    ///
    /// Validation is at construction time so accessors stay
    /// panic-free in the hot path.
    pub fn new(mutability: Vec<f64>, substitution: Vec<f64>) -> Self {
        assert_eq!(
            mutability.len(),
            S5F_NUM_CONTEXTS,
            "S5FKernel: mutability table must have {} entries (got {})",
            S5F_NUM_CONTEXTS,
            mutability.len()
        );
        assert_eq!(
            substitution.len(),
            S5F_SUBSTITUTION_LEN,
            "S5FKernel: substitution table must have {} entries (got {})",
            S5F_SUBSTITUTION_LEN,
            substitution.len()
        );

        for (i, &m) in mutability.iter().enumerate() {
            assert!(
                m.is_finite() && m >= 0.0,
                "S5FKernel: mutability[{}] must be a finite non-negative number, got {}",
                i,
                m
            );
        }
        for (i, &s) in substitution.iter().enumerate() {
            assert!(
                s.is_finite() && s >= 0.0,
                "S5FKernel: substitution[{}] must be a finite non-negative number, got {}",
                i,
                s
            );
        }

        Self {
            mutability,
            substitution,
        }
    }

    /// Encode a 5-mer (`b1 b2 b3 b4 b5`) into a context index in
    /// `[0, 1024)`. Returns `None` if any base is not A/C/G/T (or
    /// U, which is treated as T per `ir::encode_base`).
    pub fn encode_context(b1: u8, b2: u8, b3: u8, b4: u8, b5: u8) -> Option<u16> {
        let i1 = encode_base(b1)?;
        let i2 = encode_base(b2)?;
        let i3 = encode_base(b3)?;
        let i4 = encode_base(b4)?;
        let i5 = encode_base(b5)?;
        Some(
            ((i1 as u16) << 8)
                | ((i2 as u16) << 6)
                | ((i3 as u16) << 4)
                | ((i4 as u16) << 2)
                | (i5 as u16),
        )
    }

    /// Mutability for a given context. Caller must pass a valid
    /// context (`< S5F_NUM_CONTEXTS`); construction-time validation
    /// guarantees the value at that index is finite and non-negative.
    pub fn mutability(&self, context: u16) -> f64 {
        self.mutability[context as usize]
    }

    /// Substitution probability for a given context and destination
    /// base. Returns 0.0 if `dest_base` is not A/C/G/T (defensive).
    pub fn substitution(&self, context: u16, dest_base: u8) -> f64 {
        let dest_idx = match encode_base(dest_base) {
            Some(i) => i as usize,
            None => return 0.0,
        };
        self.substitution[context as usize * 4 + dest_idx]
    }

    /// All four substitution probabilities for a given context, in
    /// A/C/G/T order. Useful when building a per-context categorical
    /// distribution for the destination base.
    pub fn substitution_row(&self, context: u16) -> [f64; 4] {
        let base = context as usize * 4;
        [
            self.substitution[base],
            self.substitution[base + 1],
            self.substitution[base + 2],
            self.substitution[base + 3],
        ]
    }

    /// Read-only view of the entire mutability table. Useful for
    /// serialization or aggregate analysis.
    pub fn mutability_table(&self) -> &[f64] {
        &self.mutability
    }

    /// Read-only view of the entire substitution table.
    pub fn substitution_table(&self) -> &[f64] {
        &self.substitution
    }
}

// ──────────────────────────────────────────────────────────────────
// Tests
// ──────────────────────────────────────────────────────────────────

#[cfg(test)]
mod tests {
    use super::*;

    /// Helper: a kernel with all mutabilities = 1.0 and uniform
    /// substitution (0.25 across all destinations) for testing.
    fn make_uniform_kernel() -> S5FKernel {
        S5FKernel::new(
            vec![1.0; S5F_NUM_CONTEXTS],
            vec![0.25; S5F_SUBSTITUTION_LEN],
        )
    }

    #[test]
    #[should_panic(expected = "mutability table must have 1024 entries")]
    fn s5f_kernel_rejects_wrong_mutability_length() {
        let _ = S5FKernel::new(vec![0.5; 100], vec![0.25; S5F_SUBSTITUTION_LEN]);
    }

    #[test]
    #[should_panic(expected = "substitution table must have 4096 entries")]
    fn s5f_kernel_rejects_wrong_substitution_length() {
        let _ = S5FKernel::new(vec![0.5; S5F_NUM_CONTEXTS], vec![0.25; 100]);
    }

    #[test]
    #[should_panic(expected = "mutability[")]
    fn s5f_kernel_rejects_negative_mutability() {
        let mut m = vec![0.5; S5F_NUM_CONTEXTS];
        m[42] = -0.1;
        let _ = S5FKernel::new(m, vec![0.25; S5F_SUBSTITUTION_LEN]);
    }

    #[test]
    #[should_panic(expected = "mutability[")]
    fn s5f_kernel_rejects_nan_mutability() {
        let mut m = vec![0.5; S5F_NUM_CONTEXTS];
        m[5] = f64::NAN;
        let _ = S5FKernel::new(m, vec![0.25; S5F_SUBSTITUTION_LEN]);
    }

    #[test]
    #[should_panic(expected = "substitution[")]
    fn s5f_kernel_rejects_negative_substitution() {
        let mut s = vec![0.25; S5F_SUBSTITUTION_LEN];
        s[10] = -1.0;
        let _ = S5FKernel::new(vec![0.5; S5F_NUM_CONTEXTS], s);
    }

    #[test]
    #[should_panic(expected = "substitution[")]
    fn s5f_kernel_rejects_infinite_substitution() {
        let mut s = vec![0.25; S5F_SUBSTITUTION_LEN];
        s[0] = f64::INFINITY;
        let _ = S5FKernel::new(vec![0.5; S5F_NUM_CONTEXTS], s);
    }

    #[test]
    fn s5f_kernel_zero_values_allowed() {
        // Zero is fine — represents "this context never mutates" or
        // "this dest base is never the result of mutating this context."
        let kernel = S5FKernel::new(vec![0.0; S5F_NUM_CONTEXTS], vec![0.0; S5F_SUBSTITUTION_LEN]);
        assert_eq!(kernel.mutability(0), 0.0);
    }

    #[test]
    fn s5f_kernel_construction_round_trip_accessors() {
        let kernel = make_uniform_kernel();
        assert_eq!(kernel.mutability(0), 1.0);
        assert_eq!(kernel.mutability(1023), 1.0);
        assert_eq!(kernel.substitution(0, b'A'), 0.25);
        assert_eq!(kernel.substitution(0, b'T'), 0.25);
        assert_eq!(kernel.substitution_row(0), [0.25, 0.25, 0.25, 0.25]);
        assert_eq!(kernel.mutability_table().len(), S5F_NUM_CONTEXTS);
        assert_eq!(kernel.substitution_table().len(), S5F_SUBSTITUTION_LEN);
    }

    #[test]
    fn s5f_kernel_substitution_returns_zero_for_invalid_dest() {
        let kernel = make_uniform_kernel();
        assert_eq!(kernel.substitution(0, b'N'), 0.0);
        assert_eq!(kernel.substitution(0, b'-'), 0.0);
        assert_eq!(kernel.substitution(0, b'.'), 0.0);
    }

    #[test]
    fn s5f_kernel_substitution_is_case_insensitive() {
        let kernel = make_uniform_kernel();
        // Lowercase input should be accepted (ir::encode_base is
        // case-insensitive); same value as uppercase.
        assert_eq!(kernel.substitution(0, b'a'), 0.25);
        assert_eq!(kernel.substitution(0, b'c'), 0.25);
        assert_eq!(kernel.substitution(0, b'g'), 0.25);
        assert_eq!(kernel.substitution(0, b't'), 0.25);
    }

    #[test]
    fn s5f_kernel_encode_context_canonical_5mers() {
        // AAAAA = 0
        assert_eq!(
            S5FKernel::encode_context(b'A', b'A', b'A', b'A', b'A'),
            Some(0)
        );
        // AAAAC = 1
        assert_eq!(
            S5FKernel::encode_context(b'A', b'A', b'A', b'A', b'C'),
            Some(1)
        );
        // AAAAT = 3
        assert_eq!(
            S5FKernel::encode_context(b'A', b'A', b'A', b'A', b'T'),
            Some(3)
        );
        // AAACA = 4
        assert_eq!(
            S5FKernel::encode_context(b'A', b'A', b'A', b'C', b'A'),
            Some(4)
        );
        // TTTTT = 1023 (4^5 - 1)
        assert_eq!(
            S5FKernel::encode_context(b'T', b'T', b'T', b'T', b'T'),
            Some(1023)
        );
        // CACGT — verify by hand: C=1, A=0, C=1, G=2, T=3 → 1<<8 | 0<<6 | 1<<4 | 2<<2 | 3 = 256+16+8+3 = 283
        assert_eq!(
            S5FKernel::encode_context(b'C', b'A', b'C', b'G', b'T'),
            Some(283)
        );
    }

    #[test]
    fn s5f_kernel_encode_context_handles_lowercase_and_u() {
        // U treated as T per ir::encode_base.
        assert_eq!(
            S5FKernel::encode_context(b'a', b'a', b'a', b'a', b'a'),
            Some(0)
        );
        assert_eq!(
            S5FKernel::encode_context(b'T', b'T', b'T', b'T', b'T'),
            S5FKernel::encode_context(b'U', b'U', b'U', b'U', b'U')
        );
    }

    #[test]
    fn s5f_kernel_encode_context_rejects_n_or_gap() {
        assert_eq!(
            S5FKernel::encode_context(b'A', b'A', b'N', b'A', b'A'),
            None
        );
        assert_eq!(
            S5FKernel::encode_context(b'A', b'A', b'A', b'-', b'A'),
            None
        );
        assert_eq!(
            S5FKernel::encode_context(b'?', b'A', b'A', b'A', b'A'),
            None
        );
    }

    #[test]
    fn s5f_kernel_substitution_row_indexing_consistent_with_substitution() {
        // For a kernel with distinct values per (context, dest), the
        // row accessor should agree with the per-base accessor.
        let mut sub = vec![0.0; S5F_SUBSTITUTION_LEN];
        // Set context 5: dest A=0.1, C=0.2, G=0.3, T=0.4.
        sub[5 * 4] = 0.1;
        sub[5 * 4 + 1] = 0.2;
        sub[5 * 4 + 2] = 0.3;
        sub[5 * 4 + 3] = 0.4;
        let kernel = S5FKernel::new(vec![1.0; S5F_NUM_CONTEXTS], sub);

        assert_eq!(kernel.substitution(5, b'A'), 0.1);
        assert_eq!(kernel.substitution(5, b'C'), 0.2);
        assert_eq!(kernel.substitution(5, b'G'), 0.3);
        assert_eq!(kernel.substitution(5, b'T'), 0.4);
        assert_eq!(kernel.substitution_row(5), [0.1, 0.2, 0.3, 0.4]);
    }

    #[test]
    fn s5f_constants_match_dimensions() {
        assert_eq!(S5F_NUM_CONTEXTS, 1024);
        assert_eq!(S5F_SUBSTITUTION_LEN, 4096);
        assert_eq!(S5F_NUM_CONTEXTS * 4, S5F_SUBSTITUTION_LEN);
    }
}
