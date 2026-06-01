use super::Segment;

// ──────────────────────────────────────────────────────────────────
// complement_base — Watson-Crick single-base complement
// ──────────────────────────────────────────────────────────────────

/// Watson-Crick complement of a single nucleotide byte.
///
/// Convention mirrors the existing string-level `reverse_complement`
/// in [`airr_record/sequence.rs`](../../airr_record/sequence.rs)
/// so the per-byte and per-string complement paths stay aligned:
///
/// - `A ↔ T`, `C ↔ G` (case-preserving).
/// - `U` (RNA) maps to `A` to match `reverse_complement`.
/// - `N` / `n` are self-complement (degenerate base).
/// - **Any other byte is returned unchanged.** Defensive choice —
///   the strict alphabet validation lives on `ReferenceAlphabet`,
///   so leaking bytes (gap markers, lowercase mutated bases, etc.)
///   round-trip rather than silently coerce to `N`.
///
/// Wired into `AssembleSegmentPass` by Slice B: when a D
/// `AlleleInstance` carries
/// [`crate::assignment::SegmentOrientation::ReverseComplement`], the
/// assembler calls this helper once per retained byte before pushing
/// the nucleotide onto the pool.
pub fn complement_base(b: u8) -> u8 {
    match b {
        b'A' => b'T',
        b'T' => b'A',
        b'U' => b'A',
        b'C' => b'G',
        b'G' => b'C',
        b'a' => b't',
        b't' => b'a',
        b'u' => b'a',
        b'c' => b'g',
        b'g' => b'c',
        b'N' | b'n' => b,
        other => other,
    }
}

/// Per-nucleotide flags. Bitset over u8 — see `flag` constants.
#[derive(Copy, Clone, Eq, PartialEq, Hash, Debug, Default)]
pub struct NucFlags(u8);

impl NucFlags {
    /// Empty flag set.
    pub const fn empty() -> Self {
        Self(0)
    }

    /// Construct from a raw u8 bitfield.
    pub const fn from_bits(bits: u8) -> Self {
        Self(bits)
    }

    /// Raw bits (for serialization / introspection).
    pub const fn bits(self) -> u8 {
        self.0
    }

    /// Test whether all bits in `flags` are set.
    pub const fn contains(self, flags: NucFlags) -> bool {
        (self.0 & flags.0) == flags.0
    }

    /// Set the bits in `flags`. Returns the new flags value.
    #[must_use]
    pub const fn with(self, flags: NucFlags) -> Self {
        Self(self.0 | flags.0)
    }

    /// Clear the bits in `flags`. Returns the new flags value.
    #[must_use]
    pub const fn without(self, flags: NucFlags) -> Self {
        Self(self.0 & !flags.0)
    }
}

/// Flag constants. Each is a single-bit `NucFlags` value.
pub mod flag {
    use super::NucFlags;

    /// Nucleotide is a P-nucleotide (templated palindromic complement
    /// of the segment edge during V(D)J recombination).
    pub const P_NUC: NucFlags = NucFlags::from_bits(1 << 0);
    /// Nucleotide is an N-nucleotide (TdT-added at the junction).
    pub const N_NUC: NucFlags = NucFlags::from_bits(1 << 1);
    /// Nucleotide sits inside the junction region (V Cys -> J W/F + 3).
    pub const JUNCTION: NucFlags = NucFlags::from_bits(1 << 2);
    /// Nucleotide came from a D segment that was inverted (RSS-symmetric
    /// V(D)J inversion event, biologically real for D segments only).
    pub const INVERTED: NucFlags = NucFlags::from_bits(1 << 3);
    /// Nucleotide was inserted by a downstream indel pass — has no
    /// germline-allele provenance.
    pub const INDEL_INSERTED: NucFlags = NucFlags::from_bits(1 << 4);
}

/// Position in the source allele, or `NONE` for synthetic bases
/// (NP, P-nuc, contaminant, indel-inserted) with no germline
/// provenance.
#[derive(Copy, Clone, Eq, PartialEq, Hash, Debug)]
#[repr(transparent)]
pub struct GermlinePos(u16);

impl GermlinePos {
    /// Sentinel value meaning "no germline provenance".
    pub const NONE: Self = Self(u16::MAX);

    /// Construct a real allele position. Panics if `v == u16::MAX`
    /// — that value is reserved for `NONE`.
    pub const fn pos(v: u16) -> Self {
        assert!(
            v != u16::MAX,
            "GermlinePos::pos: u16::MAX is reserved for NONE"
        );
        Self(v)
    }

    /// Project to `Option<u16>`. `Some(v)` for a real position,
    /// `None` for `NONE`.
    pub const fn get(self) -> Option<u16> {
        if self.0 == u16::MAX {
            None
        } else {
            Some(self.0)
        }
    }

    /// Returns `true` when this is the `NONE` sentinel.
    pub const fn is_none(self) -> bool {
        self.0 == u16::MAX
    }

    /// Returns `true` when this carries a real allele position.
    pub const fn is_some(self) -> bool {
        !self.is_none()
    }
}

/// A single nucleotide in the simulation IR.
#[derive(Copy, Clone, Eq, PartialEq, Hash, Debug)]
pub struct Nucleotide {
    /// Current base — usually one of `b'A'`, `b'C'`, `b'G'`, `b'T'`
    /// plus `b'N'`/`b'n'` for ambiguous / quality-masked bases.
    pub base: u8,

    /// Wild-type base at this position — the germline reference.
    pub germline: u8,

    /// Position in the source allele, or `GermlinePos::NONE` for
    /// synthetic bases with no germline provenance.
    pub germline_pos: GermlinePos,

    /// Biological segment role.
    pub segment: Segment,

    /// Per-nucleotide flags.
    pub flags: NucFlags,
}

impl Nucleotide {
    /// Construct a germline-derived nucleotide (base == germline).
    pub const fn germline(base: u8, germline_pos: u16, segment: Segment) -> Self {
        Self {
            base,
            germline: base,
            germline_pos: GermlinePos::pos(germline_pos),
            segment,
            flags: NucFlags::empty(),
        }
    }

    /// Construct a synthetic nucleotide with no germline provenance.
    pub const fn synthetic(base: u8, segment: Segment, flags: NucFlags) -> Self {
        Self {
            base,
            germline: base,
            germline_pos: GermlinePos::NONE,
            segment,
            flags,
        }
    }

    /// Return a copy of this nucleotide with `flags` set (logical-OR
    /// with the existing flag bits). Mirrors [`NucFlags::with`] and
    /// preserves the persistent-IR convention: receiver unchanged.
    ///
    /// Slice B uses this to stamp [`flag::INVERTED`] onto every base
    /// `AssembleSegmentPass(D)` emits when the D assignment carries
    /// [`crate::assignment::SegmentOrientation::ReverseComplement`].
    #[must_use]
    pub const fn with_flags(self, flags: NucFlags) -> Self {
        Self {
            flags: self.flags.with(flags),
            ..self
        }
    }
}

/// Map a DNA base byte to its 0..=3 numeric encoding (A=0, C=1, G=2,
/// T=3). Returns `None` for any other byte.
pub const fn encode_base(b: u8) -> Option<u8> {
    match b.to_ascii_uppercase() {
        b'A' => Some(0),
        b'C' => Some(1),
        b'G' => Some(2),
        b'T' | b'U' => Some(3),
        _ => None,
    }
}

// ──────────────────────────────────────────────────────────────────
// Tests
// ──────────────────────────────────────────────────────────────────

#[cfg(test)]
mod tests {
    use super::*;

    // ── complement_base (Slice A prep for D inversion) ───────────

    #[test]
    fn complement_base_canonical_dna_pairs() {
        // A ↔ T, C ↔ G — the load-bearing rule. Self-inverse:
        // complement(complement(b)) == b for every canonical base.
        for b in [b'A', b'C', b'G', b'T'] {
            assert_eq!(complement_base(complement_base(b)), b);
        }
        assert_eq!(complement_base(b'A'), b'T');
        assert_eq!(complement_base(b'T'), b'A');
        assert_eq!(complement_base(b'C'), b'G');
        assert_eq!(complement_base(b'G'), b'C');
    }

    #[test]
    fn complement_base_preserves_lowercase_case() {
        // Mirrors the existing string-level `reverse_complement` in
        // airr_record/sequence.rs: lowercase input → lowercase
        // complement. SHM substitution traces lowercase, so an
        // inverted-D walker downstream must keep the case.
        assert_eq!(complement_base(b'a'), b't');
        assert_eq!(complement_base(b't'), b'a');
        assert_eq!(complement_base(b'c'), b'g');
        assert_eq!(complement_base(b'g'), b'c');
        // Round-trip preserves case.
        for b in [b'a', b'c', b'g', b't'] {
            assert_eq!(complement_base(complement_base(b)), b);
        }
    }

    #[test]
    fn complement_base_rna_uracil_maps_to_adenine() {
        // RNA convention: U complements to A (and vice versa,
        // when the receiver was T → A). Matches the existing
        // string-level helper. NOT self-inverse — by design,
        // since A's complement is T not U.
        assert_eq!(complement_base(b'U'), b'A');
        assert_eq!(complement_base(b'u'), b'a');
    }

    #[test]
    fn complement_base_n_self_complement_preserves_case() {
        // Degenerate base — self-complement per IUPAC convention.
        assert_eq!(complement_base(b'N'), b'N');
        assert_eq!(complement_base(b'n'), b'n');
    }

    #[test]
    fn complement_base_unknown_bytes_pass_through_unchanged() {
        // Defensive choice: unknown bytes (gaps, IUPAC ambiguity
        // codes other than N, stray ASCII) round-trip rather than
        // silently coerce to N. The strict-alphabet validator is
        // upstream — this helper does not police the input.
        assert_eq!(complement_base(b'.'), b'.');
        assert_eq!(complement_base(b'-'), b'-');
        assert_eq!(complement_base(b'R'), b'R'); // IUPAC purine
        assert_eq!(complement_base(b'X'), b'X');
        assert_eq!(complement_base(0), 0);
        assert_eq!(complement_base(255), 255);
    }

    // ── NucFlags::INVERTED scaffold (Slice A prep) ───────────────

    #[test]
    fn nuc_flags_inverted_round_trips_through_with_without() {
        // Pin that the INVERTED bit can be set, tested, cleared,
        // and re-tested without interfering with neighboring flags.
        // Slice A doesn't fire INVERTED yet; this test guards the
        // scaffold so a future bitset reshuffle doesn't break the
        // promise the design doc relies on.
        let empty = NucFlags::empty();
        assert!(!empty.contains(flag::INVERTED));

        let set = empty.with(flag::INVERTED);
        assert!(set.contains(flag::INVERTED));
        // Neighboring flags untouched by INVERTED.
        assert!(!set.contains(flag::JUNCTION));
        assert!(!set.contains(flag::N_NUC));

        let cleared = set.without(flag::INVERTED);
        assert!(!cleared.contains(flag::INVERTED));
        assert_eq!(cleared, empty);
    }

    #[test]
    fn nuc_flags_inverted_composes_with_other_flags() {
        // An inverted, junction-flagged base is a legitimate
        // combination — the design doc proposes INVERTED for
        // every byte AssembleSegmentPass emits under
        // SegmentOrientation::ReverseComplement, and JUNCTION
        // is the existing flag for junction-window bases.
        let f = NucFlags::empty()
            .with(flag::INVERTED)
            .with(flag::JUNCTION);
        assert!(f.contains(flag::INVERTED));
        assert!(f.contains(flag::JUNCTION));
        assert!(!f.contains(flag::N_NUC));
        // bits() round-trip via from_bits.
        assert_eq!(NucFlags::from_bits(f.bits()), f);
    }
}
