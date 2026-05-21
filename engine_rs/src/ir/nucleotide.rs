use super::Segment;

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
