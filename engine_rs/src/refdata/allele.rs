use crate::ir::Segment;

/// One germline allele in a reference data set.
///
/// Immutable once constructed. Per-simulation state (which allele was
/// sampled, what trim was applied, what ambiguity set the post-trim
/// retained bases project to) lives on `AlleleInstance`.
///
/// **Field discipline:**
/// - `seq` is uppercase (`b'A'`, `b'C'`, `b'G'`, `b'T'`). Mixed case
///   is not allowed — case carries semantic meaning later in the
///   simulation (lowercase = SHM-mutated, etc.) and should not leak
///   into reference data.
/// - `anchor` is `Some(pos)` when the allele has a conserved codon
///   (Cys for V, W/F for J) at position `pos` (0-indexed within
///   `seq`). `None` for anchorless / partial alleles. The
///   `AnchorPreserved` contract reads this field.
/// - `name` is the canonical allele identifier (e.g.,
///   `"IGHV1-2*01"`). `gene` is the truncation to the gene level
///   (e.g., `"IGHV1-2"`).
#[derive(Clone, Debug)]
pub struct Allele {
    pub name: String,
    pub gene: String,
    pub seq: Vec<u8>,
    pub segment: Segment,
    pub anchor: Option<u16>,
}

impl Allele {
    /// Length of the allele sequence in nucleotides.
    pub fn len(&self) -> u32 {
        self.seq.len() as u32
    }

    /// Whether the allele has a known anchor position.
    pub fn has_anchor(&self) -> bool {
        self.anchor.is_some()
    }
}
