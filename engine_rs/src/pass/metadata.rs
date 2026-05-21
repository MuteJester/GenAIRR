use crate::ir::Segment;

/// Static requirement a pass has before it can execute safely.
#[derive(Copy, Clone, Debug, Eq, PartialEq)]
pub enum PassRequirement {
    /// Pass needs reference data in `PassContext::refdata`.
    RefData,
    /// Pass needs an allele assignment for the given segment.
    AlleleAssignment(Segment),
}

/// Static effect a pass has on the simulation state.
#[derive(Copy, Clone, Debug, Eq, PartialEq)]
pub enum PassEffect {
    /// Assigns an allele id to a recombinable segment.
    AssignAllele(Segment),
    /// Updates trim metadata on an assigned allele instance.
    TrimAllele(Segment),
    /// Assembles a germline segment into the nucleotide pool/sequence.
    AssembleSegment(Segment),
    /// Appends a region to the sequence.
    AppendRegion(Segment),
    /// Appends one or more nucleotides without creating a region.
    AppendNucleotides,
    /// Edits existing nucleotide bases without changing pool length.
    EditBases,
    /// Changes pool length and shifts downstream handles/ranges.
    StructuralIndel,
}
