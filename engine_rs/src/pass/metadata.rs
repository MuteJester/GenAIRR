use crate::ir::Segment;

/// Static requirement a pass has before it can execute safely.
#[derive(Copy, Clone, Debug, Eq, PartialEq)]
pub enum PassRequirement {
    /// Pass needs reference data in `PassContext::refdata`.
    RefData,
    /// Pass needs an allele assignment for the given segment.
    AlleleAssignment(Segment),
}

/// **Static, compile-time** declaration of the category of state
/// change a pass produces. Used by the schedule analyzer for
/// ordering / dependency reasoning — *not* by runtime derived-state
/// refresh.
///
/// # Static facts, not runtime consequences
///
/// A pass declares its compile effects up-front via
/// [`crate::pass::Pass::effects`]. Two surfaces describe what a
/// pass does, and they serve different layers of the engine:
///
/// - **`PassCompileEffect`** — what a pass *claims it will do*.
///   Read at compile time by the schedule analyzer to determine
///   pass ordering, requirement satisfaction, and downstream
///   feasibility domains. Never consulted at runtime by the
///   derived-state refresh hooks.
/// - **[`crate::ir::SimulationEvent`]** — what a pass
///   *actually emitted* during execution via its internal
///   [`crate::ir::SimulationBuilder`]s. Read at runtime by
///   [`crate::live_call::LiveCallRefreshHook`] (and any future
///   derived-state consumer) to decide what to refresh.
///
/// A pass can declare effects that don't materialize as events
/// (e.g. an `EditBases` pass that no-ops because every site was
/// already at the chosen base) or emit events without declaring
/// the matching effect (a quick-and-dirty test pass that skips
/// the declaration). The schedule trusts declarations; the
/// runtime trusts events.
#[derive(Copy, Clone, Debug, Eq, PartialEq)]
pub enum PassCompileEffect {
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

/// Transitional alias for the legacy name. Reads as "the static
/// pass effect" but the type behind it has been renamed to
/// [`PassCompileEffect`] to make the static-vs-runtime distinction
/// unambiguous. New code should use `PassCompileEffect` directly;
/// this alias keeps the old `PassEffect` spelling working while
/// callers migrate.
pub type PassEffect = PassCompileEffect;
