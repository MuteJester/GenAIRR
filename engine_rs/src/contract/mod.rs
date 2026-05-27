//! Contracts — composable predicates over the simulation IR.
//!
//! ## What a contract is
//!
//! Per design doc §6 + D6 + D7, a contract is a first-class
//! predicate that asserts something about the simulation. Contracts
//! have two architectural modes:
//!
//! 1. **Verify** — given a simulation IR (and optionally the
//!    reference data), decide whether the contract holds. Returns
//!    `Ok(())` if satisfied, `Err(ContractViolation)` with a
//!    structured reason if not. Used at any point where invariants
//!    need to be checked: build-time validation, debug
//!    introspection, post-pipeline assertions.
//!
//! 2. **Filter** — given a candidate sampling action and the
//!    current state, decide whether the action is admissible
//!    *before* sampling. This is what makes constraint-aware
//!    sampling work (D6 — `respect=[productive()]`).
//!
//! ## Module layout
//!
//! - [`anchor_preserved`] — V/D/J anchor codon must remain in the
//!   retained slice after trimming.
//! - [`productive_junction_frame`] — junction length divisible by 3.
//! - [`no_stop_codon_in_junction`] — no stop codon inside the junction
//!   (with NP-base and substitution filtering).
//! - [`set`] — `ContractSet` composition.
//!
//! The trait surface itself, `ChoiceContext`, `ContractViolation`,
//! and the canonical `productive()` bundle live in this `mod.rs`.

use crate::address::ChoiceAddress;
use crate::contract::junction_stop_state::JunctionStopState as JunctionStopStateInner;
use crate::ir::{NucHandle, Segment, Simulation};
use crate::refdata::RefDataConfig;
use crate::trace::ChoiceValue;

pub mod admissible_set;
pub(crate) mod admit_mask_observer;
pub mod anchor_preserved;
pub mod junction_stop_state;
pub mod no_stop_codon_in_junction;
pub mod productive_junction_frame;
pub mod set;

#[cfg(test)]
mod mask_vs_admits_equivalence;

pub use admissible_set::{
    BaseMask, IndelEventClass, IndelKindHint, LengthSupport, TrimEnd, TrimTarget,
};

#[cfg(test)]
pub(crate) mod test_support;

pub use anchor_preserved::AnchorPreserved;
pub use junction_stop_state::JunctionStopState;
pub use no_stop_codon_in_junction::NoStopCodonInJunction;
pub use productive_junction_frame::ProductiveJunctionFrame;
pub use set::ContractSet;

// ──────────────────────────────────────────────────────────────────
// ContractViolation — structured diagnostic for a failed verify
// ──────────────────────────────────────────────────────────────────

/// One reason a contract verification failed.
///
/// Carries a stable contract identifier (`contract_name`) for
/// programmatic dispatch and a human-readable `reason` string for
/// diagnostics. The strict-mode failure shape (D7) builds on this.
#[derive(Clone, Debug, Eq, PartialEq)]
pub struct ContractViolation {
    pub contract_name: String,
    pub reason: String,
}

impl ContractViolation {
    pub fn new(contract_name: impl Into<String>, reason: impl Into<String>) -> Self {
        Self {
            contract_name: contract_name.into(),
            reason: reason.into(),
        }
    }
}

// ──────────────────────────────────────────────────────────────────
// ChoiceContext — optional execution context for candidate filtering
// ──────────────────────────────────────────────────────────────────

/// Semantic class of the candidate choice being filtered.
///
/// Addresses remain useful for trace readability and backwards
/// compatibility, but contracts should prefer this typed signal when
/// deciding whether a candidate has biological meaning.
#[derive(Copy, Clone, Debug, Eq, PartialEq)]
#[non_exhaustive]
pub enum ChoiceKind {
    /// No additional semantic class is known.
    Plain,
    /// Candidate is the destination base for a substitution at
    /// `ChoiceContext::target`.
    TargetedBaseSubstitution,
    /// Candidate is a structural insertion at `ChoiceContext::target`.
    IndelInsertion,
    /// Candidate is a structural deletion at `ChoiceContext::target`.
    IndelDeletion,
}

/// Extra context for a candidate choice being filtered by contracts.
///
/// Plain addressed choices only carry `"np.np1.bases[3]"` plus the
/// candidate value. Some contracts need bounded local execution
/// context to remain precise:
/// - NP-base filters need the current draw index and total planned
///   draw count to distinguish known future fixed bases from still
///   random future bases.
/// - Site-based transforms need the target nucleotide handle so a
///   contract can evaluate the exact post-candidate local state.
#[derive(Copy, Clone, Debug)]
#[non_exhaustive]
pub struct ChoiceContext<'a> {
    /// Typed address for the choice being filtered when it belongs
    /// to the built-in address vocabulary. Persisted traces still
    /// store the string projection; this field lets contracts stop
    /// parsing that string in hot-path predicate code.
    pub address: Option<ChoiceAddress>,
    pub draw_index: Option<u32>,
    pub draw_count: Option<u32>,
    pub target: Option<NucHandle>,
    pub kind: ChoiceKind,
    /// Precomputed junction-stop state for the active record. Set
    /// by `GenerateNPPass` before its NP-base draw loop so the
    /// `NoStopCodonInJunction` filter can use the O(1) fast path
    /// instead of rebuilding the hypothetical junction buffer per
    /// candidate. `None` outside the NP-base hot path, in which
    /// case the contract falls back to the slow rebuild path.
    pub junction_stop_state: Option<&'a JunctionStopStateInner>,
}

impl<'a> ChoiceContext<'a> {
    pub const fn none() -> Self {
        Self {
            address: None,
            draw_index: None,
            draw_count: None,
            target: None,
            kind: ChoiceKind::Plain,
            junction_stop_state: None,
        }
    }

    pub const fn indexed(draw_index: u32, draw_count: u32) -> Self {
        Self {
            address: None,
            draw_index: Some(draw_index),
            draw_count: Some(draw_count),
            target: None,
            kind: ChoiceKind::Plain,
            junction_stop_state: None,
        }
    }

    pub const fn targeted_base_substitution(
        draw_index: u32,
        draw_count: u32,
        target: NucHandle,
    ) -> Self {
        Self {
            address: None,
            draw_index: Some(draw_index),
            draw_count: Some(draw_count),
            target: Some(target),
            kind: ChoiceKind::TargetedBaseSubstitution,
            junction_stop_state: None,
        }
    }

    pub const fn indel_insertion(draw_index: u32, draw_count: u32, target: NucHandle) -> Self {
        Self {
            address: None,
            draw_index: Some(draw_index),
            draw_count: Some(draw_count),
            target: Some(target),
            kind: ChoiceKind::IndelInsertion,
            junction_stop_state: None,
        }
    }

    pub const fn indel_deletion(draw_index: u32, draw_count: u32, target: NucHandle) -> Self {
        Self {
            address: None,
            draw_index: Some(draw_index),
            draw_count: Some(draw_count),
            target: Some(target),
            kind: ChoiceKind::IndelDeletion,
            junction_stop_state: None,
        }
    }

    pub const fn indel_deletion_noop(draw_index: u32, draw_count: u32) -> Self {
        Self {
            address: None,
            draw_index: Some(draw_index),
            draw_count: Some(draw_count),
            target: None,
            kind: ChoiceKind::IndelDeletion,
            junction_stop_state: None,
        }
    }

    /// Attach a precomputed `JunctionStopState` reference to this
    /// `ChoiceContext`. The `NoStopCodonInJunction` filter consults
    /// it when present to take the O(1) fast path.
    #[must_use]
    pub fn with_junction_stop_state(mut self, state: &'a JunctionStopStateInner) -> Self {
        self.junction_stop_state = Some(state);
        self
    }

    /// Attach the typed built-in address for this candidate.
    #[must_use]
    pub const fn with_address(mut self, address: ChoiceAddress) -> Self {
        self.address = Some(address);
        self
    }

    /// Attach `address` only when this context does not already
    /// carry one. This is useful for compatibility adapters that
    /// parse a legacy string address at the boundary.
    #[must_use]
    pub const fn with_address_if_missing(mut self, address: Option<ChoiceAddress>) -> Self {
        if self.address.is_none() {
            self.address = address;
        }
        self
    }

    /// Stable string projection for diagnostics and legacy
    /// compatibility.
    pub fn address_string(self) -> Option<String> {
        self.address.map(String::from)
    }
}

// PartialEq / Eq manually so the borrowed reference doesn't have
// to participate (we compare structural fields only). Sufficient
// for existing test sites that match on draw_index / draw_count.
impl PartialEq for ChoiceContext<'_> {
    fn eq(&self, other: &Self) -> bool {
        self.draw_index == other.draw_index
            && self.address == other.address
            && self.draw_count == other.draw_count
            && self.target == other.target
            && self.kind == other.kind
            && std::ptr::eq(
                self.junction_stop_state
                    .map(|s| s as *const _)
                    .unwrap_or(std::ptr::null()),
                other
                    .junction_stop_state
                    .map(|s| s as *const _)
                    .unwrap_or(std::ptr::null()),
            )
    }
}
impl Eq for ChoiceContext<'_> {}

impl Default for ChoiceContext<'_> {
    fn default() -> Self {
        Self::none()
    }
}

// ──────────────────────────────────────────────────────────────────
// Contract trait
// ──────────────────────────────────────────────────────────────────

/// Stable semantic identity for built-in contracts.
///
/// `Contract::name()` remains the human-readable diagnostic surface.
/// The compiler uses `ContractKind` for typed dispatch so build-time
/// validation does not depend on parsing names.
#[derive(Copy, Clone, Debug, Eq, PartialEq, Hash)]
pub enum ContractKind {
    ProductiveJunctionFrame,
    NoStopCodonInJunction,
    AnchorPreserved { segment: Segment },
    Custom,
}

/// A predicate over the simulation IR.
///
/// Currently only `verify` is surfaced. Filter methods (and
/// possibly an `upstream_bound` for backward constraint
/// propagation) will land as defaulted trait methods so existing
/// contract implementations continue to compile.
pub trait Contract {
    /// Stable, human-readable identifier for this contract.
    /// Convention: hierarchical-string addresses matching D3 (e.g.,
    /// `"anchor_preserved.v"`, `"productive_junction_frame"`).
    fn name(&self) -> &str;

    /// Typed semantic identity for built-in contracts.
    ///
    /// Custom contracts inherit `Custom`. Built-ins override this so
    /// the compiled simulator can validate preconditions without
    /// downcasting trait objects or string-matching names.
    fn kind(&self) -> ContractKind {
        ContractKind::Custom
    }

    /// Verify mode: does this contract hold for `sim`?
    ///
    /// `refdata` is optional because some contracts don't need it
    /// (e.g., a contract that only inspects the assembled bases).
    /// Contracts that *do* need it should treat `None` as
    /// "insufficient data — skip" by returning `Ok(())`, not by
    /// panicking. This keeps verification safe to call in any
    /// context.
    fn verify(
        &self,
        sim: &Simulation,
        refdata: Option<&RefDataConfig>,
    ) -> Result<(), ContractViolation>;

    /// Typed candidate filter — the contract's per-candidate
    /// verdict under the v3.x constrain-before-propose API.
    ///
    /// Used by constraint-aware sampling: before a sampling pass
    /// commits to a draw, it asks every active contract whether the
    /// candidate is admissible at the current `ChoiceContext`. Only
    /// candidates all contracts admit are accepted.
    ///
    /// **Default behaviour**: `Ok(())` — "no opinion." Contracts
    /// that can usefully prune candidates at sampling time override
    /// this method. Contracts that can only check after a transform
    /// applies (e.g., `NoStopCodonInJunction` looking at codons
    /// that don't exist yet) keep the default; their constraints
    /// get enforced by the typed support hooks below where
    /// available, or by post-event / strict verification for
    /// whole-IR checks.
    ///
    /// **Returning `Err` is not a fatal failure** — it's "this
    /// specific candidate is inadmissible." The caller (the
    /// sampling pass) skips the candidate and tries another.
    fn admits_typed(
        &self,
        sim: &Simulation,
        refdata: Option<&RefDataConfig>,
        context: ChoiceContext<'_>,
        candidate: &ChoiceValue,
    ) -> Result<(), ContractViolation> {
        let _ = (sim, refdata, context, candidate);
        Ok(())
    }

    /// Post-event filter mode for structural candidates.
    ///
    /// Structural edits such as indels are best evaluated against the
    /// complete hypothetical post-event IR, because they can shift
    /// ranges, frame phases, anchors, and codon rails. The default
    /// implementation delegates to `verify(post_sim, refdata)`, giving
    /// every existing contract safe structural filtering without each
    /// contract needing a bespoke indel implementation.
    ///
    /// The typed `context` carries the trace address (via
    /// `context.address`), draw index, and kind. Contracts that need
    /// to dispatch on the address read it from there rather than a
    /// separate string parameter.
    fn admits_post_event(
        &self,
        pre_sim: &Simulation,
        post_sim: &Simulation,
        refdata: Option<&RefDataConfig>,
        context: ChoiceContext<'_>,
    ) -> Result<(), ContractViolation> {
        let _ = (pre_sim, context);
        self.verify(post_sim, refdata)
    }

    // ──────────────────────────────────────────────────────────
    // v3.0 constrain-before-propose API
    //
    // Split per-kind methods so contract authors override only the
    // decision shapes they care about. Composition lives in
    // `ContractSet`; the default implementations below mean
    // "no opinion" for an unmodified contract.
    //
    // **Invariant the contract author maintains**: a candidate `v`
    // is admitted by `admits` (the existing predicate API) iff `v`
    // is in the support returned by the matching per-kind method
    // below. Distribution-invariant tests in `tests/` should
    // verify this on every built-in contract.
    // ──────────────────────────────────────────────────────────

    /// Return the bitmask of admissible canonical bases (A/C/G/T)
    /// for a per-site substitution at the given pool handle.
    ///
    /// Default: [`BaseMask::UNCONSTRAINED`] — all four bases
    /// admissible. Contracts that can usefully prune the per-site
    /// support (e.g. `NoStopCodonInJunction` for junction
    /// positions, `AnchorPreserved` for anchor codon positions)
    /// override this.
    ///
    /// **Non-canonical writes** (lowercase, `N`, IUPAC ambiguity)
    /// do not flow through this mask — they use
    /// [`Self::admits_fixed_base_at`] instead. The 4-bit mask is
    /// the canonical hot path; the candidate check is the
    /// pinned-value path.
    fn admissible_bases_at(
        &self,
        sim: &Simulation,
        refdata: Option<&RefDataConfig>,
        site: NucHandle,
    ) -> BaseMask {
        let _ = (sim, refdata, site);
        BaseMask::UNCONSTRAINED
    }

    /// Yes/no candidate check for a pinned non-canonical write
    /// (e.g. `N` injection, lowercase quality marker that maps to
    /// the same canonical base). Defaults to `true` (admitted)
    /// for contracts with no opinion on the pinned write.
    ///
    /// Used by passes like `NCorruptionPass` whose candidate value
    /// is structurally fixed: there's no 4-base support to
    /// narrow; the contract just decides whether the specific
    /// byte at the specific site is admissible.
    fn admits_fixed_base_at(
        &self,
        sim: &Simulation,
        refdata: Option<&RefDataConfig>,
        site: NucHandle,
        byte: u8,
    ) -> bool {
        let _ = (sim, refdata, site, byte);
        true
    }

    /// Classify an indel candidate's impact on this contract.
    ///
    /// Returns one of:
    /// - [`IndelEventClass::FrameNeutral`] — no opinion / no
    ///   effect (the default).
    /// - [`IndelEventClass::FrameDelta`] — the event introduces a
    ///   ±1 frame shift the cross-slot coordinator must account
    ///   for (e.g. junction-site insertions / deletions under
    ///   `ProductiveJunctionFrame`).
    /// - [`IndelEventClass::Forbidden`] — the event is rejected
    ///   outright (e.g. deleting through the V anchor under
    ///   `AnchorPreserved`).
    ///
    /// The indel pass collects per-slot classifications, then runs
    /// a mod-3 DP over the FrameDelta values to enumerate
    /// frame-preserving tuples. The mod-3 DP handles the **frame
    /// part only**; passes must still run a final
    /// `admits_post_event` check on the sampled tuple to catch
    /// contracts whose admissibility depends on exact bases and
    /// sites (e.g. `NoStopCodonInJunction` after frame-balanced
    /// insertions can still produce a stop codon).
    fn admissible_indel_class_at(
        &self,
        sim: &Simulation,
        refdata: Option<&RefDataConfig>,
        site: u32,
        kind: IndelKindHint,
    ) -> IndelEventClass {
        let _ = (sim, refdata, site, kind);
        IndelEventClass::FrameNeutral
    }

    /// Return the support of admissible trim lengths for an
    /// end-loss-style pass targeting [`TrimTarget`].
    ///
    /// `requested_max` is the upper bound the pass would otherwise
    /// apply (its sampled length, clamped to pool length). The
    /// returned [`LengthSupport`] is the subset of `0..=requested_max`
    /// the contract admits. Default: `LengthSupport::Full(requested_max)`
    /// — full range admissible.
    ///
    /// Contracts that own segment-anchor geometry (`AnchorPreserved`)
    /// narrow this to the lengths that don't trim through the
    /// anchor codon. Contracts that own frame
    /// (`ProductiveJunctionFrame`) can narrow further to lengths
    /// that preserve the junction frame.
    fn admissible_trim_lengths(
        &self,
        sim: &Simulation,
        refdata: Option<&RefDataConfig>,
        target: TrimTarget,
        requested_max: u32,
    ) -> LengthSupport {
        let _ = (sim, refdata, target);
        LengthSupport::Full(requested_max)
    }
}

// ──────────────────────────────────────────────────────────────────
// productive() — the canonical productive-sequence bundle
// ──────────────────────────────────────────────────────────────────

/// The canonical productive-sequence contract bundle:
///
/// 1. `ProductiveJunctionFrame` — junction length divisible by 3
/// 2. `NoStopCodonInJunction` — no stops inside the junction
/// 3. `AnchorPreserved::V` — V Cys codon retained after V trim
/// 4. `AnchorPreserved::J` — J W/F codon retained after J trim
///
/// All four must hold for a sequence to be considered productive
/// in the standard biological sense. The DSL `respect=[productive()]`
/// compiles to this bundle internally.
pub fn productive() -> ContractSet {
    ContractSet::new()
        .with(Box::new(ProductiveJunctionFrame::new()))
        .with(Box::new(NoStopCodonInJunction::new()))
        .with(Box::new(AnchorPreserved::new(Segment::V)))
        .with(Box::new(AnchorPreserved::new(Segment::J)))
}

// ──────────────────────────────────────────────────────────────────
// Cross-cutting tests — productive() bundle, default admits, etc.
// ──────────────────────────────────────────────────────────────────

#[cfg(test)]
mod tests {
    use super::*;
    use crate::contract::test_support::{
        make_assembled_sim_from_refdata, make_vj_with_anchor_codons,
    };

    #[test]
    fn contract_violation_construction_round_trip() {
        let v = ContractViolation::new("test.contract", "something failed");
        assert_eq!(v.contract_name, "test.contract");
        assert_eq!(v.reason, "something failed");
    }

    // ── productive() bundle tests (D.5) ────────────────────────────

    #[test]
    fn productive_bundle_contains_four_contracts() {
        let s = productive();
        assert_eq!(s.len(), 4);
        let names: Vec<&str> = s.iter().map(|c| c.name()).collect();
        assert!(names.contains(&"productive_junction_frame"));
        assert!(names.contains(&"no_stop_codon_in_junction"));
        assert!(names.contains(&"anchor_preserved.v"));
        assert!(names.contains(&"anchor_preserved.j"));
    }

    #[test]
    fn productive_bundle_verifies_clean_sim() {
        // V anchor codon GGG (Gly), J anchor codon TTT (Phe). In-frame.
        let cfg = make_vj_with_anchor_codons(b"GGG", b"TTT");
        let sim = make_assembled_sim_from_refdata(&cfg);

        assert!(productive().verify(&sim, Some(&cfg)).is_ok());
    }

    #[test]
    fn productive_bundle_flags_stop_codon() {
        let cfg = make_vj_with_anchor_codons(b"TAA", b"GGG");
        let sim = make_assembled_sim_from_refdata(&cfg);

        let violations = productive().verify(&sim, Some(&cfg)).unwrap_err();
        let names: Vec<&str> = violations
            .iter()
            .map(|v| v.contract_name.as_str())
            .collect();
        assert!(names.contains(&"no_stop_codon_in_junction"));
    }
}
