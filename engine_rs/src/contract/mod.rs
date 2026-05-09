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
//!    need to be checked: build-time validation (D7 Phase 1),
//!    debug introspection, post-pipeline assertions.
//!
//! 2. **Filter** — given a candidate sampling action and the
//!    current state, decide whether the action is admissible
//!    *before* sampling. This is what makes constraint-aware
//!    sampling work (D6 — `respect=[productive()]`). Filter mode
//!    arrives in Phase D when contracts get wired into the
//!    `Distribution` trait surface.
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

use crate::ir::{NucHandle, Segment, Simulation};
use crate::refdata::RefDataConfig;
use crate::trace::ChoiceValue;

pub mod anchor_preserved;
pub mod no_stop_codon_in_junction;
pub mod productive_junction_frame;
pub mod set;

#[cfg(test)]
pub(crate) mod test_support;

pub use anchor_preserved::AnchorPreserved;
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
/// diagnostics. Phase D's strict-mode failure shape (D7) builds
/// on this.
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
#[derive(Copy, Clone, Debug, Eq, PartialEq)]
#[non_exhaustive]
pub struct ChoiceContext {
    pub draw_index: Option<u32>,
    pub draw_count: Option<u32>,
    pub target: Option<NucHandle>,
    pub kind: ChoiceKind,
}

impl ChoiceContext {
    pub const fn none() -> Self {
        Self {
            draw_index: None,
            draw_count: None,
            target: None,
            kind: ChoiceKind::Plain,
        }
    }

    pub const fn indexed(draw_index: u32, draw_count: u32) -> Self {
        Self {
            draw_index: Some(draw_index),
            draw_count: Some(draw_count),
            target: None,
            kind: ChoiceKind::Plain,
        }
    }

    pub const fn indexed_target(draw_index: u32, draw_count: u32, target: NucHandle) -> Self {
        Self::targeted_base_substitution(draw_index, draw_count, target)
    }

    pub const fn targeted_base_substitution(
        draw_index: u32,
        draw_count: u32,
        target: NucHandle,
    ) -> Self {
        Self {
            draw_index: Some(draw_index),
            draw_count: Some(draw_count),
            target: Some(target),
            kind: ChoiceKind::TargetedBaseSubstitution,
        }
    }

    pub const fn indel_insertion(draw_index: u32, draw_count: u32, target: NucHandle) -> Self {
        Self {
            draw_index: Some(draw_index),
            draw_count: Some(draw_count),
            target: Some(target),
            kind: ChoiceKind::IndelInsertion,
        }
    }

    pub const fn indel_deletion(draw_index: u32, draw_count: u32, target: NucHandle) -> Self {
        Self {
            draw_index: Some(draw_index),
            draw_count: Some(draw_count),
            target: Some(target),
            kind: ChoiceKind::IndelDeletion,
        }
    }

    pub const fn indel_deletion_noop(draw_index: u32, draw_count: u32) -> Self {
        Self {
            draw_index: Some(draw_index),
            draw_count: Some(draw_count),
            target: None,
            kind: ChoiceKind::IndelDeletion,
        }
    }
}

impl Default for ChoiceContext {
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
/// Phase C.9 surfaces only `verify`. Phase D will add filter
/// methods (and possibly an `upstream_bound` for backward
/// constraint propagation) as defaulted trait methods so existing
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

    /// Filter mode: would `candidate` at `address` keep this
    /// contract satisfiable given the current `sim` state?
    ///
    /// Used by constraint-aware sampling (D.6 onward): before a
    /// sampling pass commits to a draw, it asks every active
    /// contract `admits(sim, refdata, addr, candidate)` and only
    /// accepts candidates that all contracts admit.
    ///
    /// **Default behaviour**: `Ok(())` — "always allow." Contracts
    /// that can usefully prune candidates at sampling time
    /// override this method. Contracts that can only check after
    /// a transform applies (e.g., `NoStopCodonInJunction` looking
    /// at codons that don't exist yet) keep the default; their
    /// constraints get enforced via `verify` post-hoc, or by
    /// future contract-aware base-sampling passes (Phase E).
    ///
    /// **Returning `Err` is not a fatal failure** — it's "this
    /// specific candidate is inadmissible." The caller (the
    /// sampling pass via the `ContractSet` from D.5) treats it as
    /// a filter signal and tries another candidate.
    ///
    /// `address` follows the D3 hierarchical-string convention:
    /// `"trim.v_3"`, `"np.np1.length"`, `"sample_allele.v"`, etc.
    /// Contracts dispatch on prefix to handle the addresses they
    /// care about and ignore the rest.
    fn admits(
        &self,
        sim: &Simulation,
        refdata: Option<&RefDataConfig>,
        address: &str,
        candidate: &ChoiceValue,
    ) -> Result<(), ContractViolation> {
        // Defaulted: ignore inputs, always allow.
        let _ = (sim, refdata, address, candidate);
        Ok(())
    }

    /// Context-aware filter mode. Defaults to the simpler `admits`
    /// implementation so existing contracts only override this when
    /// they need execution-local metadata such as indexed draw count.
    fn admits_with_context(
        &self,
        sim: &Simulation,
        refdata: Option<&RefDataConfig>,
        address: &str,
        candidate: &ChoiceValue,
        context: ChoiceContext,
    ) -> Result<(), ContractViolation> {
        let _ = context;
        self.admits(sim, refdata, address, candidate)
    }

    /// Post-event filter mode for structural candidates.
    ///
    /// Structural edits such as indels are best evaluated against the
    /// complete hypothetical post-event IR, because they can shift
    /// ranges, frame phases, anchors, and codon rails. The default
    /// implementation delegates to `verify(post_sim, refdata)`, giving
    /// every existing contract safe structural filtering without each
    /// contract needing a bespoke indel implementation.
    fn admits_post_event(
        &self,
        pre_sim: &Simulation,
        post_sim: &Simulation,
        refdata: Option<&RefDataConfig>,
        address: &str,
        context: ChoiceContext,
    ) -> Result<(), ContractViolation> {
        let _ = (pre_sim, address, context);
        self.verify(post_sim, refdata)
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
/// (Phase F) compiles to this bundle internally.
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
        make_assembled_sim_from_refdata, make_v_anchor_at, make_vj_with_anchor_codons,
    };

    #[test]
    fn contract_violation_construction_round_trip() {
        let v = ContractViolation::new("test.contract", "something failed");
        assert_eq!(v.contract_name, "test.contract");
        assert_eq!(v.reason, "something failed");
    }

    #[test]
    fn default_admits_returns_ok_for_anchor_preserved() {
        let cfg = make_v_anchor_at(10, Some(3));
        let contract = AnchorPreserved::new(Segment::V);
        let sim = Simulation::new();

        // Default admits ignores the address + candidate and returns Ok.
        assert!(contract
            .admits(&sim, Some(&cfg), "trim.v_3", &ChoiceValue::Int(5))
            .is_ok());
        assert!(contract
            .admits(&sim, None, "np.np1.length", &ChoiceValue::Int(3))
            .is_ok());
        assert!(contract
            .admits(
                &sim,
                Some(&cfg),
                "sample_allele.v",
                &ChoiceValue::AlleleId(0)
            )
            .is_ok());
    }

    #[test]
    fn default_admits_returns_ok_for_productive_junction_frame() {
        let contract = ProductiveJunctionFrame::new();
        let sim = Simulation::new();

        // Same default behaviour for the new contract — D.6 will
        // override `admits` on this contract specifically; for now
        // it inherits the trait default.
        assert!(contract
            .admits(&sim, None, "np.np1.length", &ChoiceValue::Int(7))
            .is_ok());
    }

    #[test]
    fn no_stop_codon_admits_passes_vacuously_without_refdata() {
        let contract = NoStopCodonInJunction::new();
        let sim = Simulation::new();

        assert!(contract
            .admits(&sim, None, "np.np1.bases[0]", &ChoiceValue::Base(b'T'))
            .is_ok());
    }

    #[test]
    fn admits_works_through_box_dyn() {
        let contract: Box<dyn Contract> = Box::new(AnchorPreserved::new(Segment::V));
        let sim = Simulation::new();
        assert!(contract
            .admits(&sim, None, "trim.v_3", &ChoiceValue::Int(0))
            .is_ok());
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

    #[test]
    fn productive_bundle_admits_returns_ok_for_default_admits() {
        // None of the productive bundle's contracts override admits
        // yet (D.6 will add it). All defaults → Ok.
        let s = productive();
        let sim = Simulation::new();
        assert!(s
            .admits(&sim, None, "np.np1.length", &ChoiceValue::Int(3))
            .is_ok());
    }
}
