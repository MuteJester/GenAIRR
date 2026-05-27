//! `ContractSet` — composition of multiple contracts (D.5).

use crate::ir::{NucHandle, Simulation};
use crate::refdata::RefDataConfig;
use crate::trace::ChoiceValue;
use std::sync::Arc;

use super::{
    BaseMask, ChoiceContext, Contract, ContractKind, ContractViolation, IndelEventClass,
    IndelKindHint, LengthSupport, TrimTarget,
};

/// A set of contracts that all must hold for the simulation.
///
/// Composition is *intersection*: `verify` succeeds only if every
/// contained contract verifies; `admits` accepts a candidate only
/// if every contained contract admits it. This mirrors the
/// architectural commitment from D6 — "Multiple contract specs
/// intersect — all must hold throughout the simulation."
///
/// **`verify` collects all violations** (not first-only) so
/// diagnostic output shows the full picture of what failed.
///
/// **`admits` short-circuits** on the first inadmissible contract
/// because sampling speed matters and one violation is enough to
/// reject the candidate.
///
/// Contract objects are shared behind `Arc` so compiled simulators can
/// cheaply capture immutable contract bundles. Contract implementations
/// should remain pure predicates over their inputs.
#[derive(Clone)]
pub struct ContractSet {
    contracts: Vec<Arc<dyn Contract>>,
}

impl ContractSet {
    /// Empty set — admits everything, verifies anything.
    pub fn new() -> Self {
        Self {
            contracts: Vec::new(),
        }
    }

    /// Builder-style append. Returns `self` for chaining.
    #[must_use]
    pub fn with(mut self, contract: Box<dyn Contract>) -> Self {
        self.contracts.push(contract.into());
        self
    }

    /// In-place append. Returns `&mut self` for builder-mut chains.
    pub fn add(&mut self, contract: Box<dyn Contract>) -> &mut Self {
        self.contracts.push(contract.into());
        self
    }

    /// Number of contained contracts.
    pub fn len(&self) -> usize {
        self.contracts.len()
    }

    /// Whether the set contains zero contracts.
    pub fn is_empty(&self) -> bool {
        self.contracts.is_empty()
    }

    /// Iterate over the contained contracts.
    pub fn iter(&self) -> impl Iterator<Item = &dyn Contract> {
        self.contracts.iter().map(|c| c.as_ref())
    }

    /// Whether a built-in semantic contract kind is active.
    pub fn contains_kind(&self, kind: ContractKind) -> bool {
        self.contracts
            .iter()
            .any(|contract| contract.kind() == kind)
    }

    /// Verify every contract. Returns `Ok(())` only if all pass.
    /// On failure returns the *complete* list of violations, so
    /// callers see every contract that failed in one pass.
    pub fn verify(
        &self,
        sim: &Simulation,
        refdata: Option<&RefDataConfig>,
    ) -> Result<(), Vec<ContractViolation>> {
        let mut violations = Vec::new();
        for c in &self.contracts {
            if let Err(v) = c.verify(sim, refdata) {
                violations.push(v);
            }
        }
        if violations.is_empty() {
            Ok(())
        } else {
            Err(violations)
        }
    }

    /// Test whether `candidate` at the current `ChoiceContext` is
    /// admissible by every contract. Short-circuits on the first
    /// violator (sampling hot path) and returns its
    /// `ContractViolation`.
    pub fn admits_typed(
        &self,
        sim: &Simulation,
        refdata: Option<&RefDataConfig>,
        context: ChoiceContext<'_>,
        candidate: &ChoiceValue,
    ) -> Result<(), ContractViolation> {
        for c in &self.contracts {
            c.admits_typed(sim, refdata, context, candidate)?;
        }
        Ok(())
    }

    /// Structural-event variant of [`ContractSet::admits_typed`].
    ///
    /// Used for candidates whose admissibility depends on the complete
    /// post-event IR, such as indels that change pool length and shift
    /// downstream sequence metadata.
    pub fn admits_post_event(
        &self,
        pre_sim: &Simulation,
        post_sim: &Simulation,
        refdata: Option<&RefDataConfig>,
        context: ChoiceContext<'_>,
    ) -> Result<(), ContractViolation> {
        for c in &self.contracts {
            c.admits_post_event(pre_sim, post_sim, refdata, context)?;
        }
        Ok(())
    }

    // ──────────────────────────────────────────────────────────
    // v3.0 constrain-before-propose composition
    //
    // Each method below intersects the per-contract supports into
    // a single typed support the pass samples from. Composition is
    // intersection (D6 semantics — every contract in the bundle
    // must admit). The shape-specific reducers live on the support
    // types themselves; this layer just folds.
    // ──────────────────────────────────────────────────────────

    /// Compose every contract's per-site admissible-base mask into
    /// the bundle's intersection mask. The pass samples from
    /// `natural_per_position_kernel & mask` over the surviving
    /// bases.
    pub fn admissible_bases_at(
        &self,
        sim: &Simulation,
        refdata: Option<&RefDataConfig>,
        site: NucHandle,
    ) -> BaseMask {
        let mut acc = BaseMask::UNCONSTRAINED;
        for c in &self.contracts {
            let m = c.admissible_bases_at(sim, refdata, site);
            acc = BaseMask(acc.0 & m.0);
            if !acc.is_satisfiable() {
                // Short-circuit: no canonical base survives.
                return BaseMask::EMPTY;
            }
        }
        acc
    }

    /// Yes/no candidate check for a pinned non-canonical write
    /// (e.g. `N` injection) against every contract in the bundle.
    pub fn admits_fixed_base_at(
        &self,
        sim: &Simulation,
        refdata: Option<&RefDataConfig>,
        site: NucHandle,
        byte: u8,
    ) -> bool {
        self.contracts
            .iter()
            .all(|c| c.admits_fixed_base_at(sim, refdata, site, byte))
    }

    /// Compose every contract's indel-event classification for a
    /// candidate `(site, kind)`. Reducer semantics live on
    /// [`IndelEventClass::compose`] — Forbidden absorbs,
    /// FrameNeutral is identity, conflicting deltas collapse to
    /// Forbidden defensively.
    pub fn admissible_indel_class_at(
        &self,
        sim: &Simulation,
        refdata: Option<&RefDataConfig>,
        site: u32,
        kind: IndelKindHint,
    ) -> IndelEventClass {
        let mut acc = IndelEventClass::FrameNeutral;
        for c in &self.contracts {
            let cls = c.admissible_indel_class_at(sim, refdata, site, kind);
            acc = acc.compose(cls);
            if matches!(acc, IndelEventClass::Forbidden) {
                return IndelEventClass::Forbidden;
            }
        }
        acc
    }

    /// Compose every contract's admissible trim-length support
    /// into the bundle's intersection. Pass samples from
    /// `natural_length_distribution ∩ support` over the surviving
    /// lengths.
    pub fn admissible_trim_lengths(
        &self,
        sim: &Simulation,
        refdata: Option<&RefDataConfig>,
        target: TrimTarget,
        requested_max: u32,
    ) -> LengthSupport {
        let mut acc = LengthSupport::Full(requested_max);
        for c in &self.contracts {
            let s = c.admissible_trim_lengths(sim, refdata, target, requested_max);
            acc = acc.intersect(s);
            if matches!(acc, LengthSupport::Empty) {
                return LengthSupport::Empty;
            }
        }
        acc
    }
}

impl Default for ContractSet {
    fn default() -> Self {
        Self::new()
    }
}

#[cfg(test)]
mod tests {
    use super::super::test_support::{
        make_assembled_sim_from_refdata, make_v_anchor_at, make_vj_with_anchor_codons,
    };
    use super::super::{AnchorPreserved, NoStopCodonInJunction, ProductiveJunctionFrame};
    use super::*;
    use crate::address::ChoiceAddress;
    use crate::assignment::AlleleInstance;
    use crate::ir::Segment;
    use crate::refdata::AlleleId;

    #[test]
    fn contract_set_empty_verifies_and_admits_everything() {
        let s = ContractSet::new();
        assert!(s.is_empty());
        assert_eq!(s.len(), 0);

        let sim = Simulation::new();
        assert!(s.verify(&sim, None).is_ok());
        assert!(s
            .admits_typed(&sim, None, ChoiceContext::none(), &ChoiceValue::Int(42),)
            .is_ok());
    }

    #[test]
    fn contract_set_with_chains_contracts() {
        let s = ContractSet::new()
            .with(Box::new(AnchorPreserved::new(Segment::V)))
            .with(Box::new(ProductiveJunctionFrame::new()));
        assert_eq!(s.len(), 2);
        assert!(!s.is_empty());

        let names: Vec<&str> = s.iter().map(|c| c.name()).collect();
        assert_eq!(
            names,
            vec!["anchor_preserved.v", "productive_junction_frame"]
        );
    }

    #[test]
    fn contract_set_add_in_place() {
        let mut s = ContractSet::new();
        s.add(Box::new(AnchorPreserved::new(Segment::V)));
        s.add(Box::new(AnchorPreserved::new(Segment::J)));
        assert_eq!(s.len(), 2);
    }

    #[test]
    fn contract_set_verify_succeeds_when_all_pass() {
        let cfg = make_vj_with_anchor_codons(b"TGT", b"TGG");
        let sim = make_assembled_sim_from_refdata(&cfg);

        let s = ContractSet::new()
            .with(Box::new(AnchorPreserved::new(Segment::V)))
            .with(Box::new(ProductiveJunctionFrame::new()));
        assert!(s.verify(&sim, Some(&cfg)).is_ok());
    }

    #[test]
    fn contract_set_verify_returns_only_failing_contracts() {
        // V anchor codon TAA → NoStopCodonInJunction violates.
        // V trim_3 = 0, V trim_5 = 0 → AnchorPreserved.V passes.
        // → exactly one violation.
        let cfg = make_vj_with_anchor_codons(b"TAA", b"GGG");
        let sim = make_assembled_sim_from_refdata(&cfg);

        let s = ContractSet::new()
            .with(Box::new(AnchorPreserved::new(Segment::V)))
            .with(Box::new(NoStopCodonInJunction::new()))
            .with(Box::new(ProductiveJunctionFrame::new()));

        let violations = s.verify(&sim, Some(&cfg)).unwrap_err();
        assert_eq!(violations.len(), 1);
        assert_eq!(violations[0].contract_name, "no_stop_codon_in_junction");
    }

    #[test]
    fn contract_set_verify_collects_multiple_real_violations() {
        // Force two genuine violations: AnchorPreserved.V (trim_3=99)
        // AND a contract that always fails. Use a custom test contract.
        struct AlwaysFail;
        impl Contract for AlwaysFail {
            fn name(&self) -> &str {
                "always_fail"
            }
            fn verify(
                &self,
                _sim: &Simulation,
                _refdata: Option<&RefDataConfig>,
            ) -> Result<(), ContractViolation> {
                Err(ContractViolation::new("always_fail", "always fails"))
            }
        }

        let cfg = make_v_anchor_at(10, Some(3));
        let bad_sim = Simulation::new().with_allele_assigned(
            Segment::V,
            AlleleInstance::new(AlleleId::new(0)).with_trim_5(5),
        );

        let s = ContractSet::new()
            .with(Box::new(AnchorPreserved::new(Segment::V)))
            .with(Box::new(AlwaysFail));

        let violations = s.verify(&bad_sim, Some(&cfg)).unwrap_err();
        assert_eq!(violations.len(), 2);
        let names: Vec<&str> = violations
            .iter()
            .map(|v| v.contract_name.as_str())
            .collect();
        assert!(names.contains(&"anchor_preserved.v"));
        assert!(names.contains(&"always_fail"));
    }

    #[test]
    fn contract_set_admits_short_circuits_on_first_violator() {
        // Two contracts, both have admits_typed returning Err. The
        // set should return the first one's violation.
        struct RejectAt(&'static str);
        impl Contract for RejectAt {
            fn name(&self) -> &str {
                self.0
            }
            fn verify(
                &self,
                _sim: &Simulation,
                _refdata: Option<&RefDataConfig>,
            ) -> Result<(), ContractViolation> {
                Ok(())
            }
            fn admits_typed(
                &self,
                _sim: &Simulation,
                _refdata: Option<&RefDataConfig>,
                _context: ChoiceContext<'_>,
                _candidate: &ChoiceValue,
            ) -> Result<(), ContractViolation> {
                Err(ContractViolation::new(self.name(), "rejected"))
            }
        }

        let s = ContractSet::new()
            .with(Box::new(RejectAt("first")))
            .with(Box::new(RejectAt("second")));

        let sim = Simulation::new();
        let err = s
            .admits_typed(&sim, None, ChoiceContext::none(), &ChoiceValue::Int(0))
            .unwrap_err();
        assert_eq!(err.contract_name, "first");
    }

    #[test]
    fn contract_set_admits_succeeds_when_all_admit() {
        // Built-in productive-bundle contracts have no opinion on
        // a `trim.v_3` candidate against an empty `Simulation` (no
        // refdata, no V allele assigned) — `admits_typed` defaults
        // to `Ok` along every relevant path.
        let s = ContractSet::new()
            .with(Box::new(AnchorPreserved::new(Segment::V)))
            .with(Box::new(ProductiveJunctionFrame::new()))
            .with(Box::new(NoStopCodonInJunction::new()));

        let sim = Simulation::new();
        assert!(s
            .admits_typed(
                &sim,
                None,
                ChoiceContext::none().with_address(ChoiceAddress::Trim {
                    segment: crate::address::VdjSegment::V,
                    end: crate::assignment::TrimEnd::Three,
                }),
                &ChoiceValue::Int(0),
            )
            .is_ok());
    }
}
