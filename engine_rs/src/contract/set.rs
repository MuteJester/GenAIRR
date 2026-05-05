//! `ContractSet` — composition of multiple contracts (D.5).

use crate::ir::Simulation;
use crate::refdata::RefDataConfig;
use crate::trace::ChoiceValue;

use super::{ChoiceContext, Contract, ContractViolation};

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
pub struct ContractSet {
    contracts: Vec<Box<dyn Contract>>,
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
        self.contracts.push(contract);
        self
    }

    /// In-place append. Returns `&mut self` for builder-mut chains.
    pub fn add(&mut self, contract: Box<dyn Contract>) -> &mut Self {
        self.contracts.push(contract);
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

    /// Test whether `candidate` at `address` is admissible by every
    /// contract. Short-circuits on the first violator (sampling
    /// hot path). Returns the violator's `ContractViolation`.
    pub fn admits(
        &self,
        sim: &Simulation,
        refdata: Option<&RefDataConfig>,
        address: &str,
        candidate: &ChoiceValue,
    ) -> Result<(), ContractViolation> {
        for c in &self.contracts {
            c.admits(sim, refdata, address, candidate)?;
        }
        Ok(())
    }

    /// Context-aware variant of [`ContractSet::admits`].
    pub fn admits_with_context(
        &self,
        sim: &Simulation,
        refdata: Option<&RefDataConfig>,
        address: &str,
        candidate: &ChoiceValue,
        context: ChoiceContext,
    ) -> Result<(), ContractViolation> {
        for c in &self.contracts {
            c.admits_with_context(sim, refdata, address, candidate, context)?;
        }
        Ok(())
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
        make_assembled_sim, make_assembled_sim_from_refdata, make_v_anchor_at,
        make_vj_for_frame_test, make_vj_with_anchor_codons,
    };
    use super::super::{
        AnchorPreserved, NoStopCodonInJunction, ProductiveJunctionFrame,
    };
    use super::*;
    use crate::assignment::AlleleInstance;
    use crate::ir::Segment;
    use crate::refdata::AlleleId;

    #[test]
    fn contract_set_empty_admits_everything() {
        let s = ContractSet::new();
        assert!(s.is_empty());
        assert_eq!(s.len(), 0);

        let sim = Simulation::new();
        assert!(s.verify(&sim, None).is_ok());
        assert!(s
            .admits(&sim, None, "any.address", &ChoiceValue::Int(42))
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
        let cfg = make_vj_for_frame_test(Some(6), Some(0));
        let v_inst = AlleleInstance::new(AlleleId::new(0));
        let j_inst = AlleleInstance::new(AlleleId::new(0));
        let sim = make_assembled_sim(0, 9, 9, 6, v_inst, j_inst);

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
        // Two contracts, both have admits returning Err. The set
        // should return the first one's violation.
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
            fn admits(
                &self,
                _sim: &Simulation,
                _refdata: Option<&RefDataConfig>,
                _address: &str,
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
            .admits(&sim, None, "any.address", &ChoiceValue::Int(0))
            .unwrap_err();
        assert_eq!(err.contract_name, "first");
    }

    #[test]
    fn contract_set_admits_succeeds_when_all_admit() {
        // Default `admits` is Ok, so a set of default contracts
        // admits everything.
        let s = ContractSet::new()
            .with(Box::new(AnchorPreserved::new(Segment::V)))
            .with(Box::new(ProductiveJunctionFrame::new()))
            .with(Box::new(NoStopCodonInJunction::new()));

        let sim = Simulation::new();
        assert!(s
            .admits(&sim, None, "trim.v_3", &ChoiceValue::Int(0))
            .is_ok());
    }
}
