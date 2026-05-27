use crate::address::ChoiceAddress;
use crate::ir::{translate_codon, NucHandle, Simulation, AMINO_STOP};
use crate::junction::compute_junction;
use crate::refdata::RefDataConfig;
use crate::trace::ChoiceValue;

use super::NoStopCodonInJunction;
use crate::contract::{ChoiceContext, ChoiceKind, Contract, ContractViolation};

impl NoStopCodonInJunction {
    /// Classify a [`ChoiceContext`] as a targeted base substitution.
    ///
    /// Two signals are accepted:
    /// 1. Explicit kind tag set by the substitution helper:
    ///    [`ChoiceKind::TargetedBaseSubstitution`] — built-in
    ///    callers via [`ChoiceContext::targeted_base_substitution`]
    ///    flow through this branch.
    /// 2. Typed [`ChoiceAddress`] variant matching one of the
    ///    targeted-substitution producer addresses — picks up
    ///    legacy string callers routed through the
    ///    [`Contract::admits`] trait default, which attaches the
    ///    parsed address but leaves `kind = Plain`.
    ///
    /// Note: `CorruptQualityBase` is intentionally excluded from
    /// the address-variant list to match the pre-typed-migration
    /// behavior — quality substitutions reach the contract via
    /// the explicit kind tag from the substitution helper, not
    /// via the address-variant fallback.
    fn is_targeted_substitution_candidate(context: ChoiceContext<'_>) -> bool {
        if context.kind == ChoiceKind::TargetedBaseSubstitution {
            return true;
        }
        matches!(
            context.address,
            Some(
                ChoiceAddress::MutateUniformBase(_)
                    | ChoiceAddress::MutateS5fBase(_)
                    | ChoiceAddress::CorruptPcrBase(_)
                    | ChoiceAddress::CorruptContaminantBase(_)
            )
        )
    }

    fn pool_base_with_candidate(
        sim: &Simulation,
        handle: NucHandle,
        target: NucHandle,
        candidate_base: u8,
    ) -> Option<u8> {
        if handle == target {
            Some(candidate_base)
        } else {
            sim.pool.get(handle).map(|n| n.base)
        }
    }

    pub(super) fn admits_targeted_substitution_candidate(
        &self,
        sim: &Simulation,
        refdata: Option<&RefDataConfig>,
        candidate: &ChoiceValue,
        context: ChoiceContext<'_>,
    ) -> Result<(), ContractViolation> {
        if !Self::is_targeted_substitution_candidate(context) {
            return Ok(());
        }

        let candidate_base = match candidate {
            ChoiceValue::Base(b) => *b,
            _ => return Ok(()),
        };
        let target = match context.target {
            Some(h) => h,
            None => return Ok(()),
        };
        let refdata = match refdata {
            None => return Ok(()),
            Some(r) => r,
        };
        let junction = match compute_junction(sim, refdata) {
            None => return Ok(()),
            Some(j) => j,
        };

        if !junction.is_in_frame() {
            return Ok(());
        }

        let start = junction.start.index();
        let end = junction.end.index();
        let target_idx = target.index();
        if target_idx < start || target_idx >= end {
            return Ok(());
        }

        let mut pos = start;
        while pos + 3 <= end {
            let b1 =
                Self::pool_base_with_candidate(sim, NucHandle::new(pos), target, candidate_base);
            let b2 = Self::pool_base_with_candidate(
                sim,
                NucHandle::new(pos + 1),
                target,
                candidate_base,
            );
            let b3 = Self::pool_base_with_candidate(
                sim,
                NucHandle::new(pos + 2),
                target,
                candidate_base,
            );

            if let (Some(b1), Some(b2), Some(b3)) = (b1, b2, b3) {
                if translate_codon(b1, b2, b3) == AMINO_STOP {
                    let diagnostic_address = context.address_string().unwrap_or_default();
                    return Err(ContractViolation::new(
                        self.name(),
                        format!(
                            "candidate at {} for target {} would leave stop codon '{}{}{}' at junction position {}",
                            diagnostic_address,
                            target_idx,
                            b1 as char,
                            b2 as char,
                            b3 as char,
                            pos
                        ),
                    ));
                }
            }
            pos += 3;
        }

        Ok(())
    }
}
