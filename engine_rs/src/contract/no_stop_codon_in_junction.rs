//! `NoStopCodonInJunction` — junction codons must translate to amino acids.

use crate::ir::{NucHandle, Simulation};
use crate::refdata::RefDataConfig;
use crate::trace::ChoiceValue;

use super::{ChoiceContext, Contract, ContractKind, ContractViolation};

mod np_candidate;
mod support;
mod targeted;
mod verify;

/// Verifies that no codon inside the junction translates to a stop
/// (TAA, TAG, TGA). Walks codons from `junction.start` in steps of
/// 3 and checks each translation.
///
/// Vacuously satisfied (returns `Ok`) when:
/// - no `RefDataConfig` is provided,
/// - the junction is undefined (`compute_junction` returns `None`),
/// - the junction is out of frame — that's a different contract's
///   concern (`ProductiveJunctionFrame`); we pass on rather than
///   double-reporting.
///
/// Returns a `ContractViolation` with the first stop's position
/// and codon string when a stop is detected.
///
/// **Note:** the contract reads from `sim.pool`, not from the
/// allele reference. Mutations / NP base substitutions are
/// reflected — this is the *current* pool's stop-codon status, not
/// the germline's.
pub struct NoStopCodonInJunction;

impl NoStopCodonInJunction {
    pub fn new() -> Self {
        Self
    }
}

impl Default for NoStopCodonInJunction {
    fn default() -> Self {
        Self::new()
    }
}

impl Contract for NoStopCodonInJunction {
    fn name(&self) -> &str {
        "no_stop_codon_in_junction"
    }

    fn kind(&self) -> ContractKind {
        ContractKind::NoStopCodonInJunction
    }

    fn verify(
        &self,
        sim: &Simulation,
        refdata: Option<&RefDataConfig>,
    ) -> Result<(), ContractViolation> {
        self.verify_impl(sim, refdata)
    }

    fn admits_typed(
        &self,
        sim: &Simulation,
        refdata: Option<&RefDataConfig>,
        context: ChoiceContext<'_>,
        candidate: &ChoiceValue,
    ) -> Result<(), ContractViolation> {
        // Fast path: when the caller wired the [`JunctionStopState`]
        // precomputed admit-mask observer onto the context, the
        // per-candidate stop check is O(1) — we only re-translate
        // the codon touched by the candidate. Falls through to the
        // generic `admits_np_candidate` slow path when no state is
        // available (e.g. test harnesses that don't compile a
        // reference index).
        if let Some(state) = context.junction_stop_state {
            if let Some((segment, addr_index)) = Self::np_base_for_context(context) {
                if let ChoiceValue::Base(candidate_base) = *candidate {
                    let index = context.draw_index.unwrap_or(addr_index);
                    state.admits_np_candidate(sim, segment, index, candidate_base)?;
                    return self
                        .admits_targeted_substitution_candidate(sim, refdata, candidate, context);
                }
            }
        }

        self.admits_np_candidate(sim, refdata, candidate, context)?;
        self.admits_targeted_substitution_candidate(sim, refdata, candidate, context)
    }

    /// v3.0 constrain-before-propose: per-site admissible-base mask.
    ///
    /// Returns the bitmask of canonical bases (A/C/G/T) that, when
    /// substituted at `site`, would not introduce a stop codon
    /// anywhere in the junction.
    ///
    /// Mirrors `admits_targeted_substitution_candidate` semantics
    /// exactly for the TargetedBaseSubstitution case:
    /// - No refdata / no junction / out-of-frame junction / out-of-junction
    ///   site → all four bases admissible (the contract has no opinion).
    /// - In-junction site: walk junction codons in steps of 3, build the
    ///   hypothetical codon with each candidate substituted at `site`,
    ///   set the bit only for candidates that produce no stop.
    fn admissible_bases_at(
        &self,
        sim: &Simulation,
        refdata: Option<&RefDataConfig>,
        site: NucHandle,
    ) -> super::BaseMask {
        self.admissible_bases_at_impl(sim, refdata, site)
    }

    /// Pinned non-canonical write: walk junction codons with `byte`
    /// substituted at `site`, reject if any resulting codon
    /// translates to a stop. `translate_codon` is case-insensitive,
    /// so lowercase quality writes (`a`/`c`/`g`/`t`) and uppercase
    /// candidates yield the same verdict. Non-canonical bytes
    /// (`N`, IUPAC ambiguities) translate to `X` rather than a
    /// stop, so they are admitted by this contract.
    fn admits_fixed_base_at(
        &self,
        sim: &Simulation,
        refdata: Option<&RefDataConfig>,
        site: NucHandle,
        byte: u8,
    ) -> bool {
        self.admits_fixed_base_at_impl(sim, refdata, site, byte)
    }
}

#[cfg(test)]
mod tests;
