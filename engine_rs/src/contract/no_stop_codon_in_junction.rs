//! `NoStopCodonInJunction` — junction codons must translate to amino acids.

use crate::ir::{translate_codon, NucHandle, Simulation, AMINO_STOP};
use crate::junction::compute_junction;
use crate::refdata::RefDataConfig;
use crate::trace::ChoiceValue;

use super::{ChoiceContext, Contract, ContractKind, ContractViolation};

mod np_candidate;
mod targeted;

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
        let refdata = match refdata {
            None => return Ok(()),
            Some(r) => r,
        };

        let junction = match compute_junction(sim, refdata) {
            None => return Ok(()),
            Some(j) => j,
        };

        // Out of frame → not our concern. ProductiveJunctionFrame
        // is the contract that owns the frame check; here we'd just
        // be walking codon-shaped triples that aren't real codons.
        if !junction.is_in_frame() {
            return Ok(());
        }

        let start = junction.start.index();
        let end = junction.end.index();

        let mut pos = start;
        while pos + 3 <= end {
            // Defensive lookups — pool pointers from the junction
            // window should always be valid for a well-formed
            // simulation, but we don't want to panic on a malformed
            // input. Treat missing handles as ambiguous.
            let b1 = sim
                .pool
                .get(NucHandle::new(pos))
                .map(|n| n.base)
                .unwrap_or(b'N');
            let b2 = sim
                .pool
                .get(NucHandle::new(pos + 1))
                .map(|n| n.base)
                .unwrap_or(b'N');
            let b3 = sim
                .pool
                .get(NucHandle::new(pos + 2))
                .map(|n| n.base)
                .unwrap_or(b'N');

            let aa = translate_codon(b1, b2, b3);
            if aa == AMINO_STOP {
                return Err(ContractViolation::new(
                    self.name(),
                    format!(
                        "stop codon at junction position {}: '{}{}{}'",
                        pos, b1 as char, b2 as char, b3 as char
                    ),
                ));
            }
            pos += 3;
        }

        Ok(())
    }

    fn admits(
        &self,
        sim: &Simulation,
        refdata: Option<&RefDataConfig>,
        address: &str,
        candidate: &ChoiceValue,
    ) -> Result<(), ContractViolation> {
        self.admits_np_candidate(sim, refdata, address, candidate, ChoiceContext::none())
    }

    fn admits_with_context(
        &self,
        sim: &Simulation,
        refdata: Option<&RefDataConfig>,
        address: &str,
        candidate: &ChoiceValue,
        context: ChoiceContext<'_>,
    ) -> Result<(), ContractViolation> {
        if let Some(state) = context.junction_stop_state {
            if let Some((segment, addr_index)) = Self::parse_np_base_address(address) {
                if let ChoiceValue::Base(candidate_base) = *candidate {
                    let index = context.draw_index.unwrap_or(addr_index);
                    state.admits_np_candidate(sim, segment, index, candidate_base)?;
                    return self.admits_targeted_substitution_candidate(
                        sim, refdata, address, candidate, context,
                    );
                }
            }
        }

        self.admits_np_candidate(sim, refdata, address, candidate, context)?;
        self.admits_targeted_substitution_candidate(sim, refdata, address, candidate, context)
    }
}

#[cfg(test)]
mod tests;
