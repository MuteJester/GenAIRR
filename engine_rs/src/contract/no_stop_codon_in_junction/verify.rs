use crate::ir::{translate_codon, NucHandle, Simulation, AMINO_STOP};
use crate::junction::compute_junction;
use crate::refdata::RefDataConfig;

use super::NoStopCodonInJunction;
use crate::contract::{Contract, ContractViolation};

impl NoStopCodonInJunction {
    pub(super) fn verify_impl(
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

        // Out of frame -> not our concern. ProductiveJunctionFrame
        // owns the frame check.
        if !junction.is_in_frame() {
            return Ok(());
        }

        let start = junction.start.index();
        let end = junction.end.index();

        let mut pos = start;
        while pos + 3 <= end {
            // Defensive lookups: malformed pool pointers are treated
            // as ambiguous rather than panicking.
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
}
