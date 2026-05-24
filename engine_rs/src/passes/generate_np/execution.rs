use crate::contract::JunctionStopState;
use crate::ir::{flag, NucHandle, Nucleotide, Region, Simulation, SimulationBuilder};
use crate::pass::{PassContext, PassError};
use crate::trace::ChoiceValue;


use super::GenerateNPPass;

impl GenerateNPPass {
    pub(super) fn execute_with_sampling_mode(
        &self,
        sim: &Simulation,
        ctx: &mut PassContext,
        strict: bool,
    ) -> Result<Simulation, PassError> {
        let address = self.length_address();
        let length = self.sample_length(sim, ctx, address, strict)?;
        if strict && length < 0 {
            return Err(PassError::invalid_distribution_output(
                self.pass_name(),
                address,
                length,
                "negative_length",
            ));
        }
        if strict && length > u32::MAX as i64 {
            return Err(PassError::invalid_distribution_output(
                self.pass_name(),
                address,
                length,
                "length_exceeds_u32",
            ));
        }
        assert!(
            length >= 0,
            "GenerateNPPass({}): length distribution returned negative value {}",
            address,
            length
        );
        assert!(
            length <= u32::MAX as i64,
            "GenerateNPPass({}): length distribution returned {} > u32::MAX",
            address,
            length
        );
        ctx.trace.record(address, ChoiceValue::Int(length));

        let length = length as u32;
        let region_start = NucHandle::new(sim.pool.len() as u32);

        let junction_stop_state = match (ctx.refdata, ctx.contracts) {
            (Some(refdata), Some(_)) => {
                JunctionStopState::build(sim, refdata, self.np_segment, length)
            }
            _ => None,
        };

        // Frame phase is determined by upstream regions' cumulative
        // length mod 3 — captured here so the codon-rail observer
        // can be initialised before the push loop.
        let cumulative_len: u32 = sim.sequence.regions.iter().map(|r| r.len()).sum();
        let frame_phase = (cumulative_len % 3) as u8;

        let mut builder = SimulationBuilder::from_simulation(sim.clone());

        // Attach the productive-admit-mask observer when a
        // JunctionStopState was built. It caches the 4-bit admit
        // mask for the current NP slot via the precomputed contract
        // state; the sampler reads `builder.current_admit_mask()`
        // before each draw and bypasses `sample_filtered_result`'s
        // per-candidate contract dispatch via
        // `sample_base_with_admit_mask`.
        if let Some(ref state) = junction_stop_state {
            builder.attach_admit_mask_observer(state, self.np_segment);
        }

        for i in 0..length {
            let base_address = format!("{}[{}]", self.bases_prefix(), i);
            // Pull the current admit mask BEFORE sampling so the
            // sampler can route through the fast path. The mask
            // reflects every previously-committed NP byte (via
            // `builder.peek()` inside `current_admit_mask`).
            let admit_mask = builder.current_admit_mask();
            let base = self.sample_base(
                builder.peek(),
                ctx,
                &base_address,
                i,
                length,
                strict,
                junction_stop_state.as_ref(),
                admit_mask,
            )?;
            ctx.trace.record(base_address, ChoiceValue::Base(base));
            builder.push_nucleotide(Nucleotide::synthetic(base, self.np_segment, flag::N_NUC));
        }

        let current = builder.seal();

        let region_end = NucHandle::new(current.pool.len() as u32);
        // The per-region codon rail (amino_acids, stop_codon_positions)
        // is not maintained in the hot path — no production consumer
        // reads it. Callers that need it compute on demand via
        // `Region::with_codon_rail_recomputed(&pool)`.
        let region = Region::new(self.np_segment, region_start, region_end)
            .with_frame_phase(frame_phase);
        Ok(current.with_region_added(region))
    }
}
