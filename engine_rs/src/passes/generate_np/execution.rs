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
        let seq_start = sim.pool.len() as u32;

        let mut builder = SimulationBuilder::from_simulation(sim.clone());

        // Phase 2 extension: attach the codon-rail observer so the
        // NP region's amino_acids and stop_codon_positions accumulate
        // inline with the per-base push loop, just like the assembly
        // pass already does.
        builder.attach_codon_rail_observer(self.np_segment, seq_start, frame_phase);

        // Phase 3: attach the productive-admit-mask observer when a
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

        let sealed_rail = builder.seal_codon_rail_observer();
        let current = builder.seal();

        let region_end = NucHandle::new(current.pool.len() as u32);
        // Phase 2 extension: Region built directly from the sealed
        // codon-rail observer state — byte-identical to the
        // previous `Region::new(...).with_frame_phase(p)
        // .with_codon_rail_recomputed(&current.pool)` chain but
        // skipping the post-pass region rebuild.
        let region = Region::from_sealed_codon_rail(
            self.np_segment,
            region_start,
            region_end,
            frame_phase,
            sealed_rail,
        );
        Ok(current.with_region_added(region))
    }
}
