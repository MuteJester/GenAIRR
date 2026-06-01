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
        // Trace-injected replay (length only — base draws still go
        // through the RNG path in this tier). Consume the recorded
        // length from the cursor; the downstream validation +
        // assertion path runs unchanged, so bad recorded values
        // still surface as `InvalidDistributionOutput` rather than
        // corrupting the IR.
        let length = if let Some(cursor) = ctx.replay_cursor.as_deref_mut() {
            cursor
                .expect_int(self.length_choice_address())
                .map_err(|reason| PassError::replay(self.pass_name(), reason))?
        } else {
            self.sample_length(sim, ctx, address, strict)?
        };
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
        ctx.trace
            .record_choice(self.length_choice_address(), ChoiceValue::Int(length));

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

        // Attach an event-log observer to the base-push builder
        // when a caller-supplied sink is present, so each NP
        // `push_nucleotide` fires a `BasePushed` event into the
        // pass's `EventRecord`. Drained right before the seal
        // below (which consumes the builder).
        if ctx.event_log_sink.is_some() {
            builder.attach_event_log_observer();
        }

        // Tracks the most recently emitted NP base across the
        // per-position loop so a Markov generator can condition
        // its `support()` on it. Uniform / categorical
        // generators ignore this. `None` for position 0; the
        // recorded base from the previous iteration otherwise.
        let mut previous: Option<u8> = None;
        for i in 0..length {
            let base_choice_address = crate::address::ChoiceAddress::NpBase {
                segment: self.typed_np_segment(),
                index: i,
            };
            let base_address = base_choice_address.to_string();
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
                previous,
            )?;
            ctx.trace
                .record_choice(base_choice_address, ChoiceValue::Base(base));
            builder.push_nucleotide(Nucleotide::synthetic(base, self.np_segment, flag::N_NUC));
            // Update only after the base is accepted/recorded
            // so a mid-stream replay error doesn't leak a
            // half-committed previous-base state.
            previous = Some(base);
        }

        // Drain the per-push events into the caller-supplied sink
        // before sealing (which consumes the builder).
        if ctx.event_log_sink.is_some() {
            let captured = builder.seal_event_log_observer();
            if let Some(sink) = ctx.event_log_sink.as_deref_mut() {
                sink.extend(captured);
            }
        }
        let current = builder.seal();

        let region_end = NucHandle::new(current.pool.len() as u32);
        // Region carries no codon-rail data; the pool is the
        // authoritative source. Callers that need translation call
        // `crate::ir::compute_codon_rail(&region, &pool)` on demand.
        let region =
            Region::new(self.np_segment, region_start, region_end).with_frame_phase(frame_phase);
        // Route the region append through a fresh builder so the
        // `RegionAdded` event flows to any attached sink. The
        // base-sampling builder above was already sealed; this
        // second builder exists solely to broadcast the consequence.
        //
        // Forward the captured event into `ctx.event_log_sink` when
        // the caller supplied one — same pass-level event-
        // observability primitive used by `AssembleSegmentPass`.
        let mut region_builder = SimulationBuilder::from_simulation(current);
        if ctx.event_log_sink.is_some() {
            region_builder.attach_event_log_observer();
        }
        region_builder.add_region(region);
        if let Some(sink) = ctx.event_log_sink.as_deref_mut() {
            sink.extend(region_builder.seal_event_log_observer());
        }
        Ok(region_builder.seal())
    }
}
