use crate::address;
use crate::dist::FilteredSampleError;
use crate::ir::{NucHandle, Simulation, SimulationBuilder};
use crate::pass::{Pass, PassContext, PassError};
use crate::trace::ChoiceValue;

use super::event::WeightedS5FMutationEvent;
use super::S5FMutationPass;

impl S5FMutationPass {
    pub(super) fn execute_with_sampling_mode(
        &self,
        sim: &Simulation,
        ctx: &mut PassContext,
        strict: bool,
    ) -> Result<Simulation, PassError> {
        let count_raw = self.count_dist.sample(ctx.rng);
        if strict && count_raw < 0 {
            return Err(PassError::invalid_distribution_output(
                self.name(),
                address::MUTATE_S5F_COUNT,
                count_raw,
                "negative_count",
            ));
        }
        if strict && count_raw > u32::MAX as i64 {
            return Err(PassError::invalid_distribution_output(
                self.name(),
                address::MUTATE_S5F_COUNT,
                count_raw,
                "count_exceeds_u32",
            ));
        }
        assert!(
            count_raw >= 0,
            "S5FMutationPass: count distribution returned negative {}",
            count_raw
        );
        assert!(
            count_raw <= u32::MAX as i64,
            "S5FMutationPass: count distribution returned {} > u32::MAX",
            count_raw
        );
        ctx.trace
            .record(address::MUTATE_S5F_COUNT, ChoiceValue::Int(count_raw));

        let count = count_raw as u32;
        if count == 0 || sim.pool.len() < 5 {
            return Ok(sim.clone());
        }

        // Phase 7: route mutations through `SimulationBuilder` so
        // each `change_base` notifies attached observers.
        //
        // Phase 10: attach BOTH codon-rail (one per existing region)
        // AND walker (one per existing V/D/J region) observers.
        // After the mutation loop, the combined seal helper stages
        // walker calls in V→D→J order with monotonically increasing
        // `evidence_version`s so the post-pass
        // `PassEffect::EditBases` refresh hits the Phase-1 fast path
        // for each segment — skipping the from-scratch
        // `call_from_region` rebuild entirely.
        let mut builder = SimulationBuilder::from_simulation(sim.clone());
        builder.attach_codon_rail_observers_for_all_regions();
        if let Some(ref_index) = ctx.reference_index {
            builder.attach_dirty_signal_observer();
            // Walker observers require a SegmentRefIndex per segment.
            // When the runtime didn't compile a reference index (test
            // harness / refdata-only plans), skip walker attachment —
            // the post-pass `PassEffect::EditBases` refresh handles
            // it just as it did pre-Phase-10.
            for region in builder
                .peek()
                .sequence
                .regions
                .iter()
                .cloned()
                .collect::<Vec<_>>()
                .iter()
                .filter(|r| ref_index.get(r.segment).is_some())
            {
                let segment_index = ref_index.get(region.segment).unwrap();
                builder.attach_walker_observer_for_region(segment_index, region);
            }
        }

        for i in 0..count {
            let profile = self.build_profile(&builder.peek().pool);
            if profile.is_empty() {
                break;
            }

            let total: f64 = profile.iter().map(|(_, m)| m).sum();
            if total <= 0.0 || !total.is_finite() {
                break;
            }
            let base_address = address::mutate_s5f_base(i);

            if let Some(contracts) = ctx.contracts {
                let candidates = match self.event_candidates(builder.peek(), &profile, i, count) {
                    Ok(candidates) => candidates,
                    Err(reason) if strict => {
                        return Err(PassError::constraint_sampling(
                            self.name(),
                            &base_address,
                            reason,
                        ));
                    }
                    Err(_) => Vec::new(),
                };
                let filtered: Vec<WeightedS5FMutationEvent> = candidates
                    .into_iter()
                    .filter(|candidate| {
                        contracts
                            .admits_with_context(
                                builder.peek(),
                                ctx.refdata,
                                &base_address,
                                &ChoiceValue::Base(candidate.event.base),
                                candidate.context,
                            )
                            .is_ok()
                    })
                    .collect();

                match Self::sample_weighted_event(ctx.rng, &filtered) {
                    Ok(Some(event)) => {
                        ctx.trace.record(
                            address::mutate_s5f_site(i),
                            ChoiceValue::Int(event.site as i64),
                        );
                        ctx.trace
                            .record(base_address, ChoiceValue::Base(event.base));

                        builder.change_base(NucHandle::new(event.site), event.base);
                        continue;
                    }
                    Ok(None) if strict => {
                        return Err(PassError::constraint_sampling(
                            self.name(),
                            &base_address,
                            FilteredSampleError::EmptyAdmissibleSupport,
                        ));
                    }
                    Err(reason) if strict => {
                        return Err(PassError::constraint_sampling(
                            self.name(),
                            &base_address,
                            reason,
                        ));
                    }
                    Ok(None) | Err(_) => {}
                }
            }

            let r = ctx.rng.next_f64() * total;
            let mut cum = 0.0;
            let mut chosen_pos = profile[0].0;
            for &(pos, mu) in &profile {
                cum += mu;
                if r < cum {
                    chosen_pos = pos;
                    break;
                }
            }

            let context = match Self::context_at(&builder.peek().pool, chosen_pos) {
                Some(c) => c,
                None => continue,
            };
            let row = self.kernel.substitution_row(context);
            let candidates = Self::positive_row_candidates(row);
            let chosen_base =
                match Self::sample_weighted_base(ctx.rng, &candidates).map_err(|reason| {
                    PassError::constraint_sampling(self.name(), &base_address, reason)
                })? {
                    Some(base) => base,
                    None => continue,
                };

            ctx.trace.record(
                address::mutate_s5f_site(i),
                ChoiceValue::Int(chosen_pos as i64),
            );
            ctx.trace
                .record(base_address, ChoiceValue::Base(chosen_base));

            builder.change_base(NucHandle::new(chosen_pos), chosen_base);
        }

        // Seal observers and commit results. When a reference index
        // is available (compiled execution path), use the combined
        // helper that also stages walker calls so the post-pass
        // `PassEffect::EditBases` refresh fast-paths each segment.
        // Otherwise fall back to the codon-rail-only helper —
        // walker observers weren't attached (no segment_index to
        // borrow), so they have nothing to seal anyway.
        let mut sealed = if let Some(ref_index) = ctx.reference_index {
            builder.seal_with_committed_live_calls(ref_index)
        } else {
            builder.seal_with_committed_codon_rails()
        };

        // Phase 17: stash the mutation count on `LiveCallState` so
        // the AIRR projection reads `n_mutations` in O(1) instead of
        // trace-scanning `MUTATE_S5F_COUNT`. Preserves prior
        // semantics — we stash `count_raw` to match what the trace
        // recorded above (line 45), which is the requested count
        // rather than the actual applied count if the loop broke
        // early.
        if count_raw > 0 {
            let mut state = sealed
                .live_calls
                .as_ref()
                .map(|s| (**s).clone())
                .unwrap_or_default();
            state.mutation_count = state.mutation_count.saturating_add(count_raw as u32);
            sealed = sealed.with_live_calls(state);
        }
        Ok(sealed)
    }
}
