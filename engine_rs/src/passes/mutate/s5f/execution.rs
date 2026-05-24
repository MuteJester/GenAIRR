use crate::address;
use crate::dist::FilteredSampleError;
use crate::ir::{NucHandle, Simulation};
use crate::pass::{Pass, PassContext, PassError};
use crate::passes::mutation_transaction::MutationTransaction;
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

        // Per-step S5F sampling. Each iteration: build the per-
        // position S5F-weighted profile from the current pool, then
        // either (a) enumerate candidates + contract-filter +
        // weighted-sample when contracts are active, or (b) sample
        // the position from the profile + sample the destination
        // base from the kernel's row. All pool mutations route
        // through MutationTransaction so observer attach, contract-
        // checked vs unconditional substitution, the
        // seal_with_committed_live_calls fast path, and the
        // mutation_count sidecar all happen automatically.
        let mut tx = MutationTransaction::open(sim, ctx, self.name(), strict);

        for i in 0..count {
            let profile = self.build_profile(&tx.peek().pool);
            if profile.is_empty() {
                break;
            }

            let total: f64 = profile.iter().map(|(_, m)| m).sum();
            if total <= 0.0 || !total.is_finite() {
                break;
            }
            let base_address = address::mutate_s5f_base(i);

            // Contract-filtered candidate enumeration first; falls
            // through to the kernel-row path when contracts are
            // absent or filter exhausted.
            let (sim_ref, ctx_ref) = tx.split_borrows();
            let contracts_opt = ctx_ref.contracts;
            if let Some(contracts) = contracts_opt {
                let candidates = match self.event_candidates(sim_ref, &profile, i, count) {
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
                                sim_ref,
                                ctx_ref.refdata,
                                &base_address,
                                &ChoiceValue::Base(candidate.event.base),
                                candidate.context,
                            )
                            .is_ok()
                    })
                    .collect();

                match Self::sample_weighted_event(ctx_ref.rng, &filtered) {
                    Ok(Some(event)) => {
                        ctx_ref.trace.record(
                            address::mutate_s5f_site(i),
                            ChoiceValue::Int(event.site as i64),
                        );
                        ctx_ref.trace.record(base_address, ChoiceValue::Base(event.base));
                        tx.substitute_base_fixed(NucHandle::new(event.site), event.base)?;
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

            // Unconstrained path: kernel-row sampling.
            let r = tx.rng().next_f64() * total;
            let mut cum = 0.0;
            let mut chosen_pos = profile[0].0;
            for &(pos, mu) in &profile {
                cum += mu;
                if r < cum {
                    chosen_pos = pos;
                    break;
                }
            }

            let context = match Self::context_at(&tx.peek().pool, chosen_pos) {
                Some(c) => c,
                None => continue,
            };
            let row = self.kernel.substitution_row(context);
            let candidates = Self::positive_row_candidates(row);
            let chosen_base =
                match Self::sample_weighted_base(tx.rng(), &candidates).map_err(|reason| {
                    PassError::constraint_sampling(self.name(), &base_address, reason)
                })? {
                    Some(base) => base,
                    None => continue,
                };

            tx.trace().record(
                address::mutate_s5f_site(i),
                ChoiceValue::Int(chosen_pos as i64),
            );
            tx.trace().record(base_address, ChoiceValue::Base(chosen_base));
            tx.substitute_base_fixed(NucHandle::new(chosen_pos), chosen_base)?;
        }

        // `count_raw` matches what the trace records above — the
        // requested count, not the actual applied count if the
        // sampling loop broke early. Preserves prior semantics.
        if count_raw > 0 {
            tx.add_to_mutation_count(count_raw as u32);
        }
        tx.commit()
    }
}
