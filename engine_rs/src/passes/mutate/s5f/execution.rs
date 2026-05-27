use std::collections::HashMap;

use crate::address;
use crate::contract::BaseMask;
use crate::dist::FilteredSampleError;
use crate::ir::{NucHandle, Simulation};
use crate::pass::{Pass, PassContext, PassError};
use crate::passes::count_source::sample_validated_count;
use crate::passes::mutation_transaction::MutationTransaction;
use crate::trace::ChoiceValue;

use super::event::WeightedS5FMutationEvent;
use super::S5FMutationPass;

/// Map a canonical base byte to its index in the S5F substitution
/// row `[A, C, G, T]`. Returns `None` for non-canonical bytes —
/// S5F never writes lowercase or `N`, so any non-A/C/G/T byte in a
/// replay cursor is by definition outside the S5F natural support.
fn base_to_row_index(b: u8) -> Option<usize> {
    match b {
        b'A' => Some(0),
        b'C' => Some(1),
        b'G' => Some(2),
        b'T' => Some(3),
        _ => None,
    }
}

impl S5FMutationPass {
    pub(super) fn execute_with_sampling_mode(
        &self,
        sim: &Simulation,
        ctx: &mut PassContext,
        strict: bool,
    ) -> Result<Simulation, PassError> {
        // Validation + trace recording shared with uniform / PCR /
        // quality / N-corrupt via `sample_validated_count`. Carries
        // trace-injected replay (Tier 2) for free — the helper
        // consumes the `mutate.s5f.count` record from the cursor
        // when one is active.
        let count = sample_validated_count(
            &self.count_source,
            ctx,
            sim.pool.len() as u32,
            self.name(),
            address::ChoiceAddress::MutateS5fCount,
            strict,
        )?;

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
        // Track realized mutations — slots may be skipped in
        // permissive mode when the contract bundle narrows the
        // admissible support to empty at every site. Per the v3.0
        // contract: count applied, not requested.
        let mut applied: u32 = 0;

        for i in 0..count {
            let site_choice_address = address::ChoiceAddress::MutateS5fSite(i);
            let base_choice_address = address::ChoiceAddress::MutateS5fBase(i);
            let base_address = base_choice_address.to_string();

            // ── Trace-injected replay (Tier 2, S5F) ──────────────
            //
            // Replay consumes the recorded `(site, base)` instead
            // of weighted-sampling. The recorded pair is validated
            // through three filters in order:
            //
            //   1. **Site in range** — must be a real pool index.
            //   2. **S5F naturality** — the (site, base) must be a
            //      reachable S5F event: site has a 5-mer context
            //      and the kernel row gives `base` strictly
            //      positive weight. The fresh-RNG path never
            //      writes outside this support; replay must hold
            //      the same line.
            //   3. **Contract admissibility** — under active
            //      contracts the recorded base must be admitted by
            //      the per-site mask at the in-progress sim.
            //
            // Failing any filter is a **structured error** rather
            // than the permissive skip used by the fresh-RNG
            // contract path. That asymmetry is intentional: a
            // recorded value that the engine now rejects is a
            // diagnostic worth surfacing — the original run did
            // produce it, so silently skipping would erase the
            // mismatch.
            let consumed_pair = if tx.replay_cursor().is_some() {
                let cursor = tx
                    .replay_cursor()
                    .expect("cursor presence checked above");
                let site = cursor
                    .expect_int(site_choice_address)
                    .map_err(|reason| PassError::replay(self.name(), reason))?;
                let base = cursor
                    .expect_base(base_choice_address)
                    .map_err(|reason| PassError::replay(self.name(), reason))?;
                Some((site, base))
            } else {
                None
            };
            if let Some((site_idx, base_byte)) = consumed_pair {
                let pool_len = tx.peek().pool.len() as u32;
                if site_idx < 0 || site_idx >= pool_len as i64 {
                    return Err(PassError::invalid_plan_state(
                        self.name().to_string(),
                        format!(
                            "S5F: replayed site {} out of pool range [0, {})",
                            site_idx, pool_len
                        ),
                    ));
                }
                let site = site_idx as u32;

                // S5F naturality: the recorded `(site, base)` must
                // be a reachable event under the current pool. The
                // fresh-RNG path requires BOTH a non-zero per-site
                // mutability (so `build_profile` includes the site)
                // AND a non-zero substitution-row weight for the
                // recorded base (so `sample_weighted_base` could
                // pick it). Either factor being zero means the
                // event was unreachable by the natural distribution
                // — replay must reject.
                let row_idx = base_to_row_index(base_byte);
                let s5f_natural_weight = match row_idx {
                    Some(idx) => match Self::context_at(&tx.peek().pool, site) {
                        Some(ctx) => {
                            let mu = self.kernel.mutability(ctx);
                            let row = self.kernel.substitution_row(ctx)[idx];
                            mu * row
                        }
                        None => 0.0,
                    },
                    None => 0.0,
                };
                if !(s5f_natural_weight > 0.0) {
                    return Err(PassError::constraint_sampling(
                        self.name(),
                        &base_address,
                        FilteredSampleError::EmptyAdmissibleSupport,
                    ));
                }

                // Contract admissibility: same shape as
                // `substitute_position_constrained`'s replay branch.
                let (sim_ref, ctx_ref) = tx.split_borrows();
                if let Some(contracts) = ctx_ref.contracts {
                    let mask = contracts.admissible_bases_at(
                        sim_ref,
                        ctx_ref.refdata,
                        NucHandle::new(site),
                    );
                    if !mask.admits(base_byte) {
                        return Err(PassError::constraint_sampling(
                            self.name(),
                            &base_address,
                            FilteredSampleError::EmptyAdmissibleSupport,
                        ));
                    }
                }

                ctx_ref
                    .trace
                    .record_choice(site_choice_address, ChoiceValue::Int(site as i64));
                ctx_ref
                    .trace
                    .record_choice(base_choice_address, ChoiceValue::Base(base_byte));
                tx.force_substitute_base(NucHandle::new(site), base_byte)?;
                applied += 1;
                continue;
            }

            let profile = self.build_profile(&tx.peek().pool);
            if profile.is_empty() {
                break;
            }

            let total: f64 = profile.iter().map(|(_, m)| m).sum();
            if total <= 0.0 || !total.is_finite() {
                break;
            }

            // v3.0 constrain-before-propose path.
            //
            // For each candidate `(site, base)` produced by the
            // natural S5F event enumeration with weight
            // `mutability(site) * (row_weight / row_total)`, look up
            // the per-site admissible-base mask from the active
            // contract bundle and zero the weight of any base the
            // mask rejects. Sampling from the surviving weighted
            // support produces draws from the **natural S5F event
            // distribution restricted to the contract-admissible
            // support** — exactly the architectural invariant.
            //
            // Per-site mask lookups are cached per slot: each site
            // has 4 candidate bases enumerated by `event_candidates`,
            // so caching collapses the 4 mask queries into one.
            let (sim_ref, ctx_ref) = tx.split_borrows();
            if let Some(contracts) = ctx_ref.contracts {
                // v3.0 architectural rule: under active contracts the
                // pass must NEVER propose an action outside the
                // contract-admissible support. If the mask-filtered
                // S5F candidate list is empty, the slot is genuinely
                // empty — in strict mode that surfaces as
                // `ConstraintSampling`; in permissive mode the slot
                // is consumed as a no-op (the count slot exists in
                // the trace via `count` but no `site[i]`/`base[i]`
                // entries are recorded and no pool mutation happens).
                // Falling through to the unconstrained kernel-row
                // sampler would re-propose a contract-violating
                // (site, base) and consume extra RNG words — both
                // architectural violations.
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
                let mut site_masks: HashMap<u32, BaseMask> = HashMap::new();
                let filtered: Vec<WeightedS5FMutationEvent> = candidates
                    .into_iter()
                    .filter(|candidate| {
                        let mask = *site_masks.entry(candidate.event.site).or_insert_with(|| {
                            contracts.admissible_bases_at(
                                sim_ref,
                                ctx_ref.refdata,
                                NucHandle::new(candidate.event.site),
                            )
                        });
                        mask.admits(candidate.event.base)
                    })
                    .collect();

                match Self::sample_weighted_event(ctx_ref.rng, &filtered) {
                    Ok(Some(event)) => {
                        ctx_ref.trace.record_choice(
                            site_choice_address,
                            ChoiceValue::Int(event.site as i64),
                        );
                        ctx_ref
                            .trace
                            .record_choice(base_choice_address, ChoiceValue::Base(event.base));
                        // Candidate was already admitted by the
                        // contract mask filter above — force-write.
                        tx.force_substitute_base(NucHandle::new(event.site), event.base)?;
                        applied += 1;
                    }
                    Ok(None) | Err(_) if strict => {
                        return Err(PassError::constraint_sampling(
                            self.name(),
                            &base_address,
                            FilteredSampleError::EmptyAdmissibleSupport,
                        ));
                    }
                    Ok(None) | Err(_) => {
                        // Permissive empty-support: skip the slot.
                        // No trace entry, no RNG consumption, no
                        // pool mutation — `applied` is not
                        // incremented.
                    }
                }
                continue;
            }

            // No active contracts: unconstrained kernel-row sampling.
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

            tx.trace()
                .record_choice(site_choice_address, ChoiceValue::Int(chosen_pos as i64));
            tx.trace()
                .record_choice(base_choice_address, ChoiceValue::Base(chosen_base));
            tx.force_substitute_base(NucHandle::new(chosen_pos), chosen_base)?;
            applied += 1;
        }

        // Bump the mutation-count sidecar by *applied* mutations,
        // not the requested count. When the contract bundle narrows
        // the admissible support to empty at some slots, permissive
        // mode skips those — the AIRR `n_mutations` field should
        // reflect the realized number of base changes.
        if applied > 0 {
            tx.add_to_mutation_count(applied);
        }
        tx.commit()
    }
}
