//! v3.0 tuple sampler for the indel pass.
//!
//! ## Why a tuple sampler
//!
//! Per the architectural rule, contracts narrow the action space
//! before the pass proposes. The
//! [`ProductiveJunctionFrame`](crate::contract::ProductiveJunctionFrame)
//! contract requires the junction length to stay divisible by 3 —
//! every single indel in or around the junction changes that
//! length by ±1, so **no single-event indel is admissible** under
//! the productive bundle. Indels in productive immune repertoires
//! come in compensating tuples (e.g. one insertion + one deletion).
//!
//! ## How
//!
//! The sampler:
//! 1. At each step `k`, **rebuilds** a mod-3 reachability DP from
//!    the per-step delta-mass distribution of the *current
//!    working sim*. `dp[m][r]` is the mass of `m`-step
//!    continuations whose frame-delta sum is `≡ r (mod 3)`.
//! 2. Samples one event per step using **continuation weighting**:
//!    each candidate's effective weight is the product of its
//!    natural per-event weight and the DP mass of valid
//!    `count − k − 1`-step completions from the new accumulated
//!    delta.
//! 3. After all `count` events are sampled hypothetically, runs
//!    [`ContractSet::admits_post_event`] on the assembled
//!    post-state against the full bundle. The DP handles **frame
//!    only**; this final check catches non-frame violations the
//!    bundle still enforces (e.g.
//!    [`NoStopCodonInJunction`](crate::contract::NoStopCodonInJunction)
//!    when an inserted base completes a stop codon).
//! 4. On post-event rejection, the tuple is resampled up to
//!    [`POST_EVENT_RETRY_BUDGET`] times. On success the tuple is
//!    applied through the [`MutationTransaction`]. If every retry
//!    fails (or some step had no admissible candidates),
//!    permissive mode collapses the entire batch to `NoOp`
//!    records (the count slot exists in the trace but no
//!    structural change happens) and strict mode surfaces
//!    `PassError::ConstraintSampling`. **The retry budget reduces
//!    strict false-errors when valid support exists but the first
//!    sample missed it; it does NOT prove the support is empty
//!    when exhausted** — that's a residual approximation of the
//!    sampler.
//!
//! ## Per-step DP recompute
//!
//! Continuation mass at step `k` is computed from a DP seeded by
//! the per-step delta-mass profile of the **current working sim**
//! at step `k`, not step 0. This matches the v3.0 invariant: the
//! continuation weight reflects the natural distribution at the
//! step where the candidate is being evaluated. A candidate whose
//! commit would change the future per-step support carries the
//! same multiplier as siblings with the same mod-3 delta — the
//! approximation `mass(step k+1) ≈ mass(step k+2) ≈ ...` is the
//! only remaining bias and it's controlled (the per-step mass
//! profile changes by ±1 site per commit; for biological count
//! values of ≤ ~6 the residual error is negligible). True
//! per-candidate DP would seed from the *post-event* sim of each
//! candidate; for count ≤ 6 the cost is bearable but the
//! information gain is small and we defer that refinement.

use crate::address;
use crate::contract::{ChoiceContext, ContractSet, IndelEventClass, IndelKindHint};
use crate::dist::FilteredSampleError;
use crate::ir::{flag, NucHandle, Nucleotide, Simulation};
use crate::pass::{Pass, PassContext, PassError};
use crate::refdata::RefDataConfig;

use super::event::IndelEvent;
use super::IndelPass;

/// Maximum number of tuple-sampling attempts before falling back
/// to no-op (permissive) / surfacing an error (strict). Each
/// attempt is one full `count`-step sequential sample + one
/// final `admits_post_event` check. The budget **reduces strict
/// false-errors** when valid bundle support exists but the first
/// draw missed it (cross-event no-stop interactions, which the
/// mod-3 DP doesn't model); it does **not** prove the support is
/// empty when exhausted. A truly exhaustive strict guarantee
/// would require a proper rejection-sampling bound or
/// enumeration — out of scope for v3.0.
const POST_EVENT_RETRY_BUDGET: u32 = 16;

/// A per-event candidate evaluated against the bundle's
/// indel-class composition. `delta` ∈ {−1, 0, +1}.
#[derive(Clone, Copy, Debug)]
struct TupleCandidate {
    event: IndelEvent,
    delta: i32,
    natural_weight: f64,
}

/// Outcome of a single tuple-sampling attempt.
enum TupleSampleOutcome {
    /// A frame-balanced tuple was sampled and the assembled
    /// post-state satisfies the bundle's full admit check.
    Accepted(Vec<IndelEvent>),
    /// Some step had no admissible candidates (or the
    /// `tuple_candidates` enumeration returned `Err`). Permissive
    /// no-op; nothing to retry.
    Empty,
    /// A frame-balanced tuple was sampled, but the assembled
    /// post-state failed the bundle's `admits_post_event` check
    /// (cross-event interaction the mod-3 DP can't model). Caller
    /// should retry up to [`POST_EVENT_RETRY_BUDGET`] times.
    PostEventRejected,
}

/// Triangle of mod-3 reachability masses. `dp[m][r]` is the mass
/// of m-step continuations whose frame-delta sum is `r (mod 3)`.
type Mod3Dp = Vec<[f64; 3]>;

impl IndelPass {
    /// Sample a frame-balanced tuple of length `count`, then run
    /// the bundle's full post-event validator on the assembled
    /// post-state. Returns `Some(events)` when a valid tuple was
    /// found, `None` when permissive mode skipped the entire
    /// batch. Strict-mode failures bubble out as `PassError`.
    ///
    /// On post-event rejection (cross-event interactions that the
    /// mod-3 DP doesn't model — e.g. an insertion completing a
    /// `TAA` codon under [`NoStopCodonInJunction`]), the tuple is
    /// resampled up to [`POST_EVENT_RETRY_BUDGET`] times before
    /// declaring the sample-space exhausted. The retry budget
    /// prevents a strict-mode false-error when the bundle support
    /// is non-empty but the first sample happened to land on a
    /// cross-event violator.
    pub(super) fn sample_admissible_tuple(
        &self,
        sim: &Simulation,
        ctx: &mut PassContext,
        count: u32,
        strict: bool,
    ) -> Result<Option<Vec<IndelEvent>>, PassError> {
        let contracts = match ctx.contracts {
            Some(c) => c,
            None => return Ok(None),
        };
        let refdata = ctx.refdata;

        for _attempt in 0..POST_EVENT_RETRY_BUDGET {
            match self.sample_one_tuple(sim, ctx, contracts, refdata, count, strict)? {
                TupleSampleOutcome::Accepted(events) => return Ok(Some(events)),
                TupleSampleOutcome::Empty => return Ok(None),
                TupleSampleOutcome::PostEventRejected => {
                    // Cross-event interaction (typically a stop
                    // codon introduced by an inserted base) — the
                    // mod-3 DP doesn't model this. Resample.
                    continue;
                }
            }
        }

        // Budget exhausted. In strict mode this surfaces as
        // a constraint-sampling error; in permissive mode we
        // accept the architectural reality that the bundle's
        // *full* support couldn't be exercised within the budget
        // and collapse to no-op.
        if strict {
            return Err(PassError::constraint_sampling(
                self.name(),
                address::ChoiceAddress::CorruptIndelSite(count.saturating_sub(1)).to_string(),
                FilteredSampleError::EmptyAdmissibleSupport,
            ));
        }
        Ok(None)
    }

    /// Inner workhorse: sample a single tuple attempt. Returns
    /// `PostEventRejected` when the assembled tuple fails the
    /// bundle's `admits_post_event` check (the caller retries).
    fn sample_one_tuple(
        &self,
        sim: &Simulation,
        ctx: &mut PassContext,
        contracts: &ContractSet,
        refdata: Option<&RefDataConfig>,
        count: u32,
        strict: bool,
    ) -> Result<TupleSampleOutcome, PassError> {
        let mut events: Vec<IndelEvent> = Vec::with_capacity(count as usize);
        let mut working = sim.clone();
        let mut accumulated_mod3: i32 = 0;

        for k in 0..count {
            let remaining = (count - k - 1) as usize;

            // **Per-step DP recompute**: the continuation mass at
            // step k reflects the per-step delta-mass profile of
            // the *current* working sim. Stationarity across
            // remaining steps is the only residual approximation;
            // tied to the v3.0 invariant for sequential conditional
            // sampling with continuation weighting.
            let candidates = match self.tuple_candidates(&working, contracts, refdata) {
                Ok(cs) => cs,
                Err(reason) if strict => {
                    return Err(PassError::constraint_sampling(
                        self.name(),
                        address::ChoiceAddress::CorruptIndelSite(k).to_string(),
                        reason,
                    ));
                }
                Err(_) => Vec::new(),
            };
            let mass = step_mass(&candidates);
            let dp = compute_mod3_dp(remaining as u32, mass);

            // Effective weight at step k = natural × continuation
            // mass for the remaining-step target. `dp[remaining]`
            // is the mass vector over `count-k-1` remaining steps.
            let weighted: Vec<(IndelEvent, i32, f64)> = candidates
                .into_iter()
                .filter_map(|c| {
                    let after = accumulated_mod3 + c.delta;
                    let needed = (-after).rem_euclid(3) as usize;
                    let cont = dp[remaining][needed];
                    let eff = c.natural_weight * cont;
                    if eff.is_finite() && eff > 0.0 {
                        Some((c.event, c.delta, eff))
                    } else {
                        None
                    }
                })
                .collect();

            if weighted.is_empty() {
                if strict {
                    return Err(PassError::constraint_sampling(
                        self.name(),
                        address::ChoiceAddress::CorruptIndelSite(k).to_string(),
                        FilteredSampleError::EmptyAdmissibleSupport,
                    ));
                }
                return Ok(TupleSampleOutcome::Empty);
            }

            let total: f64 = weighted.iter().map(|(_, _, w)| *w).sum();
            if !total.is_finite() || total <= 0.0 {
                if strict {
                    return Err(PassError::constraint_sampling(
                        self.name(),
                        address::ChoiceAddress::CorruptIndelSite(k).to_string(),
                        FilteredSampleError::InvalidFilteredSupport,
                    ));
                }
                return Ok(TupleSampleOutcome::Empty);
            }
            let r = ctx.rng.next_f64() * total;
            let mut cum = 0.0;
            let mut chosen_event = weighted.last().expect("non-empty checked above").0;
            let mut chosen_delta = weighted.last().expect("non-empty checked above").1;
            for &(event, delta, w) in &weighted {
                cum += w;
                if r < cum {
                    chosen_event = event;
                    chosen_delta = delta;
                    break;
                }
            }

            // Hypothetically apply to a local working sim — final
            // commit through the TX is deferred until after the
            // full post-event validator.
            working = apply_event_hypothetical(&working, chosen_event);
            events.push(chosen_event);
            accumulated_mod3 = (accumulated_mod3 + chosen_delta).rem_euclid(3);
        }

        // Final post-event validator: the DP handled frame only.
        // The bundle's `admits_post_event` aggregates every other
        // contract's verdict on the assembled post-state (stops in
        // the junction, anchor preservation against actual bytes,
        // etc.). On rejection we signal the caller to retry — the
        // first sampled tuple may have hit a cross-event
        // interaction (e.g. a stop introduced by an inserted base)
        // even though *other* frame-balanced tuples in the support
        // satisfy the bundle.
        let validator_context =
            ChoiceContext::indel_insertion(count.saturating_sub(1), count, NucHandle::new(0));
        if contracts
            .admits_post_event(sim, &working, refdata, validator_context)
            .is_err()
        {
            return Ok(TupleSampleOutcome::PostEventRejected);
        }

        Ok(TupleSampleOutcome::Accepted(events))
    }

    /// Enumerate per-event candidates with their bundle-composed
    /// indel class + natural per-event weight. `Forbidden`
    /// candidates are filtered out.
    fn tuple_candidates(
        &self,
        sim: &Simulation,
        contracts: &ContractSet,
        refdata: Option<&RefDataConfig>,
    ) -> Result<Vec<TupleCandidate>, FilteredSampleError> {
        let pool_len = sim.pool.len() as u32;
        let mut candidates: Vec<TupleCandidate> = Vec::new();

        if self.insertion_prob > 0.0 {
            let base_support = self.base_support()?;
            let per_site = self.insertion_prob / (pool_len + 1) as f64;
            for site in 0..=pool_len {
                let class = contracts.admissible_indel_class_at(
                    sim,
                    refdata,
                    site,
                    IndelKindHint::Insertion,
                );
                if matches!(class, IndelEventClass::Forbidden) {
                    continue;
                }
                let delta = class.frame_delta().unwrap_or(0) as i32;
                for &(base, base_w) in &base_support {
                    candidates.push(TupleCandidate {
                        event: IndelEvent::Insertion { site, base },
                        delta,
                        natural_weight: per_site * base_w,
                    });
                }
            }
        }

        let deletion_prob = 1.0 - self.insertion_prob;
        if deletion_prob > 0.0 {
            if pool_len == 0 {
                // Pool-empty deletion is a no-op; classification is
                // FrameNeutral by definition (nothing to remove).
                candidates.push(TupleCandidate {
                    event: IndelEvent::Deletion { site: None },
                    delta: 0,
                    natural_weight: deletion_prob,
                });
            } else {
                let per_site = deletion_prob / pool_len as f64;
                for site in 0..pool_len {
                    let class = contracts.admissible_indel_class_at(
                        sim,
                        refdata,
                        site,
                        IndelKindHint::Deletion,
                    );
                    if matches!(class, IndelEventClass::Forbidden) {
                        continue;
                    }
                    let delta = class.frame_delta().unwrap_or(0) as i32;
                    candidates.push(TupleCandidate {
                        event: IndelEvent::Deletion { site: Some(site) },
                        delta,
                        natural_weight: per_site,
                    });
                }
            }
        }

        Ok(candidates)
    }
}

/// Aggregate per-delta mass from a candidate set. Returns an array
/// `[mass(d=0), mass(d=+1), mass(d=−1)]` keyed by `delta.rem_euclid(3)`
/// → `[0, 1, 2]`. Unnormalized — the DP propagates the mass
/// proportionally.
fn step_mass(candidates: &[TupleCandidate]) -> [f64; 3] {
    let mut mass = [0.0f64; 3];
    for c in candidates {
        let r = c.delta.rem_euclid(3) as usize;
        mass[r] += c.natural_weight;
    }
    mass
}

/// `dp[m][r]` = mass of m-step continuations whose frame-delta
/// sum is `r (mod 3)`. Rebuilt at each step of
/// [`IndelPass::sample_one_tuple`] from the per-step delta-mass
/// profile of the current working sim — the only stationarity
/// assumed is that the next `m` remaining steps share the
/// **current step's** mass profile (a strictly better
/// approximation than freezing the profile at step 0).
fn compute_mod3_dp(count: u32, per_step_mass: [f64; 3]) -> Mod3Dp {
    let count = count as usize;
    let mut dp = vec![[0.0f64; 3]; count + 1];
    // Base case: 0 steps to go ⇒ only the "already at 0" target
    // has mass 1 (multiplicative identity).
    dp[0][0] = 1.0;
    for m in 1..=count {
        for r in 0..3 {
            // dp[m][r] = sum_{d=0..3} per_step_mass[d] * dp[m-1][(r - d) mod 3]
            let mut acc = 0.0f64;
            for d in 0..3 {
                let prev_r = ((r as i32 - d as i32).rem_euclid(3)) as usize;
                acc += per_step_mass[d] * dp[m - 1][prev_r];
            }
            dp[m][r] = acc;
        }
    }
    dp
}

/// Apply an indel event to a sim without going through the TX.
/// Used by the tuple sampler to build a hypothetical chain.
/// `pub(super)` so the replay-mode validator in `replay.rs` can
/// build the same hypothetical post-state and run `admits_post_event`
/// against it before committing through the TX.
pub(super) fn apply_event_hypothetical(sim: &Simulation, event: IndelEvent) -> Simulation {
    match event {
        IndelEvent::Insertion { site, base } => {
            let segment = IndelPass::insertion_segment(sim, site);
            let nuc = Nucleotide::synthetic(base, segment, flag::INDEL_INSERTED);
            sim.with_indel_inserted(site, nuc)
        }
        IndelEvent::Deletion { site: Some(site) } => sim.with_indel_deleted(site),
        IndelEvent::Deletion { site: None } | IndelEvent::NoOp => sim.clone(),
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn mod3_dp_base_case_is_identity() {
        // 0 steps to go: only target 0 has mass 1.
        let dp = compute_mod3_dp(0, [1.0, 0.0, 0.0]);
        assert_eq!(dp[0], [1.0, 0.0, 0.0]);
    }

    #[test]
    fn mod3_dp_one_step_neutral_only_keeps_target_zero() {
        // Per-step delta is always 0 (mass concentrated at d=0).
        // After m steps, only target 0 has nonzero mass.
        let dp = compute_mod3_dp(5, [1.0, 0.0, 0.0]);
        for m in 0..=5 {
            assert!(dp[m][0] > 0.0);
            assert_eq!(dp[m][1], 0.0);
            assert_eq!(dp[m][2], 0.0);
        }
    }

    #[test]
    fn mod3_dp_two_steps_with_balanced_pm1_can_return_to_zero() {
        // Per-step mass: equal probability on +1 and -1. With 2
        // steps, P(sum ≡ 0) > 0 (via +1 then -1 or vice versa).
        let dp = compute_mod3_dp(2, [0.0, 0.5, 0.5]);
        assert!(dp[2][0] > 0.0);
    }

    #[test]
    fn mod3_dp_three_steps_with_only_plus_one_reaches_zero() {
        // Three +1 steps land at 3 ≡ 0 (mod 3).
        let dp = compute_mod3_dp(3, [0.0, 1.0, 0.0]);
        assert!(dp[3][0] > 0.0);
        // But 1 or 2 steps of +1 cannot land at 0.
        assert_eq!(dp[1][0], 0.0);
        assert_eq!(dp[2][0], 0.0);
    }

    #[test]
    fn mod3_dp_one_step_only_plus_one_lands_at_one() {
        let dp = compute_mod3_dp(1, [0.0, 1.0, 0.0]);
        assert_eq!(dp[1][0], 0.0);
        assert_eq!(dp[1][1], 1.0);
        assert_eq!(dp[1][2], 0.0);
    }
}
