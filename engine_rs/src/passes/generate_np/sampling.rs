use crate::address::ChoiceAddress;
use crate::contract::{ChoiceContext, JunctionStopState};
use crate::dist::{
    sample_base_with_admit_mask, sample_filtered_with_policy, Distribution, EmptySupport,
    FilteredSampleError,
};
use crate::ir::Simulation;
use crate::pass::{PassContext, PassError};
use crate::rng::Rng;
use crate::trace::ChoiceValue;

use super::GenerateNPPass;

/// Adapter that exposes a pre-materialised `(byte, weight)`
/// support vector as a `Distribution<Output = u8>`. Used to
/// feed a Markov / position-conditional generator's per-call
/// support through the existing generic helpers
/// ([`sample_base_with_admit_mask`] / [`sample_filtered_with_policy`])
/// without duplicating the inverse-CDF logic.
///
/// Owns the pairs by reference so each per-position call avoids
/// a second `Vec` clone — the lifetime is bounded by the
/// caller's local pair-vector.
struct SupportPairsDist<'a> {
    pairs: &'a [(u8, f64)],
}

impl<'a> SupportPairsDist<'a> {
    fn new(pairs: &'a [(u8, f64)]) -> Self {
        Self { pairs }
    }
}

impl Distribution for SupportPairsDist<'_> {
    type Output = u8;

    fn sample(&self, rng: &mut Rng) -> u8 {
        // Inverse-CDF over the pair vector — same RNG shape
        // (one `next_f64()`) as `UniformBase::sample` / the
        // generic helpers, so the unconstrained Markov path
        // consumes exactly one RNG word per base.
        let total: f64 = self.pairs.iter().map(|(_, w)| w).sum();
        // Spec-layer validation guarantees positive total; the
        // defensive assert here surfaces a malformed generator
        // before silently looping forever.
        assert!(
            total.is_finite() && total > 0.0,
            "SupportPairsDist::sample: support has non-positive total weight {total}",
        );
        let r = rng.next_f64() * total;
        let mut cum = 0.0;
        for (b, w) in self.pairs.iter() {
            cum += *w;
            if r < cum {
                return *b;
            }
        }
        // ULP fallback.
        self.pairs.last().unwrap().0
    }

    fn support(&self) -> Option<Vec<(u8, f64)>> {
        Some(self.pairs.to_vec())
    }
}

/// v3.0 empty-support policy for NP-length sampling. An NP region
/// with length 0 is a valid no-op: skipping the region entirely.
const NP_LENGTH_EMPTY_SUPPORT: EmptySupport<i64> = EmptySupport::Sentinel(0);

/// v3.0 empty-support policy for NP-base sampling. `b'N'`
/// translates to amino acid `X` (not a stop), and NP slots sit
/// between V/J anchor codons so the IUPAC sentinel never
/// violates anchor preservation either.
const NP_BASE_EMPTY_SUPPORT: EmptySupport<u8> = EmptySupport::Sentinel(b'N');

/// Resolve [`NP_BASE_EMPTY_SUPPORT`] to its sentinel byte for the
/// fast-path (admit-mask) caller. The fast path uses
/// `sample_base_with_admit_mask`, whose error shape differs from
/// `sample_filtered_result`, so it can't go through
/// `sample_filtered_with_policy` directly; resolving the const
/// here keeps the policy declaration in one place.
#[inline]
fn apply_base_empty_support_policy() -> u8 {
    match NP_BASE_EMPTY_SUPPORT {
        EmptySupport::Sentinel(s) => s,
        EmptySupport::Skip => unreachable!(
            "NP base sampler requires a sentinel — NP slots can't be skipped, the pool needs a byte"
        ),
    }
}

impl GenerateNPPass {
    /// Constraint-aware NP-length sample.
    ///
    /// Empty-support policy lives at [`NP_LENGTH_EMPTY_SUPPORT`]:
    /// strict mode surfaces `PassError::ConstraintSampling`,
    /// permissive mode emits the length-0 sentinel (the only
    /// architectural no-op for a length sampler).
    pub(super) fn sample_length(
        &self,
        sim: &Simulation,
        ctx: &mut PassContext,
        address: &'static str,
        strict: bool,
    ) -> Result<i64, PassError> {
        let refdata = ctx.refdata;
        let contracts = ctx.contracts;

        if let Some(contracts) = contracts {
            let pass_name = self.pass_name().to_string();
            let context =
                ChoiceContext::none().with_address_if_missing(ChoiceAddress::parse(address));
            let outcome = sample_filtered_with_policy(
                ctx.rng,
                self.length_dist.as_ref(),
                |&candidate: &i64| {
                    contracts
                        .admits_typed(sim, refdata, context, &ChoiceValue::Int(candidate))
                        .is_ok()
                },
                strict,
                &pass_name,
                address,
                NP_LENGTH_EMPTY_SUPPORT,
            )?;
            // The policy guarantees `Some(_)` in permissive (Skip
            // is never used here); strict mode would have errored
            // above.
            return Ok(outcome.expect("NP_LENGTH_EMPTY_SUPPORT is Sentinel(0), never Skip"));
        }
        Ok(self.length_dist.sample(ctx.rng))
    }

    /// Constraint-aware base sample.
    ///
    /// fast path: when `admit_mask` is `Some`, the caller has
    /// pre-computed (via the `ProductiveAdmitMaskObserver` attached to
    /// the `SimulationBuilder`) a 4-bit mask of which canonical bases
    /// admit at this slot. We sample directly through that mask via
    /// [`sample_base_with_admit_mask`], skipping the generic
    /// `sample_filtered_result` Vec materialisation and the
    /// per-candidate contract trait dispatch entirely. This is the
    /// architectural payoff of contracts move from
    /// "queried per candidate" to "subscribed once per slot."
    ///
    /// Slow path (no admit mask available — e.g. tests using
    /// `PassRuntime::execute_with_refdata` that doesn't compile a
    /// reference index): falls back to the
    /// `sample_filtered_result` + predicate-via-`admits_with_context`
    /// path.
    pub(super) fn sample_base(
        &self,
        sim: &Simulation,
        ctx: &mut PassContext,
        address: &str,
        index: u32,
        total_len: u32,
        strict: bool,
        junction_stop_state: Option<&JunctionStopState>,
        admit_mask: Option<u8>,
        previous: Option<u8>,
    ) -> Result<u8, PassError> {
        let refdata = ctx.refdata;
        let contracts = ctx.contracts;

        // Materialise per-position support via the generator.
        // Markov: row keyed by `previous`. Uniform / categorical
        // wrappers ignore `previous` / `position` and return the
        // pre-slice 4-way support.
        let support_pairs = self.base_generator.support(index as usize, previous);
        let base_dist = SupportPairsDist::new(&support_pairs);

        // ── Trace-injected replay (Tier 3 NP base) ───────────────
        //
        // Consume the recorded base and validate it against the same
        // support path the sampler would have followed at this slot:
        //
        //   * Fast path (admit_mask present): the recorded canonical
        //     base must be admitted by the 4-bit mask. The `N`
        //     sentinel is accepted only when the mask is empty —
        //     that's the policy-emission branch in fresh sampling.
        //   * Slow path (contracts active, no admit mask): canonical
        //     base must be admitted by `contracts.admits_typed`. `N`
        //     is accepted only when every canonical candidate in
        //     the generator's per-position support is rejected.
        //   * No-contracts path: byte must be in the generator's
        //     per-position support if enumerable.
        //
        // Any other byte (including `N` outside the empty-support
        // emission path) surfaces as `ConstraintSampling`. We never
        // force-apply a recorded value the live sampler couldn't
        // have produced.
        if let Some(cursor) = ctx.replay_cursor.as_deref_mut() {
            let parsed = ChoiceAddress::parse(address)
                .expect("NP base address is a built-in ChoiceAddress");
            let recorded = cursor
                .expect_base(parsed)
                .map_err(|reason| PassError::replay(self.pass_name(), reason))?;
            self.validate_replayed_np_base(
                recorded,
                sim,
                refdata,
                contracts,
                index,
                total_len,
                address,
                junction_stop_state,
                admit_mask,
                &support_pairs,
            )?;
            return Ok(recorded);
        }

        // fast path: admit-mask observer hands us a 4-bit
        // mask, we sample inverse-CDF directly over the admitted
        // subset of the generator's per-position support. No
        // per-candidate contract dispatch, no intermediate
        // `filtered` Vec.
        //
        // We still need the contract set to be present — if the user
        // didn't request `respect=productive()` the JunctionStopState
        // wouldn't have been built and `admit_mask` would be None
        // anyway, so contracts.is_some() is guaranteed here.
        // Fast and slow paths share the same v3 empty-support
        // policy, [`NP_BASE_EMPTY_SUPPORT`]. We resolve the
        // policy inline here for the fast path (different error
        // shape than `sample_filtered_result`) and via
        // `sample_filtered_with_policy` for the slow path.
        if let Some(mask) = admit_mask {
            match sample_base_with_admit_mask(ctx.rng, &base_dist, mask) {
                Ok(value) => return Ok(value),
                Err(reason) if strict => {
                    return Err(self.constraint_sampling_error(address, reason));
                }
                Err(_) => return Ok(apply_base_empty_support_policy()),
            }
        } else if let Some(contracts) = contracts {
            // Slow path: per-candidate contract dispatch via
            // `sample_filtered_with_policy`. Reached only when no
            // admit-mask observer is attached (test harnesses,
            // refdata-only plans).
            let mut context = ChoiceContext::indexed(index, total_len)
                .with_address_if_missing(ChoiceAddress::parse(address));
            if let Some(state) = junction_stop_state {
                context = context.with_junction_stop_state(state);
            }
            let pass_name = self.pass_name().to_string();
            let outcome = sample_filtered_with_policy(
                ctx.rng,
                &base_dist,
                |candidate: &u8| {
                    contracts
                        .admits_typed(sim, refdata, context, &ChoiceValue::Base(*candidate))
                        .is_ok()
                },
                strict,
                &pass_name,
                address,
                NP_BASE_EMPTY_SUPPORT,
            )?;
            return Ok(outcome.expect("NP_BASE_EMPTY_SUPPORT is Sentinel(b'N'), never Skip"));
        }

        Ok(base_dist.sample(ctx.rng))
    }

    pub(super) fn constraint_sampling_error(
        &self,
        address: impl Into<String>,
        reason: FilteredSampleError,
    ) -> PassError {
        PassError::constraint_sampling(self.pass_name(), address, reason)
    }

    /// Replay-mode validator for a recorded NP base. Mirrors the
    /// three branches inside [`Self::sample_base`] so the
    /// admissibility rule matches what the sampler would have
    /// applied at this slot.
    ///
    /// Returns `Ok(())` if the byte is reproducible by the live
    /// sampler, otherwise `Err(PassError::ConstraintSampling)` —
    /// never force-applies.
    #[allow(clippy::too_many_arguments)]
    fn validate_replayed_np_base(
        &self,
        recorded: u8,
        sim: &Simulation,
        refdata: Option<&crate::refdata::RefDataConfig>,
        contracts: Option<&crate::contract::ContractSet>,
        index: u32,
        total_len: u32,
        address: &str,
        junction_stop_state: Option<&JunctionStopState>,
        admit_mask: Option<u8>,
        support_pairs: &[(u8, f64)],
    ) -> Result<(), PassError> {
        let policy_sentinel = apply_base_empty_support_policy();
        let reject = || {
            Err(self.constraint_sampling_error(
                address.to_string(),
                FilteredSampleError::EmptyAdmissibleSupport,
            ))
        };

        if let Some(mask) = admit_mask {
            // Fast path: admit-mask present.
            if recorded == policy_sentinel {
                // `N` sentinel is only valid when the mask is empty
                // (the policy-emission branch).
                return if mask == 0 { Ok(()) } else { reject() };
            }
            // Canonical recorded byte: must be admitted by the mask.
            let bit_idx = match recorded {
                b'A' | b'a' => 0,
                b'C' | b'c' => 1,
                b'G' | b'g' => 2,
                b'T' | b't' => 3,
                _ => return reject(), // non-canonical, non-sentinel
            };
            if (mask >> bit_idx) & 1 == 1 {
                Ok(())
            } else {
                reject()
            }
        } else if let Some(contracts) = contracts {
            // Slow path: contracts active, admit_mask not precomputed.
            let mut context = ChoiceContext::indexed(index, total_len)
                .with_address_if_missing(ChoiceAddress::parse(address));
            if let Some(state) = junction_stop_state {
                context = context.with_junction_stop_state(state);
            }
            if recorded == policy_sentinel {
                // Valid only when every candidate in the
                // generator's per-position support is rejected
                // by contracts — same predicate
                // `sample_filtered_with_policy` uses to decide it
                // exhausted the support and must emit the sentinel.
                if support_pairs.is_empty() {
                    return Ok(());
                }
                let any_admitted = support_pairs.iter().any(|(b, _)| {
                    contracts
                        .admits_typed(sim, refdata, context, &ChoiceValue::Base(*b))
                        .is_ok()
                });
                return if any_admitted { reject() } else { Ok(()) };
            }
            // Canonical recorded byte: contracts must admit it.
            if contracts
                .admits_typed(sim, refdata, context, &ChoiceValue::Base(recorded))
                .is_ok()
            {
                Ok(())
            } else {
                reject()
            }
        } else {
            // No-contracts path: byte must be in the
            // generator's per-position support. The generator
            // always returns enumerable support (uniform /
            // categorical / Markov all do), so there's no
            // "best effort accept" fallback to consider.
            if support_pairs.iter().any(|(b, _)| *b == recorded) {
                Ok(())
            } else {
                reject()
            }
        }
    }
}
