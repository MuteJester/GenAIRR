use crate::contract::{ChoiceContext, JunctionStopState};
use crate::dist::{sample_base_with_admit_mask, sample_filtered_result, FilteredSampleError};
use crate::ir::Simulation;
use crate::pass::{PassContext, PassError};
use crate::trace::ChoiceValue;

use super::GenerateNPPass;

impl GenerateNPPass {
    /// Constraint-aware length sample.
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
            match sample_filtered_result(ctx.rng, self.length_dist.as_ref(), |&candidate: &i64| {
                contracts
                    .admits(sim, refdata, address, &ChoiceValue::Int(candidate))
                    .is_ok()
            }) {
                Ok(value) => return Ok(value),
                Err(reason) if strict => {
                    return Err(self.constraint_sampling_error(address, reason));
                }
                Err(_) => {}
            }
        }
        Ok(self.length_dist.sample(ctx.rng))
    }

    /// Constraint-aware base sample.
    ///
    /// Phase 3 fast path: when `admit_mask` is `Some`, the caller has
    /// pre-computed (via the `ProductiveAdmitMaskObserver` attached to
    /// the `SimulationBuilder`) a 4-bit mask of which canonical bases
    /// admit at this slot. We sample directly through that mask via
    /// [`sample_base_with_admit_mask`], skipping the generic
    /// `sample_filtered_result` Vec materialisation and the
    /// per-candidate contract trait dispatch entirely. This is the
    /// architectural payoff of Phase 3: contracts move from
    /// "queried per candidate" to "subscribed once per slot."
    ///
    /// Slow path (no admit mask available — e.g. tests using
    /// `PassRuntime::execute_with_refdata` that doesn't compile a
    /// reference index): falls back to the existing
    /// `sample_filtered_result` + predicate-via-`admits_with_context`
    /// path. Behavior bit-identical to pre-Phase-3.
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
    ) -> Result<u8, PassError> {
        let refdata = ctx.refdata;
        let contracts = ctx.contracts;

        // Phase 3 fast path: admit-mask observer hands us a 4-bit
        // mask, we sample inverse-CDF directly over the admitted
        // subset of `base_dist`'s support. No per-candidate contract
        // dispatch, no intermediate `filtered` Vec.
        //
        // We still need the contract set to be present — if the user
        // didn't request `respect=productive()` the JunctionStopState
        // wouldn't have been built and `admit_mask` would be None
        // anyway, so contracts.is_some() is guaranteed here.
        if let Some(mask) = admit_mask {
            match sample_base_with_admit_mask(ctx.rng, self.base_dist.as_ref(), mask) {
                Ok(value) => return Ok(value),
                Err(reason) if strict => {
                    return Err(self.constraint_sampling_error(address, reason));
                }
                Err(_) => {}
            }
        } else if let Some(contracts) = contracts {
            // Slow path: per-candidate contract dispatch via
            // `sample_filtered_result`. Reached only when no
            // admit-mask observer is attached (test harnesses,
            // refdata-only plans).
            let mut context = ChoiceContext::indexed(index, total_len);
            if let Some(state) = junction_stop_state {
                context = context.with_junction_stop_state(state);
            }
            match sample_filtered_result(ctx.rng, self.base_dist.as_ref(), |candidate: &u8| {
                contracts
                    .admits_with_context(
                        sim,
                        refdata,
                        address,
                        &ChoiceValue::Base(*candidate),
                        context,
                    )
                    .is_ok()
            }) {
                Ok(value) => return Ok(value),
                Err(reason) if strict => {
                    return Err(self.constraint_sampling_error(address, reason));
                }
                Err(_) => {}
            }
        }

        Ok(self.base_dist.sample(ctx.rng))
    }

    pub(super) fn constraint_sampling_error(
        &self,
        address: impl Into<String>,
        reason: FilteredSampleError,
    ) -> PassError {
        PassError::constraint_sampling(self.pass_name(), address, reason)
    }
}
