//! Shared helpers for contract-aware pass sampling.
//!
//! These helpers keep pass implementations focused on event order and
//! trace shape while centralizing the permissive-vs-strict behavior for
//! constrained candidate draws.

use crate::contract::ChoiceContext;
use crate::dist::{sample_filtered_result, Distribution, FilteredSampleError};
use crate::ir::{NucHandle, Simulation};
use crate::pass::{PassContext, PassError};
use crate::trace::ChoiceValue;

fn identity_base(base: u8) -> u8 {
    base
}

/// Metadata for a base draw that will replace one existing pool
/// nucleotide.
#[derive(Copy, Clone)]
pub(crate) struct TargetedBaseChoice<'a> {
    pub(crate) pass_name: &'a str,
    pub(crate) address: &'a str,
    pub(crate) index: u32,
    pub(crate) count: u32,
    pub(crate) site: NucHandle,
    pub(crate) strict: bool,
    value_transform: fn(u8) -> u8,
}

impl<'a> TargetedBaseChoice<'a> {
    pub(crate) fn new(
        pass_name: &'a str,
        address: &'a str,
        index: u32,
        count: u32,
        site: NucHandle,
        strict: bool,
    ) -> Self {
        Self {
            pass_name,
            address,
            index,
            count,
            site,
            strict,
            value_transform: identity_base,
        }
    }

    pub(crate) fn with_value_transform(mut self, value_transform: fn(u8) -> u8) -> Self {
        self.value_transform = value_transform;
        self
    }

    fn transform(&self, base: u8) -> u8 {
        (self.value_transform)(base)
    }

    fn context(&self) -> ChoiceContext<'static> {
        ChoiceContext::targeted_base_substitution(self.index, self.count, self.site)
    }
}

/// Sample a destination base for a targeted substitution.
///
/// If contracts are active, enumerable distributions are filtered
/// before a value is committed. Permissive execution falls back to the
/// legacy unconstrained draw when filtering cannot produce a value;
/// strict execution returns a structured `PassError`.
pub(crate) fn sample_targeted_base(
    sim: &Simulation,
    ctx: &mut PassContext,
    base_dist: &dyn Distribution<Output = u8>,
    choice: TargetedBaseChoice<'_>,
) -> Result<u8, PassError> {
    let refdata = ctx.refdata;
    let contracts = ctx.contracts;

    if let Some(contracts) = contracts {
        let context = choice.context();
        match sample_filtered_result(ctx.rng, base_dist, |candidate: &u8| {
            let candidate = choice.transform(*candidate);
            contracts
                .admits_with_context(
                    sim,
                    refdata,
                    choice.address,
                    &ChoiceValue::Base(candidate),
                    context,
                )
                .is_ok()
        }) {
            Ok(value) => return Ok(choice.transform(value)),
            Err(reason) if choice.strict => {
                return Err(PassError::constraint_sampling(
                    choice.pass_name,
                    choice.address,
                    reason,
                ));
            }
            Err(_) => {
                // Permissive legacy path: if support is missing or
                // filtered empty, preserve the old unconstrained draw.
            }
        }
    }

    Ok(choice.transform(base_dist.sample(ctx.rng)))
}

/// One finite structural candidate plus the IR that would result if
/// it were committed.
#[derive(Clone)]
pub(crate) struct PostEventCandidate<'a, T> {
    pub(crate) value: T,
    pub(crate) weight: f64,
    pub(crate) post_sim: Simulation,
    pub(crate) context: ChoiceContext<'a>,
}

impl<'a, T> PostEventCandidate<'a, T> {
    pub(crate) fn new(
        value: T,
        weight: f64,
        post_sim: Simulation,
        context: ChoiceContext<'a>,
    ) -> Self {
        Self {
            value,
            weight,
            post_sim,
            context,
        }
    }
}

/// Sample a structural event after filtering by whole-IR post-state.
///
/// Returns `Ok(None)` when there are no active contracts or when a
/// permissive caller should fall back to legacy unconstrained sampling.
/// Strict callers receive the same structured `PassError` shape used by
/// substitution sampling.
pub(crate) fn sample_contract_verified_event<'a, T: Clone>(
    pre_sim: &Simulation,
    ctx: &mut PassContext,
    pass_name: &str,
    address: &str,
    strict: bool,
    candidates: Result<Vec<PostEventCandidate<'a, T>>, FilteredSampleError>,
) -> Result<Option<T>, PassError> {
    let Some(contracts) = ctx.contracts else {
        return Ok(None);
    };

    let candidates = match candidates {
        Ok(candidates) => candidates,
        Err(reason) if strict => {
            return Err(PassError::constraint_sampling(pass_name, address, reason));
        }
        Err(_) => return Ok(None),
    };

    let mut filtered: Vec<PostEventCandidate<'a, T>> = candidates
        .into_iter()
        .filter(|candidate| {
            contracts
                .admits_post_event(
                    pre_sim,
                    &candidate.post_sim,
                    ctx.refdata,
                    address,
                    candidate.context,
                )
                .is_ok()
        })
        .collect();

    if filtered.is_empty() {
        if strict {
            return Err(PassError::constraint_sampling(
                pass_name,
                address,
                FilteredSampleError::EmptyAdmissibleSupport,
            ));
        }
        return Ok(None);
    }

    let total: f64 = filtered.iter().map(|candidate| candidate.weight).sum();
    if !total.is_finite() || total <= 0.0 {
        if strict {
            return Err(PassError::constraint_sampling(
                pass_name,
                address,
                FilteredSampleError::InvalidFilteredSupport,
            ));
        }
        return Ok(None);
    }

    let r = ctx.rng.next_f64() * total;
    let mut cum = 0.0;
    for candidate in &filtered {
        cum += candidate.weight;
        if r < cum {
            return Ok(Some(candidate.value.clone()));
        }
    }

    Ok(Some(
        filtered
            .pop()
            .expect("filtered candidates checked non-empty")
            .value,
    ))
}
