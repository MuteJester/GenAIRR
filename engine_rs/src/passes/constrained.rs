//! Shared helpers for contract-aware structural-event sampling.
//!
//! The substitution-pass helper (`sample_targeted_base` +
//! `TargetedBaseChoice`) is gone — all targeted base substitutions
//! now flow through [`crate::passes::mutation_transaction::MutationTransaction::substitute_base`].
//! The structural-event helper below is still used by the indel
//! pass (which evaluates each candidate by building a hypothetical
//! `post_sim` and asking contracts to verify whole-IR admissibility).

use crate::contract::ChoiceContext;
use crate::dist::FilteredSampleError;
use crate::ir::Simulation;
use crate::pass::{PassContext, PassError};

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
