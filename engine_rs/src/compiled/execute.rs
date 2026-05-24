//! Transactional execution of a compiled `Schedule`.
//!
//! `execute_transactional` is the engine's run loop. It owns the
//! trace, RNG, revision stack, and event log; each pass either commits
//! its delta or aborts the run. Per-pass derived-state refresh is
//! delegated to a slice of [`EffectHook`] implementations registered
//! on the compiled simulator — today that's the
//! [`crate::live_call::LiveCallRefreshHook`], which used to live
//! inline here as `apply_live_call_updates`. `abort` packages a
//! `PassError` together with the already-committed prefix so the
//! caller can introspect partial state.

use super::ExecutionPolicy;
use crate::contract::ContractSet;
use crate::event::{EventRecord, StateSummary, TraceSpan};
use crate::feasibility::FeasibilityContext;
use crate::ir::Simulation;
use crate::live_call::ReferenceMatchIndex;
use crate::pass::{
    EffectHook, HookContext, NodeId, Outcome, PassContext, PassError, PassPlan,
};
use crate::refdata::RefDataConfig;
use crate::rng::Rng;
use crate::trace::Trace;

#[derive(Copy, Clone)]
pub(crate) struct ExecutionInputs<'a> {
    pub(crate) plan: &'a PassPlan,
    /// Topologically-sorted execution order. `execute_transactional`
    /// iterates `sorted_order` and indexes into `plan.passes()`; the
    /// vec's positional index becomes the runtime `pass_index` in
    /// `PassContext`.
    pub(crate) sorted_order: &'a [NodeId],
    /// Effect hooks invoked after each pass commits. The compiled
    /// artifact owns the underlying `Vec<Box<dyn EffectHook>>` and
    /// hands the executor a borrowed slice of trait objects.
    pub(crate) hooks: &'a [Box<dyn EffectHook>],
    pub(crate) refdata: Option<&'a RefDataConfig>,
    pub(crate) contracts: Option<&'a ContractSet>,
    pub(crate) feasibility: Option<&'a FeasibilityContext>,
    pub(crate) reference_index: Option<&'a ReferenceMatchIndex>,
    pub(crate) policy: ExecutionPolicy,
}

#[derive(Debug)]
pub(crate) struct ExecutionAbort {
    pub(crate) error: PassError,
    /// Already-committed prefix at the point of failure. Kept private:
    /// public runners expose only the structured pass error, while
    /// internal tests can assert transaction commit boundaries.
    #[cfg_attr(not(test), allow(dead_code))]
    pub(crate) committed: Outcome,
}

pub(crate) fn execute_transactional(
    inputs: ExecutionInputs<'_>,
    initial: Simulation,
    seed: u64,
) -> Result<Outcome, ExecutionAbort> {
    let mut trace = Trace::new();
    let mut rng = Rng::new(seed);

    let mut revisions: Vec<Simulation> = Vec::with_capacity(inputs.plan.len() + 1);
    let mut pass_names: Vec<String> = Vec::with_capacity(inputs.plan.len());
    let mut events: Vec<EventRecord> = Vec::with_capacity(inputs.plan.len());

    revisions.push(initial);

    for node in inputs.sorted_order {
        let pass = &inputs.plan.passes()[node.index()];
        let prev = revisions.last().expect("non-empty by construction");
        let pass_index = pass_names.len();
        let pass_name = pass.name().to_string();
        let effects = pass.effects();
        let trace_start = trace.len();
        let pre = StateSummary::from_simulation(prev);
        let mut trace_delta = Trace::new();
        let mut ctx = PassContext {
            trace: &mut trace_delta,
            rng: &mut rng,
            pass_index,
            refdata: inputs.refdata,
            contracts: inputs.contracts,
            feasibility: inputs.feasibility,
            reference_index: inputs.reference_index,
        };

        let mut next = match inputs.policy {
            ExecutionPolicy::Permissive => pass.execute(prev, &mut ctx),
            ExecutionPolicy::Strict => match pass.execute_checked(prev, &mut ctx) {
                Ok(next) => next,
                Err(error) => {
                    return Err(abort(error, revisions, pass_names, trace, events));
                }
            },
        };
        drop(ctx);

        let hook_ctx = HookContext {
            reference_index: inputs.reference_index,
        };
        for hook in inputs.hooks {
            next = hook.apply(next, &effects, hook_ctx);
        }

        if inputs.policy == ExecutionPolicy::Strict {
            if let Some(contracts) = inputs.contracts {
                if let Err(violations) = contracts.verify(&next, inputs.refdata) {
                    return Err(abort(
                        PassError::contract_violation(pass_name, violations),
                        revisions,
                        pass_names,
                        trace,
                        events,
                    ));
                }
            }
        }

        let trace_end = trace_start + trace_delta.len();
        let post = StateSummary::from_simulation(&next);
        let event = EventRecord::pass_committed(
            pass_index,
            pass_name.clone(),
            effects,
            TraceSpan::new(trace_start, trace_end),
            pre,
            post,
        );

        trace.append_delta(trace_delta);
        pass_names.push(pass_name);
        revisions.push(next);
        events.push(event);
    }

    Ok(Outcome {
        revisions,
        pass_names,
        trace,
        events,
    })
}

fn abort(
    error: PassError,
    revisions: Vec<Simulation>,
    pass_names: Vec<String>,
    trace: Trace,
    events: Vec<EventRecord>,
) -> ExecutionAbort {
    ExecutionAbort {
        error,
        committed: Outcome {
            revisions,
            pass_names,
            trace,
            events,
        },
    }
}
