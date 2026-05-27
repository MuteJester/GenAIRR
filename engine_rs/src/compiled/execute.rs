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
use crate::ir::{Simulation, SimulationEvent};
use crate::live_call::ReferenceMatchIndex;
use crate::pass::{EffectHook, HookContext, NodeId, Outcome, PassContext, PassError, PassPlan};
use crate::refdata::RefDataConfig;
use crate::replay::TraceCursor;
use crate::rng::Rng;
use crate::trace::{ChoiceRecord, Trace};

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
    /// Pre-recorded trace to inject during execution. `Some(...)`
    /// puts the run into Option-B trace-injected replay: every
    /// sampling pass consumes records positionally from this slice
    /// instead of drawing from the RNG. After the plan completes,
    /// the cursor is asserted drained — extra trailing records
    /// surface as `PassError::Replay`.
    pub(crate) replay_records: Option<&'a [ChoiceRecord]>,
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

    // Replay-mode plumbing: when `inputs.replay_records` is `Some`,
    // build a cursor whose lifetime spans the whole run and reborrow
    // it into each `PassContext`. Passes that haven't been migrated
    // to the trace-injected path still draw from `rng` (this slice
    // wires only SampleAllelePass; other passes ignore the cursor).
    let mut cursor_opt: Option<TraceCursor> = inputs.replay_records.map(TraceCursor::new);

    let mut revisions: Vec<Simulation> = Vec::with_capacity(inputs.plan.len() + 1);
    let mut pass_names: Vec<String> = Vec::with_capacity(inputs.plan.len());
    let mut events: Vec<EventRecord> = Vec::with_capacity(inputs.plan.len());

    revisions.push(initial);

    for node in inputs.sorted_order {
        let pass = &inputs.plan.passes()[node.index()];
        let prev = revisions.last().expect("non-empty by construction");
        let pass_index = pass_names.len();
        let pass_name = pass.name().to_string();
        // `pass.effects()` returns the **static compile-time**
        // effect declarations (see [`crate::pass::PassCompileEffect`]
        // and its rename note). They flow into the schedule
        // analyzer and the `EventRecord.compile_effects` ledger.
        // They are *not* the source of truth for runtime derived-
        // state refresh — the per-pass `simulation_events` buffer
        // captured below is.
        let compile_effects = pass.effects();
        let trace_start = trace.len();
        let pre = StateSummary::from_simulation(prev);
        let mut trace_delta = Trace::new();
        // Per-pass capture buffer for `SimulationEvent`s emitted by
        // the pass's internal builders. Threaded into `PassContext`
        // so every event-emitting builder method (assign_allele,
        // update_trim, add_region, set_mutation_count, …) routes its
        // broadcast here for the duration of the pass. Drained into
        // the `EventRecord` after the pass commits + hooks run, so
        // the outcome ledger carries both the trace span (choices)
        // and the simulation-event stream (consequences) for every
        // committed pass.
        let mut simulation_events: Vec<SimulationEvent> = Vec::new();
        let mut ctx = PassContext {
            trace: &mut trace_delta,
            rng: &mut rng,
            pass_index,
            refdata: inputs.refdata,
            contracts: inputs.contracts,
            feasibility: inputs.feasibility,
            reference_index: inputs.reference_index,
            replay_cursor: cursor_opt.as_mut(),
            event_log_sink: Some(&mut simulation_events),
        };

        // Replay-mode forces the strict dispatch path even when the
        // policy is permissive. Replay errors are structural — they
        // mean "this trace doesn't match this plan" — and silently
        // panicking inside `execute().expect(...)` would erase the
        // structured diagnostic. With `replay_records: Some(_)` set,
        // every pass goes through `execute_checked` so any
        // `PassError::Replay` propagates cleanly.
        let effective_policy = if inputs.replay_records.is_some() {
            ExecutionPolicy::Strict
        } else {
            inputs.policy
        };
        let mut next = match effective_policy {
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
        // Hooks see both the static `compile_effects` (kept for
        // scheduling-aware hooks) and the runtime
        // `simulation_events` stream the pass actually emitted.
        // The events buffer is borrowed here and then moved into
        // the `EventRecord` below — no clone.
        for hook in inputs.hooks {
            next = hook.apply(next, &compile_effects, &simulation_events, hook_ctx);
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
            compile_effects,
            TraceSpan::new(trace_start, trace_end),
            pre,
            post,
            simulation_events,
        );

        trace.append_delta(trace_delta);
        pass_names.push(pass_name);
        revisions.push(next);
        events.push(event);
    }

    // Replay-mode tail check: with Tier 3 closed, every sampling
    // pass consumes from the cursor in fresh-RNG-record order, so
    // every record in the input trace must be consumed by the
    // plan. Trailing records mean the plan changed (or replay was
    // pointed at the wrong simulator); surface that as a
    // structured `PassError::Replay` so the failure is just as
    // visible as a mid-pass mismatch.
    if let Some(cursor) = cursor_opt.as_ref() {
        if let Err(reason) = cursor.assert_drained() {
            // Use the last executed pass name when available; otherwise
            // a stable sentinel that clearly isn't a real pass.
            let pass_name = pass_names
                .last()
                .cloned()
                .unwrap_or_else(|| "replay.end_of_plan".to_string());
            return Err(abort(
                PassError::replay(pass_name, reason),
                revisions,
                pass_names,
                trace,
                events,
            ));
        }
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
