//! Transactional execution of a compiled `PassPlan`.
//!
//! `execute_transactional` is the engine's run loop. It owns the
//! trace, RNG, revision stack, and event log; each pass either commits
//! its delta or aborts the run. `apply_live_call_updates` is the
//! per-pass live-caller refresh hook that keeps the V/D/J live calls
//! in sync with assembly/edit/NP/indel effects. `abort` packages a
//! `PassError` together with the already-committed prefix so the
//! caller can introspect partial state.

use super::ExecutionPolicy;
use crate::contract::ContractSet;
use crate::event::{EventRecord, StateSummary, TraceSpan};
use crate::feasibility::FeasibilityContext;
use crate::ir::{Segment, Simulation};
use crate::live_call::{with_assembled_segment_live_call, ReferenceMatchIndex};
use crate::pass::{Outcome, PassContext, PassEffect, PassError, PassPlan};
use crate::refdata::RefDataConfig;
use crate::rng::Rng;
use crate::trace::Trace;

#[derive(Copy, Clone)]
pub(super) struct ExecutionInputs<'a> {
    pub(super) plan: &'a PassPlan,
    pub(super) refdata: Option<&'a RefDataConfig>,
    pub(super) contracts: Option<&'a ContractSet>,
    pub(super) feasibility: Option<&'a FeasibilityContext>,
    pub(super) reference_index: Option<&'a ReferenceMatchIndex>,
    pub(super) policy: ExecutionPolicy,
}

#[derive(Debug)]
pub(super) struct ExecutionAbort {
    pub(super) error: PassError,
    /// Already-committed prefix at the point of failure. Kept private:
    /// public runners expose only the structured pass error, while
    /// internal tests can assert transaction commit boundaries.
    #[cfg_attr(not(test), allow(dead_code))]
    pub(super) committed: Outcome,
}

pub(super) fn execute_transactional(
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

    for pass in inputs.plan.passes() {
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

        if let Some(reference_index) = inputs.reference_index {
            next = apply_live_call_updates(next, &effects, reference_index);
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

fn apply_live_call_updates(
    mut sim: Simulation,
    effects: &[PassEffect],
    reference_index: &ReferenceMatchIndex,
) -> Simulation {
    for effect in effects {
        match effect {
            PassEffect::AssembleSegment(segment) => {
                sim = with_assembled_segment_live_call(&sim, reference_index, *segment);
                // assembling a downstream segment introduces
                // new bases that an earlier segment's right-extension
                // walker can reach into when its allele suffix happens
                // to match. Retrigger the upstream segment's refresh
                // so any cross-boundary overlap surfaces as
                // OVERLAPS_OTHER_SEGMENT on the upstream hypothesis.
                //
                // - Assembling D → refresh V (V right overlaps into D).
                // - Assembling J → refresh D (D right overlaps into J).
                //
                // No symmetric upstream hook for V right-into-J (we
                // never assemble J without D in VDJ chains, and the
                // walker would have to traverse D's region first;
                // the targeted hops stay conservative).
                match segment {
                    Segment::D => {
                        sim = with_assembled_segment_live_call(
                            &sim,
                            reference_index,
                            Segment::V,
                        );
                    }
                    Segment::J => {
                        sim = with_assembled_segment_live_call(
                            &sim,
                            reference_index,
                            Segment::D,
                        );
                    }
                    _ => {}
                }
            }
            // any base edit (SHM, uniform mutation, PCR, quality
            // / N injection, contaminant overwrite) can change which
            // alleles the assembled bases support. We do a conservative
            // refresh — every assembled V/D/J segment is recomputed from
            // the current pool. Segments without an assembled region are
            // a no-op inside `with_assembled_segment_live_call`.
            //
            // A dirty-window optimization (scanning the trace delta for
            // edited positions and skipping segments outside the dirty
            // range) was evaluated. The trace-scan cost (~5 µs/record)
            // cancelled the refresh savings — net change was 0-2
            // µs/record on typical loads, sometimes a small loss. The
            // simpler full-refresh policy is kept for clarity.
            PassEffect::EditBases => {
                for segment in [Segment::V, Segment::D, Segment::J] {
                    sim = with_assembled_segment_live_call(&sim, reference_index, segment);
                }
            }
            // an NP region appearing right-adjacent to V
            // can extend V's right boundary if the NP bases happen to
            // continue exactly into a V allele's suffix. The walker
            // inside `with_assembled_segment_live_call` does the
            // extension automatically — we just need to retrigger the
            // V refresh now that NP1 exists.
            PassEffect::AppendRegion(Segment::Np1) => {
                sim = with_assembled_segment_live_call(&sim, reference_index, Segment::V);
            }
            // NP2 appears AFTER D is assembled (the typical
            // VDJ pipeline order is `... → assemble(D) → np2 → ...`),
            // so D's right boundary cannot pick up NP2 bases at
            // assembly time. Retrigger D refresh once NP2 exists so
            // the walker can extend D's right side into NP2 if the
            // bases match.
            //
            // J left-extension does NOT need a separate hook here:
            // J is assembled after every NP region exists, so its
            // own `AssembleSegment(J)` refresh already sees the full
            // NP context.
            PassEffect::AppendRegion(Segment::Np2) => {
                sim = with_assembled_segment_live_call(&sim, reference_index, Segment::D);
            }
            // structural indels (insertions / deletions) shift
            // the pool layout under V/D/J coding regions. Indel-inserted
            // nucleotides have GermlinePos::NONE so the walker skips them;
            // deletions show up as forward jumps in `germline_pos` which
            // the walker now tolerates. Refresh every assembled V/D/J
            // segment from the post-indel pool so the live calls reflect
            // the new evidence layout.
            PassEffect::StructuralIndel => {
                for segment in [Segment::V, Segment::D, Segment::J] {
                    sim = with_assembled_segment_live_call(&sim, reference_index, segment);
                }
            }
            _ => {}
        }
    }
    sim
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
