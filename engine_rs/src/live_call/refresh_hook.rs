//! [`LiveCallRefreshHook`] — keeps the V/D/J live-call sidecar in
//! sync with each pass's emitted [`SimulationEvent`] stream.
//!
//! Architectural stance:
//!
//! - [`PassCompileEffect`] = **static compile/schedule facts**.
//!   Used by the schedule analyzer to order passes and enforce
//!   dependencies. Still on the trait surface but no longer
//!   consulted at runtime by this hook.
//! - [`SimulationEvent`] = **runtime consequences**. The hook
//!   reads the per-pass event stream the compiled executor
//!   threads in, builds a [`LiveCallRefreshPlan`], and executes
//!   the resulting steps. The biology rules (D-refreshes-V,
//!   NP1-extends-V, indel-refreshes-everything) live entirely in
//!   [`super::refresh_plan`] now; this hook is the
//!   step-interpreter.
//!
//! This swap means a pass with no event-emitting state changes
//! produces no refresh — even if its compile-effect declarations
//! claim it did. Conversely, a pass that emits a `BaseChanged`
//! without declaring [`PassCompileEffect::EditBases`] will trigger
//! the refresh anyway. Runtime consequences are the source of
//! truth.

use crate::ir::{Segment, Simulation, SimulationEvent};
use crate::pass::{EffectHook, HookContext, PassCompileEffect};

use super::refresh_plan::{LiveCallRefreshPlan, LiveCallRefreshStep};
use super::{with_assembled_segment_live_call, DirtyWindow, ReferenceMatchIndex};

/// The standard derived-state refresh for V/D/J live calls. Reads
/// the post-pass [`SimulationEvent`] stream, derives a
/// [`LiveCallRefreshPlan`], and executes its steps against the
/// committed simulation.
pub struct LiveCallRefreshHook;

impl LiveCallRefreshHook {
    pub fn new() -> Self {
        Self
    }
}

impl Default for LiveCallRefreshHook {
    fn default() -> Self {
        Self::new()
    }
}

impl EffectHook for LiveCallRefreshHook {
    fn name(&self) -> &str {
        "live_call.refresh"
    }

    fn apply(
        &self,
        mut sim: Simulation,
        _compile_effects: &[PassCompileEffect],
        events: &[SimulationEvent],
        ctx: HookContext,
    ) -> Simulation {
        let Some(reference_index) = ctx.reference_index else {
            return sim;
        };

        // Translate the pass's runtime consequence stream into an
        // ordered refresh plan. `effects` is intentionally ignored:
        // it's the static declaration, not the source of truth.
        let plan = LiveCallRefreshPlan::from_events(events);

        for step in &plan.steps {
            match step {
                LiveCallRefreshStep::Segment(segment) => {
                    sim = with_assembled_segment_live_call(&sim, reference_index, *segment);
                }
                LiveCallRefreshStep::EditedSegments => {
                    sim = refresh_segments_for_edit(sim, reference_index);
                }
                LiveCallRefreshStep::AllStructural => {
                    for &segment in Segment::assignable() {
                        sim = with_assembled_segment_live_call(&sim, reference_index, segment);
                    }
                    sim = drain_dirty_windows(sim);
                }
            }
        }
        sim
    }
}

/// Refresh V/D/J live calls after a base-edit pass, using dirty
/// windows stamped by the pass's `DirtySignalObserver` to skip
/// segments whose region doesn't overlap any dirty position.
fn refresh_segments_for_edit(
    mut sim: Simulation,
    reference_index: &ReferenceMatchIndex,
) -> Simulation {
    let dirty_windows: Vec<DirtyWindow> = sim.dirty_log.windows.clone();
    let has_dirty = !dirty_windows.is_empty();

    for &segment in Segment::assignable() {
        if has_dirty && !region_overlaps_dirty(&sim, segment, &dirty_windows) {
            continue;
        }
        sim = with_assembled_segment_live_call(&sim, reference_index, segment);
    }

    if has_dirty {
        sim = drain_dirty_windows(sim);
    }
    sim
}

/// Does the segment's assembled region overlap any of the supplied
/// dirty windows? A region with no assembled instance returns
/// `false`.
fn region_overlaps_dirty(sim: &Simulation, segment: Segment, windows: &[DirtyWindow]) -> bool {
    let Some(region) = sim.sequence.regions.iter().find(|r| r.segment == segment) else {
        return false;
    };
    let region_start = region.start.index();
    let region_end = region.end.index();
    windows
        .iter()
        .any(|w| w.start < region_end && w.end > region_start)
}

/// Clear the dirty-window log so the next pass's `EditBases`
/// dispatch starts clean. No-op when the log is already empty.
fn drain_dirty_windows(sim: Simulation) -> Simulation {
    if sim.dirty_log.is_empty() {
        return sim;
    }
    sim.with_dirty_log(crate::live_call::DirtyLog::empty())
}
