//! [`LiveCallRefreshHook`] — keeps the V/D/J live-call sidecar in
//! sync with each pass's [`PassEffect`].
//!
//! Before Stage 2 of the dep-graph scheduler migration, the
//! "which pass effect refreshes which segment" cross-segment cascade
//! lived as a hand-coded `match` inside `compiled/execute.rs`. The
//! biology rules ("assembling D refreshes V because V's
//! right-extension walker can reach into D's bases", "an indel
//! invalidates V/D/J unconditionally because pool indices shift")
//! were prose comments inside the executor.
//!
//! This hook absorbs that logic. The runtime now registers exactly
//! one `LiveCallRefreshHook` per compile and dispatches it via the
//! generic [`crate::pass::EffectHook`] trait. The executor stays
//! biology-agnostic.

use crate::ir::{Segment, Simulation};
use crate::pass::{EffectHook, HookContext, PassEffect};

use super::{
    with_assembled_segment_live_call, DirtyWindow, ReferenceMatchIndex,
};

/// The standard derived-state refresh for V/D/J live calls. Reads
/// post-pass `PassEffect`s, walks the segment-relevant ones, and
/// rebuilds the affected `SegmentLiveCall`s from current pool state.
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

    fn apply(&self, mut sim: Simulation, effects: &[PassEffect], ctx: HookContext) -> Simulation {
        let Some(reference_index) = ctx.reference_index else {
            return sim;
        };

        for effect in effects {
            match effect {
                PassEffect::AssembleSegment(segment) => {
                    sim = with_assembled_segment_live_call(&sim, reference_index, *segment);
                    // Assembling a downstream segment introduces new
                    // bases that an earlier segment's right-extension
                    // walker can reach into when its allele suffix
                    // happens to match. Retrigger the upstream
                    // segment's refresh so any cross-boundary overlap
                    // surfaces as OVERLAPS_OTHER_SEGMENT on the
                    // upstream hypothesis.
                    //
                    // - Assembling D → refresh V.
                    // - Assembling J → refresh D AND V (covers
                    //   VJ-chain V→J overlap that would otherwise be
                    //   missed if only D refreshed on J assembly).
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
                            sim = with_assembled_segment_live_call(
                                &sim,
                                reference_index,
                                Segment::V,
                            );
                        }
                        _ => {}
                    }
                }
                // Any base edit (SHM, uniform mutation, PCR, quality
                // / N injection, contaminant overwrite) can change
                // which alleles the assembled bases support. Read
                // dirty windows stamped by the pass's
                // `DirtySignalObserver` and refresh only segments
                // whose region overlaps a dirty position. When no
                // observer was attached, `dirty_windows` is empty and
                // we fall back to the conservative full V/D/J sweep.
                PassEffect::EditBases => {
                    sim = refresh_segments_for_edit(sim, reference_index);
                }
                // An NP region appearing right-adjacent to V can
                // extend V's right boundary if the NP bases happen to
                // continue exactly into a V allele's suffix.
                PassEffect::AppendRegion(Segment::Np1) => {
                    sim = with_assembled_segment_live_call(&sim, reference_index, Segment::V);
                }
                // NP2 appears AFTER D is assembled, so D's right
                // boundary cannot pick up NP2 bases at assembly time.
                // Retrigger D refresh once NP2 exists. J's
                // left-extension does NOT need a separate hook here
                // because J is assembled after every NP region
                // exists.
                PassEffect::AppendRegion(Segment::Np2) => {
                    sim = with_assembled_segment_live_call(&sim, reference_index, Segment::D);
                }
                // Structural indels shift pool layout under V/D/J.
                // Refresh all three unconditionally and drain the
                // dirty-window log so the next pass's `EditBases`
                // dispatch starts clean. Dirty-window narrowing does
                // NOT apply here because an indel anywhere shifts
                // every region with a start ≥ the indel position.
                PassEffect::StructuralIndel => {
                    for segment in [Segment::V, Segment::D, Segment::J] {
                        sim = with_assembled_segment_live_call(&sim, reference_index, segment);
                    }
                    sim = drain_dirty_windows(sim);
                }
                _ => {}
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

    for segment in [Segment::V, Segment::D, Segment::J] {
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
    let Some(region) = sim
        .sequence
        .regions
        .iter()
        .find(|r| r.segment == segment)
    else {
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
