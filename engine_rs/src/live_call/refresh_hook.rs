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
        // Guard against redundant full V/D/J re-walks: once an
        // `AllStructural` or `SegmentReplaced(_)` step has executed,
        // subsequent steps of either kind would just recompute the
        // same calls against the same final pool. Track once and
        // short-circuit.
        let mut all_segments_refreshed = false;

        for step in &plan.steps {
            match step {
                LiveCallRefreshStep::Segment(segment) => {
                    sim = with_assembled_segment_live_call(&sim, reference_index, *segment);
                }
                LiveCallRefreshStep::EditedSegments => {
                    sim = refresh_segments_for_edit(sim, reference_index);
                }
                LiveCallRefreshStep::AllStructural
                | LiveCallRefreshStep::SegmentReplaced(_) => {
                    if all_segments_refreshed {
                        continue;
                    }
                    for &segment in Segment::assignable() {
                        sim = with_assembled_segment_live_call(&sim, reference_index, segment);
                    }
                    sim = drain_dirty_windows(sim);
                    all_segments_refreshed = true;
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

#[cfg(test)]
mod tests {
    //! Receptor-revision Slice B integration tests for the
    //! `LiveCallRefreshHook`. Drives a `SimulationBuilder` through
    //! `replace_segment`, captures the emitted event stream, and
    //! verifies the hook produces a `segment_calls` sidecar
    //! identical to a from-scratch `assembled_segment_live_call`
    //! against the post-replace pool.
    //!
    //! The "AllStructural-equivalent" semantics the hook applies
    //! to [`LiveCallRefreshStep::SegmentReplaced`] are correctness
    //! tests against the cache-parity oracle, not behavioural
    //! introspection of the hook's internal control flow.

    use super::*;
    use crate::ir::{NucHandle, Nucleotide, Region, SimulationBuilder};
    use crate::live_call::{
        assembled_segment_live_call, with_assembled_segment_live_call, DirtyLog, DirtyReason,
        DirtyWindow, ReferenceMatchIndex,
    };
    use crate::refdata::{Allele, ChainType, RefDataConfig};

    fn allele(name: &str, seq: &[u8]) -> Allele {
        Allele {
            name: name.to_string(),
            gene: name.split('*').next().unwrap_or(name).to_string(),
            seq: seq.to_vec(),
            segment: Segment::V,
            anchor: None,
            functional_status: None,
            subregions: Vec::new(),
        }
    }

    fn sim_with_v_region(bases: &[u8]) -> Simulation {
        let mut sim = Simulation::new();
        for (i, &b) in bases.iter().enumerate() {
            let (next, _) = sim.with_nucleotide_pushed(Nucleotide::germline(
                b,
                i as u16,
                Segment::V,
            ));
            sim = next;
        }
        sim.with_region_added(Region::new(
            Segment::V,
            NucHandle::new(0),
            NucHandle::new(bases.len() as u32),
        ))
    }

    fn v_replacement(bases: &[u8]) -> Vec<Nucleotide> {
        bases
            .iter()
            .enumerate()
            .map(|(i, &b)| Nucleotide::germline(b, i as u16, Segment::V))
            .collect()
    }

    /// Compare the hook's installed `segment_calls` against a
    /// from-scratch oracle for every assignable segment.
    fn assert_segment_calls_match_oracle(
        sim: &Simulation,
        reference_index: &ReferenceMatchIndex,
    ) {
        let evidence_version = sim.segment_calls.version.saturating_add(1);
        for &segment in Segment::assignable() {
            let cached = sim.segment_calls.get(segment).cloned();
            let oracle = assembled_segment_live_call(
                sim,
                reference_index,
                segment,
                evidence_version,
            );
            match (&cached, &oracle) {
                (None, None) => {}
                (Some(c), Some(o)) => {
                    assert_eq!(
                        c.allele_call.to_ids(),
                        o.allele_call.to_ids(),
                        "segment {segment:?} allele tie-set diverges from oracle",
                    );
                }
                (cached, oracle) => panic!(
                    "segment {segment:?} call presence diverges: cached={:?} oracle={:?}",
                    cached.is_some(),
                    oracle.is_some(),
                ),
            }
        }
    }

    #[test]
    fn hook_refreshes_v_call_after_segment_replaced_to_distinct_allele() {
        // Build refdata with two V alleles. Start with a sim that
        // matches allele A1 byte-for-byte; replace_segment to A2's
        // bytes; the hook must update the V call to point at A2.
        let mut cfg = RefDataConfig::empty(ChainType::Vdj);
        let _a1 = cfg.v_pool.push(allele("V1*01", b"AACCGG"));
        let a2 = cfg.v_pool.push(allele("V1*02", b"TTGGCA"));
        let reference_index = ReferenceMatchIndex::build(&cfg);

        // Stage 1: assemble V matching A1; install the call.
        let mut sim = sim_with_v_region(b"AACCGG");
        sim = with_assembled_segment_live_call(&sim, &reference_index, Segment::V);
        let pre_v = sim.segment_calls.get(Segment::V).cloned().unwrap();
        assert!(
            pre_v.allele_call.iter_ids().any(|id| id != a2),
            "before replacement V should not be tagged as A2; got {:?}",
            pre_v.allele_call.to_ids(),
        );

        // Stage 2: receptor-revise V to A2's bytes via the builder.
        let mut builder = SimulationBuilder::from_simulation(sim);
        builder.attach_event_log_observer();
        let _ = builder.replace_segment(Segment::V, v_replacement(b"TTGGCA"));
        let events = builder.seal_event_log_observer();
        let post_sim = builder.seal();

        assert_eq!(events.len(), 1);
        assert!(matches!(
            events[0],
            SimulationEvent::SegmentReplaced {
                segment: Segment::V,
                ..
            }
        ));

        // Stage 3: run the hook and check parity vs. the oracle.
        let ctx = HookContext {
            reference_index: Some(&reference_index),
        };
        let hook = LiveCallRefreshHook::new();
        let refreshed = hook.apply(post_sim, &[], &events, ctx);

        assert_segment_calls_match_oracle(&refreshed, &reference_index);
        let post_v = refreshed.segment_calls.get(Segment::V).cloned().unwrap();
        assert_eq!(
            post_v.allele_call.to_ids(),
            vec![a2],
            "V call must now point at the replacement allele A2",
        );
    }

    #[test]
    fn hook_refreshes_after_length_changing_segment_replacement() {
        // Length-changing replacement (V grows from 5 → 7 bytes).
        // The new bytes correspond to allele A_long; the from-scratch
        // oracle must agree with the cached call.
        let mut cfg = RefDataConfig::empty(ChainType::Vdj);
        let _a_short = cfg.v_pool.push(allele("Vshort*01", b"AAGGC"));
        let a_long = cfg.v_pool.push(allele("Vlong*01", b"TTCCAGG"));
        let reference_index = ReferenceMatchIndex::build(&cfg);

        let sim = sim_with_v_region(b"AAGGC");
        let sim = with_assembled_segment_live_call(&sim, &reference_index, Segment::V);

        let mut builder = SimulationBuilder::from_simulation(sim);
        builder.attach_event_log_observer();
        let _ = builder.replace_segment(Segment::V, v_replacement(b"TTCCAGG"));
        let events = builder.seal_event_log_observer();
        let post_sim = builder.seal();

        // Pool grew by 2 bytes; event payload encodes that.
        assert!(post_sim.pool.len() == 7);
        match &events[0] {
            SimulationEvent::SegmentReplaced { bytes_delta, .. } => {
                assert_eq!(*bytes_delta, 2);
            }
            other => panic!("expected SegmentReplaced, got {other:?}"),
        }

        let ctx = HookContext {
            reference_index: Some(&reference_index),
        };
        let refreshed = LiveCallRefreshHook::new().apply(post_sim, &[], &events, ctx);

        assert_segment_calls_match_oracle(&refreshed, &reference_index);
        let v = refreshed.segment_calls.get(Segment::V).cloned().unwrap();
        assert_eq!(v.allele_call.to_ids(), vec![a_long]);
    }

    #[test]
    fn hook_drains_dirty_log_on_segment_replaced() {
        // The AllStructural-equivalent semantics include a dirty-log
        // drain. Plant a fake dirty window, fire SegmentReplaced,
        // and confirm the log is empty after the hook runs.
        let mut cfg = RefDataConfig::empty(ChainType::Vdj);
        let _ = cfg.v_pool.push(allele("V1*01", b"AACC"));
        let _ = cfg.v_pool.push(allele("V1*02", b"GGTT"));
        let reference_index = ReferenceMatchIndex::build(&cfg);

        let sim = sim_with_v_region(b"AACC");
        let sim = sim.with_dirty_log(DirtyLog {
            windows: vec![DirtyWindow::new(
                0,
                1,
                DirtyReason::BaseEdited { site: 0 },
            )],
        });
        assert!(!sim.dirty_log.is_empty());

        let mut builder = SimulationBuilder::from_simulation(sim);
        builder.attach_event_log_observer();
        let _ = builder.replace_segment(Segment::V, v_replacement(b"GGTT"));
        let events = builder.seal_event_log_observer();
        let post_sim = builder.seal();

        let ctx = HookContext {
            reference_index: Some(&reference_index),
        };
        let refreshed = LiveCallRefreshHook::new().apply(post_sim, &[], &events, ctx);
        assert!(
            refreshed.dirty_log.is_empty(),
            "SegmentReplaced step must drain the dirty log just like AllStructural",
        );
    }

    #[test]
    fn hook_idempotent_when_multiple_segment_replaced_events_fire() {
        // Two SegmentReplaced events in one pass (hypothetical
        // future receptor-revision sampler) collapse to one
        // refresh — the second is a no-op. The post-hook state
        // must equal the result of running the hook on a stream
        // with a single replacement event.
        let mut cfg = RefDataConfig::empty(ChainType::Vdj);
        let _ = cfg.v_pool.push(allele("V1*01", b"AACCGG"));
        let _ = cfg.v_pool.push(allele("V1*02", b"TTGGCA"));
        let reference_index = ReferenceMatchIndex::build(&cfg);

        let sim = sim_with_v_region(b"AACCGG");
        let sim = with_assembled_segment_live_call(&sim, &reference_index, Segment::V);

        let mut builder = SimulationBuilder::from_simulation(sim);
        builder.attach_event_log_observer();
        let _ = builder.replace_segment(Segment::V, v_replacement(b"AAAAAA"));
        let _ = builder.replace_segment(Segment::V, v_replacement(b"TTGGCA"));
        let events = builder.seal_event_log_observer();
        let post_sim = builder.seal();

        assert_eq!(events.len(), 2);
        let ctx = HookContext {
            reference_index: Some(&reference_index),
        };
        let refreshed = LiveCallRefreshHook::new().apply(post_sim, &[], &events, ctx);

        // The final pool reflects the *second* replacement; the
        // refreshed call must match it.
        assert_segment_calls_match_oracle(&refreshed, &reference_index);
    }
}
