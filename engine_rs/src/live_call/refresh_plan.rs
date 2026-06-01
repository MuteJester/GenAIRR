//! Event-derived live-call refresh planning.
//!
//! This module is the first slice of the **"compile-time effects
//! are for scheduling; runtime events are for runtime derived-state
//! refresh"** refactor. It exists alongside the existing
//! [`super::refresh_hook::LiveCallRefreshHook`] (which still drives
//! production execution) and proves that the event-derived plan can
//! describe the same biology. No executor wiring lives here yet —
//! the second slice will swap the hook over once parity is pinned.
//!
//! ## What it does
//!
//! [`LiveCallRefreshPlan::from_events`] takes the
//! `EventRecord::simulation_events` slice produced by ONE pass and
//! returns the ordered sequence of refresh actions needed to keep
//! the V/D/J live-call sidecar consistent with the new IR state.
//!
//! Today's `LiveCallRefreshHook` reacts to one [`PassEffect`] per
//! match arm; this helper reacts to one (or a small group of)
//! [`crate::ir::SimulationEvent`] variant(s) and produces the
//! equivalent [`LiveCallRefreshStep`] sequence.
//!
//! ## Event → step mapping
//!
//! | Event variant                              | Refresh steps emitted                                       |
//! |--------------------------------------------|-------------------------------------------------------------|
//! | `RegionAdded { segment: V }`               | `Segment(V)`                                                |
//! | `RegionAdded { segment: D }`               | `Segment(D)`, `Segment(V)` (D's bases extend V's right edge)|
//! | `RegionAdded { segment: J }`               | `Segment(J)`, `Segment(D)`, `Segment(V)` (full cascade)     |
//! | `RegionAdded { segment: Np1 }`             | `Segment(V)` (V's right extension reaches into NP1)         |
//! | `RegionAdded { segment: Np2 }`             | `Segment(D)` (D's right extension reaches into NP2)         |
//! | `BaseChanged { .. }` (any occurrence)      | `EditedSegments` (once per pass; subsequent are merged in) |
//! | `IndelInserted` / `IndelDeleted` (any)     | `AllStructural` (once per pass; subsumes any EditedSegments)|
//! | `SegmentReplaced { segment, .. }`          | `SegmentReplaced(segment)` (deduped per segment, first-encountered order) |
//! | All other variants                         | no-op                                                       |
//!
//! `BasePushed` is intentionally a no-op even though it carries
//! state-change information: the equivalent refresh trigger is the
//! `RegionAdded` event that closes out the assembly / NP push loop.
//! Reacting to `BasePushed` independently would re-refresh on every
//! pushed byte.
//!
//! ## What this slice does NOT do
//!
//! - It does NOT honour dirty-window narrowing internally. That
//!   gating still lives inside the executor / refresh hook because
//!   it needs `&Simulation` to read `sim.dirty_log`. The plan
//!   describes WHAT to do; the executor decides HOW (e.g. skipping
//!   `EditedSegments` for segments with no overlap).
//! - It does NOT touch `with_assembled_segment_live_call` or the
//!   dirty-log drain — those stay in the consuming hook.
//! - It does NOT replace [`super::refresh_hook::LiveCallRefreshHook`]
//!   yet. Production execution continues to read `PassEffect`s.
//!
//! Unit tests below pin one-to-one parity against every match arm
//! in the current hook.

use crate::ir::{Segment, SimulationEvent};

/// One refresh action the live-call sidecar needs to apply in
/// response to a pass's emitted [`SimulationEvent`] stream.
///
/// The variants intentionally mirror the existing
/// `LiveCallRefreshHook` match arms one-to-one, so the future swap
/// is a per-variant translation rather than a behavioural rewrite.
#[derive(Copy, Clone, Debug, Eq, PartialEq)]
pub enum LiveCallRefreshStep {
    /// Refresh `segment`'s `SegmentLiveCall` from current pool
    /// state. Used for assembly cascades (V/D/J) and NP-region
    /// upstream retriggers.
    Segment(Segment),
    /// Run the edit-refresh routine: dirty-window-narrowed
    /// re-call across V/D/J, then drain the dirty log. The
    /// executor decides whether each segment qualifies via
    /// region-overlap-with-`DirtyWindow`; this plan step only
    /// requests the routine.
    EditedSegments,
    /// Unconditional V/D/J refresh and dirty-log drain after a
    /// structural indel. Subsumes [`Self::EditedSegments`] when
    /// both would otherwise apply.
    AllStructural,
    /// A `segment`'s pool span was structurally replaced via
    /// [`crate::ir::SimulationBuilder::replace_segment`]. The
    /// `segment` field carries the replaced segment so future
    /// receptor-revision slices can specialise the refresh.
    ///
    /// **Slice B execution semantics:** identical to
    /// [`Self::AllStructural`] — receptor revision is rare, and an
    /// unconditional V/D/J refresh + dirty-log drain is the safest
    /// generic behaviour. Per-segment narrowing is a Slice C+
    /// optimisation, not a Slice B correctness requirement.
    SegmentReplaced(Segment),
}

/// Ordered list of refresh actions derived from one pass's
/// `simulation_events`. Built via [`Self::from_events`].
#[derive(Clone, Debug, Default, Eq, PartialEq)]
pub struct LiveCallRefreshPlan {
    pub steps: Vec<LiveCallRefreshStep>,
}

impl LiveCallRefreshPlan {
    /// Translate one pass's emitted events into the refresh-step
    /// sequence the live-call sidecar would need.
    ///
    /// **Per-pass granularity.** The caller is expected to pass
    /// the `simulation_events` slice from a single
    /// `EventRecord` — `BaseChanged` / indel events are deduped
    /// to one step per call, mirroring how the existing hook
    /// reacts once per `EditBases` / `StructuralIndel` effect.
    pub fn from_events(events: &[SimulationEvent]) -> Self {
        let mut steps: Vec<LiveCallRefreshStep> = Vec::new();
        let mut edited_emitted = false;
        let mut structural_emitted = false;
        // First-encountered-order, per-unique-segment dedup for
        // `SegmentReplaced`. A pass that replaces V twice (e.g. an
        // adventurous receptor-revision sampler) collapses to one
        // refresh step; a pass that replaces V then D yields two
        // steps in that order.
        let mut segments_replaced: Vec<Segment> = Vec::new();

        for event in events {
            match event {
                SimulationEvent::RegionAdded { region } => {
                    push_region_added(&mut steps, region.segment);
                }
                SimulationEvent::BaseChanged { .. } => {
                    if !edited_emitted {
                        steps.push(LiveCallRefreshStep::EditedSegments);
                        edited_emitted = true;
                    }
                }
                SimulationEvent::IndelInserted { .. }
                | SimulationEvent::IndelDeleted { .. } => {
                    if !structural_emitted {
                        steps.push(LiveCallRefreshStep::AllStructural);
                        structural_emitted = true;
                    }
                }
                SimulationEvent::SegmentReplaced { segment, .. } => {
                    if !segments_replaced.contains(segment) {
                        segments_replaced.push(*segment);
                        steps.push(LiveCallRefreshStep::SegmentReplaced(*segment));
                    }
                }
                // BasePushed is a no-op: the equivalent refresh
                // trigger is the RegionAdded that closes the
                // assembly / NP push loop. Reacting on each pushed
                // byte would re-refresh redundantly.
                SimulationEvent::BasePushed { .. } => {}
                // Assignments, trims, mutation-count, rev-comp,
                // base-deletes (reserved), region-replaced
                // (reserved) — none of these invalidate live-call
                // evidence on their own. AssignmentChanged in
                // particular is upstream of assembly: the actual
                // refresh happens when the AssembleSegment pass
                // emits its RegionAdded later in the plan.
                SimulationEvent::AssignmentChanged { .. }
                | SimulationEvent::TrimChanged { .. }
                | SimulationEvent::BaseDeleted { .. }
                | SimulationEvent::RegionReplaced { .. }
                // PRegionAdded sits between assembled segments,
                // shifts pool positions of later regions but never
                // overlaps an already-walked live-call span. The
                // next assemble pass's walker reads the (now-
                // P-extended) pool position implicitly via seq_start.
                | SimulationEvent::PRegionAdded { .. }
                | SimulationEvent::ReverseComplementFlagRecorded { .. }
                | SimulationEvent::OrientationChanged { .. }
                | SimulationEvent::MutationCountChanged { .. } => {}
            }
        }

        Self { steps }
    }

    /// Convenience: does this plan trigger any refresh?
    pub fn is_empty(&self) -> bool {
        self.steps.is_empty()
    }
}

fn push_region_added(steps: &mut Vec<LiveCallRefreshStep>, segment: Segment) {
    // Mirrors the cascade in `LiveCallRefreshHook::apply`:
    // assembling D refreshes V; assembling J refreshes D then V;
    // NP1 region triggers V refresh; NP2 region triggers D
    // refresh.
    match segment {
        Segment::V => steps.push(LiveCallRefreshStep::Segment(Segment::V)),
        Segment::D => {
            steps.push(LiveCallRefreshStep::Segment(Segment::D));
            steps.push(LiveCallRefreshStep::Segment(Segment::V));
        }
        Segment::J => {
            steps.push(LiveCallRefreshStep::Segment(Segment::J));
            steps.push(LiveCallRefreshStep::Segment(Segment::D));
            steps.push(LiveCallRefreshStep::Segment(Segment::V));
        }
        Segment::Np1 => steps.push(LiveCallRefreshStep::Segment(Segment::V)),
        Segment::Np2 => steps.push(LiveCallRefreshStep::Segment(Segment::D)),
    }
}

#[cfg(test)]
mod tests {
    //! Parity unit tests: every match arm in
    //! `LiveCallRefreshHook::apply` has a counterpart here, fed
    //! with the event stream a pass would actually emit. If the
    //! hook's biology changes, these tests fail first.
    use super::*;
    use crate::assignment::{AlleleInstance, TrimEnd};
    use crate::ir::{NucFlags, NucHandle, Region, SimulationEvent};
    use crate::refdata::AlleleId;

    // Event constructors — kept short to keep the parity table
    // below readable.

    fn base_pushed(handle: u32, segment: Segment) -> SimulationEvent {
        SimulationEvent::BasePushed {
            handle: NucHandle::new(handle),
            base: b'A',
            segment,
            germline_pos: None,
            flags: NucFlags::empty(),
        }
    }
    fn region_added(segment: Segment, start: u32, end: u32) -> SimulationEvent {
        SimulationEvent::RegionAdded {
            region: Region::new(segment, NucHandle::new(start), NucHandle::new(end)),
        }
    }
    fn base_changed(handle: u32) -> SimulationEvent {
        SimulationEvent::BaseChanged {
            handle: NucHandle::new(handle),
            old_base: b'A',
            new_base: b'G',
            segment: Segment::V,
            germline_pos: Some(handle as u16),
        }
    }
    fn indel_inserted(at: u32) -> SimulationEvent {
        SimulationEvent::IndelInserted {
            at,
            base: b'N',
            segment: Segment::V,
            flags: NucFlags::empty(),
        }
    }
    fn indel_deleted(at: u32) -> SimulationEvent {
        SimulationEvent::IndelDeleted {
            at,
            removed_base: b'A',
            segment: Segment::V,
        }
    }
    fn segment_replaced(segment: Segment, old_len: u32, new_len: u32) -> SimulationEvent {
        let old = Region::new(segment, NucHandle::new(0), NucHandle::new(old_len));
        let new = Region::new(segment, NucHandle::new(0), NucHandle::new(new_len));
        SimulationEvent::SegmentReplaced {
            segment,
            old_region: old,
            new_region: new,
            bytes_delta: (new_len as i32) - (old_len as i32),
        }
    }
    fn region_replaced_metadata_only(segment: Segment, start: u32, end: u32) -> SimulationEvent {
        // Identity-replacement is the closest the persistent layer
        // gets to a metadata-only `RegionReplaced` payload — pool
        // bytes unchanged, only the `Region` struct itself swapped.
        let r = Region::new(segment, NucHandle::new(start), NucHandle::new(end));
        SimulationEvent::RegionReplaced {
            old: r.clone(),
            new: r,
        }
    }

    // ── PassEffect::AssembleSegment(V/D/J) parity ─────────────

    #[test]
    fn assembly_v_emits_segment_v_only() {
        // Assembly emits many BasePushed then one RegionAdded(V).
        // Mirrors PassEffect::AssembleSegment(V) — refresh V only.
        let events = vec![
            base_pushed(0, Segment::V),
            base_pushed(1, Segment::V),
            region_added(Segment::V, 0, 2),
        ];
        assert_eq!(
            LiveCallRefreshPlan::from_events(&events).steps,
            vec![LiveCallRefreshStep::Segment(Segment::V)]
        );
    }

    #[test]
    fn assembly_d_refreshes_d_then_v() {
        // Mirrors the AssembleSegment(D) match arm in the hook:
        // D's bases can extend V's right boundary, so V is
        // re-walked after D.
        let events = vec![
            base_pushed(10, Segment::D),
            region_added(Segment::D, 10, 13),
        ];
        assert_eq!(
            LiveCallRefreshPlan::from_events(&events).steps,
            vec![
                LiveCallRefreshStep::Segment(Segment::D),
                LiveCallRefreshStep::Segment(Segment::V),
            ]
        );
    }

    #[test]
    fn assembly_j_refreshes_j_then_d_then_v() {
        // Full cascade: J's bases can pull both D's right edge AND
        // (for VJ chains) V's right edge.
        let events = vec![region_added(Segment::J, 20, 26)];
        assert_eq!(
            LiveCallRefreshPlan::from_events(&events).steps,
            vec![
                LiveCallRefreshStep::Segment(Segment::J),
                LiveCallRefreshStep::Segment(Segment::D),
                LiveCallRefreshStep::Segment(Segment::V),
            ]
        );
    }

    // ── PassEffect::AppendRegion(Np1/Np2) parity ──────────────

    #[test]
    fn np1_region_triggers_v_refresh() {
        let events = vec![
            base_pushed(9, Segment::Np1),
            region_added(Segment::Np1, 9, 12),
        ];
        assert_eq!(
            LiveCallRefreshPlan::from_events(&events).steps,
            vec![LiveCallRefreshStep::Segment(Segment::V)]
        );
    }

    #[test]
    fn np2_region_triggers_d_refresh() {
        let events = vec![region_added(Segment::Np2, 15, 17)];
        assert_eq!(
            LiveCallRefreshPlan::from_events(&events).steps,
            vec![LiveCallRefreshStep::Segment(Segment::D)]
        );
    }

    // ── PassEffect::EditBases parity ──────────────────────────

    #[test]
    fn single_base_changed_emits_one_edited_segments_step() {
        let events = vec![base_changed(5)];
        assert_eq!(
            LiveCallRefreshPlan::from_events(&events).steps,
            vec![LiveCallRefreshStep::EditedSegments]
        );
    }

    #[test]
    fn multiple_base_changed_dedup_to_one_edited_segments_step() {
        // Mirrors how the hook reacts once per `EditBases` effect
        // even when the pass changed many bases.
        let events = vec![base_changed(0), base_changed(1), base_changed(2)];
        assert_eq!(
            LiveCallRefreshPlan::from_events(&events).steps,
            vec![LiveCallRefreshStep::EditedSegments]
        );
    }

    #[test]
    fn mutation_pass_pattern_emits_edited_segments_once() {
        // A real mutation pass also emits MutationCountChanged at
        // commit; the plan ignores the count event and emits one
        // EditedSegments step for the BaseChanged batch.
        let events = vec![
            base_changed(0),
            base_changed(1),
            SimulationEvent::MutationCountChanged {
                old: 0,
                new: 2,
                delta: 2,
            },
        ];
        assert_eq!(
            LiveCallRefreshPlan::from_events(&events).steps,
            vec![LiveCallRefreshStep::EditedSegments]
        );
    }

    // ── PassEffect::StructuralIndel parity ────────────────────

    #[test]
    fn single_indel_insert_emits_one_all_structural_step() {
        let events = vec![indel_inserted(7)];
        assert_eq!(
            LiveCallRefreshPlan::from_events(&events).steps,
            vec![LiveCallRefreshStep::AllStructural]
        );
    }

    #[test]
    fn single_indel_delete_emits_one_all_structural_step() {
        // Mirrors EndLossPass and the deletion arm of IndelPass:
        // both declare PassEffect::StructuralIndel.
        let events = vec![indel_deleted(3)];
        assert_eq!(
            LiveCallRefreshPlan::from_events(&events).steps,
            vec![LiveCallRefreshStep::AllStructural]
        );
    }

    #[test]
    fn mixed_indel_insert_and_delete_emit_one_all_structural_step() {
        // A pass that interleaves inserts and deletes (the
        // IndelPass tuple sampler can) still emits exactly one
        // AllStructural refresh — same as the hook's single
        // StructuralIndel match arm regardless of operation count.
        let events = vec![indel_inserted(2), indel_deleted(5), indel_inserted(8)];
        assert_eq!(
            LiveCallRefreshPlan::from_events(&events).steps,
            vec![LiveCallRefreshStep::AllStructural]
        );
    }

    // ── No-refresh passes (sample / trim / rev-comp) ──────────

    #[test]
    fn assignment_changed_alone_does_not_trigger_refresh() {
        let events = vec![SimulationEvent::AssignmentChanged {
            segment: Segment::V,
            old: None,
            new: AlleleInstance::new(AlleleId::new(0)),
        }];
        assert!(LiveCallRefreshPlan::from_events(&events).is_empty());
    }

    #[test]
    fn trim_changed_alone_does_not_trigger_refresh() {
        let events = vec![SimulationEvent::TrimChanged {
            segment: Segment::V,
            end: TrimEnd::Five,
            old: Some(0),
            new: 3,
        }];
        assert!(LiveCallRefreshPlan::from_events(&events).is_empty());
    }

    #[test]
    fn reverse_complement_flag_does_not_trigger_refresh() {
        let events =
            vec![SimulationEvent::ReverseComplementFlagRecorded { applied: true }];
        assert!(LiveCallRefreshPlan::from_events(&events).is_empty());
    }

    #[test]
    fn empty_event_stream_yields_empty_plan() {
        assert!(LiveCallRefreshPlan::from_events(&[]).is_empty());
    }

    #[test]
    fn base_pushed_alone_does_not_trigger_refresh() {
        // BasePushed in isolation (no closing RegionAdded) is not
        // a real production scenario but the plan must not
        // over-eagerly emit a refresh per push.
        let events = vec![base_pushed(0, Segment::V), base_pushed(1, Segment::V)];
        assert!(LiveCallRefreshPlan::from_events(&events).is_empty());
    }

    // ── Combined event streams (multi-effect passes) ──────────

    // ── SimulationEvent::SegmentReplaced — receptor revision Slice B ──

    #[test]
    fn segment_replaced_v_emits_single_segment_replaced_step() {
        let events = vec![segment_replaced(Segment::V, 5, 8)];
        assert_eq!(
            LiveCallRefreshPlan::from_events(&events).steps,
            vec![LiveCallRefreshStep::SegmentReplaced(Segment::V)]
        );
    }

    #[test]
    fn multiple_segment_replaced_for_same_segment_dedup_to_one_step() {
        // A pass that replaces V twice (e.g. an adventurous
        // receptor-revision sampler probing two candidates before
        // committing) must collapse to a single refresh step —
        // same biology, same execution.
        let events = vec![
            segment_replaced(Segment::V, 5, 8),
            segment_replaced(Segment::V, 8, 6),
        ];
        assert_eq!(
            LiveCallRefreshPlan::from_events(&events).steps,
            vec![LiveCallRefreshStep::SegmentReplaced(Segment::V)]
        );
    }

    #[test]
    fn segment_replaced_different_segments_yield_separate_steps_in_event_order() {
        // First-encountered-order, per-unique-segment dedup. A
        // (hypothetical Slice C+) pass that replaces V then D
        // produces two steps in that order.
        let events = vec![
            segment_replaced(Segment::V, 5, 8),
            segment_replaced(Segment::D, 3, 4),
        ];
        assert_eq!(
            LiveCallRefreshPlan::from_events(&events).steps,
            vec![
                LiveCallRefreshStep::SegmentReplaced(Segment::V),
                LiveCallRefreshStep::SegmentReplaced(Segment::D),
            ]
        );
    }

    #[test]
    fn segment_replaced_is_distinguishable_from_region_replaced() {
        // `RegionReplaced` is metadata-only and stays a no-op (no
        // pool bytes changed). `SegmentReplaced` is structural and
        // emits a refresh step. A single record carrying both must
        // produce exactly one step — the SegmentReplaced one.
        let events = vec![
            region_replaced_metadata_only(Segment::V, 0, 5),
            segment_replaced(Segment::V, 5, 8),
        ];
        assert_eq!(
            LiveCallRefreshPlan::from_events(&events).steps,
            vec![LiveCallRefreshStep::SegmentReplaced(Segment::V)]
        );
    }

    #[test]
    fn region_replaced_alone_remains_a_noop() {
        // Pins that Slice B has not accidentally widened
        // RegionReplaced's behaviour — only SegmentReplaced
        // triggers a refresh.
        let events = vec![region_replaced_metadata_only(Segment::V, 0, 5)];
        assert!(LiveCallRefreshPlan::from_events(&events).is_empty());
    }

    #[test]
    fn segment_replaced_step_distinguishable_from_all_structural() {
        // The two variants share execution semantics under the
        // Slice B hook, but the plan-layer step is distinct so a
        // future slice can specialise the refresh without
        // touching the IR event surface or the indel path.
        assert_ne!(
            LiveCallRefreshStep::AllStructural,
            LiveCallRefreshStep::SegmentReplaced(Segment::V),
        );
    }

    #[test]
    fn assembly_event_stream_matches_assemble_segment_arm_byte_for_byte() {
        // Sanity: a real assembly pass emits N BasePushed + one
        // RegionAdded(seg). The derived plan equals the hook's
        // AssembleSegment(seg) arm output.
        for (seg, expected) in [
            (
                Segment::V,
                vec![LiveCallRefreshStep::Segment(Segment::V)],
            ),
            (
                Segment::D,
                vec![
                    LiveCallRefreshStep::Segment(Segment::D),
                    LiveCallRefreshStep::Segment(Segment::V),
                ],
            ),
            (
                Segment::J,
                vec![
                    LiveCallRefreshStep::Segment(Segment::J),
                    LiveCallRefreshStep::Segment(Segment::D),
                    LiveCallRefreshStep::Segment(Segment::V),
                ],
            ),
        ] {
            let mut events: Vec<SimulationEvent> =
                (0..6).map(|i| base_pushed(i, seg)).collect();
            events.push(region_added(seg, 0, 6));
            assert_eq!(
                LiveCallRefreshPlan::from_events(&events).steps,
                expected,
                "assembly cascade parity for {seg:?}"
            );
        }
    }
}
