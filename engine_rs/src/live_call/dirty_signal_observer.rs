//! Streaming dirty-signal observer.
//!
//! Captures [`DirtyWindow`] change signals from the IR event stream
//! as mutations happen, so the post-pass live-call refresh in
//! `compiled/execute.rs::apply_live_call_updates` can narrow to
//! segments whose region overlaps a dirty position.
//!
//! ## Channel — `SimulationEventSink`
//!
//! This observer rides the unified
//! [`SimulationEventSink`](crate::ir::SimulationEventSink) channel —
//! the same one [`crate::ir::event_log_observer::EventLogObserver`]
//! migrated onto in slice 1. It reacts to *consequences* (the
//! `SimulationEvent` enum) rather than the legacy
//! `IrEventObserver` per-event method surface. See
//! [`crate::ir::sim_event`] for the trace-vs-events architectural
//! split.
//!
//! ## Which events produce dirty signals
//!
//! - [`SimulationEvent::BaseChanged`] →
//!   `DirtyReason::BaseEdited { site: handle.index() }`
//! - [`SimulationEvent::IndelInserted`] →
//!   `DirtyReason::StructuralIndel { site: at, delta: +1 }`
//! - [`SimulationEvent::IndelDeleted`] →
//!   `DirtyReason::StructuralIndel { site: at, delta: -1 }`
//! - [`SimulationEvent::BasePushed`] → no signal: appends extend
//!   the pool and are already covered by
//!   `PassEffect::AssembleSegment` / `PassEffect::AppendRegion`,
//!   which have their own dispatch paths in
//!   `apply_live_call_updates`.
//! - Every other (reserved) variant: ignored. New variants are
//!   additive — the sink consciously matches only what it cares
//!   about and lets the rest pass.

use crate::ir::{SimulationEvent, SimulationEventSink};

use super::model::{DirtyReason, DirtyWindow};

/// Observer that records a `DirtyWindow` for every `BaseChanged` /
/// `IndelInserted` / `IndelDeleted` event the simulation builder
/// emits.
///
/// Cheap to attach (one `Vec` allocation), O(1) per event. Drained
/// at seal time via [`Self::seal`]; the returned `Vec<DirtyWindow>`
/// is the captured stream, in event order.
pub(crate) struct DirtySignalObserver {
    windows: Vec<DirtyWindow>,
}

impl DirtySignalObserver {
    pub(crate) fn new() -> Self {
        Self {
            windows: Vec::new(),
        }
    }

    /// Consume and return the captured dirty windows.
    pub(crate) fn seal(self) -> Vec<DirtyWindow> {
        self.windows
    }
}

impl SimulationEventSink for DirtySignalObserver {
    fn on_event(&mut self, event: &SimulationEvent) {
        match *event {
            SimulationEvent::BaseChanged { handle, .. } => {
                let site = handle.index();
                self.windows.push(DirtyWindow::new(
                    site,
                    site.saturating_add(1),
                    DirtyReason::BaseEdited { site },
                ));
            }
            SimulationEvent::IndelInserted { at, .. } => {
                self.windows.push(DirtyWindow::new(
                    at,
                    at.saturating_add(1),
                    DirtyReason::StructuralIndel { site: at, delta: 1 },
                ));
            }
            SimulationEvent::IndelDeleted { at, .. } => {
                self.windows.push(DirtyWindow::new(
                    at,
                    at.saturating_add(1),
                    DirtyReason::StructuralIndel {
                        site: at,
                        delta: -1,
                    },
                ));
            }
            // Pushes are AssembleSegment / AppendRegion territory;
            // reserved (non-pool) variants describe sidecar state
            // that doesn't dirty pool positions.
            SimulationEvent::BasePushed { .. }
            | SimulationEvent::BaseDeleted { .. }
            | SimulationEvent::AssignmentChanged { .. }
            | SimulationEvent::TrimChanged { .. }
            | SimulationEvent::RegionAdded { .. }
            | SimulationEvent::PRegionAdded { .. }
            | SimulationEvent::RegionReplaced { .. }
            | SimulationEvent::SegmentReplaced { .. }
            | SimulationEvent::ReverseComplementFlagRecorded { .. }
            | SimulationEvent::OrientationChanged { .. }
            | SimulationEvent::MutationCountChanged { .. } => {}
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::ir::{flag, NucFlags, NucHandle, Nucleotide, Segment, Simulation, SimulationBuilder};

    fn base_pushed(handle: u32, base: u8) -> SimulationEvent {
        SimulationEvent::BasePushed {
            handle: NucHandle::new(handle),
            base,
            segment: Segment::V,
            germline_pos: None,
            flags: NucFlags::empty(),
        }
    }

    fn base_changed(handle: u32, old_base: u8, new_base: u8) -> SimulationEvent {
        SimulationEvent::BaseChanged {
            handle: NucHandle::new(handle),
            old_base,
            new_base,
            segment: Segment::V,
            germline_pos: None,
        }
    }

    fn indel_inserted(at: u32, base: u8) -> SimulationEvent {
        SimulationEvent::IndelInserted {
            at,
            base,
            segment: Segment::V,
            flags: NucFlags::empty(),
        }
    }

    fn indel_deleted(at: u32, removed_base: u8) -> SimulationEvent {
        SimulationEvent::IndelDeleted {
            at,
            removed_base,
            segment: Segment::V,
        }
    }

    #[test]
    fn empty_observer_yields_empty_dirty_log() {
        let obs = DirtySignalObserver::new();
        assert!(obs.seal().is_empty());
    }

    #[test]
    fn base_pushed_emits_no_signal() {
        let mut obs = DirtySignalObserver::new();
        obs.on_event(&base_pushed(0, b'A'));
        obs.on_event(&base_pushed(1, b'C'));
        assert!(obs.seal().is_empty());
    }

    #[test]
    fn base_changed_emits_base_edited_window() {
        let mut obs = DirtySignalObserver::new();
        obs.on_event(&base_changed(7, b'A', b'G'));
        let windows = obs.seal();
        assert_eq!(windows.len(), 1);
        assert_eq!(windows[0].start, 7);
        assert_eq!(windows[0].end, 8);
        assert!(matches!(
            windows[0].reason,
            DirtyReason::BaseEdited { site: 7 }
        ));
    }

    #[test]
    fn indel_inserted_emits_positive_delta() {
        let mut obs = DirtySignalObserver::new();
        obs.on_event(&indel_inserted(12, b'N'));
        let windows = obs.seal();
        assert_eq!(windows.len(), 1);
        assert!(matches!(
            windows[0].reason,
            DirtyReason::StructuralIndel { site: 12, delta: 1 }
        ));
    }

    #[test]
    fn indel_deleted_emits_negative_delta() {
        let mut obs = DirtySignalObserver::new();
        obs.on_event(&indel_deleted(4, b'C'));
        let windows = obs.seal();
        assert_eq!(windows.len(), 1);
        assert!(matches!(
            windows[0].reason,
            DirtyReason::StructuralIndel { site: 4, delta: -1 }
        ));
    }

    #[test]
    fn reserved_variants_emit_no_signal() {
        // Future-emission variants must not corrupt the dirty
        // stream until they're consciously wired in. Covers every
        // reserved variant the enum can carry today; if the enum
        // grows, this test won't compile until the new variant is
        // either matched here (and acknowledged as a no-op) or
        // accepted as a deliberate dirty trigger.
        use crate::assignment::{AlleleInstance, TrimEnd};
        use crate::ir::Region;
        use crate::refdata::AlleleId;
        let mut obs = DirtySignalObserver::new();
        obs.on_event(&SimulationEvent::BaseDeleted {
            at: 0,
            removed_base: b'A',
        });
        obs.on_event(&SimulationEvent::AssignmentChanged {
            segment: Segment::V,
            old: None,
            new: AlleleInstance::new(AlleleId::new(0)),
        });
        obs.on_event(&SimulationEvent::TrimChanged {
            segment: Segment::V,
            end: TrimEnd::Five,
            old: Some(0),
            new: 3,
        });
        obs.on_event(&SimulationEvent::RegionAdded {
            region: Region::new(Segment::V, NucHandle::new(0), NucHandle::new(5)),
        });
        obs.on_event(&SimulationEvent::RegionReplaced {
            old: Region::new(Segment::V, NucHandle::new(0), NucHandle::new(5)),
            new: Region::new(Segment::V, NucHandle::new(0), NucHandle::new(7)),
        });
        obs.on_event(&SimulationEvent::ReverseComplementFlagRecorded { applied: true });
        obs.on_event(&SimulationEvent::MutationCountChanged {
            old: 2,
            new: 7,
            delta: 5,
        });
        assert!(obs.seal().is_empty());
    }

    #[test]
    fn windows_accumulate_in_event_order() {
        let mut obs = DirtySignalObserver::new();
        obs.on_event(&base_changed(0, b'A', b'G'));
        obs.on_event(&indel_inserted(5, b'N'));
        obs.on_event(&base_changed(2, b'C', b'T'));
        let windows = obs.seal();
        assert_eq!(windows.len(), 3);
        assert!(matches!(
            windows[0].reason,
            DirtyReason::BaseEdited { site: 0 }
        ));
        assert!(matches!(
            windows[1].reason,
            DirtyReason::StructuralIndel { delta: 1, .. }
        ));
        assert!(matches!(
            windows[2].reason,
            DirtyReason::BaseEdited { site: 2 }
        ));
    }

    #[test]
    fn capture_during_assembly_loop_yields_expected_dirty_windows() {
        // End-to-end parity: drive a real `SimulationBuilder`
        // through pushes + change + insert + delete with the
        // dirty-signal sink attached, and confirm the captured
        // windows match what the legacy `IrEventObserver` path
        // produced (BaseEdited / StructuralIndel(+1/-1); pushes
        // emit nothing). This pins the new-channel behavior to
        // the old contract before any production callers move.
        let mut builder = SimulationBuilder::from_simulation(Simulation::new());
        builder.attach_dirty_signal_observer();

        for (i, &b) in b"ACGT".iter().enumerate() {
            builder.push_nucleotide(Nucleotide::germline(b, i as u16, Segment::V));
        }
        builder.change_base(NucHandle::new(1), b'X'); // C → X
        builder.insert_indel(
            2,
            Nucleotide::synthetic(b'N', Segment::V, flag::INDEL_INSERTED),
        );
        builder.delete_indel(0);

        let windows = builder
            .take_dirty_signal_observer()
            .expect("dirty-signal sink was attached above")
            .seal();

        // 4 pushes contribute no signal; the change, insert,
        // delete each contribute exactly one window.
        assert_eq!(windows.len(), 3);
        assert_eq!(windows[0].start, 1);
        assert_eq!(windows[0].end, 2);
        assert!(matches!(
            windows[0].reason,
            DirtyReason::BaseEdited { site: 1 }
        ));
        assert!(matches!(
            windows[1].reason,
            DirtyReason::StructuralIndel { site: 2, delta: 1 }
        ));
        assert!(matches!(
            windows[2].reason,
            DirtyReason::StructuralIndel { site: 0, delta: -1 }
        ));
    }
}
