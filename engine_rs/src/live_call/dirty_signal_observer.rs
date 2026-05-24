//! Streaming dirty-signal observer.
//!
//! Captures `DirtyWindow` change signals from the IR event stream as
//! mutations happen, so the post-pass live-call refresh in
//! `compiled/execute.rs::apply_live_call_updates` can narrow to
//! segments whose region overlaps a dirty position.
//!
//! ## Why an observer rather than ad-hoc tracking on the builder
//!
//! Same reason every other reactive consumer in the engine is an
//! observer: it composes cleanly with attach/broadcast/seal, it's
//! opt-in per pass (the assembly path doesn't need it), and it keeps
//! the dirty-signal vocabulary co-located with the rest of the
//! live-call module.
//!
//! ## Which events produce dirty signals
//!
//! - `on_base_changed(handle, ..)` →
//!   `DirtyReason::BaseEdited { site: handle.index() }`
//! - `on_indel_inserted(at, ..)` →
//!   `DirtyReason::StructuralIndel { site: at, delta: +1 }`
//! - `on_indel_deleted(at, ..)` →
//!   `DirtyReason::StructuralIndel { site: at, delta: -1 }`
//! - `on_base_pushed` → no signal: appends extend the pool and are
//!   already covered by `PassEffect::AssembleSegment` /
//!   `PassEffect::AppendRegion`, which have their own dispatch paths
//!   in `apply_live_call_updates`.

use crate::ir::builder::IrEventObserver;
use crate::ir::{NucHandle, Nucleotide};

use super::model::{DirtyReason, DirtyWindow};

/// Observer that records a `DirtyWindow` for every `change_base` /
/// `insert_indel` / `delete_indel` event a `SimulationBuilder` emits.
///
/// Cheap to attach (one `Vec` allocation), O(1) per event. Drained at
/// seal time via [`Self::seal`]; the returned `Vec<DirtyWindow>` is
/// the captured stream, in event order.
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

impl IrEventObserver for DirtySignalObserver {
    fn on_base_pushed(&mut self, _handle: NucHandle, _n: &Nucleotide) {
        // Pushes are handled by AssembleSegment / AppendRegion
        // dispatch; we record no signal here.
    }

    fn on_base_changed(&mut self, handle: NucHandle, _old_n: &Nucleotide, _new_base: u8) {
        let site = handle.index();
        self.windows.push(DirtyWindow::new(
            site,
            site.saturating_add(1),
            DirtyReason::BaseEdited { site },
        ));
    }

    fn on_indel_inserted(&mut self, at: u32, _n: &Nucleotide) {
        self.windows.push(DirtyWindow::new(
            at,
            at.saturating_add(1),
            DirtyReason::StructuralIndel { site: at, delta: 1 },
        ));
    }

    fn on_indel_deleted(&mut self, at: u32, _removed: &Nucleotide) {
        self.windows.push(DirtyWindow::new(
            at,
            at.saturating_add(1),
            DirtyReason::StructuralIndel { site: at, delta: -1 },
        ));
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::ir::{flag, Segment};

    fn n(base: u8) -> Nucleotide {
        Nucleotide::synthetic(base, Segment::V, flag::N_NUC)
    }

    #[test]
    fn empty_observer_yields_empty_dirty_log() {
        let obs = DirtySignalObserver::new();
        assert!(obs.seal().is_empty());
    }

    #[test]
    fn base_pushed_emits_no_signal() {
        let mut obs = DirtySignalObserver::new();
        IrEventObserver::on_base_pushed(&mut obs, NucHandle::new(0), &n(b'A'));
        IrEventObserver::on_base_pushed(&mut obs, NucHandle::new(1), &n(b'C'));
        assert!(obs.seal().is_empty());
    }

    #[test]
    fn base_changed_emits_base_edited_window() {
        let mut obs = DirtySignalObserver::new();
        IrEventObserver::on_base_changed(&mut obs, NucHandle::new(7), &n(b'A'), b'G');
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
        IrEventObserver::on_indel_inserted(&mut obs, 12, &n(b'N'));
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
        IrEventObserver::on_indel_deleted(&mut obs, 4, &n(b'C'));
        let windows = obs.seal();
        assert_eq!(windows.len(), 1);
        assert!(matches!(
            windows[0].reason,
            DirtyReason::StructuralIndel { site: 4, delta: -1 }
        ));
    }

    #[test]
    fn windows_accumulate_in_event_order() {
        let mut obs = DirtySignalObserver::new();
        IrEventObserver::on_base_changed(&mut obs, NucHandle::new(0), &n(b'A'), b'G');
        IrEventObserver::on_indel_inserted(&mut obs, 5, &n(b'N'));
        IrEventObserver::on_base_changed(&mut obs, NucHandle::new(2), &n(b'C'), b'T');
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
}
