//! Phase 13: structured event-log observer.
//!
//! Captures every IR mutation event emitted by `SimulationBuilder`
//! into a flat `Vec<IrEvent>`. Replaying the resulting event stream
//! through [`replay_events`] on a copy of the pre-mutation
//! `Simulation` produces a result whose pool/sequence state is
//! bit-identical to the original post-mutation `Simulation`.
//!
//! This is the foundation for:
//!
//! - **Trace replay** — debugging tools can capture the IR event
//!   stream for one record and replay it later on a clean simulation
//!   to inspect intermediate state.
//! - **MCP audit log** — the JSON-friendly `IrEvent` enum can be
//!   serialised so the MCP server can hand the model a "what
//!   actually happened" log alongside the existing decision trace.
//! - **Sampling-bias metrics** — counting events of each kind per
//!   pass is a small downstream consumer.
//!
//! The architectural significance: this is the **fourth concrete
//! `IrEventObserver` impl** (after walker, codon rail, productive
//! admit-mask), demonstrating that the trait system supports new
//! consumers by **pure registration** — no changes to the trait
//! surface, no changes to the existing observers, no changes to the
//! builder's broadcast logic. Adding a fifth observer (e.g. a
//! per-segment freshness tracker, or a per-allele evidence
//! summariser) follows the same template.

use super::nucleotide::NucFlags;
use super::{NucHandle, Nucleotide, Segment};
use crate::ir::builder::IrEventObserver;

#[cfg(test)]
use super::nucleotide::GermlinePos;
#[cfg(test)]
use super::{Simulation, SimulationBuilder};

/// One IR mutation event captured by [`EventLogObserver`].
///
/// Variants mirror the methods on [`IrEventObserver`]. Each carries
/// just enough state to be replayable through
/// [`SimulationBuilder`]'s mutation methods without needing the
/// pre-mutation simulation as auxiliary context.
#[derive(Clone, Debug, Eq, PartialEq)]
pub enum IrEvent {
    /// Emitted by `SimulationBuilder::push_nucleotide`.
    BasePushed {
        handle: NucHandle,
        base: u8,
        segment: Segment,
        germline_pos: Option<u16>,
        flags: NucFlags,
    },
    /// Emitted by `SimulationBuilder::change_base`.
    BaseChanged {
        handle: NucHandle,
        old_base: u8,
        new_base: u8,
    },
    /// Emitted by `SimulationBuilder::insert_indel`.
    IndelInserted {
        at: u32,
        base: u8,
        segment: Segment,
        flags: NucFlags,
    },
    /// Emitted by `SimulationBuilder::delete_indel`.
    IndelDeleted {
        at: u32,
        removed_base: u8,
    },
}

/// Observer that captures every emitted event into a flat vec.
///
/// Cheap to attach (a fresh `Vec` allocation), cheap per event
/// (`Vec::push` of a small enum). Drained at seal time; the
/// returned `Vec<IrEvent>` is the captured stream.
pub struct EventLogObserver {
    events: Vec<IrEvent>,
}

impl EventLogObserver {
    pub fn new() -> Self {
        Self {
            events: Vec::new(),
        }
    }

    /// Consume and return the captured event vec.
    pub fn seal(self) -> Vec<IrEvent> {
        self.events
    }
}

impl Default for EventLogObserver {
    fn default() -> Self {
        Self::new()
    }
}

impl IrEventObserver for EventLogObserver {
    fn on_base_pushed(&mut self, handle: NucHandle, n: &Nucleotide) {
        self.events.push(IrEvent::BasePushed {
            handle,
            base: n.base,
            segment: n.segment,
            germline_pos: n.germline_pos.get(),
            flags: n.flags,
        });
    }

    fn on_base_changed(&mut self, handle: NucHandle, old_n: &Nucleotide, new_base: u8) {
        self.events.push(IrEvent::BaseChanged {
            handle,
            old_base: old_n.base,
            new_base,
        });
    }

    fn on_indel_inserted(&mut self, at: u32, n: &Nucleotide) {
        self.events.push(IrEvent::IndelInserted {
            at,
            base: n.base,
            segment: n.segment,
            flags: n.flags,
        });
    }

    fn on_indel_deleted(&mut self, at: u32, removed: &Nucleotide) {
        self.events.push(IrEvent::IndelDeleted {
            at,
            removed_base: removed.base,
        });
    }
}

/// Replay a captured event stream on top of `start` via a fresh
/// `SimulationBuilder`. Returns the post-replay `Simulation`.
///
/// Currently used by tests to validate that the `IrEvent` enum
/// captures enough state to reproduce a mutation stream byte-for-
/// byte. A future MCP audit-log consumer would promote this back
/// to `pub` once it has an actual caller.
///
/// The semantics are:
///
/// - `BasePushed` → `builder.push_nucleotide(n)` where `n` is
///   reconstructed from the event's segment / germline_pos / flags.
/// - `BaseChanged` → `builder.change_base(handle, new_base)`.
/// - `IndelInserted` → `builder.insert_indel(at, n)`.
/// - `IndelDeleted` → `builder.delete_indel(at)`.
///
/// No observers are attached during replay: the goal is to
/// reproduce the post-mutation `Simulation`, not to re-emit the
/// captured events.
///
/// This replays the *raw IR mutation stream* only — pass-side
/// effects like region creation, allele assignment, trim updates,
/// or live-call snapshots happen at pass boundaries via
/// non-builder APIs and would need their own event surface.
#[cfg(test)]
pub(crate) fn replay_events(start: Simulation, events: &[IrEvent]) -> Simulation {
    let mut builder = SimulationBuilder::from_simulation(start);
    for event in events {
        replay_one(&mut builder, event);
    }
    builder.seal()
}

#[cfg(test)]
fn replay_one(builder: &mut SimulationBuilder, event: &IrEvent) {
    match event {
        IrEvent::BasePushed {
            handle: _,
            base,
            segment,
            germline_pos,
            flags,
        } => {
            // Reconstruct the nucleotide field-by-field — neither
            // `Nucleotide::germline` nor `::synthetic` lets us set
            // germline_pos and flags independently in one call.
            // (germline == base because base-changes only happen
            // *after* push, and at push time those always agree.)
            let n = Nucleotide {
                base: *base,
                germline: *base,
                germline_pos: match germline_pos {
                    Some(g) => GermlinePos::pos(*g),
                    None => GermlinePos::NONE,
                },
                segment: *segment,
                flags: *flags,
            };
            builder.push_nucleotide(n);
        }
        IrEvent::BaseChanged {
            handle,
            old_base: _,
            new_base,
        } => {
            builder.change_base(*handle, *new_base);
        }
        IrEvent::IndelInserted {
            at,
            base,
            segment,
            flags,
        } => {
            let n = Nucleotide::synthetic(*base, *segment, *flags);
            builder.insert_indel(*at, n);
        }
        IrEvent::IndelDeleted {
            at,
            removed_base: _,
        } => {
            builder.delete_indel(*at);
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::ir::flag;

    fn n(base: u8) -> Nucleotide {
        Nucleotide::synthetic(base, Segment::V, flag::N_NUC)
    }

    #[test]
    fn empty_observer_yields_empty_log() {
        let obs = EventLogObserver::new();
        let log = obs.seal();
        assert!(log.is_empty());
    }

    #[test]
    fn captures_base_pushed_events() {
        let mut obs = EventLogObserver::new();
        IrEventObserver::on_base_pushed(&mut obs, NucHandle::new(0), &n(b'A'));
        IrEventObserver::on_base_pushed(&mut obs, NucHandle::new(1), &n(b'C'));
        let log = obs.seal();
        assert_eq!(log.len(), 2);
        match &log[0] {
            IrEvent::BasePushed { handle, base, .. } => {
                assert_eq!(*handle, NucHandle::new(0));
                assert_eq!(*base, b'A');
            }
            other => panic!("expected BasePushed, got {:?}", other),
        }
    }

    #[test]
    fn captures_all_event_kinds() {
        let mut obs = EventLogObserver::new();
        IrEventObserver::on_base_pushed(&mut obs, NucHandle::new(0), &n(b'A'));
        IrEventObserver::on_base_changed(&mut obs, NucHandle::new(0), &n(b'A'), b'G');
        IrEventObserver::on_indel_inserted(&mut obs, 1, &n(b'C'));
        IrEventObserver::on_indel_deleted(&mut obs, 0, &n(b'G'));
        let log = obs.seal();
        assert_eq!(log.len(), 4);
        assert!(matches!(log[0], IrEvent::BasePushed { .. }));
        assert!(matches!(log[1], IrEvent::BaseChanged { .. }));
        assert!(matches!(log[2], IrEvent::IndelInserted { .. }));
        assert!(matches!(log[3], IrEvent::IndelDeleted { .. }));
    }

    #[test]
    fn replay_pushes_match_original() {
        let events = vec![
            IrEvent::BasePushed {
                handle: NucHandle::new(0),
                base: b'A',
                segment: Segment::V,
                germline_pos: Some(0),
                flags: NucFlags::empty(),
            },
            IrEvent::BasePushed {
                handle: NucHandle::new(1),
                base: b'C',
                segment: Segment::V,
                germline_pos: Some(1),
                flags: NucFlags::empty(),
            },
            IrEvent::BasePushed {
                handle: NucHandle::new(2),
                base: b'G',
                segment: Segment::V,
                germline_pos: Some(2),
                flags: NucFlags::empty(),
            },
        ];
        let result = replay_events(Simulation::new(), &events);
        assert_eq!(result.pool.len(), 3);
        let bases: Vec<u8> = result.pool.as_slice().iter().map(|n| n.base).collect();
        assert_eq!(&bases, b"ACG");
    }

    #[test]
    fn capture_during_assembly_loop_then_replay_yields_same_pool() {
        // End-to-end: drive a `SimulationBuilder` through a sequence
        // of pushes + a change + an indel insert + delete, with an
        // event-log observer attached. Then replay the captured
        // stream on a fresh empty Simulation and verify the
        // resulting pool matches the original sealed simulation
        // byte-for-byte.
        // Original simulation: push 4 bases, change one, insert
        // one, delete one. Capture every event.
        let mut builder = SimulationBuilder::from_simulation(Simulation::new());
        builder.attach_event_log_observer();
        for (i, &b) in b"ACGT".iter().enumerate() {
            builder.push_nucleotide(Nucleotide::germline(b, i as u16, Segment::V));
        }
        builder.change_base(NucHandle::new(1), b'X'); // C → X
        builder.insert_indel(
            2,
            Nucleotide::synthetic(b'N', Segment::V, crate::ir::flag::INDEL_INSERTED),
        );
        builder.delete_indel(0);

        let captured = builder.seal_event_log_observer();
        let original_sim = builder.seal();

        // The captured stream is 4 pushes + 1 change + 1 insert + 1 delete = 7 events.
        assert_eq!(captured.len(), 7);

        // Replay on a fresh empty Simulation. Result's pool should
        // match the original's pool byte-for-byte.
        let replayed = replay_events(Simulation::new(), &captured);
        let original_bases: Vec<u8> =
            original_sim.pool.as_slice().iter().map(|n| n.base).collect();
        let replayed_bases: Vec<u8> =
            replayed.pool.as_slice().iter().map(|n| n.base).collect();
        assert_eq!(
            original_bases, replayed_bases,
            "replayed pool bases diverge from original"
        );
        assert_eq!(original_sim.pool.len(), replayed.pool.len());
    }

    #[test]
    fn replay_changes_match_original() {
        // Start with a sim that has bases ACGT; capture an edit
        // of pos 1 (C → T) via the observer; replay; verify.
        let mut start = Simulation::new();
        for (i, &b) in b"ACGT".iter().enumerate() {
            let (next, _) = start.with_nucleotide_pushed(Nucleotide::germline(
                b,
                i as u16,
                Segment::V,
            ));
            start = next;
        }

        let events = vec![IrEvent::BaseChanged {
            handle: NucHandle::new(1),
            old_base: b'C',
            new_base: b'T',
        }];
        let result = replay_events(start, &events);
        let bases: Vec<u8> = result.pool.as_slice().iter().map(|n| n.base).collect();
        assert_eq!(&bases, b"ATGT");
    }
}
