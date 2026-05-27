//! `EventLogObserver` — a [`SimulationEventSink`] that captures every
//! emitted [`SimulationEvent`] into a flat `Vec`.
//!
//! Originally the only state-change consumer that explicitly recorded
//! its input rather than reacting to it. With the event/observer
//! refactor (slice 1) it migrated from the legacy `IrEventObserver`
//! trait onto the unified [`SimulationEventSink`] channel — see
//! [`crate::ir::sim_event`] for the architectural split between
//! "trace = choices" and "events = consequences."
//!
//! Replaying the captured stream through [`replay_events`] on a copy
//! of the pre-mutation `Simulation` produces a result whose pool /
//! sequence state is bit-identical to the original post-mutation
//! `Simulation`. This is the foundation for:
//!
//! - **Trace replay** — debugging tools can capture the event stream
//!   for one record and replay it later on a clean simulation to
//!   inspect intermediate state.
//! - **MCP audit log** — the JSON-friendly [`SimulationEvent`] enum
//!   can be serialised so the MCP server can hand the model a "what
//!   actually happened" log alongside the existing decision trace.
//! - **Sampling-bias metrics** — counting events of each kind per
//!   pass is a small downstream consumer.

use super::sim_event::{SimulationEvent, SimulationEventSink};

#[cfg(test)]
use super::nucleotide::GermlinePos;
#[cfg(test)]
use super::Nucleotide;
#[cfg(test)]
use super::{NucHandle, Segment, Simulation, SimulationBuilder};

/// Sink that captures every emitted [`SimulationEvent`] into a flat
/// vec.
///
/// Cheap to attach (a fresh `Vec` allocation), cheap per event
/// (`Vec::push` of a small enum). Drained at seal time; the
/// returned `Vec<SimulationEvent>` is the captured stream in
/// emission order.
pub struct EventLogObserver {
    events: Vec<SimulationEvent>,
}

impl EventLogObserver {
    pub fn new() -> Self {
        Self { events: Vec::new() }
    }

    /// Consume and return the captured event vec.
    pub fn seal(self) -> Vec<SimulationEvent> {
        self.events
    }
}

impl Default for EventLogObserver {
    fn default() -> Self {
        Self::new()
    }
}

impl SimulationEventSink for EventLogObserver {
    fn on_event(&mut self, event: &SimulationEvent) {
        self.events.push(event.clone());
    }
}

/// Replay a captured event stream on top of `start` via a fresh
/// `SimulationBuilder`. Returns the post-replay `Simulation`.
///
/// Test-only today — used to validate that the [`SimulationEvent`]
/// enum captures enough state to reproduce a mutation stream
/// byte-for-byte. A future MCP audit-log consumer would promote
/// this back to `pub` once it has an actual caller.
///
/// The semantics are:
///
/// - `BasePushed` → `builder.push_nucleotide(n)` where `n` is
///   reconstructed from the event's segment / germline_pos / flags.
/// - `BaseChanged` → `builder.change_base(handle, new_base)`.
/// - `IndelInserted` → `builder.insert_indel(at, n)`.
/// - `IndelDeleted` → `builder.delete_indel(at)`.
/// - Reserved variants (`BaseDeleted`, `RegionChanged`,
///   `AssignmentChanged`, `ReverseComplemented`,
///   `MutationCountChanged`) are not yet emitted by the builder
///   and are skipped by the replay path.
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
pub(crate) fn replay_events(start: Simulation, events: &[SimulationEvent]) -> Simulation {
    let mut builder = SimulationBuilder::from_simulation(start);
    for event in events {
        replay_one(&mut builder, event);
    }
    builder.seal()
}

#[cfg(test)]
fn replay_one(builder: &mut SimulationBuilder, event: &SimulationEvent) {
    match event {
        SimulationEvent::BasePushed {
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
        SimulationEvent::BaseChanged {
            handle,
            new_base,
            ..
        } => {
            // segment / germline_pos / old_base are unchanged by
            // `change_base`; replay only needs the target handle
            // and the new byte.
            builder.change_base(*handle, *new_base);
        }
        SimulationEvent::IndelInserted {
            at,
            base,
            segment,
            flags,
        } => {
            let n = Nucleotide::synthetic(*base, *segment, *flags);
            builder.insert_indel(*at, n);
        }
        SimulationEvent::IndelDeleted {
            at,
            removed_base: _,
            segment: _,
        } => {
            builder.delete_indel(*at);
        }
        // Reserved (non-pool) variants don't drive replay. Replay
        // is, by documented policy, **pool-mutation only** — the
        // non-pool consequence variants describe sidecar state
        // (assignments, trims, regions, rev-comp, mutation
        // counters) which is reconstructed by the passes that own
        // it, not by the event stream. Skipping them here keeps
        // `replay_events` honest: if a future caller wires a
        // non-pool consumer (e.g. an audit-log rebuilder) it
        // belongs in a sibling replay function, not by widening
        // this one.
        SimulationEvent::BaseDeleted { .. }
        | SimulationEvent::AssignmentChanged { .. }
        | SimulationEvent::TrimChanged { .. }
        | SimulationEvent::RegionAdded { .. }
        | SimulationEvent::RegionReplaced { .. }
        | SimulationEvent::ReverseComplementFlagRecorded { .. }
        | SimulationEvent::MutationCountChanged { .. } => {}
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::ir::NucFlags;

    #[test]
    fn empty_observer_yields_empty_log() {
        let obs = EventLogObserver::new();
        let log = obs.seal();
        assert!(log.is_empty());
    }

    #[test]
    fn captures_base_pushed_events_via_sink_trait() {
        let mut obs = EventLogObserver::new();
        obs.on_event(&SimulationEvent::BasePushed {
            handle: NucHandle::new(0),
            base: b'A',
            segment: Segment::V,
            germline_pos: None,
            flags: NucFlags::empty(),
        });
        obs.on_event(&SimulationEvent::BasePushed {
            handle: NucHandle::new(1),
            base: b'C',
            segment: Segment::V,
            germline_pos: None,
            flags: NucFlags::empty(),
        });
        let log = obs.seal();
        assert_eq!(log.len(), 2);
        match &log[0] {
            SimulationEvent::BasePushed { handle, base, .. } => {
                assert_eq!(*handle, NucHandle::new(0));
                assert_eq!(*base, b'A');
            }
            other => panic!("expected BasePushed, got {:?}", other),
        }
    }

    #[test]
    fn captures_all_wired_event_kinds_via_sink_trait() {
        let mut obs = EventLogObserver::new();
        for ev in [
            SimulationEvent::BasePushed {
                handle: NucHandle::new(0),
                base: b'A',
                segment: Segment::V,
                germline_pos: None,
                flags: NucFlags::empty(),
            },
            SimulationEvent::BaseChanged {
                handle: NucHandle::new(0),
                old_base: b'A',
                new_base: b'G',
                segment: Segment::V,
                germline_pos: None,
            },
            SimulationEvent::IndelInserted {
                at: 1,
                base: b'C',
                segment: Segment::V,
                flags: NucFlags::empty(),
            },
            SimulationEvent::IndelDeleted {
                at: 0,
                removed_base: b'G',
                segment: Segment::V,
            },
        ] {
            obs.on_event(&ev);
        }
        let log = obs.seal();
        assert_eq!(log.len(), 4);
        assert!(matches!(log[0], SimulationEvent::BasePushed { .. }));
        assert!(matches!(log[1], SimulationEvent::BaseChanged { .. }));
        assert!(matches!(log[2], SimulationEvent::IndelInserted { .. }));
        assert!(matches!(log[3], SimulationEvent::IndelDeleted { .. }));
    }

    #[test]
    fn records_every_reserved_payload_variant_verbatim() {
        // The event log observer is the universal recorder — when
        // a future builder method starts emitting a non-pool
        // event, the recorded stream must round-trip through
        // `seal()` byte-for-byte without the recorder having to
        // be updated. This test pins that property NOW so the
        // payload shapes can't drift when emitters are added.
        use crate::assignment::{AlleleInstance, TrimEnd};
        use crate::ir::Region;
        use crate::refdata::AlleleId;
        let mut obs = EventLogObserver::new();
        let events = vec![
            SimulationEvent::BaseDeleted {
                at: 9,
                removed_base: b'T',
            },
            SimulationEvent::AssignmentChanged {
                segment: Segment::V,
                old: None,
                new: AlleleInstance::new(AlleleId::new(2)),
            },
            SimulationEvent::TrimChanged {
                segment: Segment::J,
                end: TrimEnd::Three,
                old: Some(0),
                new: 4,
            },
            SimulationEvent::RegionAdded {
                region: Region::new(Segment::V, NucHandle::new(0), NucHandle::new(10)),
            },
            SimulationEvent::RegionReplaced {
                old: Region::new(Segment::V, NucHandle::new(0), NucHandle::new(10)),
                new: Region::new(Segment::V, NucHandle::new(0), NucHandle::new(13)),
            },
            SimulationEvent::ReverseComplementFlagRecorded { applied: false },
            SimulationEvent::MutationCountChanged {
                old: 0,
                new: 7,
                delta: 7,
            },
        ];
        for ev in &events {
            obs.on_event(ev);
        }
        let log = obs.seal();
        // Equality (not just `matches!`) — exercises Eq on every
        // payload field and on the embedded `Region`.
        assert_eq!(log, events);
    }

    #[test]
    fn replay_skips_non_pool_events_by_documented_policy() {
        // `replay_events` is pool-mutation only. A captured stream
        // that interleaves non-pool variants between pool events
        // must produce the same final `Simulation` as if the non-
        // pool variants had been omitted — they're inert.
        use crate::assignment::{AlleleInstance, TrimEnd};
        use crate::ir::Region;
        use crate::refdata::AlleleId;
        let with_non_pool = vec![
            SimulationEvent::BasePushed {
                handle: NucHandle::new(0),
                base: b'A',
                segment: Segment::V,
                germline_pos: Some(0),
                flags: NucFlags::empty(),
            },
            SimulationEvent::AssignmentChanged {
                segment: Segment::V,
                old: None,
                new: AlleleInstance::new(AlleleId::new(0)),
            },
            SimulationEvent::TrimChanged {
                segment: Segment::V,
                end: TrimEnd::Five,
                old: Some(0),
                new: 1,
            },
            SimulationEvent::BasePushed {
                handle: NucHandle::new(1),
                base: b'C',
                segment: Segment::V,
                germline_pos: Some(1),
                flags: NucFlags::empty(),
            },
            SimulationEvent::RegionAdded {
                region: Region::new(Segment::V, NucHandle::new(0), NucHandle::new(2)),
            },
            SimulationEvent::ReverseComplementFlagRecorded { applied: true },
            SimulationEvent::MutationCountChanged {
                old: 0,
                new: 3,
                delta: 3,
            },
        ];
        let pool_only = vec![
            SimulationEvent::BasePushed {
                handle: NucHandle::new(0),
                base: b'A',
                segment: Segment::V,
                germline_pos: Some(0),
                flags: NucFlags::empty(),
            },
            SimulationEvent::BasePushed {
                handle: NucHandle::new(1),
                base: b'C',
                segment: Segment::V,
                germline_pos: Some(1),
                flags: NucFlags::empty(),
            },
        ];
        let a = replay_events(Simulation::new(), &with_non_pool);
        let b = replay_events(Simulation::new(), &pool_only);
        let a_bases: Vec<u8> = a.pool.as_slice().iter().map(|n| n.base).collect();
        let b_bases: Vec<u8> = b.pool.as_slice().iter().map(|n| n.base).collect();
        assert_eq!(a_bases, b_bases);
        assert_eq!(a.pool.len(), b.pool.len());
    }

    #[test]
    fn replay_pushes_match_original() {
        let events = vec![
            SimulationEvent::BasePushed {
                handle: NucHandle::new(0),
                base: b'A',
                segment: Segment::V,
                germline_pos: Some(0),
                flags: NucFlags::empty(),
            },
            SimulationEvent::BasePushed {
                handle: NucHandle::new(1),
                base: b'C',
                segment: Segment::V,
                germline_pos: Some(1),
                flags: NucFlags::empty(),
            },
            SimulationEvent::BasePushed {
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
        // End-to-end parity test: drive a `SimulationBuilder`
        // through a sequence of pushes + a change + an indel
        // insert + delete, with an event-log sink attached. Then
        // replay the captured stream on a fresh empty Simulation
        // and verify the resulting pool matches the original
        // sealed simulation byte-for-byte.
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

        // 4 pushes + 1 change + 1 insert + 1 delete = 7 events.
        assert_eq!(captured.len(), 7);

        let replayed = replay_events(Simulation::new(), &captured);
        let original_bases: Vec<u8> = original_sim
            .pool
            .as_slice()
            .iter()
            .map(|n| n.base)
            .collect();
        let replayed_bases: Vec<u8> = replayed.pool.as_slice().iter().map(|n| n.base).collect();
        assert_eq!(
            original_bases, replayed_bases,
            "replayed pool bases diverge from original",
        );
        assert_eq!(original_sim.pool.len(), replayed.pool.len());
    }

    #[test]
    fn replay_changes_match_original() {
        let mut start = Simulation::new();
        for (i, &b) in b"ACGT".iter().enumerate() {
            let (next, _) =
                start.with_nucleotide_pushed(Nucleotide::germline(b, i as u16, Segment::V));
            start = next;
        }

        let events = vec![SimulationEvent::BaseChanged {
            handle: NucHandle::new(1),
            old_base: b'C',
            new_base: b'T',
            segment: Segment::V,
            germline_pos: Some(1),
        }];
        let result = replay_events(start, &events);
        let bases: Vec<u8> = result.pool.as_slice().iter().map(|n| n.base).collect();
        assert_eq!(&bases, b"ATGT");
    }
}
