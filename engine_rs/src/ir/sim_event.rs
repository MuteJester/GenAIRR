//! Typed internal simulation-event channel.
//!
//! # Two event layers in the engine
//!
//! The engine has two distinct event surfaces — easy to confuse, so
//! the difference is spelled out here:
//!
//! 1. **[`crate::trace::Trace`] — stochastic choices.** What the user
//!    proposed or sampled at each addressed decision point. Persisted
//!    via [`crate::trace_file::TraceFile`]. The replay layer
//!    consumes it.
//!
//! 2. **[`SimulationEvent`] — state-change consequences.** What
//!    actually happened to the in-progress `Simulation` as a result
//!    of a proposed choice being applied. Emitted by
//!    [`crate::ir::SimulationBuilder`] at each mutation site;
//!    consumed by [`SimulationEventSink`] implementors (today: only
//!    the event-log observer; future: live-call walker, dirty-window
//!    refresher, codon-rail rebuilder).
//!
//! Spoken in two words:
//!
//! - **Trace** = "what was proposed/chosen?"
//! - **SimulationEvent** = "what changed?"
//!
//! Keeping these separated means every downstream consumer of
//! state-change observability talks to one coherent channel rather
//! than a scattered set of per-observer trait methods on the
//! builder. The migration to this channel is staged — see the
//! "Migration status" section below.
//!
//! # Variants
//!
//! Four variants are **emitted today** (the existing
//! `SimulationBuilder` mutation methods all funnel through them):
//!
//! - [`SimulationEvent::BasePushed`]
//! - [`SimulationEvent::BaseChanged`]
//! - [`SimulationEvent::IndelInserted`]
//! - [`SimulationEvent::IndelDeleted`]
//!
//! Seven variants are **reserved for future emission** — they have
//! a home in the enum with full payloads so consumers can match
//! exhaustively, but no builder method emits them yet. Each one
//! describes a real consequence (not a marker); when its emitting
//! code path is migrated, it carries the data a downstream sink
//! needs to react without having to re-read the simulation:
//!
//! - [`SimulationEvent::BaseDeleted`] — end-loss / primer-trim
//!   removal of a base. Semantically distinct from `IndelDeleted`
//!   (different biological provenance) even though the underlying
//!   IR primitive is the same today.
//! - [`SimulationEvent::AssignmentChanged`] — V/D/J allele
//!   assignment installed or replaced, carrying both the prior
//!   instance (if any) and the new one.
//! - [`SimulationEvent::TrimChanged`] — a V/D/J trim length was
//!   updated at one end, carrying the prior value and the new one.
//! - [`SimulationEvent::RegionAdded`] — a new structural region
//!   was appended to `Simulation::sequence`.
//! - [`SimulationEvent::RegionReplaced`] — an existing region was
//!   extended or otherwise replaced in place, carrying both the
//!   pre- and post-state regions.
//! - [`SimulationEvent::ReverseComplementFlagRecorded`] — a
//!   reverse-complement decision was committed to the simulation
//!   sidecar with `applied = true | false`.
//! - [`SimulationEvent::MutationCountChanged`] — the `n_mutations`
//!   sidecar was updated, carrying `old`, `new`, and the signed
//!   `delta` for sinks that want either the absolute count or the
//!   incremental change.
//!
//! Adding a new variant is **additive** w.r.t. the on-disk
//! [`crate::trace_file::TraceFile`] contract (that surface persists
//! `ChoiceValue`, not `SimulationEvent`), so no schema bump is
//! required.
//!
//! # Observer architecture
//!
//! Every derived-state consumer on
//! [`crate::ir::SimulationBuilder`] rides this single channel:
//!
//! - [`crate::ir::event_log_observer::EventLogObserver`] — recorder
//!   for trace replay / MCP audit.
//! - [`crate::live_call::dirty_signal_observer::DirtySignalObserver`]
//!   — refresher for segment-overlap-gated live-call updates.
//! - [`crate::live_call::walker_observer::WalkerObserverState`] —
//!   derived per-allele scoring state for the live-call walker.
//! - [`crate::contract::admit_mask_observer::ProductiveAdmitMaskObserver`]
//!   — derived NP-slot index for the productive admit-mask fast
//!   path.
//!
//! The legacy `IrEventObserver` trait was a transitional dual
//! channel used during the four-slice observer migration; it was
//! retired in the admit-mask slice, leaving `SimulationEventSink`
//! as the sole fanout surface.

use super::nucleotide::NucFlags;
use super::region::Region;
use super::{NucHandle, Segment};
use crate::assignment::{AlleleInstance, TrimEnd};

/// One typed state-change event emitted by the simulation builder.
///
/// See the module docs for the four-emitted / five-reserved split.
#[derive(Clone, Debug, Eq, PartialEq)]
pub enum SimulationEvent {
    /// A nucleotide was pushed onto the end of the pool.
    BasePushed {
        handle: NucHandle,
        base: u8,
        segment: Segment,
        germline_pos: Option<u16>,
        flags: NucFlags,
    },
    /// The base byte at `handle` was replaced. `segment` and
    /// `germline_pos` describe the nucleotide at `handle` —
    /// they are unchanged by the substitution (only the base byte
    /// flips), so the same values apply both before and after. The
    /// `flags` field is not carried because no consumer needs it.
    BaseChanged {
        handle: NucHandle,
        old_base: u8,
        new_base: u8,
        segment: Segment,
        germline_pos: Option<u16>,
    },
    /// **Reserved for future emission.** End-loss / primer-trim
    /// removal of a base. Today the underlying IR primitive is the
    /// same as `IndelDeleted`; this variant is scaffolding for the
    /// future separation when end-loss gains its own event path.
    BaseDeleted { at: u32, removed_base: u8 },
    /// A nucleotide was inserted into the pool at position `at`;
    /// every handle `>= at` shifts up by 1.
    IndelInserted {
        at: u32,
        base: u8,
        segment: Segment,
        flags: NucFlags,
    },
    /// The nucleotide at position `at` was removed from the pool;
    /// every handle `> at` shifts down by 1. `segment` records the
    /// biological segment of the removed nucleotide (symmetric with
    /// `IndelInserted.segment`) so downstream per-segment counters
    /// don't need to reach back into pre-event simulation state.
    IndelDeleted {
        at: u32,
        removed_base: u8,
        segment: Segment,
    },
    /// **Reserved for future emission.** A V/D/J allele assignment
    /// was installed or replaced. `old == None` on first install;
    /// `old == Some(prior)` when one assignment displaces another
    /// (e.g. a re-call slot replacement).
    AssignmentChanged {
        segment: Segment,
        old: Option<AlleleInstance>,
        new: AlleleInstance,
    },
    /// **Reserved for future emission.** One side's trim length on
    /// `segment` was updated. `old == None` if the assignment is
    /// being established with an initial trim and no prior value
    /// existed; otherwise the prior trim byte.
    TrimChanged {
        segment: Segment,
        end: TrimEnd,
        old: Option<u16>,
        new: u16,
    },
    /// **Reserved for future emission.** A new structural region
    /// was appended to `Simulation::sequence`. Carries the freshly-
    /// inserted `Region` so a downstream consumer that needs the
    /// segment/start/end/frame_phase doesn't have to scan the
    /// regions list.
    RegionAdded { region: Region },
    /// **Reserved for future emission.** An existing region was
    /// replaced in place (extended, frame-rotated, etc). Carries
    /// both the pre- and post-state regions so consumers can diff
    /// without re-reading the simulation.
    RegionReplaced { old: Region, new: Region },
    /// **Reserved for future emission.** A reverse-complement
    /// decision was committed to the simulation sidecar. `applied`
    /// records the resolved boolean — `true` when the strand was
    /// flipped, `false` when the sampler/contract decided not to
    /// apply.
    ReverseComplementFlagRecorded { applied: bool },
    /// **Reserved for future emission.** The `n_mutations`
    /// sidecar was updated. `old` / `new` carry the absolute
    /// values; `delta = (new as i32) - (old as i32)` is provided
    /// pre-computed so a sink that only wants the increment
    /// doesn't have to widen the operands itself.
    MutationCountChanged { old: u32, new: u32, delta: i32 },
}

/// Consumer of [`SimulationEvent`] notifications. Implementations
/// receive every event the builder emits, in emission order.
///
/// Sinks own their reaction: a recorder appends to a `Vec`, a
/// dirty-window tracker accumulates windows, a live-call walker
/// updates per-allele scores incrementally. The trait is
/// deliberately a single method so the sink doesn't have to know
/// what kinds of events exist when it's authored — it pattern-
/// matches what it cares about and ignores the rest.
pub trait SimulationEventSink {
    fn on_event(&mut self, event: &SimulationEvent);
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn event_variants_are_clone_and_eq() {
        // Pin that all five wired variants survive a clone + Eq
        // round-trip. Reserved variants get a smoke check too.
        let cases = [
            SimulationEvent::BasePushed {
                handle: NucHandle::new(7),
                base: b'A',
                segment: Segment::V,
                germline_pos: Some(12),
                flags: NucFlags::empty(),
            },
            SimulationEvent::BaseChanged {
                handle: NucHandle::new(7),
                old_base: b'A',
                new_base: b'G',
                segment: Segment::V,
                germline_pos: Some(12),
            },
            SimulationEvent::BaseDeleted {
                at: 3,
                removed_base: b'C',
            },
            SimulationEvent::IndelInserted {
                at: 4,
                base: b'T',
                segment: Segment::Np1,
                flags: NucFlags::empty(),
            },
            SimulationEvent::IndelDeleted {
                at: 5,
                removed_base: b'G',
                segment: Segment::V,
            },
            SimulationEvent::AssignmentChanged {
                segment: Segment::V,
                old: None,
                new: AlleleInstance::new(crate::refdata::AlleleId::new(2)),
            },
            SimulationEvent::TrimChanged {
                segment: Segment::J,
                end: TrimEnd::Five,
                old: Some(0),
                new: 3,
            },
            SimulationEvent::RegionAdded {
                region: Region::new(Segment::V, NucHandle::new(0), NucHandle::new(10)),
            },
            SimulationEvent::RegionReplaced {
                old: Region::new(Segment::V, NucHandle::new(0), NucHandle::new(10)),
                new: Region::new(Segment::V, NucHandle::new(0), NucHandle::new(12)),
            },
            SimulationEvent::ReverseComplementFlagRecorded { applied: true },
            SimulationEvent::MutationCountChanged {
                old: 2,
                new: 5,
                delta: 3,
            },
        ];
        for case in &cases {
            let cloned = case.clone();
            assert_eq!(&cloned, case);
        }
    }

    /// Tiny vec-based sink used as a fixture in the parity tests
    /// below — proves the trait can be implemented externally.
    #[derive(Default)]
    struct VecSink(Vec<SimulationEvent>);

    impl SimulationEventSink for VecSink {
        fn on_event(&mut self, event: &SimulationEvent) {
            self.0.push(event.clone());
        }
    }

    #[test]
    fn sink_trait_can_be_implemented_externally() {
        let mut sink = VecSink::default();
        sink.on_event(&SimulationEvent::ReverseComplementFlagRecorded { applied: true });
        sink.on_event(&SimulationEvent::BaseChanged {
            handle: NucHandle::new(0),
            old_base: b'A',
            new_base: b'C',
            segment: Segment::V,
            germline_pos: None,
        });
        assert_eq!(sink.0.len(), 2);
        assert!(matches!(
            sink.0[0],
            SimulationEvent::ReverseComplementFlagRecorded { applied: true }
        ));
        assert!(matches!(
            sink.0[1],
            SimulationEvent::BaseChanged {
                old_base: b'A',
                new_base: b'C',
                ..
            }
        ));
    }
}
