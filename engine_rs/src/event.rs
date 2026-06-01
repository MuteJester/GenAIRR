//! Committed simulation event ledger.
//!
//! Events are not a pub/sub bus and they are not emitted directly by
//! individual passes. The compiled simulator creates one event only
//! after a pass has completed, strict fences have passed, and the
//! trace delta/state revision are ready to commit atomically.
//!
//! # Compile effects vs. runtime events — the policy
//!
//! [`EventRecord`] carries both
//! [`compile_effects`](EventRecord::compile_effects) (static
//! declarations) and
//! [`simulation_events`](EventRecord::simulation_events) (runtime
//! consequences). The architecture deliberately allows these to
//! diverge: a test pass can declare an effect it doesn't materialise,
//! and an unconventional pass can emit events it didn't declare. The
//! runtime derived-state refresh (e.g.
//! [`crate::live_call::LiveCallRefreshHook`]) trusts the
//! `simulation_events` stream — not the declarations — so divergence
//! is a *legitimate* failure mode the architecture has to support.
//!
//! That flexibility has a flip side: a new built-in mutating pass
//! could declare a compile effect, mutate `Simulation` directly via
//! `sim.with_*` (bypassing the builder), and silently produce zero
//! events. Derived-state refresh would then skip the pass.
//!
//! To prevent that anti-pattern from creeping back, the test-only
//! [`check_event_emission_consistency`] helper applies a policy
//! mapping over an `EventRecord`:
//!
//! | `PassCompileEffect`     | Required event variant(s)                        |
//! |-------------------------|--------------------------------------------------|
//! | `AssignAllele(seg)`     | `AssignmentChanged { segment: seg, .. }`         |
//! | `TrimAllele(seg)`       | `TrimChanged { segment: seg, .. }`               |
//! | `AssembleSegment(seg)`  | `BasePushed` ∧ `RegionAdded { region.segment }`  |
//! | `AppendRegion(seg)`     | `BasePushed` ∧ `RegionAdded { region.segment }`  |
//! | `AppendNucleotides`     | `BasePushed`                                     |
//! | `EditBases`             | `BaseChanged`, **or** zero events (legitimate    |
//! |                         | no-op when the pass's stochastic count was 0)    |
//! | `StructuralIndel`       | `IndelInserted` ∨ `IndelDeleted`, **or** zero    |
//! |                         | events (legitimate no-op for zero-indel counts)  |
//!
//! Built-in full-stack plans are expected to satisfy this policy.
//! The check is **not** wired into the executor — it's a test-suite
//! invariant, run by the policy-conformance integration test in
//! `passes::tests`. Out-of-tree custom passes are free to diverge;
//! the policy isn't enforced globally because compile effects are
//! intentionally a declarative surface, not a contract on emission.

use crate::ir::{Simulation, SimulationEvent};
use crate::pass::PassCompileEffect;

/// Stable class of a committed event.
#[derive(Copy, Clone, Debug, Eq, PartialEq)]
pub enum EventKind {
    /// One pass committed a state transition.
    PassCommitted,
}

/// Half-open trace range covered by one committed event.
#[derive(Copy, Clone, Debug, Eq, PartialEq)]
pub struct TraceSpan {
    pub start: usize,
    pub end: usize,
}

impl TraceSpan {
    pub fn new(start: usize, end: usize) -> Self {
        assert!(
            start <= end,
            "TraceSpan start must be <= end (got {}..{})",
            start,
            end
        );
        Self { start, end }
    }

    pub fn len(&self) -> usize {
        self.end - self.start
    }

    pub fn is_empty(&self) -> bool {
        self.start == self.end
    }
}

/// Lightweight state summary recorded before and after a committed pass.
#[derive(Copy, Clone, Debug, Eq, PartialEq)]
pub struct StateSummary {
    pub pool_len: usize,
    pub region_count: usize,
    pub assigned_allele_count: usize,
}

impl StateSummary {
    pub fn from_simulation(sim: &Simulation) -> Self {
        let assigned_allele_count = sim.assignments.iter().count();

        Self {
            pool_len: sim.pool.len(),
            region_count: sim.sequence.region_count(),
            assigned_allele_count,
        }
    }
}

/// One durable event in the committed run ledger.
///
/// Three complementary surfaces per committed pass — each captures
/// a different layer of "what happened":
///
/// - [`compile_effects`](Self::compile_effects) — the **static
///   declarations** the pass made via
///   [`crate::pass::Pass::effects`]. Used by the schedule analyzer
///   for ordering. Read at compile time, not by runtime derived-
///   state refresh. See [`PassCompileEffect`] for the
///   static-vs-runtime split.
/// - [`trace_span`](Self::trace_span) — the slice of the run trace
///   recording the *choices* this pass consumed or emitted.
/// - [`simulation_events`](Self::simulation_events) — the ordered
///   list of [`SimulationEvent`]s the pass *actually emitted*
///   during commit via its internal builder(s). This is the
///   source-of-truth runtime consequence stream consumed by
///   [`crate::live_call::LiveCallRefreshHook`] and any future
///   derived-state hook.
///
/// A pass with non-empty `compile_effects` can legitimately have
/// an empty `simulation_events` (declared an edit but no-oped) or
/// vice versa — see the divergence tests in
/// `compiled::tests::live_call_edits`.
#[derive(Clone, Debug, Eq, PartialEq)]
pub struct EventRecord {
    pub pass_index: usize,
    pub pass_name: String,
    pub kind: EventKind,
    /// Static, compile-time effect declarations from
    /// [`crate::pass::Pass::effects`]. Used by the schedule
    /// analyzer; ignored by runtime derived-state refresh.
    pub compile_effects: Vec<PassCompileEffect>,
    pub trace_span: TraceSpan,
    pub pre: StateSummary,
    pub post: StateSummary,
    /// Typed state-change events the pass emitted during commit, in
    /// emission order. See [`SimulationEvent`] for the variants.
    /// This is the **runtime consequences** surface — the source of
    /// truth for derived-state refresh.
    pub simulation_events: Vec<SimulationEvent>,
}

/// Failure mode for [`check_event_emission_consistency`]. Carries
/// the offending `compile_effect` and a human-readable reason so
/// the test-suite reporter can pinpoint the divergent pass without
/// dumping the full event vector.
#[cfg(test)]
#[derive(Clone, Debug, Eq, PartialEq)]
pub(crate) struct EventEmissionPolicyError {
    pub pass_name: String,
    pub compile_effect: PassCompileEffect,
    pub reason: String,
}

/// Apply the **declared-effect / emitted-event consistency policy**
/// to a single committed pass record.
///
/// See the module-level doc for the policy table. Returns `Ok(())`
/// when every declared `compile_effect` is matched by an
/// appropriate event variant (or qualifies under the documented
/// `EditBases` / `StructuralIndel` zero-exemptions).
///
/// **Test-only.** Production code paths never call this; the
/// runtime derived-state refresh trusts events directly. This
/// helper exists so the test suite can fail a CI build when a
/// new built-in pass mutates `Simulation` directly and forgets to
/// emit corresponding events.
#[cfg(test)]
pub(crate) fn check_event_emission_consistency(
    record: &EventRecord,
) -> Result<(), EventEmissionPolicyError> {
    use crate::ir::{Segment, SimulationEvent};

    fn has_event<F>(events: &[SimulationEvent], pat: F) -> bool
    where
        F: Fn(&SimulationEvent) -> bool,
    {
        events.iter().any(pat)
    }

    let fail = |effect: PassCompileEffect, reason: &str| EventEmissionPolicyError {
        pass_name: record.pass_name.clone(),
        compile_effect: effect,
        reason: reason.to_string(),
    };

    for effect in &record.compile_effects {
        match effect {
            PassCompileEffect::AssignAllele(seg) => {
                if !has_event(&record.simulation_events, |e| {
                    matches!(e, SimulationEvent::AssignmentChanged { segment, .. } if segment == seg)
                }) {
                    return Err(fail(
                        *effect,
                        "declared AssignAllele but emitted no AssignmentChanged for that segment",
                    ));
                }
            }
            PassCompileEffect::TrimAllele(seg) => {
                if !has_event(&record.simulation_events, |e| {
                    matches!(e, SimulationEvent::TrimChanged { segment, .. } if segment == seg)
                }) {
                    return Err(fail(
                        *effect,
                        "declared TrimAllele but emitted no TrimChanged for that segment",
                    ));
                }
            }
            PassCompileEffect::AssembleSegment(seg)
            | PassCompileEffect::AppendRegion(seg) => {
                let has_push = has_event(&record.simulation_events, |e| {
                    matches!(e, SimulationEvent::BasePushed { segment, .. } if segment == seg)
                });
                let has_region = has_event(&record.simulation_events, |e| {
                    matches!(e, SimulationEvent::RegionAdded { region } if region.segment == *seg)
                });
                if !has_push {
                    return Err(fail(
                        *effect,
                        "declared region producer but emitted no BasePushed for its segment",
                    ));
                }
                if !has_region {
                    return Err(fail(
                        *effect,
                        "declared region producer but emitted no RegionAdded for its segment",
                    ));
                }
            }
            PassCompileEffect::AppendNucleotides => {
                if !has_event(&record.simulation_events, |e| {
                    matches!(e, SimulationEvent::BasePushed { .. })
                }) {
                    return Err(fail(
                        *effect,
                        "declared AppendNucleotides but emitted no BasePushed",
                    ));
                }
            }
            PassCompileEffect::EditBases => {
                // Zero-exemption: a stochastic mutation/PCR/quality
                // pass may legitimately have rolled count=0 and
                // emitted nothing. Only fail when the pass emitted
                // *other* mutating events but skipped `BaseChanged` —
                // that pattern means a writer path silently bypassed
                // the builder.
                let has_base_changed = has_event(&record.simulation_events, |e| {
                    matches!(e, SimulationEvent::BaseChanged { .. })
                });
                let has_any_mutation_event = record
                    .simulation_events
                    .iter()
                    .any(|e| !matches!(e, SimulationEvent::MutationCountChanged { .. }));
                if has_any_mutation_event && !has_base_changed {
                    return Err(fail(
                        *effect,
                        "declared EditBases and emitted mutating events but no BaseChanged \
                         — direct sim.with_* path is bypassing the builder?",
                    ));
                }
                // `(no_other_mutating_events && no_base_changed)`
                // is the legitimate zero-edit case. Allowed.
                //
                // Reference `seg = Segment::V` to keep clippy quiet
                // about the unused import path in cases where this
                // arm is the only one to fire.
                let _ = Segment::V;
            }
            PassCompileEffect::StructuralIndel => {
                // Zero-exemption: stochastic indel pass with
                // count=0 emits nothing. Mirror the EditBases
                // logic: only fail when *other* mutating events
                // exist but indel events don't.
                let has_indel = has_event(&record.simulation_events, |e| {
                    matches!(
                        e,
                        SimulationEvent::IndelInserted { .. }
                            | SimulationEvent::IndelDeleted { .. }
                    )
                });
                let has_any_mutation_event = record
                    .simulation_events
                    .iter()
                    .any(|e| !matches!(e, SimulationEvent::MutationCountChanged { .. }));
                if has_any_mutation_event && !has_indel {
                    return Err(fail(
                        *effect,
                        "declared StructuralIndel and emitted mutating events but no \
                         IndelInserted/IndelDeleted — direct sim.with_indel_* bypass?",
                    ));
                }
            }
        }
    }
    Ok(())
}

impl EventRecord {
    pub fn pass_committed(
        pass_index: usize,
        pass_name: impl Into<String>,
        compile_effects: Vec<PassCompileEffect>,
        trace_span: TraceSpan,
        pre: StateSummary,
        post: StateSummary,
        simulation_events: Vec<SimulationEvent>,
    ) -> Self {
        Self {
            pass_index,
            pass_name: pass_name.into(),
            kind: EventKind::PassCommitted,
            compile_effects,
            trace_span,
            pre,
            post,
            simulation_events,
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::ir::{NucFlags, Nucleotide, Segment};

    #[test]
    fn trace_span_len_uses_half_open_interval() {
        let span = TraceSpan::new(3, 8);
        assert_eq!(span.len(), 5);
        assert!(!span.is_empty());

        let empty = TraceSpan::new(4, 4);
        assert_eq!(empty.len(), 0);
        assert!(empty.is_empty());
    }

    #[test]
    #[should_panic(expected = "TraceSpan start must be <= end")]
    fn trace_span_rejects_inverted_range() {
        let _ = TraceSpan::new(9, 1);
    }

    #[test]
    fn state_summary_reflects_lightweight_simulation_shape() {
        let sim = Simulation::new();
        assert_eq!(
            StateSummary::from_simulation(&sim),
            StateSummary {
                pool_len: 0,
                region_count: 0,
                assigned_allele_count: 0,
            }
        );

        let (sim, _) = sim.with_nucleotide_pushed(Nucleotide::synthetic(
            b'A',
            Segment::Np1,
            NucFlags::empty(),
        ));
        assert_eq!(StateSummary::from_simulation(&sim).pool_len, 1);
    }

    #[test]
    fn event_record_has_distinct_static_and_runtime_surfaces() {
        // Naming-discipline test: a future contributor scanning
        // `EventRecord`'s fields must see at a glance that
        // `compile_effects` and `simulation_events` are two
        // different surfaces. This test pins the names + their
        // independence: an `EventRecord` can be built with one
        // populated and the other empty, in either direction.
        use crate::ir::{NucHandle, SimulationEvent};
        use crate::pass::PassCompileEffect;

        // Static-only — declared compile-effects, no runtime events.
        let static_only = EventRecord::pass_committed(
            0,
            "declared.only",
            vec![PassCompileEffect::EditBases],
            TraceSpan::new(0, 0),
            StateSummary::from_simulation(&Simulation::new()),
            StateSummary::from_simulation(&Simulation::new()),
            Vec::new(),
        );
        assert_eq!(static_only.compile_effects.len(), 1);
        assert!(static_only.simulation_events.is_empty());

        // Runtime-only — no compile-effect declaration, but a real
        // BaseChanged event was emitted.
        let runtime_only = EventRecord::pass_committed(
            1,
            "emitted.only",
            Vec::new(),
            TraceSpan::new(0, 0),
            StateSummary::from_simulation(&Simulation::new()),
            StateSummary::from_simulation(&Simulation::new()),
            vec![SimulationEvent::BaseChanged {
                handle: NucHandle::new(0),
                old_base: b'A',
                new_base: b'G',
                segment: Segment::V,
                germline_pos: None,
            }],
        );
        assert!(runtime_only.compile_effects.is_empty());
        assert_eq!(runtime_only.simulation_events.len(), 1);

        // Type-level invariant: the legacy `PassEffect` alias
        // still resolves to `PassCompileEffect`. Removes the alias
        // safely-removable-in-the-future signal: a tool can grep
        // for `PassEffect` to find pre-rename code paths.
        let _: crate::pass::PassEffect = PassCompileEffect::EditBases;
    }

    // ──────────────────────────────────────────────────────────
    // Unit tests for the event-emission policy. Construct synthetic
    // `EventRecord`s that exercise every match arm of
    // `check_event_emission_consistency`, including the
    // zero-exemption corners.
    // ──────────────────────────────────────────────────────────

    fn record(
        name: &str,
        compile_effects: Vec<PassCompileEffect>,
        events: Vec<SimulationEvent>,
    ) -> EventRecord {
        let empty = StateSummary::from_simulation(&Simulation::new());
        EventRecord::pass_committed(
            0,
            name,
            compile_effects,
            TraceSpan::new(0, 0),
            empty,
            empty,
            events,
        )
    }

    fn assignment_changed_v() -> SimulationEvent {
        use crate::assignment::AlleleInstance;
        use crate::refdata::AlleleId;
        SimulationEvent::AssignmentChanged {
            segment: Segment::V,
            old: None,
            new: AlleleInstance::new(AlleleId::new(0)),
        }
    }
    fn base_pushed_v() -> SimulationEvent {
        use crate::ir::NucHandle;
        SimulationEvent::BasePushed {
            handle: NucHandle::new(0),
            base: b'A',
            segment: Segment::V,
            germline_pos: None,
            flags: NucFlags::empty(),
        }
    }
    fn region_added_v() -> SimulationEvent {
        use crate::ir::{NucHandle, Region};
        SimulationEvent::RegionAdded {
            region: Region::new(Segment::V, NucHandle::new(0), NucHandle::new(1)),
        }
    }
    fn base_changed_v() -> SimulationEvent {
        use crate::ir::NucHandle;
        SimulationEvent::BaseChanged {
            handle: NucHandle::new(0),
            old_base: b'A',
            new_base: b'G',
            segment: Segment::V,
            germline_pos: None,
        }
    }
    fn indel_inserted() -> SimulationEvent {
        SimulationEvent::IndelInserted {
            at: 0,
            base: b'N',
            segment: Segment::V,
            flags: NucFlags::empty(),
        }
    }

    #[test]
    fn policy_passes_when_assign_allele_declaration_matches_emitted_event() {
        let r = record(
            "happy.assign",
            vec![PassCompileEffect::AssignAllele(Segment::V)],
            vec![assignment_changed_v()],
        );
        assert!(check_event_emission_consistency(&r).is_ok());
    }

    #[test]
    fn policy_fails_when_assign_allele_declared_but_no_event_emitted() {
        let r = record(
            "buggy.assign",
            vec![PassCompileEffect::AssignAllele(Segment::V)],
            vec![],
        );
        let err = check_event_emission_consistency(&r).unwrap_err();
        assert_eq!(
            err.compile_effect,
            PassCompileEffect::AssignAllele(Segment::V)
        );
        assert!(err.reason.contains("AssignmentChanged"));
    }

    #[test]
    fn policy_requires_both_base_pushed_and_region_added_for_assemble_segment() {
        // Only RegionAdded, no BasePushed → fail.
        let r = record(
            "buggy.assemble",
            vec![PassCompileEffect::AssembleSegment(Segment::V)],
            vec![region_added_v()],
        );
        let err = check_event_emission_consistency(&r).unwrap_err();
        assert!(err.reason.contains("BasePushed"));

        // Only BasePushed, no RegionAdded → fail too.
        let r = record(
            "buggy.assemble2",
            vec![PassCompileEffect::AssembleSegment(Segment::V)],
            vec![base_pushed_v()],
        );
        let err = check_event_emission_consistency(&r).unwrap_err();
        assert!(err.reason.contains("RegionAdded"));

        // Both present → ok.
        let r = record(
            "happy.assemble",
            vec![PassCompileEffect::AssembleSegment(Segment::V)],
            vec![base_pushed_v(), region_added_v()],
        );
        assert!(check_event_emission_consistency(&r).is_ok());
    }

    #[test]
    fn policy_allows_edit_bases_with_zero_events_legitimate_zero_count() {
        // EditBases declared but pass rolled count=0 and emitted
        // nothing. Should pass (zero-exemption).
        let r = record(
            "stochastic.editbases_zero",
            vec![PassCompileEffect::EditBases],
            vec![],
        );
        assert!(check_event_emission_consistency(&r).is_ok());
    }

    #[test]
    fn policy_allows_edit_bases_with_only_mutation_count_event() {
        // A future mutation pass might fire `MutationCountChanged`
        // even with zero actual edits (commit path bumps by delta=0?
        // Today guard prevents that, but the policy is conservative:
        // MutationCountChanged-only is treated as "no real mutation
        // events" and the zero-exemption still applies.
        let r = record(
            "stochastic.editbases_count_only",
            vec![PassCompileEffect::EditBases],
            vec![SimulationEvent::MutationCountChanged {
                old: 0,
                new: 0,
                delta: 0,
            }],
        );
        assert!(check_event_emission_consistency(&r).is_ok());
    }

    #[test]
    fn policy_fails_when_edit_bases_declared_and_mutating_events_present_but_no_base_changed() {
        // Pathological pass: declared EditBases, emitted an indel
        // event but no BaseChanged. Either the declaration is
        // wrong or the writer bypassed `change_base` — caller
        // intent is unclear, fail the policy.
        let r = record(
            "buggy.editbases_indel_only",
            vec![PassCompileEffect::EditBases],
            vec![indel_inserted()],
        );
        let err = check_event_emission_consistency(&r).unwrap_err();
        assert!(err.reason.contains("BaseChanged"));
    }

    #[test]
    fn policy_allows_structural_indel_with_zero_events() {
        let r = record(
            "stochastic.indel_zero",
            vec![PassCompileEffect::StructuralIndel],
            vec![],
        );
        assert!(check_event_emission_consistency(&r).is_ok());
    }

    #[test]
    fn policy_fails_when_structural_indel_declared_with_mutating_events_but_no_indel() {
        let r = record(
            "buggy.indel_only_base_change",
            vec![PassCompileEffect::StructuralIndel],
            vec![base_changed_v()],
        );
        let err = check_event_emission_consistency(&r).unwrap_err();
        assert!(err.reason.contains("IndelInserted"));
    }

    #[test]
    fn policy_passes_for_segment_replaced_with_no_declared_effects() {
        // Receptor-revision Slice B: `SegmentReplaced` events ride
        // the runtime-consequence channel without yet having a
        // matching `PassCompileEffect`. The declared/emitted policy
        // is keyed off `compile_effects`, so an empty compile list
        // means there's nothing to validate against — the runtime
        // event is accepted regardless. Pins this for Slice C: the
        // future `ReceptorRevisionPass` is free to declare no
        // compile-effect, or to declare a new one whose policy arm
        // we'd add then.
        use crate::ir::{NucHandle, Region};
        let r = record(
            "receptor_revision.no_declared_effects",
            vec![],
            vec![SimulationEvent::SegmentReplaced {
                segment: Segment::V,
                old_region: Region::new(Segment::V, NucHandle::new(0), NucHandle::new(5)),
                new_region: Region::new(Segment::V, NucHandle::new(0), NucHandle::new(8)),
                bytes_delta: 3,
            }],
        );
        assert!(check_event_emission_consistency(&r).is_ok());
    }

    #[test]
    fn policy_passes_for_no_declared_effects_regardless_of_events() {
        // RevCompPass declares no effects. Whatever it emits (or
        // doesn't), the policy has nothing to check against.
        let r = record(
            "no_effects",
            vec![],
            vec![SimulationEvent::ReverseComplementFlagRecorded { applied: true }],
        );
        assert!(check_event_emission_consistency(&r).is_ok());
    }
}
