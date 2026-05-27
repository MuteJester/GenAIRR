//! Shared test fixtures **and** test-driver helpers used across
//! multiple pass test modules.
//!
//! Compiled only under `#[cfg(test)]`. Two responsibilities:
//!
//! 1. **Fixtures** — stop-codon-filtering refdata / sim builders
//!    reusable across substitution-style pass tests.
//! 2. **Driver helpers** — `run_pass_capturing_events`,
//!    `run_pass_with_replay_records`, plus their assertion
//!    counterparts. Subsume the per-pass `run_replay` / `run_substitute`
//!    / `run_np_replay` boilerplate (~15-25 lines of `PassContext { ... }`
//!    construction each) into named helpers so new pass tests start
//!    from a one-line invocation.
//!
//! See [`crate::passes`] module docs and `docs/adding_a_pass.md`
//! for the recommended test patterns.

#![cfg(test)]
#![allow(dead_code)]

use crate::assignment::AlleleInstance;
use crate::contract::ContractSet;
use crate::dist::Distribution;
use crate::ir::{NucHandle, Nucleotide, Region, Segment, Simulation, SimulationEvent};
use crate::pass::{Outcome, Pass, PassContext, PassError};
use crate::refdata::{Allele, AlleleId, ChainType, RefDataConfig};
use crate::replay::TraceCursor;
use crate::rng::Rng;
use crate::trace::{ChoiceRecord, Trace};

/// A base distribution that always samples `A` but reports a
/// 2-element support `{A, C}`. Use to demonstrate that an active
/// productive-contract filter can divert from the default `A` to
/// the safe alternative `C` when `A` would create a stop codon.
#[derive(Clone, Debug)]
pub(crate) struct StopThenSafeMutationBaseDist;

impl Distribution for StopThenSafeMutationBaseDist {
    type Output = u8;

    fn sample(&self, _rng: &mut crate::rng::Rng) -> u8 {
        b'A'
    }

    fn support(&self) -> Option<Vec<(u8, f64)>> {
        Some(vec![(b'A', 1.0), (b'C', 1.0)])
    }
}

/// A degenerate distribution whose support is `{A}` only. When a
/// contract rejects `A` at the filter stage, strict mode must
/// surface `EmptyAdmissibleSupport`. Use to exercise the strict-mode
/// failure path.
#[derive(Clone, Debug)]
pub(crate) struct StopOnlyMutationBaseDist;

impl Distribution for StopOnlyMutationBaseDist {
    type Output = u8;

    fn sample(&self, _rng: &mut crate::rng::Rng) -> u8 {
        b'A'
    }

    fn support(&self) -> Option<Vec<(u8, f64)>> {
        Some(vec![(b'A', 1.0)])
    }
}

/// VJ fixture used by the substitution-style pass tests
/// (uniform-mutation, PCR error, contaminant). The pool already
/// holds `TAC` (V) + `TGG` (J) so position 2's neighbours are
/// `T A · T G G`. With anchors at V[0] and J[0], the junction is
/// `[0, 6)` and codon 0 (`TAC`) becomes `TAA` if position 2 is
/// substituted to `A` — the perfect stop-codon filter probe.
pub(crate) fn make_substitution_productive_vj_fixture() -> (RefDataConfig, Simulation) {
    let mut cfg = RefDataConfig::empty(ChainType::Vj);
    let _ = cfg.v_pool.push(Allele {
        name: "v_mut*01".into(),
        gene: "v_mut".into(),
        seq: b"TAC".to_vec(),
        segment: Segment::V,
        anchor: Some(0),
    });
    let _ = cfg.j_pool.push(Allele {
        name: "j_mut*01".into(),
        gene: "j_mut".into(),
        seq: b"TGG".to_vec(),
        segment: Segment::J,
        anchor: Some(0),
    });

    let mut sim = Simulation::new();
    for (i, &b) in b"TAC".iter().enumerate() {
        let (next, _) = sim.with_nucleotide_pushed(Nucleotide::germline(b, i as u16, Segment::V));
        sim = next;
    }
    let v_region = Region::new(Segment::V, NucHandle::new(0), NucHandle::new(3));
    sim = sim.with_region_added(v_region);

    for (i, &b) in b"TGG".iter().enumerate() {
        let (next, _) = sim.with_nucleotide_pushed(Nucleotide::germline(b, i as u16, Segment::J));
        sim = next;
    }
    let j_region = Region::new(Segment::J, NucHandle::new(3), NucHandle::new(6));
    sim = sim.with_region_added(j_region);

    sim = sim
        .with_allele_assigned(Segment::V, AlleleInstance::new(AlleleId::new(0)))
        .with_allele_assigned(Segment::J, AlleleInstance::new(AlleleId::new(0)));

    (cfg, sim)
}

/// VJ fixture where every junction position has a single-base anchor
/// mask (V = TGG / W, J = TGG / W). Under a base distribution whose
/// support is disjoint from `{T, G}` (for example `{A}`), every site
/// admits zero mass. Substitution-style passes use this to exercise
/// strict-mode `EmptyAdmissibleSupport` without each test module
/// carrying its own copy of the setup.
pub(crate) fn make_fully_locked_vj_fixture() -> (RefDataConfig, Simulation) {
    let mut cfg = RefDataConfig::empty(ChainType::Vj);
    let _ = cfg.v_pool.push(Allele {
        name: "v_lock*01".into(),
        gene: "v_lock".into(),
        seq: b"TGG".to_vec(),
        segment: Segment::V,
        anchor: Some(0),
    });
    let _ = cfg.j_pool.push(Allele {
        name: "j_lock*01".into(),
        gene: "j_lock".into(),
        seq: b"TGG".to_vec(),
        segment: Segment::J,
        anchor: Some(0),
    });

    let mut sim = Simulation::new();
    for (i, &b) in b"TGG".iter().enumerate() {
        let (next, _) = sim.with_nucleotide_pushed(Nucleotide::germline(b, i as u16, Segment::V));
        sim = next;
    }
    sim = sim.with_region_added(Region::new(
        Segment::V,
        NucHandle::new(0),
        NucHandle::new(3),
    ));

    for (i, &b) in b"TGG".iter().enumerate() {
        let (next, _) = sim.with_nucleotide_pushed(Nucleotide::germline(b, i as u16, Segment::J));
        sim = next;
    }
    sim = sim.with_region_added(Region::new(
        Segment::J,
        NucHandle::new(3),
        NucHandle::new(6),
    ));

    sim = sim
        .with_allele_assigned(Segment::V, AlleleInstance::new(AlleleId::new(0)))
        .with_allele_assigned(Segment::J, AlleleInstance::new(AlleleId::new(0)));

    (cfg, sim)
}

// ──────────────────────────────────────────────────────────────────
// Pass-driver helpers
//
// Subsume the recurring per-pass `run_replay` / `run_substitute` /
// `run_np_replay` boilerplate. Each one constructs a fresh
// `PassContext` with documented defaults, threads the closure-owned
// trace / rng / event buffers, and returns the captured artifacts.
// Callers pick the helper that matches their pattern:
//
// | Pattern                                | Helper                              |
// |----------------------------------------|-------------------------------------|
// | Run a pass and inspect its events      | `run_pass_capturing_events`         |
// | Run a pass under trace-injected replay | `run_pass_with_replay_records`      |
// | Drive a full plan through round-trip   | `assert_compiled_simulator_replay_round_trip` |
//
// The corresponding `assert_*` wrappers compose on top of these for
// the most common pin assertions.
// ──────────────────────────────────────────────────────────────────

/// Seed used by every test driver below. A magic number keeps test
/// failures across different pass tests reproducible: the same seed
/// is asserted at the same call site no matter which pass is being
/// driven.
pub(crate) const TEST_DRIVER_SEED: u64 = 0xc0ff_ee;

/// Captured artifacts after running a pass through one of the test
/// driver helpers below. The pass author asserts on whichever fields
/// the test cares about — the others are silently available for
/// diagnostic prints.
pub(crate) struct PassRunCapture {
    /// Post-pass simulation.
    pub sim: Simulation,
    /// The `simulation_events` stream the pass emitted via its
    /// internal builders. Empty when no builder was used or when
    /// the pass routes through a `Simulation::with_*` path
    /// (which would also trip the production-path lockdown).
    pub events: Vec<SimulationEvent>,
    /// The trace delta the pass wrote.
    pub trace: Trace,
}

/// Run `pass` through `execute_checked` with permissive defaults
/// and an event-log sink attached. The captured stream lets the
/// caller pin which `SimulationEvent` variants the pass actually
/// emitted — the load-bearing assertion for any new mutating pass.
///
/// Defaults:
/// - `pass_index = 0`, `feasibility = None`, `reference_index = None`
/// - No replay cursor (fresh-RNG path)
/// - RNG seeded with [`TEST_DRIVER_SEED`]
///
/// Returns `Err(_)` if the pass returns a structured `PassError`.
/// Callers that want to assert on the failure mode rather than the
/// success path can pattern-match on the result directly.
pub(crate) fn run_pass_capturing_events(
    pass: &dyn Pass,
    sim: &Simulation,
    refdata: Option<&RefDataConfig>,
    contracts: Option<&ContractSet>,
) -> Result<PassRunCapture, PassError> {
    let mut trace = Trace::new();
    let mut rng = Rng::new(TEST_DRIVER_SEED);
    let mut events: Vec<SimulationEvent> = Vec::new();
    let next = {
        let mut ctx = PassContext {
            trace: &mut trace,
            rng: &mut rng,
            pass_index: 0,
            refdata,
            contracts,
            feasibility: None,
            reference_index: None,
            replay_cursor: None,
            event_log_sink: Some(&mut events),
        };
        pass.execute_checked(sim, &mut ctx)?
    };
    Ok(PassRunCapture {
        sim: next,
        events,
        trace,
    })
}

/// Drive `pass` through trace-injected replay against `records`.
/// The cursor is drained as the pass consumes each recorded value
/// — passes that haven't migrated to the replay-injected path still
/// branch on `rng`, in which case the cursor stays full and the
/// recorded values flow through downstream sampling sites.
///
/// Strict policy: every replay-validation mismatch surfaces as
/// `Err(PassError::*)` so tests can pattern-match on the failure
/// mode without unwrap-panics burying the diagnostic.
pub(crate) fn run_pass_with_replay_records(
    pass: &dyn Pass,
    sim: &Simulation,
    records: Vec<ChoiceRecord>,
    refdata: Option<&RefDataConfig>,
    contracts: Option<&ContractSet>,
) -> (Result<Simulation, PassError>, Trace, Vec<SimulationEvent>) {
    let mut cursor = TraceCursor::from_owned(records);
    let mut trace = Trace::new();
    let mut rng = Rng::new(TEST_DRIVER_SEED);
    let mut events: Vec<SimulationEvent> = Vec::new();
    let result = {
        let mut ctx = PassContext {
            trace: &mut trace,
            rng: &mut rng,
            pass_index: 0,
            refdata,
            contracts,
            feasibility: None,
            reference_index: None,
            replay_cursor: Some(&mut cursor),
            event_log_sink: Some(&mut events),
        };
        pass.execute_checked(sim, &mut ctx)
    };
    (result, trace, events)
}

/// Assert at least one event in `events` matches `predicate`.
/// Failure messages dump the full captured stream so the test
/// reporter shows what *was* emitted alongside what was expected.
pub(crate) fn assert_event_matching<F>(
    events: &[SimulationEvent],
    predicate: F,
    description: &str,
) where
    F: Fn(&SimulationEvent) -> bool,
{
    assert!(
        events.iter().any(predicate),
        "Expected at least one event matching: {description}.\n\
         Captured events ({n}): {events:?}",
        n = events.len(),
    );
}

/// Assert two committed traces are record-for-record equal on
/// `(address, value)`. Diverges with a structured message naming
/// the first mismatching index plus its address — the canonical
/// diagnostic for replay-equality regressions.
pub(crate) fn assert_traces_record_equal(a: &[ChoiceRecord], b: &[ChoiceRecord], label: &str) {
    assert_eq!(
        a.len(),
        b.len(),
        "{label}: trace length mismatch ({} vs {})",
        a.len(),
        b.len()
    );
    for (i, (x, y)) in a.iter().zip(b.iter()).enumerate() {
        assert_eq!(
            x.address, y.address,
            "{label}: address divergence at record {i}: {:?} vs {:?}",
            x.address, y.address
        );
        assert_eq!(
            x.value, y.value,
            "{label}: value divergence at record {i} (addr={})",
            x.address
        );
    }
}

/// Drive an `OwnedCompiledSimulator` through one fresh run plus a
/// trace-injected replay against the same recorded trace, and
/// assert the two outcomes' traces match record-for-record. The
/// load-bearing replay-determinism check for plan-level tests.
///
/// Permissive policy is used for the fresh run; replay uses strict
/// (the standard contract — `replay_records` forces strict
/// dispatch internally regardless of policy).
pub(crate) fn assert_compiled_simulator_replay_round_trip(
    compiled: &crate::compiled::OwnedCompiledSimulator,
    seed: u64,
) {
    use crate::compiled::ExecutionPolicy;

    let original: Outcome = compiled
        .run_one(seed)
        .expect("fresh permissive run should succeed");
    let original_records: Vec<ChoiceRecord> = original.trace.choices().to_vec();
    let replayed: Outcome = compiled
        .replay_from_trace_records(&original_records, seed, ExecutionPolicy::Strict)
        .expect("replay against the same plan should succeed");
    assert_traces_record_equal(
        original.trace.choices(),
        replayed.trace.choices(),
        "compiled-simulator replay round-trip",
    );
}

// ──────────────────────────────────────────────────────────────────
// Sanity tests for the helpers themselves. Each one exercises a
// helper through a small built-in pass so a regression in the
// helper plumbing fails fast, separate from any one pass's tests.
// ──────────────────────────────────────────────────────────────────

#[cfg(test)]
mod helper_tests {
    use super::*;
    use crate::address::ChoiceAddress;
    use crate::compiled::{ExecutionPolicy, OwnedCompiledSimulator};
    use crate::dist::AllelePoolDist;
    use crate::pass::PassPlan;
    use crate::passes::sample_allele::test_support::make_test_pool;
    use crate::passes::{AssembleSegmentPass, SampleAllelePass};
    use crate::trace::ChoiceValue;

    fn vj_refdata_one_allele_each() -> RefDataConfig {
        let mut cfg = RefDataConfig::empty(ChainType::Vj);
        let _ = cfg.v_pool.push(Allele {
            name: "v*01".into(),
            gene: "v".into(),
            seq: b"AAACCCGGG".to_vec(),
            segment: Segment::V,
            anchor: Some(6),
        });
        let _ = cfg.j_pool.push(Allele {
            name: "j*01".into(),
            gene: "j".into(),
            seq: b"TTTAAA".to_vec(),
            segment: Segment::J,
            anchor: Some(0),
        });
        cfg
    }

    #[test]
    fn run_pass_capturing_events_returns_emitted_assignment_changed() {
        // SampleAllelePass with a uniform single-allele distribution
        // emits exactly one `AssignmentChanged`. Use the helper to
        // capture and inspect.
        let pool = make_test_pool(1, Segment::V);
        let pass =
            SampleAllelePass::new(Segment::V, Box::new(AllelePoolDist::uniform(&pool)));
        let cfg = vj_refdata_one_allele_each();

        let capture =
            run_pass_capturing_events(&pass, &Simulation::new(), Some(&cfg), None).unwrap();
        assert_event_matching(
            &capture.events,
            |e| matches!(e, SimulationEvent::AssignmentChanged { segment, .. } if *segment == Segment::V),
            "AssignmentChanged for V",
        );
        // Trace also captures the recorded choice.
        assert_eq!(capture.trace.len(), 1);
    }

    #[test]
    fn run_pass_with_replay_records_consumes_the_cursor() {
        // Replay a recorded AlleleId=0 through SampleAllelePass.
        // The strict replay path validates the recorded value
        // against refdata + support; success returns the new sim.
        let pool = make_test_pool(1, Segment::V);
        let pass =
            SampleAllelePass::new(Segment::V, Box::new(AllelePoolDist::uniform(&pool)));
        let cfg = vj_refdata_one_allele_each();
        let records = vec![ChoiceRecord::new(
            ChoiceAddress::SampleAllele(crate::address::VdjSegment::V).to_string(),
            ChoiceValue::AlleleId(0),
        )];
        let (result, trace, _events) =
            run_pass_with_replay_records(&pass, &Simulation::new(), records, Some(&cfg), None);
        result.expect("valid replay record");
        // Trace re-emitted by the replay path.
        assert_eq!(trace.len(), 1);
    }

    #[test]
    fn assert_traces_record_equal_passes_on_identical_traces() {
        let a = vec![ChoiceRecord::new(
            "sample_allele.v".to_string(),
            ChoiceValue::AlleleId(0),
        )];
        let b = a.clone();
        assert_traces_record_equal(&a, &b, "identical");
    }

    #[test]
    #[should_panic(expected = "address divergence")]
    fn assert_traces_record_equal_fails_on_address_divergence() {
        let a = vec![ChoiceRecord::new(
            "sample_allele.v".to_string(),
            ChoiceValue::AlleleId(0),
        )];
        let b = vec![ChoiceRecord::new(
            "sample_allele.d".to_string(),
            ChoiceValue::AlleleId(0),
        )];
        assert_traces_record_equal(&a, &b, "divergent");
    }

    #[test]
    fn assert_event_matching_succeeds_when_predicate_finds_event() {
        let events = vec![SimulationEvent::ReverseComplementFlagRecorded { applied: true }];
        assert_event_matching(
            &events,
            |e| matches!(e, SimulationEvent::ReverseComplementFlagRecorded { applied: true }),
            "rev-comp applied",
        );
    }

    #[test]
    #[should_panic(expected = "Expected at least one event matching")]
    fn assert_event_matching_fails_when_predicate_finds_nothing() {
        let events = vec![SimulationEvent::ReverseComplementFlagRecorded { applied: false }];
        assert_event_matching(
            &events,
            |e| matches!(e, SimulationEvent::ReverseComplementFlagRecorded { applied: true }),
            "rev-comp applied (true)",
        );
    }

    #[test]
    fn assert_compiled_replay_round_trip_holds_on_minimal_recombine_plan() {
        let cfg = vj_refdata_one_allele_each();
        let mut plan = PassPlan::new();
        plan.push(Box::new(SampleAllelePass::new(
            Segment::V,
            Box::new(AllelePoolDist::uniform(&cfg.v_pool)),
        )));
        plan.push(Box::new(SampleAllelePass::new(
            Segment::J,
            Box::new(AllelePoolDist::uniform(&cfg.j_pool)),
        )));
        plan.push(Box::new(AssembleSegmentPass::new(Segment::V)));
        plan.push(Box::new(AssembleSegmentPass::new(Segment::J)));
        let compiled = OwnedCompiledSimulator::compile(
            plan,
            Some(cfg),
            None,
            ExecutionPolicy::Permissive,
        )
        .expect("plan compiles");
        assert_compiled_simulator_replay_round_trip(&compiled, 1234);
    }
}
