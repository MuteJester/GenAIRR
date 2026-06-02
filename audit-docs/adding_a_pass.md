# Adding a new pass

A practical, copy-pasteable companion to
[`engine_architecture.md`](engine_architecture.md). Read this when
you want to add a new biological mechanism to the simulation
pipeline. Read the architecture doc when you want to understand
*why* the rules below exist.

> **Time budget.** A pass that follows this template typically
> lands in ~1-2 hours including tests. The CI guardrails will tell
> you immediately if you skipped a step.

---

## TL;DR — the contributor's checklist

For each new mutating pass:

- [ ] **Declare it.** Implement
  [`Pass`](../engine_rs/src/pass/traits.rs) with `name()`,
  `effects()`, `execute()`, `execute_checked()`. Add typed
  `declared_choice_patterns()` if the pass samples anything.
- [ ] **Route mutations through `SimulationBuilder`.** Use the
  event-emitting builder methods or
  [`MutationTransaction`](../engine_rs/src/passes/mutation_transaction/mod.rs).
  Never call `Simulation::with_*` from production pass code —
  the lockdown will fail your build.
- [ ] **Handle replay.** When `ctx.replay_cursor.is_some()`,
  consume the recorded value, run it through the same
  admissibility chain a fresh draw would, and only then apply.
- [ ] **Three test layers.** Event-emission test, replay-validation
  test (when applicable), policy-conformance entry. Use the
  helpers in [`passes/test_support.rs`](../engine_rs/src/passes/test_support.rs)
  to keep boilerplate low.
- [ ] **Biological coverage** (if the pass can break productivity).
  Add an entry to
  [`tests/test_productive_stress_matrix.py`](../tests/test_productive_stress_matrix.py).

---

## Minimal pass template

A skeleton you can paste into a new file under
`engine_rs/src/passes/`. Replace `MyMechanismPass` with the real
name; fill in the body where commented.

```rust
//! `MyMechanismPass` — short biological description of what this
//! pass does and which AIRR fields it influences.

use crate::address;
use crate::ir::{Simulation, SimulationBuilder};
use crate::pass::{Pass, PassContext, PassEffect, PassError};
use crate::trace::ChoiceValue;

pub struct MyMechanismPass {
    // ... pass parameters (distributions, configuration, ...)
}

impl MyMechanismPass {
    pub fn new(/* … */) -> Self {
        Self { /* … */ }
    }

    fn choice_address(&self) -> address::ChoiceAddress {
        // Add a variant to ChoiceAddress; address vocabulary is
        // schema-versioned (ADDRESS_SCHEMA_VERSION). Pin the new
        // spelling via the `frozen_address_spellings_for_choice_address_schema_v1`
        // test.
        address::ChoiceAddress::MyMechanism /* (optional segment) */
    }

    fn execute_with_sampling_mode(
        &self,
        sim: &Simulation,
        ctx: &mut PassContext,
        _strict: bool,
    ) -> Result<Simulation, PassError> {
        // 1. Replay branch: consume the recorded value and validate it
        //    through the SAME admissibility chain a fresh draw would.
        //    See SampleAllelePass::validate_replay_candidate for the
        //    canonical template.
        if ctx.replay_cursor.is_some() {
            // Consume → validate → record to ctx.trace → apply.
            todo!("replay branch");
        }

        // 2. Fresh-RNG branch: sample with the active contract /
        //    feasibility filter, record to ctx.trace, then apply via
        //    the event-emitting builder.
        //
        // For substitutions / indels:    use MutationTransaction.
        // For assignment / trim / region: use SimulationBuilder
        //                                 directly + event-log drain.
        let mut builder = SimulationBuilder::from_simulation(sim.clone());
        if ctx.event_log_sink.is_some() {
            builder.attach_event_log_observer();
        }

        // ... call builder.assign_allele / update_trim / add_region /
        // set_mutation_count / record_reverse_complement_flag ...

        if let Some(sink) = ctx.event_log_sink.as_deref_mut() {
            sink.extend(builder.seal_event_log_observer());
        }
        Ok(builder.seal())
    }
}

impl Pass for MyMechanismPass {
    fn name(&self) -> &str {
        "mechanism.my_mechanism"
    }

    fn execute(&self, sim: &Simulation, ctx: &mut PassContext) -> Simulation {
        self.execute_with_sampling_mode(sim, ctx, false)
            .expect("MyMechanismPass permissive execution must not return PassError")
    }

    fn execute_checked(
        &self,
        sim: &Simulation,
        ctx: &mut PassContext,
    ) -> Result<Simulation, PassError> {
        self.execute_with_sampling_mode(sim, ctx, true)
    }

    fn declared_choice_patterns(&self) -> Vec<address::ChoiceAddressPattern> {
        vec![address::ChoiceAddressPattern::MyMechanism]
    }

    fn effects(&self) -> Vec<PassEffect> {
        // Static, compile-time declaration. NOT consulted at
        // runtime — declare honestly what the SCHEDULE needs to
        // know, then route the mutation through the event-emitting
        // builder so the RUNTIME refresh hook sees real events.
        vec![PassEffect::EditBases /* or AssembleSegment(seg), etc. */]
    }
}
```

---

## Required test patterns

Three load-bearing tests per mutating pass. The helpers in
[`passes::test_support`](../engine_rs/src/passes/test_support.rs)
keep each one to ~10 lines.

### 1. Event emission

Prove the pass emits the consequence variants its compile effect
implies. The
[`builtin_passes_emit_events_consistent_with_declared_compile_effects`](../engine_rs/src/passes/mod.rs)
policy will fail your CI run otherwise, but a focused unit test
gives a faster signal.

```rust
use crate::passes::test_support::{
    assert_event_matching, run_pass_capturing_events,
};
use crate::ir::SimulationEvent;

#[test]
fn my_mechanism_pass_emits_base_changed() {
    let pass = MyMechanismPass::new(/* … */);
    let sim = make_test_fixture(); // your own initial sim
    let cfg = make_test_refdata(); // your own refdata

    let capture =
        run_pass_capturing_events(&pass, &sim, Some(&cfg), None).unwrap();
    assert_event_matching(
        &capture.events,
        |e| matches!(e, SimulationEvent::BaseChanged { .. }),
        "BaseChanged from MyMechanismPass",
    );
}
```

### 2. Replay-validation

If the pass samples anything, prove that an invalid recorded value
surfaces as a structured `PassError::*` — not silently applied.

```rust
use crate::passes::test_support::run_pass_with_replay_records;
use crate::pass::PassError;
use crate::trace::{ChoiceRecord, ChoiceValue};

#[test]
fn my_mechanism_pass_replay_rejects_out_of_support_value() {
    let pass = MyMechanismPass::new(/* … */);
    let cfg = make_test_refdata();
    let records = vec![ChoiceRecord::new(
        pass.choice_address().to_string(),
        ChoiceValue::Int(999_999), // outside the natural support
    )];

    let (result, _trace, _events) = run_pass_with_replay_records(
        &pass,
        &make_test_fixture(),
        records,
        Some(&cfg),
        None,
    );

    match result.unwrap_err() {
        PassError::InvalidDistributionOutput { reason, .. } => {
            assert!(reason.contains("not_in_distribution_support"));
        }
        other => panic!("expected InvalidDistributionOutput, got {other:?}"),
    }
}
```

The canonical template is
[`SampleAllelePass`'s replay-validation tests](../engine_rs/src/passes/sample_allele.rs)
— four parallel tests covering refdata-miss, support-miss,
contract-rejection, and the happy path.

### 3. Compile-effect / event consistency

Append your pass to
[`build_coverage_audit_plan`](../engine_rs/src/passes/mod.rs) so
the [event-emission policy
test](../engine_rs/src/passes/mod.rs) (`builtin_passes_emit_events_consistent_with_declared_compile_effects`)
exercises it. No new test file needed — extending the audit plan
is the canonical way to keep the policy comprehensive.

If your pass declares an effect with a stricter or new emission
contract, also add a unit case under
[`event::tests`](../engine_rs/src/event.rs) so the policy logic
itself stays comprehensive.

---

## Drive a full-plan replay round-trip

When your pass interacts with other passes (almost always), pin
that the full plan still survives trace replay end-to-end. The
helper:

```rust
use crate::passes::test_support::assert_compiled_simulator_replay_round_trip;
use crate::compiled::{ExecutionPolicy, OwnedCompiledSimulator};
use crate::pass::PassPlan;

#[test]
fn full_plan_with_my_mechanism_survives_replay() {
    let cfg = make_test_refdata();
    let mut plan = PassPlan::new();
    plan.push(Box::new(SampleAllelePass::new(/* … */)));
    plan.push(Box::new(AssembleSegmentPass::new(Segment::V)));
    plan.push(Box::new(MyMechanismPass::new(/* … */)));

    let compiled = OwnedCompiledSimulator::compile(
        plan,
        Some(cfg),
        None,
        ExecutionPolicy::Permissive,
    )
    .expect("plan compiles");

    assert_compiled_simulator_replay_round_trip(&compiled, /* seed */ 1234);
}
```

The helper runs once, captures the trace, replays through
`replay_from_trace_records`, and asserts the replayed trace
matches record-for-record.

---

## Biological coverage (if relevant)

If your pass can break productivity, add a stack to
[`tests/test_productive_stress_matrix.py`](../tests/test_productive_stress_matrix.py):

```python
def stack_my_mechanism(config: str, *, productive: bool) -> ga.Experiment:
    exp = ga.Experiment.on(config).recombine().my_mechanism(...)
    return exp.productive_only() if productive else exp


# Add to PRODUCTIVE_STACKS so the contract-holds test sweeps it.
PRODUCTIVE_STACKS.append(
    ("my_mechanism", lambda c, *, productive: stack_my_mechanism(c, productive=productive))
)

# Add to NEGATIVE_CONTROL_STACKS too if the mechanism CAN break
# productivity (i.e. without productive_only some records will be
# non-productive across 30 seeds).
```

That single appendage gets you parametrised coverage across VJ
and VDJ, with and without `productive_only`, plus the trace-replay
leg.

---

## Postcondition validator (the last step)

Every new mechanism must prove that the records it produces agree
with the engine's own truth oracle. Run the postcondition validator
over a representative stack and assert it passes:

```python
def test_my_mechanism_records_validate():
    exp = (
        ga.Experiment.on("human_igk")
        .recombine()
        .my_mechanism(...)
        .productive_only()
    )
    refdata = exp.refdata
    result = exp.run_records(n=100, seed=4242)
    report = result.validate_records(refdata)
    assert report, report.summary()
```

`report.summary()` returns a histogram of issue kinds across the
batch — if your mechanism breaks the engine-vs-projection
agreement on counters, junction, or allele calls, this test
surfaces the divergence with the failing record indices.

Use a VJ chain config (`human_igk` / `human_igl`) unless your
mechanism specifically requires D-segment behaviour — the C4
allele oracle is currently incomplete on D under non-zero trim
(known gap, follow-up slice). For VDJ-specific testing, filter
out `AlleleCallTieSetMismatch` on D as the release-tier test
does (see [`tests/test_release_validation.py`](../tests/test_release_validation.py)).

This is the **final required step** in the contributor checklist:

1. Pass emits events through `SimulationBuilder`.
2. Replay round-trip reproduces sequence + AIRR coords + events.
3. Distribution invariant tested if stochastic.
4. Biological stress matrix updated if it can break productivity.
5. **`result.validate_records(refdata)` passes** on a representative
   stack.

---

## Anti-patterns CI catches

The lockdowns and policy tests will fail with a structured
message and a pointer to the canonical fix. The full list lives
in [`engine_architecture.md` §3](engine_architecture.md#3-anti-patterns).
The short version:

| Anti-pattern | Test that catches it |
|---|---|
| `sim.with_*` in production code | [`no_pass_calls_persistent_with_star_mutators_in_production_code`](../engine_rs/src/passes/mutation_transaction/mod.rs) |
| `builder.change_base` / `insert_indel` / `delete_indel` outside `MutationTransaction` | [`no_pass_calls_low_level_builder_mutators_directly`](../engine_rs/src/passes/mutation_transaction/mod.rs) |
| Declared compile effect with no matching event | [`builtin_passes_emit_events_consistent_with_declared_compile_effects`](../engine_rs/src/passes/mod.rs) |
| Replay force-apply without re-validation | [`replay_rejects_id_*`](../engine_rs/src/passes/sample_allele.rs) per pass |
| Live-call refresh consuming `PassCompileEffect` instead of `SimulationEvent` | [`refresh_follows_events_not_effects_effect_without_event_skips_refresh`](../engine_rs/src/compiled/tests/live_call_edits.rs) |

---

## Existing passes to copy from

When in doubt, copy from the closest existing pass:

| Mechanism                       | Reference pass                                                                 |
|---------------------------------|--------------------------------------------------------------------------------|
| Sample an allele                | [`SampleAllelePass`](../engine_rs/src/passes/sample_allele.rs)                 |
| Apply a trim                    | [`TrimPass`](../engine_rs/src/passes/trim.rs)                                   |
| Assemble a germline segment     | [`AssembleSegmentPass`](../engine_rs/src/passes/assemble_segment/execution.rs) |
| Generate an NP region           | [`GenerateNPPass`](../engine_rs/src/passes/generate_np/execution.rs)           |
| Substitution-style SHM / errors | [`UniformMutationPass`](../engine_rs/src/passes/mutate/uniform.rs), [`S5FMutationPass`](../engine_rs/src/passes/mutate/s5f.rs), [`PCRErrorPass`](../engine_rs/src/passes/corrupt/pcr.rs), [`QualityErrorPass`](../engine_rs/src/passes/corrupt/quality.rs), [`NCorruptionPass`](../engine_rs/src/passes/corrupt/ncorrupt.rs), [`ContaminantPass`](../engine_rs/src/passes/corrupt/contaminant.rs) |
| Indels                          | [`IndelPass`](../engine_rs/src/passes/corrupt/indel.rs)                         |
| End-loss / primer trim          | [`EndLossPass`](../engine_rs/src/passes/corrupt/end_loss.rs)                    |
| Flag-only metadata pass         | [`RevCompPass`](../engine_rs/src/passes/corrupt/rev_comp.rs)                    |
| Single-Bool decision that mutates an assignment | [`InvertDPass`](../engine_rs/src/passes/invert_d.rs) (V(D)J D-segment inversion; emits `OrientationChanged` and commits via `SimulationBuilder::update_allele_orientation`) |

Every one of these passes:

- routes mutations through the event-emitting builder,
- validates replayed values against the admissibility chain,
- has a per-pass event-emission test,
- appears in the [coverage audit plan](../engine_rs/src/passes/mod.rs).

If your new pass looks unlike all of them, something has probably
gone off-template — re-read [`engine_architecture.md`](engine_architecture.md)
or open an issue.
