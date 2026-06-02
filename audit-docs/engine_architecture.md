# Engine architecture

This document codifies the invariants and contributor-facing
conventions of the Rust simulation kernel under
[`engine_rs/`](../engine_rs/). It exists so the load-bearing
architectural decisions live somewhere greppable, not just in
tests and conversation history.

The code enforces most of these rules — CI lockdowns, policy
checks, divergence tests — but the doc teaches the right way
before CI yells.

> **Audience.** Contributors adding a new pass, a new contract, a
> new derived-state consumer, or any code that mutates the
> [`Simulation`](../engine_rs/src/ir/simulation.rs) IR.

---

## 1. The seven invariants

The engine's correctness rests on seven layered invariants. Each
has a canonical type and a test that proves it.

### 1.1 Contracts constrain *support* before proposals

Sampling is **constrain-before-propose**, not propose-and-retry.
Every contract-aware sampler narrows its distribution's support
through a filter pass *before* drawing a candidate. A draw that
returns a value is always admissible by construction.

- **Trait:** [`Contract::admits_typed`](../engine_rs/src/contract/mod.rs)
- **Set:** [`ContractSet`](../engine_rs/src/contract/set.rs)
- **Typical site:** `sample_filtered_result(rng, dist, |candidate| contracts.admits_typed(...).is_ok())`

A permissive mode failure surfaces as the pass's *explicit
sentinel* (e.g. `b'N'` for ambiguous-base sampling, no-op skip
for indel position sampling). A strict mode failure surfaces as
[`PassError::ConstraintSampling`](../engine_rs/src/pass/error.rs).

### 1.2 Trace = choices/proposals

The trace records **what was proposed or sampled** at each
addressed decision point. It is the persisted artifact that
makes runs reproducible across rebuilds.

- **Type:** [`Trace`](../engine_rs/src/trace.rs), [`ChoiceValue`](../engine_rs/src/trace.rs)
- **Address vocabulary:** [`ChoiceAddress`](../engine_rs/src/address.rs)
- **File format:** [`TraceFile`](../engine_rs/src/trace_file.rs) — schema-versioned,
  refdata-content-hashed, address-schema-versioned

Two traces are *value-equal* iff every recorded `(address,
ChoiceValue)` pair matches in order.

### 1.3 Replay = validated proposal consumption

Trace-injected replay (Option B) consumes the recorded values
through a [`TraceCursor`](../engine_rs/src/replay.rs) at each
sampling site. The cursor's strict-positional contract means
the cursor *must* be drained when the plan completes — trailing
records are a structured `PassError::Replay`.

**Critical rule** — every replayed value must pass the same
admissibility chain a fresh draw would: refdata lookup,
distribution-support membership, contract acceptance,
feasibility acceptance. The trace supplies proposals; the engine
decides whether they apply. See [`SampleAllelePass::validate_replay_candidate`](../engine_rs/src/passes/sample_allele.rs)
for the canonical example.

**`strict=True` on replay does not re-raise permissive sentinels.**
Sentinel values (indel `site = -1` NoOp, NP length `0`, NP base
`N`, trim `0`) are themselves valid trace records — they represent
"the pass committed nothing" rather than "the pass committed
something the contract rejects." Replay consumes them verbatim and
they reproduce the original outcome. A user who wants the
strict-fresh-sampling behaviour for a recorded seed should call
`simulator.run(seed=<trace.seed>, strict=True)`, not
`replay_from_trace_file(tf, strict=True)`. See
[`docs/productive_failure_mode_audit.md`](productive_failure_mode_audit.md)
§5 and §6.2 for the failure matrix and the test pin.

### 1.4 Simulation events = consequences

[`SimulationEvent`](../engine_rs/src/ir/sim_event.rs) is the
typed runtime stream describing **what actually happened** to
the IR. Variants today:

| Variant                             | Carried payload                                      |
|-------------------------------------|------------------------------------------------------|
| `BasePushed`                        | handle, base, segment, germline_pos, flags           |
| `BaseChanged`                       | handle, old_base, new_base, segment, germline_pos    |
| `IndelInserted` / `IndelDeleted`    | at, base/removed_base, segment, flags                |
| `AssignmentChanged`                 | segment, old, new                                    |
| `TrimChanged`                       | segment, end, old, new                               |
| `RegionAdded` / `RegionReplaced`    | region or (old, new)                                 |
| `MutationCountChanged`              | old, new, delta                                      |
| `ReverseComplementFlagRecorded`     | applied                                              |
| `BaseDeleted`                       | at, removed_base (reserved)                          |

Events are emitted by [`SimulationBuilder`](../engine_rs/src/ir/builder.rs)
at every mutation site and fanned out to attached
[`SimulationEventSink`](../engine_rs/src/ir/sim_event.rs)
implementors. The compiled executor captures the per-pass event
stream into [`EventRecord.simulation_events`](../engine_rs/src/event.rs).

### 1.5 Compile effects = scheduling facts

[`PassCompileEffect`](../engine_rs/src/pass/metadata.rs) is the
**static, declarative** category of state change a pass produces.
It's read at compile time by the schedule analyzer for ordering
and dependency reasoning. It is **never** consulted at runtime
by derived-state refresh.

A transitional `pub type PassEffect = PassCompileEffect` alias
keeps the old spelling compiling; new code should use the
canonical name. The static-vs-runtime split is asymmetric on
purpose: declarations are *intent*; events are *what happened*.

### 1.6 Live-call refresh follows events

The post-pass V/D/J live-call refresh
([`LiveCallRefreshHook`](../engine_rs/src/live_call/refresh_hook.rs))
reads the pass's `simulation_events` stream — **not** the
`compile_effects` declarations — and runs the steps produced by
[`LiveCallRefreshPlan::from_events`](../engine_rs/src/live_call/refresh_plan.rs).
A pass that declares an effect but emits no event triggers no
refresh; a pass that emits an event without declaring the
matching effect still triggers refresh.

Two divergence tests pin this: [`refresh_follows_events_not_effects_effect_without_event_skips_refresh`](../engine_rs/src/compiled/tests/live_call_edits.rs)
and `..._event_without_effect_triggers_refresh`.

### 1.7 Built-in mutating passes route through `SimulationBuilder`

The only sanctioned mutation path for production pass code is the
event-emitting [`SimulationBuilder`](../engine_rs/src/ir/builder.rs)
API:

| Pass intent                  | Builder method                          |
|------------------------------|------------------------------------------|
| Assign allele                | `builder.assign_allele(segment, instance)`|
| Update trim                  | `builder.update_trim(segment, end, value)`|
| Add region                   | `builder.add_region(region)`              |
| Replace region               | `builder.replace_region(replacement)`     |
| Bump mutation count          | `builder.set_mutation_count(new_count)`   |
| Record rev-comp flag         | `builder.record_reverse_complement_flag(applied)` |
| Substitute / insert / delete a base | route through `MutationTransaction` |

Calling `Simulation::with_*` directly from production pass code
bypasses event emission, so the
[`LiveCallRefreshHook`](../engine_rs/src/live_call/refresh_hook.rs)
sees nothing and the live-call sidecar goes stale.

Two CI lockdowns make this unrepresentable:

- [`no_pass_calls_low_level_builder_mutators_directly`](../engine_rs/src/passes/mutation_transaction/mod.rs)
  — guards the low-level `change_base` / `insert_indel` / `delete_indel`
  builder methods.
- [`no_pass_calls_persistent_with_star_mutators_in_production_code`](../engine_rs/src/passes/mutation_transaction/mod.rs)
  — guards the persistent `with_allele_assigned` / `with_trim` /
  `with_region_added` / `with_region_replaced_for_segment` /
  `with_mutation_count` mutators. Test code inside `#[cfg(test)]`
  regions and under `tests/` subtrees stays free.

---

## 2. How to add a new pass

A new pass that mutates biology must touch six surfaces. The
template below is the canonical order; each step has an
existing analog you can copy from.

> **Practical companion.** [`docs/adding_a_pass.md`](adding_a_pass.md)
> is the copy-pasteable version of this section: minimal pass
> template, three required test patterns with helpers from
> [`passes::test_support`](../engine_rs/src/passes/test_support.rs),
> and a per-mechanism crib sheet of which existing pass to copy
> from. Read it when you're about to write code; come back here
> for the *why*.

### 2.1 Declare typed choice address patterns

If the pass samples anything, every sampling site needs a
`ChoiceAddress` variant in
[`engine_rs/src/address.rs`](../engine_rs/src/address.rs).

```rust
fn declared_choice_patterns(&self) -> Vec<ChoiceAddressPattern> {
    vec![ChoiceAddressPattern::MyNewPassCount, ...]
}
```

The address vocabulary is schema-versioned via
`ADDRESS_SCHEMA_VERSION` (see [`address.rs`](../engine_rs/src/address.rs));
adding a variant is additive. Pin the spelling via the
`frozen_address_spellings_for_choice_address_schema_v1` test.

### 2.2 Declare compile effects — for the scheduler only

```rust
fn effects(&self) -> Vec<PassCompileEffect> {
    vec![PassCompileEffect::EditBases]   // or AssembleSegment(seg), etc.
}
```

This is **static**. The scheduler reads it to order this pass
correctly relative to its dependencies. It is *not* read at
runtime for derived-state refresh.

### 2.3 Emit events through `SimulationBuilder`

Mutate via the builder, not `Simulation::with_*`. For
substitutions and indels, route through
[`MutationTransaction`](../engine_rs/src/passes/mutation_transaction/mod.rs).
For non-pool consequences (assignment, trim, region, mutation
count, rev-comp), call the typed builder methods listed in
[§1.7](#17-built-in-mutating-passes-route-through-simulationbuilder).

If you build your own `SimulationBuilder` instance (e.g. to
isolate a region-add step) and the caller passed
`PassContext.event_log_sink`, forward the captured events:

```rust
let mut builder = SimulationBuilder::from_simulation(sim);
if ctx.event_log_sink.is_some() {
    builder.attach_event_log_observer();
}
builder.add_region(region);
if let Some(sink) = ctx.event_log_sink.as_deref_mut() {
    sink.extend(builder.seal_event_log_observer());
}
let sim = builder.seal();
```

### 2.4 Validate replayed values

When the pass is in replay mode (`ctx.replay_cursor.is_some()`),
the cursor supplies the value. Run the **same admissibility
chain a fresh draw would** before applying. The canonical
template is [`SampleAllelePass::validate_replay_candidate`](../engine_rs/src/passes/sample_allele.rs):

1. Look up resource (refdata, etc.).
2. Check the value lies in the natural distribution's support.
3. Run `contracts.admits_typed(...)`.
4. Run `feasibility.admits(...)`.

Mismatches surface as the structured `PassError` variant that
best names *what* failed (`MissingAllele`,
`InvalidDistributionOutput { reason: ... }`, `ContractViolation`).

### 2.5 Add contract support if the pass needs to be filter-aware

If the pass should respect the productive bundle (or any future
contract bundle), wire it into the
[`sample_filtered_result`](../engine_rs/src/dist/filtered.rs)
pattern its category uses. Mutation/PCR/quality passes go through
[`MutationTransaction::substitute_position_constrained`](../engine_rs/src/passes/mutation_transaction/substitution.rs);
indel passes go through `MutationTransaction::insert_base` /
`delete_base_admitting`.

### 2.6 Add coverage

Three layers of test coverage:

- **Pass-level event-log test** — exercise the pass through
  `PassRuntime::execute_with_refdata` (or a direct
  `PassContext` invocation), capture
  `EventRecord.simulation_events`, assert the expected event
  variants fire. See [`outcome_event_records_carry_assignment_trim_and_region_events`](../engine_rs/src/passes/mod.rs)
  for the template.
- **Event-emission policy entry** — if the pass declares any
  `PassCompileEffect`, add it to [`build_coverage_audit_plan`](../engine_rs/src/passes/mod.rs)
  so [`builtin_passes_emit_events_consistent_with_declared_compile_effects`](../engine_rs/src/passes/mod.rs)
  pins the declared/emitted relationship.
- **Biological coverage** — if the pass mutates anything that
  could break the productive contract, add a stack to
  [`tests/test_productive_stress_matrix.py`](../tests/test_productive_stress_matrix.py)
  exercising it under `productive_only` and as a negative
  control.

---

## 3. Anti-patterns

The patterns below are caught by tests today. The doc lists them
so a new contributor sees the rule before tripping the CI.

### 3.1 Direct `sim.with_*` in production passes

> ❌ `sim.with_allele_assigned(...)` / `sim.with_trim(...)` /
> `sim.with_region_added(...)` / `sim.with_mutation_count(...)` /
> `sim.with_base_changed(...)` / `sim.with_indel_inserted(...)`

The runtime [`LiveCallRefreshHook`](../engine_rs/src/live_call/refresh_hook.rs)
follows events. A direct `with_*` call leaves the event stream
silent and the V/D/J live-call sidecar goes stale.

**Caught by:**
[`no_pass_calls_persistent_with_star_mutators_in_production_code`](../engine_rs/src/passes/mutation_transaction/mod.rs)
plus the existing low-level lockdown.

**Fix:** route through `SimulationBuilder` or `MutationTransaction`.

### 3.2 Replay force-apply

> ❌ A replay path that applies the recorded value without
> running the admissibility chain.

Replay does not mean "trust the trace." A recorded value can
have been invalidated by a refdata swap, an address-schema
revision, or a contract change since the trace was written.
Force-applying silently corrupts the IR; validation surfaces a
structured `PassError::Replay` or `PassError::MissingAllele`
that the caller can diagnose.

**Caught by:**
[`replay_rejects_id_missing_from_refdata`](../engine_rs/src/passes/sample_allele.rs)
and friends.

**Fix:** mirror [`SampleAllelePass::validate_replay_candidate`](../engine_rs/src/passes/sample_allele.rs).

### 3.3 Contract-violating permissive fallback

> ❌ Permissive mode falling back to an unconstrained draw.

When the contract narrows the support to empty, the canonical
permissive fallback is the pass's *explicit* sentinel — `b'N'`
for ambiguous-base, no-op skip for trim, etc. Falling back to
the unconstrained natural distribution defeats the purpose of
declaring the contract.

The one documented exception is [`SampleAllelePass::sample_allele`](../engine_rs/src/passes/sample_allele.rs)
(an empty allele pool can't sentinel-skip without panicking the
next assembly pass). That exception is spelled out in the
method-level doc and tested via the productive-only NP1 frame
filter.

### 3.4 Using `PassCompileEffect` for runtime refresh

> ❌ A new derived-state consumer that reads `effects` instead of `events`.

`PassCompileEffect` is the *intent* a pass declared at compile
time. `SimulationEvent` is the consequence stream the pass
actually emitted. A consumer reading effects sees declarations
the pass may not have honored; a consumer reading events sees
the truth.

**Caught by:**
[`refresh_follows_events_not_effects_effect_without_event_skips_refresh`](../engine_rs/src/compiled/tests/live_call_edits.rs)
plus the policy-conformance test.

**Fix:** consume `events: &[SimulationEvent]` in your hook's
`apply`. Use a [`LiveCallRefreshPlan`-style](../engine_rs/src/live_call/refresh_plan.rs)
event-to-step translator if the consumer logic is non-trivial.

### 3.5 Declaring an effect without emitting matching events

> ❌ `fn effects(&self) -> vec![PassCompileEffect::EditBases]` plus
> mutation via a path that doesn't fire `BaseChanged`.

The runtime refresh trusts events. Declaring an effect without
emitting the matching event variant silently skips refresh —
the same anti-pattern as direct `sim.with_*`, just phrased
differently.

**Caught by:**
[`builtin_passes_emit_events_consistent_with_declared_compile_effects`](../engine_rs/src/passes/mod.rs)
and the per-pass policy unit tests in [`event::tests`](../engine_rs/src/event.rs).

**Fix:** either route through the event-emitting builder, or
drop the effect declaration if the pass legitimately doesn't
produce that consequence.

---

## 4. Where the architecture lives

Quick navigation for new contributors:

| Concept                          | Canonical home                                                                 |
|----------------------------------|---------------------------------------------------------------------------------|
| Contract trait + set             | [`engine_rs/src/contract/`](../engine_rs/src/contract/)                         |
| Choice address vocabulary        | [`engine_rs/src/address.rs`](../engine_rs/src/address.rs)                       |
| Persisted trace                  | [`engine_rs/src/trace.rs`](../engine_rs/src/trace.rs), [`trace_file.rs`](../engine_rs/src/trace_file.rs) |
| Replay cursor                    | [`engine_rs/src/replay.rs`](../engine_rs/src/replay.rs)                         |
| Persistent IR                    | [`engine_rs/src/ir/`](../engine_rs/src/ir/)                                     |
| Event-emitting builder           | [`engine_rs/src/ir/builder.rs`](../engine_rs/src/ir/builder.rs)                 |
| Simulation event channel         | [`engine_rs/src/ir/sim_event.rs`](../engine_rs/src/ir/sim_event.rs)             |
| Compile effects + Pass trait     | [`engine_rs/src/pass/`](../engine_rs/src/pass/)                                 |
| Mutation transaction             | [`engine_rs/src/passes/mutation_transaction/`](../engine_rs/src/passes/mutation_transaction/) |
| Live-call refresh + plan         | [`engine_rs/src/live_call/refresh_hook.rs`](../engine_rs/src/live_call/refresh_hook.rs), [`refresh_plan.rs`](../engine_rs/src/live_call/refresh_plan.rs) |
| Event record / outcome ledger    | [`engine_rs/src/event.rs`](../engine_rs/src/event.rs)                           |
| Compiled execution loop          | [`engine_rs/src/compiled/execute.rs`](../engine_rs/src/compiled/execute.rs)     |
| Production passes                | [`engine_rs/src/passes/`](../engine_rs/src/passes/)                             |
| AIRR record projection           | [`engine_rs/src/airr_record/`](../engine_rs/src/airr_record/)                   |
| Productive-contract stress matrix| [`tests/test_productive_stress_matrix.py`](../tests/test_productive_stress_matrix.py) |

---

## 5. Glossary

- **Pass.** A scheduled unit of simulation work. Implements
  [`Pass`](../engine_rs/src/pass/traits.rs). Declares
  requirements, effects, and choice-address patterns; runs
  inside the compiled executor's transactional loop.
- **Plan.** Ordered list of passes via
  [`PassPlan`](../engine_rs/src/pass/) or the higher-level
  [`Experiment`](../src/GenAIRR/experiment.py) DSL on the
  Python side.
- **Contract.** A typed admissibility predicate over a
  candidate choice. Constraint-aware samplers narrow their
  distribution's support via contracts before drawing.
- **Effect.** Renamed [`PassCompileEffect`](../engine_rs/src/pass/metadata.rs).
  Static declaration of *category* of state change. Not a
  runtime fact.
- **Event.** [`SimulationEvent`](../engine_rs/src/ir/sim_event.rs).
  Typed runtime consequence emitted by the builder at each
  mutation site.
- **Outcome.** [`Outcome`](../engine_rs/src/pass/outcome.rs).
  The per-run result: revisions, trace, pass names, and the
  per-pass `EventRecord` ledger.
- **Live call.** The V/D/J allele-call sidecar maintained on
  each `Simulation`, refreshed by the
  [`LiveCallRefreshHook`](../engine_rs/src/live_call/refresh_hook.rs)
  after every pass that emits relevant events.
- **AIRR record.** The flat record dict produced at the Python
  boundary by [`outcome_to_airr_record`](../src/GenAIRR/_airr_record.py)
  — the user-visible projection of an `Outcome`.
