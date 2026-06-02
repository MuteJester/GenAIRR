# Architecture (Contributor)

<p class="lead">A curated entrance to GenAIRR's engine
architecture. The full corpus — 39 audit and design docs at
<code>docs/*.md</code> — is contributor-facing and stays at the
repository root. This page is the map: the mental model, the
anchor documents, the audit-first workflow, and the checklist for
adding a new mechanism. Deep-link out from here when you need the
authoritative source.</p>

!!! tip "Your learning path"
    You're at the entry of the **"I'm contributing to GenAIRR"**
    path. After the engine mental model and the audit-first
    workflow on this page, the three anchor docs at
    `audit-docs/engine_architecture.md`, `audit-docs/adding_a_pass.md`, and
    `audit-docs/validation_matrix.md` are the contributor-tone
    references for actually writing kernel code.
    [See all paths →](../learn.md)

## Who this section is for

This section is **not** for new users. It assumes you're already
running simulations and need to look under the hood:

- **Contributors** adding new passes, contracts, AIRR fields, or
  cartridge planes.
- **Advanced users** debugging an unexpected output or reviewing
  a published GenAIRR dataset against its source.
- **Maintainers** auditing one mechanism's drift across a
  release, or scoping a refactor's blast radius before touching
  Rust kernel code.

If you're learning the simulator for the first time, the
[Getting Started](../getting-started/index.md) and
[Core Concepts](../concepts/index.md) tracks are the right
entries; come back here when something below those tracks
surprises you.

## The engine mental model

Three abstractions describe how the engine runs:

1. **`Experiment` lowers to a pass plan.** Every method call on
   the fluent `Experiment` builder appends one entry to a
   typed pass plan. `Experiment.compile()` resolves cartridge
   defaults, plan signatures, and contract sets into a frozen
   `CompiledExperiment`. The plan is the engine's input — not
   the `Experiment` itself.
2. **Passes produce three streams.** Each pass run by the
   simulator emits:
   - **Trace choices** — `(address, value)` records of every
     sampled draw, the substrate for replay.
   - **Simulation events** — typed records of *consequences*
     (base changed, segment trimmed, indel inserted, etc.), the
     substrate for derived counters and validators.
   - **A final `Simulation`** — the assembled molecule plus
     per-segment claim ledger plus `SegmentLiveCall` caches.
3. **AIRR records are projections.** The dict you get from
   `result[i]` is **derived** from `Simulation` + refdata via a
   pure projection layer (Rust: `airr_record::builder`). It
   carries no state of its own; rerunning the projection on the
   same `Simulation` produces a byte-identical record.

Validators don't trust projection — they re-derive every field
they cover from the underlying events and assert agreement. That
independence is the reason the AIRR validator catches projection
bugs the projector can't see.

## Core architecture documents

Three anchor documents at `docs/` are the canonical contributor
reference. Each is contributor-tone, not user-tone, and lives at
the repository root rather than inside the docs site:

| Document | Purpose |
|---|---|
| `audit-docs/engine_architecture.md` | The seven engine invariants — contracts constrain support, trace = choices, events = consequences, live-call refresh follows events, plan signature gates replay, refdata content hash gates content drift, projection is pure. Required reading before any kernel work. |
| `audit-docs/adding_a_pass.md` | The pass-author playbook: minimal pass template, the three required test patterns (deterministic-by-trace, plan-signature-stable, event-ledger-complete), the live-call refresh hook contract. Copy-paste-ready. |
| `audit-docs/validation_matrix.md` | The navigable map: every guarantee → its audit doc → its test file → the Rust kernel invariant. Use it to find the load-bearing tests for any behaviour you're about to change. |
| `audit-docs/plan_signature_completeness_audit.md` | The replay-safety audit: every parameterised DSL surface + every cartridge-driven compile parameter mapped against plan-signature participation. The reference when you wire a new parameter. |
| `audit-docs/docs_website_audit.md` | The information architecture audit of the docs corpus itself, plus the migration roadmap to the MkDocs site you're reading now. |

## The audit-first workflow

GenAIRR's release process is **audit-first**: a mechanism gets
specified, validated against the existing engine, and pinned by
contract tests *before* implementation begins. The pattern, in
order:

1. **Audit doc.** Open `audit-docs/<topic>_audit.md`. Define the
   biology, name the v1 boundary (what's in / what's deferred),
   enumerate the invariants the implementation must preserve,
   and list the drift section's open questions.
2. **Contract pins.** Add or extend `tests/test_*_contract.py`
   with the test patterns from `adding_a_pass.md` (
   deterministic-by-trace, plan-signature-stable,
   event-ledger-complete) and any audit-specific invariants. The
   pins fail on `master` because the implementation isn't there
   yet — that's the point.
3. **Implementation slices.** Land the mechanism in small,
   independently-reviewable slices. Each slice flips one set of
   contract pins from red to green. The audit doc's invariants
   guide what each slice must hold.
4. **Release consolidation.** Once every slice has shipped, a
   release commit renames `audit-docs/<topic>_audit.md` →
   `audit-docs/<topic>_design.md` (audit → design transition), updates
   `audit-docs/validation_matrix.md` with the new row, and stamps the
   audit as resolved.

**Why this matters for biological mechanisms.** Most simulators
treat invariants as informal — "the recombination should produce
a valid junction." GenAIRR makes them executable: the validator
re-derives the junction, the contract pins enforce
counter-partition equalities, and the plan signature catches
any kwarg that would silently shift output. The audit-first
discipline keeps that invariant set explicit so a new mechanism
can't quietly break an old one.

## Where mechanisms live

The Rust crate at `engine_rs/` carries the simulation kernel;
the Python wrapper at `src/GenAIRR/` carries the DSL surface,
cartridge handling, projection helpers, and validation. By
mechanism area:

| Mechanism area | Rust kernel | Python wrapper |
|---|---|---|
| **Recombination** | `engine_rs/src/compiled/recombine.rs`, `engine_rs/src/sample/allele.rs`, `engine_rs/src/junction.rs` | `Experiment.recombine`, `restrict_alleles`, allele-usage estimator |
| **Mutation (SHM)** | `engine_rs/src/compiled/mutate.rs`, `engine_rs/src/contract/productive_*.rs` | `Experiment.mutate`, `productive_only`, `segment_rates`, `v_subregion_rates` |
| **Corruption** | `engine_rs/src/compiled/pcr_amplify.rs`, `polymerase_indels.rs`, `sequencing_errors.rs`, `ambiguous_base_calls.rs`, `end_loss.rs` | `Experiment.pcr_amplify`, `polymerase_indels`, `sequencing_errors`, `ambiguous_base_calls`, `end_loss_*prime` |
| **Recombination editing** | `engine_rs/src/compiled/invert_d.rs`, `receptor_revision.rs` | `Experiment.invert_d`, `receptor_revision` |
| **Projection (AIRR)** | `engine_rs/src/airr_record/builder.rs`, `record.rs` | `SimulationResult`, `result.to_*` exporters, projection coordinate views |
| **Reference cartridge** | `engine_rs/src/dataconfig/*` (bridged), `engine_rs/src/refdata.rs` | `DataConfig`, `ReferenceCartridgeBuilder`, `cfg.cartridge_manifest()`, `cfg.compute_checksum()` |
| **Trace + replay** | `engine_rs/src/trace_file.rs`, `engine_rs/src/replay.rs` | `compiled.simulator.trace_file_from / replay_from_trace_file / rerun_from_trace_file` (PyO3 surface) |
| **Validation** | `engine_rs/src/contract/*`, `engine_rs/src/airr_record/validate.rs` | `result.validate_records`, `validate_families`, `validate_families_with_parents` |

When following an audit doc from the table above, the `docs/`
audit will name the kernel files it touches in its §3 or §4 —
use that as the authoritative source.

## Before you add a new mechanism

The checklist new contributors run before opening a slice PR:

- [ ] **Define the biology and provenance.** What does this
      mechanism model in vivo? What counters and provenance
      fields does it surface on the AIRR record? Write it in
      `audit-docs/<topic>_audit.md` §1-2 before writing code.
- [ ] **Trace / replay choices.** Every sampling site this
      mechanism introduces needs a stable address (`engine_rs/src/address.rs`)
      and a deterministic projection from `ChoiceValue`. Write
      the address schema in the audit before lowering it.
- [ ] **Event consequences.** Decide which event types this
      mechanism emits, and ensure they're sufficient to
      reconstruct every counter the projection will read.
      Re-derivability from events is the validator's contract.
- [ ] **Validators.** Add the contract pin that re-derives the
      mechanism's counters from events and asserts equality with
      the projected record. Without this, the AIRR validator
      can't see the new mechanism.
- [ ] **Cartridge ownership** (when applicable). If the
      mechanism takes per-cartridge parameters, decide which
      `reference_models` plane owns them, whether they
      participate in the plan signature, and whether they fold
      into `refdata_content_hash`. Document the decision in the
      audit's "cartridge boundary" section.
- [ ] **Docs / user guide.** A new user-facing mechanism gets a
      user-tone guide under `site_docs/guides/`. Use the existing
      mechanism guides as templates (e.g. SHM targeting, junction
      additions, corruption + sequencing artefacts). Reference
      the audit doc only at the bottom under "Deep architecture
      notes".
- [ ] **Release-tier test.** The two-layer integrity model
      (AIRR validator + live-call cache parity) must pass on the
      new mechanism's output before release. Add a release-tier
      pin if the mechanism touches the live-call refresh hook.

The `audit-docs/adding_a_pass.md` companion turns this checklist into
copy-paste templates; use it for the actual implementation.

## Two-layer integrity model

GenAIRR's release-readiness rests on two independent layers:

| Layer | Question | API |
|---|---|---|
| **AIRR validator** | Does the projected AIRR record agree with an independent re-derivation from the `Outcome`? | `result.validate_records(refdata)` → `ValidationReport` |
| **Cache parity** | Does the cached `SegmentLiveCall` on the final `Simulation` equal a from-scratch recompute? | `outcome.check_live_call_cache_parity(refdata)` → `List[dict]` |

Both layers run in release-tier CI. The AIRR validator catches
projection / counter drift; cache parity catches live-call
refresh-hook bugs that the AIRR validator can't see because
they don't surface on the record. See
`audit-docs/validation_matrix.md`
§1.0 for the rationale.

## Deep links

For convenience, the full set of GitHub-tree links to the
audit corpus. (MkDocs strict can't validate links outside the
site source tree, so these are absolute URLs.)

### Anchor docs (linked above)

- `audit-docs/engine_architecture.md` 
- `audit-docs/adding_a_pass.md` 
- `audit-docs/validation_matrix.md` 
- `audit-docs/plan_signature_completeness_audit.md` 

### Sample audits by mechanism area

The full corpus carries 39 markdown files. A few representative
audits to illustrate scope:

- `audit-docs/allele_call_audit.md` — recombination 
- `audit-docs/junction_n_addition_audit.md` — recombination 
- `audit-docs/mutation_provenance_audit.md` — mutation 
- `audit-docs/v_region_substructure_audit.md` — mutation 
- `audit-docs/indel_provenance_audit.md` — corruption 
- `audit-docs/airr_record_validator.md` — projection 
- `audit-docs/reference_cartridge_completeness_audit.md` — cartridge 

The full inventory is browsable at the repository root in
`docs/`; a Phase 5 audit-index page will enumerate everything
with one-line summaries.
