# Clonal Plan-Split — Architecture Contract

**Status: shipped.** This audit is the canonical contract for how
`Experiment.expand_clones(...)` partitions a pipeline into an
ancestor-phase plan (run once per clone) and a descendant-phase
plan (run once per read inside the clone). It documents the
corrected behaviour after the six clonal-correctness slices that
preceded it:

- **Bug A fix** — `invert_d` post-fork silent drop, guarded.
- **Bug B fix** — `receptor_revision` post-fork silent drop, guarded.
- **Bug C fix** — `paired_end` pre-fork silent r1/r2 empty, guarded.
- **Bug D fix** — receptor-revision provenance moved from trace
  to IR (`AlleleInstance.receptor_revision_original_id`).
- **Bug E fix** — `random_strand_orientation` pre-fork silent
  `rev_comp=False`, guarded via the unified descendant-phase table.
- **Bug F fix** — `end_loss_5prime/3prime` pre-fork silent
  `end_loss_*_length=0`, guarded via the same table.

This is a re-entry of the audit deferred earlier (the initial pass
surfaced the bugs above and was paused per the "stop and report"
protocol). The companion file
[`tests/test_clonal_plan_split_contract.py`](../tests/test_clonal_plan_split_contract.py)
pins every claim here as a positive `pin_scaffold_*` test plus
guard / divergence pins. Deferred-architecture absences are pinned
with `pin_absence_*` so a future implementation slice flips them
in lockstep.

When a future slice changes how clonal lowering works, this doc +
its contract file are the change-control surface — update both,
then any other layer that needs to follow.

---

## Existing scaffolding the contract pins

Surfaces the architecture depends on. Each is pinned by a
`pin_scaffold_*` test in the companion contract file.

| Surface | Where | What it pins |
|---|---|---|
| `_ClonalForkStep` marker | [`src/GenAIRR/_pipeline_ir.py:105`](../src/GenAIRR/_pipeline_ir.py#L105) | Structural partition point; frozen dataclass with `n_clones` + `size`. |
| `Experiment.compile()` fork-split | [`src/GenAIRR/experiment.py:1691-1742`](../src/GenAIRR/experiment.py#L1691-L1742) | Splits `_steps` at the fork; returns `CompiledClonalExperiment`. |
| `CompiledClonalExperiment.run_records` orchestration | [`src/GenAIRR/_compiled.py:436-454`](../src/GenAIRR/_compiled.py#L436-L454) | Per-clone loop: `pre.run(clone_seed)` → `final_simulation()` → `post.run_from(parent_sim, desc_seed)`. |
| `_descendant_phase_step_classifier` | [`src/GenAIRR/experiment.py`](../src/GenAIRR/experiment.py) (module-level helper) | Single source of truth for "is this step descendant-phase?" Drives the unified `expand_clones` guard. |
| `Experiment._has_clonal_fork()` | [`src/GenAIRR/experiment.py`](../src/GenAIRR/experiment.py) | Used by `invert_d` / `receptor_revision` for the opposite-direction (pre-fork only) check. |
| `AlleleInstance.receptor_revision_original_id` | [`engine_rs/src/assignment.rs:94-178`](../engine_rs/src/assignment.rs#L94-L178) | IR-side provenance the Bug D fix added; `Some(id)` when revision applied, `None` otherwise. Survives parent→descendant boundary via `Simulation.assignments`. |
| `SimulationResult.parents` | [`src/GenAIRR/result.py`](../src/GenAIRR/result.py) (Slice 2) | One parent `Outcome` per clone retained for replay/inspection; flat `.outcomes` keeps only descendants. |
| `SimulationResult.validate_families` / `validate_families_with_parents` | [`src/GenAIRR/result.py`](../src/GenAIRR/result.py) | The two family-layer validators; field-only and parent-aware respectively. |

---

## 1. Current state — the partition

`Experiment.expand_clones(n_clones, per_clone)` appends a
`_ClonalForkStep` marker to `_steps`. At
`Experiment.compile()` the step list is split at the marker
into a **pre-fork** half (per-clone) and a **post-fork** half
(per-descendant); each half compiles into its own
`CompiledSimulator`. The two compiled artifacts plus the fork
metadata are bundled into a `CompiledClonalExperiment`.

`CompiledClonalExperiment.run_records(seed=S)` runs the
canonical orchestration:

```python
for clone_idx in range(self._fork.n_clones):
    clone_seed = int(seed) + clone_idx * 1_000_000
    parent = self._pre.run(seed=clone_seed, strict=strict)
    parents.append(parent)                       # Slice 2 retention
    parent_sim = parent.final_simulation()
    for desc_idx in range(self._fork.size):
        desc_seed = clone_seed + 1 + desc_idx
        desc = self._post.run_from(
            parent_sim, desc_seed, strict=strict
        )
        outcomes.append(desc)
        rec = outcome_to_airr_record(
            desc, self._refdata,
            sequence_id=f"clone{clone_idx}_desc{desc_idx}",
        )
        rec["clone_id"] = clone_idx
        rec["parent_id"] = clone_idx
        ...
```

Boundary handoff: only `parent.final_simulation()` (a `Simulation`
IR snapshot) crosses into descendants. Each descendant gets a
fresh `Trace::new()` and `Rng::new(desc_seed)` from
`execute_transactional`. The parent's trace, events, and revision
history stay on the parent `Outcome`, which Slice 2 retains in
`SimulationResult.parents` rather than dropping.

---

## 2. Ancestor-phase steps (run pre-fork, once per clone)

DSL methods whose effects must be inherited by every descendant
of a clone — recombination identity, allele provenance, structural
decisions. Lowered into the pre-fork plan; rejected at the DSL
boundary by their method-level guards if appended **after**
`expand_clones`.

| DSL method | Step type | Why ancestor-phase |
|---|---|---|
| `recombine()` | `_RecombineStep` | Establishes V/D/J + trim + NP — the family identity. |
| `invert_d(prob)` | `_InvertDStep` (inlined into `_lower_recombine`) | D orientation is recombination-time; persisted on `AlleleInstance.orientation`, inherited via assignments. |
| `receptor_revision(prob)` | `_ReceptorRevisionStep` (inlined into `_lower_recombine`) | V replacement is recombination-time; provenance persisted on `AlleleInstance.receptor_revision_original_id` (Bug D fix), inherited via assignments. |
| Productive contracts (via `productive_only()`) | Contract bundle on the compiled artifact, NOT a pipeline step | Shapes recombination-time choice admissibility; runs on every pass including pre-fork. |

### Ancestor-phase guards

- `invert_d()` rejects if `_has_clonal_fork()` is `True` —
  `"invert_d must be called before expand_clones(); D inversion
  is a recombination-time decision and must be inherited by all
  clone descendants."`
- `receptor_revision()` rejects if `_has_clonal_fork()` is
  `True` — `"receptor_revision must be called before
  expand_clones(); receptor revision is a recombination-time
  decision and must be inherited by all clone descendants."`

These are **method-level** guards (the opposite direction from the
unified descendant-phase guard at `expand_clones`). They live on
the methods themselves because each individual call needs to know
"has the fork already been placed?" before appending.

### Productive contracts — a clarification

`productive_only()` is **not** a pipeline step; it's a constraint
bundle attached to the compiled artifact via
`_build_contracts()`. The bundle runs on every pass — pre-fork
recombination AND post-fork mutation/corruption — gating
admissibility. There is no pre/post fork distinction to enforce:
the same contracts apply across the boundary. The compile-time
precondition check for `ProductiveJunctionFrame` only runs on the
pre-fork plan (per
[`engine_rs/src/compiled/analyze.rs:58-68`](../engine_rs/src/compiled/analyze.rs#L58-L68)
— `has_recombination_support` returns `false` for the post-fork
plan, so the productive-frame precondition is skipped there; the
runtime contract bundle stays active on both halves).

---

## 3. Descendant-phase steps (run post-fork, once per descendant)

DSL methods whose effects must be sampled independently per clone
member. Lowered into the post-fork plan; rejected at the DSL
boundary by the unified guard at `expand_clones` if any of these
were appended **before** the fork.

| DSL method | Step type | Why descendant-phase |
|---|---|---|
| `mutate(...)` | `_MutateStep` | SHM is a per-cell process; descendants of a clone diverge through independent mutation events. |
| `pcr_amplify(...)` | `_CorruptStep(kind="pcr")` | PCR errors are stochastic per template; each clone member's amplicon should diverge independently. |
| `sequencing_errors(...)` | `_CorruptStep(kind="quality")` | Quality-error sites are per-read. |
| `polymerase_indels(...)` | `_CorruptStep(kind="indel")` | Indel sites are per-template. |
| `end_loss_5prime(length)` | `_CorruptStep(kind="5prime_loss")` | Observation-stage end loss is per-read. Bug F's silent-drop mechanism: trace-sourced `end_loss_5_length` projects to 0 on descendants when pre-fork. |
| `end_loss_3prime(length)` | `_CorruptStep(kind="3prime_loss")` | Same as 5' — Bug F applies. |
| `ambiguous_base_calls(...)` | `_CorruptStep(kind="ns")` | N-corruption sites are per-read. |
| `random_strand_orientation(prob)` | `_CorruptStep(kind="rev_comp")` | Strand orientation is per-fragment. Bug E's silent-drop mechanism: trace-sourced `rev_comp` projects to `False` on descendants when pre-fork. |
| `paired_end(...)` | `_PairedEndStep` | r1/r2 fragmentation is per-read. Bug C's silent-drop mechanism: trace-sourced r1/r2 fields project empty on descendants when pre-fork. |

### Descendant-phase guard

The unified guard at `expand_clones(...)` scans `_steps`,
classifies each step via
`_descendant_phase_step_classifier(step)`, and rejects the first
match. The message names the offending DSL method and the
canonical fix:

- For `mutate`:
  `"mutate must be called after expand_clones(); SHM is
  descendant-specific in GenAIRR's current clonal model. Move
  mutate(...) after expand_clones(...)."`
- For the other 8 methods, uniform phrasing:
  `"{method} must be called after expand_clones(); it is
  descendant-specific and must be sampled independently for each
  clone member. Move {method}(...) after expand_clones(...)."`

### Unguarded: `contaminate(...)`

`contaminate(prob=...)` (lowered as
`_CorruptStep(kind="contaminant")`) is deliberately **not** in
the descendant-phase table. The placement semantics weren't
relitigated in the Bug E/F slice; the pin file documents the
omission so a follow-up can add the classification when the
biology is decided.

### Non-clonal use stays unguarded

The guard fires only at `expand_clones`. Non-clonal pipelines can
freely use any descendant-phase step. The biological rationale:
in a non-clonal experiment there's no clone-level identity to
violate; each record is its own ancestor-and-read.

---

## 4. Parent / child trace contents

The boundary is a `Simulation` IR snapshot; the `Trace`, the
event ledger, and the per-pass revision history live on the
parent `Outcome` and do NOT cross into the descendant. The
descendant's `execute_transactional` allocates a fresh
`Trace::new()` at the top
([`engine_rs/src/compiled/execute.rs:66`](../engine_rs/src/compiled/execute.rs#L66)).

### What lives on the parent trace

For a canonical full clonal pipeline `recombine + invert_d +
receptor_revision`, the parent's `pass_names()` lists all 14
recombination passes:

```
sample_allele.v / sample_allele.d / sample_allele.j
trim.v_3 / trim.d_5 / trim.d_3 / trim.j_5
assemble.v / generate_np.np1 / invert_d / assemble.d
generate_np.np2 / assemble.j / receptor_revision
```

`event_count() == 14`, `len(trace()) == 27` (the 27 addressed
choices the recombination machinery samples — V/D/J allele ids,
trims, NP1/NP2 lengths and bases, invert_d Bool,
receptor_revision Bool + allele + trim).

### What lives on the descendant trace

For the same pipeline plus a canonical descendant tail
`mutate(count=5) + end_loss_5prime(length=10) +
random_strand_orientation(prob=0.5) + paired_end(...)`, every
descendant's `pass_names()` lists only the 4 post-fork passes:

```
mutate.s5f / corrupt.end_loss.5 / corrupt.rev_comp / paired_end
```

`event_count() == 4`, `len(trace()) == 16` (the descendant-
phase addressed choices — SHM count + sites + replacement bases,
end-loss length, rev_comp applied flag, paired-end r1/r2/insert
lengths).

### What does NOT live on either side

- The descendant trace does **not** duplicate the parent's trace.
- The parent trace does **not** carry descendant-phase records.
- In valid public pipelines (those that pass the
  descendant-phase guard), `corrupt.rev_comp` /
  `corrupt.end_loss.*` / `paired_end.*` trace records appear
  **only on descendants** — never only on a parent. The guard
  makes the "only on parent" state structurally unreachable from
  the public DSL.

### Parent count vs. record count

`len(SimulationResult.parents) == n_clones`. `len(.outcomes) ==
n_clones × per_clone`. Parents are stored once per clone, not
duplicated per descendant.

---

## 5. Projection inheritance

The AIRR projection
([`engine_rs/src/airr_record/builder.rs`](../engine_rs/src/airr_record/builder.rs))
reads from two sources: the descendant's IR (`Simulation`) and
its trace. Field-by-field the audit pins which source each
AIRR field uses, and therefore whether it survives the boundary.

### Inherited via IR (survives the boundary)

These fields read from `sim.assignments` / `sim.pool` /
`sim.sequence` — all part of `Simulation`, all cloned to the
descendant's `revisions[0]` at fork time:

| AIRR field | Source | Survives | Reason |
|---|---|---|---|
| `d_inverted` | `sim.assignments.d.orientation.is_reverse()` | Yes | Bug A fix's headline pin — ancestor decision inherited. |
| `receptor_revision_applied` | `sim.assignments.v.receptor_revision_original_id.is_some()` | Yes (Bug D fix) | Was trace-sourced and silently dropped on descendants; now IR-sourced. |
| `original_v_call` | refdata.get(V, `receptor_revision_original_id`).name | Yes (Bug D fix) | Same as above. |
| `v_call` / `d_call` / `j_call` (truth, modulo SHM widening) | `sim.assignments.{v,d,j}.allele_id` | Yes | Truth alleles inherited; live-call tie-set can widen post-SHM. |
| `truth_v_call` / `truth_d_call` / `truth_j_call` (with `expose_provenance=True`) | Same allele ids | Yes | Pinned by `test_g8_truth_columns_via_clonal_pipeline`. |
| `junction` (length, germline coords) | `sim.sequence` regions + assignments | Yes (pre-SHM coords) | Pre-SHM junction string is invariant; post-SHM bytes diverge. |

### NOT inherited — projected from descendant trace only

These fields read from the descendant's trace, which doesn't
carry pre-fork events. In valid pipelines (where the relevant
passes run post-fork) the descendant trace carries them
correctly. In invalid pipelines, the guards make pre-fork
placement structurally unreachable.

| AIRR field | Trace source | When populated |
|---|---|---|
| `rev_comp` | `ChoiceAddress::CorruptRevCompApplied` | Only when `random_strand_orientation` ran post-fork. |
| `end_loss_5_length` / `end_loss_3_length` | `corrupt.end_loss.5` / `corrupt.end_loss.3` | Only when the corresponding `end_loss_*prime` ran post-fork. |
| `r1_sequence` / `r2_sequence` / `r1_start` / `r2_start` / etc. | `paired_end.r1_length` / `paired_end.r2_length` / `paired_end.insert_size` | Only when `paired_end` ran post-fork. |
| `read_layout` | Same paired-end addresses | Same. |

### Parent projection — observation fields absent

Projecting a parent `Outcome` directly via
`outcome_to_airr_record(parent, refdata)` produces a record where:

- `d_inverted` / `receptor_revision_applied` / `original_v_call`
  are correctly populated (IR-sourced).
- `rev_comp == False`, `end_loss_*_length == 0`,
  `r1_sequence == ""`, `r2_sequence == ""`, `read_layout` empty
  — because the parent ran no observation-stage passes.

This makes the parent projection useful for ancestor-truth
queries (what V/D/J / orientation / revision did the family
ancestor have?) without confusing it with observation data.

---

## 6. Validator posture

Three validators currently observe the clonal output. The audit
pins their posture on the canonical full clonal pipeline:

### `validate_records(refdata)` — per-record AIRR postconditions

Iterates `(outcome, record)` pairs independently. Each
descendant is re-projected from its own outcome state and the
projected record is compared with the supplied record. **Stays
per-record only**; ignores `clone_id` and `parent_id`. Pinned by
the existing `pin_validate_records_stays_per_record_only` test
in `tests/test_clonal_family_contract.py`.

On the canonical full clonal stack, every record passes.

### `validate_families()` — field-only family checks

Groups records by `clone_id` and asserts the truth fields
(`truth_v_call`, `truth_d_call`, `truth_j_call`) agree across
siblings. Dict-only — works on records-only `SimulationResult`s
loaded from TSV. Skipped silently when the truth fields aren't
present (no `expose_provenance`).

On the canonical full clonal stack, the report is ok with
`family_count == n_clones` and the right
`members_per_family` map.

### `validate_families_with_parents(refdata)` — parent-aware

Reads from `result.parents` and compares each descendant's
truth fields, `d_inverted`, `original_v_call` against the
parent's projected values. Structural checks
(`ParentsMissing`, `ParentIdMissing`, `ParentIdOutOfRange`)
also run. With `refdata=None`, only structural checks run.

On the canonical full clonal stack, the report is ok with
`family_count == n_clones`.

### Wiring into `validate_records=True`

`CompiledClonalExperiment.run_records(validate_records=True)`
runs **two** gates after the batch is built:

1. Per-record postconditions via `validate_records(refdata)`.
2. Field-only family checks via `validate_families()`.

It does **NOT** invoke `validate_families_with_parents` —
parent-aware validation is an explicit deeper diagnostic the
caller invokes manually. Pinned by
`pin_present_parent_aware_validator` in
`test_clonal_parent_contract.py`.

---

## 7. Replay determinism

Per-descendant replay is fully deterministic: given
`(pre, post, refdata, seed)`, descendant `(clone_idx, desc_idx)`
reproduces byte-for-byte. The seed-derivation formula
`clone_seed = seed + clone_idx * 1_000_000` is internal but
stable (pinned by `pin_scaffold_clone_seed_formula_is_stable`
in `test_clonal_family_contract.py`).

Per-family replay is two-pass: replay the parent via
`compiled._pre.replay_from_trace_file(...)` to get the parent
IR, then call `compiled._post.run_from(parent_sim, desc_seed)`
per descendant. This works but reaches through `_pre` / `_post`
private attributes. A future `run_family(seed, clone_idx)`
helper would expose it publicly — deferred per audit §14.

The `trace_file_from` clonal-bundle format is **not** shipped
yet — pinned absent below.

---

## 8. Edge cases the contract covers

1. **Non-clonal experiments.** All descendant-phase methods work
   unguarded; `SimulationResult.parents is None`; family
   validators return ok no-op with `family_count == 0`.

2. **Empty post-fork plan.** `recombine().expand_clones(...)`
   with no descendant-phase steps after the fork produces
   byte-identical descendants — every clone member equals the
   parent IR. Family validator accepts this as a degenerate but
   valid case.

3. **`productive_only()` + heavy SHM.** The post-fork plan can
   reject descendants under heavy SHM. Total record count
   guarantee (`n_clones × per_clone`) is structural; the
   `productive_junction_frame` precondition runs only on the
   pre-fork plan (per
   [`compiled/analyze.rs:208-221`](../engine_rs/src/compiled/analyze.rs#L208-L221)).

4. **`paired_end` + clonal.** Each descendant gets its own
   r1/r2 layout sampled from the configured distribution. The
   guard rejects pre-fork placement; the unguarded canonical
   placement runs paired-end last in the post-fork plan (after
   every biology / corruption pass).

5. **TCR + clonal.** `mutate` is BCR-only (rejected on TCR
   refdata at the DSL boundary), so a TCR clonal pipeline can't
   exercise SHM. Other descendant-phase steps are locus-agnostic.

6. **Strict-mode failure ordering.** The orchestration loop
   forwards `strict=strict` to both `pre.run` and
   `post.run_from`. A parent strict-error aborts the whole clone
   (no descendants emitted); a descendant strict-error aborts the
   batch from that point. Slice 2 retained the parent on
   successful clones; strict-aborted clones produce no parent
   entry.

7. **Contaminant placement unguarded.** `contaminate(...)` is
   not classified as descendant-phase in this slice. Pin the
   omission so the placement-policy decision is explicit.

8. **Multiple descendant-phase steps pre-fork.** The unified
   guard reports the **first** offender — the message points the
   user at the call they need to look at first. Pinned by
   `test_guard_reports_first_descendant_phase_step_encountered`.

---

## 9. Backwards compatibility

This audit does NOT introduce new public surface. The earlier
slices it documents (A/B/C/D/E/F + the unified guard) shipped
the user-visible changes:

- `SimulationResult.parents` (Slice 2).
- `record["parent_id"]` (Slice 2).
- `SimulationResult.validate_families_with_parents` (Slice 3).
- `FamilyValidationReport.parent_value` / `child_values`
  failure-dict keys (Slice 3).
- `Outcome.airr_record.receptor_revision_applied` /
  `original_v_call` flipped from trace-sourced to IR-sourced
  (Bug D).
- `Outcome.airr_record.details.source` for receptor-revision
  issues flipped from `trace:*` to
  `ir:assignments.v.receptor_revision_original_id` (Bug D).
- DSL ordering rejections at `invert_d` / `receptor_revision` /
  `expand_clones` (Bug A/B/C + E/F).

The audit's contract pins make those changes the documented
architecture; any future slice that touches clonal lowering must
update this doc + the contract file in lockstep.

---

## 10. Validator integration

(Slice 3 / Bug D landed the validator pieces; this audit only
pins the resulting posture.)

- Per-record validation: continues to fire automatically when
  `validate_records=True` is passed to the clonal
  `run_records`.
- Field-only family validation: ditto, runs after the per-
  record gate.
- Parent-aware family validation: explicit. Callers invoke
  `result.validate_families_with_parents(refdata)` manually
  when they want the deeper diagnostic.
- Cache parity: per-outcome (descendant only). The retained
  parent outcomes can be cache-parity-checked manually but the
  release workflow doesn't fold that in by default.

---

## 11. Implementation order (history)

Recorded here for cross-doc traceability; no new slice is
proposed.

1. **Slice 0 — `expand_clones` + `clone_id`** (pre-audit).
2. **Slice 1 — `validate_families` (field-only)**.
3. **Slice 2 — Parent-Outcome Read-Only Surface**
   (`SimulationResult.parents`, `parent_id` on records).
4. **Slice 3 — Parent-aware family validator**
   (`validate_families_with_parents`).
5. **Bug D fix — receptor-revision provenance in IR**
   (`AlleleInstance.receptor_revision_original_id`).
6. **Slice — DSL ordering guards for Bugs A/B/C**
   (`invert_d` / `receptor_revision` post-fork + `paired_end`
   pre-fork rejections).
7. **Slice — Unified descendant-phase guard (Bugs E/F + class)**
   (every descendant-phase step now rejects at `expand_clones`
   if placed pre-fork).
8. **This audit** — architecture contract, no implementation.

---

## 12. Test surface — what this audit pins

Mirrored in [`tests/test_clonal_plan_split_contract.py`](../tests/test_clonal_plan_split_contract.py).

### `pin_scaffold_*` — the corrected architecture

1. `_ClonalForkStep` exists as a frozen dataclass; `expand_clones`
   appends exactly one.
2. `Experiment.compile()` returns `CompiledClonalExperiment` iff
   a fork step is present.
3. The orchestration loop uses `pre.run(...)` →
   `final_simulation()` → `post.run_from(...)`.
4. `clone_seed = seed + clone_idx * 1_000_000` (re-pinned for
   cross-doc traceability).
5. **Ancestor-phase set** — `recombine`, `invert_d`,
   `receptor_revision` lower into pre-fork plan; their effects
   are inherited by descendants (`d_inverted=True`,
   `receptor_revision_applied=True`, `original_v_call`
   populated).
6. **Descendant-phase set** — every guarded method (`mutate`,
   `pcr_amplify`, `ambiguous_base_calls`, `sequencing_errors`,
   `polymerase_indels`, `end_loss_5prime`, `end_loss_3prime`,
   `random_strand_orientation`, `paired_end`) is classified by
   `_descendant_phase_step_classifier`.
7. **Ancestor-phase guards** — `invert_d` / `receptor_revision`
   after `expand_clones` rejects.
8. **Descendant-phase guards** — each descendant-phase method
   before `expand_clones` rejects at `expand_clones`.
9. **Non-clonal use unaffected** — all descendant-phase methods
   compile and run in non-clonal pipelines.
10. **Parent trace contents** — for a canonical configuration,
    `parent.pass_names()` contains `invert_d` and
    `receptor_revision`; trace and event counts match the
    configuration.
11. **Descendant trace contents** — for a canonical
    configuration, descendant `pass_names()` contains exactly
    the post-fork passes and nothing more.
12. **Descendant trace doesn't duplicate parent's** — strict
    inequality: descendant trace length < non-clonal equivalent.
13. **Parent count = clone count** — `len(.parents) == n_clones`,
    `len(.outcomes) == n_clones * per_clone`.
14. **IR-sourced provenance survives boundary** — `d_inverted`,
    `receptor_revision_applied`, `original_v_call` on every
    descendant match the parent's projected values.
15. **Observation fields absent on parent projection** —
    `rev_comp=False`, `end_loss_*_length=0`,
    `r1_sequence==""`, `r2_sequence==""` on the parent's AIRR
    projection.
16. **Descendant-phase divergence** — per-descendant
    independence (e.g. `rev_comp` differs across siblings at
    `prob=0.5`).
17. **All three validators pass** on the canonical full clonal
    stack.
18. **`validate_records=True` runs field-only family, not
    parent-aware** — pinned via source-inspection of the
    clonal `run_records` body.

### `pin_absence_*` — deferred architecture still pinned absent

19. `pin_absence_no_clonalfamily_aggregate` — no `ClonalFamily`
    / `FamilyOutcome` type exported from the public package.
20. `pin_absence_no_clonal_trace_file_bundle` — no
    `trace_file_from_family` / `clonal_trace_file_from` helper
    on `CompiledClonalExperiment`.
21. `pin_absence_no_pre_shm_junction_validator` — no
    `validate_pre_shm_junction_invariance` /
    `validate_families_pre_shm_junction` method.
22. `pin_absence_no_mutation_distance_aggregator` — no
    `mutations_per_clone` / `family_shm_summary` /
    `clone_distance_matrix` accessor on `SimulationResult`.

### Doc anchor

23. The audit doc continues to exist and references this
    contract file; the 14-section structure is intact.

---

## 13. Out of scope

Documented here so a future contributor doesn't accidentally
expand the contract.

- **`contaminate` classification.** Deferred — pin its omission,
  decide placement in a future slice.
- **`ClonalFamily` aggregate type.** Slice 4+ per the
  parent-outcome audit; pinned absent here.
- **Clonal trace-file bundle.** Slice 5+ per the parent-outcome
  audit; pinned absent.
- **Pre-SHM junction invariance validator.** Requires either a
  `FamilyRecord` projection or a parent-derived junction column;
  pinned absent.
- **Mutation-distance distribution aggregator.** Same — needs
  parent sequence projection or per-clone SHM-event grouping.
- **Lineage trees / GC selection / time evolution.** All
  compose against the `ClonalFamily` aggregate; out of scope
  until that lands.
- **Rust-side `Outcome.parent_id` back-pointer.** Pinned absent
  in `test_clonal_parent_contract.py`; not relitigated here.
- **`run_family(seed, clone_idx)` public helper.** Deferred per
  `test_clonal_parent_contract.py`.

---

## 14. Summary table

| Concern | Answer |
|---|---|
| Ancestor phase | `recombine`, `invert_d`, `receptor_revision`; productive contracts apply across the boundary. |
| Descendant phase | `mutate`, `pcr_amplify`, `sequencing_errors`, `polymerase_indels`, `ambiguous_base_calls`, `end_loss_5prime`, `end_loss_3prime`, `random_strand_orientation`, `paired_end`. |
| Ancestor-phase guards | Method-level on `invert_d` / `receptor_revision`. |
| Descendant-phase guards | Unified table at `expand_clones`; first offender named. |
| Non-clonal use | All descendant-phase methods work freely. |
| Parent trace | Recombination-time choices including `invert_d` / `receptor_revision` when configured. |
| Descendant trace | Descendant-phase choices only — never duplicates parent. |
| Pre-fork observation records | Cannot exist in valid public pipelines (guarded). |
| IR-sourced provenance | `d_inverted`, `receptor_revision_applied`, `original_v_call` survive boundary; descendants inherit them. |
| Observation projection fields | `rev_comp`, `end_loss_*_length`, paired-end fields — descendant-only; absent on parent projection. |
| Parent / record counts | `len(parents) == n_clones`, `len(outcomes) == n_clones * per_clone`. |
| `validate_records(refdata)` | Per-record only; ignores clonal structure. |
| `validate_families()` | Field-level; passes on canonical full stack. |
| `validate_families_with_parents(refdata)` | Parent-aware; passes on canonical full stack. |
| `validate_records=True` wiring | Per-record + field-only family; NOT parent-aware. |
| Deferred architecture | `ClonalFamily` aggregate, clonal trace-file bundle, pre-SHM junction validator, mutation-distance aggregator — all pinned absent. |

The clonal plan-split is now a **typed contract**: ancestor-phase
goes pre-fork, descendant-phase goes post-fork, the DSL rejects
misorder at the boundary, projection sources are pinned per field,
and validator posture is documented. Future biology slices
(lineage, selection, time evolution) compose against this
architecture instead of re-discovering it.
