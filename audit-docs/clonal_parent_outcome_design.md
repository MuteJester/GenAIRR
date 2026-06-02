# Clonal Parent Outcome — Pre-Implementation Audit

**Status: audit only.** Pins today's parent-outcome lifecycle and
the architecture gaps any future "expose the ancestor" slice will
need to close. No implementation is proposed in this pass; the
deliverable is the shared vocabulary + the contract pins so a later
slice can land against an audited baseline.

Companion to
[`tests/test_clonal_parent_contract.py`](../tests/test_clonal_parent_contract.py)
which freezes today's behaviour (`pin_scaffold_*`) and the gaps the
audit identifies (`pin_absence_*`). When an implementation slice
lands, the `pin_absence_*` tests flip to `pin_*_present` /
`pin_*_used` in lockstep with the slice.

Sibling to [`docs/clonal_family_design.md`](clonal_family_design.md)
(the higher-level family/lineage audit) which already pinned the
flat-records-plus-`clone_id` shape. This audit zooms in on the
**parent-outcome lifecycle** specifically — what gets built, what
gets retained, what gets thrown away, and where the boundary
sits — so the next architectural slice (whatever shape it takes)
has a precise object to design against.

This audit is intentionally narrow. The family audit answered
"do we need a real family concept?" (yes); this one answers "where
exactly does the parent live and die today, and which surface
should expose it?"

---

## Existing scaffolding the audit relies on

Surfaces today's engine provides, each pinned by a `pin_scaffold_*`
test in the companion contract file.

| Existing surface | Where | What it gives a parent-outcome audit |
|---|---|---|
| `Simulation` IR (the parent → descendant carrier) | [`engine_rs/src/ir/simulation.rs:18-44`](../engine_rs/src/ir/simulation.rs#L18-L44) | The exact shape of what crosses the parent → descendant boundary today: `pool`, `sequence`, `assignments`, `segment_calls`, `dirty_log`, `mutation_count`. Trace and events are **not** in this struct. |
| `Outcome` (the per-run aggregate) | [`engine_rs/src/pass/outcome.rs:7-34`](../engine_rs/src/pass/outcome.rs#L7-L34) | The Rust-side run aggregate: `revisions: Vec<Simulation>`, `pass_names: Vec<String>`, `trace: Trace`, `events: Vec<EventRecord>`. The parent run produces one of these; the descendant run produces another, independent one. |
| `execute_transactional` | [`engine_rs/src/compiled/execute.rs:61-80`](../engine_rs/src/compiled/execute.rs#L61-L80) | Where the descendant's trace is allocated fresh: `let mut trace = Trace::new()`. Confirms the parent trace is structurally inaccessible from the descendant's outcome. |
| `CompiledSimulator::run_one_from_with_policy` | [`engine_rs/src/compiled/mod.rs:301-309`](../engine_rs/src/compiled/mod.rs#L301-L309) | The Rust entry point that takes a caller-supplied `initial: Simulation` and runs the plan from it. Parent → descendant data flow goes through this method, parameterised by `(parent.final_simulation(), desc_seed)`. |
| `Outcome.final_simulation()` | [`engine_rs/src/pass/outcome.rs:16-20`](../engine_rs/src/pass/outcome.rs#L16-L20) | The accessor the orchestration loop calls on the parent to obtain the IR snapshot the descendants start from. Returns `&Simulation` — a reference into `revisions.last()`. |
| `CompiledClonalExperiment.run_records` orchestration loop | [`src/GenAIRR/_compiled.py:380-391`](../src/GenAIRR/_compiled.py#L380-L391) and [`:436-454`](../src/GenAIRR/_compiled.py#L436-L454) | The Python-side orchestration where the parent `Outcome` is materialized, has `.final_simulation()` extracted, and is then **dropped** at the end of the per-clone iteration. The transient-parent lifecycle is structural here, not in Rust. |
| `_ClonalForkStep` partitioning at compile time | [`src/GenAIRR/experiment.py:1691-1742`](../src/GenAIRR/experiment.py#L1691-L1742) | Splits the step list into pre-fork (per-clone) and post-fork (per-descendant) plans. Pinned by the existing clonal-family contract. The audit references this for "which passes run before descendants branch." |
| `validate_families` (record-only family checker) | [`src/GenAIRR/result.py`](../src/GenAIRR/result.py) — Slice 1 of the family audit | Today's family-layer validator. Field-only; can't see anything the parent's `Simulation` carried that didn't make it onto descendant AIRR records. |

---

## 1. Current state — where is the ancestor today?

### 1.1 In Python orchestration

`CompiledClonalExperiment.run_records` runs this loop (paraphrasing
the implementation at
[`_compiled.py:436-454`](../src/GenAIRR/_compiled.py#L436-L454)):

```python
for clone_idx in range(self._fork.n_clones):
    clone_seed = int(seed) + clone_idx * 1_000_000
    parent = self._pre.run(seed=clone_seed, strict=strict)   # PyOutcome
    parent_sim = parent.final_simulation()                   # PySimulation
    for desc_idx in range(self._fork.size):
        desc_seed = clone_seed + 1 + desc_idx
        desc = self._post.run_from(parent_sim, desc_seed, strict=strict)
        outcomes.append(desc)
        # ...record projection + clone_id tag...
    # `parent` goes out of scope here at the end of the clone_idx iteration.
```

The parent's lifetime is **exactly the duration of one outer-loop
iteration**. It exists as a Python local, has its `final_simulation()`
extracted (a clone of the last `Simulation` revision), is fed to every
descendant via `run_from`, and is dropped at the top of the next
iteration. It is never:

- Stored on `self`.
- Returned to the caller.
- Inserted into the descendant-`outcomes` list.
- Referenced by any descendant `Outcome`.
- Persisted to a trace file.

### 1.2 In Rust execution

`run_one_from_with_policy` calls `execute_transactional`
([`execute.rs:61-80`](../engine_rs/src/compiled/execute.rs#L61-L80)),
which allocates fresh state at the top of the function:

```rust
let mut trace = Trace::new();
let mut rng = Rng::new(seed);
let mut revisions: Vec<Simulation> = Vec::with_capacity(inputs.plan.len() + 1);
revisions.push(initial);
```

The descendant gets a brand-new `Trace`, a brand-new `Rng` seeded
from `desc_seed`, and a fresh `revisions` vector whose **only
parent-derived entry** is `initial` — the parent's
`final_simulation()`, cloned. From this point on, parent and
descendant share no mutable state and no addressable identity.

### 1.3 The "one parent per clone" claim, confirmed

The user's question 1.1: *"Is there one parent simulation per
clone?"* Yes, exactly one. Each outer-loop iteration drives one
`pre.run(clone_seed)` call. There is no batching, no parent
sharing across clone indices, and no parent reuse. Different
`clone_idx` values produce different `clone_seed` values
(formula: `seed + clone_idx * 1_000_000`) and therefore different
parent `Outcome`s.

### 1.4 The "parent trace retained anywhere" question, answered

Question 1.2: *"Is the parent trace retained anywhere?"* No.
The parent `PyOutcome` carries `.trace()` while it is alive in
the Python local `parent`, but the orchestration loop never
copies, stores, or projects that trace. When the Python local is
rebound at the next iteration, the GC reclaims the entire parent
`PyOutcome` — including the trace, the revision history, and the
event ledger. None of it survives.

### 1.5 The "which passes run before descendants branch" question

Question 1.3: *"Which passes run before descendants branch?"*
Exactly the passes in `Experiment._steps[:fork_idx]` — by
default, just `recombine()` (which lowers to the canonical
`V → NP1 → D → NP2 → J` sequence plus any inlined `invert_d` /
`receptor_revision` probabilities — see
[`_compile.py:_extract_invert_d_prob`](../src/GenAIRR/_compile.py#L94-L132)
and
[`_extract_receptor_revision_prob`](../src/GenAIRR/_compile.py#L135-L180)).
Any DSL step appended **between** `.recombine()` and
`.expand_clones(...)` also runs pre-fork. The DSL doesn't enforce
"only recombine before the fork"; technically a user could call:

```python
exp = (Experiment.on("human_igh")
       .recombine()
       .mutate(count=10)            # ← would be per-clone, not per-descendant
       .expand_clones(n_clones=5, per_clone=4)
       .pcr_amplify(count=2))
```

The mutate step would then be **shared by every descendant** of
each clone, which is biologically wrong for SHM (SHM is meant to
be per-descendant). The DSL doesn't reject this today; the
compile-time fork-split silently accepts whatever the user put
pre-fork. This is a separate, smaller architecture gap — see
§8 edge case 4 below.

### Gap pin

`pin_absence_parent_outcome_unreachable_from_descendant` — given
a descendant `Outcome`, there is no API to ask "what was the
parent's trace / events / revision history?" The information
exists during orchestration and is structurally inaccessible
afterwards.

---

## 2. Q2 — What is actually shared?

The user's second question, answered against the precise
shape of the `Simulation` struct.

### What crosses the parent → descendant boundary

`parent.final_simulation()` returns a `PySimulation` wrapping a
cloned `crate::ir::Simulation`
([`ir/simulation.rs:18-44`](../engine_rs/src/ir/simulation.rs#L18-L44)).
The fields, and what each carries to the descendant:

| Field | Type | Carries to descendant? | Biological meaning |
|---|---|---|---|
| `pool` | `NucleotidePool` | Yes — full pool with all NP / V / D / J bytes | The complete assembled-sequence arena (every byte the V+NP1+D+NP2+J recombination put down) |
| `sequence` | `Arc<Sequence>` | Yes — the assembled region-mapped sequence | V / NP1 / D / NP2 / J coords + per-region cargo |
| `assignments` | `AlleleAssignments` | Yes — full V/D/J/C allele instances | Truth alleles + trims + (post-receptor-revision) `original_v_call` + (post-invert) `d_inverted` |
| `segment_calls` | `Arc<SegmentCalls>` | Yes — live-call evidence at fork time | Walker state representing what an aligner would call on the assembled IR |
| `dirty_log` | `Arc<DirtyLog>` | Yes — derived-state machinery | Per-pass dirty-window log used by `LiveCallRefreshHook` |
| `mutation_count` | `u32` | Yes — starts at parent's value (typically 0) | Running SHM counter; only ever non-zero post-fork in canonical pipelines |

### What therefore **is** family-invariant by construction

Because all of the above are cloned to every descendant at fork
time, every descendant starts with the **identical** V/D/J
assignments, trim values, NP bytes, junction germline coordinates,
receptor-revision provenance, and D-inversion flag. The user's
itemised list maps directly:

- **`assignments`** — fully shared. Every descendant has the same
  `assignments.{v,d,j,c}` truth ids.
- **`trims`** — fully shared. `assignments.{v,d,j}.{trim_5,trim_3}`
  identical across descendants of a clone.
- **`sequence regions`** — fully shared *at fork time*. Region
  coordinates and contents are identical; post-fork passes can
  modify them per descendant (SHM rewrites bytes in place, indels
  shift coordinates, end-loss truncates regions).
- **`junction`** — fully shared *pre-SHM*. The post-SHM junction
  (the `junction` AIRR field) is per-descendant because each
  descendant's SHM independently rewrites junction bytes; only
  the **pre-SHM** junction is invariant, and today the pre-SHM
  junction is **not surfaced anywhere** (gap pinned in the
  family audit).
- **`D inversion`** — fully shared. `assignments.d.orientation`
  is set pre-fork; the post-fork plan can't (and doesn't) modify
  it.
- **`receptor revision`** — fully shared. Receptor revision is
  inlined into `recombine` and runs pre-fork; every descendant
  has the same `assignments.v.original_allele_id`.
- **`paired-end layout`** — **not** shared in canonical
  pipelines: `.paired_end(...)` is a post-fork step and produces
  per-descendant r1/r2 cuts. If a user manages to put
  `.paired_end(...)` **before** `.expand_clones(...)`, the layout
  decisions get baked into the parent `Simulation` and every
  descendant inherits the same r1/r2 layout — almost certainly
  not what the user wanted. See §8 edge case 4.

### What is structurally **not** shared at the Rust level

The parent's `Outcome` fields **not in `Simulation`** are
structurally inaccessible to descendants:

- `Outcome.revisions: Vec<Simulation>` — only `revisions.last()`
  crosses the boundary (as `initial`); the per-pass history is
  dropped.
- `Outcome.pass_names: Vec<String>` — dropped.
- `Outcome.trace: Trace` — dropped. The descendant gets
  `Trace::new()`.
- `Outcome.events: Vec<EventRecord>` — dropped.

This is by design — `Simulation` is the state-only carrier; the
event-and-trace ledger lives on `Outcome` — but the descendant has
no way to refer back to the parent's events even if it wanted to.

### Gap pin

`pin_scaffold_simulation_carries_shared_recombination_state` —
the descendant's `revisions[0]` (its initial IR, accessible via
`revision(0)`) equals the parent's `final_simulation()` field-by-
field for the shared fields above. Pinning this protects the
invariant that the boundary handoff is *whole IR*, not a partial
subset.

`pin_absence_parent_trace_and_events_dropped_at_fork` — the
descendant's `outcome.event_count()` reflects only post-fork
passes; the parent's events do not appear. The descendant's
`outcome.trace()` length is strictly less than a non-clonal run
of the same overall pipeline (`recombine + mutate`) on a fresh
initial IR.

---

## 3. Q3 — What is descendant-specific?

The dual of §2. Per-descendant divergence is driven entirely by
the post-fork plan running with its own seed.

### Post-fork passes that produce per-descendant divergence

| Pass family | DSL surface | What diverges per descendant |
|---|---|---|
| SHM substitutions | `.mutate(rate=…)` / `.mutate(count=…)` | Number of substitutions, sites, replacement bases |
| PCR errors | `.pcr_amplify(count=…)` | PCR substitution sites + identities |
| Sequencing-quality errors | `.sequencing_errors(...)` | Quality-error sites |
| Polymerase indels | `.polymerase_indels(...)` | Indel sites, lengths, insertion/deletion choice |
| Ambiguous `N` corruption | `.ambiguous_base_calls(...)` | N-corruption positions |
| End-loss / strand orientation | `.corrupt_*` family | End-loss start/end clips, RC flip |
| Paired-end read layout | `.paired_end(...)` | r1/r2 lengths, insert-size sample |

### Ordering edge — `random_strand_orientation`

`random_strand_orientation` exists in the corruption family.
It's a sampling pass: it flips the strand of the assembled
sequence with some probability. Whether the flip is shared
across a clone or per-descendant depends entirely on whether
the DSL step is placed before or after `expand_clones`. Today's
`mcp_server.py:290-293` comment names the canonical order as
"recombine → restrict_alleles → expand_clones → mutate"; strand
orientation is not in that list, which makes its canonical
position ambiguous. A user who wants "every descendant has the
same strand orientation as the family ancestor" puts it
**pre-fork**; a user who wants "each descendant independently
flipped" puts it **post-fork**. The DSL doesn't surface a
preferred default — see §8 edge case 4.

### Gap pin

`pin_scaffold_descendant_seqs_differ_under_post_fork_passes` —
under any post-fork pass that introduces stochasticity (SHM,
PCR, indels), the descendants of a clone produce distinct
`sequence` strings. This is already pinned by the existing
[`test_g5_clonal_descendants_share_junction_diverge_via_mutate`](../tests/test_experiment.py#L1815);
we re-pin here at the audit-doc layer for §3 traceability.

`pin_absence_no_per_descendant_shm_event_aggregation` — there
is no API surfacing "for each descendant of clone C, the list of
SHM substitution sites." A consumer can read each descendant's
`outcome.events()` and filter for mutate events, but there is no
family-grouped view.

---

## 4. Q4 — What object should eventually exist?

This is the architectural design space. The audit names the
candidate shapes, weighs each, and recommends the one a future
implementation slice should aim at. **No commitment is made in
this slice;** the recommendation is for the next-slice author.

### Candidate A — `Outcome.parent_trace` (additive, minimal)

Add an `Optional<Trace>` field to `Outcome` that the orchestration
loop populates with the parent's trace before each descendant is
emitted. Symmetrically, `Outcome.parent_revisions` /
`Outcome.parent_events` if needed.

**Pros:**
- Smallest API surface change.
- Backwards-compatible — existing consumers ignore the new field.
- Zero new types.

**Cons:**
- **Massive duplication.** Every descendant in a clone carries a
  full copy of the parent's trace; for a 1k-descendant clone the
  parent trace is held 1000 times.
- Conceptually wrong: the parent isn't a per-descendant property;
  it's a clone-level property. Putting it on `Outcome` overstates
  its scope.
- Family-level operations (e.g. "give me the parent's projected
  AIRR record once") still need a group-by-parent step.

**Verdict:** rejected as the long-term shape, but acceptable as
a **Slice 2 read-only stop-gap** if the parent is also stored
once and pointed at via `Arc`. The duplication concern goes
away with `Arc<Trace>` / `Arc<Outcome>` sharing.

### Candidate B — `ParentOutcome` exposed via `SimulationResult.parents`

Materialize the parent `Outcome` once per clone, expose them as
a parallel list. Each descendant gains a `parent_id: int` field
indexing into `result.parents`.

**Shape:**
```python
class SimulationResult:
    @property
    def parents(self) -> Optional[List["Outcome"]]: ...   # None if non-clonal
    # `outcomes` unchanged (one per descendant)
    # `records[i]["clone_id"]` continues to identify the family
```

```rust
#[pyclass]
struct PyOutcome {
    // existing fields...
    parent_id: Option<u32>,   // index into SimulationResult.parents
}
```

**Pros:**
- Each parent stored exactly once — no duplication.
- Symmetric with the existing `result.outcomes` flat list:
  callers who want the parent for a record do
  `result.parents[record["clone_id"]]`.
- Doesn't introduce a new aggregate type; reuses `Outcome`.
- The parent's full `revisions` + `trace` + `events` are
  accessible for replay / debugging / lineage analysis.
- No biology coupling — purely a provenance/replay surface.

**Cons:**
- Adds a Python-side concept (`SimulationResult.parents`) that
  the Rust side doesn't model. The wiring lives in
  `CompiledClonalExperiment.run_records`.
- A future `FamilyRecord` / `FamilyOutcome` aggregate (Candidate
  D) would supersede this with a cleaner structure.

**Verdict:** **recommended** as the Slice 2 shape. It is the
minimum-viable read-only surface that:
- Exposes the ancestor without duplicating it.
- Lets `validate_families` (or a future `validate_families_with_parent`)
  compare descendants against the parent's truth fields, the
  parent's pre-SHM junction (projected once), and the parent's
  receptor-revision provenance.
- Lets users replay just the parent half of a clone (call
  `compiled._pre.replay_from_trace_file(...)` against the
  parent's trace).

### Candidate C — `ClonalFamily` aggregate type

Introduce a new Python type bundling `(parent_outcome,
descendant_outcomes, family_metadata)`. `SimulationResult.families:
List[ClonalFamily]` becomes the primary clonal view; the existing
`outcomes` / `records` flat lists are derived views.

**Pros:**
- Clean aggregate matches the biology (one family = one
  recombination ancestor + K descendants).
- Validator code reads naturally: `for family in result.families: ...`
- Future biology (selection, time evolution) hangs off
  `ClonalFamily` instead of stacking `clone_id`-keyed dicts.

**Cons:**
- Larger surface change — two new types (`ClonalFamily` +
  `FamilyMetadata`) plus the `families` accessor.
- Risks fragmenting "I want the flat record list" callers from
  "I want the family view" callers; needs careful design to
  keep flat-list access cheap.
- More disruption to existing tests / consumers; the
  `SimulationResult.outcomes` flat list semantics needs explicit
  preservation.

**Verdict:** **recommended as Slice 3+**, after Candidate B
ships. Slice 2 (Candidate B) makes parents addressable; Slice 3+
(Candidate C) builds the aggregate type biology will hang off of.

### Candidate D — `FamilyRecord` AIRR projection

Project the parent IR to an AIRR record dict — a `FamilyRecord`
that carries truth alleles, pre-SHM junction, pre-SHM junction_aa,
CDR3 germline coords, etc. — and expose it alongside descendant
records.

**Pros:**
- Lets downstream consumers (AIRR-strict TSV writers, clonotype
  clusterers, etc.) treat the family ancestor as a first-class
  AIRR row.
- Closes the family audit's §3 gap on pre-SHM junction.

**Cons:**
- The parent IR has no observed sequence (assembly is "germline-
  frame" — the assembled sequence is the germline-faithful
  V+NP1+D+NP2+J string). Several AIRR fields (`productive`,
  `sequence_alignment`, mutation counts) are ill-defined for a
  family ancestor.
- Needs a separate projection path that handles "no observed
  read" semantics, distinct from `outcome_to_airr_record`.

**Verdict:** **deferred to Slice 4+**. Builds on top of Candidate
B/C and requires its own design pass.

### Candidate E — `Outcome.parent_id` only (no parent storage)

Just add a `parent_id: int` field on descendant `Outcome`s; the
parent itself is not retained. This is purely a marker.

**Pros:**
- Zero memory overhead.
- Cleanly groups descendants without committing to a storage shape.

**Cons:**
- The `parent_id` integer doesn't point at anything — it's a
  re-encoding of `clone_id`. No new information.
- Doesn't help mutation-distance / pre-SHM-junction checks; the
  parent IR is still unreachable.

**Verdict:** **rejected**. Adds noise without providing the
parent itself; equivalent to today's `clone_id` tag.

### Recommendation summary

| Slice | Shape | Why |
|---|---|---|
| Slice 2 (immediate next) | **Candidate B**: `SimulationResult.parents` + `Outcome.parent_id` | Minimum-viable read-only ancestor surface. No new aggregate type, no biology coupling. |
| Slice 3 | **Candidate C**: `ClonalFamily` aggregate | Biology-friendly view. Builds on Slice 2. |
| Slice 4 | **Candidate D**: `FamilyRecord` projection | Pre-SHM junction + truth alleles as a first-class AIRR row. |
| Rejected | Candidates A and E | A duplicates, E adds noise. |

### Gap pin

`pin_absence_no_parents_attr_on_simulation_result` — Slice 2 not
landed yet. `SimulationResult` exposes only `.records` and
`.outcomes`.

`pin_absence_no_parent_id_field_on_outcome` — Slice 2 back-pointer
not landed yet.

(These are sharper restatements of pins the family audit already
flagged; this audit's contract file restates them for the
parent-outcome design lockstep.)

---

## 5. Q5 — Replay model

The user's question: "parent trace + child traces? one flattened
trace with branch markers? trace file schema implications?"

### Today

A single `Outcome` carries a single `Trace`. The descendant's
trace contains only post-fork sampling sites; the parent's
trace contains only pre-fork sampling sites. There is **no
flattened trace with branch markers** — branches live in the
Python orchestration loop, not in the trace format.

`CompiledSimulator.trace_file_from` ([`engine_rs/src/python/compiled.rs:197`](../engine_rs/src/python/compiled.rs#L197))
bundles `(plan_signature, refdata_signature, seed, trace)` for
**one** outcome's simulator. A clonal-batch trace file would
need either:

- Two plan signatures (`pre` + `post`) + two seeds
  (`clone_seed` + `desc_seed`) + two traces — replay the
  parent + replay the descendant from the parent's IR.
- Or a flattened trace where the post-fork sampling sites are
  appended to the pre-fork sites with a `fork()` marker between,
  and `replay_from_trace_file` knows to switch plans at the
  marker.

### Recommendation: parent trace + child trace, no flattening

The **two-trace** model maps cleanly to today's data flow:
- Parent trace = `pre.run(clone_seed).trace()`
- Descendant trace = `post.run_from(parent_sim, desc_seed).trace()`

A clonal trace file (Slice 4+ scope) would bundle:
```
{
  "pre_plan_signature": "...",
  "pre_refdata_signature": "...",
  "clone_seed": <u64>,
  "parent_trace": <Trace>,
  "post_plan_signature": "...",
  "desc_seed": <u64>,
  "descendant_trace": <Trace>,
}
```

Replay: `pre.replay_from_trace(parent_trace, clone_seed)` →
`parent_sim` → `post.replay_from_trace_with_initial(descendant_trace,
parent_sim, desc_seed)`.

A **flattened trace with branch markers** is structurally
appealing (single bundle per descendant) but introduces a new
trace-format concept (`fork()` markers) and forces the replay
machinery to switch plans mid-trace — neither benefit justifies
the cost. Reject flattening.

### Per-clone vs per-descendant trace file scope

The per-clone half (parent trace) is replayable once for the
whole family; the per-descendant halves are independent. A
trace file that wants to round-trip "give me family C" would
bundle one parent trace + K descendant traces; "give me
descendant (C, D)" would bundle one parent trace + one descendant
trace. The format is the same; only the descendant-trace count
varies.

This is out of scope for any near-term slice but inheriting the
two-trace model now (Slice 2's parent-outcome accessor) makes it
trivial to add later — the parent's `outcome.trace()` is already
there to bundle.

### Gap pin

`pin_absence_no_clonal_trace_file_format` — neither
`compiled._pre.trace_file_from` nor `compiled._post.trace_file_from`
exposes a clonal-aware bundling method. Trace files today are
per-`Outcome` and don't encode the parent → descendant relation.

---

## 6. Q6 — Validator model

The family audit's Slice 1 (`validate_families`) is field-level:
it groups records by `clone_id` and compares `truth_v_call` /
`truth_d_call` / `truth_j_call` for invariance. **It can't compare
records against the parent**, because the parent isn't there.

### What a parent-aware validator unlocks

Once Slice 2 (Candidate B) lands, the validator surface grows to:

1. **Junction divergence detection.** The current `validate_families`
   can't enforce junction invariance because the descendant's
   `junction` field is post-SHM. With the parent IR accessible, the
   validator can extract the **pre-SHM junction** from the parent's
   `Simulation` (the bytes between CDR3 anchors at fork time) and
   assert every descendant's `junction` is consistent with it modulo
   SHM substitutions. Concretely: descendant `junction` length must
   equal parent pre-SHM junction length (SHM doesn't insert/delete);
   per-base divergence must match the descendant's SHM event count.

2. **Mutation-distance distribution.** For each clone, compute
   `hamming(descendant.sequence, parent.sequence)` against the
   parent's projected sequence. Assert the distribution is
   consistent with the configured mutate rate / count. Today this
   check is unrunnable.

3. **`original_v_call` invariance.** Receptor-revision provenance
   lives in the parent's `assignments.v.original_allele_id`.
   Already projected onto each descendant by
   `outcome_to_airr_record`, but a parent-aware validator can
   double-check "every descendant agrees with the parent's value"
   instead of just "descendants agree with each other."

4. **`d_inverted` invariance.** Same as `original_v_call`.

5. **Truth-allele invariance with provenance.** Today's
   `validate_families` compares descendant truth alleles to each
   other. A parent-aware validator can compare them to the parent
   directly — useful when a future bug widens descendants'
   `assignments` past what the parent set.

### Sketch — `validate_families_with_parents()`

A future Slice 2+ surface would look like:

```python
class SimulationResult:
    def validate_families_with_parents(
        self, refdata: Optional[Any] = None
    ) -> "FamilyValidationReport":
        """Parent-aware family validation. Requires the result
        carries `.parents` (Slice 2). Raises if .parents is None.
        """
```

The existing `validate_families()` stays as the dict-only fallback
for results-from-TSV; the parent-aware version is the release-tier
gate. `validate_records=True` on a clonal `run_records` would
prefer the parent-aware version when parents are available.

### What this audit does *not* commit to

The validator model above is sketched, not specified. Slice 2
(parent accessor) ships **without** changing the validator; Slice
3+ adds the parent-aware checks. Splitting the surfaces keeps
each slice diffable and reviewable.

### Gap pin

`pin_absence_no_parent_aware_validator` —
`validate_families_with_parents` does not exist. The field-only
`validate_families` (Slice 1) is the only family-layer gate.

---

## 7. Q7 — Implementation-readiness recommendation

The user's prompt: "First implementation slice after audit should
be read-only: expose parent metadata/trace without changing
simulation semantics. Avoid changing branching behavior until we
can validate it."

The audit endorses this. Concretely:

### Recommended Slice 2 — Parent-Outcome Read-Only Surface

**Scope:**
1. Modify `CompiledClonalExperiment.run_records` to retain the
   parent `Outcome` per clone in a `parents: List[Outcome]` list.
2. Add `SimulationResult.parents: Optional[List[Outcome]]`
   property — `None` for non-clonal results, populated for
   clonal.
3. Stamp `parent_id: int` on each descendant's projected AIRR
   record (mirroring the existing `clone_id` tag — `parent_id ==
   clone_id` today; we keep both because `parent_id` is the
   indexing semantic and `clone_id` is the family-identity
   semantic).
4. Reserve but **do not implement** `Outcome.parent_id` on the
   Rust side — the back-pointer needs more design (see §11
   Backwards compatibility).

**Out of scope:**
- No validator changes.
- No new aggregate type (`ClonalFamily`, `FamilyRecord`).
- No trace-file format changes.
- No simulation-semantics changes — the orchestration loop
  produces byte-identical descendant outcomes; only the parent
  retention is new.

**Why this scope:**
- Pure additive — `SimulationResult.parents` defaults to `None`
  for everything that doesn't opt in.
- No biology coupling — the parent is exposed exactly as it
  exists today, just with an extra reference keeping it alive.
- Byte-reproducibility of descendants stays guaranteed
  (orchestration produces the same outcomes; we just don't
  drop the parent at the end of the iteration).

### Why **not** Slice 3 first

`ClonalFamily` would be the cleaner long-term shape but:
- Larger diff, more surface area to review.
- Forces a decision on whether `outcomes` / `records` stay flat
  or get nested under families — that decision should be made
  after a few weeks of using Slice 2 in practice, not upfront.
- Slice 2's `parents` list is trivially upgradeable to "the
  `parent` field of a `ClonalFamily`" later; the reverse is
  harder.

### Why **not** validator changes in Slice 2

Each slice should land **one** kind of change:
- Slice 2 = parent accessor (data flow).
- Slice 3+ = parent-aware validator (semantics).

Combining them risks the validator finding a parent-data bug and
the slice author conflating "the data flow is wrong" with "the
validator is too strict." Splitting keeps the diff causal.

---

## 8. Edge cases the implementation slice must handle

1. **Non-clonal experiment.** `SimulationResult.parents` returns
   `None`. Existing per-record APIs unchanged. `validate_families`
   continues to return ok-no-op (already covered by Slice 1).

2. **Empty post-fork plan.** `expand_clones()` immediately
   followed by `run_records()` (no `.mutate` / `.pcr_amplify` /
   etc.): every descendant is byte-identical to the parent. The
   parent retention is still valid (the parent IR exists, even
   when descendants equal it). Slice 2 must not assume "if
   parents == descendants we don't need parents."

3. **Productive-only with heavy SHM.** Under `productive_only()`
   + high mutate, the post-fork pass can repeatedly reject
   descendants. The total record count guarantee
   (`n_clones * per_clone`) is structural; the parent retention
   guarantee must not depend on the descendant count.

4. **Steps placed pre-fork that should be post-fork (and vice
   versa).** The DSL accepts:
   ```python
   exp.recombine().mutate(count=10).expand_clones(...)   # mutate pre-fork
   exp.recombine().expand_clones(...).receptor_revision(...)  # revision post-fork
   ```
   Neither raises today. The first produces "shared SHM across
   all descendants of a clone" (biologically wrong); the second
   produces "per-descendant receptor revision" (biologically
   wrong; receptor revision is a B-cell recombination-time
   decision, one per ancestor). **Slice 2 does not fix this.**
   The audit flags it as a separate DSL-validation gap that
   should ship with the validator slice (Slice 3+), not the
   read-only accessor slice.

5. **`random_strand_orientation` ordering.** Same shape as #4 —
   pre-fork vs post-fork placement changes biology. Out of scope
   for Slice 2.

6. **Paired-end + clonal.** When `.paired_end(...)` is in the
   post-fork plan (canonical), each descendant gets its own r1/r2
   layout. Slice 2 retains the parent IR (pre-paired-end) and the
   K descendant outcomes (post-paired-end); there is no separate
   paired-end-aware parent representation needed.

7. **TCR + clonal.** TCR loci accept `expand_clones`; SHM is
   rejected on TCR so `mutate` can't be in the post-fork plan
   for TCR. Slice 2 is locus-agnostic; the parent retention
   logic doesn't branch on locus.

8. **Strict mode for parent vs descendant runs.** The
   orchestration loop forwards `strict=strict` to **both**
   `pre.run(clone_seed, strict=strict)` and `post.run_from(...,
   strict=strict)`. If the parent strict-errors, the loop aborts
   that whole clone (no descendants emitted); if a descendant
   strict-errors, the loop aborts the whole batch (no further
   descendants from that clone or others). Slice 2 must preserve
   this failure ordering — a parent strict-error must not be
   caught + papered over to keep emitting descendants.

---

## 9. Replay determinism

(Section structure mirrors the family audit §9 for
cross-doc traceability.)

### Per-descendant replay — sufficient today

Given `(pre, post, refdata, seed)`, the descendant
`(clone_idx, desc_idx)` reproduces byte-for-byte. The
descendant trace + descendant seed are enough to replay the
descendant **from the parent IR**, but they are NOT enough to
reconstruct the parent IR — that requires re-running the parent
with `clone_seed`.

### Per-family replay — two-pass, manual today

Today's replay path: `compiled._pre.replay_from_trace_file(parent_tf)`
gives back the parent outcome → `compiled._post.run_from(parent_sim,
desc_seed)` re-derives descendant D. **But:**
- The parent trace file isn't produced today (no orchestration
  step calls `compiled._pre.trace_file_from(parent, clone_seed)`).
- A consumer who wants to do this must call `compiled._pre.run(clone_seed)`
  themselves and then `compiled._pre.trace_file_from(...)`.

Slice 2 fixes the *first* gap by retaining the parent on the
result. Slice 4+ would fix the *second* gap by adding clonal
trace-file bundling.

### Cross-version replay

The seed-derivation formula `clone_seed = seed + clone_idx *
1_000_000` is internal but visible by behavioural side-effect.
A future Slice 2 (or later) should consider exposing it as
`CompiledClonalExperiment.clone_seed_for(seed, clone_idx)` so
consumers don't reinvent it. This audit doesn't commit to the
helper; the existing pin
(`pin_scaffold_clone_seed_formula_is_stable` in the family
audit's contract file) covers the current scheme.

### Gap pin

`pin_absence_no_parent_trace_file_emission` — the clonal
orchestration doesn't emit parent trace files. Slice 2 doesn't
need to fix this; Slice 4+ does.

---

## 10. Validator integration

(Slice 2 ships **no validator changes**; this section sketches
what becomes possible after Slice 2 + 3.)

### Slice 2 (read-only accessor) — no validator hook

After Slice 2:
- `SimulationResult.parents` is populated for clonal runs.
- `validate_families()` continues to be field-only — it doesn't
  read `.parents`.
- `validate_records=True` continues to run per-record gate +
  field-only family gate, in that order.

### Slice 3 (parent-aware validator) — new hook

After Slice 3:
- New `validate_families_with_parents(refdata)` returns the
  same `FamilyValidationReport` shape but adds parent-aware
  issue kinds: `JunctionDivergesFromParentPreSHM`,
  `MutationDistanceOutOfSpec`, `OriginalVCallDivergesFromParent`,
  `DInvertedDivergesFromParent`.
- `validate_records=True` calls
  `validate_families_with_parents` instead of
  `validate_families` when `.parents` is populated.
- Field-only `validate_families` stays as a fallback for results
  loaded from TSV (no parents available).

### Cache parity interaction

`Outcome.check_live_call_cache_parity` is per-outcome and runs
on each descendant outcome independently. Slice 2 doesn't change
that — the new `.parents` entries are also `Outcome`s, but cache
parity isn't asked of the parent in the standard release
workflow today (the parent's live-call state is mostly
"germline" — V/D/J freshly assembled and untouched by SHM).
A future Slice 3+ could add `parent.check_live_call_cache_parity()`
to the release gate; out of scope here.

---

## 11. Backwards compatibility

Slice 2 must not break:

- `SimulationResult.records` — same flat list, same fields.
- `SimulationResult.outcomes` — same flat list of descendant
  outcomes. **The parent outcomes do NOT appear here**;
  they live on the separate `.parents` list. Conflating them
  would break any consumer who does `len(result.outcomes) ==
  n_clones * per_clone`.
- `Outcome` Rust type — no new required fields. Adding
  `parent_id: Option<u32>` is acceptable as an optional field
  but should default to `None` and only be set on descendant
  outcomes in a clonal run; non-clonal outcomes keep
  `parent_id = None`. **Recommendation: defer the Rust-side
  field to Slice 3+; Slice 2 lives entirely in Python with the
  `parent_id` written onto the AIRR record dict only.**
- `validate_records=True` semantics — unchanged. Slice 2 doesn't
  touch the validator hook.
- `clone_id` integer tag — unchanged; remains in `[0, n_clones)`.
- All `pin_scaffold_*` tests in this audit + the family audit
  must continue to pass.

The one place a careful contributor must watch: extending
`SimulationResult` with `__slots__ = (..., "_parents", ...)`. The
current slots are `("_records", "_outcomes")`. Adding a slot is
a transparent change for consumers but is a breaking change for
any subclass that defined its own slots — we don't have any in
the codebase today, but the pin
`pin_scaffold_simulationresult_slots_documented_for_extension`
will flag the addition so reviewers notice.

### Gap pin

`pin_scaffold_simulationresult_slots_documented_for_extension` —
the current slots are documented and any new slot must be
recorded.

---

## 12. Implementation order

Same shape as the family audit's §12, refined for the parent
accessor.

1. **Slice 2 — Parent-Outcome Read-Only Surface** (recommended
   next). Python-only:
   - Retain parents in `CompiledClonalExperiment.run_records`.
   - Expose `SimulationResult.parents`.
   - Stamp `parent_id: int` (= `clone_id`) on each AIRR record.
   - No validator changes, no Rust changes, no aggregate type.
   - **Closes the audit's main architectural gap** ("there is no
     addressable parent today") in a minimal-diff slice.

2. **Slice 3 — Parent-Aware Validator.** Python-only:
   - Add `validate_families_with_parents(refdata)`.
   - Wire it into `validate_records=True` for clonal runs (taking
     precedence over the field-only `validate_families` when
     `.parents` is populated).
   - Add the parent-aware issue kinds enumerated in §6.

3. **Slice 4 — `ClonalFamily` aggregate type.** Larger surface:
   - New `ClonalFamily` Python type bundling `(parent, descendants,
     metadata)`.
   - `SimulationResult.families: List[ClonalFamily]`.
   - Migrate the validator to iterate families.
   - The flat `records` / `outcomes` lists become derived views.

4. **Slice 5 — Clonal trace files.** Extend `trace_file_from` for
   the clonal path; bundle parent + descendant traces in a
   single file. Out of scope until per-family replay demand
   surfaces.

5. **Slice 6+ — Biology on top.** Lineage trees, germinal-center
   selection, antigen-specific weighting. All compose against
   Slice 4's `ClonalFamily`.

---

## 13. Test surface — what this audit pins

Mirrored in [`tests/test_clonal_parent_contract.py`](../tests/test_clonal_parent_contract.py).

### `pin_scaffold_*` — existing parent-lifecycle surfaces

1. `pin_scaffold_outcome_struct_shape_today` —
   `Outcome.revisions` / `pass_names` / `trace` / `events` are
   the four fields; no `parent_id` / `parent_trace` field.
2. `pin_scaffold_simulation_struct_shape_today` — `Simulation`
   carries `pool` / `sequence` / `assignments` / `segment_calls`
   / `dirty_log` / `mutation_count`; no trace, no events.
3. `pin_scaffold_run_one_from_with_policy_takes_initial_sim` —
   the Rust signature: `run_one_from_with_policy(initial:
   Simulation, seed, policy) -> Outcome`.
4. `pin_scaffold_execute_transactional_allocates_fresh_trace` —
   the descendant trace is allocated at the top of
   `execute_transactional` (regex on the source: `let mut trace
   = Trace::new()`).
5. `pin_scaffold_orchestration_uses_pre_run_then_post_run_from` —
   the Python loop signature inside
   `CompiledClonalExperiment.run_records` (regex on the source).
6. `pin_scaffold_simulation_carries_shared_recombination_state` —
   the descendant's `revision(0)` equals the parent's
   `final_simulation()` on the shared fields (assignment ids,
   trim values, junction germline coords).
7. `pin_scaffold_descendant_seqs_differ_under_post_fork_passes` —
   re-pin of the existing `test_g5_*` invariant for §3
   traceability.
8. `pin_scaffold_clonal_truth_calls_stable_within_clone_under_normal_fixtures` —
   under the standard IGH fixture (no extreme SHM), truth alleles
   within a clone are byte-identical. Pin the "audit baseline"
   behaviour.
9. `pin_scaffold_simulationresult_slots_documented_for_extension` —
   the current slots are exactly `("_records", "_outcomes")`.
   Adding `_parents` must be recorded in lockstep.
10. `pin_scaffold_validate_families_today_is_field_only` —
    `validate_families` doesn't reference `.parents` /
    `parent_id`; the Slice 1 surface stays narrow.

### `pin_absence_*` — gaps Slice 2+ closes

1. `pin_absence_no_parents_attr_on_simulation_result` —
   `SimulationResult` has no `.parents` attribute.
2. `pin_absence_no_parent_id_field_on_outcome` — `Outcome` has no
   `parent_id` field (Rust or Python).
3. `pin_absence_no_parent_id_on_airr_records` — record dicts
   carry `clone_id` but not `parent_id`.
4. `pin_absence_parent_outcome_unreachable_from_descendant` —
   given a descendant `Outcome`, there is no API path to recover
   the parent `Outcome` / parent `Trace` / parent revisions.
5. `pin_absence_parent_trace_and_events_dropped_at_fork` — the
   descendant's `outcome.event_count()` is strictly less than a
   non-clonal run of `recombine + same post-fork plan` from a
   fresh empty IR (the difference equals the recombine-pass
   events the parent ran).
6. `pin_absence_no_parent_aware_validator` —
   `validate_families_with_parents` does not exist.
7. `pin_absence_no_clonal_trace_file_format` —
   `trace_file_from` is per-`Outcome`; no clonal-aware bundler.
8. `pin_absence_no_per_descendant_shm_event_aggregation` — no
   API surfaces a per-clone or per-family view of SHM events.

### Flip behaviour when slices land

When **Slice 2** lands:
- `pin_absence_no_parents_attr_on_simulation_result` →
  `pin_present_parents_attr_on_clonal_result`.
- `pin_absence_no_parent_id_on_airr_records` →
  `pin_present_parent_id_on_descendant_records`.
- `pin_absence_parent_outcome_unreachable_from_descendant`
  partially flips (parent reachable via `result.parents[r["clone_id"]]`).
- `pin_scaffold_simulationresult_slots_documented_for_extension`
  gets updated to add `_parents` to the documented slot list.

When **Slice 3** lands:
- `pin_absence_no_parent_aware_validator` →
  `pin_present_parent_aware_validator`.
- `pin_scaffold_validate_families_today_is_field_only`
  gets updated to note both paths exist.

When **Slice 4** lands (clonal trace files):
- `pin_absence_no_clonal_trace_file_format` →
  `pin_present_clonal_trace_file_format`.

The `pin_scaffold_*` tests on `Outcome` / `Simulation` shape stay
unchanged unless those underlying types change — which Slice 2
deliberately does not.

---

## 14. Out of scope

Documented here so a future Slice 2 author doesn't accidentally
expand the implementation.

- **No simulation-semantics changes.** Slice 2 retains the
  parent; it does not change how parents or descendants are
  produced. Byte-reproducibility of every descendant must hold.
- **No Rust-side `Outcome.parent_id` field.** Defer to Slice 3+.
  Slice 2 lives entirely in Python; the back-pointer goes on
  the AIRR record dict (`record["parent_id"]`) — which is
  trivially the same value as `record["clone_id"]` today but
  carries a different semantic.
- **No new aggregate type.** `ClonalFamily` is Slice 4+.
- **No validator changes.** Field-only `validate_families`
  remains the only family-layer gate after Slice 2.
- **No trace-file format change.** Clonal trace files are Slice
  5+.
- **No DSL ordering enforcement.** The pre-fork / post-fork step
  ordering gaps (§8 #4, #5) are flagged but not fixed.
- **No new `clone_seed_for(seed, clone_idx)` helper.** Convenient
  but out of scope; the formula stays internal-but-stable.
- **No biology layered on parents.** Selection, time evolution,
  lineage topology — all Slice 6+.

---

## Summary table

| Question | Answer |
|---|---|
| Where is the ancestral recombination generated today? | In `CompiledClonalExperiment.run_records`, exactly one `pre.run(clone_seed)` per `clone_idx`. The parent `Outcome` is a Python local; the orchestration loop calls `.final_simulation()` to extract the IR and drops the parent at the next iteration. The parent trace is never copied, stored, or projected. |
| What is actually shared? | The full `Simulation` IR: `pool`, `sequence`, `assignments`, `segment_calls`, `dirty_log`, `mutation_count`. By extension: truth alleles, trims, NP bytes, junction germline coords, `d_inverted`, `original_v_call`. **Not** shared: trace, events, revision history. |
| What is descendant-specific? | Everything the post-fork plan introduces: SHM substitutions, PCR errors, sequencing errors, indels, N corruption, end-loss, paired-end r1/r2, random strand orientation (depending on pre-fork vs post-fork placement). |
| What object should eventually exist? | **Slice 2: `SimulationResult.parents: Optional[List[Outcome]]`** + `parent_id` on AIRR records. Slice 3: `validate_families_with_parents`. Slice 4: `ClonalFamily` aggregate. Slice 5: `FamilyRecord` projection. Reject `Outcome.parent_trace` (duplicates) and `Outcome.parent_id` alone (no new info). |
| Replay model? | **Parent trace + child trace, no flattening.** A clonal trace file (Slice 5+) bundles two plan signatures + two seeds + two traces. Rejected: flattened trace with branch markers. |
| Validator model? | Slice 1's `validate_families` stays field-only. Slice 3 adds `validate_families_with_parents` for pre-SHM junction, mutation-distance distribution, parent-truth comparison. Splitting the validators keeps the field-only path available for results-from-TSV. |
| First implementation slice after this audit? | **Slice 2 — Parent-Outcome Read-Only Surface**, Python-only, no Rust changes, no validator changes, no aggregate type. Retain parents + expose `.parents` + stamp `parent_id` on records. Smallest diff that closes the addressability gap. |

The audit's verdict: today's clonal path produces a real ancestor
and then throws it away. Closing the gap is one Python-only
slice that retains and exposes what already exists — no new
simulation semantics, no new types, no biology coupling.
