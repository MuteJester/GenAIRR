# Clonal Family / Lineage — Pre-Implementation Audit

**Status: audit only.** Pins today's behaviour and the architecture
gaps a future "real lineage" slice will have to close. No
implementation is proposed in this pass; the deliverable is the
shared vocabulary + the contract pins so a later slice can land
against an audited baseline.

Companion to
[`tests/test_clonal_family_contract.py`](../tests/test_clonal_family_contract.py)
which freezes today's behaviour (`pin_scaffold_*`) and the gaps the
audit identifies (`pin_absence_*`). When an implementation slice
lands, the `pin_absence_*` tests flip to `pin_*_present` /
`pin_*_used` in lockstep with the slice.

The clonal-expansion surface — `expand_clones(n_clones, per_clone)`
+ `CompiledClonalExperiment` + the `clone_id` tag — has shipped and
the existing G5 tests
([`tests/test_experiment.py:1784-1891`](../tests/test_experiment.py#L1784-L1891))
already pin: total record count, "identical without post-fork
passes", "share junction backbone, diverge under SHM", double-call
rejection, explicit-`n` reconciliation. **This audit does not
relitigate that surface.** It asks the next-layer question: do the
records returned by clonal expansion model a *family*, or do they
model `n_clones × per_clone` independent records that happen to
share a parent IR by construction?

---

## Existing scaffolding the audit relies on

Surfaces today's engine provides, which any future implementation
slice will compose against. Each is pinned by a `pin_scaffold_*`
test in the companion contract file.

| Existing surface | Where | What it gives a family-layer audit |
|---|---|---|
| `_ClonalForkStep` marker | [`src/GenAIRR/_pipeline_ir.py:105`](../src/GenAIRR/_pipeline_ir.py#L105) | Frozen dataclass with `n_clones` + `size`. Pinned partitioning point in the step list — the audit's "pre-fork plan" / "post-fork plan" terminology is grounded here. |
| `Experiment.compile()` fork split | [`src/GenAIRR/experiment.py:1691-1742`](../src/GenAIRR/experiment.py#L1691-L1742) | Detects the fork, splits `self._steps`, builds two `CompiledSimulator`s, returns `CompiledClonalExperiment`. The "pre-fork = per-clone, post-fork = per-descendant" partitioning is structural, not a runtime branch. |
| `CompiledClonalExperiment.run_records` | [`src/GenAIRR/_compiled.py:393-453`](../src/GenAIRR/_compiled.py#L393-L453) | The parent→descendants orchestration loop: `pre.run(clone_seed)` produces a parent `Outcome`, `parent.final_simulation()` is shared with `post.run_from(parent_sim, desc_seed)` for each descendant. |
| `CompiledSimulator.run_from` (Rust) | [`engine_rs/src/python/compiled.rs:166-183`](../engine_rs/src/python/compiled.rs#L166-L183) | The Rust-side entry point that takes a parent `Simulation` IR and runs the post-fork plan from it. Docstring already names "clonal-family generation" semantics. |
| Post-fork contract bundle inheritance | [`engine_rs/src/compiled/analyze.rs:58-68`](../engine_rs/src/compiled/analyze.rs#L58-L68) and [`:207-221`](../engine_rs/src/compiled/analyze.rs#L207-L221) | `has_recombination_support()` lets the post-fork plan skip productive-frame compile-time preconditions (no recombination facts available) while keeping runtime contract filtering live. |
| Per-descendant `clone_id` tag | [`src/GenAIRR/_compiled.py:444`](../src/GenAIRR/_compiled.py#L444) | The single "family identity" surface visible to downstream consumers today: an integer in `[0, n_clones)` written onto each AIRR record dict. |
| `expose_provenance` truth columns | [`src/GenAIRR/_compiled.py:445-446`](../src/GenAIRR/_compiled.py#L445-L446) and `_inject_truth_columns` | Adds `truth_v_call/d_call/j_call`. Already verified family-invariant by [`test_g8_truth_columns_via_clonal_pipeline`](../tests/test_experiment.py#L1707). |
| `SimulationResult.validate_records` | [`src/GenAIRR/result.py:374-431`](../src/GenAIRR/result.py#L374-L431) | Per-record postcondition validator; iterates `(outcome, record)` pairs independently. No family-level pass exists. |
| `Outcome.trace` | [`engine_rs/src/python/outcome.rs:85-88`](../engine_rs/src/python/outcome.rs#L85-L88) | Addressed-choice trace for replay. Each descendant carries only its own post-fork trace; the pre-fork parent trace is discarded after `final_simulation()` extracts the IR. |

---

## 1. Current state — what shipping clonal expansion does today

`Experiment.expand_clones(n_clones=N, per_clone=K)` appends a
`_ClonalForkStep(n_clones=N, size=K)` marker to the step list.
`Experiment.compile()` detects the marker, splits the step list at
its index, builds a pre-fork `CompiledSimulator` and a post-fork
`CompiledSimulator`, and returns `CompiledClonalExperiment` instead
of `CompiledExperiment`.

`CompiledClonalExperiment.run_records(seed=S)` runs the orchestration
loop:

```python
for clone_idx in range(self._fork.n_clones):
    clone_seed = int(seed) + clone_idx * 1_000_000
    parent = self._pre.run(seed=clone_seed, strict=strict)
    parent_sim = parent.final_simulation()
    for desc_idx in range(self._fork.size):
        desc_seed = clone_seed + 1 + desc_idx
        desc = self._post.run_from(parent_sim, desc_seed, strict=strict)
        outcomes.append(desc)
        rec = outcome_to_airr_record(desc, self._refdata,
                                     sequence_id=f"clone{clone_idx}_desc{desc_idx}")
        rec["clone_id"] = clone_idx
        ...
```

The pre-fork `Outcome` is **discarded** after its
`final_simulation()` is extracted; the post-fork descendants are
each materialized as full independent `Outcome`s (their own trace,
their own event ledger, their own final `Simulation`). The returned
`SimulationResult` is a flat list of `n_clones * per_clone` records
+ `outcomes`, plus a `clone_id` int per record. **There is no family
object, no parent record, no parent outcome, no pre-fork trace
preservation, and no family-level validator.**

What ships today is a *deterministic clonal-expansion record
generator*. What it does *not* ship is a first-class family concept.

---

## 2. Q1 — What does the current clonal path actually share?

The user's first audit question, answered concretely against the
code above.

### Shared by construction (via `parent_sim`):

| Thing | Where it lives in `parent_sim` | Visible on descendant record? |
|---|---|---|
| V/D/J allele identity | `Simulation.assignments.{v,d,j}.allele_id` | Yes, as `v_call/d_call/j_call` (modulo SHM-driven tie growth) |
| `original_v_call` (receptor-revision provenance) | `Simulation.assignments.v.original_allele_id` | Yes, projected by `outcome_to_airr_record` |
| V/D/J trim amounts (`trim_5`, `trim_3`) | `Simulation.assignments.{v,d,j}.{trim_5,trim_3}` | Indirectly via coordinate fields |
| NP1 / NP2 bytes | `Simulation.pool[…]` for the np region | Indirectly via the assembled `sequence` (modulo SHM) |
| Junction (pre-SHM) | Pool bytes between `cdr3_start_pos` / `cdr3_end_pos` | Yes, as `junction` (post-SHM modifies it) |
| Receptor-revision applied flag | `Simulation.events.receptor_revision_*` events | Yes, as `original_v_call` divergence |
| D inversion (`d_inverted`) | `Simulation.assignments.d.orientation` | Yes, as the `d_inverted` AIRR field |
| Assembled IR (pool + region map) | `Simulation` IR itself | Yes, every descendant starts from this IR before post-fork passes |

### Not shared — diverges per descendant via post-fork passes:

| Thing | Driven by | Per-descendant divergence today |
|---|---|---|
| SHM substitutions (count, sites, identities) | `mutate(rate=…)` or `mutate(count=…)` | Independent per descendant |
| PCR errors | `pcr_amplify(count=…)` | Independent per descendant |
| Sequencing-quality errors | `sequencing_errors(...)` | Independent per descendant |
| Polymerase indels | `polymerase_indels(...)` | Independent per descendant |
| Ambiguous `N` bases | `ambiguous_base_calls(...)` | Independent per descendant |
| End-loss / random-strand-orientation | `corrupt_*` family | Independent per descendant |
| Paired-end r1/r2 fragmentation | `paired_end(...)` | Independent per descendant |
| `productive` / `complete_vdj` | Computed at projection time after post-fork passes | Can differ per descendant under heavy SHM |

### Lost between parent and descendants today:

| Thing | Where it existed | Why it's lost |
|---|---|---|
| Parent recombination trace | `parent.trace` (a `PyTrace`) | The orchestration loop calls `parent_sim = parent.final_simulation()` and then drops `parent`. **No descendant carries the pre-fork trace.** |
| Parent event ledger | `parent.events` | Same — `final_simulation()` returns the IR, not the events that produced it. |
| Parent-level `clone_seed` provenance | `seed + clone_idx * 1_000_000` formula | Lives in Python orchestration code; not stamped on any descendant. A consumer who wants "which seed produced clone 0?" has to recompute it. |
| Parent's "would-have-been" AIRR projection | Never materialized | The pre-fork run never gets projected to an AIRR record; there is no "row 0" for the family ancestor. |

The architectural shape today: **the parent is a transient IR
snapshot**. The family is reconstructable *only* because every
descendant's IR contains the post-recombination state — but that
reconstruction is implicit and partial (e.g. the parent's
NP-generation event identities are gone).

### Gap pin

`pin_absence_parent_trace_not_carried_on_descendants` — every
descendant's `outcome.trace` length equals the number of post-fork
sampling sites, never the parent's (pre-fork + post-fork) total.

---

## 3. Q2 — Which AIRR fields should be family-invariant?

A *family-invariant* field is one whose value should be identical
across every descendant of a clone, by biological necessity. A
*descendant-variant* field is one whose value is allowed (and
expected) to differ within a clone.

The biological model: every descendant of a clone shares the same
V(D)J recombination event — the same V, D, J allele picks, the
same trims, the same NP bytes, the same junction coordinates *in
the germline frame of reference*. SHM, sequencing artefacts, and
library-prep passes are the only sources of within-clone divergence.

### Proposed family-invariant fields

| AIRR field | Source | Why invariant |
|---|---|---|
| `v_call` (truth, not live) | Pre-fork V allele assignment | Same recombination event |
| `d_call` (truth, not live) | Pre-fork D allele assignment | Same recombination event |
| `j_call` (truth, not live) | Pre-fork J allele assignment | Same recombination event |
| `truth_v_call` / `truth_d_call` / `truth_j_call` | `expose_provenance=True` | Already verified by [`test_g8_truth_columns_via_clonal_pipeline`](../tests/test_experiment.py#L1707) |
| `original_v_call` | Pre-fork receptor-revision provenance | Receptor revision happens in the pre-fork pass; the "what V did we start with?" value is identical for every descendant |
| `d_inverted` | Pre-fork D orientation flag | D inversion is a pre-fork decision |
| Junction in germline frame | Pre-fork pool bytes between CDR3 anchors | SHM may modify *observed* junction bases but the germline-frame coordinates and the *pre-SHM* junction string are family-invariant |
| `junction_aa` *evaluated pre-SHM* | Translation of pre-SHM junction | If we projected the parent IR, this would be the canonical clonotype protein sequence |
| `cdr3_start` / `cdr3_end` in germline coords | Pre-fork anchor positions | Same recombination event, same anchors |
| `clone_id` | Orchestration-loop assignment | Tautological — the tag itself defines the family |

### Caveats — what is *not* invariant under heavy SHM today

The live-call (`v_call` / `d_call` / `j_call` columns as projected)
**can** diverge across descendants under heavy SHM even though the
truth allele is invariant: the tie-set widens when SHM erases the
distinguishing positions. [`test_g5_clonal_descendants_share_
junction_diverge_via_mutate`](../tests/test_experiment.py#L1815)
already pins a soft form of this: it asserts only that the
**non-empty** v_call sets agree, allowing the heavy-SHM cases where
some descendants get a wider tie set or a tie-list-of-equals. A
proper family validator must distinguish:

- **Truth-allele invariance** — must hold (provenance).
- **Live-call invariance** — does not hold under heavy SHM; the
  family validator should check `truth_v_call in v_call.split(",")`
  for every descendant, not equality of `v_call` strings.

### Junction-as-string invariance is subtle

`junction` (the AIRR field) is the *observed* junction in the
final sequence — post-SHM, post-PCR, post-N-corruption. It is
**not family-invariant**. What is invariant is the *pre-SHM*
junction, which today exists only inside `parent_sim` and is never
projected to a record. A family validator that wants to check
"shared junction-pre-SHM" needs either:

1. The parent IR materialized as a record (not done today), or
2. A way to reconstruct the pre-SHM junction from each descendant
   by reversing the SHM events in its trace.

Option 2 is brittle (per-descendant reconstruction, requires
complete SHM event coverage). Option 1 is clean but needs a
`FamilyRecord` / parent-projection concept (see Q6).

### Gap pin

`pin_absence_pre_shm_junction_not_projected` — no AIRR field today
exposes the pre-SHM junction; the only way to read it is via the
parent IR snapshot, which is not retained on any returned object.

---

## 4. Q3 — Which fields should vary per clone member?

The dual of §3. A *per-member-variant* field is one whose value
should differ across descendants of a clone, by biological
necessity — checking that they don't *all* match would be a real
test that SHM / library prep / sequencing is actually doing
something.

### Proposed per-clone-variant fields

| AIRR field | Source | Why variant |
|---|---|---|
| `sequence` | Post-SHM, post-corruption final pool | SHM and library-prep passes are independent per descendant |
| `sequence_alignment` | Sequence + alignment string | Tracks `sequence` divergence |
| `v_sequence_alignment` / `d_sequence_alignment` / `j_sequence_alignment` | Per-segment alignment | Tracks `sequence` divergence |
| Mutation count / rate per record | SHM events | Each descendant draws independent SHM events |
| `productive` | Computed from final `sequence` | Heavy SHM can break productive frame on some descendants |
| `complete_vdj` | Computed | End-loss + corruptions can clip some descendants |
| `n_count` / ambiguous-base markers | Post-corruption | Independent per descendant |
| `r1_sequence` / `r2_sequence` (paired-end) | Paired-end pass | r1/r2 fragmentation is per descendant |
| `read_layout` choice variability | Paired-end pass | If `random_strand_orientation` is in the pipeline, per-descendant |
| Indel events / indel-coordinate fields | `polymerase_indels` | Independent per descendant |
| `sequence_id` | `f"clone{c}_desc{d}"` formatter | Tautologically distinct per descendant |

### Caveat — descendant SHM can reach "saturation"

Under `mutate(rate=0.3)` or `mutate(count=100)`, descendants of a
small clone can saturate to near-uncorrelated sequences. A
family-level "mutation-distance distribution" check (Q5) must use
the **truth allele** as the reference point, not the parent-IR
junction string, so that a descendant whose post-SHM live-call
widens to a tie-set still contributes a sensible distance value.

### Gap pin

`pin_absence_per_descendant_mutation_distance_not_aggregated` —
today no API surfaces "for each descendant, distance to parent
sequence" or "for clone C, the distribution of distances". A
consumer has to write the aggregation themselves from the records.

---

## 5. Q4 — Is the trace model sufficient?

The user's most architecturally loaded question.

### What the trace model looks like today

- **Pre-fork run.** `pre.run(clone_seed)` produces a parent
  `Outcome` with its own `PyTrace`. The trace contains the
  recombination addressed choices: V/D/J allele samples, trim
  draws, NP-length and NP-base samples, anchor decisions.
- **`final_simulation()`.** Extracts the IR (assembled pool, region
  map, assignments). The IR contains the *state* the trace
  produced, but not the trace itself. After this call, the parent's
  trace is unreachable from the descendant.
- **Post-fork run.** `post.run_from(parent_sim, desc_seed)`
  produces a descendant `Outcome` whose `PyTrace` only contains the
  post-fork sampling sites (SHM substitutions, PCR errors, indel
  sites, etc.). The descendant's trace **cannot** be used to
  replay the recombination — those choices were already baked into
  the IR it inherited.

### What's sufficient

For per-descendant replay of post-fork passes, the model is
sufficient: given `parent_sim`, `desc_seed`, and `post`, you can
re-derive the descendant exactly (this is what the existing G5
"identical-without-post-fork-passes" test relies on).

### What's not sufficient

1. **Per-family replay is two-pass.** To replay clone C
   descendant D, a consumer needs both `clone_seed = seed + C *
   1_000_000` to reconstruct the parent IR (via `pre.run`) and
   `desc_seed = clone_seed + 1 + D` to reconstruct the descendant
   (via `post.run_from`). The descendant's `trace_file` alone is
   insufficient — it carries the post-fork plan signature but not
   the pre-fork plan signature.

2. **Trace files don't round-trip a family.** Today's
   `trace_file_from` ([`engine_rs/src/python/compiled.rs:197`](../engine_rs/src/python/compiled.rs#L197))
   bundles `(plan, refdata, seed, trace)` for a single outcome's
   simulator. A clonal-family trace file would need *two* plan
   signatures (`pre`, `post`) plus both seeds. **No clonal trace
   file format exists.**

3. **Parent provenance not addressable.** Consumers can't ask "what
   are the unique parent recombination identities in this batch?"
   without re-running every `pre.run(clone_seed)`. There is no
   `ParentOutcome` materialized per clone.

4. **No "parent-replay" entry point.** `CompiledClonalExperiment`
   has no public method to ask "give me just the parent IR / parent
   outcome / parent record for clone C." A consumer who wants the
   parent representation must reach through `._pre` and call
   `.run(clone_seed)` themselves, recomputing the seed formula.

### The "parent trace + child traces" alternative

A future design could materialize one `ParentOutcome` per clone
(holding the pre-fork trace + event ledger) plus K
`DescendantOutcome`s per clone (holding only post-fork traces),
with descendants carrying a back-pointer to their parent. This
would:

- Make family replay a single bundle (`parent + descendant traces`).
- Make parent-only projection possible (e.g. "give me the AIRR
  record of the family ancestor").
- Make per-descendant SHM-event analysis cleaner (the pre-SHM
  baseline is right there on the parent).

The cost: structural change to `Outcome` (or a new type), and
non-trivial wire changes through `SimulationResult` and
`validate_records`. This audit recommends but does not implement.

### Gap pin

- `pin_absence_no_parent_outcome_per_clone` — there is no
  parent-`Outcome` materialized per clone; the parent IR is
  transient.
- `pin_absence_descendant_trace_excludes_parent_events` — a
  descendant's `outcome.trace` does not include pre-fork addressed
  choices.

---

## 6. Q5 — How should family validation work?

`SimulationResult.validate_records` iterates `(outcome, record)`
pairs independently — it has no awareness of clonal structure.
[`pin_absence_validate_records_is_per_record_only`] in the contract
file freezes this. A family-level validator is the audit's primary
proposed extension.

### What a family validator must check

Given the §3/§4 split, a family validator should assert, per
`clone_id` group:

1. **Shared truth allele.** Every descendant's
   `truth_v_call`, `truth_d_call`, `truth_j_call` matches the
   group's first value. (Already true today, tested manually in
   [`test_g8_truth_columns_via_clonal_pipeline`](../tests/test_experiment.py#L1707)
   but not exposed as a built-in check.)

2. **Shared `original_v_call`.** Receptor-revision provenance is
   pre-fork; every descendant must agree.

3. **Shared `d_inverted` flag.** D inversion is pre-fork.

4. **Shared CDR3 germline coordinates.** Anchor positions in
   germline frame are pre-fork.

5. **Truth allele in every descendant's live-call tie set.** For
   every descendant `r`, `r["truth_v_call"] in r["v_call"].split(",")`
   (the live caller never drops the truth allele under SHM, by the
   tie-set design of the score-and-tie caller).

6. **Mutation distance distribution looks like SHM.** Per clone,
   the distribution of `hamming(descendant.sequence,
   parent.sequence)` should be consistent with the configured
   mutate rate / count. Today there is no parent sequence
   available — the validator would have to reconstruct it from any
   descendant by reversing its SHM trace. Cleaner: project the
   parent IR (Q6).

7. **Within-clone sequence non-identity under SHM.** If `mutate`
   or any per-descendant pass is in the post-fork plan, every
   descendant pair within a clone should differ (heavy-SHM
   collisions are vanishingly rare). Today this is checked in
   isolation by [`test_g5_clonal_descendants_share_junction_diverge_
   via_mutate`](../tests/test_experiment.py#L1815) but it is a
   single test, not a release-tier postcondition.

### What a family validator must *not* check

- **`v_call` / `d_call` / `j_call` string equality.** Under heavy
  SHM the live caller widens tie sets per descendant; literal
  equality would surface noise.
- **`junction` (observed) equality.** Post-SHM junctions diverge.
  The pre-SHM junction is the right invariant (and today is not
  projectable).
- **`productive` equality.** Per-descendant by design.

### Default behaviour proposal

If a family validator slice lands, it should slot into the
existing `validate_records=True` kwarg:

- `exp.run_records(seed=0, validate_records=True)` on a
  `CompiledClonalExperiment` runs both per-record validation and
  the family-level group checks.
- A new `RecordValidationFailedError` issue kind
  (`family_invariant_violated`) carries the `clone_id` and the
  divergent field name.
- Default stays `False` for performance and back-compat (matches
  the projection-validator slice that already shipped).

### Gap pin

- `pin_absence_no_family_level_validator` — `validate_records`
  ignores `clone_id`. The contract test demonstrates the gap by
  asserting that a synthetic divergence of `truth_v_call` across
  descendants of a clone (impossible in practice, but constructible
  by mutating the returned records) is **not** caught by
  `validate_records(refdata)`.

---

## 7. Q6 — Architecture debt: do we need `ClonalOutcome` / `FamilyRecord`?

The user's headline question: "if clonal simulation currently
reuses APIs meant for independent records, this audit will tell us
whether we need a real `ClonalOutcome` / `FamilyRecord` concept
before adding more biology on top."

### The answer

**Yes — before any further biology lands on clonal expansion, we
need a first-class family concept.**

### Why

The current shape — flat `List[Outcome]` + flat
`SimulationResult` records + a `clone_id` int tag — is **load-
bearing exactly to the point of "n_clones × per_clone independent
records with shared truth alleles"**. Everything beyond that point
requires either ad-hoc per-call reconstruction (write your own
group-by-`clone_id`) or external recomputation (re-run `pre.run`
with the right seed). Concrete consequences of the missing family
concept:

1. **The parent recombination is unobservable.** No `Outcome`, no
   AIRR record, no trace, no addressable identity. Anything that
   wants "the family ancestor" — lineage trees, germinal-center
   selection coefficients, mutation-distance baselines —  has to
   recompute it.

2. **The "family" identity is a free-floating int.** `clone_id` is
   just a position in the orchestration loop. It has no stable
   identity tied to refdata or to a hash of the recombination
   choices. Two batches with the same seed will produce the same
   `clone_id → V/D/J` mapping by coincidence of the seed scheme,
   not by design.

3. **Validation can't address the family.** `validate_records`
   iterates pairs; there is no second-pass group hook.

4. **Trace files don't bundle families.** `trace_file_from` is
   per-`Outcome`; a clonal-batch trace round-trip is not possible.

5. **The orchestration loop is in Python.** All seed accounting
   (`clone_seed = seed + clone_idx * 1_000_000`), all `clone_id`
   tagging, and all parent-IR sharing is Python code in
   `CompiledClonalExperiment.run_records`. The Rust side is unaware
   of family structure. Anything new (selection, time evolution,
   parent projection, etc.) that wants to be efficient or part of
   the pass framework would need a Rust counterpart.

6. **The post-fork plan can't see family-level facts.** Today the
   post-fork compile (`analyze.rs:has_recombination_support()`)
   has to skip the productive-frame precondition because it sees
   no recombination support — it has no addressable concept of
   "the parent IR carries those facts." If a future pass needed to
   look at parent assignments at compile-time (e.g. "set SHM rate
   higher in CDR than FR" with germline-position resolution
   referenced from the parent), it would have to invent the
   addressability from scratch.

7. **Composability with `paired_end` / `random_strand_orientation`
   etc. is already at the edge.** These passes treat each
   descendant as independent. That happens to be correct
   biologically (paired-end fragmentation is per-read), but the
   pattern reinforces "descendants are independent records that
   happen to share an IR" rather than "descendants are reads
   sampled from a family."

### What "first-class family" should mean

A minimal, non-presumptuous proposal — not specified for
implementation here, only sketched so the audit's contract pins
have a target:

- **`FamilyOutcome` / `ClonalFamily`** — a Rust-side type holding
  the parent `Simulation` (post-recombination IR), the parent
  trace, the parent event ledger, and a list of descendant
  `Outcome`s.
- **`FamilyRecord`** — a Python-side projection: a dict of family-
  invariant AIRR fields (`truth_v_call`, `truth_d_call`,
  `truth_j_call`, `d_inverted`, `original_v_call`, pre-SHM
  `junction`, pre-SHM `junction_aa`, CDR3 germline coords) plus a
  list of descendant records.
- **`SimulationResult.families`** — `Optional[List[FamilyRecord]]`;
  `None` for non-clonal results, populated by
  `CompiledClonalExperiment.run_records`.
- **Family-level `validate_families`** — runs §6 group checks.
- **`Outcome.parent_id`** — `Optional[int]` back-pointer from a
  descendant to its `FamilyOutcome` in the same batch.

The current `clone_id` int + flat record list **should remain**
as the "thin tag" view. Downstream consumers who want the family
structure use `.families`; consumers who want a flat per-read
table keep using the existing API. No back-compat break.

### Cost

Non-trivial: at least one new Rust type, changes to
`SimulationResult`, a new validator pass, and an opt-in
projection of the parent IR (which needs an AIRR-record shape
that supports "no observed sequence" — the parent is a germline-
frame entity). This audit deliberately does not commit to a slice
shape.

### Gap pin

- `pin_absence_no_familyrecord_type` — no `FamilyRecord` /
  `ClonalFamily` / `FamilyOutcome` type exists in the public
  package.
- `pin_absence_simulationresult_has_no_families_attr` —
  `SimulationResult` exposes `.records` and `.outcomes`; no
  `.families`.
- `pin_absence_outcome_has_no_parent_id` — `Outcome` has no
  `parent_id` / `family_id` field.

---

## 8. Edge cases the implementation slice must handle

A non-exhaustive list of cases any future family-layer
implementation must address. These inform the contract pins but
no slice is proposed.

1. **Non-clonal experiments.** `Experiment` without
   `expand_clones()` must still work; `SimulationResult.families`
   stays `None`; `validate_records` keeps its current per-record
   semantics.

2. **Empty post-fork plan.** `expand_clones()` immediately followed
   by `run_records()` with nothing after (the
   "identical-without-post-fork-passes" case): every descendant is
   identical to the parent. The family validator must accept this
   as the trivial case, not flag "no within-clone divergence" as a
   bug.

3. **Productive-only with heavy SHM.** Under `productive_only()` +
   high `mutate(rate)`, the post-fork pass can repeatedly reject
   descendants. The family-level total record count guarantee
   today is structural (`n_clones * per_clone`); the family
   validator must not assume `len(family.descendants) ==
   per_clone` until that interaction is audited.

4. **Receptor revision in the post-fork plan.** Today receptor
   revision is inlined into `recombine` (see
   [`docs/receptor_revision_design.md`](receptor_revision_design.md))
   and is therefore always pre-fork. A future contributor who
   accidentally lowers a receptor-revision step into the post-fork
   plan would produce per-descendant `original_v_call` divergence,
   which the family validator should catch.

5. **D inversion in the post-fork plan.** Same as above —
   `invert_d` is structural and pre-fork. A misordered pipeline
   that pushed `invert_d` into the post-fork plan would produce
   per-descendant `d_inverted` divergence; family validator
   catches it.

6. **`expose_provenance=True` interaction.** The truth columns
   are already family-invariant (test_g8 pins this). A family-
   level projection should reuse `_inject_truth_columns` rather
   than re-derive.

7. **Seed collision across `clone_seed = seed + clone_idx *
   1_000_000`.** Today this formula assumes `n_clones * 1_000_000
   < u64::MAX`, which is fine in practice but is an undocumented
   precondition. A future implementation slice should either pin
   it explicitly or replace the formula with `(seed, clone_idx)`
   addressability.

8. **`paired_end` + clonal.** The release-tier paired-end test in
   [`test_paired_end.py`](../tests/test_paired_end.py) doesn't
   currently cross with `expand_clones`. The family validator
   must not assume one record per descendant in a paired-end
   world (each descendant still maps to a single record today,
   but the audit flags this as a future check).

---

## 9. Replay determinism

The audit's stance:

- **Per-descendant replay** is already deterministic. Given
  `(pre, post, refdata, seed)`, descendant
  `(clone_idx, desc_idx)` reproduces byte-for-byte across runs.
  Existing test infrastructure (sha256-of-records) covers this.

- **Per-family replay** is two-pass deterministic but not
  encapsulated in a single object. A consumer who wants "give me
  family C from this batch" today has to:
  1. Reconstruct `clone_seed = seed + C * 1_000_000`
  2. Call `compiled._pre.run(clone_seed)` to get the parent
  3. Loop K times with `desc_seed = clone_seed + 1 + i` and
     `compiled._post.run_from(parent_sim, desc_seed)`.

  This works but reaches through private `_pre` / `_post`
  attributes. A future slice should expose a public method like
  `CompiledClonalExperiment.run_family(seed, clone_idx)` returning
  `(parent_outcome, descendant_outcomes)`.

- **Family replay across batches** is not addressable. A
  `trace_file` for a single descendant exists but does not encode
  the parent's plan signature; a "family trace file" doesn't
  exist. A `FamilyRecord`-bearing format would be the natural
  place to add it.

### Gap pin

`pin_absence_no_run_family_public_method` —
`CompiledClonalExperiment` has no `run_family(seed, clone_idx)`
entry point.

---

## 10. Validator integration

The audit's recommended integration shape — not implemented here,
but specified so the family-validator slice has a target.

### Two layers, mirroring projection validator

The shipped `validate_records=True` kwarg already mirrors this
shape for per-record postcondition validation. The family-layer
addition:

- **Per-record layer (today, unchanged).** Every descendant runs
  through `outcome.validate_record(refdata, sequence_id=...)`.
  Family-level concerns are out of scope at this layer.

- **Family layer (proposed).** A second pass over the records
  grouped by `clone_id` runs §6's invariants. Issues are reported
  with a new `RecordValidationFailedError` family of issue kinds:
  `family_truth_v_divergent`, `family_truth_d_divergent`,
  `family_truth_j_divergent`, `family_d_inverted_divergent`,
  `family_original_v_call_divergent`,
  `family_descendants_identical_under_shm` (degenerate-SHM
  warning), `family_truth_not_in_live_call_tie_set` (per
  descendant).

### Trigger surface

Reuse `validate_records=True` on `CompiledClonalExperiment.run_records`
(already accepts the kwarg, today only does per-record validation).
No new kwarg; the family-layer pass activates automatically when
the result has clonal structure.

### Output shape

`ValidationReport.failures` already carries dicts; add an
optional `clone_id` key to the failure dict. A failure on a
family invariant uses `record_index = None` and `clone_id = C`;
a per-record failure uses `record_index = i` and omits or
duplicates `clone_id = records[i]["clone_id"]`.

### Cache parity interaction

`Outcome.check_live_call_cache_parity` is per-outcome. It does
not change at the family layer — every descendant continues to
run the cache-parity check independently if a release workflow
opts in. The audit recommends *not* layering a "family cache
parity" check (the per-outcome check already covers the live-
call surface).

---

## 11. Backwards compatibility

A family-layer slice should aim for zero break to existing
consumers. The proposed shape preserves:

- `SimulationResult.records` — same flat list, same fields, plus
  the existing `clone_id` int.
- `SimulationResult.outcomes` — same flat list.
- `clone_id` integer — still in `[0, n_clones)`, still tagged at
  projection time.
- `expand_clones(n_clones, per_clone)` — same constructor.
- `validate_records=True` default `False` — unchanged.
- All G5 tests continue to pass.

The additions (`SimulationResult.families`, `Outcome.parent_id`,
`FamilyRecord`, `FamilyOutcome`, family-validator failures) are
new surfaces; existing call sites that don't reference them are
unaffected.

The one place a careful contributor must watch: the seed scheme
(`clone_seed = seed + clone_idx * 1_000_000`) is currently
internal but visible by behavioural side-effect. A future slice
that changes the scheme breaks byte-reproducibility across
versions. The audit recommends documenting the scheme as
**internal but stable** (don't change it without a deprecation
cycle) and exposing `run_family(seed, clone_idx)` so consumers
don't reinvent it.

### Gap pin

`pin_scaffold_clone_seed_formula_is_stable` — pin the current
seed scheme. If a future slice changes it, this test surfaces it
in CI before the change reaches a release.

---

## 12. Implementation order

This audit deliberately does not propose an implementation slice.
A future slice (or sequence of slices) might shape up like:

1. **Slice 1 — family-level validator (Python-only).** Add
   `_validate_family_invariants(records, refdata)` that groups by
   `clone_id`, runs §6 invariants 1-5 (no parent IR needed —
   uses `truth_*` columns + live-call tie-set check). Hook into
   `validate_records=True` when records carry `clone_id`. New
   failure kinds in `RecordValidationFailedError`. No Rust changes.
   No `FamilyRecord` yet. Most user benefit per token cost; lowest
   architectural risk. **Recommended as the first slice.**

2. **Slice 2 — parent outcome materialization.** Retain the
   parent `Outcome` per clone inside `CompiledClonalExperiment`
   and expose `run_family(seed, clone_idx)`. Optionally surface
   `SimulationResult.parents` as a flat list of parent outcomes
   (one per clone). No `FamilyRecord` yet — consumers compose
   manually. Modest Rust changes (no new types yet).

3. **Slice 3 — `FamilyRecord` projection.** Add the
   `FamilyRecord` type, populate `SimulationResult.families`,
   add `Outcome.parent_id` back-pointer. This unlocks
   germline-frame projections of the parent (pre-SHM junction,
   pre-SHM junction_aa, etc.) and the §6.6 mutation-distance
   distribution check.

4. **Slice 4 — clonal trace files.** Extend `trace_file_from`
   for `CompiledClonalExperiment`, producing a family-bundle
   format with both plan signatures + both seeds. Defer until
   demand surfaces.

5. **Slice 5+ — biology on top.** Lineage trees, germinal-center
   selection, antigen-specific weighting, time-evolved clones.
   All of these compose against Slice 3's `FamilyRecord` /
   `FamilyOutcome` surface; attempting them on the flat
   `clone_id` tag surface would repeat the architecture debt
   this audit pinpoints.

---

## 13. Test surface — what the audit pins

Mirrored in [`tests/test_clonal_family_contract.py`](../tests/test_clonal_family_contract.py).

### `pin_scaffold_*` — existing surfaces the audit relies on

1. `_ClonalForkStep` exists as a frozen dataclass with `n_clones`
   and `size` fields.
2. `Experiment.compile()` returns `CompiledClonalExperiment` when
   a fork step is present; `CompiledExperiment` otherwise.
3. `CompiledClonalExperiment` exposes `n_clones`, `size`,
   `total_records`, `refdata`, `run`, `run_records`.
4. The parent-IR sharing path uses `pre.run(...)` →
   `.final_simulation()` → `post.run_from(...)`.
5. The orchestration loop sets `clone_id` in `[0, n_clones)` on
   every record.
6. `clone_seed = seed + clone_idx * 1_000_000` is the current
   scheme — if this changes, the byte-reproducibility contract
   breaks.
7. `CompiledSimulator.run_from` exists on the Rust side and is
   docstring-named "clonal expansion."

### `pin_absence_*` — gaps the audit identifies

1. `pin_absence_parent_trace_not_carried_on_descendants` — every
   descendant `outcome.trace` carries only post-fork events.
2. `pin_absence_no_parent_outcome_per_clone` — no
   parent-`Outcome` is exposed; the parent IR is transient.
3. `pin_absence_pre_shm_junction_not_projected` — no AIRR field
   exposes the pre-SHM junction.
4. `pin_absence_per_descendant_mutation_distance_not_aggregated`
    — no API surfaces a per-clone mutation-distance distribution.
5. `pin_absence_validate_records_is_per_record_only` —
   `validate_records` ignores `clone_id`; a synthetic
   `truth_v_call` divergence across clone members is not caught.
6. `pin_absence_no_family_level_validator` — `_validate_family_
   invariants` does not exist.
7. `pin_absence_no_familyrecord_type` — `FamilyRecord` /
   `ClonalFamily` / `FamilyOutcome` types not in the public
   package.
8. `pin_absence_simulationresult_has_no_families_attr` — no
   `.families` attribute on `SimulationResult`.
9. `pin_absence_outcome_has_no_parent_id` — no `parent_id` /
   `family_id` field on `Outcome`.
10. `pin_absence_no_run_family_public_method` —
    `CompiledClonalExperiment` has no `run_family(seed, clone_idx)`.

### Flip behaviour when implementation slices land

When Slice 1 lands, `pin_absence_no_family_level_validator` and
`pin_absence_validate_records_is_per_record_only` flip to
`pin_family_validator_catches_truth_v_divergence` (and similar).
When Slice 3 lands, `pin_absence_no_familyrecord_type` and
`pin_absence_simulationresult_has_no_families_attr` flip to
`pin_family_record_projects_pre_shm_junction` (and similar).

---

## 14. Out of scope

Documented here so a future contributor doesn't accidentally
expand any implementation slice.

- **Germinal-center / selection biology.** Out of scope until at
  least Slice 3 lands. Stacking selection on the flat-record
  surface would repeat the architecture debt.
- **Lineage tree shape (branching, time evolution).** Out of
  scope; would require a tree-of-`FamilyOutcome`s extension well
  past Slice 3.
- **Antigen-specific allele weighting.** Composes against today's
  allele-weight surface (`restrict_alleles` /
  `_AlleleWeightsStep`); does not require family infrastructure.
  Out of scope for this audit because it doesn't touch the
  family layer.
- **Cross-clone shared SHM hotspots / kernel sharing.** All SHM
  passes today are per-descendant by construction; cross-clone
  sharing would be a kernel-design change orthogonal to the
  family concept.
- **TCR clonal expansion.** The DSL accepts `expand_clones` for
  TCR loci already (since `mutate` is the only step gated to
  BCR); no biology change is implied. Out of scope.
- **Mixed BCR/TCR families.** Biologically nonsensical; not on
  the table.
- **`SimulationResult` loaded from TSV regaining family
  structure.** TSV round-trip loses `outcomes`; family
  reconstruction from records-only is out of scope for the same
  reason `validate_records` already rejects records-only results.

---

## Summary table

| Question | Answer |
|---|---|
| What does the current clonal path share? | Parent IR (V/D/J + trim + NP + assembled pool) — transmitted as a `Simulation` snapshot, then discarded. Parent trace + ledger are not retained. |
| Which AIRR fields should be family-invariant? | Truth alleles (`truth_v/d/j_call`), `original_v_call`, `d_inverted`, CDR3 germline coords, pre-SHM junction. Live `v_call` is only invariant in the "truth in tie set" sense under heavy SHM. |
| Which fields should vary per clone member? | Final `sequence` (post-SHM/PCR/N/indel), per-segment alignment columns, `productive`, `complete_vdj`, paired-end r1/r2, indel coordinates. |
| Is the trace model sufficient? | For per-descendant replay yes; for family replay it is two-pass and not encapsulated. No `parent_id`, no clonal-trace-file format, no `run_family` entry point. |
| How should family validation work? | Group records by `clone_id`, assert truth-allele / d_inverted / original_v_call invariance, assert truth-in-live-tie-set, assert within-clone divergence under SHM. Surface via the existing `validate_records=True` kwarg. |
| Do we need `FamilyRecord` / `ClonalOutcome`? | **Yes**, before any further biology lands. Today's flat-list-of-independent-records-plus-`clone_id` model is load-bearing exactly to the current capability frontier; lineage / selection / time evolution all require addressable family identity. |
| First implementation slice? | A Python-only family-level validator hook into `validate_records=True`. Lowest risk, highest user benefit per token, no Rust changes. |

The audit's verdict: clonal expansion ships a correct
*expansion generator*. It does not yet ship a *family concept*.
Closing the gap is a structural slice, not an incremental one.
