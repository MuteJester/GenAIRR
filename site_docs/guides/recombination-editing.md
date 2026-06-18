# D inversion and receptor revision

<p class="lead">Two advanced recombination-stage mechanisms that
edit the parent's V(D)J truth: <strong>D inversion</strong> flips
the D allele into reverse-complement orientation; <strong>receptor
revision</strong> replaces V after the initial recombination with
a different germline allele. Both are heavy-chain biology, both
must be configured before clonal fork, and both leave structured
provenance on every AIRR record so the downstream analysis can
tell what happened.</p>

## What this guide covers

D inversion and receptor revision belong together because they
share four properties:

- Both are **recombination-stage / ancestor-phase** decisions -
  every descendant of a clone inherits the parent's edited state.
- Both are **heavy-chain (VDJ) only**; the DSL rejects them on
  light chains and TCR-alpha/gamma at chain time.
- Both leave **dedicated AIRR provenance fields** (`d_inverted`
  for inversion; `receptor_revision_applied` + `original_v_call`
  for revision) so downstream pipelines can branch on what
  happened biologically.
- Both are **validator-checked**: `validate_records` independently
  re-derives the provenance fields from the simulation IR and
  surfaces any divergence.

Neither is on by default - each requires an explicit
probability-keyword opt-in.

## D inversion

In real B-cell biology, the D segment is occasionally recombined
in reverse-complement orientation rather than its canonical
forward orientation. GenAIRR models that with a per-record
Bernoulli draw: with probability `prob`, the D allele's bytes are
written into the pool reverse-complemented; with probability
`1 - prob`, forward orientation is committed.

```python
import GenAIRR as ga

result = (
    ga.Experiment.on("human_igh")
      .recombine()
      .invert_d(prob=0.05)
      .run_records(n=100, seed=1)
)

# How many records drew the inverted D?
print(sum(1 for r in result.records if r["d_inverted"]))
```

### API

```python
.invert_d(*, prob: float = 0.05) -> Experiment
```

- **Keyword-only**, single argument.
- **VDJ-only.** Raises `ValueError: invert_d is only valid for
  VDJ chains (current chain_type='vj')` at chain time on VJ
  cartridges.
- **Single call per pipeline.** A second `.invert_d(...)` raises
  `ValueError: invert_d already configured on this experiment;
  v1 accepts at most one inversion step per pipeline`.
- **Must come before clonal forks** (ancestor-phase). See
  [Clonal placement](#clonal-placement) below.

### The `d_inverted` AIRR field

```python
rec["d_inverted"]    # bool - True when this record's D was reverse-complemented
```

| Pipeline + outcome | `d_inverted` value |
|---|---|
| VJ cartridge (any pipeline) | `False` on every record |
| VDJ cartridge, no `.invert_d(...)` in pipeline | `False` on every record |
| VDJ cartridge, `.invert_d(prob=0.05)`, draw fired False | `False` |
| VDJ cartridge, `.invert_d(prob=0.05)`, draw fired True | `True` |

The field is sourced from the final IR, not from the trace - the
post-pipeline D orientation is treated as canonical ground
truth.

### Coordinates and CIGAR under inversion

A crucial detail: when `d_inverted == True`, the **`d_germline_start`
/ `d_germline_end` coordinates remain in the original (forward)
allele orientation**. They don't track the inverted orientation
that's in the assembled pool. `d_cigar` continues to use only
canonical M/I/D ops (never X), tallied over the inverted pool
bytes vs the forward allele bytes - so under inversion the
match-rate `d_identity` drops, but the M-op length still equals
the retained slice's allele-coordinate span.

This convention means:

- `d_germline_start/end` are stable forward-allele coordinates
  regardless of orientation.
- `d_identity` reflects the alignment evidence (low identity is
  expected under inversion).
- `d_inverted` is the single source of truth for "did inversion
  happen" - don't try to infer it from `d_identity` or coordinate
  patterns.

### D allele calls are orientation-aware

The live-call walker that produces `d_call` is orientation-aware:
it complements the observed pool byte on the fly before scoring
against the reference index. So inverted records still get a
meaningful `d_call` - the same forward allele id, just scored
through a single-complement pass instead of via a parallel
inverted index. The orientation lives on the AlleleInstance,
not on the call.

## Receptor revision

In real B-cell biology, the V segment can be replaced after
initial recombination - receptor editing in the bone marrow.
GenAIRR models that with a per-record Bernoulli draw at the
post-recombine point: with probability `prob`, the originally
sampled V is replaced with a different germline V allele; with
probability `1 - prob`, the original V is kept.

```python
result = (
    ga.Experiment.on("human_igh")
      .recombine()
      .receptor_revision(prob=0.05)
      .run_records(n=100, seed=1)
)

print(sum(1 for r in result.records if r["receptor_revision_applied"]))
```

### API

```python
.receptor_revision(*, prob: float = 0.05) -> Experiment
```

- **Keyword-only**, single argument.
- **VDJ-only.** Raises `ValueError: receptor_revision is only
  valid for VDJ chains` at chain time on VJ cartridges.
- **Single call per pipeline.** A second call raises
  `ValueError: receptor_revision already configured on this
  experiment`.
- **Must come before clonal forks** (ancestor-phase).
- The compile path also checks that the cartridge has a V pool
  in refdata - required because the replacement allele draws from
  it.

### Provenance fields

```python
rec["receptor_revision_applied"]    # bool - did revision fire on this record?
rec["original_v_call"]              # str - the pre-revision V identity (when applied)
rec["v_call"]                       # str - the POST-revision V evidence-call
```

| Field | When revision NOT applied | When revision applied |
|---|---|---|
| `receptor_revision_applied` | `False` | `True` |
| `original_v_call` | `""` (empty string) | The pre-revision V allele name |
| `v_call` | The recombination-time V evidence-call | The post-revision V evidence-call |

Two facts to internalise about the `v_call` / `original_v_call`
distinction:

- **`v_call` always reports the post-revision identity.** When
  revision fires, the V slot in the IR is rewritten to point at
  the replacement allele; `v_call` reads from that slot and so
  carries the new V. The pre-revision identity survives only on
  `original_v_call`.
- **`original_v_call` is empty when revision DIDN'T fire.** Not
  `None`, not equal to `v_call` - empty string. The empty
  sentinel is deliberate: it reads cleaner than `original_v_call
  == v_call` because the latter would force you to cross-check
  `receptor_revision_applied` to disambiguate "not revised" from
  "revised but happened to land on the same allele".

## Using them together

Both mechanisms compose freely in the recombination-stage block:

```python
result = (
    ga.Experiment.on("human_igh")
      .recombine()
      .invert_d(prob=0.05)
      .receptor_revision(prob=0.02)
      .run_records(n=100, seed=1)
)

# Inverted, revised, or both per record
import pandas as pd
df = pd.DataFrame(result.records)
df.groupby(["d_inverted", "receptor_revision_applied"]).size()
```

The two are independent - a single record can carry both
`d_inverted=True` and `receptor_revision_applied=True`. The
schedule analyser orders them per their pass requirements
(`invert_d` runs before `assemble.D`; `receptor_revision` runs
after V assembly).

## Clonal placement

Both methods must be configured BEFORE a clonal fork
(`clonal_lineage`, `clonal_repertoire`, or legacy `expand_clones`) -
they're recombination-time decisions that the family inherits.
The DSL enforces this at chain time with two symmetric guards:

```python
# WRONG - invert_d after a clonal fork raises ValueError
ga.Experiment.on("human_igh") \
   .recombine() \
   .clonal_repertoire(n_clones=10, max_size=20) \
   .invert_d(prob=0.05)
# ValueError: invert_d must be called before the clonal fork;
# D inversion is a recombination-time decision and must be
# inherited by all clone descendants. Move the invert_d(...)
# call before clonal_lineage(...), clonal_repertoire(...), or expand_clones(...).
```

Symmetric message for `receptor_revision` placed after
a clonal fork. The right order is:

```python
result = (
    ga.Experiment.on("human_igh")
      .recombine()
      .invert_d(prob=0.05)              # ancestor phase
      .receptor_revision(prob=0.02)     # ancestor phase
      .clonal_repertoire(n_clones=10, max_size=20)
      .mutate(model="s5f", rate=0.03)   # descendant phase
      .run_records(seed=1)
)
```

Each clone inherits the parent's `d_inverted` value and
`receptor_revision_applied` / `original_v_call`. Every descendant
of the same clone carries the same triple of provenance fields -
that's precisely the family-truth invariance `validate_families`
checks (via the `truth_*_call` fields when provenance exposure is
on; the parent-aware validator also compares `d_inverted` and
`original_v_call` against the parent Outcome).

See [Clonal simulation overview](clonal-families.md) for the full
ancestor-vs-descendant phase discipline.

## Validation and replay

The postcondition validator re-derives both mechanisms' fields
from the IR, independently of the AIRR projection path. Three
issue kinds can fire:

- **`DInvertedMismatch { reported, expected }`** - `d_inverted`
  on the record disagrees with the IR's D orientation. Re-derived
  by inspecting `Simulation.assignments[D].orientation`.
- **`ReceptorRevisionAppliedMismatch { reported, expected }`** -
  the applied flag disagrees with whether the IR's V slot carries
  a `receptor_revision_original_id`. Re-derived from the slot,
  not from the trace.
- **`OriginalVCallMismatch { reported, expected }`** - the
  reported `original_v_call` doesn't match what `refdata.get(V,
  receptor_revision_original_id).name` resolves to.

Two notable things about the validator coverage:

- **No chain-type gate.** The validator runs the checks on every
  record. VJ records simply have `expected_d_inverted=False`,
  `expected_receptor_revision_applied=False`,
  `expected_original_v_call=""`, so they pass by default.
- **The walker is D-orientation-aware**, so unlike the
  random-strand-orientation rev-comp case (which makes the C4
  allele-call oracle skip), D inversion does NOT skip any
  validator checks. Inverted records validate fully.

### Replay

Both mechanisms record their decisions on stable trace addresses
so replay reproduces them deterministically:

| Mechanism | Trace address(es) |
|---|---|
| D inversion | `sample_allele.d.inverted` (Bool, always recorded) |
| Receptor revision | `receptor_revision.applied` (Bool, always recorded) |
| Receptor revision (when applied) | `receptor_revision.v_allele` (AlleleId) + `receptor_revision.v_trim_3` (Int) |

`replay_from_trace_file(...)` consumes these values verbatim;
`rerun_from_trace_file(...)` re-runs the samplers from the same
seed. Plan-signature and refdata-content-hash gates fire first
in both cases, so changing `prob` between record-time and
replay-time surfaces immediately as a plan-signature mismatch.
See the [Validation hub](../validation/index.md#trace-and-replay)
for the full replay surface.

## Common mistakes

A handful of issues that show up repeatedly with these two
mechanisms.

**Calling `.invert_d()` on a VJ cartridge.** Raises `ValueError:
invert_d is only valid for VDJ chains (current
chain_type='vj')` at chain time. The same applies to
`.receptor_revision()`. Neither mechanism exists for light chains
or TCR-alpha/gamma - recombination biology only models them in
heavy chains.

**Expecting `original_v_call` to be non-empty when revision
DIDN'T fire.** It's the empty string `""`, not `None`, not
`v_call`. The empty sentinel makes "revision didn't fire" a
clean check: `if rec["original_v_call"]: ...` reads as "if we
have a pre-revision identity". Don't try to compute
`original_v_call == v_call`; that's a different question
(answered by `receptor_revision_applied` combined with whether
the replacement allele happened to land on the same V).

**Treating `d_cigar` coordinates as reversed.** Under
`d_inverted=True`, the CIGAR ops still use the forward allele's
coordinate space - `d_germline_start/end` are forward-allele
coordinates. The match-rate `d_identity` drops because the
inverted pool bytes don't match the forward allele bytes, but
the coordinate system stays stable. If your downstream pipeline
re-RCs the bytes before comparing, do it on
`rec["sequence"][d_sequence_start:d_sequence_end]` - don't
remap the coordinates.

**Putting `.invert_d()` or `.receptor_revision()` after a clonal
fork.** The DSL rejects this immediately at chain
time with the symmetric message above. Both decisions must be
shared across every descendant of a clonal family; the API will
not let you accidentally split them.

## Where to go next

- **[Recombination and junction biology](recombination-junction.md)**
  what `.recombine()` produces in the first place, before
  either of these mechanisms run.
- **[Clonal simulation overview](clonal-families.md)** - the ancestor-vs-
  descendant phase discipline both mechanisms participate in.
- **[Validation & reproducibility](../validation/index.md)** -
  the validator's issue catalogue (including the three issue
  kinds named above) and the replay model.
- **[SHM and mutation targeting](shm-targeting.md)** - the
  biology stage that runs *after* recombination edits; SHM
  accumulates on the post-revision V and the inverted D.
- For the engine-side mechanics + audit details, see
  [`docs/d_inversion_design.md`](https://github.com/MuteJester/GenAIRR/blob/master/docs/d_inversion_design.md)
  and
  [`docs/receptor_revision_design.md`](https://github.com/MuteJester/GenAIRR/blob/master/docs/receptor_revision_design.md).
