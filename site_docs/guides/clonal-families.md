# Generate clonal families

<p class="lead">A clonal family is one parent recombination plus
many descendants that share the parent's V(D)J truth but diverge
through somatic hypermutation, library prep, and sequencing
artefacts. <code>expand_clones</code> is the DSL marker that turns
a flat pipeline into a clonal one — and once it's in the chain,
the API rejects misordered steps so you can't accidentally collapse
SHM diversity or split a recombination decision across descendants.</p>

## What clonal simulation means

In real biology, a clonal family is the lineage descended from a
single B cell whose V(D)J recombination has happened. Every cell
in that lineage shares the same V/D/J allele choices, the same
trims, the same NP bases, the same D orientation — those are
*recombination-time* decisions, frozen at fork. What diverges
between cousins is somatic hypermutation (each B cell mutates
independently), library prep artefacts (PCR errors don't cross
between reads), and sequencing geometry (R1/R2 windows are per-read).

GenAIRR models that with a single DSL marker, `expand_clones`. The
methods you append BEFORE it run *once per clone* (and that
outcome is shared across every descendant). The methods you
append AFTER run *once per descendant* (independent draws per
read).

## A minimal clonal example

```python
import GenAIRR as ga

result = (
    ga.Experiment.on("human_igh")
      .recombine()                                          # ancestor phase
      .expand_clones(n_clones=5, per_clone=10)              # fork
      .mutate(model="s5f", rate=0.02)                       # descendant phase
      .run_records(seed=1)
)

print(len(result))            # 50  = 5 clones × 10 descendants
print(result[0]["clone_id"], result[0]["parent_id"])    # 0 0
print(result[10]["clone_id"])                            # 1
```

`expand_clones(n_clones=5, per_clone=10)` produces exactly 5 × 10 =
50 records. Note that `run_records(seed=1)` doesn't take an `n`
argument — the clonal pipeline computes its own total. You can
pass `n=50` as a consistency check, but any other value (including
`n=100`) raises `ValueError` at run time. The pipeline knows how
many records it produces; you don't override it.

## Ancestor vs descendant phase

The DSL partitions steps into two phases:

**Ancestor-phase** — runs once per clonal parent. Use for any
mechanism whose decision must be inherited by every descendant
of a family.

| Pass | Why ancestor |
|---|---|
| `.recombine()` | V/D/J allele choices + trims + NP define the family's identity |
| `.invert_d(prob=...)` | D orientation is a recombination-time event; the family shares it |
| `.receptor_revision(prob=...)` | V replacement happens once during B-cell development; the family shares the post-revision V |

**Descendant-phase** — runs independently per descendant. Use for
mechanisms that should vary within a family.

| Pass | Why descendant |
|---|---|
| `.mutate(...)` | Each memory B cell mutates independently |
| `.pcr_amplify(...)` | PCR errors are read-specific |
| `.polymerase_indels(...)` | Same — per-read library artefacts |
| `.ambiguous_base_calls(...)` | Per-read N-injection |
| `.sequencing_errors(...)` | Per-read quality artefacts |
| `.end_loss_5prime(...)`, `.end_loss_3prime(...)` | Per-read 5'/3' adapter loss |
| `.random_strand_orientation(...)` | Per-read strand decision |
| `.paired_end(...)` | R1/R2 windows are per-read |

The DSL enforces ordering at chain time — not at compile time, not
silently at run time. If you call a descendant-phase method
*before* `expand_clones`, the next call to `expand_clones` raises:

```text
ValueError: mutate must be called after expand_clones();
SHM is descendant-specific in GenAIRR's current clonal model.
Move mutate(...) after expand_clones(...).
```

If you call an ancestor-phase method (`invert_d` / `receptor_revision`)
*after* `expand_clones`, the offending method itself raises with
the symmetric message:

```text
ValueError: invert_d must be called before expand_clones();
D inversion is a recombination-time decision and must be inherited
by all clone descendants. Move the invert_d(...) call before
expand_clones(...).
```

A second call to `expand_clones` raises `expand_clones() can only
be called once per pipeline`.

In practice the rule of thumb is: anything biology fixed once per
cell goes before the fork, anything happening per read goes after.
The DSL catches violations the moment you write them.

## Reading clone IDs and parents

Every descendant record carries two integer fields that wire it
back to its clonal family:

```python
rec = result[0]

rec["clone_id"]    # 0  — family identity; every descendant of clone 0 carries 0
rec["parent_id"]   # 0  — addressing index into result.parents
```

**Today `clone_id == parent_id` by construction** — the two are
stamped separately because they carry distinct semantics
(`clone_id` is the family identity, `parent_id` is the lookup
index into the parents collection), so a future change to the
addressing scheme can move one without retrofitting the other.
Treat them as a pair; address downstream joins by `clone_id`.

The parent outcomes live separately on the result:

```python
result.parents                       # list of length n_clones
len(result.parents)                  # 5

parent_outcome = result.parents[rec["parent_id"]]
# Outcome carrying the pre-fork plan's full trace + event ledger
parent_outcome.final_simulation()    # post-recombine IR
parent_outcome.trace()               # pre-fork addressed-choice trace
parent_outcome.events()              # pre-fork event ledger
```

The flat `result.outcomes` list continues to hold exactly one
entry per descendant record (length `n_clones × per_clone`);
parents are exposed on the separate `.parents` collection so
clonal consumers get extra information without changing the
descendant-list shape.

!!! info "Non-clonal results"
    `result.parents` is `None` when `expand_clones` is NOT in the
    pipeline. Record dicts from non-clonal runs also have
    `record.get("clone_id") is None` and `record.get("parent_id")
    is None`. Code that may receive either shape can safely use
    `result.parents or []` and `record.get("clone_id")` without
    branching on clonal-ness.

## Validation

The validation hub covers the layers in depth; here's the clonal
slice. Three calls compose for the strongest gate:

```python
# 1. Per-record AIRR consistency.
record_report = result.validate_records(refdata)
assert record_report, record_report.summary()

# 2. Within-family invariants (no refdata required).
family_report = result.validate_families()
assert family_report, family_report.summary()

# 3. Parent-aware checks: descendant truth fields agree with their parent.
full_report = result.validate_families_with_parents(refdata)
assert full_report, full_report.summary()
```

What each adds:

- **`validate_families()`** groups records by `clone_id` and asserts
  the recombination-time truth fields (`truth_v_call`, `truth_d_call`,
  `truth_j_call`) are invariant across every descendant of a clone.
  No refdata required; works on records-only results (e.g.
  round-tripped from TSV). On a non-clonal batch (no record carries
  a non-null `clone_id`) it's a safe no-op and returns ok with
  `family_count == 0` — you don't need to branch on whether the
  pipeline was clonal.
- **`validate_families_with_parents(refdata)`** adds structural
  checks (`ParentsMissing` / `ParentIdMissing` / `ParentIdOutOfRange`)
  and value-comparison checks against the actual parent `Outcome`
  on `result.parents` (`ParentDInvertedMismatch`,
  `ParentOriginalVCallMismatch`, `ParentTruthVCallMismatch`,
  etc.). It does NOT re-run the within-family checks
  `validate_families` does; it's a sibling validator, not a
  superset.

The runtime opt-in `run_records(..., validate_records=True)` runs
the per-record validator *and* the field-level family checks — but
NOT the parent-aware checks. If you want
`validate_families_with_parents` in CI, call it explicitly:

```python
result = exp.run_records(seed=42, validate_records=True)
assert result.validate_families_with_parents(refdata), "parent mismatch"
```

## Common workflows

A few patterns that come up repeatedly with clonal output.

**Clone-level SHM benchmark.** Hold the parent V(D)J fixed and
let SHM diverge:

```python
result = (
    ga.Experiment.on("human_igh")
      .recombine()
      .expand_clones(n_clones=100, per_clone=50)
      .mutate(model="s5f", rate=0.05)
      .run_records(seed=42)
)

# Per-clone SHM distribution:
import pandas as pd
df = pd.DataFrame(result.records)
df.groupby("clone_id")["n_mutations"].describe()
```

**Parent / descendant comparison.** Pull the parent IR for each
descendant and compare against the post-SHM sequence:

```python
for rec in result.records:
    parent = result.parents[rec["parent_id"]]
    parent_seq = parent.final_simulation().bases()
    # rec["sequence"] is the post-SHM descendant; parent_seq is
    # the pre-fork assembled IR.
```

**Export records with clone IDs.** Every export format carries the
`clone_id` and `parent_id` columns as-is — they ship with the
record dict:

```python
result.to_tsv("repertoire.tsv")            # clone_id + parent_id columns included
df = result.to_dataframe()                 # same in the DataFrame
result.to_fasta("seqs.fa", prefix="seq")   # headers include sequence_id; clone IDs in TSV
```

**Paired FASTQ from clonal output.** `paired_end` is descendant-
phase, so it composes with `expand_clones` naturally:

```python
result = (
    ga.Experiment.on("human_igh")
      .recombine()
      .expand_clones(n_clones=50, per_clone=20)
      .mutate(model="s5f", rate=0.05)
      .paired_end(r1_length=150, insert_size=300)
      .run_records(seed=1)
)

result.to_paired_fastq("reads_R1.fastq", "reads_R2.fastq")
# 1,000 R1 records + 1,000 R2 records, sequence_ids match
```

## Common mistakes

A handful of issues that show up with clonal pipelines.

**Calling `invert_d()` after `expand_clones()`.** D orientation is
a recombination-time decision; it must be inherited by every
descendant. `invert_d` rejects this immediately at chain time with
"D inversion is a recombination-time decision and must be
inherited by all clone descendants." Move it before
`expand_clones`. Same goes for `receptor_revision`.

**Putting `.paired_end()` before `expand_clones()`.** R1/R2 windows
are per-read; placing `paired_end` in the ancestor phase would
share one R1/R2 window across the whole family. `expand_clones`
scans the prior steps when called and rejects descendant-phase
methods that appear before it: "paired_end must be called after
expand_clones(); it is descendant-specific and must be sampled
independently for each clone member."

**Expecting child traces to include the full parent trace.** A
descendant `Outcome.trace()` carries only the post-fork plan's
addressed choices — SHM substitutions, PCR errors, indels, etc.
The pre-fork plan's trace (recombination choices, NP bases, D
inversion) lives on `result.parents[i].trace()` instead. They're
two separate addressing namespaces because the pre-fork plan ran
once per parent and the post-fork plan ran once per descendant.

**Expecting mutation-distance aggregation fields today.** The
clonal validator's per-clone mutation-distance distribution and
pre-SHM junction invariance checks are deliberately deferred
(future-slice scope). Today's `validate_families` checks the
recombination-time truth fields are invariant within a clone; it
does NOT verify that descendant SHM counts cluster around a
biologically plausible distribution. If you need that, compute it
yourself from `df.groupby("clone_id")["n_mutations"]`.

## Where to go next

- **[Validation & reproducibility](../validation/index.md)** — the
  hub explaining the three validation layers and how the runtime
  opt-in composes with `validate_families`.
- **[SHM and mutation targeting](shm-targeting.md)** — what runs
  per descendant in the SHM pass.
- **[Export the results](../getting-started/export-results.md)** —
  TSV / DataFrame / paired FASTQ formats; `clone_id` ships on
  every record.
- **[The Experiment builder](experiment-builder.md)** — the full
  control-panel page including the canonical ancestor-phase /
  descendant-phase rule.
- For the engine-side mechanics of the plan split and the
  pre-fork / post-fork compile path, see the contributor audit
  [`docs/clonal_plan_split_design.md`](https://github.com/MuteJester/GenAIRR/blob/master/docs/clonal_plan_split_design.md).
