# `validate_records` — AIRR record postcondition validation

<p class="lead">Every AIRR record GenAIRR emits is supposed to be
internally consistent with the outcome that produced it and with the
reference data the engine ran against. <code>validate_records</code>
is the single call that proves it — re-deriving every reported field
from upstream sources and surfacing any divergence as a structured
report.</p>

!!! info "Section hub"
    This page covers the per-record validator in depth. For the
    full validation + reproducibility model (family validation,
    trace files, replay, strict mode, recommended workflows), see
    [**Validation & reproducibility**](index.md).

## What validation answers

> Given a finished `SimulationResult` and the `RefDataConfig` that
> drove it, is every reported field on every record internally
> consistent and biologically derivable?

The validator does not check whether your simulation is biologically
*realistic* — that's a question about your cartridge and your
parameters. It checks whether the engine's output is *self-consistent*:
that the `v_call` matches what an independent walker would assign,
that `n_mutations` matches the event ledger, that `junction_length`
matches `len(junction)`, that the `productive` flag corresponds to a
real in-frame-no-stop-anchors-preserved evaluation, and so on.

## Quick start

After any `run_records(...)` call, run the validator over the
returned `SimulationResult`:

```python
import GenAIRR as ga

refdata = ga.dataconfig_to_refdata(ga.HUMAN_IGH_OGRDB)

result = (
    ga.Experiment.on("human_igh")
      .recombine()
      .productive_only()
      .run_records(n=100, seed=1)
)

report = result.validate_records(refdata)
assert report                     # ValidationReport is truthy iff every record passed
print(report.summary())           # Histogram of issue kinds (empty when ok)
```

`refdata` is the cartridge in its engine-side `RefDataConfig` shape.
When you bound the `Experiment` to a string like `"human_igh"` or a
`DataConfig`, the bridge produces a `RefDataConfig` for the engine;
you can obtain the same object via `ga.dataconfig_to_refdata(cfg)`
on the Python side, or read it off
`exp.compile().simulator.refdata` after a compile.

## Runtime validation

For workloads where you want the validator to run on every batch
without an explicit second call, opt in at simulation time:

```python
result = (
    ga.Experiment.on("human_igh")
      .recombine()
      .run_records(n=100, seed=1, validate_records=True)
)
# The validator already ran. result.validate_records(refdata) re-runs
# if you call it explicitly; the runtime flag is just a convenience.
```

The runtime flag is off by default — validation is intentionally an
*opt-in* discipline so production simulation loops that have already
qualified their parameters don't pay the per-record overhead.

## Reading a `ValidationReport`

`validate_records` returns a `ValidationReport`. The dataclass is
designed to be both immediately readable and CI-assertable:

```python
report = result.validate_records(refdata)

report.ok            # True / False — short-circuit gate for CI
report.count         # int — total records validated
report.failures      # list[dict] — one entry per failing record
report.summary()     # str — single-line histogram of issue kinds
bool(report)         # equivalent to report.ok

assert report, report.summary()   # the canonical CI one-liner
```

When `ok` is `False`, each entry of `failures` carries:

```python
{
    "record_index": 17,
    "sequence_id": "seq17",
    "issues": [
        {"kind": "JunctionLength", "details": {...}},
        {"kind": "VCallMismatch",  "details": {...}},
    ],
}
```

`issues` is a list because a single record can fail multiple
invariants — the validator never short-circuits on the first
problem.

## Common issue categories

The validator catalogues divergences into five families. Each family
groups together the checks that share the same upstream source, so
when you see one fire you can localise the suspect code path
quickly.

### Sequence + coordinate consistency

- `sequence_length` matches `len(record.sequence)`.
- `record.sequence` matches the simulation pool (case-insensitive
  because sequencing-error corruption lowercases bytes).
- V / D / J `*_sequence_start/end` and `*_germline_start/end` pairs
  are well-ordered (`start ≥ 0`, `end ≥ start`).
- CIGAR strings parse with only canonical M / I / D / S / N / P /
  X / = ops.
- CIGAR query span matches the segment's reported sequence span.

### Counter provenance

Every counter on the record is re-derived from a canonical source:

| Counter | Source |
|---|---|
| `n_mutations` | `Simulation.mutation_count` — set by S5F / Uniform at seal time |
| `n_v_mutations` / `n_d_mutations` / `n_j_mutations` / `n_np_mutations` | The same SHM events, partitioned by carried segment |
| `n_pcr_errors` / `n_quality_errors` | Trace addresses on the corruption passes |
| `n_indels` / `n_v_indels` / `n_d_indels` / `n_j_indels` | `IndelInserted` + `IndelDeleted` events from the indel pass |
| `end_loss_5_length` / `end_loss_3_length` | Trace addresses on the end-loss pass |

A `MutationCountSumMismatch` issue means the per-segment counters
don't add up to `n_mutations` — almost always a sign that a
mechanism added events to a segment the partition didn't anticipate.

### Junction + productivity

- `junction` content matches the recomputed pool slice from the V
  anchor through the J anchor + 3.
- `junction_length` matches `len(junction)`.
- `vj_in_frame` matches `junction_length % 3 == 0`.
- `stop_codon` matches `junction_has_stop(junction)` (only meaningful
  when in-frame).
- `productive` matches the full triad: in-frame ∧ no junction stop ∧
  V Cys preserved ∧ J W-or-F preserved.

When `productive=False`, the issue payload also names *which*
predicate fired (`OutOfFrame` / `JunctionStopCodon` /
`VAnchorAaChanged` / `JAnchorAaChanged`), so a "non-productive
storm" is immediately diagnosable.

### Allele calls

The validator rescores `v_call` / `d_call` / `j_call` against an
independent walker — same matching rules, different implementation
path. Mismatches surface the engine-reported and oracle-reported
tie-sets side-by-side so you can see whether the engine picked
the wrong allele or just ordered the tie-set differently.

For rev-comp records the C4 oracle is skipped by design (the
post-projection flip moves the bytes out from under the oracle's
reference index). The other four categories still run.

### Paired-end layout

When `Experiment.paired_end(...)` is part of the pipeline, the
record carries eight extra fields (`r1_sequence`, `r2_sequence`,
`r1_start/end`, `r2_start/end`, `insert_size`, `read_layout`). The
validator re-derives them from the trace and confirms the R1/R2
windows are sliced correctly, R2 is correctly reverse-complemented,
and the insert size matches.

## Family validation

For workloads using `expand_clones(...)` to generate clonal
families, two sibling validators check family-level consistency:

```python
result = (
    ga.Experiment.on("human_igh")
      .recombine()
      .expand_clones(n_clones=50, per_clone=20)
      .mutate(rate=0.05)
      .run_records(n=1000, seed=42)
)

family_report = result.validate_families(refdata)
assert family_report, family_report.summary()
```

`validate_families` checks within-family invariants: every
descendant of a clone shares its V(D)J recombination, every
descendant's `clone_id` matches the parent, no descendant has a
SHM count below the parent's, and the per-clone size matches what
`expand_clones` was asked for.

When you also want to check that each family's *parent* record
itself validates against the cartridge (in addition to family
invariants), use the parent-aware variant:

```python
report = result.validate_families_with_parents(refdata)
```

This is the most expensive of the three — it runs the per-record
validator on every parent before checking family invariants — but
it's the strongest gate. Use it in release-tier CI.

## What `validate_records` does NOT do

A few things the validator deliberately does not check, so you know
what to reach for when you need them:

- **It is not a biological truth oracle for real data.** The
  validator answers "is this record internally consistent with how
  the engine claims it was produced". It can't tell you whether your
  parameters are biologically realistic. For that you compare the
  output against published distributions or use the
  `audit-realism` recipe.
- **It is not automatic unless you opt in.** Either pass
  `validate_records=True` to `run_records`, or call
  `result.validate_records(refdata)` explicitly. Production loops
  that have already qualified their pipeline don't pay the cost.
- **It is per-record, not per-batch.** Distribution invariants
  ("V-gene usage is approximately uniform") are statistical
  guarantees over a batch; `validate_records` is point-wise. Use the
  distribution-invariant test suite for batch-level checks.
- **Cache parity is a separate, internal layer.** GenAIRR's
  release-tier CI runs a second integrity check called
  `check_live_call_cache_parity` on every outcome, comparing the
  cached `SegmentLiveCall` against a from-scratch recompute. That
  layer is internal to engine maintenance — users rarely need it
  unless they're filing a bug against the cache. Both layers are
  documented in the [two-layer integrity model](two-layer-model.md).

---

## Deep architecture notes

The validator's implementation lives in
[`engine_rs/src/airr_record/validate.rs`](https://github.com/MuteJester/GenAIRR/blob/master/engine_rs/src/airr_record/validate.rs)
and is exposed to Python through `SimulationResult.validate_records`.
For the full check catalogue, the §5A/§5B silent invariants the
validator pins, and the empirical sweep that drove the original
landing, see the contributor audit at
[`docs/airr_record_validator.md`](https://github.com/MuteJester/GenAIRR/blob/master/docs/airr_record_validator.md).
For the engine-wide validation matrix (every guarantee → audit doc
→ test file → Rust kernel mapping), see
[`docs/validation_matrix.md`](https://github.com/MuteJester/GenAIRR/blob/master/docs/validation_matrix.md).
