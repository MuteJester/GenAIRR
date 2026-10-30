# `validate_records` - AIRR record postcondition validation

<p class="lead">Every AIRR record GenAIRR emits is supposed to be
internally consistent with the outcome that produced it and with the
reference data the engine ran against. <code>validate_records</code>
checks supported fields against those upstream sources and surfaces
detected inconsistencies as a structured report.</p>

!!! info "Section hub"
    This page covers the per-record validator in depth. For the
    full validation + reproducibility model (family validation,
    trace files, replay, strict mode, recommended workflows), see
    [**Validation & reproducibility**](index.md).

## What validation answers

> Given a finished `SimulationResult` and the `RefDataConfig` that
> drove it, do the implemented consistency checks agree with the
> reported records?

The validator does not check whether your simulation is biologically
*realistic* - that's a question about your cartridge and your
parameters. It checks whether the engine's output is *self-consistent*:
that `v_call` matches the generating assignment and `recalled_v_call`
matches independent scoring of the observed bases,
that `n_mutations` matches the event ledger, that `junction_length`
matches `len(junction)`, that the `productive` flag corresponds to a
coding-frame, stop-codon and anchor checks, and so on. Passing these
checks is not a universal guarantee of biological validity.

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

The runtime flag is off by default - validation is intentionally an
*opt-in* discipline so production simulation loops that have already
qualified their parameters don't pay the per-record overhead.

## Reading a `ValidationReport`

`validate_records` returns a `ValidationReport`. The dataclass is
designed to be both immediately readable and CI-assertable:

```python
report = result.validate_records(refdata)

report.ok            # True / False - short-circuit gate for CI
report.count         # int - total records validated
report.failures      # list[dict] - one entry per failing record
report.summary()     # str - single-line histogram of issue kinds
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
invariants - the validator never short-circuits on the first
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
| `n_mutations` | `Simulation.mutation_count` - set by S5F / Uniform at seal time |
| `n_v_mutations` / `n_d_mutations` / `n_j_mutations` / `n_np_mutations` | The same SHM events, partitioned by carried segment |
| `n_pcr_errors` / `n_quality_errors` | Each pass's own `BaseChanged` events, counting only writes that changed a base |
| `n_indels` / `n_v_indels` / `n_d_indels` / `n_j_indels` | `IndelInserted` + `IndelDeleted` events from the indel pass |
| `end_loss_5_length` / `end_loss_3_length` | Trace addresses on the end-loss pass |

A `MutationCountSumMismatch` issue means the per-segment counters
don't add up to `n_mutations` - almost always a sign that a
mechanism added events to a segment the partition didn't anticipate.

### Junction + productivity

- `junction` content matches the recomputed pool slice from the V
  anchor through the J anchor + 3.
- `junction_length` matches `len(junction)`.
- `junction_aa` matches translation of complete codons in the physical
  junction, including when the rearrangement is out of frame.
- `vj_in_frame` requires a junction length divisible by three and
  alignment of its start with the retained V coding frame.
- `stop_codon` checks the retained coding sequence in that frame,
  including sequence outside the junction.
- `productive` also requires the retained J tail to preserve the
  coding frame, no coding stop, V Cys and J W-or-F anchors.

A productivity mismatch includes the deciding predicate
(`OutOfFrame`, `CodingStopCodon`, `VAnchorAaChanged`, or
`JAnchorAaChanged`).

### Allele calls

The validator checks `v_call` / `d_call` / `j_call` against the final
generating assignments. It separately rescores `recalled_v_call` /
`recalled_d_call` / `recalled_j_call` using the same matching rules through
an independent scoring path. In genotype runs, candidate alleles come from
the donor's two haplotypes. Re-call mismatches report both ordered tie sets.
The scoring oracle abstains when no germline evidence remains.

These checks also run for reverse-strand records: annotations and the
internal sequence used for scoring remain in coding orientation.

### Paired-end layout

When `Experiment.paired_end(...)` is part of the pipeline, the
record carries eight extra fields (`r1_sequence`, `r2_sequence`,
`r1_start/end`, `r2_start/end`, `insert_size`, `read_layout`). The
validator re-derives them from the trace and confirms the R1/R2
windows are sliced correctly, R2 is correctly reverse-complemented,
and the insert size matches.

## Family validation

For workloads that stamp `clone_id` (`clonal_lineage(...)`,
`clonal_repertoire(...)`, or legacy `expand_clones(...)`), family
validation checks clone-level consistency:

```python
result = (
    ga.Experiment.on("human_igh")
      .recombine()
      .clonal_repertoire(n_clones=50, max_size=100)
      .mutate(rate=0.01)
      .run_records(seed=42, expose_provenance=True)
)

family_report = result.validate_families()
assert family_report, family_report.summary()
```

`validate_families` is records-only. It checks that a clonal batch is
not mixed with non-clonal records and, when `truth_v_call`,
`truth_d_call`, and `truth_j_call` are present, that those
recombination-time truth calls are invariant within each `clone_id`
group.

For legacy `expand_clones(...)`, you can also compare each
descendant against its actual parent `Outcome`:

```python
report = result.validate_families_with_parents(refdata)
```

`clonal_repertoire` and `clonal_lineage` do not expose
`result.parents`, so `validate_families_with_parents` is not their
primary validator. For `clonal_lineage`, validate the tree objects
directly with `tree.validate()` when topology matters.

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
  layer is internal to engine maintenance - users rarely need it
  unless they're filing a bug against the cache. Both layers are
  documented in the [architecture guide](../architecture/index.md#the-two-layer-integrity-model).

---

## Deep architecture notes

The validator's implementation lives in
[`engine_rs/src/airr_record/validate.rs`](https://github.com/MuteJester/GenAIRR/blob/master/engine_rs/src/airr_record/validate.rs)
and is exposed to Python through `SimulationResult.validate_records`.
For the full check catalogue, the §5A/§5B silent invariants the
validator pins, and the empirical sweep that drove the original
landing, see the contributor audit at
[`audit-docs/airr_record_validator.md`](https://github.com/MuteJester/GenAIRR/blob/master/audit-docs/airr_record_validator.md).
For the engine-wide validation matrix (every guarantee → audit doc
→ test file → Rust kernel mapping), see
[`audit-docs/validation_matrix.md`](https://github.com/MuteJester/GenAIRR/blob/master/audit-docs/validation_matrix.md).
