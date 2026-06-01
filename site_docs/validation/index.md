# Validation & reproducibility

<p class="lead">GenAIRR's outputs are rich — 50+ AIRR fields per
record, full event ledgers per outcome, recoverable traces per
seed. Validation answers one question across all of it: does what
the record reports agree with what the engine actually did? This
page is the hub for the three user-facing validation layers, the
replay surface that makes runs durable, and the recommended
workflows for development and CI.</p>

!!! tip "Your learning path"
    You're at the hub of the **"I want reproducible / validated
    output"** path. Two focused deep dives plug in here:
    [`validate_records`](validate-records.md) is the per-record
    output-correctness gate, and
    [Trace, replay, reproducibility](../guides/trace-replay.md)
    is the durable-replay surface. Start with this page for
    orientation, then pick the deep dive that matches your task.
    [See all paths →](../learn.md)

## Why validation exists

The engine drives every AIRR field from internal state — the
persistent IR, the event ledger, and the per-draw trace. The
projection that maps that state into the record runs *outside*
the engine, so a bug in the projection (or in any code path that
feeds the projection — live-call caching, V-subregion attribution,
counter aggregation) could silently produce a record whose fields
don't match the outcome that produced them. Validation closes that
loop: independently re-derive every reported field from the
upstream source, and surface any divergence as a structured
report.

If you ever ask "does this record's `v_call` match the allele the
engine actually sampled?", "do `n_v + n_d + n_j + n_np` add up to
`n_mutations`?", "is `productive` consistent with the actual
junction translation?" — validation is the call that answers,
without you writing the cross-check yourself.

## The three validation layers

GenAIRR exposes three user-facing validation surfaces. All of
them are opt-in — production loops that have qualified their
pipeline don't pay the per-record overhead by default.

### Record validation

```python
report = result.validate_records(refdata)
assert report, report.summary()
```

The point-wise validator. Re-derives every field on every record
from the upstream source and reports per-record divergences in a
`ValidationReport`. This is the gate every release-tier CI run
should carry.

→ See [the dedicated guide](validate-records.md) for the issue
catalogue, the `ValidationReport` API, and the canonical CI
one-liner.

### Family validation

For workloads using `expand_clones(...)`, two sibling validators
check family-level invariants:

```python
family_report = result.validate_families()
assert family_report, family_report.summary()

# Stronger gate — runs per-record validation on every parent first.
full_report = result.validate_families_with_parents(refdata)
assert full_report, full_report.summary()
```

`validate_families` checks within-family invariants alone — every
descendant of a clone shares its V(D)J recombination, every
descendant's `clone_id` matches the parent's, no descendant has a
SHM count below the parent's, and the per-clone size matches what
`expand_clones` was asked for. No refdata required.

`validate_families_with_parents(refdata)` adds full per-record
validation on every parent before checking family invariants —
it's the strongest gate. Use it in release-tier CI when the
pipeline ships clonal output.

### Runtime opt-in

For workloads where you want the validator to run on every batch
without an explicit second call:

```python
result = exp.run_records(n=100, seed=1, validate_records=True)
```

The validator runs inline after each batch; the result is the same
as calling `result.validate_records(refdata)` explicitly. Off by
default — opt in during development; leave off in production
hot loops.

## What each layer catches

| Issue category | Caught by |
|---|---|
| Sequence + coordinate consistency (`sequence` matches pool; `*_start/end` well-ordered) | `validate_records` |
| CIGAR ops (canonical M/I/D/S/N/P/X/=) + query span | `validate_records` |
| Counter provenance (`n_mutations`, per-segment partition, `n_pcr_errors`, indel counts) | `validate_records` |
| Allele calls (`v_call` / `d_call` / `j_call` matches an independent walker) | `validate_records` |
| Junction + productive triad (junction content, `vj_in_frame`, `stop_codon`, `productive` predicate identification) | `validate_records` |
| Paired-end geometry (R1/R2 windows, R2 reverse-complement, `insert_size`) | `validate_records` (when `paired_end()` is in the pipeline) |
| Family size, parent/descendant `clone_id` match, descendant `n_mutations ≥ parent's` | `validate_families` |
| Each family's *parent* validates as a standalone record | `validate_families_with_parents` |

The three layers are independent — `validate_records` doesn't
require a clonal structure; `validate_families` doesn't require
a refdata. They compose for the strongest possible gate:

```python
report = result.validate_records(refdata)
assert report, report.summary()
families = result.validate_families_with_parents(refdata)
assert families, families.summary()
```

## Trace and replay

!!! tip "Full guide"
    See **[Trace, replay, and reproducibility](../guides/trace-replay.md)**
    for the deep dive: seed vs trace, replay vs rerun, every
    failure mode, strict-mode interactions, and recommended
    workflows for debugging, regression tests, and reproducible
    examples. The summary below covers the essentials.

GenAIRR's reproducibility model rests on two facts:

- **Same seed + same plan + same cartridge → byte-identical output.**
  Across runs, machines, and platforms. The `seed=` argument on
  `run_records(...)` is the canonical surface; `n` records use
  seeds `[seed, seed+1, …, seed+n-1]`, so batches stitch together
  if you offset the starting seed.
- **Every random draw is recorded on the outcome's trace.** Each
  draw lives at a stable hierarchical address (e.g.
  `"sample_allele.v"`, `"np.np1.length"`, `"np.np1.bases[3]"`).
  You can inspect the trace, dump it to disk, and replay it later.

For most users, the seed argument is the whole reproducibility
story. The trace+replay surface matters when:

- You want a *durable* replay artifact that captures the cartridge
  identity, engine version, and DSL signature — not just the seed.
- You want to verify, weeks later, that a recorded run *still*
  reproduces against the current code + cartridge.
- You're filing a bug and want a self-contained artifact a
  maintainer can reproduce against.

### Saving a trace file

A `TraceFile` bundles a single outcome's recorded trace together
with the plan signature, refdata signature, refdata content
hash, engine version, and producing seed:

```python
compiled = exp.compile()
outcome = compiled.simulator.run(seed=42)

trace_file = compiled.simulator.trace_file_from(outcome, seed=42)
trace_file.write_to("run-42.trace.json")
```

Read it back later with:

```python
from GenAIRR._engine import TraceFile

trace_file = TraceFile.read_from("run-42.trace.json")
trace_file.to_json()           # round-trip the JSON
trace_file.seed                # 42
trace_file.engine_version      # "X.Y.Z"
trace_file.schema_version      # int
```

### Replaying a trace file

Two complementary replay paths on the simulator:

- **`replay_from_trace_file(trace_file)`** — consumes the recorded
  values verbatim at every sampling slot. The trace becomes the
  source of randomness, not the RNG. This is the strongest
  reproducibility gate: byte-identical output even if the seed's
  RNG sequence drifts in a future version.
- **`rerun_from_trace_file(trace_file)`** — re-runs the sampler
  from the trace's recorded seed. The trace acts as a signature
  bundle (plan + refdata gates) rather than as a value source.
  Useful when you want a fresh draw against the same configured
  pipeline.

```python
compiled = exp.compile()
outcome = compiled.simulator.replay_from_trace_file(trace_file)
```

### Mismatch errors

Replay is gated on three signatures, each of which fires a
`ValueError` if it disagrees with the trace:

| Gate | When it fires | Message prefix |
|---|---|---|
| **Plan signature** | DSL chain changed (different passes / different rates / different kwargs that fold into the signature) | `"pass plan signature mismatch"` |
| **Refdata signature** | Cartridge structure changed (different catalogue / rules) | `"refdata signature mismatch"` |
| **Refdata content hash** | Cartridge bytes changed (curation, V-subregion annotation, allele content) | `"refdata content hash mismatch"` (replay only) |

These fail loudly, before any choices are consumed. If your
replay errors with `"refdata content hash mismatch"`, the trace
was produced against a different cartridge (rules / identity /
curation may differ); load the original cartridge to replay.

## Strict vs permissive

Every `run_records` (and `run`, `stream`, `stream_records`) takes a
`strict=False` keyword that controls what happens when a sampler
runs out of admissible candidates at sample time. This is rare,
but it's the canonical failure mode of an unsatisfiable plan
(e.g. a productive constraint that admits no junction under your
NP-length distribution).

```python
# Default — permissive
result = exp.run_records(n=100, seed=0)

# Strict — fail loud on empty admissible support
result = exp.run_records(n=100, seed=0, strict=True)
```

The two modes diverge only when admissible support is empty:

| Mode | When admissible support is empty | What you see |
|---|---|---|
| `strict=False` (default) | Falls back to a documented sentinel value — indel site `-1`, NP length `0`, NP base `N`, trim `0`; SHM substitution skips the slot. | Execution continues; the record may end up non-productive at that site. |
| `strict=True` | Raises `StrictSamplingError` immediately. | The call fails with `(pass_name, address, reason)` args naming the failing site. |

`StrictSamplingError` is **not** a `ValueError` subclass — `except
ValueError` will not catch it. Catch it by name from the top-level
import:

```python
import GenAIRR as ga

try:
    result = exp.run_records(n=10, seed=42, strict=True)
except ga.StrictSamplingError as e:
    pass_name, address, reason = e.args
    print(f"{pass_name} couldn't satisfy the contract at {address}: {reason}")
```

The `reason` field is a stable lowercase code — common values:
`"empty_admissible_support"`, `"support_unavailable"`,
`"missing_allele.V.42"`, `"contract_violation.NoStopCodonInJunction"`.

The recommended posture: **leave strict off in production**, where
the sentinel fallback keeps your batch flowing; **turn strict on
during cartridge / DSL development**, where you want the
unsatisfiable plan to surface immediately rather than silently.

Note that strict mode applies only to *fresh* sampling — replay
consumes recorded values verbatim, so it never re-evaluates
contract admissibility. To force strict-fresh semantics on a
recorded trace, call `simulator.run(seed=<original_seed>,
strict=True)` instead of `replay_from_trace_file`.

## Recommended workflows

The three layers + the two reproducibility surfaces compose into
three workflow patterns covering the common cases.

### Development

Catch problems immediately. Runtime validation surfaces drift
between your DSL changes and the engine's expectations; strict
mode surfaces unsatisfiable plans:

```python
result = exp.run_records(
    n=100,
    seed=1,
    validate_records=True,
    strict=True,
)
```

If the run succeeds and the inline validation passes, the
pipeline is qualified.

### Continuous integration

Run the validator as an explicit gate, surfacing a structured
report you can keep as a build artifact:

```python
result = exp.run_records(n=1000, seed=0)

report = result.validate_records(refdata)
assert report, report.summary()
```

A failing report dumps the per-record issues to the test runner;
the `failures` list is JSON-serialisable for downstream tooling.

### Clonal output

Stack record + family validation when the pipeline ships clonal
families:

```python
result = (
    exp
    .recombine()
    .expand_clones(n_clones=50, per_clone=20)
    .mutate(model="s5f", rate=0.05)
    .run_records(n=1000, seed=42)
)

assert result.validate_records(refdata), "AIRR record divergence"
assert result.validate_families_with_parents(refdata), "Family invariant divergence"
```

The two layers catch different things and don't substitute for
each other.

## When validation is not enough

GenAIRR's validation answers *is this record internally
consistent with how the engine claims it was produced* — not *is
this output biologically realistic*. A pipeline can validate
clean and still produce simulations that don't match the
biology you're targeting. Three places the validators don't help:

- **Are my cartridge's empirical distributions right?** The
  cartridge's `cartridge_manifest()` block is the canonical
  provenance source for cartridge content. The build report
  attached to a `ReferenceCartridgeBuilder`-produced cartridge
  captures every estimator's inputs and inferred distributions.
- **Are my output distributions calibrated to a reference dataset?**
  Use the distribution-invariant test suite or compare your
  simulated marginals against the dataset you're benchmarking
  against. The audit-realism workflow pattern is the natural
  starting point.
- **Does my simulation match a specific aligner's expectations?**
  GenAIRR's calls come from an independent walker; an aligner
  with different scoring rules will produce different calls on
  the same sequence. That's not a validator problem; it's an
  aligner-comparison problem.

For input-model provenance specifically:

```python
manifest = cfg.cartridge_manifest()
print(manifest["models"]["allele_usage"])     # what's authored on the cartridge
print(manifest["hashes"]["data_config_checksum"])  # canonical content hash

if cfg.build_report is not None:
    for stage in cfg.build_report.stages:
        print(stage["stage"], stage["inputs"])  # what every estimator saw
```

## Where to go next

- **[`validate_records`](validate-records.md)** — the full guide
  to the per-record validator: API, the five issue categories,
  reading a `ValidationReport`.
- **[Your first AIRR record](../getting-started/first-airr-record.md)**
  — the field catalogue the validator checks against.
- **[Reference cartridge](../concepts/reference-cartridge.md)** —
  the cartridge model the validator gates against and the
  manifest that documents input-model provenance.
- **[The Experiment builder](../guides/experiment-builder.md)** —
  how the pipeline composes and where `validate_records=True`
  and `strict=True` sit in the call surface.
