# Trace, replay, and reproducibility

<p class="lead">A trace is the durable record of every choice the
engine sampled during one simulation. With a trace on disk you can
replay a run byte-for-byte, rerun the same plan from the same seed,
or hand a colleague a reproducible example without shipping
gigabytes of FASTQ. This guide walks through the surface and the
failure modes you actually hit in practice.</p>

## Seed vs trace

GenAIRR has two layers of reproducibility, and they answer
different questions.

**Seed reproducibility.** Run the same code, against the same
config, with the same seed — and the engine's RNG produces the
same sequence of draws. Seeds are tiny, immortal, and require
the *exact* same environment to mean anything. Change a single
default in your `Experiment` and the same seed produces a
different outcome.

**Trace reproducibility.** A trace is the **record** of choices
the engine made during one run. It includes the seed, but it
also includes the resolved plan signature, the cartridge
signature, the content hash, and the actual value at every
address the engine touched. With a trace you can replay an
outcome even if the environment has shifted slightly — as long
as the plan and refdata signatures match, every recorded value
is consumed verbatim.

The mental model:

| You want to ... | Use |
|---|---|
| Hand someone "run this and you get my batch" | **Seed**, when code is shared |
| Reconstruct one specific record from a paper / bug report | **Trace** |
| Run a regression test that survives a non-breaking refactor | **Trace** |
| Verify a code change didn't perturb sampling | **Seed** before / **trace** after |

A seed is a promise about RNG output. A trace is a record of
sampling. Use both.

## Run once, save a trace

`compile()` once, run with a seed, and ask the simulator for the
trace that produced the resulting outcome:

```python
import GenAIRR as ga
from GenAIRR._engine import TraceFile

exp = ga.Experiment.on("human_igh").recombine().mutate(rate=0.02)
compiled = exp.compile()

outcome = compiled.simulator.run(seed=1)
trace_file = compiled.simulator.trace_file_from(outcome, seed=1)
trace_file.write_to("example.trace.json")
```

A few things to know:

- **`TraceFile` is an engine-level export.** It lives at
  `GenAIRR._engine.TraceFile`, not at top-level `ga.TraceFile`.
  The PyO3 type is stable; the placement on the public namespace
  is intentional — traces are a *contract* with the engine, not
  with the Python wrapper.
- **`trace_file_from(outcome, seed)`** packages one outcome + its
  seed into a trace. Pass the same seed you used for the run — the
  trace's `seed` field carries it forward for inspection.
- **`write_to(path)`** writes the trace as JSON. The file is plain
  text, gzip-able, diffable, and small (one entry per sampled
  choice plus the four signature fields below).

## Replay exactly

Read the trace back and replay it against any compiled experiment
whose plan and refdata signatures match:

```python
trace_file = TraceFile.read_from("example.trace.json")
replayed = compiled.simulator.replay_from_trace_file(
    trace_file,
    strict=False,
)
```

`replay_from_trace_file` consumes the recorded values verbatim.
The simulator walks the plan, hits every address in the same
order, and reads the corresponding value out of the trace instead
of asking the RNG. The returned `Outcome` is byte-identical to
the original (modulo non-determinism in code paths that aren't
sampling — there shouldn't be any).

## Rerun from trace

`rerun_from_trace_file` is the seed-based sibling of replay:

```python
rerun = compiled.simulator.rerun_from_trace_file(
    trace_file,
    strict=False,
)
```

The difference matters:

| Operation | What the engine does at each address |
|---|---|
| `replay_from_trace_file` | Reads the recorded value verbatim from the trace |
| `rerun_from_trace_file` | Re-draws from the RNG using the trace's seed and plan |

Both walk the **same** address sequence (gated by the same plan
signature). Replay reproduces the *exact* outcome from the trace.
Rerun reproduces a *fresh* outcome that should match the original
when the environment hasn't drifted — and exposes any drift
loudly when it has, since the re-drawn values differ from the
trace.

A practical mental rule:

- **Replay** answers "what did happen?"
- **Rerun** answers "what would happen now?"

A regression test that runs both and asserts they match is the
strongest reproducibility check you can author with one trace.

## What the trace contains

A serialised trace carries eight fields:

| Field | Purpose |
|---|---|
| `schema_version` | The trace-file schema version |
| `engine_version` | The engine version that wrote the trace |
| `seed` | The seed passed to the original `run(...)` |
| `pass_plan_signature` | Hash of the resolved pass plan (every method-call signature on the `Experiment`) |
| `refdata_signature` | Identity hash of the refdata (cartridge identity + catalogue shape) |
| `address_schema_version` | The schema version of the choice-address encoding |
| `refdata_content_hash` | Content hash of every plane the trace depends on |
| `trace` | The ordered list of `(address, value)` choice records |

What's deliberately **not** in the trace:

- The full reference cartridge bytes (only its signature + content hash)
- The assembled sequences from the original outcome
- The compiled pass results (these are derived during replay)
- The contract set used at compile time
- The execution policy (strict / permissive) — that's a replay-time argument

The trace is small because the cartridge is its provenance, not
its payload. You ship the trace + the cartridge name; the engine
re-derives everything else from the recorded choices.

## What can make replay fail

Replay raises `ValueError` (never `KeyError`, never a custom
exception) with one of these distinct shapes:

### Validation phase — checked before any choices are consumed

| Failure mode | When it fires |
|---|---|
| `plan_signature_mismatch` | The compiled experiment's plan hash doesn't match `pass_plan_signature` |
| `refdata_signature_mismatch` | The cartridge's identity hash doesn't match `refdata_signature` |
| `refdata_content_hash_mismatch` | A plane referenced by the trace has changed bytes |

These three are the load-bearing gates. A `plan_signature_mismatch`
means the pipeline differs — somebody added an `.invert_d()`
between recording and replaying, or changed a `mutate(rate=...)`
value. A `refdata_signature_mismatch` means the cartridge identity
differs — you're trying to replay an `HUMAN_IGH_OGRDB` trace
against `HUMAN_IGH_EXTENDED`. A `refdata_content_hash_mismatch`
is the subtle one: the cartridge name is the same but a plane
referenced by the trace has changed — usually because someone
re-estimated a model and overwrote the cartridge in place.

### Loading phase — checked when the trace file is parsed

| Failure mode | When it fires |
|---|---|
| `address_schema_version_mismatch` | The trace was written against an older or newer address encoding |

This fires before replay even starts. Old traces eventually become
unreadable when the address schema bumps; the engine version field
in the file tells you which release wrote it.

### Execution phase — checked while choices are consumed

| Failure mode | When it fires |
|---|---|
| `trace_exhausted` | The replay walked past the last recorded choice |
| `address_mismatch` | The next address the engine wants doesn't match the next recorded address |
| `value_kind_mismatch` | The recorded value type doesn't match what the address expects (e.g. an integer where a base draw is required) |
| `unused_trailing_records` | Replay finished but the trace still had entries left |

These four don't usually fire under normal use — they catch
*real* engine-level corruption (a hand-edited trace, a refactor
that changed address ordering without bumping the schema). The
three validation-phase gates fire far more often in day-to-day
work.

## Strict mode and replay

`strict=True` on a fresh `run(...)` raises `StrictSamplingError`
when an admissible-support gate fires — for example, no
productive-safe NP composition exists for the current draws.

When a trace was **recorded** under strict mode:

- The strict gate either fired (in which case the run never
  produced an outcome to trace) or didn't (in which case the
  trace records only the successful choices, including any
  recovery sentinels).
- **Replay** consumes those recorded values verbatim. The
  `strict` flag on `replay_from_trace_file` doesn't change what
  values get returned — replay always returns what was recorded.
  The flag controls how the replay engine handles trailing /
  exhausted conditions described above.
- **Rerun** does re-sample, so a `rerun_from_trace_file(...,
  strict=True)` can *newly* fire `StrictSamplingError` even if
  the original trace ran clean, because the RNG drew a different
  value at one address.

The contract: replay returns what was; rerun returns what would
be. Strict mode changes "what would be" but not "what was."

## Recommended workflows

A few patterns that come up in practice.

### Debugging one surprising record

Save a trace alongside the AIRR output:

```python
result = exp.run_records(n=1000, seed=42)
compiled = exp.compile()
outcome = compiled.simulator.run(seed=42)
trace_file = compiled.simulator.trace_file_from(outcome, seed=42)
trace_file.write_to("debug-record-42.trace.json")
```

You can hand that one JSON file plus the cartridge name to a
colleague and they can replay it on their machine — no
multi-gigabyte FASTQ exchange required.

### Sharing reproducible examples

Ship traces with bug reports. The maintainer reads the trace,
inspects choices at each address, and reproduces the bug
deterministically.

### Regression tests

A trace pinned at one engine version + cartridge revision is the
strongest regression test you can write. The test loads the
trace, replays it, and asserts the outcome matches a stored
snapshot. Any code path change that touches sampling fires one
of the three validation gates loudly:

```python
def test_replay_unchanged_record():
    trace = TraceFile.read_from("tests/fixtures/canonical-record.trace.json")
    out = compiled.simulator.replay_from_trace_file(trace)
    assert out.assembled_sequence == EXPECTED_SEQUENCE
```

When the test fails with `plan_signature_mismatch`, your refactor
moved sampling around. When it fails with
`refdata_content_hash_mismatch`, somebody re-estimated the
cartridge. Either way the failure mode tells you what changed.

### Storing traces alongside AIRR output

If you publish a dataset, write `record_id → trace_path` next to
the AIRR table. Readers who want to reconstruct a specific record
load the trace, replay it against the published cartridge, and
get the exact same row back. The trace is small enough that
shipping one per record is feasible for small / curated batches;
for large datasets, ship traces only for the records you've
annotated.

## Common mistakes

A handful of issues that show up repeatedly with the trace
surface.

**Changing a parameter and expecting replay to work.** Any change
to the `Experiment` pipeline that affects the plan signature
breaks replay — adding a pass, removing a pass, changing
`mutate(rate=0.02)` to `mutate(rate=0.03)`, even toggling
`productive_only()`. The `plan_signature_mismatch` failure is
loud on purpose. Use **rerun** if you want to see what would
change; use **replay** when you want to verify nothing did.

**Using a different cartridge.** Replaying a `HUMAN_IGH_OGRDB`
trace against `HUMAN_IGH_EXTENDED` fires
`refdata_signature_mismatch`. The two cartridges have different
catalogues, different allele names at the same index — replay
would silently corrupt. The gate catches it.

**Confusing replay with rerun.** Replay consumes recorded values;
rerun re-draws. Calling `replay_from_trace_file` and expecting
"the same plan but with a different RNG path" is a category
error — pass the original trace to `rerun_from_trace_file` if you
want fresh draws under the same plan.

**Expecting the trace to contain the full sequence or cartridge.**
It doesn't. The trace is a sequence of `(address, value)` choices
plus four signature fields. To reproduce an outcome you need the
trace AND the cartridge that the trace's `refdata_signature`
identifies. Ship them together or rely on the canonical bundled
cartridges (`ga.HUMAN_IGH_OGRDB`, etc.) for shareable examples.

**Hand-editing a trace JSON file.** The trace is human-readable
but not human-editable. Changing a value will fire
`value_kind_mismatch` or `address_mismatch` during execution;
changing a signature will fire one of the validation gates. If
you want to perturb a record, run a fresh simulation with a new
parameter and trace *that*.

## Where to go next

- **[Validation hub](../validation/index.md)** — the broader
  picture of GenAIRR's reproducibility and validation guarantees.
- **[`validate_records`](../validation/validate-records.md)** —
  per-record AIRR-output gate that pairs cleanly with
  replay-pinned regression tests.
- **[The Experiment builder](experiment-builder.md)** — the
  pipeline whose plan signature replay's first gate checks.
- **[Reference cartridge](../concepts/reference-cartridge.md)** —
  the four-plane model whose identity + content hashes ride into
  the trace's validation gates.
