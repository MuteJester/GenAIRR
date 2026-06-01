# Estimate cartridge models from real data

<p class="lead">A reference cartridge ships with biology baked in — but
the empirical models inside (allele usage, trim distributions, NP
lengths, NP base composition, P-nucleotide lengths) come from
data. GenAIRR's estimator surface lets you turn a real AIRR
rearrangement table into a cartridge whose simulator output
matches the empirical distributions in that table. This guide is
the focused workflow for that loop.</p>

## What estimators do

An estimator reads a stream of AIRR-shaped rearrangement records
and writes one typed empirical model onto the cartridge's
`reference_models` plane:

```text
records  ──estimator──►  cfg.reference_models.<plane>
```

The output is still a normal `DataConfig`. The cartridge's
identity (`cfg.identity`), allele catalogue (`cfg.alleles`),
rules (`cfg.rules`), and any models you didn't estimate stay
exactly as they were. Estimators are **typed-plane mutators**:
they touch one plane each, leave the rest alone, and stamp an
auditable trail on the build report.

## Which estimators ship today

Five estimators are available in v1:

| Estimator | Plane it populates | Required columns | Optional columns |
|---|---|---|---|
| `estimate_allele_usage` | `reference_models.allele_usage` | `v_call`, `j_call` | `d_call` (required on VDJ) |
| `estimate_trim_distributions` | `reference_models.trims` | `v_trim_3`, `j_trim_5` | `d_trim_5`, `d_trim_3` (required on VDJ) |
| `estimate_np_length_distributions` | `reference_models.np_lengths` | `np1_length` | `np2_length` (required on VDJ) |
| `estimate_np_base_model` | `reference_models.np_bases` | `np1` | `np2` (required on VDJ) |
| `estimate_p_nucleotide_lengths` | `reference_models.p_nucleotide_lengths` | `p_v_3_length`, `p_j_5_length` | `p_d_5_length`, `p_d_3_length` (required on VDJ) |

!!! warning "`estimate_shm_rates` is deferred"
    SHM-model estimation is not in v1. Calling
    `builder.estimate_shm_rates(...)` raises `AttributeError`. SHM
    kernels still ship from the cartridge author's chosen S5F
    table; for now you can swap them at build time but not estimate
    them from records.

## Input record requirements

Estimators consume one of three input shapes:

- **List of dicts**: `[{"v_call": "IGHV3-23*01", ...}, ...]`
- **AIRR TSV path**: `"my_repertoire.tsv"` (tab-delimited)
- **Open text file handle**: an iterable yielding TSV lines

The expected column names by estimator:

| Estimator | Columns it reads |
|---|---|
| `estimate_allele_usage` | `v_call`, `d_call`, `j_call` |
| `estimate_trim_distributions` | `v_trim_3`, `d_trim_5`, `d_trim_3`, `j_trim_5` |
| `estimate_np_length_distributions` | `np1_length`, `np2_length` |
| `estimate_np_base_model` | `np1`, `np2` |
| `estimate_p_nucleotide_lengths` | `p_v_3_length`, `p_d_5_length`, `p_d_3_length`, `p_j_5_length` |

Columns are read **per-field**, not per-row: a row missing one
column for one estimator is rejected only for that estimator's
stage, not for the whole pipeline. Run all five estimators in
sequence against the same record set — each takes what it needs
and reports independently on the rest.

VJ-chain cartridges (kappa, lambda) silently ignore the D-segment
columns; VDJ-chain cartridges require them.

## One complete estimator chain

The full builder workflow, with every estimator wired up:

```python
import GenAIRR as ga

records = [...]              # list[dict] of AIRR rearrangements

builder = (
    ga.ReferenceCartridgeBuilder
    .from_fasta(
        v_fasta="v.fasta",
        d_fasta="d.fasta",
        j_fasta="j.fasta",
        chain_type="BCR_HEAVY",
    )
    .infer_identity(species="HUMAN", locus="IGH")
    .estimate_allele_usage(records)
    .estimate_trim_distributions(records)
    .estimate_np_length_distributions(records)
    .estimate_np_base_model(records, kind="markov", pseudocount=1.0)
    .estimate_p_nucleotide_lengths(records)
)

cfg = builder.build()
report = builder.report()
```

Every call returns the same `ReferenceCartridgeBuilder` so the
chain stays fluent. Order doesn't matter for correctness — each
estimator writes its own plane — but the build report records the
calls in the order they ran, so put them in the order a future
reader will find easiest to follow.

`build()` runs cartridge finalisation (alphabet, identity, rules
defaulting), stamps a `schema_sha256`, and returns the
`DataConfig`. After `build()`, the builder is frozen; further
estimator calls raise.

## Ambiguity and skipped rows

Estimators **don't crash on imperfect data**. Every malformed,
ambiguous, or unknown-allele row is captured in
`report.rejected` with the reason it was dropped, so you can
audit data quality without losing the rest of the corpus.

The recognised rejection reasons:

| Reason | What triggers it |
|---|---|
| `missing_required_column` | The field was empty / missing on this row |
| `unknown_allele` | The allele name isn't in the cartridge's catalogue |
| `missing_d_call_on_vdj` | VDJ chain row with no `d_call` |
| `malformed_trim_value` | Trim column couldn't parse as integer |
| `malformed_length_value` | NP / P length column couldn't parse as integer |
| `negative_trim_value` | Trim value parsed but was negative |
| `negative_length_value` | NP / P length parsed but was negative |
| `noncanonical_base` | NP sequence contained characters outside `{A,C,G,T}` |

Each rejection entry is a dict carrying at minimum `stage`,
`row_index`, and `reason`. Stage-specific entries add fields
like `column`, `segment`, `allele_name`, `value`, or
`unknown_chars`. Inspect them post-build:

```python
report = builder.report()
for entry in report.rejected:
    if entry["reason"] == "unknown_allele":
        print(entry["stage"], entry["segment"], entry["allele_name"])
```

### `estimate_allele_usage` ambiguity policy

The `ambiguous` keyword chooses how comma-separated tie sets
(`"IGHV3-23*01,IGHV3-23*04"`) are handled:

| Policy | Behaviour |
|---|---|
| `"fractional"` (default) | Splits credit evenly across the known alleles in the tie set. Unknown alleles inside a tie set surface as `unknown_allele` rejection entries, one per segment per row, without dropping the row if at least one allele in the set is known. |
| `"truth_first"` | Collapses to the first allele in the tie set before resolution — closest to "treat the call as a single allele". |
| `"reject"` | Drops the entire row if any segment has multiple alleles. Useful when you want to estimate from confident calls only. |

### P-nucleotide naive input detector

`estimate_p_nucleotide_lengths` emits a stage-level warning when
≥95 % of contributing rows reported zero for a given key. That
threshold is the **P-naive detection signal**: most external AIRR
tools don't report palindromic insertions, so a tsv from MiXCR /
IMGT / AIRR-tools that doesn't have those columns will trip the
warning and you'll know the estimated distribution is essentially
"always zero" rather than meaningful biology.

## How estimated models affect simulation

The empirical models you write to the cartridge feed the engine's
recombination passes in the obvious way:

| Plane | What it drives |
|---|---|
| `allele_usage` | Recombination's V / D / J segment-pick weights |
| `trims` | Per-segment trim distributions at recombine time |
| `np_lengths` | NP1 / NP2 insertion length draws |
| `np_bases` | NP1 / NP2 base composition (first-base + Markov transitions on `kind="markov"`) |
| `p_nucleotide_lengths` | Palindromic-insertion length draws |

Two important composition rules:

- **Explicit `Experiment` kwargs always override the cartridge.**
  When you pass `.trim(v_3=..., d_5=..., d_3=..., j_5=...)` to an
  experiment, those overrides win — the estimated cartridge
  distribution is the *default*, not a floor.
- **The plan signature folds these distributions** where applicable,
  so replaying a trace against a cartridge whose estimated
  distributions have drifted will fail cleanly with a plan-signature
  mismatch. The one documented soft gap is in allele-usage
  weights — see `docs/replay_protocol.md` for the precise
  envelope.

The cartridge itself doesn't expose run-time sampling code; the
engine consumes the typed planes through `dataconfig_to_refdata`
and translates them into compiled recombination draws at
`Experiment.compile()` time.

## Inspecting the result

Two readable handles after `build()`:

### The build report

```python
report = builder.report()
report.stages              # list[dict] — one entry per builder call
report.warnings            # list[str]  — build-finalisation warnings
report.rejected            # list[dict] — per-row + per-allele drops
report.manifest_snapshot   # dict | None — cfg.cartridge_manifest() at build time
report.checksum_at_build_time  # str | None — schema_sha256 stamped on cfg
report.to_dict()           # JSON-clean dict for CI artifacts
```

`stages` carries an `inputs` and `inferred` block per stage. For
estimator stages you'll find counts of `records_seen`,
`records_used`, `replaced` (idempotency marker), and per-key
diagnostic dicts.

### The cartridge manifest

```python
manifest = cfg.cartridge_manifest()
manifest["models"]
# {
#   "has_reference_models": True,
#   "allele_usage": {"v": 320, "d": 30, "j": 6, ...},
#   "trim_models": {"v_trim_3": {"support_size": 18, ...}, ...},
#   "np_length_models": {"np1": {...}, "np2": {...}},
#   "np_base_models": {"first_base": {...}, "transition": {...}},
#   "p_nucleotide_models": {"p_v_3_length": {...}, ...},
#   "shm": {"v_subregion_support": {"available": True, ...}, ...},
#   "trim_keys": ["v_trim_3", "d_trim_5", ...],
#   "np_length_keys": ["np1", "np2"],
#   "legacy_trim_dicts_present": False,
#   "legacy_np_lengths_present": False,
# }
```

`manifest["models"]` is the canonical place to verify what got
estimated. Each top-level key is `None` when its plane is
empty, so `manifest["models"]["allele_usage"] is None` after a
build that skipped `estimate_allele_usage` is normal.

## Common mistakes

A handful of issues that show up repeatedly with the estimator
surface.

**Expecting legacy single-dict fields to auto-lift.** The old
`cfg.trim_dicts`-style fields don't drive the engine anymore.
The estimators write to `cfg.reference_models.<plane>`. If
you carry forward a hand-authored cartridge with legacy fields,
`manifest["models"]["legacy_trim_dicts_present"]` flags it.

**Calling `estimate_p_nucleotide_lengths` on a generic AIRR
table.** External tools usually don't populate `p_*_length`
columns at all. The estimator's 95 % zero-fraction warning is
the signal that you're estimating off P-naive input — at that
point the model is "always-zero", not biology. Either source the
records from a P-aware tool or skip this estimator.

**Expecting the SHM-rate estimator today.** It's deferred to a
follow-up slice. Calling `builder.estimate_shm_rates(...)` raises
`AttributeError`. Use the cartridge's bundled S5F kernel and tune
the rate via the `Experiment.mutate(rate=..., model=...)` knob —
see [SHM and mutation targeting](shm-targeting.md).

**Forgetting `replace=False` when you want single-call
discipline.** All estimators default to `replace=True`: calling
the same estimator twice silently overwrites the first result
and stamps `replaced=True` on the stage entry. That's the
ergonomic default for iterative work. When you want a build
script that *asserts* you only ran each estimator once — for CI
hygiene, or to catch a duplicated method call in a longer chain —
pass `replace=False` on every estimator call. The second call
then raises `ValueError` instead of overwriting.

**Markov NP-base estimation with `pseudocount=0` on tiny
inputs.** If a from-base never appears in the observed sequences,
the engine can't form a transition row for it, and `ValueError`
fires before spec construction. Either use the default
`kind="markov", pseudocount=1.0`, or drop to
`kind="empirical_first_base"` if the from-base coverage isn't
representative.

## Where to go next

- **[Build a reference cartridge](build-reference-cartridge.md)**
  — the broader builder workflow these estimators plug into.
- **[Reference cartridge](../concepts/reference-cartridge.md)** —
  the four-plane conceptual model and why the empirical-models
  plane is typed.
- **[Junction N/P additions](junction-additions.md)** — how the
  NP-base and P-length models affect what the engine samples.
- **[Validation hub](../validation/index.md)** — replay and
  plan-signature guarantees that ride on top of an estimated
  cartridge.
