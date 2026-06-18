# Inspect a cartridge: manifest and build report

<p class="lead">Once you have a cartridge - built, estimated, and
ready to simulate - two inspection surfaces let you audit it.
The <strong>manifest</strong> describes the cartridge's current
state. The <strong>build report</strong> describes how that state
was produced. Together they let you verify identity, catalogue
coverage, model availability, and reproducibility hashes - before
you run a simulation, before you ship a dataset, and inside CI.</p>

## Why inspect a cartridge

Three moments where this guide pays off:

- **Before simulating.** Confirm the identity is what you think
  it is (`HUMAN_IGH_OGRDB` vs `HUMAN_IGH_EXTENDED` is the kind
  of mismatch that propagates silently into 100 % of your
  output). Confirm the empirical planes you need are actually
  attached.
- **Before sharing.** Stamp the manifest and the build report
  next to a published dataset. Readers reproducing the dataset
  can verify their cartridge is byte-identical to yours via the
  hashes, and can audit which empirical planes were estimated
  vs hand-authored.
- **In CI.** Pin invariants ("V catalogue has ≥ 300 alleles",
  "no errors in the manifest", "schema checksum stable") in
  unit tests so cartridge drift fails the build instead of
  silently shifting downstream output.
- **When reproducing a dataset-matched model.** If the goal is
  to match an external repertoire's empirical distributions,
  the manifest tells you which planes are estimated from records
  vs falling through to defaults - i.e. where the cartridge
  *isn't* matched yet.

## The two inspection surfaces

```python
manifest = cfg.cartridge_manifest()           # current state
report   = cfg.build_report                   # how it was built
# (after .build(): identical to builder.report())
```

| Surface | What it answers |
|---|---|
| `cfg.cartridge_manifest()` | What is in this cartridge **now**? |
| `cfg.build_report` / `builder.report()` | How was this cartridge **produced**? |

A cartridge can mutate after build - you can pop a model plane,
swap a rules block, or hand-author a field. The manifest reflects
those edits; the build report doesn't. The pair tells you both
the present state and the immutable history of how you got here.

## Reading the manifest

`cfg.cartridge_manifest()` returns a single nested dict with
nine top-level keys:

```python
manifest = cfg.cartridge_manifest()
sorted(manifest)
# ['catalogue', 'curation', 'dropped_allele_fields', 'errors',
#  'hashes', 'identity', 'models', 'orphan_dataconfig_fields',
#  'rules', 'schema_version']
```

The remaining sections walk through each.

### `identity`

```python
manifest["identity"]
# {
#   "name": "HUMAN_IGH_OGRDB",
#   "species": "HUMAN",
#   "locus": "IGH",
#   "reference_set": "OGRDB",
#   "source": "ogrdb-2024-Q3"
# }
```

`name`, `species`, and `reference_set` come from the cartridge
itself. `locus` and `source` are populated through the engine
bridge (`dataconfig_to_refdata`) - they reflect what the engine
actually sees, not what the Python-side dataclass holds.

### `catalogue`

```python
manifest["catalogue"]
# {
#   "v_count": 320,
#   "d_count": 30,
#   "j_count": 6,
#   "c_count": 0,
#   "functional_status_counts": {"functional": 285, "ORF": 25, "pseudogene": 10}
# }
```

Counts come from the Python-side allele lists. Use these as the
first sanity check on any custom-built cartridge - a count of
zero on a segment your chain needs is the most common failure
mode after building from FASTA.

### `rules`

```python
manifest["rules"]
# {
#   "has_explicit_rules": True,
#   "allowed_bases": ["A", "C", "G", "T"],
#   "v_anchor": {"expected_aa": "C", "required": True},
#   "j_anchor": {"expected_aa": "W", "required": True}
# }
```

`has_explicit_rules` is `True` when a `ReferenceRulesSpec` is
attached (via `with_rules` or directly on `cfg.reference_rules`).
The `allowed_bases`, `v_anchor`, and `j_anchor` entries are
populated through the engine bridge - so they reflect the
*effective* rules including any cartridge-shipped defaults.

### `models`

The model plane surfaces every empirical block the cartridge
carries:

```python
manifest["models"]
# {
#   "has_reference_models": True,
#   "shm": {...},
#   "allele_usage": {...},
#   "trim_models": {...},
#   "np_length_models": {...},
#   "np_base_models": {...},
#   "p_nucleotide_models": {...},
#   "np_length_keys": ["NP1", "NP2"],
#   "trim_keys": ["v_trim_3", "d_trim_5", "d_trim_3", "j_trim_5"],
#   "legacy_np_lengths_present": False,
#   "legacy_trim_dicts_present": False
# }
```

Each typed block carries:

- `in_plan_signature` - does this plane affect the plan
  signature? (replay safety surface)
- `in_content_hash` - does this plane affect the refdata
  content hash? (cartridge identity surface)
- A model-specific availability marker (`available`,
  `length_keys`, `models`, etc.)
- `legacy_*_present` - does the cartridge carry an older
  pre-typed-plane field for this surface that didn't auto-lift?

The two `legacy_*_present` flags at the top of the block expose
the pre-typed-plane fields (`NP_lengths`, `trim_dicts`). They
don't drive the engine - only the typed
`reference_models.<plane>` paths do - but the flags tell you
whether a stale dict is sitting on the cartridge.

For the SHM block specifically, `shm.v_subregion_support.available`
is the load-bearing flag for whether `v_subregion_rates` on
`.mutate(...)` is allowed (see
[SHM and mutation targeting](shm-targeting.md)).

### `curation`

```python
manifest["curation"]
# {
#   "source_tag": "ogrdb-2024-Q3|curated:functional_status:functional",
#   "policies": ["functional_status:functional"]
# }
```

The `source_tag` is the unmodified bridged identity source.
`policies` is the parsed list of curation policies applied via
`Experiment.curate_refdata(...)` - empty when the cartridge is
uncurated. The encoding is parsed from `|curated:` markers in
the source tag; see
[The Experiment builder](experiment-builder.md) for the curation
surface.

### `hashes`

```python
manifest["hashes"]
# {
#   "data_config_checksum": "sha256:e3b0c44...",
#   "refdata_content_hash": "sha256:6a1f8b2..."
# }
```

- **`data_config_checksum`** is `cfg.compute_checksum()` - the
  sha256 of the pickled `DataConfig` with a few transient fields
  zeroed (build report and `schema_sha256` itself, so the field
  can't include itself in its own hash).
- **`refdata_content_hash`** is the engine-side content hash -
  the same value that gates replay (see [Trace, replay, and
  reproducibility](trace-replay.md)). It hashes the cartridge
  bytes the engine consumes; changing a plane that affects
  refdata flips it.

The two hashes serve different purposes: `data_config_checksum`
identifies the Python-side dataclass (useful for "did anything
about this cartridge change?"); `refdata_content_hash` identifies
the engine-side material (useful for "does this trace replay?").

### `dropped_allele_fields` and `orphan_dataconfig_fields`

```python
manifest["dropped_allele_fields"]
# ["pseudo_allele_flag", "...other audited drops..."]
manifest["orphan_dataconfig_fields"]
# ["NP_lengths", "NP_transitions", "trim_dicts", "..."]
```

These two lists name documented completeness gaps from the
manifest audit - fields the manifest doesn't surface (because
they're not used by the engine) but that an inspecting tool
should know exist. Reading them tells you whether the cartridge
carries pre-typed-plane orphans that didn't auto-lift.

### `errors`

```python
manifest.get("errors")
# Either absent (clean manifest) or a list of strings
```

The `errors` key **is not present** in a clean manifest. It only
appears when the manifest builder encountered a recoverable
problem reading the bridge or a plane. Treat `manifest.get("errors")`
truthy as a CI gate worth failing on - it means part of the
manifest is incomplete and downstream consumers may be reading
stale data.

## Reading the build report

`builder.report()` (or `cfg.build_report` after `.build()`) returns
a `CartridgeBuildReport` dataclass:

```python
report = cfg.build_report
report.stages                   # list[dict] - one entry per builder stage
report.warnings                 # list[str]  - finalisation warnings
report.rejected                 # list[dict] - per-row + per-allele drops
report.manifest_snapshot        # dict | None - manifest at build time
report.checksum_at_build_time   # str  | None - schema_sha256 stamped on cfg
report.to_dict()                # JSON-clean dict for CI artifacts
```

### `stages`

One entry per builder method that ran. Each carries:

- `stage` - the method name (e.g. `"from_fasta"`,
  `"estimate_allele_usage"`)
- `inputs` - the kwargs the user passed
- `inferred` - what the stage computed (counts, model summaries,
  fallback flags)
- `warnings` - stage-level warnings (e.g. P-naive zero-fraction)
- `replaced` - `True` when an estimator overwrote a prior typed
  plane (idempotency marker; see
  [Estimate cartridge models from real data](estimate-cartridge-models.md))

### `warnings`

Build-finalisation warnings - distinct from the per-stage warnings
inside `stages[*].warnings`. These are emitted by `.build()` itself
during plane normalisation / identity resolution.

### `rejected`

The per-row + per-allele drop list. Each entry has `stage`,
`row_index` (or equivalent), and `reason`. Reasons documented
explicitly in
[Estimate cartridge models from real data](estimate-cartridge-models.md#ambiguity-and-skipped-rows).

### `manifest_snapshot` and `checksum_at_build_time`

These two stamp the cartridge's state **at the moment `.build()`
returned**. If you later mutate the cartridge (rare - usually
discouraged), the live manifest will drift but the snapshot stays
pinned to the build-time state.

### JSON export

```python
import json
report_dict = report.to_dict()
json.dumps(report_dict, indent=2)
```

`to_dict()` returns a JSON-clean shape (no numpy types, no
non-stringable keys). Stamp it into a CI artifact next to the
built cartridge.

## Using it in CI

The canonical pattern:

```python
manifest = cfg.cartridge_manifest()
report   = cfg.build_report

# Identity assertions
assert manifest["identity"]["species"] == "HUMAN"
assert manifest["identity"]["locus"] == "IGH"
assert manifest["identity"]["reference_set"] == "OGRDB"

# Catalogue thresholds
assert manifest["catalogue"]["v_count"] >= 300
assert manifest["catalogue"]["d_count"] >= 25
assert manifest["catalogue"]["j_count"] >= 5

# Models - every plane you need must be attached
assert manifest["models"]["allele_usage"]["available"]
assert manifest["models"]["trim_models"]["available"]
assert manifest["models"]["np_length_models"]["available"]
assert manifest["models"]["np_base_models"]["models"]

# No deferred errors
assert not manifest.get("errors")

# No rejected rows above an acceptable threshold (records-driven build)
assert len(report.rejected) < 0.01 * len(input_records)

# Pin the content hash if you ship traces with this cartridge
assert manifest["hashes"]["refdata_content_hash"] == "sha256:6a1f8b2..."
```

Use the assertions you care about; ignore the ones that drift
legitimately between releases.

## Common warning patterns

Five warnings that surface repeatedly in practice.

**Missing gapped V sequences.** A `with_models` or
`infer_v_subregions` call ran against a cartridge whose V
catalogue lacks IMGT gapped sequences. `v_subregion_rates` won't
be available on `.mutate(...)`; flag-check
`manifest["models"]["shm"]["v_subregion_support"]["available"]`.

**Curatable anchors.** The cartridge ships with anchor positions
that fail strict validation but can be auto-curated by the engine.
Surfaces as a warning under `.allow_curatable_refdata()` or
`.curate_refdata(policy)`; refusing the curation surfaces as a
build-time error instead.

**P-naive input warning.** From `estimate_p_nucleotide_lengths`:
≥ 95 % of input rows reported zero on a `p_*_length` field. The
estimated distribution is "always-zero" - not biology. See
[Estimate cartridge models from real data](estimate-cartridge-models.md#p-nucleotide-naive-input-detector).

**Unknown alleles in estimator records.** Per-row entries in
`report.rejected` with `reason: "unknown_allele"` and a
`segment` / `allele_name` pair. The records reference allele
names that aren't in the FASTA catalogue. Either the cartridge
should be rebuilt with a broader catalogue or the records should
be filtered to known alleles.

**Zero coverage for a model plane.** An estimator stage that
ran on zero rows (all rows rejected). The cartridge silently
falls back to defaults for that plane; the build report's stage
entry stamps `records_used: 0`. Pin the count in CI to catch
this.

## Hashes and reproducibility

Three reproducibility surfaces live in different places:

| Hash | Where it lives | What it identifies |
|---|---|---|
| `cfg.compute_checksum()` | DataConfig method | Python-side dataclass identity |
| `refdata_content_hash` | manifest.hashes / TraceFile | Engine-side material identity |
| `pass_plan_signature` | TraceFile only | Resolved pipeline-stage signature |

**`cfg.compute_checksum()`** answers "did this `DataConfig`
change?" It hashes the pickle bytes with build-report and
schema-checksum fields zeroed. Useful for caching cartridge
loads and detecting accidental mutation.

**`refdata_content_hash`** answers "does this trace replay against
this cartridge?" It hashes the engine-consumed material - the
same hash that gates `replay_from_trace_file` and
`rerun_from_trace_file`. Use this for trace pinning. See
[Trace, replay, and reproducibility](trace-replay.md).

**`pass_plan_signature`** lives only on the trace file, not on
the cartridge. It hashes the *pipeline* - every method call on
the `Experiment`. Changes to the pipeline flip the plan
signature; changes to the cartridge content flip the refdata
content hash. The two are independent.

What none of these prove on their own:

- A matching `data_config_checksum` doesn't prove engine-side
  equality (the engine may translate the dataclass through a
  curation pass that drops fields).
- A matching `refdata_content_hash` doesn't prove the Python
  dataclass is bit-identical (curation can produce identical
  engine material from different dataclass states).
- A matching `pass_plan_signature` doesn't prove the cartridge
  is the same (different cartridges can pin the same plan).

The triple - `(plan_signature, refdata_signature,
refdata_content_hash)` - is what replay actually checks.

## Common mistakes

A handful of issues that show up repeatedly with the inspection
surface.

**Treating the manifest as build history.** The manifest
describes the cartridge's **current** state. If you pop a plane
or hand-edit a field, the manifest reflects the edit and there's
no record on the manifest itself that the cartridge was ever
different. Read the build report alongside it whenever the
question is "how did this cartridge get this way?"

**Treating the build report as current state after mutating
`cfg`.** The build report is **frozen at `.build()`**. If you
mutate `cfg.reference_models` afterwards, the manifest will
change but `cfg.build_report.manifest_snapshot` will still
reflect the pre-mutation state. The snapshot is the historical
record, not a live view.

**Ignoring `manifest.get("errors")`.** The `errors` key is
absent by default and only appears when the manifest builder
encountered a recoverable problem. A truthy `errors` list means
part of the manifest is incomplete - downstream consumers
reading the affected blocks may be reading stale data. Fail CI
on it.

**Expecting legacy fields to auto-lift.** `NP_lengths`,
`NP_transitions`, `NP_first_bases`, `p_nucleotide_length_probs`,
and the per-segment `trim_dicts` all exist on legacy cartridges
but **don't** drive the modern engine. The
`legacy_*_present` flags in the manifest tell you whether
they're sitting on the cartridge as orphans; auto-lift is
deliberately disabled because it would silently change output
bytes. Estimate the typed planes explicitly instead - see
[Estimate cartridge models from real data](estimate-cartridge-models.md).

**Pinning `data_config_checksum` instead of
`refdata_content_hash`.** If your downstream test verifies trace
replay, pin `refdata_content_hash`. If it verifies "did anybody
touch this dataclass", pin `compute_checksum()`. The two answer
different questions; pick the one that matches what your test
actually needs.

## Where to go next

- **[Build a reference cartridge](build-reference-cartridge.md)**
  the builder workflow these inspection surfaces describe.
- **[Estimate cartridge models from real data](estimate-cartridge-models.md)**
  the estimator workflow that drives the `rejected` /
  `warnings` entries in the build report.
- **[Reference cartridge](../concepts/reference-cartridge.md)**
  the four-plane conceptual model behind the `models` block.
- **[Trace, replay, and reproducibility](trace-replay.md)** - how
  `refdata_content_hash` gates trace replay.
