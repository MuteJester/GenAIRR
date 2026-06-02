# FASTQ / Read Export — Audit + Slice Shipped

**Status: Slice shipped.** `SimulationResult.to_paired_fastq(
r1_path, r2_path, *, quality="illumina", overwrite=False,
**quality_kwargs)` landed alongside 13 implementation tests
([`tests/test_to_paired_fastq.py`](../tests/test_to_paired_fastq.py))
and pin flips in two contract files. Python-only; zero engine
changes; no validator issue kinds, no manifest block, no trace
addresses. The audit body below is preserved for traceability;
sections marked **[Shipped]** describe how the recommendations
actually landed.

The slice adds a single export method to `SimulationResult`
that writes the per-record R1 and R2 sequences as two
synchronized FASTQ files with `@{sequence_id}/1` and
`@{sequence_id}/2` headers, reusing the existing pluggable
quality models (`constant` / `illumina`) from `_qmodel.py`
verbatim.

The audit is **export-layer only**. Nothing in this slice
touches engine semantics, the IR, the trace, validator issue
kinds, manifest blocks, or the simulation passes. FASTQ writing
is pure projection over already-projected AIRR record fields —
the same boundary the existing `to_fasta` / `to_fastq` /
`to_csv` / `to_tsv` exporters live behind.

Companion to
[`tests/test_fastq_export_contract.py`](../tests/test_fastq_export_contract.py)
which freezes today's surfaces (`pin_scaffold_*`) and the gaps
the implementation slice would close (`pin_absence_*`).

**Pre-flight finding (Q3 below): clean yes.** The paired-end
AIRR fields ship populated whenever `.paired_end(...)` ran;
R2 is already reverse-complemented at projection time (the
existing validator enforces this); no per-base quality model
exists on the engine side, so v1 reuses the same constant /
illumina quality models the single-end exporter already
shipped. No engine changes are needed.

---

## Existing scaffolding the contract pins

| Surface | Where | What it pins |
|---|---|---|
| `SimulationResult.to_dataframe` / `.to_csv` / `.to_tsv` | [`src/GenAIRR/result.py`](../src/GenAIRR/result.py) | Tabular export delegating to `csv.DictWriter` via `_write_delimited`. Uses the canonical column order from `_DEFAULT_COLUMN_ORDER`. The new FASTQ exporter does not extend this path — it writes its own format. |
| `SimulationResult.to_fasta(path, *, prefix="seq")` | [`src/GenAIRR/result.py:1066-1076`](../src/GenAIRR/result.py#L1066-L1076) | Single-end FASTA export. Header format: `">{prefix}{i}\|v_call=...\|j_call=..."`. Body is `rec.get("sequence", "")` verbatim. The new paired-FASTQ surface mirrors this naming convention. |
| `SimulationResult.to_fastq(path, *, quality="illumina", prefix="seq", **quality_kwargs)` | [`src/GenAIRR/result.py:1078-1133`](../src/GenAIRR/result.py#L1078-L1133) | Single-end FASTQ export. Header `"@{prefix}{i}\|v_call=...\|j_call=..."`. Body uppercased. Quality string Phred+33 from one of the two pluggable models in `_qmodel.py`. The new paired-FASTQ surface reuses this kwarg vocabulary verbatim. |
| `_qmodel.ConstantQualityModel` / `IlluminaQualityModel` | [`src/GenAIRR/_qmodel.py`](../src/GenAIRR/_qmodel.py) | The two pluggable quality models. `ConstantQualityModel(q=30, n_q=2)` is the canonical "all bases Q30, N bases Q2". `IlluminaQualityModel(peak_q, start_q, end_q, ramp_len, tail_len, low_q, n_q)` is the smoothed trapezoid. Both consumed by `to_fastq` via `resolve_quality_model(name, **kwargs)`. |
| Paired-end AIRR fields | [`engine_rs/src/airr_record/record.rs:269-301`](../engine_rs/src/airr_record/record.rs#L269-L301) | `read_layout: String`, `r1_sequence: String`, `r2_sequence: String`, `r1_start/r1_end/r2_start/r2_end: Option<i64>`, `insert_size: i64`. Always present on every record; defaults are `""` / `None` / `0` when `.paired_end()` didn't run. `read_layout == "paired_end"` is the canonical "paired-end fields are real" signal. |
| R2-orientation invariant | [`docs/paired_end_design.md`](paired_end_design.md) §7 + validator's `PairedEndWindowMismatch { side: R2 }` ([`engine_rs/src/airr_record/validate.rs`](../engine_rs/src/airr_record/validate.rs)) | `r2_sequence` is **already** the reverse complement of `sequence[r2_start:r2_end]` — populated at projection time. The FASTQ writer writes it verbatim; applying a second RC would corrupt the read. |
| `insert_size == r2_end` invariant | [`docs/paired_end_design.md`](paired_end_design.md) §8 | The insert size equals the R2 window's exclusive end position from the 5' molecule end. NOT an overlap / gap count. The FASTQ writer does not consult `insert_size` for record body content; it's metadata. |
| `sequence_id` populated on every record | [`engine_rs/src/airr_record/record.rs:11-12`](../engine_rs/src/airr_record/record.rs#L11-L12) + [`src/GenAIRR/result.py`](../src/GenAIRR/result.py) | Assigned by `from_outcomes` as `f"{id_prefix}{i}"` (default `"seq0"`, `"seq1"`, …). Never empty after a successful run. The new paired-FASTQ writer uses this for read names. |
| `n_quality_errors` is the only quality surface | [`engine_rs/src/passes/corrupt/quality.rs`](../engine_rs/src/passes/corrupt/quality.rs) + [`engine_rs/src/airr_record/record.rs:93`](../engine_rs/src/airr_record/record.rs#L93) | The quality pass records `(error_site[i], error_base[i])` pairs in the trace but **does not** project per-position Phred scores anywhere. v1 constant-quality FASTQ does not collide with any engine surface. |
| Existing absence pin for `to_paired_fastq` | [`tests/test_paired_end_contract.py:125-141`](../tests/test_paired_end_contract.py#L125-L141) | `test_pin_absence_no_single_end_or_to_paired_fastq_aliases` already flags the surface as deferred. The implementation slice flips this pin in lockstep with this audit's pins. |

---

## 1. Q1 — Current output surfaces

The audit confirms five export surfaces today:

```python
result.records              # list[dict[str, Any]] — the live records
result.to_dataframe(*, airr_strict=False)        # pandas DataFrame
result.to_csv(path, *, airr_strict=False)         # csv.DictWriter
result.to_tsv(path, *, airr_strict=False)         # csv.DictWriter, tab
result.to_fasta(path, *, prefix="seq")            # 2-line FASTA
result.to_fastq(path, *, quality="illumina",
                prefix="seq", **quality_kwargs)   # 4-line FASTQ, single-end
```

`to_csv` / `to_tsv` both delegate to `_write_delimited`, which
uses `csv.DictWriter(extrasaction="ignore")` and serialises
`None` as the empty string (not literal `"None"`).

The new surface this audit specifies is `to_paired_fastq`. The
existing five surfaces stay untouched.

### Pinned

- `pin_scaffold_existing_export_surfaces_present` — all five
  shipped exporters exist on `SimulationResult` today.
- `pin_scaffold_to_fastq_uses_constant_or_illumina_models` —
  the existing single-end FASTQ takes `quality="illumina"` /
  `"constant"` only.

---

## 2. Q2 — Shape decision: single-end vs paired-end

The user's recommendation:

> Support both, but paired-end only when `read_layout ==
> "paired_end"` and R1/R2 fields exist.

**Endorsed.** Single-end already shipped. The new slice adds
paired-end **as a separate method** so:

- The single-end signature stays clean (no mode discriminator
  on `to_fastq`).
- The paired-end method's signature takes **two paths**
  (`r1_path`, `r2_path`), making the user intent unambiguous
  at the call site.
- A future "interleaved FASTQ" mode (single file alternating
  R1/R2 records) can land as a new method without rewriting
  either of the two existing ones.

### Recommended API

```python
result.to_paired_fastq(
    r1_path: str,
    r2_path: str,
    *,
    quality: str = "illumina",
    prefix: str = "seq",
    **quality_kwargs,
) -> None
```

Same kwarg vocabulary as the single-end `to_fastq` — the user's
muscle memory carries over. Default `quality="illumina"`
matches the single-end default.

### What about `to_fastq_records()` (record-generator form)?

The user's brief mentioned a "return strings" form as an
alternative to writer methods. **The audit recommends NOT
adding a record-generator surface in v1.**

Reasons:

- The two existing exporters (`to_fasta`, `to_fastq`) both
  write to disk; a record-generator surface would break the
  pattern.
- Tests can read the file back via `path.read_text(...)`. The
  existing single-end FASTQ test suite uses this pattern
  exclusively (`tests/test_experiment.py` lines 1502+).
- A user who wants in-memory FASTQ bytes can `open(path)` with
  a `tempfile.NamedTemporaryFile` or `io.StringIO` if they
  patch `open()`. The added API surface doesn't pay for
  itself.

### Pinned

- `pin_absence_no_to_paired_fastq_method` — carried from
  `tests/test_paired_end_contract.py:125-141`.
- `pin_absence_no_record_generator_fastq_form` — no
  `to_fastq_records()` / `to_paired_fastq_records()` exists.

---

## 3. Q3 — API shape (detailed)

### Method signature

```python
def to_paired_fastq(
    self,
    r1_path: str,
    r2_path: str,
    *,
    quality: str = "illumina",
    prefix: str = "seq",
    **quality_kwargs,
) -> None:
```

### Header format

Two read records per AIRR record. R1 and R2 share the same
`{prefix}{i}` body plus an Illumina-style `/1` and `/2` suffix:

```
R1 file (r1_path):
@seq0/1
<r1_sequence (uppercase)>
+
<phred+33 quality string of length len(r1_sequence)>
@seq1/1
…

R2 file (r2_path):
@seq0/2
<r2_sequence (uppercase, already RC)>
+
<phred+33 quality string of length len(r2_sequence)>
@seq1/2
…
```

The user's brief recommends `@{sequence_id}/1` / `@{sequence_id}/2`
— **endorsed**. Simpler than the full Illumina header
(`@instrument:run:flowcell:lane:tile:x:y 1:N:0:CGATGT`),
unambiguous for pairing, and survives every downstream tool's
read-name parser the audit checked (BWA, STAR, samtools, the
Picard tools).

The header **deliberately omits** the `|v_call=…|j_call=…`
metadata fragment that the single-end exporters include.
Rationale:

- Paired-end FASTQ is consumed by alignment pipelines that
  *expect* short, stable read names. The `|`-pipe convention
  some BWA versions tolerate is not universally portable;
  STAR < 2.7 splits on `|` and corrupts the rest.
- Users who want AIRR metadata round-trip after alignment use
  `to_csv` / `to_tsv` and join on `sequence_id`.
- The `{prefix}{i}` body matches the single-end FASTQ's
  default prefix exactly, so a user mixing single-end + paired
  files in the same downstream pipeline gets consistent
  identifiers.

### Quality string

Same pluggable model surface as the single-end exporter — the
v1 slice reuses `_qmodel.resolve_quality_model(name,
**quality_kwargs)` verbatim, then calls
`model.quality_array(r1_sequence)` for R1 and
`model.quality_array(r2_sequence)` for R2 independently. Each
call gets the read-length-appropriate quality string.

**Important:** the `IlluminaQualityModel`'s trapezoid shape is
parameterised by the read length internally — calling it
separately on R1 (e.g. 150 bp) and R2 (e.g. 150 bp) produces
two separate ramp-up + tail-down shapes, NOT a single shape
truncated at the insert size. This is the correct Illumina-style
behaviour: both reads have their own quality profile starting
from position 0 of the read.

### Pinned

- `pin_absence_no_to_paired_fastq_method` (carry).
- `pin_scaffold_quality_models_take_sequence_only` — both
  `ConstantQualityModel.quality_array(seq)` and
  `IlluminaQualityModel.quality_array(seq)` take the sequence
  string as their only argument; they don't read any AIRR
  field directly. The new paired-FASTQ slice calls each model
  twice per record (R1 + R2) and the model has no per-read
  state.

---

## 4. Q4 — Quality strings

### Today's quality data

The engine emits:

- `n_quality_errors` (integer count of `corrupt.quality` base
  changes per record).
- `corrupt.quality.error_site[i]` (pool position, in the
  trace).
- `corrupt.quality.error_base[i]` (lowercase destination base,
  in the trace).

There is **no per-position Phred score** anywhere — no
`quality_string: String` on `AirrRecord`, no `Vec<u8>` quality
array on the IR, no quality address in the trace. The engine's
quality model is a count + position + base substitution; the
Phred string is invented at FASTQ-write time by the Python
quality model.

### Recommended discipline

The v1 slice **does not** add a per-base quality model.
Reasons:

- Inventing a sequencer-aware quality model is a separate
  biology slice with its own audit — instrument profiles, GC
  ramps, error-rate-by-cycle, etc.
- The existing `to_fastq` already supports two pluggable
  models (constant + illumina) and tests pin both their
  behaviour. v1 reuses them verbatim.
- A user who wants paired-end FASTQ with realistic per-base
  quality can author their own `QualityModel` and pass it via
  `**quality_kwargs` — the resolve surface accepts it. (Or
  more practically — they can post-process the output FASTQ
  through `simseq`, `wgsim`, `art_illumina`, etc.)

### Composition with corruption-stage events

The single-end `to_fastq` routes **lowercase bases** (the
GenAIRR corruption-marker convention) to the `low_q`
parameter. This means a `n_quality_errors > 0` record carries
its quality-degraded positions through to the FASTQ as
lower-Phred bases:

```python
to_fastq(path, quality="constant", q=30, low_q=10)
# → uppercase bases get Q30, lowercase bases get Q10
```

The paired-FASTQ slice inherits this behaviour for free —
both R1 and R2 windows are substrings of the full
(case-preserved) `sequence`, so a lowercase base in the full
sequence stays lowercase in the window and routes to `low_q`
in the FASTQ output.

### Pinned

- `pin_absence_no_per_base_quality_field_on_airr_record` —
  no `quality_string` / `phred_array` / `r1_quality` /
  `r2_quality` field exists on `AirrRecord`.
- `pin_scaffold_quality_pass_records_count_only` — the
  `corrupt.quality` pass records the integer count + per-event
  (site, base) pairs in the trace, but does NOT emit a
  per-position Phred score.
- `pin_scaffold_constant_and_illumina_models_exist` — the
  two pluggable models the v1 slice reuses.

---

## 5. Q5 — Read names

The user's recommendation:

> Use `sequence_id`. For paired-end, append `/1` and `/2`.

**Endorsed exactly as written.** Header format:

```
@{sequence_id}/1
@{sequence_id}/2
```

### Why not full Illumina header

A "real" Illumina header looks like:

```
@HISEQ:1234:H7AT2BCXX:1:1101:1234:5678 1:N:0:CGATGT
```

The seven-field colon-separated tuple is mostly random
metadata that bcl2fastq writes and downstream tools mostly
ignore. The interesting suffix `1:N:0:CGATGT` carries the
read direction (1/2), the "not passing filter" flag, the
control bits, and the index sequence. None of these have a
GenAIRR analogue (no flow cell, no lane, no index — every
record is a single in-silico molecule).

The `/1` / `/2` suffix is the **older, universally-portable**
read-direction marker. Pre-Illumina-Casava-1.8 tools all
understand it; current tools (BWA / STAR / samtools / Picard)
all accept it. Adopting it is a no-future-regret choice.

### Why `sequence_id` rather than `{prefix}{i}`

Two paths to consider:

- **`{prefix}{i}`** (mirrors `to_fastq` / `to_fasta`). Stable
  but doesn't survive a `result.records` reorder.
- **`sequence_id`** (the AIRR record's own identifier).
  Survives any record-list reordering; trivially joinable
  against `result.to_csv(...)`.

**Audit recommends `sequence_id`**. The single-end exporters
should arguably also migrate, but that's a separate
ergonomics change with backwards-compat implications. v1
paired-end gets it right from day one.

### Header format for both reads

Final shape:

```
R1:  @{sequence_id}/1\n
     {r1_sequence.upper()}\n
     +\n
     {q_string_r1}\n
R2:  @{sequence_id}/2\n
     {r2_sequence.upper()}\n
     +\n
     {q_string_r2}\n
```

No `|v_call=…|j_call=…` metadata — see §3 above for
rationale.

### Pinned

- `pin_scaffold_sequence_id_is_string_field_on_airr_record`
  — the field exists and is always populated.

---

## 6. Q6 — Validation

### What v1 rejects

The v1 slice rejects three classes of misconfiguration at
write time:

1. **Records with `read_layout != "paired_end"`** — the
   record doesn't carry R1/R2 windows. Raise `ValueError`
   with a hint pointing at `.paired_end(...)`.
2. **Empty `r1_sequence` or `r2_sequence`** on a paired-layout
   record — defensive; shouldn't happen if the projection
   ran cleanly, but a future regression in the builder would
   produce empty windows that the FASTQ writer must not
   silently emit (a downstream aligner would fail on a
   zero-length read).
3. **Quality-string length mismatch** — `len(q_array_r1) !=
   len(r1_sequence)` indicates a buggy custom quality model.
   Mirror the existing single-end `to_fastq` check at
   [`result.py:1124-1128`](../src/GenAIRR/result.py#L1124-L1128).

### What v1 does NOT do

- **Does NOT call `validate_records(refdata)` automatically.**
  The user opts into AIRR validation as a separate step. FASTQ
  export is a pure projection over the record dict; it has no
  refdata dependency.
- **Does NOT validate the R2-is-RC invariant.** That's the
  AIRR validator's job (`PairedEndWindowMismatch { side: R2
  }`); the FASTQ writer takes the record's `r2_sequence` at
  face value. If a downstream consumer hand-edited the record,
  `validate_records(refdata)` will fire the issue before the
  FASTQ writer would even need to.
- **Does NOT check `insert_size`.** The FASTQ output doesn't
  consume insert_size; that field is record-level metadata
  for downstream aligner pairs.

### Error message conventions

Mirror the existing `to_fastq` error shape:

```python
raise ValueError(
    f"to_paired_fastq: record {i} ({rec.get('sequence_id', '?')!r}) "
    f"has read_layout={rec.get('read_layout', '')!r} — paired-end "
    "FASTQ export requires read_layout='paired_end'. Run "
    "Experiment.paired_end(r1_length=…, insert_size=…) on the "
    "experiment before exporting."
)
```

### Pinned

- `pin_absence_no_paired_fastq_export_validator_kinds` — the
  validator does NOT need new issue kinds (FASTQ export
  doesn't surface validator issues; it raises Python errors
  directly).

---

## 7. Q7 — Interaction with `rev_comp`, `end_loss`, paired-end

### `rev_comp`

The `apply_rev_comp_projection` projection step runs **before**
the paired-end window selection. So when `RevCompPass` fires
`applied=True`:

- `rec["sequence"]` is already the reverse-complemented full
  molecule.
- `rec["r1_sequence"]` is the forward window of that already-
  flipped molecule.
- `rec["r2_sequence"]` is the RC of the 3' window of that
  already-flipped molecule.

The FASTQ writer takes the strings verbatim. No
double-flip / under-flip / interaction-with-rev_comp logic
needed.

### `end_loss`

End-loss truncates the assembled sequence before paired-end
window selection. So R1/R2 windows are always within the
post-end-loss `sequence`. The FASTQ writer doesn't need to
know about end-loss at all.

### Paired-end window math

The user's brief says explicitly: **"FASTQ uses already-projected
`sequence`, `r1_sequence`, `r2_sequence`. It must not recompute
layout."**

Confirmed. The writer reads the three fields directly. It
does NOT consult `r1_start` / `r1_end` / `r2_start` / `r2_end`
/ `insert_size` to slice `sequence` itself. The projection
kernel already did that; double-slicing in the writer would
re-introduce the R2-orientation bug class that the
`PairedEndWindowMismatch` validator was designed to catch.

### Composition with `n_quality_errors`

A record that went through `sequencing_errors()` has lowercase
bases scattered through `r1_sequence` and `r2_sequence`. The
existing `to_fastq` routes lowercase bases to `low_q` (default
10). The paired-FASTQ slice inherits this behaviour because
`model.quality_array(seq)` already runs on case-sensitive
input.

### Pinned

- `pin_scaffold_paired_end_fields_default_when_disabled` —
  on a record without `.paired_end()`, all eight fields hold
  their documented defaults (`read_layout=""`, `r1_sequence=""`,
  …, `insert_size=0`).

---

## 8. Edge cases the implementation slice must handle

| Case | Expected behaviour |
|---|---|
| No `.paired_end()` ran on the experiment | Every record has `read_layout=""` → `to_paired_fastq` raises `ValueError` on the first record. |
| Mixed batch — some records paired-end, some not | First non-paired record raises. v1 does NOT silently skip; a partial output file is worse than a clean error. |
| `quality="constant", q=20, low_q=5` | Both R1 and R2 quality strings emit Q20 for uppercase, Q5 for lowercase, Q2 for `N`. |
| `quality="illumina"` with default kwargs | Both R1 and R2 get the smoothed trapezoid independently. |
| Custom quality model via `**quality_kwargs` | Forwarded to `resolve_quality_model`. |
| `r1_sequence` contains an `N` base (e.g. from `ambiguous_base_calls()`) | The `N` gets the model's `n_q` value (default 2). Same as single-end. |
| Both files written to the same path | Caller error; `open(path, "w")` overwrites the first file. v1 does not detect this — let the OS handle it. |
| Empty result (no records) | Both files written empty (zero bytes). Caller-detectable via filesize / record count. |
| `result` is `None` or has no records | `to_paired_fastq` requires a populated result; this is the same precondition as `to_fasta` / `to_fastq` / `to_csv`. |
| Trace replay then FASTQ export | Replays produce identical AIRR records → identical FASTQ output by construction. No special handling needed. |

---

## 9. Performance

### Per-record cost

Two string writes per record (R1 line + R2 line), one
`model.quality_array(seq)` call per read. Same shape as the
single-end exporter. No new architectural cost.

### File I/O pattern

Two open file handles instead of one. Write paths are
sequential; no buffering / `os.fsync` discipline needed
beyond Python's default. The single-end exporter doesn't
explicitly flush either.

### Pinned

- `pin_scaffold_existing_to_fastq_is_record_serial` — the
  single-end exporter iterates records once and writes
  sequentially; the paired exporter does the same.

---

## 10. Trace / replay / cartridge identity

**Zero impact on any of the three.** FASTQ writing is pure
projection over already-projected AIRR records:

- No new trace addresses.
- Plan signature unchanged.
- `refdata_content_hash` unchanged.
- No new manifest block.

### Pinned

- `pin_scaffold_no_new_trace_addresses_under_fastq_export` —
  source-level pin that no `Export*` / `Fastq*` choice
  addresses exist in `address.rs`.

---

## 11. Manifest extension

The v1 slice does **NOT** add a manifest block. FASTQ export
isn't a cartridge-level capability; every cartridge can export
to FASTQ if `.paired_end()` ran. A `fastq_export_support`
block would carry zero per-cartridge information.

### Pinned

- `pin_absence_no_fastq_export_support_in_manifest` —
  manifest doesn't carry the block.

---

## 12. Implementation order **[Shipped]**

A single self-contained slice. Three sub-steps, all landed:

### Step 1 — `to_paired_fastq` method

- Add to `SimulationResult` in [`src/GenAIRR/result.py`](../src/GenAIRR/result.py)
  immediately after `to_fastq` (line 1133).
- Signature per §3: `(r1_path, r2_path, *, quality="illumina",
  prefix="seq", **quality_kwargs)`.
- Loop over `self._records`. For each record:
  - Validate `read_layout == "paired_end"` (raise on
    mismatch).
  - Pull `r1_sequence`, `r2_sequence` (raise on either
    empty).
  - Compute `q_array_r1 = model.quality_array(r1_sequence)`
    and the R2 equivalent.
  - Assert length parity.
  - Write four lines to each file:
    `@{sequence_id}/1\n` / `{r1_sequence.upper()}\n` / `+\n` /
    `{phred_to_ascii(q_array_r1)}\n` (and `/2` mirror for
    R2).

Cost estimate: ~50 lines.

### Step 2 — Tests

Mirror the existing single-end FASTQ test discipline in
[`tests/test_experiment.py`](../tests/test_experiment.py)
lines 1502+. New file
`tests/test_to_paired_fastq.py` (or inline in
`test_experiment.py` next to the single-end block).

Tests:

- `test_to_paired_fastq_writes_r1_and_r2_in_canonical_format`
  — 4-line-per-record per file, header format, sequence
  uppercase, quality length == sequence length.
- `test_to_paired_fastq_r2_is_already_rc` — the R2 file's body
  equals `rec["r2_sequence"].upper()` verbatim (no double
  RC).
- `test_to_paired_fastq_rejects_non_paired_layout` — clear
  error message when `read_layout != "paired_end"`.
- `test_to_paired_fastq_rejects_empty_r1` /
  `_rejects_empty_r2` — defensive empty-string guards.
- `test_to_paired_fastq_constant_quality_writes_q_string` —
  constant Q30, mixed-case bases route to low_q.
- `test_to_paired_fastq_illumina_quality_shape` — trapezoid
  shape on R1 and R2 independently.
- `test_to_paired_fastq_with_n_bases_get_low_quality` —
  N bases get `n_q` (default 2).
- `test_to_paired_fastq_round_trips_via_pyfastq_parser` (if
  Bio.SeqIO available — opt-in via importorskip).
- `test_to_paired_fastq_empty_result_writes_empty_files` —
  zero records produces zero-byte files.

Cost estimate: ~150 lines.

### Step 3 — Pin flips

- Flip `tests/test_paired_end_contract.py::test_pin_absence_no_single_end_or_to_paired_fastq_aliases`
  → `pin_present_to_paired_fastq_method`.
- Flip the four `pin_absence_*` entries in
  `tests/test_fastq_export_contract.py` (this audit's
  companion).

Cost estimate: ~30 lines.

**Total estimate:** ~230 lines, all Python. Zero Rust
changes.

### Why a single slice

Same argument as the V-subregion counters slice: the export
extension has no internal phase boundary. The method, tests,
and pin flips all land together as one mechanical change to
the result module + one new test file.

---

## 13. Test surface — what this audit pins

Mirrored in
[`tests/test_fastq_export_contract.py`](../tests/test_fastq_export_contract.py).

### `pin_scaffold_*` — pre-existing surfaces the slice builds on

1. `SimulationResult.to_fasta` exists with the documented
   signature.
2. `SimulationResult.to_fastq` exists with `quality` /
   `prefix` / `**quality_kwargs`.
3. `_qmodel.ConstantQualityModel` and `IlluminaQualityModel`
   both expose `quality_array(seq)` returning a Phred-int list
   the same length as the input sequence.
4. `_qmodel.resolve_quality_model("constant" | "illumina",
   **kwargs)` and `phred_to_ascii(q_array)` are the canonical
   helpers the new exporter reuses.
5. The paired-end AIRR fields (`read_layout`, `r1_sequence`,
   `r2_sequence`, `r1_start/end`, `r2_start/end`,
   `insert_size`) are present and documented on every record.
6. `read_layout == "paired_end"` is the canonical "paired-end
   fields are real" sentinel.
7. `r2_sequence` is already reverse-complemented at
   projection time (validator's `PairedEndWindowMismatch` enforces
   it); the new writer outputs it verbatim.
8. `sequence_id` is a string field on every record, never
   empty after a successful run.
9. `n_quality_errors` and `corrupt.quality.error_site[i]` /
   `error_base[i]` are the only quality-related surfaces;
   no per-position Phred field exists.

### `pin_absence_*` — gaps the implementation slice closes

10. No `SimulationResult.to_paired_fastq` method.
11. No `to_fastq_records()` / `to_paired_fastq_records()`
    record-generator surface.
12. No `r1_quality` / `r2_quality` / `quality_string` field
    on `AirrRecord`.
13. No `fastq_export_support` block in the cartridge
    manifest.

### Doc anchor

14. The audit doc exists and references the contract file;
    structure intact.

---

## 14. Out of scope

Documented here so a future implementer doesn't accidentally
expand the work.

- **Interleaved FASTQ** (R1 and R2 alternating in one file).
  Trivial follow-up; not in v1.
- **Compressed FASTQ output** (`.fastq.gz`). Trivial wrapper
  with `gzip.open(...)`; not in v1.
- **Per-base error model FASTQ** — inventing a sequencer-aware
  per-position quality model. Separate biology slice.
- **BCL / FASTA-with-quality** export. Out of scope; users
  who need these can author their own writer.
- **CLI command** for FASTQ export (`genairr to-paired-fastq
  --r1 …`). The MCP audit is a separate slice.
- **Full Illumina-format read header** (`@HISEQ:1234:…`). The
  `/1` / `/2` suffix is the v1 standard; full headers are a
  follow-up.
- **Reading FASTQ back** — `from_fastq(path)`. The engine is
  one-way (in-silico → output); reading external FASTQ for
  benchmarking is a separate "benchmark ingest" slice.
- **Quality model that reads engine-side per-base data**.
  Would require an event-payload change; out of scope.

---

## 15. Summary table

| Concern | Post-slice state | Status |
|---|---|---|
| Single-end FASTQ export | `to_fastq(path, quality, prefix, **q_kwargs)` ships | **Unchanged** |
| Paired-end FASTQ export | `to_paired_fastq(r1_path, r2_path, *, quality="illumina", overwrite=False, **quality_kwargs)` | **Shipped** |
| Quality model surface | `ConstantQualityModel` + `IlluminaQualityModel` via `resolve_quality_model("name", **kwargs)` reused verbatim | **Shipped, no new model** |
| Read names | Paired-end uses `@{sequence_id}/1` and `@{sequence_id}/2` — no metadata pipe (tool portability) | **Shipped** |
| R2 orientation | `r2_sequence` is already RC of `sequence[r2_start:r2_end]` at projection time (enforced by `PairedEndWindowMismatch`); writer outputs verbatim | **Shipped, no second RC** |
| Validation discipline | Refuses non-paired records (`ValueError`), empty R1/R2 (`ValueError`), quality-length mismatch (`RuntimeError`), existing-path-without-overwrite (`FileExistsError`) | **Shipped** |
| AIRR validator interaction | FASTQ writer does NOT call `validate_records(refdata)` — pure projection | **Shipped, no change** |
| Trace addresses | No new addresses | **Shipped, no change** |
| Plan signature | No change | **Shipped, no change** |
| `refdata_content_hash` | No change | **Shipped, no change** |
| Manifest block | None — FASTQ export isn't a cartridge capability | **Shipped, no manifest block** |
| Engine-side per-base quality | None — v1 uses constant / Illumina models only | **Shipped, no engine surface** |
| Interaction with `rev_comp` / `end_loss` / `paired_end` | Writer reads the three projected fields directly; never re-slices | **Shipped, composes correctly** |
| Overwrite guard | `overwrite=False` default refuses to clobber; `overwrite=True` allows replacement; guard fires when ONLY ONE of the two paths exists | **Shipped** |
| Performance | Two string writes + two `quality_array` calls per record | **Negligible — confirmed by test pass times** |
| Pre-flight bugs found | **None.** | — |

The FASTQ export slice shipped clean. The paired-end projection
plumbing was already in place from earlier slices; the quality
model surface was already in place from the single-end FASTQ
exporter; the file-write pattern was already in place from
`to_fasta` / `to_csv` / `to_tsv`. The slice added one method,
13 tests, and three pin flips — ~250 lines of Python, zero
Rust changes, zero validator changes.

Follow-up slices remain explicitly out of scope per §14
(compressed `.fastq.gz`, interleaved single-file output, CLI
command, per-base sequencer error model, BCL / FASTA-with-
quality export, full Illumina-format read headers,
`from_fastq(path)` ingest).
