# Paired-end reads and FASTQ export

<p class="lead">Paired-end is a projection — it adds R1/R2 read
windows over the final assembled molecule without changing what
the engine simulated. One DSL method enables it; eight AIRR
fields carry the geometry; one writer dumps R1.fastq and
R2.fastq files synchronised by <code>sequence_id</code>.</p>

## What paired-end means in GenAIRR

Paired-end is a **read-layout layer**, not a recombination or
biology mechanism. The simulated receptor truth is already
complete by the time `.paired_end(...)` runs — `paired_end`
just decides where to slice R1 and R2 windows over the final
assembled molecule and stamps the windows onto the AIRR record.

Three consequences follow from that framing:

- **There's still one record per molecule.** Paired-end doesn't
  produce two AIRR rows; it produces one row carrying both R1
  and R2 fields. The two reads share `sequence_id`, V/D/J calls,
  junction, productive flag, mutation counters — everything.
- **R2 is already reverse-complemented at projection time.** The
  `r2_sequence` field carries the RC sequence, sliced from the
  3' window of the insert. The FASTQ writer outputs it verbatim
  — no second flip.
- **The molecule is the simulated truth, not the reads.** All the
  validators run against the molecule. Slicing the R1/R2 windows
  is downstream of biology; the molecule's `sequence`, `v_call`,
  `junction`, etc. are the canonical ground truth regardless of
  the read layout.

## A minimal paired-end simulation

```python
import GenAIRR as ga

result = (
    ga.Experiment.on("human_igh")
      .recombine()
      .paired_end(r1_length=150, insert_size=300)
      .run_records(n=10, seed=1)
)

rec = result[0]
print(rec["read_layout"])           # 'paired_end'
print(rec["r1_sequence"][:40])      # first 40 bp of R1
print(rec["r2_sequence"][:40])      # first 40 bp of R2 (already RC)
print(rec["insert_size"])           # 300
print(rec["r1_start"], rec["r1_end"], rec["r2_start"], rec["r2_end"])
```

Three things to know about `.paired_end(...)`:

- **`r1_length` is required**; `r2_length` defaults to `r1_length`.
- **`insert_size` is required**. Each length argument accepts an
  `int` (fixed), a `(low, high)` tuple (uniform), or a
  `[(value, weight), ...]` list (empirical distribution) — the
  same shape every length-or-count argument uses elsewhere in the
  DSL.
- The pass runs late in the pipeline (after end-loss, after strand
  flips) so its R1/R2 windows reflect every prior corruption.

## Read fields in AIRR records

Eight AIRR fields carry the read layout:

```python
rec["read_layout"]      # 'single' by default; 'paired_end' after .paired_end(...)
rec["r1_sequence"]      # forward window over the insert
rec["r2_sequence"]      # reverse-complement of the 3' window
rec["r1_start"]         # 0-based start of R1 over the molecule
rec["r1_end"]           # exclusive end of R1
rec["r2_start"]         # 0-based start of R2's source window (pre-RC)
rec["r2_end"]           # exclusive end of R2's source window (pre-RC)
rec["insert_size"]      # drawn insert size for this molecule
```

When `.paired_end(...)` is NOT in the pipeline, every field is
present on the record dict for schema stability but carries its
default:

| Field | Default (no paired-end) |
|---|---|
| `read_layout` | `'single'` |
| `r1_sequence`, `r2_sequence` | `''` (empty string) |
| `r1_start`, `r1_end`, `r2_start`, `r2_end` | `0` |
| `insert_size` | `0` |

Code that may receive either shape can safely check
`rec["read_layout"] == "paired_end"` before reading the windows.

## Export paired FASTQ

`SimulationResult.to_paired_fastq(...)` writes two synchronised
FASTQ files:

```python
result.to_paired_fastq("reads_R1.fastq", "reads_R2.fastq")
```

The writer produces one R1 entry and one R2 entry per AIRR
record. Headers use the universally-portable suffix convention:

```text
@seq0/1            ← R1.fastq header for record sequence_id=seq0
gaggtgcagctg...
+
EFFGGGHHIIII...

@seq0/2            ← R2.fastq header for same record (note /2)
ccaggcaccggg...
+
HHIIIJJJKKKK...
```

No `|`-pipe metadata in the headers (some aligners — STAR < 2.7
in particular — split on `|` and lose context). The
`sequence_id` joins back to the AIRR TSV row by row.

### Quality models

Two quality models ship; they reuse the same plumbing as the
single-end `to_fastq` writer:

```python
result.to_paired_fastq(
    "R1.fastq", "R2.fastq",
    quality="illumina",        # default — trapezoid shape, per read
)

result.to_paired_fastq(
    "R1.fastq", "R2.fastq",
    quality="constant", q=30,  # uniform Q-score
)
```

`illumina` produces realistic ramp-up + tail-down quality
profiles. Each read is scored independently — R1 and R2 each get
their own ramp from position 0, matching what real sequencer
output looks like. `constant` is the simpler choice for testing.

### Overwrite behavior

By default the writer refuses to clobber existing files:

```python
result.to_paired_fastq("R1.fastq", "R2.fastq")               # raises if either exists
result.to_paired_fastq("R1.fastq", "R2.fastq", overwrite=True)  # replaces both
```

The pre-existence check covers BOTH files. If one of the two
already exists and `overwrite=False`, the writer raises before
touching either — you won't end up with a stale R1 paired to a
fresh R2.

## Interaction with strand and end-loss

`paired_end` runs after the biology and library-prep passes, so a
few interactions are worth knowing:

- **End-loss happens before read layout.** When
  `.end_loss_5prime(...)` and `.end_loss_3prime(...)` are in the
  pipeline, they shorten the molecule first; `paired_end` then
  slices R1/R2 windows over what's left. A 300-bp insert request
  against a 290-bp post-end-loss molecule clips to what's
  available — no error, just a shorter insert.
- **Random strand orientation flips before R1/R2 windows.** When
  `.random_strand_orientation(prob=0.5)` is in the pipeline and
  the record draws a strand flip, the molecule's sequence is
  reverse-complemented and then R1 / R2 are sliced from the
  flipped molecule. The record's `rev_comp` flag tells you which
  way the strand went; `r1_sequence` / `r2_sequence` reflect the
  post-flip orientation.
- **The writer doesn't recompute or re-RC.** Whatever bytes are
  in `r1_sequence` and `r2_sequence` at AIRR-projection time go
  to the FASTQ file verbatim. R2's reverse-complement is applied
  at projection time, exactly once; the writer never applies a
  second flip.

## Clonal workflows

`paired_end` is a descendant-phase pass — R1/R2 windows are
per-read, so each emitted clone member gets its own independent layout
draw. It is fully supported after legacy `expand_clones(...)`, and
accepted after `clonal_repertoire(...)` with one important caveat:
`clonal_repertoire` collapses identical assembled sequences into one
record with `duplicate_count`, and FASTQ export does not expand that
abundance back into multiple read pairs. `clonal_lineage(...)` does
not support paired-end projection yet. If you put `.paired_end(...)`
before a flat clonal fork, the DSL raises at chain time. The right
order for a collapsed clonal-repertoire record surface:

```python
result = (
    ga.Experiment.on("human_igh")
      .recombine()
      .clonal_repertoire(n_clones=50, max_size=100)
      .mutate(model="s5f", rate=0.05)
      .paired_end(r1_length=150, insert_size=300)
      .run_records(seed=1)
)

result.to_paired_fastq("reads_R1.fastq", "reads_R2.fastq")
# Each emitted read has its own R1 / R2 windows, possibly drawn from
# a different insert_size; identical reads may have collapsed into
# duplicate_count before export.
```

See [Clonal simulation overview](clonal-families.md) for the clonal model
chooser and phase rules in full.

## Validation

`validate_records` checks paired-end geometry as part of the
per-record postcondition catalogue. When `.paired_end(...)` is
in the pipeline, every record is checked for:

- `r1_sequence` and `r2_sequence` length match `r1_end - r1_start`
  and `r2_end - r2_start`.
- The R1 window content matches the molecule's forward slice over
  `[r1_start, r1_end)`.
- The R2 window content matches the reverse-complement of the
  molecule's slice over `[r2_start, r2_end)`.
- `insert_size` equals the drawn insert distance.

When `.paired_end(...)` is NOT in the pipeline, the validator
short-circuits the paired-end checks — `read_layout == "single"`
records pass the eight-field defaults without issue.

The FASTQ writer's job is different: it checks the *shape* of
the output (synchronised file lengths, correct header format,
quality length matches sequence length) but not the *biology*
(no V/D/J cross-checking, no junction validation). For biological
validity, run `validate_records` on the `SimulationResult` before
exporting:

```python
report = result.validate_records(refdata)
assert report, report.summary()
result.to_paired_fastq("R1.fastq", "R2.fastq")
```

See the [validation hub](../validation/index.md) for the full
discipline.

## Common mistakes

A handful of issues that show up repeatedly with paired-end.

**Calling `.paired_end()` before a flat clonal fork.** R1/R2 windows
are per-read, so `paired_end` is descendant-phase. Move
`.paired_end(...)` after `clonal_repertoire(...)` or legacy
`expand_clones(...)`. For exact per-copy paired FASTQ depth today,
use legacy `expand_clones`; `clonal_repertoire` is abundance-collapsed
and exposes copy number through `duplicate_count`. `clonal_lineage(...)`
currently rejects `.paired_end(...)` even after the fork; paired-end
lineage output is a future addition.

**Expecting two AIRR rows per molecule.** There's still one row
per record — paired-end is a layered projection, not a record
multiplier. Both R1 and R2 fields live on the same row. If you
want a row per read for downstream joining, derive it yourself
from the DataFrame (`df.melt(...)`-style); GenAIRR doesn't do the
explode for you.

**Reverse-complementing R2 again on the way out.** `r2_sequence`
is already reverse-complemented at AIRR-projection time. Reading
it and applying another `reverse_complement()` would flip it
back to the forward orientation and the FASTQ would no longer
match what a sequencer actually emits. Trust the field as
written; the validator's `PairedEndWindowMismatch { side: R2 }`
check enforces it.

**Expecting FASTQ export without `.paired_end()` in the
pipeline.** `to_paired_fastq` raises on records whose
`read_layout == "single"` (or whose `r1_sequence` / `r2_sequence`
are empty). For single-end FASTQ output, use `to_fastq(...)`
instead.

## Where to go next

- **[Export the results](../getting-started/export-results.md)** —
  every export format (TSV / CSV / FASTA / FASTQ / DataFrame /
  paired FASTQ) in one place.
- **[Validation & reproducibility](../validation/index.md)** — the
  validator that checks paired-end geometry, plus the
  reproducibility model.
- **[Clonal simulation overview](clonal-families.md)** — the clonal
  model chooser and why `paired_end` must come after flat clonal
  forks.
- **[The Experiment builder](experiment-builder.md)** — where
  `paired_end` sits in the full DSL pipeline.
