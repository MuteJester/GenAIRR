# Export the results

<p class="lead">A <code>SimulationResult</code> can write itself to
every common AIRR / sequencing format with a single method call.
Most exports are pure-stdlib; only <code>to_dataframe</code> needs
pandas.</p>

## TSV (AIRR-format)

```python
import GenAIRR as ga

result = ga.Experiment.on("human_igh").recombine().run_records(n=1000, seed=42)

result.to_tsv("repertoire.tsv")
```

Writes the canonical AIRR-spec TSV — one row per record, 50+
columns covering sequence, V/D/J calls, junction, productive flag,
mutation counters, paired-end fields (when present), and arbitrary
metadata fields added by `with_metadata(...)`. The column order is
stable across releases.

Pass `airr_strict=True` to convert coordinate fields to the
1-based-inclusive convention the AIRR-C schema mandates (GenAIRR
emits 0-based-exclusive coordinates by default for downstream
consistency with Python slicing).

## CSV

```python
result.to_csv("repertoire.csv")
```

Identical schema to TSV, comma-separated. Accepts the same
`airr_strict=True` flag.

## FASTA

```python
result.to_fasta("sequences.fasta")
```

Writes one record per sequence with the `sequence_id` plus
`v_call` and `j_call` in the header:

```text
>seq0 IGHVF10-G38*04 IGHJ2*01
gaggtgcagctggtg...
```

Pass `prefix="myseq"` to change the per-record ID prefix.

### Blinded handoff — strip truth columns

FASTA is the natural format when the next stage is an aligner
that should NOT see the ground-truth columns — e.g. an aligner
benchmark or a tool comparison where leaking `v_call` /
`truth_v_call` would defeat the point. `to_fasta` already
emits only `sequence_id` + the assembled sequence, so it's
implicitly blinded.

For a blinded TSV, two options:

```python
# Option A — write a regular TSV, then drop truth columns
import pandas as pd

result.to_tsv("panel.tsv")
df = pd.read_csv("panel.tsv", sep="\t")
df.drop(
    columns=[c for c in df.columns if c.startswith("truth_")],
    inplace=True,
)
df.to_csv("panel_blinded.tsv", sep="\t", index=False)

# Option B — don't expose provenance in the first place
result = (
    ga.Experiment.on("human_igh")
       .recombine()
       .run_records(n=10_000, seed=42)   # no expose_provenance=True
)
result.to_tsv("panel_blinded.tsv")
```

Option B is the default — `truth_*` columns are only emitted
when you pass `expose_provenance=True` to `run_records(...)`.
Option A is the recipe for *post hoc* blinding when the same
batch needs to be shared in both forms.

## FASTQ

```python
result.to_fastq("sequences.fastq")
```

Writes single-end FASTQ with the same `sequence_id` headers and
synthetic quality scores. Two quality models ship:

```python
result.to_fastq("sequences.fastq", quality="illumina")             # default — trapezoid shape
result.to_fastq("sequences.fastq", quality="constant", q=30)       # uniform Q-score
```

`illumina` produces realistic ramp-up + tail-down quality profiles
that mirror real reads. `constant` is convenient for testing.

## Paired-end FASTQ

When the pipeline includes `Experiment.paired_end(r1_length=...,
r2_length=..., insert_size=...)`, each record carries R1 / R2
windows. `to_paired_fastq` writes them to two synchronised files:

```python
result = (
    ga.Experiment.on("human_igh")
      .recombine()
      .paired_end(r1_length=150, insert_size=300)
      .run_records(n=1000, seed=1)
)

result.to_paired_fastq("reads_R1.fastq", "reads_R2.fastq")
```

Headers use the universally-portable `/1` and `/2` suffix
(`@seq0/1` and `@seq0/2`) — no `|`-pipe metadata that some aligners
split on. R2 is already reverse-complemented at projection time;
the writer doesn't apply a second flip.

The default `overwrite=False` refuses to clobber existing files;
pass `overwrite=True` to replace. Quality models match the
single-end writer — `quality="illumina"` (default) or
`quality="constant"` with `q=...`.

## pandas DataFrame

```python
df = result.to_dataframe()
```

Returns one row per record with every AIRR field as a column.
Requires `pip install GenAIRR[all]` (or `pip install pandas` on
its own). Accepts the same `airr_strict=True` flag as the TSV /
CSV writers.

```python
df["productive"].mean()                    # productive fraction
df.groupby("v_call")["v_identity"].mean()  # per-allele mean identity
df["n_mutations"].describe()               # mutation-count summary
```

For downstream analyses that fit comfortably in memory, the
DataFrame is the most direct entry point.

## How `sequence_id` joins everything

Every export uses the same `sequence_id` field as the join key:

- **TSV / CSV / DataFrame** → `sequence_id` is one of the columns.
- **FASTA / FASTQ / paired-end FASTQ** → `sequence_id` is the
  header (with `/1` `/2` suffixes for paired-end).

So a downstream analysis that aligns the FASTQ and then joins the
result back to the TSV's per-record metadata uses `sequence_id`
as the foreign key without any extra wiring. Custom metadata
stamped via `with_metadata(experiment_id=..., tissue=...)` rides
on the TSV/DataFrame; the FASTQ files only carry sequence + quality.

---

## Next step

→ Pick the workflow that matches your goal:

- **[Validate AIRR records](../validation/validate-records.md)** —
  confirm every field is internally consistent before downstream
  analysis.
- **[The simulation pipeline concept](../concepts/pipeline.md)** —
  the mental model behind the DSL and what each pass actually does.
- **[Productive sampling guide](../guides/productive.md)** — what
  "productive" means biologically and how the constraint composes
  with other passes.
- **[Build your own cartridge](../concepts/reference-cartridge.md)**
  — author a custom reference cartridge from FASTA.
