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

Writes one row per record with GenAIRR's canonical columns
covering sequence, V/D/J calls, junction, productive flag,
mutation counters, paired-end fields (when present), and arbitrary
metadata fields added by `with_metadata(...)`.

`to_tsv()` defaults to `airr_strict=True`: biological coordinate
starts become 1-based, ends are inclusive, booleans are written as
`T` / `F`, and per-segment CIGARs gain leading `S` and `N` runs for
query and reference offsets when nonzero. Missing values are blank.
This conversion does not validate custom metadata against an AIRR
schema.

Python record dicts use 0-based, half-open intervals and local
per-segment CIGARs. Pass `airr_strict=False` to retain that
representation in TSV.

Biological query coordinates refer to the coding-oriented sequence:
`sequence` when `rev_comp=False`, or its reverse complement when
`rev_comp=True`. Paired-read windows (`r1_start/end`, `r2_start/end`)
always use 0-based, half-open coordinates on the observed `sequence`,
including in strict exports.

## CSV

```python
result.to_csv("repertoire.csv")
```

Uses the same columns as TSV, comma-separated, but defaults to
`airr_strict=False`. Pass `airr_strict=True` for the coordinate,
CIGAR, and `T` / `F` boolean conversion described above.

## FASTA

```python
result.to_fasta("sequences.fasta")
```

Writes one record per sequence with the `sequence_id` plus
`v_call` and `j_call` in the header:

```text
>seq0|v_call=IGHVF6-G23*03|j_call=IGHJ1*01
gaggtgcagctggtg...
```

Pass `prefix="myseq"` to change the per-record ID prefix.

### Blinded handoff - strip truth columns

When the next stage is a tool you are scoring, it must not see the
generating calls. The exporters differ in what they carry:

| Exporter | Carries truth? | What ends up in the header / row |
|---|---|---|
| `to_paired_fastq(r1, r2)` | No | `@seq0/1` and `@seq0/2` - `sequence_id` only |
| `to_fasta(path)` | **Yes** | `>seq0|v_call=...|j_call=...` |
| `to_fastq(path)` | **Yes** | `@seq0|v_call=...|j_call=...` |
| `to_tsv` / `to_csv` / `to_dataframe` | **Yes** | `v_call`, `d_call`, `j_call`, every counter, plus `truth_*` under `expose_provenance=True` |

`to_fasta` and `to_fastq` annotate their headers on purpose, which is
convenient when you are the consumer and wrong when the consumer is a
tool under test. For a blinded panel, write the records yourself -
`sequence_id` is all the downstream join needs:

```python
import GenAIRR as ga

result = (
    ga.Experiment.on("human_igh")
      .recombine()
      .productive_only()
      .run_records(n=10_000, seed=42)
)

# Blinded FASTA - sequence_id and bases, nothing else.
with open("panel_blinded.fasta", "w") as fh:
    for rec in result:
        fh.write(f">{rec['sequence_id']}\n{rec['sequence'].upper()}\n")
```

`.upper()` is deliberate: the corruption passes lowercase the bases
they touch, and some aligners read case as a soft-mask. Normalise on
the way out and keep the original `sequence` in your truth table.

For the truth side, two options:

```python
# Option A - write a regular TSV, then drop truth columns
import pandas as pd

result.to_tsv("panel.tsv")
df = pd.read_csv("panel.tsv", sep="\t")
df.drop(
    columns=[c for c in df.columns if c.startswith("truth_")],
    inplace=True,
)
df.to_csv("panel_blinded.tsv", sep="\t", index=False)

# Option B - don't expose provenance in the first place
result = (
    ga.Experiment.on("human_igh")
       .recombine()
       .run_records(n=10_000, seed=42)   # no expose_provenance=True
)
result.to_tsv("panel_blinded.tsv")
```

Option B is the default - `truth_*` columns are only emitted when you
pass `expose_provenance=True` to `run_records(...)`. Option A is the
recipe for *post hoc* blinding when the same run has to produce both
the panel and the answer key.

Either way, note what Option B does not do: it blinds the `truth_*`
aliases, not the record. `v_call`, `d_call` and `j_call` are always the
generating alleles. A panel for a tool you are scoring should be FASTA
or FASTQ written as above, not a TSV with columns removed one at a
time.

### Keeping the truth table beside the panel

The usual shape is two files joined on `sequence_id`:

```python
result = (
    ga.Experiment.on("human_igh")
      .recombine()
      .productive_only()
      .mutate(model="s5f", rate=0.03)
      .run_records(n=10_000, seed=42, expose_provenance=True)
)

# 1. the blinded input the tool sees
with open("panel.fasta", "w") as fh:
    for rec in result:
        fh.write(f">{rec['sequence_id']}\n{rec['sequence'].upper()}\n")

# 2. the answer key you keep
result.to_tsv("panel_truth.tsv")
```

Run the tool on `panel.fasta`, then join its output back on
`sequence_id` - see
[How `sequence_id` joins everything](#how-sequence_id-joins-everything).
Two things to settle before you score:

- **Allele naming.** On `"human_igh"` the truth calls are ASC cluster
  names (`IGHVF6-G23*03`), while an aligner against an IMGT reference
  emits IMGT names (`IGHV4-31*10`). Translate one side using
  `cfg.asc_tables["V"]`, or simulate on an IMGT-named cartridge
  instead:

    ```python
    from GenAIRR.data import HUMAN_IGH_IMGT

    result = (
        ga.Experiment.on(HUMAN_IGH_IMGT)
          .curate_refdata("functional_anchors_only")
          .recombine()
          .run_records(n=10_000, seed=42)
    )
    ```

- **Which truth column.** `v_call` is the allele the engine
  *generated*. `recalled_v_call` is what GenAIRR's own walker infers
  from the mutated sequence, and it may be a comma-separated tie.
  Score against `v_call` to ask "did the tool recover the generating
  allele"; against `recalled_v_call` to ask "did it agree with an
  independent caller on the same evidence". The two diverge exactly
  where mutation has made the generating allele unrecoverable, which
  is the interesting part of the benchmark. See
  [Truth and evidence-based re-calls](first-airr-record.md#truth-and-evidence-based-re-calls).

## FASTQ

```python
result.to_fastq("sequences.fastq")
```

Writes single-end FASTQ with the same `sequence_id` headers and
synthetic quality scores. Two quality models ship:

```python
result.to_fastq("sequences.fastq", quality="illumina")             # default - trapezoid shape
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
(`@seq0/1` and `@seq0/2`) - no `|`-pipe metadata that some aligners
split on. R2 is already reverse-complemented at projection time;
the writer doesn't apply a second flip.

The default `overwrite=False` refuses to clobber existing files;
pass `overwrite=True` to replace. Quality models match the
single-end writer - `quality="illumina"` (default) or
`quality="constant"` with `q=...`.

## pandas DataFrame

```python
df = result.to_dataframe()
```

Returns one row per record with every AIRR field as a column.
Requires `pip install GenAIRR[all]` (or `pip install pandas` on
its own). Defaults to `airr_strict=False`; `airr_strict=True`
converts biological coordinates and CIGAR placement offsets.
DataFrame boolean values remain booleans in either mode.

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

- **[Validate AIRR records](../validation/validate-records.md)** -
  confirm every field is internally consistent before downstream
  analysis.
- **[The Experiment builder](../guides/experiment-builder.md)** -
  the mental model behind the DSL and what each pass actually does.
- **[Productive sampling guide](../guides/productive.md)** - what
  "productive" means biologically and how the constraint composes
  with other passes.
- **[Build your own cartridge](../concepts/reference-cartridge.md)**
  author a custom reference cartridge from FASTA.
