# SimulationResult

<p class="lead"><code>SimulationResult</code> is the output wrapper
returned by <code>Experiment.run_records(...)</code>. It holds the
list of AIRR record dicts, the underlying engine <code>Outcome</code>
objects (for trace / replay / validation), legacy parent
<code>Outcome</code>s when <code>expand_clones</code> produced fixed-size
families, and lineage trees when <code>clonal_lineage</code> produced BCR
tree output.
Treat it as a list-like view of records plus the typed validators
and export helpers below.</p>

## Common methods

The eight methods you'll reach for in real pipelines:

| Method | Purpose |
|---|---|
| `.validate_records(refdata)` | Per-record AIRR-output correctness gate |
| `.validate_families()` | Clonal family consistency gate (groups by `clone_id`) |
| `.validate_families_with_parents(refdata=None)` | Parent-aware family validator |
| `.to_dataframe(*, airr_strict=False)` | Return a pandas DataFrame with canonical AIRR columns |
| `.to_tsv(path, *, airr_strict=False)` | Write AIRR-style TSV with the canonical header |
| `.to_csv(path, *, airr_strict=False)` | Write CSV (sibling of `to_tsv`) |
| `.to_fasta(path, *, prefix="seq")` | Write assembled sequences as FASTA |
| `.to_fastq(...)` / `.to_paired_fastq(...)` | Write FASTQ; paired-end requires `read_layout="paired_end"` |

`clonal_repertoire` returns ordinary `SimulationResult` records with
`clone_id` and `duplicate_count`. `clonal_lineage` returns
`SimulationResultWithLineages`, a subclass that adds `.lineage_trees`.
Legacy `expand_clones` returns `SimulationResult` with `.parents`.

### FASTQ exports (prose only)

The two FASTQ-emitting methods are documented in
[Export the results](../getting-started/export-results.md) and
[Paired-end reads and FASTQ](../guides/paired-end-fastq.md), and
are intentionally omitted from the generated block below because
their `**quality_kwargs` parameter is untyped (griffe rejects it
under strict mode). Their public signatures:

```python
result.to_fastq(
    path: str,
    *,
    quality: str = "illumina",       # "illumina" (trapezoid) or "constant"
    prefix: str = "seq",
    **quality_kwargs,                # see paired-end-fastq.md
) -> None

result.to_paired_fastq(
    r1_path: str,
    r2_path: str,
    *,
    quality: str = "illumina",
    overwrite: bool = False,
    **quality_kwargs,
) -> None
```

`to_paired_fastq` requires the experiment to have included
`.paired_end(r1_length=..., r2_length=..., insert_size=...)` -
otherwise it raises. It also refuses to overwrite an existing
output file unless you pass `overwrite=True`.

## Class reference

::: GenAIRR.result.SimulationResult
    options:
      show_source: false
      show_root_heading: true
      heading_level: 3
      members:
        - records
        - outcomes
        - parents
        - validate_records
        - validate_families
        - validate_families_with_parents
        - to_dataframe
        - to_tsv
        - to_csv
        - to_fasta
        - from_outcomes
